[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  [gmg]
    type = AnnularMeshGenerator
    dmin = 0
    dmax = 90
    nr = ${nr}
    nt = ${nt}
    rmin = 1
    rmax = 4
  []
  [block1]
    type = ParsedSubdomainMeshGenerator
    input = 'gmg'
    block_id = 1
    combinatorial_geometry = 'r:=sqrt(x^2+y^2);r<2'
  []
  [block2]
    type = ParsedSubdomainMeshGenerator
    input = 'block1'
    block_id = 2
    combinatorial_geometry = 'r:=sqrt(x^2+y^2);r>=2'
  []
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
[]

[AuxVariables]
  [stress_rr]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_tt]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [stress_rr]
    type = ADRankTwoAux
    variable = 'stress_rr'
    boundary = dmin
    rank_two_tensor = 'stress'
    index_i = 0
    index_j = 0
  []
  [stress_tt]
    type = ADRankTwoAux
    variable = 'stress_tt'
    boundary = dmin
    rank_two_tensor = 'stress'
    index_i = 1
    index_j = 1
  []
[]

[Kernels]
  [solid_x]
    type = ADStressDivergenceTensors
    variable = 'disp_x'
    component = 0
    use_displaced_mesh = false
  []
  [solid_y]
    type = ADStressDivergenceTensors
    variable = 'disp_y'
    component = 1
    use_displaced_mesh = false
  []
[]

[BCs]
  [Pressure]
    [inner]
      boundary = 'rmin'
      function = '10*t'
    []
    [outer]
      boundary = 'rmax'
      function = '5*t'
    []
  []
  [symmetry_y]
    type = DirichletBC
    variable = 'disp_y'
    value = 0
    boundary = 'dmin'
  []
  [symmetry_x]
    type = DirichletBC
    variable = 'disp_x'
    value = 0
    boundary = 'dmax'
  []
[]

[Materials]
  [elasticity_tensor_1]
    type = ADComputeIsotropicElasticityTensor
    shear_modulus = 1e3
    poissons_ratio = 0.3
    block = 1
  []
  [elasticity_tensor_2]
    type = ADComputeIsotropicElasticityTensor
    shear_modulus = 2e3
    poissons_ratio = 0.3
    block = 2
  []
  [strain]
    type = ADComputeSmallStrain
  []
  [stress]
    type = ADComputeLinearElasticStress
  []
[]

[VectorPostprocessors]
  [u_r]
    type = NodalValueSampler
    variable = 'disp_x'
    boundary = dmin
    sort_by = x
  []
  [stress_rr]
    type = LineValueSampler
    variable = 'stress_rr'
    num_points = 1000
    start_point = '1 0 0'
    end_point = '4 0 0'
    sort_by = x
  []
  [stress_tt]
    type = LineValueSampler
    variable = 'stress_tt'
    num_points = 1000
    start_point = '1 0 0'
    end_point = '4 0 0'
    sort_by = x
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist'
  automatic_scaling = true

  # controls for nonlinear iterations
  nl_rel_tol = 1e-06
  nl_abs_tol = 1e-08

  dt = 1
  end_time = 10
[]

[Outputs]
  print_linear_converged_reason = false
  print_nonlinear_converged_reason = false
  print_linear_residuals = false
  [exodus]
    type = Exodus
    file_base = 'output/visualize_nr_${nr}_nt_${nt}'
    execute_on = 'FINAL'
  []
  [csv]
    type = CSV
    file_base = 'output/data_nr_${nr}_nt_${nt}'
    execute_on = 'FINAL'
  []
[]
