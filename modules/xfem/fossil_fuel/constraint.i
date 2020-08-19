# This test is for two layer materials with different youngs modulus
# and different eigenstrains
# The global stress is determined by switching the stress based on level set values
# The material interface is marked by a level set function
# The two layer materials are glued together

T_ref = 753
CTE_oxide = 1e-2
E_metal = 1.93e5
E_oxide = 2.72e5
nu_metal = 0.3
nu_oxide = 0.3

[GlobalParams]
  displacements = 'disp_r disp_z'
  # decomposition_method = EigenSolution
[]

[Problem]
  type = FEProblem
  coord_type = RZ
[]

[Mesh]
  [generated_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    elem_type = QUAD4
    nx = 44
    ny = 40
    xmin = 13.8
    xmax = 25
    ymin = 0
    ymax = 10
  []
  [oxide]
    type = SubdomainBoundingBoxGenerator
    input = 'generated_mesh'
    block_id = 1
    bottom_left = '13.8 0 0'
    top_right = '15 10 0'
  []
  [metal]
    type = SubdomainBoundingBoxGenerator
    input = 'oxide'
    block_id = 2
    bottom_left = '15 0 0'
    top_right = '25 10 0'
  []
[]

[Variables]
  [disp_r]
  []
  [disp_z]
  []
[]

[AuxVariables]
  [temperature]
    [InitialCondition]
      type = ConstantIC
      value = ${T_ref}
    []
  []
  [stress]
    order = CONSTANT
    family = MONOMIAL
  []
  [pid]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [temperature]
    type = FunctionAux
    variable = temperature
    function = 'if(t<10,${T_ref}+t,${T_ref}+10)'
  []
  [stress]
    type = ADRankTwoScalarAux
    variable = 'stress'
    rank_two_tensor = stress
    scalar_type = 'VonMisesStress'
  []
  [pid]
    type = ProcessorIDAux
    variable = 'pid'
  []
[]

[Kernels]
  [solid_r]
    type = ADStressDivergenceRZTensors
    variable = 'disp_r'
    component = 0
    use_displaced_mesh = true
  []
  [solid_z]
    type = ADStressDivergenceRZTensors
    variable = 'disp_z'
    component = 1
    use_displaced_mesh = true
  []
[]

[Constraints]
  [disp_z_fixed]
    type = EqualValueBoundaryConstraint
    variable = 'disp_z'
    secondary = top
    formulation = penalty
    penalty = 1e10
  []
[]

[BCs]
  [fixed_z]
    type = DirichletBC
    variable = 'disp_z'
    boundary = 'bottom'
    value = 0.0
  []
[]

[Materials]
  # oxide
  [elasticity_tensor_oxide]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = ${E_oxide}
    poissons_ratio = ${nu_oxide}
    block = 1
  []
  [eigenstrain_oxide]
    type = ADComputeThermalExpansionEigenstrain
    stress_free_temperature = ${T_ref}
    temperature = temperature
    thermal_expansion_coeff = ${CTE_oxide}
    eigenstrain_name = 'eigenstrain_oxide'
    block = 1
  []
  [strain_oxide]
    type = ADComputeAxisymmetricRZFiniteStrain
    eigenstrain_names = 'eigenstrain_oxide'
    block = 1
  []
  [stress_oxide]
    type = ADComputeFiniteStrainElasticStress
    block = 1
  []
  # metal
  [elasticity_tensor_metal]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = ${E_metal}
    poissons_ratio = ${nu_metal}
    block = 2
  []
  [strain_metal]
    type = ADComputeAxisymmetricRZFiniteStrain
    block = 2
  []
  [stress_metal]
    type = ADComputeFiniteStrainElasticStress
    block = 2
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  # petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter'
  # petsc_options_value = '201                hypre    boomeramg      8'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist'
  automatic_scaling = true
  # line_search = 'none'

  # controls for nonlinear iterations
  nl_rel_tol = 1e-06
  nl_abs_tol = 1e-08

  # time control
  dt = 1
  end_time = 10

  max_xfem_update = 1
[]

[Outputs]
  print_linear_converged_reason = false
  print_nonlinear_converged_reason = false
  print_linear_residuals = false
  [exodus]
    type = Exodus
    file_base = 'output/constraint'
  []
[]
