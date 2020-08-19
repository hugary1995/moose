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
nx = 2

[GlobalParams]
  displacements = 'disp_r disp_z'
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
    nx = ${nx}
    ny = 1
    xmin = 1
    xmax = 3
    ymin = 0
    ymax = 1
  []
  [oxide]
    type = SubdomainBoundingBoxGenerator
    input = generated_mesh
    block_id = 1
    block_name = oxide
    bottom_left = '1 0 0'
    top_right = '2 1 0'
  []
  [metal]
    type = SubdomainBoundingBoxGenerator
    input = oxide
    block_id = 2
    block_name = metal
    bottom_left = '2 0 0'
    top_right = '3 1 0'
  []
  [fixed]
    type = BoundingBoxNodeSetGenerator
    input = metal
    new_boundary = pin
    bottom_left = '3 0 0'
    top_right = '3 0 0'
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
[]

[AuxKernels]
  [temperature]
    type = FunctionAux
    variable = temperature
    function = '${T_ref}+t'
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

[BCs]
  [fixed_r]
    type = DirichletBC
    variable = 'disp_r'
    boundary = 'right'
    value = 0.0
  []
  [fixed_z]
    type = DirichletBC
    variable = 'disp_z'
    boundary = 'top bottom'
    value = 0.0
  []
[]

[Materials]
  # oxide
  [elasticity_tensor_oxide]
    type = ADComputeIsotropicElasticityTensor
    block = oxide
    youngs_modulus = ${E_oxide}
    poissons_ratio = ${nu_oxide}
  []
  [eigenstrain_oxide]
    type = ADComputeThermalExpansionEigenstrain
    block = oxide
    stress_free_temperature = ${T_ref}
    temperature = temperature
    thermal_expansion_coeff = ${CTE_oxide}
    eigenstrain_name = 'eigenstrain_oxide'
  []
  [strain_oxide]
    type = ADComputeAxisymmetricRZFiniteStrain
    block = oxide
    eigenstrain_names = 'eigenstrain_oxide'
  []
  [stress_oxide]
    type = ADComputeFiniteStrainElasticStress
    block = oxide
  []
  # metal
  [elasticity_tensor_metal]
    type = ADComputeIsotropicElasticityTensor
    block = metal
    youngs_modulus = ${E_metal}
    poissons_ratio = ${nu_metal}
  []
  [strain_metal]
    type = ADComputeAxisymmetricRZFiniteStrain
    block = metal
  []
  [stress_metal]
    type = ADComputeFiniteStrainElasticStress
    block = metal
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist'
  automatic_scaling = true
  # line_search = 'none'

  # controls for nonlinear iterations
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12

  # time control
  # dt = 8640.0
  dt = 1
  end_time = 30
  # end_time = 5184000 #60 days
  # end_time = 7776000 #90 days
  # end_time = 15552000 #180 days
[]

[Postprocessors]
  [disp_r]
    type = PointValue
    variable = 'disp_r'
    point = '1 0 0'
  []
[]

[Outputs]
  print_linear_converged_reason = false
  print_nonlinear_converged_reason = false
  print_linear_residuals = false
  [exodus]
    type = Exodus
    file_base = 'output/visualize_rz'
  []
[]
