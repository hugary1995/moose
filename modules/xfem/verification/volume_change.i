T_ref = 753
CTE_oxide = 1e-2
E_oxide = 2.72e5
nu_oxide = 0.3

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [generated_mesh]
    type = GeneratedMeshGenerator
    dim = 3
  []
  [fixed]
    type = BoundingBoxNodeSetGenerator
    input = generated_mesh
    new_boundary = pin
    bottom_left = '0 0 0'
    top_right = '0 0 0'
  []
[]

[Variables]
  [disp_x]
  []
  [disp_y]
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
  [solid_x]
    type = ADStressDivergenceTensors
    variable = 'disp_x'
    component = 0
    use_displaced_mesh = true
  []
  [solid_y]
    type = ADStressDivergenceTensors
    variable = 'disp_y'
    component = 1
    use_displaced_mesh = true
  []
  [solid_z]
    type = ADStressDivergenceTensors
    variable = 'disp_z'
    component = 2
    use_displaced_mesh = true
  []
[]

[BCs]
  [fixed_x]
    type = DirichletBC
    variable = 'disp_x'
    boundary = 'pin'
    value = 0.0
  []
  [fixed_y]
    type = DirichletBC
    variable = 'disp_y'
    boundary = 'pin'
    value = 0.0
  []
  [fixed_z]
    type = DirichletBC
    variable = 'disp_z'
    boundary = 'pin'
    value = 0.0
  []
[]

[Materials]
  # oxide
  [elasticity_tensor_oxide]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = ${E_oxide}
    poissons_ratio = ${nu_oxide}
  []
  [eigenstrain_oxide]
    type = ADComputeThermalExpansionEigenstrain
    stress_free_temperature = ${T_ref}
    temperature = temperature
    thermal_expansion_coeff = ${CTE_oxide}
    eigenstrain_name = 'eigenstrain_oxide'
  []
  [strain_oxide]
    type = ADComputeFiniteStrain
    eigenstrain_names = 'eigenstrain_oxide'
  []
  [stress_oxide]
    type = ADComputeFiniteStrainElasticStress
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
  [volume]
    type = VolumePostprocessor
    use_displaced_mesh = true
  []
[]

[Outputs]
  print_linear_converged_reason = false
  print_nonlinear_converged_reason = false
  print_linear_residuals = false
  csv = true
[]
