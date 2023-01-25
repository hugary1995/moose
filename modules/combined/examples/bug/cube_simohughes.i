T = 300

[Mesh]
  [Mesh]
    type = GeneratedMeshGenerator
    dim = 3
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  stabilize_strain = true
  large_kinematics = true
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
  [T]
    initial_condition = ${T}
  []
  [stress_yy]
    order = CONSTANT
    family = MONOMIAL
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = cauchy_stress
      index_i = 1
      index_j = 1
      execute_on = 'INITIAL TIMESTEP_END'
    []
  []
  [strain_yy]
    order = CONSTANT
    family = MONOMIAL
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = total_strain
      index_i = 1
      index_j = 1
      execute_on = 'INITIAL TIMESTEP_END'
    []
  []
[]

[Kernels]
  [sdx]
    type = TotalLagrangianStressDivergence
    variable = disp_x
    component = 0
    temperature = T
  []
  [sdy]
    type = TotalLagrangianStressDivergence
    variable = disp_y
    component = 1
    temperature = T
  []
  [sdz]
    type = TotalLagrangianStressDivergence
    variable = disp_z
    component = 2
    temperature = T
  []
[]

[Materials]
  [MMPDS_Elasticity_tensor]
    # Compute the varying elasticity tensor.
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 77223.9
    poissons_ratio = 0.31
  []
  [Thermal_Strain]
    # Thermal expansion.
    type = ComputeThermalExpansionEigenstrain
    thermal_expansion_coeff = 9.81e-6
    temperature = T
    eigenstrain_name = eigenstrain
    stress_free_temperature = 1100
  []
  [Strain]
    type = ComputeLagrangianStrain
    eigenstrain_names = "eigenstrain"
  []
  [flow_stress]
    type = DerivativeParsedMaterial
    f_name = flow_stress
    function = '((-5.598*10^-4)*T+1.17)*(869.6294291160563+(272.9450547421583*((effective_plastic_strain+0.0001)^0.21157910828593482)))'
    material_property_names = 'effective_plastic_strain'
    additional_derivative_symbols = 'effective_plastic_strain'
    args = 'T'
    derivative_order = 2
    compute = false
  []
  [plasticStress]
    type = ComputeSimoHughesJ2PlasticityStress
    flow_stress_material = flow_stress
  []
[]

[BCs]
  [xfix]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0
  []
  [yfix]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
  []
  [zfix]
    type = DirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0
  []
  [ydisp]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 'top'
    function = 't'
    preset = false
  []
[]

[Executioner]
  type = Transient

  dt = 1e-4
  end_time = 0.1
  automatic_scaling = true

  solve_type = 'NEWTON'

  nl_rel_tol = 1e-08
  nl_abs_tol = 5e-13

  petsc_options_iname = '-pc_type '
  petsc_options_value = 'lu'
  line_search = NONE

  [Predictor]
    type = SimplePredictor
    scale = 1
  []
[]

[Postprocessors]
  [stress]
    type = ElementAverageValue
    variable = stress_yy
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [strain]
    type = ElementAverageValue
    variable = strain_yy
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Outputs]
  file_base = 'simohughes_T_${T}'
  csv = true
[]
