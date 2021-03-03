[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[XFEM]
  output_cut_plane = true
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 50
    ny = 50
  []
  [left]
    type = SubdomainBoundingBoxGenerator
    input = 'gen'
    block_id = 0
    bottom_left = '0 0 0'
    top_right = '0.5 1 1'
  []
  [right]
    type = SubdomainBoundingBoxGenerator
    input = 'left'
    block_id = 1
    bottom_left = '0.5 0 0'
    top_right = '1 1 1'
  []
[]

[UserObjects]
  [cut]
    type = LevelSetCutUserObject
    level_set_var = phi
    negative_id = 0
    positive_id = 1
    heal_always = true
  []
  [moving_circle]
    type = GeometricCutElementSubdomainModifier
    cut = cut
    apply_initial_conditions = false
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Functions]
  [moving_circle]
    type = ParsedFunction
    value = 'if(t<0.1,x-0.013,x-0.013-(t-0.1))'
  []
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
[]

[AuxVariables]
  [phi]
  []
  [temp]
  []
  [stress]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Kernels]
  [solid_x]
    type = StressDivergenceTensors
    variable = 'disp_x'
    component = 0
    use_displaced_mesh = true
    displacements = 'disp_x disp_y'
  []
  [solid_y]
    type = StressDivergenceTensors
    variable = 'disp_y'
    component = 1
    use_displaced_mesh = true
    displacements = 'disp_x disp_y'
  []
[]

[AuxKernels]
  [phi]
    type = FunctionAux
    variable = phi
    function = moving_circle
    execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
  []
  [temp]
    type = FunctionAux
    variable = temp
    function = 'if(t<0.1,t,0.1)'
  []
  [stress]
    type = RankTwoScalarAux
    variable = stress
    rank_two_tensor = stress
    scalar_type = VonMisesStress
  []
[]

[BCs]
  [fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = right
    value = 0
  []
  [fix_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'right bottom'
    value = 0
  []
[]

[Constraints]
  [continuity_x]
    type = XFEMSingleVariableConstraint
    variable = disp_x
    geometric_cut_userobject = cut
    use_penalty = true
    alpha = 1e10
  []
  [continuity_y]
    type = XFEMSingleVariableConstraint
    variable = disp_y
    geometric_cut_userobject = cut
    use_penalty = true
    alpha = 1e10
  []
[]

[Materials]
  [C]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 2.1e5
    poissons_ratio = 0.3
  []
  [eigenstrain]
    type = ComputeThermalExpansionEigenstrain
    temperature = temp
    stress_free_temperature = 0
    thermal_expansion_coeff = 1
    eigenstrain_name = e0
    block = 0
  []
  [strain_1]
    type = ComputeFiniteStrain
    eigenstrain_names = e0
    displacements = 'disp_x disp_y'
    block = 0
  []
  [strain_2]
    type = ComputeFiniteStrain
    displacements = 'disp_x disp_y'
    block = 1
  []
  [stress]
    type = ComputeFiniteStrainElasticStress
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  automatic_scaling = true

  nl_rel_tol = 1e-08
  nl_abs_tol = 1e-10
  dt = 0.01
  num_steps = 30

  max_xfem_update = 1
[]

[Outputs]
  print_linear_residuals = false
  exodus = true
[]
