[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[XFEM]
  output_cut_plane = true
[]

[UserObjects]
  [cut]
    type = LevelSetCutUserObject
    level_set_var = phi
    heal_always = true
  []
[]

[Mesh]
  [generated_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 40
    ny = 40
  []
[]

[AuxVariables]
  [phi]
  []
[]

[AuxKernels]
  [phi]
    type = FunctionAux
    variable = phi
    function = 'if(t<0.3,x-0.713,x-0.713+(t-0.3))'
  []
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
[]

[AuxVariables]
  [maximum_principal_strain]
    order = CONSTANT
    family = MONOMIAL
  []
  [maximum_principal_stress]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Kernels]
  [TensorMechanics]
    displacements = 'disp_x disp_y'
  []
[]

[AuxKernels]
  [maximum_principal_strain]
    type = RankTwoScalarAux
    variable = maximum_principal_strain
    rank_two_tensor = elastic_strain
    scalar_type = MaxPrincipal
  []
  [maximum_principal_stress]
    type = RankTwoScalarAux
    variable = maximum_principal_stress
    rank_two_tensor = stress
    scalar_type = MaxPrincipal
  []
[]

[Constraints]
  [dispx_constraint]
    type = XFEMSingleVariableConstraint
    variable = disp_x
    geometric_cut_userobject = 'cut'
    alpha = 1e10
    use_penalty = true
    use_displaced_mesh = false
  []
  [dispy_constraint]
    type = XFEMSingleVariableConstraint
    variable = disp_y
    geometric_cut_userobject = 'cut'
    alpha = 1e10
    use_penalty = true
    use_displaced_mesh = false
  []
[]

[BCs]
  [disp_x]
    type = FunctionDirichletBC
    variable = disp_x
    function = 'if(t<0.3,t,0.3)'
    boundary = right
  []
  [disp_x_fix]
    type = DirichletBC
    variable = disp_x
    value = 0
    boundary = left
  []
  [disp_y_fix]
    type = DirichletBC
    variable = disp_y
    value = 0
    boundary = 'left right'
  []
[]

[Materials]
  [elasticity_tensor_A]
    type = ComputeIsotropicElasticityTensor
    base_name = A
    youngs_modulus = 2.1e5
    poissons_ratio = 0.3
  []
  [strain_A]
    type = ComputeFiniteStrain
    base_name = A
  []
  [creep_A]
    type = PowerLawCreepStressUpdate
    base_name = A
    activation_energy = 0
    coefficient = 1e-24
    n_exponent = 4
  []
  [stress_A]
    type = ComputeMultipleInelasticStress
    base_name = A
    inelastic_models = 'creep_A'
  []
  [elasticity_tensor_B]
    type = ComputeIsotropicElasticityTensor
    base_name = B
    youngs_modulus = 2.1e5
    poissons_ratio = 0.24
  []
  [strain_B]
    type = ComputeFiniteStrain
    base_name = B
  []
  [creep_B]
    type = PowerLawCreepStressUpdate
    base_name = B
    activation_energy = 0
    coefficient = 3e-24
    n_exponent = 4
  []
  [stress_B]
    type = ComputeMultipleInelasticStress
    base_name = B
    inelastic_models = 'creep_B'
  []
  [switching_creep]
    type = GeometricCutSwitchingMaterialReal
    geometric_cut_userobject = cut
    prop_name = effective_creep_strain
    base_name_keys = '0 1'
    base_name_vals = 'A B'
    outputs = 'exodus'
  []
  [switching_strain]
    type = GeometricCutSwitchingMaterialRankTwoTensor
    geometric_cut_userobject = cut
    prop_name = elastic_strain
    base_name_keys = '0 1'
    base_name_vals = 'A B'
  []
  [switching_stress]
    type = GeometricCutSwitchingMaterialRankTwoTensor
    geometric_cut_userobject = cut
    prop_name = stress
    base_name_keys = '0 1'
    base_name_vals = 'A B'
  []
  [switching_dstress_dstrain]
    type = GeometricCutSwitchingMaterialRankFourTensor
    geometric_cut_userobject = cut
    prop_name = Jacobian_mult
    base_name_keys = '0 1'
    base_name_vals = 'A B'
  []
[]

[Executioner]
  type = Transient

  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'
  automatic_scaling = true

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10

  dt = 0.01
  end_time = 0.6

  max_xfem_update = 1
[]

[Outputs]
  print_linear_converged_reason = false
  print_nonlinear_converged_reason = false
  print_linear_residuals = false
  exodus = true
[]
