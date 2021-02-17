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
    negative_id = 1
    positive_id = 2
    heal_always = true
  []
  [moving_interface]
    type = ChangeElementSubdomainByGeometricCut
    cut = cut
    origin_subdomain_id = 1
    target_subdomain_id = 2
    moving_boundary_name = 'moving_interface'
    include_domain_boundary = false
    reversible = true
    initialize_solution_on_move = false
    execute_on = 'TIMESTEP_BEGIN TIMESTEP_END'
  []
[]

[Mesh]
  [generated_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 40
    ny = 40
  []
  [left]
    type = SubdomainBoundingBoxGenerator
    input = 'generated_mesh'
    block_id = 1
    bottom_left = '0 0 0'
    top_right = '0.8 1 1'
  []
  [right]
    type = SubdomainBoundingBoxGenerator
    input = 'left'
    block_id = 2
    bottom_left = '0.8 0 0'
    top_right = '1 1 1'
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
    rank_two_tensor = total_strain
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
    block = 1
    youngs_modulus = 2.1e5
    poissons_ratio = 0.3
  []
  [elasticity_tensor_B]
    type = ComputeIsotropicElasticityTensor
    block = 2
    youngs_modulus = 2.1e5
    poissons_ratio = 0.24
  []
  [strain_A]
    type = ComputeFiniteStrain
    block = 1
  []
  [strain_B]
    type = ComputeFiniteStrain
    block = 2
  []
  [creep_A]
    type = PowerLawCreepStressUpdate
    block = 1
    activation_energy = 0
    coefficient = 1e-24
    n_exponent = 4
  []
  [creep_B]
    type = PowerLawCreepStressUpdate
    block = 2
    activation_energy = 0
    coefficient = 3e-24
    n_exponent = 4
  []
  [stress_A]
    type = ComputeMultipleInelasticStress
    block = 1
    inelastic_models = 'creep_A'
  []
  [stress_B]
    type = ComputeMultipleInelasticStress
    block = 2
    inelastic_models = 'creep_B'
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
