# This test is for two layer materials with different youngs modulus
# and different eigenstrains
# The global stress is determined by switching the stress based on level set values
# The material interface is marked by a level set function
# The two layer materials are glued together

[GlobalParams]
  order = FIRST
  family = LAGRANGE
[]

[XFEM]
  qrule = volfrac
  output_cut_plane = true
[]

[UserObjects]
  [./level_set_cut_uo]
    type = LevelSetCutUserObject
    level_set_var = oxide_interface
    heal_always = true
  [../]
[]

[Mesh]
  displacements = 'disp_x disp_y'
  [./generated_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 100 #600
    ny = 600 #600
    xmin = 0.0
    xmax = 1e-4# originally 1mm thick
    ymin = -2e-4
    ymax =  2e-4 # originally 4mm long
    elem_type = QUAD4
  [../]
  [./fixed]
    type = BoundingBoxNodeSetGenerator
    new_boundary = 100
    bottom_left = '0 0 0'
    top_right = '0 0 0'
    input = generated_mesh
  [../]
  [./pivot]
    type = BoundingBoxNodeSetGenerator
    new_boundary = 101
    bottom_left = '1e-4 0 0'
    top_right = '1e-4 0 0'
    input = fixed
  [../]
[]

[AuxVariables]
  [./oxide_interface]
    order = FIRST
    family = LAGRANGE
  [../]
  [./temperature]
    initial_condition = 801.0
  [../]
[]

[AuxKernels]
  [./ls_function]
    type = FunctionAux
    variable = oxide_interface
    function = estimated_oxide_growth
  [../]
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
[]

[Functions]
  [./estimated_oxide_growth] #represents the oxide layer growth into the metal
    type = ParsedFunction
    value = 'x-0.288997e-2*sqrt(t*4.32243e-14)' #0.288931 * sqrt(t*3.26e-20*exp(0.0176*T)) [cm] at 801K
  [../]
[]

[AuxVariables]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./von_mises_stress]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./oxide_strain_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./oxide_strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./oxide_strain_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./metal_strain_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./metal_strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./metal_strain_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y'
    use_automatic_differentiation = true
    use_displaced_mesh = true
  [../]
[]

[AuxKernels]
  [./stress_xx]
    type = ADRankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
    variable = stress_xx
  [../]
  [./stress_yy]
    type = ADRankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    variable = stress_yy
  [../]
  [./stress_xy]
    type = ADRankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
    variable = stress_xy
  [../]
  [./von_mises_stress]
    type = ADRankTwoScalarAux
    rank_two_tensor = stress
    scalar_type = VonMisesStress
    variable = von_mises_stress
  [../]
  [./oxide_strain_xx]
    type = ADRankTwoAux
    rank_two_tensor = oxide_total_strain
    index_i = 0
    index_j = 0
    variable = oxide_strain_xx
  [../]
  [./oxide_strain_yy]
    type = ADRankTwoAux
    rank_two_tensor = oxide_total_strain
    index_i = 1
    index_j = 1
    variable = oxide_strain_yy
  [../]
  [./oxide_strain_xy]
    type = ADRankTwoAux
    rank_two_tensor = oxide_total_strain
    index_i = 0
    index_j = 1
    variable = oxide_strain_xy
  [../]
  [./metal_strain_xx]
    type = ADRankTwoAux
    rank_two_tensor = metal_total_strain
    index_i = 0
    index_j = 0
    variable = metal_strain_xx
  [../]
  [./metal_strain_yy]
    type = ADRankTwoAux
    rank_two_tensor = metal_total_strain
    index_i = 1
    index_j = 1
    variable = metal_strain_yy
  [../]
  [./metal_strain_xy]
    type = ADRankTwoAux
    rank_two_tensor = metal_total_strain
    index_i = 0
    index_j = 1
    variable = metal_strain_xy
  [../]
[]

[Constraints]
  [./dispx_constraint]
    type = XFEMSingleVariableConstraint
    use_displaced_mesh = false
    variable = disp_x
    alpha = 1e20 #1e8 #1e18
    use_penalty = true
    geometric_cut_userobject = 'level_set_cut_uo'
  [../]
  [./dispy_constraint]
    type = XFEMSingleVariableConstraint
    use_displaced_mesh = false
    variable = disp_y
    alpha = 1e20 #1e8 #1e18
    use_penalty = true
    geometric_cut_userobject = 'level_set_cut_uo'
  [../]
[]

[BCs]
  [./fixed_x_pivot]
    type = ADDirichletBC
    boundary = 101
    variable = disp_x
    value = 0.0
  [../]
  [./fixed_y_pivot]
    type = ADDirichletBC
    boundary = 101
    variable = disp_y
    value = 0.0
  [../]
  [./fixed_symmetry]
    type = ADDirichletBC
    boundary = 100
    variable = disp_y
    value = 0.0
  [../]
[]

[Materials]
  [./elasticity_tensor_oxide]
    type = ADComputeIsotropicElasticityTensor
    base_name = oxide
    youngs_modulus = 2.72e11 ##272 GPa P. Tortorelli. Mechanical properties of chromia scales. Journal de Physique IV Colloque, 1993
    poissons_ratio = 0.3
  [../]
  [./strain_oxide]
    type = ADComputeFiniteStrain
    base_name = oxide
    displacements = 'disp_x disp_y'
    eigenstrain_names = 'pseudo_eigenstrain'
  [../]
  [./stress_oxide]
    type = ADComputeFiniteStrainElasticStress
    base_name = oxide
  [../]
  # [./eigenstrain_oxide]
  #   type = ComputeVariableEigenstrain
  #   args = oxide_interface
  #   base_name = oxide
  #   eigen_base = '0.0262 0.0262 0.0262 0 0 0' #10% eigenstrain approximation
  #   eigenstrain_name = 'pseudo_eigenstrain'
  # [../]
  [./eigenstrain_oxide]
    type = ADComputeThermalExpansionEigenstrain
    stress_free_temperature = 701.0
    temperature = temperature
    thermal_expansion_coeff = 2.62e-6 ##stand-in so that I can see the switch in the eigenstrains
    base_name = oxide
    eigenstrain_name = 'pseudo_eigenstrain'
  [../]
  [./elasticity_tensor_metal]
    type = ADComputeIsotropicElasticityTensor
    base_name = metal
    youngs_modulus = 1.93e11 ## 193 GPa
    poissons_ratio = 0.3
  [../]
  [./strain_metal]
    type = ADComputeFiniteStrain
    base_name = metal
    displacements = 'disp_x disp_y'
    # eigenstrain_names = 'no_eigenstrain'
  [../]
  [./stress_metal]
    type = ADComputeMultipleInelasticStress
    inelastic_models = rom_stress_prediction
    base_name = metal
  [../]
  [./rom_stress_prediction]
    type = P91LAROMANCEStressUpdate
    temperature = temperature
    stress_input_window_failure = IGNORE
    cell_input_window_failure = WARN
    wall_input_window_failure = WARN
    temperature_input_window_failure = WARN
    environment_input_window_failure = WARN #This is the phase fraction, below
    environmental_factor = mx_phase_fraction
    initial_cell_dislocation_density = 4.0e12 #bounds 1e12 to 6e12
    initial_wall_dislocation_density = 1.0e13 # bounds 6e12 1.8e13
  [../]
  [./mx_phase_fraction]
    type = ADGenericConstantMaterial
    prop_names = mx_phase_fraction
    prop_values = 8.0e-3  #precipitation bounds: 6e-3, 1e-1
  [../]
  # [./eigenstrain_metal]
  #   type = ComputeThermalExpansionEigenstrain
  #   stress_free_temperature = 800.0
  #   temperature = temperature
  #   thermal_expansion_coeff = 0.0 ## faking it out so that I can see the switch in the eigenstrains
  #   base_name = B
  #   eigenstrain_name = 'no_eigenstrain'
  # [../]
  [./combined_stress]
    type = ADLevelSetBiMaterialRankTwo
    levelset_positive_base = 'metal'
    levelset_negative_base = 'oxide'
    level_set_var = oxide_interface
    prop_name = stress
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter'
  petsc_options_value = '201                hypre    boomeramg      8'
  automatic_scaling = false #Don't use this option because PETSc doesn't like it
  line_search = 'none'

# controls for linear iterations
  l_max_its = 50
  l_tol = 1e-4

# controls for nonlinear iterations
  nl_max_its = 100
  nl_rel_tol = 5e-4
  nl_abs_tol = 1e-8 #1e-20

# time control
  start_time = 0.0
 dt = 86400.0 #86400.0
  dtmin = 1.0
  # dtmax = 10000.0
  #num_steps = 10
  end_time = 5184000 #60 days
  # end_time = 7776000 #90 days
  # end_time = 15552000 #180 days

  max_xfem_update = 1
[]

[Outputs]
  exodus = true
  execute_on = timestep_end
  csv = true
  perf_graph = true
  checkpoint = true
  [./console]
    type = Console
    output_linear = true
  [../]
[]

[Postprocessors]
  [./oxide_position]
    type = 'FunctionValuePostprocessor'
    function = estimated_oxide_growth
  [../]
[]
