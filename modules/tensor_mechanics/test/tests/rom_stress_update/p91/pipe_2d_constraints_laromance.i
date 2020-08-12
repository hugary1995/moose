## Version of a 2D line-type model through the wall thickness of a pipe
## Does include the constraint in the y-direction

[GlobalParams]
  order = SECOND
  family = LAGRANGE
  displacements = 'disp_x disp_y'
[]

[Problem]
  coord_type = RZ
[]

[Mesh]
  [pipe]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0.0478536 #1.884 in inches, LAROMANCE assumes meshes in m
    xmax = 0.0529336 #2.084 in inches, LAROMANCE assumes meshes in m
    ymin = 0
    ymax = 0.001016 #0.04 in inches, LAROMANCE assumes meshes in m
    nx = 5
    ny = 1
    elem_type = quad9
  []
[]

[Modules/TensorMechanics/Master]
  [all]
    strain = FINITE
    add_variables = true
    generate_output = 'max_principal_stress max_principal_strain vonmises_stress
                      strain_xx strain_yy strain_zz
                      stress_xx stress_yy stress_zz'
    use_automatic_differentiation = true
  []
[]

[AuxVariables]
  [temperature]
    initial_condition = 873
  []
  [sint]
    order = CONSTANT
    family = MONOMIAL
  []
  [effective_inelastic_strain]
    order = FIRST
    family = MONOMIAL
  []
  [cell_dislocations]
    order = FIRST
    family = MONOMIAL
  []
  [wall_dislocations]
    order = FIRST
    family = MONOMIAL
  []
[]

[Functions]
  [inner_pressure]
   type = PiecewiseLinear
   x = '0      2e8'  #time
   y = '1.0e6 1.0e6'  #pressure in Pa
   # y = '0.0 1.0 1.0'  #pressure in MPa
   scale_factor = 5.0
  []
  [timefunc]
    type=ParsedFunction
    value = t
  []
[]

[AuxKernels]
  [sint]
    type = ADRankTwoScalarAux
    rank_two_tensor = stress
    variable = sint
    scalar_type = stressIntensity
    execute_on = timestep_end
  []
  [effective_inelastic_strain]
    type = ADMaterialRealAux
    variable = effective_inelastic_strain
    property = effective_creep_strain
  []
  [cell_dislocations]
    type = ADMaterialRealAux
    variable = cell_dislocations
    property = cell_dislocations
  []
  [wall_dislocations]
    type = ADMaterialRealAux
    variable = wall_dislocations
    property = wall_dislocations
  []
[]

[BCs]
  [fixBottom]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
    preset = true
  []
  [Pressure]
    [inside]
      boundary = left
      function = inner_pressure
    []
    [top]
      boundary = top
      function = inner_pressure
      factor = -4.4726
    []
  []
[]

[Constraints]
  [top]
    type = EqualValueBoundaryConstraint
    variable = disp_y
    secondary = top #sideset on top boundary
    penalty = 1e+14
    formulation = penalty
  []
[]

[Materials]
  [elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 1.68e11 #at 600C, from neml
    poissons_ratio = 0.31
  []
  [stress]
    type = ADComputeMultipleInelasticStress
    inelastic_models = rom_stress_prediction
  []
  [rom_stress_prediction]
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
  []
  [mx_phase_fraction]
    type = ADGenericConstantMaterial
    prop_names = mx_phase_fraction
    prop_values = 8.0e-3  #precipitation bounds: 6e-3, 1e-1
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-ksp_type  -pc_type  -pc_factor_mat_solver_package'
  petsc_options_value = ' preonly    lu        superlu_dist'
  automatic_scaling = true
  line_search = 'none'

  l_max_its = 20
  l_tol = 1e-4
  nl_max_its = 10
  nl_abs_tol = 1e-16
  nl_rel_tol = 1e-4
  end_time = 6.3e8 #roughly 20 years
  dtmin = 0.1
  dtmax = 5e6 #1e8 lets it get to a state with dislocations below the input limit

  [Predictor]
    type = SimplePredictor
    scale = 1.0
  []
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1 ## This model requires a tiny timestep at the onset for the first 10s
    iteration_window = 4
    optimal_iterations = 12
    time_t = ' 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0 5.1 5.2 5.3 5.4 5.5 5.6 5.7 5.8 5.9 6.0 6.1 6.2 6.3 6.4 6.5 6.6 6.7 6.8 6.9 7.0 7.1 7.2 7.3 7.4 7.5 7.6 7.7 7.8 7.9 8.0 8.1 8.2 8.3 8.4 8.5 8.6 8.7 8.8 8.9 9.0 9.1 9.2 9.3 9.4 9.5 9.6 9.7 9.8 9.9 10.0'
    time_dt = '0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1'
  []
[]

[Postprocessors]
  [./num_lin_it]
    type = NumLinearIterations
  [../]
  [./num_nonlin_it]
    type = NumNonlinearIterations
  [../]
  [./tot_lin_it]
    type = CumulativeValuePostprocessor
    postprocessor = num_lin_it
  [../]
  [./tot_nonlin_it]
    type = CumulativeValuePostprocessor
    postprocessor = num_nonlin_it
  [../]


  [./end_time]
    type = 'FunctionValuePostprocessor'
    function = timefunc
  [../]

  [./max_stress_xx]
    type = ElementExtremeValue
    variable = stress_xx
    value_type = max
  [../]
  [./max_stress_yy]
    type = ElementExtremeValue
    variable = stress_yy
    value_type = max
  [../]
  [./max_stress_zz]
    type = ElementExtremeValue
    variable = stress_zz
    value_type = max
  [../]

  [./max_strain_xx]
    type = ElementExtremeValue
    variable = strain_xx
    value_type = max
  [../]
  [./max_strain_yy]
    type = ElementExtremeValue
    variable = strain_yy
    value_type = max
  [../]
  [./max_strain_zz]
    type = ElementExtremeValue
    variable = strain_zz
    value_type = max
  [../]

  [./max_principal_strain]
    type = ElementExtremeValue
    variable = max_principal_strain
    value_type = max
  [../]
  [./max_principal_stress]
    type = ElementExtremeValue
    variable = max_principal_stress
    value_type = max
  [../]

  [./effective_strain_avg]
    type = ElementAverageValue
    variable = effective_inelastic_strain
  [../]
  [./max_effective_strain]
    type = ElementExtremeValue
    variable = effective_inelastic_strain
    value_type = max
  [../]
  [./vonmises_stress_avg]
    type = ElementAverageValue
    variable = vonmises_stress
  [../]
  [./max_vonmises_stress]
    type = ElementExtremeValue
    variable = vonmises_stress
    value_type = max
  [../]

  [./temperature]
    type = ElementAverageValue
    variable = temperature
  [../]
  [./cell_dislocations]
    type = ElementAverageValue
    variable = cell_dislocations
  [../]
  [./max_cell_dislocations]
    type = ElementExtremeValue
    variable = cell_dislocations
    value_type = max
  [../]
  [./min_cell_dislocations]
    type = ElementExtremeValue
    variable = cell_dislocations
    value_type = min
  [../]
  [./wall_dislocations]
    type = ElementAverageValue
    variable = wall_dislocations
  [../]
  [./max_wall_dislocations]
    type = ElementExtremeValue
    variable = wall_dislocations
  [../]
  [./min_wall_dislocations]
    type = ElementExtremeValue
    variable = wall_dislocations
    value_type = min
  [../]
[]

[UserObjects]
  [max_strain]
    type = Terminator
    expression = 'max_principal_strain > 0.01'
  []
[]

[Outputs]
  exodus = true
  csv = true
  perf_graph = true
[]
