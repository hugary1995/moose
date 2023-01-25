T = 300

[Mesh]
  [Mesh]
    type = GeneratedMeshGenerator
    dim = 3
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  volumetric_locking_correction = true
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
      type = ADRankTwoAux
      rank_two_tensor = stress
      index_i = 1
      index_j = 1
      execute_on = 'INITIAL TIMESTEP_END'
    []
  []
  [strain_yy]
    order = CONSTANT
    family = MONOMIAL
    [AuxKernel]
      type = ADRankTwoAux
      rank_two_tensor = total_strain
      index_i = 1
      index_j = 1
      execute_on = 'INITIAL TIMESTEP_END'
    []
  []
[]

[Kernels]
  [sdx]
    type = ADStressDivergenceTensors
    component = 0
    variable = disp_x
  []
  [sdy]
    type = ADStressDivergenceTensors
    component = 1
    variable = disp_y
  []
  [sdz]
    type = ADStressDivergenceTensors
    component = 2
    variable = disp_z
  []
[]

[Materials]
  [MMPDS_Elasticity_tensor]
    # Compute the varying elasticity tensor.
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 77223.9
    poissons_ratio = 0.31
  []
  [Plasticity_Stress_Update]
    type = ADTemperatureDependentHardeningStressUpdate
    hardening_functions = 'hf300 hf2000'
    temperatures = "300 2000"
    temperature = T
  []
  [Stress]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'Plasticity_Stress_Update'
  []
  [Thermal_Strain]
    # Thermal expansion.
    type = ADComputeThermalExpansionEigenstrain
    thermal_expansion_coeff = 9.81e-6
    temperature = T
    eigenstrain_name = eigenstrain
    stress_free_temperature = 1100
  []
  [Strain]
    type = ADComputeFiniteStrain
    eigenstrain_names = "eigenstrain"
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

[Functions]
  [hf300]
    type = PiecewiseLinear
    x = '0	0.002	0.006	0.012	0.022	0.038	0.052'
    y = '910	937.68736	965.2664	979.05592	992.84544	1006.63496	1013.52972'
  []
  [hf2000]
    type = PiecewiseLinear
    x = '0	0.002	0.006	0.012	0.022	0.038	0.052'
    y = '45.5	46.884368	48.26332	48.952796	49.642272	50.331748	50.676486'
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
  file_base = 'isotropic_T_${T}'
  csv = true
[]
