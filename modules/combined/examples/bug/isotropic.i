[Mesh]
  [Mesh]
    type = GeneratedMeshGenerator
    nx = 50
    ny = 50
    nz = 50
    dim = 3
    xmax = 200
    ymax = 200
    zmax = 200
  []
  [Outer_Surface]
    type = SideSetsAroundSubdomainGenerator
    input = Mesh
    new_boundary = 'OuterSurface'
    block = 0
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]
[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [disp_z]
  []
  [T]
  []
[]
[AuxVariables]
  [vonMises]
    order = CONSTANT
    family = MONOMIAL
    [AuxKernel]
      type = ADRankTwoScalarAux
      rank_two_tensor = stress
      variable = vonMises
      scalar_type = VonMisesStress
      execute_on = 'TIMESTEP_END'
    []
  []
  [totalStrain]
    order = CONSTANT
    family = MONOMIAL
    [AuxKernel]
      type = ADRankTwoScalarAux
      scalar_type = L2norm
      variable = totalStrain
      execute_on = 'TIMESTEP_END'
      rank_two_tensor = total_strain
    []
  []
  [TRef]
    order = FIRST
    family = LAGRANGE
  []
  [plasticStrain]
    order = CONSTANT
    family = MONOMIAL
    [AuxKernel]
      type = ADMaterialRealAux
      property = 'effective_plastic_strain'
      variable = plasticStrain
      execute_on = 'TIMESTEP_END'
    []
  []
  [eigen_xx]
    order = CONSTANT
    family = MONOMIAL
    [AuxKernel]
      type = ADRankTwoAux
      variable = eigen_xx
      rank_two_tensor = eigenstrain
      index_i = 0
      index_j = 0
    []
  []
  [eigen_yy]
    order = CONSTANT
    family = MONOMIAL
    [AuxKernel]
      type = ADRankTwoAux
      variable = eigen_yy
      rank_two_tensor = eigenstrain
      index_i = 1
      index_j = 1
    []
  []
  [eigen_zz]
    order = CONSTANT
    family = MONOMIAL
    [AuxKernel]
      type = ADRankTwoAux
      variable = eigen_zz
      rank_two_tensor = eigenstrain
      index_i = 2
      index_j = 2
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
  [heatConduction]
    type = ADHeatConduction
    variable = T
    use_displaced_mesh = true
  []
  [heatConductionTimeDerivative]
    type = ADHeatConductionTimeDerivative
    variable = T
    use_displaced_mesh = true
  []
[]

[Materials]
  [Density]
    type = ADDensity
    density = 4.43E-06
  []
  [invariant_Properties]
    # Poissons ratio is not a function of temperature for this material
    type = ADGenericConstantMaterial
    prop_names = 'poissons_ratio youngs_modulus thermal_expansion thermal_conductivity specific_heat'
    prop_values = '0.31 77223.9 9.81e-6 0.013494 774.558'
  []
  [MMPDS_Elasticity_tensor]
    # Compute the varying elasticity tensor.
    type = ADComputeVariableIsotropicElasticityTensor
    youngs_modulus = 'youngs_modulus'
    poissons_ratio = 'poissons_ratio'
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
    type = ADComputeVariableThermalExpansionEigenstrain
    thermal_expansion_coeff = thermal_expansion
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
  [Fixed_X]
    type = ADDirichletBC
    boundary = 'left'
    value = 0
    variable = disp_x
  []
  [Fixed_Y]
    type = ADDirichletBC
    boundary = 'left'
    value = 0
    variable = disp_y
  []
  [Fixed_Z]
    type = ADDirichletBC
    boundary = 'left'
    value = 0
    variable = disp_z
  []
  [left_Temp]
    type = ADDirichletBC
    boundary = 'left'
    variable = T
    value = 2000
  []
  [Tip_Temp]
    type = ADDirichletBC
    boundary = 'right'
    variable = T
    value = 300
  []
  [Radiative_Cooling]
    type = ADFunctionRadiativeBC
    variable = T
    boundary = 'OuterSurface'
    Tinfinity = 300.0
    stefan_boltzmann_constant = 5.67e-14
    emissivity_function = 0.6
    use_displaced_mesh = true
  []
  [Convective_Cooling]
    type = ADConvectiveHeatFluxBC
    variable = T
    boundary = 'OuterSurface'
    T_infinity = 300.0
    heat_transfer_coefficient = 0.00002
    use_displaced_mesh = true
  []
[]

[ICs]
  [initialT]
    type = ConstantIC
    variable = T
    value = 1100
  []
  [referenceT]
    type = ConstantIC
    variable = TRef
    value = 1100
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

  dt = 15
  end_time = 315
  automatic_scaling = true

  solve_type = 'NEWTON'

  nl_rel_tol = 1e-08
  nl_abs_tol = 5e-13

  petsc_options_iname = '-pc_type  -pc_hypre_type -pc_hypre_boomeramg_strong_threshold  -pc_hypre_boomeramg_coarsen_type  -pc_hypre_boomeramg_interp_type'
  petsc_options_value = 'hypre    boomeramg 0.8   PMIS  ext+i'
  line_search = NONE
  [Predictor]
    type = SimplePredictor
    scale = 1
  []
[]
[Postprocessors]
  [Stress]
    type = PointValue
    point = '${measurePoint0} ${measurePoint1} ${measurePoint2}'
    variable = vonMises
  []
  [Total_Strain]
    type = PointValue
    point = '${measurePoint0} ${measurePoint1} ${measurePoint2}'
    variable = totalStrain
  []
  [Plastic_Strain]
    type = PointValue
    point = '${measurePoint0} ${measurePoint1} ${measurePoint2}'
    variable = plasticStrain
  []
  [Eigen_xx]
    type = PointValue
    point = '${measurePoint0} ${measurePoint1} ${measurePoint2}'
    variable = eigen_xx
  []
  [Eigen_yy]
    type = PointValue
    point = '${measurePoint0} ${measurePoint1} ${measurePoint2}'
    variable = eigen_yy
  []
  [Eigen_zz]
    type = PointValue
    point = '${measurePoint0} ${measurePoint1} ${measurePoint2}'
    variable = eigen_zz
  []
[]
[Outputs]
  [Exodus]
    type = Exodus
  []
  [Performance]
    type = PerfGraphOutput
  []
[]
