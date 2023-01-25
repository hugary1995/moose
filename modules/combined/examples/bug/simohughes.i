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
  large_kinematics = true
  eigenstrain_names = "eigenstrain"
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
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoScalarAux
      rank_two_tensor = cauchy_stress
      variable = vonmises
      scalar_type = VonMisesStress
      execute_on = timestep_end
    []
  []
  [plasticStrain]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = MaterialRealAux
      property = effective_plastic_strain
      variable = plasticStrain
      execute_on = 'TIMESTEP_END'
    []
  []
  [totalStrain]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoScalarAux
      variable = totalStrain
      scalar_type = L2norm
      rank_two_tensor = total_strain
      execute_on = 'TIMESTEP_END'
    []
  []
  [TRef]
    family = LAGRANGE
    order = FIRST
  []
  [eigen_xx]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      variable = eigen_xx
      rank_two_tensor = eigenstrain
      index_i = 0
      index_j = 0
    []
  []
  [eigen_yy]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      variable = eigen_yy
      rank_two_tensor = eigenstrain
      index_i = 1
      index_j = 1
    []
  []
  [eigen_zz]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      variable = eigen_zz
      rank_two_tensor = eigenstrain
      index_i = 2
      index_j = 2
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
  [heatConduction]
    type = HeatConduction
    variable = T
    use_displaced_mesh = true
  []
  [heatConductionTimeDiv]
    type = HeatConductionTimeDerivative
    variable = T
    use_displaced_mesh = true
  []
[]
[Materials]
  [invariant_Properties]
    # Poissons ratio is not a function of temperature for this material
    type = GenericConstantMaterial
    prop_names = 'thermal_expansion thermal_conductivity specific_heat'
    prop_values = '9.81e-6 0.013494 774.558'
  []
  [elasticity_tensor]
    # Compute the varying elasticity tensor.
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 77223.9
    poissons_ratio = 0.31
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
  [Strain]
    type = ComputeLagrangianStrain
    eigenstrain_names = "eigenstrain"
  []
  [MMPDS_Density]
    type = Density
    density = 4.43E-06
  []
  [Thermal_Strain]
    # Thermal expansion.
    type = ComputeVariableThermalExpansionEigenstrain
    thermal_expansion_coeff = thermal_expansion
    temperature = T
    eigenstrain_name = eigenstrain
    stress_free_temperature = 1100
  []
[]
[BCs]
  [Fixed_X]
    type = DirichletBC
    boundary = 'left'
    value = 0
    variable = disp_x
  []
  [Fixed_Y]
    type = DirichletBC
    boundary = 'left'
    value = 0
    variable = disp_y
  []
  [Fixed_Z]
    type = DirichletBC
    boundary = 'left'
    value = 0
    variable = disp_z
  []
  [left_Temp]
    type = DirichletBC
    boundary = 'left'
    variable = T
    value = 2000
  []
  [Right_Temp]
    type = DirichletBC
    boundary = 'right'
    variable = T
    value = 300
  []
  [Radiative_Cooling]
    type = FunctionRadiativeBC
    variable = T
    boundary = 'OuterSurface'
    Tinfinity = 300.0
    stefan_boltzmann_constant = 5.67e-14
    emissivity_function = 0.6
    use_displaced_mesh = true
  []
  [Convective_Cooling]
    type = ConvectiveHeatFluxBC
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
  [TotalStrain]
    type = PointValue
    point = '${measurePoint0} ${measurePoint1} ${measurePoint2}'
    variable = totalStrain
  []
  [PlasticStrain]
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
