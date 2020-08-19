a = 13.8
b = 14
c = 25
ymax = 0.96
E_oxide = 4e4
E_metal = 2e5
nu_oxide = 0.3
nu_metal = 0.3
CTE_oxide = 2.1e-4
CTE_metal = 1.4e-4
T_ref = '${fparse 580+273}'
T_fin = '${fparse 80+273}'
end_time = 50
A = 0.02
lambda = 0.08
e_A = '${fparse A/2}'

[GlobalParams]
  displacements = 'disp_r disp_z'
[]

[Problem]
  type = FEProblem
  coord_type = RZ
[]

[Mesh]
  [oxide]
    type = GeneratedMeshGenerator
    dim = 2
    elem_type = QUAD4
    nx = '${fparse ceil((14.03-a)/e_A)-1}'
    ny = 160
    xmin = ${a}
    xmax = 14.03
    ymin = 0
    ymax = ${ymax}
  []
  [metal]
    type = GeneratedMeshGenerator
    dim = 2
    elem_type = QUAD4
    nx = 80
    ny = 10
    xmin = 14.03
    xmax = ${c}
    ymin = 0
    ymax = ${ymax}
  []
  [stitch]
    type = StitchedMeshGenerator
    inputs = 'oxide metal'
    stitch_boundaries_pairs = 'right left'
  []
[]

[Variables]
  [disp_r]
  []
  [disp_z]
  []
  [c]
  []
[]

[AuxVariables]
  [interface]
  []
  [temperature]
    [InitialCondition]
      type = ConstantIC
      value = ${T_ref}
    []
  []
  [stress_rr]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_zz]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_tt]
    order = CONSTANT
    family = MONOMIAL
  []
  [bounds_dummy]
  []
[]

[XFEM]
  qrule = volfrac
  output_cut_plane = true
  minimum_weight_multiplier = 0.2
[]

[UserObjects]
  [cut]
    type = LevelSetCutUserObject
    level_set_var = interface
    heal_always = false
  []
[]

[AuxKernels]
  [interface]
    type = FunctionAux
    variable = interface
    function = 'x-${b}-${A}*sin(2*pi/${lambda}*y)'
  []
  [temperature]
    type = FunctionAux
    variable = temperature
    function = '${T_ref}+t*(${T_fin}-${T_ref})/${end_time}'
  []
  [stress_rr]
    type = RankTwoAux
    variable = 'stress_rr'
    rank_two_tensor = 'stress'
    index_i = 0
    index_j = 0
  []
  [stress_zz]
    type = RankTwoAux
    variable = 'stress_zz'
    rank_two_tensor = 'stress'
    index_i = 1
    index_j = 1
  []
  [stress_tt]
    type = RankTwoAux
    variable = 'stress_tt'
    rank_two_tensor = 'stress'
    index_i = 2
    index_j = 2
  []
[]

[Bounds]
  [c_lb]
    type = VariableOldValueBoundsAux
    variable = 'bounds_dummy'
    bounded_variable = 'c'
    bound_type = lower
  []
  [c_ub]
    type = ConstantBoundsAux
    variable = 'bounds_dummy'
    bounded_variable = 'c'
    bound_type = upper
    bound_value = 1
  []
[]

[Kernels]
  [solid_r]
    type = StressDivergenceRZTensors
    variable = 'disp_r'
    component = 0
    use_displaced_mesh = true
  []
  [solid_z]
    type = StressDivergenceRZTensors
    variable = 'disp_z'
    component = 1
    use_displaced_mesh = true
  []
  [offr]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_r
    component = 0
    c = c
  []
  [offz]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_z
    component = 1
    c = c
  []
  [ac_interface]
    type = ACInterface
    variable = 'c'
    variable_L = false
    mob_name = L
    kappa_name = kappa_op
  []
  [ac]
    type = AllenCahn
    variable = 'c'
    f_name = F
    mob_name = L
  []
[]

[Constraints]
  [disp_r_constraint]
    type = XFEMSingleVariableConstraintThreshold
    variable = 'disp_r'
    indicator_variable = 'c'
    threshold = 0.9
    geometric_cut_userobject = 'cut'
    alpha = 1e10
    use_displaced_mesh = true
  []
  [disp_z_constraint]
    type = XFEMSingleVariableConstraintThreshold
    variable = 'disp_z'
    indicator_variable = 'c'
    threshold = 0.9
    geometric_cut_userobject = 'cut'
    alpha = 1e10
    use_displaced_mesh = true
  []
  # [disp_z_fixed]
  #   type = EqualValueBoundaryConstraint
  #   variable = 'disp_z'
  #   secondary = top
  #   formulation = penalty
  #   penalty = 1e11
  #   # primary = 6836
  # []
[]

[BCs]
  [fixed_z]
    type = DirichletBC
    variable = 'disp_z'
    boundary = 'bottom top'
    value = 0.0
  []
  [c_top]
    type = DirichletBC
    variable = 'c'
    boundary = 'top'
    value = 0.0
  []
[]

[Materials]
  # fracture
  [degradation]
    type = DerivativeParsedMaterial
    f_name = degradation
    args = 'c'
    # function = '(1-c)^2'
    material_property_names = 'psic l gc_prop'
    function = 'M:=3*gc_prop/8/l; (1-c)^2/((1-c)^2+M/psic*c*(1+c))'
    derivative_order = 2
  []
  [define_mobility]
    type = ParsedMaterial
    f_name = L
    material_property_names = 'gc_prop'
    function = '1/gc_prop'
  []
  [define_kappa]
    type = ParsedMaterial
    f_name = kappa_op
    material_property_names = 'gc_prop l'
    # function = 'oxide_gc_prop*oxide_l'
    function = '0.75*gc_prop*l'
  []
  [pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'psic gc_prop l oxide_bulk_modulus oxide_shear_modulus'
    prop_values = '0.4 0.05 0.015 ${fparse E_oxide/3/(1-2*nu_oxide)} ${fparse '
                  'E_oxide/2/(1+nu_oxide)}'
  []
  [pfdummy]
    type = GenericConstantMaterial
    prop_names = 'metal_elastic_energy dmetal_elastic_energy/dc d^2metal_elastic_energy/dc^2'
    prop_values = '0 0 0'
  []
  [local_fracture_energy]
    type = DerivativeParsedMaterial
    f_name = local_fracture_energy
    args = 'c'
    material_property_names = 'gc_prop l'
    function = 'c*0.375*gc_prop/l'
    derivative_order = 2
  []
  [fracture_driving_energy]
    type = DerivativeSumMaterial
    f_name = F
    args = c
    sum_materials = 'elastic_energy local_fracture_energy'
    derivative_order = 2
  []
  # oxide
  [eigenstrain_oxide]
    type = ComputeThermalExpansionEigenstrain
    base_name = oxide
    stress_free_temperature = ${T_ref}
    temperature = temperature
    thermal_expansion_coeff = ${CTE_oxide}
    eigenstrain_name = 'thermal_eigenstrain'
  []
  [elasticity_tensor_oxide]
    type = ComputeIsotropicElasticityTensor
    base_name = oxide
    youngs_modulus = ${E_oxide}
    poissons_ratio = ${nu_oxide}
  []
  [strain_oxide]
    type = ComputeAxisymmetricRZSmallStrain
    base_name = oxide
    eigenstrain_names = 'thermal_eigenstrain'
  []
  [stress_oxide]
    type = ComputeLinearElasticPFFractureStress
    c = 'c'
    E_name = oxide_elastic_energy
    I_name = local_fracture_energy
    base_name = oxide
    use_snes_vi_solver = true
    use_current_history_variable = false
    decomposition_type = none
  []
  # [strain_oxide]
  #   type = ADComputeDeformationGradient
  #   base_name = oxide
  #   eigenstrain_names = 'thermal_eigenstrain'
  # []
  # [stress_oxide]
  #   type = ADComputeHyperelasticPFFractureStress
  #   base_name = oxide
  #   c = 'c'
  #   history_maximum = false
  # []
  # metal
  [elasticity_tensor_metal]
    type = ComputeIsotropicElasticityTensor
    base_name = metal
    youngs_modulus = ${E_metal}
    poissons_ratio = ${nu_metal}
  []
  [eigenstrain_metal]
    type = ComputeThermalExpansionEigenstrain
    base_name = metal
    stress_free_temperature = ${T_ref}
    temperature = temperature
    thermal_expansion_coeff = ${CTE_metal}
    eigenstrain_name = 'thermal_eigenstrain'
  []
  [strain_metal]
    type = ComputeAxisymmetricRZSmallStrain
    base_name = metal
    eigenstrain_names = 'thermal_eigenstrain'
  []
  [stress_metal]
    type = ComputeLinearElasticStress
    base_name = metal
  []
  # bimaterial
  [combined_stress]
    type = LevelSetBiMaterialRankTwo
    levelset_positive_base = 'metal'
    levelset_negative_base = 'oxide'
    level_set_var = 'interface'
    prop_name = stress
  []
  [combined_jacobian]
    type = LevelSetBiMaterialRankFour
    levelset_positive_base = 'metal'
    levelset_negative_base = 'oxide'
    level_set_var = 'interface'
    prop_name = Jacobian_mult
  []
  [combined_elastic_energy]
    type = LevelSetBiMaterialRealDerivatives
    levelset_positive_base = 'metal'
    levelset_negative_base = 'oxide'
    level_set_var = 'interface'
    prop_name = elastic_energy
    args = c
    outputs = exodus
    output_properties = 'elastic_energy'
  []
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type'
  petsc_options_value = 'lu       superlu_dist                  vinewtonrsls'
  # petsc_options_iname = '-pc_type -sub_pc_type -ksp_max_it -ksp_gmres_restart -sub_pc_factor_levels '
  #                       '-snes_type'
  # petsc_options_value = 'asm      ilu          1000        200                0                     '
  #                       'vinewtonrsls'
  automatic_scaling = true

  # line_search = none

  nl_rel_tol = 1e-06
  nl_abs_tol = 1e-08

  dt = 0.1
  end_time = ${end_time}

  max_xfem_update = 1
[]

[Outputs]
  print_linear_converged_reason = false
  print_nonlinear_converged_reason = false
  print_linear_residuals = false
  [exodus]
    type = Exodus
    file_base = 'output/test'
  []
[]
