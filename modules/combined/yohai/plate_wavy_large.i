a = 13.8
b = 14
c = 25
ymax = 0.96
E_oxide = 1.9e5
E_metal = 1.9e5
nu_oxide = 0.3
nu_metal = 0.3
CTE_oxide = 2.1e-5
CTE_metal = 1.4e-5
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
    nx = '${fparse ceil((14.17-a)/e_A)-1}'
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
  minimum_weight_multiplier = 5e-2
[]

[UserObjects]
  [cut]
    type = LevelSetCutUserObject
    level_set_var = interface
    heal_always = true
  []
[]

[AuxKernels]
  [interface]
    type = FunctionAux
    variable = interface
    function = 'x-${b}-${A}*sin(2*pi/${lambda}*y)-0.01*y'
  []
  [temperature]
    type = FunctionAux
    variable = temperature
    function = '${T_ref}+t*(${T_fin}-${T_ref})/${end_time}'
  []
  [stress_rr]
    type = ADRankTwoAux
    variable = 'stress_rr'
    rank_two_tensor = 'stress'
    index_i = 0
    index_j = 0
  []
  [stress_zz]
    type = ADRankTwoAux
    variable = 'stress_zz'
    rank_two_tensor = 'stress'
    index_i = 1
    index_j = 1
  []
  [stress_tt]
    type = ADRankTwoAux
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
    bound_value = 0.99
  []
[]

[Kernels]
  [solid_r]
    type = ADStressDivergenceRZTensors
    variable = 'disp_r'
    component = 0
    use_displaced_mesh = true
  []
  [solid_z]
    type = ADStressDivergenceRZTensors
    variable = 'disp_z'
    component = 1
    use_displaced_mesh = true
  []
  [ac_interface]
    type = ADACInterface
    variable = 'c'
    variable_L = false
    mob_name = oxide_L
    kappa_name = oxide_kappa_op
  []
  [ac]
    type = ADAllenCahn
    variable = 'c'
    f_name = F
    mob_name = oxide_L
  []
[]

[Constraints]
  [disp_r_constraint]
    type = XFEMSingleVariableConstraint
    variable = 'disp_r'
    # indicator_variable = 'c'
    # threshold = 0.25
    geometric_cut_userobject = 'cut'
    alpha = 1e10
    use_displaced_mesh = true
  []
  [disp_z_constraint]
    type = XFEMSingleVariableConstraint
    variable = 'disp_z'
    # indicator_variable = 'c'
    # threshold = 0.25
    geometric_cut_userobject = 'cut'
    alpha = 1e10
    use_displaced_mesh = true
  []
  [disp_z_fixed]
    type = EqualValueBoundaryConstraint
    variable = 'disp_z'
    secondary = top
    formulation = penalty
    penalty = 1e11
    primary = 6836
  []
[]

[BCs]
  [fixed_z]
    type = DirichletBC
    variable = 'disp_z'
    boundary = 'bottom'
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
    type = ADDerivativeParsedMaterial
    f_name = oxide_degradation
    args = 'c'
    # function = '(1-c)^2'
    material_property_names = 'oxide_psic oxide_l oxide_gc_prop'
    function = 'M:=3*oxide_gc_prop/8/oxide_l; (1-c)^2/((1-c)^2+M/oxide_psic*c*(1+c))'
    derivative_order = 2
  []
  [define_mobility]
    type = ADParsedMaterial
    f_name = oxide_L
    material_property_names = 'oxide_gc_prop'
    function = '1/oxide_gc_prop'
  []
  [define_kappa]
    type = ADParsedMaterial
    f_name = oxide_kappa_op
    material_property_names = 'oxide_gc_prop oxide_l'
    # function = 'oxide_gc_prop*oxide_l'
    function = '0.75*oxide_gc_prop*oxide_l'
  []
  [pfbulkmat]
    type = ADGenericConstantMaterial
    prop_names = 'oxide_psic oxide_gc_prop oxide_l oxide_bulk_modulus oxide_shear_modulus'
    prop_values = '0.2 0.05 0.0132 ${fparse E_oxide/3/(1-2*nu_oxide)} ${fparse '
                  'E_oxide/2/(1+nu_oxide)}'
  []
  [pfdummy]
    type = ADGenericConstantMaterial
    prop_names = 'metal_elastic_energy dmetal_elastic_energy/dc d^2metal_elastic_energy/dc^2'
    prop_values = '0 0 0'
  []
  [local_fracture_energy]
    type = ADDerivativeParsedMaterial
    f_name = local_fracture_energy
    args = 'c'
    material_property_names = 'oxide_gc_prop oxide_l'
    function = 'c*0.375*oxide_gc_prop/oxide_l'
    derivative_order = 2
  []
  [fracture_driving_energy]
    type = ADDerivativeSumMaterial
    f_name = F
    args = c
    sum_materials = 'elastic_energy local_fracture_energy'
    derivative_order = 2
  []
  # oxide
  [eigenstrain_oxide]
    type = ADComputeThermalExpansionEigenstrain
    base_name = oxide
    stress_free_temperature = ${T_ref}
    temperature = temperature
    thermal_expansion_coeff = ${CTE_oxide}
    eigenstrain_name = 'thermal_eigenstrain'
  []
  [strain_oxide]
    type = ADComputeDeformationGradient
    base_name = oxide
    eigenstrain_names = 'thermal_eigenstrain'
  []
  [stress_oxide]
    type = ADComputeHyperelasticPFFractureStress
    base_name = oxide
    c = 'c'
    history_maximum = false
  []
  # metal
  [elasticity_tensor_metal]
    type = ADComputeIsotropicElasticityTensor
    base_name = metal
    youngs_modulus = ${E_metal}
    poissons_ratio = ${nu_metal}
  []
  [eigenstrain_metal]
    type = ADComputeThermalExpansionEigenstrain
    base_name = metal
    stress_free_temperature = ${T_ref}
    temperature = temperature
    thermal_expansion_coeff = ${CTE_metal}
    eigenstrain_name = 'thermal_eigenstrain'
  []
  [strain_metal]
    type = ADComputeAxisymmetricRZFiniteStrain
    base_name = metal
    eigenstrain_names = 'thermal_eigenstrain'
    decomposition_method = TaylorExpansion
  []
  [stress_metal]
    type = ADComputeFiniteStrainElasticStress
    base_name = metal
  []
  # bimaterial
  [combined_stress]
    type = ADLevelSetBiMaterialRankTwo
    levelset_positive_base = 'metal'
    levelset_negative_base = 'oxide'
    level_set_var = 'interface'
    prop_name = stress
  []
  [combined_elastic_energy]
    type = ADLevelSetBiMaterialRealDerivatives
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
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type'
  petsc_options_value = 'lu       mumps                         vinewtonrsls'
  # petsc_options_iname = '-pc_type -sub_pc_type -ksp_max_it -ksp_gmres_restart -sub_pc_factor_levels '
  #                       '-snes_type'
  # petsc_options_value = 'asm      ilu          1000        200                0                     '
  #                       'vinewtonrsls'
  automatic_scaling = true

  # line_search = none

  nl_rel_tol = 1e-06
  nl_abs_tol = 1e-08

  dt = 0.5
  end_time = ${end_time}

  max_xfem_update = 1
[]

[Outputs]
  print_linear_converged_reason = false
  print_nonlinear_converged_reason = false
  print_linear_residuals = false
  [exodus]
    type = Exodus
    file_base = 'output/plate_wavy_large'
  []
[]
