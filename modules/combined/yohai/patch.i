E_oxide = 1.9e2
E_metal = 1.9e5
nu_oxide = 0.3
nu_metal = 0.3
CTE_oxide = 2.1e-4
CTE_metal = 1.4e-5
T_ref = '${fparse 580+273}'
T_fin = '${fparse 80+273}'
end_time = 50

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
    nx = 2
    ny = 2
    xmin = 10
    xmax = 30
    ymin = 0
    ymax = 1
  []
[]

[Variables]
  [disp_r]
  []
  [disp_z]
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
  [c]
  []
[]

[XFEM]
  qrule = volfrac
  output_cut_plane = true
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
    function = 'x-19.25-y'
  []
  [temperature]
    type = FunctionAux
    variable = temperature
    function = '${T_ref}+t*(${T_fin}-${T_ref})/${end_time}'
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
[]

[Constraints]
  [disp_r_constraint]
    type = XFEMSingleVariableConstraint
    variable = 'disp_r'
    geometric_cut_userobject = 'cut'
    alpha = 1e8
    use_penalty = true
    use_displaced_mesh = true
  []
  [disp_z_constraint]
    type = XFEMSingleVariableConstraint
    variable = 'disp_z'
    geometric_cut_userobject = 'cut'
    alpha = 1e8
    use_penalty = true
    use_displaced_mesh = true
  []
  [disp_z_fixed]
    type = EqualValueBoundaryConstraint
    variable = 'disp_z'
    secondary = top
    formulation = penalty
    penalty = 1e11
  []
[]

[BCs]
  [fixed_z]
    type = DirichletBC
    variable = 'disp_z'
    boundary = 'bottom'
    value = 0.0
  []
[]

[Materials]
  # fracture
  [degradation]
    type = ADDerivativeParsedMaterial
    f_name = oxide_degradation
    args = 'c'
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
    function = '0.75*oxide_gc_prop*oxide_l'
  []
  [pfbulkmat]
    type = ADGenericConstantMaterial
    prop_names = 'oxide_psic oxide_gc_prop oxide_l oxide_bulk_modulus oxide_shear_modulus'
    prop_values = '0.2 0.02 0.01 ${fparse E_oxide/3/(1-2*nu_oxide)} ${fparse E_oxide/2/(1+nu_oxide)}'
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
    function = 'c*3*oxide_gc_prop/8/oxide_l'
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
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist                 '
  automatic_scaling = true

  nl_rel_tol = 1e-06
  nl_abs_tol = 1e-08

  dt = 1
  end_time = ${end_time}

  max_xfem_update = 1

  abort_on_solve_fail = true
[]

[Outputs]
  print_linear_converged_reason = false
  print_nonlinear_converged_reason = false
  print_linear_residuals = false
  [exodus]
    type = Exodus
    file_base = 'output/patch'
  []
[]
