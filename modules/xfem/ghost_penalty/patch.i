E_oxide = 1.9e2
E_metal = 1.9e5
nu_oxide = 0.3
nu_metal = 0.3
CTE_oxide = 2.1e-2
CTE_metal = 1.4e-3
T_ref = '${fparse 580+273}'
T_fin = '${fparse 80+273}'
end_time = 50
xc = 1.5

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
    xmin = 1
    xmax = 3
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
    function = 'x-${xc}'
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

[DGKernels]
  # [ghost_penalty]
  #   type = ADGhostPenalty
  #   variable = 'disp_z'
  #   eta = 1e5
  #   alpha = 10
  #   use_displaced_mesh = true
  # []
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
  # oxide
  [elasticity_tensor_oxide]
    type = ADComputeIsotropicElasticityTensor
    base_name = oxide
    youngs_modulus = ${E_oxide}
    poissons_ratio = ${nu_oxide}
  []
  [eigenstrain_oxide]
    type = ADComputeThermalExpansionEigenstrain
    base_name = oxide
    stress_free_temperature = ${T_ref}
    temperature = temperature
    thermal_expansion_coeff = ${CTE_oxide}
    eigenstrain_name = 'thermal_eigenstrain'
  []
  [strain_oxide]
    type = ADComputeAxisymmetricRZSmallStrain
    base_name = oxide
    eigenstrain_names = 'thermal_eigenstrain'
  []
  [stress_oxide]
    type = ADComputeLinearElasticStress
    base_name = oxide
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
    type = ADComputeAxisymmetricRZSmallStrain
    base_name = metal
    eigenstrain_names = 'thermal_eigenstrain'
  []
  [stress_metal]
    type = ADComputeLinearElasticStress
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
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options = '-pc_svd_monitor'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'svd'
  automatic_scaling = true

  nl_rel_tol = 1e-06
  nl_abs_tol = 1e-08

  dt = 1
  num_steps = 1
  # end_time = ${end_time}

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
