nx = 79
ny = 80
a = 13.8
b = 14
c = 25
ymax = 0.24
E_oxide = 1.9e5
E_metal = 1.9e5
nu_oxide = 0.3
nu_metal = 0.3
CTE_oxide = 2.1e-5
CTE_metal = 1.4e-5
T_ref = '${fparse 580+273}'
T_fin = '${fparse 480+273}'
end_time = 10

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
    nx = ${nx}
    ny = ${ny}
    xmin = ${a}
    xmax = 14.1
    ymin = 0
    ymax = ${ymax}
  []
  [metal]
    type = GeneratedMeshGenerator
    dim = 2
    elem_type = QUAD4
    nx = ${nx}
    ny = ${ny}
    xmin = 14.1
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
[]

[XFEM]
  qrule = moment_fitting
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
    function = 'x-${b}'
    # function = 'if(t<10,x-14-0.1*sin(4*y),x-14-0.1*sin(4*y)-0.04781*(t-10))'
  []
  [temperature]
    type = FunctionAux
    variable = temperature
    function = '${T_ref} + t * (${T_fin} - ${T_ref}) / ${end_time}'
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
    alpha = 1e10
    use_penalty = true
    use_displaced_mesh = true
  []
  [disp_z_constraint]
    type = XFEMSingleVariableConstraint
    variable = 'disp_z'
    geometric_cut_userobject = 'cut'
    alpha = 1e10
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
    type = ADComputeAxisymmetricRZFiniteStrain
    base_name = oxide
    eigenstrain_names = 'thermal_eigenstrain'
  []
  [stress_oxide]
    type = ADComputeFiniteStrainElasticStress
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
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist'
  automatic_scaling = true

  # controls for nonlinear iterations
  nl_rel_tol = 1e-08
  nl_abs_tol = 1e-10

  # time control
  dt = 2
  end_time = ${end_time}

  max_xfem_update = 1
[]

[Outputs]
  print_linear_converged_reason = false
  print_nonlinear_converged_reason = false
  print_linear_residuals = false
  [exodus]
    type = Exodus
    file_base = 'output/plate'
    execute_on = 'FINAL'
  []
[]
