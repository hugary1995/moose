# geometry
a = 13.8
b = 14
c = 25
ymax = 0.24

# oxide parameters
E_oxide = 1.2e5
nu_oxide = 0.24
CTE_oxide = 2.1e-5
kappa_oxide = 3

# metal parameters
E_metal = 1.9e5
nu_metal = 0.3
CTE_metal = 1.4e-5
kappa_metal = 30

# boundary condition
h_gas = 0.1
h_steam = 2.8
T_ref = '${fparse 580+273}'
T_gas = '${fparse 1100+273}'
T_steam = '${fparse 530+273}'
p_gas = 0.2
p_steam = 17

[GlobalParams]
  displacements = 'disp_r disp_z'
[]

[Problem]
  type = FEProblem
  coord_type = RZ
[]

[Mesh]
  [refined]
    type = GeneratedMeshGenerator
    dim = 2
    elem_type = QUAD4
    nx = 40
    ny = 10
    xmin = ${a}
    xmax = 14.03
    ymin = 0
    ymax = ${ymax}
  []
  [coarse]
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
    inputs = 'refined coarse'
    stitch_boundaries_pairs = 'right left'
  []
[]

[Variables]
  [disp_r]
  []
  [disp_z]
  []
  [temp]
    [InitialCondition]
      type = ConstantIC
      value = ${T_ref}
    []
  []
[]

[AuxVariables]
  [interface]
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
  qrule = volfrac
  output_cut_plane = true
  minimum_weight_multiplier = 0.00
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
    function = 'x-${b}'
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
    use_displaced_mesh = true
    use_penalty = true
  []
  [disp_z_constraint]
    type = XFEMSingleVariableConstraint
    variable = 'disp_z'
    geometric_cut_userobject = 'cut'
    alpha = 1e10
    use_displaced_mesh = true
    use_penalty = true
  []
  [disp_z_fixed]
    type = EqualValueBoundaryConstraint
    variable = 'disp_z'
    secondary = top
    formulation = penalty
    penalty = 1e11
    # primary = 6836
  []
[]

[BCs]
  [fixed_z]
    type = DirichletBC
    variable = 'disp_z'
    boundary = 'bottom'
    value = 0.0
  []
  [Pressure]
    [left]
      boundary = 'left'
      factor = ${p_steam}
      use_automatic_differentiation = true
    []
    [right]
      boundary = 'right'
      factor = ${p_gas}
      use_automatic_differentiation = true
    []
  []
[]

[Materials]
  # oxide
  [eigenstrain_oxide]
    type = ADComputeThermalExpansionEigenstrain
    base_name = oxide
    stress_free_temperature = ${T_ref}
    temperature = temp
    thermal_expansion_coeff = ${CTE_oxide}
    eigenstrain_name = 'thermal_eigenstrain'
  []
  [elasticity_tensor_oxide]
    type = ADComputeIsotropicElasticityTensor
    base_name = oxide
    youngs_modulus = ${E_oxide}
    poissons_ratio = ${nu_oxide}
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
    temperature = temp
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
  petsc_options_value = 'lu       superlu_dist                 '
  automatic_scaling = true

  # line_search = none

  nl_rel_tol = 1e-06
  nl_abs_tol = 1e-08

  num_steps = 1

  max_xfem_update = 1
[]

[Outputs]
  print_linear_converged_reason = false
  print_nonlinear_converged_reason = false
  print_linear_residuals = false
  [exodus]
    type = Exodus
    file_base = 'output/elastic'
  []
[]
