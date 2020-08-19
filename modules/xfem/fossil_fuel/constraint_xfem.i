# This test is for two layer materials with different youngs modulus
# and different eigenstrains
# The global stress is determined by switching the stress based on level set values
# The material interface is marked by a level set function
# The two layer materials are glued together

T_ref = 753
CTE_oxide = 1e-2
E_metal = 1.93e5
E_oxide = 2.72e5
nu_metal = 0.3
nu_oxide = 0.3

[GlobalParams]
  displacements = 'disp_r disp_z'
  # decomposition_method = EigenSolution
[]

[Problem]
  type = FEProblem
  coord_type = RZ
[]

[Mesh]
  [generated_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    elem_type = QUAD4
    nx = 88
    ny = 80
    xmin = 13.8
    xmax = 25
    ymin = 0
    ymax = 10
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
  [stress]
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
    function = 'if(t<10,x-14-0.1*sin(4*y),x-14-0.1*sin(4*y)-0.04781*(t-10))'
  []
  [temperature]
    type = FunctionAux
    variable = temperature
    function = 'if(t<10,${T_ref}+t,${T_ref}+10)'
  []
  [stress]
    type = ADRankTwoAux
    variable = 'stress'
    rank_two_tensor = 'stress'
    index_i = 1
    index_j = 1
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
    penalty = 1e10
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
    eigenstrain_name = 'eigenstrain_oxide'
  []
  [strain_oxide]
    type = ADComputeAxisymmetricRZSmallStrain
    base_name = oxide
    eigenstrain_names = 'eigenstrain_oxide'
  []
  [stress_oxide]
    type = ADComputeLinearElasticStress
    base_name = oxide
  []
  # bimaterial
  [combined_stress]
    type = ADLevelSetBiMaterialRankTwo
    levelset_positive_base = 'metal'
    levelset_negative_base = 'oxide'
    level_set_var = 'interface'
    prop_name = stress
  []
  # metal
  [elasticity_tensor_metal]
    type = ADComputeIsotropicElasticityTensor
    base_name = metal
    youngs_modulus = ${E_metal}
    poissons_ratio = ${nu_metal}
  []
  [strain_metal]
    type = ADComputeAxisymmetricRZSmallStrain
    base_name = metal
  []
  [stress_metal]
    type = ADComputeLinearElasticStress
    base_name = metal
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  # petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter'
  # petsc_options_value = '201                hypre    boomeramg      8'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist'
  automatic_scaling = true
  # line_search = 'none'

  # controls for nonlinear iterations
  nl_rel_tol = 1e-06
  nl_abs_tol = 1e-08

  # time control
  dt = 1
  end_time = 210

  max_xfem_update = 1
[]

[Outputs]
  print_linear_converged_reason = false
  print_nonlinear_converged_reason = false
  print_linear_residuals = false
  [exodus]
    type = Exodus
    file_base = 'output/constraint_xfem'
  []
[]
