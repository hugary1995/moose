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
nx = 3
heal_always = true

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  [generated_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    elem_type = QUAD4
    nx = ${nx}
    ny = 1
    xmin = 1
    xmax = 3
    ymin = 0
    ymax = 1
  []
  [fixed]
    type = BoundingBoxNodeSetGenerator
    input = generated_mesh
    new_boundary = pin
    bottom_left = '3 0 0'
    top_right = '3 0 0'
  []
[]

[Variables]
  [disp_x]
  []
  [disp_y]
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
    heal_always = ${heal_always}
  []
[]

[AuxKernels]
  [interface]
    type = FunctionAux
    variable = interface
    function = 'x-2'
  []
  [temperature]
    type = FunctionAux
    variable = temperature
    function = '${T_ref}+t'
  []
[]

[Kernels]
  [solid_x]
    type = ADStressDivergenceTensors
    variable = 'disp_x'
    component = 0
    use_displaced_mesh = true
  []
  [solid_y]
    type = ADStressDivergenceTensors
    variable = 'disp_y'
    component = 1
    use_displaced_mesh = true
  []
[]

[BCs]
  [fixed_x]
    type = DirichletBC
    variable = 'disp_x'
    boundary = 'left right'
    value = 0.0
  []
  [fixed_y]
    type = DirichletBC
    variable = 'disp_y'
    boundary = 'top bottom'
    value = 0.0
  []
[]

# [Constraints]
#   [disp_x_constraint]
#     type = XFEMSingleVariableConstraint
#     variable = 'disp_x'
#     geometric_cut_userobject = 'cut'
#     alpha = 1e10
#     use_penalty = true
#     use_displaced_mesh = false
#   []
#   [disp_y_constraint]
#     type = XFEMSingleVariableConstraint
#     variable = 'disp_y'
#     geometric_cut_userobject = 'cut'
#     alpha = 1e10
#     use_penalty = true
#     use_displaced_mesh = false
#   []
# []

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
    type = ADComputeSmallStrain
    base_name = oxide
    eigenstrain_names = 'eigenstrain_oxide'
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
  [strain_metal]
    type = ADComputeSmallStrain
    base_name = metal
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
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist'
  automatic_scaling = true
  line_search = 'none'

  # controls for nonlinear iterations
  nl_rel_tol = 1e-08
  nl_abs_tol = 1e-10

  # time control
  # dt = 8640.0
  dt = 1
  end_time = 1
  # end_time = 5184000 #60 days
  # end_time = 7776000 #90 days
  # end_time = 15552000 #180 days

  max_xfem_update = 2
[]

[Postprocessors]
  [disp_x]
    type = PointValue
    variable = 'disp_x'
    point = '1 0 0'
  []
[]

[Outputs]
  print_linear_converged_reason = false
  print_nonlinear_converged_reason = false
  print_linear_residuals = false
[]
