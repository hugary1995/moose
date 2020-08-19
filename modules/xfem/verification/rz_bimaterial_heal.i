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

[GlobalParams]
  displacements = 'disp_r disp_z'
  decomposition_method = EigenSolution
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
    function = 'x-2'
  []
  [temperature]
    type = FunctionAux
    variable = temperature
    function = '${T_ref}+t'
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
    use_displaced_mesh = false
  []
  [disp_z_constraint]
    type = XFEMSingleVariableConstraint
    variable = 'disp_z'
    geometric_cut_userobject = 'cut'
    alpha = 1e10
    use_penalty = true
    use_displaced_mesh = false
  []
[]

[BCs]
  [fixed_r]
    type = DirichletBC
    variable = 'disp_r'
    boundary = 'right'
    value = 0.0
  []
  [fixed_z]
    type = DirichletBC
    variable = 'disp_z'
    boundary = 'top bottom'
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
    type = ADComputeAxisymmetricRZFiniteStrain
    base_name = oxide
    eigenstrain_names = 'eigenstrain_oxide'
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
  [strain_metal]
    type = ADComputeAxisymmetricRZFiniteStrain
    base_name = metal
  []
  [stress_metal]
    type = ADComputeFiniteStrainElasticStress
    base_name = metal
  []
  # bimaterial
  [combined_stress]
    type = ADLevelSetMultiRankTwoTensorMaterial
    level_set_vars = 'interface'
    base_name_keys = '+     -'
    base_name_vals = 'metal oxide'
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
  dt = 1
  end_time = 30

  max_xfem_update = 1
[]

[Postprocessors]
  [disp_r]
    type = PointValue
    variable = 'disp_r'
    point = '1 0 0'
  []
[]

[Outputs]
  print_linear_converged_reason = false
  print_nonlinear_converged_reason = false
  print_linear_residuals = false
  [exodus]
    type = Exodus
    file_base = 'output/visualize_rz_bimaterial_heal'
  []
[]
