nr = 10
nt = 20
a = 13.8
b = 14
c = 25
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
  displacements = 'disp_x disp_y'
[]

[Mesh]
  [oxide]
    type = AnnularMeshGenerator
    dmin = 0
    dmax = 90
    nr = ${nr}
    nt = ${nt}
    rmin = ${a}
    rmax = 14.2
  []
  [metal]
    type = AnnularMeshGenerator
    dmin = 0
    dmax = 90
    nr = ${nr}
    nt = ${nt}
    rmin = 14.2
    rmax = ${c}
  []
  [stitch]
    type = StitchedMeshGenerator
    inputs = 'oxide metal'
    stitch_boundaries_pairs = 'rmax rmin'
  []
  # [amg]
  #   type = AnnularMeshGenerator
  #   dmin = 0
  #   dmax = 90
  #   nr = ${nr}
  #   nt = ${nt}
  #   rmin = ${a}
  #   rmax = '${c}'
  # []
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

[Variables]
  [disp_x]
  []
  [disp_y]
  []
[]

[AuxVariables]
  [T]
  []
  [interface]
  []
  [stress_xx]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_yy]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [T]
    type = FunctionAux
    variable = 'T'
    function = '${T_ref} + t * (${T_fin} - ${T_ref}) / ${end_time}'
  []
  [interface]
    type = FunctionAux
    variable = interface
    function = 'sqrt(x^2+y^2)-${b}'
    # function = 'if(t<10,sqrt(x^2+y^2)-2,sqrt(x^2+y^2)-2-0.1001*(t-10))'
  []
  [stress_xx]
    type = ADRankTwoAux
    variable = 'stress_xx'
    rank_two_tensor = 'stress'
    index_i = 0
    index_j = 0
  []
  [stress_yy]
    type = ADRankTwoAux
    variable = 'stress_yy'
    rank_two_tensor = 'stress'
    index_i = 1
    index_j = 1
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

[Constraints]
  [disp_x_constraint]
    type = XFEMSingleVariableConstraint
    variable = 'disp_x'
    geometric_cut_userobject = 'cut'
    alpha = 1e10
    use_penalty = true
    use_displaced_mesh = true
  []
  [disp_y_constraint]
    type = XFEMSingleVariableConstraint
    variable = 'disp_y'
    geometric_cut_userobject = 'cut'
    alpha = 1e10
    use_penalty = true
    use_displaced_mesh = true
  []
[]

[BCs]
  [symmetry_y]
    type = DirichletBC
    variable = 'disp_y'
    value = 0
    boundary = 'dmin'
  []
  [symmetry_x]
    type = DirichletBC
    variable = 'disp_x'
    value = 0
    boundary = 'dmax'
  []
[]

[Materials]
  [eigenstrain_oxide]
    type = ADComputeThermalExpansionEigenstrain
    eigenstrain_name = 'thermal_eigenstrain'
    stress_free_temperature = ${T_ref}
    temperature = 'T'
    thermal_expansion_coeff = ${CTE_oxide}
    base_name = 'oxide'
  []
  [elasticity_tensor_oxide]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = ${E_oxide}
    poissons_ratio = ${nu_oxide}
    base_name = oxide
  []
  [strain_oxide]
    type = ADComputeFiniteStrain
    eigenstrain_names = 'thermal_eigenstrain'
    base_name = oxide
  []
  [stress_oxide]
    type = ADComputeFiniteStrainElasticStress
    base_name = oxide
  []
  [eigenstrain_metal]
    type = ADComputeThermalExpansionEigenstrain
    eigenstrain_name = 'thermal_eigenstrain'
    stress_free_temperature = ${T_ref}
    temperature = 'T'
    thermal_expansion_coeff = ${CTE_metal}
    base_name = 'metal'
  []
  [elasticity_tensor_metal]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = ${E_metal}
    poissons_ratio = ${nu_metal}
    base_name = metal
  []
  [strain_metal]
    type = ADComputeFiniteStrain
    eigenstrain_names = 'thermal_eigenstrain'
    base_name = metal
  []
  [stress_metal]
    type = ADComputeFiniteStrainElasticStress
    base_name = metal
  []
  # bimaterial
  [combined_stress]
    type = ADLevelSetBiMaterialRankTwo
    levelset_negative_base = 'oxide'
    levelset_positive_base = 'metal'
    level_set_var = 'interface'
    prop_name = stress
  []
[]

[VectorPostprocessors]
  [stress_rr]
    type = LineValueSampler
    variable = 'stress_xx'
    num_points = '${fparse nr*5}'
    start_point = '${a} 0 0'
    end_point = '${c} 0 0'
    sort_by = x
  []
  [stress_tt]
    type = LineValueSampler
    variable = 'stress_yy'
    num_points = '${fparse nr*5}'
    start_point = '${a} 0 0'
    end_point = '${c} 0 0'
    sort_by = x
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist'
  automatic_scaling = true

  # controls for nonlinear iterations
  nl_rel_tol = 1e-06
  nl_abs_tol = 1e-08

  dt = 1
  end_time = ${end_time}

  max_xfem_update = 1
[]

[Outputs]
  print_linear_converged_reason = false
  print_nonlinear_converged_reason = false
  print_linear_residuals = false
  [exodus]
    type = Exodus
    file_base = 'output/visualize_nr_${nr}_nt_${nt}'
    execute_on = 'FINAL'
  []
  [csv]
    type = CSV
    file_base = 'output/data_nr_${nr}_nt_${nt}'
    execute_on = 'FINAL'
  []
[]
