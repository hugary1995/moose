strain_rate = 8.33e-5
strain_max = 0.5
nsteps_max = 200
nsteps_min = 20000

[GlobalParams]
  large_kinematics = true
  constraint_types = 'stress strain strain stress stress strain stress stress strain'
  macro_gradient = hvar
  stabilize_strain = true
[]

[Problem]
  #kernel_coverage_check = true
[]

[Mesh]
  [msh]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 4
    ny = 4
    nz = 4
  []
  [matrix]
    type = SubdomainBoundingBoxGenerator
    input = msh
    bottom_left = '0 0 0'
    top_right = '1 1 1'
    block_id = 0
    block_name = matrix
  []
  [particle]
    type = SubdomainBoundingBoxGenerator
    input = matrix
    bottom_left = '0 0 0'
    top_right = '0.5 0.5 0.5'
    block_id = 1
    block_name = particle
  []
  [secondary_left]
    type = LowerDBlockFromSidesetGenerator
    input = particle
    sidesets = 'left'
    new_block_id = 11
    new_block_name = 'secondary_left'
  []
  [primary_right]
    type = LowerDBlockFromSidesetGenerator
    input = secondary_left
    sidesets = 'right'
    new_block_id = '10'
    new_block_name = 'primary_right'
  []
  [secondary_bottom]
    type = LowerDBlockFromSidesetGenerator
    input = primary_right
    sidesets = 'bottom'
    new_block_id = '21'
    new_block_name = 'secondary_bottom'
  []
  [primary_top]
    type = LowerDBlockFromSidesetGenerator
    input = secondary_bottom
    sidesets = 'top'
    new_block_id = '20'
    new_block_name = 'primary_top'
  []
  [secondary_back]
    type = LowerDBlockFromSidesetGenerator
    input = primary_top
    sidesets = 'back'
    new_block_id = '31'
    new_block_name = 'secondary_back'
  []
  [primary_front]
    type = LowerDBlockFromSidesetGenerator
    input = secondary_back
    sidesets = 'front'
    new_block_id = '30'
    new_block_name = 'primary_front'
  []
  [fix_all]
    type = ExtraNodesetGenerator
    input = primary_front
    coord = '1 0 0'
    new_boundary = 'fix_all'
  []
  [fix_xy]
    type = ExtraNodesetGenerator
    input = fix_all
    coord = '0 0 0'
    new_boundary = 'fix_xy'
  []
  [fix_x]
    type = ExtraNodesetGenerator
    input = fix_xy
    coord = '1 1 0'
    new_boundary = 'fix_x'
  []

[]

[Functions]
  [stretch]
    type = PiecewiseLinear
    x = '0 ${fparse strain_max / strain_rate}'
    y = '0 ${fparse strain_max}'
  []
[]

[Variables]
  [disp_x]
    block = 'matrix particle'
  []
  [disp_y]
    block = 'matrix particle'
  []
  [disp_z]
    block = 'matrix particle'
  []
  [hvar]
    family = SCALAR
    order = NINTH
  []
[]

[AuxVariables]
  [stress_00]
    order = CONSTANT
    family = MONOMIAL
    block = 'matrix particle'
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = pk1_stress
      index_i = 0
      index_j = 0
      block = 'matrix particle'
    []
  []
[]

[Kernels]
  [sdx]
    type = HomogenizedTotalLagrangianStressDivergence
    variable = disp_x
    component = 0
    displacements = 'disp_x disp_y disp_z'
    block = 'matrix particle'
  []
  [sdy]
    type = HomogenizedTotalLagrangianStressDivergence
    variable = disp_y
    component = 1
    displacements = 'disp_x disp_y disp_z'
    block = 'matrix particle'
  []
  [sdz]
    type = HomogenizedTotalLagrangianStressDivergence
    variable = disp_z
    component = 2
    displacements = 'disp_x disp_y disp_z'
    block = 'matrix particle'
  []
[]

[ScalarKernels]
  [enforce]
    type = HomogenizationConstraintScalarKernel
    variable = hvar
    integrator = integrator
    ndim = 3
  []
[]

[UserObjects]
  [integrator]
    type = HomogenizationConstraintIntegral
    targets = '0 0 0 0 0 0 0 0 stretch'
    displacements = 'disp_x disp_y disp_z'
    block = 'particle matrix'
    execute_on = 'INITIAL LINEAR NONLINEAR'
  []
[]

[BCs]
  [fix_all_x]
    type = DirichletBC
    value = 0.0
    boundary = fix_all
    variable = disp_x
  []
  [fix_all_y]
    type = DirichletBC
    value = 0.0
    boundary = fix_all
    variable = disp_y
  []
  [fix_all_z]
    type = DirichletBC
    value = 0.0
    boundary = fix_all
    variable = disp_z
  []
  [fix_xy_x]
    type = DirichletBC
    value = 0.0
    boundary = fix_xy
    variable = disp_x
  []
  [fix_xy_y]
    type = DirichletBC
    value = 0.0
    boundary = fix_xy
    variable = disp_y
  []
  [fix_x_x]
    type = DirichletBC
    value = 0.0
    boundary = fix_x
    variable = disp_x
  []
[]

[Constraints]
  [ev_x]
    type = PenaltyEqualValueConstraint
    primary_boundary = right
    secondary_boundary = left
    primary_subdomain = '10'
    secondary_subdomain = '11'
    secondary_variable = disp_x
    penalty_value = 1.0e10
    quadrature = SECOND
  []
  [ev_y]
    type = PenaltyEqualValueConstraint
    primary_boundary = top
    secondary_boundary = bottom
    primary_subdomain = '20'
    secondary_subdomain = '21'
    secondary_variable = disp_y
    penalty_value = 1.0e10
    quadrature = SECOND
  []
  [ev_z]
    type = PenaltyEqualValueConstraint
    primary_boundary = front
    secondary_boundary = back
    primary_subdomain = '30'
    secondary_subdomain = '31'
    secondary_variable = disp_z
    penalty_value = 1.0e10
    quadrature = SECOND
  []
[]

[Materials]
  [homogenization_gradient]
    type = ComputeHomogenizedLagrangianStrain
    displacements = 'disp_x disp_y disp_z'
    block = 'matrix particle'
  []
  [strain]
    type = ComputeLagrangianStrain
    displacements = 'disp_x disp_y disp_z'
    homogenization_gradient_names = 'homogenization_gradient'
    block = 'matrix particle'
  []
  [stress]
    type = ComputeLagrangianLinearElasticStress
    block = 'matrix particle'
  []
  [C1]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 2e5
    poissons_ratio = 0.3
    block = 'matrix'
  []
  [C2]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 9e5
    poissons_ratio = 0.2
    block = 'particle'
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options = '-snes_converged_reason -ksp_converged_reason -pc_svd_monitor'
  # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -ksp_gmres_restart '
  #                       '-pc_hypre_boomeramg_strong_threshold -pc_hypre_boomeramg_interp_type '
  #                       '-pc_hypre_boomeramg_coarsen_type -pc_hypre_boomeramg_agg_nl '
  #                       '-pc_hypre_boomeramg_agg_num_paths -pc_hypre_boomeramg_truncfactor '
  #                       '-pc_factor_shift_amount'
  # petsc_options_value = 'hypre boomeramg 200 0.7 ext+i PMIS 4 2 0.4 1e-15'

  petsc_options_iname = '-pc_type -ksp_grmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  petsc_options_value = 'asm      31                  preonly       lu           2'

  # automatic_scaling = true
  l_max_its = 150
  l_tol = 1e-6
  nl_max_its = 100
  nl_rel_tol = 1e-4
  nl_abs_tol = 1e-6
  nl_forced_its = 1

  line_search = 'basic'
  end_time = '${fparse strain_max / strain_rate}'
  dtmin = '${fparse strain_max / strain_rate / nsteps_min}'
  dtmax = '${fparse strain_max / strain_rate / nsteps_max}'
  num_steps = 10
[]

[Outputs]
  exodus = true
  print_linear_residuals = false
[]

[Postprocessors]
  [stress_00]
    type = ElementAverageValue
    variable = stress_00
    block = 'matrix particle'
  []
[]
