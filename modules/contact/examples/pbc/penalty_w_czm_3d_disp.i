strain_rate = 8.33e-5
strain_max = 0.5
nsteps_max = 200
nsteps_min = 20000

# Bilinear mixed mode model parameters for fracture in particle phase
p_penalty_stiffness = 1e6
p_GI_C = 473
p_GII_C = 397
p_normal_strength = 1630
p_shear_strength = 934
p_eta = 1.0

[Modules/TensorMechanics/CohesiveZoneMaster]
  [czm]
    strain = FINITE
    boundary = 'matrix_particle fracture'
  []
[]

[Materials]
  [czm_3dc_particle_particle]
    type = BiLinearMixedModeTraction
    penalty_stiffness = ${p_penalty_stiffness}
    GI_c = ${p_GI_C}
    GII_c = ${p_GII_C}
    normal_strength = ${p_normal_strength}
    shear_strength = ${p_shear_strength}
    eta = ${p_eta}
    mixed_mode_criterion = POWER_LAW
    lag_mode_mixity = true
    lag_displacement_jump = true
    viscosity = 1e-6
    boundary = 'matrix_particle fracture'
  []
[]

[GlobalParams]
  large_kinematics = true
  stabilize_strain = true
  displacements = 'disp_x disp_y disp_z'
[]

[Problem]
  kernel_coverage_check = false
[]

[Mesh]
  use_displaced_mesh = false
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
    top_right = '0.25 0.25 0.5'
    block_id = 1
    block_name = particle
  []
  [split]
    type = BreakMeshByBlockGenerator
    input = particle
    split_interface = true
  []
  [explode]
    type = ExplodeMeshGenerator
    input = split
    subdomains = 1
    interface_name = fracture
  []
  [secondary_left]
    type = LowerDBlockFromSidesetGenerator
    input = explode
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
  [fix_xyz]
    type = ExtraNodesetGenerator
    input = primary_top
    coord = '0.5 0.5 0.5'
    new_boundary = 'fix_xyz'
  []
  [fix_xy]
    type = ExtraNodesetGenerator
    input = fix_xyz
    coord = '0.75 0.75 0.75'
    new_boundary = 'fix_xy'
  []
  [fix_x]
    type = ExtraNodesetGenerator
    input = fix_xy
    coord = '0.5 0.5 0.75'
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
  []
  [disp_y]
  []
  [disp_z]
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
    type = TotalLagrangianStressDivergence
    variable = disp_x
    component = 0
    block = 'matrix particle'
  []
  [sdy]
    type = TotalLagrangianStressDivergence
    variable = disp_y
    component = 1
    block = 'matrix particle'
  []
  [sdz]
    type = TotalLagrangianStressDivergence
    variable = disp_z
    component = 2
    block = 'matrix particle'
  []
[]

[BCs]
  [x]
    type = DirichletBC
    value = 0.0
    boundary = 'back'
    variable = disp_x
  []
  [y]
    type = DirichletBC
    value = 0.0
    boundary = 'back'
    variable = disp_y
  []
  [z]
    type = DirichletBC
    value = 0.0
    boundary = 'back'
    variable = disp_z
  []
  [disp]
    type = FunctionDirichletBC
    function = stretch
    boundary = 'front'
    variable = disp_z
  []
[]

[Constraints]
  [ev_xx]
    type = PenaltyEqualValueConstraint
    primary_boundary = right
    secondary_boundary = left
    primary_subdomain = primary_right
    secondary_subdomain = secondary_left
    secondary_variable = disp_x
    penalty_value = 1e10
  []
  [ev_xy]
    type = PenaltyEqualValueConstraint
    primary_boundary = top
    secondary_boundary = bottom
    primary_subdomain = primary_top
    secondary_subdomain = secondary_bottom
    secondary_variable = disp_x
    penalty_value = 1e10
  []
  [ev_yx]
    type = PenaltyEqualValueConstraint
    primary_boundary = right
    secondary_boundary = left
    primary_subdomain = primary_right
    secondary_subdomain = secondary_left
    secondary_variable = disp_y
    penalty_value = 1e10
  []
  [ev_yy]
    type = PenaltyEqualValueConstraint
    primary_boundary = top
    secondary_boundary = bottom
    primary_subdomain = primary_top
    secondary_subdomain = secondary_bottom
    secondary_variable = disp_y
    penalty_value = 1e10
  []
  [ev_zx]
    type = PenaltyEqualValueConstraint
    primary_boundary = right
    secondary_boundary = left
    primary_subdomain = primary_right
    secondary_subdomain = secondary_left
    secondary_variable = disp_z
    penalty_value = 1e10
  []
  [ev_zy]
    type = PenaltyEqualValueConstraint
    primary_boundary = top
    secondary_boundary = bottom
    primary_subdomain = primary_top
    secondary_subdomain = secondary_bottom
    secondary_variable = disp_z
    penalty_value = 1e10
  []
[]

[Materials]
  [strain]
    type = ComputeLagrangianStrain
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
  # petsc_options_value = 'hypre boomeramg 200 0.7 ext+i PMIS 4 2 0.4 1e-10'

  # petsc_options_iname = '-pc_type -ksp_grmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap '
  #                       '-sub_pc_factor_shift_type -sub_pc_factor_shift_amout'
  # petsc_options_value = 'asm      31                  preonly       lu           2               '
  #                       'NONZERO                   1e-10'

  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amout'
  petsc_options_value = 'svd      NONZERO               1e-10'

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
