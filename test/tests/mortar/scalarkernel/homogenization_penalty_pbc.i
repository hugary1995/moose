[GlobalParams]
  strain_constraint = false
[]

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 8
    ny = 8
    nz = 8
  []
  [matrix]
    type = SubdomainBoundingBoxGenerator
    input = gmg
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
  [pin]
    type = ParsedGenerateSideset
    input = particle
    new_sideset_name = pin
    combinatorial_geometry = 'x>0.499 & x<0.501'
  []
  [secondary]
    type = LowerDBlockFromSidesetGenerator
    input = pin
    sidesets = left
    new_block_id = 11
    new_block_name = secondary
  []
  [primary]
    type = LowerDBlockFromSidesetGenerator
    input = secondary
    sidesets = right
    new_block_id = 12
    new_block_name = primary
  []
[]

[Variables]
  [u]
    block = 'matrix particle'
  []
  [h]
    order = FIRST
    family = SCALAR
  []
[]

[BCs]
  [pin]
    type = DirichletBC
    variable = u
    value = 0
    boundary = pin
  []
[]

[Constraints]
  [mortar]
    type = PenaltyEqualValueConstraint
    secondary_variable = u
    primary_boundary = right
    secondary_boundary = left
    primary_subdomain = 12
    secondary_subdomain = 11
    penalty_value = 1e10
  []
[]

[Kernels]
  [diffusion]
    type = HomogenizedDiffusion
    variable = u
    scalar_variable = h
    diffusivity = D
    block = 'matrix particle'
  []
  [body_force]
    type = BodyForce
    variable = u
    function = '-2'
    block = 'matrix particle'
  []
[]

[Materials]
  [D_matrix]
    type = GenericConstantMaterial
    prop_names = 'D'
    prop_values = '1e5'
    block = matrix
  []
  [D_particle]
    type = GenericConstantMaterial
    prop_names = 'D'
    prop_values = '2e5'
    block = particle
  []
[]

[UserObjects]
  [constraint]
    type = IntegralConstraint
    variable = u
    scalar_variable = h
    target = 2
    diffusivity = D
    block = 'matrix particle'
    execute_on = 'INITIAL LINEAR NONLINEAR'
  []
[]

[ScalarKernels]
  [constraint]
    type = IntegralConstraintScalarKernel
    variable = h
    integrator = constraint
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options = '-snes_converged_reason -ksp_converged_reason -pc_svd_monitor'
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'svd      NONZERO               1e-15'
  num_steps = 1
[]

[Outputs]
  exodus = true
[]
