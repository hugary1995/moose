[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 2
    ny = 1
  []
  [pin]
    type = ParsedGenerateSideset
    input = gmg
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
    block = 0
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
    penalty_value = 10
  []
[]

[Kernels]
  [diffusion]
    type = Diffusion
    variable = u
    block = 0
  []
  [body_force]
    type = BodyForce
    variable = u
    function = '-2'
    block = 0
  []
[]

[UserObjects]
  [constraint]
    type = IntegralConstraint
    variable = u
    scalar_variable = h
    target = 2
    block = 0
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
  petsc_options_iname = '-pc_type -mat_view'
  petsc_options_value = 'lu       ::ascii_matlab'
  num_steps = 1
[]

[Outputs]
  exodus = true
[]
