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
[]

[Variables]
  [u]
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
  [Periodic]
    [x]
      variable = u
      auto_direction = x
    []
  []
[]

[Kernels]
  [diffusion]
    type = Diffusion
    variable = u
  []
  [body_force]
    type = BodyForce
    variable = u
    function = '-2'
  []
[]

[UserObjects]
  [constraint]
    type = IntegralConstraint
    variable = u
    scalar_variable = h
    target = 2
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
