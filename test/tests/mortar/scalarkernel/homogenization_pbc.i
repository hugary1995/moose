[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 4
    ny = 4
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
    top_right = '0.5 0.5 1'
    block_id = 1
    block_name = particle
  []
  [pin]
    type = ParsedGenerateSideset
    input = particle
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
    type = MatDiffusion
    variable = u
    diffusivity = D
  []
  [body_force]
    type = BodyForce
    variable = u
    function = '-2'
  []
[]

[Materials]
  [D_matrix]
    type = GenericConstantMaterial
    prop_names = 'D'
    prop_values = '1'
    block = matrix
  []
  [D_particle]
    type = GenericConstantMaterial
    prop_names = 'D'
    prop_values = '2'
    block = particle
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
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu      '
  num_steps = 1
[]

[Outputs]
  exodus = true
[]
