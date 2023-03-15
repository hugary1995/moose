[Problem]
  error_on_jacobian_nonzero_reallocation = true
[]

[Mesh]
  [base]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 4
    ny = 4
    nz = 4
  []
[]

[Variables]
  [u]
  []
[]

[Kernels]
  [diff]
    type = Diffusion
    variable = u
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = u
    boundary = left
    value = 0
  []
  [right]
    type = FunctionNeumannBC
    variable = u
    boundary = right
    function = y
  []
[]

[Constraints]
  [u]
    type = EqualValueBoundaryConstraint
    variable = u
    secondary = right
    penalty = 1000
  []
[]

[Executioner]
  type = Transient
  solve_type = 'newton'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  num_steps = 1
  automatic_scaling = false
  line_search = none
[]

[Outputs]
  exodus = true
[]

