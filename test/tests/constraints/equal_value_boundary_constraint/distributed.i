[Problem]
  error_on_jacobian_nonzero_reallocation = true
[]

[Mesh]
  [base]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 10
    ny = 10
    nz = 10
  []
[]

[Variables]
  [u]
  []
[]

[Kernels]
  [diff_u]
    type = Diffusion
    variable = u
  []
[]

[BCs]
  [left_u]
    type = DirichletBC
    variable = u
    boundary = 'left'
    value = 0
  []
  [right_u]
    type = DirichletBC
    variable = u
    boundary = 'right'
    value = 1
  []
[]

[Constraints]
  [u]
    type = EqualValueBoundaryConstraint
    variable = u
    secondary = 'top'
    penalty = 1e3
  []
[]

[Executioner]
  type = Transient
  solve_type = 'newton'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  num_steps = 1
[]

