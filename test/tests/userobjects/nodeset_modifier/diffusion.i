[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 10
    ny = 10
  []
  [fix]
    type = BoundingBoxNodeSetGenerator
    input = gen
    new_boundary = fix
    bottom_left = '0 0 0'
    top_right = '0.5 0 0'
  []
[]

[Variables]
  [u]
  []
[]

[AuxVariables]
  [bnd_level_set]
  []
[]

[AuxKernels]
  [bnd_level_set]
    type = ParsedAux
    variable = bnd_level_set
    function = 'x-0.5+t'
    use_xyzt = true
  []
[]

[UserObjects]
  [modify_nodeset]
    type = VariableThresholdNodesetModifier
    variable = bnd_level_set
    boundary = fix
    from_boundary = fix
    to_boundary = bottom
    criterion_type = ABOVE
    threshold = 0
    execute_on = TIMESTEP_END
  []
[]

[Kernels]
  [diffusion]
    type = Diffusion
    variable = u
  []
[]

[BCs]
  [bottom_fix]
    type = DirichletBC
    variable = u
    boundary = fix
    value = 1
  []
  [top_fix]
    type = DirichletBC
    variable = u
    boundary = top
    value = 2
  []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  dt = 0.1
  end_time = 0.5

  nl_abs_tol = 1e-08
  nl_rel_tol = 1e-10
[]

[Outputs]
  exodus = true
[]
