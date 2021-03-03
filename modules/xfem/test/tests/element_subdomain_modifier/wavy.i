[XFEM]
  output_cut_plane = true
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 40
    ny = 40
  []
  [left]
    type = SubdomainBoundingBoxGenerator
    input = 'gen'
    block_id = 1
    bottom_left = '0 0 0'
    top_right = '0.5 1 1'
  []
  [right]
    type = SubdomainBoundingBoxGenerator
    input = 'left'
    block_id = 2
    bottom_left = '0.5 0 0'
    top_right = '1 1 1'
  []
  [sidesets]
    type = SideSetsAroundSubdomainGenerator
    input = right
    normal = '1 0 0'
    block = 1
    new_boundary = 'moving_interface'
  []
[]

[UserObjects]
  [cut]
    type = LevelSetCutUserObject
    level_set_var = phi
    negative_id = 1
    positive_id = 2
    heal_always = true
  []
  [moving_interface]
    type = GeometricCutElementSubdomainModifier
    cut = cut
    moving_boundary_name = moving_interface
    apply_initial_conditions = false
    execute_on = 'TIMESTEP_BEGIN TIMESTEP_END'
  []
[]

[Functions]
  [wavy_interface]
    type = ParsedFunction
    value = 'if(t<0.8, x-0.113-t+0.05*sin(5*pi*(y+t))+0.025*cos(8*pi*(y+2*t)), '
            'x-0.113+t-1.6+0.05*sin(5*pi*(y+t))+0.025*cos(8*pi*(y+2*t)))'
  []
[]

[Variables]
  [u]
  []
[]

[AuxVariables]
  [phi]
    [InitialCondition]
      type = FunctionIC
      function = wavy_interface
    []
  []
[]

[AuxKernels]
  [phi]
    type = FunctionAux
    variable = phi
    function = wavy_interface
    execute_on = 'TIMESTEP_BEGIN'
  []
[]

[Kernels]
  [diff]
    type = MatDiffusion
    diffusivity = 'diffusivity'
    variable = 'u'
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = u
    boundary = 'left'
    value = 0
  []
  [right]
    type = DirichletBC
    variable = u
    boundary = 'right'
    value = 1
  []
[]

[Constraints]
  [u_continuity]
    type = XFEMSingleVariableConstraint
    variable = u
    geometric_cut_userobject = cut
    alpha = 100
    use_penalty = true
  []
[]

[Materials]
  [diffusivity1]
    type = GenericConstantMaterial
    prop_names = 'diffusivity'
    prop_values = '4'
    block = 1
    outputs = exodus
  []
  [diffusivity2]
    type = GenericConstantMaterial
    prop_names = 'diffusivity'
    prop_values = '1'
    block = 2
    outputs = exodus
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-06
  nl_abs_tol = 1e-08

  dt = 0.01
  end_time = 1.6
  max_xfem_update = 1
[]

[Outputs]
  print_linear_converged_reason = false
  print_nonlinear_converged_reason = false
  print_linear_residuals = false
  exodus = true
[]
