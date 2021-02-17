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
    top_right = '0.1 1 1'
  []
  [right]
    type = SubdomainBoundingBoxGenerator
    input = 'left'
    block_id = 2
    bottom_left = '0.1 0 0'
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
    type = ChangeElementSubdomainByGeometricCut
    cut = cut
    origin_subdomain_id = 2
    target_subdomain_id = 1
    moving_boundary_name = 'moving_interface'
    include_domain_boundary = true
    reversible = true
    initialize_solution_on_move = false
    execute_on = 'INITIAL TIMESTEP_BEGIN'
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
[]

[Constraints]
  [u_continuity]
    type = XFEMEqualValueAtInterface
    variable = u
    geometric_cut_userobject = cut
    alpha = 1e6
    value = 1
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
  automatic_scaling = true

  dt = 0.01
  end_time = 0.8
  max_xfem_update = 1
[]

[Outputs]
  print_linear_converged_reason = false
  print_nonlinear_converged_reason = false
  print_linear_residuals = false
  exodus = true
[]
