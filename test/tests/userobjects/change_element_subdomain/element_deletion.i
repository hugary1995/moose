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
    top_right = '0.9 1 1'
  []
  [right]
    type = SubdomainBoundingBoxGenerator
    input = 'left'
    block_id = 2
    bottom_left = '0.9 0 0'
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
  [moving_interface]
    type = ChangeElementSubdomainByCoupledVariable
    coupled_var = 'phi'
    activate_value = 0
    origin_subdomain_id = 1
    target_subdomain_id = 2
    moving_boundary_name = 'moving_interface'
    include_domain_boundary = false
    reversible = false
    initialize_solution_on_move = false
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
[]

[Variables]
  [u]
  []
[]

[Functions]
  [wavy_interface]
    type = ParsedFunction
    value = 'if(t<0.8, x-0.1-t+0.05*sin(5*pi*(y+t))+0.025*cos(8*pi*(y+2*t)), '
            'x-0.1+t-1.6+0.05*sin(5*pi*(y+t))+0.025*cos(8*pi*(y+2*t)))'
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
  [interface]
    type = DirichletBC
    variable = u
    boundary = 'moving_interface'
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

  dt = 0.01
  start_time = 0.8
  end_time = 1.6
[]

[Outputs]
  print_linear_converged_reason = false
  print_nonlinear_converged_reason = false
  print_linear_residuals = false
  exodus = true
[]
