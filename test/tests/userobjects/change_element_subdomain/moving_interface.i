[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 100
    ny = 100
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
  [moving_interface]
    type = ChangeElementSubdomainByCoupledVariable
    coupled_var = 'phi'
    activate_value = 0
    origin_subdomain_id = 1
    target_subdomain_id = 2
    moving_boundary_name = 'moving_interface'
    include_domain_boundary = false
    reversible = true
    initialize_solution_on_move = false
    execute_on = 'TIMESTEP_BEGIN'
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
  [right]
    type = DirichletBC
    variable = u
    boundary = 'right'
    value = 1
  []
[]

[Materials]
  [stateful1]
    type = StatefulMaterial
    outputs = exodus
    initial_diffusivity = 1
    multiplier = 1.01
    block = 1
  []
  [stateful2]
    type = StatefulMaterial
    outputs = exodus
    initial_diffusivity = 1
    multiplier = 1.01
    block = 2
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

  # num_steps = 5
[]

[Outputs]
  print_linear_converged_reason = false
  print_nonlinear_converged_reason = false
  print_linear_residuals = false
  exodus = true
[]
