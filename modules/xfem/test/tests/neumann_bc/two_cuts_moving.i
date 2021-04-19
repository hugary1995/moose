N = 4

[Problem]
  kernel_coverage_check = false
  material_coverage_check = false
[]

[XFEM]
  qrule = volfrac
  output_cut_plane = true
[]

[UserObjects]
  [cut1]
    type = LevelSetCutUserObject
    level_set_var = phi1
    negative_id = 1
    positive_id = 33
    heal_always = true
  []
  [cut2]
    type = LevelSetCutUserObject
    level_set_var = phi2
    negative_id = 5
    positive_id = 7
    heal_always = true
  []
  [combo]
    type = ComboCutUserObject
    geometric_cut_userobjects = 'cut1 cut2'
    cut_subdomain_combinations = '1 5 33 5 33 7'
    cut_subdomains = '1 2 3'
    heal_always = false
    execute_on = NONE
  []
  [esm]
    type = CutElementSubdomainModifier
    geometric_cut_userobject = combo
    apply_initial_conditions = false
  []
[]

[Mesh]
  [square]
    type = GeneratedMeshGenerator
    dim = 2
    nx = ${N}
    ny = 1
  []
  [left]
    type = SubdomainBoundingBoxGenerator
    input = square
    block_id = 1
    bottom_left = '0 0 0'
    top_right = '0.3 1 0'
  []
  [middle]
    type = SubdomainBoundingBoxGenerator
    input = left
    block_id = 2
    bottom_left = '0.3 0 0'
    top_right = '0.6 1 0'
  []
  [right]
    type = SubdomainBoundingBoxGenerator
    input = middle
    block_id = 3
    bottom_left = '0.6 0 0'
    top_right = '1 1 0'
  []
[]

[Functions]
  [solution]
    type = ParsedFunction
    value = '3*x^3+2*x^2+5*x-10'
  []
  [flux]
    type = ParsedFunction
    value = '9*x^2+4*x+5'
  []
  [b]
    type = ParsedFunction
    value = '-18*x-4'
  []
[]

[Variables]
  [u]
    block = '1 2'
  []
[]

[AuxVariables]
  [phi1]
  []
  [phi2]
  []
  [solution]
  []
[]

[AuxKernels]
  [phi1]
    type = FunctionAux
    variable = phi1
    function = 'x-0.313'
  []
  [phi2]
    type = FunctionAux
    variable = phi2
    function = 'x-0.713-t'
  []
  [solution]
    type = FunctionAux
    variable = solution
    function = solution
  []
[]

[Kernels]
  [diff]
    type = MatDiffusion
    variable = u
    diffusivity = D
    block = '1 2'
  []
  [b]
    type = BodyForce
    variable = u
    function = b
    block = '1 2'
  []
[]

[BCs]
  [u_fix]
    type = FunctionDirichletBC
    variable = u
    boundary = 'left'
    function = 'solution'
  []
[]

[DiracKernels]
  [u_flux]
    type = XFEMFunctionNeumannBC
    variable = u
    geometric_cut_userobject = cut2
    function = flux
    block = 2
  []
[]

[Constraints]
  [continuity]
    type = XFEMSingleVariableConstraint
    variable = u
    geometric_cut_userobject = cut1
    use_penalty = true
  []
[]

[Materials]
  [diffusivity]
    type = GenericConstantMaterial
    prop_names = 'D'
    prop_values = '1'
    block = '1 2'
  []
[]

[Postprocessors]
  [error]
    type = ElementL2Error
    variable = u
    function = solution
    block = '1 2'
  []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  dt = 0.1
  num_steps = 2
  nl_rel_tol = 1e-08
  nl_abs_tol = 1e-10

  max_xfem_update = 1
  abort_on_solve_fail = true
[]

[Outputs]
  exodus = true
[]
