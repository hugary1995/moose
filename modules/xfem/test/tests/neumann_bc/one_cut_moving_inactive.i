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
  [cut]
    type = LevelSetCutUserObject
    level_set_var = phi
    negative_id = 1
    positive_id = 33
    heal_always = true
  []
  [esm]
    type = CutElementSubdomainModifier
    geometric_cut_userobject = cut
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
  [A]
    type = SubdomainBoundingBoxGenerator
    input = square
    block_id = 1
    bottom_left = '0 0 0'
    top_right = '0.5 1 0'
  []
  [B]
    type = SubdomainBoundingBoxGenerator
    input = A
    block_id = 33
    bottom_left = '0.5 0 0'
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
    block = 1
  []
[]

[AuxVariables]
  [phi]
  []
  [solution]
  []
[]

[AuxKernels]
  [phi]
    type = FunctionAux
    variable = phi
    function = 'x-0.213-t'
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
    block = 1
  []
  [b]
    type = BodyForce
    variable = u
    function = b
    block = 1
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
    geometric_cut_userobject = cut
    function = flux
    block = 1
  []
[]

[Materials]
  [diffusivity_B]
    type = GenericConstantMaterial
    prop_names = 'D'
    prop_values = '1'
    block = 1
  []
[]

[Postprocessors]
  [error]
    type = ElementL2Error
    variable = u
    function = solution
    block = 1
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
