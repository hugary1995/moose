n = 99
penalty = 10
v = 0.01

[XFEM]
  output_cut_plane = true
[]

[UserObjects]
  [level_set_cut_uo]
    type = LevelSetCutUserObject
    level_set_var = ls
    heal_always = true
  []
  [esm]
    type = CutElementSubdomainModifier
    geometric_cut_userobject = level_set_cut_uo
    apply_initial_conditions = false
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Mesh]
  [generated_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = ${n}
    ny = 1
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
    elem_type = QUAD4
  []
  [left]
    type = SubdomainBoundingBoxGenerator
    input = generated_mesh
    block_id = 1
    bottom_left = '0 0 0'
    top_right = '0.5 1 0'
  []
  [right]
    type = SubdomainBoundingBoxGenerator
    input = left
    block_id = 0
    bottom_left = '0.5 0 0'
    top_right = '1 1 0'
  []
[]

[Variables]
  [u]
  []
[]

[Functions]
  [u_mms]
    type = ParsedFunction
    value = '(0.25-(x-0.5)^2)*t'
  []
[]

[AuxVariables]
  [u_mms]
  []
  [ls]
  []
  [D]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [u_mms]
    type = FunctionAux
    variable = u_mms
    function = u_mms
  []
  [ls]
    type = FunctionAux
    variable = ls
    function = '${v}*t-x'
  []
  [D]
    type = MaterialRealAux
    variable = D
    property = D
    block = 1
    execute_on = 'TIMESTEP_END'
  []
[]

[Kernels]
  [diff]
    type = MatDiffusion
    variable = u
    diffusivity = D
  []
  [b1]
    type = BodyForce
    variable = u
    function = '(2*t+2-(4*x-1)/${v})*t'
    block = 1
  []
  [b0]
    type = BodyForce
    variable = u
    function = '2*t'
    block = 0
  []
[]

[Constraints]
  [u_continuity]
    type = XFEMSingleVariableConstraint
    variable = u
    alpha = ${penalty}
    geometric_cut_userobject = 'level_set_cut_uo'
    use_displaced_mesh = false
    use_penalty = true
  []
[]

[BCs]
  [zero]
    type = DirichletBC
    variable = u
    boundary = 'left right'
    value = 0
  []
[]

[Materials]
  [D1]
    type = RampedStatefulProperty
    property_name = D
    block = 1
  []
  [D0]
    type = GenericConstantMaterial
    prop_names = 'D'
    prop_values = '1'
    block = 0
  []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  automatic_scaling = true

  nl_rel_tol = 1e-06
  nl_abs_tol = 1e-08

  dt = '${fparse 100/102}'
  end_time = 100

  max_xfem_update = 1
[]

[Outputs]
  print_linear_residuals = false
  [exodus]
    type = Exodus
    file_base = 'mms_n_${n}/out'
  []
[]
