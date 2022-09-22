a = 250

E = 6.88e4
nu = 0.33

sigma_y = 320
n = 5
ep0 = 0.01

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  stabilize_strain = true
  large_kinematics = true
[]

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 3
    xmax = ${a}
    ymax = ${a}
    zmax = ${a}
  []
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [disp_z]
  []
[]

[AuxVariables]
  [rx]
  []
  [ry]
  []
  [rz]
  []
[]

[Kernels]
  [solid_x]
    type = TotalLagrangianStressDivergence
    variable = disp_x
    component = 0
    save_in = rx
  []
  [solid_y]
    type = TotalLagrangianStressDivergence
    variable = disp_y
    component = 1
    save_in = ry
  []
  [solid_z]
    type = TotalLagrangianStressDivergence
    variable = disp_z
    component = 2
    save_in = rz
  []
[]

[Materials]
  [elastic_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = ${E}
    poissons_ratio = ${nu}
  []
  [compute_strain]
    type = ComputeLagrangianStrain
  []
  [flow_stress]
    type = DerivativeParsedMaterial
    f_name = flow_stress
    function = '${sigma_y}*(1+effective_plastic_strain/${ep0})^(1/${n})'
    material_property_names = 'effective_plastic_strain'
    additional_derivative_symbols = 'effective_plastic_strain'
    derivative_order = 2
    compute = false
  []
  [compute_stress]
    type = ComputeSimoHughesJ2PlasticityStress
    flow_stress_material = flow_stress
  []
  [psip]
    type = ParsedMaterial
    f_name = plastic_energy
    function = '${n}*${sigma_y}*${ep0}/(${n}+1)*((1+effective_plastic_strain/${ep0})^(1/${n}+1)-1)'
    material_property_names = 'effective_plastic_strain'
  []
[]

[BCs]
  [xfix]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0
  []
  [yfix]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
  []
  [zfix]
    type = DirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0
  []
  [ydisp]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 'top'
    function = 't'
    preset = false
  []
[]

[Postprocessors]
  [ep]
    type = ElementAverageMaterialProperty
    mat_prop = effective_plastic_strain
  []
  [W]
    type = ExternalWork
    forces = 'rx ry rz'
    displacements = 'disp_x disp_y disp_z'
  []
  [psie]
    type = ElementIntegralMaterialProperty
    mat_prop = elastic_energy
  []
  [psip]
    type = ElementIntegralMaterialProperty
    mat_prop = plastic_energy
  []
  [psi]
    type = ParsedPostprocessor
    function = 'psie+psip'
    pp_names = 'psie psip'
  []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu      '

  line_search = none

  nl_rel_tol = 1e-08
  nl_abs_tol = 1e-10
  nl_max_its = 50

  dt = '${fparse 0.0001 * a}'
  end_time = '${fparse 0.1 * a}'

  automatic_scaling = true

  abort_on_solve_fail = true
[]

[Outputs]
  file_base = stress_deformation
  print_linear_residuals = false
  csv = true
  exodus = true
[]
