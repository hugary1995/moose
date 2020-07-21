
[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  type = PeridynamicsMesh
  horizon_number = 3
  cracks_start = '0.25 0.5 0'
  cracks_end = '0.75 0.5 0'

  [./gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 8
    ny = 8
  [../]
  [./gpd]
    type = MeshGeneratorPD
    input = gmg
    retain_fe_mesh = false
  [../]
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
[]

[AuxVariables]
  [./critical_stretch]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]

[AuxKernels]
  [./bond_status]
    type = StretchBasedFailureCriterionPD
    critical_variable = critical_stretch
    variable = bond_status
  [../]
[]

[ICs]
  [./critical_stretch]
    type = ConstantIC
    variable = critical_stretch
    value = 0.001
  [../]
[]

[BCs]
  [./left_x]
    type = DirichletBC
    variable = disp_x
    boundary = 1003
    value = 0.0
  [../]
  [./top_y]
    type = DirichletBC
    variable = disp_y
    boundary = 1002
    value = 0.0
  [../]
  [./bottom_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 1000
    function = '-0.001*t'
  [../]

  [./rbm_x]
    type = RBMPresetOldValuePD
    variable = disp_x
    boundary = 999
  [../]
  [./rbm_y]
    type = RBMPresetOldValuePD
    variable = disp_y
    boundary = 999
  [../]
[]

[Modules/Peridynamics/Mechanics/Master]
  [./all]
    formulation = BOND
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 2e5
    poissons_ratio = 0.33
  [../]
  [./force_density]
    type = ComputeSmallStrainConstantHorizonMaterialBPD
  [../]
[]

[Postprocessors]
  [./bond_status_updated_times]
    type = BondStatusConvergedPostprocessorPD
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  line_search = none
  start_time = 0
  dt = 0.5
  end_time = 1

  picard_max_its = 5
  accept_on_max_picard_iteration = true
  picard_custom_pp = bond_status_updated_times
  custom_abs_tol = 2
  disable_picard_residual_norm_check = true
[]

[Outputs]
  file_base = 2D_bond_status_convergence_BPD
  exodus = true
[]
