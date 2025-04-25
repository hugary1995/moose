!include 'expdyn.i'

[NEML2]
  input = '../models/elasticity.i'
  [all]
    executor_name = 'linear_elasticity'
    model = 'model'
    verbose = true
    moose_input_kernels = 'strain'
    # To run on GPU, pass --libtorch-device=cuda on the command line
    # and uncomment the following lines:
    # device = 'cuda'
    # output_device = 'cuda'
  []
[]

[UserObjects]
  [assembly]
    type = NEML2Assembly
  []
  [fem]
    type = NEML2FEMInterpolation
    assembly = 'assembly'
  []
  [strain]
    type = NEML2SmallStrain
    assembly = 'assembly'
    fem = 'fem'
    to_neml2 = 'forces/E'
  []
  [residual]
    type = NEML2StressDivergence
    assembly = 'assembly'
    fem = 'fem'
    executor = 'linear_elasticity'
    stress = 'state/S'
    residual = 'NONTIME'
  []
[]

[Executioner]
  type = Transient

  [TimeIntegrator]
    type = NEML2CentralDifference
    mass_matrix_tag = 'mass'
    use_constant_mass = true
    second_order_vars = 'disp_x disp_y disp_z'
    assembly = 'assembly'
    fem = 'fem'
  []

  start_time = 0.0
  num_steps = 30
  dt = 0.02
[]
