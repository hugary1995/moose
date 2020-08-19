[Mesh]
  type = GeneratedMesh
  dim = 1
[]

[Variables]
  [c]
  []
[]

[Kernels]
  [diff]
    type = MatDiffusion
    variable = 'c'
    diffusivity = 'dF/dc'
  []
[]

[Materials]
  [dF_dc]
    type = DerivativeParsedMaterial
    f_name = F
    args = 'c'
    function = '1'
    # function = 'c' # this works
  []
[]

[Executioner]
  type = Steady
[]
