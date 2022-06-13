//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VectorizedMaterialFake.h"

registerMooseObject("TensorMechanicsApp", VectorizedMaterialFake);

InputParameters
VectorizedMaterialFake::validParams()
{
  InputParameters params = VectorizedMaterialBase::validParams();
  return params;
}

VectorizedMaterialFake::VectorizedMaterialFake(const InputParameters & parameters)
  : VectorizedMaterialBase(parameters)
{
}

void
VectorizedMaterialFake::GPUCalls()
{
  // This is merely for testing purposes
  for (auto elem : make_range(_mesh.maxElemId()))
    for (auto qp : make_range(_material_data->nQPoints()))
    {
      const RankFourTensor C = get<RankFourTensor>(_input, "elasticity_tensor", elem, qp);
      const RankTwoTensor E = get<RankTwoTensor>(_input, "mechanical_strain", elem, qp);
      fill<RankTwoTensor>(_output, "small_stress", elem, qp, C * E);
      fill<RankFourTensor>(_output, "small_jacobian", elem, qp, C);
    }
}
