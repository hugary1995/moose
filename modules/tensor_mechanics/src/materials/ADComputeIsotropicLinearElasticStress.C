//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADComputeIsotropicLinearElasticStress.h"

registerMooseObject("TensorMechanicsApp", ADComputeIsotropicLinearElasticStress);

InputParameters
ADComputeIsotropicLinearElasticStress::validParams()
{
  InputParameters params = ADComputeStressBase::validParams();
  params.addClassDescription("Compute stress for small strains");
  return params;
}

ADComputeIsotropicLinearElasticStress::ADComputeIsotropicLinearElasticStress(
    const InputParameters & parameters)
  : ADComputeStressBase(parameters),
    _C(getADMaterialPropertyByName<SymmetricIsotropicRankFourTensor>(_base_name +
                                                                     "elasticity_tensor"))
{
}

void
ADComputeIsotropicLinearElasticStress::computeQpStress()
{
  _stress[_qp] = _C[_qp] * _mechanical_strain[_qp];

  // Assign value for elastic strain, which is equal to the mechanical strain
  _elastic_strain[_qp] = _mechanical_strain[_qp];
}
