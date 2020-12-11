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
    _K(getADMaterialPropertyByName<Real>(_base_name + "bulk_modulus")),
    _G(getADMaterialPropertyByName<Real>(_base_name + "shear_modulus"))
{
}

void
ADComputeIsotropicLinearElasticStress::computeQpStress()
{
  _stress[_qp].fillFromInputVector({_K[_qp] * _mechanical_strain[_qp].trace()});
  _stress[_qp] += 2.0 * _G[_qp] * _mechanical_strain[_qp].deviatoric();

  // Assign value for elastic strain, which is equal to the mechanical strain
  _elastic_strain[_qp] = _mechanical_strain[_qp];
}
