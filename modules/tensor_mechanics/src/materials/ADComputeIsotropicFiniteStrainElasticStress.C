//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADComputeIsotropicFiniteStrainElasticStress.h"

registerMooseObject("TensorMechanicsApp", ADComputeIsotropicFiniteStrainElasticStress);

InputParameters
ADComputeIsotropicFiniteStrainElasticStress::validParams()
{
  InputParameters params = ADComputeStressBase::validParams();
  params.addClassDescription("Compute stress using elasticity for finite strains");
  return params;
}

ADComputeIsotropicFiniteStrainElasticStress::ADComputeIsotropicFiniteStrainElasticStress(
    const InputParameters & parameters)
  : ADComputeStressBase(parameters),
    GuaranteeConsumer(this),
    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(
        getADMaterialProperty<SymmetricIsotropicRankFourTensor>(_elasticity_tensor_name)),
    _strain_increment(getADMaterialPropertyByName<RankTwoTensor>(_base_name + "strain_increment")),
    _rotation_increment(
        getADMaterialPropertyByName<RankTwoTensor>(_base_name + "rotation_increment")),
    _stress_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "stress")),
    _elastic_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "elastic_strain"))
{
}

void
ADComputeIsotropicFiniteStrainElasticStress::initialSetup()
{
  if (!hasGuaranteedMaterialProperty(_elasticity_tensor_name, Guarantee::ISOTROPIC))
    mooseError("ADComputeIsotropicFiniteStrainElasticStress can only be used with elasticity "
               "tensor materials "
               "that guarantee isotropic tensors.");
}

void
ADComputeIsotropicFiniteStrainElasticStress::computeQpStress()
{
  // Calculate the stress in the intermediate configuration
  ADRankTwoTensor intermediate_stress;

  intermediate_stress =
      _elasticity_tensor[_qp] * (_strain_increment[_qp] + _elastic_strain_old[_qp]);

  // Rotate the stress state to the current configuration
  _stress[_qp] =
      _rotation_increment[_qp] * intermediate_stress * _rotation_increment[_qp].transpose();

  // Assign value for elastic strain, which is equal to the mechanical strain
  _elastic_strain[_qp] = _mechanical_strain[_qp];
}
