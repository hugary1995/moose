//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADRadialReturnCreepStressIsotropicUpdateBase.h"
#include "RankTwoTensor.h"

InputParameters
ADRadialReturnCreepStressIsotropicUpdateBase::validParams()
{
  InputParameters params = ADRadialReturnStressIsotropicUpdate::validParams();
  params.set<std::string>("effective_inelastic_strain_name") = "effective_creep_strain";
  return params;
}

ADRadialReturnCreepStressIsotropicUpdateBase::ADRadialReturnCreepStressIsotropicUpdateBase(
    const InputParameters & parameters)
  : ADRadialReturnStressIsotropicUpdate(parameters),
    _creep_strain(declareADProperty<RankTwoTensor>(_base_name + "creep_strain")),
    _creep_strain_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "creep_strain"))
{
}

void
ADRadialReturnCreepStressIsotropicUpdateBase::initQpStatefulProperties()
{
  _creep_strain[_qp].zero();

  ADRadialReturnStressIsotropicUpdate::initQpStatefulProperties();
}

void
ADRadialReturnCreepStressIsotropicUpdateBase::propagateQpStatefulProperties()
{
  _creep_strain[_qp] = _creep_strain_old[_qp];

  propagateQpStatefulPropertiesRadialReturn();
}

void
ADRadialReturnCreepStressIsotropicUpdateBase::computeStressFinalize(
    const ADRankTwoTensor & plastic_strain_increment)
{
  _creep_strain[_qp] = _creep_strain_old[_qp] + plastic_strain_increment;
}
