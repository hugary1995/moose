//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADComputeDegradedThermalExpansionEigenstrain.h"

registerMooseObject("TensorMechanicsApp", ADComputeDegradedThermalExpansionEigenstrain);

InputParameters
ADComputeDegradedThermalExpansionEigenstrain::validParams()
{
  InputParameters params = ADComputeThermalExpansionEigenstrainBase::validParams();
  params.addClassDescription("Computes eigenstrain due to thermal expansion "
                             "with a constant coefficient");
  params.addRequiredParam<Real>("thermal_expansion_coeff", "Thermal expansion coefficient");
  params.addParam<MaterialPropertyName>(
      "D_name", "degradation", "Name of material property for energetic degradation function.");
  return params;
}

ADComputeDegradedThermalExpansionEigenstrain::ADComputeDegradedThermalExpansionEigenstrain(
    const InputParameters & parameters)
  : ADComputeThermalExpansionEigenstrainBase(parameters),
    _thermal_expansion_coeff(getParam<Real>("thermal_expansion_coeff")),
    _D(getADMaterialProperty<Real>(_base_name + getParam<MaterialPropertyName>("D_name")))
{
}

void
ADComputeDegradedThermalExpansionEigenstrain::computeThermalStrain(ADReal & thermal_strain)
{
  thermal_strain =
      _D[_qp] * _thermal_expansion_coeff * (_temperature[_qp] - _stress_free_temperature[_qp]);
}
