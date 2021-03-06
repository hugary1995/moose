//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeThermalExpansionEigenstrainBase.h"
#include "RankTwoTensor.h"

InputParameters
ComputeThermalExpansionEigenstrainBase::validParams()
{
  InputParameters params = ComputeEigenstrainBase::validParams();
  params.addCoupledVar("temperature", "Coupled temperature");
  params.addRequiredCoupledVar("stress_free_temperature",
                               "Reference temperature at which there is no "
                               "thermal expansion for thermal eigenstrain "
                               "calculation");
  return params;
}

ComputeThermalExpansionEigenstrainBase::ComputeThermalExpansionEigenstrainBase(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<ComputeEigenstrainBase>(parameters),
    _temperature(coupledValue("temperature")),
    _deigenstrain_dT(declarePropertyDerivative<RankTwoTensor>(_eigenstrain_name,
                                                              getVar("temperature", 0)->name())),
    _stress_free_temperature(coupledValue("stress_free_temperature"))
{
}

void
ComputeThermalExpansionEigenstrainBase::computeQpEigenstrain()
{
  Real thermal_strain = 0.0;

  // instantaneous_cte is just the derivative of thermal_strain with respect to temperature
  Real instantaneous_cte = 0.0;

  computeThermalStrain(thermal_strain, instantaneous_cte);

  _eigenstrain[_qp].zero();
  _eigenstrain[_qp].addIa(thermal_strain);

  _deigenstrain_dT[_qp].zero();
  _deigenstrain_dT[_qp].addIa(instantaneous_cte);
}
