//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// Moose Includes
#include "XFEMInterfaceKernel.h"

class XFEMInterfaceConvectiveHeatFluxBC : public XFEMInterfaceKernel
{
public:
  static InputParameters validParams();

  XFEMInterfaceConvectiveHeatFluxBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /// Far-field temperature variable
  const MaterialProperty<Real> & _T_infinity;

  /// Convective heat transfer coefficient
  const MaterialProperty<Real> & _htc;

  /// Derivative of convective heat transfer coefficient with respect to temperature
  const MaterialProperty<Real> & _htc_dT;
};
