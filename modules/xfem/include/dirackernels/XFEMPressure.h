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
#include "XFEMIntegratedBC.h"
#include "Function.h"

class XFEMPressure : public XFEMIntegratedBC
{
public:
  static InputParameters validParams();

  XFEMPressure(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  /// Component of the normal traction
  const int _component;

  /// Factor to be multiplied to the pressure
  const Real _factor;

  /// Function to be multiplied to the pressure
  const Function * const _function;
};
