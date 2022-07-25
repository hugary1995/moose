//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ScalarKernel.h"

#include "IntegralConstraint.h"

class IntegralConstraintScalarKernel : public ScalarKernel
{
public:
  static InputParameters validParams();

  IntegralConstraintScalarKernel(const InputParameters & parameters);

  virtual void reinit() override {}

  virtual void computeResidual();

  virtual void computeJacobian();

protected:
  const IntegralConstraint & _integrator;

  const Real & _residual;

  const Real & _jacobian;
};
