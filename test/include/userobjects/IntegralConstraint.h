//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementUserObject.h"

#include "Function.h"

class IntegralConstraint : public ElementUserObject
{
public:
  static InputParameters validParams();

  IntegralConstraint(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void threadJoin(const UserObject & y) override;
  virtual void finalize() override;

  virtual const Real & getResidual() const;
  virtual const Real & getJacobian() const;

protected:
  const VariableGradient & _grad_u;

  const VariableValue & _macro_gradient;

  const VariablePhiGradient & _grad_phi;

  const Function & _target;

  unsigned int _qp;

  Real _residual;

  Real _jacobian;
};
