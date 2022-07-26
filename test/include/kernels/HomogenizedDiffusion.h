//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MatDiffusion.h"

class HomogenizedDiffusion : public MatDiffusion
{
public:
  static InputParameters validParams();

  HomogenizedDiffusion(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual void computeOffDiagJacobianScalar(unsigned int jvar) override;

  const unsigned int _macro_gradient_num;

  const VariableValue & _h;

  const bool _strain_constraint;
};
