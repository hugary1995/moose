//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Diffusion.h"

class HomogenizedDiffusion : public Diffusion
{
public:
  static InputParameters validParams();

  HomogenizedDiffusion(const InputParameters & parameters);

protected:
  virtual void computeOffDiagJacobianScalar(unsigned int jvar) override;

  const unsigned int _macro_gradient_num;
};
