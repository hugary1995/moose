//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "VectorizedMaterialBase.h"

// This is a faked vectorized material for testing purposes
class VectorizedMaterialFake : public VectorizedMaterialBase
{
public:
  static InputParameters validParams();

  VectorizedMaterialFake(const InputParameters & parameters);

protected:
  virtual void GPUCalls() override;
};
