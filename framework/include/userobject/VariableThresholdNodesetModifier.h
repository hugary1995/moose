//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ThresholdNodesetModifier.h"
#include "MooseVariableInterface.h"

class VariableThresholdNodesetModifier : public ThresholdNodesetModifier,
                                         public MooseVariableInterface<Real>
{
public:
  static InputParameters validParams();

  VariableThresholdNodesetModifier(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

private:
  /// The coupled variable used in the criterion
  const VariableValue & _u;
};
