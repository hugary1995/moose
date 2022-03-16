//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "NodesetModifier.h"
#include "MooseEnum.h"

class ThresholdNodesetModifier : public NodesetModifier
{
public:
  static InputParameters validParams();

  ThresholdNodesetModifier(const InputParameters & parameters);

protected:
  virtual bool shouldModify() override;

  /// Compute the value used in the criterion
  virtual Real computeValue() = 0;

private:
  /// Threshold to modify the element subdomain ID
  const Real _threshold;

  /// Criterion type
  const enum class CriterionType { Below, Equal, Above } _criterion_type;
};
