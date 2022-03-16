//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ThresholdNodesetModifier.h"

InputParameters
ThresholdNodesetModifier::validParams()
{
  InputParameters params = NodesetModifier::validParams();
  params.addRequiredParam<Real>(
      "threshold", "The value above (or below) which to change the boundary ID of a node");
  params.addParam<MooseEnum>("criterion_type",
                             MooseEnum("BELOW EQUAL ABOVE", "ABOVE"),
                             "Criterion to use for the threshold");
  return params;
}

ThresholdNodesetModifier::ThresholdNodesetModifier(const InputParameters & parameters)
  : NodesetModifier(parameters),
    _threshold(getParam<Real>("threshold")),
    _criterion_type(getParam<MooseEnum>("criterion_type").getEnum<CriterionType>())
{
}

bool
ThresholdNodesetModifier::shouldModify()
{
  Real value = computeValue();

  switch (_criterion_type)
  {
    case CriterionType::Equal:
      return MooseUtils::absoluteFuzzyEqual(value - _threshold, 0);

    case CriterionType::Below:
      return value < _threshold;

    case CriterionType::Above:
      return value > _threshold;
  }

  mooseError("Internal error. Unsupported criterion type");
  return false;
}
