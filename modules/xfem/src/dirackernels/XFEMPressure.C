//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "XFEMPressure.h"

registerMooseObject("XFEMApp", XFEMPressure);

InputParameters
XFEMPressure::validParams()
{
  InputParameters params = XFEMInterfaceKernel::validParams();
  params.addClassDescription("Applies a pressure on an interface cut by XFEM.");
  params.addRequiredParam<unsigned int>("component", "The component for the pressure");
  params.addParam<Real>("factor", 1.0, "The magnitude to use in computing the pressure");
  params.addParam<FunctionName>("function", "The function that describes the pressure");
  return params;
}

XFEMPressure::XFEMPressure(const InputParameters & parameters)
  : XFEMInterfaceKernel(parameters),
    _component(getParam<unsigned int>("component")),
    _factor(getParam<Real>("factor")),
    _function(isParamValid("function") ? &getFunction("function") : NULL)
{
}

Real
XFEMPressure::computeQpResidual()
{
  Real factor = _factor;

  if (_function)
    factor *= _function->value(_t, _current_point);

  return -_test[_i][_qp] * _interface_normal(_component) * factor;
}
