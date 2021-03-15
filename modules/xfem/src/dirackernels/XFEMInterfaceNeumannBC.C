//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "XFEMInterfaceNeumannBC.h"

registerMooseObject("XFEMApp", XFEMInterfaceNeumannBC);

InputParameters
XFEMInterfaceNeumannBC::validParams()
{
  InputParameters params = XFEMInterfaceKernel::validParams();
  params.addClassDescription("A constant Neumann BC on an XFEM interface.");
  params.addRequiredParam<Real>("value", "Value of the BC");
  return params;
}

XFEMInterfaceNeumannBC::XFEMInterfaceNeumannBC(const InputParameters & parameters)
  : XFEMInterfaceKernel(parameters), _v(getParam<Real>("value"))
{
}

Real
XFEMInterfaceNeumannBC::computeQpResidual()
{
  return -_test[_i][_qp] * _v;
}

Real
XFEMInterfaceNeumannBC::computeQpJacobian()
{
  return 0;
}
