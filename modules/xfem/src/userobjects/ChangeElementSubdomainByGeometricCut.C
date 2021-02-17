//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ChangeElementSubdomainByGeometricCut.h"

registerMooseObject("XFEMApp", ChangeElementSubdomainByGeometricCut);

InputParameters
ChangeElementSubdomainByGeometricCut::validParams()
{
  InputParameters params = ChangeElementSubdomainUserObjectBase::validParams();

  params.addRequiredParam<UserObjectName>(
      "cut",
      "The geometric cut userobject that determines the GeometricCutSubdomainID for each element.");

  return params;
}

ChangeElementSubdomainByGeometricCut::ChangeElementSubdomainByGeometricCut(
    const InputParameters & parameters)
  : ChangeElementSubdomainUserObjectBase(parameters),
    _cut(&getUserObject<GeometricCutUserObject>("cut"))
{
  _xfem = MooseSharedNamespace::dynamic_pointer_cast<XFEM>(_fe_problem.getXFEM());
  if (_xfem == nullptr)
    mooseError("Problem casting to XFEM in GeometricCutUserObject");
}

bool
ChangeElementSubdomainByGeometricCut::shouldChangeSubdomainFromOriginToTarget()
{
  return _xfem->getGeometricCutSubdomainID(_current_elem, _cut) == targetSubdomainID();
}

bool
ChangeElementSubdomainByGeometricCut::shouldChangeSubdomainFromTargetToOrigin()
{
  return _xfem->getGeometricCutSubdomainID(_current_elem, _cut) == originSubdomainID();
}
