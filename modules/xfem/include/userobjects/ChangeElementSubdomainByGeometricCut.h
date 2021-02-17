//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ChangeElementSubdomainUserObjectBase.h"
#include "XFEM.h"

class ChangeElementSubdomainByGeometricCut : public ChangeElementSubdomainUserObjectBase
{
public:
  static InputParameters validParams();

  ChangeElementSubdomainByGeometricCut(const InputParameters & parameters);

  virtual bool shouldChangeSubdomainFromOriginToTarget() override;

  virtual bool shouldChangeSubdomainFromTargetToOrigin() override;

protected:
  const GeometricCutUserObject * _cut;

  /// Pointer to the XFEM controller object
  std::shared_ptr<XFEM> _xfem;
};
