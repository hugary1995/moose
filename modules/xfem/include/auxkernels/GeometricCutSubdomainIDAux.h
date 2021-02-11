//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"
#include "GeometricCutUserObject.h"

class XFEM;

class GeometricCutSubdomainIDAux : public AuxKernel
{
public:
  static InputParameters validParams();

  GeometricCutSubdomainIDAux(const InputParameters & parameters);

protected:
  virtual Real computeValue();

private:
  std::shared_ptr<XFEM> _xfem;
  const GeometricCutUserObject * _cut;
};
