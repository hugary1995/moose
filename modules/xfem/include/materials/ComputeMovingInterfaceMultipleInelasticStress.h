//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ComputeMultipleInelasticStress.h"
#include "XFEM.h"

class ComputeMovingInterfaceMultipleInelasticStress : public ComputeMultipleInelasticStress
{
public:
  static InputParameters validParams();

  ComputeMovingInterfaceMultipleInelasticStress(const InputParameters & parameters);

protected:
  virtual void computeQpStress() override;

private:
  bool isElemPhysical(const Elem *) const;

  const GeometricCutUserObject * _cut;
  const std::vector<GeometricCutSubdomainID> _physical_cut_subdomains;
  /// Pointer to the XFEM controller object
  std::shared_ptr<XFEM> _xfem;
};
