//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADComputeStressBase.h"

/**
 * ADComputeIsotropicLinearElasticStress computes the stress following linear elasticity theory
 * (small strains)
 */
class ADComputeIsotropicLinearElasticStress : public ADComputeStressBase
{
public:
  static InputParameters validParams();

  ADComputeIsotropicLinearElasticStress(const InputParameters & parameters);

protected:
  virtual void computeQpStress() override;

  const ADMaterialProperty<Real> & _K;
  const ADMaterialProperty<Real> & _G;
};
