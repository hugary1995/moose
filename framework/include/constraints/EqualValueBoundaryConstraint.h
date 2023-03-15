//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "NodalConstraint.h"

class EqualValueBoundaryConstraint : public NodalConstraint
{
public:
  static InputParameters validParams();

  EqualValueBoundaryConstraint(const InputParameters & parameters);

protected:
  /**
   * Computes the residual for the current secondary node
   */
  virtual Real computeQpResidual(Moose::ConstraintType type) override;

  /**
   * Computes the jacobian for the constraint
   */
  virtual Real computeQpJacobian(Moose::ConstraintJacobianType type) override;

  // Penalty if constraint is not satisfied
  Real _penalty;
};
