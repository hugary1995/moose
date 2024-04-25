//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Action.h"

/**
 * The action associated with the common area under
 * [Physics/SolidMechanics/QuasiStatic] and
 * [Physics/SolidMechanics/Dynamic]
 *
 * This action serves to
 *   1. Check the existence of sub-blocks under QuasiStatic and Dynamic
 *   2. Allow sub-block actions to pull in common parameters through the action warehouse
 */
class SolidMechanicsPhysicsCommon : public Action
{
public:
  static InputParameters validParams();

  SolidMechanicsPhysicsCommon(const InputParameters & parameters);

  virtual void act() override;
};
