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
 * Base class for merging parameters under the common area into the sub-block.
 *
 * This object serves as the intermediate class so that sub-block actions derived from it can
 * directly work with the *union* of the input parameters under the common area and under each of
 * the sub-block.
 */
class SolidMechanicsPhysicsSubBlock : public Action
{
public:
  SolidMechanicsPhysicsSubBlock(const InputParameters & params);

protected:
  /// Pull in parameters under the common area
  virtual void applyCommonParameters() final;

  /// Modify input parameters before they are retrieved
  virtual void modifyParameters();
};
