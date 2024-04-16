//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SolidMechanicsPhysicsCommon.h"
#include "SolidMechanicsPhysicsBase.h"
#include "ActionWarehouse.h"

registerMooseAction("SolidMechanicsApp", SolidMechanicsPhysicsCommon, "meta_action");

InputParameters
SolidMechanicsPhysicsCommon::validParams()
{
  InputParameters params = SolidMechanicsPhysicsBase::validParams();
  params.addClassDescription("Store common solid mechanics parameters");
  return params;
}

SolidMechanicsPhysicsCommon::SolidMechanicsPhysicsCommon(const InputParameters & parameters)
  : Action(parameters)
{
}

void
SolidMechanicsPhysicsCommon::act()
{
  // Check if sub-blocks block are found which will use the common parameters.
  auto actions = _awh.getActions<SolidMechanicsPhysicsBase>();
  if (actions.empty())
    mooseWarning("Common parameters for the solid mechanics physics are supplied, but not used in ",
                 parameters().blockLocation());
}
