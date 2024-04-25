//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SolidMechanicsPhysicsSubBlock.h"

InputParameters
SolidMechanicsPhysicsSubBlock::validParams()
{
  InputParameters params = SolidMechanicsPhysicsCommon::validParams();
  return params;
}

SolidMechanicsPhysicsSubBlock::SolidMechanicsPhysicsSubBlock(const InputParameters & parameters)
  : Action(parameters)
{
  applyCommonParameters();
  modifyParameters();
}

void
SolidMechanicsPhysicsSubBlock::applyCommonParameters()
{
  auto action = _awh.getActions<SolidMechanicsPhysicsCommon>();
  mooseAssert(!action.empty(), "SolidMechanicsPhysicsCommon not found");

  if (action.size() != 1)
    mooseError("Duplicate SolidMechanicsPhysicsCommon actions not allowed. This could happen if "
               "the input file contains both [QuasiStatic] and [Dynamic] sub-blocks under the "
               "[SolidMechanics] physics block.");

  const_cast<InputParameters *>(&parameters())->applyParameters(action[0]->parameters());
}

void
SolidMechanicsPhysicsSubBlock::modifyParameters()
{
  // Set strain_base_name to base_name if nothing is specified by the user
  if (!isParamSetByUser("strain_base_name"))
    params.set<std::string>("strain_base_name") = getParam<std::string>("base_name");

  // Set default incremental based on strain
  if (!isParamSetByUser("incremental"))
  {
    const std::map<Strain, bool> incremental_default = {{"Small", false}, {"Finite", true}};
    params.set<bool>("incremental") = incremental_default.at(getParam<MooseEnum>("strain"));
  }
}
