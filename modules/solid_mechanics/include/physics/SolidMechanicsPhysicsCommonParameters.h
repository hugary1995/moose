//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InputParameters.h"
#include "SolidMechanicsPropertyRegistry.h"

/**
 * Store parameters common to [Physics/SolidMechanics/QuasiStatic/*] and
 * [Physics/SolidMechanics/Dynamic/*].
 *
 * These parameters can either be specified under
 * [Physics/SolidMechanics/QuasiStatic], [Physics/SolidMechanics/Dynamic], or under any of their
 * sub-blocks.
 *
 * IMPORTANT: This object does not inherit from Action -- it merely stores the parameters.
 */
class SolidMechanicsPhysicsCommonParameters
{
public:
  static InputParameters validParams();

  static MultiMooseEnum outputProperties();
  static MultiMooseEnum materialOutputOrders();
  static MultiMooseEnum materialOutputFamilies();

  static bool isPrefix(const std::string & key);

  static const std::map<std::string, SolidMechanics::PropertyRegistryEntry>
      output_property_restriction;
};
