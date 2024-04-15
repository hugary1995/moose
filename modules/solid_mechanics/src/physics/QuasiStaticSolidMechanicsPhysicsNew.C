//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "QuasiStaticSolidMechanicsPhysicsNew.h"

registerMooseAction("SolidMechanicsApp", QuasiStaticSolidMechanicsPhysicsNew, "meta_action");
registerMooseAction("SolidMechanicsApp",
                    QuasiStaticSolidMechanicsPhysicsNew,
                    "setup_mesh_complete");
registerMooseAction("SolidMechanicsApp",
                    QuasiStaticSolidMechanicsPhysicsNew,
                    "validate_coordinate_systems");
registerMooseAction("SolidMechanicsApp", QuasiStaticSolidMechanicsPhysicsNew, "add_variable");
registerMooseAction("SolidMechanicsApp", QuasiStaticSolidMechanicsPhysicsNew, "add_aux_variable");
registerMooseAction("SolidMechanicsApp", QuasiStaticSolidMechanicsPhysicsNew, "add_kernel");
registerMooseAction("SolidMechanicsApp", QuasiStaticSolidMechanicsPhysicsNew, "add_aux_kernel");
registerMooseAction("SolidMechanicsApp", QuasiStaticSolidMechanicsPhysicsNew, "add_material");
registerMooseAction("SolidMechanicsApp",
                    QuasiStaticSolidMechanicsPhysicsNew,
                    "add_master_action_material");
registerMooseAction("SolidMechanicsApp", QuasiStaticSolidMechanicsPhysicsNew, "add_scalar_kernel");
registerMooseAction("SolidMechanicsApp", QuasiStaticSolidMechanicsPhysicsNew, "add_user_object");

InputParameters
QuasiStaticSolidMechanicsPhysicsNew::validParams()
{
  InputParameters params = QuasiStaticSolidMechanicsPhysicsBase::validParams();
  params.addClassDescription("Set up stress divergence kernels with coordinate system aware logic");
  return params;
}

QuasiStaticSolidMechanicsPhysicsNew::QuasiStaticSolidMechanicsPhysicsNew(
    const InputParameters & params)
  : QuasiStaticSolidMechanicsPhysicsBase(params)
{
}

void
QuasiStaticSolidMechanicsPhysicsNew::act()
{
  std::cout << "In QuasiStaticSolidMechanicsPhysicsNew::act()" << std::endl;
}
