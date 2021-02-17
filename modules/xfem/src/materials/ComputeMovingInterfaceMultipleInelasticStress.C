//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeMovingInterfaceMultipleInelasticStress.h"

registerMooseObject("XFEMApp", ComputeMovingInterfaceMultipleInelasticStress);

InputParameters
ComputeMovingInterfaceMultipleInelasticStress::validParams()
{
  InputParameters params = ComputeMultipleInelasticStress::validParams();

  params.addRequiredParam<UserObjectName>(
      "geometric_cut_userobject",
      "The geometric cut userobject that provides the moving interface definition.");
  params.addRequiredParam<std::vector<GeometricCutSubdomainID>>("physical_cut_subdomains",
                                                                "Physical cut subdomains");
  return params;
}

ComputeMovingInterfaceMultipleInelasticStress::ComputeMovingInterfaceMultipleInelasticStress(
    const InputParameters & parameters)
  : ComputeMultipleInelasticStress(parameters),
    _cut(&getUserObject<GeometricCutUserObject>("geometric_cut_userobject")),
    _physical_cut_subdomains(
        getParam<std::vector<GeometricCutSubdomainID>>("physical_cut_subdomains"))
{
  _xfem = MooseSharedNamespace::dynamic_pointer_cast<XFEM>(_fe_problem.getXFEM());
  if (_xfem == nullptr)
    mooseError("Problem casting to XFEM in GeometricCutUserObject");
}

bool
ComputeMovingInterfaceMultipleInelasticStress::isElemPhysical(const Elem * e) const
{
  if (std::find(_physical_cut_subdomains.begin(),
                _physical_cut_subdomains.end(),
                _xfem->getGeometricCutSubdomainID(e, _cut)) != _physical_cut_subdomains.end())
    return true;
  return false;
}

void
ComputeMovingInterfaceMultipleInelasticStress::computeQpStress()
{
  if (isElemPhysical(_current_elem))
    ComputeMultipleInelasticStress::computeQpStress();
  else
  {
    ComputeMultipleInelasticStress::initQpStatefulProperties();
    // for (unsigned i_rmm = 0; i_rmm < _num_models; ++i_rmm)
    //   _models[i_rmm]->propagateQpStatefulProperties();
  }
}
