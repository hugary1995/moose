//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADGhostPenalty.h"

// MOOSE includes
#include "MooseVariableFE.h"

#include "libmesh/utility.h"

registerMooseObject("XFEMApp", ADGhostPenalty);

InputParameters
ADGhostPenalty::validParams()
{
  InputParameters params = ADDGKernel::validParams();
  params.addClassDescription("DG kernel for ghost penalty");
  params.addParam<MaterialPropertyName>(
      "eta", 1., "The stiffness (or thermal conductivity or viscosity) coefficient.");
  params.addParam<Real>("alpha", 1., "The penalty parameter");
  return params;
}

ADGhostPenalty::ADGhostPenalty(const InputParameters & parameters)
  : ADDGKernel(parameters),
    _eta(getADMaterialProperty<Real>("eta")),
    _alpha(getParam<Real>("alpha"))
{
}

void
ADGhostPenalty::computeResidual()
{
  std::cout << "current element: " << _current_elem->id() << ", side: " << _current_side
            << std::endl;
  ADDGKernel::computeResidual();
}

ADReal
ADGhostPenalty::computeQpResidual(Moose::DGResidualType type)
{
  ADReal r = 0;

  const unsigned int elem_b_order = _var.order();
  const Real h_elem =
      this->_current_elem_volume / this->_current_side_volume * 1. / Utility::pow<2>(elem_b_order);

  ADReal u_flux_jump = _grad_u[_qp] * _normals[_qp] - _grad_u_neighbor[_qp] * _normals[_qp];
  ADReal test_flux_jump =
      _grad_test[_i][_qp] * _normals[_qp] - _grad_test_neighbor[_i][_qp] * _normals[_qp];

  switch (type)
  {
    case Moose::Element:
      r += _eta[_qp] * _alpha * h_elem * u_flux_jump * test_flux_jump;
      break;

    case Moose::Neighbor:
      r += _eta[_qp] * _alpha * h_elem * u_flux_jump * test_flux_jump;
      break;
  }

  return r;
}
