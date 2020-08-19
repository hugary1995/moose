//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "XFEMSingleVariableConstraintThreshold.h"

// MOOSE includes
#include "Assembly.h"
#include "ElementPairInfo.h"
#include "FEProblem.h"
#include "GeometricCutUserObject.h"
#include "XFEM.h"
#include "Function.h"

#include "libmesh/quadrature.h"

registerMooseObject("XFEMApp", XFEMSingleVariableConstraintThreshold);

InputParameters
XFEMSingleVariableConstraintThreshold::validParams()
{
  InputParameters params = ElemElemConstraint::validParams();
  params.addParam<Real>("alpha",
                        100,
                        "Stabilization parameter in Nitsche's formulation and penalty factor "
                        "in the Penalty Method. In Nitsche's formulation this should be as "
                        "small as possible while the method is still stable; while in the "
                        "Penalty Method you want this to be quite large (e.g. 1e6).");
  params.addParam<FunctionName>("jump", 0, "Jump at the interface. Can be a Real or FunctionName.");
  params.addParam<FunctionName>(
      "jump_flux", 0, "Flux jump at the interface. Can be a Real or FunctionName.");
  params.addParam<UserObjectName>(
      "geometric_cut_userobject",
      "Name of GeometricCutUserObject associated with this constraint.");
  params.addRequiredParam<NonlinearVariableName>(
      "indicator_variable", "variable to indicate when to deactivate the constraint");
  params.addRequiredParam<Real>(
      "threshold", "when indicator_variable >= threshold, the constraint will be deactivated");
  params.addClassDescription("Enforce constraints on the value or flux associated with a variable "
                             "at an XFEM interface.");
  return params;
}

XFEMSingleVariableConstraintThreshold::XFEMSingleVariableConstraintThreshold(
    const InputParameters & parameters)
  : ElemElemConstraint(parameters),
    _alpha(getParam<Real>("alpha")),
    _jump(getFunction("jump")),
    _jump_flux(getFunction("jump_flux")),
    _indicator_var(_sys.getFieldVariable<Real>(
        _tid, parameters.get<NonlinearVariableName>("indicator_variable"))),
    _indicator(_indicator_var.sln()),
    _indicator_neighbor(_indicator_var.slnNeighbor()),
    _threshold(getParam<Real>("threshold"))
{
  _xfem = std::dynamic_pointer_cast<XFEM>(_fe_problem.getXFEM());
  if (_xfem == nullptr)
    mooseError("Problem casting to XFEM in XFEMSingleVariableConstraintThreshold");

  const UserObject * uo =
      &(_fe_problem.getUserObjectBase(getParam<UserObjectName>("geometric_cut_userobject")));

  if (dynamic_cast<const GeometricCutUserObject *>(uo) == nullptr)
    mooseError(
        "UserObject casting to GeometricCutUserObject in XFEMSingleVariableConstraintThreshold");

  _interface_id = _xfem->getGeometricCutID(dynamic_cast<const GeometricCutUserObject *>(uo));
}

XFEMSingleVariableConstraintThreshold::~XFEMSingleVariableConstraintThreshold() {}

void
XFEMSingleVariableConstraintThreshold::reinitConstraintQuadrature(
    const ElementPairInfo & element_pair_info)
{
  _interface_normal = element_pair_info._elem1_normal;
  ElemElemConstraint::reinitConstraintQuadrature(element_pair_info);
}

Real
XFEMSingleVariableConstraintThreshold::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;

  if (_indicator[_qp] >= _threshold || _indicator_neighbor[_qp] >= _threshold)
    return r;

  switch (type)
  {
    case Moose::Element:
      r += 0.5 * _test[_i][_qp] * _jump_flux.value(_t, _u[_qp]);
      r += _alpha * (_u[_qp] - _u_neighbor[_qp] - _jump.value(_t, _u[_qp])) * _test[_i][_qp];
      break;

    case Moose::Neighbor:
      r += 0.5 * _test_neighbor[_i][_qp] * _jump_flux.value(_t, _u_neighbor[_qp]);
      r -= _alpha * (_u[_qp] - _u_neighbor[_qp] - _jump.value(_t, _u_neighbor[_qp])) *
           _test_neighbor[_i][_qp];
      break;
  }
  return r;
}

Real
XFEMSingleVariableConstraintThreshold::computeQpJacobian(Moose::DGJacobianType type)
{
  Real r = 0;

  if (_indicator[_qp] >= _threshold || _indicator_neighbor[_qp] >= _threshold)
    return r;

  switch (type)
  {
    case Moose::ElementElement:
      r += _alpha * _phi[_j][_qp] * _test[_i][_qp];
      break;

    case Moose::ElementNeighbor:
      r -= _alpha * _phi_neighbor[_j][_qp] * _test[_i][_qp];
      break;

    case Moose::NeighborElement:
      r -= _alpha * _phi[_j][_qp] * _test_neighbor[_i][_qp];
      break;

    case Moose::NeighborNeighbor:
      r += _alpha * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp];
      break;
  }

  return r;
}
