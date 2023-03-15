//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// MOOSE includes
#include "EqualValueBoundaryConstraint.h"
#include "MooseMesh.h"

registerMooseObject("MooseApp", EqualValueBoundaryConstraint);

InputParameters
EqualValueBoundaryConstraint::validParams()
{
  InputParameters params = NodalConstraint::validParams();
  params.addClassDescription(
      "Constraint for enforcing that variables on each side of a boundary are equivalent.");
  params.addParam<dof_id_type>(
      "primary",
      std::numeric_limits<dof_id_type>::max(),
      "The ID of the primary node. If no ID is provided, first node of secondary set is chosen.");
  params.addParam<std::vector<unsigned int>>("secondary_node_ids", "The IDs of the secondary node");
  params.addParam<BoundaryName>(
      "secondary", "NaN", "The boundary ID associated with the secondary side");
  params.addRequiredParam<Real>("penalty", "The penalty used for the boundary term");

  params.addRelationshipManager("AugmentSparsityOnBoundary",
                                Moose::RelationshipManagerType::COUPLING,
                                [](const InputParameters & obj_params, InputParameters & rm_params)
                                {
                                  rm_params.set<bool>("use_displaced_mesh") =
                                      obj_params.get<bool>("use_displaced_mesh");
                                  rm_params.set<BoundaryName>("secondary_boundary") =
                                      obj_params.get<BoundaryName>("secondary");
                                  rm_params.set<dof_id_type>("primary_node") =
                                      obj_params.get<dof_id_type>("primary");
                                });

  return params;
}

EqualValueBoundaryConstraint::EqualValueBoundaryConstraint(const InputParameters & parameters)
  : NodalConstraint(parameters), _penalty(getParam<Real>("penalty"))
{
}

Real
EqualValueBoundaryConstraint::computeQpResidual(Moose::ConstraintType type)
{
  switch (type)
  {
    case Moose::Secondary:
      return (_u_secondary[_i] - _u_primary[_j]) * _penalty;
    case Moose::Primary:
      return (_u_primary[_j] - _u_secondary[_i]) * _penalty;
  }
  return 0.;
}

Real
EqualValueBoundaryConstraint::computeQpJacobian(Moose::ConstraintJacobianType type)
{
  switch (type)
  {
    case Moose::SecondarySecondary:
      return _penalty;
    case Moose::SecondaryPrimary:
      return -_penalty;
    case Moose::PrimaryPrimary:
      return _penalty;
    case Moose::PrimarySecondary:
      return -_penalty;
    default:
      mooseError("Unsupported type");
      break;
  }
  return 0.;
}
