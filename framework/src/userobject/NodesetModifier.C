//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NodesetModifier.h"
#include "DisplacedProblem.h"

#include "libmesh/parallel_algebra.h"
#include "libmesh/parallel.h"
#include "libmesh/dof_map.h"
#include "libmesh/remote_elem.h"
#include "libmesh/parallel_ghost_sync.h"

InputParameters
NodesetModifier::validParams()
{
  InputParameters params = NodalUserObject::validParams();
  params.addClassDescription("Modify the nodeset a node belongs to.");
  params.addRequiredParam<BoundaryName>("from_boundary", "The original boundary name");
  params.addRequiredParam<BoundaryName>("to_boundary", "The boundary name to set to");
  return params;
}

NodesetModifier::NodesetModifier(const InputParameters & parameters)
  : NodalUserObject(parameters),
    _from_boundary_name(getParam<BoundaryName>("to_boundary")),
    _from_boundary_id(_mesh.getBoundaryID(_from_boundary_name)),
    _to_boundary_name(getParam<BoundaryName>("to_boundary"))
{
  // Make sure the restricted boundaries contains the from_boundary
  const auto & boundary_ids = boundaryIDs();
  if (std::find(boundary_ids.begin(), boundary_ids.end(), _from_boundary_id) != boundary_ids.end())
    mooseError("boundary should include from_boundary");

  // Get the ID of the to_boundary. Create to_boundary if it doesn't exist yet.
  const std::vector<BoundaryID> to_boundary_ids = _mesh.getBoundaryIDs({{_to_boundary_name}}, true);
  mooseAssert(to_boundary_ids.size() == 1, "Expect exactly one boundary ID.");
  _to_boundary_id = to_boundary_ids[0];
}

void
NodesetModifier::initialize()
{
  _nodes_to_modify.clear();
}

void
NodesetModifier::execute()
{
  // First, get the boundary ids of the current node
  BoundaryInfo & bnd_info = _mesh.getMesh().get_boundary_info();
  std::vector<BoundaryID> current_boundary_ids;
  bnd_info.boundary_ids(_current_node, current_boundary_ids);

  // If we want to modify the boundary id of the current node AND
  // the current node is on the from_boundary AND
  // the current node isn't already on the to_boundary
  if (shouldModify() &&
      std::find(current_boundary_ids.begin(), current_boundary_ids.end(), _from_boundary_id) !=
          current_boundary_ids.end())
    _nodes_to_modify.push_back(_current_node);
}

void
NodesetModifier::finalize()
{
  auto & bnd_info = _mesh.getMesh().get_boundary_info();
  for (const Node * n : _nodes_to_modify)
  {
    bnd_info.remove_node(n, _from_boundary_id);
    bnd_info.add_node(n, _to_boundary_id);
  }
  _mesh.update();

  // Reinit equation systems
  _fe_problem.meshChanged();
}
