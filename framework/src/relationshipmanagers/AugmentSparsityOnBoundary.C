//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// App includes
#include "AugmentSparsityOnBoundary.h"

// libMesh includes
#include "libmesh/elem.h"
#include "libmesh/mesh_base.h"
#include "libmesh/boundary_info.h"

registerMooseObject("MooseApp", AugmentSparsityOnBoundary);

using namespace libMesh;

InputParameters
AugmentSparsityOnBoundary::validParams()
{
  InputParameters params = RelationshipManager::validParams();
  params.addParam<dof_id_type>(
      "primary_node",
      std::numeric_limits<dof_id_type>::max(),
      "The ID of the primary node. If no ID is provided, first node of secondary set is chosen.");
  params.addRequiredParam<BoundaryName>("secondary_boundary",
                                        "The boundary ID associated with the secondary side");
  params.set<bool>("attach_geometric_early") = false;
  return params;
}

AugmentSparsityOnBoundary::AugmentSparsityOnBoundary(const InputParameters & params)
  : RelationshipManager(params),
    _primary_node_id_in(getParam<dof_id_type>("primary_node")),
    _secondary_node_set_in(getParam<BoundaryName>("secondary_boundary")),
    _primary_node_id(std::numeric_limits<dof_id_type>::max())
{
}

AugmentSparsityOnBoundary::AugmentSparsityOnBoundary(const AugmentSparsityOnBoundary & other)
  : RelationshipManager(other),
    _primary_node_id_in(other._primary_node_id_in),
    _secondary_node_set_in(other._secondary_node_set_in),
    _primary_node_id(other._primary_node_id)
{
}

void
AugmentSparsityOnBoundary::internalInit()
{
  if (_primary_node_id_in != std::numeric_limits<dof_id_type>::max())
  {
    _primary_node_id = _primary_node_id_in;
    return;
  }

  const auto secondary_node_ids =
      _moose_mesh->getNodeList(_moose_mesh->getBoundaryID(_secondary_node_set_in));
  const auto in = std::min_element(secondary_node_ids.begin(), secondary_node_ids.end());
  _primary_node_id = (in == secondary_node_ids.end()) ? DofObject::invalid_id : *in;
  _communicator.min(_primary_node_id);
}

bool
AugmentSparsityOnBoundary::containsNode(const MeshBase::const_element_iterator & range_begin,
                                        const MeshBase::const_element_iterator & range_end,
                                        const dof_id_type id) const
{
  for (const auto & elem : as_range(range_begin, range_end))
    if (elem->on_boundary())
      for (const auto & node : elem->node_ref_range())
        if (node.id() == id)
          return true;
  return false;
}

std::string
AugmentSparsityOnBoundary::getInfo() const
{
  std::ostringstream oss;
  oss << "AugmentSparsityOnBoundary";
  return oss.str();
}

void
AugmentSparsityOnBoundary::operator()(const MeshBase::const_element_iterator & range_begin,
                                      const MeshBase::const_element_iterator & range_end,
                                      const processor_id_type p,
                                      map_type & coupled_elements)
{
  // Don't do anything if the mesh hasn't been set up
  if (!_moose_mesh->getMeshPtr())
    return;

  // If the range doesn't contain the primary node, there is no coupling
  if (!containsNode(range_begin, range_end, _primary_node_id))
    return;

  // Couple the elements that contain node(s) on the secondary node set
  const BoundaryInfo & binfo = _mesh->get_boundary_info();
  const auto secondary = _moose_mesh->getBoundaryID(_secondary_node_set_in);
  for (const auto elem : _mesh->active_element_ptr_range())
  {
    // Don't bother adding elements already on processor p
    if (elem->processor_id() == p)
      continue;

    for (auto side : elem->side_index_range())
      if (binfo.has_boundary_id(elem, side, secondary))
        coupled_elements.insert(std::make_pair(elem, _null_mat));
  }
}

bool
AugmentSparsityOnBoundary::operator>=(const RelationshipManager & other) const
{
  if (auto asob = dynamic_cast<const AugmentSparsityOnBoundary *>(&other))
    if (_secondary_node_set_in >= asob->_secondary_node_set_in && baseGreaterEqual(*asob))
      return true;

  return false;
}
