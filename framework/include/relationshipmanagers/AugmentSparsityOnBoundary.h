//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// App includes
#include "AutomaticMortarGeneration.h"
#include "RelationshipManager.h"

// libMesh includes
#include "libmesh/mesh_base.h"

using libMesh::boundary_id_type;
using libMesh::CouplingMatrix;
using libMesh::Elem;
using libMesh::GhostingFunctor;
using libMesh::MeshBase;
using libMesh::processor_id_type;

class AugmentSparsityOnBoundary : public RelationshipManager
{
public:
  AugmentSparsityOnBoundary(const InputParameters &);

  AugmentSparsityOnBoundary(const AugmentSparsityOnBoundary & others);

  static InputParameters validParams();

  virtual void operator()(const MeshBase::const_element_iterator & range_begin,
                          const MeshBase::const_element_iterator & range_end,
                          processor_id_type p,
                          map_type & coupled_elements) override;

  virtual std::unique_ptr<GhostingFunctor> clone() const override
  {
    return std::make_unique<AugmentSparsityOnBoundary>(*this);
  }

  std::string getInfo() const override;

  virtual bool operator>=(const RelationshipManager & other) const override;

protected:
  virtual void internalInit() override;

  bool containsNode(const MeshBase::const_element_iterator & range_begin,
                    const MeshBase::const_element_iterator & range_end,
                    const dof_id_type id) const;

  const dof_id_type _primary_node_id_in;

  const BoundaryName _secondary_node_set_in;

  dof_id_type _primary_node_id;

  /// null matrix for generating full variable coupling
  const CouplingMatrix * const _null_mat = nullptr;
};
