//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementUserObject.h"
#include "NonlinearSystemBase.h"
#include "AuxiliarySystem.h"

class ElementSubdomainModifier : public ElementUserObject
{
public:
  static InputParameters validParams();

  ElementSubdomainModifier(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void threadJoin(const UserObject & /*uo*/) override{};
  virtual void finalize() override;

protected:
  virtual SubdomainID computeSubdomainID() = 0;

  BoundaryID movingBoundaryID() const
  {
    if (!_moving_boundary_specified)
      mooseError("no moving boundary specified");
    return _moving_boundary_id;
  }

  const BoundaryName movingBoundaryName() const
  {
    if (!_moving_boundary_specified)
      mooseError("no moving boundary specified");
    return _moving_boundary_name;
  }

  ConstElemRange & movedElemsRange() const { return *_moved_elems_range.get(); }

  ConstBndNodeRange & movedBndNodesRange() const { return *_moved_bnd_nodes_range.get(); }

  DisplacedProblem * _displaced_problem;

private:
  void setMovingBoundaryName(MooseMesh & mesh);

  void updateBoundaryInfo(MooseMesh & mesh, const std::vector<const Elem *> & moved_elems);

  void recordNodeIdsOnElemSide(const Elem * elem,
                               const unsigned short int side,
                               std::set<dof_id_type> & node_ids);

  void pushBoundarySideInfo(
      MooseMesh & mesh,
      std::unordered_map<processor_id_type, std::vector<std::pair<dof_id_type, unsigned int>>> &
          elems_to_push);

  void pushBoundaryNodeInfo(
      MooseMesh & mesh,
      std::unordered_map<processor_id_type, std::vector<dof_id_type>> & nodes_to_push);

  void buildMovedElemsRange();
  void buildMovedBndNodesRange();

  /**
   * Initialize solutions for the nodes
   */
  void setSolutionForMovedNodes(SystemBase & sys);

  std::vector<const Elem *> _moved_elems;
  std::vector<const Elem *> _moved_displaced_elems;
  std::set<dof_id_type> _moved_nodes;

  std::unique_ptr<ConstElemRange> _moved_elems_range;
  std::unique_ptr<ConstBndNodeRange> _moved_bnd_nodes_range;

  const bool _apply_ic;

  /// Whether a moving boundary name is provided
  const bool _moving_boundary_specified;

  /// The name of the moving boundary
  BoundaryName _moving_boundary_name;

  /// The Id of the moving boundary
  BoundaryID _moving_boundary_id;
};
