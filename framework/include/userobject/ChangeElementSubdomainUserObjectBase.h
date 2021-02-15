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

class ChangeElementSubdomainUserObjectBase : public ElementUserObject
{
public:
  static InputParameters validParams();

  ChangeElementSubdomainUserObjectBase(const InputParameters & parameters);

  void initialize() override{};
  void execute() override;
  void threadJoin(const UserObject & /*uo*/) override{};
  void finalize() override;

protected:
  virtual bool shouldChangeSubdomainFromOriginToTarget() = 0;
  virtual bool shouldChangeSubdomainFromTargetToOrigin() = 0;

  subdomain_id_type targetSubdomainID() const { return _target_subdomain_id; }

  subdomain_id_type originSubdomainID() const { return _origin_subdomain_id; }

  BoundaryID movingBoundary() const { return _moving_boundary_id; }

  /**
   * get ranges for use with threading.
   */
  ConstElemRange & movedElemsRange() const { return *_moved_elems_range.get(); }
  ConstBndNodeRange & movedBndNodesRange() const { return *_moved_bnd_nodes_range.get(); }

  std::shared_ptr<DisplacedProblem> _displaced_problem;

private:
  void setMovingBoundaryName(MooseMesh & mesh);

  void updateBoundaryInfo(MooseMesh & mesh, const std::vector<Elem *> & moved_elems);

  void pushBoundarySideInfo(
      MooseMesh & mesh,
      std::unordered_map<processor_id_type, std::vector<std::pair<dof_id_type, unsigned int>>> &
          elems_to_push);

  void pushBoundaryNodeInfo(
      MooseMesh & mesh,
      std::unordered_map<processor_id_type, std::vector<dof_id_type>> & nodes_to_push);

  /**
   * build ranges for use with threading.
   */
  void buildMovedElemsRange();
  void buildMovedBndNodesRange();

  void recordNodeIdsOnElemSide(const Elem * elem,
                               const unsigned short int side,
                               std::set<dof_id_type> & node_ids);

  /**
   * Initialize solutions for the nodes
   */
  void setSolutionForMovedNodes(SystemBase & sys);

  /// activate subdomain ID
  const subdomain_id_type _target_subdomain_id;

  /// inactivate subdomain ID (the subdomain that you want to keep the same)
  const subdomain_id_type _origin_subdomain_id;

  std::vector<Elem *> _moved_elems;
  std::vector<Elem *> _moved_displaced_elems;
  std::set<dof_id_type> _moved_nodes;

  enum Direction
  {
    forward,
    backward
  };
  std::unordered_map<Elem *, Direction> _elem_to_direction_map;

  /**
   * Ranges for use with threading.
   */
  std::unique_ptr<ConstElemRange> _moved_elems_range;
  std::unique_ptr<ConstBndNodeRange> _moved_bnd_nodes_range;

  /// The name of the moving boundary
  const BoundaryName _moving_boundary_name;

  /// The Id of the moving boundary
  BoundaryID _moving_boundary_id;

  const bool _include_domain_boundary;

  const bool _reversible;

  const bool _init_solution_on_move;
};
