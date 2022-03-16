//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "NodalUserObject.h"

class NodesetModifier : public NodalUserObject
{
public:
  static InputParameters validParams();

  NodesetModifier(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void threadJoin(const UserObject & /*uo*/) override {}
  virtual void finalize() override;

protected:
  // Compute the subdomain ID of the current element
  virtual bool shouldModify() = 0;

private:
  // Nodes to be modified
  std::vector<const Node *> _nodes_to_modify;

  // The original name and ID of the boundary
  // @{
  BoundaryName _from_boundary_name;
  BoundaryID _from_boundary_id;
  // @}

  // The name and ID of the boundary to set to
  // @{
  BoundaryName _to_boundary_name;
  BoundaryID _to_boundary_id;
  // @}
};
