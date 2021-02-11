//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeometricCutUserObject.h"

class ComboCutUserObject : public GeometricCutUserObject
{
public:
  static InputParameters validParams();

  ComboCutUserObject(const InputParameters & parameters);

  virtual bool cutElementByGeometry(const Elem * elem,
                                    std::vector<Xfem::CutEdge> & cut_edges,
                                    std::vector<Xfem::CutNode> & cut_nodes,
                                    Real time) const override;
  virtual bool cutElementByGeometry(const Elem * elem,
                                    std::vector<Xfem::CutFace> & cut_faces,
                                    Real time) const override;

  virtual bool cutFragmentByGeometry(std::vector<std::vector<Point>> & frag_edges,
                                     std::vector<Xfem::CutEdge> & cut_edges,
                                     Real time) const override;
  virtual bool cutFragmentByGeometry(std::vector<std::vector<Point>> & frag_faces,
                                     std::vector<Xfem::CutFace> & cut_faces,
                                     Real time) const override;

  virtual const std::vector<Point>
  getCrackFrontPoints(unsigned int num_crack_front_points) const override;

  virtual GeometricCutSubdomainID getCutSubdomainID(const Node * node) const override;

protected:
private:
  void buildMap();
  const std::vector<UserObjectName> _cut_names;
  unsigned int _num_cuts;
  std::vector<const GeometricCutUserObject *> _cuts;
  std::vector<GeometricCutSubdomainID> _keys;
  std::vector<GeometricCutSubdomainID> _vals;
  std::map<std::vector<GeometricCutSubdomainID>, GeometricCutSubdomainID> _combo_ids;
};
