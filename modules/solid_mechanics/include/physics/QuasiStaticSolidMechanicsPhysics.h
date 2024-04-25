//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "SolidMechanicsPhysicsBase.h"
#include "libmesh/point.h"

/**
 * The action associated with [Physics/SolidMechanics/QuasiStatic/*]
 */
class QuasiStaticSolidMechanicsPhysics : public SolidMechanicsPhysicsBase
{
public:
  static InputParameters validParams();

  QuasiStaticSolidMechanicsPhysicsBase(const InputParameters & params);

  virtual void act();

protected:
  virtual void modifyParameters() override;
  void mergeMultiMooseEnum(const std::string & dest, const std::string & src);

  void actSubdomainChecks();
  void actOutputGeneration();
  void actEigenstrainNames();
  void actOutputMatProp();
  void actGatherActionParameters();
  void actStrain();

  virtual std::string getKernelType();
  virtual InputParameters getKernelParameters(std::string type);

  Moose::CoordinateSystemType _coord_system;

  /// if this vector is not empty the variables, kernels and materials are restricted to these subdomains
  std::vector<SubdomainName> _subdomain_names;

  /// set generated from the passed in vector of subdomain names
  std::set<SubdomainID> _subdomain_ids;

  /// set generated from the combined block restrictions of all SolidMechanics/Master action blocks
  std::set<SubdomainID> _subdomain_id_union;
};
