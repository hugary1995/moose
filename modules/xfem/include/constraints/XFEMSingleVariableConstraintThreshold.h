//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "ElemElemConstraint.h"
#include "MooseMesh.h"
#include "SystemBase.h"

// Forward Declarations

class XFEM;
class Function;

class XFEMSingleVariableConstraintThreshold : public ElemElemConstraint
{
public:
  static InputParameters validParams();

  XFEMSingleVariableConstraintThreshold(const InputParameters & parameters);
  virtual ~XFEMSingleVariableConstraintThreshold();

protected:
  virtual void reinitConstraintQuadrature(const ElementPairInfo & element_pair_info) override;

  virtual Real computeQpResidual(Moose::DGResidualType type) override;

  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;

  /// Vector normal to the internal interface
  Point _interface_normal;

  /// Stabilization parameter in Nitsche's formulation and penalty factor in the
  /// Penalty Method
  Real _alpha;

  /// Change in variable value at the interface
  const Function & _jump;

  /// Change in flux of variable value at the interface
  const Function & _jump_flux;

  MooseVariable & _indicator_var;
  const VariableValue & _indicator;
  const VariableValue & _indicator_neighbor;
  Real _threshold;

  /// Pointer to the XFEM controller object
  std::shared_ptr<XFEM> _xfem;
};
