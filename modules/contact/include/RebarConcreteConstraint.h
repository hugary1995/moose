//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef RebarConcreteConstraint_H
#define RebarConcreteConstraint_H

// MOOSE includes
#include "NodeElemConstraint.h"
#include "ContactMaster.h"

// Forward Declarations
class RebarConcreteConstraint;
class ContactLineSearchBase;

template <>
InputParameters validParams<RebarConcreteConstraint>();

/**
 * A RebarConcreteConstraint forces the value of a variable to be the same on both sides of an
 * interface.
 */
class RebarConcreteConstraint : public NodeElemConstraint
{
public:
  RebarConcreteConstraint(const InputParameters & parameters);

  virtual void timestepSetup() override;
  virtual void jacobianSetup() override;
  virtual void residualEnd() override;

  virtual Real computeQpSlaveValue() override;

  virtual Real computeQpResidual(Moose::ConstraintType type) override;

  /**
   * Computes the jacobian for the current element.
   */
  virtual void computeJacobian() override;

  /**
   * Compute off-diagonal Jacobian entries
   * @param jvar The index of the coupled variable
   */
  virtual void computeOffDiagJacobian(unsigned int jvar) override;

  virtual Real computeQpJacobian(Moose::ConstraintJacobianType type) override;

  /**
   * Compute off-diagonal Jacobian entries
   * @param type The type of coupling
   * @param jvar The index of the coupled variable
   */
  virtual Real computeQpOffDiagJacobian(Moose::ConstraintJacobianType type,
                                        unsigned int jvar) override;

  /**
   * Get the dof indices of the nodes connected to the slave node for a specific variable
   * @param var_num The number of the variable for which dof indices are gathered
   * @return bool indicating whether the coupled variable is one of the displacement variables
   */
  virtual void getConnectedDofIndices(unsigned int var_num) override;

  /**
   * Determine whether the coupled variable is one of the displacement variables,
   * and find its component
   * @param var_num The number of the variable to be checked
   * @param component The component index computed in this routine
   * @return bool indicating whether the coupled variable is one of the displacement variables
   */
  bool getCoupledVarComponent(unsigned int var_num, unsigned int & component);

  virtual bool addCouplingEntriesToJacobian() override { return _master_slave_jacobian; }

  bool shouldApply(Node * node);
  void computeContactForce(Node * node, bool update_contact_set);

protected:
  MooseSharedPointer<DisplacedProblem> _displaced_problem;
  FEProblem & _fe_problem;

  const unsigned int _component;
  ContactModel _model;
  const ContactFormulation _formulation;

  const Real _penalty;

  NumericVector<Number> & _residual_copy;
  //  std::map<Point, PenetrationInfo *> _point_to_info;

  const unsigned int _mesh_dimension;

  std::vector<unsigned int> _vars;


  /// Whether to include coupling between the master and slave nodes in the Jacobian
  const bool _master_slave_jacobian;
  /// Whether to include coupling terms with the nodes connected to the slave nodes in the Jacobian
  const bool _connected_slave_nodes_jacobian;
  /// Whether to include coupling terms with non-displacement variables in the Jacobian
  const bool _non_displacement_vars_jacobian;

  static Threads::spin_mutex _contact_set_mutex;
  RealVectorValue _contact_force;
};

#endif
