//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// MOOSE includes
#include "RebarConcreteConstraint.h"
#include "FEProblem.h"
#include "DisplacedProblem.h"
#include "AuxiliarySystem.h"
#include "SystemBase.h"
#include "Assembly.h"
#include "MooseMesh.h"
#include "Executioner.h"
#include "AddVariableAction.h"

#include "libmesh/string_to_enum.h"
#include "libmesh/sparse_matrix.h"

registerMooseObject("ContactApp", RebarConcreteConstraint);

template <>
InputParameters
validParams<RebarConcreteConstraint>()
{
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());

  InputParameters params = validParams<NodeElemConstraint>();
  params.addRequiredParam<unsigned int>("component",
                                        "An integer corresponding to the direction "
                                        "the variable this kernel acts in. (0 for x, "
                                        "1 for y, 2 for z)");

  params.addCoupledVar("disp_x", "The x displacement");
  params.addCoupledVar("disp_y", "The y displacement");
  params.addCoupledVar("disp_z", "The z displacement");

  params.addCoupledVar(
      "displacements",
      "The displacements appropriate for the simulation geometry and coordinate system");

  //only supports glued for now
  params.addParam<std::string>("model", "glued", "The contact model to use");
  //only supports kinematic(default) for now
  params.addParam<std::string>("formulation", "default", "The contact formulation");

  params.addParam<Real>(
      "penalty",
      1e8,
      "The penalty to apply.  This can vary depending on the stiffness of your materials");
  params.addParam<MooseEnum>("order", orders, "The finite element order");

  params.set<bool>("use_displaced_mesh") = false;

  return params;
}

RebarConcreteConstraint::RebarConcreteConstraint(const InputParameters & parameters)
  : NodeElemConstraint(parameters),
    _fe_problem(*parameters.get<FEProblem *>("_fe_problem")),
    _component(getParam<unsigned int>("component")),
    _model(ContactMaster::contactModel(getParam<std::string>("model"))),
    _formulation(ContactMaster::contactFormulation(getParam<std::string>("formulation"))),
    _penalty(getParam<Real>("penalty")),
    _residual_copy(_sys.residualGhosted()),
    _mesh_dimension(_mesh.dimension()),
    _vars(3, libMesh::invalid_uint)
{
  // std::cout << "In RebarConcreteConstraint()" << std::endl;
  _overwrite_slave_residual = false;

  if (isParamValid("displacements"))
  {
    // modern parameter scheme for displacements
    for (unsigned int i = 0; i < coupledComponents("displacements"); ++i)
      _vars[i] = coupled("displacements", i);
  }
  else
  {
    // Legacy parameter scheme for displacements
    if (isParamValid("disp_x"))
      _vars[0] = coupled("disp_x");
    if (isParamValid("disp_y"))
      _vars[1] = coupled("disp_y");
    if (isParamValid("disp_z"))
      _vars[2] = coupled("disp_z");

    mooseDeprecated("use the `displacements` parameter rather than the `disp_*` parameters (those "
                    "will go away with the deprecation of the Solid Mechanics module).");
  }
}

void
RebarConcreteConstraint::timestepSetup()
{

}

void
RebarConcreteConstraint::jacobianSetup()
{

}

bool
RebarConcreteConstraint::shouldApply()
{
  // std::cout << "          In RebarConcreteConstraint::shouldApply" << std::endl;
  //currently always assume in contact
  bool in_contact = true;

  //currently does not support updating contact set
  //bool is_nonlinear = _fe_problem.computingNonlinearResid();

  // This computes the contact force once per constraint, rather than once per quad point
  // and for both master and slave cases.
  if (_component == 0)
    computeContactForce();

  // std::cout << "          Out RebarConcreteConstraint::shouldApply" << std::endl;
  return in_contact;
}

void
RebarConcreteConstraint::computeContactForce()
{
  // std::cout << "          In RebarConcreteConstraint::computeContactForce()" << std::endl;
  //const Node * node = pinfo->_node;
  const Node * node = _current_node;
  // std::cout << "               current node id: " << node->id() << std::endl;

  // Build up residual vector
  RealVectorValue res_vec;
  for (unsigned int i = 0; i < _mesh_dimension; ++i)
  {
    dof_id_type dof_number = node->dof_number(0, _vars[i], 0);
    res_vec(i) = _residual_copy(dof_number);
  }

  switch (_model)
  {
    case CM_GLUED:
      switch (_formulation)
      {
        case CF_KINEMATIC:
          _contact_force =  -res_vec;
          break;

        default:
          mooseError("Invalid contact formulation");
          break;
      }
      break;

    default:
      mooseError("Invalid or unavailable contact model");
      break;
  }
  // std::cout << "          Out RebarConcreteConstraint::computeContactForce()" << std::endl;
}

Real
RebarConcreteConstraint::computeQpSlaveValue()
{
  return _u_slave[_qp];
}

Real
RebarConcreteConstraint::computeQpResidual(Moose::ConstraintType type)
{
  // std::cout << "          In RebarConcreteConstraint::computeQpResidual()" << std::endl;
  Real resid = _contact_force(_component);
  std::cout << "                    component: " << _component << std::endl;
  switch (type)
  {
    case Moose::Slave:
      if (_formulation == CF_KINEMATIC)
      {
        // std::cout << "               _u_slave[_qp] = " << _u_slave[_qp] << std::endl;
        // std::cout << "               _u_master[_qp] = " << _u_master[_qp] << std::endl;
        RealVectorValue distance_vec(_u_slave[_qp] - _u_master[_qp]);
        RealVectorValue pen_force(_penalty * distance_vec);
        if (_model == CM_GLUED)
          resid += pen_force(_component);
      }
      // std::cout << "          Out RebarConcreteConstraint::computeQpResidual(Slave)" << std::endl;
      std::cout << "                    resid = " << resid << std::endl;
      std::cout << "                    _test_slave" << "[" << _i << "][" << _qp << "] = " << _test_slave[_i][_qp] << std::endl;
      return _test_slave[_i][_qp] * resid;

    case Moose::Master:
      // std::cout << "          Out RebarConcreteConstraint::computeQpResidual(Master)" << std::endl;
      std::cout << "                    resid = " << resid << std::endl;
      std::cout << "                    _test_master" << "[" << _i << "][" << _qp << "] = " << _test_master[_i][_qp] << std::endl;
      return _test_master[_i][_qp] * -resid;
  }

  return 0.0;
}

Real
RebarConcreteConstraint::computeQpJacobian(Moose::ConstraintJacobianType type)
{
  const Real penalty = _penalty;

  switch (type)
  {
    case Moose::SlaveSlave:
      switch (_model)
      {
        case CM_GLUED:
          switch (_formulation)
          {
            case CF_KINEMATIC:
            {
              const Real curr_jac = (*_jacobian)(_current_node->dof_number(0, _vars[_component], 0),
                                                 _connected_dof_indices[_j]);
              return -curr_jac + _phi_slave[_j][_qp] * penalty * _test_slave[_i][_qp];
            }

            default:
              mooseError("Invalid contact formulation");
          }
        default:
          mooseError("Invalid or unavailable contact model");
      }

    case Moose::SlaveMaster:
      switch (_model)
      {
        case CM_GLUED:
          switch (_formulation)
          {
            case CF_KINEMATIC:
            {
              Node * curr_master_node = _current_master->get_node(_j);
              const Real curr_jac =
                  (*_jacobian)(_current_node->dof_number(0, _vars[_component], 0),
                               curr_master_node->dof_number(0, _vars[_component], 0));
              return -curr_jac - _phi_master[_j][_qp] * penalty * _test_slave[_i][_qp];
            }

            default:
              mooseError("Invalid contact formulation");
          }

        default:
          mooseError("Invalid or unavailable contact model");
      }

    case Moose::MasterSlave:
      switch (_model)
      {
        case CM_GLUED:
          switch (_formulation)
          {
            case CF_KINEMATIC:
            {
              const Real slave_jac = (*_jacobian)(
                  _current_node->dof_number(0, _vars[_component], 0), _connected_dof_indices[_j]);
              return slave_jac * _test_master[_i][_qp];
            }

            default:
              mooseError("Invalid contact formulation");
          }

        default:
          mooseError("Invalid or unavailable contact model");
      }

    case Moose::MasterMaster:
      switch (_model)
      {
        case CM_GLUED:
          switch (_formulation)
          {
            case CF_KINEMATIC:
              return 0.0;

            default:
              mooseError("Invalid contact formulation");
          }

        default:
          mooseError("Invalid or unavailable contact model");
      }
  }

  return 0.0;
}

Real
RebarConcreteConstraint::computeQpOffDiagJacobian(Moose::ConstraintJacobianType type,
                                                      unsigned int jvar)
{

  switch (type)
  {
    case Moose::SlaveSlave:
      switch (_model)
      {
        case CM_GLUED:
        {
          const Real curr_jac = (*_jacobian)(_current_node->dof_number(0, _vars[_component], 0),
                                             _connected_dof_indices[_j]);
          return -curr_jac;
        }
        default:
          mooseError("Invalid or unavailable contact model");
      }

    case Moose::SlaveMaster:
      switch (_model)
      {
        case CM_GLUED:
          return 0;

        default:
          mooseError("Invalid or unavailable contact model");
      }

    case Moose::MasterSlave:
      switch (_model)
      {
        case CM_GLUED:
          switch (_formulation)
          {
            case CF_KINEMATIC:
            {
              const Real slave_jac = (*_jacobian)(
                  _current_node->dof_number(0, _vars[_component], 0), _connected_dof_indices[_j]);
              return slave_jac * _test_master[_i][_qp];
            }

            default:
              mooseError("Invalid contact formulation");
          }

        default:
          mooseError("Invalid or unavailable contact model");
      }

    case Moose::MasterMaster:
      switch (_model)
      {
        case CM_GLUED:
            return 0.0;

        default:
          mooseError("Invalid or unavailable contact model");
      }
  }

  return 0.0;
}

void
RebarConcreteConstraint::computeJacobian()
{
  // std::cout << "          In RebarConcreteConstraint::computeJacobian()" << std::endl;
  getConnectedDofIndices(_var.number());

  // std::cout << "               _phi_master.size() = " << _phi_master.size() << std::endl;
  // std::cout << "               _phi_slave.size() = " << _phi_slave.size() << std::endl;

  DenseMatrix<Number> & Knn =
      _assembly.jacobianBlockNeighbor(Moose::NeighborNeighbor, _master_var.number(), _var.number());

  _Kee.resize(_test_slave.size(), _connected_dof_indices.size());

  for (_i = 0; _i < _test_slave.size(); _i++)
    // Loop over the connected dof indices so we can get all the jacobian contributions
    for (_j = 0; _j < _connected_dof_indices.size(); _j++)
      _Kee(_i, _j) += computeQpJacobian(Moose::SlaveSlave);

  DenseMatrix<Number> & Ken =
      _assembly.jacobianBlockNeighbor(Moose::ElementNeighbor, _var.number(), _var.number());
  if (Ken.m() && Ken.n())
    for (_i = 0; _i < _test_slave.size(); _i++)
      for (_j = 0; _j < _phi_master.size(); _j++)
        Ken(_i, _j) += computeQpJacobian(Moose::SlaveMaster);

  _Kne.resize(_test_master.size(), _connected_dof_indices.size());
  for (_i = 0; _i < _test_master.size(); _i++)
    // Loop over the connected dof indices so we can get all the jacobian contributions
    for (_j = 0; _j < _connected_dof_indices.size(); _j++)
      _Kne(_i, _j) += computeQpJacobian(Moose::MasterSlave);

  if (Knn.m() && Knn.n())
    for (_i = 0; _i < _test_master.size(); _i++)
      for (_j = 0; _j < _phi_master.size(); _j++)
        Knn(_i, _j) += computeQpJacobian(Moose::MasterMaster);

  // std::cout << "          Out RebarConcreteConstraint::computeJacobian()" << std::endl;
}

void
RebarConcreteConstraint::computeOffDiagJacobian(unsigned int jvar)
{
  getConnectedDofIndices(jvar);

  _Kee.resize(_test_slave.size(), _connected_dof_indices.size());

  DenseMatrix<Number> & Knn =
      _assembly.jacobianBlockNeighbor(Moose::NeighborNeighbor, _master_var.number(), jvar);

  for (_i = 0; _i < _test_slave.size(); _i++)
    // Loop over the connected dof indices so we can get all the jacobian contributions
    for (_j = 0; _j < _connected_dof_indices.size(); _j++)
      _Kee(_i, _j) += computeQpOffDiagJacobian(Moose::SlaveSlave, jvar);

  DenseMatrix<Number> & Ken =
      _assembly.jacobianBlockNeighbor(Moose::ElementNeighbor, _var.number(), jvar);
  for (_i = 0; _i < _test_slave.size(); _i++)
    for (_j = 0; _j < _phi_master.size(); _j++)
      Ken(_i, _j) += computeQpOffDiagJacobian(Moose::SlaveMaster, jvar);

  _Kne.resize(_test_master.size(), _connected_dof_indices.size());
  if (_Kne.m() && _Kne.n())
    for (_i = 0; _i < _test_master.size(); _i++)
      // Loop over the connected dof indices so we can get all the jacobian contributions
      for (_j = 0; _j < _connected_dof_indices.size(); _j++)
        _Kne(_i, _j) += computeQpOffDiagJacobian(Moose::MasterSlave, jvar);

  for (_i = 0; _i < _test_master.size(); _i++)
    for (_j = 0; _j < _phi_master.size(); _j++)
      Knn(_i, _j) += computeQpOffDiagJacobian(Moose::MasterMaster, jvar);
}

void
RebarConcreteConstraint::getConnectedDofIndices(unsigned int var_num)
{
  // std::cout << "          In RebarConcreteConstraint::getConnectedDofIndices()" << std::endl;
  unsigned int component;
  if (getCoupledVarComponent(var_num, component))
    NodeElemConstraint::getConnectedDofIndices(var_num);

  _phi_slave.resize(_connected_dof_indices.size());

  dof_id_type current_node_var_dof_index = _sys.getVariable(0, var_num).nodalDofIndex();

  // Fill up _phi_slave so that it is 1 when j corresponds to the dof associated with this node
  // and 0 for every other dof
  // This corresponds to evaluating all of the connected shape functions at _this_ node
  _qp = 0;
  for (unsigned int j = 0; j < _connected_dof_indices.size(); j++)
  {
    _phi_slave[j].resize(1);

    if (_connected_dof_indices[j] == current_node_var_dof_index)
      _phi_slave[j][_qp] = 1.0;
    else
      _phi_slave[j][_qp] = 0.0;
  }
  // std::cout << "          Out RebarConcreteConstraint::getConnectedDofIndices()" << std::endl;
}

bool
RebarConcreteConstraint::getCoupledVarComponent(unsigned int var_num, unsigned int & component)
{
  component = std::numeric_limits<unsigned int>::max();
  bool coupled_var_is_disp_var = false;
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  {
    if (var_num == _vars[i])
    {
      coupled_var_is_disp_var = true;
      component = i;
      break;
    }
  }
  return coupled_var_is_disp_var;
}

void
RebarConcreteConstraint::residualEnd()
{

}
