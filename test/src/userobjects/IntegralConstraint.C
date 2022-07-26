//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "IntegralConstraint.h"

registerMooseObject("MooseTestApp", IntegralConstraint);

InputParameters
IntegralConstraint::validParams()
{
  InputParameters params = ElementUserObject::validParams();
  params.addClassDescription("This object computes the residual and Jacobian to constrain a "
                             "variable gradient to a given target value.");
  params.addRequiredCoupledVar("variable", "Variable to constrain");
  params.addRequiredCoupledVar("scalar_variable", "Scalar variable providing the macro gradient");
  params.addRequiredParam<FunctionName>("target", "Function giving the target to hit");
  params.addRequiredParam<MaterialPropertyName>("diffusivity", "Diffusivity");
  params.addRequiredParam<bool>("strain_constraint",
                                "Whether to constrain the strain or the stress");
  return params;
}

IntegralConstraint::IntegralConstraint(const InputParameters & parameters)
  : ElementUserObject(parameters),
    _grad_u(coupledGradient("variable")),
    _macro_gradient(coupledScalarValue("scalar_variable")),
    _grad_phi(_assembly.gradPhi()),
    _target(getFunction("target")),
    _D(getMaterialProperty<Real>("diffusivity")),
    _strain_constraint(getParam<bool>("strain_constraint"))
{
}

void
IntegralConstraint::initialize()
{
  _residual = 0;
  _jacobian = 0;
}

void
IntegralConstraint::execute()
{
  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
  {
    Real strain = _grad_u[_qp](0) + _macro_gradient[0];
    if (_strain_constraint)
    {
      _residual += (strain - _target.value(_t, _q_point[_qp])) * _JxW[_qp] * _coord[_qp];
      _jacobian += 1.0 * _JxW[_qp] * _coord[_qp];
    }
    else
    {
      _residual += (_D[_qp] * strain - _target.value(_t, _q_point[_qp])) * _JxW[_qp] * _coord[_qp];
      _jacobian += _D[_qp] * _JxW[_qp] * _coord[_qp];
    }
  }
}

void
IntegralConstraint::threadJoin(const UserObject & y)
{
  const IntegralConstraint & other = static_cast<const IntegralConstraint &>(y);
  _residual += other._residual;
  _jacobian += other._jacobian;
}

void
IntegralConstraint::finalize()
{
  gatherSum(_residual);
  gatherSum(_jacobian);
}

const Real &
IntegralConstraint::getResidual() const
{
  return _residual;
}

const Real &
IntegralConstraint::getJacobian() const
{
  return _jacobian;
}