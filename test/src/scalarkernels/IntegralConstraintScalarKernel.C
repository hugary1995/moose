//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "IntegralConstraintScalarKernel.h"

#include "Assembly.h"

registerMooseObject("MooseTestApp", IntegralConstraintScalarKernel);

InputParameters
IntegralConstraintScalarKernel::validParams()
{
  InputParameters params = ScalarKernel::validParams();
  params.addRequiredParam<UserObjectName>("integrator",
                                          "The integrator user object doing the "
                                          "element calculations.");

  return params;
}

IntegralConstraintScalarKernel::IntegralConstraintScalarKernel(const InputParameters & parameters)
  : ScalarKernel(parameters),
    _integrator(getUserObject<IntegralConstraint>("integrator")),
    _residual(_integrator.getResidual()),
    _jacobian(_integrator.getJacobian())
{
}

void
IntegralConstraintScalarKernel::computeResidual()
{
  DenseVector<Number> & re = _assembly.residualBlock(_var.number());
  re(0) += _residual;
}

void
IntegralConstraintScalarKernel::computeJacobian()
{
  DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
  ke(0, 0) += _jacobian;
}
