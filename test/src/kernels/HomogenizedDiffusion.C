//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HomogenizedDiffusion.h"

registerMooseObject("MooseTestApp", HomogenizedDiffusion);

InputParameters
HomogenizedDiffusion::validParams()
{
  InputParameters params = MatDiffusion::validParams();
  params.addRequiredCoupledVar("scalar_variable", "Scalar variable providing the macro gradient");
  return params;
}

HomogenizedDiffusion::HomogenizedDiffusion(const InputParameters & parameters)
  : MatDiffusion(parameters),
    _macro_gradient_num(coupledScalar("scalar_variable")),
    _h(coupledScalarValue("scalar_variable"))
{
}

Real
HomogenizedDiffusion::computeQpResidual()
{
  RealVectorValue macro_gradient(_h[0], 0, 0);
  return _D[_qp] * (_grad_u[_qp] + macro_gradient) * _grad_test[_i][_qp];
}

Real
HomogenizedDiffusion::computeQpJacobian()
{
  return _D[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
}

void
HomogenizedDiffusion::computeOffDiagJacobianScalar(unsigned int jvar)
{
  if (jvar == _macro_gradient_num)
  {
    DenseMatrix<Number> & ken = _assembly.jacobianBlock(_var.number(), jvar);
    DenseMatrix<Number> & kne = _assembly.jacobianBlock(jvar, _var.number());

    for (_i = 0; _i < _test.size(); _i++)
    {
      _j = _i;
      for (_qp = 0; _qp < _qrule->n_points(); _qp++)
      {
        ken(_i, 0) += _grad_test[_i][_qp](0) * _D[_qp] * _JxW[_qp] * _coord[_qp];
        kne(0, _i) += _grad_phi[_j][_qp](0) * _JxW[_qp] * _coord[_qp];
      }
    }
  }
}