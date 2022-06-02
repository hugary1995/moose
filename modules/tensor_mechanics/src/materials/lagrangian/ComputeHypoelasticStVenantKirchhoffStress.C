//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeHypoelasticStVenantKirchhoffStress.h"

registerMooseObject("TensorMechanicsApp", ComputeHypoelasticStVenantKirchhoffStress);

InputParameters
ComputeHypoelasticStVenantKirchhoffStress::validParams()
{
  InputParameters params = ComputeLagrangianObjectiveStress::validParams();
  params.addClassDescription(
      "Calculate a small strain elastic stress that is equivalent to the hyperelastic St. "
      "Venant-Kirchhoff model if integrated using the Truesdell rate.");

  params.addParam<MaterialPropertyName>(
      "elasticity_tensor", "elasticity_tensor", "The name of the elasticity tensor.");

  return params;
}

ComputeHypoelasticStVenantKirchhoffStress::ComputeHypoelasticStVenantKirchhoffStress(
    const InputParameters & parameters)
  : ComputeLagrangianObjectiveStress(parameters),
    _elasticity_tensor(getMaterialProperty<RankFourTensor>(
        getParam<MaterialPropertyName>(_base_name + "elasticity_tensor"))),
    _def_grad(getMaterialProperty<RankTwoTensor>(_base_name + "deformation_gradient"))
{
}

void
ComputeHypoelasticStVenantKirchhoffStress::computeQpSmallStress()
{
  // If small kinematics, it falls back to the grade-zero hypoelasticity
  if (!_large_kinematics)
  {
    _small_stress[_qp] = _elasticity_tensor[_qp] * _mechanical_strain[_qp];
    _small_jacobian[_qp] = _elasticity_tensor[_qp];
    return;
  }

  // Large kinematics:
  // First push forward the elasticity tensor
  const RankTwoTensor F = _def_grad[_qp];
  const RankTwoTensor Ft = F.transpose();
  const Real J = F.det();
  const RankFourTensor C0 = _elasticity_tensor[_qp];
  const RankFourTensor C = F.mixedProductIkJl(F) * C0 * Ft.mixedProductIkJl(Ft) / J;

  // Update the small stress
  const RankTwoTensor dD = _strain_increment[_qp];
  const RankTwoTensor dS = C * dD;
  _small_stress[_qp] = _small_stress_old[_qp] + dS;

  if (_fe_problem.currentlyComputingJacobian())
  {
    // Derivative of the deformation gradient w.r.t. the strain increment depends on kinmatics:
    const RankFourTensor Isym(RankFourTensor::initIdentitySymmetricFour);
    const RankTwoTensor dL = RankTwoTensor::Identity() - _inv_df[_qp];
    const RankFourTensor dFddL = _inv_df[_qp].inverse().mixedProductIkJl(Ft);

    // Compute the small Jacobian
    _small_jacobian[_qp] =
        C - dS.outerProduct(_inv_def_grad[_qp].transpose().initialContraction(dFddL));
    for (auto i : make_range(3))
      for (auto j : make_range(3))
        for (auto r : make_range(3))
          for (auto s : make_range(3))
            for (auto m : make_range(3))
              for (auto n : make_range(3))
                for (auto p : make_range(3))
                  for (auto q : make_range(3))
                    for (auto k : make_range(3))
                      for (auto l : make_range(3))
                        _small_jacobian[_qp](i, j, r, s) +=
                            (dFddL(i, m, r, s) * F(j, n) * F(k, p) * F(l, q) +
                             F(i, m) * dFddL(j, n, r, s) * F(k, p) * F(l, q) +
                             F(i, m) * F(j, n) * dFddL(k, p, r, s) * F(l, q) +
                             F(i, m) * F(j, n) * F(k, p) * dFddL(l, q, r, s)) *
                            C0(m, n, p, q) * dL(k, l) / J;
  }
}
