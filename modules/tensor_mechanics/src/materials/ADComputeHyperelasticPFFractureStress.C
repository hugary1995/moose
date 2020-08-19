//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADComputeHyperelasticPFFractureStress.h"

registerADMooseObject("TensorMechanicsApp", ADComputeHyperelasticPFFractureStress);

InputParameters
ADComputeHyperelasticPFFractureStress::validParams()
{
  InputParameters params = Material::validParams();
  params.addParam<std::string>("base_name",
                               "Optional parameter that allows the user to define "
                               "multiple mechanics material systems on the same "
                               "block, i.e. for multiple phases");
  params.suppressParameter<bool>("use_displaced_mesh");
  params.addRequiredCoupledVar("c", "Name of damage variable");
  params.addParam<MaterialPropertyName>(
      "E_name", "elastic_energy", "Name of material property for elastic energy");
  params.addParam<MaterialPropertyName>(
      "D_name", "degradation", "Name of material property for energetic degradation function.");
  params.addParam<MaterialPropertyName>(
      "bulk_modulus_name", "bulk_modulus", "linearized bulk modulus");
  params.addParam<MaterialPropertyName>(
      "shear_modulus_name", "shear_modulus", "linearized shear modulus");
  params.addParam<bool>("history_maximum", false, "Use the history maximum of the elastic energy");
  return params;
}

ADComputeHyperelasticPFFractureStress::ADComputeHyperelasticPFFractureStress(
    const InputParameters & parameters)
  : Material(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _stress(declareADProperty<RankTwoTensor>(_base_name + "stress")),
    _deformation_gradient(
        getADMaterialProperty<RankTwoTensor>(_base_name + "deformation_gradient")),
    _E(declareADProperty<Real>(_base_name + getParam<MaterialPropertyName>("E_name"))),
    _dE_dc(declareADProperty<Real>(derivativePropertyNameFirst(
        _base_name + getParam<MaterialPropertyName>("E_name"), this->getVar("c", 0)->name()))),
    _D(getADMaterialProperty<Real>(_base_name + getParam<MaterialPropertyName>("D_name"))),
    _dD_dc(getADMaterialProperty<Real>(derivativePropertyNameFirst(
        _base_name + getParam<MaterialPropertyName>("D_name"), this->getVar("c", 0)->name()))),
    _K(getADMaterialProperty<Real>(_base_name +
                                   getParam<MaterialPropertyName>("bulk_modulus_name"))),
    _G(getADMaterialProperty<Real>(_base_name +
                                   getParam<MaterialPropertyName>("shear_modulus_name"))),
    _E_active(
        declareADProperty<Real>(_base_name + getParam<MaterialPropertyName>("E_name") + "_active")),
    _E_active_old(getMaterialPropertyOld<Real>(
        _base_name + getParam<MaterialPropertyName>("E_name") + "_active")),
    _hist_max(declareADProperty<Real>(_base_name + getParam<MaterialPropertyName>("E_name") +
                                      "_hist_max")),
    _hist_max_old(getMaterialPropertyOld<Real>(
        _base_name + getParam<MaterialPropertyName>("E_name") + "_hist_max")),
    _hist(getParam<bool>("history_maximum"))
{
}

void
ADComputeHyperelasticPFFractureStress::initQpStatefulProperties()
{
  _E_active[_qp] = 0;
  _hist_max[_qp] = 0;
}

void
ADComputeHyperelasticPFFractureStress::computeQpProperties()
{
  ADRankTwoTensor F = _deformation_gradient[_qp];
  ADRankTwoTensor C = F.transpose() * F;
  ADReal J = F.det();
  ADReal J2 = C.det();

  // dilatational potential
  // U = 0.5 * K (0.5 * (J^2 - 1) - lnJ)
  // Ua = H(J - 1) * U
  // Ui = H(1 - J) * U
  // dU_dJ = 0.5 * K (J - 1/J)
  ADReal U = 0.25 * _K[_qp] * (J2 - std::log(J2) - 1.0);
  ADReal Ua = J > 1.0 ? U : 0.0;
  ADReal Ui = J > 1.0 ? 0.0 : U;
  ADReal dU_dJ = 0.5 * _K[_qp] * (J - 1.0 / J);

  // isochoric potential
  // W = 0.5 * G * (J^(-2/3) tr(C) - 3)
  ADReal W = 0.5 * _G[_qp] * (C.trace() / std::cbrt(J2) - 3.0);

  // elastic potential
  // E = D * (Ua + W) + Ui
  _E_active[_qp] = Ua + W;
  if (_E_active[_qp] > _hist_max_old[_qp])
    _hist_max[_qp] = _E_active[_qp];
  else
    _hist_max[_qp] = _hist_max_old[_qp];
  _E[_qp] = _D[_qp] * (Ua + W) + Ui;
  _dE_dc[_qp] = _dD_dc[_qp] * (_hist ? _hist_max_old[_qp] : _E_active_old[_qp]);

  // PK2 stress
  ADRankTwoTensor I2(ADRankTwoTensor::initIdentity);
  ADRankTwoTensor S_vol = J * dU_dJ * C.inverse();
  ADRankTwoTensor S_dev = _G[_qp] / std::cbrt(J2) * (I2 - C.trace() * C.inverse() / 3.0);
  ADRankTwoTensor S = J > 1.0 ? _D[_qp] * (S_vol + S_dev) : S_vol + _D[_qp] * S_dev;

  // cauchy stress
  _stress[_qp] = F * S * F.transpose() / J;

  // ADReal lambda = _K[_qp] - 2 * _G[_qp] / 3;
  //
  // ADRankTwoTensor F = _deformation_gradient[_qp];
  // ADRankTwoTensor b = F * F.transpose();
  // ADReal J = F.det();
  //
  // // dilatational potential
  // // Ua = 0.5 * lambda * <ln(J)>+ ^ 2
  // // Ui = 0.5 * lambda * <ln(J)>- ^ 2
  // ADReal lnJ_pos = J > 1 ? std::log(J) : 0;
  // ADReal lnJ_neg = J > 1 ? 0 : std::log(J);
  // ADReal Ua = 0.5 * lambda * lnJ_pos * lnJ_pos;
  // ADReal Ui = 0.5 * lambda * lnJ_neg * lnJ_neg;
  //
  // // isochoric potential
  // // Wa = G * sum(<e_i>+ ^ 2)
  // // Wi = G * sum(<e_i>- ^ 2)
  // // e_1 = 0.5 * ln(d_i)
  // // d_i are eigenvalues of b
  // ADRankTwoTensor eigvecs, etens, ha, hi;
  // std::vector<ADReal> eigvals;
  // b.symmetricEigenvaluesEigenvectors(eigvals, eigvecs);
  // for (unsigned int i = 0; i < eigvals.size(); i++)
  // {
  //   etens.vectorOuterProduct(eigvecs.column(i), eigvecs.column(i));
  //   if (eigvals[i] > 1)
  //     ha += 0.5 * std::log(eigvals[i]) * etens;
  //   else
  //     hi += 0.5 * std::log(eigvals[i]) * etens;
  // }
  // ADReal Wa = _G[_qp] * ha.doubleContraction(ha);
  // ADReal Wi = _G[_qp] * hi.doubleContraction(hi);
  //
  // // elastic potential
  // // E = D * (Ua + Wa) + Ui + Wi
  // _E_active[_qp] = Ua + Wa;
  // if (_E_active[_qp] > _hist_max_old[_qp])
  //   _hist_max[_qp] = _E_active[_qp];
  // else
  //   _hist_max[_qp] = _hist_max_old[_qp];
  // _E[_qp] = _D[_qp] * (Ua + Wa) + Ui + Wi;
  // _dE_dc[_qp] = _dD_dc[_qp] * (_hist ? _hist_max_old[_qp] : _E_active_old[_qp]);
  //
  // // Kirchhoff stress
  // ADRankTwoTensor tau_a = lambda * lnJ_pos * (ha + hi) + 2 * _G[_qp] * ha;
  // ADRankTwoTensor tau_i = lambda * lnJ_neg * (ha + hi) + 2 * _G[_qp] * hi;
  // ADRankTwoTensor tau = _D[_qp] * tau_a + tau_i;
  //
  // // cauchy stress
  // _stress[_qp] = tau / J;
}
