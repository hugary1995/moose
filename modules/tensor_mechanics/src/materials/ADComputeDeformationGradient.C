//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ADComputeDeformationGradient.h"

registerADMooseObject("TensorMechanicsApp", ADComputeDeformationGradient);

InputParameters
ADComputeDeformationGradient::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Compute the deformation gradient");
  params.addRequiredCoupledVar(
      "displacements",
      "The displacements appropriate for the simulation geometry and coordinate system");
  params.addParam<std::string>("base_name",
                               "Optional parameter that allows the user to define "
                               "multiple mechanics material systems on the same "
                               "block, i.e. for multiple phases");
  params.addParam<std::vector<MaterialPropertyName>>(
      "eigenstrain_names", "List of eigenstrains to be applied in this strain calculation");
  params.suppressParameter<bool>("use_displaced_mesh");
  return params;
}

ADComputeDeformationGradient::ADComputeDeformationGradient(const InputParameters & parameters)
  : Material(parameters),
    _ndisp(coupledComponents("displacements")),
    _disp(3),
    _grad_disp(3),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _F(declareADProperty<RankTwoTensor>(_base_name + "deformation_gradient")),
    _eigenstrain_names(getParam<std::vector<MaterialPropertyName>>("eigenstrain_names")),
    _eigenstrains(_eigenstrain_names.size())
{
  for (unsigned int i = 0; i < _eigenstrains.size(); ++i)
  {
    _eigenstrain_names[i] = _base_name + _eigenstrain_names[i];
    _eigenstrains[i] = &getADMaterialProperty<RankTwoTensor>(_eigenstrain_names[i]);
  }

  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    _disp[i] = &adCoupledValue("displacements", i);
    _grad_disp[i] = &adCoupledGradient("displacements", i);
  }

  // set unused dimensions to zero
  for (unsigned i = _ndisp; i < 3; ++i)
  {
    _disp[i] = &_ad_zero;
    _grad_disp[i] = &_ad_grad_zero;
  }
}

ADReal
ADComputeDeformationGradient::computeQpOutOfPlaneGradDisp()
{
  if (getBlockCoordSystem() == Moose::COORD_RZ)
    if (!MooseUtils::absoluteFuzzyEqual(_q_point[_qp](0), 0.0))
      return (*_disp[0])[_qp] / _q_point[_qp](0);

  return 0.0;
}

void
ADComputeDeformationGradient::computeQpProperties()
{
  ADRankTwoTensor A((*_grad_disp[0])[_qp], (*_grad_disp[1])[_qp], (*_grad_disp[2])[_qp]);
  if (_ndisp == 2)
    A(2, 2) = computeQpOutOfPlaneGradDisp();

  _F[_qp] = A;
  _F[_qp].addIa(1.0);

  for (auto es : _eigenstrains)
  {
    ADRankTwoTensor eigvecs, E;
    std::vector<ADReal> eigvals;
    (*es)[_qp].symmetricEigenvaluesEigenvectors(eigvals, eigvecs);
    for (unsigned int i = 0; i < eigvals.size(); i++)
      eigvals[i] = std::exp(-eigvals[i]);
    E.fillFromInputVector(eigvals);
    E = eigvecs * E * eigvecs.transpose();
    _F[_qp] = _F[_qp] * E;
  }
}
