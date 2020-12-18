//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeSymmetricIsotropicElasticityTensor.h"

registerMooseObject("TensorMechanicsApp", ComputeSymmetricIsotropicElasticityTensor);
registerMooseObject("TensorMechanicsApp", ADComputeSymmetricIsotropicElasticityTensor);

template <bool is_ad>
InputParameters
ComputeSymmetricIsotropicElasticityTensorTempl<is_ad>::validParams()
{
  InputParameters params = Material::validParams();
  params.addParam<std::string>("base_name",
                               "Optional parameter that allows the user to define "
                               "multiple mechanics material systems on the same "
                               "block, i.e. for multiple phases");
  params.addRequiredParam<Real>("bulk_modulus", "The bulk modulus for the material.");
  params.addRequiredParam<Real>("shear_modulus", "The shear modulus of the material.");
  return params;
}

template <bool is_ad>
ComputeSymmetricIsotropicElasticityTensorTempl<
    is_ad>::ComputeSymmetricIsotropicElasticityTensorTempl(const InputParameters & parameters)
  : Material(parameters),
    GuaranteeProvider(this),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(
        declareGenericProperty<SymmetricIsotropicRankFourTensor, is_ad>(_elasticity_tensor_name)),
    _K(getParam<Real>("bulk_modulus")),
    _G(getParam<Real>("shear_modulus"))
{
  issueGuarantee(_elasticity_tensor_name, Guarantee::ISOTROPIC);
  _Cijkl(0) = 3. * _K;
  _Cijkl(1) = 2. * _G;
}

template <bool is_ad>
void
ComputeSymmetricIsotropicElasticityTensorTempl<is_ad>::computeQpProperties()
{
  _elasticity_tensor[_qp] = _Cijkl;
}

template class ComputeSymmetricIsotropicElasticityTensorTempl<false>;
template class ComputeSymmetricIsotropicElasticityTensorTempl<true>;
