//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
#include "SymmetricIsotropicRankFourTensor.h"
#include "GuaranteeProvider.h"

template <bool is_ad>
class ComputeSymmetricIsotropicElasticityTensorTempl : public Material, public GuaranteeProvider
{
public:
  static InputParameters validParams();

  ComputeSymmetricIsotropicElasticityTensorTempl(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Base name of the material system
  const std::string _base_name;

  std::string _elasticity_tensor_name;

  GenericMaterialProperty<SymmetricIsotropicRankFourTensor, is_ad> & _elasticity_tensor;

  Real _K;
  Real _G;

  SymmetricIsotropicRankFourTensor _Cijkl;
};

typedef ComputeSymmetricIsotropicElasticityTensorTempl<false>
    ComputeSymmetricIsotropicElasticityTensor;
typedef ComputeSymmetricIsotropicElasticityTensorTempl<true>
    ADComputeSymmetricIsotropicElasticityTensor;
