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
#include "DerivativeMaterialPropertyNameInterface.h"

class ADComputeHyperelasticPFFractureStress : public Material,
                                              public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  ADComputeHyperelasticPFFractureStress(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  const std::string _base_name;

  ADMaterialProperty<RankTwoTensor> & _stress;

  const ADMaterialProperty<RankTwoTensor> & _deformation_gradient;

  /// Material property for elastic energy
  ADMaterialProperty<Real> & _E;

  /// Derivative of elastic energy w.r.t damage variable
  ADMaterialProperty<Real> & _dE_dc;

  /// Material property for energetic degradation function
  const ADMaterialProperty<Real> & _D;

  /// Derivative of degradation function w.r.t damage variable
  const ADMaterialProperty<Real> & _dD_dc;

  const ADMaterialProperty<Real> & _K;
  const ADMaterialProperty<Real> & _G;

  ADMaterialProperty<Real> & _E_active;
  const MaterialProperty<Real> & _E_active_old;

  ADMaterialProperty<Real> & _hist_max;
  const MaterialProperty<Real> & _hist_max_old;
  const bool _hist;
};
