//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LevelSetBiMaterialBase.h"
#include "DerivativeMaterialPropertyNameInterface.h"

/**
 * Compute a Real material property for bi-materials problem (consisting of two
 * different materials) defined by a level set function
 */
template <bool is_ad>
class LevelSetBiMaterialRealDerivativesTempl : public LevelSetBiMaterialBase,
                                               public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  LevelSetBiMaterialRealDerivativesTempl(const InputParameters & parameters);

protected:
  virtual void assignQpPropertiesForLevelSetPositive() override;
  virtual void assignQpPropertiesForLevelSetNegative() override;

  std::string _arg_name;

  /// Global Real material property (switch bi-material diffusion coefficient based on level set values)
  GenericMaterialProperty<Real, is_ad> & _prop;
  GenericMaterialProperty<Real, is_ad> & _dprop_dc;

  /// Real Material properties for the two separate materials in the bi-material system
  std::pair<const GenericMaterialProperty<Real, is_ad> &,
            const GenericMaterialProperty<Real, is_ad> &>
      _biprop;
  std::pair<const GenericMaterialProperty<Real, is_ad> &,
            const GenericMaterialProperty<Real, is_ad> &>
      _dbiprop_dc;
};

typedef LevelSetBiMaterialRealDerivativesTempl<false> LevelSetBiMaterialRealDerivatives;
typedef LevelSetBiMaterialRealDerivativesTempl<true> ADLevelSetBiMaterialRealDerivatives;
