//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LevelSetBiMaterialRealDerivatives.h"

registerMooseObject("XFEMApp", LevelSetBiMaterialRealDerivatives);
registerMooseObject("XFEMApp", ADLevelSetBiMaterialRealDerivatives);

template <bool is_ad>
InputParameters
LevelSetBiMaterialRealDerivativesTempl<is_ad>::validParams()
{
  InputParameters params = LevelSetBiMaterialBase::validParams();
  params.addClassDescription(
      "Compute a Real material property for bi-materials problem (consisting of two "
      "different materials) defined by a level set function.");
  params.addRequiredCoupledVar("args", "function arguments");
  return params;
}

template <bool is_ad>
LevelSetBiMaterialRealDerivativesTempl<is_ad>::LevelSetBiMaterialRealDerivativesTempl(
    const InputParameters & parameters)
  : LevelSetBiMaterialBase(parameters),
    _arg_name(this->getVar("args", 0)->name()),
    _prop(declareGenericProperty<Real, is_ad>(_base_name + _prop_name)),
    _dprop_dc(declareGenericProperty<Real, is_ad>(
        derivativePropertyNameFirst(_base_name + _prop_name, _arg_name))),
    _biprop(getGenericMaterialProperty<Real, is_ad>(_positive_base_name + _prop_name),
            getGenericMaterialProperty<Real, is_ad>(_negative_base_name + _prop_name)),
    _dbiprop_dc(getGenericMaterialProperty<Real, is_ad>(
                    derivativePropertyNameFirst(_positive_base_name + _prop_name, _arg_name)),
                getGenericMaterialProperty<Real, is_ad>(
                    derivativePropertyNameFirst(_negative_base_name + _prop_name, _arg_name)))
{
}

template <bool is_ad>
void
LevelSetBiMaterialRealDerivativesTempl<is_ad>::assignQpPropertiesForLevelSetPositive()
{
  _prop[_qp] = _biprop.first[_qp];
  _dprop_dc[_qp] = _dbiprop_dc.first[_qp];
}

template <bool is_ad>
void
LevelSetBiMaterialRealDerivativesTempl<is_ad>::assignQpPropertiesForLevelSetNegative()
{
  _prop[_qp] = _biprop.second[_qp];
  _dprop_dc[_qp] = _dbiprop_dc.second[_qp];
}

template class LevelSetBiMaterialRealDerivativesTempl<false>;
template class LevelSetBiMaterialRealDerivativesTempl<true>;
