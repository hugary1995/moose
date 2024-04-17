//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Action.h"
#include "Factory.h"
#include "FEProblem.h"
#include "PropertyRegistry.h"

class SolidMechanicsPhysicsBase : public Action
{
public:
  static InputParameters validParams();

  SolidMechanicsPhysicsBase(const InputParameters & params);

  static MultiMooseEnum outputProperties();
  static MultiMooseEnum materialOutputOrders();
  static MultiMooseEnum materialOutputFamilies();

protected:
  void applyCommonParameters();
  void mergeMultiMooseEnum(const std::string & dest, const std::string & src);
  std::string getOutputPropertyName(const std::string & out, const std::string & key) const;

  template <typename... ParamType, typename... ValueType>
  bool
  addOutputMatPropRankTwo(const std::string & out,
                          const std::string & obj_type,
                          const std::vector<std::string> & param_names,
                          const std::map<std::string, std::tuple<ValueType...>> & param_val_table,
                          const std::string & basename);

  ///@{ Data for output generation
  static const std::vector<SolidMechanics::PropertyRegistry> _output_properties;
  static const std::map<std::string, std::tuple<std::string>> _output_r2_invariant_params;
  static const std::map<std::string, std::tuple<>> _output_r2_directional_params;
  static const std::map<std::string, std::tuple<std::string>> _output_r2_cylindrical_params;
  static const std::map<std::string, std::tuple<std::string>> _output_r2_spherical_params;
  static const std::map<std::string, std::tuple<unsigned int, unsigned int>>
      _output_r2_cartesian_params;
  static const std::map<
      std::string,
      std::tuple<bool, SolidMechanics::Rank, SolidMechanics::Symmetry, SolidMechanics::Type>>
      _output_property_restriction;
  ///@}

  const bool _use_ad;

private:
  template <size_t I, typename... ParamType, typename... ValueType>
  void applyOutputMatPropParameters(InputParameters & params,
                                    const std::vector<std::string> & param_names,
                                    const std::tuple<ValueType...> & param_vals) const;
};

template <typename... ParamType, typename... ValueType>
bool
SolidMechanicsPhysicsBase::addOutputMatPropRankTwo(
    const std::string & out,
    const std::string & obj_type,
    const std::vector<std::string> & param_names,
    const std::map<std::string, std::tuple<ValueType...>> & param_val_table,
    const std::string & basename)
{
  // Check if the output is prefixed/suffixed by any of the supported prefix/suffix
  for (auto & [key, param_vals] : param_val_table)
  {
    auto prefix = std::get<0>(_output_property_restriction.at(key));
    if (prefix && !std::equal(key.begin(), key.end(), out.begin()))
      continue;
    if (!prefix && !std::equal(key.rbegin(), key.rend(), out.rbegin()))
      continue;

    // Map from property alias to the material property name declared by the object.
    // For example, this maps from "strain" (specified by the user) to "total_strain" (declared in
    // the strain calculators).
    auto prop_name = getOutputPropertyName(out, key);

    // Retrieve valid parameters and apply action parameters where applicable
    auto ad_obj_type = (_use_ad ? "AD" : "") + obj_type;
    auto params = _factory.getValidParams(ad_obj_type);
    params.applyParameters(parameters());

    // These two properties are common
    params.template set<MaterialPropertyName>("rank_two_tensor") = basename + prop_name;
    params.template set<MaterialPropertyName>("property_name") = basename + out;

    // Apply additional parameter values specific to each obj type.
    // For example, (AD)RankTwoInvariant requires "invariant",
    //              (AD)RankTwoCartesianComponent requires "index_i" and "index_j".
    applyOutputMatPropParameters<0, ParamType...>(params, param_names, param_vals);

    // Finally, add the output helper material to the problem
    _problem->addMaterial(ad_obj_type, basename + out + '_' + name(), params);
    return true;
  }

  return false;
}

template <size_t I, typename... ParamType, typename... ValueType>
void
SolidMechanicsPhysicsBase::applyOutputMatPropParameters(
    InputParameters & params,
    const std::vector<std::string> & param_names,
    const std::tuple<ValueType...> & param_vals) const
{
  static_assert(sizeof...(ParamType) == sizeof...(ValueType));

  if constexpr (sizeof...(ParamType) > I)
  {
    typedef typename std::tuple_element<I, std::tuple<ParamType...>>::type T;
    params.set<T>(param_names[I]) = std::get<I>(param_vals);
    applyOutputMatPropParameters<I + 1, ParamType...>(params, param_names, param_vals);
  }
}
