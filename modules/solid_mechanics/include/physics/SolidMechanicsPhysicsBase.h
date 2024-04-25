//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "SolidMechanicsPhysicsSubBlock.h"
#include "Factory.h"
#include "FEProblem.h"

class SolidMechanicsPhysicsBase : public SolidMechanicsPhysicsSubBlock
{
public:
  static InputParameters validParams();

  SolidMechanicsPhysicsBase(const InputParameters & params);

protected:
  /// Whether the displaced mesh should be used.
  virtual bool useDisplacedMesh() const = 0;

  std::string getOutputPropertyName(const std::string & out, const std::string & key) const;

  template <typename... ParamType, typename... ValueType>
  bool
  addOutputMatPropRankTwo(const std::string & out,
                          const std::string & obj_type,
                          const std::vector<std::string> & param_names,
                          const std::map<std::string, std::tuple<ValueType...>> & param_val_table,
                          const std::string & basename);

  ///@{ Data for output generation
  static const std::map<std::string, std::tuple<std::string>> _output_r2_invariant_params;
  static const std::map<std::string, std::tuple<>> _output_r2_directional_params;
  static const std::map<std::string, std::tuple<std::string>> _output_r2_cylindrical_params;
  static const std::map<std::string, std::tuple<std::string>> _output_r2_spherical_params;
  static const std::map<std::string, std::tuple<unsigned int, unsigned int>>
      _output_r2_cartesian_params;
  ///@}

  /// Displacement variable names
  std::vector<VariableName> _displacements;

  /// Number of displacement variables
  unsigned int _ndisp;

  /// Scaling for each of the displacement variable
  std::vector<Real> _scaling;

  ///@{ Residual
  const std::vector<AuxVariableName> _save_in;
  const std::vector<AuxVariableName> _diag_save_in;
  ///@}

  /// Strain formulation
  const MooseEnum _strain;

  /// Incremental or total material update
  const bool _incremental;

  /// Whether to use the displaced mesh
  const bool _use_displaced_mesh;

  /// Planar formulation
  const MooseEnum _planar_formulation;

  /// Out-of-plane direction
  const MooseEnum _out_of_plane_direction;

  /// Base name for the current block (appended with '_' if non-empty)
  const std::string _base_name;

  /// Output materials to generate scalar output quantities
  MultiMooseEnum _generate_output;
  MultiMooseEnum _material_output_order;
  MultiMooseEnum _material_output_family;

  /// automatically gather names of eigenstrain tensors provided by simulation objects
  const bool _auto_eigenstrain;

  std::vector<MaterialPropertyName> _eigenstrain_names;

  /// Verbosity
  const bool _verbose;

  /// Whether to use automatic differentiation
  const bool _use_ad;

private:
  template <size_t I, typename... ParamType, typename... ValueType>
  void applyOutputMatPropParameters(InputParameters & params,
                                    const std::vector<std::string> & param_names,
                                    const std::tuple<ValueType...> & param_vals) const;

  static const std::map<Strain, bool> _incremental_default;
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
    auto prefix = SolidMechanicsPhysicsBase::isPrefix(key);
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
