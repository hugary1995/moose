//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SolidMechanicsPhysicsBase.h"
#include "SolidMechanicsPhysicsCommon.h"
#include "ActionWarehouse.h"
#include "AddAuxVariableAction.h"
#include "ComputeFiniteStrain.h"
#include "MooseApp.h"
#include "InputParameterWarehouse.h"

using namespace SolidMechanics;

// Map aux variable name prefixes to enum option used in RankTwoInvariant
// clang-format off
const std::map<std::string, std::tuple<std::string>>
    SolidMechanicsPhysicsBase::_output_r2_invariant_params = {
        {"vonmises_",      {"VonMisesStress"}},
        {"hydrostatic_",   {"Hydrostatic"}},
        {"l2norm_",        {"L2norm"}},
        {"volumetric_",    {"VolumetricStrain"}},
        {"firstinv_",      {"FirstInvariant"}},
        {"secondinv_",     {"SecondInvariant"}},
        {"thirdinv_",      {"ThirdInvariant"}},
        {"triaxiality_",   {"TriaxialityStress"}},
        {"maxshear_",      {"MaxShear"}},
        {"intensity_",     {"StressIntensity"}},
        {"max_principal_", {"MaxPrincipal"}},
        {"mid_principal_", {"MidPrincipal"}},
        {"min_principal_", {"MinPrincipal"}}
    };
// clang-format on

// Map aux variable name prefixes to enum option used in RankTwoDirectionalComponent
// clang-format off
const std::map<std::string, std::tuple<>>
    SolidMechanicsPhysicsBase::_output_r2_directional_params = {
        {"directional_",  {}}
    };
// clang-format on

// Map aux variable name prefixes to enum option used in RankTwoCylindricalComponent
// clang-format off
const std::map<std::string, std::tuple<std::string>>
    SolidMechanicsPhysicsBase::_output_r2_cylindrical_params = {
        {"axial_",  {"AxialStress"}},
        {"hoop_",   {"HoopStress"}},
        {"radial_", {"RadialStress"}}
    };
// clang-format on

// Map aux variable name prefixes to enum option used in RankTwoSphericalComponent
// clang-format off
const std::map<std::string, std::tuple<std::string>>
    SolidMechanicsPhysicsBase::_output_r2_spherical_params = {
        {"hoop_",   {"HoopStress"}},
        {"radial_", {"RadialStress"}}
    };
// clang-format on

// Map aux variable name suffixes to indices used in RankTwoCartesianComponent
// clang-format off
const std::map<std::string, std::tuple<unsigned int, unsigned int>>
    SolidMechanicsPhysicsBase::_output_r2_cartesian_params = {
        {"_xx", {0, 0}}, {"_xy", {0, 1}}, {"_xz", {0, 2}},
        {"_yx", {1, 0}}, {"_yy", {1, 1}}, {"_yz", {1, 2}},
        {"_zx", {2, 0}}, {"_zy", {2, 1}}, {"_zz", {2, 2}}
    };
// clang-format on

// Default choice on whether to perform incremental material update
const std::map<Strain, bool> _incremental_default = {{Strain::Small, false},
                                                     {Strain::Finite, true}};

InputParameters
SolidMechanicsPhysicsBase::validParams()
{
  InputParameters params = SolidMechanicsPhysicsSubBlock::validParams();
  return params;
}

SolidMechanicsPhysicsBase::SolidMechanicsPhysicsBase(const InputParameters & parameters)
  : SolidMechanicsPhysicsSubBlock(parameters),
    _displacements(getParam<std::vector<VariableName>>("displacements")),
    _ndisp(_displacements.size()),
    _scaling(getParam<std::vector<Real>>("scaling")),
    _save_in(getParam<std::vector<AuxVariableName>>("save_in")),
    _diag_save_in(getParam<std::vector<AuxVariableName>>("diag_save_in")),
    _strain(getParam<MooseEnum>("strain")),
    _incremental(getParam<bool>("incremental")),
    _use_displaced_mesh(useDisplacedMesh()),
    _planar_formulation(getParam<MooseEnum>("planar_formulation")),
    _out_of_plane_direction(getParam<MooseEnum>("out_of_plane_direction")),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _generate_output(getParam<MultiMooseEnum>("generate_output")),
    _material_output_order(getParam<MultiMooseEnum>("material_output_order")),
    _material_output_family(getParam<MultiMooseEnum>("material_output_family")),
    _auto_eigenstrain(getParam<bool>("automatic_eigenstrain_names")),
    _eigenstrain_names(getParam<std::vector<MaterialPropertyName>>("eigenstrain_names")),
    _verbose(getParam<bool>("verbose")),
    _use_ad(getParam<bool>("use_automatic_differentiation"))
{
  // Error checks
  checkVLC();
  checkScaling();
  checkSaveIn();
  checkPlanarFormulation();
  checkOutputOrder();
  checkOutputFamily();
}

void
SolidMechanicsPhysicsBase::checkVLC() const
{
  // Error if volumetric locking correction is true for 1D problems
  if (_ndisp == 1 && getParam<bool>("volumetric_locking_correction"))
    mooseError("Volumetric locking correction should be set to false for 1D problems.");
}

void
SolidMechanicsPhysicsBase::checkScaling()
{
  if (!isParamSetByUser("scaling"))
    continue;

  if (!getParam<bool>("add_variables"))
    paramError("scaling",
               "The scaling parameter has no effect unless add_variables is set to true. Did you "
               "mean to set 'add_variables = true'?");

  // check number of scaling numbers
  if (_scaling.size() != 1 && _scaling.size() != _ndisp)
    paramError("scaling",
               "The number of variable scaling numbers must be: 1 to use the same scaling for all "
               "displacement variables or the same size as the number of displacement variables.");

  // For only one order, make all orders the same magnitude
  if (_scaling.size() == 1)
    _scaling = std::vector<Real>(_ndisp, _scaling[0]);
}

void
SolidMechanicsPhysicsBase::checkSaveIn() const
{
  if (_save_in.size() != 0 && _save_in.size() != _ndisp)
    paramError("save_in",
               "Number of save_in variables should equal to the number of displacement variables ",
               _ndisp);

  if (_diag_save_in.size() != 0 && _diag_save_in.size() != _ndisp)
    paramError(
        "diag_save_in",
        "Number of diag_save_in variables should equal to the number of displacement variables ",
        _ndisp);
}

void
SolidMechanicsPhysicsBase::checkPlanarFormulation() const
{
  if (_planar_formulation == "NONE")
  {
    if (isParamValid("scalar_out_of_plane_strain"))
      paramError("scalar_out_of_plane_strain",
                 "scalar_out_of_plane_strain has no effect because planar_formulation=NONE");
    if (isParamValid("out_of_plane_strain"))
      paramError("out_of_plane_strain",
                 "out_of_plane_strain is ignored because planar_formulation=NONE");
    if (isParamSetByUser("out_of_plane_direction"))
      paramError("out_of_plane_direction",
                 "out_of_plane_direction is ignored because planar_formulation=NONE");
  }
  else if (_planar_formulation == "WEAK_PLANE_STRESS")
  {
    if (isParamValid("scalar_out_of_plane_strain"))
      paramError(
          "scalar_out_of_plane_strain",
          "scalar_out_of_plane_strain has no effect because planar_formulation=WEAK_PLANE_STRESS");
    if (!isParamValid("out_of_plane_strain"))
      paramError("out_of_plane_strain",
                 "out_of_plane_strain is required with planar_formulation=WEAK_PLANE_STRESS");
  }
  else if (_planar_formulation == "PLANE_STRAIN")
  {
    if (isParamValid("scalar_out_of_plane_strain"))
      paramError(
          "scalar_out_of_plane_strain",
          "scalar_out_of_plane_strain has no effect because planar_formulation=PLANE_STRAIN");
    if (isParamValid("out_of_plane_strain"))
      paramError("out_of_plane_strain",
                 "out_of_plane_strain has no effect because planar_formulation=PLANE_STRAIN");
    if (isParamSetByUser("out_of_plane_direction"))
      paramError("out_of_plane_direction",
                 "out_of_plane_direction has no effect because planar_formulation=PLANE_STRAIN");
  }
  else if (_planar_formulation == "GENERALIZED_PLANE_STRAIN")
  {
    if (!isParamValid("scalar_out_of_plane_strain"))
      paramError("scalar_out_of_plane_strain",
                 "scalar_out_of_plane_strain is required with "
                 "planar_formulation=GENERALIZED_PLANE_STRAIN");
    if (isParamValid("out_of_plane_strain"))
      paramError(
          "out_of_plane_strain",
          "out_of_plane_strain has no effect because planar_formulation=GENERALIZED_PLANE_STRAIN");
  }
}

void
SolidMechanicsPhysicsBase::checkOutputOrder()
{
  // Ensure material output order and family vectors are same size as generate output

  // check number of supplied orders and families
  if (_material_output_order.size() > 1 && _material_output_order.size() < _generate_output.size())
    paramError("material_output_order",
               "The number of orders assigned to material outputs must be: 0 to be assigned "
               "CONSTANT; 1 to assign all outputs the same value, or the same size as the number "
               "of generate outputs listed.");

  // if no value was provided, chose the default CONSTANT
  if (_material_output_order.size() == 0)
    _material_output_order.push_back("CONSTANT");

  // For only one order, make all orders the same magnitude
  if (_material_output_order.size() == 1)
    _material_output_order =
        std::vector<std::string>(_generate_output.size(), _material_output_order[0]);
}

void
SolidMechanicsPhysicsBase::checkOutputFamily()
{
  if (_material_output_family.size() > 1 &&
      _material_output_family.size() < _generate_output.size())
    paramError("material_output_family",
               "The number of families assigned to material outputs must be: 0 to be assigned "
               "MONOMIAL; 1 to assign all outputs the same value, or the same size as the number "
               "of generate outputs listed.");

  // if no value was provided, chose the default MONOMIAL
  if (_material_output_family.size() == 0)
    _material_output_family.push_back("MONOMIAL");

  // For only one family, make all families that value
  if (_material_output_family.size() == 1)
    _material_output_family =
        std::vector<std::string>(_generate_output.size(), _material_output_family[0]);
}

std::string
SolidMechanicsPhysicsBase::getOutputPropertyName(const std::string & out,
                                                 const std::string & key) const
{
  auto prefix = SolidMechanicsPhysicsBase::isPrefix(key);
  auto alias = prefix ? out.substr(key.length())                    // remove prefix
                      : out.substr(0, out.length() - key.length()); // remove suffix
  auto candidates =
      PropertyRegistry::query().filter(_output_property_restriction.at(key)).alias(alias).get();

  // If no candidate is found, treat the alias as the name of a custom tensor property
  if (candidates.empty())
    return alias;

  mooseAssert(candidates.size() == 1,
              "Internal error: Requested output yields multiple candidates.");
  return candidates[0].name;
}
