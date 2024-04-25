//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SolidMechanicsPhysicsCommonParameters.h"
#include "ComputeFiniteStrain.h"

using namespace SolidMechanics;

// Map aux variable name prefixes/suffixes to applicable tensor queries
// clang-format off
const std::map<std::string, PropertyRegistryEntry>
    SolidMechanicsPhysicsBase::_output_property_restriction = {
      // RankTwoInvariant
      {"vonmises_",      {"", "", Rank::TWO, Symmetry::ANY,                                      Type::STRESS_LIKE}},
      {"hydrostatic_",   {"", "", Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"l2norm_",        {"", "", Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"volumetric_",    {"", "", Rank::TWO, Symmetry::ANY,                                      Type::STRAIN_LIKE}},
      {"firstinv_",      {"", "", Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"secondinv_",     {"", "", Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"thirdinv_",      {"", "", Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"triaxiality_",   {"", "", Rank::TWO, Symmetry::ANY,                                      Type::STRESS_LIKE}},
      {"maxshear_",      {"", "", Rank::TWO, Symmetry::SYMMETRIC | Symmetry::POSSIBLY_SYMMETRIC, Type::ANY}},
      {"intensity_",     {"", "", Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"max_principal_", {"", "", Rank::TWO, Symmetry::SYMMETRIC | Symmetry::POSSIBLY_SYMMETRIC, Type::ANY}},
      {"mid_principal_", {"", "", Rank::TWO, Symmetry::SYMMETRIC | Symmetry::POSSIBLY_SYMMETRIC, Type::ANY}},
      {"min_principal_", {"", "", Rank::TWO, Symmetry::SYMMETRIC | Symmetry::POSSIBLY_SYMMETRIC, Type::ANY}},
      // RankTwoDirectionalComponent
      {"directional_",   {"", "", Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      // RankTwoCartesianComponent
      {"_xx",            {"", "", Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"_xy",            {"", "", Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"_xz",            {"", "", Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"_yx",            {"", "", Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"_yy",            {"", "", Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"_yz",            {"", "", Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"_zx",            {"", "", Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"_zy",            {"", "", Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"_zz",            {"", "", Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      // RankTwo[Cylindrical,Spherical]Component
      {"axial_",         {"", "", Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"hoop_",          {"", "", Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"radial_",        {"", "", Rank::TWO, Symmetry::ANY,                                      Type::ANY}}
    };
// clang-format on

InputParameters
SolidMechanicsPhysicsCommonParameters::validParams()
{
  /////////////////////////////////////////////////////////////////////////////
  // Parameters common to the old and the new systems
  /////////////////////////////////////////////////////////////////////////////
  params.addRequiredParam<std::vector<VariableName>>(
      "displacements", "The nonlinear displacement variables for the problem");
  params.addParam<std::vector<VariableName>>("temperature", "The temperature");
  params.addParam<bool>("add_variables", false, "Add the displacement variables");

  MooseEnum strainType("SMALL FINITE", "SMALL");
  params.addParam<MooseEnum>("strain", strainType, "Strain formulation");

  params.addParam<std::string>("base_name", "Material property base name");
  params.addParam<std::string>(
      "strain_base_name",
      "The base name used for the strain. If not provided, it will be set equal to base_name");
  params.addParam<bool>(
      "volumetric_locking_correction", false, "Flag to correct volumetric locking");
  params.addParam<std::vector<MaterialPropertyName>>(
      "eigenstrain_names", {}, "List of eigenstrains to be applied in this strain calculation");
  params.addParam<bool>("use_automatic_differentiation",
                        false,
                        "Flag to use automatic differentiation (AD) objects when possible");
  params.addParam<bool>("automatic_eigenstrain_names",
                        false,
                        "Collects all material eigenstrains and passes to required strain "
                        "calculator within TMA internally.");

  // Advanced
  params.addParam<std::vector<AuxVariableName>>("save_in", {}, "The displacement residuals");
  params.addParam<std::vector<AuxVariableName>>(
      "diag_save_in", {}, "The displacement diagonal preconditioner terms");
  params.addParam<std::vector<TagName>>(
      "extra_vector_tags",
      "The tag names for extra vectors that residual data should be saved into");
  params.addParam<std::vector<TagName>>("absolute_value_vector_tags",
                                        "The tag names for extra vectors that the absolute value "
                                        "of the residual should be accumulated into");
  params.addParam<std::vector<Real>>(
      "scaling",
      {},
      "The scaling to apply to the displacement variables. If one number is provided, the same "
      "scaling is applied to all displacement variables. If multiple scaling numbers are "
      "specified, the scaling numbers are applied to each of the displacement variable.");

  params.addParamNamesToGroup(
      "save_in diag_save_in extra_vector_tags absolute_value_vector_tags scaling", "Advanced");

  // Planar Formulation
  MooseEnum planarFormulationType("NONE WEAK_PLANE_STRESS PLANE_STRAIN GENERALIZED_PLANE_STRAIN",
                                  "NONE");
  params.addParam<MooseEnum>(
      "planar_formulation", planarFormulationType, "Out-of-plane stress/strain formulation");
  params.addParam<VariableName>("scalar_out_of_plane_strain",
                                "Scalar variable for the out-of-plane strain (in y "
                                "direction for 1D Axisymmetric or in z direction for 2D "
                                "Cartesian problems)");
  params.addParam<VariableName>("out_of_plane_strain",
                                "Variable for the out-of-plane strain for plane stress models");
  MooseEnum outOfPlaneDirection("x y z", "z");
  params.addParam<MooseEnum>(
      "out_of_plane_direction", outOfPlaneDirection, "The direction of the out-of-plane strain.");
  params.addParamNamesToGroup(
      "planar_formulation scalar_out_of_plane_strain out_of_plane_strain out_of_plane_direction",
      "Planar formulation");

  // Output
  params.addParam<MultiMooseEnum>("generate_output",
                                  SolidMechanicsPhysicsCommonParameters::outputProperties(),
                                  "Add scalar quantity output for stress and/or strain");
  params.addParam<MultiMooseEnum>(
      "material_output_order",
      SolidMechanicsPhysicsCommonParameters::materialOutputOrders(),
      "Specifies the order of the FE shape function to use for this variable.");
  params.addParam<MultiMooseEnum>(
      "material_output_family",
      SolidMechanicsPhysicsCommonParameters::materialOutputFamilies(),
      "Specifies the family of FE shape functions to use for this variable.");
  params.addParam<Point>(
      "cylindrical_axis_point1",
      "Starting point for direction of axis of rotation for cylindrical stress/strain.");
  params.addParam<Point>(
      "cylindrical_axis_point2",
      "Ending point for direction of axis of rotation for cylindrical stress/strain.");
  params.addParam<Point>("spherical_center_point",
                         "Center point of the spherical coordinate system.");
  params.addParam<Point>("direction", "Direction stress/strain is calculated in");
  params.addParamNamesToGroup(
      "generate_output material_output_order material_output_family "
      "cylindrical_axis_point1 cylindrical_axis_point2 spherical_center_point direction",
      "Output");
  params.addParam<bool>("verbose", false, "Display extra information.");

  /////////////////////////////////////////////////////////////////////////////
  // Parameters specific to the old system
  /////////////////////////////////////////////////////////////////////////////
  params.addParam<bool>("incremental",
                        "Use incremental or total strain (if not explicitly specified this "
                        "defaults to incremental for finite strain and total for small strain)");
  params.addParam<bool>(
      "use_finite_deform_jacobian", false, "Jacobian for corrotational finite strain");

  params.addParam<MaterialPropertyName>(
      "global_strain",
      "Name of the global strain material to be applied in this strain calculation. "
      "The global strain tensor is constant over the whole domain and allows visualization "
      "of the deformed shape with the periodic BC");
  params.addParam<MooseEnum>("decomposition_method",
                             ComputeFiniteStrain::decompositionType(),
                             "Methods to calculate the finite strain and rotation increments");
  params.addDeprecatedParam<FunctionName>(
      "out_of_plane_pressure",
      "Function used to prescribe pressure (applied toward the body) in the out-of-plane direction "
      "(y for 1D Axisymmetric or z for 2D Cartesian problems)",
      "This has been replaced by 'out_of_plane_pressure_function'");
  params.addParam<FunctionName>(
      "out_of_plane_pressure_function",
      "Function used to prescribe pressure (applied toward the body) in the out-of-plane direction "
      "(y for 1D Axisymmetric or z for 2D Cartesian problems)");
  params.addParam<Real>(
      "pressure_factor",
      "Scale factor applied to prescribed out-of-plane pressure (both material and function)");
  params.addParam<MaterialPropertyName>("out_of_plane_pressure_material",
                                        "0",
                                        "Material used to prescribe pressure (applied toward the "
                                        "body) in the out-of-plane direction");
  params.addParamNamesToGroup("out_of_plane_pressure out_of_plane_pressure_material "
                              "out_of_plane_pressure_function pressure_factor",
                              "Planar formulation");

  /////////////////////////////////////////////////////////////////////////////
  // Parameters specific to the new system
  /////////////////////////////////////////////////////////////////////////////

  return params;
}

MultiMooseEnum
SolidMechanicsPhysicsCommonParameters::outputProperties()
{
  std::string options = "";

  // Rank-0 (scalar-valued) properties
  for (auto prop : PropertyRegistry::query().rank(Rank::ZERO).get())
    options += prop.alias + ' ';

  // Higher rank properties will need one additional layer of postprocessing so that we can write
  // them into scalar-valued aux variables.
  // The type of postprocessing is specified by the prefix/suffix.
  for (auto & [key, restriction] : output_property_restriction)
    for (auto prop : PropertyRegistry::query().filter(restriction).get())
    {
      auto prefix = isPrefix(key);
      options += prefix ? key + prop.alias + ' ' : prop.alias + key + ' ';
    }

  return MultiMooseEnum(options, "", /*allow_out_of_range=*/true);
}

MultiMooseEnum
SolidMechanicsPhysicsCommonParameters::materialOutputOrders()
{
  return AddAuxVariableAction::getAuxVariableOrders().getRawNames();
}

MultiMooseEnum
SolidMechanicsPhysicsCommonParameters::materialOutputFamilies()
{
  return MultiMooseEnum("MONOMIAL LAGRANGE");
}

bool
SolidMechanicsPhysicsBase::isPrefix(const std::string & key)
{
  mooseAssert(key.front() == '_' || key.back() == '_',
              "Key in output property table must begin or end with character '_'.");

  if (key.front() == '_')
    return false;

  return true;
}
