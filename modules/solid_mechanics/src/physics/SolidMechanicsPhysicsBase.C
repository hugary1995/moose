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

using namespace SolidMechanicsMaterialProperty;

// Register material properties name for outputting purposes.
// clang-format off
std::vector<Registry> SolidMechanicsPhysicsBase::_output_properties = {
  {"total_strain",             "strain",                   Rank::TWO,  Symmetry::SYMMETRIC,          Type::STRAIN_LIKE},
  {"mechanical_strain",        "mechanical_strain",        Rank::TWO,  Symmetry::POSSIBLY_SYMMETRIC, Type::STRAIN_LIKE},
  {"stress",                   "stress",                   Rank::TWO,  Symmetry::POSSIBLY_SYMMETRIC, Type::STRESS_LIKE},
  {"cauchy_stress",            "cauchy_stress",            Rank::TWO,  Symmetry::POSSIBLY_SYMMETRIC, Type::STRESS_LIKE},
  {"deformation_gradient",     "deformation_gradient",     Rank::TWO,  Symmetry::NONE,               Type::STRAIN_LIKE},
  {"pk1_stress",               "pk1_stress",               Rank::TWO,  Symmetry::NONE,               Type::STRESS_LIKE},
  {"pk2_stress",               "pk2_stress",               Rank::TWO,  Symmetry::POSSIBLY_SYMMETRIC, Type::STRESS_LIKE},
  {"small_stress",             "small_stress",             Rank::TWO,  Symmetry::POSSIBLY_SYMMETRIC, Type::STRESS_LIKE},
  {"elastic_strain",           "elastic_strain",           Rank::TWO,  Symmetry::POSSIBLY_SYMMETRIC, Type::STRAIN_LIKE},
  {"plastic_strain",           "plastic_strain",           Rank::TWO,  Symmetry::POSSIBLY_SYMMETRIC, Type::STRAIN_LIKE},
  {"creep_strain",             "creep_strain",             Rank::TWO,  Symmetry::POSSIBLY_SYMMETRIC, Type::STRAIN_LIKE},
  {"effective_plastic_strain", "effective_plastic_strain", Rank::ZERO, Symmetry::NONE,               Type::STRAIN_LIKE},
  {"effective_creep_strain",   "effective_creep_strain",   Rank::ZERO, Symmetry::NONE,               Type::STRAIN_LIKE}
};
// clang-format on

// Map aux variable name prefixes to applicable tensor queries
// clang-format off
std::map<std::string, std::tuple<bool, Rank, Symmetry, Type>>
    SolidMechanicsPhysicsBase::_output_property_restriction = {
      // RankTwoInvariant
      {"vonmises",      {true,  Rank::TWO, Symmetry::ANY,                                      Type::STRESS_LIKE}},
      {"hydrostatic",   {true,  Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"l2norm",        {true,  Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"volumetric",    {true,  Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"firstinv",      {true,  Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"secondinv",     {true,  Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"thirdinv",      {true,  Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"triaxiality",   {true,  Rank::TWO, Symmetry::ANY,                                      Type::STRESS_LIKE}},
      {"maxshear",      {true,  Rank::TWO, Symmetry::SYMMETRIC | Symmetry::POSSIBLY_SYMMETRIC, Type::ANY}},
      {"intensity",     {true,  Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"max_principal", {true,  Rank::TWO, Symmetry::SYMMETRIC | Symmetry::POSSIBLY_SYMMETRIC, Type::ANY}},
      {"mid_principal", {true,  Rank::TWO, Symmetry::SYMMETRIC | Symmetry::POSSIBLY_SYMMETRIC, Type::ANY}},
      {"min_principal", {true,  Rank::TWO, Symmetry::SYMMETRIC | Symmetry::POSSIBLY_SYMMETRIC, Type::ANY}},
      // RankTwoDirectionalComponent
      {"directional",   {true,  Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      // RankTwoCartesianComponent
      {"xx",            {false, Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"xy",            {false, Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"xz",            {false, Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"yx",            {false, Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"yy",            {false, Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"yz",            {false, Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"zx",            {false, Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"zy",            {false, Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"zz",            {false, Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      // RankTwo[Cylindrical,Spherical]Component
      {"axial",         {true,  Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"hoop",          {true,  Rank::TWO, Symmetry::ANY,                                      Type::ANY}},
      {"radial",        {true,  Rank::TWO, Symmetry::ANY,                                      Type::ANY}}
    };
// clang-format on

// Map aux variable name prefixes to enum option used in RankTwoInvariant
// clang-format off
const std::map<std::string, std::tuple<MooseEnum>> SolidMechanicsPhysicsBase::_output_r2_invariant_params = {
    {"vonmises",      {"VonMisesStress"}},
    {"hydrostatic",   {"Hydrostatic"}},
    {"l2norm",        {"L2norm"}},
    {"volumetric",    {"VolumetricStrain"}},
    {"firstinv",      {"FirstInvariant"}},
    {"secondinv",     {"SecondInvariant"}},
    {"thirdinv",      {"ThirdInvariant"}},
    {"triaxiality",   {"TriaxialityStress"}},
    {"maxshear",      {"MaxShear"}},
    {"intensity",     {"StressIntensity"}},
    {"max_principal", {"MaxPrincipal"}},
    {"mid_principal", {"MidPrincipal"}},
    {"min_principal", {"MinPrincipal"}}
};
// clang-format on

// Map aux variable name prefixes to enum option used in RankTwoDirectionalComponent
// clang-format off
const std::map<std::string, std::tuple<MooseEnum>> SolidMechanicsPhysicsBase::_output_r2_directional_params = {
    {"directional",  {}}
};
// clang-format on

// Map aux variable name prefixes to enum option used in RankTwoCylindricalComponent
// clang-format off
const std::map<std::string, std::tuple<MooseEnum>> SolidMechanicsPhysicsBase::_output_r2_cylindrical_params = {
    {"axial",  {"AxialStress"}},
    {"hoop",   {"HoopStress"}},
    {"radial", {"RadialStress"}}
};
// clang-format on

// Map aux variable name prefixes to enum option used in RankTwoSphericalComponent
// clang-format off
const std::map<std::string, std::tuple<MooseEnum>> SolidMechanicsPhysicsBase::_output_r2_spherical_params = {
    {"hoop",   {"HoopStress"}},
    {"radial", {"RadialStress"}}
};
// clang-format on

// Map aux variable name suffixes to indices used in RankTwoCartesianComponent
// clang-format off
const std::map<std::string, std::tuple<unsigned int, unsigned int>> SolidMechanicsPhysicsBase::_output_cartesian_params = {
    {"xx", {0, 0}}, {"xy", {0, 1}}, {"xz", {0, 2}},
    {"yx", {1, 0}}, {"yy", {1, 1}}, {"yz", {1, 2}},
    {"zx", {2, 0}}, {"zy", {2, 1}}, {"zz", {2, 2}}
};
// clang-format on

InputParameters
SolidMechanicsPhysicsBase::validParams()
{
  InputParameters params = Action::validParams();

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
  params.addParam<bool>(
      "volumetric_locking_correction", false, "Flag to correct volumetric locking");
  params.addParam<std::vector<MaterialPropertyName>>(
      "eigenstrain_names", {}, "List of eigenstrains to be applied in this strain calculation");
  params.addParam<bool>("use_automatic_differentiation",
                        false,
                        "Flag to use automatic differentiation (AD) objects when possible");

  // Advanced
  params.addParam<std::vector<AuxVariableName>>("save_in", {}, "The displacement residuals");
  params.addParam<std::vector<AuxVariableName>>(
      "diag_save_in", {}, "The displacement diagonal preconditioner terms");
  params.addParamNamesToGroup("save_in diag_save_in", "Advanced");

  // Planar Formulation
  MooseEnum planarFormulationType(
      "Symmetry::NONE WEAK_PLANE_STRESS PLANE_STRAIN GENERALIZED_PLANE_STRAIN", "Symmetry::NONE");
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
                                  SolidMechanicsPhysicsBase::outputProperties(),
                                  "Add scalar quantity output for stress and/or strain");
  params.addParam<MultiMooseEnum>(
      "material_output_order",
      SolidMechanicsPhysicsBase::materialOutputOrders(),
      "Specifies the order of the FE shape function to use for this variable.");
  params.addParam<MultiMooseEnum>(
      "material_output_family",
      SolidMechanicsPhysicsBase::materialOutputFamilies(),
      "Specifies the family of FE shape functions to use for this variable.");
  params.addParamNamesToGroup("generate_output material_output_order material_output_family",
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

SolidMechanicsPhysicsBase::SolidMechanicsPhysicsBase(const InputParameters & parameters)
  : Action(parameters), _use_ad(getParam<bool>("use_automatic_differentiation"))
{
  // Merge SolidMechanicsPhysicsCommon parameters into the parameters of this action
  applyCommonParameters();

  // Append additional output parameters
  mergeMultiMooseEnum("generate_output", "additional_generate_output");
  mergeMultiMooseEnum("material_output_order", "additional_material_output_order");
  mergeMultiMooseEnum("material_output_family", "additional_material_output_family");
}

void
SolidMechanicsPhysicsBase::applyCommonParameters()
{
  auto action = _awh.getActions<SolidMechanicsPhysicsCommon>();
  mooseAssert(!action.empty(), "SolidMechanicsPhysicsCommon not found");

  if (action.size() != 1)
    mooseError("Duplicate SolidMechanicsPhysicsCommon actions not allowed. This could happen if "
               "the input file contains both [QuasiStatic] and [Dynamic] sub-blocks under the "
               "[SolidMechanics] physics block.");

  const_cast<InputParameters *>(&parameters())->applyParameters(action[0]->parameters());
}

void
SolidMechanicsPhysicsBase::mergeMultiMooseEnum(const std::string & dest, const std::string & src)
{
  if (!isParamValid(src))
    return;

  auto dest_enum = getParam<MultiMooseEnum>(dest);
  dest_enum.push_back(getParam<MultiMooseEnum>(src));
  const_cast<InputParameters *>(&parameters())->set<MultiMooseEnum>(dest) = dest_enum;
}

MultiMooseEnum
SolidMechanicsPhysicsBase::outputProperties()
{
  std::string options = "";

  // Rank-0 (scalar-valued) properties
  for (auto prop : query(_output_properties).rank(Rank::ZERO).get())
    options += prop.alias;

  // Higher rank properties will need one additional layer of postprocessing so that we can write
  // them into scalar-valued aux variables.
  // The type of postprocessing is specified by the prefix/suffix.
  for (auto & [key, restriction] : _output_property_restriction)
  {
    auto [prefix, rank, symmetry, type] = restriction;
    for (auto prop : query(_output_properties).rank(rank).symmetry(symmetry).type(type).get())
      options += prefix ? key + '_' + prop.alias : prop.alias + '_' + key;
  }

  return options;
}

MultiMooseEnum
SolidMechanicsPhysicsBase::materialOutputOrders()
{
  return AddAuxVariableAction::getAuxVariableOrders();
}

MultiMooseEnum
SolidMechanicsPhysicsBase::materialOutputFamilies()
{
  return "MONOMIAL LAGRANGE";
}

std::string
SolidMechanicsMaterialProperty::getOutputPropertyName(const std::string & out,
                                                      const std::string & key) const
{
  auto [prefix, rank, symmetry, type] = _output_property_restriction[key];
  auto alias = prefix ? out.substr(key.length() + 1) : out.substr(out.length() - key.length() - 1);
  auto candidates =
      query(_output_properties).alias(alias).rank(q_rank).symmetry(q_symmetry).type(q_type).get();

  // If no candidate is found, treat the alias as the name of a custom tensor property
  if (candidates.empty())
    prop_name = q_alias;
  else
  {
    mooseAssert(candidates.size() == 1,
                "Internal error: Requested output yields multiple candidates.");
    prop_name = candidates[0].name;
  }
}
