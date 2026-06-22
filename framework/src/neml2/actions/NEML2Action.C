//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NEML2Action.h"
#include "NEML2ActionCommon.h"
#include "FEProblem.h"
#include "Factory.h"
#include "NEML2Utils.h"
#include "InputParameterWarehouse.h"

#ifdef NEML2_ENABLED
#include <algorithm>

// parseLag / lagName / contains are shared with NEML2ModelExecutor and live in NEML2Utils so the
// NEML2 unity build does not see duplicate definitions.
using NEML2Utils::contains;
using NEML2Utils::lagName;
using NEML2Utils::parseLag;

namespace
{
// Map a NEML2 variable base shape to the corresponding MOOSE C++ tensor type string. Returns an
// empty string if the shape is not (yet) supported. A derivative block's type is obtained by
// passing the concatenation of the output and input base shapes.
std::string
shapeToMooseType(const std::vector<int64_t> & s)
{
  if (s.empty())
    return "Real";
  if (s == std::vector<int64_t>{3})
    return "RealVectorValue";
  if (s == std::vector<int64_t>{6})
    return "SymmetricRankTwoTensor";
  if (s == std::vector<int64_t>{3, 3})
    return "RankTwoTensor";
  if (s == std::vector<int64_t>{6, 6})
    return "SymmetricRankFourTensor";
  if (s == std::vector<int64_t>{3, 3, 3, 3})
    return "RankFourTensor";
  return "";
}

// Build a name -> base-shape lookup for the model's inputs / outputs.
std::map<std::string, std::vector<int64_t>>
shapeMap(const std::vector<std::string> & names, const std::vector<std::vector<int64_t>> & shapes)
{
  std::map<std::string, std::vector<int64_t>> m;
  for (std::size_t i = 0; i < names.size(); ++i)
    m[names[i]] = shapes[i];
  return m;
}
}
#endif

registerMooseAction("MooseApp", NEML2Action, "parse_neml2");
registerMooseAction("MooseApp", NEML2Action, "add_material");
registerMooseAction("MooseApp", NEML2Action, "add_user_object");

InputParameters
NEML2Action::validParams()
{
  InputParameters params = NEML2ActionCommon::commonParams();
  params.addClassDescription("Set up the NEML2 material model");
  params.addParam<std::string>("executor_name",
                               "Name of the NEML2ModelExecutor user object. The default name is "
                               "'neml2_<model-name>_<block-name>' where <model-name> is the NEML2 "
                               "model's name, and <block-name> is this action sub-block's name.");
  params.addParam<std::string>(
      "batch_index_generator_name",
      "Name of the NEML2BatchIndexGenerator user object. The default name is "
      "'neml2_index_<model-name>_<block-name>' where <model-name> is the NEML2 model's name, and "
      "<block-name> is this action sub-block's name.");
  params.addParam<std::vector<SubdomainName>>(
      "block", {}, "List of blocks (subdomains) where the material model is defined");
  return params;
}

NEML2Action::NEML2Action(const InputParameters & params)
  : Action(params),
    _executor_name(isParamValid("executor_name")
                       ? getParam<std::string>("executor_name")
                       : "neml2_" + getParam<std::string>("model") + "_" + name()),
    _idx_generator_name(isParamValid("batch_index_generator_name")
                            ? getParam<std::string>("batch_index_generator_name")
                            : "neml2_index_" + getParam<std::string>("model") + "_" + name()),
    _block(getParam<std::vector<SubdomainName>>("block")),
    _skip_input_variables(getParam<std::vector<std::string>>("skip_input_variables"))
{
  NEML2Utils::assertNEML2Enabled();

  // Apply parameters under the common area, i.e., under [NEML2]
  const auto & all_params = _app.getInputParameterWarehouse().getInputParameters();
  auto & sub_block_params = *(all_params.find(uniqueActionName())->second.get());
  const auto & common_action = getCommonAction();
  sub_block_params.applyParameters(common_action.parameters());

  // Set up optional output variable initialization
  auto init_vars = getParam<std::vector<MaterialPropertyName>>("initialize_outputs");
  auto init_vals = getParam<std::vector<MaterialPropertyName>>("initialize_output_values");
  if (init_vars.size() != init_vals.size())
    paramError("initialize_outputs",
               "initialize_outputs should have the same length as initialize_output_values");
  for (auto i : index_range(init_vars))
    _initialize_output_values[init_vars[i]] = init_vals[i];

  // Set up additional outputs for each requested material property
  auto outputs = getParam<std::vector<MaterialPropertyName>>("export_outputs");
  auto output_targets = getParam<std::vector<std::vector<OutputName>>>("export_output_targets");
  if (outputs.size() != output_targets.size())
    paramError("export_outputs",
               "export_outputs should have the same length as export_output_targets");
  for (auto i : index_range(outputs))
    _export_output_targets[outputs[i]] = output_targets[i];

#ifdef NEML2_ENABLED
  // File name and CLI args
  _fname = getParam<DataFileName>("input");
  _cli_args = getParam<std::vector<std::string>>("cli_args");

  // Load the model via the eager runtime, used here only for introspection during setup. The
  // 'cli_args' parameter is currently not forwarded (NEML2 v3's eager load_model takes none).
  _model = std::make_unique<neml2::eager::Model>(std::string(_fname), getParam<std::string>("model"));
#endif
}

const NEML2ActionCommon &
NEML2Action::getCommonAction() const
{
  auto common_block = _awh.getActions<NEML2ActionCommon>();
  mooseAssert(common_block.size() == 1, "There must exist one and only one common NEML2 action.");
  return *common_block[0];
}

#ifndef NEML2_ENABLED

void
NEML2Action::act()
{
}

#else

static std::string
obscureObjectName(const std::string & name,
                  const std::string & prefix,
                  const std::string & suffix,
                  const std::string & block)
{
  return "__" + prefix + "_" + name + "_" + suffix + "_" + block + "__";
}

void
NEML2Action::act()
{
  if (_current_task == "parse_neml2")
  {
    if (_app.parameters().have_parameter<bool>("parse_neml2_only"))
      if (!_app.parameters().get<bool>("parse_neml2_only"))
        return;
    printSummary();
  }

  if (_current_task == "add_user_object")
  {
    setupInputMappings(*_model);
    setupOutputMappings(*_model);
    setupDerivativeMappings(*_model);

    printSummary();

    // Create and register a MOOSEToNEML2 gatherer user object, returning its name
    auto addGatherer = [&](const std::string & moose_name,
                           const std::string & neml2_name,
                           NEML2Utils::MOOSEIOType moose_type,
                           const std::string & moose_tensor_type,
                           const std::string & suffix,
                           const std::string & type_prefix = "")
    {
      auto obj_name = obscureObjectName(moose_name, "moose_to_neml2", suffix, name());
      auto obj_type = "MOOSE" + type_prefix + moose_tensor_type + "ToNEML2";
      auto obj_params = _factory.getValidParams(obj_type);
      obj_params.set<std::string>("from_moose") = moose_name;
      obj_params.set<std::string>("to_neml2") = neml2_name;
      obj_params.set<MooseEnum>("quantity_type").assign(static_cast<int>(moose_type));
      obj_params.set<std::vector<SubdomainName>>("block") = _block;
      _problem->addUserObject(obj_type, obj_name, obj_params);
      return obj_name;
    };

    // MOOSEToNEML2 input gatherers. The NEML2 target name carries the lag suffix (var~N); the
    // MOOSE source is the un-lagged base name, and old (lag 1) values are read with the "Old"
    // gatherer variant.
    std::vector<UserObjectName> gatherers;
    for (const auto & input : _inputs)
      gatherers.push_back(addGatherer(input.name,
                                      lagName(input.name, input.history_order),
                                      input.moose_type,
                                      input.moose_tensor_type,
                                      std::to_string(input.history_order),
                                      input.history_order == 1 ? "Old" : ""));

    // Additional NEML2Kernels that provide input data
    for (const auto & kernel_name : getParam<std::vector<std::string>>("input_kernels"))
    {
      if (!contains(_model->input_names(), kernel_name))
        paramError("input_kernels",
                   "The NEML2 kernel ",
                   kernel_name,
                   " name does not match any NEML2 input variable.");
      gatherers.push_back(kernel_name);
    }

    // The index generator UO
    {
      auto type = "NEML2BatchIndexGenerator";
      auto params = _factory.getValidParams(type);
      params.applyParameters(parameters());
      _problem->addUserObject(type, _idx_generator_name, params);
    }

    // The Executor UO
    {
      auto type = "NEML2ModelExecutor";
      auto params = _factory.getValidParams(type);
      params.applyParameters(parameters());
      params.set<UserObjectName>("batch_index_generator") = _idx_generator_name;
      params.set<std::vector<UserObjectName>>("gatherers") = gatherers;
      _problem->addUserObject(type, _executor_name, params);
    }
  }

  if (_current_task == "add_material")
  {
    // Create and register a NEML2ToMOOSE material property retriever; `extra` lets each caller
    // add the bits that are unique to outputs vs. derivatives.
    auto addRetriever = [&](const std::string & moose_name,
                            const std::string & neml2_var,
                            const std::string & moose_tensor_type,
                            auto && extra)
    {
      auto obj_name = obscureObjectName(moose_name, "neml2_to_moose", "", name());
      auto obj_type = "NEML2ToMOOSE" + moose_tensor_type + "MaterialProperty";
      auto obj_params = _factory.getValidParams(obj_type);
      obj_params.set<UserObjectName>("neml2_executor") = _executor_name;
      obj_params.set<MaterialPropertyName>("to_moose") = moose_name;
      obj_params.set<std::string>("from_neml2") = neml2_var;
      obj_params.set<std::vector<SubdomainName>>("block") = _block;
      if (_export_output_targets.count(moose_name))
        obj_params.set<std::vector<OutputName>>("outputs") = _export_output_targets[moose_name];
      extra(obj_params);
      _problem->addMaterial(obj_type, obj_name, obj_params);
    };

    // NEML2ToMOOSE output retrievers
    for (const auto & output : _outputs)
    {
      if (output.moose_type != NEML2Utils::MOOSEIOType::MATERIAL)
        paramError("moose_output_types",
                   "Unsupported type corresponding to the moose output ",
                   output.name);
      addRetriever(output.name,
                   output.name,
                   output.moose_tensor_type,
                   [&](InputParameters & p)
                   {
                     if (_initialize_output_values.count(output.name))
                       p.set<MaterialPropertyName>("moose_material_property_init") =
                           _initialize_output_values[output.name];
                   });
    }

    // NEML2ToMOOSE derivative retrievers
    for (const auto & deriv : _derivs)
      addRetriever(deriv.name,
                   deriv.y,
                   deriv.moose_tensor_type,
                   [&](InputParameters & p)
                   { p.set<std::string>("neml2_input_derivative") = deriv.x; });
  }
}

NEML2Utils::MOOSEIOType
NEML2Action::inferMOOSEIOType(const std::string & name, bool is_scalar) const
{
  // only scalar-typed variables can come from a scalar variable, function, or field variable
  if (is_scalar)
  {
    bool is_time = _problem->isTransient() && (name == "t" || name == "time");
    bool has_scalar = _problem->hasScalarVariable(name);
    bool has_func = _problem->hasFunction(name);
    bool has_var = _problem->hasVariable(name);
    if (int(is_time) + int(has_scalar) + int(has_func) + int(has_var) > 1)
      mooseError("Trying to infer MOOSE data type for NEML2 variable ",
                 name,
                 ". The name matches multiple types (",
                 (is_time ? "time " : ""),
                 (has_scalar ? "scalar variable " : ""),
                 (has_func ? "function " : ""),
                 (has_var ? "variable " : ""),
                 "). To avoid ambiguity, please explicitly specify the type in the NEML2 action.");
    if (is_time)
      return NEML2Utils::MOOSEIOType::TIME;
    if (has_scalar)
      return NEML2Utils::MOOSEIOType::SCALAR;
    if (has_func)
      return NEML2Utils::MOOSEIOType::FUNCTION;
    if (has_var)
      return NEML2Utils::MOOSEIOType::VARIABLE;
    // if neither function nor variable exists, let's assume it's a material property
    // note that we can't explicitly check if a material property with the given name exists,
    // because materials are added _after_ user objects (see Moose.C)
  }

  // non-scalar can only come from material properties
  return NEML2Utils::MOOSEIOType::MATERIAL;
}

void
NEML2Action::setupInputMappings(const neml2::eager::Model & model)
{
  const auto & kernels = getParam<std::vector<std::string>>("input_kernels");
  const auto in_shapes = shapeMap(model.input_names(), model.input_base_shapes());

  // Default mapping
  for (const auto & vname : model.input_names())
  {
    const auto [base, order] = parseLag(vname);

    // user requested to skip (by base name)
    if (contains(_skip_input_variables, base))
      continue;

    // skip if the input is directly provided by a custom kernel (matched by full NEML2 name)
    if (contains(kernels, vname))
      continue;

    // skip stateful (old) inputs when managing state advance on the compute device
    if (getParam<bool>("manage_state_advance") && order > 0)
      continue;

    const auto & base_shape = in_shapes.at(vname);
    const auto mtype = shapeToMooseType(base_shape);
    if (mtype.empty())
      mooseError("NEML2 input variable ", vname, " has a base shape not yet mapped to a MOOSE type.");

    _inputs.push_back({base,
                       inferMOOSEIOType(base, base_shape.empty()),
                       mtype,
                       static_cast<std::size_t>(order)});
  }

  // User-specified mapping (overrides default mapping)
  const auto [input_types, inputs] =
      getInputParameterMapping<NEML2Utils::MOOSEIOType, std::string>("input_types", "inputs");

  for (auto i : index_range(inputs))
  {
    const auto base_i = parseLag(inputs[i]).first;
    // Check if the input variable also appears in skip_input_variables
    if (contains(_skip_input_variables, base_i) || contains(_skip_input_variables, inputs[i]))
      paramError("skip_input_variables",
                 "The input variable ",
                 inputs[i],
                 " is listed in skip_input_variables, but it also appears in inputs. Please "
                 "remove it from either list.");
    // Check if the input variable exists in the NEML2 model
    if (!contains(model.input_names(), inputs[i]))
      paramError("inputs", "The neml2 input variable ", inputs[i], " does not exist.");
    // Check if the input variable is already gathered by a custom kernel
    if (contains(kernels, inputs[i]))
      paramError("inputs",
                 "The input variable ",
                 inputs[i],
                 " is listed in inputs, but it also appears in input_kernels. Please "
                 "remove it from either list.");
    // Check if the input variable is stateful and manage_state_advance is true
    if (getParam<bool>("manage_state_advance") && parseLag(inputs[i]).second > 0)
      paramError("inputs",
                 "The input variable ",
                 inputs[i],
                 " is listed in inputs, but it is an old (stateful) variable and "
                 "manage_state_advance is true. Please remove it from inputs, or set "
                 "manage_state_advance to false.");

    // Get the existing mapping for this neml2 input variable and override its MOOSE type
    for (auto & input : _inputs)
      if (input.name == base_i)
      {
        input.moose_type = input_types[i];
        break;
      }
  }
}

void
NEML2Action::setupOutputMappings(const neml2::eager::Model & model)
{
  if (!getParam<bool>("auto_output"))
    return;

  const auto out_shapes = shapeMap(model.output_names(), model.output_base_shapes());

  for (const auto & oname : model.output_names())
  {
    const auto mtype = shapeToMooseType(out_shapes.at(oname));
    if (mtype.empty())
      mooseError("NEML2 output variable ", oname, " has a base shape not yet mapped to a MOOSE type.");
    _outputs.push_back({oname, NEML2Utils::MOOSEIOType::MATERIAL, mtype, /*history_order=*/0});
  }
}

void
NEML2Action::setupDerivativeMappings(const neml2::eager::Model & model)
{
  const auto derivs = getParam<std::vector<std::vector<std::string>>>("derivatives");
  const auto in_shapes = shapeMap(model.input_names(), model.input_base_shapes());
  const auto out_shapes = shapeMap(model.output_names(), model.output_base_shapes());

  for (auto i : index_range(derivs))
  {
    if (derivs[i].size() != 2)
      paramError("derivatives", "The length of each pair in derivatives must be 2.");

    const auto & y = derivs[i][0];
    const auto & x = derivs[i][1];
    if (!contains(model.output_names(), y))
      paramError("derivatives", "The NEML2 output variable ", y, " does not exist.");
    if (!contains(model.input_names(), x))
      paramError("derivatives", "The NEML2 input variable ", x, " does not exist.");

    // The derivative block's base shape is the concatenation of the output and input base shapes.
    std::vector<int64_t> dshape = out_shapes.at(y);
    const auto & xs = in_shapes.at(x);
    dshape.insert(dshape.end(), xs.begin(), xs.end());
    const auto mtype = shapeToMooseType(dshape);
    if (mtype.empty())
      mooseError("The NEML2 derivative of ",
                 y,
                 " with respect to ",
                 x,
                 " has a base shape not yet mapped to a MOOSE type.");

    const auto deriv_name = derivativePropertyNameFirst(y, x);
    _derivs.push_back({deriv_name, y, x, mtype});
  }
}

void
NEML2Action::printSummary() const
{
  if (!_app.parameters().have_parameter<bool>("parse_neml2_only"))
    return;

  // Save formatting of the output stream so that we can restore it later
  const auto flags = _console.flags();

  // Default width for the summary
  const int width = 79;

  _console << std::endl;
  _console << COLOR_CYAN << std::setw(width) << std::setfill('*') << std::left
           << "NEML2 MATERIAL MODEL SUMMARY BEGIN " << std::setfill(' ') << COLOR_DEFAULT
           << std::endl;

  // Metadata
  _console << "NEML2 input file location: " << fname() << std::endl;
  _console << "NEML2 action path:         " << parameters().blockFullpath() << std::endl;

  // List inputs and outputs of the model
  _console << COLOR_CYAN << std::setw(width) << std::setfill('-') << std::left
           << "Material model structure " << std::setfill(' ') << COLOR_DEFAULT << std::endl;
  _console << "NEML2 model inputs:" << std::endl;
  for (const auto & iname : _model->input_names())
    _console << "  - " << iname << std::endl;
  _console << "NEML2 model outputs:" << std::endl;
  for (const auto & oname : _model->output_names())
    _console << "  - " << oname << std::endl;

  // List transfer between MOOSE and NEML2
  if (!_app.parameters().get<bool>("parse_neml2_only"))
  {
    _console << COLOR_CYAN << std::setw(width) << std::setfill('-') << std::left
             << "Transfer between MOOSE and NEML2 " << std::setfill(' ') << COLOR_DEFAULT
             << std::endl;

    _console << "MOOSE --> NEML2" << std::endl;

    // List input transfer, MOOSE -> NEML2
    for (const auto & input : _inputs)
      _console << "  - " << (input.history_order > 0 ? ("(old) " + input.name) : input.name) << " ("
               << NEML2Utils::stringify(input.moose_type) << ")" << std::endl;

    _console << "MOOSE <-- NEML2" << std::endl;

    // List output transfer, NEML2 -> MOOSE
    for (const auto & output : _outputs)
      _console << "  - " << output.name << std::endl;

    // List derivative transfer, NEML2 -> MOOSE
    for (const auto & deriv : _derivs)
      _console << "  - " << deriv.name << std::endl;
  }

  _console << COLOR_CYAN << std::setw(width) << std::setfill('*') << std::left
           << "NEML2 MATERIAL MODEL SUMMARY END " << std::setfill(' ') << COLOR_DEFAULT << std::endl
           << std::endl;

  // Restore the formatting of the output stream
  _console.flags(flags);
}

std::size_t
NEML2Action::getLongestMOOSEName() const
{
  std::size_t max_moose_name_length = 0;
  for (const auto & input : _inputs)
  {
    auto n = input.name.size();
    if (input.history_order > 0)
      n += 6; // 6 is the length of "(old) "
    max_moose_name_length = std::max(max_moose_name_length, n);
  }
  for (const auto & output : _outputs)
    max_moose_name_length = std::max(max_moose_name_length, output.name.size());
  for (const auto & deriv : _derivs)
    max_moose_name_length = std::max(max_moose_name_length, deriv.name.size());
  return max_moose_name_length;
}
#endif // NEML2_ENABLED
