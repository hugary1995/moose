//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include <thread>
#include <utility>
#include <tuple>
#include "NEML2Utils.h"
#include "InputParameters.h"

#ifdef NEML2_ENABLED
#include "NEML2ModelHandle.h"
#endif

/**
 * Interface class to provide common input parameters, members, and methods for MOOSEObjects that
 * use NEML2 models.
 */
template <class T>
class NEML2ModelInterface : public T
{
public:
  static InputParameters validParams();

  template <typename... P>
  NEML2ModelInterface(const InputParameters & params, P &&... args);

#ifdef NEML2_ENABLED

protected:
  /**
   * Validate the NEML2 material model. Note that the developer is responsible for calling this
   * method at the appropriate times, for example, at initialSetup().
   */
  virtual void validateModel() const;

  /// Get the NEML2 model handle (cpp-aoti by default, cpp-eager when eager=true)
  const NEML2ModelHandle & model() const { return *_model; }

  /// Non-const access to the NEML2 model handle, for runtime parameter writes (set_parameter)
  NEML2ModelHandle & model() { return *_model; }

  /// Get the target compute device
  const at::Device & device() const { return _device; }

  /// Get the target output device
  const at::Device & output_device() const { return _output_device; }

private:
  /// Resolve the 'load' DataFileName list (whose relative paths InputParameters has already
  /// resolved against the input file + data search path) to plain path strings for the eager
  /// runtime's external-extension loader.
  static std::vector<std::string> loadExtensionPaths(const InputParameters & params)
  {
    const auto files = params.get<std::vector<DataFileName>>("load");
    return std::vector<std::string>(files.begin(), files.end());
  }

  /// The device on which to evaluate the NEML2 model
  const at::Device _device;
  /// The device on which to store the outputs
  const at::Device _output_device;
  /// The NEML2 material model handle. Defaults to the cpp-aoti runtime (loads an ahead-of-time-
  /// compiled artifact specified by 'meta'); the cpp-eager runtime (embeds a CPython interpreter
  /// and loads the model from the source 'input' file) is opt-in via the 'eager' parameter.
  std::unique_ptr<NEML2ModelHandle> _model;

#endif // NEML2_ENABLED
};

template <class T>
InputParameters
NEML2ModelInterface<T>::validParams()
{
  InputParameters params = T::validParams();
  params.addParam<bool>(
      "eager",
      false,
      "Use the cpp-eager NEML2 runtime, which embeds a Python interpreter and loads the model "
      "from the source 'input' file at runtime. Default false -- use the cpp-aoti runtime, which "
      "loads an ahead-of-time-compiled artifact (no Python at runtime) from the compile stub "
      "given by 'input'. The 'load' parameter is only valid with eager=true.");
  params.addParam<DataFileName>(
      "input",
      "Path to the NEML2 '.i' file: the source model file for the cpp-eager runtime, or the "
      "'<model>_aoti.i' stub produced by neml2-compile for the cpp-aoti runtime.");
  params.addParam<std::vector<std::string>>(
      "cli_args",
      {},
      "Additional command line arguments to use when parsing the NEML2 input file.");
  params.addParam<std::vector<DataFileName>>(
      "load",
      {},
      "External Python extension files (paths to .py files or package directories, resolved "
      "relative to the input file and the application's data search path) imported into the "
      "embedded interpreter before the NEML2 model is built. Importing them registers any "
      "@register_neml2_object types they define -- e.g. NEML2 models hosted inside a MOOSE app -- "
      "so the NEML2 input file can reference them. Mirrors the neml2 CLI --load flag. Only valid "
      "with the cpp-eager runtime (eager=true).");
  params.addParam<std::string>(
      "model",
      "",
      "Name of the NEML2 model, i.e., the string inside the brackets [] in the NEML2 input file "
      "that corresponds to the model you want to use.");
  params.addParam<std::string>(
      "device",
      "Device on which to evaluate the NEML2 model. The string supplied must follow the following "
      "schema: (cpu|cuda)[:<device-index>] where cpu or cuda specifies the device type, and "
      ":<device-index> optionally specifies a device index. For example, device='cpu' sets the "
      "target compute device to be CPU, and device='cuda:1' sets the target compute device to be "
      "CUDA with device ID 1. If not specified, default to the compute device specified via the "
      "command line argument --compute-device.");
  params.addParam<std::string>(
      "output_device",
      "Similar to the 'device' parameter, this parameter specifies the device on which to store "
      "the outputs. Default to be the same as 'device'.");

  return params;
}

#ifndef NEML2_ENABLED

template <class T>
template <typename... P>
NEML2ModelInterface<T>::NEML2ModelInterface(const InputParameters & params, P &&... args)
  : T(params, args...)
{
}

#else

template <class T>
template <typename... P>
NEML2ModelInterface<T>::NEML2ModelInterface(const InputParameters & params, P &&... args)
  : T(params, args...),
    _device(params.isParamValid("device") ? at::Device(params.get<std::string>("device"))
                                          : this->getMooseApp().getLibtorchDevice()),
    _output_device(params.isParamValid("output_device")
                       ? at::Device(params.get<std::string>("output_device"))
                       : _device)
{
  // The cpp-eager runtime loads the model from the source .i (embedding CPython); the cpp-aoti
  // runtime (default) loads the ahead-of-time-compiled artifact named by 'meta'. The 'cli_args'
  // parameter is currently not forwarded (NEML2 v3's load takes no extra parse arguments).
  const bool eager = params.get<bool>("eager");
  const auto load = loadExtensionPaths(params);

  // 'load' imports Python extensions into the embedded interpreter -- meaningful only for cpp-eager.
  if (!eager && !load.empty())
    this->paramError("load",
                     "The 'load' parameter imports Python extensions into the embedded "
                     "interpreter and is only valid with the cpp-eager runtime. Set eager=true to "
                     "use it, or remove it for the cpp-aoti runtime.");

  // 'input' is the source .i (cpp-eager) or the compiled-artifact stub .i (cpp-aoti).
  if (!params.isParamValid("input"))
    this->paramError("input", "'input' (the NEML2 source or AOTI stub '.i') is required.");

  _model = makeNEML2ModelHandle(eager,
                                std::string(params.get<DataFileName>("input")),
                                params.get<std::string>("model"),
                                _device,
                                load);
}

template <class T>
void
NEML2ModelInterface<T>::validateModel() const
{
  // NEML2 v3 provides no C++ diagnose() equivalent; constructing the eager model above already
  // parses and validates it. Kept as a hook for derived classes.
}

#endif // NEML2_ENABLED
