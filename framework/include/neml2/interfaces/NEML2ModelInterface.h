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
#include "neml2/csrc/eager/Model.h"
#include "neml2/csrc/eager/load_model.h"
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

  /// Get the NEML2 model (eager runtime)
  const neml2::eager::Model & model() const { return _model; }

  /// Non-const access to the NEML2 model, for runtime parameter writes (set_parameter)
  neml2::eager::Model & model() { return _model; }

  /// Get the target compute device
  const at::Device & device() const { return _device; }

  /// Get the target output device
  const at::Device & output_device() const { return _output_device; }

private:
  /// The device on which to evaluate the NEML2 model
  const at::Device _device;
  /// The device on which to store the outputs
  const at::Device _output_device;
  /// The NEML2 material model. The eager runtime embeds a CPython interpreter and loads the
  /// model directly from the input file -- no ahead-of-time compilation required.
  neml2::eager::Model _model;

#endif // NEML2_ENABLED
};

template <class T>
InputParameters
NEML2ModelInterface<T>::validParams()
{
  InputParameters params = T::validParams();
  params.addParam<DataFileName>("input",
                                "Path to the NEML2 input file containing the NEML2 model(s).");
  params.addParam<std::vector<std::string>>(
      "cli_args",
      {},
      "Additional command line arguments to use when parsing the NEML2 input file.");
  params.addParam<std::vector<std::string>>(
      "load",
      {},
      "External Python extension modules imported (into the embedded interpreter) before the NEML2 "
      "model is built: file paths to .py files or package directories, or importable dotted module "
      "names. Importing them registers any @register_neml2_object types they define -- e.g. NEML2 "
      "models hosted inside a MOOSE app -- so the NEML2 input file can reference them. Mirrors the "
      "neml2 CLI --load flag.");
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
                       : _device),
    // The eager runtime loads the model directly from the .i file (embedding CPython) and pins
    // it to the requested device. The 'cli_args' parameter is currently not forwarded (NEML2 v3's
    // eager load_model takes no extra parse arguments).
    _model(std::string(params.get<DataFileName>("input")),
           params.get<std::string>("model"),
           _device,
           params.get<std::vector<std::string>>("load"))
{
}

template <class T>
void
NEML2ModelInterface<T>::validateModel() const
{
  // NEML2 v3 provides no C++ diagnose() equivalent; constructing the eager model above already
  // parses and validates it. Kept as a hook for derived classes.
}

#endif // NEML2_ENABLED
