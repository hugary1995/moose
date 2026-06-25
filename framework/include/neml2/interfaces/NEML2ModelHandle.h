//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "NEML2Utils.h"

#ifdef NEML2_ENABLED

#include <memory>
#include <string>
#include <map>
#include <vector>
#include <utility>

#include "neml2/csrc/eager/Model.h"
#include "neml2/csrc/eager/load_model.h"
#include "neml2/csrc/dispatchers/factory.h"
#include "neml2/csrc/dispatchers/DispatchedModel.h"
#include "neml2/csrc/dispatchers/SimpleScheduler.h"

/**
 * Runtime-agnostic handle over a NEML2 v3 model.
 *
 * NEML2 v3 ships two C++ runtimes with parallel APIs but no common base class:
 * `neml2::eager::Model` (embeds CPython, loads a model from a source `.i`) and the cpp-aoti
 * `neml2::aoti::DispatchedModel` (loads ahead-of-time-compiled `.pt2` + `_meta.json` artifacts
 * via `neml2::aoti::load_model`, Python-free, dispatched to a device by a scheduler). This handle
 * exposes the single surface the MOOSE integration needs -- introspection (names + base shapes)
 * and evaluation (jacobian / param_jacobian / set_parameter) -- so both the setup-time
 * introspection in NEML2Action and the runtime evaluation in NEML2ModelExecutor use one
 * consistent API regardless of which runtime backs the model.
 *
 * Construct one with ::makeNEML2ModelHandle.
 */
class NEML2ModelHandle
{
public:
  virtual ~NEML2ModelHandle() = default;

  /// @name Introspection, in graph-call order (params keyed by qualified name)
  ///@{
  virtual const std::vector<std::string> & input_names() const = 0;
  virtual const std::vector<std::string> & output_names() const = 0;
  virtual const std::vector<std::vector<int64_t>> & input_base_shapes() const = 0;
  virtual const std::vector<std::vector<int64_t>> & output_base_shapes() const = 0;
  virtual const std::map<std::string, std::vector<int64_t>> & parameter_base_shapes() const = 0;
  ///@}

  /// @name Evaluation
  ///@{
  /// Evaluate the model and its full input Jacobian.
  virtual std::pair<std::map<std::string, at::Tensor>, neml2::aoti::VariablePairJacobian>
  jacobian(const std::map<std::string, at::Tensor> & inputs) const = 0;
  /// Evaluate the model and its parameter Jacobian (d output / d parameter).
  virtual std::pair<std::map<std::string, at::Tensor>, neml2::aoti::VariablePairJacobian>
  param_jacobian(const std::map<std::string, at::Tensor> & inputs) const = 0;
  /// Replace a (runtime-flexible) model parameter's value.
  virtual void set_parameter(const std::string & name, const at::Tensor & value) = 0;
  ///@}
};

/**
 * Handle wrapping the cpp-eager runtime (embeds CPython; loads from a source `.i`).
 */
class EagerModelHandle : public NEML2ModelHandle
{
public:
  EagerModelHandle(const std::string & input_file,
                   const std::string & model_name,
                   const at::Device & device,
                   const std::vector<std::string> & load)
    : _m(input_file, model_name, device, load)
  {
  }

  const std::vector<std::string> & input_names() const override { return _m.input_names(); }
  const std::vector<std::string> & output_names() const override { return _m.output_names(); }
  const std::vector<std::vector<int64_t>> & input_base_shapes() const override
  {
    return _m.input_base_shapes();
  }
  const std::vector<std::vector<int64_t>> & output_base_shapes() const override
  {
    return _m.output_base_shapes();
  }
  const std::map<std::string, std::vector<int64_t>> & parameter_base_shapes() const override
  {
    return _m.parameter_base_shapes();
  }

  std::pair<std::map<std::string, at::Tensor>, neml2::aoti::VariablePairJacobian>
  jacobian(const std::map<std::string, at::Tensor> & inputs) const override
  {
    return _m.jacobian(inputs);
  }
  std::pair<std::map<std::string, at::Tensor>, neml2::aoti::VariablePairJacobian>
  param_jacobian(const std::map<std::string, at::Tensor> & inputs) const override
  {
    return _m.param_jacobian(inputs);
  }
  void set_parameter(const std::string & name, const at::Tensor & value) override
  {
    _m.set_parameter(name, value);
  }

private:
  neml2::eager::Model _m;
};

/**
 * Handle wrapping the cpp-aoti runtime: a `neml2::aoti::DispatchedModel` loaded from a
 * `neml2-compile` stub `.i` via `neml2::aoti::load_model`. A `SimpleScheduler` (whole batch in
 * one chunk) pins the workload to the requested device.
 */
class AOTIModelHandle : public NEML2ModelHandle
{
public:
  AOTIModelHandle(const std::string & stub_file,
                  const std::string & model_name,
                  const at::Device & device)
    : _m(load(stub_file, model_name, device))
  {
  }

  const std::vector<std::string> & input_names() const override { return _m.input_names(); }
  const std::vector<std::string> & output_names() const override { return _m.output_names(); }
  const std::vector<std::vector<int64_t>> & input_base_shapes() const override
  {
    return _m.input_base_shapes();
  }
  const std::vector<std::vector<int64_t>> & output_base_shapes() const override
  {
    return _m.output_base_shapes();
  }
  const std::map<std::string, std::vector<int64_t>> & parameter_base_shapes() const override
  {
    return _m.parameter_base_shapes();
  }

  std::pair<std::map<std::string, at::Tensor>, neml2::aoti::VariablePairJacobian>
  jacobian(const std::map<std::string, at::Tensor> & inputs) const override
  {
    return _m.jacobian(inputs);
  }
  std::pair<std::map<std::string, at::Tensor>, neml2::aoti::VariablePairJacobian>
  param_jacobian(const std::map<std::string, at::Tensor> & inputs) const override
  {
    return _m.param_jacobian(inputs);
  }
  void set_parameter(const std::string & name, const at::Tensor & value) override
  {
    _m.set_parameter(name, value);
  }

private:
  /// Load the dispatched aoti model, pinning the whole batch to `device` in one chunk.
  static neml2::aoti::DispatchedModel
  load(const std::string & stub_file, const std::string & model_name, const at::Device & device)
  {
    neml2::aoti::SimpleScheduler::Config config;
    config.device = device.str();
    config.batch_size = 0; // run the whole batch at once (no chunking)
    return neml2::aoti::load_model(
        stub_file, model_name, std::make_shared<neml2::aoti::SimpleScheduler>(config));
  }

  neml2::aoti::DispatchedModel _m;
};

/**
 * Construct a NEML2 model handle for the requested runtime.
 *
 * @param eager   When true, build a cpp-eager model from the source `.i`; when false, build a
 *                cpp-aoti dispatched model from the compiled-artifact stub `.i`.
 * @param input   Path to the NEML2 `.i`: the source file for cpp-eager, the `neml2-compile` stub
 *                for cpp-aoti.
 * @param model   Name of the model in the `.i`.
 * @param device  Compute device. cpp-eager pins the model to it; cpp-aoti dispatches to it via a
 *                SimpleScheduler (the device must name an artifact subfolder).
 * @param load    External Python extension paths to import before building (cpp-eager only).
 */
inline std::unique_ptr<NEML2ModelHandle>
makeNEML2ModelHandle(bool eager,
                     const std::string & input,
                     const std::string & model,
                     const at::Device & device,
                     const std::vector<std::string> & load)
{
  if (eager)
    return std::make_unique<EagerModelHandle>(input, model, device, load);
  return std::make_unique<AOTIModelHandle>(input, model, device);
}

#endif // NEML2_ENABLED
