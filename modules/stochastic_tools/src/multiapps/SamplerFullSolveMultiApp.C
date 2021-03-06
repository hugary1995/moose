//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// StochasticTools includes
#include "SamplerFullSolveMultiApp.h"
#include "Sampler.h"
#include "StochasticToolsTransfer.h"

registerMooseObject("StochasticToolsApp", SamplerFullSolveMultiApp);

InputParameters
SamplerFullSolveMultiApp::validParams()
{
  InputParameters params = FullSolveMultiApp::validParams();
  params += SamplerInterface::validParams();
  params.addClassDescription(
      "Creates a full-solve type sub-application for each row of each Sampler matrix.");
  params.addParam<SamplerName>("sampler", "The Sampler object to utilize for creating MultiApps.");
  params.suppressParameter<std::vector<Point>>("positions");
  params.suppressParameter<bool>("output_in_position");
  params.suppressParameter<std::vector<FileName>>("positions_file");
  params.suppressParameter<Real>("move_time");
  params.suppressParameter<std::vector<Point>>("move_positions");
  params.suppressParameter<std::vector<unsigned int>>("move_apps");
  params.set<bool>("use_positions") = false;

  MooseEnum modes("normal=0 batch-reset=1 batch-restore=2", "normal");
  params.addParam<MooseEnum>(
      "mode",
      modes,
      "The operation mode, 'normal' creates one sub-application for each row in the Sampler and "
      "'batch' creates one sub-application for each processor and re-executes for each row.");
  return params;
}

SamplerFullSolveMultiApp::SamplerFullSolveMultiApp(const InputParameters & parameters)
  : FullSolveMultiApp(parameters),
    SamplerInterface(this),
    _sampler(SamplerInterface::getSampler("sampler")),
    _mode(getParam<MooseEnum>("mode").getEnum<StochasticTools::MultiAppMode>()),
    _local_batch_app_index(0),
    _solved_once(false),
    _perf_solve_step(registerTimedSection("solveStep", 1)),
    _perf_solve_batch_step(registerTimedSection("solveStepBatch", 1)),
    _perf_command_line_args(registerTimedSection("getCommandLineArgsParamHelper", 4))
{
  if (_mode == StochasticTools::MultiAppMode::BATCH_RESET ||
      _mode == StochasticTools::MultiAppMode::BATCH_RESTORE)
  {
    if (n_processors() > _sampler.getNumberOfRows())
      paramError(
          "mode",
          "There appears to be more available processors (",
          n_processors(),
          ") than samples (",
          _sampler.getNumberOfRows(),
          "), this is not supported in "
          "batch mode. Consider switching to \'normal\' to allow multiple processors per sample.");
    init(n_processors());
  }
  else
    init(_sampler.getNumberOfRows());

  _number_of_sampler_rows = _sampler.getNumberOfRows();
}

void SamplerFullSolveMultiApp::preTransfer(Real /*dt*/, Real /*target_time*/)
{
  if (_mode == StochasticTools::MultiAppMode::NORMAL)
  {
    // Reinitialize MultiApp size
    const auto num_rows = _sampler.getNumberOfRows();
    if (num_rows != _number_of_sampler_rows)
    {
      init(_sampler.getNumberOfRows());
      _number_of_sampler_rows = num_rows;
    }

    // Reinitialize app to original state prior to solve, if a solve has occured
    if (_solved_once)
      initialSetup();
  }
}

bool
SamplerFullSolveMultiApp::solveStep(Real dt, Real target_time, bool auto_advance)
{
  TIME_SECTION(_perf_solve_step);

  mooseAssert(_my_num_apps, _sampler.getNumberOfLocalRows());

  bool last_solve_converged = true;

  if (_mode == StochasticTools::MultiAppMode::BATCH_RESET ||
      _mode == StochasticTools::MultiAppMode::BATCH_RESTORE)
    last_solve_converged = solveStepBatch(dt, target_time, auto_advance);
  else
    last_solve_converged = FullSolveMultiApp::solveStep(dt, target_time, auto_advance);

  _solved_once = true;

  return last_solve_converged;
}

bool
SamplerFullSolveMultiApp::solveStepBatch(Real dt, Real target_time, bool auto_advance)
{
  TIME_SECTION(_perf_solve_batch_step);

  // Value to return
  bool last_solve_converged = true;

  // Reinitialize app to original state prior to solve, if a solve has occured
  if (_solved_once)
    initialSetup();

  // List of active relevant Transfer objects
  std::vector<std::shared_ptr<StochasticToolsTransfer>> to_transfers =
      getActiveStochasticToolsTransfers(MultiAppTransfer::TO_MULTIAPP);
  std::vector<std::shared_ptr<StochasticToolsTransfer>> from_transfers =
      getActiveStochasticToolsTransfers(MultiAppTransfer::FROM_MULTIAPP);

  // Initialize to/from transfers
  for (auto transfer : to_transfers)
    transfer->initializeToMultiapp();

  for (auto transfer : from_transfers)
    transfer->initializeFromMultiapp();

  if (_mode == StochasticTools::MultiAppMode::BATCH_RESTORE)
    backup();

  // Perform batch MultiApp solves
  _local_batch_app_index = 0;
  for (dof_id_type i = _sampler.getLocalRowBegin(); i < _sampler.getLocalRowEnd(); ++i)
  {
    for (auto & transfer : to_transfers)
    {
      transfer->setGlobalMultiAppIndex(i);
      transfer->executeToMultiapp();
    }

    last_solve_converged = FullSolveMultiApp::solveStep(dt, target_time, auto_advance);

    for (auto & transfer : from_transfers)
    {
      transfer->setGlobalMultiAppIndex(i);
      transfer->executeFromMultiapp();
    }

    if (i < _sampler.getLocalRowEnd() - 1)
    {
      if (_mode == StochasticTools::MultiAppMode::BATCH_RESTORE)
        restore();
      else
      {
        // The app is being reset for the next loop, thus the batch index must be indexed as such
        _local_batch_app_index = i + 1;
        resetApp(_local_batch_app_index, target_time);
        initialSetup();
      }
    }
  }

  // Finalize to/from transfers
  for (auto transfer : to_transfers)
    transfer->finalizeToMultiapp();
  for (auto transfer : from_transfers)
    transfer->finalizeFromMultiapp();

  return last_solve_converged;
}

std::vector<std::shared_ptr<StochasticToolsTransfer>>
SamplerFullSolveMultiApp::getActiveStochasticToolsTransfers(Transfer::DIRECTION direction)
{
  std::vector<std::shared_ptr<StochasticToolsTransfer>> output;
  const ExecuteMooseObjectWarehouse<Transfer> & warehouse =
      _fe_problem.getMultiAppTransferWarehouse(direction);
  for (std::shared_ptr<Transfer> transfer : warehouse.getActiveObjects())
  {
    auto ptr = std::dynamic_pointer_cast<StochasticToolsTransfer>(transfer);
    if (ptr && ptr->getMultiApp().get() == this)
      output.push_back(ptr);
  }
  return output;
}

std::string
SamplerFullSolveMultiApp::getCommandLineArgsParamHelper(unsigned int local_app)
{
  TIME_SECTION(_perf_command_line_args);

  // Since we only store param_names in cli_args, we need to find the values for each param from
  // sampler data and combine them to get full command line option strings.
  std::vector<Real> row = _sampler.getNextLocalRow();

  std::ostringstream oss;
  const std::vector<std::string> & cli_args_name =
      MooseUtils::split(FullSolveMultiApp::getCommandLineArgsParamHelper(local_app), ";");

  for (dof_id_type col = 0; col < _sampler.getNumberOfCols(); ++col)
  {
    if (col > 0)
      oss << ";";
    oss << cli_args_name[col] << "=" << Moose::stringify(row[col]);
  }
  return oss.str();
}
