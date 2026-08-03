//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TransientExecutor.h"
#include "FEProblemBase.h"
#include "NonlinearSystemBase.h"
#include "AuxiliarySystem.h"

#include <limits>
#include <algorithm>

registerMooseObject("MooseApp", TransientExecutor);

InputParameters
TransientExecutor::validParams()
{
  InputParameters params = Executor::validParams();
  params.addClassDescription(
      "A transient executor: drives inner executors (e.g. the outer field-split SPIN "
      "NewtonSNESExecutor) over pseudo-time steps, advancing state between converged steps so "
      "old solutions (d_old for irreversibility) are populated. Fixed target dt with cutback on "
      "non-convergence; no TimeStepper object (suitable for quasi-static problems).");
  params.addRequiredParam<std::vector<ExecutorName>>(
      "inner_executors", "The executor(s) to run each time step (e.g. the outer SPIN executor).");
  params.addParam<Real>("start_time", 0.0, "Simulation start time.");
  params.addParam<Real>("end_time", 1.0e30, "Simulation end time.");
  params.addParam<Real>("dt", 1.0, "Target (nominal) time step.");
  params.addParam<unsigned int>(
      "num_steps", std::numeric_limits<unsigned int>::max(), "Maximum number of time steps.");
  params.addParam<Real>("dtmin", 1.0e-12, "Minimum dt; the run fails if cutback goes below this.");
  params.addParam<Real>("dtmax", 1.0e30, "Maximum dt.");
  params.addParam<Real>(
      "growth_factor", 2.0, "Factor to grow dt by after a converged step (capped at dt).");
  params.addParam<Real>(
      "cutback_factor", 0.5, "Factor to shrink dt by after a non-converged step.");
  params.addParam<Real>(
      "timestep_tolerance", 1.0e-12, "Tolerance when comparing the current time against end_time.");
  return params;
}

TransientExecutor::TransientExecutor(const InputParameters & params)
  : Executor(params),
    _fe_problem(*params.getCheckedPointerParam<FEProblemBase *>("_fe_problem_base")),
    _start_time(getParam<Real>("start_time")),
    _end_time(getParam<Real>("end_time")),
    _dt_nominal(getParam<Real>("dt")),
    _num_steps(getParam<unsigned int>("num_steps")),
    _dtmin(getParam<Real>("dtmin")),
    _dtmax(getParam<Real>("dtmax")),
    _growth_factor(getParam<Real>("growth_factor")),
    _cutback_factor(getParam<Real>("cutback_factor")),
    _timestep_tolerance(getParam<Real>("timestep_tolerance"))
{
  for (const auto & name : getParam<std::vector<ExecutorName>>("inner_executors"))
    _inner_executors.push_back(&getExecutorByName<Executor>(name));

  _fe_problem.transient(true);

  // Ensure EVERY solver system (and the aux system) carries an old-solution state
  // (state 1, Time). Otherwise advanceState()/copySolutionsBackwards() and, on
  // cutback, restoreSolutions() fail for a system that has no other consumer of old
  // solutions -- e.g. the quasi-static displacement system has no time-derivative or
  // old-value kernel, so its old-solution vector would never be allocated. Requested
  // here in the ctor (before initialSetup) so the vectors are allocated at init.
  for (std::size_t i = 0; i < _fe_problem.numNonlinearSystems(); ++i)
    _fe_problem.getNonlinearSystemBase(static_cast<unsigned int>(i)).needSolutionState(1);
  _fe_problem.getAuxiliarySystem().needSolutionState(1);
}

Executor::Result
TransientExecutor::run()
{
  auto result = newResult();

  _fe_problem.initialSetup();

  // --- initialize time state ---
  _fe_problem.timeStep() = 0;
  _fe_problem.time() = _start_time;
  _fe_problem.timeOld() = _start_time;
  _fe_problem.dt() = _dt_nominal;
  _fe_problem.dtOld() = _dt_nominal;

  // Output the initial condition and seed old == current so the first step's
  // slnOld() is well defined.
  _fe_problem.outputStep(EXEC_INITIAL);
  _fe_problem.copySolutionsBackwards();

  Real time = _start_time;    // last committed (converged) time
  int t_step = 0;             // last committed step count
  Real dt = _dt_nominal;      // desired dt for the next step
  bool advance_next = false;  // advance state only between converged steps
  bool aborted = false;

  // NB: cast num_steps to long, not int -- the default UINT_MAX cast to int is -1.
  while (t_step < static_cast<long>(_num_steps) && time < _end_time - _timestep_tolerance)
  {
    // Populate old solutions from the previous *converged* step. Not done before
    // the very first step (copySolutionsBackwards handled that) nor on a cutback
    // retry (advance_next stays false there).
    if (advance_next)
    {
      _fe_problem.advanceState();
      advance_next = false;
    }

    // Determine the dt to attempt, clamped to dtmax and end_time.
    Real step_dt = std::min(dt, _dtmax);
    if (time + step_dt > _end_time)
      step_dt = _end_time - time;
    if (step_dt < _dtmin)
    {
      result.fail("dt below dtmin before solve");
      aborted = true;
      break;
    }

    // Set trial time state for this step.
    _fe_problem.dtOld() = _fe_problem.dt();
    _fe_problem.dt() = step_dt;
    _fe_problem.timeOld() = time;
    _fe_problem.time() = time + step_dt;
    _fe_problem.timeStep() = t_step + 1;

    _console << "\nTime Step " << t_step + 1 << ", time = " << time + step_dt
             << ", dt = " << step_dt << std::endl;

    _fe_problem.timestepSetup();
    _fe_problem.onTimestepBegin();
    _fe_problem.execute(EXEC_TIMESTEP_BEGIN);

    // ---- one outer (SPIN) solve of the coupled system ----
    bool step_converged = true;
    for (auto * const inner : _inner_executors)
    {
      auto r = inner->exec();
      if (!r.convergedAll())
        step_converged = false;
    }

    if (!step_converged)
    {
      // Cutback: restore the previous converged solution, roll time/step back,
      // shrink dt and retry the same step. Do NOT advance state on the retry.
      _fe_problem.restoreSolutions();
      _fe_problem.time() = time;
      _fe_problem.timeStep() = t_step;
      dt = step_dt * _cutback_factor;
      _console << " Solve did not converge; cutting back to dt = " << dt << std::endl;
      if (dt < _dtmin)
      {
        result.fail("cutback dt below dtmin");
        aborted = true;
        break;
      }
      continue;
    }

    // ---- step accepted ----
    time += step_dt;
    ++t_step;
    _fe_problem.onTimestepEnd();
    _fe_problem.execute(EXEC_TIMESTEP_END);
    _fe_problem.outputStep(EXEC_TIMESTEP_END);
    advance_next = true;

    // Grow dt back toward the target for the next step.
    dt = std::min(step_dt * _growth_factor, _dt_nominal);
  }

  _fe_problem.execute(EXEC_FINAL);
  _fe_problem.outputStep(EXEC_FINAL);
  _fe_problem.postExecute();

  if (!aborted)
    result.pass("transient completed", /*overwrite=*/true);

  return result;
}
