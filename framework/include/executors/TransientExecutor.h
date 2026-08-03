//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Executor.h"
#include <vector>

class FEProblemBase;

/**
 * A transient executor for the Executor system.
 *
 * Drives its inner executors (e.g. the outer field-split SPIN NewtonSNESExecutor)
 * over a sequence of (pseudo-)time steps. Between converged steps it advances the
 * problem state so that old solutions (`_u_old`, `_var.slnOld()`) are populated --
 * which is what quasi-static phase-field irreversibility (penalty kernel reading
 * d_old) and any time-derivative terms rely on.
 *
 * Modeled on SteadyExecutor + the legacy TransientBase time loop, but self
 * contained: it does NOT use a TimeStepper object (which is tightly coupled to the
 * legacy Transient executioner). Time stepping is a simple fixed target dt with
 * geometric cutback on non-convergence and geometric growth back toward the target
 * after a converged step. Suitable for quasi-static problems (no TimeIntegrator
 * required); TimeIntegrators, if present, are still driven via advanceState() /
 * onTimestepBegin() / the inner solve.
 */
class TransientExecutor : public Executor
{
public:
  static InputParameters validParams();
  TransientExecutor(const InputParameters & params);
  virtual Result run() override;

protected:
  /// The finite element problem providing solve/state/output routines
  FEProblemBase & _fe_problem;

  /// Inner executors run each timestep (e.g. the outer SPIN NewtonSNESExecutor)
  std::vector<Executor *> _inner_executors;

  /// Simulation start time
  const Real _start_time;
  /// Simulation end time
  const Real _end_time;
  /// Target (nominal) time step
  const Real _dt_nominal;
  /// Maximum number of time steps
  const unsigned int _num_steps;
  /// Minimum allowable dt before the run is declared failed
  const Real _dtmin;
  /// Maximum allowable dt
  const Real _dtmax;
  /// Factor to grow dt by after a converged step (capped at the target dt)
  const Real _growth_factor;
  /// Factor to shrink dt by after a non-converged step (cutback)
  const Real _cutback_factor;
  /// Tolerance used when comparing the current time against end_time
  const Real _timestep_tolerance;
};
