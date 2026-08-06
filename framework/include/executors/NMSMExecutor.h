//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "SNESNPCExecutor.h"
#include "MooseEnum.h"

class NewtonSNESExecutor;

/**
 * Executor implementing nonlinear block Gauss-Seidel, also called MSPIN for Multiplicative Schwarz
 * Preconditioned Inexact Newton when used as a preconditioner for an outer Newton solver, by
 * sweeping over a set of sub-SNES executors via a SNESSHELL.
 *
 * sweep_type = multiplicative           forward sweep only (1..N)
 * sweep_type = symmetric_multiplicative forward then backward sweep (1..N, N-1..1)
 */
class NMSMExecutor : public SNESNPCExecutor
{
public:
  static InputParameters validParams();
  NMSMExecutor(const InputParameters & params);
  virtual ~NMSMExecutor();

  virtual Result run() override;
  virtual PetscErrorCode applyBA(Mat A, Vec X, Vec Y) override;

  /// Arm the sub-executor whose nonlinear system is \p sys_num to do a REDUCED (bound-constrained)
  /// solve on the next sweep, holding \p frozen (that system's global DOF indices) at d_old. Used by
  /// the parent coupled PDAS executor so the phase-field block sub-solve respects the shared active set.
  void armBoundedSubSolve(unsigned int sys_num,
                          const std::vector<PetscInt> & frozen,
                          const std::vector<PetscReal> & vals,
                          bool bounded_box = true,
                          bool restrict_assembly = false);
  /// Undo armBoundedSubSolve on all sub-executors (restore normal stock sub-solves).
  void disarmBoundedSubSolve();

  /// NEPIN nonlinear elimination: restrict the sweep to the single block whose system number is
  /// sys_num (the bad-set pf sub-solve); skip the other ("good") blocks and the backward sweep.
  void setNepinOnly(unsigned int sys_num) { _nepin_only_sys = static_cast<int>(sys_num); }
  /// Restore the full multiplicative sweep.
  void clearNepinOnly() { _nepin_only_sys = -1; }

protected:
  virtual void setupSNES() override;

private:
  std::vector<NewtonSNESExecutor *> _sub_snes;
  const MooseEnum _sweep_type;
  /// NEPIN: system number of the sole block to sweep (-1 = full multiplicative sweep, the default).
  int _nepin_only_sys = -1;
  Vec _block_residual = nullptr;
  Vec _block_update = nullptr;

  /// On-demand (SPIN_VERBOSE) one-line summary of a block sub-solve: its iteration count and final
  /// residual (from the sub-system's MOOSE nonlinear-solve stats), indented to nest under the outer
  /// SNES function-norm monitor. \p back tags the backward half of a symmetric sweep.
  void logSubSolve(NewtonSNESExecutor * sub, bool back);

  static PetscErrorCode shellSolveCallback(SNES snes, Vec x);
  PetscErrorCode applyBlockUpdate(Mat A, Vec rhs, Vec Y, PetscInt i);
};
