//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "SNESExecutor.h"
#include "MooseTypes.h"

#include <vector>
#include <map>
#include <memory>
#include <utility>
#include <cstdlib>

namespace libMesh
{
template <typename>
class PetscMatrix;
template <typename>
class NumericVector;
class Elem;
class Node;
}

/**
 * Executor implementing a Newton-type outer solve (SNESEWTONLS) for one or more libMesh
 * nonlinear systems.
 */
class NewtonSNESExecutor : public SNESExecutor
{
public:
  static InputParameters validParams();
  NewtonSNESExecutor(const InputParameters & params);
  virtual ~NewtonSNESExecutor();

  virtual Result run() override;
  virtual SNES getSNES() override;

  /// Global nonlinear-system numbers this executor solves (for a block sub-executor, its single system).
  const std::vector<unsigned int> & nlSysNums() const { return _nl_sys_nums; }

  /**
   * @returns the libMesh system corresponding to the index \p i
   */
  libMesh::System & getSystem();

protected:
  virtual void setupSNES() override;

private:
  /// The nonlinear systems this SNES targets
  std::vector<unsigned int> _nl_sys_nums;

  /// Name of the postprocessor giving the total potential energy Psi(x); when set, the outer SPIN
  /// line search uses the energy merit (paper App. B/C) instead of the residual-norm merit.
  const PostprocessorName _energy_pp_name;
  /// Whether an energy postprocessor was provided (selects the line-search merit)
  const bool _has_energy_pp;

  /// If true (requires _has_energy_pp), the energy-merit path uses a ratio-based trust region
  /// instead of the App. B/C strong-Wolfe line search.
  const bool _use_trust_region;

  /// Trust-region subproblem solver: true = SPIN-preconditioned Steihaug-Toint truncated CG over
  /// the full space (handles indefinite Hessians); false = 1-D model along the single SPIN direction.
  const bool _tr_steihaug;

  /// True when a per-iteration NPC "sweep" is applied inside the Steihaug TR (i.e. npc != mono):
  /// either a full MSPIN field sweep (npc=mspin) or a NEPIN hard-set elimination (npc=mspin_pne/fne).
  /// Resolved in the ctor from the npc enum (or the deprecated tr_npc bool). false = MONO.
  bool _tr_npc;

  // --- NEPIN nonlinear-elimination layer (INB-NE; Liu-Hwang-Luo-Cai-Keyes, SISC 2022) ----------
  /// The tr_npc "sweep" is a nonlinear ELIMINATION of the hard (strongly nonlinear) set instead of the full
  /// MSPIN field sweep. _nepin = any NE (partial or full). _nepin_full = the coupled front-band variant
  /// (npc=mspin_fne): the displacement solve is ALSO restricted to the front band (both fields move
  /// together); _nepin && !_nepin_full = the d-only partial variant (npc=mspin_pne): disp is left to
  /// the outer TR. Resolved in the ctor from the npc enum.
  bool _nepin;
  bool _nepin_full;
  /// Process-zone band defining the hard (eliminated) damage set: nepin_band_lo < d < 1 - nepin_band_hi.
  /// Everything else (far-field d~0, broken d~1, irreversibly-frozen front) is "easy" and pinned.
  const Real _nepin_lo;
  const Real _nepin_hi;
  /// (FULL only) displacement hard-set relative threshold: a disp DOF is eliminated (freed on the band)
  /// iff |R_u,i| > nepin_u_rtol * max_i|R_u,i| -- the front force-imbalance selects the coupled front.
  const Real _nepin_u_rtol;
  /// NEPIN[FULL] cost option: restrict each block sub-solve's FE assembly to the elements incident to
  /// its free set (crack-front band). vinewtonrsls already solves only the free reduced space, so this
  /// makes assembly band-sized too, reproducing the full-assembly iterates bit-for-bit. Default true.
  const bool _nepin_restrict_assembly;
  /// Damage easy (frozen) DOFs + current values -- the complement of the damage band, handed to
  /// armReducedSolve (which pins these and frees the rest in [d_old, 1]). Filled by computeHardSet.
  std::vector<PetscInt> _nepin_frozen;
  std::vector<PetscReal> _nepin_frozen_val;
  /// (FULL only) displacement easy (frozen far-field) DOFs + current values -- pinned so only the front
  /// band of disp is free. The frozen ring is the Dirichlet anchor that removes the rigid-body modes of
  /// the free patch (this is NEPIN's eq 2.24 pinning of the easy set). Filled by computeHardSet.
  std::vector<PetscInt> _nepin_frozen_u;
  std::vector<PetscReal> _nepin_frozen_u_val;
  /// Displacement system (the non-bounded coupled system) global number + local index in this
  /// executor's VecNest. Resolved alongside _bounded_sys_* for the FULL disp elimination.
  unsigned int _disp_sys_num = libMesh::invalid_uint;
  unsigned int _disp_sys_local = libMesh::invalid_uint;

  // --- NEPIN[FULL] band-restricted assembly (nepin_restrict_assembly) ---------------------------
  /// Backing store for the band-restricted assembly range: the locally-owned elements incident to a
  /// block's free set. Filled per block sub-solve in reducedNewtonSolve. Must persist across the solve
  /// (libMesh ConstElemRange stores iterators into this, not a copy).
  std::vector<const libMesh::Elem *> _band_elems;

  // --- PDAS bound-constraint layer (design (b), Heister-Wheeler-Wick Alg. 3.2) ------------------
  /// If true, enforce the irreversibility bound d_old <= d <= 1 on _bounded_var_name via a
  /// primal-dual active set merged into the outer solve, applied consistently to the coupled operator
  /// AND every field-split preconditioner block.
  const bool _bounds;
  /// Name of the bounded (phase-field) variable; its owning nonlinear system carries the constraint.
  const NonlinearVariableName _bounded_var_name;
  /// Primal-dual active-set complementarity constant c > 0. Retained for input compatibility (unused
  /// by the projected-gradient active set).
  const Real _pdas_c;
  /// Multiplier dead-band: a bounded DOF at its bound joins the active set only if |B^{-1} R| exceeds
  /// this. Suppresses roundoff-level active-set chatter in the unstressed pre-damage region.
  const Real _pdas_lambda_tol;
  /// Name of the matrix tag holding the (consistent) mass matrix on the bounded variable.
  const TagName _mass_tag_name;
  /// On-demand verbose diagnostics (env SPIN_VERBOSE): print the one-line-per-category solver traces
  /// ([PDAS], [PDAS sweep], [spinTRcg], [gradPsi conv], NPC-sweep headers). Default off -> only the
  /// outer SNES function norm is logged.
  const bool _verbose = (std::getenv("SPIN_VERBOSE") != nullptr);

  /// Resolved (once, lazily at first solve): global nonlinear-system number and variable number of
  /// the bounded variable, its local index within this executor's VecNest, and the mass matrix tag.
  unsigned int _bounded_sys_num = libMesh::invalid_uint;
  unsigned int _bounded_var_num = libMesh::invalid_uint;
  unsigned int _bounded_sys_local = libMesh::invalid_uint;
  TagID _mass_tag = Moose::INVALID_TAG_ID;

  /// Cached B^{-1}: reciprocal of the lumped (row-summed) mass diagonal on the bounded variable's
  /// system. Computed once (geometry-only) and reused every iteration for the multiplier estimate
  /// lambda ~ B^{-1} R. Layout matches the bounded variable's system solution vector.
  std::unique_ptr<libMesh::NumericVector<libMesh::Number>> _binv;
  /// Whether _binv has been built yet.
  bool _mass_built = false;

  /// Assemble the mass matrix on the bounded variable once, lump it (row sum), invert, and cache the
  /// result in _binv. Also resolves _bounded_sys_num/_bounded_var_num/_bounded_sys_local/_mass_tag.
  void buildLumpedMassInverse();

  /// PDAS active set on the bounded variable's DOFs, recomputed each outer iteration (block-global
  /// indices in the bounded system's index space, this rank's owned dofs). _active_d_prev holds the
  /// previous iteration's set for the HWW active-set-unchanged convergence test.
  std::vector<PetscInt> _active_d_dofs;
  std::vector<PetscInt> _active_d_prev;
  /// Freeze value for each entry of _active_d_dofs: the bound the DOF is pinned to (d_old for a
  /// lower-active DOF, 1 for an upper-active DOF). Aligned index-for-index with _active_d_dofs.
  std::vector<PetscReal> _active_val;
  /// Whether the active set changed on the last computeActiveSet (HWW converges only when it does not
  /// change AND ||R|| < tol). Starts true so a solve cannot "converge" before any set is computed.
  bool _active_set_changed = true;
  /// Mask over the bounded variable's system (1.0 on inactive DOFs, 0.0 on active), used to zero the
  /// active columns of the (u,d) off-diagonal block (MatDiagonalScale) and the active residual entries.
  Vec _active_mask = nullptr;
  /// Previous accepted step's increment on the bounded variable (delta d^{k-1}) for the criterion;
  /// reset to zero at the start of each solve.
  Vec _pdas_prev_dstep = nullptr;
  /// Scratch for the HIK merit safeguard on the NPC sweep (bounds=true, tr_npc): _x_presweep holds the
  /// pre-sweep iterate X, _sweep_dir holds NPC(X) - X so the sweep can be backtracked on the KKT merit.
  Vec _x_presweep = nullptr;
  Vec _sweep_dir = nullptr;
  /// Scratch holding the reduced residual R_I(X) across the trust-region trial loop (bounds + energy
  /// merit), so the decoupled KKT-progress test can evaluate ||R_I(W)|| (which overwrites _r_plain) and
  /// then restore R_I(X) for the next Steihaug re-solve.
  Vec _r_stash = nullptr;

  /// Recompute the PDAS active set from R (= _r_plain), the current iterate \p X, d_old, and the
  /// previous step. Updates _active_d_dofs, _active_mask, and _active_set_changed. Active criterion
  /// (raccoon d >= d_old convention): (B^{-1})_ii R_i + c (d_old - d^k - delta d^{k-1})_i > 0.
  void computeActiveSet(Vec X);
  /// NEPIN: identify the hard (eliminated) damage set -- the process-zone band on the current iterate
  /// \p X -- and fill _nepin_frozen / _nepin_frozen_val with its complement (the easy DOFs pinned at
  /// their current value) for armReducedSolve. Solving only the hard band is the nonlinear elimination.
  void computeHardSet(Vec X);
  /// Dynamic Dirichlet elimination of the active set from the coupled Hessian _mat_nest and the
  /// residual _r_plain: zero active rows+cols (diag 1) of the d-diagonal block, active rows of the
  /// (d,u) block, active cols of the (u,d) block, and active entries of R -- so the reduced step has
  /// delta d = 0 on the active set. Must run after assembleCoupledJacobian and before steihaugTRS.
  PetscErrorCode applyActiveSetElimination();
  /// Project the bounded variable of the coupled iterate \p X onto the active-set freeze values
  /// (d_old on lower-active DOFs, 1 on upper-active DOFs) -- the feasibility projection of a
  /// bound-constrained active-set step. Run after computeActiveSet, before assembling at X.
  void clampActiveToBounds(Vec X);
  /// Project a trial coupled iterate \p W onto the box d_old <= d <= 1 (clip the whole bounded-variable
  /// block). Applied to EVERY trial in the TR and LS globalizations so a bound-constrained step cannot
  /// overshoot. No-op unless _bounds. Shared PDAS treatment across globalizations.
  void pdasProjectTrial(Vec W);
  /// Record the accepted step's bounded-variable increment (delta d^k) into _pdas_prev_dstep for the
  /// next iteration's active-set criterion. No-op unless _bounds.
  void pdasRecordStep(Vec step);

public:
  /// Arm this (single-system) sub-executor to solve the REDUCED subsystem on its next run(): a Newton
  /// solve of its nonlinear system with the DOFs in \p frozen held at the values in \p vals (their
  /// rows/cols eliminated). Called by the parent coupled executor (via NMSMExecutor) before an NPC
  /// sweep so the phase-field block sub-solve respects the shared two-sided PDAS active set. \p frozen
  /// are global DOF indices in this executor's (single) nonlinear system; \p vals the bound each is
  /// pinned to (d_old or 1), aligned with \p frozen.
  /// \p bounded_box selects the FREE-DOF box: true = the damage irreversibility box [d_old, 1] (pf
  /// block); false = unbounded [-inf, +inf] (displacement block, NEPIN[FULL]) so the free DOFs solve
  /// as plain Newton while the frozen ring is pinned (a Dirichlet-anchored reduced solve).
  /// \p restrict_assembly (NEPIN[FULL], nepin_restrict_assembly): restrict the block's FE assembly to the
  /// elements incident to its free (non-frozen) set, so residual/Jacobian cost scales with the band.
  void armReducedSolve(const std::vector<PetscInt> & frozen,
                       const std::vector<PetscReal> & vals,
                       bool bounded_box = true,
                       bool restrict_assembly = false);
  /// Disarm the reduced solve (restore the normal stock sub-solve).
  void disarmReducedSolve();

private:
  /// Whether run() should do a reduced (frozen-set) Newton solve instead of the stock _fe_problem
  /// solve, and the frozen DOF set (this system's global indices).
  bool _reduced_solve_armed = false;
  /// Free-DOF box for the reduced solve: true = [d_old, 1] (pf), false = unbounded (disp). See above.
  bool _reduced_bounded_box = true;
  /// NEPIN[SUBDOMAIN]: restrict this block's FE assembly to the band incident to its free set.
  bool _reduced_restrict_assembly = false;
  std::vector<PetscInt> _frozen_dofs;
  /// Freeze value for each _frozen_dofs entry (bound the DOF is pinned to; aligned with _frozen_dofs).
  std::vector<PetscReal> _frozen_vals;

  /// Newton solve of this single nonlinear system with _frozen_dofs held at d_old (their rows/cols
  /// eliminated from the Jacobian and residual). Used for the bound-constrained NPC phase-field block.
  void reducedNewtonSolve();

  /// Trust-region radius (P-norm for Steihaug, solution-norm for the 1-D variant), persistent across
  /// outer iterations of a solve and reset at the start of each solve. Negative => uninitialized.
  Real _tr_radius = -1.0;

  /// Block-Jacobi field-split preconditioner solves (A00^{-1}, A11^{-1}) for the Steihaug TRS.
  KSP _pc_ksp0 = nullptr;
  KSP _pc_ksp1 = nullptr;
  /// CG scratch VecNests for the Steihaug TRS (lazily duplicated from _vec_func).
  Vec _cg_p = nullptr;
  Vec _cg_r = nullptr;
  Vec _cg_y = nullptr;
  Vec _cg_d = nullptr;
  Vec _cg_Ad = nullptr;
  Vec _cg_Pp = nullptr;
  Vec _cg_Pd = nullptr;

  /// PETSc MatNest holding all (i,j) Jacobian sub-blocks.
  Mat _mat_nest = nullptr;

  /// KSP on the full Jacobian _mat_nest for the App. C full-J inexact-Newton fallback direction
  /// (energy-merit line search only). Created lazily when the energy merit is in use.
  KSP _fullJ_ksp = nullptr;

  /// Dedicated VecNest holding the plain (unpreconditioned) residual R = grad Psi, assembled by
  /// computePlainResidual(). The outer SNES function vector F is a SEPARATE vector holding the
  /// preconditioned residual F_SPIN, so the energy line search must use _r_plain for grad Psi.
  Vec _r_plain = nullptr;

  /// Stash of the reduced R(X) for the RESIDUAL-merit trust region, preserved across the ratio-test
  /// retry loop: steihaugTRS re-reads _r_plain=R(X) on every rejected re-solve, but the per-trial
  /// residual eval (computeMerit(W)) overwrites _r_plain with R(W). Allocated only for the
  /// residual-merit Steihaug path (TR steihaug without an energy postprocessor).
  Vec _r_base = nullptr;

  /// ||R|| = ||grad Psi|| cached from the last line-search call, and its value at the start of the
  /// current outer solve. The energy-merit path converges on ||grad Psi|| (the nonlinear
  /// preconditioner pre-minimizes the energy, so ||x - NPC(x)|| plateaus while grad Psi -> 0).
  Real _r_norm = -1.0;
  Real _r0_norm = -1.0;

  struct MatData
  {
    std::unique_ptr<libMesh::PetscMatrix<libMesh::Number>> mat;
    TagID tag;
  };

  /// Off-diagonal Jacobian blocks J_ij (i != j).
  std::map<std::pair<unsigned int, unsigned int>, MatData> _off_diag_mats;

  // Multi-system helpers -----------------------------------------------------------

  void allocateOffDiagMats();
  void assembleOffDiagJacobian();
  void buildMatNest();

  /// Apply the SPD block-Jacobi (additive) field-split preconditioner: out = diag(A00,A11)^{-1} in.
  PetscErrorCode applyBlockJacobi(Vec in, Vec out);
  /// Apply the SPD symmetric block Gauss-Seidel preconditioner (the symmetrized MSPIN; uses the
  /// coupling, much stronger than block-Jacobi): out = P_SGS^{-1} in with
  /// P_SGS = (D+L) D^{-1} (D+U), D=diag(A00,A11), L=[[0,0],[A10,0]], U=L^T. Applied as
  /// (D+U)^{-1} D (D+L)^{-1} = forward block sweep, block-diagonal scale, backward block sweep.
  PetscErrorCode applyBlockSGS(Vec in, Vec out);
  /// Solve the trust-region subproblem  min g^T p + 1/2 p^T A p  s.t.  ||p||_P <= Delta  by
  /// preconditioned Steihaug-Toint truncated CG (A = _mat_nest, P = block-Jacobi). Returns the step
  /// in _cg_p, the model reduction pred = -(g^T p + 1/2 p^T A p), and whether the step is on the
  /// trust-region boundary (radius hit or negative curvature). \p eps is the inner forcing tolerance.
  PetscErrorCode
  steihaugTRS(Vec g, Real Delta, Real eps, Real & pred, bool & on_boundary, PetscInt & cg_its);

  /// Scatter a solution VecNest's sub-vecs into the libMesh nonlinear system solution vectors
  /// (and update() them). Used to place a trial outer iterate before evaluating a merit.
  void scatterToSystems(Vec x);
  /// Total potential energy Psi at the outer iterate \p x: scatters \p x, runs EXEC_LINESEARCH so
  /// the energy postprocessor re-integrates, and returns its value. Requires _has_energy_pp.
  Real computeEnergy(Vec x);
  /// Merit at the outer iterate \p x for the trust-region ratio test: the total potential energy
  /// Psi(x) if an energy postprocessor is set, else the residual-norm merit 1/2||R(x)||^2 (reduced by
  /// the active-set mask when _bounds). NOTE: computeEnergy leaves _r_plain untouched, but the
  /// residual branch OVERWRITES _r_plain with R(x); callers preserve R(X) via _r_base across the TR
  /// retry loop (see spinLineSearch's Steihaug block).
  Real computeMerit(Vec x);
  /// KKT-residual merit 1/2||mask (.) R(x)||^2 on the reduced (inactive) space of the CURRENT active-set
  /// mask. Independent of _has_energy_pp: this is the semismooth-Newton merit for the bound-constrained
  /// KKT system, used to globalize the NPC sweep (the energy Psi can fall while the KKT residual rises,
  /// so Psi is the wrong merit for the active-set safeguard). OVERWRITES _r_plain with R(x).
  Real reducedResidualMerit(Vec x);
  /// Assemble the plain (unpreconditioned) residual F(x) at \p x into _vec_func (== the outer SNES
  /// function vector). This is grad Psi(x), used for the energy line search's slope and App. C.
  void computePlainResidual(Vec x);

  /// Cubic energy-Armijo backtracking along direction \p Y from base iterate \p X with initial slope
  /// \p slope (= grad Psi . (-Y) < 0). On return \p W = X - lambda_out*Y for the accepted step.
  /// Returns true if sufficient decrease was achieved. Used by each App. B/C direction tier.
  bool energyBacktrack(Vec X,
                       Vec Y,
                       Vec W,
                       Real psiX,
                       Real slope,
                       Real minlambda,
                       PetscInt max_it,
                       Real & lambda_out,
                       Real & psiW_out);

  /// Strong-Wolfe line search (paper App. B: cubic bracketing + zoom, c1=1e-4, c2=0.9) on the energy
  /// merit along direction \p Y from \p X with phi'(0) = \p phi0p (= grad Psi . (-Y) < 0). phi(alpha)
  /// = Psi(X - alpha*Y) via computeEnergy; phi'(alpha) = -R(X-alpha*Y).Y via computePlainResidual.
  /// The curvature condition |phi'(alpha)| <= c2|phi'(0)| forbids the tiny steps Armijo alone would
  /// accept. On return \p W = X - lambda_out*Y. Returns true if an acceptable step was found.
  bool strongWolfe(Vec X,
                   Vec Y,
                   Vec W,
                   Real psiX,
                   Real phi0p,
                   Real minlambda,
                   PetscInt max_it,
                   Real & lambda_out,
                   Real & psiW_out);

  static PetscErrorCode outerResidualCallback(SNES snes, Vec x, Vec f, void * ctx);
  static PetscErrorCode outerJacobianCallback(SNES snes, Vec x, Mat A, Mat P, void * ctx);

  /// Assemble the coupled Jacobian _mat_nest = grad^2 Psi at the outer iterate \p x: scatter x into
  /// the systems, compute each diagonal block J_ii (+ hanging-node constraint enforcement), assemble
  /// the off-diagonal coupling blocks, finalize _mat_nest, and rebind the per-field sub-KSP and the
  /// Steihaug field-split preconditioner operators to the fresh diagonal blocks. Shared by
  /// outerJacobianCallback (non-B modes) and the (B) NPC-warm-started Steihaug branch (which needs
  /// the Hessian re-assembled at NPC(x)).
  PetscErrorCode assembleCoupledJacobian(Vec x);

  //
  // shell machinery for an outer SNES
  //

  /// Shell matrix representing multiplication of the linearized nonlinear
  /// preconditioner and the full Jacobian. I think of it as representing B*A,
  /// where B is what PETSc typically uses to denote the preconditioner (*not*
  /// what MOOSE calls the preconditioning matrix P, but something more like
  /// P^{-1}). The nonlinear preconditioner supplies the exact block application
  /// used by this shell matrix.
  Mat _jac_shell = nullptr;

  /**
   * Shell routine associated with the _jac_shell that will forward to the nonlinear preconditioner
   * for applying the preconditioning at the linearized operator level and for applying our global
   * MatNest Jacobian
   */
  static PetscErrorCode shellMatMult(Mat m, Vec X, Vec Y);

  /**
   * Custom shell line search for the outer SPIN solve. The outer SNES runs with
   * SNES_FUNCTION_PRECONDITIONED + a left NPC, so its residual is the preconditioned
   * SPIN residual F_SPIN(x) = x - NPC(x). The stock line search evaluates trial points
   * with the plain (unpreconditioned) residual, giving an inconsistent merit that fails
   * to backtrack once the coupling turns on. This implements cubic backtracking with the
   * CONSISTENT merit 1/2||F_SPIN(W)||^2 (via SNESApplyNPC at each trial) and Armijo
   * sufficient decrease (paper App. B, c1 = 1e-4).
   */
  static PetscErrorCode spinLineSearch(SNESLineSearch ls, void * ctx);

  /**
   * Custom outer convergence test for the energy-merit path: converges on ||R|| = ||grad Psi||
   * (cached from the line search) rather than the preconditioned residual ||x - NPC(x)||. The
   * nonlinear preconditioner pre-minimizes the energy block-wise, so the fixed-point residual
   * plateaus while grad Psi -> 0; the energy line search drives grad Psi to zero, so that is the
   * consistent, physically meaningful convergence measure.
   */
  static PetscErrorCode outerConvergenceTest(SNES snes,
                                             PetscInt it,
                                             PetscReal xnorm,
                                             PetscReal ynorm,
                                             PetscReal fnorm,
                                             SNESConvergedReason * reason,
                                             void * ctx);
};
