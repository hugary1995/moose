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

namespace libMesh
{
template <typename>
class PetscMatrix;
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

  /// Nonlinear-preconditioner usage inside the Steihaug trust region (only meaningful when
  /// _use_trust_region && _tr_steihaug). false = (A) NO NPC: plain coupled Newton-TR on grad^2 Psi
  /// (monolithic-family). true = (B) NPC EVERY outer iteration: apply one multiplicative-Schwarz
  /// sweep X <- NPC(X), then a Steihaug-TR step on grad^2 Psi(NPC(X)) -- the "TR version of MSPIN",
  /// same nonlinear-preconditioning cadence as the App. C line search but TR globalization.
  const bool _tr_npc;

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
