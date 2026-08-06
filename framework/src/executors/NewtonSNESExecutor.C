//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NewtonSNESExecutor.h"
#include "Convergence.h"
#include "FEProblem.h"
#include "NonlinearSystem.h"
#include "NonlinearSystemBase.h"
#include "SNESNPCExecutor.h"
#include "NMSMExecutor.h"
#include "NodalBCBase.h"
#include "MooseVariableBase.h"
#include "MooseMesh.h"
#include "BndNode.h"
#include "Moose.h"

#include "libmesh/libmesh.h"
#include "libmesh/petsc_solver_exception.h"
#include "libmesh/implicit_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/mesh_base.h"
#include "libmesh/elem.h"
#include "libmesh/dof_map.h"
#include "libmesh/remote_elem.h"

#include <petscsnes.h>
#include <petscsys.h>

#include <array>
#include <set>
#include <string>

registerMooseObject("MooseApp", NewtonSNESExecutor);

namespace
{
// libMesh nonlinear-solver bounds callback for the bound-constrained NPC phase-field sub-solve. It
// copies the system's (MOOSE-owned) lower_bound/upper_bound vectors -- filled by the executor with
// the shared PDAS set -- into the SNES's VI bound vectors XL/XU. Using MOOSE/libMesh vectors here (no
// raw PETSc, no shared Vec ownership) avoids the double-free that SNESVISetVariableBounds caused.
void
pdasCopyBounds(libMesh::NumericVector<libMesh::Number> & XL,
               libMesh::NumericVector<libMesh::Number> & XU,
               libMesh::NonlinearImplicitSystem & sys)
{
  const auto & lb = sys.get_vector("lower_bound");
  const auto & ub = sys.get_vector("upper_bound");
  for (auto i = XL.first_local_index(); i < XL.last_local_index(); ++i)
  {
    XL.set(i, lb(i));
    XU.set(i, ub(i));
  }
  XL.close();
  XU.close();
}
}

InputParameters
NewtonSNESExecutor::validParams()
{
  InputParameters params = SNESExecutor::validParams();
  params.addClassDescription(
      "Newton-type outer solver executor (SNESEWTONLS). "
      "Delegates to _fe_problem.solve() for a single nonlinear system, or builds a combined "
      "outer SNES with VecNest/MatNest when multiple systems are present (nl_preconditioning "
      "required in that path).");
  params += Moose::PetscSupport::kspRelatedParams();
  params += Moose::PetscSupport::newtonKrylovParams();
  params.addRequiredParam<std::vector<NonlinearSystemName>>(
      "nonlinear_system_names", "Name of the nonlinear systems this executor targets.");
  params.addParam<std::vector<ConvergenceName>>(
      "convergence_names",
      "Convergence object name(s) for each nonlinear system. If provided, must have the same "
      "length as 'nonlinear_system_names'.");
  params.addParam<PostprocessorName>(
      "energy_postprocessor",
      "Postprocessor giving the total potential energy Psi(x). If set, the outer SPIN line search "
      "uses the energy as its merit (paper App. B/C: it is C^1 across the penalty kink and admits "
      "a descent-direction fallback). If unset, the line search uses the preconditioned residual "
      "norm merit 1/2||x - NPC(x)||^2.");
  params.addParam<bool>(
      "use_trust_region",
      false,
      "If true (requires energy_postprocessor), globalize the outer SPIN solve with a ratio-based "
      "TRUST REGION on the energy merit instead of the App. B/C strong-Wolfe line search. The "
      "actual/predicted energy-reduction ratio accepts/shrinks the step and adapts the radius, "
      "degrading gracefully at crack propagation rather than spiraling the time step to dtmin. See "
      "'trust_region_solver' for how the trust-region subproblem is solved.");
  params.addParam<MooseEnum>(
      "trust_region_solver",
      MooseEnum("steihaug tr_linesearch", "steihaug"),
      "Trust-region subproblem solver (when use_trust_region=true): 'steihaug' = SPIN-preconditioned "
      "Steihaug-Toint truncated CG over the full space (explores the Krylov subspace, handles "
      "indefinite Hessians via a negative-curvature-to-boundary exit); 'tr_linesearch' = a 1-D "
      "quadratic energy model along the single SPIN direction (cheaper, but stalls when that "
      "direction is nearly orthogonal to grad(Psi), e.g. during crack propagation).");
  params.addParam<bool>(
      "tr_npc",
      false,
      "Nonlinear-preconditioner usage inside the Steihaug trust region (use_trust_region=true, "
      "steihaug). false = (A) NO nonlinear preconditioning: plain coupled Newton trust region on "
      "grad^2 Psi (monolithic-family). true = (B) apply one multiplicative-Schwarz NPC sweep "
      "X <- NPC(X) EVERY outer iteration, then take the Steihaug-TR step at NPC(X) -- the "
      "'trust-region version of MSPIN' (same nonlinear-preconditioning cadence as the App. C line "
      "search, but TR globalization).");
  params.addParam<MooseEnum>(
      "npc",
      MooseEnum("mono mspin mspin_pne mspin_fne", "mono"),
      "Nonlinear preconditioner / per-iteration sweep applied inside the Steihaug trust region "
      "(the coupled-solver FAMILY). 'mono' = no NPC (plain coupled Newton-TR). 'mspin' = one full "
      "multiplicative-Schwarz field sweep X<-NPC(X) every outer iteration (TR version of MSPIN). "
      "'mspin_pne' = NEPIN[PARTIAL]: partial nonlinear elimination -- the sweep eliminates only the "
      "hard DAMAGE band (nepin_band_lo < d < 1-nepin_band_hi), freezing all other damage DOFs at their "
      "current value and solving the band alone (pf only); displacement is left to the outer TR. "
      "'mspin_fne' = NEPIN[FULL]: full nonlinear elimination -- BOTH the damage band and the coupled "
      "displacement front (|R_u|>nepin_u_rtol*max|R_u|) are eliminated together (block Gauss-Seidel of "
      "vinewtonrsls reduced-space sub-solves, far-field frozen). By default (nepin_restrict_assembly) "
      "each block sub-solve restricts its FE assembly to the crack-front band, so cost scales with the "
      "band -- identical numerics to the full-assembly variant. If npc is unset, it is derived from "
      "the deprecated 'tr_npc' bool (true->mspin, false->mono).");
  params.addParam<bool>(
      "nepin_restrict_assembly",
      true,
      "NEPIN[FULL] (npc=mspin_fne): restrict each block sub-solve's FE assembly to the elements "
      "incident to its free set (crack-front band), so residual/Jacobian cost scales with the band "
      "rather than the whole mesh. The vinewtonrsls reduced-space solve is unchanged, so this is a pure "
      "cost optimization -- the iterates are bit-identical to the full-assembly path. Set false to "
      "assemble the full domain (reference / for profiling comparisons).");
  params.addParam<Real>("nepin_band_lo",
                        1e-2,
                        "NEPIN hard-set lower threshold: damage DOFs with d > nepin_band_lo (and below "
                        "the upper threshold) are eliminated. Below this they are far-field 'easy'.");
  params.addParam<Real>("nepin_band_hi",
                        1e-2,
                        "NEPIN hard-set upper threshold: damage DOFs with d < 1 - nepin_band_hi (and "
                        "above the lower threshold) are eliminated. Above this they are broken 'easy'.");
  params.addParam<Real>("nepin_u_rtol",
                        1e-2,
                        "NEPIN[FULL] displacement hard-set relative threshold: a disp DOF is eliminated "
                        "(freed on the front band) iff |R_u,i| > nepin_u_rtol * max_i|R_u,i|. Larger "
                        "-> narrower disp band. Unused for mono/mspin/mspin_pne.");
  params.addParam<bool>(
      "bounds",
      false,
      "Enforce the irreversibility bound constraint d_old <= d <= 1 on 'bounded_variable' via a "
      "primal-dual active set (Heister-Wheeler-Wick Alg. 3.2) merged into the outer solve. The same "
      "active set is applied to the coupled operator and every field-split preconditioner block, so "
      "the nonlinear preconditioner is consistent with the constraint.");
  params.addParam<NonlinearVariableName>(
      "bounded_variable",
      "The bounded (phase-field) variable when bounds=true; its owning nonlinear system carries the "
      "irreversibility constraint d_old <= d <= 1.");
  params.addParam<Real>(
      "pdas_c",
      1.0,
      "Primal-dual active-set complementarity constant c > 0 (bounds=true). Retained for input "
      "compatibility; the projected-gradient active set uses the multiplier sign + dead-band "
      "'pdas_lambda_tol', not c.");
  params.addParam<Real>(
      "pdas_lambda_tol",
      1e-8,
      "Multiplier dead-band (bounds=true): a bounded DOF at its bound joins the active set only if its "
      "lumped-mass multiplier estimate |B^{-1} R| exceeds this tolerance. Excludes numerically "
      "indifferent DOFs (multiplier ~ roundoff, e.g. the unstressed pre-damage region) that would "
      "otherwise oscillate in/out of the set and inflate the outer-iteration count.");
  params.addParam<TagName>(
      "mass_matrix_tag",
      "mass",
      "Name of the matrix tag holding the (consistent) mass matrix on 'bounded_variable' (bounds="
      "true), used to form the lumped-mass multiplier scaling B^{-1}. Provide a MassMatrix kernel on "
      "the bounded variable targeting this tag plus '[Problem] extra_tag_matrices'.");
  return params;
}

NewtonSNESExecutor::NewtonSNESExecutor(const InputParameters & params)
  : SNESExecutor(params),
    _energy_pp_name(isParamValid("energy_postprocessor")
                        ? getParam<PostprocessorName>("energy_postprocessor")
                        : PostprocessorName("")),
    _has_energy_pp(isParamValid("energy_postprocessor")),
    _use_trust_region(getParam<bool>("use_trust_region")),
    _tr_steihaug(getParam<MooseEnum>("trust_region_solver") == "steihaug"),
    _nepin_lo(getParam<Real>("nepin_band_lo")),
    _nepin_hi(getParam<Real>("nepin_band_hi")),
    _nepin_u_rtol(getParam<Real>("nepin_u_rtol")),
    _nepin_restrict_assembly(getParam<bool>("nepin_restrict_assembly")),
    _bounds(getParam<bool>("bounds")),
    _bounded_var_name(_bounds ? getParam<NonlinearVariableName>("bounded_variable")
                              : NonlinearVariableName("")),
    _pdas_c(getParam<Real>("pdas_c")),
    _pdas_lambda_tol(getParam<Real>("pdas_lambda_tol")),
    _mass_tag_name(getParam<TagName>("mass_matrix_tag"))
{
  // Resolve the coupled-solver FAMILY (npc enum) into the internal sweep flags. The enum is the modern
  // single mutually-exclusive selector; if unset, fall back to the deprecated tr_npc bool so existing
  // inputs (mono_tr.i / mspin_tr.i) keep working. NEPIN is a hard-set-restricted MSPIN sweep, so it
  // lives on the same axis as mono/mspin -- not a separate composable knob.
  MooseEnum npc = getParam<MooseEnum>("npc");
  if (!isParamSetByUser("npc") && isParamSetByUser("tr_npc"))
    npc = getParam<bool>("tr_npc") ? "mspin" : "mono";
  _tr_npc = (npc != "mono");
  _nepin = (npc == "mspin_pne" || npc == "mspin_fne");
  _nepin_full = (npc == "mspin_fne");

  // The trust region can globalize either the energy merit Psi (requires energy_postprocessor) or, for
  // the Steihaug variant, the residual merit 1/2||R||^2 (no energy postprocessor needed -- used where
  // R != grad Psi, e.g. irreversible mixed-mode CZM). The 1-D (non-steihaug) TR is energy-only.
  if (_use_trust_region && !_tr_steihaug && !_has_energy_pp)
    paramError("use_trust_region",
               "the 1-D (non-steihaug) trust region requires energy_postprocessor (energy merit). Use "
               "trust_region_solver=steihaug for the residual-merit trust region.");
  if (_bounds && !isParamValid("bounded_variable"))
    paramError("bounded_variable", "bounds=true requires 'bounded_variable'.");
  if (_bounds && _pdas_c <= 0.0)
    paramError("pdas_c", "must be > 0.");
  // NEPIN reuses the bound-constrained NPC pf sub-solve (vinewtonrsls + the shared PDAS box) to solve
  // the free band, so it requires the same bounds=true + steihaug TR path as MSPIN.
  if (_nepin && !(_bounds && _use_trust_region && _tr_steihaug))
    paramError("npc",
               "npc=mspin_pne/mspin_fne (NEPIN) requires bounds=true, use_trust_region=true, "
               "trust_region_solver=steihaug (same path as npc=mspin).");
  // I don't actually know if this is possible with the parser like if the user passes an empty
  // string
  const auto & nl_sys_names = getParam<std::vector<NonlinearSystemName>>("nonlinear_system_names");
  if (nl_sys_names.empty())
    paramError("nonlinear_system_names", "Empty string passed?");
  for (const auto & nl_sys_name : nl_sys_names)
    _nl_sys_nums.push_back(_fe_problem.nlSysNum(nl_sys_name));

  // For the bound-constrained PDAS path, add lower_bound/upper_bound vectors to the target systems
  // now (before EquationSystems init). The NPC phase-field sub-solve fills them with the shared active
  // set and its bounds callback (pdasCopyBounds) feeds them to vinewtonrsls -- MOOSE's native VI path.
  // PARALLEL (non-ghosted): VI bounds only need owned entries. GHOSTED corrupts in parallel because
  // this ctor runs before the mesh is partitioned, so the ghost pattern is built wrong.
  if (_bounds)
    for (const auto nl : _nl_sys_nums)
    {
      auto & s = _fe_problem.getNonlinearSystemBase(nl);
      s.addVector("lower_bound", false, libMesh::PARALLEL);
      s.addVector("upper_bound", false, libMesh::PARALLEL);
    }

  //
  // Store PETSc options
  //

  if (_nl_sys_nums.size() > 1)
    Moose::PetscSupport::storePetscOptions(_fe_problem, this->name() + "_", *this);
  else
  {
    const auto sys_num = _nl_sys_nums[0];
    const auto & sys = _fe_problem.getNonlinearSystemBase(sys_num);
    Moose::PetscSupport::storePetscOptions(_fe_problem, sys.prefix(), *this);
    auto & solver_params = _fe_problem.solverParams(sys_num);
    solver_params._prefix = sys.prefix();
    solver_params._solver_sys_num = sys_num;
  }
  Moose::PetscSupport::setESLinearSolverParams(_fe_problem.es(), *this);

  //
  // Setup matrices
  //

  if (_nl_sys_nums.size() > 0)
  {
    // Need full Jacobian in order to actually compute cross-system coupling
    _fe_problem.setCoupling(Moose::COUPLING_FULL);
    for (const auto i : index_range(_nl_sys_nums))
      for (const auto j : index_range(_nl_sys_nums))
      {
        if (i == j)
          continue;
        const TagID tag =
            _fe_problem.addMatrixTag("NPC_J_" + std::to_string(i) + "_" + std::to_string(j));
        _off_diag_mats[{i, j}] = {nullptr, tag};
      }
  }

  //
  // Setup convergence objects
  //

  if (isParamValid("convergence_names"))
  {
    const auto & conv_names = getParam<std::vector<ConvergenceName>>("convergence_names");
    if (conv_names.size() != nl_sys_names.size())
      paramError("convergence_names",
                 "Must have the same length as 'nonlinear_system_names' (",
                 nl_sys_names.size(),
                 ")");
    for (const auto i : index_range(nl_sys_names))
      _fe_problem.setNonlinearConvergence(nl_sys_names[i], conv_names[i]);
  }
}

NewtonSNESExecutor::~NewtonSNESExecutor()
{
  if (_r_plain)
    PetscCallAbort(this->comm().get(), VecDestroy(&_r_plain));
  if (_r_base)
    PetscCallAbort(this->comm().get(), VecDestroy(&_r_base));
  if (_fullJ_ksp)
    PetscCallAbort(this->comm().get(), KSPDestroy(&_fullJ_ksp));
  if (_pc_ksp0)
    PetscCallAbort(this->comm().get(), KSPDestroy(&_pc_ksp0));
  if (_pc_ksp1)
    PetscCallAbort(this->comm().get(), KSPDestroy(&_pc_ksp1));
  for (Vec * v : {&_cg_p, &_cg_r, &_cg_y, &_cg_d, &_cg_Ad, &_cg_Pp, &_cg_Pd})
    if (*v)
      PetscCallAbort(this->comm().get(), VecDestroy(v));
  for (Vec * v : {&_active_mask, &_pdas_prev_dstep, &_x_presweep, &_sweep_dir, &_r_stash})
    if (*v)
      PetscCallAbort(this->comm().get(), VecDestroy(v));
  if (_jac_shell)
    PetscCallAbort(this->comm().get(), MatDestroy(&_jac_shell));
  if (_mat_nest)
    PetscCallAbort(this->comm().get(), MatDestroy(&_mat_nest));
}

void
NewtonSNESExecutor::setupSNES()
{
  const auto n_sys = _nl_sys_nums.size();

  if (n_sys == 1)
    // We don't need to create a new SNES with Nest data structures. We'll just be leveraging the
    // libMesh solver's SNES
    return;

  // Build VecNest for solution and residual.
  // _vec_sol uses independent (non-aliasing) sub-vecs so the outer Newton iterate is not
  // corrupted when the NPC sub-solve updates the libmesh solution vectors.
  std::vector<Vec> ind_sol_vecs(n_sys);
  std::vector<Vec> rhs_vecs(n_sys);
  for (unsigned int i = 0; i < n_sys; ++i)
  {
    auto & sys_i = _fe_problem.getNonlinearSystem(i);
    Vec libmesh_sol =
        cast_ptr<libMesh::PetscVector<libMesh::Number> *>(sys_i.system().solution.get())->vec();
    LibmeshPetscCallA(this->comm().get(), VecDuplicate(libmesh_sol, &ind_sol_vecs[i]));
    LibmeshPetscCallA(this->comm().get(), VecCopy(libmesh_sol, ind_sol_vecs[i]));
    rhs_vecs[i] = cast_ptr<libMesh::PetscVector<libMesh::Number> *>(&sys_i.RHS())->vec();
  }
  LibmeshPetscCallA(
      this->comm().get(),
      VecCreateNest(this->comm().get(), n_sys, nullptr, ind_sol_vecs.data(), &_vec_sol));
  for (unsigned int i = 0; i < n_sys; ++i)
    LibmeshPetscCallA(this->comm().get(), VecDestroy(&ind_sol_vecs[i]));
  LibmeshPetscCallA(this->comm().get(),
                    VecCreateNest(this->comm().get(), n_sys, nullptr, rhs_vecs.data(), &_vec_func));

  allocateOffDiagMats();
  buildMatNest();

  LibmeshPetscCallA(this->comm().get(), SNESCreate(this->comm().get(), &_snes));
  // Make sure to set type early because otherwise SNESSetFromOptions at the end of this function
  // will call SNESSetType which will forward to SNESCreate_NEWTONLS which will overwrite our NPC
  // default side choice
  LibmeshPetscCallA(this->comm().get(), SNESSetType(_snes, SNESNEWTONLS));
  LibmeshPetscCallA(this->comm().get(),
                    SNESSetFunction(_snes, _vec_func, outerResidualCallback, this));

  PetscInt M, N, m, n;
  LibmeshPetscCallA(this->comm().get(), MatGetSize(_mat_nest, &M, &N));
  LibmeshPetscCallA(this->comm().get(), MatGetLocalSize(_mat_nest, &m, &n));
  LibmeshPetscCallA(this->comm().get(),
                    MatCreateShell(this->comm().get(), m, n, M, N, this, &_jac_shell));
  LibmeshPetscCallA(this->comm().get(),
                    MatShellSetOperation(_jac_shell, MATOP_MULT, (PetscErrorCodeFn *)shellMatMult));

  // For the Steihaug trust-region variants -- (A) no-NPC and (B) NPC-every-iteration -- the shell
  // line search computes its own step from the coupled Hessian _mat_nest; the outer Krylov direction
  // is discarded. So the outer SNES Jacobian operator is the raw _mat_nest (never the applyBA shell),
  // the outer function is the PLAIN residual (no automatic NPC sweep -- (B) applies the NPC itself in
  // the line search), and the outer KSP is made a trivial PREONLY solve below. The App. B/C line
  // search paths keep the SPIN wiring (shell operator A = M^{-1}A, preconditioned fixed-point residual).
  const bool tr_coupled = _use_trust_region && _tr_steihaug;

  LibmeshPetscCallA(
      this->comm().get(),
      SNESSetJacobian(
          _snes, tr_coupled ? _mat_nest : _jac_shell, _mat_nest, outerJacobianCallback, this));

  // Outer KSP needs no linear preconditioner: all preconditioning is at the NL level.
  KSP ksp;
  PC pc;
  LibmeshPetscCallA(this->comm().get(), SNESGetKSP(_snes, &ksp));
  LibmeshPetscCallA(this->comm().get(), KSPGetPC(ksp, &pc));
  LibmeshPetscCallA(this->comm().get(), PCSetType(pc, PCNONE));

  if (!tr_coupled)
  {
    // Drive ||x - NPC(x)|| to zero (ASPIN-style fixed-point convergence).
    LibmeshPetscCallA(this->comm().get(), SNESSetNPCSide(_snes, PC_LEFT));
    LibmeshPetscCallA(this->comm().get(),
                      SNESSetFunctionType(_snes, SNES_FUNCTION_PRECONDITIONED));
  }

  // Set options prefix and set from options
  LibmeshPetscCallA(this->comm().get(), SNESSetOptionsPrefix(_snes, (this->name() + "_").c_str()));
  LibmeshPetscCallA(this->comm().get(), SNESSetFromOptions(_snes));

  // For the coupled trust-region variants make the outer Krylov a trivial PREONLY solve -- the shell
  // line search discards the outer direction and steps via the Steihaug TRS. Done AFTER
  // SNESSetFromOptions so it wins over any -<name>_ksp_type option (e.g. gmres) in the input.
  if (tr_coupled)
    LibmeshPetscCallA(this->comm().get(), KSPSetType(ksp, KSPPREONLY));

  // Install a custom shell line search for the preconditioned SPIN residual. Done AFTER
  // SNESSetFromOptions so it wins over any -<name>_snes_linesearch_type option. The stock line
  // search evaluates trial points with the plain residual (SNESComputeFunction), which is
  // inconsistent with the preconditioned residual F_SPIN = x - NPC(x) that this outer SNES
  // actually drives to zero, and consequently cannot backtrack once the fields couple strongly.
  // spinLineSearch uses SNESApplyNPC at every trial for a consistent merit (see its definition).
  SNESLineSearch ls;
  LibmeshPetscCallA(this->comm().get(), SNESGetLineSearch(_snes, &ls));
  LibmeshPetscCallA(this->comm().get(), SNESLineSearchSetType(ls, SNESLINESEARCHSHELL));
  LibmeshPetscCallA(this->comm().get(),
                    SNESLineSearchShellSetApply(ls, &NewtonSNESExecutor::spinLineSearch, this));

  // Dedicated plain-residual VecNest R: grad Psi (energy merit) or the coupled nonlinear residual
  // (residual merit) / its reduced form (PDAS). Read by computePlainResidual, computeActiveSet,
  // applyActiveSetElimination, steihaugTRS, and the energy line search -- allocate for ANY path that
  // needs it, not just the energy merit.
  const bool tr_steihaug_path = _use_trust_region && _tr_steihaug;
  if (_has_energy_pp || tr_steihaug_path || _bounds)
    LibmeshPetscCallA(this->comm().get(), VecDuplicate(_vec_func, &_r_plain));

  // Full-Jacobian KSP for the App. C inexact-Newton fallback direction (ENERGY merit only; the
  // Steihaug TR path never uses it). Fieldsplit preconditioner; tunable via -<name>_fullj_* options.
  if (_has_energy_pp)
  {
    PC pc_fullj;
    LibmeshPetscCallA(this->comm().get(), KSPCreate(this->comm().get(), &_fullJ_ksp));
    LibmeshPetscCallA(this->comm().get(),
                      KSPSetOptionsPrefix(_fullJ_ksp, (this->name() + "_fullj_").c_str()));
    LibmeshPetscCallA(this->comm().get(), KSPSetType(_fullJ_ksp, KSPGMRES));
    LibmeshPetscCallA(this->comm().get(), KSPGetPC(_fullJ_ksp, &pc_fullj));
    LibmeshPetscCallA(this->comm().get(), PCSetType(pc_fullj, PCFIELDSPLIT));
    LibmeshPetscCallA(this->comm().get(),
                      KSPSetTolerances(_fullJ_ksp, 1e-3, PETSC_DEFAULT, PETSC_DEFAULT, 200));
    LibmeshPetscCallA(this->comm().get(), KSPSetFromOptions(_fullJ_ksp));
  }

  // Outer convergence on the COUPLED residual ||R|| (grad Psi for the energy merit; the coupled
  // nonlinear residual for the residual merit) with the HWW active-set-unchanged gate -- NOT the
  // preconditioned residual ||x - NPC(x)||. The paper's Algorithm 3 while-condition is
  // ||F(U,C)|| >= eps_rel ||F(U0,C0)|| (F = the coupled/monolithic residual, eq. 14 = grad Psi), and
  // Table A.7 documents -snes_atol/-snes_rtol as "(coupled residual)". Every trust-region path
  // converges on ||R|| cached in _r_norm/_r0_norm by spinLineSearch (for the energy merit the NPC
  // pre-minimizes block-wise so ||x - NPC(x)|| plateaus while ||R|| -> 0).
  if (_has_energy_pp || tr_steihaug_path)
    LibmeshPetscCallA(this->comm().get(),
                      SNESSetConvergenceTest(_snes, outerConvergenceTest, this, nullptr));

  // Steihaug-Toint TRS machinery (energy OR residual merit): block field-split preconditioner on the
  // two SPD diagonal blocks A00 (disp), A11 (pf) + CG scratch VecNests. Block PC default LU/MUMPS;
  // tunable via -<name>_trpc0_* / -<name>_trpc1_*. Operators (re)bound to the current blocks each solve.
  if (tr_steihaug_path)
  {
    Mat A00, A11;
    LibmeshPetscCallA(this->comm().get(), MatNestGetSubMat(_mat_nest, 0, 0, &A00));
    LibmeshPetscCallA(this->comm().get(), MatNestGetSubMat(_mat_nest, 1, 1, &A11));
    const std::array<std::pair<KSP *, Mat>, 2> blocks{{{&_pc_ksp0, A00}, {&_pc_ksp1, A11}}};
    for (std::size_t b = 0; b < blocks.size(); ++b)
    {
      KSP * ksp = blocks[b].first;
      PC pc_b;
      LibmeshPetscCallA(this->comm().get(), KSPCreate(this->comm().get(), ksp));
      LibmeshPetscCallA(
          this->comm().get(),
          KSPSetOptionsPrefix(*ksp, (this->name() + "_trpc" + std::to_string(b) + "_").c_str()));
      LibmeshPetscCallA(this->comm().get(), KSPSetType(*ksp, KSPPREONLY));
      LibmeshPetscCallA(this->comm().get(), KSPGetPC(*ksp, &pc_b));
      LibmeshPetscCallA(this->comm().get(), PCSetType(pc_b, PCLU));
      // Parallel-capable direct solver (petsc's built-in LU is serial-only). Use MUMPS, NOT
      // SuperLU_DIST: SuperLU_DIST frees its process-grid MPI communicator twice at teardown
      // (once in MatDestroy_SuperLU_DIST -> superlu_gridexit, again via its comm-keyval delete
      // callback at PetscFinalize), which double-frees the OpenMPI communicator and corrupts the
      // heap ("corrupted double-linked list") at exit. MUMPS manages its comm cleanly.
      LibmeshPetscCallA(this->comm().get(), PCFactorSetMatSolverType(pc_b, MATSOLVERMUMPS));
      LibmeshPetscCallA(this->comm().get(),
                        KSPSetOperators(*ksp, blocks[b].second, blocks[b].second));
      LibmeshPetscCallA(this->comm().get(), KSPSetFromOptions(*ksp));
    }
    // CG scratch VecNests (same block layout as the residual/solution nest).
    for (Vec * v : {&_cg_p, &_cg_r, &_cg_y, &_cg_d, &_cg_Ad, &_cg_Pp, &_cg_Pd})
      LibmeshPetscCallA(this->comm().get(), VecDuplicate(_vec_func, v));
    // Residual-merit TR only: preserve R(X) across the ratio-test retry loop (see spinLineSearch).
    if (!_has_energy_pp)
      LibmeshPetscCallA(this->comm().get(), VecDuplicate(_vec_func, &_r_base));
  }

  _snes_setup_done = true;
}

Executor::Result
NewtonSNESExecutor::run()
{
  auto & result = newResult();

  static const std::string solve_converged_msg = "Solve converged";
  static const std::string solve_didnt_converge_msg = "Solve failed to converge";

  if (_nl_sys_nums.size() == 1)
  {
    const auto nl_sys_num = _nl_sys_nums[0];

    // Bound-constrained NPC block: solve the reduced subsystem (frozen DOFs held at d_old) instead of
    // the stock nonlinear solve. Armed by the parent coupled executor before an NPC sweep.
    if (_reduced_solve_armed)
    {
      reducedNewtonSolve();
      result.pass(solve_converged_msg);
      return result;
    }

    // Wire the nonlinear preconditioner if we have it
    if (_npc_executor)
      LibmeshPetscCallA(this->comm().get(),
                        SNESSetNPC(_fe_problem.getNonlinearSystem(nl_sys_num).getSNES(),
                                   _npc_executor->getSNES()));

    _fe_problem.solve(nl_sys_num);

    PetscInt iter;
    LibmeshPetscCallA(
        this->comm().get(),
        SNESGetIterationNumber(_fe_problem.getNonlinearSystem(nl_sys_num).getSNES(), &iter));
    const auto & conv_name = _fe_problem.getNonlinearConvergenceNames()[nl_sys_num];
    const auto status = _fe_problem.getConvergence(conv_name).checkConvergence(iter);
    if (status == Convergence::MooseConvergenceStatus::CONVERGED)
      result.pass(solve_converged_msg);
    else
      result.fail(solve_didnt_converge_msg);
    return result;
  }

  //
  // Multi-system: combined outer Newton with VecNest/MatNest.
  //

  if (!_npc_executor)
    mooseError("NewtonSNESExecutor: multiple nonlinear systems currently require "
               "'nl_preconditioning' to be set.");

  // Ensure the PETSc options database is populated before SNESSetFromOptions() runs inside
  // setupSNES(). Without this, petscSetOptions() called from the first FEProblemBase::solve()
  // during the NPC sub-solve would clear and re-insert options with used=false, causing PETSc to
  // report options like -outer_snes_monitor as unused at finalization even though they are active.
  _fe_problem.insertPetscOptionsIfNeeded();

  // SNES setup also calls buildMatNest()
  if (!_snes_setup_done)
    setupSNES();
  else
    buildMatNest();

  // Build the cached lumped-mass inverse once for the bound-constrained (PDAS) layer.
  if (_bounds && !_mass_built)
    buildLumpedMassInverse();

  // Attach the nonlinear preconditioner for the App. B/C line-search paths. NOT for the coupled
  // trust-region variants (A) and (B): (A) uses no NPC at all; (B) applies the NPC sweep EXPLICITLY
  // in the line search via SNESSolve on the NMSM shell (see spinLineSearch). Attaching it would make
  // PETSc's NEWTONLS auto-apply the NPC every iteration, doubling the sweep (measured 2x).
  const bool tr_coupled = _use_trust_region && _tr_steihaug;
  if (!tr_coupled)
    LibmeshPetscCallA(this->comm().get(), SNESSetNPC(_snes, _npc_executor->getSNES()));

  // For an energy-based line search, the outer iterate must satisfy the current step's Dirichlet
  // boundary conditions from iteration 0 -- otherwise Psi(X_0) omits the boundary-condition energy
  // (e.g. Psi(0)=0 at an undeformed start) and no step can reduce the (non-negative) energy below
  // it, so the line search fails immediately. Preset the nodal BCs into each system's solution and
  // re-seed the outer iterate from it. (The residual-norm merit does not need this, since
  // ||x - NPC(x)|| -> 0 at the solution regardless of how the BC is reached.)
  if (_has_energy_pp)
    for (const auto i : index_range(_nl_sys_nums))
    {
      auto & sys_i = _fe_problem.getNonlinearSystemBase(static_cast<unsigned int>(i));
      sys_i.setInitialSolution();
      Vec libmesh_sol =
          cast_ptr<libMesh::PetscVector<libMesh::Number> *>(sys_i.system().solution.get())->vec();
      Vec sub;
      LibmeshPetscCallA(this->comm().get(), VecNestGetSubVec(_vec_sol, i, &sub));
      LibmeshPetscCallA(this->comm().get(), VecCopy(libmesh_sol, sub));
    }

  // Reset the ||grad Psi|| convergence baseline for this solve (energy-merit path).
  _r_norm = -1.0;
  _r0_norm = -1.0;
  // Reset the trust-region radius; it is initialized to ||Y|| on the first line search of the solve.
  _tr_radius = -1.0;
  // Reset the PDAS state for this solve: the previous d-increment is zero (initial guess d = d_old),
  // and the active set starts unknown so the first iteration is never spuriously "converged".
  if (_bounds)
  {
    _active_d_prev.clear();
    _active_set_changed = true;
    if (_pdas_prev_dstep)
      LibmeshPetscCallA(this->comm().get(), VecSet(_pdas_prev_dstep, 0.0));
  }

  LibmeshPetscCallA(this->comm().get(), SNESSolve(_snes, nullptr, _vec_sol));

  SNESConvergedReason reason;
  LibmeshPetscCallA(this->comm().get(), SNESGetConvergedReason(_snes, &reason));
  if (reason > 0)
    result.pass(solve_converged_msg);
  else
    result.fail(solve_didnt_converge_msg);
  return result;
}

void
NewtonSNESExecutor::buildLumpedMassInverse()
{
  // Resolve the bounded variable's system/var numbers (variables exist by the first solve, unlike at
  // executor-construction time). getVariable searches all systems, so it works with the field split.
  const auto & var = _fe_problem.getVariable(/*tid=*/0,
                                             _bounded_var_name,
                                             Moose::VarKindType::VAR_SOLVER,
                                             Moose::VarFieldType::VAR_FIELD_ANY);
  _bounded_sys_num = var.sys().number();
  _bounded_var_num = var.number();
  _bounded_sys_local = libMesh::invalid_uint;
  for (const auto i : index_range(_nl_sys_nums))
    if (_nl_sys_nums[i] == _bounded_sys_num)
      _bounded_sys_local = static_cast<unsigned int>(i);
  if (_bounded_sys_local == libMesh::invalid_uint)
    mooseError("NewtonSNESExecutor: bounded_variable '",
               _bounded_var_name,
               "' is not in any of this executor's nonlinear systems.");

  // NEPIN[FULL]/[SUBDOMAIN] additionally eliminate the displacement front, so resolve the (single)
  // non-bounded coupled system as the "disp" system. Both assume a two-system disp+pf coupling.
  if (_nepin_full)
  {
    for (const auto i : index_range(_nl_sys_nums))
      if (_nl_sys_nums[i] != _bounded_sys_num)
      {
        _disp_sys_local = static_cast<unsigned int>(i);
        _disp_sys_num = _nl_sys_nums[i];
      }
    if (_disp_sys_local == libMesh::invalid_uint || _nl_sys_nums.size() != 2)
      mooseError("NewtonSNESExecutor: npc=mspin_fne expects exactly two coupled systems "
                 "(displacement + bounded phase field); found ",
                 _nl_sys_nums.size(),
                 ".");
  }

  _mass_tag = _fe_problem.getMatrixTagID(_mass_tag_name);
  auto & dsys = _fe_problem.getNonlinearSystemBase(_bounded_sys_num);
  if (!dsys.hasMatrix(_mass_tag))
    mooseError("NewtonSNESExecutor: no mass matrix on tag '",
               _mass_tag_name,
               "'. Add a MassMatrix kernel on '",
               _bounded_var_name,
               "' with matrix_tags='",
               _mass_tag_name,
               "' and [Problem] extra_tag_matrices='",
               _mass_tag_name,
               "'.");

  libMesh::SparseMatrix<libMesh::Number> & M = dsys.getMatrix(_mass_tag);

  // Assemble the mass matrix ONCE (geometry-only; constant across the whole run).
  _fe_problem.setCurrentNonlinearSystem(_bounded_sys_num);
  _fe_problem.computeJacobianTag(*dsys.system().current_local_solution, M, _mass_tag);

  // Lump by row sum (= integral of N_i for a partition-of-unity basis; the mass system holds only
  // the bounded variable, so every diagonal is nonzero), then invert to get the B^{-1} diagonal.
  _binv = dsys.solution().zero_clone();
  auto ones = dsys.solution().zero_clone();
  *ones = 1.0;
  ones->close();
  M.vector_mult(*_binv, *ones); // _binv = M * 1 = row sums = lumped mass
  _binv->close();

  const Real mmin = _binv->min();
  const Real mmax = _binv->max();

  _binv->reciprocal(); // -> B^{-1} diagonal
  _binv->close();

  _mass_built = true;

  _console << "[PDAS] lumped mass on '" << _bounded_var_name << "' (nl sys " << _bounded_sys_num
           << ", local index " << _bounded_sys_local << "): lumped-mass min=" << mmin
           << " max=" << mmax << ", c=" << _pdas_c << std::endl;
}

void
NewtonSNESExecutor::computeActiveSet(Vec X)
{
  const auto comm = this->comm().get();
  auto & dsys = _fe_problem.getNonlinearSystemBase(_bounded_sys_num);
  const libMesh::NumericVector<libMesh::Number> & d_old = dsys.solutionOld();

  // Bounded-variable block of the coupled residual (R = grad Psi) and current iterate.
  Vec Rd, Xd;
  LibmeshPetscCallA(comm, VecNestGetSubVec(_r_plain, _bounded_sys_local, &Rd));
  LibmeshPetscCallA(comm, VecNestGetSubVec(X, _bounded_sys_local, &Xd));

  // Lazy-allocate the mask and previous-step vectors (bounded-system layout, matching Rd).
  if (!_active_mask)
    LibmeshPetscCallA(comm, VecDuplicate(Rd, &_active_mask));
  if (!_pdas_prev_dstep)
  {
    LibmeshPetscCallA(comm, VecDuplicate(Rd, &_pdas_prev_dstep));
    LibmeshPetscCallA(comm, VecSet(_pdas_prev_dstep, 0.0));
  }

  PetscInt lo, hi;
  LibmeshPetscCallA(comm, VecGetOwnershipRange(Rd, &lo, &hi));

  const PetscScalar *r_arr, *x_arr, *ds_arr;
  PetscScalar * mask_arr;
  LibmeshPetscCallA(comm, VecGetArrayRead(Rd, &r_arr));
  LibmeshPetscCallA(comm, VecGetArrayRead(Xd, &x_arr));
  LibmeshPetscCallA(comm, VecGetArrayRead(_pdas_prev_dstep, &ds_arr));
  LibmeshPetscCallA(comm, VecGetArray(_active_mask, &mask_arr));

  _active_d_dofs.clear();
  _active_val.clear();
  unsigned int n_lower = 0, n_upper = 0;
  const PetscInt n = hi - lo;
  // Projected-gradient (min-map) active set, following PETSc's vinewtonrsls (SNESVIGetActiveSetIS,
  // vi.c): a bounded DOF joins the active set from its current POSITION and multiplier SIGN, with NO
  // predictor delta d^{k-1} and NO complementarity constant c (the HWW predictor drove a period-2
  // active-set limit cycle). The multiplier lambda = (B^{-1})_ii R_i is lumped-mass scaled, so the
  // dead-band _pdas_lambda_tol is a mesh-independent pointwise driving force: DOFs whose multiplier is
  // numerically indifferent (~roundoff, e.g. the unstressed pre-damage region) stay inactive rather
  // than oscillating in/out of the set.
  //   lower-active: d <= d_old + ztol  AND  lambda >  +lambda_tol  -> freeze at d_old
  //   upper-active: d >= 1     - ztol  AND  lambda <  -lambda_tol  -> freeze at 1
  const PetscReal ztol = 1e-8;
  (void)ds_arr; // no predictor in the projected-gradient criterion (ds_arr consumed by Restore below)
  for (PetscInt j = 0; j < n; ++j)
  {
    const libMesh::dof_id_type gdof = static_cast<libMesh::dof_id_type>(lo + j);
    const PetscReal di = PetscRealPart(x_arr[j]);
    const PetscReal lambda = (*_binv)(gdof)*PetscRealPart(r_arr[j]); // pointwise multiplier estimate
    const PetscReal dlo = d_old(gdof);
    const PetscReal dup = 1.0; // physical upper bound on damage
    if (di <= dlo + ztol && lambda > _pdas_lambda_tol)
    {
      _active_d_dofs.push_back(lo + j);
      _active_val.push_back(dlo);
      mask_arr[j] = 0.0;
      ++n_lower;
    }
    else if (di >= dup - ztol && lambda < -_pdas_lambda_tol)
    {
      _active_d_dofs.push_back(lo + j);
      _active_val.push_back(dup);
      mask_arr[j] = 0.0;
      ++n_upper;
    }
    else
      mask_arr[j] = 1.0;
  }

  LibmeshPetscCallA(comm, VecRestoreArrayRead(Rd, &r_arr));
  LibmeshPetscCallA(comm, VecRestoreArrayRead(Xd, &x_arr));
  LibmeshPetscCallA(comm, VecRestoreArrayRead(_pdas_prev_dstep, &ds_arr));
  LibmeshPetscCallA(comm, VecRestoreArray(_active_mask, &mask_arr));

  // Global change flag: the active set has "changed" if it changed on ANY rank (so the HWW
  // convergence decision in outerConvergenceTest is identical on all ranks, as PETSc requires).
  unsigned int changed = (_active_d_dofs != _active_d_prev) ? 1u : 0u;
  this->comm().max(changed);
  _active_set_changed = (changed != 0u);
  _active_d_prev = _active_d_dofs;

  if (_verbose)
  {
    this->comm().sum(n_lower);
    this->comm().sum(n_upper);
    _console << "    [PDAS] active=" << (n_lower + n_upper) << " / " << dsys.system().n_dofs()
             << " (lower=" << n_lower << " upper=" << n_upper << ")  changed=" << _active_set_changed
             << std::endl;
  }
}

void
NewtonSNESExecutor::computeHardSet(Vec X)
{
  // NEPIN nonlinear-elimination set. The "hard" (strongly nonlinear) DOFs are the phase-field process
  // zone -- the moving damage band nepin_lo < d < 1 - nepin_hi where g(d) and the u-d coupling vary
  // steeply. Those get ELIMINATED (left free in the reduced pf sub-solve, box [d_old, 1]); every other
  // damage DOF (far-field d~0, broken d~1, irreversibly-frozen front) is "easy" and pinned at its
  // current value. armReducedSolve pins the DOFs it is handed and frees the rest, so we store the EASY
  // (complement) set here. For NEPIN[PARTIAL] displacement is not touched -- it stays "easy" and the
  // outer TR moves it. For NEPIN[FULL] we ALSO build the displacement front set below (both fields move
  // together on the band), which is the variant that actually cuts outer iterations.
  const auto comm = this->comm().get();
  Vec Xd;
  LibmeshPetscCallA(comm, VecNestGetSubVec(X, _bounded_sys_local, &Xd));

  PetscInt lo, hi;
  LibmeshPetscCallA(comm, VecGetOwnershipRange(Xd, &lo, &hi));
  const PetscScalar * x_arr;
  LibmeshPetscCallA(comm, VecGetArrayRead(Xd, &x_arr));

  _nepin_frozen.clear();
  _nepin_frozen_val.clear();
  unsigned int n_bad = 0;
  const PetscInt n = hi - lo;
  const PetscReal d_hi = 1.0 - _nepin_hi;
  for (PetscInt j = 0; j < n; ++j)
  {
    const PetscReal di = PetscRealPart(x_arr[j]);
    if (di > _nepin_lo && di < d_hi)
      ++n_bad; // process zone -> eliminate (free in the reduced pf sub-solve)
    else
    {
      _nepin_frozen.push_back(lo + j); // easy -> pin at current value
      _nepin_frozen_val.push_back(di);
    }
  }
  LibmeshPetscCallA(comm, VecRestoreArrayRead(Xd, &x_arr));

  // NEPIN[FULL]: the displacement front. The hard disp DOFs are where the coupled residual (force
  // imbalance from the changing g(d)) is largest -- |R_u,i| > nepin_u_rtol * max|R_u| -- which is the
  // moving crack front. Freeze the far-field disp (its complement) at current so the reduced disp solve
  // moves only the band; the frozen ring supplies the Dirichlet anchor that removes rigid-body modes.
  unsigned int n_bad_u = 0;
  _nepin_frozen_u.clear();
  _nepin_frozen_u_val.clear();
  if (_nepin_full)
  {
    Vec Ru, Xu;
    LibmeshPetscCallA(comm, VecNestGetSubVec(_r_plain, _disp_sys_local, &Ru));
    LibmeshPetscCallA(comm, VecNestGetSubVec(X, _disp_sys_local, &Xu));
    PetscReal rmax;
    LibmeshPetscCallA(comm, VecNorm(Ru, NORM_INFINITY, &rmax));
    const PetscReal thr = _nepin_u_rtol * rmax;

    PetscInt ulo, uhi;
    LibmeshPetscCallA(comm, VecGetOwnershipRange(Ru, &ulo, &uhi));
    const PetscScalar *ru_arr, *xu_arr;
    LibmeshPetscCallA(comm, VecGetArrayRead(Ru, &ru_arr));
    LibmeshPetscCallA(comm, VecGetArrayRead(Xu, &xu_arr));
    const PetscInt nu = uhi - ulo;
    for (PetscInt j = 0; j < nu; ++j)
    {
      const PetscReal ri = PetscAbsReal(PetscRealPart(ru_arr[j]));
      if (rmax > 0.0 && ri > thr)
        ++n_bad_u; // front DOF -> eliminate (free, unbounded box)
      else
      {
        _nepin_frozen_u.push_back(ulo + j); // far-field -> pin at current value
        _nepin_frozen_u_val.push_back(PetscRealPart(xu_arr[j]));
      }
    }
    LibmeshPetscCallA(comm, VecRestoreArrayRead(Ru, &ru_arr));
    LibmeshPetscCallA(comm, VecRestoreArrayRead(Xu, &xu_arr));
  }

  if (_verbose)
  {
    this->comm().sum(n_bad);
    auto & dsys = _fe_problem.getNonlinearSystemBase(_bounded_sys_num);
    _console << "    [NEPIN] hard d(eliminated)=" << n_bad << " / " << dsys.system().n_dofs()
             << " (band " << _nepin_lo << " < d < " << d_hi << ")";
    if (_nepin_full)
    {
      this->comm().sum(n_bad_u);
      auto & usys = _fe_problem.getNonlinearSystemBase(_disp_sys_num);
      _console << "  |  hard u(eliminated)=" << n_bad_u << " / " << usys.system().n_dofs()
               << " (|R_u|>" << _nepin_u_rtol << "*max)";
    }
    _console << std::endl;
  }
}

PetscErrorCode
NewtonSNESExecutor::applyActiveSetElimination()
{
  PetscFunctionBeginUser;
  const PetscInt i_d = static_cast<PetscInt>(_bounded_sys_local);
  const PetscInt n_active = static_cast<PetscInt>(_active_d_dofs.size());
  const PetscInt * idx = _active_d_dofs.empty() ? nullptr : _active_d_dofs.data();

  // (d,d) diagonal block: zero active rows AND columns, unit diagonal -> delta d = 0 on the active
  // set (dynamic Dirichlet elimination). Modifying this block (the pf system matrix) bumps its state
  // so the block preconditioner (_pc_ksp1 LU / block-SGS) refactors the reduced block downstream.
  Mat Add;
  PetscCall(MatNestGetSubMat(_mat_nest, i_d, i_d, &Add));
  PetscCall(MatZeroRowsColumns(Add, n_active, idx, 1.0, nullptr, nullptr));

  // Off-diagonal coupling blocks (restore symmetry, HWW step 4).
  for (const auto j : index_range(_nl_sys_nums))
  {
    if (static_cast<PetscInt>(j) == i_d)
      continue;
    // (d,u): zero the active d ROWS -- the frozen d-equations no longer couple to u.
    Mat Adu;
    PetscCall(MatNestGetSubMat(_mat_nest, i_d, static_cast<PetscInt>(j), &Adu));
    PetscCall(MatZeroRows(Adu, n_active, idx, 0.0, nullptr, nullptr));
    // (u,d): zero the active d COLUMNS -- column-scale by the inactive mask (no MatZeroColumns in
    // PETSc); keeps the coupled operator symmetric for the Steihaug CG.
    Mat Aud;
    PetscCall(MatNestGetSubMat(_mat_nest, static_cast<PetscInt>(j), i_d, &Aud));
    PetscCall(MatDiagonalScale(Aud, nullptr, _active_mask));
  }

  // Zero the active residual entries so the reduced gradient drives delta d = 0 there.
  Vec Rd;
  PetscCall(VecNestGetSubVec(_r_plain, i_d, &Rd));
  PetscCall(VecPointwiseMult(Rd, Rd, _active_mask));

  PetscFunctionReturn(PETSC_SUCCESS);
}

void
NewtonSNESExecutor::clampActiveToBounds(Vec X)
{
  if (_active_d_dofs.empty())
    return;
  const auto comm = this->comm().get();
  Vec Xd;
  LibmeshPetscCallA(comm, VecNestGetSubVec(X, _bounded_sys_local, &Xd));
  PetscInt lo, hi;
  LibmeshPetscCallA(comm, VecGetOwnershipRange(Xd, &lo, &hi));
  PetscScalar * xa;
  LibmeshPetscCallA(comm, VecGetArray(Xd, &xa));
  for (std::size_t i = 0; i < _active_d_dofs.size(); ++i)
  {
    const PetscInt g = _active_d_dofs[i];
    if (g >= lo && g < hi)
      xa[g - lo] = _active_val[i];
  }
  LibmeshPetscCallA(comm, VecRestoreArray(Xd, &xa));
}

void
NewtonSNESExecutor::armReducedSolve(const std::vector<PetscInt> & frozen,
                                    const std::vector<PetscReal> & vals,
                                    bool bounded_box,
                                    bool restrict_assembly)
{
  _reduced_solve_armed = true;
  _reduced_bounded_box = bounded_box;
  _reduced_restrict_assembly = restrict_assembly;
  _frozen_dofs = frozen;
  _frozen_vals = vals;
}

void
NewtonSNESExecutor::disarmReducedSolve()
{
  _reduced_solve_armed = false;
  _reduced_restrict_assembly = false;
  _frozen_dofs.clear();
  _frozen_vals.clear();
}

void
NewtonSNESExecutor::reducedNewtonSolve()
{
  const auto sys_num = _nl_sys_nums[0];
  auto & sys = _fe_problem.getNonlinearSystemBase(sys_num);
  auto & lm = _fe_problem.getNonlinearSystem(sys_num).sys();

  // Install our bounds callback on this system's nonlinear solver once (replaces MOOSE's default
  // compute_bounds for the VI sub-solve). It copies lower_bound/upper_bound (filled below) into the
  // SNES's VI bound vectors -- MOOSE/libMesh only, no raw PETSc, no shared Vec ownership.
  if (lm.nonlinear_solver->bounds != pdasCopyBounds)
    lm.nonlinear_solver->bounds = pdasCopyBounds;

  // Fill the reduced-solve VI bounds into the system's MOOSE-owned bound vectors. FREE DOFs get the
  // field-appropriate box: the bounded (pf) block uses the damage irreversibility box [d_old, 1]; the
  // displacement block (NEPIN[FULL]) is unbounded [-inf, +inf], so its free DOFs solve as plain Newton.
  // FROZEN DOFs are pinned (lower = upper = freeze value) -> dynamic Dirichlet: the pf block reuses the
  // outer active/inactive partition exactly; the disp block anchors the free front patch against RBM.
  auto & lb = sys.getVector("lower_bound");
  auto & ub = sys.getVector("upper_bound");
  if (_reduced_bounded_box)
  {
    const libMesh::NumericVector<libMesh::Number> & d_old = sys.solutionOld();
    for (auto i = lb.first_local_index(); i < lb.last_local_index(); ++i)
    {
      lb.set(i, d_old(i));
      ub.set(i, 1.0);
    }
  }
  else
    for (auto i = lb.first_local_index(); i < lb.last_local_index(); ++i)
    {
      lb.set(i, PETSC_NINFINITY);
      ub.set(i, PETSC_INFINITY);
    }
  for (std::size_t k = 0; k < _frozen_dofs.size(); ++k)
  {
    const auto g = static_cast<libMesh::dof_id_type>(_frozen_dofs[k]);
    lb.set(g, _frozen_vals[k]);
    ub.set(g, _frozen_vals[k]);
  }
  lb.close();
  ub.close();

  // NEPIN[FULL] band-restricted assembly (nepin_restrict_assembly): restrict the FE assembly to the
  // locally-owned elements incident to this block's FREE set (owned DOFs not frozen), so the
  // residual/Jacobian cost scales with the crack-front
  // band, not the whole mesh. vinewtonrsls already solves only the free (inactive) reduced space, so the
  // assembly is the sole remaining full-domain cost -- this removes it. All-gather the free set (the mesh
  // is replicated) so a free DOF's incident elements owned by another rank are assembled there and the
  // parallel assembly sums the contributions -> free rows are complete; the frozen ring is VI-pinned.
  if (_reduced_restrict_assembly)
  {
    std::set<PetscInt> frozen(_frozen_dofs.begin(), _frozen_dofs.end());
    const libMesh::DofMap & dm = lm.get_dof_map();
    std::vector<PetscInt> gfree;
    for (auto g = dm.first_dof(); g < dm.end_dof(); ++g)
      if (!frozen.count(static_cast<PetscInt>(g)))
        gfree.push_back(static_cast<PetscInt>(g));
    this->comm().allgather(gfree);
    std::set<PetscInt> band(gfree.begin(), gfree.end());
    libMesh::MeshBase & mesh = _fe_problem.mesh().getMesh();
    _band_elems.clear();
    std::vector<libMesh::dof_id_type> di;
    for (const libMesh::Elem * e : mesh.active_local_element_ptr_range())
    {
      dm.dof_indices(e, di);
      for (const auto g : di)
        if (band.count(static_cast<PetscInt>(g)))
        {
          _band_elems.push_back(e);
          break;
        }
    }
    libMesh::ConstElemRange band_range(&_band_elems);
    _fe_problem.setCurrentAlgebraicElementRange(&band_range); // deep-copies iterators into _band_elems

    if (_verbose)
    {
      unsigned int nb = static_cast<unsigned int>(_band_elems.size()), nf = static_cast<unsigned int>(gfree.size());
      this->comm().sum(nb);
      _console << "    [NEPIN band] sys " << sys_num << " band elems=" << nb << " free=" << nf << std::endl;
    }
  }

  // Stock nonlinear solve of the block under the VI bounds (vinewtonrsls via the sub-executor petsc
  // options; correct re-entrant assembly, unlike a hand-rolled Newton). Assembly is band-restricted
  // above when nepin_restrict_assembly; the reduced-space linear solve is unchanged (bit-identical).
  _fe_problem.solve(sys_num);

  if (_reduced_restrict_assembly)
    _fe_problem.setCurrentAlgebraicElementRange(nullptr); // restore full-domain assembly
}

void
NewtonSNESExecutor::allocateOffDiagMats()
{
  for (const auto i : index_range(_nl_sys_nums))
  {
    auto & sys_i = _fe_problem.getNonlinearSystemBase(_nl_sys_nums[i]);
    for (const auto j : index_range(_nl_sys_nums))
    {
      if (i == j)
        continue;

      auto & sys_j = _fe_problem.getNonlinearSystemBase(_nl_sys_nums[j]);

      const auto m = sys_i.system().n_dofs();
      const auto n = sys_j.system().n_dofs();
      const auto m_l = sys_i.system().n_local_dofs();
      const auto n_l = sys_j.system().n_local_dofs();

      auto mat = std::make_unique<libMesh::PetscMatrix<libMesh::Number>>(_fe_problem.comm());
      // We'll allow the matrix to be hash table assembled
      mat->init_without_preallocation(m, n, m_l, n_l, 1);
      mat->finish_initialization();
      // Allow new nonzeros: the AMR column-constraint distribution in assembleOffDiagJacobian adds
      // parent-column entries (J_ij(:,parent) += w*J_ij(:,hanging)) that need not be in the raw pattern.
      LibmeshPetscCallA(_fe_problem.comm().get(),
                        MatSetOption(mat->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));

      const std::pair<unsigned int, unsigned int> key{i, j};
      auto & [our_mat, tag] = libmesh_map_find(_off_diag_mats, key);
      our_mat = std::move(mat);
    }
  }
}

void
NewtonSNESExecutor::buildMatNest()
{
  const unsigned int n_sys = _fe_problem.numNonlinearSystems();

  if (_mat_nest)
  {
    LibmeshPetscCallA(_fe_problem.comm().get(), MatDestroy(&_mat_nest));
    _mat_nest = nullptr;
  }

  std::vector<Mat> sub_mats(n_sys * n_sys, nullptr);
  for (const auto i : make_range(n_sys))
  {
    auto & sys_i = _fe_problem.getNonlinearSystem(i);
    auto & J_ii =
        static_cast<libMesh::PetscMatrix<libMesh::Number> &>(sys_i.sys().get_system_matrix());
    sub_mats[i * n_sys + i] = J_ii.mat();

    for (unsigned int j = 0; j < n_sys; ++j)
    {
      if (i == j)
        continue;
      const std::pair<unsigned int, unsigned int> key{i, j};
      sub_mats[i * n_sys + j] = libmesh_map_find(_off_diag_mats, key).mat->mat();
    }
  }

  LibmeshPetscCallA(
      _fe_problem.comm().get(),
      MatCreateNest(
          _fe_problem.comm().get(), n_sys, nullptr, n_sys, nullptr, sub_mats.data(), &_mat_nest));
}

void
NewtonSNESExecutor::assembleOffDiagJacobian()
{
  const auto n_sys = _nl_sys_nums.size();
  for (const auto i : make_range(n_sys))
  {
    const auto nl_sys_i_num = _nl_sys_nums[i];
    auto & sys_i = _fe_problem.getNonlinearSystemBase(nl_sys_i_num);
    for (const auto j : make_range(n_sys))
    {
      if (i == j)
        // We already computed and will re-use the sub-SNES diagonal Jacobian blocks
        continue;

      const auto nl_sys_j_num = _nl_sys_nums[j];
      _fe_problem.setJacobianBlockContext(nl_sys_i_num, nl_sys_j_num);
      // NB: we deliberately do NOT force the element-level constraint path here (e.g. via
      // setCurrentNonlinearSystem(i)) -- MOOSE's constrain_element_matrix uses a single (row-system)
      // DofMap and throws "Column too large" on the cross-system columns (col in system j exceeds
      // system i's dof range). That single-DofMap limitation is idaholab/moose#33478. We instead apply
      // the row-system constraints to the assembled global block below (post-processing).
      const std::pair<unsigned int, unsigned int> key{i, j};
      auto & [J_ij, tag] = libmesh_map_find(_off_diag_mats, key);
      _fe_problem.computeJacobianTag(*sys_i.system().current_local_solution, *J_ij, tag);

      // Make the off-diagonal block consistent with the AMR hanging-node constraints, matching the
      // FORM used by dof_map.enforce_constraints_on_jacobian() on the diagonal blocks: a hanging DOF h
      // of the ROW system i keeps its DOF but has its equation replaced by the interpolation
      // constraint (u_h - sum_p w_p u_p = 0), which has NO system-j dependence, so row h of the
      // off-diagonal block J_ij must be ZERO. No parent-row redistribution is applied (that would be
      // the C^T A C form; enforce_constraints_on_jacobian instead sets the constraint row, so the
      // residual/diagonal are NOT redistributed -- redistributing the off-diagonal makes it
      // inconsistent and degrades the SPIN direction, verified). System-j hanging COLUMNS need no
      // treatment for the SPIN direction: Y on those DOFs is fixed by the constraint row of the
      // diagonal block A_jj, so the raw column already yields the correct A_ij*Y (verified: distributing
      // system-j columns leaves both the SPIN direction and the FD Jacobian check unchanged). A fully
      // consistent cross-system block (P_i^T A_ij P_j, both DofMaps) would need a two-DofMap constraint
      // API that neither MOOSE nor libMesh provides; that general framework fix is idaholab/moose#33478.
      auto & dof_map_i = sys_i.system().get_dof_map();
      std::vector<PetscInt> crows;
      for (libMesh::dof_id_type d = dof_map_i.first_dof(); d < dof_map_i.end_dof(); ++d)
        if (dof_map_i.is_constrained_dof(d))
          crows.push_back(static_cast<PetscInt>(d));

      // ALSO zero the row-system's NODAL-BC (Dirichlet) rows. These are applied via the nodal-BC path
      // (setResidual/constrainJacobianRow on the diagonal block), NOT libMesh DofMap constraints, so
      // is_constrained_dof() above MISSES them -- and on a static conforming mesh (no hanging nodes)
      // that leaves crows empty, so the off-diagonal block keeps its Dirichlet rows. Then the coupled
      // Jacobian's constrained row is [I | A_ij[bc,:]] instead of a clean identity, and the Newton/CG
      // step drifts the constrained DOFs by -A_ij[bc,:]*p_j (proportional to the OTHER field's step:
      // ~0 in the elastic regime, large at damage nucleation -> the constrained disp DOFs blow up and
      // the solve stalls). Enumerate them from the nodal-BC warehouse of the row system so the coupled
      // Dirichlet row is a pure identity (row zero in every off-diagonal block).
      const auto sysnum = sys_i.system().number();
      const auto & nodal_bcs = sys_i.getNodalBCWarehouse();
      if (nodal_bcs.hasActiveObjects())
      {
        ConstBndNodeRange & bnd_nodes = *_fe_problem.mesh().getBoundaryNodeRange();
        for (const auto & bnode : bnd_nodes)
        {
          const BoundaryID bid = bnode->_bnd_id;
          const libMesh::Node * node = bnode->_node;
          if (node->processor_id() != this->processor_id())
            continue;
          if (!nodal_bcs.hasActiveBoundaryObjects(bid, 0))
            continue;
          for (const auto & bc : nodal_bcs.getActiveBoundaryObjects(bid, 0))
          {
            const auto varnum = bc->variable().number();
            const auto ncomp = node->n_comp(sysnum, varnum);
            for (unsigned int comp = 0; comp < ncomp; ++comp)
              crows.push_back(static_cast<PetscInt>(node->dof_number(sysnum, varnum, comp)));
          }
        }
      }

      // MatZeroRows is COLLECTIVE: every rank must call it (with its own local constrained rows, which
      // may be zero), else ranks with no constrained DOFs skip it and the others deadlock.
      LibmeshPetscCallA(
          _fe_problem.comm().get(),
          MatZeroRows(
              J_ij->mat(), crows.size(), crows.empty() ? nullptr : crows.data(), 0.0, nullptr, nullptr));
    }
  }
}

PetscErrorCode
NewtonSNESExecutor::outerResidualCallback(SNES /*snes*/, Vec /*x*/, Vec /*f*/, void * ctx)
{
  PetscFunctionBegin;
  auto * ex = static_cast<NewtonSNESExecutor *>(ctx);
  const unsigned int n_sys = ex->_fe_problem.numNonlinearSystems();

  // x and f are VecNests whose sub-Vecs are the per-system solution and RHS Vecs.
  for (unsigned int i = 0; i < n_sys; ++i)
  {
    auto & sys_i = ex->_fe_problem.getNonlinearSystemBase(i);
    ex->_fe_problem.setCurrentNonlinearSystem(i);
    sys_i.computeResidualTag(sys_i.RHS(), sys_i.residualVectorTag());
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
NewtonSNESExecutor::shellMatMult(Mat m, Vec X, Vec Y)
{
  void * ctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(m, &ctx));
  auto * ex = static_cast<NewtonSNESExecutor *>(ctx);
  mooseAssert(ex->_npc_executor,
              "Should not be using a shell matrix without a nonlinear preconditioner");
  PetscCall(ex->_npc_executor->applyBA(ex->_mat_nest, X, Y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

void
NewtonSNESExecutor::scatterToSystems(Vec x)
{
  PetscInt n_sub;
  LibmeshPetscCallA(this->comm().get(), VecNestGetSize(x, &n_sub));
  for (PetscInt i = 0; i < n_sub; ++i)
  {
    Vec sub_x;
    LibmeshPetscCallA(this->comm().get(), VecNestGetSubVec(x, i, &sub_x));
    auto & lm_sys = _fe_problem.getNonlinearSystemBase(static_cast<unsigned int>(i)).system();
    Vec sol_i = cast_ptr<libMesh::PetscVector<libMesh::Number> *>(lm_sys.solution.get())->vec();
    LibmeshPetscCallA(this->comm().get(), VecCopy(sub_x, sol_i));
    lm_sys.update();
  }
}

Real
NewtonSNESExecutor::computeEnergy(Vec x)
{
  scatterToSystems(x);
  // Re-integrate the total-potential-energy postprocessor at the scattered iterate. EXEC_LINESEARCH
  // carries only the energy postprocessor, so this recomputes its material dependency chain
  // (strain -> elasticity -> psie -> psi_total) at x without other exec-flag side effects.
  _fe_problem.execute(EXEC_LINESEARCH);
  return _fe_problem.getPostprocessorValueByName(_energy_pp_name);
}

Real
NewtonSNESExecutor::computeMerit(Vec x)
{
  // Energy merit: total potential energy Psi(x) (leaves _r_plain untouched).
  if (_has_energy_pp)
    return computeEnergy(x);

  // Residual merit 1/2||R(x)||^2. Used where R != grad Psi (e.g. irreversible mixed-mode CZM), so no
  // energy postprocessor exists. Reduced by the active-set mask so the merit is on the same space as
  // the step and the base merit 1/2||reduced R(X)||^2.
  return reducedResidualMerit(x);
}

Real
NewtonSNESExecutor::reducedResidualMerit(Vec x)
{
  // KKT-residual merit 1/2||mask (.) R(x)||^2 on the current active-set mask. computePlainResidual
  // OVERWRITES _r_plain with R(x); zero the active d-entries (mirroring applyActiveSetElimination) so
  // the merit lives on the same REDUCED (inactive) space as the step.
  computePlainResidual(x);
  if (_bounds && _active_mask)
  {
    Vec Rd;
    LibmeshPetscCallA(this->comm().get(), VecNestGetSubVec(_r_plain, _bounded_sys_local, &Rd));
    LibmeshPetscCallA(this->comm().get(), VecPointwiseMult(Rd, Rd, _active_mask));
  }
  PetscReal n;
  LibmeshPetscCallA(this->comm().get(), VecNorm(_r_plain, NORM_2, &n));
  return 0.5 * n * n;
}

void
NewtonSNESExecutor::pdasProjectTrial(Vec W)
{
  // Feasibility projection: clip the bounded variable's block of the trial iterate W onto [d_old, 1].
  // Shared by every globalization (TR / LS) so a bound-constrained step cannot overshoot.
  if (!_bounds)
    return;
  const auto comm = this->comm().get();
  auto & dsys = _fe_problem.getNonlinearSystemBase(_bounded_sys_num);
  const libMesh::NumericVector<libMesh::Number> & d_old = dsys.solutionOld();
  Vec Wd;
  LibmeshPetscCallA(comm, VecNestGetSubVec(W, _bounded_sys_local, &Wd));
  PetscInt wlo, whi;
  LibmeshPetscCallA(comm, VecGetOwnershipRange(Wd, &wlo, &whi));
  PetscScalar * wa;
  LibmeshPetscCallA(comm, VecGetArray(Wd, &wa));
  for (PetscInt i = wlo; i < whi; ++i)
  {
    const PetscReal dlo = d_old(static_cast<libMesh::dof_id_type>(i));
    const PetscReal v = PetscRealPart(wa[i - wlo]);
    wa[i - wlo] = PetscMax(dlo, PetscMin(v, 1.0));
  }
  LibmeshPetscCallA(comm, VecRestoreArray(Wd, &wa));
}

void
NewtonSNESExecutor::pdasRecordStep(Vec step)
{
  // Cache the accepted step's bounded-variable increment (delta d^k) for the next active-set criterion.
  if (!_bounds)
    return;
  Vec sd;
  LibmeshPetscCallA(this->comm().get(), VecNestGetSubVec(step, _bounded_sys_local, &sd));
  LibmeshPetscCallA(this->comm().get(), VecCopy(sd, _pdas_prev_dstep));
}

void
NewtonSNESExecutor::computePlainResidual(Vec x)
{
  scatterToSystems(x);
  const unsigned int n_sys = _fe_problem.numNonlinearSystems();

  // Repoint each system's current-solution pointer at the scattered iterate (the residual assembly
  // reads _current_solution, repointed only by SolverSystem::setSolution()). Do this for ALL
  // systems first so cross-system coupling terms see the correct state.
  for (unsigned int i = 0; i < n_sys; ++i)
  {
    auto & sys_i = _fe_problem.getNonlinearSystemBase(i);
    sys_i.setSolution(*sys_i.system().current_local_solution);
  }
  for (unsigned int i = 0; i < n_sys; ++i)
  {
    auto & sys_i = _fe_problem.getNonlinearSystemBase(i);
    _fe_problem.setCurrentNonlinearSystem(i);

    // FEProblem-level residual assembly, mirroring the normal solve path (and the working
    // outerJacobianCallback, which uses _fe_problem.computeJacobian). Driving the bare system-level
    // sys_i.computeResidualTags() instead produces an identically-zero residual: it SKIPS
    // FEProblemBase::computeResidualTags' pre-assembly setup, and in particular never calls
    // setCurrentResidualVectorTags(). FEProblemBase::addResidual() distributes each element's
    // accumulated residual into exactly the vectors named by currentResidualVectorTags(); with that
    // set empty (cleared by the preceding Jacobian eval's resetState()), every kernel contribution
    // is scattered into no vector, so _Re_non_time stays 0 and residual += *_Re_non_time gives 0.
    // We call computeResidualInternal() (not the computeResidual(soln,residual,nl) wrapper) because
    // only the wrapper touches the shared _fe_vector_tags member + mooseAssert(empty()) that is not
    // reentrant inside the outer SNES line search; computeResidualInternal/computeResidualTags take
    // the tag set by argument and are reentrancy-safe. Associate RHS to the residual tag BEFORE
    // selecting tags so residualVectorTag() is included (required for the residual += _Re_non_time
    // step). Do NOT additionally invoke the bare system-level path here: its element loop caches
    // residual contributions that a later good assembly then double-counts (observed as ratio 0.5).
    sys_i.associateVectorToTag(sys_i.RHS(), sys_i.residualVectorTag());
    std::set<TagID> local_tags;
    _fe_problem.selectVectorTagsFromSystem(
        sys_i, _fe_problem.getVectorTags(Moose::VECTOR_TAG_RESIDUAL), local_tags);
    _fe_problem.computeResidualInternal(
        *sys_i.system().current_local_solution, sys_i.RHS(), local_tags);

    // Copy the freshly assembled residual into _r_plain. The outer SNES function vector F is
    // recomputed (SNESApplyNPC) before it is used for convergence, so temporarily overwriting
    // sys.RHS() here is safe; the energy line search reads grad Psi from _r_plain.
    Vec rhs_i = cast_ptr<libMesh::PetscVector<libMesh::Number> *>(&sys_i.RHS())->vec();
    Vec sub;
    LibmeshPetscCallA(this->comm().get(),
                      VecNestGetSubVec(_r_plain, static_cast<PetscInt>(i), &sub));
    LibmeshPetscCallA(this->comm().get(), VecCopy(rhs_i, sub));
  }
}

bool
NewtonSNESExecutor::energyBacktrack(Vec X,
                                    Vec Y,
                                    Vec W,
                                    Real psiX,
                                    Real slope,
                                    Real minlambda,
                                    PetscInt max_it,
                                    Real & lambda_out,
                                    Real & psiW_out)
{
  const Real c1 = 1e-4; // sufficient-decrease constant (paper)
  // Energy sufficient-decrease cushion. The SPIN direction targets F_SPIN = 0, not the energy
  // minimum; near the solution its Newton step (alpha=1) reduces ||F_precond|| (verified to machine
  // precision) but raises Psi by a hair (residual NPC inexactness -> the F_SPIN fixed point differs
  // from the energy min at ~1e-8 of |Psi|). A strict Armijo test rejects that step and clips
  // alpha->0, freezing ||F_precond|| ~1.5 orders above tolerance. The cushion ignores energy changes
  // below a meaningful fraction of the total energy, so the SPIN step is accepted whenever it does
  // not INCREASE the energy meaningfully -- the correct globalization for a preconditioned-residual
  // Newton direction. (A catastrophic energy increase, e.g. an unstable jump at crack nucleation,
  // still exceeds the cushion and triggers the App. C fallback.)
  const Real energy_tol = 1e-8 * PetscAbsReal(psiX);
  Real lambda = 1.0, lambdaprev = 1.0, psiprev = psiX, psiW = psiX;
  bool ok = false;
  for (PetscInt it = 0; it <= max_it; ++it)
  {
    LibmeshPetscCallA(this->comm().get(), VecWAXPY(W, -lambda, Y, X)); // W = X - lambda*Y
    psiW = computeEnergy(W);

    if (!PetscIsInfOrNanReal(psiW) && psiW <= psiX + c1 * lambda * slope + energy_tol) // energy Armijo
    {
      ok = true;
      break;
    }
    if (lambda <= minlambda)
      break;

    // Interpolate a shorter step from the energy values (quadratic first, cubic thereafter).
    Real lambdatemp;
    if (it == 0 || PetscIsInfOrNanReal(psiW))
      lambdatemp = -slope * lambda * lambda / (2.0 * (psiW - psiX - lambda * slope));
    else
    {
      const Real t1 = psiW - psiX - lambda * slope;
      const Real t2 = psiprev - psiX - lambdaprev * slope;
      const Real a =
          (t1 / (lambda * lambda) - t2 / (lambdaprev * lambdaprev)) / (lambda - lambdaprev);
      const Real b =
          (-lambdaprev * t1 / (lambda * lambda) + lambda * t2 / (lambdaprev * lambdaprev)) /
          (lambda - lambdaprev);
      Real d = b * b - 3.0 * a * slope;
      if (d < 0.0)
        d = 0.0;
      lambdatemp = (a == 0.0) ? -slope / (2.0 * b) : (-b + PetscSqrtReal(d)) / (3.0 * a);
    }
    if (PetscIsInfOrNanReal(lambdatemp))
      lambdatemp = 0.5 * lambda;

    lambdaprev = lambda;
    psiprev = psiW;
    lambda = PetscClipInterval(lambdatemp, 0.1 * lambda, 0.5 * lambda);
  }
  lambda_out = lambda;
  psiW_out = psiW;
  return ok;
}

bool
NewtonSNESExecutor::strongWolfe(Vec X,
                                Vec Y,
                                Vec W,
                                Real psiX,
                                Real phi0p,
                                Real minlambda,
                                PetscInt /*max_it*/,
                                Real & lambda_out,
                                Real & psiW_out)
{
  // Strong Wolfe (Nocedal & Wright Alg. 3.5 bracketing + 3.6 zoom), paper App. B: c1=1e-4, c2=0.9.
  // phi(a) = Psi(X - a*Y); phi'(a) = -R(X - a*Y).Y (grad Psi = R). phi0p = phi'(0) < 0 (descent).
  const Real c1 = 1e-4, c2 = 0.9;
  const Real phi0 = psiX;
  const Real amax = 16.0; // cap so we do not chase overly long steps
  const PetscInt maxit = 40;

  // Evaluate phi(a) (sets W = X - a*Y) and, when needed, phi'(a) from _r_plain.
  auto eval_phi = [&](Real a) -> Real
  {
    LibmeshPetscCallA(this->comm().get(), VecWAXPY(W, -a, Y, X));
    return computeEnergy(W);
  };
  auto eval_phip = [&]() -> Real // W must currently hold X - a*Y
  {
    computePlainResidual(W);
    PetscScalar RdotY;
    LibmeshPetscCallA(this->comm().get(), VecDot(_r_plain, Y, &RdotY));
    return -PetscRealPart(RdotY);
  };

  auto accept = [&](Real a, Real phi_a) -> bool
  {
    lambda_out = a;
    psiW_out = phi_a;
    LibmeshPetscCallA(this->comm().get(), VecWAXPY(W, -a, Y, X)); // leave W = X - a*Y
    return true;
  };

  // ---- zoom(a_lo, a_hi): a_lo satisfies Armijo and has the lower phi ----
  auto zoom =
      [&](Real a_lo, Real a_hi, Real phi_lo, Real /*phi_hi*/, Real & a_out, Real & phi_out) -> bool
  {
    for (PetscInt j = 0; j < maxit; ++j)
    {
      const Real a_j = 0.5 * (a_lo + a_hi); // bisection (robust)
      const Real phi_j = eval_phi(a_j);
      if (PetscIsInfOrNanReal(phi_j) || phi_j > phi0 + c1 * a_j * phi0p || phi_j >= phi_lo)
        a_hi = a_j;
      else
      {
        const Real phip_j = eval_phip();
        if (PetscAbsReal(phip_j) <= -c2 * phi0p) // strong Wolfe satisfied
        {
          a_out = a_j;
          phi_out = phi_j;
          return true;
        }
        if (phip_j * (a_hi - a_lo) >= 0.0)
          a_hi = a_lo;
        a_lo = a_j;
        phi_lo = phi_j;
      }
      if (PetscAbsReal(a_hi - a_lo) < minlambda)
        break;
    }
    // Fall back to a_lo (Armijo-satisfying, lower phi).
    a_out = a_lo;
    phi_out = eval_phi(a_lo);
    return (a_lo > 0.0 && !PetscIsInfOrNanReal(phi_out) && phi_out <= phi0 + c1 * a_lo * phi0p);
  };

  // ---- bracketing ----
  Real a_prev = 0.0, phi_prev = phi0, a = 1.0;
  for (PetscInt i = 0; i < maxit; ++i)
  {
    const Real phi_a = eval_phi(a);
    if (PetscIsInfOrNanReal(phi_a) || phi_a > phi0 + c1 * a * phi0p || (i > 0 && phi_a >= phi_prev))
    {
      Real a_out, phi_out;
      const bool ok = zoom(a_prev, a, phi_prev, phi_a, a_out, phi_out);
      return ok ? accept(a_out, phi_out) : (lambda_out = a_out, psiW_out = phi_out, false);
    }
    const Real phip_a = eval_phip(); // W currently holds X - a*Y
    if (PetscAbsReal(phip_a) <= -c2 * phi0p)
      return accept(a, phi_a); // strong Wolfe satisfied at a
    if (phip_a >= 0.0)
    {
      Real a_out, phi_out;
      const bool ok = zoom(a, a_prev, phi_a, phi_prev, a_out, phi_out);
      return ok ? accept(a_out, phi_out) : (lambda_out = a_out, psiW_out = phi_out, false);
    }
    a_prev = a;
    phi_prev = phi_a;
    a = PetscMin(2.0 * a, amax);
    if (a_prev >= amax) // already at the cap and still increasing: accept the capped step
      return accept(a_prev, phi_prev);
  }
  lambda_out = a_prev;
  psiW_out = phi_prev;
  return false;
}

PetscErrorCode
NewtonSNESExecutor::outerConvergenceTest(SNES snes,
                                         PetscInt it,
                                         PetscReal /*xnorm*/,
                                         PetscReal /*ynorm*/,
                                         PetscReal /*fnorm*/,
                                         SNESConvergedReason * reason,
                                         void * ctx)
{
  PetscFunctionBeginUser;
  auto * ex = static_cast<NewtonSNESExecutor *>(ctx);
  *reason = SNES_CONVERGED_ITERATING;

  PetscReal atol, rtol, stol;
  PetscInt maxit, maxf;
  PetscCall(SNESGetTolerances(snes, &atol, &rtol, &stol, &maxit, &maxf));

  // Converge on ||R|| = ||grad Psi|| cached from the line search (one outer iteration in arrears,
  // which is inconsequential). _r_norm/_r0_norm are set once the first line search of this solve has
  // run; at it==0 (before any line search) we simply keep iterating. For the bound-constrained PDAS
  // path (HWW Remark 3.3) convergence additionally requires the active set to have stopped changing:
  // ||R|| here is the reduced (inactive) residual, which can be tiny while the active set still moves.
  if (it > 0 && ex->_r_norm >= 0.0 && ex->_r0_norm > 0.0)
  {
    // (A) Absolute KKT satisfaction ALWAYS counts: a reduced (inactive) residual below atol means the
    // stationarity+complementarity residual is at the noise floor, so any residual-level active-set
    // churn (borderline DOFs with sub-tolerance multipliers) is meaningless -- accept regardless of
    // _active_set_changed. This clears the spurious "converged residual but active set still flipping"
    // stall (reduced ||R|| ~ 1e-10 blocked by a 16-DOF flip, then a degenerate rho rejection).
    if (ex->_r_norm < atol)
      *reason = SNES_CONVERGED_FNORM_ABS;
    // Relative convergence still requires a settled active set (HWW Remark 3.3): a modest relative drop
    // with the set still moving is not yet a PDAS solution.
    else if ((!ex->_bounds || !ex->_active_set_changed) && ex->_r_norm < rtol * ex->_r0_norm)
      *reason = SNES_CONVERGED_FNORM_RELATIVE;
  }
  if (*reason == SNES_CONVERGED_ITERATING && it >= maxit)
    *reason = SNES_DIVERGED_MAX_IT;

  if (ex->_verbose)
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
                          "    [gradPsi conv] it=%d ||R||=%.4e rtol*||R0||=%.4e atol=%.4e reason=%d\n",
                          (int)it,
                          (double)(ex->_r_norm),
                          (double)(rtol * ex->_r0_norm),
                          (double)atol,
                          (int)*reason));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
NewtonSNESExecutor::spinLineSearch(SNESLineSearch ls, void * ctx)
{
  PetscFunctionBeginUser;
  auto * ex = static_cast<NewtonSNESExecutor *>(ctx);

  SNES snes;
  Vec X, F, Y, W, G;
  PetscReal fnorm, ynorm, minlambda;
  PetscInt max_it;
  PetscCall(SNESLineSearchGetSNES(ls, &snes));
  // X = base iterate, F = F_SPIN(X) = X - NPC(X); Y = search direction (solving J_SPIN Y = F);
  // W, G are scratch (trial solution / trial residual).
  PetscCall(SNESLineSearchGetVecs(ls, &X, &F, &Y, &W, &G));
  PetscCall(SNESLineSearchGetNorms(ls, nullptr, &fnorm, &ynorm));
  PetscCall(
      SNESLineSearchGetTolerances(ls, &minlambda, nullptr, nullptr, nullptr, nullptr, &max_it));

  // Enter this block for the energy-merit line search (App. B/C) AND for the Steihaug trust region
  // regardless of merit -- the TR shares the same prefix (tr_npc sweep, plain-residual assembly, PDAS
  // elimination) and the residual-merit TR reuses the whole Steihaug machinery below, swapping only
  // the ratio-test merit (energy Psi -> 1/2||R||^2). The 1-D TR and App-C line search (further down)
  // are energy-only and are unreachable here without _has_energy_pp (the Steihaug block returns).
  if (ex->_has_energy_pp || (ex->_use_trust_region && ex->_tr_steihaug))
  {
    // ============ energy-merit line search (App. B + C) / trust-region shared prefix ============
    // Merit is the total potential energy Psi (C^1 across the penalty kink), evaluated at the OUTER
    // iterate W = X - alpha*Y. grad Psi = plain residual R = _r_plain. App. C selects the search
    // direction UP FRONT by the ascent test (paper): use the SPIN direction if it is a descent
    // direction for Psi; else the full-J inexact-Newton direction J^{-1} R; else (Cauchy-scaled)
    // steepest descent -R. A SINGLE strong-Wolfe line search (App. B, c1=1e-4, c2=0.9) then
    // globalizes the chosen, guaranteed-descent direction. (Selecting by ascent -- rather than
    // trying each tier and falling through on line-search failure -- guarantees a descent direction,
    // so the line search cannot fail on all tiers.)

    // (B) NPC-every-iteration trust region: apply ONE multiplicative-Schwarz NPC sweep and ADOPT the
    // swept iterate NPC(X) as the outer iterate, then re-assemble grad^2 Psi there. We invoke the NMSM
    // shell SNES DIRECTLY (SNESSolve sweeps the VecNest in place, X <- NPC(X)) rather than
    // SNESApplyNPC -- the NPC is deliberately NOT attached to the outer SNES (see run()), so PETSc's
    // NEWTONLS does not also auto-apply it (which doubled the sweep). The residual/energy below are
    // then taken at NPC(X) and the Steihaug TRS steps from there -- the "trust-region version of
    // MSPIN" (NPC cadence of the App. C line search, TR globalization). The Jacobian callback skipped
    // the raw-X assembly for this mode, so assemble here.
    if (ex->_use_trust_region && ex->_tr_steihaug && ex->_tr_npc)
    {
      if (ex->_bounds)
      {
        // Bound-constrained (B): compute the PDAS active set at the PRE-SWEEP iterate X, then arm the
        // phase-field sub-executor to solve only the INACTIVE subsystem (active DOFs held at d_old).
        // This makes the NPC sweep bound-consistent (no unconstrained d excursion). The same set is
        // eliminated from the coupled operator below (computeActiveSet is not called again this
        // iteration -- see the guarded block after the sweep).
        ex->computePlainResidual(X);
        ex->computeActiveSet(X); // for the coupled Steihaug step below (all MSPIN/NEPIN modes)
        if (ex->_nepin)
          ex->computeHardSet(X); // NEPIN: process-zone band -> easy-set frozen lists (per field)
        PetscReal r_pre;
        PetscCall(VecNorm(ex->_r_plain, NORM_2, &r_pre)); // full coupled residual BEFORE the sweep
        auto * const nmsm = dynamic_cast<NMSMExecutor *>(ex->_npc_executor);
        if (!nmsm)
          ex->mooseError("bounds=true with tr_npc requires an NMSMExecutor nonlinear preconditioner "
                         "(nl_preconditioning).");

        // ---- HIK merit safeguard on the multiplicative NPC sweep --------------------------------
        // The sweep is a semismooth-Newton block step armed by the active set A(X). Under rapid crack
        // propagation A(X) can freeze DOFs that should move, so the swept iterate NPC(X) carries a
        // LARGER KKT residual than X; adopting it unconditionally (outside the trust region's rho
        // test) sustains a limit cycle (the active set oscillates, reduced ||R|| orbits at ~1e-2).
        // Globalize the sweep by a backtracking line search on the KKT-residual merit
        // phi = 1/2||mask (.) R||^2 (NOT the energy Psi: the block minimization lowers Psi while
        // raising ||R||, so Psi would not see the oscillation). Take the full sweep, then accept the
        // largest t in {1, 1/2, ...} with phi(X_pre + t (NPC(X)-X_pre)) <= phi(X_pre); t=0 rejects a
        // pure-ascent sweep and steps from X_pre (the rho-guarded TR below then makes the progress).
        // Monotone phi across outer iterations -> no limit cycle. In the healthy regime the full
        // sweep reduces phi and t=1 is taken with a single extra residual evaluation.
        if (!ex->_x_presweep)
          PetscCall(VecDuplicate(X, &ex->_x_presweep));
        if (!ex->_sweep_dir)
          PetscCall(VecDuplicate(X, &ex->_sweep_dir));
        const PetscReal phi0 = ex->reducedResidualMerit(X); // 1/2||mask (.) R(X_pre)||^2
        PetscCall(VecCopy(X, ex->_x_presweep));             // stash pre-sweep iterate

        // Arm the sub-solves according to the sweep's SCOPE (npc family):
        //  MSPIN      : freeze the active set on pf (d=d_old there), solve the whole inactive damage
        //               field, and sweep every block (disp then pf).
        //  mspin_pne  : freeze the EASY damage DOFs at current so only the process-zone band is free,
        //               and restrict the sweep to the pf block alone (disp left to the outer TR).
        //  mspin_fne  : additionally freeze the far-field displacement (unbounded box on the free
        //               front band) and run the full multiplicative sweep (disp-band -> pf-band), so
        //               u and d move together on the front -- the coupled elimination.
        // Arm the block sub-solves for the elimination sweep (npc family):
        //  mspin      : freeze the active set on pf; solve the whole inactive damage field; sweep all blocks.
        //  mspin_pne  : freeze the easy damage DOFs -> only the process-zone band is free; pf block only.
        //  mspin_fne  : additionally free the disp front (far-field frozen); full sweep (both fields).
        //               By default (nepin_restrict_assembly) each block sub-solve restricts its FE
        //               ASSEMBLY to the elements incident to its free set -- band-sized assembly,
        //               otherwise identical (same vinewtonrsls reduced-space solve, bit-identical iterates).
        const bool restrict_asm = ex->_nepin_full && ex->_nepin_restrict_assembly;
        if (ex->_nepin)
        {
          nmsm->armBoundedSubSolve(ex->_bounded_sys_num,
                                   ex->_nepin_frozen,
                                   ex->_nepin_frozen_val,
                                   /*bounded_box=*/true,
                                   restrict_asm);
          if (ex->_nepin_full)
            nmsm->armBoundedSubSolve(ex->_disp_sys_num,
                                     ex->_nepin_frozen_u,
                                     ex->_nepin_frozen_u_val,
                                     /*bounded_box=*/false,
                                     restrict_asm);
          else
            nmsm->setNepinOnly(ex->_bounded_sys_num); // PARTIAL: pf-only, no disp sweep
        }
        else
          nmsm->armBoundedSubSolve(ex->_bounded_sys_num, ex->_active_d_dofs, ex->_active_val);
        PetscCall(SNESSolve(ex->_npc_executor->getSNES(), nullptr, X)); // X <- NPC(X)/E(X)
        if (ex->_nepin && !ex->_nepin_full)
          nmsm->clearNepinOnly();
        nmsm->disarmBoundedSubSolve();
        PetscCall(VecWAXPY(ex->_sweep_dir, -1.0, ex->_x_presweep, X)); // sweep_dir = NPC(X) - X_pre

        PetscReal phi = ex->reducedResidualMerit(X); // phi at the full sweep (t = 1)
        PetscReal t = 1.0;
        PetscInt bt = 0;
        const PetscInt max_bt = 12;
        while (phi > phi0 && bt < max_bt) // sweep raised the KKT merit -> backtrack
        {
          t *= 0.5;
          ++bt;
          PetscCall(VecWAXPY(X, t, ex->_sweep_dir, ex->_x_presweep)); // X = X_pre + t (NPC(X)-X_pre)
          phi = ex->reducedResidualMerit(X);
        }
        if (phi > phi0) // even the shortest tried step ascends -> reject the sweep entirely
        {
          PetscCall(VecCopy(ex->_x_presweep, X));
          t = 0.0;
        }

        if (ex->_verbose)
        {
          ex->computePlainResidual(X); // (recomputed unconditionally at line ~1505; here only for print)
          PetscReal r_post;
          PetscCall(VecNorm(ex->_r_plain, NORM_2, &r_post)); // full coupled residual AFTER safeguard
          PetscCall(PetscPrintf(PETSC_COMM_WORLD,
                                "    [PDAS sweep] full||R|| pre=%.4e post=%.4e (ratio %.2f) t=%.3g bt=%d\n",
                                (double)r_pre,
                                (double)r_post,
                                (double)(r_post / PetscMax(r_pre, 1e-300)),
                                (double)t,
                                (int)bt));
        }
      }
      else
        PetscCall(SNESSolve(ex->_npc_executor->getSNES(), nullptr, X)); // X <- NPC(X) (one sweep)
      PetscCall(ex->assembleCoupledJacobian(X)); // grad^2 Psi at NPC(X)
    }

    ex->computePlainResidual(X); // _r_plain <- R(X) = grad Psi(X)  (F holds F_SPIN, do not use it)

    // PDAS bound constraint: eliminate the active set (dynamic Dirichlet) from the coupled Hessian and
    // residual so the TR/CG step below has delta d = 0 on the active set. Done BEFORE the norm so ||R||
    // is the reduced (inactive) residual (HWW Remark 3.4), and BEFORE steihaugTRS so the block
    // preconditioners refactor the reduced blocks. The set is computed here for the non-NPC paths; the
    // (B) tr_npc path already computed it at the pre-sweep iterate (used to freeze the NPC pf sweep).
    if (ex->_bounds)
    {
      const bool tr_npc_path = ex->_use_trust_region && ex->_tr_steihaug && ex->_tr_npc;
      if (!tr_npc_path)
        ex->computeActiveSet(X);
      ex->clampActiveToBounds(X); // project active DOFs onto their bound (d_old lower, 1 upper)
      // Diagnostic: the bounded variable's range on the current iterate (should stay >= d_old) and
      // the full (pre-elimination) coupled residual. If d drops below the lower bound, a DOF that
      // should be frozen is free -> the step is not respecting the constraint.
      if (ex->_verbose)
      {
        Vec Xd;
        PetscReal dmn, dmx, rfull;
        PetscCall(VecNestGetSubVec(X, ex->_bounded_sys_local, &Xd));
        PetscCall(VecMin(Xd, nullptr, &dmn));
        PetscCall(VecMax(Xd, nullptr, &dmx));
        PetscCall(VecNorm(ex->_r_plain, NORM_2, &rfull));
        PetscCall(PetscPrintf(PETSC_COMM_WORLD,
                              "    [PDAS diag] d in [%.4e, %.4e]  full||R||=%.4e\n",
                              (double)dmn,
                              (double)dmx,
                              (double)rfull));
      }
      PetscCall(ex->applyActiveSetElimination());
    }

    PetscReal gradPsiNorm;
    PetscCall(VecNorm(ex->_r_plain, NORM_2, &gradPsiNorm)); // ||reduced R|| (= ||grad Psi|| for energy)
    // Cache for outerConvergenceTest (the paper converges on this coupled residual). _r0_norm is the
    // baseline set at the first line search of the solve; _r_norm is the latest value.
    ex->_r_norm = gradPsiNorm;
    if (ex->_r0_norm < 0.0)
      ex->_r0_norm = gradPsiNorm;
    // Base merit for the trust-region ratio test: total potential energy Psi(X) if provided, else the
    // residual merit 1/2||R(X)||^2 (reduced; reuses gradPsiNorm -- no extra assembly). The trial
    // point uses computeMerit(W).
    const PetscReal psiX =
        ex->_has_energy_pp ? ex->computeEnergy(X) : 0.5 * gradPsiNorm * gradPsiNorm;
    // Residual-merit TR: stash the reduced R(X) so it survives the ratio-test retry loop (steihaugTRS
    // re-reads _r_plain=R(X) each rejected re-solve; the per-trial computeMerit(W) overwrites it).
    if (!ex->_has_energy_pp && ex->_r_base)
      PetscCall(VecCopy(ex->_r_plain, ex->_r_base));

    if (ex->_has_energy_pp && std::getenv("SPIN_AUDIT"))
    { // ===== SPIN INFRA AUDIT: is _mat_nest the SYMMETRIC, correct coupled Jacobian dR/dx? =====
      // (0) FD GRADIENT check: is g = _r_plain actually grad(Psi)? The trust-region ratio test rests
      // ENTIRELY on this: rho = ared/pred with ared = Psi(X)-Psi(X+p), pred = -(g^T p + 1/2 p^T A p).
      // If R != grad(Psi) even on the active subspace, rho can never approach 1 and the model is junk.
      // Test: (Psi(X + eps g) - Psi(X))/eps  ?=  g^T g  (directional derivative of Psi along g).
      {
        Vec gcopy;
        PetscCall(VecDuplicate(X, &gcopy));
        PetscCall(VecCopy(ex->_r_plain, gcopy)); // gcopy = g = R(X)
        PetscReal gnrm, xnrm;
        PetscCall(VecNorm(gcopy, NORM_2, &gnrm));
        PetscCall(VecNorm(X, NORM_2, &xnrm));
        const PetscReal epsg = 1e-7 * (xnrm + 1.0) / PetscMax(gnrm, 1e-300);
        PetscCall(VecWAXPY(W, epsg, gcopy, X));         // W = X + eps g
        const PetscReal psiWg = ex->computeEnergy(W);   // NOTE: clobbers _r_plain with R(W)
        const PetscReal fdg = (psiWg - psiX) / epsg;    // ~ g^T g if R = grad Psi
        const PetscReal gTg = gnrm * gnrm;
        ex->computePlainResidual(X);                    // restore _r_plain = R(X)
        PetscCall(PetscPrintf(PETSC_COMM_WORLD,
                              "  [spinGrad] dPsi/dg=%.6e  g^Tg=%.6e  ratio=%.4f (=1 iff R=gradPsi)\n",
                              (double)fdg,
                              (double)gTg,
                              (double)(fdg / PetscMax(gTg, 1e-300))));
        PetscCall(VecDestroy(&gcopy));
      }
      Vec v, w, Av, Aw, R0, FDv;
      PetscCall(VecDuplicate(X, &v));
      PetscCall(VecDuplicate(X, &w));
      PetscCall(VecDuplicate(X, &Av));
      PetscCall(VecDuplicate(X, &Aw));
      PetscCall(VecDuplicate(X, &R0));
      PetscCall(VecDuplicate(X, &FDv));
      PetscRandom rng;
      PetscCall(PetscRandomCreate(PETSC_COMM_WORLD, &rng));
      PetscCall(PetscRandomSetFromOptions(rng));
      PetscCall(VecSetRandom(v, rng));
      PetscCall(VecSetRandom(w, rng));
      // (1) SYMMETRY: w^T A v vs v^T A w. A = grad^2 Psi must be symmetric (A01 = A10^T); the
      // Steihaug CG is only valid for a symmetric operator, so an asymmetric _mat_nest is a bug.
      PetscCall(MatMult(ex->_mat_nest, v, Av));
      PetscCall(MatMult(ex->_mat_nest, w, Aw));
      PetscScalar wAv, vAw;
      PetscCall(VecDot(w, Av, &wAv));
      PetscCall(VecDot(v, Aw, &vAw));
      const PetscReal sym =
          PetscAbsReal(PetscRealPart(wAv - vAw)) /
          (PetscAbsReal(PetscRealPart(wAv)) + PetscAbsReal(PetscRealPart(vAw)) + 1e-300);
      // Does the asymmetry ACTIVATE along an active-subspace vector? R = _r_plain has R_bc = 0 (the
      // BCs are satisfied/preset), so if ||(A - A^T) R|| / ||A R|| is small, the (Dirichlet-BC)
      // asymmetry lives only on constrained rows/cols and is BENIGN for the solve (p_bc = 0);
      // if it is O(sym), the asymmetry genuinely corrupts the Newton/CG step.
      PetscReal asym_act = -1.0;
      {
        // R and AR both have _bc = 0 (BCs satisfied), so this stays on the active subspace and
        // EXCLUDES the benign Dirichlet-BC row asymmetry: t1 = R^T A(AR) = R^T A^2 R, t2 = (AR)^T(AR)
        // = R^T A^T A R. Equal iff A is symmetric on the active subspace (where the CG/energy solve
        // actually operates). If ~0 the 12% is purely benign BC; if O(0.1) it corrupts the step.
        Vec AR, AAR;
        PetscCall(VecDuplicate(X, &AR));
        PetscCall(VecDuplicate(X, &AAR));
        PetscCall(MatMult(ex->_mat_nest, ex->_r_plain, AR)); // AR = A R
        PetscCall(MatMult(ex->_mat_nest, AR, AAR));          // AAR = A (A R)
        PetscScalar t1, t2;
        PetscCall(VecDot(ex->_r_plain, AAR, &t1)); // R^T A^2 R
        PetscCall(VecDot(AR, AR, &t2));            // R^T A^T A R
        asym_act = PetscAbsReal(PetscRealPart(t1 - t2)) /
                   (PetscAbsReal(PetscRealPart(t1)) + PetscAbsReal(PetscRealPart(t2)) + 1e-300);
        PetscCall(VecDestroy(&AR));
        PetscCall(VecDestroy(&AAR));
      }
      // (2) FD JACOBIAN along a RANDOM v: A v vs (R(X+eps v) - R(X))/eps (exercises the FULL coupled
      // operator incl. off-diagonals, unlike the v=R check). Per block: 0=disp, 1=pf.
      PetscCall(VecCopy(ex->_r_plain, R0)); // R(X)
      PetscReal xn, vn;
      PetscCall(VecNorm(X, NORM_2, &xn));
      PetscCall(VecNorm(v, NORM_2, &vn));
      const PetscReal eps = 1e-7 * (xn + 1.0) / vn;
      PetscCall(VecWAXPY(W, eps, v, X));
      ex->computePlainResidual(W); // _r_plain = R(X + eps v)
      PetscCall(VecWAXPY(FDv, -1.0, R0, ex->_r_plain));
      PetscCall(VecScale(FDv, 1.0 / eps)); // FDv = dR/dx . v
      PetscCall(VecAXPY(Av, -1.0, FDv));   // Av = A v - dR/dx v
      PetscReal nfd, ea, e0 = -1.0, e1 = -1.0;
      PetscCall(VecNorm(FDv, NORM_2, &nfd));
      PetscCall(VecNorm(Av, NORM_2, &ea));
      {
        Vec a0, a1, f0, f1;
        PetscReal na0, na1, nf0, nf1;
        PetscCall(VecNestGetSubVec(Av, 0, &a0));
        PetscCall(VecNestGetSubVec(Av, 1, &a1));
        PetscCall(VecNestGetSubVec(FDv, 0, &f0));
        PetscCall(VecNestGetSubVec(FDv, 1, &f1));
        PetscCall(VecNorm(a0, NORM_2, &na0));
        PetscCall(VecNorm(a1, NORM_2, &na1));
        PetscCall(VecNorm(f0, NORM_2, &nf0));
        PetscCall(VecNorm(f1, NORM_2, &nf1));
        e0 = (nf0 > 0) ? na0 / nf0 : na0;
        e1 = (nf1 > 0) ? na1 / nf1 : na1;
      }
      // (3) PER-BLOCK symmetry: locate the asymmetry (coupling A01=A10^T already verified; check the
      // diagonal blocks A00 (disp) and A11 (pf)).
      PetscReal s00 = -1.0, s11 = -1.0;
      // A00: vd=(v0,0), wd=(w0,0) -> A00 v0 = (A vd)_disp, A00 w0 = (A wd)_disp.
      PetscCall(VecCopy(v, Av));
      PetscCall(VecCopy(w, Aw));
      {
        Vec p;
        PetscCall(VecNestGetSubVec(Av, 1, &p));
        PetscCall(VecSet(p, 0.0)); // Av = (v0,0)
        PetscCall(VecNestGetSubVec(Aw, 1, &p));
        PetscCall(VecSet(p, 0.0)); // Aw = (w0,0)
      }
      PetscCall(MatMult(ex->_mat_nest, Av, R0));  // R0  = A (v0,0)
      PetscCall(MatMult(ex->_mat_nest, Aw, FDv)); // FDv = A (w0,0)
      {
        Vec r0d, fdd, v0, w0;
        PetscScalar a, b;
        PetscCall(VecNestGetSubVec(R0, 0, &r0d));  // A00 v0
        PetscCall(VecNestGetSubVec(FDv, 0, &fdd)); // A00 w0
        PetscCall(VecNestGetSubVec(v, 0, &v0));
        PetscCall(VecNestGetSubVec(w, 0, &w0));
        PetscCall(VecDot(w0, r0d, &a)); // w0^T A00 v0
        PetscCall(VecDot(v0, fdd, &b)); // v0^T A00 w0
        s00 = PetscAbsReal(PetscRealPart(a - b)) /
              (PetscAbsReal(PetscRealPart(a)) + PetscAbsReal(PetscRealPart(b)) + 1e-300);
      }
      // A11 + coupling: vp=(0,v1), wp=(0,w1).
      PetscCall(VecCopy(v, Av));
      PetscCall(VecCopy(w, Aw));
      {
        Vec p;
        PetscCall(VecNestGetSubVec(Av, 0, &p));
        PetscCall(VecSet(p, 0.0)); // Av = (0,v1)
        PetscCall(VecNestGetSubVec(Aw, 0, &p));
        PetscCall(VecSet(p, 0.0)); // Aw = (0,w1)
      }
      PetscCall(MatMult(ex->_mat_nest, Av, R0));  // R0  = A (0,v1) -> (R0)_pf=A11 v1, (R0)_disp=A01 v1
      PetscCall(MatMult(ex->_mat_nest, Aw, FDv)); // FDv = A (0,w1)
      {
        Vec r0p, fdp, v1v, w1v;
        PetscScalar a, b;
        PetscCall(VecNestGetSubVec(R0, 1, &r0p));  // A11 v1
        PetscCall(VecNestGetSubVec(FDv, 1, &fdp)); // A11 w1
        PetscCall(VecNestGetSubVec(v, 1, &v1v));
        PetscCall(VecNestGetSubVec(w, 1, &w1v));
        PetscCall(VecDot(w1v, r0p, &a)); // w1^T A11 v1
        PetscCall(VecDot(v1v, fdp, &b)); // v1^T A11 w1
        s11 = PetscAbsReal(PetscRealPart(a - b)) /
              (PetscAbsReal(PetscRealPart(a)) + PetscAbsReal(PetscRealPart(b)) + 1e-300);
      }
      ex->computePlainResidual(X); // restore _r_plain = R(X)
      PetscCall(PetscPrintf(
          PETSC_COMM_WORLD,
          "  [spinAudit] ||gradPsi||=%.3e sym(A)=%.3e asymActive(R)=%.3e | s(A00)=%.3e s(A11)=%.3e | "
          "FD=%.3e (disp=%.3e pf=%.3e)\n",
          (double)gradPsiNorm,
          (double)sym,
          (double)asym_act,
          (double)s00,
          (double)s11,
          (double)(nfd > 0 ? ea / nfd : ea),
          (double)e0,
          (double)e1));
      PetscCall(PetscRandomDestroy(&rng));
      PetscCall(VecDestroy(&v));
      PetscCall(VecDestroy(&w));
      PetscCall(VecDestroy(&Av));
      PetscCall(VecDestroy(&Aw));
      PetscCall(VecDestroy(&R0));
      PetscCall(VecDestroy(&FDv));
    }

    if (ex->_has_energy_pp && std::getenv("SPIN_DIAG"))
    { // ===== SPIN-DIRECTION DIAGNOSTIC (opt-in; does an extra full-J solve per line search) =====
      // Exact energy-Newton direction: solve A z = R accurately (z = A^{-1}R) into G.
      PetscCall(KSPSetOperators(ex->_fullJ_ksp, ex->_mat_nest, ex->_mat_nest));
      PetscCall(KSPSetTolerances(ex->_fullJ_ksp, 1e-11, 1e-14, PETSC_DEFAULT, 2000)); // tight ref for cos
      PetscCall(KSPSolve(ex->_fullJ_ksp, ex->_r_plain, G)); // G = A^{-1} R (exact energy-Newton dir)
      KSPConvergedReason kr_full;
      PetscInt kits_full;
      PetscCall(KSPGetConvergedReason(ex->_fullJ_ksp, &kr_full));
      PetscCall(KSPGetIterationNumber(ex->_fullJ_ksp, &kits_full));
      PetscScalar RAiR;
      PetscCall(VecDot(ex->_r_plain, G, &RAiR)); // kappa_inv = R.A^{-1}R (>0 => -A^{-1}R descent)
      PetscReal gnorm_ex, ynorm_in;
      PetscCall(VecNorm(G, NORM_2, &gnorm_ex)); // ||A^{-1}R||
      PetscCall(VecNorm(Y, NORM_2, &ynorm_in)); // ||Y_spin||
      PetscScalar YdotG;
      PetscCall(VecDot(Y, G, &YdotG)); // Y . A^{-1}R
      const PetscReal cosYG =
          (ynorm_in > 0 && gnorm_ex > 0) ? PetscRealPart(YdotG) / (ynorm_in * gnorm_ex) : 0.0;
      // Outer Krylov that produced Y: did it actually converge (M^{-1}A) Y = F_SPIN ?
      KSP outer_ksp;
      KSPConvergedReason kr_out;
      PetscInt kits_out;
      PetscReal krnorm_out;
      PetscCall(SNESGetKSP(snes, &outer_ksp));
      PetscCall(KSPGetConvergedReason(outer_ksp, &kr_out));
      PetscCall(KSPGetIterationNumber(outer_ksp, &kits_out));
      PetscCall(KSPGetResidualNorm(outer_ksp, &krnorm_out));
      // RHS-consistency check: is F_SPIN = M^{-1} R ?  Equivalent: M F_SPIN = R with the
      // block-lower-triangular M = [[A00,0],[A10,A11]] (uses only the trusted _mat_nest blocks).
      PetscReal mfr_rel = -1.0;
      {
        Mat A00, A10, A11;
        Vec F0, F1, W0, W1;
        PetscCall(MatNestGetSubMat(ex->_mat_nest, 0, 0, &A00));
        PetscCall(MatNestGetSubMat(ex->_mat_nest, 1, 0, &A10));
        PetscCall(MatNestGetSubMat(ex->_mat_nest, 1, 1, &A11));
        if (A00 && A10 && A11)
        {
          PetscCall(VecNestGetSubVec(F, 0, &F0));
          PetscCall(VecNestGetSubVec(F, 1, &F1));
          PetscCall(VecNestGetSubVec(W, 0, &W0));
          PetscCall(VecNestGetSubVec(W, 1, &W1));
          PetscCall(MatMult(A00, F0, W0));        // W0 = A00 F0
          PetscCall(MatMult(A10, F0, W1));        // W1 = A10 F0
          PetscCall(MatMultAdd(A11, F1, W1, W1)); // W1 += A11 F1  => W = M F_SPIN
          PetscCall(VecAXPY(W, -1.0, ex->_r_plain)); // W = M F_SPIN - R
          PetscReal mfr;
          PetscCall(VecNorm(W, NORM_2, &mfr));
          mfr_rel = (gradPsiNorm > 0) ? mfr / gradPsiNorm : mfr;
        }
      }
      // OPERATOR-consistency check: is the shell truly M^{-1}A ?  Test M*(shell*R) == A*R.
      // If applyBA's block solves use a STALE factorization M_stale != M (KSPSetOperators skipped on
      // unchanged pointer), then M*(M_stale^{-1} A R) != A R and this ratio is large.
      PetscReal shell_rel = -1.0;
      {
        Mat A00, A10, A11;
        PetscCall(MatNestGetSubMat(ex->_mat_nest, 0, 0, &A00));
        PetscCall(MatNestGetSubMat(ex->_mat_nest, 1, 0, &A10));
        PetscCall(MatNestGetSubMat(ex->_mat_nest, 1, 1, &A11));
        if (A00 && A10 && A11)
        {
          Vec sR, AR, MsR;
          PetscCall(VecDuplicate(ex->_r_plain, &sR));
          PetscCall(VecDuplicate(ex->_r_plain, &AR));
          PetscCall(VecDuplicate(ex->_r_plain, &MsR));
          PetscCall(MatMult(ex->_jac_shell, ex->_r_plain, sR)); // sR = shell*R = M_stale^{-1} A R
          PetscCall(MatMult(ex->_mat_nest, ex->_r_plain, AR));   // AR = A R
          Vec s0, s1, m0, m1;
          PetscCall(VecNestGetSubVec(sR, 0, &s0));
          PetscCall(VecNestGetSubVec(sR, 1, &s1));
          PetscCall(VecNestGetSubVec(MsR, 0, &m0));
          PetscCall(VecNestGetSubVec(MsR, 1, &m1));
          PetscCall(MatMult(A00, s0, m0));        // m0 = A00 s0
          PetscCall(MatMult(A10, s0, m1));        // m1 = A10 s0
          PetscCall(MatMultAdd(A11, s1, m1, m1)); // m1 += A11 s1 => MsR = M*(shell*R)
          PetscCall(VecAXPY(MsR, -1.0, AR));      // MsR = M*shell*R - A*R
          PetscReal e_shell, nAR;
          PetscCall(VecNorm(MsR, NORM_2, &e_shell));
          PetscCall(VecNorm(AR, NORM_2, &nAR));
          shell_rel = (nAR > 0) ? e_shell / nAR : e_shell;
          PetscCall(VecDestroy(&sR));
          PetscCall(VecDestroy(&AR));
          PetscCall(VecDestroy(&MsR));
        }
      }
      PetscCall(PetscPrintf(
          PETSC_COMM_WORLD,
          "  [spinDir] ||R||=%.3e ||F_SPIN||=%.3e | Y:||%.3e|| exact||%.3e|| cos=%.4f | "
          "MFspin-R=%.3e MshR-AR=%.3e | outerKSP its=%d reason=%d | fullJ(refG) its=%d reason=%d\n",
          (double)gradPsiNorm,
          (double)fnorm,
          (double)ynorm_in,
          (double)gnorm_ex,
          (double)cosYG,
          (double)mfr_rel,
          (double)shell_rel,
          (int)kits_out,
          (int)kr_out,
          (int)kits_full,
          (int)kr_full));
    }

    if (ex->_use_trust_region && ex->_tr_steihaug)
    {
      // ========== SPIN-preconditioned Steihaug-Toint TRUST REGION (full-space) ==========
      // Solve min g^T p + 1/2 p^T A p s.t. ||p||_P <= Delta by preconditioned truncated CG (see
      // the SPIN manuscript, spin repo docs/alg/mspin_tr_pdas), accept on the energy-reduction ratio,
      // adapt the radius. Unlike
      // the 1-D variant, CG explores the whole Krylov subspace, so it finds a descent step even when
      // the single SPIN direction is nearly orthogonal to grad(Psi) (the propagation regime), and it
      // handles the indefinite (softening) Hessian via a negative-curvature-to-boundary exit.
      const PetscReal eta1 = 0.1, eta2 = 0.75, gshrink = 0.25, gexpand = 2.0, Dmax = 1.0e8;
      // Eisenstat-Walker inner forcing (tightened): CG explores the Krylov subspace rather than
      // stopping after one preconditioned step.
      const PetscReal eps =
          PetscMin(0.1, PetscSqrtReal(gradPsiNorm / PetscMax(ex->_r0_norm, 1e-300)));

      // Near-convergence threshold (couple of x the outer atol). There the quadratic energy model
      // breaks down -- pred (~||p||^2) drops below the energy's cubic term (~||p||^3), so the ratio
      // test is unreliable (NOT noise: refuted; genuine higher-order/local non-convexity). The SGS
      // step is still a valid (Newton-like) descent direction, so we trust it there and let the outer
      // solve converge on ||grad Psi||. Away from convergence the ratio test governs fully.
      PetscReal atol_o, rtol_o, stol_o;
      PetscInt maxit_o, maxf_o;
      PetscCall(SNESGetTolerances(snes, &atol_o, &rtol_o, &stol_o, &maxit_o, &maxf_o));
      const PetscReal conv_thresh = 10.0 * PetscMax(atol_o, rtol_o * ex->_r0_norm);

      // First use: Delta0 = ||P^{-1}g||_P = sqrt(g^T P^{-1} g) (natural field-split step scale).
      if (ex->_tr_radius <= 0.0)
      {
        PetscCall(ex->applyBlockSGS(ex->_r_plain, ex->_cg_y)); // y = P^{-1} g
        PetscScalar gy;
        PetscCall(VecDot(ex->_r_plain, ex->_cg_y, &gy));
        ex->_tr_radius = PetscSqrtReal(PetscMax(PetscRealPart(gy), 1e-300));
      }

      bool tr_accepted = false, on_bnd = false, accepted_via_kkt = false;
      PetscReal rho = 0.0, pred = 0.0;
      PetscInt cg_its = 0;
      // Radius floor for FAILURE: once the radius is a tiny fraction of its initial value and the
      // ratio is still bad, the model cannot guide a step here (abrupt softening at nucleation) ->
      // fail so the driver cuts dt, rather than accepting meaningless (near-zero) steps forever.
      const PetscReal rad_fail = 1e-8 * ex->_tr_radius;
      // Decoupled acceptance (bounds + energy merit): stash the reduced residual R_I(X) so the
      // KKT-progress test can evaluate ||R_I(W)|| per trial (which overwrites _r_plain) and then
      // restore R_I(X) for the next TRS re-solve. ||R_I(X)|| itself is the scalar gradPsiNorm.
      if (ex->_bounds && ex->_has_energy_pp)
      {
        if (!ex->_r_stash)
          PetscCall(VecDuplicate(ex->_r_plain, &ex->_r_stash));
        PetscCall(VecCopy(ex->_r_plain, ex->_r_stash));
      }
      for (PetscInt tr_it = 0; tr_it <= max_it; ++tr_it)
      {
        // Residual-merit TR: restore R(X) into _r_plain for this (re-)solve -- the previous trial's
        // computeMerit(W) overwrote it, and steihaugTRS re-reads _r_plain=R(X).
        if (!ex->_has_energy_pp && ex->_r_base)
          PetscCall(VecCopy(ex->_r_base, ex->_r_plain));
        PetscCall(
            ex->steihaugTRS(ex->_r_plain, ex->_tr_radius, eps, pred, on_bnd, cg_its)); // step -> _cg_p
        PetscCall(VecWAXPY(W, 1.0, ex->_cg_p, X)); // W = X + p  (p solves A p = -g, a descent step)
        // Feasibility projection of the trial step: clip d onto [d_old, 1] so a bound-constrained
        // step cannot overshoot (the ratio test then sees the FEASIBLE energy/residual).
        ex->pdasProjectTrial(W);
        // Residual merit: override the energy-model pred (from steihaugTRS) with the residual-model
        // reduction 1/2(||R(X)||^2 - ||R(X)+A p||^2), formed BEFORE the trial residual eval overwrites
        // _r_plain. A p = _cg_Ad (left by steihaugTRS); R(X) = _r_base; G is free scratch here.
        if (!ex->_has_energy_pp)
        {
          PetscCall(VecWAXPY(G, 1.0, ex->_cg_Ad, ex->_r_base)); // G = R(X) + A p
          PetscReal rApnorm;
          PetscCall(VecNorm(G, NORM_2, &rApnorm));
          pred = 0.5 * (gradPsiNorm * gradPsiNorm - rApnorm * rApnorm);
        }
        const PetscReal psiW = ex->computeMerit(W); // Psi(W) (energy) or 1/2||R(W)||^2 (residual)
        const PetscReal ared = psiX - psiW;
        rho = (pred > 0.0 && !PetscIsInfOrNanReal(ared)) ? ared / pred : -1.0;
        // Accept on a good ratio, OR when we are in the local (Newton) basin: once ||grad Psi|| has
        // dropped ~2 orders from this step's start (or is within a few x the outer atol), the
        // quadratic ENERGY model is unreliable -- pred sinks below the energy's higher-order term, so
        // the ratio goes erratic -- but the strong SGS step is a valid (near-Newton) descent
        // direction, so trust it and let the outer solve converge on ||grad Psi||. The threshold is
        // RELATIVE (load-independent): an absolute one lags the breakdown floor, which rises with
        // load. Away from the basin (e.g. the residual JUMP at nucleation) the ratio test governs
        // fully, so a genuinely bad model (abrupt softening) still rejects -> shrink -> cutback.
        // near_conv trusts the SGS step without the ratio test when close to a solution. For the
        // bound-constrained path this must NOT fire while the active set is still changing (the
        // eliminated residual can look tiny mid-transition), or it accepts overshooting steps.
        const bool near_conv =
            (!ex->_bounds || !ex->_active_set_changed) &&
            (gradPsiNorm < conv_thresh || gradPsiNorm < 1e-2 * PetscMax(ex->_r0_norm, 1e-300));
        // Decoupled KKT-progress clause (bounds): the energy Psi drives the smooth descent via rho, but
        // near the bound-constrained minimizer Psi is flat while the reduced (inactive) residual still
        // resolves the complementarity/active set. Accept if ||R_I(W)|| < ||R_I(X)|| = gradPsiNorm even
        // when rho abstains, so the active-set resolution is not vetoed by a flat-energy ratio. This
        // ||R_I||-decrease is itself the overshoot safeguard the old !active_set_changed gate provided.
        // ||R_I(W)|| via reducedResidualMerit(W): energy path clobbers _r_plain -> restore from stash.
        bool kkt_progress = false;
        if (ex->_bounds)
        {
          PetscReal kktW;
          if (ex->_has_energy_pp)
          {
            const PetscReal mW = ex->reducedResidualMerit(W); // 1/2||R_I(W)||^2 (overwrites _r_plain)
            kktW = PetscSqrtReal(2.0 * PetscMax(mW, 0.0));
            PetscCall(VecCopy(ex->_r_stash, ex->_r_plain)); // restore reduced R(X) for the next TRS
          }
          else
            kktW = PetscSqrtReal(2.0 * PetscMax(psiW, 0.0)); // residual merit: psiW = 1/2||R_I(W)||^2
          kkt_progress = (kktW < gradPsiNorm);
        }
        const bool accept = (rho >= eta1 || near_conv || kkt_progress);
        accepted_via_kkt = accept && rho < eta1 && !near_conv; // accepted purely on KKT-residual progress
        // Per-TRIAL trace: one line per Steihaug trial, so rejected trials are visible and the residual/
        // energy eval count for this outer iteration = number of these lines. acc=1 marks the accepted one.
        if (ex->_verbose)
        {
          PetscReal pnorm_t;
          PetscCall(VecNorm(ex->_cg_p, NORM_2, &pnorm_t));
          PetscCall(PetscPrintf(
              PETSC_COMM_WORLD,
              "    [spinTRcg] rho=%.3e pred=%.3e radius=%.3e ||p||=%.3e cg=%d onBnd=%d "
              "||gradPsi||=%.6e kkt=%d acc=%d\n",
              (double)rho, (double)pred, (double)ex->_tr_radius, (double)pnorm_t,
              (int)cg_its, (int)on_bnd, (double)gradPsiNorm, (int)accepted_via_kkt, (int)accept));
        }
        if (accept)
        {
          if (rho > eta2 && on_bnd) // very good AND radius-limited -> expand
            ex->_tr_radius = PetscMin(gexpand * ex->_tr_radius, Dmax);
          tr_accepted = true;
          break;
        }
        // Reject: shrink the radius and re-solve. Log the rejection + the shrink on its own line.
        const PetscReal rad_before = ex->_tr_radius;
        ex->_tr_radius *= gshrink;
        const bool tr_fail = (ex->_tr_radius < rad_fail || ex->_tr_radius < 1e-14);
        if (ex->_verbose)
          PetscCall(PetscPrintf(PETSC_COMM_WORLD,
                                "    [spinTR reject] rho=%.3e < eta1=%.2f  radius %.3e -> %.3e%s\n",
                                (double)rho,
                                (double)eta1,
                                (double)rad_before,
                                (double)ex->_tr_radius,
                                tr_fail ? "  (radius floor: fail -> dt cutback)" : ""));
        if (tr_fail)
          break; // radius collapsed far from convergence -> fail this solve (dt cutback, not a hang)
      }

      PetscReal pnorm;
      PetscCall(VecNorm(ex->_cg_p, NORM_2, &pnorm)); // ||p|| for SNESLineSearchSetNorms below
      // (per-trial [spinTRcg] and [spinTR reject] are logged inside the retry loop above)

      if (tr_accepted)
      {
        // Commit X <- X + p. In the Steihaug path the outer SNES function F_SPIN and its Krylov
        // direction Y are NOT used (convergence is on ||grad Psi|| via outerConvergenceTest; the step
        // comes from the TRS; the block-Jacobi preconditioner, not the multiplicative NPC, is the
        // preconditioner). So we deliberately do NOT run SNESApplyNPC here -- that extra
        // multiplicative sweep is a full inner disp+pf Newton solve (the pf one oscillating on the
        // penalty kink) whose result is discarded. Report ||grad Psi|| as the line-search fnorm for
        // the monitor; F is stale but only ever feeds the discarded Krylov RHS.
        PetscReal xnorm;
        PetscCall(VecCopy(W, X));
        PetscCall(VecNorm(X, NORM_2, &xnorm));
        // Record the accepted d-increment (delta d^k) for the next iteration's active-set criterion.
        ex->pdasRecordStep(ex->_cg_p);
        PetscCall(SNESLineSearchSetNorms(ls, xnorm, gradPsiNorm, pnorm));
        PetscCall(SNESLineSearchSetLambda(ls, 1.0));
        PetscCall(SNESLineSearchSetReason(ls, SNES_LINESEARCH_SUCCEEDED));
      }
      else
        PetscCall(SNESLineSearchSetReason(ls, SNES_LINESEARCH_FAILED_REDUCT));

      PetscFunctionReturn(PETSC_SUCCESS);
    }

    if (ex->_use_trust_region)
    {
      // ================= ratio-based TRUST REGION on the energy Psi =================
      // 1-D trust region restricted to the SPIN direction Y. Quadratic energy model along the ray
      // X - alpha*D:  m(alpha) = Psi - alpha*(g.D) + 1/2*alpha^2*(D^T A D), with g = grad Psi = R,
      // A = _mat_nest (the coupled Hessian at X, assembled by outerJacobianCallback). The step is the
      // model minimizer clipped to the radius; accept/shrink on the actual/predicted reduction ratio;
      // adapt the radius. If Y is not a descent direction for Psi (g.Y <= 0, near instability where A
      // is indefinite), fall back to steepest descent D = g (self-scaled by the model, Cauchy-like).
      // For the exact SPIN direction Y ~ A^{-1}g the model minimizer is alpha ~ 1 (the full SPIN =
      // energy-Newton step); when the model is poor (nonlinear/propagation region) the ratio shrinks
      // the step, so the solve makes small monotone-energy progress instead of spiraling dt to dtmin.
      const PetscReal eta1 = 0.1;      // accept if rho >= eta1
      const PetscReal eta2 = 0.75;     // expand radius if rho > eta2 and step hit the boundary
      const PetscReal gshrink = 0.25;  // radius shrink factor on rejection
      const PetscReal gexpand = 2.0;   // radius expansion factor
      const PetscReal Dmax = 1.0e8;    // max radius (solution norm)

      // Outer convergence tolerances, for the near-convergence bypass below.
      PetscReal atol_o, rtol_o, stol_o;
      PetscInt maxit_o, maxf_o;
      PetscCall(SNESGetTolerances(snes, &atol_o, &rtol_o, &stol_o, &maxit_o, &maxf_o));
      const PetscReal conv_thresh = 10.0 * PetscMax(atol_o, rtol_o * ex->_r0_norm);

      Vec D = Y;                 // step direction (default: SPIN direction)
      PetscReal dnorm = ynorm;   // ||D||
      PetscScalar gDs;
      PetscCall(VecDot(ex->_r_plain, Y, &gDs));
      PetscReal gD = PetscRealPart(gDs); // g.D  (>0 => -D is a descent direction for Psi)

      // Near-convergence last mile: Psi is flat/noisy and the SPIN direction is ~ the energy-Newton
      // direction (cos ~ 1), so the ratio test is dominated by energy-eval noise (the SPIN step
      // reduces ||F_SPIN|| while nudging Psi within noise -> a meaningless negative ratio that would
      // collapse the radius). Take the full SPIN step; the outer solve converges on ||gradPsi||,
      // which the SPIN step drives to zero. In the HARD region ||gradPsi|| >> conv_thresh, so the
      // ratio-test globalization below applies fully.
      if (gradPsiNorm < conv_thresh && gD > 0.0)
      {
        PetscReal xnorm, fnorm_new;
        PetscCall(VecWAXPY(W, -1.0, Y, X)); // W = X - Y (full SPIN step)
        PetscCall(VecCopy(W, X));
        PetscCall(SNESApplyNPC(snes, W, nullptr, F));
        PetscCall(VecNorm(X, NORM_2, &xnorm));
        PetscCall(VecNorm(F, NORM_2, &fnorm_new));
        PetscCall(SNESLineSearchSetNorms(ls, xnorm, fnorm_new, ynorm));
        PetscCall(SNESLineSearchSetLambda(ls, 1.0));
        PetscCall(SNESLineSearchSetReason(ls, SNES_LINESEARCH_SUCCEEDED));
        PetscCall(PetscPrintf(PETSC_COMM_WORLD,
                              "  [spinTR] near-conv full SPIN step ||gradPsi||=%.6e (<%.3e)\n",
                              (double)gradPsiNorm,
                              (double)conv_thresh));
        PetscFunctionReturn(PETSC_SUCCESS);
      }

      if (gD <= 0.0)
      {
        D = ex->_r_plain;               // steepest descent fallback (X - alpha*g reduces Psi)
        gD = gradPsiNorm * gradPsiNorm; // g.D = ||g||^2 > 0
        dnorm = gradPsiNorm;
      }
      PetscCall(MatMult(ex->_mat_nest, D, G)); // G = A D
      PetscScalar cs;
      PetscCall(VecDot(D, G, &cs));
      const PetscReal curv = PetscRealPart(cs); // D^T A D (model curvature)

      if (ex->_tr_radius <= 0.0)
        ex->_tr_radius = (dnorm > 0.0) ? dnorm : 1.0; // first use: allow ~one full step (alpha_max~1)

      bool tr_accepted = false;
      PetscReal alpha = 0.0, psiW_tr = psiX, rho = 0.0;
      for (PetscInt it = 0; it <= max_it; ++it)
      {
        const PetscReal amax = ex->_tr_radius / dnorm;
        // 1-D model minimizer of m(alpha) clipped to [0, amax]
        PetscReal a = (curv > 0.0) ? PetscMin(gD / curv, amax) : amax;
        if (a <= 0.0)
          a = amax;
        PetscCall(VecWAXPY(W, -a, D, X)); // W = X - a D
        psiW_tr = ex->computeEnergy(W);
        const PetscReal pred = a * gD - 0.5 * a * a * curv; // model reduction m(0)-m(a) (>0)
        const PetscReal ared = psiX - psiW_tr;              // actual reduction
        rho = (pred > 0.0 && !PetscIsInfOrNanReal(ared)) ? ared / pred : -1.0;
        if (rho >= eta1)
        {
          alpha = a;
          if (rho > eta2 && a >= amax * (1.0 - 1e-10)) // very good AND at the boundary -> expand
            ex->_tr_radius = PetscMin(gexpand * ex->_tr_radius, Dmax);
          tr_accepted = true;
          break;
        }
        // reject: shrink the radius below the failed step and retry the same direction
        ex->_tr_radius = gshrink * a * dnorm;
        if (ex->_tr_radius < minlambda * dnorm)
          break; // radius collapsed -> fail this outer solve (triggers a dt cutback)
      }

      PetscCall(PetscPrintf(
          PETSC_COMM_WORLD,
          "  [spinTR] rho=%.3e alpha=%.3e radius=%.3e curv=%.3e ||gradPsi||=%.6e accepted=%d\n",
          (double)rho,
          (double)alpha,
          (double)ex->_tr_radius,
          (double)curv,
          (double)gradPsiNorm,
          (int)tr_accepted));

      if (tr_accepted)
      {
        PetscReal xnorm, fnorm_new;
        PetscCall(VecCopy(W, X));                        // X <- X - alpha*D
        PetscCall(SNESApplyNPC(snes, W, nullptr, F));    // F <- F_SPIN(X) for the convergence test
        PetscCall(VecNorm(X, NORM_2, &xnorm));
        PetscCall(VecNorm(F, NORM_2, &fnorm_new));
        PetscCall(SNESLineSearchSetNorms(ls, xnorm, fnorm_new, alpha * dnorm));
        PetscCall(SNESLineSearchSetLambda(ls, alpha));
        PetscCall(SNESLineSearchSetReason(ls, SNES_LINESEARCH_SUCCEEDED));
      }
      else
        PetscCall(SNESLineSearchSetReason(ls, SNES_LINESEARCH_FAILED_REDUCT));

      PetscFunctionReturn(PETSC_SUCCESS);
    }

    PetscScalar dot;
    PetscReal slope, lambda = 0.0, psiW = psiX;
    int tier = 1;
    bool accepted = false;

    // --- App. C direction selection: try the SPIN direction first, fall back only if it GENUINELY
    // increases the energy. This is the paper's App. C, but the "is the SPIN direction acceptable"
    // test is an actual (noise-cushioned) energy check rather than the raw sign of slope = -R.Y --
    // which near convergence is dominated by rounding noise (the SPIN direction becomes nearly
    // orthogonal to grad Psi, |cos(R,Y)| ~ 1e-6, since it targets F_SPIN=0 where the energy is already
    // ~minimal). A raw-sign test escalates on that noise and abandons the excellent SPIN direction
    // (which drives ||F_precond|| to machine precision) for the full-J Newton, which stalls when the
    // energy is this flat.
    const PetscReal energy_tol = 1e-8 * PetscAbsReal(psiX); // energy sufficient-decrease cushion (see energyBacktrack)
    PetscCall(VecDot(ex->_r_plain, Y, &dot));                // Y = SPIN direction from the outer KSP
    slope = -PetscRealPart(dot);                             // phi'(0) = -R . Y
    if (slope < 0.0) // SPIN is a descent direction for Psi -> cushioned energy-Armijo backtrack
      accepted = ex->energyBacktrack(X, Y, W, psiX, slope, minlambda, max_it, lambda, psiW);
    else // SPIN is (noise-level or genuine) ascent: try the full SPIN step, accept if energy is flat
    {
      PetscCall(VecWAXPY(W, -1.0, Y, X)); // W = X - Y
      psiW = ex->computeEnergy(W);
      if (!PetscIsInfOrNanReal(psiW) && psiW <= psiX + energy_tol)
      {
        lambda = 1.0;
        accepted = true; // energy flat within noise -> the SPIN Newton step still reduces ||F_precond||
      }
    }

    // --- App. C fallback: only if the SPIN direction genuinely increases the energy. ---
    if (!accepted)
    {
      // full-J inexact Newton J^{-1} R (the energy's Newton step), globalized by strong Wolfe.
      PetscCall(KSPSetOperators(ex->_fullJ_ksp, ex->_mat_nest, ex->_mat_nest));
      PetscCall(KSPSolve(ex->_fullJ_ksp, ex->_r_plain, Y)); // Y = J^{-1} R
      PetscCall(VecDot(ex->_r_plain, Y, &dot));
      slope = -PetscRealPart(dot);
      tier = 2;
      if (slope >= 0.0) // Newton direction also ascent -> Cauchy-scaled steepest descent (tier 3)
      {
        PetscScalar RR, RJR;
        PetscCall(VecDot(ex->_r_plain, ex->_r_plain, &RR));
        PetscCall(MatMult(ex->_mat_nest, ex->_r_plain, W)); // W = J R  (scratch; LS resets it)
        PetscCall(VecDot(ex->_r_plain, W, &RJR));
        // Cauchy step length (1-D minimizer of the quadratic model along -R): keeps alpha=1 well
        // scaled, since raw grad Psi = R is in force units, not solution units.
        const PetscReal aC =
            (PetscRealPart(RJR) > 0.0) ? PetscRealPart(RR) / PetscRealPart(RJR) : 1.0;
        PetscCall(VecCopy(ex->_r_plain, Y));
        PetscCall(VecScale(Y, aC)); // Y = aC * R
        PetscCall(VecDot(ex->_r_plain, Y, &dot));
        slope = -PetscRealPart(dot); // = -aC*||R||^2 < 0
        tier = 3;
      }
      accepted = ex->strongWolfe(X, Y, W, psiX, slope, minlambda, max_it, lambda, psiW);
    }

    // Outer-Krylov iteration count that produced the SPIN direction Y (J_SPIN Y = F_SPIN via applyBA).
    // This is the per-outer-iteration SPIN-direction Krylov cost; sum over the run for the benchmark.
    // (Cheap: reads the count PETSc already computed before the line search; no extra solve.)
    KSP outer_ksp;
    PetscInt kits_out = 0;
    PetscCall(SNESGetKSP(snes, &outer_ksp));
    PetscCall(KSPGetIterationNumber(outer_ksp, &kits_out));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
                          "  [spinLS] tier=%d alpha=%.4e slope=%.4e gradPsi=%.6e outerKSP=%d accepted=%d\n",
                          tier,
                          (double)lambda,
                          (double)slope,
                          (double)gradPsiNorm,
                          (int)kits_out,
                          (int)accepted));

    if (accepted)
    {
      // Commit: outer iterate X <- W; set F <- F_SPIN(W) for the outer convergence test and next
      // Jacobian. SNESApplyNPC also leaves the libMesh state at NPC(W) -- where
      // outerJacobianCallback linearizes -- so no extra bookkeeping is needed.
      PetscReal xnorm, fnorm_new;
      PetscCall(VecNorm(Y, NORM_2, &ynorm));
      PetscCall(VecCopy(W, X));
      PetscCall(SNESApplyNPC(snes, W, nullptr, F));
      PetscCall(VecNorm(X, NORM_2, &xnorm));
      PetscCall(VecNorm(F, NORM_2, &fnorm_new));
      PetscCall(SNESLineSearchSetNorms(ls, xnorm, fnorm_new, ynorm));
      PetscCall(SNESLineSearchSetLambda(ls, lambda));
      PetscCall(SNESLineSearchSetReason(ls, SNES_LINESEARCH_SUCCEEDED));
    }
    else
      PetscCall(SNESLineSearchSetReason(ls, SNES_LINESEARCH_FAILED_REDUCT));

    PetscFunctionReturn(PETSC_SUCCESS);
  }

  // ================= residual-norm merit (default: 1/2||x - NPC(x)||^2) =================
  const PetscReal alpha = 1e-4;               // App. B sufficient-decrease constant (paper c1)
  const PetscReal f = 0.5 * fnorm * fnorm;    // merit at the base iterate

  // Initial slope phi'(0) = -F^T (J_SPIN Y), using the shell operator that generated Y (W is free
  // scratch until the backtracking loop overwrites it). For an accurate SPIN direction J_SPIN Y = F
  // so this is -||F||^2 < 0 (a descent direction for the merit); force-negative as a safety net if
  // an inexact/indefinite operator makes it non-negative.
  PetscCall(MatMult(ex->_jac_shell, Y, W));
  PetscScalar FtJY;
  PetscCall(VecDot(F, W, &FtJY));
  PetscReal initslope = -PetscRealPart(FtJY);
  if (initslope >= 0.0)
    initslope = -(PetscAbsReal(PetscRealPart(FtJY)) + PETSC_MACHINE_EPSILON);

  PetscReal lambda = 1.0, lambdaprev = 1.0, gprev = f, gnorm = fnorm, g = f;
  PetscBool success = PETSC_FALSE;
  for (PetscInt it = 0; it <= max_it; ++it)
  {
    // Trial point W = X - lambda*Y and CONSISTENT merit via the preconditioned residual:
    // G = F_SPIN(W) = W - NPC(W). This is the crux -- each trial costs one inner NPC sub-solve.
    PetscCall(VecWAXPY(W, -lambda, Y, X));
    PetscCall(SNESApplyNPC(snes, W, nullptr, G));
    PetscCall(VecNorm(G, NORM_2, &gnorm));
    g = 0.5 * gnorm * gnorm;

    if (PetscIsInfOrNanReal(g))
    {
      // Non-finite merit: shrink hard and retry (do not attempt interpolation).
      if (lambda <= minlambda)
        break;
      lambdaprev = lambda;
      gprev = g;
      lambda = PetscMax(0.1 * lambda, minlambda);
      continue;
    }

    if (g <= f + lambda * alpha * initslope) // Armijo sufficient decrease
    {
      success = PETSC_TRUE;
      break;
    }

    if (lambda <= minlambda)
      break;

    // Interpolate a shorter step (quadratic on the first backtrack, cubic thereafter), mirroring
    // PETSc's SNESLineSearchApply_BT, then clip to [0.1, 0.5]*lambda.
    PetscReal lambdatemp;
    if (it == 0)
      lambdatemp = -initslope * lambda * lambda / (2.0 * (g - f - lambda * initslope));
    else
    {
      const PetscReal t1 = g - f - lambda * initslope;
      const PetscReal t2 = gprev - f - lambdaprev * initslope;
      const PetscReal a =
          (t1 / (lambda * lambda) - t2 / (lambdaprev * lambdaprev)) / (lambda - lambdaprev);
      const PetscReal b =
          (-lambdaprev * t1 / (lambda * lambda) + lambda * t2 / (lambdaprev * lambdaprev)) /
          (lambda - lambdaprev);
      PetscReal d = b * b - 3.0 * a * initslope;
      if (d < 0.0)
        d = 0.0;
      if (a == 0.0)
        lambdatemp = -initslope / (2.0 * b);
      else
        lambdatemp = (-b + PetscSqrtReal(d)) / (3.0 * a);
    }
    if (PetscIsInfOrNanReal(lambdatemp))
      lambdatemp = 0.5 * lambda;

    lambdaprev = lambda;
    gprev = g;
    lambda = PetscClipInterval(lambdatemp, 0.1 * lambda, 0.5 * lambda);
  }

  if (success)
  {
    // Commit the accepted trial: X <- W, F <- F_SPIN(W); set consistent norms/lambda/reason.
    PetscReal xnorm;
    PetscCall(VecCopy(W, X));
    PetscCall(VecCopy(G, F));
    PetscCall(VecNorm(X, NORM_2, &xnorm));
    PetscCall(SNESLineSearchSetNorms(ls, xnorm, gnorm, ynorm));
    PetscCall(SNESLineSearchSetLambda(ls, lambda));
    PetscCall(SNESLineSearchSetReason(ls, SNES_LINESEARCH_SUCCEEDED));
  }
  else
    PetscCall(SNESLineSearchSetReason(ls, SNES_LINESEARCH_FAILED_REDUCT));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
NewtonSNESExecutor::applyBlockJacobi(Vec in, Vec out)
{
  PetscFunctionBegin;
  Vec in0, in1, out0, out1;
  PetscCall(VecNestGetSubVec(in, 0, &in0));
  PetscCall(VecNestGetSubVec(in, 1, &in1));
  PetscCall(VecNestGetSubVec(out, 0, &out0));
  PetscCall(VecNestGetSubVec(out, 1, &out1));
  PetscCall(KSPSolve(_pc_ksp0, in0, out0)); // out0 = A00^{-1} in0
  PetscCall(KSPSolve(_pc_ksp1, in1, out1)); // out1 = A11^{-1} in1
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
NewtonSNESExecutor::applyBlockSGS(Vec in, Vec out)
{
  // out = P_SGS^{-1} in = (D+U)^{-1} D (D+L)^{-1} in, with D=diag(A00,A11), L=[[0,0],[A10,0]], U=L^T.
  //   forward  : z0 = A00^{-1} in0;      z1 = A11^{-1} (in1 - A10 z0)   -> out
  //   scale    : w0 = A00 z0;            w1 = A11 z1                    -> _cg_Pd (scratch)
  //   backward : y1 = A11^{-1} w1;       y0 = A00^{-1} (w0 - A01 y1)    -> out
  PetscFunctionBegin;
  Mat A00, A01, A10, A11;
  PetscCall(MatNestGetSubMat(_mat_nest, 0, 0, &A00));
  PetscCall(MatNestGetSubMat(_mat_nest, 0, 1, &A01));
  PetscCall(MatNestGetSubMat(_mat_nest, 1, 0, &A10));
  PetscCall(MatNestGetSubMat(_mat_nest, 1, 1, &A11));
  Vec in0, in1, out0, out1, t0, t1;
  PetscCall(VecNestGetSubVec(in, 0, &in0));
  PetscCall(VecNestGetSubVec(in, 1, &in1));
  PetscCall(VecNestGetSubVec(out, 0, &out0));
  PetscCall(VecNestGetSubVec(out, 1, &out1));
  PetscCall(VecNestGetSubVec(_cg_Pd, 0, &t0)); // scratch (w0, then reused)
  PetscCall(VecNestGetSubVec(_cg_Pd, 1, &t1)); // scratch (in1 - A10 z0, w1, A01 y1)

  // forward sweep (D+L)^{-1}
  PetscCall(KSPSolve(_pc_ksp0, in0, out0));  // z0 = A00^{-1} in0
  PetscCall(MatMult(A10, out0, t1));         // t1 = A10 z0
  PetscCall(VecAYPX(t1, -1.0, in1));         // t1 = in1 - A10 z0
  PetscCall(KSPSolve(_pc_ksp1, t1, out1));   // z1 = A11^{-1} (in1 - A10 z0)
  // block-diagonal scale D
  PetscCall(MatMult(A00, out0, t0));         // t0 = w0 = A00 z0 (disp)
  PetscCall(MatMult(A11, out1, t1));         // t1 = w1 = A11 z1 (pf)
  // backward sweep (D+U)^{-1}. A01 maps pf->disp, so A01*y1 is disp-sized: reuse out0 (z0 is done).
  PetscCall(KSPSolve(_pc_ksp1, t1, out1));   // y1 = A11^{-1} w1 (pf)
  PetscCall(MatMult(A01, out1, out0));       // out0 = A01 y1 (disp)
  PetscCall(VecAYPX(out0, -1.0, t0));        // out0 = w0 - A01 y1 (disp)
  PetscCall(KSPSolve(_pc_ksp0, out0, t0));   // t0 = A00^{-1} (w0 - A01 y1) = y0
  PetscCall(VecCopy(t0, out0));              // out0 = y0 (disp)
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
NewtonSNESExecutor::steihaugTRS(
    Vec g, Real Delta, Real eps, Real & pred, bool & on_boundary, PetscInt & cg_its)
{
  // Preconditioned Steihaug-Toint truncated CG on A p = -g within ||p||_P <= Delta, P = symmetric
  // block Gauss-Seidel. The P-norm quantities (||p||_P^2, <p,d>_P, ||d||_P^2) are tracked by scalar
  // recurrences (no forward P-apply needed) using the CG orthogonalities r_j _|_ d_i, r_j _|_ y_i
  // (i<j): ||d_j||_P^2 = zeta_j + beta_j^2 ||d_{j-1}||_P^2, <p_j,d_j>_P = beta_j(<p_{j-1},d_{j-1}>_P
  // + alpha_{j-1}||d_{j-1}||_P^2), with zeta_j = r_j^T P^{-1} r_j. Step returned in _cg_p.
  // See the SPIN manuscript: spin repo docs/alg/mspin_tr_pdas.
  PetscFunctionBegin;
  // Bind the block preconditioner to the current diagonal blocks (buildMatNest recreates _mat_nest).
  Mat A00, A11;
  PetscCall(MatNestGetSubMat(_mat_nest, 0, 0, &A00));
  PetscCall(MatNestGetSubMat(_mat_nest, 1, 1, &A11));
  PetscCall(KSPSetOperators(_pc_ksp0, A00, A00));
  PetscCall(KSPSetOperators(_pc_ksp1, A11, A11));

  Vec p = _cg_p, r = _cg_r, y = _cg_y, d = _cg_d, Ad = _cg_Ad;
  const PetscReal Delta2 = Delta * Delta;

  PetscCall(VecSet(p, 0.0));
  PetscCall(VecCopy(g, r));          // r = A p + g = g at p = 0
  PetscCall(applyBlockSGS(r, y));    // y = P^{-1} r
  PetscScalar ry_s;
  PetscCall(VecDot(r, y, &ry_s));    // zeta_0 = r^T P^{-1} r  (>= 0)
  PetscReal ry = PetscRealPart(ry_s);
  PetscReal gnorm;
  PetscCall(VecNorm(g, NORM_2, &gnorm));
  PetscCall(VecCopy(y, d));
  PetscCall(VecScale(d, -1.0));      // d = -y

  on_boundary = false;
  cg_its = 0;
  PetscReal np = 0.0;  // ||p||_P^2
  PetscReal nd = ry;   // ||d||_P^2 (= zeta_0 for d_0 = -y_0)
  PetscReal pdip = 0.0; // <p,d>_P (=0 at p=0)
  const PetscInt maxcg = 100;
  for (PetscInt j = 0; j < maxcg; ++j)
  {
    cg_its = j + 1;
    PetscCall(MatMult(_mat_nest, d, Ad)); // Ad = A d
    PetscScalar dAd_s;
    PetscCall(VecDot(d, Ad, &dAd_s));
    const PetscReal kappa = PetscRealPart(dAd_s); // curvature d^T A d

    if (kappa <= 0.0) // (a) negative curvature: model unbounded along d -> go to the boundary
    {
      const PetscReal tau =
          (nd > 0.0) ? (-pdip + PetscSqrtReal(PetscMax(pdip * pdip + nd * (Delta2 - np), 0.0))) / nd
                     : 0.0;
      PetscCall(VecAXPY(p, tau, d));
      on_boundary = true;
      break;
    }

    const PetscReal alpha = ry / kappa;
    const PetscReal npnew = np + 2.0 * alpha * pdip + alpha * alpha * nd; // ||p + alpha d||_P^2
    if (npnew >= Delta2) // (b) boundary hit
    {
      const PetscReal tau =
          (nd > 0.0) ? (-pdip + PetscSqrtReal(PetscMax(pdip * pdip + nd * (Delta2 - np), 0.0))) / nd
                     : 0.0;
      PetscCall(VecAXPY(p, tau, d));
      on_boundary = true;
      break;
    }

    PetscCall(VecAXPY(p, alpha, d)); // p += alpha d
    np = npnew;
    PetscCall(VecAXPY(r, alpha, Ad)); // r += alpha A d
    PetscReal rnorm;
    PetscCall(VecNorm(r, NORM_2, &rnorm));
    if (rnorm <= eps * gnorm) // (c) interior convergence
      break;
    PetscCall(applyBlockSGS(r, y)); // y = P^{-1} r
    PetscScalar ry_new_s;
    PetscCall(VecDot(r, y, &ry_new_s));
    const PetscReal ry_new = PetscRealPart(ry_new_s); // zeta_{j+1}
    const PetscReal beta = ry_new / ry;
    // scalar P-norm recurrences for d_{j+1} = -y_{j+1} + beta d_j (use OLD nd/pdip)
    pdip = beta * (pdip + alpha * nd); // <p_{j+1}, d_{j+1}>_P
    nd = ry_new + beta * beta * nd;    // ||d_{j+1}||_P^2
    ry = ry_new;
    PetscCall(VecScale(d, beta));   // d = beta d
    PetscCall(VecAXPY(d, -1.0, y)); // d = -y + beta d
  }

  // Model reduction pred = -(g^T p + 1/2 p^T A p) > 0 (Steihaug guarantees model descent).
  PetscCall(MatMult(_mat_nest, p, Ad)); // Ad reused as A p
  PetscScalar gp_s, pAp_s;
  PetscCall(VecDot(g, p, &gp_s));
  PetscCall(VecDot(p, Ad, &pAp_s));
  pred = -(PetscRealPart(gp_s) + 0.5 * PetscRealPart(pAp_s));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
NewtonSNESExecutor::assembleCoupledJacobian(Vec x)
{
  KSP sub_ksp;

  PetscFunctionBegin;
  // Assemble the coupled Jacobian _mat_nest = grad^2 Psi at the OUTER ITERATE x = (U^k, C^k),
  // matching the paper (Algorithm 3: J^k <- F'(U^k, C^k)). The residual eval (SNESApplyNPC) leaves
  // the libMesh state at NPC(x), so WITHOUT this scatter computeJacobian would linearize at
  // NPC(x) != x. (In the (B) NPC-every-iteration path the caller passes x = NPC(x) deliberately, so
  // the Hessian and the residual are both taken at the swept point.)
  scatterToSystems(x);
  const unsigned int n_sys_scatter = _fe_problem.numNonlinearSystems();
  for (unsigned int i = 0; i < n_sys_scatter; ++i)
  {
    auto & sys_i = _fe_problem.getNonlinearSystemBase(i);
    sys_i.setSolution(*sys_i.system().current_local_solution);
  }

  for (const auto nl_sys_num : _nl_sys_nums)
  {
    auto & moose_sys = _fe_problem.getNonlinearSystem(nl_sys_num);
    auto & lm_sys = moose_sys.sys();
    auto & J_ii = lm_sys.get_system_matrix();
    // current_local_solution now holds the outer iterate x (scattered above).
    _fe_problem.computeJacobian(*lm_sys.current_local_solution, J_ii, nl_sys_num);
    auto & dof_map = lm_sys.get_dof_map();
    if (dof_map.n_constrained_dofs())
    {
      // Even though our solution was constrained we stupidly don't apply our asymmetric constraints
      // at the element level so the Jacobian will be wrong without this call below
      dof_map.enforce_constraints_on_jacobian(lm_sys, &J_ii);
      J_ii.close();
    }
  }

  assembleOffDiagJacobian();

  PetscCall(MatAssemblyBegin(_mat_nest, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(_mat_nest, MAT_FINAL_ASSEMBLY));

  // Update per-field KSP operators from the freshly assembled diagonal blocks so that
  // shellMatMult (and the Steihaug field-split preconditioner) reuse the current factorizations.
  for (const auto i : index_range(_nl_sys_nums))
  {
    Mat J_ii;
    PetscCall(MatNestGetSubMat(_mat_nest, i, i, &J_ii));
    auto sub_snes = _fe_problem.getNonlinearSystem(_nl_sys_nums[i]).getSNES();
    PetscCall(SNESGetKSP(sub_snes, &sub_ksp));
    PetscCall(KSPSetOperators(sub_ksp, J_ii, J_ii));
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
NewtonSNESExecutor::outerJacobianCallback(SNES /*snes*/, Vec x, Mat /*A*/, Mat /*P*/, void * ctx)
{
  PetscFunctionBegin;
  auto * ex = static_cast<NewtonSNESExecutor *>(ctx);

  // (B) NPC-every-iteration TR (_tr_npc): the Steihaug line search sweeps X <- NPC(X) and RE-assembles
  // grad^2 Psi at NPC(X) itself, so skip this (raw-x) assembly to avoid a redundant Jacobian per
  // outer iteration. All other paths assemble the coupled Jacobian at the outer iterate here.
  if (ex->_use_trust_region && ex->_tr_steihaug && ex->_tr_npc)
    PetscFunctionReturn(PETSC_SUCCESS);

  PetscCall(ex->assembleCoupledJacobian(x));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SNES
NewtonSNESExecutor::getSNES()
{
  if (_nl_sys_nums.size() == 1)
    // We cannot rely on caching this during some setup routine because libMesh destroys its SNES at
    // the end of each nonlinear solve. For similar reasons we cannot have _snes point to this
    // because then at destruction time we'll attempt to destroy through a pointer to garbage since
    // the SNES was already destroyed at the end of the solve
    return _fe_problem.getNonlinearSystem(_nl_sys_nums[0]).getSNES();

  // For multiple systems _snes is our PETSc-owned object; use the base-class lazy-setup path
  return SNESExecutor::getSNES();
}

System &
NewtonSNESExecutor::getSystem()
{
  if (_nl_sys_nums.size() != 1)
    mooseError("Ambiguous call to getSystem()");

  return _fe_problem.getNonlinearSystem(_nl_sys_nums[0]).system();
}
