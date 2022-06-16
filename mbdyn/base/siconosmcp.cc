/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2022
 *
 * Pierangelo Masarati	<masarati@aero.polimi.it>
 * Paolo Mantegazza	<mantegazza@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 *
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

 /*
 AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
        Copyright (C) 2022(-2022) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
  */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#ifdef USE_SICONOS
#include <numerics/SolverOptions.h>
#include <numerics/MCP_cst.h>
#include <numerics/NumericsMatrix.h>
#include <numerics/NonSmoothDrivers.h>
#include <numerics/MixedComplementarityProblem.h>
#include <numerics/Newton_methods.h>

#include "siconoswrap.h"
#include "siconosmcp.h"
#include "solver.h"

SiconosMCPSolver::SiconosMCPSolver(const NonlinearSolverTestOptions& options,
                                   const LineSearchParameters& oLineSearch)
     :NonlinearSolver(options),
      LineSearchParameters(oLineSearch),
      pRes(nullptr),
      pAbsRes(nullptr),
      pSol(nullptr),
      pJac(nullptr),
      pNLP(nullptr),
      pSolver(nullptr),
      pMCP(nullptr),
      pOptions(nullptr),
      iJacPrev(-1),
      iIterCurr(-1)
{
}

SiconosMCPSolver::~SiconosMCPSolver()
{
     Cleanup();
}

void SiconosMCPSolver::Solve(const NonlinearProblem *pNLPCurr,
                             Solver *pS,
                             const integer iMaxIter,
                             const doublereal& Tol,
                             integer& iIterCnt,
                             doublereal& dErr,
                             const doublereal& SolTol,
                             doublereal& dSolErr)
{
     if (pNLP != pNLPCurr) {
          Attach(pNLPCurr, pS);
     }

     pOptions->iparam[SICONOS_IPARAM_MAX_ITER] = iMaxIter;
     pOptions->dparam[SICONOS_DPARAM_TOL] = Tol;

     pSol->Reset();
     W.Reset();
     iJacPrev = TotJac;
     iIterCurr = 0;

     int status = mcp_driver(pMCP, pSol->pdGetVec(), W.pdGetVec(), pOptions);

     iIterCnt = pOptions->iparam[SICONOS_IPARAM_ITER_DONE];
     dErr = pOptions->dparam[SICONOS_DPARAM_RESIDU];

     ASSERT(iIterCnt == iIterCurr);

     if (status != 0) {
          throw NoConvergence(MBDYN_EXCEPT_ARGS);
     }
}

void SiconosMCPSolver::Cleanup()
{
     if (pMCP) {
          pMCP->nabla_Fmcp = nullptr;
          mixedComplementarityProblem_free(pMCP);
     }

     pMCP = nullptr;

     if (pOptions) {
          solver_options_delete(pOptions);
          free(pOptions);
     }

     pOptions = nullptr;
}

void SiconosMCPSolver::Attach(const NonlinearProblem* pNLPCurr, Solver* pS)
{
     Cleanup();

     pNLP = pNLPCurr;
     pSolver = pS;
     pMCP = mixedComplementarityProblem_new();
     pMCP->compute_Fmcp = &compute_Fmcp;
     pMCP->compute_nabla_Fmcp = &compute_nabla_Fmcp;
     pMCP->env = this;
     pOptions = solver_options_create(SICONOS_MCP_NEWTON_MIN_FBLSA);
     pOptions->callback = static_cast<Callback*>(malloc(sizeof(Callback)));
     pOptions->callback->collectStatsIteration = &collectStatsIteration;
     pOptions->callback->env = this;
     pOptions->iparam[SICONOS_IPARAM_STOPPING_CRITERION] = SICONOS_STOPPING_CRITERION_USER_ROUTINE;
     pOptions->dparam[SICONOS_DPARAM_LSA_ALPHA_MIN] = dLambdaMin;

     auto* const pSM = dynamic_cast<SiconosDenseSolutionManager*>(pS->pGetSolutionManager());

     if (!pSM) {
          silent_cerr("For the siconos nonlinear solver the linear solver type must be \"siconos dense\"\n");
          throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
     }

     pRes = pSM->pResHdl();
     pAbsRes = pGetResTest()->GetAbsRes();
     pSol = pSM->pSolHdl();
     pJac = pSM->pMatHdl();

     ASSERT(pJac->pGetMatrix()->storageType == NM_DENSE);

     const DataManager* const pDM = pSolver->pGetDataManager();

     const integer iNumDofsTot = pDM->iGetNumDofs();

     ASSERT(iNumDofsTot == pJac->iGetNumCols());
     ASSERT(iNumDofsTot == pRes->iGetSize());
     ASSERT(iNumDofsTot == pSol->iGetSize());
     ASSERT(pAbsRes == nullptr || iNumDofsTot == pAbsRes->iGetSize());

     std::array<integer, DofOrder::LASTCOMP> rgNumDofs{0};

     for (integer i = 1; i <= iNumDofsTot; ++i) {
          ++rgNumDofs[pDM->GetCompType(i)];
     }

     pMCP->n1 = rgNumDofs[DofOrder::EQUALITY];
     pMCP->n2 = rgNumDofs[DofOrder::COMPLEMENTARY];

     std::vector<integer> rgRowMap(iNumDofsTot);

     std::array<integer, DofOrder::LASTCOMP> rgCurrIndex{0, rgNumDofs[DofOrder::EQUALITY]};

     for (integer i = 1; i <= iNumDofsTot; ++i) {
          rgRowMap[i - 1] = ++rgCurrIndex[pDM->GetCompType(i)];
     }

     pSM->pGetRowMap()->SetIndex(std::move(rgRowMap));

     pMCP->nabla_Fmcp = pJac->pGetMatrix();

     ASSERT(pSol->iGetSize() == pDM->iGetNumDofs());
     ASSERT(pMCP->n1 + pMCP->n2 == pDM->iGetNumDofs());
     ASSERT(pSol->iGetSize() == pJac->iGetNumCols());
     ASSERT(pRes->iGetSize() == pJac->iGetNumRows());

     W.SetRowMap(pSol->pGetRowMap());
     XPrev.SetRowMap(pSol->pGetRowMap());
     DeltaX.SetRowMap(pSol->pGetRowMap());

     W.ResizeReset(pSol->iGetSize());

     if (pSol->iGetSize() != XPrev.iGetSize()) {
          XPrev.ResizeReset(pSol->iGetSize()); // Preserve previous state!
     }

     DeltaX.ResizeReset(pSol->iGetSize());
}

void SiconosMCPSolver::compute_Fmcp(void *env, int n, doublereal *z, doublereal *F)
{
     SiconosMCPSolver* const pSiconosMCP = static_cast<SiconosMCPSolver*>(env);

     ASSERT(n == pSiconosMCP->pSol->iGetSize());

     const SiconosVectorHandler zhdl(n, z, pSiconosMCP->pSol->pGetRowMap());
     SiconosVectorHandler Fhdl(n, F, pSiconosMCP->pRes->pGetRowMap());

     pSiconosMCP->DeltaX.ScalarAddMul(zhdl, pSiconosMCP->XPrev, -1.); // Convert to incremental solution

     pSiconosMCP->XPrev = zhdl;

     pSiconosMCP->pNLP->Update(&pSiconosMCP->DeltaX);

     Fhdl.Reset();

     if (pSiconosMCP->pAbsRes) {
          pSiconosMCP->pAbsRes->Reset();
     }

     pSiconosMCP->pNLP->Residual(&Fhdl, pSiconosMCP->pAbsRes);

     if (pSiconosMCP->outputRes()) {
          pSiconosMCP->pSolver->PrintResidual(*pSiconosMCP->pRes, pSiconosMCP->iIterCurr);
     }

     Fhdl *= -1;
}

void SiconosMCPSolver::compute_nabla_Fmcp(void *env, int n, doublereal *z, NumericsMatrix *F)
{
     SiconosMCPSolver* const pSiconosMCP = static_cast<SiconosMCPSolver*>(env);

     ASSERT(F == pSiconosMCP->pMCP->nabla_Fmcp);

     pSiconosMCP->pNLP->Jacobian(pSiconosMCP->pJac);

     pSiconosMCP->TotJac++;

     if (pSiconosMCP->outputJac()) {
          silent_cout("Jacobian:" << '\n');

          if (silent_out) {
                pSiconosMCP->pJac->Print(std::cout, MatrixHandler::MAT_PRINT_TRIPLET);
          }
     }
}

void SiconosMCPSolver::collectStatsIteration(void *env, int size, double *reaction, double *velocity, double error, void *extra_data)
{
     SiconosMCPSolver* const pSiconosMCP = static_cast<SiconosMCPSolver*>(env);

     ++pSiconosMCP->iIterCurr;

     if (pSiconosMCP->outputIters()) {
          silent_cout("\tIteration(" << pSiconosMCP->iIterCurr << ") " << error);

          if (pSiconosMCP->TotJac > pSiconosMCP->iJacPrev) {
               silent_cout(" J");
               pSiconosMCP->iJacPrev = pSiconosMCP->TotJac;
          }

          silent_cout("\n");
     }

     // allow to bail out in case of multiple CTRL^C
     if (mbdyn_stop_at_end_of_iteration()) {
          throw ErrInterrupted(MBDYN_EXCEPT_ARGS);
     }
}

#endif
