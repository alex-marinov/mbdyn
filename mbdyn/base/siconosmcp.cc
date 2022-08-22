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

// This code provides interfaces to INRIA's Siconos library
// https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos/index.html
// https://github.com/siconos/siconos

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
                                   const LineSearchParameters& oLineSearch,
                                   integer iSolverId)
     :NonlinearSolver(options),
      LineSearchParameters(oLineSearch),
      pIndexMap(nullptr),
      pRes(nullptr),
      pSol(nullptr),
      pJac(nullptr),
      pNLP(nullptr),
      pSolver(nullptr),
      pMCP(nullptr),
      pOptions(nullptr),
      iJacPrev(-1),
      iIterCurr(-1),
      iSolverId(iSolverId)
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

     const DataManager* const pDM = pS->pGetDataManager();

     const VectorHandler& X = *pDM->GetpXCurr();
     //const VectorHandler& XP = *pDM->GetpXPCurr();
     const integer iNumDof = pSol->iGetSize();

     ASSERT(X.iGetSize() == iNumDof);
     ASSERT(XP.iGetSize() == iNumDof);

     pSol->Reset();
     zPrev.Reset();

     for (integer i = 1; i <= iNumDof; ++i) {
          if (pDM->GetEqualityType(i) == DofOrder::INEQUALITY) {
               zPrev(i) = (*pSol)(i) = X(i);
          }
     }

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
     pOptions = solver_options_create(iSolverId);
     pOptions->callback = static_cast<Callback*>(malloc(sizeof(Callback)));
     pOptions->callback->collectStatsIteration = &collectStatsIteration;
     pOptions->callback->env = this;
     pOptions->dparam[SICONOS_DPARAM_LSA_ALPHA_MIN] = dLambdaMin;

     const bool bResTest = pGetResTest()->GetType() != NonlinearSolverTest::NONE;
     const bool bSolTest = pGetSolTest()->GetType() != NonlinearSolverTest::NONE;

     SICONOS_STOPPING_CRITERION eCriterion = SICONOS_STOPPING_CRITERION_RESIDU;
     
     if (bResTest && bSolTest) {
          eCriterion = SICONOS_STOPPING_CRITERION_RESIDU_AND_STATIONARITY;
     } else if (!bResTest && bSolTest) {
          eCriterion = SICONOS_STOPPING_CRITERION_STATIONARITY;
     }

     pOptions->iparam[SICONOS_IPARAM_STOPPING_CRITERION] = eCriterion;
     
     auto* params = static_cast<newton_LSA_param*>(calloc(1, sizeof(newton_LSA_param)));

     pOptions->solverParameters = params;
     params->rho = dMcpRho;
     params->p = dMcpP;
     params->sigma = dMcpSigma;
     params->check_dir_quality = true;
     
     auto* const pSM = dynamic_cast<SiconosDenseSolutionManager*>(pS->pGetSolutionManager());

     if (!pSM) {
          silent_cerr("Use linear solver \"siconos dense\" for nonlinear solver \"siconos mcp\"!\n");
          throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
     }

     pRes = pSM->pResHdl();
     pSol = pSM->pSolHdl();
     pJac = pSM->pMatHdl();

     ASSERT(pJac->pGetMatrix()->storageType == NM_DENSE);

     const DataManager* const pDM = pSolver->pGetDataManager();

     const integer iNumDofsTot = pDM->iGetNumDofs();

     ASSERT(iNumDofsTot == pJac->iGetNumCols());
     ASSERT(iNumDofsTot == pRes->iGetSize());
     ASSERT(iNumDofsTot == pSol->iGetSize());

     std::array<integer, DofOrder::LASTEQUALITY> rgNumDofs{0};

     for (integer i = 1; i <= iNumDofsTot; ++i) {
          ++rgNumDofs[pDM->GetEqualityType(i)];
     }

     pMCP->n1 = rgNumDofs[DofOrder::EQUALITY];
     pMCP->n2 = rgNumDofs[DofOrder::INEQUALITY];

     if (!(pMCP->n1 > 0 && pMCP->n2 > 0)) {
          silent_cerr("Warning: nonlinear solver \"siconos mcp\" requires at least one equation and one inequality!\n");
     }

     std::vector<integer> rgRowMap(iNumDofsTot);

     std::array<integer, DofOrder::LASTEQUALITY> rgCurrIndex{0, pMCP->n1};

     for (integer i = 1; i <= iNumDofsTot; ++i) {
          rgRowMap[i - 1] = ++rgCurrIndex[pDM->GetEqualityType(i)];
     }

     pIndexMap = pSM->pGetIndexMap();

     pIndexMap->SetIndex(std::move(rgRowMap));

     pMCP->nabla_Fmcp = pJac->pGetMatrix();

     ASSERT(pSol->iGetSize() == pDM->iGetNumDofs());
     ASSERT(pMCP->n1 + pMCP->n2 == pDM->iGetNumDofs());
     ASSERT(pSol->iGetSize() == pJac->iGetNumCols());
     ASSERT(pRes->iGetSize() == pJac->iGetNumRows());

     W.SetIndexMap(pIndexMap);
     W.ResizeReset(pSol->iGetSize());

     zPrev.SetIndexMap(pIndexMap);
     zDelta.SetIndexMap(pIndexMap);
     zPrev.ResizeReset(pSol->iGetSize());
     zDelta.ResizeReset(pSol->iGetSize());
}

void SiconosMCPSolver::compute_Fmcp(void *env, int n, doublereal *z, doublereal *F)
{
     DEBUGCOUTFNAME("compute_Fmcp");

#ifdef DEBUG
     for (int i = 0; i < n; ++i) {
          DEBUGCOUT("z(" << i + 1 << ")=" << z[i] << "\n");
     }
#endif

     SiconosMCPSolver* const pSiconosMCP = static_cast<SiconosMCPSolver*>(env);

     ASSERT(n == pSiconosMCP->pSol->iGetSize());

     const SiconosVectorHandler zhdl(n, z, pSiconosMCP->pIndexMap);
     SiconosVectorHandler Fhdl(n, F, pSiconosMCP->pIndexMap);

     pSiconosMCP->zDelta.ScalarAddMul(zhdl, pSiconosMCP->zPrev, -1.);
     pSiconosMCP->zPrev = zhdl;

     pSiconosMCP->pNLP->Update(&pSiconosMCP->zDelta);

     Fhdl.Reset();

     pSiconosMCP->pNLP->Residual(&Fhdl);

     if (pSiconosMCP->outputSol()) {
          //DataManager* pDM = pSiconosMCP->pSolver->pGetDataManager();
          //const VectorHandler& X = *pDM->GetpXCurr();
          //const VectorHandler& XP = *pDM->GetpXPCurr();

          pSiconosMCP->pSolver->PrintSolution(zhdl, pSiconosMCP->iIterCurr);
     }

     if (pSiconosMCP->outputRes()) {
          pSiconosMCP->pSolver->PrintResidual(Fhdl, pSiconosMCP->iIterCurr);
     }
}

void SiconosMCPSolver::compute_nabla_Fmcp(void *env, int n, doublereal *z, NumericsMatrix *F)
{
     DEBUGCOUTFNAME("compute_nabla_Fmcp");

#ifdef DEBUG
     for (int i = 0; i < n; ++i) {
          DEBUGCOUT("z(" << i + 1 << ")=" << z[i] << "\n");
     }
#endif

     SiconosMCPSolver* const pSiconosMCP = static_cast<SiconosMCPSolver*>(env);

     ASSERT(F == pSiconosMCP->pMCP->nabla_Fmcp);

     pSiconosMCP->pNLP->Jacobian(pSiconosMCP->pJac);

     NM_scal(-1., pSiconosMCP->pJac->pGetMatrix()); // Allow us to state the complementarity condition with correct sign

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

SiconosMCPNewton::SiconosMCPNewton(const NonlinearSolverTestOptions& options, const LineSearchParameters& oLineSearch)
     :SiconosMCPSolver(options, oLineSearch, SICONOS_MCP_NEWTON_FB_FBLSA)
{
}

SiconosMCPNewton::~SiconosMCPNewton()
{
}

SiconosMCPNewtonMin::SiconosMCPNewtonMin(const NonlinearSolverTestOptions& options, const LineSearchParameters& oLineSearch)
     :SiconosMCPSolver(options, oLineSearch, SICONOS_MCP_NEWTON_MIN_FBLSA)
{
}

SiconosMCPNewtonMin::~SiconosMCPNewtonMin()
{
}

#endif
