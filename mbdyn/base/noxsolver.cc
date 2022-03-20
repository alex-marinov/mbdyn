/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2022
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 * Paolo Mantegazza     <mantegazza@aero.polimi.it>
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

#include "mbconfig.h"

#ifdef USE_TRILINOS
#include "solver.h"
#include "noxsolver.h"
#include "output.h"

#undef HAVE_BLAS // FIXME: conflicting declaration
#undef HAVE_BOOL

#include <Teuchos_ParameterList.hpp>
#include <Epetra_Operator.h>
#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>
#include <NOX.H>
#include <NOX_Epetra.H>
#include <NOX_Solver_Generic.H>
#include <NOX_Epetra_Vector.H>
#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Preconditioner.H>
#include <NOX_Epetra_LinearSystem_AztecOO.H>

class NoxStatusTest: public NOX::StatusTest::Generic {
public:
     friend class NoxNonlinearSolver;
     NoxStatusTest();
     virtual ~NoxStatusTest();
     virtual NOX::StatusTest::StatusType
     checkStatus(const NOX::Solver::Generic& problem,
                 NOX::StatusTest::CheckType checkType) override;
     virtual NOX::StatusTest::StatusType getStatus() const override;
protected:
     NOX::StatusTest::StatusType eStatus;
};

class NoxResidualTest: public NoxStatusTest {
public:
     friend class NoxNonlinearSolver;

     NoxResidualTest();
     virtual ~NoxResidualTest();
     virtual std::ostream& print(std::ostream& stream, int indent) const override;

private:
     doublereal dErrRes, dErrResDiff;
     doublereal dTolRes;
};

class NoxSolutionTest: public NoxStatusTest {
public:
     friend class NoxNonlinearSolver;

     NoxSolutionTest();
     virtual ~NoxSolutionTest();
     virtual std::ostream& print(std::ostream& stream, int indent) const override;
private:
     doublereal dErrSol;
     doublereal dTolSol;
};

class NoxNonlinearSolver : public NonlinearSolver,
                           private NoxSolverParameters,
                           private NOX::Epetra::Interface::Required,
                           private NOX::Epetra::Interface::Jacobian,
                           private NOX::Epetra::Interface::Preconditioner,
                           private Epetra_Operator
{
public:
     NoxNonlinearSolver(const NonlinearSolverTestOptions& oSolverOpt,
                        const NoxSolverParameters& oParam);

     ~NoxNonlinearSolver();

     virtual void
     Solve(const NonlinearProblem *pNLP,
           Solver *pS,
           const integer iMaxIter,
           const doublereal& Tol,
           integer& iIterCnt,
           doublereal& dErr,
           const doublereal& SolTol,
           doublereal& dSolErr) override;

private:
     struct CPUTimeGuard
     {
          explicit CPUTimeGuard(const NoxNonlinearSolver& oSolver, CPUTimeType eType)
               :oWatch(oSolver, eType) {
               oWatch.Tic();
          }

          ~CPUTimeGuard() {
               oWatch.Toc();
          }

          CPUStopWatch oWatch;
     };

     void Attach(Solver* pSolver, const NonlinearProblem* pNLP);
     bool Residual(const VectorHandler* pSol, VectorHandler* pRes);
     void Jacobian(const VectorHandler* pSol, MatrixHandler* pJac);
     virtual bool computeF(const Epetra_Vector &x, Epetra_Vector &f, const FillType fillFlag) override;
     virtual bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac) override;
     virtual bool computePreconditioner(const Epetra_Vector &x, Epetra_Operator &M, Teuchos::ParameterList *precParams) override;
     virtual int SetUseTranspose(bool UseTranspose) override;
     virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override;
     virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override;
     virtual double NormInf() const override;
     virtual const char* Label() const override;
     virtual bool UseTranspose() const override;
     virtual bool HasNormInf() const override;
     virtual const Epetra_Comm& Comm() const override;
     virtual const Epetra_Map& OperatorDomainMap() const override;
     virtual const Epetra_Map& OperatorRangeMap() const override;

     void BuildSolver(integer iMaxIter);
     void OutputIteration(integer iIterCnt, integer iTotJacPrev) const;
     Teuchos::RCP<NOX::Solver::Generic> pNonlinearSolver;
     Teuchos::RCP<NOX::Epetra::Vector> pSolutionView;
     Teuchos::ParameterList oSolverParam;
     const NonlinearProblem* pNonlinearProblem;
     Solver* pSolver;
     SolutionManager* pSolutionManager;
     mutable MyVectorHandler DeltaX, XPrev, TmpRes;
     NoxResidualTest oResTest;
     NoxSolutionTest oSolTest;
     bool bUseTranspose;
     bool bUpdateJacobian;
     Epetra_SerialComm oComm;
     Epetra_Map oMap;
};

NoxSolverParameters::NoxSolverParameters()
     :CommonNonlinearSolverParam(ALGORITHM_LINESEARCH_BASED | JACOBIAN_NEWTON, 0, false),
      dNewtonKrylovPerturbation(1e-3),
      dWrmsRelTol(1e-2),
      dWrmsAbsTol(1e-8),
      dTolLinSol(1e-10),
      iMaxIterLinSol(1000)
{
}

NoxStatusTest::NoxStatusTest()
     :eStatus(NOX::StatusTest::Unevaluated)
{
}

NoxStatusTest::~NoxStatusTest()
{
}

NOX::StatusTest::StatusType
NoxStatusTest::checkStatus(const NOX::Solver::Generic& problem,
                           NOX::StatusTest::CheckType checkType)
{
     return eStatus;
}

NOX::StatusTest::StatusType NoxStatusTest::getStatus() const
{
     return eStatus;
}

NoxResidualTest::NoxResidualTest()
     :dErrRes(-1.),
      dErrResDiff(-1.),
      dTolRes(-2.)
{
}

NoxResidualTest::~NoxResidualTest()
{
}

std::ostream& NoxResidualTest::print(std::ostream& stream, int indent) const
{
     for (int j = 0; j < indent; j ++) {
          stream << ' ';
     }

     stream << eStatus;

     stream << "F-Norm = " << NOX::Utils::sciformat(dErrRes, 3)
            << " < " << NOX::Utils::sciformat(dTolRes, 3) << "\n";

     return stream;
}

NoxSolutionTest::NoxSolutionTest()
     :dErrSol(-1.),
      dTolSol(-2.)
{
}

NoxSolutionTest::~NoxSolutionTest()
{
}

std::ostream& NoxSolutionTest::print(std::ostream& stream, int indent) const
{
     for (int j = 0; j < indent; j ++) {
          stream << ' ';
     }

     stream << eStatus;

     stream << "X-Norm = " << NOX::Utils::sciformat(dErrSol, 3)
            << " < " << NOX::Utils::sciformat(dTolSol, 3) << "\n";

     return stream;
}

NoxNonlinearSolver::NoxNonlinearSolver(const NonlinearSolverTestOptions& oSolverOpt,
                                       const NoxSolverParameters& oParam)
     :NonlinearSolver(oSolverOpt),
      NoxSolverParameters(oParam),
      pNonlinearProblem(nullptr),
      pSolver(nullptr),
      pSolutionManager(nullptr),
      bUseTranspose(false),
      bUpdateJacobian(true),
      oMap(Size, 1, oComm)
{
}

NoxNonlinearSolver::~NoxNonlinearSolver()
{
}

void
NoxNonlinearSolver::Solve(const NonlinearProblem *pNLP,
                          Solver *pS,
                          const integer iMaxIter,
                          const doublereal& dTolRes,
                          integer& iIterCnt,
                          doublereal& dResErr,
                          const doublereal& dTolSol,
                          doublereal& dSolErr)
{
     oResTest.dTolRes = dTolRes;
     oSolTest.dTolSol = dTolSol;

     Attach(pS, pNLP);

     if (!pNonlinearSolver.get()) {
          BuildSolver(iMaxIter);
     }

     pNonlinearSolver->reset(*pSolutionView);

     iIterCnt = 0;

     for (;;) {
          oResTest.dErrRes = oResTest.dErrResDiff = oSolTest.dErrSol = -1.;
          oResTest.eStatus = oSolTest.eStatus = NOX::StatusTest::Unevaluated;

          const integer iTotJacPrev = TotJac;

          NOX::StatusTest::StatusType solvStatus = pNonlinearSolver->step();

          ++iIterCnt;

          OutputIteration(iIterCnt, iTotJacPrev);

          if (solvStatus == NOX::StatusTest::Converged) {
               break;
          }

          if (solvStatus == NOX::StatusTest::Failed) {
               throw NoConvergence(MBDYN_EXCEPT_ARGS);
          }

          // allow to bail out in case of multiple CTRL^C
          if (mbdyn_stop_at_end_of_iteration()) {
               throw ErrInterrupted(MBDYN_EXCEPT_ARGS);
          }
     }

     dResErr = oResTest.dErrRes;
}

void NoxNonlinearSolver::BuildSolver(const integer iMaxIter)
{
     oMap = Epetra_Map(Size, 1, oComm);

     Teuchos::ParameterList& oPrintParam(oSolverParam.sublist("Printing"));
     Teuchos::ParameterList& oDirectionParam = oSolverParam.sublist("Direction");
     Teuchos::ParameterList& oNewtonParam = oDirectionParam.sublist("Newton");
     Teuchos::ParameterList& oLinSolParam = oNewtonParam.sublist("Linear Solver");
     Teuchos::ParameterList& oLineSearchParam = oSolverParam.sublist("Line Search");

     int iSolverOutput = 0x0;

     if (outputIters()) {
          if (uFlags & VERBOSE_MODE) {
               iSolverOutput |= NOX::Utils::Warning;
          }

          if (uFlags & PRINT_CONVERGENCE_INFO) {
               iSolverOutput |= NOX::Utils::OuterIteration |
                                NOX::Utils::InnerIteration |
                                NOX::Utils::Parameters |
                                NOX::Utils::Details |
                                NOX::Utils::OuterIterationStatusTest |
                                NOX::Utils::TestDetails;
          }
     }

     oPrintParam.set("Output Information", iSolverOutput);

     static constexpr char szNonlinearSolver[] = "Nonlinear Solver";

     if (uFlags & ALGORITHM_LINESEARCH_BASED) {
          oSolverParam.set(szNonlinearSolver, "Line Search Based");
          oLineSearchParam.set("Method", "Backtrack");
     } else if (uFlags & ALGORITHM_TRUST_REGION_BASED) {
          oSolverParam.set(szNonlinearSolver, "Trust Region Based");
     } else if (uFlags & ALGORITHM_INEXACT_TRUST_REGION_BASED) {
          oSolverParam.set(szNonlinearSolver, "Inexact Trust Region Based");
     } else if (uFlags & ALGORITHM_TENSOR_BASED) {
          oSolverParam.set(szNonlinearSolver, "Tensor Based");
          oLineSearchParam.set("Method", "Curvilinear");
     }
     
     oDirectionParam.set("Method", "Newton");

     oLinSolParam.set("Aztec Solver", "GMRES");
     oLinSolParam.set("Max Iterations", iMaxIterLinSol);
     oLinSolParam.set("Tolerance", dTolLinSol);

     oLinSolParam.set("Preconditioner", "User Defined");
     oLinSolParam.set("Preconditioner Reuse Policy", "Reuse");
     oLinSolParam.set("Max Age Of Prec", iIterationsBeforeAssembly);

     Teuchos::RCP<Epetra_Vector> pSolution = Teuchos::rcp(new Epetra_Vector(oMap, true));
     pSolutionView = Teuchos::rcp(new NOX::Epetra::Vector{pSolution, NOX::Epetra::Vector::CreateView});

     NOX::Epetra::Interface::Required& oResidualInt = *this;
     Teuchos::RCP<NOX::Epetra::Interface::Required> pResidualInt{Teuchos::rcpFromRef(oResidualInt)};
     NOX::Epetra::Interface::Preconditioner& oPrecondInt = *this;
     Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> pPrecondInt{Teuchos::rcpFromRef(oPrecondInt)};
     Epetra_Operator& oPrecondOper = *this;
     Teuchos::RCP<Epetra_Operator> pPrecondOper{Teuchos::rcpFromRef(oPrecondOper)};

     Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> pLinearSystem;

     if ((uFlags & JACOBIAN_NEWTON)) {
          NOX::Epetra::Interface::Jacobian& oJacobianInt = *this;
          Teuchos::RCP<NOX::Epetra::Interface::Jacobian> pJacobianInt{Teuchos::rcpFromRef(oJacobianInt)};
          pLinearSystem = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(oPrintParam,
                                                                            oLinSolParam,
                                                                            pJacobianInt,
                                                                            pPrecondOper,
                                                                            pPrecondInt,
                                                                            pPrecondOper,
                                                                            *pSolutionView));
     } else {
          Teuchos::RCP<NOX::Epetra::MatrixFree> pMatrixFree{new NOX::Epetra::MatrixFree(oPrintParam, pResidualInt, *pSolutionView)};

          pMatrixFree->setLambda(dNewtonKrylovPerturbation);

          Teuchos::RCP<NOX::Epetra::Interface::Jacobian> pJacobianInt(pMatrixFree);
          pLinearSystem = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(oPrintParam,
                                                                            oLinSolParam,
                                                                            pJacobianInt,
                                                                            pMatrixFree,
                                                                            pPrecondInt,
                                                                            pPrecondOper,
                                                                            *pSolutionView));
     }

     Teuchos::RCP<NOX::Epetra::Group> grpPtr =
          Teuchos::rcp(new NOX::Epetra::Group(oPrintParam,
                                              pResidualInt,
                                              *pSolutionView,
                                              pLinearSystem));

     Teuchos::RCP<NOX::StatusTest::Combo> converged =
          Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

     converged->addStatusTest(Teuchos::rcpFromRef(oResTest));

     if (oSolTest.dTolSol > 0.) {
          converged->addStatusTest(Teuchos::rcpFromRef(oSolTest));
     }

     if (dWrmsRelTol > 0. && dWrmsAbsTol > 0.) {
          Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
               Teuchos::rcp(new NOX::StatusTest::NormWRMS(dWrmsRelTol, dWrmsAbsTol));

          converged->addStatusTest(wrms);
     }

     Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
          Teuchos::rcp(new NOX::StatusTest::MaxIters(iMaxIter));

     Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
          Teuchos::rcp(new NOX::StatusTest::FiniteValue);

     Teuchos::RCP<NOX::StatusTest::Combo> pCombCriteria =
          Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

     Teuchos::RCP<Teuchos::ParameterList> pSolverParam =
          Teuchos::rcpFromRef(oSolverParam);

     pCombCriteria->addStatusTest(fv);
     pCombCriteria->addStatusTest(converged);
     pCombCriteria->addStatusTest(maxiters);

     pNonlinearSolver = NOX::Solver::buildSolver(grpPtr, pCombCriteria, pSolverParam);
}

void NoxNonlinearSolver::OutputIteration(const integer iIterCnt, const integer iTotJacPrev) const
{

#ifdef USE_MPI
     if (!bParallel || MBDynComm.Get_rank() == 0)
#endif
     {
          if (!(outputIters() || outputSolverConditionNumber())) {
               return;
          }

          silent_cout("\tIteration(" << iIterCnt << ") " << oResTest.dErrRes);

          if (TotJac > iTotJacPrev) {
               silent_cout(" J");

               if (outputSolverConditionNumber()) {
                    silent_cout(" cond=");

                    doublereal dCond;

                    if (pSolutionManager->bGetConditionNumber(dCond)) {
                         silent_cout(dCond);

                         if (outputSolverConditionStat()) {
                              AddCond(dCond);
                              silent_cout(" " << dGetCondMin() << " " << dGetCondMax() << " " << dGetCondAvg());
                         }
                    } else {
                         silent_cout("NA");
                    }
               }
          }

          if (outputCPUTime()) {
               typedef std::chrono::duration<float, std::ratio<1, 1> > FloatSec;
               auto flags = std::cout.flags();
               auto prec = std::cout.precision();

               std::cout.setf(std::ios::scientific);
               std::cout.precision(2);

               silent_cout(" CPU:" << FloatSec(dGetTimeCPU(CPU_RESIDUAL)).count()
                           << '+' << FloatSec(dGetTimeCPU(CPU_JACOBIAN)).count()
                           << '+' << FloatSec(dGetTimeCPU(CPU_LINEAR_SOLVER)).count());

               std::cout.flags(flags);
               std::cout.precision(prec);
          }

          silent_cout("\n\t\tSolErr " << oSolTest.dErrSol << std::endl);
     }
}

void NoxNonlinearSolver::Attach(Solver* pS, const NonlinearProblem* pNLP)
{
     pSolver = pS;
     pSolutionManager = pSolver->pGetSolutionManager();

     if (pNLP != pNonlinearProblem) {
          DEBUGCERR("Resetting nonlinear solver\n");
          pNonlinearSolver.reset();
          bUpdateJacobian = true;
          ResetCond();
     }

     if (!bKeepJacAcrossSteps) {
          bUpdateJacobian = true;
     }

     pNonlinearProblem = pNLP;

     VectorHandler* const pSol = pSolutionManager->pSolHdl();

     Size = pSol->iGetSize();
     pSol->Reset();
     DeltaX.ResizeReset(Size);
     XPrev.ResizeReset(Size);
     TmpRes.ResizeReset(Size);
}

bool NoxNonlinearSolver::Residual(const VectorHandler* const pSol, VectorHandler* const pRes)
{
     CPUTimeGuard oCPUTimeRes(*this, CPU_RESIDUAL);

     DeltaX.ScalarAddMul(*pSol, XPrev, -1.); // Convert to incremental solution

     bool bSolutionUpdate = false;

     for (integer i = 1; i <= Size; ++i) {
          if (DeltaX.dGetCoef(i)) {
               bSolutionUpdate = true;
               break;
          }
     }

     if (bSolutionUpdate) {
          bUpdateJacobian = true;
     }

     XPrev = *pSol;

     pNonlinearProblem->Update(&DeltaX);

     pRes->Reset();

     VectorHandler* const pAbsRes = pGetResTest()->GetAbsRes();

     if (pAbsRes) {
          pAbsRes->Reset();
     }

     pNonlinearProblem->Residual(pRes, pAbsRes);

     return bSolutionUpdate;
}

void NoxNonlinearSolver::Jacobian(const VectorHandler* const pSol, MatrixHandler* const pJac)
{
     ASSERT(pJac != nullptr);

     CPUTimeGuard oCPUTimeJac(*this, CPU_JACOBIAN);

     pNonlinearProblem->Jacobian(pJac);
     pJac->PacMat();

#ifdef USE_MPI
     if (!bParallel || MBDynComm.Get_rank() == 0)
#endif
     {
          if (outputJac()) {
               silent_cout("Jacobian:" << '\n');

               if (silent_out) {
                    pJac->Print(std::cout, MatrixHandler::MAT_PRINT_TRIPLET);
               }
          }
     }

     ++TotJac;
}

bool NoxNonlinearSolver::computeF(const Epetra_Vector& x, Epetra_Vector& F, FillType)
{
     ASSERT(pNonlinearProblem != nullptr);
     ASSERT(pSolutionManager != nullptr);

     ASSERT(x.GlobalLength() == Size);
     ASSERT(F.GlobalLength() == Size);

     const MyVectorHandler oSol(Size, x.Values());
     MyVectorHandler oRes(Size, F.Values());

     const bool bSolutionUpdate = Residual(&oSol, &oRes);

     oSolTest.eStatus = MakeSolTest(pSolver, DeltaX, oSolTest.dTolSol, oSolTest.dErrSol)
          ? NOX::StatusTest::Converged
          : NOX::StatusTest::Unconverged;

     oResTest.eStatus = MakeResTest(pSolver, pNonlinearProblem, oRes, oResTest.dTolRes, oResTest.dErrRes, oResTest.dErrResDiff)
          ? NOX::StatusTest::Converged
          : NOX::StatusTest::Unconverged;

     pSolver->CheckTimeStepLimit(oResTest.dErrRes, oResTest.dErrResDiff);

     if (bSolutionUpdate && outputSol()) {
          pSolver->PrintSolution(DeltaX, pNonlinearSolver->getNumIterations());
     }

     if (outputRes()) {
          pSolver->PrintResidual(oRes, pNonlinearSolver->getNumIterations());
     }

     oRes *= -1.; // According to conventions assumed by NOX solver

     return true;
}

bool NoxNonlinearSolver::computeJacobian(const Epetra_Vector& x, Epetra_Operator& J)
{
     ASSERT(pNonlinearProblem != nullptr);
     ASSERT(pSolutionManager != nullptr);

     ASSERT(x.GlobalLength() == Size);
     ASSERT(J.GlobalLength() == Size);

     const MyVectorHandler oSol(Size, x.Values());

     Residual(&oSol, &TmpRes); // By convention AssRes must be called always before AssJac

     if (bUpdateJacobian) {
          pSolutionManager->MatrReset();
          Jacobian(&oSol, pSolutionManager->pMatHdl());
          bUpdateJacobian = false;
     }

     return true;
}

bool NoxNonlinearSolver::computePreconditioner(const Epetra_Vector &x, Epetra_Operator &M, Teuchos::ParameterList *precParams)
{
     return computeJacobian(x, M);
}

int NoxNonlinearSolver::SetUseTranspose(bool UseTranspose)
{
     bUseTranspose = UseTranspose;

     return 0;
}

int NoxNonlinearSolver::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
     CPUTimeGuard oCPUTimeJac(*this, CPU_JACOBIAN);

     const MatrixHandler* const pJacMat = pSolutionManager->pMatHdl();

     ASSERT(X.GlobalLength() == Size);
     ASSERT(Y.GlobalLength() == Size);
     ASSERT(pJacMat->iGetSize() == Size);

     const MyVectorHandler XVec(Size, X.Values());
     MyVectorHandler YVec(Size, Y.Values());

     if (bUseTranspose) {
          pJacMat->MatTVecMul(YVec, XVec);
     } else {
          pJacMat->MatVecMul(YVec, XVec);
     }

     return 0;
}

int NoxNonlinearSolver::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
     CPUTimeGuard oCPULinearSolver(*this, CPU_LINEAR_SOLVER);

     VectorHandler* const pResVec = pSolutionManager->pResHdl();
     const VectorHandler* const pSolVec = pSolutionManager->pSolHdl();

     ASSERT(X.GlobalLength() == Size);
     ASSERT(Y.GlobalLength() == Size);
     ASSERT(Size == pResVec->iGetSize());
     ASSERT(Size == pSolVec->iGetSize());

     std::copy(X.Values(), X.Values() + Size, pResVec->pdGetVec()); // FIXME: use pSolutionManager->pdSetResVec(X.Values()) instead!

     if (bUseTranspose) {
          pSolutionManager->SolveT();
     } else {
          pSolutionManager->Solve();
     }

     std::copy(pSolVec->pdGetVec(), pSolVec->pdGetVec() + Size, Y.Values());

     return 0;
}

double NoxNonlinearSolver::NormInf() const
{
     return -1.;
}

const char* NoxNonlinearSolver::Label() const
{
     static constexpr char szLabel[] = "Jac^-1";

     return szLabel;
}

bool NoxNonlinearSolver::UseTranspose() const
{
     return bUseTranspose;
}

bool NoxNonlinearSolver::HasNormInf() const
{
     return false;
}

const Epetra_Comm& NoxNonlinearSolver::Comm() const
{
     return oComm;
}

const Epetra_Map& NoxNonlinearSolver::OperatorDomainMap() const
{
     return oMap;
}

const Epetra_Map& NoxNonlinearSolver::OperatorRangeMap() const
{
     return oMap;
}

NonlinearSolver*
pAllocateNoxNonlinearSolver(const NonlinearSolverTestOptions& oSolverOpt,
                            const NoxSolverParameters& oParam)
{
     NoxNonlinearSolver* pNLS = nullptr;

     SAFENEWWITHCONSTRUCTOR(pNLS,
                            NoxNonlinearSolver,
                            NoxNonlinearSolver(oSolverOpt,
                                               oParam));

     return pNLS;
}
#endif
