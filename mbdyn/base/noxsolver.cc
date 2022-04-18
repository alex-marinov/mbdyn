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
#ifdef USE_MPI
#include "mbcomm.h"
#endif

#undef HAVE_BLAS // FIXME: conflicting declaration
#undef HAVE_BOOL

#include <Epetra_config.h>
#include <Epetra_Operator.h>
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <Epetra_Map.h>
#include <NOX.H>
#include <NOX_Epetra.H>
#include <NOX_Solver_Generic.H>
#include <NOX_Solver_LineSearchBased.H>
#include <NOX_Epetra_Vector.H>
#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Preconditioner.H>
#include <NOX_Epetra_LinearSystem_AztecOO.H>
#include <NOX_Abstract_PrePostOperator.H>
#include <Teuchos_ParameterList.hpp>

#ifdef DEBUG_JACOBIAN
#include "sp_gradient_spmh.h"
#endif

NoxSolverParameters::NoxSolverParameters()
     :CommonNonlinearSolverParam(ALGORITHM_LINESEARCH_BASED |
                                 JACOBIAN_NEWTON |
                                 DIRECTION_NEWTON |
                                 FORCING_TERM_CONSTANT |
                                 LINESEARCH_BACKTRACK |
                                 LINEAR_SOLVER_GMRES |
                                 RECOVERY_STEP_TYPE_CONST,
                                 0,
                                 false),
#ifndef USE_SPARSE_AUTODIFF
      dNewtonKrylovPerturbation(1e-3),
#endif
      dWrmsRelTol(0.),
      dWrmsAbsTol(0.),
      dTolLinSol(1e-10),
      dMinStep(1e-12),
      dRecoveryStep(1.),
      iMaxIterLinSol(1000),
      iKrylovSubSpaceSize(300),
      iMaxIterLineSearch(200)
{
}

namespace {
     class NoxNonlinearSolver;

     class NoxStatusTest: public NOX::StatusTest::Generic {
     public:
          explicit NoxStatusTest(NoxNonlinearSolver& oNoxSolver);
          virtual ~NoxStatusTest();
          virtual NOX::StatusTest::StatusType getStatus() const override;
          inline void Reset();
     protected:
          NoxNonlinearSolver& oNoxSolver;
          NOX::StatusTest::StatusType eStatus;
     };

     class NoxResidualTest: public NoxStatusTest {
     public:
          explicit NoxResidualTest(NoxNonlinearSolver& oNoxSolver);
          virtual ~NoxResidualTest();
          virtual NOX::StatusTest::StatusType
          checkStatus(const NOX::Solver::Generic& problem,
                      NOX::StatusTest::CheckType checkType) override;
          virtual std::ostream& print(std::ostream& stream, int indent) const override;
          inline void Reset();
          inline void SetTolerance(doublereal dTol);
          inline doublereal dGetTest() const;
          inline doublereal dGetTestDiff() const;

     private:
          doublereal dErrRes, dErrResDiff;
          doublereal dTolRes;
     };

     class NoxSolutionTest: public NoxStatusTest {
     public:
          explicit NoxSolutionTest(NoxNonlinearSolver& oNoxSolver);
          virtual ~NoxSolutionTest();
          virtual NOX::StatusTest::StatusType
          checkStatus(const NOX::Solver::Generic& problem,
                      NOX::StatusTest::CheckType checkType) override;
          virtual std::ostream& print(std::ostream& stream, int indent) const override;
          inline void Reset();
          inline void SetTolerance(doublereal dTol);
          inline doublereal dGetTolerance() const;
          inline doublereal dGetTest() const;

     private:
          doublereal dErrSol;
          doublereal dTolSol;
     };

#ifdef USE_SPARSE_AUTODIFF
     class NoxMatrixFreeJacOper: public Epetra_Operator, public NOX::Epetra::Interface::Jacobian {
          friend NoxNonlinearSolver;

          explicit NoxMatrixFreeJacOper(NoxNonlinearSolver& oNoxSolver);
          virtual ~NoxMatrixFreeJacOper();

     public:
          virtual int SetUseTranspose(bool UseTranspose) override;
          virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override;
          virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override;
          virtual double NormInf() const override;
          virtual const char * Label() const override;
          virtual bool UseTranspose() const override;
          virtual bool HasNormInf() const override;
          virtual const Epetra_Comm& Comm() const override;
          virtual const Epetra_Map& OperatorDomainMap() const override;
          virtual const Epetra_Map& OperatorRangeMap() const override;
          virtual bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac) override;

     private:
          NoxNonlinearSolver& oNoxSolver;
#ifdef DEBUG_JACOBIAN
          SpGradientSparseMatrixHandler* pA;
          mutable MyVectorHandler AX;
#endif
     };
#endif

     class NoxNonlinearSolver : public NonlinearSolver,
                                private NoxSolverParameters,
                                private NOX::Epetra::Interface::Required,
                                private NOX::Epetra::Interface::Jacobian,
                                private NOX::Epetra::Interface::Preconditioner,
                                private Epetra_Operator,
                                private NOX::Abstract::PrePostOperator
     {
     public:
#ifdef USE_SPARSE_AUTODIFF
          friend NoxMatrixFreeJacOper;
#endif
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

          bool MakeSolTest(const VectorHandler& XPrev,
                           const VectorHandler& XCurr,
                           const doublereal& dTol,
                           doublereal& dTest);

          bool MakeResTest(const VectorHandler& oResVec,
                           const doublereal& dTol,
                           doublereal& dTest,
                           doublereal& dTestDiff);
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
          virtual void runPreIterate(const NOX::Solver::Generic& solver) override;
          virtual void runPostIterate(const NOX::Solver::Generic& solver) override;
          virtual void runPreSolve(const NOX::Solver::Generic& solver) override;
          virtual void runPostSolve(const NOX::Solver::Generic& solver) override;
          virtual void runPreSolutionUpdate(const NOX::Abstract::Vector& update, const NOX::Solver::Generic& solver) override;
          virtual void runPreLineSearch(const NOX::Solver::Generic& solver) override;
          virtual void runPostLineSearch(const NOX::Solver::Generic& solver) override;

          void BuildSolver(integer iMaxIter);
          void OutputIteration(integer iIterCnt, bool bJacobian) const;
          Teuchos::RCP<NOX::Solver::Generic> pNonlinearSolver;
          Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> pLinearSystem;
          Teuchos::RCP<NOX::Epetra::Vector> pSolutionView;
          Teuchos::ParameterList oSolverParam;
          Teuchos::ParameterList* pLinSolParam;
          const NonlinearProblem* pNonlinearProblem;
          Solver* pSolver;
          SolutionManager* pSolutionManager;
          mutable MyVectorHandler DeltaX, XPrev, TmpRes;
          NoxResidualTest oResTest;
          NoxSolutionTest oSolTest;
#ifdef USE_SPARSE_AUTODIFF
          NoxMatrixFreeJacOper oJacobianOperator;
#endif
          bool bUseTranspose;
          bool bUpdateJacobian;
#ifdef HAVE_MPI
          Epetra_MpiComm oComm;
#else
          Epetra_SerialComm oComm;
#endif
          Epetra_Map oMap;
          bool bInDerivativeSolver;
          bool bInLineSearch;
          mutable bool bRecomputePrecond;
     };

     NoxStatusTest::NoxStatusTest(NoxNonlinearSolver& oNoxSolver)
          :oNoxSolver(oNoxSolver),
           eStatus(NOX::StatusTest::Unevaluated)
     {
     }

     NoxStatusTest::~NoxStatusTest()
     {
     }

     void NoxStatusTest::Reset()
     {
          eStatus = NOX::StatusTest::Unevaluated;
     }

     NOX::StatusTest::StatusType NoxStatusTest::getStatus() const
     {
          return eStatus;
     }

     NoxResidualTest::NoxResidualTest(NoxNonlinearSolver& oNoxSolver)
          :NoxStatusTest(oNoxSolver),
           dErrRes(-1.),
           dErrResDiff(-1.),
           dTolRes(-2.)
     {
     }

     NoxResidualTest::~NoxResidualTest()
     {
     }

     NOX::StatusTest::StatusType
     NoxResidualTest::checkStatus(const NOX::Solver::Generic& problem,
                                  NOX::StatusTest::CheckType checkType)
     {
          if (checkType == NOX::StatusTest::None) {
               eStatus = NOX::StatusTest::Unevaluated;
               return eStatus;
          }

          const NOX::Abstract::Group& grp = problem.getSolutionGroup();

          if (!grp.isF()) {
               eStatus = NOX::StatusTest::Unevaluated;
               return eStatus;
          }

          const NOX::Abstract::Vector& FA = grp.getF();
          const NOX::Epetra::Vector& FE = dynamic_cast<const NOX::Epetra::Vector&>(FA);
          const Epetra_Vector& F = FE.getEpetraVector();
          const MyVectorHandler oResVec(F.GlobalLength(), F.Values());

          eStatus = oNoxSolver.MakeResTest(oResVec, dTolRes, dErrRes, dErrResDiff)
               ? NOX::StatusTest::Converged
               : NOX::StatusTest::Unconverged;

          return eStatus;
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

     void NoxResidualTest::Reset()
     {
          NoxStatusTest::Reset();

          dErrRes = dErrResDiff = -1.;
     }

     void NoxResidualTest::SetTolerance(doublereal dTol)
     {
          ASSERT(dTol >= 0.);

          dTolRes = dTol;
     }

     doublereal NoxResidualTest::dGetTest() const
     {
          ASSERT(eStatus != NOX::StatusTest::Unevaluated);
          ASSERT(dErrRes >= 0.);

          return dErrRes;
     }

     doublereal NoxResidualTest::dGetTestDiff() const
     {
          ASSERT(eStatus != NOX::StatusTest::Unevaluated);
          ASSERT(dErrResDiff >= 0.);

          return dErrResDiff;
     }

     NoxSolutionTest::NoxSolutionTest(NoxNonlinearSolver& oNoxSolver)
          :NoxStatusTest(oNoxSolver),
           dErrSol(-1.),
           dTolSol(-2.)
     {
     }

     NoxSolutionTest::~NoxSolutionTest()
     {
     }

     NOX::StatusTest::StatusType
     NoxSolutionTest::checkStatus(const NOX::Solver::Generic& problem,
                                  NOX::StatusTest::CheckType checkType)
     {
          if (checkType == NOX::StatusTest::None)
          {
               eStatus = NOX::StatusTest::Unevaluated;
               return eStatus;
          }

          // On the first iteration, the old and current solution are the same so
          // we should return the test as unconverged until there is a valid
          // old solution (i.e. the number of iterations is greater than zero).
          int niters = problem.getNumIterations();

          if (niters == 0)
          {
               eStatus = NOX::StatusTest::Unconverged;
               return eStatus;
          }

          // Check that F exists!
          if (!problem.getSolutionGroup().isF())
          {
               eStatus = NOX::StatusTest::Unconverged;
               return eStatus;
          }

          const NOX::Abstract::Vector& XPrevA = problem.getPreviousSolutionGroup().getX();
          const NOX::Abstract::Vector& XCurrA = problem.getSolutionGroup().getX();

          const NOX::Epetra::Vector& XPrevNE = dynamic_cast<const NOX::Epetra::Vector&>(XPrevA);
          const NOX::Epetra::Vector& XCurrNE = dynamic_cast<const NOX::Epetra::Vector&>(XCurrA);

          const Epetra_Vector& XPrevE = XPrevNE.getEpetraVector();
          const Epetra_Vector& XCurrE = XCurrNE.getEpetraVector();

          const MyVectorHandler XPrev(XPrevE.GlobalLength(), XPrevE.Values());
          const MyVectorHandler XCurr(XCurrE.GlobalLength(), XCurrE.Values());

          eStatus = oNoxSolver.MakeSolTest(XPrev, XCurr, dTolSol, dErrSol)
               ? NOX::StatusTest::Converged
               : NOX::StatusTest::Unconverged;

          return eStatus;
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

     void NoxSolutionTest::Reset()
     {
          NoxStatusTest::Reset();
          dErrSol = -1.;
     }

     void NoxSolutionTest::SetTolerance(doublereal dTol)
     {
          ASSERT(dTol >= 0.);
          dTolSol = dTol;
     }

     doublereal NoxSolutionTest::dGetTolerance() const
     {
          ASSERT(dTolSol >= 0.);
          return dTolSol;
     }

     doublereal NoxSolutionTest::dGetTest() const
     {
          ASSERT(eStatus != NOX::StatusTest::Unevaluated);
          ASSERT(dErrSol >= 0.);

          return dErrSol;
     }

#ifdef USE_SPARSE_AUTODIFF
     NoxMatrixFreeJacOper::NoxMatrixFreeJacOper(NoxNonlinearSolver& oNoxSolver)
          :oNoxSolver(oNoxSolver)
#ifdef DEBUG_JACOBIAN
          ,pA(nullptr)
#endif
     {
     }

     NoxMatrixFreeJacOper::~NoxMatrixFreeJacOper()
     {
#ifdef DEBUG_JACOBIAN
          if (pA) {
               SAFEDELETE(pA);
          }
#endif
     }

     int NoxMatrixFreeJacOper::SetUseTranspose(bool UseTranspose)
     {
          if (UseTranspose) {
               throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
          }

          return 0;
     }

     int NoxMatrixFreeJacOper::Apply(const Epetra_MultiVector& x, Epetra_MultiVector& y) const
     {
          NoxNonlinearSolver::CPUTimeGuard oCPUTimeRes(oNoxSolver, NoxNonlinearSolver::CPU_JACOBIAN);

          ASSERT(oNoxSolver.Size > 0);
          ASSERT(oNoxSolver.Size == y.GlobalLength());
          ASSERT(oNoxSolver.Size == x.GlobalLength());

          MyVectorHandler Jac(y.GlobalLength(), y.Values());
          const MyVectorHandler X(x.GlobalLength(), x.Values());

          oNoxSolver.pNonlinearProblem->Jacobian(&Jac, &X);

#ifdef DEBUG_JACOBIAN
          for (integer i = 1; i <= X.iGetSize(); ++ i) {
               DEBUGCERR("X(" << i << ")=" << X(i) << "\n");
          }

          for (integer i = 1; i <= Jac.iGetSize(); ++i) {
               DEBUGCERR("Jac(" << i << ") = " << Jac(i) << "\n");
          }

          for (integer i = 1; i <= AX.iGetSize(); ++i) {
               DEBUGCERR("AX(" << i << ") = " << AX(i) << "\n");
          }

          pA->MatVecMul(AX, X);

          const doublereal dTolRel = sqrt(std::numeric_limits<doublereal>::epsilon());
          const doublereal dTolAbs = sqrt(std::numeric_limits<doublereal>::epsilon());
          const doublereal dNormAX = AX.Norm();

          for (integer i = 1; i <= X.iGetSize(); ++i) {
               if (std::fabs(AX(i) - Jac(i)) > dTolRel + dTolAbs * dNormAX) {
                    DEBUGCERR("Error: AX(" << i << ")=" << AX(i) << " Jac(" << i << ")=" << Jac(i) << "\n");
               }
          }

          Jac = AX;
#endif
          return 0;
     }

     int NoxMatrixFreeJacOper::ApplyInverse(const Epetra_MultiVector& x, Epetra_MultiVector& y) const
     {
          throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
     }

     double NoxMatrixFreeJacOper::NormInf() const
     {
          throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
     }

     const char* NoxMatrixFreeJacOper::Label() const
     {
          static constexpr char szLabel[] = "Jac";
          return szLabel;
     }

     bool NoxMatrixFreeJacOper::UseTranspose() const
     {
          return false;
     }

     bool NoxMatrixFreeJacOper::HasNormInf() const
     {
          return false;
     }

     const Epetra_Comm & NoxMatrixFreeJacOper::Comm() const
     {
          return oNoxSolver.oComm;
     }

     const Epetra_Map & NoxMatrixFreeJacOper::OperatorDomainMap() const
     {
          return oNoxSolver.oMap;
     }

     const Epetra_Map & NoxMatrixFreeJacOper::OperatorRangeMap() const
     {
          return oNoxSolver.oMap;
     }

     bool NoxMatrixFreeJacOper::computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac)
     {
          DEBUGCERR("computeJacobian()\n");

          const MyVectorHandler oSol(oNoxSolver.Size, x.Values());

          oNoxSolver.Residual(&oSol, &oNoxSolver.TmpRes); // By convention AssRes must be called always before AssJac

#ifdef DEBUG_JACOBIAN
          AX.Resize(oNoxSolver.Size);

          if (!pA) {
               SAFENEWWITHCONSTRUCTOR(pA,
                                      SpGradientSparseMatrixHandler,
                                      SpGradientSparseMatrixHandler(x.GlobalLength(), x.GlobalLength()));
          }

          MyVectorHandler Y(oNoxSolver.Size), JY(oNoxSolver.Size), JYRef(oNoxSolver.Size);

          for (integer i = 1; i <= oNoxSolver.Size; ++i) {
               Y.Reset();
               JY.Reset();
               Y(i) = 1.;

               DEBUGCERR("Assemble Jacobian vector product for degree of freedom number " << i << "\n");

               oNoxSolver.pNonlinearProblem->Jacobian(&JY, &Y);

               DEBUGCERR("Assemble sparse Jacobian matrix ...\n");

               oNoxSolver.pNonlinearProblem->Jacobian(pA);

               pA->PacMat();

               pA->MatVecMul(JYRef, Y);

               const doublereal dTol = (1. + JYRef.Norm()) * pow(std::numeric_limits<doublereal>::epsilon(), 0.5);

               DEBUGCERR("Check Jacobian vector product " << i << "\n");

               for (integer j = 1; j <= oNoxSolver.Size; ++j) {
                    if (fabs(JYRef(j) - JY(j)) > dTol) {
                         DEBUGCERR("Jacobian check: " << i << " JYRef(" << j << ")=" << JYRef(j) << ", JY(" << j << ")=" << JY(j) << "\n");
                    }
               }

               DEBUGCERR("End of check Jacobian vector product " << i << "\n");
          }
#endif
          return true;
     }
#endif

     NoxNonlinearSolver::NoxNonlinearSolver(const NonlinearSolverTestOptions& oSolverOpt,
                                            const NoxSolverParameters& oParam)
          :NonlinearSolver(oSolverOpt),
           NoxSolverParameters(oParam),
           pLinSolParam(nullptr),
           pNonlinearProblem(nullptr),
           pSolver(nullptr),
           pSolutionManager(nullptr),
           oResTest(*this),
           oSolTest(*this),
           oJacobianOperator(*this),
           bUseTranspose(false),
           bUpdateJacobian(true),
#ifdef HAVE_MPI
           oComm(MBDynComm),
#endif
           oMap(Size, 1, oComm),
           bInDerivativeSolver(true),
           bInLineSearch(false),
           bRecomputePrecond(false)
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
          DEBUGCERR("Solve()\n");

          bInLineSearch = false;
          SetNonlinearSolverHint(LINESEARCH_ITERATION_CURR, 0);
          SetNonlinearSolverHint(LINESEARCH_LAMBDA_CURR, 1.);
          
          oResTest.SetTolerance(dTolRes);
          oSolTest.SetTolerance(dTolSol);

          Attach(pS, pNLP);

          if (!pNonlinearSolver.get()) {
               BuildSolver(iMaxIter);
          }

          pNonlinearSolver->reset(*pSolutionView);
        
          iIterCnt = 0;
          bRecomputePrecond = !bKeepJacAcrossSteps;
          
          for (;;) {
               oResTest.Reset();
               oSolTest.Reset();

               const integer TotJacPrev = TotJac;

               DEBUGCERR("Nonlinear solver step(" << iIterCnt << ")\n");

               NOX::StatusTest::StatusType solvStatus = pNonlinearSolver->step();

               if (oResTest.getStatus() == NOX::StatusTest::Unevaluated) {
                    oResTest.checkStatus(*pNonlinearSolver, NOX::StatusTest::CheckType::Complete);
               }

               ++iIterCnt;

               OutputIteration(iIterCnt, TotJac > TotJacPrev);
               pSolver->CheckTimeStepLimit(oResTest.dGetTest(), oResTest.dGetTestDiff());

               dResErr = oResTest.dGetTest();
               dSolErr = oSolTest.getStatus() != NOX::StatusTest::Unevaluated
                    ? oSolTest.dGetTest()
                    : 0.;

               if (solvStatus == NOX::StatusTest::Converged) {
                    break;
               }

               if (iIterCnt > iMaxIter || solvStatus == NOX::StatusTest::Failed) {
                    throw NoConvergence(MBDYN_EXCEPT_ARGS);
               }

               // allow to bail out in case of multiple CTRL^C
               if (mbdyn_stop_at_end_of_iteration()) {
                    throw ErrInterrupted(MBDYN_EXCEPT_ARGS);
               }
          }
     }

     bool NoxNonlinearSolver::MakeSolTest(const VectorHandler& XPrev,
                                          const VectorHandler& XCurr,
                                          const doublereal& dTol,
                                          doublereal& dTest)
     {
          DeltaX.ScalarAddMul(XCurr, XPrev, -1.);

          return NonlinearSolver::MakeSolTest(pSolver, DeltaX, dTol, dTest);
     }

     bool NoxNonlinearSolver::MakeResTest(const VectorHandler& oResVec,
                                          const doublereal& dTol,
                                          doublereal& dTest,
                                          doublereal& dTestDiff)
     {
          bool bStatus = NonlinearSolver::MakeResTest(pSolver, pNonlinearProblem, oResVec, dTol, dTest, dTestDiff);

          return bStatus;
     }

     void NoxNonlinearSolver::BuildSolver(const integer iMaxIter)
     {
          oMap = Epetra_Map(Size, 1, oComm);

          Teuchos::ParameterList& oPrintParam(oSolverParam.sublist("Printing"));
          Teuchos::ParameterList& oDirectionParam = oSolverParam.sublist("Direction");
          Teuchos::ParameterList& oNewtonParam = oDirectionParam.sublist("Newton");
          Teuchos::ParameterList& oLinSolParam = oNewtonParam.sublist("Linear Solver");
          
          pLinSolParam = &oLinSolParam;
          
          NOX::Abstract::PrePostOperator& oPrePost = *this;
          Teuchos::RCP<NOX::Abstract::PrePostOperator> pPrePost = Teuchos::rcpFromRef(oPrePost);
          oSolverParam.sublist("Solver Options").set("User Defined Pre/Post Operator", pPrePost);

          int iSolverOutput = 0x0;

          if (outputIters()) {
               if (uFlags & VERBOSE_MODE) {
                    iSolverOutput |= NOX::Utils::Warning;
               }

               if (uFlags & PRINT_CONVERGENCE_INFO) {
                    iSolverOutput |= NOX::Utils::OuterIteration |
                         NOX::Utils::InnerIteration |
                         NOX::Utils::LinearSolverDetails |
                         NOX::Utils::StepperIteration |
                         NOX::Utils::StepperDetails |
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
               Teuchos::ParameterList& oLineSearchParam = oSolverParam.sublist("Line Search");

               std::string strLineSearchMethod;
               
               if (uFlags & LINESEARCH_BACKTRACK) {
                    strLineSearchMethod = "Backtrack";
               } else if (uFlags & LINESEARCH_POLYNOMIAL) {
                    strLineSearchMethod = "Polynomial";
               } else if (uFlags & LINESEARCH_MORE_THUENTE) {
                    strLineSearchMethod = "More'-Thuente";
               } else {
                    strLineSearchMethod = "Full Step";
               }

               oLineSearchParam.set("Method", strLineSearchMethod);
               Teuchos::ParameterList& oLineSearchMethod = oLineSearchParam.sublist(strLineSearchMethod);
               oLineSearchMethod.set("Max Iters", iMaxIterLineSearch);
               oLineSearchMethod.set("Minimum Step", dMinStep);
               oLineSearchMethod.set("Recovery Step", dRecoveryStep);

               if (uFlags & RECOVERY_STEP_TYPE_CONST) {
                    oLineSearchMethod.set("Recovery Step Type", "Constant");
               } else if (uFlags & RECOVERY_STEP_TYPE_LAST_STEP) {
                    oLineSearchMethod.set("Recovery Step Type", "Last Computed Step");
               }
          } else if (uFlags & ALGORITHM_TRUST_REGION_BASED) {
               oSolverParam.set(szNonlinearSolver, "Trust Region Based");
          } else if (uFlags & ALGORITHM_INEXACT_TRUST_REGION_BASED) {
               oSolverParam.set(szNonlinearSolver, "Inexact Trust Region Based");
          } else if (uFlags & ALGORITHM_TENSOR_BASED) {
               oSolverParam.set(szNonlinearSolver, "Tensor Based");
               Teuchos::ParameterList& oLineSearchParam = oSolverParam.sublist("Line Search");
               oLineSearchParam.set("Method", "Curvilinear");
               Teuchos::ParameterList& oCurvilinearParam = oLineSearchParam.sublist("Curvilinear");
               oCurvilinearParam.set("Minimum Step", dMinStep);
          }

          if (uFlags & DIRECTION_NEWTON) {
               oDirectionParam.set("Method", "Newton");

               if (uFlags & FORCING_TERM_CONSTANT) {
                    oNewtonParam.set("Forcing Term Method", "Constant");
               } else if (uFlags & FORCING_TERM_TYPE1) {
                    oNewtonParam.set("Forcing Term Method", "Type 1");
               } else if (uFlags & FORCING_TERM_TYPE2) {
                    oNewtonParam.set("Forcing Term Method", "Type 2");
               }
          } else if (uFlags & DIRECTION_STEEPEST_DESCENT) {
               oDirectionParam.set("Method", "Steepest Descent");
          } else if (uFlags & DIRECTION_NONLINEAR_CG) {
               oDirectionParam.set("Method", "NonlinearCG");
          } else if (uFlags & DIRECTION_BROYDEN) {
               oDirectionParam.set("Method", "Broyden");
               Teuchos::ParameterList& oBroyden = oDirectionParam.sublist("Broyden");
               oBroyden.set("Restart Frequency", iIterationsBeforeAssembly);

               if (uFlags & FORCING_TERM_CONSTANT) {
                    oBroyden.set("Forcing Term Method", "Constant");
               } else if (uFlags & FORCING_TERM_TYPE1) {
                    oBroyden.set("Forcing Term Method", "Type 1");
               } else if (uFlags & FORCING_TERM_TYPE2) {
                    oBroyden.set("Forcing Term Method", "Type 2");
               }
          }

          if (uFlags & LINEAR_SOLVER_GMRES) {
               oLinSolParam.set("Aztec Solver", "GMRES");
          } else if (uFlags & LINEAR_SOLVER_CG) {
               oLinSolParam.set("Aztec Solver", "CG");
          } else if (uFlags & LINEAR_SOLVER_CGS) {
               oLinSolParam.set("Aztec Solver", "CGS");
          } else if (uFlags & LINEAR_SOLVER_TFQMR) {
               oLinSolParam.set("Aztec Solver", "TFQMR");
          } else if (uFlags & LINEAR_SOLVER_BICGSTAB) {
               oLinSolParam.set("Aztec Solver", "BiCGStab");
          }

          oLinSolParam.set("Max Iterations", iMaxIterLinSol);
          oLinSolParam.set("Tolerance", dTolLinSol);

          oLinSolParam.set("Preconditioner", "User Defined");
          oLinSolParam.set("Preconditioner Reuse Policy", "Reuse");
          oLinSolParam.set("Max Age Of Prec", iIterationsBeforeAssembly);
          oLinSolParam.set("Size of Krylov Subspace", std::min(Size, iKrylovSubSpaceSize));
          oLinSolParam.set("Use Preconditioner as Solver", (uFlags & USE_PRECOND_AS_SOLVER) != 0);
          
          if (uFlags & PRINT_CONVERGENCE_INFO) {
               oLinSolParam.set("Output Frequency", 1);
          }
          
          Teuchos::RCP<Epetra_Vector> pSolution = Teuchos::rcp(new Epetra_Vector(oMap, true));
          pSolutionView = Teuchos::rcp(new NOX::Epetra::Vector{pSolution, NOX::Epetra::Vector::CreateView});

          NOX::Epetra::Interface::Required& oResidualInt = *this;
          Teuchos::RCP<NOX::Epetra::Interface::Required> pResidualInt{Teuchos::rcpFromRef(oResidualInt)};
          NOX::Epetra::Interface::Preconditioner& oPrecondInt = *this;
          Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> pPrecondInt{Teuchos::rcpFromRef(oPrecondInt)};
          Epetra_Operator& oPrecondOper = *this;
          Teuchos::RCP<Epetra_Operator> pPrecondOper{Teuchos::rcpFromRef(oPrecondOper)};

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
#ifdef USE_SPARSE_AUTODIFF
               Teuchos::RCP<NOX::Epetra::Interface::Jacobian> pJacobianInt{Teuchos::rcpFromRef(oJacobianOperator)};
               Teuchos::RCP<Epetra_Operator> pJacobianOper{Teuchos::rcpFromRef(oJacobianOperator)};
#else
               Teuchos::RCP<NOX::Epetra::MatrixFree> pJacobianOper{new NOX::Epetra::MatrixFree(oPrintParam, pResidualInt, *pSolutionView)};

               pJacobianOper->setLambda(dNewtonKrylovPerturbation);

               Teuchos::RCP<NOX::Epetra::Interface::Jacobian> pJacobianInt(pJacobianOper);
#endif
               pLinearSystem = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(oPrintParam,
                                                                                 oLinSolParam,
                                                                                 pJacobianInt,
                                                                                 pJacobianOper,
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

          if (oSolTest.dGetTolerance() > 0.) {
               converged->addStatusTest(Teuchos::rcpFromRef(oSolTest));
          }

          // FIXME: In some cases the derivative solver cannot converge with WRMS test enabled.
          if (dWrmsRelTol > 0. && dWrmsAbsTol > 0. && !bInDerivativeSolver) {
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

     void NoxNonlinearSolver::OutputIteration(integer iIterCnt, bool bJacobian) const
     {

#ifdef USE_MPI
          if (!bParallel || MBDynComm.Get_rank() == 0)
#endif
          {
               if (!(outputIters() || outputSolverConditionNumber())) {
                    return;
               }

               silent_cout("\tIteration(" << iIterCnt << ") " << oResTest.dGetTest());

               if (bJacobian) {
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

               if (oSolTest.getStatus() != NOX::StatusTest::Unevaluated) {
                    silent_cout("\n\t\tSolErr " << oSolTest.dGetTest());
               }

               silent_cout('\n');
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

               if (pNonlinearProblem) {
                    bInDerivativeSolver = false;
               }
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
          DEBUGCERR("computeF()\n");

          ASSERT(pNonlinearProblem != nullptr);
          ASSERT(pSolutionManager != nullptr);

          ASSERT(x.GlobalLength() == Size);
          ASSERT(F.GlobalLength() == Size);

          const MyVectorHandler oSol(Size, x.Values());
          MyVectorHandler oRes(Size, F.Values());

          if (bInLineSearch) {
               auto& oLineSearch = dynamic_cast<const NOX::Solver::LineSearchBased&>(*pNonlinearSolver);
               
               SetNonlinearSolverHint(LINESEARCH_LAMBDA_CURR, oLineSearch.getStepSize());
               
               DEBUGCERR("line search iteration " << GetNonlinearSolverHint(LINESEARCH_ITERATION_CURR)
                         << ": lambda=" << oLineSearch.getStepSize() << "\n");
          }
          
          const bool bSolutionUpdate = Residual(&oSol, &oRes);

          if (bSolutionUpdate && outputSol()) {
               pSolver->PrintSolution(DeltaX, pNonlinearSolver->getNumIterations());
          }

          if (outputRes()) {
               pSolver->PrintResidual(oRes, pNonlinearSolver->getNumIterations());
          }

          oRes *= -1.; // According to conventions assumed by NOX solver

          if (bInLineSearch) {
               integer iIterCurr = GetNonlinearSolverHint(LINESEARCH_ITERATION_CURR);
               
               SetNonlinearSolverHint(LINESEARCH_ITERATION_CURR, iIterCurr + 1);
          }
          
          return true;
     }

     bool NoxNonlinearSolver::computeJacobian(const Epetra_Vector& x, Epetra_Operator& J)
     {
          ASSERT(pNonlinearProblem != nullptr);
          ASSERT(pSolutionManager != nullptr);

          ASSERT(x.GlobalLength() == Size);

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
          DEBUGCERR("computePreconditioner()\n");

          return computeJacobian(x, M);
     }

     int NoxNonlinearSolver::SetUseTranspose(bool UseTranspose)
     {
          DEBUGCERR("SetUseTranspose(" << UseTranspose << ")\n");
          bUseTranspose = UseTranspose;

          return 0;
     }

     int NoxNonlinearSolver::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
     {
          DEBUGCERR("Apply()\n");

          CPUTimeGuard oCPUTimeJac(*this, CPU_JACOBIAN);

          const MatrixHandler* const pJacMat = pSolutionManager->pMatHdl();

          ASSERT(X.GlobalLength() == Size);
          ASSERT(Y.GlobalLength() == Size);
          ASSERT(pJacMat->iGetNumRows() == Size);
          ASSERT(pJacMat->iGetNumCols() == Size);

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
          if (bRecomputePrecond) {
               DEBUGCERR("Recompute preconditioner for the first iteration ...\n");

               const NOX::Abstract::Group& oGroupA = pNonlinearSolver->getSolutionGroup();
               const NOX::Abstract::Vector& oSolVecA = oGroupA.getX();
               const NOX::Epetra::Vector& oSolVecE = dynamic_cast<const NOX::Epetra::Vector&>(oSolVecA);
               
               pLinearSystem->recomputePreconditioner(oSolVecE, *pLinSolParam);
               
               bRecomputePrecond = false;
          }
          
          DEBUGCERR("ApplyInverse()\n");

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
          DEBUGCERR("NormInf()\n");

          return -1.;
     }

     const char* NoxNonlinearSolver::Label() const
     {
          DEBUGCERR("Label()\n");

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

     void NoxNonlinearSolver::runPreIterate(const NOX::Solver::Generic& solver)
     {
          DEBUGCERR("runPreIterate()\n");
     }

     void NoxNonlinearSolver::runPostIterate(const NOX::Solver::Generic& solver)
     {
          DEBUGCERR("runPostIterate()\n");
     }

     void NoxNonlinearSolver::runPreSolve(const NOX::Solver::Generic& solver)
     {
          DEBUGCERR("runPreSolve()\n");
     }

     void NoxNonlinearSolver::runPostSolve(const NOX::Solver::Generic& solver)
     {
          DEBUGCERR("runPostSolve()\n");
     }

     void NoxNonlinearSolver::runPreSolutionUpdate(const NOX::Abstract::Vector& update, const NOX::Solver::Generic& solver)
     {
          DEBUGCERR("runPreSolutionUpdate()\n");
     }

     void NoxNonlinearSolver::runPreLineSearch(const NOX::Solver::Generic& solver)
     {
          DEBUGCERR("runPreLineSearch()\n");
          
          bInLineSearch = true;

          SetNonlinearSolverHint(LINESEARCH_LAMBDA_CURR, 1.);
          SetNonlinearSolverHint(LINESEARCH_ITERATION_CURR, 0);
     }

     void NoxNonlinearSolver::runPostLineSearch(const NOX::Solver::Generic& solver)
     {
          DEBUGCERR("runPostLineSearch()\n");
          
          bInLineSearch = false;
          
          SetNonlinearSolverHint(LINESEARCH_LAMBDA_CURR, 1.);
          SetNonlinearSolverHint(LINESEARCH_ITERATION_CURR, 0);
     }
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
