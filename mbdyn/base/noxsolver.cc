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
#include "solman.h"
#include "solver.h"
#include "noxsolver.h"
#include "output.h"
#ifdef USE_MPI
#include "mbcomm.h"
#endif

// FIXME: Compiler fails with -Werror if HAVE_BLAS is redefined inside a Trilinos header
#define HAVE_BLAS_SAVE HAVE_BLAS
#define HAVE_BOOL_SAVE HAVE_BOOL
#undef HAVE_BLAS
#undef HAVE_BOOL

#include <Epetra_config.h>
#include <Epetra_Operator.h>
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <Epetra_Map.h>
#include <AztecOO.h>
#include <NOX.H>
#include <NOX_Epetra.H>
#include <NOX_Solver_Generic.H>
#include <NOX_Solver_LineSearchBased.H>
#include <NOX_Epetra_Vector.H>
#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Preconditioner.H>
#include <NOX_Epetra_LinearSystem.H>
#include <NOX_Abstract_PrePostOperator.H>
#include <Teuchos_ParameterList.hpp>

#undef HAVE_BLAS
#undef HAVE_BOOL
#define HAVE_BLAS HAVE_BLAS_SAVE
#define HAVE_BOOL HAVE_BOOL_SAVE
#undef HAVE_BLAS_SAVE
#undef HAVE_BOOL_SAVE

#ifdef DEBUG_JACOBIAN
#include "sp_gradient_spmh.h"
#endif

NoxSolverParameters::NoxSolverParameters()
     :CommonNonlinearSolverParam(SOLVER_LINESEARCH_BASED |
                                 JACOBIAN_NEWTON |
                                 DIRECTION_NEWTON |
                                 FORCING_TERM_CONSTANT |
                                 LINESEARCH_BACKTRACK |
                                 LINEAR_SOLVER_GMRES |
                                 RECOVERY_STEP_TYPE_CONST,
                                 0,
                                 false),
      dWrmsRelTol(0.),
      dWrmsAbsTol(0.),
      dTolLinSol(1e-10),
      dMinStep(1e-12),
      dRecoveryStep(1.),
      dForcingTermMinTol(1e-6),
      dForcingTermMaxTol(1e-2),
      dForcingTermAlpha(1.5),
      dForcingTermGamma(0.9),
      iMaxIterLinSol(1000),
      iKrylovSubSpaceSize(300),
      iMaxIterLineSearch(200),
      iInnerIterBeforeAssembly(std::numeric_limits<integer>::max())
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

     class NoxNonlinearSolver : public NonlinearSolver,
                                private NoxSolverParameters,
                                private NOX::Epetra::Interface::Required,
                                private NOX::Epetra::Interface::Jacobian,
                                private NOX::Epetra::Interface::Preconditioner,
                                private Epetra_Operator,
                                private NOX::Abstract::PrePostOperator,
                                private NOX::Epetra::LinearSystem
     {
     public:
          friend NoxMatrixFreeJacOper;

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
          void Residual(const VectorHandler* pSol, VectorHandler* pRes);
          void Jacobian();
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

          virtual bool
          applyJacobian(const NOX::Epetra::Vector& input,
                        NOX::Epetra::Vector& result) const override;

          virtual bool
          applyJacobianTranspose(const NOX::Epetra::Vector& input,
                                 NOX::Epetra::Vector& result) const override;

          virtual bool
          applyJacobianInverse(Teuchos::ParameterList &params,
                               const NOX::Epetra::Vector &input,
                               NOX::Epetra::Vector &result) override;

          virtual bool
          applyRightPreconditioning(bool useTranspose,
                                    Teuchos::ParameterList& params,
                                    const NOX::Epetra::Vector& input,
                                    NOX::Epetra::Vector& result) const override;

          virtual Teuchos::RCP<NOX::Epetra::Scaling>
          getScaling() override;

          virtual void
          resetScaling(const Teuchos::RCP<NOX::Epetra::Scaling>& s) override;

          virtual bool
          computeJacobian(const NOX::Epetra::Vector& x) override;

          virtual bool
          createPreconditioner(const NOX::Epetra::Vector& x,
                               Teuchos::ParameterList& p,
                               bool recomputeGraph) const override;

          virtual bool
          destroyPreconditioner() const override;

          virtual bool
          recomputePreconditioner(const NOX::Epetra::Vector& x,
                                  Teuchos::ParameterList& linearSolverParams) const override;

          virtual PreconditionerReusePolicyType
          getPreconditionerPolicy(bool advanceReuseCounter=true) override;

          virtual bool
          isPreconditionerConstructed() const override;

          virtual Teuchos::RCP<const Epetra_Operator>
          getJacobianOperator() const override;

          virtual Teuchos::RCP<Epetra_Operator>
          getJacobianOperator() override;

          virtual Teuchos::RCP<const Epetra_Operator>
          getGeneratedPrecOperator() const override;

          virtual Teuchos::RCP<Epetra_Operator>
          getGeneratedPrecOperator() override;

          virtual void
          setJacobianOperatorForSolve(const Teuchos::RCP<const Epetra_Operator>& solveJacOp) override;

          virtual bool
          hasPreconditioner() const;

          virtual void
          setPrecOperatorForSolve(const Teuchos::RCP<const Epetra_Operator>& solvePrecOp) override;

          inline void ResetPrecondReuse() const;
          inline void ForcePrecondRebuild() const;

          void BuildSolver(integer iMaxIter);
          void OutputIteration(integer iIterCnt, bool bJacobian) const;
          Teuchos::RCP<NOX::Solver::Generic> pNonlinearSolver;
          Teuchos::RCP<NOX::Epetra::Vector> pSolutionView;
          Teuchos::ParameterList oSolverParam;
          const NonlinearProblem* pNonlinearProblem;
          Solver* pSolver;
          SolutionManager* pSolutionManager;
          mutable MyVectorHandler DeltaX, XPrev, TmpRes;
          NoxResidualTest oResTest;
          NoxSolutionTest oSolTest;
          NoxMatrixFreeJacOper oMatFreeJacOper;
          NOX::Epetra::Interface::Jacobian* pJacInt;
          Teuchos::RCP<Epetra_Operator> pJacOper;
          NOX::Epetra::Interface::Preconditioner* pPrecInt;
          mutable bool bUseTranspose;
          bool bUpdateJacobian;
#ifdef HAVE_MPI
          Epetra_MpiComm oComm;
#else
          Epetra_SerialComm oComm;
#endif
          Epetra_Map oMap;
          bool bInDerivativeSolver;
          bool bInLineSearch;
          mutable integer iPrecInnerIterCnt;
          mutable integer iPrecInnerIterCntTot;
          mutable integer iInnerIterCntTot;
          AztecOO oIterativeLinSol;
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

          // According to "Numerical recipes in C the art of scientific computing" / William H. Press [et al.]. – 2nd ed.
          const doublereal dFirstResFact = problem.getNumIterations() == 0 ? 1e-2 : 1.;

          eStatus = oNoxSolver.MakeResTest(oResVec, dFirstResFact * dTolRes, dErrRes, dErrResDiff)
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

#ifdef DEBUG
          if(dErrResDiff < 0.) {
               DEBUGCERR("Warning: dErrResDiff was not evaluated!\n");
          }
#endif
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
                    DEBUGCERR("Jacobian check failed: AX(" << i << ")=" << AX(i) << " Jac(" << i << ")=" << Jac(i) << "\n");
                    ASSERT(0);
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

          DEBUGCERR("Assemble sparse Jacobian matrix ...\n");

          oNoxSolver.pNonlinearProblem->Jacobian(pA);

          pA->PacMat();

#if DEBUG_JACOBIAN >= 2
          MyVectorHandler Y(oNoxSolver.Size), JY(oNoxSolver.Size), JYRef(oNoxSolver.Size);

          for (integer i = 1; i <= oNoxSolver.Size; ++i) {
               Y.Reset();
               JY.Reset();
               Y(i) = 1.;

               DEBUGCERR("Assemble Jacobian vector product for degree of freedom number " << i << "\n");

               oNoxSolver.pNonlinearProblem->Jacobian(&JY, &Y);

               pA->MatVecMul(JYRef, Y);

               const doublereal dTol = (1. + JYRef.Norm()) * pow(std::numeric_limits<doublereal>::epsilon(), 0.5);

               DEBUGCERR("Check Jacobian vector product " << i << "\n");

               for (integer j = 1; j <= oNoxSolver.Size; ++j) {
                    if (fabs(JYRef(j) - JY(j)) > dTol) {
                         DEBUGCERR("Jacobian check failed: " << i << " JYRef(" << j << ")=" << JYRef(j) << ", JY(" << j << ")=" << JY(j) << "\n");
                         ASSERT(0);
                    }
               }

               DEBUGCERR("End of check Jacobian vector product " << i << "\n");
          }
#endif
#endif
          return true;
     }

     NoxNonlinearSolver::NoxNonlinearSolver(const NonlinearSolverTestOptions& oSolverOpt,
                                            const NoxSolverParameters& oParam)
          :NonlinearSolver(oSolverOpt),
           NoxSolverParameters(oParam),
           pNonlinearProblem(nullptr),
           pSolver(nullptr),
           pSolutionManager(nullptr),
           oResTest(*this),
           oSolTest(*this),
           oMatFreeJacOper(*this),
           pJacInt(nullptr),
           pJacOper(nullptr),
           pPrecInt(nullptr),
           bUseTranspose(false),
           bUpdateJacobian(true),
#ifdef HAVE_MPI
           oComm(MBDynComm),
#endif
           oMap(Size, 1, oComm),
           bInDerivativeSolver(true),
           bInLineSearch(false),
           iInnerIterCntTot(0)
     {
          ForcePrecondRebuild();
     }

     NoxNonlinearSolver::~NoxNonlinearSolver()
     {
          silent_cerr("total inner iterations: " << iInnerIterCntTot << "\n");
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

          if (!bKeepJacAcrossSteps) {
               ForcePrecondRebuild();
          }

          iIterCnt = 0;

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

          oNewtonParam.set("Forcing Term Minimum Tolerance", dForcingTermMinTol);
          oNewtonParam.set("Forcing Term Maximum Tolerance", dForcingTermMaxTol);
          oNewtonParam.set("Forcing Term Alpha", dForcingTermAlpha);
          oNewtonParam.set("Forcing Term Gamma", dForcingTermGamma);

          NOX::Abstract::PrePostOperator& oPrePost = *this;
          Teuchos::RCP<NOX::Abstract::PrePostOperator> pPrePost = Teuchos::rcpFromRef(oPrePost);
          oSolverParam.sublist("Solver Options").set("User Defined Pre/Post Operator", pPrePost);

          int iSolverOutput = 0;

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

          if (uFlags & SOLVER_LINESEARCH_BASED) {
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

               if (uFlags & SUFFICIENT_DEC_COND_ARMIJO_GOLDSTEIN) {
                    oLineSearchParam.sublist(strLineSearchMethod).set("Sufficient Decrease Condition", "Armijo-Goldstein");
               } else if (uFlags & SUFFICIENT_DEC_COND_ARED_PRED) {
                    oLineSearchParam.sublist(strLineSearchMethod).set("Sufficient Decrease Condition", "Ared/Pred");
               }
               
               Teuchos::ParameterList& oLineSearchMethod = oLineSearchParam.sublist(strLineSearchMethod);
               oLineSearchMethod.set("Max Iters", iMaxIterLineSearch);
               oLineSearchMethod.set("Minimum Step", dMinStep);
               oLineSearchMethod.set("Recovery Step", dRecoveryStep);

               if (uFlags & RECOVERY_STEP_TYPE_CONST) {
                    oLineSearchMethod.set("Recovery Step Type", "Constant");
               } else if (uFlags & RECOVERY_STEP_TYPE_LAST_STEP) {
                    oLineSearchMethod.set("Recovery Step Type", "Last Computed Step");
               }
          } else if (uFlags & SOLVER_TRUST_REGION_BASED) {
               oSolverParam.set(szNonlinearSolver, "Trust Region Based");
          } else if (uFlags & SOLVER_INEXACT_TRUST_REGION_BASED) {
               oSolverParam.set(szNonlinearSolver, "Inexact Trust Region Based");
          } else if (uFlags & SOLVER_TENSOR_BASED) {
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
               oIterativeLinSol.SetAztecOption(AZ_solver, AZ_gmres);
          } else if (uFlags & LINEAR_SOLVER_CG) {
               oIterativeLinSol.SetAztecOption(AZ_solver, AZ_cg);
          } else if (uFlags & LINEAR_SOLVER_CGS) {
               oIterativeLinSol.SetAztecOption(AZ_solver, AZ_cgs);
          } else if (uFlags & LINEAR_SOLVER_TFQMR) {
               oIterativeLinSol.SetAztecOption(AZ_solver, AZ_tfqmr);
          } else if (uFlags & LINEAR_SOLVER_BICGSTAB) {
               oIterativeLinSol.SetAztecOption(AZ_solver, AZ_bicgstab);
          }

          oIterativeLinSol.SetAztecOption(AZ_kspace, std::min(Size, std::min(iMaxIterLinSol, iKrylovSubSpaceSize)));
          oIterativeLinSol.SetAztecOption(AZ_output, (uFlags & PRINT_CONVERGENCE_INFO) ? 1 : 0);
          oIterativeLinSol.SetOutputStream(std::cout);
          oIterativeLinSol.SetErrorStream(std::cerr);

          Teuchos::RCP<Epetra_Vector> pSolution = Teuchos::rcp(new Epetra_Vector(oMap, true));
          pSolutionView = Teuchos::rcp(new NOX::Epetra::Vector{pSolution, NOX::Epetra::Vector::CreateView});

          NOX::Epetra::Interface::Required& oResidualInt = *this;
          Teuchos::RCP<NOX::Epetra::Interface::Required> pResidualInt{Teuchos::rcpFromRef(oResidualInt)};

          if (uFlags & JACOBIAN_NEWTON_KRYLOV) {
               pJacOper = Teuchos::rcpFromRef(oMatFreeJacOper);
               pJacInt = &oMatFreeJacOper;
          } else {
               Epetra_Operator& oJacOper = *this;
               pJacOper = Teuchos::rcpFromRef(oJacOper);
               pJacInt = this;
          }

          pPrecInt = this;
          oIterativeLinSol.SetUserOperator(pJacOper.get());
          oIterativeLinSol.SetPrecOperator(this);

          NOX::Epetra::LinearSystem& oLinSys = *this;

          Teuchos::RCP<NOX::Epetra::Group> grpPtr =
               Teuchos::rcp(new NOX::Epetra::Group(oPrintParam,
                                                   pResidualInt,
                                                   *pSolutionView,
                                                   Teuchos::rcpFromRef(oLinSys)));

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

          ForcePrecondRebuild();

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

     void NoxNonlinearSolver::Residual(const VectorHandler* const pSol, VectorHandler* const pRes)
     {
          CPUTimeGuard oCPUTimeRes(*this, CPU_RESIDUAL);

          DeltaX.ScalarAddMul(*pSol, XPrev, -1.); // Convert to incremental solution

          XPrev = *pSol;

          pNonlinearProblem->Update(&DeltaX);

          pRes->Reset();

          VectorHandler* const pAbsRes = pGetResTest()->GetAbsRes();

          if (pAbsRes) {
               pAbsRes->Reset();
          }

          try {
               pNonlinearProblem->Residual(pRes, pAbsRes);
          } catch (const SolutionDataManager::ChangedEquationStructure& oErr) {
               DEBUGCERR("Caught exception change equation structure ...\n");
               
               if (bHonorJacRequest) {
                    DEBUGCERR("Force update of preconditioner ...\n");
                    
                    ForcePrecondRebuild();
               }
          }
     }

     void NoxNonlinearSolver::Jacobian()
     {
          CPUTimeGuard oCPUTimeJac(*this, CPU_JACOBIAN);

          bool bDone = false;
          MatrixHandler* pJac = nullptr;

          do {
               try {
                    pJac = pSolutionManager->pMatHdl();

                    pNonlinearProblem->Jacobian(pJac);

                    pJac->PacMat(); // Needed for Epetra_CrsMatrix only

                    bDone = true;
               } catch (const MatrixHandler::ErrRebuildMatrix& oErr) {
                    silent_cout("NoxNonlinearSolver: "
                                "rebuilding matrix...\n");
                    pSolutionManager->MatrInitialize();
               }
          } while (!bDone);

          ASSERT(pJac != nullptr);

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

          bUpdateJacobian = true;

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

          Residual(&oSol, &oRes);

          if (outputSol()) {
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

          if (bUpdateJacobian) {
               const MyVectorHandler oSol(Size, x.Values());
               Residual(&oSol, &TmpRes); // By convention AssRes must be called always before AssJac
               pSolutionManager->MatrReset();
               Jacobian();
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
          DEBUGCERR("ApplyInverse()\n");

          ++iPrecInnerIterCnt;
          ++iInnerIterCntTot;

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

     bool
     NoxNonlinearSolver::applyJacobian(const NOX::Epetra::Vector& input,
                                       NOX::Epetra::Vector& result) const
     {
          const bool bPrevUseTranspose = bUseTranspose;

          pJacOper->SetUseTranspose(false);

          integer status = pJacOper->Apply(input.getEpetraVector(), result.getEpetraVector());

          pJacOper->SetUseTranspose(bPrevUseTranspose);

          return status == 0;
     }

     bool
     NoxNonlinearSolver::applyJacobianTranspose(const NOX::Epetra::Vector& input,
                                                NOX::Epetra::Vector& result) const
     {
          const bool bPrevUseTranspose = pJacOper->UseTranspose();

          pJacOper->SetUseTranspose(true);

          integer status = pJacOper->Apply(input.getEpetraVector(), result.getEpetraVector());

          pJacOper->SetUseTranspose(bPrevUseTranspose);

          return status == 0;
     }

     bool
     NoxNonlinearSolver::applyJacobianInverse(Teuchos::ParameterList &params,
                                              const NOX::Epetra::Vector &input,
                                              NOX::Epetra::Vector &result)
     {
          const bool bPrevUseTranspose = UseTranspose();

          SetUseTranspose(false);

          if (uFlags & USE_PRECOND_AS_SOLVER) {
               integer status = ApplyInverse(input.getEpetraVector(), result.getEpetraVector());

               SetUseTranspose(bPrevUseTranspose);

               return status == 0;
          }

          oIterativeLinSol.SetLHS(&result.getEpetraVector());
          oIterativeLinSol.SetRHS(const_cast<Epetra_Vector*>(&input.getEpetraVector()));

          integer iMaxIter = params.get("Max Iterations", iMaxIterLinSol);
          doublereal dTol = params.get("Tolerance", dTolLinSol);

          iPrecInnerIterCnt = 0;

          integer status;

          bool bPrecReuse = iPrecInnerIterCntTot > 1;

          for (;;) {
               status = oIterativeLinSol.Iterate(iMaxIter, dTol);

               if (status == 0) {
                    DEBUGCERR("Linear solver converged\n");
                    break;
               }

               if (!bPrecReuse) {
                    DEBUGCERR("Linear solver failed to converge even with up to date preconditioner\n");
                    break;
               }

               DEBUGCERR("Linear solver failed to converge with non up to date preconditioner\n");
               DEBUGCERR("Preconditioner will be rebuild\n");

               const NOX::Abstract::Vector& XCurrA = pNonlinearSolver->getSolutionGroup().getX();
               const NOX::Epetra::Vector& XCurrE = dynamic_cast<const NOX::Epetra::Vector&>(XCurrA);

               recomputePreconditioner(XCurrE, oSolverParam);

               // Avoid starting the iterative solution with a residual too close to convergence
               // because it could cause a "loss of precision error"
               oIterativeLinSol.GetRHS()->PutScalar(0.);

               bPrecReuse = false;
          }

          oIterativeLinSol.UnsetLHSRHS();

          SetUseTranspose(bPrevUseTranspose);

          return status == 0;
     }

     bool
     NoxNonlinearSolver::applyRightPreconditioning(bool useTranspose,
                                                   Teuchos::ParameterList& params,
                                                   const NOX::Epetra::Vector& input,
                                                   NOX::Epetra::Vector& result) const
     {
          bool bPrevUseTranspose = UseTranspose();

          const_cast<NoxNonlinearSolver*>(this)->SetUseTranspose(useTranspose);

          integer status = ApplyInverse(input.getEpetraVector(), result.getEpetraVector());

          const_cast<NoxNonlinearSolver*>(this)->SetUseTranspose(bPrevUseTranspose);

          return status == 0;
     }

     Teuchos::RCP<NOX::Epetra::Scaling>
     NoxNonlinearSolver::getScaling()
     {
          return Teuchos::null;
     }

     void
     NoxNonlinearSolver::resetScaling(const Teuchos::RCP<NOX::Epetra::Scaling>& s)
     {
     }

     bool
     NoxNonlinearSolver::computeJacobian(const NOX::Epetra::Vector& x)
     {
          return pJacInt->computeJacobian(x.getEpetraVector(), *pJacOper);
     }

     bool
     NoxNonlinearSolver::createPreconditioner(const NOX::Epetra::Vector& x,
                                              Teuchos::ParameterList& p,
                                              bool recomputeGraph) const
     {
          return recomputePreconditioner(x, p);
     }

     bool
     NoxNonlinearSolver::destroyPreconditioner() const
     {
          return true;
     }

     bool
     NoxNonlinearSolver::recomputePreconditioner(const NOX::Epetra::Vector& x,
                                                 Teuchos::ParameterList& linearSolverParams) const
     {
          bool bStatus = pPrecInt->computePreconditioner(x.getEpetraVector(),
                                                         *oIterativeLinSol.GetPrecOperator(),
                                                         &linearSolverParams);

          if (bStatus) {
               ResetPrecondReuse();
          }

          return bStatus;
     }

     NOX::Epetra::LinearSystem::PreconditionerReusePolicyType
     NoxNonlinearSolver::getPreconditionerPolicy(bool advanceReuseCounter)
     {
          if (iPrecInnerIterCnt >= iInnerIterBeforeAssembly || iPrecInnerIterCntTot >= iIterationsBeforeAssembly) {
               return NOX::Epetra::LinearSystem::PRPT_RECOMPUTE;
          }

          if (advanceReuseCounter) {
               ++iPrecInnerIterCntTot;
          }

          return NOX::Epetra::LinearSystem::PRPT_REUSE;
     }

     bool
     NoxNonlinearSolver::isPreconditionerConstructed() const
     {
          return pPrecInt != nullptr;
     }

     Teuchos::RCP<const Epetra_Operator>
     NoxNonlinearSolver::getJacobianOperator() const
     {
          return pJacOper;
     }

     Teuchos::RCP<Epetra_Operator>
     NoxNonlinearSolver::getJacobianOperator()
     {
          return pJacOper;
     }

     Teuchos::RCP<Epetra_Operator>
     NoxNonlinearSolver::getGeneratedPrecOperator()
     {
          Epetra_Operator& oJac = *this;
          return Teuchos::rcpFromRef(oJac);
     }

     Teuchos::RCP<const Epetra_Operator>
     NoxNonlinearSolver::getGeneratedPrecOperator() const
     {
          const Epetra_Operator& oJac = *this;
          return Teuchos::rcpFromRef(oJac);
     }

     void
     NoxNonlinearSolver::setJacobianOperatorForSolve(const Teuchos::RCP<const Epetra_Operator>& solveJacOp)
     {
          DEBUGCERR("setJacobianOperatorForSolve()\n");

          ASSERT(solveJacOp.get() == this || solveJacOp.get() == &oMatFreeJacOper);
     }

     bool
     NoxNonlinearSolver::hasPreconditioner() const
     {
          return true;
     }

     void
     NoxNonlinearSolver::setPrecOperatorForSolve(const Teuchos::RCP<const Epetra_Operator>& solvePrecOp)
     {
          DEBUGCERR("setPrecOperatorForSolve()\n");

          ASSERT(solvePrecOp.get() == this);
     }

     void NoxNonlinearSolver::ResetPrecondReuse() const
     {
          iPrecInnerIterCnt = iPrecInnerIterCntTot = 0;
     }

     void NoxNonlinearSolver::ForcePrecondRebuild() const
     {
          iPrecInnerIterCnt = iPrecInnerIterCntTot = std::numeric_limits<integer>::max();
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
