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
#include "ls.h"
#include "linsol.h"
#ifdef USE_MPI
#include "mbcomm.h"
#endif

// FIXME: Compiler fails with -Werror if HAVE_BLAS is redefined inside a Trilinos header
#define HAVE_BLAS_SAVE HAVE_BLAS
#define HAVE_BOOL_SAVE HAVE_BOOL
#undef HAVE_BLAS
#undef HAVE_BOOL

#include "aztecoowrap.h"
#include "epetravh.h"
#include "epetraspmh.h"
#ifdef USE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <AztecOO.h>
#include <Amesos.h>
#include <Amesos_BaseSolver.h>

#undef HAVE_BLAS
#undef HAVE_BOOL
#define HAVE_BLAS HAVE_BLAS_SAVE
#define HAVE_BOOL HAVE_BOOL_SAVE
#undef HAVE_BLAS_SAVE
#undef HAVE_BOOL_SAVE

class AmesosSolver {
public:
     AmesosSolver(Epetra_RowMatrix* pA, Epetra_Vector* pX, Epetra_Vector* pB, unsigned uFlags);

     int Solve();

     void SetLHS(Epetra_MultiVector* pLhs) {
          oProblem.SetLHS(pLhs);
     }

     void SetRHS(Epetra_MultiVector* pRhs) {
          oProblem.SetRHS(pRhs);
     }

     int SymbolicFactorization() {
          return pSolver->SymbolicFactorization();
     }

     int NumericFactorization() {
          return pSolver->NumericFactorization();
     }

     bool SetUseTranspose(bool bTrans) {
          return pSolver->SetUseTranspose(bTrans);
     }

     bool UseTranspose() const {
          return pSolver->UseTranspose();
     }

     void MatrReset() {
          bRebuildNumeric = true;
     }

     void MatrInitialize() {
          bRebuildSymbolic = true;
          bRebuildNumeric = true;
     }

     static const char* GetSolverName(unsigned uFlags);

private:
     Epetra_LinearProblem oProblem;
     Teuchos::RCP<Amesos_BaseSolver> pSolver;
     bool bRebuildSymbolic, bRebuildNumeric;
};

class AmesosPreconditioner: public Epetra_Operator {
public:
     AmesosPreconditioner(unsigned uPrecondFlag, const Teuchos::RCP<Epetra_RowMatrix>& pOperator);
     virtual ~AmesosPreconditioner();

     const Teuchos::RCP<Epetra_RowMatrix>& GetOperator() const;

     void MatrReset() {
          oSolver.MatrReset();
     }

     void MatrInitialize() {
          oSolver.MatrInitialize();
     }

     virtual int SetUseTranspose(bool UseTranspose) override;

     virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override;

     virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override;

     virtual double NormInf() const override;

     virtual const char* Label() const override;

     virtual bool UseTranspose() const override;

     virtual bool HasNormInf() const override;

     virtual const Epetra_Comm & Comm() const override;

     virtual const Epetra_Map & OperatorDomainMap() const override;

     virtual const Epetra_Map & OperatorRangeMap() const override;

private:
     Teuchos::RCP<Epetra_RowMatrix> pOperator;
     mutable Epetra_Vector oRhs;
     mutable AmesosSolver oSolver;
};

class EpetraLinearSystem: public SolutionManager {
public:

     explicit EpetraLinearSystem(
#ifdef USE_MPI
          MPI::Intracomm& Comm,
#endif
          integer Dim);

     virtual ~EpetraLinearSystem();

#ifdef DEBUG
     virtual void IsValid() const override;
#endif

     virtual MatrixHandler* pMatHdl() const override;

     virtual VectorHandler* pResHdl() const override;

     virtual VectorHandler* pSolHdl() const override;

     virtual bool bGetConditionNumber(doublereal& dCond) const override;

protected:
#ifdef USE_MPI
     Epetra_MpiComm Comm;
#else
     Epetra_SerialComm Comm;
#endif
     mutable EpetraVectorHandler x, b;
     mutable EpetraSparseMatrixHandler A;
};

class AmesosSolutionManager: public EpetraLinearSystem {
public:
     AmesosSolutionManager(
#ifdef USE_MPI
          MPI::Intracomm& oComm,
#endif
          integer Dim,
          integer iVerbose,
          unsigned uFlags);

     virtual ~AmesosSolutionManager();
     virtual void MatrReset() override;
     virtual void MatrInitialize() override;
     virtual void Solve() override;

private:
     mutable AmesosSolver oSolver;
};

class AztecOOSolutionManager: public EpetraLinearSystem {
public:
     AztecOOSolutionManager(
#ifdef USE_MPI
          MPI::Intracomm& oComm,
#endif
          integer Dim,
          integer iMaxIter,
          doublereal dTol,
          integer iVerbose,
          unsigned uPrecondFlag);

     virtual ~AztecOOSolutionManager();
     virtual void Solve() override;
     virtual void MatrReset() override;

protected:
     Epetra_LinearProblem oProblem;
     AztecOO oSolver;

private:
     const integer iMaxIter;
     const doublereal dTol;
};

class AztecOOPrecondSolutionManager: public AztecOOSolutionManager {
public:
     AztecOOPrecondSolutionManager(
#ifdef USE_MPI
          MPI::Intracomm& oComm,
#endif
          integer Dim,
          integer iMaxIter,
          doublereal dTol,
          integer iVerbose,
          unsigned uPrecondFlag);

     virtual ~AztecOOPrecondSolutionManager();
     virtual void MatrReset() override;
     virtual void MatrInitialize() override;

private:
     AmesosPreconditioner oPrecond;
};

AmesosSolver::AmesosSolver(Epetra_RowMatrix* pA, Epetra_Vector* pX, Epetra_Vector* pB, unsigned uFlags)
     :oProblem(pA, pX, pB),
      pSolver(Teuchos::rcp(Amesos().Create(GetSolverName(uFlags), oProblem))),
      bRebuildSymbolic(true),
      bRebuildNumeric(true)
{
     if (!pSolver.get()) {
          throw Teuchos::ExceptionBase("Amesos failed to create linear solver interface");
     }
}

int AmesosSolver::Solve()
{
     DEBUGCERR("AmesosSolver::Solve()\n");

     integer ierr = 0;

     do {
          if (bRebuildSymbolic) {
               DEBUGCERR("AmesosSolver::SymbolicFactorization()\n");

               ierr = pSolver->SymbolicFactorization();

               if (ierr != 0) {
                    return ierr;
               }

               bRebuildSymbolic = false;
          }

          if (bRebuildNumeric) {
               DEBUGCERR("AmesosSolver::NumericFactorization()\n");

               ierr = pSolver->NumericFactorization();

               if (ierr != 0) {
                    if (!bRebuildSymbolic) {
                         bRebuildSymbolic = true;
                         continue;
                    }

                    return ierr;
               }

               bRebuildNumeric = false;
          }
     } while (ierr);

     ierr = pSolver->Solve();

     DEBUGCERR("AmesosSolver::Solve returned with status " << ierr << "\n");

     return ierr;
}

const char* AmesosSolver::GetSolverName(unsigned uFlags)
{
     switch (uFlags) {
     case LinSol::SOLVER_FLAGS_ALLOWS_PRECOND_UMFPACK:
          return "Umfpack";

     case LinSol::SOLVER_FLAGS_ALLOWS_PRECOND_KLU:
          return "Klu";

     case LinSol::SOLVER_FLAGS_ALLOWS_PRECOND_LAPACK:
          return "Lapack";

     default:
          throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
     }
}

AmesosPreconditioner::AmesosPreconditioner(unsigned uFlags, const Teuchos::RCP<Epetra_RowMatrix>& pOperator)
     :pOperator(pOperator),
      oRhs(pOperator->Map(), false),
      oSolver(pOperator.get(), nullptr, &oRhs, uFlags)
{
     oSolver.SetUseTranspose(pOperator->UseTranspose());
}

AmesosPreconditioner::~AmesosPreconditioner()
{
}

const Teuchos::RCP<Epetra_RowMatrix>& AmesosPreconditioner::GetOperator() const
{
     return pOperator;
}

int AmesosPreconditioner::SetUseTranspose(bool UseTranspose)
{
     int ierr = pOperator->SetUseTranspose(UseTranspose);

     if (ierr != 0) {
          return ierr;
     }

     return oSolver.SetUseTranspose(UseTranspose);
}

int AmesosPreconditioner::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
     return pOperator->Apply(X, Y);
}

int AmesosPreconditioner::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
     ASSERT(X.GlobalLength() == pOperator->NumGlobalCols());
     ASSERT(Y.GlobalLength() == pOperator->NumGlobalRows());
     ASSERT(X.GlobalLength() == oRhs.GlobalLength());

     // AztecOO requires that the same object may be passed for X and Y
     std::copy(X.Values(), X.Values() + X.GlobalLength(), oRhs.Values());

     oSolver.SetLHS(&Y);

     int ierr = oSolver.Solve();

     oSolver.SetLHS(nullptr);

     return ierr;
}

double AmesosPreconditioner::NormInf() const
{
     return -1.;
}

const char* AmesosPreconditioner::Label() const
{
     return pOperator->Label();
}

bool AmesosPreconditioner::UseTranspose() const
{
     ASSERT(oSolver.UseTranspose() == pOperator->UseTranspose());

     return pOperator->UseTranspose();
}

bool AmesosPreconditioner::HasNormInf() const
{
     return false;
}

const Epetra_Comm & AmesosPreconditioner::Comm() const
{
     return pOperator->Comm();
}

const Epetra_Map & AmesosPreconditioner::OperatorDomainMap() const
{
     return pOperator->OperatorDomainMap();
}

const Epetra_Map & AmesosPreconditioner::OperatorRangeMap() const
{
     return pOperator->OperatorRangeMap();
}

EpetraLinearSystem::EpetraLinearSystem(
#ifdef USE_MPI
     MPI::Intracomm& oComm,
#endif
     integer Dim)
     :
#ifdef USE_MPI
     Comm(oComm),
#endif
     x(Dim, Comm),
     b(Dim, Comm),
     A(Dim, Dim, 1, Comm)
{
}

EpetraLinearSystem::~EpetraLinearSystem()
{
}

#ifdef DEBUG
void EpetraLinearSystem::IsValid() const
{
     A.IsValid();
     x.IsValid();
     b.IsValid();
}
#endif

MatrixHandler* EpetraLinearSystem::pMatHdl() const
{
     return &A;
}

VectorHandler* EpetraLinearSystem::pResHdl() const
{
     return &b;
}

VectorHandler* EpetraLinearSystem::pSolHdl() const
{
     return &x;
}

bool EpetraLinearSystem::bGetConditionNumber(doublereal& dCond) const
{
     return false;
}

AmesosSolutionManager::AmesosSolutionManager(
#ifdef USE_MPI
     MPI::Intracomm& oComm,
#endif
     integer Dim,
     integer iVerbose,
     unsigned uFlags)
     :EpetraLinearSystem(
#ifdef USE_MPI
          oComm,
#endif
          Dim),
      oSolver(A.pGetEpetraCrsMatrix(),
              x.pGetEpetraVector(),
              b.pGetEpetraVector(),
              uFlags)
{
}

AmesosSolutionManager::~AmesosSolutionManager()
{
}

void AmesosSolutionManager::Solve()
{
     DEBUGCERR("AmesosSolutionManager::Solve()\n");

     int ierr = oSolver.Solve();

     if (ierr != 0) {
          silent_cerr("Amesos error: solution failed with status " << ierr << "\n");
          throw LinearSolver::ErrFactor(-1, MBDYN_EXCEPT_ARGS);
     }
}

void AmesosSolutionManager::MatrReset()
{
     oSolver.MatrReset();
}

void AmesosSolutionManager::MatrInitialize()
{
     oSolver.MatrInitialize();
}

AztecOOSolutionManager::AztecOOSolutionManager(
#ifdef USE_MPI
     MPI::Intracomm& oComm,
#endif
     integer Dim,
     integer iMaxIter,
     doublereal dTol,
     integer iVerbose,
     unsigned uPrecondFlag)
     :EpetraLinearSystem(
#ifdef USE_MPI
          oComm,
#endif
          Dim),
      oProblem(A.pGetEpetraCrsMatrix(), x.pGetEpetraVector(), b.pGetEpetraVector()),
      oSolver(oProblem),
      iMaxIter(iMaxIter),
      dTol(dTol)
{
     oSolver.SetAztecOption(AZ_output, iVerbose);
     oSolver.SetAztecOption(AZ_kspace, iMaxIter);
}

AztecOOSolutionManager::~AztecOOSolutionManager()
{
}

void AztecOOSolutionManager::Solve()
{
     integer ierr = oSolver.Iterate(iMaxIter, dTol);

#ifdef DEBUG
     MyVectorHandler Ax(b.iGetSize());

     A.MatVecMul(Ax, x);

     doublereal Axmb = 0.;
     doublereal Axpb = 0.;

     for (integer i = 1; i <= b.iGetSize(); ++i) {
          Axmb += std::pow(Ax(i) - b(i), 2);
          Axpb += std::pow(Ax(i) + b(i), 2);
     }

     DEBUGCERR("||A * x - b|| / ||A * x + b|| = " << sqrt(Axmb / Axpb) << "\n");
#endif

     if (ierr != 0) {
          silent_cerr("AztecOO error: iterative solution did not converge\n");
          throw LinearSolver::ErrFactor(-1, MBDYN_EXCEPT_ARGS);
     }
}

void AztecOOSolutionManager::MatrReset()
{
     // Must override because it is pure virtual
}

AztecOOPrecondSolutionManager::AztecOOPrecondSolutionManager(
#ifdef USE_MPI
     MPI::Intracomm& oComm,
#endif
     integer Dim,
     integer iMaxIter,
     doublereal dTol,
     integer iVerbose,
     unsigned uPrecondFlag)
     :AztecOOSolutionManager(
#ifdef USE_MPI
          oComm,
#endif
          Dim,
          iMaxIter,
          dTol,
          iVerbose,
          uPrecondFlag),
      oPrecond(uPrecondFlag, Teuchos::rcpFromRef(*A.pGetEpetraCrsMatrix()))
{
     oSolver.SetPrecOperator(&oPrecond);
}

AztecOOPrecondSolutionManager::~AztecOOPrecondSolutionManager()
{
}

void AztecOOPrecondSolutionManager::MatrReset()
{
     oPrecond.MatrReset();
}

void AztecOOPrecondSolutionManager::MatrInitialize()
{
     oPrecond.MatrInitialize();
}

SolutionManager*
pAllocateAztecOOSolutionManager(
#ifdef USE_MPI
     MPI::Intracomm& oComm,
#endif
     integer iNLD,
     integer iMaxIter,
     doublereal dTolRes,
     integer iVerbose,
     unsigned uSolverFlags)
{
     SolutionManager* pCurrSM = nullptr;

     unsigned uPrecondFlag = uSolverFlags & LinSol::SOLVER_FLAGS_PRECOND_MASK;

     switch (uPrecondFlag) {
     case LinSol::SOLVER_FLAGS_ALLOWS_PRECOND_ILUT:
          SAFENEWWITHCONSTRUCTOR(pCurrSM,
                                 AztecOOSolutionManager,
                                 AztecOOSolutionManager(
#ifdef USE_MPI
                                      oComm,
#endif
                                      iNLD,
                                      iMaxIter,
                                      dTolRes,
                                      iVerbose,
                                      uPrecondFlag));
          break;

     default:
          SAFENEWWITHCONSTRUCTOR(pCurrSM,
                                 AztecOOPrecondSolutionManager,
                                 AztecOOPrecondSolutionManager(
#ifdef USE_MPI
                                      oComm,
#endif
                                      iNLD,
                                      iMaxIter,
                                      dTolRes,
                                      iVerbose,
                                      uPrecondFlag));
     }

     return pCurrSM;
}

SolutionManager*
pAllocateAmesosSolutionManager(
#ifdef USE_MPI
     MPI::Intracomm& oComm,
#endif
     integer iNLD,
     integer iVerbose,
     unsigned uSolverFlags)
{
     SolutionManager* pCurrSM = nullptr;

     unsigned uPrecondFlag = uSolverFlags & LinSol::SOLVER_FLAGS_PRECOND_MASK;

     SAFENEWWITHCONSTRUCTOR(pCurrSM,
                            AmesosSolutionManager,
                            AmesosSolutionManager(
#ifdef USE_MPI
                                 oComm,
#endif
                                 iNLD,
                                 iVerbose,
                                 uPrecondFlag));
     return pCurrSM;
}

#endif
