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

class AmesosPreconditioner: public Epetra_Operator {
public:
     AmesosPreconditioner(unsigned uPrecondFlag, const Teuchos::RCP<Epetra_RowMatrix>& pOperator);     
     virtual ~AmesosPreconditioner();

     const Teuchos::RCP<Epetra_RowMatrix>& GetOperator() const;

     int SymbolicFactorization();
     
     int NumericFactorization();
     
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
     static const char* GetSolverName(unsigned uPrecondFlag);
     
     Teuchos::RCP<Epetra_RowMatrix> pOperator;
     mutable Epetra_Vector oRhs;
     mutable Epetra_LinearProblem oProblem;
     Teuchos::RCP<Amesos_BaseSolver> pSolver;
};

class EpetraLinearSystem: public SolutionManager {
public:

     explicit EpetraLinearSystem(
#ifdef USE_MPI
          MPI::Intracomm& Comm,
#endif
          integer Dim);
     
     virtual ~EpetraLinearSystem(void);

#ifdef DEBUG
     virtual void IsValid(void) const override;
#endif

     virtual void MatrReset(void) override;

     virtual MatrixHandler* pMatHdl(void) const override;

     virtual VectorHandler* pResHdl(void) const override;

     virtual VectorHandler* pSolHdl(void) const override;

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

     virtual ~AztecOOSolutionManager(void);
     virtual void Solve() override;
     
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
     virtual void Solve() override;     
     virtual void MatrReset() override;
     virtual void MatrInitialize() override;

private:
     AmesosPreconditioner oPrecond;
     bool bRebuildSymbolic, bRebuildNumeric;
};

AmesosPreconditioner::AmesosPreconditioner(unsigned uPrecondFlag, const Teuchos::RCP<Epetra_RowMatrix>& pOperator)
     :pOperator(pOperator),
      oRhs(pOperator->Map(), false),
      oProblem(pOperator.get(), nullptr, &oRhs),
      pSolver(Teuchos::rcp(Amesos().Create(GetSolverName(uPrecondFlag), oProblem)))
{          
     if (!pSolver.get()) {
          throw Teuchos::ExceptionBase("Amesos failed to create linear solver interface");
     }
          
     pSolver->SetUseTranspose(pOperator->UseTranspose());
}
     
AmesosPreconditioner::~AmesosPreconditioner()
{
}

const Teuchos::RCP<Epetra_RowMatrix>& AmesosPreconditioner::GetOperator() const
{
     return pOperator;
}

int AmesosPreconditioner::SymbolicFactorization()
{
     return pSolver->SymbolicFactorization();
}

int AmesosPreconditioner::NumericFactorization()
{
     return pSolver->NumericFactorization();
}
     
int AmesosPreconditioner::SetUseTranspose(bool UseTranspose) 
{
     int ierr = pOperator->SetUseTranspose(UseTranspose);

     if (ierr != 0) {
          return ierr;
     }
          
     return pSolver->SetUseTranspose(UseTranspose);
}
     
int AmesosPreconditioner::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
     return pOperator->Apply(X, Y);
}
     
int AmesosPreconditioner::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
     ASSERT(X.GlobalLength() == pOperator->NumGlobalCols());
     ASSERT(Y.GlobalLength() == pOperator->NumGlobalRows());

     // AztecOO requires that the same object may be passed for X and Y
     std::copy(X.Values(), X.Values() + X.GlobalLength(), oRhs.Values());
          
     oProblem.SetLHS(&Y);

     int ierr = pSolver->Solve();
          
     oProblem.SetLHS(nullptr);

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
     ASSERT(pSolver->UseTranspose() == pOperator->UseTranspose());
          
     return pSolver->UseTranspose();
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

const char* AmesosPreconditioner::GetSolverName(unsigned uPrecondFlag)
{
     switch (uPrecondFlag) {
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

EpetraLinearSystem::~EpetraLinearSystem(void)
{
}

#ifdef DEBUG
void EpetraLinearSystem::IsValid(void) const
{
     A.IsValid();
     x.IsValid();
     b.IsValid();
}
#endif

void EpetraLinearSystem::MatrReset(void)
{
     A.Reset();
}

MatrixHandler* EpetraLinearSystem::pMatHdl(void) const
{
     return &A;
}

VectorHandler* EpetraLinearSystem::pResHdl(void) const
{
     return &b;
}

VectorHandler* EpetraLinearSystem::pSolHdl(void) const
{
     return &x;
}

bool EpetraLinearSystem::bGetConditionNumber(doublereal& dCond) const
{
     return false;
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
}

AztecOOSolutionManager::~AztecOOSolutionManager(void)
{
}

void AztecOOSolutionManager::Solve(void)
{
     A.PacMat();
    
     integer ierr = oSolver.Iterate(iMaxIter, dTol);

     if (ierr != 0) {
          silent_cerr("AztecOO error: iterative solution did not converge\n");
          throw LinearSolver::ErrFactor(-1, MBDYN_EXCEPT_ARGS);
     }
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
      oPrecond(uPrecondFlag, Teuchos::rcpFromRef(*A.pGetEpetraCrsMatrix())),
      bRebuildSymbolic(true),
      bRebuildNumeric(true)
{
     oSolver.SetPrecOperator(&oPrecond);
}

AztecOOPrecondSolutionManager::~AztecOOPrecondSolutionManager(void)
{
}

void AztecOOPrecondSolutionManager::Solve(void)
{
     A.PacMat();

     integer ierr = 0;

     do {
          if (bRebuildSymbolic) {
               ierr = oPrecond.SymbolicFactorization();

               if (ierr != 0) {
                    silent_cerr("Amesos error: symbolic factorization failed\n");               
                    throw LinearSolver::ErrFactor(-1, MBDYN_EXCEPT_ARGS);
               }
          }

          if (bRebuildNumeric) {
               ierr = oPrecond.NumericFactorization();

               if (ierr != 0) {
                    if (!bRebuildSymbolic) {
                         bRebuildSymbolic = true;
                    } else {
                         silent_cerr("Amesos error: numeric factorization failed\n");
                         throw LinearSolver::ErrFactor(-1, MBDYN_EXCEPT_ARGS);                         
                    }
               }
          }
     } while (ierr);

     bRebuildNumeric = false;     
     bRebuildSymbolic = false;
     
     AztecOOSolutionManager::Solve();
}

void AztecOOPrecondSolutionManager::MatrReset()
{
     bRebuildNumeric = true;
}
     
void AztecOOPrecondSolutionManager::MatrInitialize()
{
     bRebuildSymbolic = true;
     bRebuildNumeric = true;
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

#endif
