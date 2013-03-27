/* $Header$ */
/* 
 * HmFe (C) is a FEM analysis code. 
 *
 * Copyright (C) 1996-2013
 *
 * Marco Morandini  <morandini@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */
/* December 2001 
 * Modified to add a Sparse matrix in row form and to implement methods
 * to be used in the parallel MBDyn Solver.
 *
 * Copyright (C) 2001-2013
 *
 * Giuseppe Quaranta  <quaranta@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *      
 */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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
 * Umfpack is used by permission; please read its Copyright,
 * License and Availability note.
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#ifdef USE_UMFPACK
#include "solman.h"
#include "spmapmh.h"
#include "ccmh.h"
#include "dirccmh.h"
#include "umfpackwrap.h"
#include "dgeequ.h"
#include <cstring>

/* in some cases, int and int32_t differ */

#ifdef USE_UMFPACK_LONG

#define UMFPACKWRAP_defaults 		umfpack_dl_defaults
#define UMFPACKWRAP_free_symbolic 	umfpack_dl_free_symbolic
#define UMFPACKWRAP_free_numeric 	umfpack_dl_free_numeric

#define UMFPACKWRAP_symbolic(size, app, aip, axp, sym, ctrl, info) \
	umfpack_dl_symbolic(size, size, app, aip, axp, sym, ctrl, info)

#define UMFPACKWRAP_report_info 	umfpack_dl_report_info
#define UMFPACKWRAP_report_status 	umfpack_dl_report_status
#define UMFPACKWRAP_numeric 		umfpack_dl_numeric
#define UMFPACKWRAP_solve 		umfpack_dl_solve

#else // ! USE_UMFPACK_LONG

#define UMFPACKWRAP_defaults 		umfpack_di_defaults
#define UMFPACKWRAP_free_symbolic 	umfpack_di_free_symbolic
#define UMFPACKWRAP_free_numeric 	umfpack_di_free_numeric

#define UMFPACKWRAP_symbolic(size, app, aip, axp, sym, ctrl, info) \
	umfpack_di_symbolic(size, size, app, aip, axp, sym, ctrl, info)

#define UMFPACKWRAP_report_info 	umfpack_di_report_info
#define UMFPACKWRAP_report_status 	umfpack_di_report_status
#define UMFPACKWRAP_numeric 		umfpack_di_numeric
#define UMFPACKWRAP_solve 		umfpack_di_solve

#endif // ! USE_UMFPACK_LONG

/* required factorization type (A * x = b) */
#define SYS_VALUE 			UMFPACK_A
#define SYS_VALUET 			UMFPACK_Aat

/* UmfpackSolver - begin */
	
UmfpackSolver::UmfpackSolver(const integer &size,
	const doublereal &dPivot,
	const doublereal &dDropTolerance,
	const unsigned blockSize,
	Scale scale,
	integer iMaxIter)
: LinearSolver(0),
iSize(size),
Axp(0),
Aip(0),
App(0),
Symbolic(0),
Numeric(0),
bHaveCond(false)
{
	// silence static analyzers
	memset(&Info[0], 0, sizeof(Info));
	UMFPACKWRAP_defaults(Control);

	if (dPivot != -1. && (dPivot >= 0. && dPivot <= 1.)) {
		/*
		 * 1.0: true partial pivoting
		 * 0.0: treated as 1.0
		 * 
		 * default: 0.1
		 */
		Control[UMFPACK_PIVOT_TOLERANCE] = dPivot;
	}

	if (dDropTolerance != 0.) {
#ifdef UMFPACK_DROPTOL
		ASSERT(dDropTolerance > 0.);
		Control[UMFPACK_DROPTOL] = dDropTolerance;
#endif
	}

	if (blockSize > 0) {
		Control[UMFPACK_BLOCK_SIZE] = blockSize;
	}

	if (scale != SCALE_UNDEF) {
		Control[UMFPACK_SCALE] = scale;
	}

	if (iMaxIter >= 0) {
		Control[UMFPACK_IRSTEP] = iMaxIter;
	}
}

UmfpackSolver::~UmfpackSolver(void)
{
	UMFPACKWRAP_free_symbolic(&Symbolic);
	ASSERT(Symbolic == 0);

	UMFPACKWRAP_free_numeric(&Numeric);
	ASSERT(Numeric == 0);
}

void
UmfpackSolver::Reset(void)
{
	bHaveCond = false;

	if (Numeric) {
		UMFPACKWRAP_free_numeric(&Numeric);
		ASSERT(Numeric == 0);
	}

	bHasBeenReset = true;
}

void
UmfpackSolver::Solve(void) const
{
	Solve(false);
}

void
UmfpackSolver::SolveT(void) const
{
	Solve(true);
}

void
UmfpackSolver::Solve(bool bTranspose) const
{
	if (bHasBeenReset) {
      		const_cast<UmfpackSolver *>(this)->Factor();
      		bHasBeenReset = false;
	}
		
	int status;

	/*
	 * NOTE: Axp, Aip, App should have been set by * MakeCompactForm()
	 */
		
#ifdef UMFPACK_REPORT
	doublereal t = t_iniz;
#endif /* UMFPACK_REPORT */

	status = UMFPACKWRAP_solve((bTranspose ? SYS_VALUET : SYS_VALUE),
			App, Aip, Axp,
			LinearSolver::pdSol, LinearSolver::pdRhs, 
			Numeric, Control, Info);

	if (status != UMFPACK_OK) {
		UMFPACKWRAP_report_info(Control, Info);
		UMFPACKWRAP_report_status(Control, status) ;
		silent_cerr("UMFPACKWRAP_solve failed" << std::endl);
		
		/* de-allocate memory */
		UMFPACKWRAP_free_numeric(&Numeric);
		ASSERT(Numeric == 0);

		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (Control[UMFPACK_IRSTEP] > 0 && Info[UMFPACK_IR_TAKEN] >= Control[UMFPACK_IRSTEP]) {
		silent_cerr("warning: UMFPACK_IR_TAKEN = " << Info[UMFPACK_IR_TAKEN]
		             << " >= UMFPACK_IRSTEP = " << Control[UMFPACK_IRSTEP] <<  std::endl
		             << "\tThe flag \"max iterations\" of the linear solver is too small or the condition number of the Jacobian matrix is too high" << std::endl);
	}

#ifdef UMFPACK_REPORT
	UMFPACKWRAP_report_info(Control, Info);
	doublereal t1 = umfpack_timer() - t;	/* ?!? not used! */
#endif /* UMFPACK_REPORT */
}

void
UmfpackSolver::Factor(void)
{
	int status;

	/*
	 * NOTE: Axp, Aip, App should have been set by * MakeCompactForm()
	 */
		
	if (Symbolic == 0 && !bPrepareSymbolic()) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

#ifdef UMFPACK_REPORT
	UMFPACKWRAP_report_symbolic ("Symbolic factorization of A",
			Symbolic, Control) ;
	UMFPACKWRAP_report_info(Control, Info);
	doublereal t1 = umfpack_timer() - t;	/* ?!? not used! */
#endif /* UMFPACK_REPORT */

	status = UMFPACKWRAP_numeric(App, Aip, Axp, Symbolic, 
			&Numeric, Control, Info);
	if (status == UMFPACK_ERROR_different_pattern) {
		UMFPACKWRAP_free_symbolic(&Symbolic);
		if (!bPrepareSymbolic()) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		status = UMFPACKWRAP_numeric(App, Aip, Axp, Symbolic, 
				&Numeric, Control, Info);
	}

	if (status != UMFPACK_OK) {
		UMFPACKWRAP_report_info(Control, Info);
		UMFPACKWRAP_report_status(Control, status);
		silent_cerr("UMFPACKWRAP_numeric failed" << std::endl);

		/* de-allocate memory */
		UMFPACKWRAP_free_symbolic(&Symbolic);
		UMFPACKWRAP_free_numeric(&Numeric);
		ASSERT(Numeric == 0);

		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	} else {
		bHaveCond = true;
	}
		
#ifdef UMFPACK_REPORT
	UMFPACKWRAP_report_numeric ("Numeric factorization of A",
			Numeric, Control);
	UMFPACKWRAP_report_info(Control, Info);
	t1 = umfpack_timer() - t;	/* ?!? not used! */
#endif /* UMFPACK_REPORT */
}

void
UmfpackSolver::MakeCompactForm(SparseMatrixHandler& mh,
		std::vector<doublereal>& Ax,
		std::vector<integer>& Ai,
		std::vector<integer>& Ac,
		std::vector<integer>& Ap) const
{
	if (!bHasBeenReset) {
		return;
	}
	
	mh.MakeCompressedColumnForm(Ax, Ai, Ap, 0);

	Axp = &(Ax[0]);
	Aip = &(Ai[0]);
	App = &(Ap[0]);
}

bool 
UmfpackSolver::bPrepareSymbolic(void)
{
	int status;

	status = UMFPACKWRAP_symbolic(iSize, App, Aip, Axp,
			&Symbolic, Control, Info);
	if (status != UMFPACK_OK) {
		UMFPACKWRAP_report_info(Control, Info) ;
		UMFPACKWRAP_report_status(Control, status);
		silent_cerr("UMFPACKWRAP_symbolic failed" << std::endl);

		/* de-allocate memory */
		UMFPACKWRAP_free_symbolic(&Symbolic);
		ASSERT(Symbolic == 0);

		return false;
	}

	return true;
}

bool UmfpackSolver::bGetConditionNumber(doublereal& dCond)
{
	if (!bHaveCond) {
		return false;
	}

	dCond = 1. / Info[UMFPACK_RCOND];

	return true;
}

/* UmfpackSolver - end */

/* UmfpackSparseSolutionManager - begin */

UmfpackSparseSolutionManager::UmfpackSparseSolutionManager(integer Dim,
		doublereal dPivot,
		doublereal dDropTolerance,
		const unsigned blockSize,
		UmfpackSolver::Scale scale,
		integer iMaxIter,
		ScaleWhen ms)
: A(Dim),
x(Dim),
b(Dim),
xVH(Dim, &x[0]),
bVH(Dim, &b[0]),
ms(ms)
{
	SAFENEWWITHCONSTRUCTOR(pLS, UmfpackSolver,
			UmfpackSolver(Dim, dPivot, dDropTolerance, blockSize, scale, iMaxIter));

	(void)pLS->pdSetResVec(&b[0]);
	(void)pLS->pdSetSolVec(&x[0]);
	pLS->SetSolutionManager(this);

	if (scale != UmfpackSolver::SCALE_NONE && ms != NEVER) {
		ASSERT(0); // avoid double scaling
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

UmfpackSparseSolutionManager::~UmfpackSparseSolutionManager(void) 
{
	NO_OP;
}

void
UmfpackSparseSolutionManager::MatrReset(void)
{
	pLS->Reset();
}

void
UmfpackSparseSolutionManager::MakeCompressedColumnForm(void)
{
	ScaleMatrixAndRightHandSide<SpMapMatrixHandler>(A);

	pLS->MakeCompactForm(A, Ax, Ai, Adummy, Ap);
}

template <class MH>
void UmfpackSparseSolutionManager::ScaleMatrixAndRightHandSide(MH &mh)
{
#if defined(DEBUG)
	static bool bPrintMat = true;
#endif

	if (ms != SolutionManager::NEVER) {
		if (pLS->bReset()) {
			// FIXME: if matrix is CC or Dir and the matrix was regenerated,
			// the scaling needs to be recomputed
			if (msr.empty() || ms == SolutionManager::ALWAYS) {
#if defined(DEBUG)
				if (bPrintMat) {
					silent_cout("% matrix before scaling:\n");
					silent_cout("Apre=[");
					for (typename MH::const_iterator i = mh.begin(); i != mh.end(); ++i) {
						if (i->dCoef != 0.) {
							silent_cout(i->iRow + 1 << ", " << i->iCol + 1 << ", " << i->dCoef << ";\n");
						}
					}
					silent_cout("];\n");
				}
#endif
				// (re)compute
				doublereal rowcnd = -1., colcnd = -1., amax = -1.;
				dgeequ<MH>(mh, msr, msc, rowcnd, colcnd, amax);

#if defined(DEBUG)
				if (amax < std::numeric_limits<doublereal>::epsilon()
					|| amax > 1./std::numeric_limits<doublereal>::epsilon())
				{
					silent_cerr("Warning: The matrix should be scaled\n");
				}

				if (colcnd >= 0.1) {
					silent_cerr("Warning: it is not worth scaling by C\n");
				}

				if (rowcnd >= 0.1
					&& amax >= std::numeric_limits<doublereal>::epsilon()
					&& amax <= 1./std::numeric_limits<doublereal>::epsilon())
				{
					silent_cerr("Warning: it is not worth scaling by R\n");
				}

				if (bPrintMat) {
					silent_cout("% scale factors:\n");
					silent_cout("msr=[");
					for (std::vector<doublereal>::const_iterator i = msr.begin(); i < msr.end(); ++i) {
						silent_cout(*i << ";\n");
					}

					silent_cout("];\n");
					silent_cout("msc=[");
					for (std::vector<doublereal>::const_iterator i = msc.begin(); i < msc.end(); ++i) {
						silent_cout(*i << ";\n");
					}
					silent_cout("];\n");
				}
#endif
			}

			// in any case scale matrix and right-hand-side
			dgeequ_scale<MH>(mh, msr, msc);

#if defined(DEBUG)
			if (bPrintMat) {
				silent_cout("% matrix after scaling:\n");
				silent_cout("Apost=[");
				for (typename MH::const_iterator i = mh.begin(); i != mh.end(); ++i) {
					if (i->dCoef != 0.) {
						silent_cout(i->iRow + 1 << ", " << i->iCol + 1 << ", " << i->dCoef << ";\n");
					}
				}
				silent_cout("];\n");
			}

			bPrintMat = false;
#endif
		}
		dgeequ_scale(bVH, &msr[0]);
	}
}

void UmfpackSparseSolutionManager::ScaleSolution(void)
{
	if (ms != SolutionManager::NEVER) {
		// scale solution
		dgeequ_scale(xVH, &msc[0]);
	}
}

/* Risolve il sistema Fattorizzazione + Backward Substitution */
void
UmfpackSparseSolutionManager::Solve(void)
{
	MakeCompressedColumnForm();

	pLS->Solve();

	ScaleSolution();
}

/* Risolve il sistema (trasposto) Fattorizzazione + Backward Substitution */
void
UmfpackSparseSolutionManager::SolveT(void)
{
	MakeCompressedColumnForm();

	pLS->SolveT();

	ScaleSolution();
}

/* Rende disponibile l'handler per la matrice */
MatrixHandler*
UmfpackSparseSolutionManager::pMatHdl(void) const
{
	return &A;
}

/* Rende disponibile l'handler per il termine noto */
MyVectorHandler*
UmfpackSparseSolutionManager::pResHdl(void) const
{
	return &bVH;
}

/* Rende disponibile l'handler per la soluzione */
MyVectorHandler*
UmfpackSparseSolutionManager::pSolHdl(void) const
{
	return &xVH;
}

/* UmfpackSparseSolutionManager - end */

template <class CC>
UmfpackSparseCCSolutionManager<CC>::UmfpackSparseCCSolutionManager(integer Dim,
		doublereal dPivot,
		doublereal dDropTolerance,
		const unsigned& blockSize,
		UmfpackSolver::Scale scale,
		integer iMaxIter,
		ScaleWhen ms)
: UmfpackSparseSolutionManager(Dim, dPivot, dDropTolerance, blockSize, scale, iMaxIter, ms),
CCReady(false),
Ac(0)
{
	NO_OP;
}

template <class CC>
UmfpackSparseCCSolutionManager<CC>::~UmfpackSparseCCSolutionManager(void) 
{
	if (Ac) {
		SAFEDELETE(Ac);
	}
}

template <class CC>
void
UmfpackSparseCCSolutionManager<CC>::MatrReset(void)
{
	pLS->Reset();
}

template <class CC>
void
UmfpackSparseCCSolutionManager<CC>::MakeCompressedColumnForm(void)
{
	if (!CCReady) {
		pLS->MakeCompactForm(A, Ax, Ai, Adummy, Ap);

		if (Ac == 0) {
			SAFENEWWITHCONSTRUCTOR(Ac, CC, CC(Ax, Ai, Ap));
		}

		CCReady = true;
	}

	ScaleMatrixAndRightHandSide<CC>(*dynamic_cast<CC *>(Ac));
}

/* Inizializzatore "speciale" */
template <class CC>
void
UmfpackSparseCCSolutionManager<CC>::MatrInitialize()
{
	CCReady = false;

	MatrReset();
}
	
/* Rende disponibile l'handler per la matrice */
template <class CC>
MatrixHandler*
UmfpackSparseCCSolutionManager<CC>::pMatHdl(void) const
{
	if (!CCReady) {
		return &A;
	}

	ASSERT(Ac != 0);
	return Ac;
}

template class UmfpackSparseCCSolutionManager<CColMatrixHandler<0> >;
template class UmfpackSparseCCSolutionManager<DirCColMatrixHandler<0> >;

/* UmfpackSparseCCSolutionManager - end */

#endif /* USE_UMFPACK */

