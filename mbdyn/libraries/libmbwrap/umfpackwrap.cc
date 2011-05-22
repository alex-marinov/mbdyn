/* $Header$ */
/* 
 * HmFe (C) is a FEM analysis code. 
 *
 * Copyright (C) 1996-2011
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
 * Copyright (C) 2001-2011
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
 * Copyright (C) 1996-2011
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
	const unsigned blockSize)
: LinearSolver(0),
iSize(size),
Axp(0),
Aip(0),
App(0),
Symbolic(0),
Numeric(0)
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
      		((UmfpackSolver *)this)->Factor();
      		bHasBeenReset = false;
	}
		
	int status;

	/*
	 * NOTE: Axp, Aip, App should have been set by * MakeCompactForm()
	 */
		
#ifdef UMFPACK_REPORT
	doublereal t = t_iniz;
#endif /* UMFPACK_REPORT */

	Control[UMFPACK_IRSTEP] = 0;
	status = UMFPACKWRAP_solve((bTranspose ? SYS_VALUET : SYS_VALUE),
			App, Aip, Axp,
			LinearSolver::pdSol, LinearSolver::pdRhs, 
			Numeric, Control, Info);
	if (status != UMFPACK_OK) {
		UMFPACKWRAP_report_info(Control, Info) ;
		UMFPACKWRAP_report_status(Control, status) ;
		silent_cerr("UMFPACKWRAP_solve failed" << std::endl);
		
		/* de-allocate memory */
		UMFPACKWRAP_free_numeric(&Numeric);
		ASSERT(Numeric == 0);

		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
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

/* UmfpackSolver - end */

/* UmfpackSparseSolutionManager - begin */

UmfpackSparseSolutionManager::UmfpackSparseSolutionManager(integer Dim,
		doublereal dPivot,
		doublereal dDropTolerance,
		const unsigned blockSize)
: A(Dim),
x(Dim),
b(Dim),
xVH(Dim, &x[0]),
bVH(Dim, &b[0])
{
	SAFENEWWITHCONSTRUCTOR(pLS, UmfpackSolver,
			UmfpackSolver(Dim, dPivot, dDropTolerance, blockSize));

	(void)pLS->pdSetResVec(&b[0]);
	(void)pLS->pdSetSolVec(&x[0]);
	pLS->SetSolutionManager(this);
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
	pLS->MakeCompactForm(A, Ax, Ai, Adummy, Ap);
}

/* Risolve il sistema Fattorizzazione + Bacward Substitution */
void
UmfpackSparseSolutionManager::Solve(void)
{
	MakeCompressedColumnForm();
	pLS->Solve();
}

/* Risolve il sistema (trasposto) Fattorizzazione + Bacward Substitution */
void
UmfpackSparseSolutionManager::SolveT(void)
{
	MakeCompressedColumnForm();
	pLS->SolveT();
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
		const unsigned& blockSize)
: UmfpackSparseSolutionManager(Dim, dPivot, dDropTolerance, blockSize),
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

