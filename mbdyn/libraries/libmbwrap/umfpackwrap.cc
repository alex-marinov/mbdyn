/* 
 * HmFe (C) is a FEM analysis code. 
 *
 * Copyright (C) 1996-2004
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
 * Copyright (C) 1996-2004
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
 * Copyright (C) 1996-2004
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#ifdef USE_UMFPACK
#include "solman.h"
#include "spmapmh.h"
#include "ccmh.h"
#include "umfpackwrap.h"

#ifdef HAVE_UMFPACK4

/*
 * Wrap of the required functions; (mostly) uses the calling convention 
 * of umfpack 3.0 so the matrix must be square
 */

#define UMFPACKWRAP_defaults 		umfpack_di_defaults
#define UMFPACKWRAP_free_symbolic 	umfpack_di_free_symbolic
#define UMFPACKWRAP_free_numeric 	umfpack_di_free_numeric

#ifdef HAVE_UMFPACK4_1
#define UMFPACKWRAP_symbolic(size, app, aip, axp, sym, ctrl, info) \
	umfpack_di_symbolic(size, size, app, aip, axp, sym, ctrl, info)
#else /* !HAVE_UMFPACK4_1 */
#define UMFPACKWRAP_symbolic(size, app, aip, axp, sym, ctrl, info) \
	umfpack_di_symbolic(size, size, app, aip, sym, ctrl, info)
#endif /* !HAVE_UMFPACK4_1 */

#define UMFPACKWRAP_report_info 	umfpack_di_report_info
#define UMFPACKWRAP_report_status 	umfpack_di_report_status
#define UMFPACKWRAP_numeric 		umfpack_di_numeric
#define UMFPACKWRAP_solve 		umfpack_di_solve

/* required factorization type (A * x = b) */
#define SYS_VALUE 			UMFPACK_A

#else /* !HAVE_UMFPACK4 */

/*
 * Wrap of the required functions; uses the calling convention 
 * of umfpack 3.0
 */
#define UMFPACKWRAP_defaults 		umfpack_defaults
#define UMFPACKWRAP_free_symbolic 	umfpack_free_symbolic
#define UMFPACKWRAP_free_numeric 	umfpack_free_numeric
#define UMFPACKWRAP_symbolic(size, app, aip, axp, sym, ctrl, info) \
	umfpack_symbolic(size, app, aip, sym, ctrl, info)
#define UMFPACKWRAP_report_info 	umfpack_report_info
#define UMFPACKWRAP_report_status 	umfpack_report_status
#define UMFPACKWRAP_numeric 		umfpack_numeric
#define UMFPACKWRAP_solve 		umfpack_solve

/* required factorization type (A * x = b) */
#define SYS_VALUE 			"Ax=b"

#endif /* HAVE_UMFPACK4 */

/* UmfpackSolver - begin */
	
UmfpackSolver::UmfpackSolver(const int &size, const doublereal &dPivot)
: iSize(size),
Axp(0),
Aip(0),
App(0),
Symbolic(0),
Numeric(0)
{
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
}

UmfpackSolver::~UmfpackSolver(void)
{
	UMFPACKWRAP_free_symbolic(&Symbolic);
	ASSERT(Symbolic == 0);

	UMFPACKWRAP_free_numeric(&Numeric);
	ASSERT(Numeric == 0);
}

void
UmfpackSolver::Init(void)
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
	status = UMFPACKWRAP_solve(SYS_VALUE,
			App, Aip, Axp, pdSol, pdRhs, 
			Numeric, Control, Info);
	if (status != UMFPACK_OK) {
		UMFPACKWRAP_report_info(Control, Info) ;
		UMFPACKWRAP_report_status(Control, status) ;
		silent_cerr("UMFPACKWRAP_solve failed" << std::endl);
		
		/* de-allocate memory */
		UMFPACKWRAP_free_numeric(&Numeric);
		ASSERT(Numeric == 0);

		THROW(ErrGeneric());
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
		THROW(ErrGeneric());
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
			THROW(ErrGeneric());
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

		THROW(ErrGeneric());
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
		std::vector<int>& Ai,
		std::vector<int>& Ac,
		std::vector<int>& Ap) const
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
		integer dummy, doublereal dPivot)
: A(Dim),
xVH(0),
bVH(0), 
x(Dim),
b(Dim),
pLS(0)
{
	SAFENEWWITHCONSTRUCTOR(pLS, UmfpackSolver, UmfpackSolver(Dim, dPivot));
	
	pLS->ChangeResPoint(&(b[0]));
	pLS->ChangeSolPoint(&(x[0]));
	pLS->SetSolutionManager(this);

	SAFENEWWITHCONSTRUCTOR(xVH, MyVectorHandler, 
			MyVectorHandler(Dim, &(x[0])));
	SAFENEWWITHCONSTRUCTOR(bVH, MyVectorHandler, 
			MyVectorHandler(Dim, &(b[0])));
}

UmfpackSparseSolutionManager::~UmfpackSparseSolutionManager(void) 
{
	SAFEDELETE(xVH);
	SAFEDELETE(bVH);
}

void
UmfpackSparseSolutionManager::MatrReset(const doublereal& d)
{
	A.Reset(d);
}

void
UmfpackSparseSolutionManager::MakeCompressedColumnForm(void)
{
	pLS->MakeCompactForm(A, Ax, Ai, Adummy, Ap);
}

void
UmfpackSparseSolutionManager::MatrInit(const doublereal& d)
{
	MatrReset(d);
	pLS->Init();
}

/* Risolve il sistema  Fattorizzazione + Bacward Substitution*/
void
UmfpackSparseSolutionManager::Solve(void)
{
	MakeCompressedColumnForm();
	pLS->Solve();
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
	return bVH;
}

/* Rende disponibile l'handler per la soluzione */
MyVectorHandler*
UmfpackSparseSolutionManager::pSolHdl(void) const
{
	return xVH;
}

/* UmfpackSparseSolutionManager - end */

/* UmfpackSparseCCSolutionManager - begin */

UmfpackSparseCCSolutionManager::UmfpackSparseCCSolutionManager(integer Dim,
		integer dummy, doublereal dPivot)
: UmfpackSparseSolutionManager(Dim, dummy, dPivot),
CCReady(false),
Ac(0)
{
	NO_OP;
}

UmfpackSparseCCSolutionManager::~UmfpackSparseCCSolutionManager(void) 
{
	if (Ac) {
		SAFEDELETE(Ac);
	}
}

void
UmfpackSparseCCSolutionManager::MatrReset(const doublereal& d)
{
	if (!CCReady) {
		A.Reset(d);
	} else {
		Ac->Reset(d);
	}
}

/* Risolve il sistema  Fattorizzazione + Bacward Substitution*/
void
UmfpackSparseCCSolutionManager::MakeCompressedColumnForm(void)
{
	if (!CCReady) {
		pLS->MakeCompactForm(A, Ax, Ai, Adummy, Ap);

		ASSERT(Ac == 0);

		typedef CColMatrixHandler<0> CC;
		
		SAFENEWWITHCONSTRUCTOR(Ac, CC, CC(Ax, Ai, Ap));

		CCReady = true;
	}
}

/* Inizializzatore "speciale" */
void
UmfpackSparseCCSolutionManager::MatrInitialize(const doublereal& d)
{
	SAFEDELETE(Ac);
	Ac = 0;

	CCReady = false;

	MatrInit();
}
	
/* Rende disponibile l'handler per la matrice */
MatrixHandler*
UmfpackSparseCCSolutionManager::pMatHdl(void) const
{
	if (!CCReady) {
		return &A;
	}

	ASSERT(Ac != 0);
	return Ac;
}

/* UmfpackSparseCCSolutionManager - end */

#endif /* USE_UMFPACK */

