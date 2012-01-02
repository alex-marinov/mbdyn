/* $Header$ */
/* 
 * HmFe (C) is a FEM analysis code. 
 *
 * Copyright (C) 1996-2012
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
 * Copyright (C) 2001-2012
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
 * Copyright (C) 1996-2012
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
 * Wsmp is used by permission; please read its Copyright,
 * License and Availability note.
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#ifdef USE_WSMP
#include "solman.h"
#include "spmapmh.h"
#include "ccmh.h"
#include "dirccmh.h"
#include "wsmpwrap.h"

extern "C" {
extern void wgsmp_(int *N,
		int *IA,
		int *JA,
		double *AVALS,
		double *B,
		int *LDB,
		int *NRHS,
		double *RMISC,
		int *IPARAM,
		double *DPARAM);
	extern void wsetmaxthrds_(int *);
	extern void wsmp_clear_();
}


/* WsmpSolver - begin */
	
WsmpSolver::WsmpSolver(const integer &size, const doublereal &dPivot,
		const unsigned blockSize, const unsigned nt)
: LinearSolver(0),
iSize(size),
Axp(0),
Aip(0),
App(0),
ldb(size),
nrhs(1),
rmisc(0),
Symbolic(false)
{
	int tnt = nt;
	wsetmaxthrds_(&tnt);
	iparm[0] = 0;
	iparm[1] = 0;
	iparm[2] = 0;
	
	wgsmp_(&iSize, 0, 0, 0, 0, &ldb, &nrhs, rmisc, iparm, dparm);

	/* CSC format */
	iparm[3] = 1;
	/* 0 index */
	iparm[4] = 0;
	/* pivot */

	if (dPivot != -1. && (dPivot >= 0. && dPivot <= 1.)) {
		/*
		 * 1.0: true partial pivoting
		 * 0.0: treated as 1.0
		 * 
		 * default: 0.1
		 */
		dparm[10] = dPivot;
		dparm[21] = dPivot;
	}

	if (blockSize > 0) {
		iparm[31] = blockSize;
	}
}

WsmpSolver::~WsmpSolver(void)
{
	wsmp_clear_();
}

void
WsmpSolver::Reset(void)
{
	bHasBeenReset = true;
}

void
WsmpSolver::Solve(void) const
{
	if (bHasBeenReset) {
      		((WsmpSolver *)this)->Factor();
      		bHasBeenReset = false;
	}
		
	int status;

	/*
	 * NOTE: Axp, Aip, App should have been set by * MakeCompactForm()
	 */
		
	iparm[1] = 3;
	iparm[2] = 3;
	
	int tsize = iSize;
	
	wgsmp_(&tsize, App, Aip, Axp, LinearSolver::pdRhs, &ldb, &nrhs, rmisc, iparm, dparm);
	
	status = iparm[63];
	
	if (status != 0) {
		silent_cerr("Wsmp back substitution failed" << std::endl);
		
		/* de-allocate memory */
		wsmp_clear_();

		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	
	iparm[1] = 4;
	iparm[2] = 4;
	
	wgsmp_(&tsize, App, Aip, Axp, LinearSolver::pdRhs, &ldb, &nrhs, rmisc, iparm, dparm);
	
	status = iparm[63];
	
	if (status != 0) {
		silent_cerr("Wsmp iterative refinement failed" << std::endl);
		
		/* de-allocate memory */
		wsmp_clear_();

		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	
}

void
WsmpSolver::Factor(void)
{
	int status;

	/*
	 * NOTE: Axp, Aip, App should have been set by * MakeCompactForm()
	 */
		
	if (!Symbolic && !bPrepareSymbolic()) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	iparm[1] = 2;
	iparm[2] = 2;

	int tsize = iSize;
	
	wgsmp_(&tsize, App, Aip, Axp, LinearSolver::pdRhs, &ldb, &nrhs, rmisc, iparm, dparm);
	
	status = iparm[63];

	if (status != 0) {
		silent_cerr("Wsmp factorization failed" << std::endl);
		
		/* de-allocate memory */
		wsmp_clear_();

		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

void
WsmpSolver::MakeCompactForm(SparseMatrixHandler& mh,
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
WsmpSolver::bPrepareSymbolic(void)
{
	int status;

	iparm[1] = 1;
	iparm[2] = 1;

	int tsize = iSize;
	
	wgsmp_(&tsize, App, Aip, Axp, LinearSolver::pdRhs, &ldb, &nrhs, rmisc, iparm, dparm);
	
	status = iparm[63];

	if (status != 0) {
		silent_cerr("Wsmp factorization failed" << std::endl);
		
		/* de-allocate memory */
		wsmp_clear_();

		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	Symbolic = true;

	return true;
}

/* WsmpSolver - end */

/* WsmpSparseSolutionManager - begin */

WsmpSparseSolutionManager::WsmpSparseSolutionManager(integer Dim,
		doublereal dPivot,
		const unsigned blockSize,
		const unsigned nt)
: A(Dim),
xb(Dim),
xbVH(Dim, &xb[0])
{
	SAFENEWWITHCONSTRUCTOR(pLS, WsmpSolver,
			WsmpSolver(Dim, dPivot, blockSize, nt));

	(void)pLS->pdSetResVec(&xb[0]);
	(void)pLS->pdSetSolVec(&xb[0]);
	pLS->SetSolutionManager(this);
}

WsmpSparseSolutionManager::~WsmpSparseSolutionManager(void) 
{
	NO_OP;
}

void
WsmpSparseSolutionManager::MatrReset(void)
{
	pLS->Reset();
}

void
WsmpSparseSolutionManager::MakeCompressedColumnForm(void)
{
	pLS->MakeCompactForm(A, Ax, Ai, Adummy, Ap);
}

/* Risolve il sistema  Fattorizzazione + Bacward Substitution*/
void
WsmpSparseSolutionManager::Solve(void)
{
	MakeCompressedColumnForm();
	pLS->Solve();
}

/* Rende disponibile l'handler per la matrice */
MatrixHandler*
WsmpSparseSolutionManager::pMatHdl(void) const
{
	return &A;
}

/* Rende disponibile l'handler per il termine noto */
MyVectorHandler*
WsmpSparseSolutionManager::pResHdl(void) const
{
	return &xbVH;
}

/* Rende disponibile l'handler per la soluzione */
MyVectorHandler*
WsmpSparseSolutionManager::pSolHdl(void) const
{
	return &xbVH;
}

/* WsmpSparseSolutionManager - end */

template <class CC>
WsmpSparseCCSolutionManager<CC>::WsmpSparseCCSolutionManager(integer Dim,
		doublereal dPivot,
		const unsigned& blockSize,
		const unsigned nt)
: WsmpSparseSolutionManager(Dim, dPivot, blockSize, nt),
CCReady(false),
Ac(0)
{
	NO_OP;
}

template <class CC>
WsmpSparseCCSolutionManager<CC>::~WsmpSparseCCSolutionManager(void) 
{
	if (Ac) {
		SAFEDELETE(Ac);
	}
}

template <class CC>
void
WsmpSparseCCSolutionManager<CC>::MatrReset(void)
{
	pLS->Reset();
}

template <class CC>
void
WsmpSparseCCSolutionManager<CC>::MakeCompressedColumnForm(void)
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
WsmpSparseCCSolutionManager<CC>::MatrInitialize()
{
	CCReady = false;

	MatrReset();
}
	
/* Rende disponibile l'handler per la matrice */
template <class CC>
MatrixHandler*
WsmpSparseCCSolutionManager<CC>::pMatHdl(void) const
{
	if (!CCReady) {
		return &A;
	}

	ASSERT(Ac != 0);
	return Ac;
}

template class WsmpSparseCCSolutionManager<CColMatrixHandler<0> >;
template class WsmpSparseCCSolutionManager<DirCColMatrixHandler<0> >;

/* WsmpSparseCCSolutionManager - end */

#endif /* USE_WSMP */

