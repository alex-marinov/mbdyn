/* 
 * HmFe (C) is a FEM analysis code. 
 *
 * Copyright (C) 1996-2005
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
 * Copyright (C) 2001-2005
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
 * Copyright (C) 1996-2005
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
 * Lapack is used by permission; please read its Copyright,
 * License and Availability note.
 */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#ifdef USE_LAPACK
#include "solman.h"
#include "lapackwrap.h"

extern "C" {
int dgetrf_(integer *N, integer *N2, doublereal *A, integer *LDA,
		integer *IPIV, integer *INFO);
int dgetrs_(char *MODE, integer *N, integer *NRHS, doublereal *A,
		integer *LDA, integer *IPIV, doublereal *B, integer *LDB,
		integer *INFO);
}

/* LapackSolver - begin */
	
LapackSolver::LapackSolver(const integer &size, const doublereal &dPivot,
		doublereal *pa, doublereal *pb)
: LinearSolver(0),
iSize(size),
pA(pa),
pB(pb),
piIPIV(0)
{
	ASSERT(pA);

	SAFENEWARR(piIPIV, integer, iSize);
}

LapackSolver::~LapackSolver(void)
{
	if (piIPIV) {
		SAFEDELETEARR(piIPIV);
	}
}

void
LapackSolver::Reset(void)
{
	bHasBeenReset = true;
}

void
LapackSolver::Solve(void) const
{
	if (bHasBeenReset) {
      		((LapackSolver *)this)->Factor();
      		bHasBeenReset = false;
	}

	integer	iNRHS = 1, iINFO = 0;
	integer iN = iSize;

	dgetrs_("No transpose", &iN, &iNRHS, pA, &iN, piIPIV, pB, &iN, &iINFO);
}

void
LapackSolver::Factor(void)
{
	integer	iINFO = 0;

	dgetrf_(&iSize, &iSize, pA, &iSize, piIPIV, &iINFO);
}

/* LapackSolver - end */

/* LapackSolutionManager - begin */

LapackSolutionManager::LapackSolutionManager(integer Dim, doublereal dPivot)
: A(Dim),
VH(Dim)
{
	SAFENEWWITHCONSTRUCTOR(pLS, LapackSolver,
			LapackSolver(Dim, dPivot, A.pdGetMat(), VH.pdGetVec()));

	(void)pLS->ChangeResPoint(VH.pdGetVec());
	(void)pLS->ChangeSolPoint(VH.pdGetVec());

	pLS->SetSolutionManager(this);
}

LapackSolutionManager::~LapackSolutionManager(void) 
{
	NO_OP;
}

void
LapackSolutionManager::MatrReset(void)
{
	pLS->Reset();
}

/* Risolve il sistema  Fattorizzazione + Bacward Substitution*/
void
LapackSolutionManager::Solve(void)
{
	pLS->Solve();
}

/* Rende disponibile l'handler per la matrice */
MatrixHandler*
LapackSolutionManager::pMatHdl(void) const
{
	return &A;
}

/* Rende disponibile l'handler per il termine noto */
MyVectorHandler*
LapackSolutionManager::pResHdl(void) const
{
	return &VH;
}

/* Rende disponibile l'handler per la soluzione */
MyVectorHandler*
LapackSolutionManager::pSolHdl(void) const
{
	return &VH;
}

/* LapackSolutionManager - end */

#endif /* USE_LAPACK */

