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
 * The Naive Solver is copyright (C) 2004 by
 * Paolo Mantegazza <mantegazza@aero.polimi.it>
 */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <sys/types.h>
#include <unistd.h>
#include <signal.h>

#include "spmh.h"
#include "spmapmh.h"
#include "naivewrap.h"

extern "C" int
naivfct(void *mat_a, integer neq, integer *vet_nzr,
		void *mat_ri, integer *vet_nzc, void *mat_ci,
		integer ncd, integer *vet_piv, doublereal minpiv);
extern "C" void
naivslv(void *mat_a, integer neq, integer *vet_nzc,
		void *mat_ci, integer ncd, doublereal *vet_rhs,
		integer *vet_piv);

/* NaiveSolver - begin */
	
NaiveSolver::NaiveSolver(const integer &size)
: LinearSolver(0),
iSize(size),
Axp(0),
Arowp(0),
Acolp(0),
nzr(size),
nzc(size),
piv(size)
{
	NO_OP;
}

NaiveSolver::~NaiveSolver(void)
{
	NO_OP;
}

void
NaiveSolver::Init(void)
{
	bHasBeenReset = true;
}

void
NaiveSolver::Solve(void) const
{
	if (bHasBeenReset) {
      		((NaiveSolver *)this)->Factor();
      		bHasBeenReset = false;
	}

	naivslv(Axp, iSize, &nzc[0], Acolp, iSize,
			LinearSolver::pdRhs, &piv[0]);
}

void
NaiveSolver::Factor(void)
{
	int		rc;
	doublereal	minpiv = 1e-9;
	
	rc = naivfct(Axp, iSize, &nzr[0], Arowp, &nzc[0], Acolp,
			iSize, &piv[0], minpiv);

#define ENULCOL  -1
#define ENOPIV   -2

	switch (rc) {
	case ENULCOL:
		silent_cerr("NaiveSolver: ENULCOL" << std::endl);
		THROW(ErrGeneric());
		break;

	case ENOPIV:
		silent_cerr("NaiveSolver: ENOPIV" << std::endl);
		THROW(ErrGeneric());
		break;

	case 0:
		break;

	default:
		THROW(ErrGeneric());
		break;
	}
}

void
NaiveSolver::MakeCompactForm(SparseMatrixHandler& mh,
		std::vector<doublereal>& Ax,
		std::vector<integer>& Arow,
		std::vector<integer>& Acol,
		std::vector<integer>& Admy) const
{
	if (!bHasBeenReset) {
		return;
	}
	
	mh.MakeNaiveForm(Ax, Arow, Acol, nzr, nzc);

	Axp = &Ax[0];
	Arowp = &Arow[0];
	Acolp = &Acol[0];
}

/* NaiveSolver - end */

/* NaiveSparseSolutionManager - begin */

NaiveSparseSolutionManager::NaiveSparseSolutionManager(integer Dim)
: A(Dim),
VH(Dim)
{
	SAFENEWWITHCONSTRUCTOR(pLS, NaiveSolver, NaiveSolver(Dim));

	(void)pLS->ChangeResPoint(VH.pdGetVec());
	(void)pLS->ChangeSolPoint(VH.pdGetVec());
	pLS->SetSolutionManager(this);
}

NaiveSparseSolutionManager::~NaiveSparseSolutionManager(void) 
{
	NO_OP;
}

void
NaiveSparseSolutionManager::MatrReset(const doublereal& d)
{
	A.Reset(d);
}

void
NaiveSparseSolutionManager::MakeNaiveForm(void)
{
	pLS->MakeCompactForm(A, Ax, Arow, Acol, Admy);
}

void
NaiveSparseSolutionManager::MatrInit(const doublereal& d)
{
	MatrReset(d);
	pLS->Init();
}

/* Risolve il sistema  Fattorizzazione + Bacward Substitution*/
void
NaiveSparseSolutionManager::Solve(void)
{
	MakeNaiveForm();
	pLS->Solve();
}

/* Rende disponibile l'handler per la matrice */
MatrixHandler*
NaiveSparseSolutionManager::pMatHdl(void) const
{
	return &A;
}

/* Rende disponibile l'handler per il termine noto */
MyVectorHandler*
NaiveSparseSolutionManager::pResHdl(void) const
{
	return &VH;
}

/* Rende disponibile l'handler per la soluzione */
MyVectorHandler*
NaiveSparseSolutionManager::pSolHdl(void) const
{
	return &VH;
}

/* NaiveSparseSolutionManager - end */

