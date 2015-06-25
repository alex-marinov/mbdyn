/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

/* solution manager */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <string.h>	/* for memset() */

#include <iostream>
#include <iomanip>

#include "solman.h"
#include "ls.h"

/* LinearSolver - begin */

LinearSolver::LinearSolver(SolutionManager *psm)
: pSM(psm), bHasBeenReset(true), pdRhs(0), pdSol(0)
{
	NO_OP;
}

LinearSolver::~LinearSolver(void)
{
	NO_OP;
}

#ifdef DEBUG
void
LinearSolver::IsValid(void) const
{
	ASSERT(pSM);
	ASSERT(pdRhs);
	ASSERT(pdSol);
}
#endif /* DEBUG */

void
LinearSolver::Reset(void)
{
	bHasBeenReset = true;
}

void
LinearSolver::SolveT(void) const
{
	silent_cerr("LinearSolver::SolveT() not supported" << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

void
LinearSolver::SetSolutionManager(SolutionManager *psm)
{
	pSM = psm;
}

/* ritorna il puntatore al vettore del residuo */
doublereal *
LinearSolver::pdGetResVec(void) const
{
	return pdRhs;
}

/* sposta il puntatore al vettore del residuo */
doublereal *
LinearSolver::pdSetResVec(doublereal* pd)
{
	doublereal *p = pdRhs;

	pdRhs = pd;

	return p;
}

/* ritorna il puntatore al vettore della soluzione */
doublereal *
LinearSolver::pdGetSolVec(void) const
{
	return pdSol;
}

/* sposta il puntatore al vettore della soluzione */
doublereal *
LinearSolver::pdSetSolVec(doublereal* pd)
{
	doublereal *p = pdSol;

	pdSol = pd;

	return p;
}

void
LinearSolver::MakeCompactForm(SparseMatrixHandler& mh,
		std::vector<doublereal>& Ax,
		std::vector<integer>& Ar,
		std::vector<integer>& Ac,
		std::vector<integer>& Ap) const
{
	NO_OP;
}

bool LinearSolver::bGetConditionNumber(doublereal& dCond)
{
	return false; // true means that the condition number was returned in dCond
}

/* LinearSolver - end */

