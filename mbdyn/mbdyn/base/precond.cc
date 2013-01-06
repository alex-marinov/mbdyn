/* $Header$ */
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
  *
  * Copyright (C) 2003-2013
  * Giuseppe Quaranta	<quaranta@aero.polimi.it>
  *
  * classi che implementano la risoluzione del sistema nonlineare 
  */
  
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
  
#include "precond_.h"  

Preconditioner::~Preconditioner(void)
{
	NO_OP;
}
	
FullJacobianPr::~FullJacobianPr(void)
{
	NO_OP;
}


void
FullJacobianPr::Precond(VectorHandler& b, VectorHandler& x, 
		SolutionManager* pSM) const
{
	/*
	 * To be conservative, restore the existing vectors when done :)
	 */
	doublereal *pdr = 0, *pds = 0;

	/*
	 * FIXME: what if they're null, but in general need be different?
	 *
	 * better do a

	 bool bChangePointers(doublereal *newRhs, doublereal *newSol,
	 		doublereal *&oldRhs, doublereal *&oldSol)
	
	 * that, if both are needed does what expected, otherwise
	 * copies rhs in sol and uses only sol
	 */
	if (pSM->pSolHdl() != pSM->pResHdl()) {	
		pdr = pSM->pdSetResVec(b.pdGetVec());
		pds = pSM->pdSetSolVec(x.pdGetVec());

	} else {
		x = b;
		pdr = pSM->pdSetResVec(x.pdGetVec());
	}

	pSM->Solve();

	(void)pSM->pdSetResVec(pdr);
	if (pds != 0 && pds != pdr) {
		(void)pSM->pdSetSolVec(pds);
	}
}

