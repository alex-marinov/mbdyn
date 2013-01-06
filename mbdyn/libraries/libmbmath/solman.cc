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

/* solution manager */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <string.h>	/* for memset() */

#include <iostream>
#include <iomanip>

#include "solman.h"
#include "matvec3.h"
#include "ls.h"

/* SolutionDataManager - begin */

SolutionDataManager::~SolutionDataManager(void)
{
	NO_OP;
}

/* SolutionDataManager - end */


/* SolutionManager - begin */

SolutionManager::SolutionManager(void)
: pLS(0)
{
	NO_OP;
}

SolutionManager::~SolutionManager(void)
{
   	if (pLS != NULL) {	
      		SAFEDELETE(pLS);
   	}
}

/* Inizializzatore "speciale" */
void
SolutionManager::MatrInitialize(void)
{
	MatrReset();
}

void
SolutionManager::SolveT(void)
{
	silent_cerr("SolutionManager::SolveT() not supported" << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

void
SolutionManager::LinkToSolution(VectorHandler& XCurr,
	VectorHandler& XPrimeCurr)
{
	NO_OP;
}

/* sposta il puntatore al vettore del residuo */
doublereal *
SolutionManager::pdSetResVec(doublereal* pd)
{
	ASSERT(pLS);
	return pLS->pdSetResVec(pd);
}
   
/* sposta il puntatore al vettore della soluzione */
doublereal *
SolutionManager::pdSetSolVec(doublereal* pd)
{
	ASSERT(pLS);
	return pLS->pdSetSolVec(pd);
}

bool SolutionManager::bGetConditionNumber(doublereal& dCond) const
{
	return pLS->bGetConditionNumber(dCond);
}

/* SolutionManager - end */

