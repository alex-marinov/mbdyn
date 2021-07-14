/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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
 * Copyright (C) 1999-2017
 * Giuseppe Quaranta     <quaranta@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 */

/* Schur Matrix Handler */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "schurmh.h"

/* dimensioni */
integer
SchurMatrixHandler::iGetNumRows(void) const
{
	return LSize + ISize;
}

integer
SchurMatrixHandler::iGetNumCols(void) const
{
	return LSize + ISize;
}

/* Restituisce l'handler alla matrice B */
MatrixHandler*
SchurMatrixHandler::GetBMat(void)
{
	return pB;
}

/* Setta l'handler alla matrice B */
void
SchurMatrixHandler::SetBMat(MatrixHandler* pBM)
{
	pB = pBM;
}

/* Restituisce il puntatore alla matrice C
 * che e' un vettore doublereal contenente le righe in successione */
doublereal*
SchurMatrixHandler::GetCMat(void)
{
	return pdC;
}

SchurMatrixHandler::SchurMatrixHandler(int LocSize, int IntSize,
		MatrixHandler* pBM, integer* pGlobToLoc, doublereal* pdEv)
: LSize(LocSize),
ISize(IntSize),
pB(pBM),
pE(0),
pdE(0),
pF(0),
pC(0),
pdC(0),
pGTL(pGlobToLoc),
bExtpdE(pdEv != 0 ? true : false)
{
	if (!bExtpdE) {
		SAFENEWARR(pdE, doublereal, LSize*ISize);
		pdEv = pdE;
	}

	SAFENEWWITHCONSTRUCTOR(pE,
				MyVectorHandler,
				MyVectorHandler(LSize*ISize, pdEv));
	SAFENEWARR(pdC, doublereal, ISize*ISize);
	SAFENEWWITHCONSTRUCTOR(pC,
				MyVectorHandler,
				MyVectorHandler(ISize*ISize, pdC));
	SAFENEWWITHCONSTRUCTOR(pF,
				SpMapMatrixHandler,
				SpMapMatrixHandler(ISize, LSize));

#ifdef DEBUG
	IsValid();
#endif /* DEBUG */
}

SchurMatrixHandler::~SchurMatrixHandler(void)
{
	if (pE != 0) {
		SAFEDELETE(pE);
	}
	if (pC != 0) {
		SAFEDELETE(pC);
	}
	if (pF != 0) {
		SAFEDELETE(pF);
	}
	if (!bExtpdE && pdE != 0) {
		SAFEDELETEARR(pdE);
	}
	if (pdC != 0) {
		SAFEDELETEARR(pdC);
	}
}

#ifdef DEBUG
/* Usata per il debug */
void
SchurMatrixHandler::IsValid(void) const
{
	ASSERT(LSize >0);
	ASSERT(ISize >0);
	ASSERT(pB != 0);
	ASSERT(pE != 0);
	ASSERT(pF != 0);
	ASSERT(pC != 0);
	ASSERT(pGTL != 0);
}
#endif /* DEBUG */

/* Ridimensiona la matrice */
void
SchurMatrixHandler::Resize(integer, integer)
{
	silent_cerr("cannot resize a SchurMatrixHandler" << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

/* SchurMatrixHandler - End */

/* SchurVectorHandler - Start */

SchurVectorHandler::SchurVectorHandler(int LocSize, int IntSize,
		VectorHandler* pLocVec, integer* pGlobToLoc)
: LSize(LocSize),
ISize(IntSize),
pLV(pLocVec),
pIV(0),
bExtpIV(false),
pGTL(pGlobToLoc)
{
	SAFENEWWITHCONSTRUCTOR(pIV,
			MyVectorHandler,
			MyVectorHandler(IntSize));
}

SchurVectorHandler::SchurVectorHandler(int LocSize, int IntSize,
		VectorHandler* pLocV, VectorHandler* pIntV,
		integer* pGlobToLoc)
: LSize(LocSize),
ISize(IntSize),
pLV(pLocV),
pIV(pIntV),
bExtpIV(pIntV != 0 ? true : false),
pGTL(pGlobToLoc)
{
	NO_OP;
}


SchurVectorHandler::~SchurVectorHandler(void)
{
	if (!bExtpIV && pIV != 0) {
		SAFEDELETE(pIV);
	}
}

#ifdef DEBUG
/* Usata per il debug */
void
SchurVectorHandler::IsValid(void) const
{
	NO_OP;
}
#endif /* DEBUG */

/* SchurVectorHandler - End */


/* SchurMatrixHandlerUm - Begin */

SchurMatrixHandlerUm::SchurMatrixHandlerUm(int LocSize, int IntSize,
		MatrixHandler* pBM,
		integer* pGlobToLoc)
: SchurMatrixHandler(LocSize, IntSize, pBM, pGlobToLoc, 0),
pdEs(0),
pEs(0),
Eflag(1)
{
	SAFENEWARR(pdEs, doublereal, LSize*(ISize + 1));
	pdE = pdEs + LSize;
	pE->Attach(LSize*ISize, pdE, LSize*ISize);
	SAFENEWWITHCONSTRUCTOR(pEs,
			MyVectorHandler,
			MyVectorHandler(LSize*ISize, pdEs));
}

SchurMatrixHandlerUm::~SchurMatrixHandlerUm(void)
{
	if (pE != 0) {
		SAFEDELETE(pE);
	}

	if (pdEs != 0) {
		SAFEDELETEARR(pdEs);
	}

	pdE = 0;
}

/* SchurMatrixHandlerUm - End */

