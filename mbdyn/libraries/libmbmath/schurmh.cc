/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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
 * Copyright (C)1996-2001 
 * Giuseppe Quaranta     <quaranta@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 */

/* Schur Matrix Handler */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <schurmh.h>

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

/* Restituisce il puntatore alla alla matrice C 
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
pE(NULL),
pdE(NULL),
pF(NULL),
pC(NULL),
pdC(NULL), 
pGTL(pGlobToLoc),
extpdE(pdEv == NULL ? true : false)
{
	if (pdEv == NULL) {
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
	if (pE != NULL) {
		SAFEDELETE(pE);
	}
	if (pC != NULL) {	
		SAFEDELETE(pC);
	}
	if (pF != NULL) {
		SAFEDELETE(pF);
	}
	if (extpdE) {
		if (pdE != NULL) {
			SAFEDELETEARR(pdE);
		}
	}
	if (pdC != NULL) {
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
	ASSERT(pB != NULL);
	ASSERT(pE != NULL);
	ASSERT(pF != NULL);
	ASSERT(pC != NULL);
	ASSERT(pGTL != NULL);
}
#endif /* DEBUG */

/* SchurMatrixHandler - End */

/* SchurVectorHandler - Start */

SchurVectorHandler::SchurVectorHandler(int LocSize, int IntSize,
		VectorHandler* pLocVec, integer* pGlobToLoc)
: LSize(LocSize),
ISize(IntSize),
pLV(pLocVec),
pIV(NULL),
pIntVec(NULL),
pGTL(pGlobToLoc) 
{
	SAFENEWARR(pIntVec, doublereal, IntSize); 
	SAFENEWWITHCONSTRUCTOR(pIV,
			MyVectorHandler,
			MyVectorHandler(IntSize, pIntVec));
}

SchurVectorHandler::SchurVectorHandler(int LocSize, int IntSize,
		VectorHandler* pLocV, VectorHandler* pIntV,
		integer* pGlobToLoc)
: LSize(LocSize),
ISize(IntSize),
pLV(pLocV),
pIV(pIntV),
pIntVec(NULL),
pGTL(pGlobToLoc) 
{
	NO_OP;
}
     

SchurVectorHandler::~SchurVectorHandler(void)
{
	if (pIntVec != NULL) {
      		SAFEDELETEARR(pIntVec);
	}
	if (pIV != NULL) {
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
: SchurMatrixHandler(LocSize, IntSize, pBM, pGlobToLoc, NULL), 
pdEs(NULL),
pEs(NULL),
Eflag(1)
{ 
	SAFENEWARR(pdEs, doublereal, LSize*(ISize+1));
	pdE = pdEs+LSize;
	pE->Attach(LSize*ISize,pdE,LSize*ISize);
	SAFENEWWITHCONSTRUCTOR(pEs,
			MyVectorHandler,
			MyVectorHandler((LSize*ISize), pdEs));
}

SchurMatrixHandlerUm::~SchurMatrixHandlerUm(void)
{
	if (pE != NULL) {
		SAFEDELETE(pE);
	}
	if (pdEs != NULL) {
		SAFEDELETEARR(pdEs);
	};
}

/* SchurMatrixHandlerUm - End*/

