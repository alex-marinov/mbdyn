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
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
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

/*****************************************************************************
 *                                                                           *
 *                            HARWLIB C++ WRAPPER                            *
 *                                                                           *
 *****************************************************************************/

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <sparsemh.h>

/* SparseData - begin */

/* Inizializza i vettori con il valore -(iSize+1) ecc. */
void
SparseData::ResetVec(void) {
	ASSERT(iMaxSize > 0);
	ASSERT(iCurSize > 0);
	ASSERT(ppiTable != NULL);
	ASSERT(ppiKeys != NULL);
	ASSERT(*ppiTable != NULL);
	ASSERT(*ppiKeys != NULL);
#ifdef DEBUG_MEMMANAGER
	ASSERT(SMmm.fIsPointerToBlock((void*)*ppiTable));
	ASSERT(SMmm.fIsPointerToBlock((void*)*ppiKeys));
#endif /* DEBUG_MEMMANAGER */

	integer iSize = iCurSize;
	__FC_DECL__(kd01a)(&iSize, *ppiTable, *ppiKeys);
};

/*
 * Costruttore.
 * iSize -   dimensione degli spazi di lavoro
 * piTable - puntatore a punt. ad un array di interi di dim. iSize,
 *           che viene usato come tabella di stato del campo,
 * piKeys  - puntatore a punt. ad un array delle stesse dimensioni del
 *           precedente, che viene usato per contenere le keys
 */
SparseData::SparseData(integer iSize, 
		       integer** ppiTmpTable,
		       integer** ppiTmpKeys)
: iMaxSize(iSize),
iCurSize(iSize),
iNumItem(0),
ppiTable(ppiTmpTable),
ppiKeys(ppiTmpKeys) {
	ASSERT(iCurSize > 0);
	ASSERT(ppiTable != NULL);
	ASSERT(ppiKeys != NULL);
	ASSERT(*ppiTable != NULL);
	ASSERT(*ppiKeys != NULL);
	
#ifdef DEBUG_MEMMANAGER
	ASSERT(SMmm.fIsPointerToBlock((void*)*ppiTable));
	ASSERT(SMmm.fIsPointerToBlock((void*)*ppiKeys));
#endif
};

/* Distruttore (banalissimo: non c'e' nulla da distruggere) */
SparseData::~SparseData(void) {
	NO_OP;
};

/* SparseData - end */


/* SparseMatrixHandler - begin */

#ifdef DEBUG_MEMMANAGER
clMemMan MHmm("SparseMatrixHandler");
#endif

SparseMatrixHandler::SparseMatrixHandler(integer iMSize, 
				         integer** ppiTmpRow,
					 integer** ppiTmpCol,
					 doublereal** ppdTmpMat,
					 integer iWSSize)
: iWorkSpaceSize(iWSSize),
iCurSize(iWSSize), 
iNumItem(0),
iMatSize(iMSize),
pSD(NULL),
ppiRow(ppiTmpRow),
ppiCol(ppiTmpCol), 
ppdMat(ppdTmpMat)
{
#ifdef DEBUG
   	IsValid();
#endif

   	SAFENEWWITHCONSTRUCTOR(pSD,
			       SparseData,
			       SparseData(iCurSize, ppiRow, ppiCol), 
			       MHmm);
}


SparseMatrixHandler::~SparseMatrixHandler(void) 
{
#ifdef DEBUG
   	IsValid();
#endif
   
   	if (pSD != NULL) {
      		SAFEDELETE(pSD, MHmm);
   	}
}

integer 
SparseMatrixHandler::PacMat(void)
{
   	doublereal* pdMatNew = *ppdMat;
   	integer* piRowNew = *ppiRow;
   	integer* piColNew = *ppiCol;
   
   	doublereal* pdMatOld = *ppdMat;
   	integer* piRowOld = *ppiRow;
   	integer* piColOld = *ppiCol;
   
   	integer iEmpty = -(iCurSize+1);
   	integer iCount = 0;
   
   	for ( ;
	     piRowOld < *ppiRow+this->iCurSize; 
	     pdMatOld++, piRowOld++, piColOld++) {
      		if (*(piRowOld) != iEmpty) {
	 		SparseData::uPacVec uPV;
	 		uPV.iInt = *piColOld;
	 		*pdMatNew++ = *pdMatOld;
	 		*piRowNew++ = integer(uPV.sRC.ir);
	 		*piColNew++ = integer(uPV.sRC.ic);	     
	 		iCount++;
      		}
   	}
	
   	return iCount;
}

void 
SparseMatrixHandler::IsValid(void) const
{   
   	ASSERT(iMatSize > 0);
   	ASSERT(iWorkSpaceSize > 0);
   	ASSERT(ppiRow != NULL);
   	ASSERT(ppiCol != NULL);
   	ASSERT(ppdMat != NULL);
   	ASSERT(*ppiRow != NULL);
   	ASSERT(*ppiCol != NULL);
   	ASSERT(*ppdMat != NULL);
   
#ifdef DEBUG_MEMMANAGER
   	ASSERT(SMmm.fIsValid((void*)*ppiRow, iWorkSpaceSize*sizeof(integer)));
   	ASSERT(SMmm.fIsValid((void*)*ppiCol, iWorkSpaceSize*sizeof(integer)));
   	ASSERT(SMmm.fIsValid((void*)*ppdMat, 
			     iWorkSpaceSize*sizeof(doublereal)));
#endif
}

void
SparseMatrixHandler::Init(const doublereal& dResetVal)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   
   	ASSERT(pSD != NULL);
#ifdef DEBUG_MEMMANAGER
   	ASSERT(MHmm.fIsPointerToBlock((void*)pSD));
#endif /* DEBUG_MEMMANAGER */
   
   	pSD->ResetVec();
   	doublereal* pdTmp = *ppdMat;
   	while ( pdTmp < *ppdMat+iCurSize ) {
      		*pdTmp++ = dResetVal; 
   	}
}

/* Prepara i vettori e la matrice per il solutore */
integer
SparseMatrixHandler::iPacVec(void)
{
#ifdef DEBUG
	IsValid();
#endif  

	ASSERT(pSD->iCurSize > 0);
	
	doublereal* pdMatNew = *ppdMat;
	integer* piRowNew = *ppiRow;
	integer* piColNew = *ppiCol;
	
	doublereal* pdMatOld = *ppdMat;
	integer* piRowOld = *ppiRow;
	integer* piColOld = *ppiCol;
	
	integer iEmpty = -(iCurSize+1);
	integer iCount = 0;              
	
	for ( ;
	     piRowOld < *ppiRow+iCurSize;
	     pdMatOld++, piRowOld++, piColOld++) {	     
	     	if (*piRowOld != iEmpty) {
			SparseData::uPacVec uPV;
			uPV.iInt = *piColOld;		
			*pdMatNew++ = *pdMatOld;
			*piRowNew++ = integer(uPV.sRC.ir);
			*piColNew++ = integer(uPV.sRC.ic);
			iCount++;
		}
	}
	
	return iCount;
}

/* SparseMatrixHandler - end */

