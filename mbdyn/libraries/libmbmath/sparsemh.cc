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
	ASSERT(defaultMemoryManager.fIsPointerToBlock(*ppiTable));
	ASSERT(defaultMemoryManager.fIsPointerToBlock(*ppiKeys));
#endif /* DEBUG_MEMMANAGER */

	/*
	 * Note: sets the arrays to IEMPTY = - (iCurSize + 1)
	 * so size MUST NOT change unless right before calling 
	 * this function
	 */

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
	ASSERT(defaultMemoryManager.fIsPointerToBlock(*ppiTable));
	ASSERT(defaultMemoryManager.fIsPointerToBlock(*ppiKeys));
#endif /* DEBUG_MEMMANAGER */
};

/* Distruttore (banalissimo: non c'e' nulla da distruggere) */
SparseData::~SparseData(void) {
	NO_OP;
};

bool 
SparseData::SetCurSize(integer i)
{
	if ( i < 1 || i > iMaxSize ) {
		return false;
	}

	iCurSize = i;

	return true;
}

/* SparseData - end */


/* SparseMatrixHandler - begin */

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
ppdMat(ppdTmpMat),
pdMatm1(NULL)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */

	pdMatm1 = *ppdMat - 1;

   	SAFENEWWITHCONSTRUCTOR(pSD, SparseData,
			       SparseData(iCurSize, ppiRow, ppiCol));
}


SparseMatrixHandler::~SparseMatrixHandler(void) 
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   
   	if (pSD != NULL) {
      		SAFEDELETE(pSD);
   	}
}

integer 
SparseMatrixHandler::PacMat(void)
{
	ASSERT(iCurSize > 0);

#if 0
   	doublereal* pdMatNew = *ppdMat;
   	integer* piRowNew = *ppiRow;
   	integer* piColNew = *ppiCol;
   
   	doublereal* pdMatOld = *ppdMat;
   	integer* piRowOld = *ppiRow;
   	integer* piColOld = *ppiCol;
#endif
   
   	integer iEmpty = -(iCurSize+1);
	integer iNew = 0;

#if 0
   	for ( ;
	     piRowOld < *ppiRow+this->iCurSize; 
	     pdMatOld++, piRowOld++, piColOld++) {
      		if (*(piRowOld) != iEmpty) {
	 		SparseData::uPacVec uPV;
	 		uPV.iInt = *piColOld;
	 		*pdMatNew++ = *pdMatOld;
	 		*piRowNew++ = integer(uPV.sRC.ir);
	 		*piColNew++ = integer(uPV.sRC.ic);	     
	 		iNew++;
      		}
   	}
#endif
	
	for (integer iOld = 0; iOld < iCurSize; iOld++) {
		if ((*ppiRow)[iOld] != iEmpty) {
			SparseData::uPacVec uPV;
			uPV.iInt = (*ppiCol)[iOld];
			(*ppdMat)[iNew] = (*ppdMat)[iOld];
			(*ppiRow)[iNew] = integer(uPV.sRC.ir);
			(*ppiCol)[iNew] = integer(uPV.sRC.ic);
			iNew++;
		}
	}


   	return iNew;
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
   	ASSERT(pdMatm1 != NULL);
   
#ifdef DEBUG_MEMMANAGER
   	ASSERT(defaultMemoryManager.fIsValid(*ppiRow, 
				iWorkSpaceSize*sizeof(integer)));
   	ASSERT(defaultMemoryManager.fIsValid(*ppiCol, 
				iWorkSpaceSize*sizeof(integer)));
   	ASSERT(defaultMemoryManager.fIsValid(*ppdMat, 
				iWorkSpaceSize*sizeof(doublereal)));
#endif /* DEBUG_MEMMANAGER */
}

void
SparseMatrixHandler::Init(const doublereal& dResetVal)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   
   	ASSERT(pSD != NULL);
#ifdef DEBUG_MEMMANAGER
   	ASSERT(defaultMemoryManager.fIsPointerToBlock(pSD));
#endif /* DEBUG_MEMMANAGER */
   
   	pSD->ResetVec();
#ifndef SPARSE_MATRIX_NO_RESET
#ifdef HAVE_MEMSET
	if (dResetVal == 0.) {
		/*
		 * Should be optimized on many architectures ...
		 * 
		 * 0.5 % gain on PIII!!!
		 */
		memset(*ppdMat, 0, sizeof(doublereal)*iCurSize);

	} else {
#endif /* HAVE_MEMSET */
		for ( integer cnt = 0; cnt < iCurSize; cnt++ ) {
			(*ppdMat)[cnt] = dResetVal;
		}
#ifdef HAVE_MEMSET
	}
#endif /* HAVE_MEMSET */
#endif /* !SPARSE_MATRIX_NO_RESET */
}

integer
SparseMatrixHandler::iGetCurSize(void) const
{
	return iCurSize;
}

bool
SparseMatrixHandler::SetCurSize(integer i)
{
	if (i < 1 || i > iWorkSpaceSize) {
		return false;
	}

	if (!pSD->SetCurSize(i)) {
		return false;
	}

	iCurSize = i;

	return true;
}

/* Prepara i vettori e la matrice per il solutore */
integer
SparseMatrixHandler::iPacVec(void)
{
#ifdef DEBUG
	IsValid();
#endif  /* DEBUG */

	ASSERT(iCurSize > 0);

#if 0	
	doublereal* pdMatNew = *ppdMat;
	integer* piRowNew = *ppiRow;
	integer* piColNew = *ppiCol;
	
	doublereal* pdMatOld = *ppdMat;
	integer* piRowOld = *ppiRow;
	integer* piColOld = *ppiCol;
#endif
	
	integer iEmpty = -(iCurSize+1);
	integer iNew = 0;

#if 0
	for ( ;
	     piRowOld < *ppiRow+iCurSize;
	     pdMatOld++, piRowOld++, piColOld++) {	     
	     	if (*piRowOld != iEmpty) {
			SparseData::uPacVec uPV;
			uPV.iInt = *piColOld;		
			*pdMatNew++ = *pdMatOld;
			*piRowNew++ = integer(uPV.sRC.ir);
			*piColNew++ = integer(uPV.sRC.ic);
			iNew++;
		}
	}
#endif
	
	for (integer iOld = 0; iOld < iCurSize; iOld++) {
		if ((*ppiRow)[iOld] != iEmpty) {
			SparseData::uPacVec uPV;
			uPV.iInt = (*ppiCol)[iOld];
			(*ppdMat)[iNew] = (*ppdMat)[iOld];
			(*ppiRow)[iNew] = integer(uPV.sRC.ir);
			(*ppiCol)[iNew] = integer(uPV.sRC.ic);
			iNew++;
		}
	}

	return iNew;
}

/* SparseMatrixHandler - end */

