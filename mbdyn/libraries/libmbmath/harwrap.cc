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

#include <harwrap.h>

/* HSMatrixHandler - begin: code */

#ifdef DEBUG_MEMMANAGER
clMemMan MHmm("HSMatrixHandler");
#endif

HSMatrixHandler::HSMatrixHandler(integer iMSize, 
				 integer** ppiTmpRow,
				 integer** ppiTmpCol,
				 doublereal** ppdTmpMat,
				 integer iWSSize)
: iWorkSpaceSize(iWSSize),
iCurSize(iWSSize), 
iNumItem(0), iMatSize(iMSize),
pHS(NULL), dZero(0.),
ppiRow(ppiTmpRow),
ppiCol(ppiTmpCol), 
ppdMat(ppdTmpMat)
{
#ifdef DEBUG
   IsValid();
#endif

   SAFENEWWITHCONSTRUCTOR(pHS,
			  HarwellSparse,
			  HarwellSparse(iCurSize, ppiRow, ppiCol), 
			  MHmm);
}


HSMatrixHandler::~HSMatrixHandler(void) 
{
#ifdef DEBUG
   IsValid();
#endif
   
   if (pHS != NULL) {
      SAFEDELETE(pHS, MHmm);
   }
}

integer HSMatrixHandler::PacMat(void)
{
  //  ASSERT(pMH->iCurSize > 0);   
   
   doublereal* pdMatNew = *ppdMat;
   integer* piRowNew = *ppiRow;
   integer* piColNew = *ppiCol;
   
   doublereal* pdMatOld = *ppdMat;
   integer* piRowOld = *ppiRow;
   integer* piColOld = *ppiCol;
   
   integer iEmpty = -(this->pHS->iCurSize+1);
   integer iCount = 0;
   
   for (; piRowOld < *ppiRow+this->iCurSize; pdMatOld++, piRowOld++, piColOld++) {
      if (*(piRowOld) != iEmpty) {
	 union uPacVec uPV;
	 uPV.iInt = *piColOld;
	 *pdMatNew++ = *pdMatOld;
	 *piRowNew++ = integer(uPV.sRC.ir);
	 *piColNew++ = integer(uPV.sRC.ic);	     
	 iCount++;
      }
   }
   return iCount;
}

void HSMatrixHandler::IsValid(void) const
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
   ASSERT(SMmm.fIsValid((void*)*ppdMat, iWorkSpaceSize*sizeof(doublereal)));
#endif
   
}


void HSMatrixHandler::Init(const doublereal& dResetVal)
{
#ifdef DEBUG
   IsValid();
#endif
   
   ASSERT(pHS != NULL);
#ifdef DEBUG_MEMMANAGER
   ASSERT(MHmm.fIsPointerToBlock((void*)pHS));
#endif
   
   pHS->ResetVec();
   doublereal* pdTmp = *ppdMat;
   while ( pdTmp < *ppdMat+iCurSize ) {
      *pdTmp++ = dResetVal; 
   }  
}

/* HSMatrixHandler - end */






/* HarwellLUSolver - begin */

#ifdef DEBUG_MEMMANAGER
clMemMan LUmm(sLUClassName);
#endif

/* HarwellLUSolver - end */








/* HSLUSolutionManager - begin: code */

#ifdef DEBUG_MEMMANAGER
clMemMan SMmm("HSLUSolutionManager");
#endif

/* Costruttore */

HSLUSolutionManager::HSLUSolutionManager(integer iSize, 
					 integer iWorkSpaceSize,
					 const doublereal& dPivotFactor) :
iMatMaxSize(iSize),
iMatSize(iSize), 
piRow(NULL), piCol(NULL), 
pdMat(NULL), pdVec(NULL),
pMH(NULL), pVH(NULL), pLU(NULL),
fHasBeenReset(1)
{
   ASSERT(iSize > 0);
   ASSERT((dPivotFactor >= 0.0) && (dPivotFactor <= 1.0));

   /* Valore di default */
   if (iWorkSpaceSize == 0) {
      iWorkSpaceSize = iSize*iSize;
   }
      
   /* Alloca arrays */
   SAFENEWARR(piRow, integer, iWorkSpaceSize, SMmm);
   SAFENEWARR(piCol, integer, iWorkSpaceSize, SMmm);
   SAFENEWARR(pdMat, doublereal, iWorkSpaceSize, SMmm);
   SAFENEWARR(pdVec, doublereal, iMatSize, SMmm);
   
   /* Alloca handlers ecc. */
   SAFENEWWITHCONSTRUCTOR(pMH, 
			  HSMatrixHandler,
			  HSMatrixHandler(iMatSize, &piRow, &piCol, &pdMat, 
					  iWorkSpaceSize), 
			  SMmm);
   SAFENEWWITHCONSTRUCTOR(pVH,
			  MyVectorHandler,
			  MyVectorHandler(iMatSize, pdVec), SMmm);
   SAFENEWWITHCONSTRUCTOR(pLU, 
			  HarwellLUSolver,
			  HarwellLUSolver(iMatSize, iWorkSpaceSize, &piRow, &piCol, 
					  &pdMat, pdVec, dPivotFactor), SMmm);
   
#ifdef DEBUG
   IsValid();
#endif      
}


/* Distruttore; verificare la distruzione degli oggetti piu' complessi */

HSLUSolutionManager::~HSLUSolutionManager(void)
{
#ifdef DEBUG
   IsValid();
#endif      
   
   /* Dealloca oggetti strani */
   if (pLU != NULL) {	
      SAFEDELETE(pLU, SMmm);
   }   
   if (pVH != NULL) {      
      SAFEDELETE(pVH, SMmm);
   }   
   if (pMH != NULL) {
      SAFEDELETE(pMH, SMmm);
   }      
   
   /* Dealloca arrays */
   if (pdVec != NULL) {	
      SAFEDELETEARR(pdVec, SMmm);
   }   
   if (pdMat != NULL) {	
      SAFEDELETEARR(pdMat, SMmm);
   }   
   if (piCol != NULL) {	
      SAFEDELETEARR(piCol, SMmm);
   }   
   if (piRow != NULL) {	
      SAFEDELETEARR(piRow, SMmm);
   }   
}


/* Test di validita' del manager */

void HSLUSolutionManager::IsValid(void) const
{   
   ASSERT(iMatMaxSize > 0);
   ASSERT(iMatSize > 0);
   ASSERT(pMH != NULL);
   ASSERT(pdMat != NULL);
   ASSERT(piRow != NULL);
   ASSERT(piCol != NULL);
   
#ifdef DEBUG_MEMMANAGER
   ASSERT(SMmm.fIsPointerToBlock((void*)piRow));
   ASSERT(SMmm.fIsPointerToBlock((void*)piCol));
   ASSERT(SMmm.fIsPointerToBlock((void*)pdMat));
   ASSERT(SMmm.fIsPointerToBlock((void*)pMH));
   ASSERT(SMmm.fIsPointerToBlock((void*)pVH));
   ASSERT(SMmm.fIsPointerToBlock((void*)pLU));
#endif
   
   ASSERT((pMH->IsValid(), 1));
   ASSERT((pVH->IsValid(), 1));
   ASSERT((pLU->IsValid(), 1));
}

/* Prepara i vettori e la matrice per il solutore */

void HSLUSolutionManager::PacVec(void)
{
#ifdef DEBUG
   IsValid();
#endif
   
   /*
   ASSERT(pMH != NULL);
   ASSERT(pdMat != NULL);
   ASSERT(piRow != NULL);
   ASSERT(piCol != NULL);
#ifdef DEBUG_MEMMANAGER
   ASSERT(SMmm.fIsPointerToBlock((void*)pMH));
   ASSERT(SMmm.fIsPointerToBlock((void*)piRow));
   ASSERT(SMmm.fIsPointerToBlock((void*)piCol));
   ASSERT(SMmm.fIsPointerToBlock((void*)pdMat));
#endif
    */
   
   ASSERT(pMH->iCurSize > 0);   
   ASSERT(fHasBeenReset == 1);
   
   doublereal* pdMatNew = pdMat;
   integer* piRowNew = piRow;
   integer* piColNew = piCol;
   
   doublereal* pdMatOld = pdMat;
   integer* piRowOld = piRow;
   integer* piColOld = piCol;
   
   integer iEmpty = -(pMH->pHS->iCurSize+1);
   integer iCount = 0;
   
   for (; piRowOld < piRow+pMH->iCurSize; pdMatOld++, piRowOld++, piColOld++) {
      if (*(piRowOld) != iEmpty) {
	 union uPacVec uPV;
	 uPV.iInt = *piColOld;
	 *pdMatNew++ = *pdMatOld;
	 *piRowNew++ = integer(uPV.sRC.ir);
	 *piColNew++ = integer(uPV.sRC.ic);	     
	 iCount++;
      }
   }
   
   pLU->iNonZeroes = iCount;
}


/* Inizializza il gestore delle matrici */

void HSLUSolutionManager::MatrInit(const doublereal& dResetVal)
{
#ifdef DEBUG
   IsValid();
#endif
   
   /*
   ASSERT(pMH != NULL);
#ifdef DEBUG_MEMMANAGER
   ASSERT(SMmm.fIsPointerToBlock((void*)pMH));
#endif
    */
   
   pMH->Init(dResetVal);
   fHasBeenReset = flag(1);
}


/* Risolve il problema */

void HSLUSolutionManager::Solve(void)
{
#ifdef DEBUG
   IsValid();
#endif

   /*
   ASSERT(iMatSize > 0);
   ASSERT(piRow != NULL);
   ASSERT(piCol != NULL);
   ASSERT(pdMat != NULL);
   ASSERT(pdVec != NULL);
   ASSERT(pMH != NULL);
   ASSERT(pVH != NULL);
   ASSERT(pLU != NULL);
#ifdef DEBUG_MEMMANAGER
   ASSERT(SMmm.fIsPointerToBlock((void*)piRow));
   ASSERT(SMmm.fIsPointerToBlock((void*)piCol));
   ASSERT(SMmm.fIsPointerToBlock((void*)pdMat));
   ASSERT(SMmm.fIsPointerToBlock((void*)pdVec));
   ASSERT(SMmm.fIsPointerToBlock((void*)pMH));
   ASSERT(SMmm.fIsPointerToBlock((void*)pVH));
   ASSERT(SMmm.fIsPointerToBlock((void*)pLU));
#endif
    */
   
   if (fHasBeenReset == 1) {
      this->PacVec();
      fHasBeenReset = flag(0);
      flag fReturnFlag = pLU->fLUFactor();
      if (fReturnFlag < 0) {	 
	 THROW(HSLUSolutionManager::ErrGeneric());
      }
   }
   pLU->Solve();
}

/* HSLUSolutionManager - end */
