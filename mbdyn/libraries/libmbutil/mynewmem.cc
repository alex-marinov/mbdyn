/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 * Paolo Mantegazza     <mantegazza@aero.polimi.it>
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

	Allocazione dinamica di memoria controllata con ASSERT e Memory Manager

	Scritto da
	Pierangelo Masarati
	il 05/04/1997

 *****************************************************************************/

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#ifdef DEBUG

#include <string.h>
#include <iostream>
#include <iomanip>

#include "myassert.h"
#include "mynewmem.h"


/* Funzioni usate anche senza memory manager */
void 
_Safenew(const char *file, int line, int flag)
{
   	std::cout.flush();
   	if (flag == 0) {
      		std::cerr << std::endl 
			<< "SAFENEW fault: NULL return pointer in file " 
			<< file << " at line " << line << std::endl;
   	} else if (flag == 1) {
      		std::cerr << std::endl 
			<< "SAFENEWARR fault: NULL return pointer in file " 
			<< file << " at line " << line << std::endl;
   	}
}


void 
_Safenewfill(void *pv, size_t size, char fill)
{
   	ASSERT(pv);
   	ASSERT(size);
   
   	char* pb = (char*)pv;
   	while (size--) {
      		*pb++ = fill;
   	}
}

#ifdef DEBUG_MEMMANAGER

clMemMan defaultMemoryManager("Default");

/* Funzioni proprie private del memory manager */
clMemMan::stList *
clMemMan::pstFindElem(const void* pvToFind) const
{
   	ASSERT(pvToFind);
   
   	stList *pstL = pstRoot;
   	ASSERT(pstL);
   
   	while (pstL->pstNext && (pstL->pstNext->stMB.pv <= pvToFind)) {
      		pstL = pstL->pstNext;
      		if (pstL->stMB.pv == pvToFind) {
	 		if (pstL->stMB.eSt == ALLOCATED) { 
	    			return pstL; 
	 		}
      		}
   	}
   	CERR << std::endl << "clMemMan " << sName << " error: pointer " 
     		<< (void*)pvToFind << " not found in pstFindElem()" << std::endl;
   	return NULL;
}

clMemMan::stList *
clMemMan::pstFindPrev(const void *pvToFindPrev) const
{
   	ASSERT(pvToFindPrev);
   
   	stList *pstL = pstRoot;
   	stList *pstN = NULL;
   	ASSERT(pstL);
   
   	while (pstL->pstNext && (pstL->pstNext->stMB.pv <= pvToFindPrev)) {
      		pstN = pstL->pstNext;
      		if (pstN->stMB.pv == pvToFindPrev) {
	 		if (pstN->stMB.eSt == ALLOCATED) {
	    			return pstL; 
	 		}
      		}
      		pstL = pstN;
   	}
   
   	CERR << std::endl << "clMemMan " << sName << " error: pointer " 
     		<< (void*)pvToFindPrev << " not found in pstFindPrev()" 
		<< std::endl;
   	return NULL;
}

/* enum eRemoveMode { RELEASE, DELBUTKEEPLOG, DELBUTNOTRELEASE }; */
void 
clMemMan::_remove(const void *pvToRemove, clMemMan::eRemoveMode eMode, flag fArr, flag fFill)
{
   	ASSERT(pvToRemove);
   
   	stList *pstL = pstFindPrev(pvToRemove);
   	ASSERT(pstL);
   
   	if (!pstL) {
      		CERR << std::endl << "clMemMan " << sName << " warning: pointer " 
			<< (void*)pvToRemove;
      		if (fArr) {
	 		std::cerr << " to array";
      		}
      		std::cerr << " not found in _remove()" << std::endl;
      		throw clMemMan::ErrNotFound();
   	}
   
   	stList *pstN = pstL->pstNext;
   	ASSERT(pstN);
   
   	if (fFill) {
      		ASSERT(fArr && pstN->stMB.fArr);
      		_Safenewfill(pstN->stMB.pv, pstN->stMB.size, cDebugFree);
   	}
   
   	switch (eMode) {
    	case RELEASE:
	   	/* caso di cancellazione totale */
	   	pstL->pstNext = pstN->pstNext;
	   	delete pstN;
	   	break;

    	case DELBUTKEEPLOG:
	   	/* Cancellazione della memoria mantenendo il registro */
	   	pstN->stMB.eSt = FREED;
	   	break;

    	case DELBUTNOTRELEASE:
	   	/* Eliminazione del riferimento senza cancellazione memoria */
	   	pstN->stMB.eSt = FREEDBUTNOTRELEASED;
	   	break;
	}
}

/* Funzioni proprie pubbliche del memory manager */
clMemMan::clMemMan(char *sNameIn)
: pstRoot(NULL), sName(NULL)
{
   	if (sNameIn) {
      		sName = new char[strlen(sNameIn)+1];
      		strcpy(sName, sNameIn);
   	}
   
   	pstRoot = new stList(stMemBlock());
}

clMemMan::~clMemMan(void)
{
   	stList* pstL = pstRoot;
   	stList* pstP = NULL;
   	ASSERT(pstL);
   
   	while (pstL) {
      		pstP = pstL;
      		pstL = pstL->pstNext;
      		delete pstP;
   	}
   
   	if (sName) { 
      		delete[] sName; 
   	}
}

flag 
clMemMan::fIsBlock(const void *pvBlock, size_t sizeBlock) const
{
   	ASSERT(pvBlock);
   	ASSERT(sizeBlock);
   
   	stList *pstL = pstFindElem(pvBlock);
   	ASSERT(pstL);
   
   	if (pstL && (pstL->stMB.size == sizeBlock)) {
      		ASSERT(pstL->stMB.eSt != UNKNOWN);
      		if (pstL->stMB.eSt == ALLOCATED) { 
	 		return 1; 
      		}
   	}
	
   	return 0;
}

flag 
clMemMan::fIsPointerToBlock(const void *pvBlock) const
{
   	ASSERT(pvBlock);
   
   	stList *pstL = pstFindElem(pvBlock);
   	ASSERT(pstL);

   	if (pstL) {
      		ASSERT(pstL->stMB.eSt != UNKNOWN);
      		if (pstL->stMB.eSt == ALLOCATED) { 
	 		return 1; 
      		}
   	}
	
   	return 0;
}

flag 
clMemMan::fIsValid(const void *pvValid, size_t sizeValid) const
{
   	ASSERT(pvValid);
   	ASSERT(sizeValid);
   
   	stList *pstL = pstRoot;
   	ASSERT(pstL);

   	while (pstL->pstNext) {
      		pstL = pstL->pstNext;
      		flag fCond1 = (pstL->stMB.pv <= pvValid);
      		flag fCond2 = (((void*)pstL->stMB.pv+pstL->stMB.size) 
			>= ((void*)pvValid+sizeValid));
      
      		if (fCond1 && fCond2) {
	 		ASSERT(pstL->stMB.eSt != UNKNOWN);
	 		if (pstL->stMB.eSt == ALLOCATED) { 
	    			return 1; 
	 		}
      		}
   	}
	
   	return 0;
}

size_t 
clMemMan::sizeOfBlock(const void* pvSizeOf) const
{
   	ASSERT(pvSizeOf);
   	stList *pstL = pstFindElem(pvSizeOf);
   	ASSERT(pstL);
   
   	if (pstL) {
      		return pstL->stMB.size; 
   	}
	
   	return 0;
}


flag 
clMemMan::fIsArray(const void *pvIsArray) const
{
   	ASSERT(pvIsArray);
   	stList *pstL = pstFindElem(pvIsArray);
   	ASSERT(pstL);
   
   	if (pstL) { 
      		return pstL->stMB.fArr; 
   	}
   	return 0;
}

eStatus 
clMemMan::eBlockStatus(const void *pvBStatus) const
{
   	ASSERT(pvBStatus);
   	stList *pstL = pstFindElem(pvBStatus);
   	ASSERT(pstL);
   
   	if (pstL) { 
      		return pstL->stMB.eSt; 
   	}
   	return UNKNOWN;
}

void 
clMemMan::ClearRefs(void)
{
   	stList *pstL = pstRoot;
   	ASSERT(pstL);
   
   	while(pstL->pstNext) {
      		pstL = pstL->pstNext;
      		pstL->stMB.fRef = 0;
   	}
}

void 
clMemMan::PutRef(const void *pvRef)
{
   	ASSERT(pvRef);
   
   	stList *pstL = pstFindElem(pvRef);
   	ASSERT(pstL);
   
   	if (pstL) { 
      		pstL->stMB.fRef = 1; 
   	}
}

flag 
clMemMan::fIsRefd(const void *pvIsRefd) const
{
   	ASSERT(pvIsRefd);
   
   	stList *pstL = pstFindElem(pvIsRefd);
   	ASSERT(pstL);
   
   	if (pstL && (pstL->stMB.eSt == ALLOCATED) && (pstL->stMB.fRef == 1)) {
      		return 1; 
   	}
	
   	return 0;
}

/* enum eStatus { UNKNOWN, ALLOCATED, FREED, FREEDBUTNOTRELEASED }; */
std::ostream& 
clMemMan::DumpRef(std::ostream& rout) const
{
   	rout << "Memory Manager 1.0";
   	if (sName) { 
		rout << ": " << sName; 
	}
   	rout << std::endl;
   
   	rout << "Ref'd blocks:" << std::endl;
   
   	stList *pstL = pstRoot;
   	ASSERT(pstL);
   
   	int iCount = 0;
   	while (pstL->pstNext) {
      		pstL = pstL->pstNext;
      		ASSERT(pstL->stMB.eSt != UNKNOWN);
      
      		iCount++;
      		if (pstL->stMB.fRef) {
	 		rout << "Block " << setw(4) << iCount << ':' << std::endl
	   			<< "pointer " << pstL->stMB.pv
	   			<< ", size " << pstL->stMB.size 
	   			<< ", status: ";
	 
	 		if (pstL->stMB.eSt == ALLOCATED) {
	    			rout << "ALLOCATED" << std::endl;
	 		} else if (pstL->stMB.eSt == FREED) {
	    			rout << "FREED" << std::endl;
	 		} else if (pstL->stMB.eSt == FREEDBUTNOTRELEASED) {
	    			rout << "FREEDBUTNOTRELEASED" << std::endl;
	 		}
      		}
   	}
   
   	rout << "Unref'd blocks:" << std::endl;
   
   	pstL = pstRoot;
   	iCount = 0;
   	while (pstL->pstNext) {
      		pstL = pstL->pstNext;
      		ASSERT(pstL->stMB.eSt != UNKNOWN);
      
      		iCount++;
      		if (!pstL->stMB.fRef) {
	 		rout << "Block " << setw(4) << iCount << ':' << std::endl
	   			<< "pointer " << pstL->stMB.pv
	   			<< ", size " << pstL->stMB.size
	   			<< ", status: ";
	 
	 		if (pstL->stMB.eSt == ALLOCATED) {
	    			rout << "ALLOCATED" << std::endl;
	 		} else if (pstL->stMB.eSt == FREED) {
	    			rout << "FREED" << std::endl;
	 		} else if (pstL->stMB.eSt == FREEDBUTNOTRELEASED) {
	    			rout << "FREEDBUTNOTRELEASED" << std::endl;
	 		}
      		}
   	}
   
   	return rout;
}

void 
clMemMan::add(const void *pvIn, size_t sizeIn, flag fArr)
{
   	ASSERT(pvIn);
   	ASSERT(sizeIn);
   	ASSERT(fArr == 0 || fArr == 1);
   
   	stList *pstL = pstRoot;
   	stList *pstN = NULL;
   	ASSERT(pstL);
   
   	while (pstL->pstNext) {
      		pstN = pstL->pstNext;
      		if ((pstN->stMB.pv == pvIn) && (pstN->stMB.eSt == ALLOCATED)) {
	 		CERR << std::endl << "clMemMan" << sName 
	   			<< " error: block pointed by " 
	   			<< (void*)pvIn << ", size " << sizeIn << std::endl
	   			<< "is already defined. Previous size is " 
	   			<< pstN->stMB.size << std::endl;
	 		return;
      		}
		
      		if (pstN->stMB.pv >= pvIn) {
	 		break; 
      		}
      		pstL = pstN;
   	}
   
   	pstN = new stList(stMemBlock((void*)pvIn, sizeIn, ALLOCATED, fArr));
   	ASSERT(pstN);
   
   	if (pstN == NULL) {
      		CERR << std::endl << "clMemMan " << sName 
			<< ": error in allocation in add()" << std::endl;
      		throw ErrMemory(MBDYN_EXCEPT_ARGS);
   	}
   
   	pstN->pstNext = pstL->pstNext;
   	pstL->pstNext = pstN;
}

/* Operatore friend del memory manager */
std::ostream& 
operator << (std::ostream& rout, const clMemMan& rm)
{
   	rout << "Memory Manager 1.0";
   	if (rm.sName) { 
      		rout << ": " << rm.sName; 
   	}
   	rout << std::endl;
   
   	clMemMan::stList *pstL = rm.pstRoot;
   
   	int iCount = 0;
   	while (pstL->pstNext) {
      		pstL = pstL->pstNext;
      		ASSERT(pstL->stMB.eSt != UNKNOWN);
	
      		rout << "Block " << setw(4) << (++iCount) << ':' << std::endl
			<< "pointer " << pstL->stMB.pv
			<< ", size " << pstL->stMB.size
			<< ", status: ";
      
      		if (pstL->stMB.eSt == ALLOCATED) {
	 		rout << "ALLOCATED" << std::endl;
      		} else if (pstL->stMB.eSt == FREED) {
	 		rout << "FREED" << std::endl;
      		} else if (pstL->stMB.eSt == FREEDBUTNOTRELEASED) {
	 		rout << "FREEDBUTNOTRELEASED" << std::endl;
      		}
   	}
   
   	return rout;
}

#endif /* DEBUG_MEMMANAGER */

#endif /* DEBUG */

