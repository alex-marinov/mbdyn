/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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

/* dof manager */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "dataman.h"

/* DataManager - begin */

void DataManager::DofManager(void)
{
   DummyDofOwner.iFirstIndex = 0;
   DummyDofOwner.iNumDofs = 0;
   
   /* Resetta la struttura statica */
   for(int i = 0; i < DofOwner::LASTDOFTYPE; i++) {
      DofData[i].pFirstDofOwner = NULL;
      DofData[i].iNum = 0;
      DofData[i].iSize = 0;
      DofData[i].dDefScale = 1.;
   }   
}


void DataManager::DofManagerDestructor(void)
{
   DEBUGCOUTFNAME("DataManager::DofManagerDestructor");
   
   /* Distrugge le strutture dinamiche */
   ASSERT(pDofOwners != NULL);   
   if(pDofOwners != NULL) {	
      DEBUGLCOUT(MYDEBUG_INIT, "deleting dof owners structure" << std::endl);
      SAFEDELETEARR(pDofOwners);
   }
   
   if(pDofs != NULL) {	
      DEBUGLCOUT(MYDEBUG_INIT, "deleting dofs structure" << std::endl);
      SAFEDELETEARR(pDofs);
   }   
}

doublereal
DataManager::dGetDefaultScale(DofOwner::Type t) const
{
	return DofData[t].dDefScale;
}

void DataManager::DofDataInit(void)
{
   /* struttura dei DofOwner */  
   
   /* Calcola il numero totale di DofOwner */
   for (int iCnt = 0; iCnt < DofOwner::LASTDOFTYPE; iCnt++) {      
      iTotDofOwners += DofData[iCnt].iNum;
   }   

   DEBUGLCOUT(MYDEBUG_INIT, "iTotDofOwners = " << iTotDofOwners << std::endl);
	
   /* Crea la struttura dinamica dei DofOwner */
   if (iTotDofOwners > 0) {	     
      SAFENEWARRNOFILL(pDofOwners, DofOwner, iTotDofOwners);
      
      /* Resetta la struttura dinamica dei DofOwner */
      for (int iCnt = 0; iCnt < iTotDofOwners; iCnt++) {
	 pDofOwners[iCnt].iFirstIndex = 0;
	 pDofOwners[iCnt].iNumDofs = 0;
      }
      
      /* Inizializza la struttura dinamica dei DofOwner
       * con il numero di Dof di ognuno */
      DofData[0].pFirstDofOwner = pDofOwners;
      for (int iType = 0; iType < DofOwner::LASTDOFTYPE - 1; iType++) {
	 DofData[iType + 1].pFirstDofOwner =
	   DofData[iType].pFirstDofOwner+
	   DofData[iType].iNum;

	 for (int iDof = 0; iDof < DofData[iType].iNum; iDof++) {
	    DofData[iType].pFirstDofOwner[iDof].SetScale(dGetDefaultScale(DofOwner::Type(iType)));
	 }
      }

   } else {
      /* Se non sono definiti DofOwners, la simulazione non ha senso,
       * quindi il programma termina */
      silent_cerr("warning, no dof owners are defined" << std::endl);
      throw NoErr(MBDYN_EXCEPT_ARGS);
   }
}

void DataManager::DofInit(void)
{  
   if (iTotDofOwners > 0) {	
      
      /* Di ogni DofOwner setta il primo indice
       * e calcola il numero totale di Dof */
      DofOwner* pTmp = pDofOwners;  /* viene fatto scorrere 
						 * sulla struttura dei 
						 * DofOwners */
      integer iIndex = 0;    /* contatore dei Dof */
      integer iNumDofs = 0;  /* numero di dof di un owner */
      for(int iCnt = 1; 
	  pTmp < pDofOwners+iTotDofOwners; 
	  iCnt++, pTmp++) {
	 if((iNumDofs = pTmp->iNumDofs) > 0) {
	    pTmp->iFirstIndex = iIndex;
	    iIndex += iNumDofs;
	 } else {
	    pTmp->iFirstIndex = -1;
	    DEBUGCERR("warning, item " << iCnt << " has 0 dofs" << std::endl);
	 }
      }
      
      iTotDofs = iIndex;
	
      DEBUGLCOUT(MYDEBUG_INIT, "iTotDofs = " << iTotDofs << std::endl);
   } else {
      DEBUGCERR("");
      silent_cerr("no dof owners are defined" << std::endl);
      
      throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
   }	   	
   
   	
   /* Crea la struttura dinamica dei Dof */
   if(iTotDofs > 0) {	
      SAFENEWARRNOFILL(pDofs, Dof, iTotDofs);
      
      /* Inizializza l'iteratore sui Dof */
      DofIter.Init(pDofs, iTotDofs);
      
      /* Inizializza la struttura dinamica dei Dof */
      Dof* pTmp = pDofs;
      integer iIndex = pDofOwners[0].iFirstIndex;
      while(pTmp < pDofs+iTotDofs) {
	 pTmp->iIndex = iIndex++;
	 pTmp->Order = DofOrder::DIFFERENTIAL;
	 pTmp++;
      }	
   } else {
      DEBUGCERR("");
      silent_cerr("no dofs are defined" << std::endl);
      
      throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
   }	   	
}  

void
DataManager::SetScale(VectorHandler& XScale) const
{
	for (integer iCnt = 0; iCnt < iTotDofOwners; iCnt++) {
		integer iFirstIndex = pDofOwners[iCnt].iFirstIndex;
		unsigned int iNumDofs = pDofOwners[iCnt].iNumDofs;
		doublereal dScale = pDofOwners[iCnt].dScale;

		for (unsigned int iDof = 1; iDof <= iNumDofs; iDof++) {
			XScale.PutCoef(iFirstIndex + iDof, dScale);
		}
	}
}

/* DataManager - end */
