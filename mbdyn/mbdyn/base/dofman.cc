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

/* dof manager */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <dataman.h>

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
   }   
}


void DataManager::DofManagerDestructor(void)
{
   DEBUGCOUTFNAME("DataManager::DofManagerDestructor");
   
   /* Distrugge le strutture dinamiche */
   ASSERT(pDofOwners != NULL);   
   if(pDofOwners != NULL) {	
      DEBUGLCOUT(MYDEBUG_INIT, "deleting dof owners structure" << endl);
      SAFEDELETEARR(pDofOwners, DMmm);
   }
   
   // ASSERT(pDofs != NULL);  // it could haven't been allocated
   if(pDofs != NULL) {	
      DEBUGLCOUT(MYDEBUG_INIT, "deleting dofs structure" << endl);
      SAFEDELETEARR(pDofs, DMmm);
   }   
}


void DataManager::DofDataInit(void)
{
   /* struttura dei DofOwner */  
   
   /* Calcola il numero totale di DofOwner */
   for (int iCnt = 0; iCnt < DofOwner::LASTDOFTYPE; iCnt++) {      
      iTotDofOwners += DofData[iCnt].iNum;
   }   

   DEBUGLCOUT(MYDEBUG_INIT, "iTotDofOwners = " << iTotDofOwners << endl);
	
   /* Crea la struttura dinamica dei DofOwner */
   if (iTotDofOwners > 0) {	     
      SAFENEWARR(pDofOwners, DofOwner, iTotDofOwners, DMmm);
      
      /* Resetta la struttura dinamica dei DofOwner */
      DofOwner* pTmp = pDofOwners;
      while (pTmp < pDofOwners+iTotDofOwners) {
	 pTmp->iFirstIndex = 0;
	 pTmp->iNumDofs = 0;
	 pTmp++;
      }
      
      /* Inizializza la struttura dinamica dei DofOwner
       * con il numero di Dof di ognuno */
      DofData[0].pFirstDofOwner = pDofOwners;
      for (int iCnt = 1; iCnt < DofOwner::LASTDOFTYPE; iCnt++) {
	 DofData[iCnt].pFirstDofOwner =
	   DofData[iCnt-1].pFirstDofOwner+
	   DofData[iCnt-1].iNum;
      }
   } else {
      /* Se non sono definiti DofOwners, la simulazione non ha senso,
       * quindi il programma termina */
      DEBUGCERR("");
      cerr << "warning, no dof owners are defined; program is aborting ..." 
	<< endl;
#ifdef USE_EXCEPTIONS      
      throw NoErr();
#else
      exit(EXIT_SUCCESS);
#endif      
   }	   
}

void DataManager::DofInit(void)
{  
   if( iTotDofOwners > 0) {	
      
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
	    DEBUGCERR("warning, item " << iCnt << " has 0 dofs" << endl);
	 }
      }
      
      iTotDofs = iIndex;
	
      DEBUGLCOUT(MYDEBUG_INIT, "iTotDofs = " << iTotDofs << endl);
   } else {
      DEBUGCERR("");
      cerr << "no dof owners are defined; aborting ..." << endl;
      
      THROW(DataManager::ErrGeneric());
   }	   	
   
   	
   /* Crea la struttura dinamica dei Dof */
   if(iTotDofs > 0) {	
      SAFENEWARR(pDofs, Dof, iTotDofs, DMmm);
      
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
      cerr << "no dofs are defined, aborting ..." << endl;
      
      THROW(DataManager::ErrGeneric());
   }	   	
}  

/* DataManager - end */
