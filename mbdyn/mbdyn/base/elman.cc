/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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

/* ElementManager */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <dataman.h>
#include <search.h>

/* DataManager - begin */

const char sEMClassName[] = "DataManager";

/* costruttore: resetta i dati */
void DataManager::ElemManager(void)
{
   /* Reset della struttura ElemData */
   for(int i = 0; i < Elem::LASTELEMTYPE; i++) {
      ElemData[i].ppFirstElem = NULL;
      ElemData[i].iNum = 0;
      ElemData[i].DofOwnerType = DofOwner::UNKNOWN;
      ElemData[i].fIsUnique = flag(0);
      ElemData[i].fToBeUsedInAssembly = flag(0);
      ElemData[i].fGeneratesInertialForces = flag(0);
      ElemData[i].fUsesAirProperties = flag(0);
      ElemData[i].fDefaultOut = fDefaultOut; /* Da "output.h" */
      ElemData[i].OutFile = OutputHandler::UNKNOWN; /* Da "output.h" */
   }
   
   /* Se un tipo possiede Dof propri aggiungere qui il tipo di Dof */
   ElemData[Elem::JOINT].DofOwnerType = DofOwner::JOINT;
   ElemData[Elem::ELECTRIC].DofOwnerType = DofOwner::ELECTRIC;
   ElemData[Elem::ELECTRICBULK].DofOwnerType = DofOwner::ELECTRICBULK;
   ElemData[Elem::GENEL].DofOwnerType = DofOwner::GENEL;
   ElemData[Elem::ROTOR].DofOwnerType = DofOwner::ROTOR;
   ElemData[Elem::AEROMODAL].DofOwnerType = DofOwner::AEROMODAL;
   ElemData[Elem::HYDRAULIC].DofOwnerType = DofOwner::HYDRAULIC;
   ElemData[Elem::LOADABLE].DofOwnerType = DofOwner::LOADABLE;
   
   /* Se un tipo scrive su un file di output, aggiungere qui il tipo di file */
   ElemData[Elem::AUTOMATICSTRUCTURAL].OutFile = OutputHandler::INERTIA;
   ElemData[Elem::JOINT].OutFile = OutputHandler::JOINTS;
   ElemData[Elem::FORCE].OutFile = OutputHandler::FORCES;
   ElemData[Elem::BEAM].OutFile = OutputHandler::BEAMS;
   ElemData[Elem::ROTOR].OutFile = OutputHandler::ROTORS;   
   ElemData[Elem::AEROMODAL].OutFile = OutputHandler::AEROMODALS;   
   ElemData[Elem::AERODYNAMIC].OutFile = OutputHandler::AERODYNAMIC;
   ElemData[Elem::HYDRAULIC].OutFile = OutputHandler::HYDRAULIC;
   ElemData[Elem::LOADABLE].OutFile = OutputHandler::LOADABLE;
   ElemData[Elem::GENEL].OutFile = OutputHandler::GENELS;
   ElemData[Elem::EXTERNAL].OutFile = OutputHandler::EXTERNALS;
   
   /* Tabella delle derivazioni */
   ElemData[Elem::AUTOMATICSTRUCTURAL].iDerivation = ELEM;
   ElemData[Elem::GRAVITY].iDerivation = ELEM;
   ElemData[Elem::BODY].iDerivation = ELEM | GRAVITYOWNER | INITIALASSEMBLY;
   ElemData[Elem::JOINT].iDerivation = ELEM | DOFOWNER | INITIALASSEMBLY;
   ElemData[Elem::GENEL].iDerivation = ELEM | DOFOWNER;
   ElemData[Elem::FORCE].iDerivation = ELEM | INITIALASSEMBLY;
   ElemData[Elem::BEAM].iDerivation = ELEM | GRAVITYOWNER | INITIALASSEMBLY;
   ElemData[Elem::PLATE].iDerivation = ELEM | INITIALASSEMBLY;
   ElemData[Elem::AIRPROPERTIES].iDerivation = ELEM | INITIALASSEMBLY;
   ElemData[Elem::ROTOR].iDerivation = ELEM | DOFOWNER | AIRPROPOWNER;
   ElemData[Elem::AEROMODAL].iDerivation = ELEM | DOFOWNER | AIRPROPOWNER;
   ElemData[Elem::AERODYNAMIC].iDerivation = ELEM | AIRPROPOWNER |INITIALASSEMBLY;
   ElemData[Elem::ELECTRIC].iDerivation = ELEM | DOFOWNER;
   ElemData[Elem::ELECTRICBULK].iDerivation = ELEM | DOFOWNER;
   ElemData[Elem::BULK].iDerivation = ELEM;
   ElemData[Elem::LOADABLE].iDerivation = ELEM | GRAVITYOWNER | INITIALASSEMBLY | DOFOWNER;   
   ElemData[Elem::DRIVEN].iDerivation = ELEM;
   ElemData[Elem::AEROMODAL].iDerivation = ELEM | AIRPROPOWNER;

   /* Aggiungere qui il flag di elemento unico */
   ElemData[Elem::GRAVITY].fIsUnique = flag(1);
   ElemData[Elem::AIRPROPERTIES].fIsUnique = flag(1);

   /* Aggiungere qui se un tipo deve essere usato di default 
    * nell'assemblaggio iniziale */
   ElemData[Elem::JOINT].fToBeUsedInAssembly = flag(1);
   ElemData[Elem::BEAM].fToBeUsedInAssembly = flag(1);  
   
   /* Aggiungere qui se un tipo genera forze d'inerzia e quindi
    * deve essere collegato all'elemento accelerazione di gravita' */
   ElemData[Elem::BODY].fGeneratesInertialForces = flag(1);
   ElemData[Elem::LOADABLE].fGeneratesInertialForces = flag(1);
   
   /* Aggiungere qui se un tipo usa le proprieta' dell'aria e quindi
    * deve essere collegato all'elemento proprieta' dell'aria */
   ElemData[Elem::ROTOR].fUsesAirProperties = flag(1);
   ElemData[Elem::AEROMODAL].fUsesAirProperties = flag(1);
   ElemData[Elem::AERODYNAMIC].fUsesAirProperties = flag(1);
   ElemData[Elem::LOADABLE].fUsesAirProperties = flag(1);
   ElemData[Elem::EXTERNAL].fUsesAirProperties = flag(1);
   
   
   /* Reset della struttura DriveData */
   for(int i = 0; i < Drive::LASTDRIVETYPE; i++) {
      DriveData[i].ppFirstDrive = NULL;
      DriveData[i].iNum = 0;
   }
}


/* distruttore */
void DataManager::ElemManagerDestructor(void)
{
   DEBUGCOUT("Entering DataManager::ElemManagerDestructor()" << std::endl);
   
   /* Distruzione matrici di lavoro per assemblaggio */ 
   if(pWorkMatB != NULL) {
      DEBUGCOUT("deleting assembly structure, SubMatrix B" << std::endl);
      SAFEDELETE(pWorkMatB);
   }   
   
   if(pWorkMatA != NULL) {
      DEBUGCOUT("deleting assembly structure, SubMatrix A" << std::endl);
      SAFEDELETE(pWorkMatA);
   }   
   
   
   if(pdWorkMat != NULL) {
      DEBUGCOUT("deleting assembly structure, double workspace" << std::endl);
      SAFEDELETEARR(pdWorkMat);
   }   

   // ASSERT(piWorkIndex != NULL);

   if(piWorkIndex != NULL) {
      DEBUGCOUT("deleting assembly structure, integer workspace" << std::endl);
      SAFEDELETEARR(piWorkIndex);
   }   

   /* Distruzione elementi */
   ASSERT(ppElems != NULL);

   if (ppElems != NULL) {
#if defined(USE_ELEM_ITER)
      /* Usa l'iteratore built-in */
       Elem* p = NULL;
       if (ElemIter.fGetFirst(p)) {
	  do {
	     ASSERT(p != NULL);
	     if(p != NULL) {		  
		DEBUGCOUT("deleting element " << p->GetLabel()
			  << ", type " << psElemNames[p->GetElemType()]
			  << std::endl);
		SAFEDELETE(p);
	     }
	  } while (ElemIter.fGetNext(p));
       }
#else /* !USE_ELEM_ITER */
      Elem** pp = ppElems;
      while(pp < ppElems+iTotElem) {
	 ASSERT(*pp != NULL);
	 if(*pp != NULL) {		  
	    DEBUGCOUT("deleting element " << (*pp)->GetLabel() 
		      << ", type " << psElemNames[(*pp)->GetElemType()]
		      << std::endl);
	    SAFEDELETE(*pp);		  
	 }	     
	 pp++;
      }
#endif /* !USE_ELEM_ITER */
      DEBUGCOUT("deleting elements structure" << std::endl);
      SAFEDELETEARR(ppElems);
   }   
     
   /* Distruzione drivers */
   if(ppDrive != NULL) {
      Drive** pp = ppDrive;
      while(pp < ppDrive+iTotDrive) {
	 if(*pp != NULL) {
	    DEBUGCOUT("deleting driver " << (*pp)->GetLabel() << ", type " 
		      << psDriveNames[(*pp)->GetDriveType()] << std::endl);
	    SAFEDELETE(*pp);
	 }	     
	 pp++;
      }
      
      DEBUGCOUT("deleting drivers structure" << std::endl);
      SAFEDELETEARR(ppDrive);
   }
      
}


/* Inizializza la struttura dei dati degli elementi ed alloca
 * l'array degli elementi */
void DataManager::ElemDataInit(void)
{
   /* struttura degli elementi */
   for (int iCnt = 0; iCnt < Elem::LASTELEMTYPE; iCnt++) {	   
      iTotElem += ElemData[iCnt].iNum;
   }	
   
   DEBUGCOUT("iTotElem = " << iTotElem << std::endl);
   
   if (iTotElem > 0) {	     
      SAFENEWARR(ppElems, Elem*, iTotElem);
      
      /* Inizializza l'iteratore degli elementi usato all'interno
       * dell'ElemManager */
      ElemIter.Init(ppElems, iTotElem);

      Elem** ppTmp = ppElems;
      while (ppTmp < ppElems+iTotElem) {	      
	 *ppTmp++ = NULL;
      }
      
      ElemData[0].ppFirstElem = ppElems;
      for (int iCnt = 0; iCnt < Elem::LASTELEMTYPE-1; iCnt++) {	      
	 ElemData[iCnt+1].ppFirstElem =
	   ElemData[iCnt].ppFirstElem+ElemData[iCnt].iNum;
      }
      
   } else {
      std::cerr << "warning, no elements are defined" << std::endl;
   }	

		
   /* struttura dei drivers */
   for (int iCnt = 0; iCnt < Drive::LASTDRIVETYPE; iCnt++) {	   
      iTotDrive += DriveData[iCnt].iNum;
   }	
   
   DEBUGCOUT("iTotDrive = " << iTotDrive << std::endl);
   
   if (iTotDrive > 0) {
      SAFENEWARR(ppDrive, Drive*, iTotDrive);

      /* Puo' essere sostituito con un opportuno iteratore:
       VecIter<Drive*> DriveIter(ppDrive, iTotDrive);
       Drive* pTmp = NULL;
       if (DriveIter.fGetFirst(pTmp)) {	       
	  do {
	     pTmp = NULL;
	  } while (DriveIter.fGetNext(pTmp));
       }
       */
      
      Drive** ppTmp = ppDrive;
      while (ppTmp < ppDrive+iTotDrive) {	      
	 *ppTmp++ = NULL;
      }	   
      
      DriveData[0].ppFirstDrive = ppDrive;
      for (int iCnt = 0; iCnt < Drive::LASTDRIVETYPE-1; iCnt++) {	      
	 DriveData[iCnt+1].ppFirstDrive =
	   DriveData[iCnt].ppFirstDrive +
	   DriveData[iCnt].iNum;
      }
      
   } else {
      DEBUGCERR("warning, no drivers are defined" << std::endl);
   }	   
}


/* Prepara per l'assemblaggio */
void DataManager::ElemAssInit(void)
{
   ASSERT(iMaxWorkNumRows > 0);
   ASSERT(iMaxWorkNumCols > 0);
   
   /* Per evitare problemi, alloco tanto spazio quanto necessario per
    * scrivere in modo sparso la matrice piu' grande
    * 
    * Nota: le dimensioni sono state moltiplicate per due per
    * poter creare due matrici (in quanto la seconda e' 
    * richiesta per gli autovalori)
    */
   iWorkIntSize = 2*iMaxWorkNumRows*iMaxWorkNumCols;
   iWorkDoubleSize = iMaxWorkNumRows*iMaxWorkNumCols;

   if (iWorkIntSize > 0) {
      SAFENEWARR(piWorkIndex, integer, 2*iWorkIntSize);
      SAFENEWARR(pdWorkMat, doublereal, 2*iWorkDoubleSize);      
      
      /* SubMatrixHandlers */
      SAFENEWWITHCONSTRUCTOR(pWorkMatA,
			     VariableSubMatrixHandler,
			     VariableSubMatrixHandler(iWorkIntSize,
						      iWorkDoubleSize,
						      piWorkIndex,
						      pdWorkMat));
      
      SAFENEWWITHCONSTRUCTOR(pWorkMatB,
			     VariableSubMatrixHandler,
			     VariableSubMatrixHandler(iWorkIntSize,
						      iWorkDoubleSize,
						      piWorkIndex+iWorkIntSize, 
						      pdWorkMat+iWorkDoubleSize));
      
      DEBUGCOUT("Creating working matrices: work int size = " << iWorkIntSize 
		<< ", work double size = " << iWorkDoubleSize << std::endl
		<< "work matrices are " << iMaxWorkNumRows 
		<< " x " << iMaxWorkNumCols << std::endl);
      
   } else {
      std::cerr << "warning, null size of working matrix" << std::endl;
   }      
}


/* Assemblaggio dello jacobiano.
 * Di questa routine e' molto importante l'efficienza, quindi vanno valutate
 * correttamente le opzioni. */

void DataManager::AssJac(MatrixHandler& JacHdl, doublereal dCoef)
{
   DEBUGCOUT("Entering DataManager::AssJac()" << std::endl);

   ASSERT(pWorkMatA != NULL);   
   ASSERT(ppElems != NULL);

#if !defined(USE_ELEM_ITER)
   /* ciclo sugli elementi */
   for(Elem** ppTmpEl = ppElems; ppTmpEl < ppElems+iTotElem; ppTmpEl++) {
      /* Nuova versione, piu' compatta.
       * La funzione propria AssJac, comune a tutti gli elementi,
       * scrive nella WorkMat (passata come reference) il contributo
       * dell'elemento allo jacobiano e restituisce un reference 
       * alla workmat stessa, che viene quindi sommata allo jacobiano.
       * Ogni elemento deve provvedere al resizing della WorkMat e al
       * suo reset ove occorra */
      
      
      
      
      /* il SubMatrixHandler e' stato modificato in modo da essere
       * in grado di trasformarsi agevolmente da Full a Sparse e quindi 
       * viene gestito in modo automatico, e del tutto trasparente, 
       * il passaggio da un modo all'altro. L'elemento switcha la matrice
       * nel modo che ritiene opportuno; l'operatore += capisce di quale
       * matrice si sta occupando ed agisce di conseguenza.
       */
      
      /* Con VariableSubMatrixHandler */
      JacHdl += (*ppTmpEl)->AssJac(*pWorkMatA, dCoef, *pXCurr, *pXPrimeCurr);
   }	   
#else /* USE_ELEM_ITER */
   
   /* Versione con iteratore: */
    Elem* pTmpEl = NULL;
    if (ElemIter.fGetFirst(pTmpEl)) {
       do {	
	  JacHdl += pTmpEl->AssJac(*pWorkMatA, dCoef, *pXCurr, *pXPrimeCurr);
       } while (ElemIter.fGetNext(pTmpEl));    
    }   
#endif /* USE_ELEM_ITER */
}


void DataManager::AssEig(MatrixHandler& A_Hdl, MatrixHandler& B_Hdl)
{
   DEBUGCOUT("Entering DataManager::AssEig()" << std::endl);
      
   ASSERT(pWorkMatA != NULL);
   ASSERT(pWorkMatB != NULL);
   ASSERT(ppElems != NULL);
         
#if !defined(USE_ELEM_ITER)
   /* ciclo sugli elementi */
   for(Elem** ppTmpEl = ppElems; ppTmpEl < ppElems+iTotElem; ppTmpEl++) {
      /* Nuova versione, piu' compatta.
       * La funzione propria AssJac, comune a tutti gli elementi,
       * scrive nella WorkMat (passata come reference) il contributo
       * dell'elemento allo jacobiano e restituisce un reference 
       * alla workmat stessa, che viene quindi sommata allo jacobiano.
       * Ogni elemento deve provvedere al resizing della WorkMat e al
       * suo reset ove occorra */
      
      
      
      
      /* il SubMatrixHandler e' stato modificato in modo da essere
       * in grado di trasformarsi agevolmente da Full a Sparse e quindi 
       * viene gestito in modo automatico, e del tutto trasparente, 
       * il passaggio da un modo all'altro. L'elemento switcha la matrice
       * nel modo che ritiene opportuno; l'operatore += capisce di quale
       * matrice si sta occupando ed agisce di conseguenza.
       */
      
      /* Con VariableSubMatrixHandler */
      (*ppTmpEl)->AssMats(*pWorkMatA, *pWorkMatB, *pXCurr, *pXPrimeCurr);
      A_Hdl += *pWorkMatA;
      B_Hdl -= *pWorkMatB;
   }	   
#else /* USE_ELEM_ITER */   
   
   /* Versione con iteratore: */
    Elem* pTmpEl = NULL;
    if (ElemIter.fGetFirst(pTmpEl)) {
       do {		 
	  pTmpEl->AssMats(*pWorkMatA, *pWorkMatB, *pXCurr, *pXPrimeCurr);
	  A_Hdl += *pWorkMatA;
	  B_Hdl += *pWorkMatB;
       } while (ElemIter.fGetNext(pTmpEl));    
    }   
#endif /* USE_ELEM_ITER */
}


/* Assemblaggio del residuo */

void DataManager::AssRes(VectorHandler& ResHdl, doublereal dCoef)
{
   DEBUGCOUT("Entering AssRes()" << std::endl);

   /* Vedi quanto scritto per lo jacobiano */
   
   ASSERT(iWorkIntSize >= iWorkDoubleSize);
   MySubVectorHandler WorkVec(iWorkDoubleSize, piWorkIndex, pdWorkMat);
      
#if !defined(USE_ELEM_ITER)
   for (Elem** ppTmpEl = ppElems; ppTmpEl < ppElems+iTotElem; ppTmpEl++) {
      ResHdl += (*ppTmpEl)->AssRes(WorkVec, dCoef, *pXCurr, *pXPrimeCurr);
   }
#else /* USE_ELEM_ITER */
   /* Versione con iteratore: */
    Elem* pTmpEl = NULL;
    if (ElemIter.fGetFirst(pTmpEl)) {       
       do {
	  DEBUGCOUT(psElemNames[pTmpEl->GetElemType()]
		    << "(" << pTmpEl->GetLabel() << ")" << std::endl);
#if 0
	  std::cout << psElemNames[pTmpEl->GetElemType()]
	    << "(" << pTmpEl->GetLabel() << ")" << std::endl;
#endif
	  ResHdl += pTmpEl->AssRes(WorkVec, dCoef, *pXCurr, *pXPrimeCurr);
       } while (ElemIter.fGetNext(pTmpEl));
    }
#endif /* USE_ELEM_ITER */
}


void DataManager::ElemOutput(OutputHandler& OH) const
{
#if !defined(USE_ELEM_ITER)
   for (Elem** ppTmpEl = ppElems; ppTmpEl < ppElems+iTotElem; ppTmpEl++) {      
      (*ppTmpEl)->Output(OH);
   }   
#else /* USE_ELEM_ITER */
   /* Versione con iteratore: */
    Elem* pTmpEl = NULL;
    if (ElemIter.fGetFirst(pTmpEl)) {       
       do {	
	  pTmpEl->Output(OH);
       } while (ElemIter.fGetNext(pTmpEl));
    }
#endif /* USE_ELEM_ITER */
}


void
DataManager::ElemOutput(
		OutputHandler& OH,
		const VectorHandler& X,
		const VectorHandler& XP
		) const
{
#if !defined(USE_ELEM_ITER)
   for (Elem** ppTmpEl = ppElems; ppTmpEl < ppElems+iTotElem; ppTmpEl++) {      
      (*ppTmpEl)->Output(OH, X, XP);
   }   
#else /* USE_ELEM_ITER */
   /* Versione con iteratore: */
    Elem* pTmpEl = NULL;
    if (ElemIter.fGetFirst(pTmpEl)) {       
       do {	
	  pTmpEl->Output(OH, X, XP);
       } while (ElemIter.fGetNext(pTmpEl));
    }
#endif /* USE_ELEM_ITER */
}


void
DataManager::ElemOutput_pch(
		std::ostream& pch
		) const
{
#if !defined(USE_ELEM_ITER)
   for (Elem** ppTmpEl = ppElems; ppTmpEl < ppElems+iTotElem; ppTmpEl++) {      
      (*ppTmpEl)->Output_pch(pch);
   }   
#else /* USE_ELEM_ITER */
   /* Versione con iteratore: */
    Elem* pTmpEl = NULL;
    if (ElemIter.fGetFirst(pTmpEl)) {       
       do {	
	  pTmpEl->Output_pch(pch);
       } while (ElemIter.fGetNext(pTmpEl));
    }
#endif /* USE_ELEM_ITER */
}


void
DataManager::ElemOutput_f06(
		std::ostream& f06,
		const VectorHandler& X
		) const
{
#if !defined(USE_ELEM_ITER)
   for (Elem** ppTmpEl = ppElems; ppTmpEl < ppElems+iTotElem; ppTmpEl++) {      
      (*ppTmpEl)->Output_f06(f06, X);
   }   
#else /* USE_ELEM_ITER */
   /* Versione con iteratore: */
    Elem* pTmpEl = NULL;
    if (ElemIter.fGetFirst(pTmpEl)) {       
       do {	
	  pTmpEl->Output_f06(f06, X);
       } while (ElemIter.fGetNext(pTmpEl));
    }
#endif /* USE_ELEM_ITER */
}


void
DataManager::ElemOutput_f06(
		std::ostream& f06,
		const VectorHandler& Xr,
		const VectorHandler& Xi
		) const
{
#if !defined(USE_ELEM_ITER)
   for (Elem** ppTmpEl = ppElems; ppTmpEl < ppElems+iTotElem; ppTmpEl++) {      
      (*ppTmpEl)->Output_f06(f06, Xr, Xi);
   }   
#else /* USE_ELEM_ITER */
   /* Versione con iteratore: */
    Elem* pTmpEl = NULL;
    if (ElemIter.fGetFirst(pTmpEl)) {       
       do {	
	  pTmpEl->Output_f06(f06, Xr, Xi);
       } while (ElemIter.fGetNext(pTmpEl));
    }
#endif /* USE_ELEM_ITER */
}


/* cerca un elemento qualsiasi */
void* DataManager::pFindElem(Elem::Type Typ, unsigned int uL) const
{
   ASSERT(ElemData[Typ].ppFirstElem != NULL);
   ASSERT(ElemData[Typ].iNum > 0);
   ASSERT(uL > 0);
   
   Elem* p = pLabelSearch(ElemData[Typ].ppFirstElem, ElemData[Typ].iNum, uL);
      
   if( p == NULL ) {
      return NULL;
   }
   
   ASSERT(p->pGetElem() != NULL);
   return p->pGetElem();
}


/* cerca un elemento qualsiasi */
Elem** DataManager::ppFindElem(Elem::Type Typ, unsigned int uL) const
{
   ASSERT(ElemData[Typ].ppFirstElem != NULL);
   ASSERT(ElemData[Typ].iNum > 0);
   ASSERT(uL > 0);
   
   int i = LabelSearch(ElemData[Typ].ppFirstElem, ElemData[Typ].iNum, uL);
      
   if (i < 0) {
      return NULL;
   }

   return ElemData[Typ].ppFirstElem+i;
}


/* cerca un elemento qualsiasi */
void* DataManager::pFindElem(Elem::Type Typ, unsigned int uL,
			     unsigned int iDeriv) const
{
   ASSERT(ElemData[Typ].ppFirstElem != NULL);
   ASSERT(ElemData[Typ].iNum > 0);
   ASSERT(uL > 0);
   ASSERT(iDeriv == int(ELEM) || ElemData[Typ].iDerivation & iDeriv);
      
   Elem* p = pLabelSearch(ElemData[Typ].ppFirstElem, ElemData[Typ].iNum, uL);
      
   if( p == NULL ) {
      return NULL;
   }
   
   return pChooseElem(p, iDeriv);
}


/* Usata dalle due funzioni precedenti */
void* DataManager::pChooseElem(Elem* p, unsigned int iDeriv) const
{
   ASSERT(p != NULL);

   switch(iDeriv) {
    case ELEM:
      ASSERT(p->pGetElem() != NULL);
      return p->pGetElem();
      
    case DOFOWNER:
      ASSERT(p->pGetElemWithDofs() != NULL);
      return p->pGetElemWithDofs();
      
    case GRAVITYOWNER:
      ASSERT(p->pGetElemGravityOwner() != NULL);
      return p->pGetElemGravityOwner();
      
    case AIRPROPOWNER:
      ASSERT(p->pGetAerodynamicElem() != NULL);
      return p->pGetAerodynamicElem();

    case INITIALASSEMBLY:
      ASSERT(p->pGetInitialAssemblyElem() != NULL);
      return p->pGetInitialAssemblyElem();
   }
   
   /* default */
   return NULL;
}


/* cerca un drive qualsiasi */
void* DataManager::pFindDrive(Drive::Type Typ, unsigned int uL) const
{
   ASSERT(DriveData[Typ].ppFirstDrive != NULL);
   ASSERT(DriveData[Typ].iNum > 0);
   ASSERT(uL > 0);
   
   Drive* p = pLabelSearch(DriveData[Typ].ppFirstDrive, DriveData[Typ].iNum, uL);

   if (p == NULL) {
      return NULL;
   }

   return p;
}


flag DataManager::fGetDefaultOutputFlag(const Elem::Type& t) const
{
   return ElemData[t].fDefaultOut;
}

/* DataManager - end */


/* InitialAssemblyIterator - begin */

InitialAssemblyIterator::
InitialAssemblyIterator(const DataManager::ElemDataStructure
			(*pED)[Elem::LASTELEMTYPE])
: pElemData(pED), 
FirstType(Elem::UNKNOWN), ppFirst(NULL), 
CurrType(Elem::UNKNOWN), ppCurr(NULL)
{
   int iCnt = 0;
   while ((*pElemData)[iCnt].fToBeUsedInAssembly == 0
	 || (*pElemData)[iCnt].iNum == 0) {
      if (++iCnt >= Elem::LASTELEMTYPE) {	 
	 break;
      }      
   }   
	   
   ASSERT(iCnt < Elem::LASTELEMTYPE);   
   ASSERT((*pElemData)[iCnt].ppFirstElem != NULL);
   
   ((Elem**&)ppFirst) = (*pElemData)[iCnt].ppFirstElem;
   ppCurr = (Elem**)ppFirst;
   (Elem::Type&)FirstType = CurrType = Elem::Type(iCnt);     
}

InitialAssemblyElem* InitialAssemblyIterator::GetFirst(void) const
{
   (Elem::Type&)CurrType = FirstType;
   ((Elem**&)ppCurr) = (Elem**)ppFirst;

   /* La variabile temporanea e' necessaria per il debug. */
   InitialAssemblyElem* p = (*ppCurr)->pGetInitialAssemblyElem();
   
#ifdef DEBUG   
   ASSERT(p != NULL);
   if (p == NULL) {
      std::cerr << "warning, element " << (*ppCurr)->GetLabel() 
	<< " is not subject to initial assembly" << std::endl;
   }
#endif
      
   return p;
}


InitialAssemblyElem* InitialAssemblyIterator::GetNext(void) const
{
   ((Elem**&)ppCurr)++;
   if (ppCurr >= (*pElemData)[CurrType].ppFirstElem
      +(*pElemData)[CurrType].iNum) {
      int iCnt = int(CurrType);
      
      do {
	 if (++iCnt >= Elem::LASTELEMTYPE) {	    
	    return NULL;
	 }	 
      } while ((*pElemData)[iCnt].fToBeUsedInAssembly == 0
	      || (*pElemData)[iCnt].iNum == 0);	
	
      ASSERT((*pElemData)[iCnt].ppFirstElem != NULL);
      (Elem::Type&)CurrType = Elem::Type(iCnt);	
      (Elem**&)ppCurr = (*pElemData)[iCnt].ppFirstElem;
      
      /* La variabile temporanea e' necessaria per il debug. */
      InitialAssemblyElem* p = (*ppCurr)->pGetInitialAssemblyElem();
      
#ifdef DEBUG   
      ASSERT(p != NULL);
      if (p == NULL) {
	 std::cerr << "warning, element " << (*ppCurr)->GetLabel() 
	   << " is not subject to initial assembly" << std::endl;
      }
#endif
      
      return p;
   }
   
   /* La variabile temporanea e' necessaria per il debug. */
   InitialAssemblyElem* p = (*ppCurr)->pGetInitialAssemblyElem();
   
#ifdef DEBUG   
   ASSERT(p != NULL);
   if (p == NULL) {
      std::cerr << "warning, element " << (*ppCurr)->GetLabel() 
	<< " is not subject to initial assembly" << std::endl;
   }
#endif
      
   return p;
}

/* InitialAssemblyIterator - end */

