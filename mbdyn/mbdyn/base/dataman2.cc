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

/* DataManager - 
 * continua qui perche' il file dataman.cc sta diventando lungo */

#include <mbconfig.h>

#include <dataman.h>
#include <dataman_.h>

#include <gravity.h>
#include <memmans.h>
#include <harwrap.h>
#ifdef USE_MESCHACH
#include <mschwrap.h>
#endif /* USE_MESCHACH */

/* DataManager - continue */

/* Setta il valore della variabile Time nel DataManager 
 * usato dal metodo numerico all'inizio di ogni step temporale */

void DataManager::SetTime(doublereal dTime)
{   
   /* Setta la variabile Time nella tabella dei simboli */   
   ASSERT(pTime != NULL);
   pTime->SetVal(dTime);

   DEBUGLCOUT(MYDEBUG_INIT|MYDEBUG_ASSEMBLY, 
	      "GlobalSymbolTable:" << endl
	      << GlobalSymbolTable << endl);
   
   /* Setta il tempo nel DriveHandler */
   DrvHdl.SetTime(dTime);
   
   /* serve i drive pending */
   for (int iType = 0; iType < DriveType::LASTDRIVETYPE; iType++) {
      for (unsigned int iCnt = 0; iCnt < DriveData[iType].iNum; iCnt++) {
	 DriveData[iType].ppFirstDrive[iCnt]->ServePending();
      }
   }
      
} /* End of DataManager::SetTime() */


/* Collega il DataManager ed il DriveHandler alla soluzione */
void DataManager::LinkToSolution(const VectorHandler& XCurr, 
				 const VectorHandler& XPrimeCurr)
{
   (VectorHandler*&)pXCurr = (VectorHandler*)&XCurr;
   (VectorHandler*&)pXPrimeCurr = (VectorHandler*)&XPrimeCurr;
   DrvHdl.LinkToSolution(XCurr, XPrimeCurr);
}  


/* Inizializzatore dei dof di ogni elemento */
void DataManager::DofOwnerInit(void)
{
   DEBUGCOUTFNAME("DataManager::DofOwnerInit");
   ASSERT(pDofs != NULL);
   ASSERT(ppNodes != NULL);
   ASSERT(ppElems != NULL);
      
   /* per ogni nodo */
   Node** ppNd = ppNodes;
   while (ppNd < ppNodes+iTotNodes) {	
      ASSERT(*ppNd != NULL);
      DEBUGLCOUT(MYDEBUG_INIT|MYDEBUG_ASSEMBLY, 
		 "Node type " << (*ppNd)->GetNodeType()
		 << " (" << psNodeNames[(*ppNd)->GetNodeType()] 
		 << "(" << (*ppNd)->GetLabel() << "))" << endl);
      
      unsigned int iNumDof;
      
      /* chiede al nodo quanti dof possiede */
      if ((iNumDof = (*ppNd)->iGetNumDof()) > 0) {	  	     
	 /* si fa passare il primo Dof */
	 Dof* pDf = pDofs+(*ppNd)->iGetFirstIndex();
	 
	 DEBUGLCOUT(MYDEBUG_INIT|MYDEBUG_ASSEMBLY,
		    psNodeNames[(*ppNd)->GetNodeType()] 
		    << "(" << (*ppNd)->GetLabel()
		    << "): first dof = " << pDf->iIndex+1 << endl);
	 
	 /* per ogni Dof, chiede al nodo di che tipo e' e lo 
	  * setta nel DofOwner */
	 for (unsigned int iCnt = 0; iCnt < iNumDof; iCnt++) {
	    /* (Dof*)->Order = */
	    (pDf+iCnt)->Order = (*ppNd)->SetDof(iCnt);
	 }
      }	
      ppNd++;
   }
   
   
   /* per ogni elemento */
#if !defined(USE_ELEM_ITER)
   Elem** ppLastElem = ppElems+iTotElem;
   for (Elem** ppEl = ppElems; ppEl < ppLastElem; ppEl++) {	
      ASSERT(*ppEl != NULL);
      DEBUGLCOUT(MYDEBUG_INIT|MYDEBUG_ASSEMBLY, 
		 "Elem type " << (*ppEl)->GetElemType() 
		 << " (" << psElemNames[(*ppEl)->GetElemType()] 
		 << "(" << (*ppEl)->GetLabel() << "))" << endl);
      
      unsigned int iNumDof;
	
      /* chiede all'elemento quanti dof possiede */
      if ((iNumDof = (*ppEl)->iGetNumDof()) > 0) {
	 ElemWithDofs* pElWD = (ElemWithDofs*)(*ppEl)->pGetElemWithDofs();
	 ASSERT(pElWD != NULL);
	 
	 /* si fa passare il DofOwner */
	 Dof* pDf = pDofs+pElWD->iGetFirstIndex();
	 
	 DEBUGLCOUT(MYDEBUG_INIT|MYDEBUG_ASSEMBLY, 
		    "element " << pElWD->GetLabel()
		    << ": first dof = " << pDf->iIndex+1 << endl);
	 
	 /* per ogni Dof, chiede all'elemento di che tipo e' e lo 
	  * setta nel DofOwner */
	 for (unsigned int iCnt = 0; iCnt < iNumDof; iCnt++) {
	    /* (Dof*)->Order = */
	    (pDf+iCnt)->Order = pElWD->SetDof(iCnt);
	 }
      }	
      ppEl++;
   }
#else /* USE_ELEM_ITER */
   Elem* pEl = NULL;
   if(ElemIter.fGetFirst(pEl)) {
      do {
	 ASSERT(pEl != NULL);
	 DEBUGLCOUT(MYDEBUG_INIT|MYDEBUG_ASSEMBLY,
		    "Elem type " << pEl->GetElemType() 
		    << " (" << psElemNames[pEl->GetElemType()] 
		    << "(" << pEl->GetLabel() << "))" << endl);
	 
	 unsigned int iNumDof;
	
	 /*
	 pEl = pEl->pGet();
	  */
	 
	 /* chiede all'elemento quanti dof possiede */
	 if ((iNumDof = pEl->iGetNumDof()) > 0) {
	    ElemWithDofs* pElWD = (ElemWithDofs*)pEl->pGetElemWithDofs();
	    ASSERT(pElWD != NULL);
	    
	    /* si fa passare il DofOwner */
	    Dof* pDf = pDofs+pElWD->iGetFirstIndex();
	    
	    DEBUGLCOUT(MYDEBUG_INIT|MYDEBUG_ASSEMBLY,
		       psElemNames[pEl->GetElemType()] 
		       << "(" << pElWD->GetLabel()
		       << "): first dof = " << pDf->iIndex+1 << endl);
	    
	    /* per ogni Dof, chiede all'elemento di che tipo e' e lo 
	     * setta nel DofOwner */
	    for (unsigned int iCnt = 0; iCnt < iNumDof; iCnt++) {
	       /* (Dof*)->Order = */
	       pDf[iCnt].Order = pElWD->SetDof(iCnt);
	    }
	 }
      } while (ElemIter.fGetNext(pEl));
   }
#endif /* USE_ELEM_ITER */
} /* End of DataManager::DofOwnerInit() */


/* Inizializzazione della struttura dei dof 
 * per l'assemblaggio iniziale dei vincoli */
#if defined(USE_STRUCT_NODES)
void DataManager::InitialJointAssembly(void)
{
   /* Costruisce la struttura temporanea dei Dof */
   
   ASSERTMSG(DofData[DofType::JOINT].iNum > 0, 
	     "Warning, no joints are defined; You shouldn't have reached this point");
   ASSERT(DofData[DofType::STRUCTURALNODE].iNum > 0);

   /* Nodi strutturali: mette gli indici ai DofOwner */
   StructNode** ppFirstNode = 
     (StructNode**)NodeData[NodeType::STRUCTURAL].ppFirstNode;
   integer iNumNodes = NodeData[NodeType::STRUCTURAL].iNum;
   
   StructNode** ppNode = ppFirstNode;     
   DofOwner* pTmp = DofData[DofType::STRUCTURALNODE].pFirstDofOwner;
   
   integer iIndex = 0;    /* Indice dei gradi di liberta' */
   unsigned int iNumDofs = 0;  /* numero di dof di un owner */
   for(int iCnt = 1; 
       pTmp < DofData[DofType::STRUCTURALNODE].pFirstDofOwner+
       DofData[DofType::STRUCTURALNODE].iNum;
       iCnt++, pTmp++, ppNode++) {
      iNumDofs = pTmp->iNumDofs = (*ppNode)->iGetInitialNumDof();
      if(iNumDofs > 0) {
	 DEBUGLCOUT(MYDEBUG_INIT|MYDEBUG_ASSEMBLY,
		    "Node " << (*ppNode)->GetLabel() 
		    << ": first index = " << iIndex+1 << endl);
	 pTmp->iFirstIndex = iIndex;
	 iIndex += iNumDofs;
      }
      else {
	 DEBUGCERR("");
	 cerr << "warning, item " << iCnt << " has 0 dofs" << endl;
      }
   }     

   
   /* Elementi: mette gli indici agli eventuali DofOwner */   
   for(int iCnt1 = 0; iCnt1 < ElemType::LASTELEMTYPE; iCnt1++) {
      /* Pre ogni tipo di elemento */
      if(ElemData[iCnt1].fToBeUsedInAssembly && ElemData[iCnt1].iNum > 0) {
	 /* Se deve essere usato nell'assemblaggio e ne sono definiti */
	 
	 /* Tipo di dof dell'elemento corrente */
	 DofType::Type CurrDofType = 
	   ElemData[iCnt1].DofOwnerType;
	 
	 if(CurrDofType != DofType::UNKNOWN) {
	    
	    /* Puntatore al primo DofOwner */
	    pTmp = DofData[CurrDofType].pFirstDofOwner;
	 	 
	    Elem** ppFirstEl = ElemData[iCnt1].ppFirstElem;
	    integer iNumEls = ElemData[iCnt1].iNum;
	    ASSERT(DofData[CurrDofType].iNum == iNumEls);
	  
	    /* Iterazione sugli Elem */
	    Elem** ppEl = ppFirstEl;	     
	    for(int iCnt = 0;
		pTmp < DofData[CurrDofType].pFirstDofOwner+iNumEls;
		pTmp++, ppEl++, iCnt++) {
	       
	       ASSERT((*ppEl)->pGetInitialAssemblyElem() != NULL);
	       iNumDofs = 
		 (*ppEl)->pGetInitialAssemblyElem()->iGetInitialNumDof();
	       if((pTmp->iNumDofs = iNumDofs) > 0) {
		  DEBUGLCOUT(MYDEBUG_INIT|MYDEBUG_ASSEMBLY,
			     "Elem " << (*ppEl)->GetLabel()
			     << " of type \"" << psElemNames[iCnt1]
			     << "\": first index = "
			     << iIndex+1 << endl);
		  pTmp->iFirstIndex = iIndex;
		  iIndex += iNumDofs;
	       }
	       else {
		  DEBUGCERR("");
		  cerr << "warning, item " << iCnt << " has 0 dofs" << endl;
	       }	
	    }
	 }	 
      }
   }
   

   /* Numero totale di Dof durante l'assemblaggio iniziale */
   integer iInitialNumDofs = iIndex;
   
   
   /* Trova le massime dimensioni del workspace 
    * per l'assemblaggio iniziale */
   integer iMaxRows = 0;   
   integer iMaxCols = 0;

   
   InitialAssemblyIterator IAIter(&ElemData);
   InitialAssemblyElem* pEl = IAIter.GetFirst();
   ASSERT(pEl != NULL);
   while(pEl != NULL) {
      integer iCurrRows = 0;
      integer iCurrCols = 0;
      pEl->InitialWorkSpaceDim(&iCurrRows, &iCurrCols);
      if(iCurrRows > iMaxRows) {	 
	 iMaxRows = iCurrRows;
      }      
      if(iCurrCols > iMaxCols) {	 
	 iMaxCols = iCurrCols;
      }      
      pEl = IAIter.GetNext();       
   }



   /* Alla fine, i DofOwner di nodi e joint contengono gli indici giusti per 
    * l'assemblaggio iniziale. Corrispondono a:
    * - per ogni nodo:
    *   - posizione x
    *   - parametri di rotazione g
    *   - velocita' xP
    *   - velocita' angolare omega
    * - per ogni joint:
    *   - se vincolo in posizione, reazione e sua derivata
    *   - se vincolo in velocita', reazione.
    * - per vincoli misti si hanno reazioni ed eventualmente loro derivate 
    *   in base al tipo */

   
   /* Creazione e costruzione array Dof */   
   SAFENEWARR(pDofs, Dof, iInitialNumDofs, DMmm);
   
   iIndex = 0;
   for(Dof* pTmpDof = pDofs; 
       pTmpDof < pDofs+iInitialNumDofs; pTmpDof++) {
      pTmpDof->iIndex = iIndex++;
   }


   /* Alloca spazio per l'assemblaggio */
   integer iIntDim = iMaxRows*iMaxCols*2;
   integer* piWI = NULL;
   SAFENEWARR(piWI, integer, iIntDim, SMmm);
   
   integer iDoubleDim = iMaxRows*iMaxCols;
   doublereal* pdWD = NULL;
   SAFENEWARR(pdWD, doublereal, iDoubleDim, SMmm);   
   
   
   /* Ciclo di iterazioni fino a convergenza */
   
   /* Crea la struttura di calcolo */
   SolutionManager* pSM = NULL;
#ifdef USE_MESCHACH
   if (1) {
      SAFENEWWITHCONSTRUCTOR(pSM,
			     MeschachSparseLUSolutionManager,
			     MeschachSparseLUSolutionManager(iInitialNumDofs, 
							     0, 1.),
			     SMmm);
   } else {
#endif
      SAFENEWWITHCONSTRUCTOR(pSM,
			     HSLUSolutionManager,
			     HSLUSolutionManager(iInitialNumDofs, 0, 1.),
			     SMmm);
#ifdef USE_MESCHACH
   }
#endif
      
   /* Crea il vettore con lo stato del sistema durante l'assemblaggio */
   doublereal* pdX = NULL;
   SAFENEWARR(pdX, doublereal, iInitialNumDofs, SMmm);      
   
#ifdef DEBUG_MEMMANAGER
   DEBUGLCOUT(MYDEBUG_MEM|MYDEBUG_ASSEMBLY, 
	      "After initialisation in InitialJointAssembly" << endl
	      << SMmm << endl
	      // << MHmm << endl
	      << LUmm << endl
	      << DMmm << endl);
#endif   
   
   MyVectorHandler X(iInitialNumDofs, pdX);
   X.Reset(0.);
   
   /* Linka il DriveHandler al vettore soluzione */
   LinkToSolution(X, X);

   /* Setta i valori iniziali dei gradi di liberta' dei nodi strutturali
    * durante l'assemblaggio iniziale */
   for(StructNode** ppTmpNode = ppFirstNode;
       ppTmpNode < ppFirstNode+iNumNodes; ppTmpNode++) {	
      (*ppTmpNode)->SetInitialValue(X);
   }   
     
   
   /* Setta i valori iniziali dei gradi di liberta' dei vincoli
    * durante l'assemblaggio iniziale */
   for(int iCnt1 = 0; iCnt1 < ElemType::LASTELEMTYPE; iCnt1++) {
      /* Pre ogni tipo di elemento */
      if(ElemData[iCnt1].DofOwnerType != DofType::UNKNOWN &&
	 ElemData[iCnt1].fToBeUsedInAssembly &&
	 ElemData[iCnt1].iNum > 0) {
	 
	 Elem** ppFirstEl = ElemData[iCnt1].ppFirstElem;
	 integer iNumEl = ElemData[iCnt1].iNum;
	 for(Elem** ppTmpEl = ppFirstEl;
	     ppTmpEl < ppFirstEl+iNumEl; ppTmpEl++) {
	    
	    ASSERT((*ppTmpEl)->pGetElemWithDofs() != NULL);
	    (*ppTmpEl)->pGetElemWithDofs()->SetInitialValue(X);
	 }
      }	
   }   
   
   
   /* Vettore di lavoro */
   VectorHandler* pResHdl = pSM->pResHdl();
   MySubVectorHandler WorkVec(iDoubleDim, piWI, pdWD);
   
   /* Matrice di lavoro */
   MatrixHandler* pMatHdl = pSM->pMatHdl();
   VariableSubMatrixHandler WorkMat(iIntDim, iDoubleDim, piWI, pdWD);

   /* Soluzione */
   VectorHandler* pSolHdl = pSM->pSolHdl();
   
   /* Ciclo di assemblaggio */
   integer iNumIter = 0;
   while(++iNumIter) {	
      /* Assemblo il residuo */
      pResHdl->Reset(0.);
      
      /* Contributo dei nodi */	
      for(StructNode** ppTmpNode = ppFirstNode; 
	  ppTmpNode < ppFirstNode+iNumNodes; ppTmpNode++) {
	 integer iFirstIndex = ((*ppTmpNode)->pGetDofOwner())->iFirstIndex;
	 
	 /* Nuova feature: ogni nodo ha la sua stiffness */
	 doublereal dPosStiff = (*ppTmpNode)->dGetPositionStiffness();
	 doublereal dVelStiff = (*ppTmpNode)->dGetVelocityStiffness();
	 
	 /* Posizione: k*Delta_x = k(x_0-x) + F */
	 Vec3 TmpVec = (*ppTmpNode)->GetXPrev()-(*ppTmpNode)->GetXCurr();
	 pResHdl->Add(iFirstIndex+1, TmpVec*dPosStiff);
	 
	 /* Rotazione: k*Delta_g = -k*g(R_Delta) + M */
	 Mat3x3 R0 = (*ppTmpNode)->GetRPrev();
	 Mat3x3 RDelta = (*ppTmpNode)->GetRCurr()*R0.Transpose();
	 TmpVec = -gparam(RDelta);
	 pResHdl->Add(iFirstIndex+4, TmpVec*dPosStiff);
	 
	 /* Velocita': k*Delta_v = k*(v0-Delta_v) + F */
	 TmpVec = (*ppTmpNode)->GetVPrev()-(*ppTmpNode)->GetVCurr();
	 pResHdl->Add(iFirstIndex+7, TmpVec*dVelStiff);
	 
	 /* Velocita' angolare: k*(Delta_w+(R_Delta*w0)/\Delta_g) = 
	  *                                    k*(R_Delta*w0-w) + M */
	 Vec3 wPrev((*ppTmpNode)->GetWPrev());
	 Vec3 wCurr((*ppTmpNode)->GetWCurr());
	 
	 
	 if((*ppTmpNode)->fOmegaRotates()) {
	    TmpVec = RDelta*wPrev-wCurr; // con questa la velocita' angolare e' solidale con il nodo
	 } else {	 
	    TmpVec = wPrev-wCurr; // con questa la velocita' angolare e' solidale col riferimento assoluto
	 }
	 
	 pResHdl->Add(iFirstIndex+10, TmpVec*dVelStiff);
	 
      }	
       	
      /* Elementi (con iteratore): */
      pEl = IAIter.GetFirst();
      while(pEl != NULL) {
	 *pResHdl += pEl->InitialAssRes(WorkVec, X);
	 pEl = IAIter.GetNext();
      }	
      
	
#ifdef DEBUG
	if (DEBUG_LEVEL_MATCH(MYDEBUG_ASSEMBLY|MYDEBUG_RESIDUAL)) {
	   /* Output del residuo */
	   cout << "Residual:" << endl;
	   for(int iTmpCnt = 1; iTmpCnt <= iInitialNumDofs; iTmpCnt++) {
	      cout << "Dof" << setw(4) << iTmpCnt << ": "
		<< pResHdl->dGetCoef(iTmpCnt) << endl;
	   }
	}
#endif
	
	
      /* Eseguo il test di convergenza; se e' positivo, esco */
      doublereal dTest = pResHdl->Dot()/(1.+X.Dot());
      if(!isfinite(dTest)) {
	 cerr << "Assembly diverged; aborting ..." << endl;

	 THROW(DataManager::ErrAssemblyDiverged());
      }
      dTest = sqrt(dTest);
      
      DEBUGLCOUT(MYDEBUG_ASSEMBLY, "Iteration: " << iNumIter 
		 << ", Test: " << dTest 
		 << " (Toll = " << dInitialAssemblyToll << ") " << endl);
	
      /* Se la tolleranza e' raggiunta, esce dal ciclo */
      if(dTest <= dInitialAssemblyToll) {
	 DEBUGLCOUT(MYDEBUG_ASSEMBLY, "Initial assembly performed successfully in "
		    << iNumIter << " iterations" << endl);
	 goto endofcycle;
      }
	
	
	
      /* Se ho raggiunto il numero massimo di iterazioni */
      if(iNumIter > iMaxInitialIterations) {
	 cerr
	   << "Initial assembly iterations reached maximum number "
	   << iMaxInitialIterations << "; aborting ..." << endl;

	 THROW(DataManager::ErrAssemblyMaxIters());
      }
      
	
      /* Assemblo lo jacobiano e risolvo */
      pSM->MatrInit(0.);
	
      /* Contributo dei nodi */	
      for(StructNode** ppTmpNode = ppFirstNode; 
	  ppTmpNode < ppFirstNode+iNumNodes; ppTmpNode++) {
	 integer iFirstIndex = ((*ppTmpNode)->pGetDofOwner())->iFirstIndex;
	 
	 /* Nuova feature: ogni nodo ha la sua stiffness */
	 doublereal dPosStiff = (*ppTmpNode)->dGetPositionStiffness();
	 doublereal dVelStiff = (*ppTmpNode)->dGetVelocityStiffness();
	 
	 for(int iCnt = 1; iCnt<= 6; iCnt++) {
	    /* Posizione, rotazione */
	    integer iTmp = iFirstIndex+iCnt;
	    pMatHdl->fPutCoef(iTmp, iTmp, dPosStiff);
	    
	    /* Velocita', velocita' angolare */
	    iTmp += 6;
	    pMatHdl->fPutCoef(iTmp, iTmp, dVelStiff);
	 }
	 

	 if((*ppTmpNode)->fOmegaRotates()) {
	    // con questi la velocita' angolare e' solidale con il nodo
	    
	    /* Velocita' angolare - termine di rotazione: R_Delta*w0/\ */
	    Mat3x3 R0 = (*ppTmpNode)->GetRPrev();
	    Mat3x3 R = (*ppTmpNode)->GetRCurr();
	    Vec3 W0 = (*ppTmpNode)->GetWPrev();
	    Vec3 TmpVec = R*(R0.Transpose()*(W0*dVelStiff));
	    
	    
	    /* W1 in m(3, 2), -W1 in m(2, 3) */
	    doublereal d = TmpVec.dGet(1);
	    pMatHdl->fPutCoef(iFirstIndex+12, iFirstIndex+5, d);
	    pMatHdl->fPutCoef(iFirstIndex+11, iFirstIndex+6, -d);
	    
	    /* W2 in m(1, 3), -W2 in m(3, 1) */
	    d = TmpVec.dGet(2);
	    pMatHdl->fPutCoef(iFirstIndex+10, iFirstIndex+6, d);
	    pMatHdl->fPutCoef(iFirstIndex+12, iFirstIndex+4, -d);
	    
	    /* W3 in m(2, 1), -W3 in m(1, 2) */
	    d = TmpVec.dGet(3);
	    pMatHdl->fPutCoef(iFirstIndex+11, iFirstIndex+4, d);
	    pMatHdl->fPutCoef(iFirstIndex+10, iFirstIndex+5, -d);
	    
	 } // con questi la velocita' angolare e' solidale con il nodo
	 
      }	
      		
      /* Contributo degli elementi */
      pEl = IAIter.GetFirst();
      while(pEl != NULL) {
	 *pMatHdl += pEl->InitialAssJac(WorkMat, X);
	 pEl = IAIter.GetNext();
      }
      


	
      /* Fattorizza e risolve con jacobiano e residuo appena calcolati */
      pSM->Solve();
      
#ifdef DEBUG
      if (DEBUG_LEVEL_MATCH(MYDEBUG_ASSEMBLY|MYDEBUG_RESIDUAL)) {
	 /* Output della soluzione */
	 cout << "Solution:" << endl;
	 for(int iTmpCnt = 1; iTmpCnt <= iInitialNumDofs; iTmpCnt++) {	   
	    cout << "Dof" << setw(4) << iTmpCnt << ": " 
	      << pSolHdl->dGetCoef(iTmpCnt) << endl;
	 }
      }
#endif
	
	
      /* Aggiorno la soluzione */
      X += *pSolHdl;
	
      /* Correggo i nodi */
      for(StructNode** ppTmpNode = ppFirstNode;
	  ppTmpNode < ppFirstNode+iNumNodes; ppTmpNode++) {
	 (*ppTmpNode)->InitialUpdate(X);
      }		     
   }       
   
   
   
endofcycle:   
      
   /* Distrugge il workspace */
   ASSERT(pdWD != NULL);
   if(pdWD != NULL) {	
      SAFEDELETEARR(pdWD, SMmm);
   }
   
   ASSERT(piWI != NULL);
   if(piWI != NULL) {
      SAFEDELETEARR(piWI, SMmm);
   }
      
   
   /* Distrugge il vettore soluzione */
   ASSERT(pdX != NULL);
   if(pdX != NULL) {	
      SAFEDELETEARR(pdX, SMmm);
   }   

   
   /* Resetta e distrugge la struttura temporanea dei Dof */

   /* Elementi: rimette a posto il numero di Dof propri dei vincoli */
   for(int iCnt1 = 0; iCnt1 < ElemType::LASTELEMTYPE; iCnt1++) {
      /* Per ogni tipo di elemento */
      if(ElemData[iCnt1].DofOwnerType != DofType::UNKNOWN &&
	 ElemData[iCnt1].fToBeUsedInAssembly &&
	 ElemData[iCnt1].iNum > 0) {
	 /* Se possiede dofs, se deve essere usato nell'assemblaggio
	  * e se ne sono presenti */
	 
	 /* Tipo di dof dell'elemento corrente */
	 DofType::Type CurrDofType = 
	   ElemData[iCnt1].DofOwnerType;
	 /* Puntatore al primo DofOwner */
	 pTmp = DofData[CurrDofType].pFirstDofOwner;
	 
	 
	 /* Puntatore al primo Elem */
	 Elem** ppFirstEl = ElemData[iCnt1].ppFirstElem;
	 /* Numero di Elem (== al numero di DofOwner) */
	 integer iNumEls = DofData[CurrDofType].iNum;
	 
	 /* Iterazione sugli Elem */
	 Elem** ppEl = ppFirstEl;
	 for(; pTmp < DofData[CurrDofType].pFirstDofOwner+iNumEls;
	     pTmp++, ppEl++) {
	    pTmp->iNumDofs = (*ppEl)->iGetNumDof();
	 }
      }
   }              
   
   
   /* Dealloca il vettore dei Dof */
   ASSERT(pDofs != NULL);
   if(pDofs != NULL) {	
      SAFEDELETEARR((Dof*&)pDofs, DMmm);
   }   
   
} /* End of InitialJointAssembly */
#endif // USE_STRUCT_NODES

/* Aggiorna i DofOwner con il numero di dofs dell'elemento */

void DataManager::DofOwnerSet(void)
{
   DEBUGCOUTFNAME("DataManager::DofOwnerSet");
   /* Setta i DofOwner dei nodi */
   Node** ppTmpNode = ppNodes;
   for (; ppTmpNode < ppNodes+iTotNodes; ppTmpNode++) {            
      DofOwner* pDO = (DofOwner*)(*ppTmpNode)->pGetDofOwner();
      pDO->iNumDofs = (*ppTmpNode)->iGetNumDof();
   }
     
   /* Setta i DofOwner degli elementi (chi li possiede) */
   for (int iCnt = 0; iCnt < ElemType::LASTELEMTYPE; iCnt++) {
      DofType::Type DT = ElemData[iCnt].DofOwnerType;
      if (DT != DofType::UNKNOWN) {
	 DEBUGLCOUT(MYDEBUG_INIT, "Elem type " << iCnt 
		    << " (" << psElemNames[iCnt] << ")" << endl);
	 
	 Elem** ppFirstEl = ElemData[iCnt].ppFirstElem;
	 for (Elem** ppTmp = ppFirstEl; 
	      ppTmp < ppFirstEl+ElemData[iCnt].iNum; 
	      ppTmp++) {
	    ASSERT((*ppTmp)->pGetElemWithDofs() != NULL);
	    ElemWithDofs* pTmp = (*ppTmp)->pGetElemWithDofs();
	    
	    DEBUGLCOUT(MYDEBUG_INIT, "    " << psElemNames[pTmp->GetElemType()]
		       << "(" << pTmp->GetLabel() << ")" << endl);
	    
	    DofOwner* pDO = (DofOwner*)pTmp->pGetDofOwner();
	    pDO->iNumDofs = pTmp->iGetNumDof();
	    DEBUGLCOUT(MYDEBUG_INIT, "    num dofs: " << pDO->iNumDofs << endl);
	 }
      }
   }   
} /* end of DofOwnerSet() */



void DataManager::SetValue(VectorHandler& X, VectorHandler& XP)
{
   /* Nodi */
   for (Node** ppNode = ppNodes; ppNode < ppNodes+iTotNodes; ppNode++) {
      ASSERT(*ppNode != NULL);
      (*ppNode)->SetValue(X, XP);
   }

   /* Elementi */
#if !defined(USE_ELEM_ITER)
   Elem** ppLastElem = ppElems+iTotElem;
   for (Elem** ppEl = ppElems; ppEl < ppLastElem; ppEl++) {
      ASSERT(*ppEl != NULL);
      (*ppEl)->SetValue(X, XP);
   }        
#else /* USE_ELEM_ITER */
   /* Versione con iteratore: */
    Elem* pEl = NULL;
    if(ElemIter.fGetFirst(pEl)) {       
       do {
	  pEl->SetValue(X, XP);
       } while (ElemIter.fGetNext(pEl));
    }
#endif /* USE_ELEM_ITER */
} /* End of SetValue */


/* Output dati */
void DataManager::Output(void) const
{ 
   /* Dati dei nodi */
   NodeOutput((OutputHandler&)OutHdl);
   
   /* Dati degli elementi */
   ElemOutput((OutputHandler&)OutHdl);
   
   /* Nota: il casting di OutHdl e' necessario in quanto la funzione propria
    * <void DataManager::Output(void) const> e' dichiarata, appunto, <const>.
    * Questo fa si' che un oggetto proprio della classe DataManager sia
    * implicitamente definito come <const> agli occhi della funzione.
    * Dal momento che le funzioni 
    * <void NodeManager::Output(OutputHandler&) const> e
    * <void ElemManager::Output(OutputHandler&) const> ricevono come argomento
    * un oggetto di tipo <OutputHandler&> che non e' <const> in quanto su di
    * esso si scrive, il casting e' necessario per spiegare alla funzione
    * <void DataManager::Output(void) const> che le funzioni invocate
    * modificano si' l'<OutputHandler> passato loro, ma solo nel modo 
    * consentito e quindi la sua dichiarazione come funzione <const> e'
    * dovuta al fatto che i dati propri non vengono modificati in modo 
    * incontrollabile */   

   
   /* Restart condizionato */
   switch(RestartEvery) {
    case ITERATIONS: {
       if (++((integer&)iCurrRestartIter) == iRestartIterations) {
	  (integer&)iCurrRestartIter = 0;
	  ((DataManager*)this)->MakeRestart();
       }
       break;
    }
    case TIME: {
       ASSERT(pTime != NULL);
       if (pTime->GetVal().GetReal()-dLastRestartTime >= dRestartTime) {
	  (doublereal&)dLastRestartTime = pTime->GetVal().GetReal();
	  ((DataManager*)this)->MakeRestart();
       }	   
      break;
    }
    default:
      break;
   }
   
   /* Se richiesto, esegue l'output delle condizioni iniziali*/
   if (fAdamsOutput()) {
      ((integer&)iAdamsOutputBlock)++;
      AdamsResOutput(iAdamsOutputBlock, "DYNAMIC", "MBDyn");
   }
   
   
}


void DataManager::BeforePredict(VectorHandler& X, VectorHandler& XP,
				VectorHandler& XPrev, 
				VectorHandler& XPPrev) const
{
   Node** ppLastNode = ppNodes+iTotNodes;
   for (Node** ppTmp = ppNodes; ppTmp < ppLastNode; ppTmp++) {
      ASSERT(*ppTmp != NULL);
      (*ppTmp)->BeforePredict(X, XP, XPrev, XPPrev);
   }   

#if !defined(USE_ELEM_ITER)
   Elem** ppLastElem = ppElems+iTotElem;
   for (Elem** ppTmp = ppElems; ppTmp < ppLastElem; ppTmp++) {
      ASSERT(*ppTmp != NULL);
      (*ppTmp)->BeforePredict(X, XP, XPrev, XPPrev);
   }
#else /* USE_ELEM_ITER */
   /* Versione con iteratore: */
    Elem* pEl = NULL;
    if(((VecIter<Elem*>&)ElemIter).fGetFirst(pEl)) {       
       do {
	  pEl->BeforePredict(X, XP, XPrev, XPPrev);
       } while (((VecIter<Elem*>&)ElemIter).fGetNext(pEl));
    }
#endif /* USE_ELEM_ITER */
}


void DataManager::AfterPredict(void) const
{   
   Node** ppLastNode = ppNodes+iTotNodes;
   for (Node** ppTmp = ppNodes; ppTmp < ppLastNode; ppTmp++) {
      ASSERT(*ppTmp != NULL);
      (*ppTmp)->AfterPredict(*(VectorHandler*)pXCurr,
			     *(VectorHandler*)pXPrimeCurr);
   }

#if !defined(USE_ELEM_ITER)
   Elem** ppLastElem = ppElems+iTotElem;
   for (Elem** ppTmp = ppElems; ppTmp < ppLastElem; ppTmp++) {
      ASSERT(*ppTmp != NULL);
      (*ppTmp)->AfterPredict(*(VectorHandler*)pXCurr,
			     *(VectorHandler*)pXPrimeCurr);
   }   
#else /* USE_ELEM_ITER */
   /* Versione con iteratore: */
    Elem* pEl = NULL;
    if(((VecIter<Elem*>&)ElemIter).fGetFirst(pEl)) {       
       do {
	  pEl->AfterPredict(*(VectorHandler*)pXCurr,
			    *(VectorHandler*)pXPrimeCurr);
       } while (((VecIter<Elem*>&)ElemIter).fGetNext(pEl));
    }
#endif /* USE_ELEM_ITER */
}


void DataManager::Update(void) const
{   
   Node** ppLastNode = ppNodes+iTotNodes;
   for (Node** ppTmp = ppNodes; ppTmp < ppLastNode; ppTmp++) {
      ASSERT(*ppTmp != NULL);
      (*ppTmp)->Update(*pXCurr, *pXPrimeCurr);
   }   

#if !defined(USE_ELEM_ITER)
   Elem** ppLastElem = ppElems+iTotElem;
   for (Elem** ppTmp = ppElems; ppTmp < ppLastElem; ppTmp++) {
      ASSERT(*ppTmp != NULL);
      (*ppTmp)->Update(*pXCurr, *pXPrimeCurr);
   }   
#else /* USE_ELEM_ITER */
   /* Versione con iteratore: */
    Elem* pEl = NULL;
    if(((VecIter<Elem*>&)ElemIter).fGetFirst(pEl)) {       
       do {
	  pEl->Update(*pXCurr, *pXPrimeCurr);
       } while (((VecIter<Elem*>&)ElemIter).fGetNext(pEl));
    }
#endif /* USE_ELEM_ITER */
}


void DataManager::DerivativesUpdate(void) const
{  
   Node** ppLastNode = ppNodes+iTotNodes;
   for (Node** ppTmp = ppNodes; ppTmp < ppLastNode; ppTmp++) {
      ASSERT(*ppTmp != NULL);
      if ((*ppTmp)->GetNodeType() == NodeType::STRUCTURAL) {
	 (*(StructNode**)ppTmp)->DerivativesUpdate(*pXCurr, *pXPrimeCurr);
      } else {	 
	 (*ppTmp)->Update(*pXCurr, *pXPrimeCurr);
      }
   }
   
#if !defined(USE_ELEM_ITER)
   Elem** ppLastElem = ppElems+iTotElem;
   for (Elem** ppTmp = ppElems; ppTmp < ppLastElem; ppTmp++) {
      ASSERT(*ppTmp != NULL);
      (*ppTmp)->Update(*pXCurr, *pXPrimeCurr);
   }
#else /* USE_ELEM_ITER */
   /* Versione con iteratore: */
    Elem* pEl = NULL;
    if(((VecIter<Elem*>&)ElemIter).fGetFirst(pEl)) {       
       do {
	  pEl->Update(*pXCurr, *pXPrimeCurr);
       } while (((VecIter<Elem*>&)ElemIter).fGetNext(pEl));
    }
#endif /* USE_ELEM_ITER */
}

/* DataManager - end */
