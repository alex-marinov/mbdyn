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

/* DataManager - 
 * continua qui perche' il file dataman.cc sta diventando lungo */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <ac/math.h>

#include <dataman.h>
#include <dataman_.h>

#include <gravity.h>
#include <harwrap.h>
#include <mschwrap.h>
#include <y12wrap.h>
#include <umfpackwrap.h>

/* DataManager - continue */

/* Setta il valore della variabile Time nel DataManager 
 * usato dal metodo numerico all'inizio di ogni step temporale */

void DataManager::SetTime(doublereal dTime)
{   
   /* Setta la variabile Time nella tabella dei simboli */   
   ASSERT(pTime != NULL);
   pTime->SetVal(dTime);

   DEBUGLCOUT(MYDEBUG_INIT|MYDEBUG_ASSEMBLY, 
	      "GlobalSymbolTable:" << std::endl
	      << GlobalSymbolTable << std::endl);
   
   /* Setta il tempo nel DriveHandler */
   DrvHdl.SetTime(dTime);
   
   /* serve i drive pending */
   for (int iType = 0; iType < Drive::LASTDRIVETYPE; iType++) {
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
      
   /* per ogni nodo strutturale */
   if (NodeData[Node::STRUCTURAL].iNum > 0) {

      /*
       * used by POD stuff: if any, output the list of the first dof (minus 1)
       * of each structural node, so it's easy to get the struct node values
       * in MATLAB: given a vector "X" with all the states, and a vector
       * "v" with the first dof of each structural node, then the x coordinate
       * is X(v+1) and so forth
       */
      OutHdl.Log() << "struct node dofs:";
      
      StructNode** ppNd = (StructNode **)NodeData[Node::STRUCTURAL].ppFirstNode;
      for (unsigned long i = 0; i < NodeData[Node::STRUCTURAL].iNum; i++) {
         ASSERT(ppNd[i] != NULL);
	 if (ppNd[i]->GetStructNodeType() == StructNode::DUMMY) {
	    continue;
	 }
	 OutHdl.Log() << " " << ppNd[i]->iGetFirstIndex();
      }

      OutHdl.Log() << std::endl;
   }
   
   /* per ogni nodo */
   Node** ppNd = ppNodes;
   while (ppNd < ppNodes+iTotNodes) {	
      ASSERT(*ppNd != NULL);
      DEBUGLCOUT(MYDEBUG_INIT|MYDEBUG_ASSEMBLY, 
		 "Node type " << (*ppNd)->GetNodeType()
		 << " (" << psNodeNames[(*ppNd)->GetNodeType()] 
		 << "(" << (*ppNd)->GetLabel() << "))" << std::endl);
      
      unsigned int iNumDof;
      
      /* chiede al nodo quanti dof possiede */
      if ((iNumDof = (*ppNd)->iGetNumDof()) > 0) {	  	     
	 /* si fa passare il primo Dof */
	 Dof* pDf = pDofs+(*ppNd)->iGetFirstIndex();
	
#ifdef DEBUG
	 DEBUGLCOUT(MYDEBUG_INIT|MYDEBUG_ASSEMBLY,
		    psNodeNames[(*ppNd)->GetNodeType()] 
		    << "(" << (*ppNd)->GetLabel()
		    << "): first dof = " << pDf->iIndex+1 << std::endl);
#else /* !DEBUG */
	 if (!silent_output && fPrintDofStats) {
	    std::cout << psNodeNames[(*ppNd)->GetNodeType()]
	      << "(" << (*ppNd)->GetLabel()
	      << "): first dof = " << pDf->iIndex+1;
	    if ((*ppNd)->iGetNumDof() > 1) {
	       std::cout << "->" << pDf->iIndex+(*ppNd)->iGetNumDof();
	    }
	    std::cout << std::endl;
	 }
#endif /* !DEBUG */
	 
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
		 << "(" << (*ppEl)->GetLabel() << "))" << std::endl);
      
      unsigned int iNumDof;
	
      /* chiede all'elemento quanti dof possiede */
      if ((iNumDof = (*ppEl)->iGetNumDof()) > 0) {
	 ElemWithDofs* pElWD = (ElemWithDofs*)(*ppEl)->pGetElemWithDofs();
	 ASSERT(pElWD != NULL);
	 
	 /* si fa passare il DofOwner */
	 Dof* pDf = pDofs+pElWD->iGetFirstIndex();

#ifdef DEBUG
	 DEBUGLCOUT(MYDEBUG_INIT|MYDEBUG_ASSEMBLY, 
		    pdElemNames[pElWD->GetElemType()]
		    << "(" << pElWD->GetLabel()
		    << "): first dof = " << pDf->iIndex+1 << std::endl);
#else /* !DEBUG */
	 if (!silent_output && fPrintDofStats) {
	    std::cout << pdElemNames[pElWD->GetElemType()] 
	      << "(" << pElWD->GetLabel()
	      << "): first dof = " << pDf->iIndex+1 << std::endl
#endif /* !DEBUG */
	 
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
   if (ElemIter.fGetFirst(pEl)) {
      do {
	 ASSERT(pEl != NULL);
	 DEBUGLCOUT(MYDEBUG_INIT|MYDEBUG_ASSEMBLY,
		    "Elem type " << pEl->GetElemType() 
		    << " (" << psElemNames[pEl->GetElemType()] 
		    << "(" << pEl->GetLabel() << "))" << std::endl);
	 
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
	    
#ifdef DEBUG
	    DEBUGLCOUT(MYDEBUG_INIT|MYDEBUG_ASSEMBLY,
		       psElemNames[pEl->GetElemType()] 
		       << "(" << pElWD->GetLabel()
		       << "): first dof = " << pDf->iIndex+1 << std::endl);
#else /* !DEBUG */
	    if (!silent_output && fPrintDofStats) {
	       std::cout << psElemNames[pEl->GetElemType()]
		 << "(" << pElWD->GetLabel()
		 << "): first dof = " << pDf->iIndex+1;
	       if (pElWD->iGetNumDof() > 1) {
		  std::cout << "->" << pDf->iIndex+pElWD->iGetNumDof();
	       }
	       std::cout << std::endl;
	    }
#endif /* !DEBUG */
	    
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
   
   ASSERTMSG(DofData[DofOwner::JOINT].iNum > 0, 
	     "Warning, no joints are defined; You shouldn't have reached this point");
   ASSERT(DofData[DofOwner::STRUCTURALNODE].iNum > 0);

   /* Nodi strutturali: mette gli indici ai DofOwner */
   StructNode** ppFirstNode = 
     (StructNode**)NodeData[Node::STRUCTURAL].ppFirstNode;
   integer iNumNodes = NodeData[Node::STRUCTURAL].iNum;
   
   StructNode** ppNode = ppFirstNode;     
   DofOwner* pTmp = DofData[DofOwner::STRUCTURALNODE].pFirstDofOwner;
   
   integer iIndex = 0;    /* Indice dei gradi di liberta' */
   unsigned int iNumDofs = 0;  /* numero di dof di un owner */
   for (int iCnt = 1; 
        pTmp < DofData[DofOwner::STRUCTURALNODE].pFirstDofOwner+
        DofData[DofOwner::STRUCTURALNODE].iNum;
        iCnt++, pTmp++, ppNode++) {
      iNumDofs = pTmp->iNumDofs = (*ppNode)->iGetInitialNumDof();
      if (iNumDofs > 0) {
	 DEBUGLCOUT(MYDEBUG_INIT|MYDEBUG_ASSEMBLY,
		    psNodeNames[(*ppNode)->GetNodeType()] 
		    << "(" << (*ppNode)->GetLabel() 
		    << "): first index = " << iIndex+1 << std::endl);
	 pTmp->iFirstIndex = iIndex;
	 iIndex += iNumDofs;
      } else {
	 pedantic_cerr(psNodeNames[(*ppNode)->GetNodeType()] << "(" << iCnt 
			 << ") has 0 dofs" << std::endl);
      }
   }     

   
   /* Elementi: mette gli indici agli eventuali DofOwner */   
   for (int iCnt1 = 0; iCnt1 < Elem::LASTELEMTYPE; iCnt1++) {
      /* Pre ogni tipo di elemento */
      if (ElemData[iCnt1].fToBeUsedInAssembly && ElemData[iCnt1].iNum > 0) {
	 /* Se deve essere usato nell'assemblaggio e ne sono definiti */
	 
	 /* Tipo di dof dell'elemento corrente */
	 DofOwner::Type CurrDofType = 
	   ElemData[iCnt1].DofOwnerType;
	 
	 if (CurrDofType != DofOwner::UNKNOWN) {
	    
	    /* Puntatore al primo DofOwner */
	    pTmp = DofData[CurrDofType].pFirstDofOwner;
	 	 
	    Elem** ppFirstEl = ElemData[iCnt1].ppFirstElem;
	    integer iNumEls = ElemData[iCnt1].iNum;
	    ASSERT(DofData[CurrDofType].iNum == iNumEls);
	  
	    /* Iterazione sugli Elem */
	    Elem** ppEl = ppFirstEl;	     
	    for (int iCnt = 0;
		 pTmp < DofData[CurrDofType].pFirstDofOwner+iNumEls;
		 pTmp++, ppEl++, iCnt++) {
	       
	       ASSERT((*ppEl)->pGetInitialAssemblyElem() != NULL);
	       iNumDofs = 
		 (*ppEl)->pGetInitialAssemblyElem()->iGetInitialNumDof();
	       if ((pTmp->iNumDofs = iNumDofs) > 0) {
		  DEBUGLCOUT(MYDEBUG_INIT|MYDEBUG_ASSEMBLY,
			     psElemNames[iCnt1] << "(" << (*ppEl)->GetLabel()
			     << "): first index = " << iIndex+1 << std::endl);
		  pTmp->iFirstIndex = iIndex;
		  iIndex += iNumDofs;
	       } else {
		  pedantic_cerr(psElemNames[iCnt1] 
				  << "(" << (*ppEl)->GetLabel() 
				  << ") has 0 dofs" << std::endl);
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
   while (pEl != NULL) {
      integer iCurrRows = 0;
      integer iCurrCols = 0;
      pEl->InitialWorkSpaceDim(&iCurrRows, &iCurrCols);
      if (iCurrRows > iMaxRows) {	 
	 iMaxRows = iCurrRows;
      }      
      if (iCurrCols > iMaxCols) {	 
	 iMaxCols = iCurrCols;
      }      
      pEl = IAIter.GetNext();       
   }



   /*
    * Alla fine, i DofOwner di nodi e joint contengono gli indici giusti per 
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
    *   in base al tipo
    */

   
   /* Creazione e costruzione array Dof */   
   SAFENEWARR(pDofs, Dof, iInitialNumDofs);
   
   iIndex = 0;
   for (Dof* pTmpDof = pDofs; 
        pTmpDof < pDofs+iInitialNumDofs; pTmpDof++) {
      pTmpDof->iIndex = iIndex++;
   }


   /* Alloca spazio per l'assemblaggio */
   integer iIntDim = iMaxRows*iMaxCols*2;
   integer* piWI = NULL;
   SAFENEWARR(piWI, integer, iIntDim);
   
   integer iDoubleDim = iMaxRows*iMaxCols;
   doublereal* pdWD = NULL;
   SAFENEWARR(pdWD, doublereal, iDoubleDim);   
   
   
   /* Ciclo di iterazioni fino a convergenza */
   
   /* Crea la struttura di calcolo */
   SolutionManager* pSM = NULL;
#if defined(USE_Y12)
   SAFENEWWITHCONSTRUCTOR(pSM,
		   Y12SparseLUSolutionManager,
		   Y12SparseLUSolutionManager(iInitialNumDofs, 0, 1.));
#elif defined(USE_UMFPACK)
   SAFENEWWITHCONSTRUCTOR(pSM,
		   UmfpackSparseLUSolutionManager,
		   UmfpackSparseLUSolutionManager(iInitialNumDofs, 0, 1.));
#elif defined(USE_HARWELL)
   SAFENEWWITHCONSTRUCTOR(pSM,
		   HarwellSparseLUSolutionManager,
		   HarwellSparseLUSolutionManager(iInitialNumDofs, 0, 1.));
#elif defined(USE_MESCHACH)
   SAFENEWWITHCONSTRUCTOR(pSM,
		   MeschachSparseLUSolutionManager,
		   MeschachSparseLUSolutionManager(iInitialNumDofs, 0, 1.));
#else
#error "no solver available"
#endif
    
   /* Crea il vettore con lo stato del sistema durante l'assemblaggio */
   doublereal* pdX = NULL;
   SAFENEWARR(pdX, doublereal, iInitialNumDofs);
   
#ifdef DEBUG_MEMMANAGER
   DEBUGLCOUT(MYDEBUG_MEM|MYDEBUG_ASSEMBLY, 
	      "After initialisation in InitialJointAssembly" << std::endl
	      << defaultMemoryManager << std::endl);
#endif /* DEBUG_MEMMANAGER */
   
   MyVectorHandler X(iInitialNumDofs, pdX);
   X.Reset(0.);
   
   /* Linka il DriveHandler al vettore soluzione */
   LinkToSolution(X, X);

   /* Setta i valori iniziali dei gradi di liberta' dei nodi strutturali
    * durante l'assemblaggio iniziale */
   for (StructNode** ppTmpNode = ppFirstNode;
        ppTmpNode < ppFirstNode+iNumNodes; ppTmpNode++) {	
      (*ppTmpNode)->SetInitialValue(X);
   }   
     
   
   /* Setta i valori iniziali dei gradi di liberta' dei vincoli
    * durante l'assemblaggio iniziale */
   for (int iCnt1 = 0; iCnt1 < Elem::LASTELEMTYPE; iCnt1++) {
      /* Pre ogni tipo di elemento */
      if (ElemData[iCnt1].DofOwnerType != DofOwner::UNKNOWN &&
	  ElemData[iCnt1].fToBeUsedInAssembly &&
	  ElemData[iCnt1].iNum > 0) {
	 
	 Elem** ppFirstEl = ElemData[iCnt1].ppFirstElem;
	 integer iNumEl = ElemData[iCnt1].iNum;
	 for (Elem** ppTmpEl = ppFirstEl;
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
   while (++iNumIter) {	
      /* Assemblo il residuo */
      pResHdl->Reset(0.);
      
      /* Contributo dei nodi */	
      for (StructNode** ppTmpNode = ppFirstNode; 
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
	 
	 
	 if ((*ppTmpNode)->fOmegaRotates()) {
	    TmpVec = RDelta*wPrev-wCurr; // con questa la velocita' angolare e' solidale con il nodo
	 } else {	 
	    TmpVec = wPrev-wCurr; // con questa la velocita' angolare e' solidale col riferimento assoluto
	 }
	 
	 pResHdl->Add(iFirstIndex+10, TmpVec*dVelStiff);
	 
      }	
       	
      /* Elementi (con iteratore): */
      pEl = IAIter.GetFirst();
      while (pEl != NULL) {
	 *pResHdl += pEl->InitialAssRes(WorkVec, X);
	 pEl = IAIter.GetNext();
      }	
      
	
#ifdef DEBUG
	if (DEBUG_LEVEL_MATCH(MYDEBUG_ASSEMBLY|MYDEBUG_RESIDUAL)) {
	   /* Output del residuo */
	   std::cout << "Residual:" << std::endl;
	   for (int iTmpCnt = 1; iTmpCnt <= iInitialNumDofs; iTmpCnt++) {
	      std::cout << "Dof" << std::setw(4) << iTmpCnt << ": "
		<< pResHdl->dGetCoef(iTmpCnt) << std::endl;
	   }
	}
#endif
	
	
      /* Eseguo il test di convergenza; se e' positivo, esco */
      doublereal dTest = pResHdl->Dot()/(1.+X.Dot());
      if (!isfinite(dTest)) {
	 std::cerr << "Assembly diverged; aborting ..." << std::endl;

	 THROW(DataManager::ErrAssemblyDiverged());
      }
      dTest = sqrt(dTest);
      
      DEBUGLCOUT(MYDEBUG_ASSEMBLY, "Iteration: " << iNumIter 
		 << ", Test: " << dTest 
		 << " (Toll = " << dInitialAssemblyToll << ") " << std::endl);
	
      /* Se la tolleranza e' raggiunta, esce dal ciclo */
      if (dTest <= dInitialAssemblyToll) {
	 DEBUGLCOUT(MYDEBUG_ASSEMBLY, "Initial assembly performed successfully in "
		    << iNumIter << " iterations" << std::endl);
	 goto endofcycle;
      }
	
	
	
      /* Se ho raggiunto il numero massimo di iterazioni */
      if (iNumIter > iMaxInitialIterations) {
	 std::cerr
	   << "Initial assembly iterations reached maximum number "
	   << iMaxInitialIterations << "; aborting ..." << std::endl;

	 THROW(DataManager::ErrAssemblyMaxIters());
      }
      
	
      /* Assemblo lo jacobiano e risolvo */
      pSM->MatrInit(0.);
	
      /* Contributo dei nodi */	
      for (StructNode** ppTmpNode = ppFirstNode; 
	  ppTmpNode < ppFirstNode+iNumNodes; ppTmpNode++) {
	 integer iFirstIndex = ((*ppTmpNode)->pGetDofOwner())->iFirstIndex;
	 
	 /* Nuova feature: ogni nodo ha la sua stiffness */
	 doublereal dPosStiff = (*ppTmpNode)->dGetPositionStiffness();
	 doublereal dVelStiff = (*ppTmpNode)->dGetVelocityStiffness();
	 
	 for (int iCnt = 1; iCnt <= 6; iCnt++) {
	    /* Posizione, rotazione */
	    integer iTmp = iFirstIndex+iCnt;
	    pMatHdl->fPutCoef(iTmp, iTmp, dPosStiff);
	    
	    /* Velocita', velocita' angolare */
	    iTmp += 6;
	    pMatHdl->fPutCoef(iTmp, iTmp, dVelStiff);
	 }
	 

	 if ((*ppTmpNode)->fOmegaRotates()) {
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
      while (pEl != NULL) {
	 *pMatHdl += pEl->InitialAssJac(WorkMat, X);
	 pEl = IAIter.GetNext();
      }
      


	
      /* Fattorizza e risolve con jacobiano e residuo appena calcolati */
      pSM->Solve(0.);
      
#ifdef DEBUG
      if (DEBUG_LEVEL_MATCH(MYDEBUG_ASSEMBLY|MYDEBUG_RESIDUAL)) {
	 /* Output della soluzione */
	 std::cout << "Solution:" << std::endl;
	 for (int iTmpCnt = 1; iTmpCnt <= iInitialNumDofs; iTmpCnt++) {	   
	    std::cout << "Dof" << std::setw(4) << iTmpCnt << ": " 
	      << pSolHdl->dGetCoef(iTmpCnt) << std::endl;
	 }
      }
#endif
	
	
      /* Aggiorno la soluzione */
      X += *pSolHdl;
	
      /* Correggo i nodi */
      for (StructNode** ppTmpNode = ppFirstNode;
	  ppTmpNode < ppFirstNode+iNumNodes; ppTmpNode++) {
	 (*ppTmpNode)->InitialUpdate(X);
      }
   }       
   
   
   
endofcycle:   
      
   /* Distrugge il workspace */
   ASSERT(pdWD != NULL);
   if (pdWD != NULL) {	
      SAFEDELETEARR(pdWD);
   }
   
   ASSERT(piWI != NULL);
   if (piWI != NULL) {
      SAFEDELETEARR(piWI);
   }
      
   
   /* Distrugge il vettore soluzione */
   ASSERT(pdX != NULL);
   if (pdX != NULL) {	
      SAFEDELETEARR(pdX);
   }   

   
   /* Resetta e distrugge la struttura temporanea dei Dof */

   /* Elementi: rimette a posto il numero di Dof propri dei vincoli */
   for (int iCnt1 = 0; iCnt1 < Elem::LASTELEMTYPE; iCnt1++) {
      /* Per ogni tipo di elemento */
      if (ElemData[iCnt1].DofOwnerType != DofOwner::UNKNOWN &&
	  ElemData[iCnt1].fToBeUsedInAssembly &&
	  ElemData[iCnt1].iNum > 0) {
	 /* Se possiede dofs, se deve essere usato nell'assemblaggio
	  * e se ne sono presenti */
	 
	 /* Tipo di dof dell'elemento corrente */
	 DofOwner::Type CurrDofType = 
	   ElemData[iCnt1].DofOwnerType;
	 /* Puntatore al primo DofOwner */
	 pTmp = DofData[CurrDofType].pFirstDofOwner;
	 
	 
	 /* Puntatore al primo Elem */
	 Elem** ppFirstEl = ElemData[iCnt1].ppFirstElem;
	 /* Numero di Elem (== al numero di DofOwner) */
	 integer iNumEls = DofData[CurrDofType].iNum;
	 
	 /* Iterazione sugli Elem */
	 Elem** ppEl = ppFirstEl;
	 for (; pTmp < DofData[CurrDofType].pFirstDofOwner+iNumEls;
	      pTmp++, ppEl++) {
	    pTmp->iNumDofs = (*ppEl)->iGetNumDof();
	 }
      }
   }              
   
   
   /* Dealloca il vettore dei Dof */
   ASSERT(pDofs != NULL);
   if (pDofs != NULL) {	
      SAFEDELETEARR((Dof*&)pDofs);
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
   for (int iCnt = 0; iCnt < Elem::LASTELEMTYPE; iCnt++) {
      DofOwner::Type DT = ElemData[iCnt].DofOwnerType;
      if (DT != DofOwner::UNKNOWN) {
	 DEBUGLCOUT(MYDEBUG_INIT, "Elem type " << iCnt 
		    << " (" << psElemNames[iCnt] << ")" << std::endl);
	 
	 Elem** ppFirstEl = ElemData[iCnt].ppFirstElem;
	 for (Elem** ppTmp = ppFirstEl; 
	      ppTmp < ppFirstEl+ElemData[iCnt].iNum; 
	      ppTmp++) {
	    ASSERT((*ppTmp)->pGetElemWithDofs() != NULL);
	    ElemWithDofs* pTmp = (*ppTmp)->pGetElemWithDofs();
	    
	    DEBUGLCOUT(MYDEBUG_INIT, "    " << psElemNames[pTmp->GetElemType()]
		       << "(" << pTmp->GetLabel() << ")" << std::endl);
	    
	    DofOwner* pDO = (DofOwner*)pTmp->pGetDofOwner();
	    pDO->iNumDofs = pTmp->iGetNumDof();
	    DEBUGLCOUT(MYDEBUG_INIT, "    num dofs: " << pDO->iNumDofs << std::endl);
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
    if (ElemIter.fGetFirst(pEl)) {       
       do {
	  pEl->SetValue(X, XP);
       } while (ElemIter.fGetNext(pEl));
    }
#endif /* USE_ELEM_ITER */
} /* End of SetValue */


/* Output dati */
void DataManager::Output(bool force) const
{ 
   /* output only at multiples of iOutputFrequency */
   if ((!force) && (iOutputCount++%iOutputFrequency)) {
      return;
   }

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


/* Output dati */
void
DataManager::Output(const VectorHandler& X, const VectorHandler& XP) const
{
   /* Dati dei nodi */
   NodeOutput((OutputHandler&)OutHdl, X, XP);
   
   /* Dati degli elementi */
   ElemOutput((OutputHandler&)OutHdl, X, XP);
};


/* Output dati pch */
void
DataManager::Output_pch(std::ostream& pch) const
{
   /* Dati dei nodi */
   NodeOutput_pch(pch);
   
   /* Dati degli elementi */
   ElemOutput_pch(pch);
};


/* Output dati f06 */
void
DataManager::Output_f06(std::ostream& f06, const VectorHandler& X) const
{
   /* Dati dei nodi */
   NodeOutput_f06(f06, X);
   
   /* Dati degli elementi */
   ElemOutput_f06(f06, X);
};


/* Output dati f06 */
void
DataManager::Output_f06(std::ostream& f06, const VectorHandler& Xr, const VectorHandler& Xi) const
{
   /* Dati dei nodi */
   NodeOutput_f06(f06, Xr, Xi);
   
   /* Dati degli elementi */
   ElemOutput_f06(f06, Xr, Xi);
};


/* OpenDX output (tbi) */
void DataManager::Output_OpenDX(std::ostream& dx, const VectorHandler& Xr, const VectorHandler& Xi) const
{
	NO_OP;
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
    if (ElemIter.fGetFirst(pEl)) {       
       do {
	  pEl->BeforePredict(X, XP, XPrev, XPPrev);
       } while (ElemIter.fGetNext(pEl));
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
    if (ElemIter.fGetFirst(pEl)) {       
       do {
	  pEl->AfterPredict(*(VectorHandler*)pXCurr,
			    *(VectorHandler*)pXPrimeCurr);
       } while (ElemIter.fGetNext(pEl));
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
    if (ElemIter.fGetFirst(pEl)) {       
       do {
	  pEl->Update(*pXCurr, *pXPrimeCurr);
       } while (ElemIter.fGetNext(pEl));
    }
#endif /* USE_ELEM_ITER */
}


void DataManager::AfterConvergence(void) const
{   
   Node** ppLastNode = ppNodes+iTotNodes;
   for (Node** ppTmp = ppNodes; ppTmp < ppLastNode; ppTmp++) {
      ASSERT(*ppTmp != NULL);
      (*ppTmp)->AfterConvergence(*(VectorHandler*)pXCurr,
			     *(VectorHandler*)pXPrimeCurr);
   }

#if !defined(USE_ELEM_ITER)
   Elem** ppLastElem = ppElems+iTotElem;
   for (Elem** ppTmp = ppElems; ppTmp < ppLastElem; ppTmp++) {
      ASSERT(*ppTmp != NULL);
      (*ppTmp)->AfterConvergence(*(VectorHandler*)pXCurr,
			     *(VectorHandler*)pXPrimeCurr);
   }   
#else /* USE_ELEM_ITER */
   /* Versione con iteratore: */
    Elem* pEl = NULL;
    if (ElemIter.fGetFirst(pEl)) {       
       do {
	  pEl->AfterConvergence(*(VectorHandler*)pXCurr,
			    *(VectorHandler*)pXPrimeCurr);
       } while (ElemIter.fGetNext(pEl));
    }
#endif /* USE_ELEM_ITER */
}


void DataManager::DerivativesUpdate(void) const
{  
   Node** ppLastNode = ppNodes+iTotNodes;
   for (Node** ppTmp = ppNodes; ppTmp < ppLastNode; ppTmp++) {
      ASSERT(*ppTmp != NULL);
      if ((*ppTmp)->GetNodeType() == Node::STRUCTURAL) {
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
    if (ElemIter.fGetFirst(pEl)) {       
       do {
	  pEl->Update(*pXCurr, *pXPrimeCurr);
       } while (ElemIter.fGetNext(pEl));
    }
#endif /* USE_ELEM_ITER */
}

/* DataManager - end */
