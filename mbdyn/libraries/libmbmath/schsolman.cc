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

/*  SchurMatrixHandler 
    SchurVectorHandler
    SchurSolutionManager */
    
/* 
 * Copyright 1999-2000 Giuseppe Quaranta <quaranta@aero.polimi.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#ifdef USE_MPI

#include <schsolman.h>
#ifdef USE_MYSLEEP
#include <mysleep.h>
#endif /* USE_MYSLEEP */
/* 
 * Le occorrenze della routine 'mysleep' vanno commentate quando 
 * si opera su di una macchina dove le comunicazioni sono efficienti
 * e non ci sono altri utenti operanti che esauriscono le risorse.
 * In caso contrario la si utilizza per ottenere delle misure significative
 * di tempi di CPU, anche se i tempi solari ottenibili non sono proporzionati 
 */ 

#ifdef MPI_PROFILING 
extern "C" {
#include <mpe.h>
#include <stdio.h>
};
#endif /* MPI_PROFILING */

#undef min
#undef max
#include <algorithm>

SchurSolutionManager::SchurSolutionManager (integer iSize, 
					    integer* pLocalDofs,
					    int iDim1,
					    integer* pInterfDofs,
					    int iDim2, 
					    integer iWorkSize,
					    const doublereal& dPivotFactor)
: iPrbmSize(iSize),
SolvCommSize(0),
pLocDofs(pLocalDofs),
LocVecDim(iDim1),
pIntDofs(pInterfDofs),
IntVecDim(iDim2),
pDispl(NULL),
pRecvDim(NULL),
pDofsRecvdList(NULL),
pSchurDofs(NULL),
iSchurIntDim(0),
pGlbToLoc(NULL),
pSchGlbToLoc(NULL),
pBlockLenght(NULL), 
pTypeDsp(NULL), 
pBuffer(NULL),
ppNewTypes(NULL),
pBMH(NULL), 
pdEMH(NULL), 
pdFMH(NULL), 
pdCMH(NULL), 
prVH(NULL),
pgVH(NULL),
pdgVec(NULL), 
pMH(NULL), 
pRVH(NULL),
pSchMH(NULL), 
pSchVH(NULL),
pGSReq(NULL), 
pGRReq(NULL),  
iWorkSpaceSize(iWorkSize), 
fNewMatrix(0) 
{
  DEBUGCOUT("Entering SchurSolutionManager::SchurSolutionManager()"
		  << endl);
  
  ASSERT(iPrbmSize > 0);
  ASSERT(pLocDofs != NULL);
  
  /* inizializzo il communicator */
  SolvComm = MPI::COMM_WORLD.Dup();
  SolvCommSize = SolvComm.Get_size();
  MyRank = SolvComm.Get_rank();
  
#ifdef DEBUG
  if (SolvCommSize != 1) {
    ASSERT(pIntDofs != NULL);
  }
#endif /* DEBUG */
  
  DEBUGCOUT("Solution Communicator Size: " << SolvCommSize << endl);
  
  /* detrmino le dimensioni locali e il vettore di trasferimento Globale -> Locale;
     sul master (myrank == 0) la matrice C è quella di tutte le interfaccie */
  
  InitializeComm();


//   ASSERT((dPivotFactor >= 0.0) && (dPivotFactor <= 1.0));
// 
// 
//   /* Valore di default */
//   if (iWorkSpaceSize != 0) {
//   	iWorkSpaceSize = static_cast<integer>(sqrt(iWorkSpaceSize) + 1);
//   } else {	
//         iWorkSpaceSize = LocVecDim;
//   }

  /* matrice C e vettore g*/
  SAFENEWARR(pCMH, doublereal, IntVecDim*IntVecDim);
  SAFENEWARR(pdgVec, doublereal, IntVecDim);
  
  /* Costruisce il solutore UMFPACK3 locale */
  DEBUGCOUT(" Using UMFPACK3 math library as local solver" << endl);  
  SAFENEWWITHCONSTRUCTOR(pLocalSM, Umfpack3SparseLUSolutionManager,
			   Umfpack3SparseLUSolutionManager(LocVecDim));			    
   pBMH = pLocalSM->pMatHdl();
   prVH = pLocalSM->pResHdl();

  if(!MyRank) {
     	DEBUGCOUT(" Using UMFPACK3 math library as interface solver" << endl);  
  	/* Costruisce il solutore UMFPACK3 di interfaccia */
  	SAFENEWWITHCONSTRUCTOR(pInterfSM, Umfpack3SparseLUSolutionManager,
				   Umfpack3SparseLUSolutionManager(iSchurIntDim));			    
   	pSchMH = pLocalSM->pMatHdl();
   	pSchVH = pLocalSM->pResHdl();
   }	

  /* Definisco le matrici globali */
   SAFENEWWITHCONSTRUCTOR(pMH, SchurMatrixHandler, 
			 SchurMatrixHandler(LocVecDim, IntVecDim, pBMH, 
					    pdCMat, pGlbToLoc));
  
  pEMH = pMH->GetEMat();
  pFMH = pMH->GetFMat();
  
  SAFENEWWITHCONSTRUCTOR(pRVH, SchurVectorHandler, 
			 SchurVectorHandler(LocVecDim, IntVecDim, prVH,
					    pdgVec, pGlbToLoc));
  
  pgVH = pRVH->GetI();
     
  /* Creazione di NewType per la trasmissione del vett soluzione calcolato */
  if(!MyRank) {
    SAFENEWARR(pBlockLenght, int, pDispl[SolvCommSize]);
    InitializeList(pBlockLenght, pDispl[SolvCommSize], 1);
    SAFENEWARR(pTypeDsp, MPI::Aint, pDispl[SolvCommSize]);
    
    MPI::Aint DispTmp = MPI::Get_address(pSchVH->pdGetVec());
    for (int i=0; i < pDispl[SolvCommSize]; i++) {
      pTypeDsp[i] = MPI::Get_address(pSchVH->pdGetVec() + pSchGlbToLoc[pDofsRecvdList[i]] - 1) - DispTmp;
    }
    
    SAFENEWARR(ppNewTypes, MPI::Datatype*, SolvCommSize);
    MPI::Aint* pActualDispl = pTypeDsp;
    for (int i=0; i < SolvCommSize; i++) {
      ppNewTypes[i] = NULL;
      SAFENEWWITHCONSTRUCTOR(ppNewTypes[i], MPI::Datatype, 
		      MPI::Datatype(MPI::DOUBLE.Create_hindexed(pRecvDim[i], 
				      pBlockLenght, pActualDispl)));
      ppNewTypes[i]->Commit();
      pActualDispl += pRecvDim[i];
    }
  }
  
  cout << "Local dimension on process " << MyRank << ": " << LocVecDim << endl;
  if(!MyRank) {
    cout << "Interface dimension: " << iSchurIntDim << endl;
  }
  
#ifdef DEBUG
  if(!MyRank) {
    cout << "Interface Dofs " <<endl;
    for (int i=0; i < iSchurIntDim; i++) { 
      cout << pSchurDofs[i] << "  ";
    }
    cout << endl;
  }
  IsValid();
#endif /* DEBUG */
}

SchurSolutionManager::~SchurSolutionManager(void)
{
  DEBUGCOUT("Entering SchurSolutionManager::~SchurSolutionManager()" << endl);

#ifdef DEBUG
  IsValid();
#endif /* DEBUG */

#if 0
  if(pLocDofs != NULL) {
    SAFEDELETEARR(pLocDofs);
  }
  if(pIntDofs != NULL) {
    SAFEDELETEARR(pIntDofs);
  } 
  if(pRecvDim != NULL) {
    SAFEDELETEARR(pRecvDim);
  }
  if(pDispl != NULL) {
    SAFEDELETEARR(pDispl);
  }
  if(pDofsRecvdList != NULL) {
    SAFEDELETEARR(pDofsRecvdList);
  }
  if(pSchurDofs != NULL) {
    SAFEDELETEARR(pSchurDofs);
  }
  if(pGlbToLoc != NULL) {
    SAFEDELETEARR(pGlbToLoc);
  }
  if(pSchGlbToLoc != NULL) {
    SAFEDELETEARR(pSchGlbToLoc);
  }
  if(pBlockLenght != NULL) {
    SAFEDELETEARR(pBlockLenght);
  }
   if(pTypeDsp != NULL) {
    SAFEDELETEARR(pTypeDsp);
  }
  if(pBuffer != NULL) {
    SAFEDELETEARR(pBuffer);
  }
  if(ppNewTypes != NULL) {
    for(int i=0; i < SolvCommSize; i++) {
       ppNewTypes[i]->Free();
      SAFEDELETE(ppNewTypes[i]);
    }
   SAFEDELETEARR(ppNewTypes);
  }
  if(pBMH != NULL) {
      SAFEDELETE(pBMH);
  }
  if(pdCMat != NULL) {
    SAFEDELETEARR(pdCMat);
  }
  if(pdgVec != NULL) {
    SAFEDELETEARR(pdgVec);
  }
  if(pMH != NULL) {
	SAFEDELETE(pMH);
  }
  if(pRVH != NULL) {
	SAFEDELETE(pRVH);
  }
  if(pSchSM != NULL) {
      SAFEDELETEARR(pSchSM);
  }
  if(pSchMH != NULL) {
    SAFEDELETEARR(pSchMH);
  }
  if(pSchVH != NULL) {
    SAFEDELETEARR(pSchVH);
  }
#endif /* 0 */
}


void SchurSolutionManager::IsValid(void) const
{
  NO_OP;
}

/* Inizializza il gestore delle matrici */

void SchurSolutionManager::MatrInit(const doublereal& dResetVal)
{
#ifdef DEBUG
  IsValid();
#endif /* DEBUG */
  pMH->Init(dResetVal);
  fNewMatrix = flag(1);
}

/* Risolve i blocchi */

void SchurSolutionManager::Solve(void)
{
  DEBUGCOUT("Entering SchurSolutionManager::Solve()" << endl);
#ifdef DEBUG
  IsValid();
#endif /* DEBUG */

#ifdef MPI_PROFILING 
  MPE_Log_event(31, 0, "start");
#endif /* MPI_PROFILING */   
  
  if(fNewMatrix) {
     pLocalSM->Solve();
  } 
  
#ifdef MPI_PROFILING 
  MPE_Log_event(32, 0, "end");
#endif /* MPI_PROFILING */   
     
  /* prodotto F * r */
  pFMH->MatVecDecMul(pgVH, prVH);
  
  /* Inizializza Trasmissione di g */
#ifdef MPI_PROFILING
   MPE_Log_event(13, 0, "start");
   MPE_Log_send(0, G_TAG, IntVecDim);
#endif /* MPI_PROFILING */   
   pGSReq[0] = SolvComm.Isend(pdgVec, IntVecDim, MPI::DOUBLE, 0, G_TAG);
  
   if(!MyRank) {
#ifdef MPI_PROFILING
     MPE_Log_event(19, 0, "start");
#endif /* MPI_PROFILING */   
     for( int i=0; i < SolvCommSize; i++){
       pGRReq[i] = SolvComm.Irecv(pBuffer + pDispl[i], pRecvDim[i], MPI::DOUBLE, i, G_TAG);
       
     }
   }
  
   if (fNewMatrix) {
     /* Calcolo di E' */ 
     for (int i=0; i < IntVecDim; i++) {
       pEMH->GetCol(i, pLocalSM->pdRhs);
       pLocalSM->BackSub();
     }
   }
   
   /* verifica completamento ricezioni e trasmissione g*/
   flag SentFlag;
   while (1) {
     SentFlag = pGSReq->Test();
     if (SentFlag) {
#ifdef MPI_PROFILING
      MPE_Log_event(14, 0, "end");
#endif /* MPI_PROFILING */       
      break;
#ifdef USE_MYSLEEP
     } else {
       mysleep(150);
#endif /* USE_MYSLEEP */
     }
   }
   
  if(!MyRank) {
    flag RecvFlag;
    while (1) {
      RecvFlag = MPI::Request::Testall(SolvCommSize, pGRReq);
      if (RecvFlag) {
#ifdef MPI_PROFILING
        MPE_Log_event(20, 0, "end");
        for( int i=0; i < SolvCommSize; i++){
          MPE_Log_receive(i, G_TAG, pRecvDim[i]);
        }
#endif /* MPI_PROFILING */
  	break;
#ifdef USE_MYSLEEP
      } else {
	mysleep(1500);
#endif /* USE_MYSLEEP */
      }
    }
  }
  
  /* Assemblaggio del vettore dell incognite g sulla macchina master */  
  if (!MyRank) {
    /* Assembla pSVH */
    for (int i=0; i < pDispl[SolvCommSize]; i++) {
      pSchVH->fIncCoef(pSchGlbToLoc[pDofsRecvdList[i]], pBuffer[i]);
    }
  }    
  
  
  /* assembla la matrice di schur se è stata appena assemblata una nuova matrice di partenza */
  
  if (fNewMatrix) {
    AssSchur();
    fNewMatrix = flag(0);
  }
  
  if ((!MyRank) && (SolvCommSize != 1)) {
#ifdef MPI_PROFILING 
    MPE_Log_event(35, 0, "start");
#endif /* MPI_PROFILING */   

  /* risoluzione schur pb sul master */
    pSchSM->Solve();
#ifdef MPI_PROFILING 
    MPE_Log_event(36, 0, "end");
#endif /* MPI_PROFILING */       
  }
  
  /* invia la soluzione calcolata */  
  if(!MyRank) {
#ifdef MPI_PROFILING
    MPE_Log_event(13, 0, "start");
#endif /* MPI_PROFILING */ 
    for( int i=0; i < SolvCommSize; i++){
      pGRReq[i] = SolvComm.Isend(pSchVH->pdGetVec(), 1, *ppNewTypes[i], i, G_TAG);	 
#ifdef MPI_PROFILING
      MPE_Log_send(i, G_TAG, pRecvDim[i]);  
#endif /* MPI_PROFILING */     
    }	
  }
    
#ifdef MPI_PROFILING
   MPE_Log_event(19, 0, "start");
#endif /* MPI_PROFILING */   

  pGSReq[0] = SolvComm.Irecv(pdgVec, IntVecDim, MPI::DOUBLE, 0, G_TAG);
  
  /* Verifica completamento trasmissioni */
  if(!MyRank) {
    flag SentFlag;
    while (1) {
      SentFlag = MPI::Request::Testall(SolvCommSize, pGRReq);
      if (SentFlag) {
#ifdef MPI_PROFILING
        MPE_Log_event(14, 0, "end");
#endif /* MPI_PROFILING */
	break;
#ifdef USE_MYSLEEP
      } else {
	mysleep(150);
#endif /* USE_MYSLEEP */
      }
    }
  }
  
  flag RecvFlag;
  while (1) {
    RecvFlag = pGSReq->Test();
    if (RecvFlag) {
#ifdef MPI_PROFILING
    MPE_Log_event(20, 0, "end");
    MPE_Log_receive(0, G_TAG, IntVecDim);
#endif /* MPI_PROFILING */ 
      break;
#ifdef USE_MYSLEEP
    } else {
      mysleep(1500);
#endif /* USE_MYSLEEP */
    }
  }
  
   
  /* calcolo la soluzione corretta per x */
  pEMH->MatVecDecMul(prVH, pgVH); 

}



void SchurSolutionManager::AssSchur(void) 
{
#ifdef MPI_PROFILING
   MPE_Log_event(41, 0, "start");
#endif /* MPI_PROFILING */ 
  
  DEBUGCOUT("Entering SchurSolutionManager::AssSchur()" << endl);
  
#ifdef DEBUG
  IsValid();
#endif /* DEBUG */
  SpMapMatrixHandler Mtmp(IntVecDim, IntVecDim);  
  pFMH->MatMatMul(&Mtmp,pEMH);
  
 /*******************************CONTROLLA**************/ 
  
  /* Calcola le Schur locali */
  for(int j=0; j < IntVecDim; j++) {
    int iColc = j * IntVecDim;
    for( int i=0; i < IntVecDim; i++) {
     pdCMat[i + iColc] -=  Mtmp(i,j);
      }
    }
  }
 
  /* Trasmette le Schur locali */
#ifdef MPI_PROFILING
   MPE_Log_event(13, 0, "start");
   MPE_Log_send(0, S_TAG, IntVecDim*IntVecDim);
#endif /* MPI_PROFILING */   
  
  pGSReq[0] = SolvComm.Isend(pdCMat, IntVecDim*IntVecDim, MPI::DOUBLE, 0, S_TAG);
  int iOffset = 0;
  if(!MyRank) { 
  #ifdef MPI_PROFILING
   MPE_Log_event(19, 0, "start");
#endif /* MPI_PROFILING */  
    for( int i=0; i < SolvCommSize; i++){
      pGRReq[i] = SolvComm.Irecv(pBuffer + iOffset, (pRecvDim[i]*pRecvDim[i]), MPI::DOUBLE, i, S_TAG);
      iOffset += (pRecvDim[i]*pRecvDim[i]);
    }
  }

  /* verifica completamento ricezioni e trasmissione */
  flag SentFlag;
  while (1) {
    SentFlag = pGSReq->Test();
    if (SentFlag) {
#ifdef MPI_PROFILING
      MPE_Log_event(14, 0, "end");
#endif /* MPI_PROFILING */        
      break;
#ifdef USE_MYSLEEP
    } else {
      mysleep(150);
#endif /* USE_MYSLEEP */
    }
  }
  
  
  if(!MyRank) {
    flag RecvFlag;
    while (1) {
      RecvFlag = MPI::Request::Testall(SolvCommSize, pGRReq);
      if (RecvFlag) {
#ifdef MPI_PROFILING
        MPE_Log_event(20, 0, "end");
        for( int i=0; i < SolvCommSize; i++){
          MPE_Log_receive(i, S_TAG, pRecvDim[i]*pRecvDim[i]);
        }
#endif /* MPI_PROFILING */      
	break;
#ifdef USE_MYSLEEP
      } else {
	mysleep(1500);
#endif /* USE_MYSLEEP */
      }
    }
    
    
    /* Assembla Schur matrix */

    
    pSchSM->MatrInit();
    iOffset = 0;
    for (int i=0; i < SolvCommSize; i++) {
      for(int  j=0; j < pRecvDim[i]; j++) {
	int iColx = j * pRecvDim[i];
	for( int k=0; k < pRecvDim[i]; k++) {
	  pSchMH->fIncCoef(pSchGlbToLoc[pDofsRecvdList[pDispl[i]+k]], pSchGlbToLoc[pDofsRecvdList[pDispl[i]+j]], pBuffer[iOffset + k + iColx]);
	}
      }
      iOffset += (pRecvDim[i]*pRecvDim[i]);
    }
  }
#ifdef MPI_PROFILING
  MPE_Log_event(42, 0, "end");
#endif /* MPI_PROFILING */ 
}



void SchurSolutionManager::InitializeComm(void)
{
  DEBUGCOUT("Entering SchurSolutionManager::InitializeComm()" << endl);

  if (!MyRank) {
    SAFENEWARR(pRecvDim, int, SolvCommSize);
    SAFENEWARR(pDispl, int, SolvCommSize+1);
    pDispl[0] = 0;
  } 
  
  SolvComm.Gather(&IntVecDim, 1, MPI::INT, pRecvDim, 1, MPI::INT, 0);
  
  if (!MyRank) {
    for(int i=1; i <= SolvCommSize; i++){
      pDispl[i] = pDispl[i-1] + pRecvDim[i-1];
    }
  }
  
 
  if (!MyRank && pDispl[SolvCommSize] != 0) {
    SAFENEWARR(pSchurDofs, integer, pDispl[SolvCommSize]);
    SAFENEWARR(pDofsRecvdList, integer, pDispl[SolvCommSize]);
  }
  
  SolvComm.Gatherv(pIntDofs, IntVecDim, MPI::INT, pDofsRecvdList, pRecvDim, pDispl, MPI::INT, 0);
  
  if (!MyRank) { 
    for(int i=0; i < pDispl[SolvCommSize]; i++) {
      pSchurDofs[i] = pDofsRecvdList[i];
    }
    
    /*ordino gli indici dei residui locali*/
    std::sort(pSchurDofs, pSchurDofs + pDispl[SolvCommSize]);
    /* elimino le ripetizioni */
    integer* p = std::unique(pSchurDofs, pSchurDofs + pDispl[SolvCommSize]);
    /* dimensione effettiva del residuo locale */
    iSchurIntDim = p - pSchurDofs;
  }
  
  
  /* Vettore di trasformazione locale globale */
  SAFENEWARR(pGlbToLoc, integer, iPrbmSize+1);
  SAFENEWARR(pSchGlbToLoc, integer, iPrbmSize+1);
  for (int i=0; i < iPrbmSize+1; i++) {
    pGlbToLoc[i] = 0;
    pSchGlbToLoc[i] = 0;
  }
  
  
  for (int i=0; i < LocVecDim; i++) {
    pGlbToLoc[pLocDofs[i]] = i + 1;
  }
  
  /* per distinguerli gli indici dei nodi di interfaccia sono negativi */
  for (int i=0; i < IntVecDim; i++) {
    pGlbToLoc[pIntDofs[i]] = -(i + 1);
  }
  
  /* Global to local per la matrice di schur */
  if (!MyRank) {
    for (int i=0; i < iSchurIntDim; i++) {
      pSchGlbToLoc[pSchurDofs[i]] = i+1;
    }
  }


  /* creo i buffer per la ricezione e trasmissione dei messaggi */
 if (SolvCommSize !=1) {
   if (!MyRank) {
     /* i  messaggi + grandi sono le trasmissioni delle matrici di schur */
     integer iTmpTot = 0;
    
     for (int i=0; i < SolvCommSize; i++) {
       iTmpTot += pRecvDim[i] * pRecvDim[i];
     }
     /* buffer di ricezione */
     SAFENEWARR(pBuffer, doublereal, iTmpTot);
   }
   else{
     /* il messaggi + grandi sono le ricezioni dei valori di interfaccia */
     SAFENEWARR(pBuffer, doublereal, IntVecDim);
   }
  }
 
  if (!MyRank){
    SAFENEWARR(pGSReq, MPI::Request, SolvCommSize);
    SAFENEWARR(pGRReq, MPI::Request, SolvCommSize);
    for (int i=0; i < SolvCommSize; i++) { 
    }
  } else {
    SAFENEWARR(pGSReq, MPI::Request, 1);
    SAFENEWARR(pGRReq, MPI::Request, 1);
  }


}


void SchurSolutionManager::StartExchInt(void)
{
  DEBUGCOUT("Entering SchurSolutionManager::StartExchInt()" << endl);
  
  if (SolvCommSize > 1) {
    /* Inizializza Trasmissione di g */
    
#ifdef MPI_PROFILING
    MPE_Log_event(13, 0, "start");
    MPE_Log_send(0, G_TAG, IntVecDim);
#endif /* MPI_PROFILING */      
    pGSReq[0] = SolvComm.Isend(pdgVec, IntVecDim, MPI::DOUBLE, 0, G_TAG);
    
    if(!MyRank) {
#ifdef MPI_PROFILING
      MPE_Log_event(19, 0, "start");
#endif /* MPI_PROFILING */     
      for( int i=0; i < SolvCommSize; i++){
	pGRReq[i] = SolvComm.Irecv(pBuffer + pDispl[i], pRecvDim[i], MPI::DOUBLE, i, G_TAG);
      }
      pSchVH->Reset();
    }
  }
}

void SchurSolutionManager::ComplExchInt(doublereal& dR, doublereal& dXP)
{
  DEBUGCOUT("Entering SchurSolutionManager::ComplExchInt()" << endl);
  if (!MyRank) {
    }
   if (SolvCommSize > 1) {
    /* verifica completamento ricezioni e trasmissione */
      flag SentFlag = 0;
    while (1) {
      SentFlag = pGSReq->Test();
      if (SentFlag) {
#ifdef MPI_PROFILING
	MPE_Log_event(14, 0, "end");
#endif /* MPI_PROFILING */       
	break;
#ifdef USE_MYSLEEP
      } else {
	mysleep(150);
#endif /* USE_MYSLEEP */
      }
    }
    
    if(!MyRank) {
      flag RecvFlag;
      while (1) {
	RecvFlag = MPI::Request::Testall(SolvCommSize, pGRReq);
	if (RecvFlag) {
#ifdef MPI_PROFILING
	  MPE_Log_event(20, 0, "end");
	  for( int i=0; i < SolvCommSize; i++){
	    MPE_Log_receive(i, G_TAG, pRecvDim[i]);
	  }
#endif /* MPI_PROFILING */      	
	  break;
#ifdef USE_MYSLEEP
	} else {
	  mysleep(1500);
#endif /* USE_MYSLEEP */
	}
      }
    }
    
    if (!MyRank) {
      /* Assembla pSVH */ 
      for (int i=0; i < pDispl[SolvCommSize]; i++) {
	pSchVH->fIncCoef(pSchGlbToLoc[pDofsRecvdList[i]], pBuffer[i]);
      }
      
      for (int iCntp1 = 1; iCntp1 <= iSchurIntDim; iCntp1++) {
	doublereal d = pSchVH->dGetCoef(iCntp1);
	dR += d * d;
      }
    }
    
    
    doublereal dTmp[2] = {0, 0};
    dTmp[0] = dR;
    dTmp[1] = dXP;
    
    SolvComm.Send(dTmp, 2, MPI::DOUBLE, 0, G_TAG);
    
    if(!MyRank) {
      for( int i=0; i < SolvCommSize; i++){
	pGRReq[i] = SolvComm.Irecv(pBuffer + 2*i, 2, MPI::DOUBLE, i, G_TAG);
      }
    }
    
    if(!MyRank) {
      flag RecvFlag;
      while (1) {
	RecvFlag = MPI::Request::Testall(SolvCommSize, pGRReq);
	if (RecvFlag) {
	  break;
#ifdef USE_MYSLEEP
	} else {
	  mysleep(1500);
#endif /* USE_MYSLEEP */
	}
      }
    }
    if(!MyRank) {
      for( int i=1; i < SolvCommSize; i++){
	dTmp[0] += pBuffer[2*i];
	dTmp[1] += pBuffer[2*i+1];
      }
      
      for( int i=0; i < SolvCommSize; i++){
	SolvComm.Send(dTmp, 2, MPI::DOUBLE, i, G_TAG);
      }
    } 
    
    pGRReq[0] = SolvComm.Irecv(dTmp, 2, MPI::DOUBLE, 0, G_TAG);
    
    while (1) {
      SentFlag = pGRReq->Test();
      if (SentFlag) {
	break;
#ifdef USE_MYSLEEP
      } else {
	mysleep(150);
#endif /* USE_MYSLEEP */
      }
    }

    dR = dTmp[0];
    dXP = dTmp[1];
  }
}
/* SchurSolutionManager - End */
  /* Inizializza le varie liste */
void InitializeList(int* list, integer dim, integer  value)
{ 
  for (int i=0; i <= dim-1; i++) {
    list[i] = value;
  }
}  


#endif /* USE_MPI */

