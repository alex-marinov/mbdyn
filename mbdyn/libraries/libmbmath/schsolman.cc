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

/*  SchurMatrixHandler 
    SchurVectorHandler
    SchurSolutionManager */
    
/* 
 * Copyright 1999-2000 Giuseppe Quaranta <giuquaranta@tiscalinet.it>
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

/*  BroySchMatrixHandler - Start */

SchurMatrixHandler::SchurMatrixHandler(int LocSize, int IntSize,
				       MatrixHandler* pMatB,
				       doublereal* pMatE,
				       doublereal* pMatF,
				       doublereal* pMatC,
				       integer* pGlobToLoc)
: LSize(LocSize),
ISize(IntSize),
pB(pMatB),
pE(pMatE),
pF(pMatF),
pC(pMatC),
pGTL(pGlobToLoc),
dZero(0.)
  
{ 
#ifdef DEBUG
  IsValid();
#endif  
}

SchurMatrixHandler::~SchurMatrixHandler(void)
{
 NO_OP;
}

void SchurMatrixHandler::IsValid(void) const
{
  ASSERT(LSize >0);
  ASSERT(ISize >0);
  ASSERT(pB != NULL);
  ASSERT(pE != NULL);
  ASSERT(pF != NULL);
  ASSERT(pC != NULL);
  ASSERT(pGTL != NULL);
}

/* SchurMatrixHandler - End */


/* SchurVectorHandler - Start */

SchurVectorHandler::SchurVectorHandler(int LocSize, int IntSize,
				       doublereal* pLocVec,
				       doublereal* pIntVec,
				       integer* pGlobToLoc)
:  LSize(LocSize),
ISize(IntSize),
pLV(pLocVec),
pIV(pIntVec),
pGTL(pGlobToLoc),
dZero(0.)

{
#ifdef DEBUG
  IsValid();
#endif 
}

SchurVectorHandler::~SchurVectorHandler(void)
{
NO_OP;
}

void SchurVectorHandler::IsValid(void) const
{
  ASSERT(LSize >0);
  ASSERT(ISize >0);
  ASSERT(pLV != NULL);
  ASSERT(pIV != NULL);
  ASSERT(pGTL != NULL);
}

/* SchurSolutionManager - Start */

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
pBlockLenght(NULL), pTypeDsp(NULL), pBuffer(NULL), ppNewTypes(NULL),
piBRow(NULL), piBCol(NULL), pdBMat(NULL), pdrVec(NULL), 
pBMH(NULL), pLU(NULL), pdEMat(NULL), pdFMat(NULL), pdCMat(NULL), 
pdgVec(NULL), 
pMH(NULL), pRVH(NULL),
pSchSM(NULL), pSchMH(NULL), pSchVH(NULL),
pGSReq(NULL), pGRReq(NULL),  
iWorkSpaceSize(0), fNewMatrix(0) 
{
  const char sFuncName[] = "BroySchSolutionManager::BroySchSolutionManager()";
 
 DEBUGCOUT("Entering " << sFuncName << endl);
  
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
#endif
  
  DEBUGCOUT("Solution Communicator Size: " << SolvCommSize << endl);
  
  /* detrmino le dimensioni locali e il vettore di trasferimento Globale -> Locale;
     sul master (myrank == 0) la matrice C è quella di tutte le interfaccie */
  
  InitializeComm();
  

  ASSERT((dPivotFactor >= 0.0) && (dPivotFactor <= 1.0));


  /* Valore di default */
  
  iWorkSpaceSize = LocVecDim;

  /* Alloca arrays per la matrice B e il vettore r */
  SAFENEWARR(pdBMat, doublereal, LocVecDim * LocVecDim, SMmm);
  SAFENEWARR(piBRow, integer, iWorkSpaceSize*iWorkSpaceSize, SMmm);
  SAFENEWARR(piBCol, integer, iWorkSpaceSize*iWorkSpaceSize, SMmm);
  SAFENEWARR(pdrVec, doublereal, iWorkSpaceSize, SMmm);

  /* Alloca arrays per il solutore */
  
  
  SAFENEWWITHCONSTRUCTOR(pBMH, 
                         SparseMatrixHandler, 
			 SparseMatrixHandler(iWorkSpaceSize, &piBRow, 
			 		     &piBCol, &pdBMat,
 					     iWorkSpaceSize*iWorkSpaceSize),
 			 SMmm);
  
  SAFENEWWITHCONSTRUCTOR(pLU, 
			 HarwellLUSolver,
			 HarwellLUSolver(iWorkSpaceSize, iWorkSpaceSize*iWorkSpaceSize, &piBRow, &piBCol, 
					 &pdBMat, pdrVec, dPivotFactor), SMmm);
  

  /* Alloca la matrice E  per colonne */
  SAFENEWARR(pdEMat, doublereal, LocVecDim * IntVecDim, SMmm);
  
  /* Alloca la matrice F */
  SAFENEWARR(pdFMat, doublereal, IntVecDim * LocVecDim, SMmm);

  /* matrice C e vettore g*/
  SAFENEWARR(pdCMat, doublereal, IntVecDim*IntVecDim, SMmm);
  SAFENEWARR(pdgVec, doublereal, IntVecDim, SMmm);
  
  /* matrice di Shur sul master */
  if(!MyRank) {
     SAFENEWWITHCONSTRUCTOR(pSchSM, MeschachSparseLUSolutionManager,
			    MeschachSparseLUSolutionManager(iSchurIntDim, iSchurIntDim),
			    SMmm);
     pSchMH = pSchSM->pMatHdl();
     pSchVH = pSchSM->pResHdl();
  }
  
  
  
  /* Hadlers per la matrice per il residuo e per la soluzione */
  SAFENEWWITHCONSTRUCTOR(pMH, SchurMatrixHandler, 
			 SchurMatrixHandler(LocVecDim, IntVecDim, pBMH, 
					    pdEMat, pdFMat, pdCMat, pGlbToLoc),
			 SMmm);
  
  SAFENEWWITHCONSTRUCTOR(pRVH, SchurVectorHandler, 
			 SchurVectorHandler(LocVecDim, IntVecDim, pdrVec,
					    pdgVec, pGlbToLoc),
			 SMmm);
  
  /* Creazione di NewType per la trasmissione del vett soluzione calcolato */
  if(!MyRank) {
    SAFENEWARR(pBlockLenght, int, pDispl[SolvCommSize], SMmm);
    InitializeList(pBlockLenght, pDispl[SolvCommSize], 1);
    SAFENEWARR(pTypeDsp, MPI::Aint, pDispl[SolvCommSize] , SMmm);
    
    MPI::Aint DispTmp = MPI::Get_address(pSchVH->pdGetVec());
    for (int i=0; i < pDispl[SolvCommSize]; i++) {
      pTypeDsp[i] = MPI::Get_address(pSchVH->pdGetVec() + pSchGlbToLoc[pDofsRecvdList[i]] - 1) - DispTmp;
    }
    
    SAFENEWARR(ppNewTypes, MPI::Datatype*, SolvCommSize, SMmm);
    MPI::Aint* pActualDispl = pTypeDsp;
    for (int i=0; i < SolvCommSize; i++) {
      ppNewTypes[i] = NULL;
      SAFENEWWITHCONSTRUCTOR(ppNewTypes[i], MPI::Datatype, MPI::Datatype(MPI::DOUBLE.Create_hindexed(pRecvDim[i], pBlockLenght, pActualDispl)), SMmm);
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
#endif      
}

SchurSolutionManager::~SchurSolutionManager(void)
{

 const char sFuncName[] = "SchurSolutionManager::~SchurSolutionManager()";
   
  DEBUGCOUT("Entering " << sFuncName << endl);

#ifdef DEBUG
  IsValid();
#endif      
 /* FIXME : la cancellazione di questi arrays da alcuni problemi
	    in chiusura. Per ora viene commentata .
  */
  /*
  if(pLocDofs != NULL) {
    SAFEDELETEARR(pLocDofs, SMmm);
  }
  if(pIntDofs != NULL) {
    SAFEDELETEARR(pIntDofs, SMmm);
  } 
  if(pRecvDim != NULL) {
    SAFEDELETEARR(pRecvDim, SMmm);
  }
  if(pDispl != NULL) {
    SAFEDELETEARR(pDispl, SMmm);
  }
  if(pDofsRecvdList != NULL) {
    SAFEDELETEARR(pDofsRecvdList, SMmm);
  }
  if(pSchurDofs != NULL) {
    SAFEDELETEARR(pSchurDofs, SMmm);
  }
  if(pGlbToLoc != NULL) {
    SAFEDELETEARR(pGlbToLoc, SMmm);
  }
  if(pSchGlbToLoc != NULL) {
    SAFEDELETEARR(pSchGlbToLoc, SMmm);
  }
  if(pBlockLenght != NULL) {
    SAFEDELETEARR(pBlockLenght, SMmm);
  }
   if(pTypeDsp != NULL) {
    SAFEDELETEARR(pTypeDsp, SMmm);
  }
  if(pBuffer != NULL) {
    SAFEDELETEARR(pBuffer, SMmm);
  }
  if(ppNewTypes != NULL) {
    for(int i=0; i < SolvCommSize; i++) {
       ppNewTypes[i]->Free();
      SAFEDELETE(ppNewTypes[i], SMmm); 	
    }
   SAFEDELETEARR(ppNewTypes, SMmm);
  }
  if(piBRow != NULL) {
    SAFEDELETEARR(piBRow, SMmm);
  }
  if(piBCol != NULL) {
    SAFEDELETEARR(piBCol, SMmm);
  }
  if(pdBMat != NULL) {
    SAFEDELETEARR(pdBMat, SMmm);
  }
  if(pdrVec != NULL) {
    SAFEDELETEARR(pdrVec, SMmm);
  }
  if(pBMH != NULL) {
      SAFEDELETE(pBMH, SMmm);
  }
  if(pLU != NULL) {
    SAFEDELETE(pLU, SMmm);
  }
  if(pdEMat != NULL) {
    SAFEDELETEARR(pdEMat, SMmm);
  }
  if(pdFMat != NULL) {
    SAFEDELETEARR(pdFMat, SMmm);
  }
  if(pdCMat != NULL) {
    SAFEDELETEARR(pdCMat, SMmm);
  }
  if(pdgVec != NULL) {
    SAFEDELETEARR(pdgVec, SMmm);
  }
  if(pMH != NULL) {
	SAFEDELETE(pMH, SMmm);
  }
  if(pRVH != NULL) {
	SAFEDELETE(pRVH, SMmm);
  }
  if(pSchSM != NULL) {
      SAFEDELETEARR(pSchSM, SMmm);
  }
  if(pSchMH != NULL) {
    SAFEDELETEARR(pSchMH, SMmm);
  }
  if(pSchVH != NULL) {
    SAFEDELETEARR(pSchVH, SMmm);
  }
  if(pGSReq != NULL) {
    SAFEDELETEARR(pGSReq, SMmm);
  }
  if(pGRReq != NULL) {
    SAFEDELETEARR(pGRReq, SMmm);
  }
 */
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
#endif
  pMH->Init(dResetVal);
  fNewMatrix = flag(1);
}

/* Risolve i blocchi */

void SchurSolutionManager::Solve(void)
{
  const char sFuncName[] = "SchurSolutionManager::Solve()";
  
  DEBUGCOUT("Entering " << sFuncName << endl);
#ifdef DEBUG
  IsValid();
#endif

#ifdef MPI_PROFILING 
  MPE_Log_event(31, 0, "start");
#endif /* MPI_PROFILING */   
  
  if(fNewMatrix) {
    integer iCount = pBMH->PacMat();
    pLU->iNonZeroes = iCount;
    flag fReturnFlag = pLU->fLUFactor();
    if (fReturnFlag < 0) {	 
      THROW(HSLUSolutionManager::ErrGeneric());
    }
  } 
  
  pLU->pdRhs = pdrVec; 
  
  pLU->Solve();
#ifdef MPI_PROFILING 
  MPE_Log_event(32, 0, "end");
#endif /* MPI_PROFILING */   
  
  InitializeList(pdgVec, IntVecDim, 0.);
  
  /* prodotto F * r */
  for (int j=0; j < LocVecDim; j++) {
    int iColx = j * IntVecDim;
    for (int i=0; i < IntVecDim; i++) {
      if(pdFMat[i + iColx] != 0) {
	pdgVec[i] -= pdFMat[i + iColx] * pdrVec[j];
      }
    }
  } 
  
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
       pLU->pdRhs = pdEMat + i*LocVecDim; 
       pLU->Solve();
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
  
   
  /* risoluzione schur pb sul master */
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
  for( int j=0; j < IntVecDim; j++) {  
    int iColx = j * LocVecDim;
    for (int i=0; i < LocVecDim; i++) {
      if (pdEMat[i + iColx] != 0) {
	pdrVec[i] -= pdEMat[i + iColx] * pdgVec[j];
      }
    }
  }
}



void SchurSolutionManager::AssSchur(void) 
{
#ifdef MPI_PROFILING
   MPE_Log_event(41, 0, "start");
#endif /* MPI_PROFILING */ 
  const char sFuncName[] = "SchurSolutionManager::AssSchur()";
  
  DEBUGCOUT("Entering " << sFuncName << endl);
  
#ifdef DEBUG
  IsValid();
#endif
  
  /* Calcola le Schur locali */
  for(int j=0; j < IntVecDim; j++) {
    int iColc = j * IntVecDim;
    int iCole = j * LocVecDim;
    for (int k=0; k < LocVecDim; k++) {
      int iColf = k * IntVecDim;
      for( int i=0; i < IntVecDim; i++) {
	pdCMat[i + iColc] -=  pdFMat[i + iColf] * pdEMat[iCole + k];
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

    
    //InitializeList(pdSchMat, iSchurIntDim * iSchurIntDim, 0.);
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
  if (!MyRank) {
    SAFENEWARR(pRecvDim, int, SolvCommSize, SMmm);
    SAFENEWARR(pDispl, int, SolvCommSize+1, SMmm);
    pDispl[0] = 0;
  } 
  
  SolvComm.Gather(&IntVecDim, 1, MPI::INT, pRecvDim, 1, MPI::INT, 0);
  
  if (!MyRank) {
    for(int i=1; i <= SolvCommSize; i++){
      pDispl[i] = pDispl[i-1] + pRecvDim[i-1];
    }
  }
  
 
  if (!MyRank && pDispl[SolvCommSize] != 0) {
    SAFENEWARR(pSchurDofs, integer, pDispl[SolvCommSize], SMmm);
    SAFENEWARR(pDofsRecvdList, integer, pDispl[SolvCommSize], SMmm);
  }
  
  SolvComm.Gatherv(pIntDofs, IntVecDim, MPI::INT, pDofsRecvdList, pRecvDim, pDispl, MPI::INT, 0);
  
  if (!MyRank) { 
    for(int i=0; i < pDispl[SolvCommSize]; i++) {
      pSchurDofs[i] = pDofsRecvdList[i];
    }
    
    /*ordino gli indici dei residulocali*/
    sort(pSchurDofs, pSchurDofs + pDispl[SolvCommSize]);
    /* elimino le ripetizioni */
    integer* p = unique(pSchurDofs, pSchurDofs + pDispl[SolvCommSize]);
    /* dimensione effettiva del residuo locale */
    iSchurIntDim = p - pSchurDofs;
  }
  
  
  /* Vettore di trasformazione locale globale */
  SAFENEWARR(pGlbToLoc, integer, iPrbmSize+1, SMmm);
  SAFENEWARR(pSchGlbToLoc, integer, iPrbmSize+1, SMmm);
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
     SAFENEWARR(pBuffer, doublereal, iTmpTot, SMmm);
   }
   else{
     /* il messaggi + grandi sono le ricezioni dei valori di interfaccia */
     SAFENEWARR(pBuffer, doublereal, IntVecDim, SMmm);
   }
  }
 

  if (!MyRank){
    SAFENEWARR(pGSReq, MPI::Request, SolvCommSize, SMmm);
    SAFENEWARR(pGRReq, MPI::Request, SolvCommSize, SMmm);
  } else {
    SAFENEWARR(pGSReq, MPI::Request, 1, SMmm);
    SAFENEWARR(pGRReq, MPI::Request, 1, SMmm);
  }
}


void SchurSolutionManager::StartExchInt(void)
{
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

void InitializeList(doublereal* list, integer dim, doublereal  value)
{ 
  for (int i=0; i <= dim-1; i++) {
    list[i] = value;
  }
}

#endif /* USE_MPI */

