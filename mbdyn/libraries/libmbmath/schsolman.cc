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

#include <harwrap.h>
#include <mschwrap.h>
#include <y12wrap.h>
#include <umfpackwrap.h>

#undef min
#undef max
#include <algorithm>

template <class S> SchurSolutionManager::SchurSolutionManager (integer iSize, 
					    integer* pLocalDofs,
					    int iDim1,
					    integer* pInterfDofs,
					    int iDim2,
					    SolutionManager* pLSM,
					    S* pISM,
					    integer iWorkSize,
					    const doublereal& dPivotFactor)
: 
SolvCommSize(0),
iPrbmSize(iSize),
pLocDofs(pLocalDofs),
pIntDofs(pInterfDofs),
iLocVecDim(iDim1),
iIntVecDim(iDim2),
pRecvDim(NULL),
pDispl(NULL),
pDofsRecvdList(NULL),
pSchurDofs(NULL),
iSchurIntDim(0),
pGlbToLoc(NULL),
pSchGlbToLoc(NULL),
pBlockLenght(NULL), 
pTypeDsp(NULL), 
pBuffer(NULL),
ppNewTypes(NULL),
iWorkSpaceSize(iWorkSize), 
pMH(NULL), 
pRVH(NULL),
pSolVH(NULL),
pSchMH(NULL), 
pSchVH(NULL),
pSolSchVH(NULL),
pBMH(NULL), 
pdCM(NULL), 
prVH(NULL),
pgVH(NULL),
pSolrVH(NULL),
pGSReq(NULL), 
pGRReq(NULL),
pLocalSM(pLSM),
pInterSM(pISM),  
fNewMatrix(0) 
{
  DEBUGCOUT("Entering SchurSolutionManager::SchurSolutionManager()"
		  << std::endl);
  
  ASSERT(iPrbmSize > 0);
  ASSERT(pLocDofs != NULL);
  
  /* inizializzo il communicator */
  SolvComm = MBDynComm.Dup();
  SolvCommSize = SolvComm.Get_size();
  MyRank = SolvComm.Get_rank();
  
#ifdef DEBUG
  if (SolvCommSize != 1) {
    	ASSERT(pIntDofs != NULL);
  }
#endif /* DEBUG */
  
  DEBUGCOUT("Solution Communicator Size: " << SolvCommSize << std::endl);
  
  /* determina le dimensioni locali e il vettore di trasferimento Globale -> Locale;
     sul master (myrank == 0) la matrice C è quella di tutte le interfaccie */
  
  InitializeComm();

   /* utilizza iWorkSpaceSize come un coefficiente moltiplicativo */
   integer IntiWorkSpaceSize = iWorkSpaceSize/(iPrbmSize*iPrbmSize)* (iSchurIntDim*iSchurIntDim);
   		
   if(!MyRank) {
  	/* solutore di interfaccia */				
	SAFENEWWITHCONSTRUCTOR(pInterSM,
				S,
				S(iSchurIntDim,
				IntiWorkSpaceSize,
				dPivotFactor));
   }
   /* estrae i puntatori alle matrici e ai vettori creati dai solutori per 
    * l'assemblaggio degli SchurHandlers */
   pBMH    = pLocalSM->pMatHdl();
   prVH    = pLocalSM->pResHdl();
   pSolrVH = pLocalSM->pSolHdl();
   
			
   if(!MyRank) {
   	pSchMH    = pInterSM->pMatHdl();
   	pSchVH    = pInterSM->pResHdl();
	pSolSchVH = pInterSM->pSolHdl();
   }	

  /* Definisce le matrici globali */
  if (prVH->pdGetVec() == pSolrVH->pdGetVec()) {
  	SAFENEWWITHCONSTRUCTOR(pMH, SchurMatrixHandler, 
			 SchurMatrixHandler(iLocVecDim, 
			 		iIntVecDim, pBMH,pGlbToLoc));
  } else {
  	SAFENEWWITHCONSTRUCTOR(pMH, SchurMatrixHandlerUm, 
			 SchurMatrixHandlerUm(iLocVecDim, 
			 		iIntVecDim, pBMH,pGlbToLoc));
  }
  
  SAFENEWWITHCONSTRUCTOR(pRVH, SchurVectorHandler, 
			 SchurVectorHandler(iLocVecDim, 
			 		iIntVecDim, prVH, pGlbToLoc));
  pdCM = pMH->GetCMat();
  pgVH = pRVH->GetIVec();

  SAFENEWWITHCONSTRUCTOR(pSolVH, SchurVectorHandler, 
			 SchurVectorHandler(iLocVecDim, 
			 		iIntVecDim, pSolrVH, pgVH, pGlbToLoc));
     
  /* Creazione di NewType per la trasmissione del vett soluzione calcolato */
  if(!MyRank) {
    	SAFENEWARR(pBlockLenght, int, pDispl[SolvCommSize]);
    	InitializeList(pBlockLenght, pDispl[SolvCommSize], 1);
    	SAFENEWARR(pTypeDsp, MPI::Aint, pDispl[SolvCommSize]);
    
    	MPI::Aint DispTmp = MPI::Get_address(pSolSchVH->pdGetVec());
    	for (int i=0; i < pDispl[SolvCommSize]; i++) {
      		pTypeDsp[i] = MPI::Get_address(pSolSchVH->pdGetVec() + pSchGlbToLoc[pDofsRecvdList[i]] - 1) - DispTmp;
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
  
  std::cout << "Local dimension on process " << MyRank << ": " << iLocVecDim << std::endl;
  if(!MyRank) {
    	std::cout << "Interface dimension: " << iSchurIntDim << std::endl;
  }
  
#ifdef DEBUG
  if(!MyRank) {
    	std::cout << "Interface Dofs " <<std::endl;
    	for (int i=0; i < iSchurIntDim; i++) { 
      		std::cout << pSchurDofs[i] << "  ";
    	}
    	std::cout << std::endl;
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
  if(pMH != NULL) {
	SAFEDELETE(pMH);
  }
  if(pRVH != NULL) {
	SAFEDELETE(pRVH);
  }
  if(pSolVH != NULL) {
	SAFEDELETE(pRVH);
  }
  if(pInterSM != NULL) {
	SAFEDELETE(pInterSM);
  }
  if(pGSReq != NULL) {
	SAFEDELETEARR(pGSReq);
  }
  if(pGRReq != NULL) {
	SAFEDELETEARR(pGRReq);
  }
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
	
	pLocalSM->MatrInit(dResetVal);
	pMH->MatEFCInit(dResetVal);
	if(!MyRank) {
		pInterSM->MatrInit(dResetVal);
	}
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

	/* Fattorizzazione matrice B */
	pLocalSM->Solve();
#ifdef MPI_PROFILING 
  	MPE_Log_event(32, 0, "end");
#endif /* MPI_PROFILING */   
     
	/* prodotto g = g - F*r */
	pMH->CompNewg(*pgVH, *pSolrVH);
#ifdef MPI_PROFILING
   	MPE_Log_event(13, 0, "start");
   	MPE_Log_send(0, G_TAG, iIntVecDim);
#endif /* MPI_PROFILING */   

  	/* Inizializza Trasmissione di g */
   	pGSReq[0] = SolvComm.Isend(pgVH->pdGetVec(), iIntVecDim, MPI::DOUBLE, 0, G_TAG);
  
   	if(!MyRank) {
#ifdef MPI_PROFILING
     		MPE_Log_event(19, 0, "start");
#endif /* MPI_PROFILING */   
     		for( int i=0; i < SolvCommSize; i++){
       			pGRReq[i] = SolvComm.Irecv(pBuffer + pDispl[i], pRecvDim[i], MPI::DOUBLE, i, G_TAG);
     		}
   	}
	
     /* Calcolo di E'. Va fatto solo se  La matrice è stata rifattorizzata*/   
   	if (fNewMatrix) {
     		for (int i=0; i < iIntVecDim; i++) {
   			pLocalSM->ChangeResPoint(pMH->GetECol(i));
			pLocalSM->ChangeSolPoint(pMH->GetEColSol(i));			
       			/* fa solo la back Substion perche' 
			 * e' stato già lanciato il solve al precedentemente */
			pLocalSM->Solve();
     		}
	}
   	pLocalSM->ChangeResPoint(prVH->pdGetVec());
	pLocalSM->ChangeSolPoint(pSolrVH->pdGetVec());			
	
   	/* verifica completamento ricezioni e trasmissione vettore g*/
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
  
  	/* Assemblaggio del vettore delle incognite g sulla macchina master */  
  	if (!MyRank) {
    		/* Assembla pSchVH */
		pSchVH->Reset(0.);
    		for (int i=0; i < pDispl[SolvCommSize]; i++) {
      			pSchVH->fIncCoef(pSchGlbToLoc[pDofsRecvdList[i]], pBuffer[i]);
    		}
  	}      
  
  	/* assembla le matrici di schur locali se è stata 
	 * appena assemblata una nuova matrice di partenza */
  	if (fNewMatrix) {
    		AssSchur();
    		fNewMatrix = flag(0);
  	}
  	/* risoluzione schur pb sul master */
  	if ((!MyRank) && (SolvCommSize != 1)) {
#ifdef MPI_PROFILING 
    		MPE_Log_event(35, 0, "start");
#endif /* MPI_PROFILING */   
    		pInterSM->Solve();
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
      			pGRReq[i] = SolvComm.Isend(pSolSchVH->pdGetVec(), 1, *ppNewTypes[i], i, G_TAG);	 
#ifdef MPI_PROFILING
      			MPE_Log_send(i, G_TAG, pRecvDim[i]);  
#endif /* MPI_PROFILING */     
    		}	
  	}
    
#ifdef MPI_PROFILING
   	MPE_Log_event(19, 0, "start");
#endif /* MPI_PROFILING */   

  	pGSReq[0] = SolvComm.Irecv(pgVH->pdGetVec(), iIntVecDim, MPI::DOUBLE, 0, G_TAG);
  
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
    			MPE_Log_receive(0, G_TAG, iIntVecDim);
#endif /* MPI_PROFILING */ 
      			break;
#ifdef USE_MYSLEEP
    		} else {
      			mysleep(1500);
#endif /* USE_MYSLEEP */
    		}
  	}   
  	/* calcolo la soluzione corretta per x = Solr - E'*g */
	pMH->CompNewf(*pSolrVH, *pgVH);
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
	pMH->CompLocSchur();

  	/* Trasmette le Schur locali */
#ifdef MPI_PROFILING
   	MPE_Log_event(13, 0, "start");
   	MPE_Log_send(0, S_TAG, iIntVecDim*iIntVecDim);
#endif /* MPI_PROFILING */   
  
  	pGSReq[0] = SolvComm.Isend(pdCM, iIntVecDim*iIntVecDim, MPI::DOUBLE, 0, S_TAG);
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
;
	}
#ifdef MPI_PROFILING
  	MPE_Log_event(42, 0, "end");
#endif /* MPI_PROFILING */ 
}



void SchurSolutionManager::InitializeComm(void)
{
  DEBUGCOUT("Entering SchurSolutionManager::InitializeComm()" << endl);

   /* il master riceve informazioni sulle dimensioni dei vettori scambiati */
  if (!MyRank) {
  	SAFENEWARR(pRecvDim, int, SolvCommSize);
    	SAFENEWARR(pDispl, int, SolvCommSize+1);
    	pDispl[0] = 0;
  } 
   /* trasmiussione dimensioni interfacce locali */
  SolvComm.Gather(&iIntVecDim, 1, MPI::INT, pRecvDim, 1, MPI::INT, 0);
  if (!MyRank) {
    	for(int i=1; i <= SolvCommSize; i++){
      		pDispl[i] = pDispl[i-1] + pRecvDim[i-1];
    	}
  }
  
 
  if (!MyRank && pDispl[SolvCommSize] != 0) {
    	SAFENEWARR(pSchurDofs, integer, pDispl[SolvCommSize]);
    	SAFENEWARR(pDofsRecvdList, integer, pDispl[SolvCommSize]);
  }
  
  SolvComm.Gatherv(pIntDofs, iIntVecDim, MPI::INT, pDofsRecvdList, pRecvDim, pDispl, MPI::INT, 0);
  
  if (!MyRank) { 
    	for(int i=0; i < pDispl[SolvCommSize]; i++) {
      		pSchurDofs[i] = pDofsRecvdList[i];
    	}
    	/*ordina gli indici dei residui locali*/
    	std::sort(pSchurDofs, pSchurDofs + pDispl[SolvCommSize]);
    	/* elimina le ripetizioni */
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
  
  
  for (int i=0; i < iLocVecDim; i++) {
    	pGlbToLoc[pLocDofs[i]] = i + 1;
  }
  
  /* per distinguerli gli indici dei nodi di interfaccia sono negativi */
  for (int i=0; i < iIntVecDim; i++) {
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
     		/* i  messaggi + grandi sono legati alla ricezione  delle matrici di schur */
     		integer iTmpTot = 0;
    		for (int i=0; i < SolvCommSize; i++) {
       			iTmpTot += pRecvDim[i] * pRecvDim[i];
     		}
     		/* buffer di ricezione */
     		SAFENEWARR(pBuffer, doublereal, iTmpTot);
   	} else {
     		/* il messaggi + grandi sono le ricezioni dei valori di interfaccia */
     		SAFENEWARR(pBuffer, doublereal, iIntVecDim);
   	}
  }
 
  if (!MyRank){
    	SAFENEWARR(pGSReq, MPI::Request, SolvCommSize);
    	SAFENEWARR(pGRReq, MPI::Request, SolvCommSize);
  } else {
    	SAFENEWARR(pGSReq, MPI::Request, 1);
    	SAFENEWARR(pGRReq, MPI::Request, 1);
  }


}


void
SchurSolutionManager::StartExchInt(void)
{
  	DEBUGCOUT("Entering SchurSolutionManager::StartExchInt()" << endl);
  
  	if (SolvCommSize > 1) {
    	/* Inizializza Trasmissione di g */
    
#ifdef MPI_PROFILING
    		MPE_Log_event(13, 0, "start");
    		MPE_Log_send(0, G_TAG, iIntVecDim);
#endif /* MPI_PROFILING */      
    		pGSReq[0] = SolvComm.Isend(pgVH->pdGetVec(), iIntVecDim, MPI::DOUBLE, 0, G_TAG);
    
    		if(!MyRank) {
#ifdef MPI_PROFILING
      			MPE_Log_event(19, 0, "start");
#endif /* MPI_PROFILING */     
      			for (int i = 0; i < SolvCommSize; i++) {
				pGRReq[i] = SolvComm.Irecv(pBuffer + pDispl[i],
						pRecvDim[i], MPI::DOUBLE, 
						i, G_TAG);
      			}
    		}
  	}
}

void
SchurSolutionManager::ComplExchInt(doublereal& dR, doublereal& dXP)
{
  	DEBUGCOUT("Entering SchurSolutionManager::ComplExchInt()" << endl);
   	if (SolvCommSize > 1) {
    		doublereal dTmp[2] = {0, 0};
    		dTmp[0] = dR;
    		dTmp[1] = dXP;
    		SolvComm.Send(dTmp, 2, MPI::DOUBLE, 0, G_TAG+100);
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
      			/* Assembla pSVH */ 
      			pSchVH->Reset(0.);
      			for (int i=0; i < pDispl[SolvCommSize]; i++) {
				pSchVH->fIncCoef(pSchGlbToLoc[pDofsRecvdList[i]], pBuffer[i]);
      			}
      			for (int iCntp1 = 1; iCntp1 <= iSchurIntDim; iCntp1++) {
				doublereal d = pSchVH->dGetCoef(iCntp1);
				dR += d * d;
      			}
    		}
    		if(!MyRank) {
      			for( int i=0; i < SolvCommSize; i++){
				pGRReq[i] = SolvComm.Irecv(pBuffer + 2*i, 2, MPI::DOUBLE, i, G_TAG+100);
      			}
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

#ifdef USE_Y12
template SchurSolutionManager::SchurSolutionManager(integer iSize, 
					    integer* pLocalDofs,
					    int iDim1,
					    integer* pInterfDofs,
					    int iDim2,
					    SolutionManager* pLSM,
					    Y12SparseLUSolutionManager* pISM,
					    integer iWorkSize,
					    const doublereal& dPivotFactor);
#endif //USE_Y12
#ifdef USE_MESCHACH					    
template SchurSolutionManager::SchurSolutionManager(integer iSize, 
					    integer* pLocalDofs,
					    int iDim1,
					    integer* pInterfDofs,
					    int iDim2,
					    SolutionManager* pLSM,
					    MeschachSparseLUSolutionManager* pISM,
					    integer iWorkSize,
					    const doublereal& dPivotFactor); 
#endif /* USE_MESCHACH */
#ifdef USE_UMFPACK
template SchurSolutionManager::SchurSolutionManager(integer iSize, 
					    integer* pLocalDofs,
					    int iDim1,
					    integer* pInterfDofs,
					    int iDim2,
					    SolutionManager* pLSM,
					    UmfpackSparseLUSolutionManager* pISM,
					    integer iWorkSize,
					    const doublereal& dPivotFactor); 
#endif /* USE_UMFPACK */
#endif /* USE_MPI */


