/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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
 * Copyright 1999-2015 Giuseppe Quaranta <quaranta@aero.polimi.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <iomanip>

#ifdef USE_MPI
#include "solman.h"
#include "parser.h"
#include "linsol.h"
#include "schsolman.h"
#include "mysleep.h"

/* FIXME: we must move the test to libmbmath */
#include "nonlin.h"

/* NOTE: define to use Wait/Waitall instead of Test/Testall
 * Apparently, this results in far better performances,
 * so we might want to extend it to all other communications */
#define USE_MPI_WAIT

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

SchurSolutionManager::SchurSolutionManager(integer iSize,
		integer iBlocks,
		integer* pLocalDofs,
		int iDim1,
		integer* pInterfDofs,
		int iDim2,
		SolutionManager* pLSM,
		LinSol &ls)
:
SolvCommSize(0),
iPrbmSize(iSize),
iPrbmBlocks(iBlocks),
iBlkSize(iSize),
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
ppNewTypes(NULL),
pBuffer(NULL),
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
pInterSM(0),
bNewMatrix(false)
{
	DEBUGCOUT("Entering SchurSolutionManager::SchurSolutionManager()"
      			<< std::endl);

	ASSERT(iPrbmSize > 0);
	ASSERT(pLocDofs != NULL);

	/* inizializzo il communicator */
	SolvComm = MBDynComm.Dup();
	SolvCommSize = SolvComm.Get_size();
	MyRank = SolvComm.Get_rank();

	if (SolvCommSize <= 1) {
		silent_cerr("SchurSolutionManager: "
			"invalid communicator size " << SolvCommSize
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

#ifdef DEBUG
	ASSERT(pIntDofs != NULL);
#endif /* DEBUG */

	DEBUGCOUT("Solution Communicator Size: " << SolvCommSize << std::endl);

	/* determina le dimensioni locali e il vettore
	 * di trasferimento Globale -> Locale;
         * sul master (MyRank == 0) la matrice C e' quella
	 * di tutte le interfacce */

	/* initializes iSchurIntDim among other stuff... */
	InitializeComm();

	/* utilizza iWorkSpaceSize come un coefficiente moltiplicativo */

	/* NOTE: iWorkSpaceSize is likely to be 0, since now most of the
	 * linear solvers don't need any; this is not a big deal, resulting
	 * in GetSolutionManager() being called with a 0 iWorkSpaceSize arg */
	iWorkSpaceSize = ls.iGetWorkSpaceSize();
	integer iIntWorkSpaceSize = iWorkSpaceSize*(iSchurIntDim*iSchurIntDim)/(iPrbmSize*iPrbmSize);

	if (MyRank == 0) {
		pInterSM = ls.GetSolutionManager(iSchurIntDim,
				iIntWorkSpaceSize);
	}

	/* estrae i puntatori alle matrici e ai vettori creati
	 * dai solutori per l'assemblaggio degli SchurHandlers */
	pBMH    = pLocalSM->pMatHdl();
	prVH    = pLocalSM->pResHdl();
	pSolrVH = pLocalSM->pSolHdl();

	if (MyRank == 0) {
		pSchMH    = pInterSM->pMatHdl();
		pSchVH    = pInterSM->pResHdl();
		pSolSchVH = pInterSM->pSolHdl();
	}

	/* Definisce le matrici globali */
	if (prVH->pdGetVec() == pSolrVH->pdGetVec()) {
		SAFENEWWITHCONSTRUCTOR(pMH, SchurMatrixHandler,
				SchurMatrixHandler(iLocVecDim,
					iIntVecDim, pBMH, pGlbToLoc));

	} else {
		SAFENEWWITHCONSTRUCTOR(pMH, SchurMatrixHandlerUm,
       				SchurMatrixHandlerUm(iLocVecDim,
					iIntVecDim, pBMH, pGlbToLoc));
	}

	SAFENEWWITHCONSTRUCTOR(pRVH, SchurVectorHandler,
			SchurVectorHandler(iLocVecDim,
				iIntVecDim, prVH, pGlbToLoc));

      	pdCM = pMH->GetCMat();

	/* NOTE: from my understanding, pgVH is the interface portion
	 * of the pRVH (the SchurVectorhandler of the residual), which
	 * is allocated internally by the constructor; it is passed
	 * to pSolVH (the SchurVectorhandler of the solution) as well,
	 * to act as interface portion of the solution, so the residual
	 * and the solution share the same memory area. */
	pgVH = pRVH->GetIVec();

	SAFENEWWITHCONSTRUCTOR(pSolVH, SchurVectorHandler,
			SchurVectorHandler(iLocVecDim,
				iIntVecDim, pSolrVH, pgVH, pGlbToLoc));

	/* Creazione di NewType per la trasmissione
	 * del vettore soluzione calcolato */
	if (MyRank == 0) {
		SAFENEWARR(pBlockLenght, int, pDispl[SolvCommSize]);
		InitializeList(pBlockLenght, pDispl[SolvCommSize], 1);
		SAFENEWARR(pTypeDsp, MPI::Aint, pDispl[SolvCommSize]);
		
		doublereal *pdm1 = pSolSchVH->pdGetVec() - 1;
		MPI::Aint DispTmp = MPI::Get_address(pdm1 + 1);
		for (int i = 0; i < pDispl[SolvCommSize]; i++) {
			int j = i%iBlkSize;
			int blk = i/iBlkSize;
			pTypeDsp[i] = MPI::Get_address(&pdm1[pSchGlbToLoc[pDofsRecvdList[j] + blk*iBlkSize]]) - DispTmp;
		}

		SAFENEWARR(ppNewTypes, MPI::Datatype*, SolvCommSize);
		MPI::Aint* pCurrDispl = pTypeDsp;

		for (int i = 0; i < SolvCommSize; i++) {
			ppNewTypes[i] = NULL;
			SAFENEWWITHCONSTRUCTOR(ppNewTypes[i], MPI::Datatype,
					MPI::Datatype(MPI::DOUBLE.Create_hindexed(pRecvDim[i],
							pBlockLenght, pCurrDispl)));
			ppNewTypes[i]->Commit();
			pCurrDispl += pRecvDim[i];
		}
	}

	silent_cout("Local dimension on process " << MyRank
		<< ": " << iLocVecDim << std::endl);

	if (MyRank == 0) {
		silent_cout("Interface dimension: " << iSchurIntDim
			<< std::endl);
	}

#ifdef DEBUG
	if (MyRank == 0) {
		silent_cout("Interface Dofs " << std::endl);
		for (int i = 0; i < iSchurIntDim; i++) {
			silent_cout(pSchurDofs[i] << "  ");
		}
		silent_cout(std::endl);
	}

	IsValid();
#endif /* DEBUG */
}

SchurSolutionManager::~SchurSolutionManager(void)
{
	DEBUGCOUT("Entering SchurSolutionManager::~SchurSolutionManager()" 
			<< endl);

#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	if (pRecvDim != NULL) {
		SAFEDELETEARR(pRecvDim);
	}

	if (pDispl != NULL) {
		SAFEDELETEARR(pDispl);
	}

	if (pDofsRecvdList != NULL) {
		SAFEDELETEARR(pDofsRecvdList);
	}

	if (pSchurDofs != NULL) {
		SAFEDELETEARR(pSchurDofs);
	}

	if (pGlbToLoc != NULL) {
		SAFEDELETEARR(pGlbToLoc);
	}

	if (pSchGlbToLoc != NULL) {
		SAFEDELETEARR(pSchGlbToLoc);
	}

	if (pBlockLenght != NULL) {
		SAFEDELETEARR(pBlockLenght);
	}

	if (pTypeDsp != NULL) {
		SAFEDELETEARR(pTypeDsp);
	}

	if (pBuffer != NULL) {
		SAFEDELETEARR(pBuffer);
	}

	if (ppNewTypes != NULL) {
		for (int i = 0; i < SolvCommSize; i++) {
			ppNewTypes[i]->Free();
			SAFEDELETE(ppNewTypes[i]);
		}

		SAFEDELETEARR(ppNewTypes);
	}

	if (pMH != NULL) {
		SAFEDELETE(pMH);
	}

	if (pRVH != NULL) {
		SAFEDELETE(pRVH);
	}

	if (pSolVH != NULL) {
		SAFEDELETE(pSolVH);
	}

	if (pInterSM != NULL) {
		SAFEDELETE(pInterSM);
	}

	if (pGSReq != NULL) {
		SAFEDELETEARR(pGSReq);
	}

	if (pGRReq != NULL) {
		SAFEDELETEARR(pGRReq);
	}
}

#ifdef DEBUG
void
SchurSolutionManager::IsValid(void) const
{
	NO_OP;
}
#endif /* DEBUG */

/* Inizializza il gestore delle matrici */

void
SchurSolutionManager::MatrReset(void)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	pLocalSM->MatrReset();
	pMH->MatEFCReset();
	if (MyRank == 0) {
		pInterSM->MatrReset();
	}
	bNewMatrix = true;
}

/* Inizializzatore "speciale" */
void
SchurSolutionManager::MatrInitialize()
{
	pLocalSM->MatrInitialize();
	pMH->MatEFCReset();
	if (MyRank == 0) {
		pInterSM->MatrInitialize();
	}
	bNewMatrix = true;
}
	
/* Risolve i blocchi */

void
SchurSolutionManager::Solve(void)
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

	/* NOTE: we need to set the matrix handler after each solution,
	 * because when the local SM is using an optimized sparse matrix,
	 * it may switch from a generic representation, e.g. SpMapMH,
	 * to a compressed representation, e.g. CCMH or DirMH
	 */
	pMH->SetBMat(pLocalSM->pMatHdl());

#ifdef MPI_PROFILING
  	MPE_Log_event(32, 0, "end");
#endif /* MPI_PROFILING */

	/* prodotto g = g - F*r */
	pMH->CompNewg(*pgVH, *pSolrVH);

#ifdef MPI_PROFILING
   	MPE_Log_event(13, 0, "start");
   	MPE_Log_send(0, G_TAG, iIntVecDim);
#endif /* MPI_PROFILING */

	/* on master node... */
	if (MyRank == 0) {
#ifdef MPI_PROFILING
     		MPE_Log_event(19, 0, "start");
#endif /* MPI_PROFILING */

		/* receive from remote only... */
     		for (int i = 1; i < SolvCommSize; i++){
       			pGRReq[i] = SolvComm.Irecv(&pBuffer[pDispl[i]],
					pRecvDim[i], MPI::DOUBLE, i, G_TAG);
     		}

		/* emulate transmission...*/
		(void)memmove(&pBuffer[pDispl[0]], pgVH->pdGetVec(),
				sizeof(double)*pRecvDim[0]);

	/* on remote nodes... */
	} else {
 	 	/* Inizializza Trasmissione di g */
   		pGSReq[0] = SolvComm.Isend(pgVH->pdGetVec(), iIntVecDim,
				MPI::DOUBLE, 0, G_TAG);
   	}

	/* Calcolo di E'. Va fatto solo se la matrice e' stata rifattorizzata */
   	if (bNewMatrix) {
     		for (int i = 0; i < iIntVecDim; i++) {
   			pLocalSM->pdSetResVec(pMH->GetECol(i));
			pLocalSM->pdSetSolVec(pMH->GetEColSol(i));
       			/* fa solo la back substitution perche'
			 * la fattorizzazione e' gia' stata lanciata */
			pLocalSM->Solve();
     		}
	}

   	pLocalSM->pdSetResVec(prVH->pdGetVec());
	pLocalSM->pdSetSolVec(pSolrVH->pdGetVec());

	/* on master node... */
	if (MyRank == 0) {
#ifdef USE_MPI_WAIT
		MPI::Request::Waitall(SolvCommSize - 1, &pGRReq[1]);
#else /* ! USE_MPI_WAIT */
		while (true) {
			/* testing remote only... */
      			if (MPI::Request::Testall(SolvCommSize - 1, &pGRReq[1])) 
			{
#ifdef MPI_PROFILING
        			MPE_Log_event(20, 0, "end");
        			for (int i = 1; i < SolvCommSize; i++){
          				MPE_Log_receive(i, G_TAG, pRecvDim[i]);
        			}
#endif /* MPI_PROFILING */
  				break;
      			}
			MYSLEEP(1500);
    		}
#endif /* ! USE_MPI_WAIT */

  		/* Assemblaggio del vettore delle incognite g 
		 * sulla macchina master */

#ifdef DEBUG_SCHUR
		std::cout << "# Solve 19:";
		for (int i = 0; i < pDispl[SolvCommSize]; i++) {
			std::cout << " " << pBuffer[i];
		}
		std::cout << std::endl;
#endif /* DEBUG_SCHUR */

    		/* Assembla pSchVH */
		pSchVH->Reset();
		if (iPrbmBlocks == 1) {
    			for (int i = 0; i < pDispl[SolvCommSize]; i++) {
     	 			pSchVH->IncCoef(pSchGlbToLoc[pDofsRecvdList[i]], pBuffer[i]);
			}

		} else {
    			for (int i = 0; i < pDispl[SolvCommSize]; i++) {
				int j = i%iBlkSize;
				int blk = i/iBlkSize;

     	 			pSchVH->IncCoef(pSchGlbToLoc[pDofsRecvdList[j] + blk*iBlkSize], pBuffer[i]);
			}
    		}

	/* on remote nodes... */
   	} else {
#ifdef USE_MPI_WAIT
     		pGSReq[0].Wait();
#else /* ! USE_MPI_WAIT */
 	  	/* verifica completamento ricezioni e trasmissione vettore g */
   		while (true) {
     			if (pGSReq[0].Test()) {
#ifdef MPI_PROFILING
      				MPE_Log_event(14, 0, "end");
#endif /* MPI_PROFILING */
				break;
     			}
       			MYSLEEP(150);
		}
#endif /* ! USE_MPI_WAIT */
  	}

  	/* assembla le matrici di schur locali se e' stata
	 * appena assemblata una nuova matrice di partenza */
  	if (bNewMatrix) {
    		AssSchur();
    		bNewMatrix = false;
  	}

	/* on master node... */
  	if (MyRank == 0) {
  		/* risoluzione schur pb sul master */

#ifdef MPI_PROFILING
    		MPE_Log_event(35, 0, "start");
#endif /* MPI_PROFILING */

    		pInterSM->Solve();

#ifdef MPI_PROFILING
    		MPE_Log_event(36, 0, "end");
#endif /* MPI_PROFILING */

	  	/* invia la soluzione calcolata */
#ifdef MPI_PROFILING
    		MPE_Log_event(13, 0, "start");
#endif /* MPI_PROFILING */
		doublereal *pd = pSolSchVH->pdGetVec();
    		for (int i = 1; i < SolvCommSize; i++) {
			/* sending to remote only... */
      			pGRReq[i] = SolvComm.Isend(pd, 1,
					*ppNewTypes[i], i, G_TAG);
#ifdef MPI_PROFILING
      			MPE_Log_send(i, G_TAG, pRecvDim[i]);
#endif /* MPI_PROFILING */
    		}

		/* emulate local transmission... */
		memmove(pgVH->pdGetVec(), pd, sizeof(double)*pRecvDim[0]);

#ifdef USE_MPI_WAIT
      		MPI::Request::Waitall(SolvCommSize - 1, &pGRReq[1]);
#else /* ! USE_MPI_WAIT */
  		/* Verifica completamento trasmissioni */
    		while (true) {
			/* testing remote only... */
      			if (MPI::Request::Testall(SolvCommSize - 1, &pGRReq[1]))
			{
#ifdef MPI_PROFILING
        			MPE_Log_event(14, 0, "end");
#endif /* MPI_PROFILING */
				break;
      			}
			MYSLEEP(150);
    		}
#endif /* ! USE_MPI_WAIT */

	/* on other nodes... */	
  	} else {
#ifdef MPI_PROFILING
		MPE_Log_event(19, 0, "start");
#endif /* MPI_PROFILING */

  		pGSReq[0] = SolvComm.Irecv(pgVH->pdGetVec(), iIntVecDim,
				MPI::DOUBLE, 0, G_TAG);

#ifdef USE_MPI_WAIT
    		pGSReq[0].Wait();
#else /* ! USE_MPI_WAIT */
 	 	while (true) {
    			if (pGSReq[0].Test()) {
#ifdef MPI_PROFILING
    				MPE_Log_event(20, 0, "end");
    				MPE_Log_receive(0, G_TAG, iIntVecDim);
#endif /* MPI_PROFILING */
      				break;
    			}
      			MYSLEEP(1500);
		}
#endif /* ! USE_MPI_WAIT */
  	}

  	/* calcolo la soluzione corretta per x = Solr - E'*g */
	pMH->CompNewf(*pSolrVH, *pgVH);

#ifdef DEBUG_SCHUR
	std::cout << "# Solve 20a:";
	for (int i = 1; i <= iIntVecDim; i++) {
		std::cout << " " << pgVH->operator()(i);
	}
	std::cout << std::endl;
	std::cout << "# Solve 20b:";
	for (int i = 1; i <= iLocVecDim; i++) {
		std::cout << " " << pSolrVH->operator()(i);
	}
	std::cout << std::endl;
#endif /* DEBUG_SCHUR */
}

void
SchurSolutionManager::AssSchur(void)
{
#ifdef MPI_PROFILING
   	MPE_Log_event(41, 0, "start");
#endif /* MPI_PROFILING */

  	DEBUGCOUT("Entering SchurSolutionManager::AssSchur()" << endl);

#ifdef DEBUG
  	IsValid();
#endif /* DEBUG */

	pMH->CompLocSchur();

	/* on local node... */
  	if (MyRank == 0) {
#ifdef MPI_PROFILING
   		MPE_Log_event(19, 0, "start");
#endif /* MPI_PROFILING */

		/* remote nodes contribution... */
  		int iOffset = pRecvDim[0]*pRecvDim[0];
    		for (int i = 1; i < SolvCommSize; i++) {
			int	iShift = pRecvDim[i]*pRecvDim[i];

      			pGRReq[i] = SolvComm.Irecv(&pBuffer[iOffset],
					iShift, MPI::DOUBLE, i, S_TAG);

      			iOffset += iShift;
    		}

		/* local node contribution... */
		(void)memmove(&pBuffer[pDispl[0]], pdCM,
			      sizeof(double)*pRecvDim[0]*pRecvDim[0]);

#ifdef USE_MPI_WAIT
      		MPI::Request::Waitall(SolvCommSize - 1, &pGRReq[1]);
#else /* ! USE_MPI_WAIT */
    		while (true) {
      			if (MPI::Request::Testall(SolvCommSize - 1, &pGRReq[1]))
			{
#ifdef MPI_PROFILING
        			MPE_Log_event(20, 0, "end");
        			for (int i = 1; i < SolvCommSize; i++){
          				MPE_Log_receive(i, S_TAG, pRecvDim[i]*pRecvDim[i]);
        			}
#endif /* MPI_PROFILING */
				break;
      			}
			MYSLEEP(1500);
    		}
#endif /* ! USE_MPI_WAIT */

    		/* Assembla Schur matrix */
		pSchMH->Reset();
    		doublereal *pbmd = &pBuffer[-pDispl[0]];
    		for (int i = 0; i < SolvCommSize; i++) {
			if (iPrbmBlocks == 1) {
      				for (int j = pDispl[i]; j < pDispl[i + 1]; j++) {
					int idxC = pSchGlbToLoc[pDofsRecvdList[j]];

					for (int k = pDispl[i]; k < pDispl[i + 1]; k++) {
						int idxR = pSchGlbToLoc[pDofsRecvdList[k]];

						pSchMH->IncCoef(idxR, idxC, pbmd[k]);
					}

      					pbmd += pRecvDim[i];
				}

			} else {
      				for (int j = pDispl[i]; j < pDispl[i + 1]; j++) {

					int ic = j%iBlkSize;
					int bc = j/iBlkSize;
					int idxC = pSchGlbToLoc[pDofsRecvdList[ic] + bc*iBlkSize];
	
					for (int k = pDispl[i]; k < pDispl[i + 1]; k++) {
						int ir = k%iBlkSize;
						int br = k/iBlkSize;
						int idxR = pSchGlbToLoc[pDofsRecvdList[ir] + br*iBlkSize];

						pSchMH->IncCoef(idxR, idxC, pbmd[k]);
					}

      					pbmd += pRecvDim[i];
				}
      			}
	
			/* NOTE: this is required to allow the innermost loop
			 * on k to be between pDispl[i] and pDispl[i + 1];
			 * as a consequence, pbmd must be decreased
			 * by pDispl[i] at each i.
			 * Note that, in general, pDispl[i + 1] - pDispl[i]
			 * is equal to pRecvDim[i] */
			pbmd -= pDispl[i + 1] - pDispl[i];
    		}

	/* on other nodes... */
	} else {
	  	/* Trasmette le Schur locali */
#ifdef 	MPI_PROFILING
		MPE_Log_event(13, 0, "start");
		MPE_Log_send(0, S_TAG, iIntVecDim*iIntVecDim);
#endif /* MPI_PROFILING */

  		pGSReq[0] = SolvComm.Isend(pdCM, iIntVecDim*iIntVecDim,
				MPI::DOUBLE, 0, S_TAG);

#ifdef USE_MPI_WAIT
    		pGSReq[0].Wait();
#else /* !USE_MPI_WAIT */
	  	/* verifica completamento ricezioni e trasmissione */
  		while (true) {
    			if (pGSReq[0].Test()) {
#ifdef MPI_PROFILING
      				MPE_Log_event(14, 0, "end");
#endif /* MPI_PROFILING */
      				break;
    			}
      			MYSLEEP(150);
  		}
#endif /* !USE_MPI_WAIT */
	}

#ifdef MPI_PROFILING
	MPE_Log_event(42, 0, "end");
#endif /* MPI_PROFILING */
}

void
SchurSolutionManager::InitializeComm(void)
{
	DEBUGCOUT("Entering SchurSolutionManager::InitializeComm()" << endl);

	/* il master riceve informazioni sulle dimensioni 
	 * dei vettori scambiati */
	if (MyRank == 0) {
		SAFENEWARR(pRecvDim, int, SolvCommSize);
		SAFENEWARR(pDispl, int, SolvCommSize + 1);
		pDispl[0] = 0;
	}

	/* trasmissione dimensioni interfacce locali */
	SolvComm.Gather(&iIntVecDim, 1, MPI::INT, pRecvDim, 1, MPI::INT, 0);
	if (MyRank == 0) {
		for (int i = 0; i < SolvCommSize; i++){
			pDispl[i + 1] = pDispl[i] + pRecvDim[i];
		}
	}

	if (MyRank == 0 && pDispl[SolvCommSize] == 0) {
		silent_cerr("SchurSolutionManager::InitializeComm(): "
				"empty problem" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (MyRank == 0) {
		SAFENEWARR(pSchurDofs, integer, pDispl[SolvCommSize]);
		SAFENEWARR(pDofsRecvdList, integer, pDispl[SolvCommSize]);
	}

	SolvComm.Gatherv(pIntDofs, iIntVecDim, MPI::INT,
			pDofsRecvdList, pRecvDim, pDispl, MPI::INT, 0);

	if (MyRank == 0) {
		for (int i = 0; i < pDispl[SolvCommSize]; i++) {
			pSchurDofs[i] = pDofsRecvdList[i];
		}

		/* ordina gli indici dei residui locali */
		std::sort(pSchurDofs, pSchurDofs + pDispl[SolvCommSize]);

		/* elimina le ripetizioni */
		integer* p = std::unique(pSchurDofs,
				pSchurDofs + pDispl[SolvCommSize]);

		/* dimensione effettiva del residuo locale */
		iSchurIntDim = p - pSchurDofs;
	}

	/* Vettore di trasformazione locale-globale */
	SAFENEWARR(pGlbToLoc, integer, iPrbmSize*iPrbmBlocks + 1);
	SAFENEWARR(pSchGlbToLoc, integer, iPrbmSize*iPrbmBlocks + 1);
	for (int i = 0; i < iPrbmSize*iPrbmBlocks + 1; i++) {
		pGlbToLoc[i] = 0;
		pSchGlbToLoc[i] = 0;
	}

	for (int i = 0; i < iLocVecDim; i++) {
		for (int j = 0; j < iPrbmBlocks; j++) {
			pGlbToLoc[pLocDofs[i] + j*iBlkSize] = i + j*iBlkSize + 1;
		}
	}

	/* per distinguerli gli indici dei nodi di interfaccia sono negativi */
	for (int i = 0; i < iIntVecDim; i++) {
		for (int j = 0; j < iPrbmBlocks; j++) {
			pGlbToLoc[pIntDofs[i] + j*iBlkSize] = -(i  + j*iBlkSize + 1);
		}
	}

	/* Global to local per la matrice di schur */
	if (MyRank == 0) {
		for (int i = 0; i < iSchurIntDim; i++) {
			for (int j = 0; j < iPrbmBlocks; j++) {
				pSchGlbToLoc[pSchurDofs[i]+ j*iBlkSize] = i + j*iBlkSize + 1;
			}
		}
	}

	iLocVecDim *= iPrbmBlocks;
	iIntVecDim *= iPrbmBlocks;
	iPrbmSize *= iPrbmBlocks;
	iSchurIntDim *= iPrbmBlocks;

	/* creo i buffer per la ricezione e trasmissione dei messaggi */
	if (MyRank == 0) {
		/* i  messaggi + grandi sono legati alla ricezione
		 * delle matrici di schur */
		integer iTmpTot = 0;
		for (int i = 0; i < SolvCommSize; i++) {
			pRecvDim[i] = pRecvDim[i]*iPrbmBlocks;
			iTmpTot += pRecvDim[i]*pRecvDim[i];
		}

		for(int i = 0; i < SolvCommSize; i++){
			pDispl[i + 1] = pDispl[i] + pRecvDim[i];
		}

		/* buffer di ricezione */
		SAFENEWARR(pBuffer, doublereal, iTmpTot);

	} else {
		/* il messaggi + grandi sono le ricezioni 
		 * dei valori di interfaccia */
		SAFENEWARR(pBuffer, doublereal, iIntVecDim);
	}

	if (MyRank == 0){
		SAFENEWARRNOFILL(pGSReq, MPI::Request, SolvCommSize);
		SAFENEWARRNOFILL(pGRReq, MPI::Request, SolvCommSize);

	} else {
		SAFENEWARRNOFILL(pGSReq, MPI::Request, 1);
		SAFENEWARRNOFILL(pGRReq, MPI::Request, 1);
	}
}

/* sposta il puntatore al vettore del residuo */
doublereal *
SchurSolutionManager::pdSetResVec(doublereal* pRes)
{
	silent_cerr("SchurSolutionManager::pdSetResVec(): "
		"you should not be here!! "
		"Aborting..." << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

/* sposta il puntatore al vettore del residuo */
doublereal *
SchurSolutionManager::pdSetSolVec(doublereal* pSol)
{
	silent_cerr("SchurSolutionManager::pdSetSolVec(): "
		"you should not be here!! "
		"Aborting..." << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

/* Rende disponibile l'handler per la matrice */
MatrixHandler*
SchurSolutionManager::pMatHdl(void) const
{
	ASSERT(pMH != NULL);
	return pMH;
}

/* Rende disponibile l'handler per il termine noto */
VectorHandler*
SchurSolutionManager::pResHdl(void) const
{
	ASSERT(pRVH != NULL);
	return pRVH;
}

/* Rende disponibile l'handler per la soluzione */
VectorHandler*
SchurSolutionManager::pSolHdl(void) const
{
	ASSERT(pSolVH != NULL);
	return pSolVH;
}

void
SchurSolutionManager::StartExchIntRes(void)
{
  	DEBUGCOUT("Entering SchurSolutionManager::StartExchIntRes()" << endl);

    	/* Inizializza Trasmissione di g */

#ifdef MPI_PROFILING
   	MPE_Log_event(13, 0, "start");
    	MPE_Log_send(0, G_TAG, iIntVecDim);
#endif /* MPI_PROFILING */

    	if (MyRank == 0) {
#ifdef MPI_PROFILING
      		MPE_Log_event(19, 0, "start");
#endif /* MPI_PROFILING */

		/* collect remote contributions... */
      		for (int i = 1; i < SolvCommSize; i++) {
			pGRReq[i] = SolvComm.Irecv(&pBuffer[pDispl[i]],
					pRecvDim[i], MPI::DOUBLE, i, G_TAG);
      		}

		/* set up local contribution... */
		(void)memmove(&pBuffer[pDispl[0]], pgVH->pdGetVec(),
			      sizeof(double)*pRecvDim[0]);

    	} else {
    		pGSReq[0] = SolvComm.Isend(pgVH->pdGetVec(), iIntVecDim,
				MPI::DOUBLE, 0, G_TAG);
	}
}

void
SchurSolutionManager::ComplExchIntRes(doublereal& dRes,
		const NonlinearSolverTest* t)
{
  	DEBUGCOUT("Entering SchurSolutionManager::ComplExchIntRes()" << endl);

	/* NOTE: right now, all we transmit is the partial result
	 * of the test, as computed by the caller of this function;
	 * the master node is in charge of computing the contribution
	 * of the interface residual */
	const size_t DBLMSGSIZE = 1;
	doublereal d[DBLMSGSIZE];
	d[0] = dRes;
	
    	if (MyRank == 0) {
#ifdef USE_MPI_WAIT
      		MPI::Request::Waitall(SolvCommSize - 1, &pGRReq[1]);
#else /* ! USE_MPI_WAIT */
      		while (true) {
      			if (MPI::Request::Testall(SolvCommSize - 1, &pGRReq[1]))
			{
#ifdef MPI_PROFILING
  				MPE_Log_event(20, 0, "end");
  				for (int i = 1; i < SolvCommSize; i++){
    					MPE_Log_receive(i, G_TAG, pRecvDim[i]);
	  			}
#endif /* MPI_PROFILING */
  				break;
			}
  			MYSLEEP(1500);
		}
#endif /* ! USE_MPI_WAIT */

		/* Assembla pSVH */
      		pSchVH->Reset();
		if (iPrbmBlocks == 1) {
      			for (int i = 0; i < pDispl[SolvCommSize]; i++) {
				pSchVH->IncCoef(pSchGlbToLoc[pDofsRecvdList[i]], pBuffer[i]);
			}

		} else {
      			for (int i = 0; i < pDispl[SolvCommSize]; i++) {
				int j = i%iBlkSize;
				int blk = i/iBlkSize;

				pSchVH->IncCoef(pSchGlbToLoc[pDofsRecvdList[j] + blk*iBlkSize],
						pBuffer[i]);
			}
      		}

		/* interface contribution to error */
		for (int iCntp1 = 1; iCntp1 <= iSchurIntDim; iCntp1++) {
			t->TestOne(d[0], *pSchVH, iCntp1);
      		}

      		for (int i = 1; i < SolvCommSize; i++){
			pGRReq[i] = SolvComm.Irecv(&pBuffer[DBLMSGSIZE*i],
					DBLMSGSIZE, MPI::DOUBLE,
					i, G_TAG + 100);
      		}

#ifdef USE_MPI_WAIT
      		MPI::Request::Waitall(SolvCommSize - 1, &pGRReq[1]);
#else /* ! USE_MPI_WAIT */
      		while (true) {
      			if (MPI::Request::Testall(SolvCommSize - 1, &pGRReq[1]))
			{
				break;
			}
  			MYSLEEP(1500);
		}
#endif /* ! USE_MPI_WAIT */

		for (int i = 1; i < SolvCommSize; i++) {
			t->TestMerge(d[0], pBuffer[DBLMSGSIZE*i]);
      		}

      		for (int i = 1; i < SolvCommSize; i++){
			SolvComm.Send(d, DBLMSGSIZE, MPI::DOUBLE, i, G_TAG);
      		}

    	} else {
		/* send the residual test value */
		SolvComm.Send(d, DBLMSGSIZE, MPI::DOUBLE, 0, G_TAG + 100);

#ifdef USE_MPI_WAIT
      		pGSReq[0].Wait();
#else /* ! USE_MPI_WAIT */
   	 	/* verifica completamento ricezioni e trasmissione */
    		while (true) {
      			if (pGSReq[0].Test()) {
#ifdef MPI_PROFILING
				MPE_Log_event(14, 0, "end");
#endif /* MPI_PROFILING */
				break;
      			}
			MYSLEEP(150);
    		}
#endif /* ! USE_MPI_WAIT */

		/* wait for the cumulative residual test */
		pGRReq[0] = SolvComm.Irecv(d, DBLMSGSIZE, MPI::DOUBLE, 
				0, G_TAG);

#ifdef USE_MPI_WAIT
      		pGRReq[0].Wait();
#else /* ! USE_MPI_WAIT */
   	 	while (true) {
      			if (pGRReq[0].Test()) {
				break;
      			}
			MYSLEEP(150);
   		}
#endif /* ! USE_MPI_WAIT */
	}

	dRes = d[0];
}

void
SchurSolutionManager::StartExchIntSol(void)
{
  	DEBUGCOUT("Entering SchurSolutionManager::StartExchIntSol()" << endl);

	/* FIXME: the solution interface should have alread been exchanged
	 * during the solution process, so we don't start anything... */
	return;
}

void
SchurSolutionManager::ComplExchIntSol(doublereal& dSol,
		const NonlinearSolverTest* t)
{
  	DEBUGCOUT("Entering SchurSolutionManager::ComplExchIntSol()" << endl);

	/* right now, all we transmit is the partial result of the test,
	 * as computed by the caller of this function */
	const size_t DBLMSGSIZE = 1;
	doublereal d[DBLMSGSIZE];
	d[0] = dSol;
	
    	if (MyRank == 0) {
		/* Note: the interface contribution should have already
		 * been transmitted during Solve(), and stored in pgVH. */

		/* interface contribution to error */
		for (int iCntp1 = 1; iCntp1 <= iSchurIntDim; iCntp1++) {
			t->TestOne(d[0], *pgVH, iCntp1);
      		}

      		for (int i = 1; i < SolvCommSize; i++){
			pGRReq[i] = SolvComm.Irecv(&pBuffer[DBLMSGSIZE*i],
					DBLMSGSIZE, MPI::DOUBLE,
					i, G_TAG + 100);
      		}

#ifdef USE_MPI_WAIT
      		MPI::Request::Waitall(SolvCommSize - 1, &pGRReq[1]);
#else /* ! USE_MPI_WAIT */
      		while (true) {
      			if (MPI::Request::Testall(SolvCommSize - 1, &pGRReq[1]))
			{
				break;
			}
  			MYSLEEP(1500);
		}
#endif /* ! USE_MPI_WAIT */

		for (int i = 1; i < SolvCommSize; i++) {
			t->TestMerge(d[0], pBuffer[DBLMSGSIZE*i]);
      		}

      		for (int i = 1; i < SolvCommSize; i++){
			SolvComm.Send(d, DBLMSGSIZE, MPI::DOUBLE, i, G_TAG);
      		}

    	} else {
		/* send the residual test value */
		SolvComm.Send(d, DBLMSGSIZE, MPI::DOUBLE, 0, G_TAG + 100);

#ifdef USE_MPI_WAIT
      		pGSReq[0].Wait();
#else /* ! USE_MPI_WAIT */
   	 	/* verifica completamento ricezioni e trasmissione */
    		while (true) {
      			if (pGSReq[0].Test()) {
#ifdef MPI_PROFILING
				MPE_Log_event(14, 0, "end");
#endif /* MPI_PROFILING */
				break;
      			}
			MYSLEEP(150);
    		}
#endif /* ! USE_MPI_WAIT */

		/* wait for the cumulative residual test */
		pGRReq[0] = SolvComm.Irecv(d, DBLMSGSIZE, MPI::DOUBLE, 
				0, G_TAG);

#ifdef USE_MPI_WAIT
      		pGRReq[0].Wait();
#else /* ! USE_MPI_WAIT */
   	 	while (true) {
      			if (pGRReq[0].Test()) {
				break;
      			}
			MYSLEEP(150);
   		}
#endif /* ! USE_MPI_WAIT */
	}

	dSol = d[0];
}

/* SchurSolutionManager - End */

#endif /* USE_MPI */

