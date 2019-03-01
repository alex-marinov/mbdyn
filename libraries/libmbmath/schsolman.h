/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

/* Solutore Schur */
 
/* 
 * Copyright 1999-2017 Giuseppe Quaranta <quaranta@aero.polimi.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */

#ifndef SCHSOLMAN_H
#define SCHSOLMAN_H

#include "myassert.h"
#include "mynewmem.h"
#include "except.h"
#include "solman.h"
#ifdef USE_MPI
#include "ac/mpi.h"
#endif /* USE_MPI */
#include "schurmh.h"

/* look ahead */
class NonlinearSolverTest;

template <class T>
inline void
InitializeList(T* list, integer dim, T value)
{
	for (int i = 0; i < dim; i++) {
		list[i] = value;
	}
} 

/* SchurMatrixHandler - Begin */

const int G_TAG = 100;
const int S_TAG = 200; 

class SchurSolutionManager : public SolutionManager {
public: 
	class ErrGeneric : public MBDynErrBase {
	public:
		ErrGeneric(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};

   
private:

	/* communicator per lo scambio di messaggi */
#ifdef USE_MPI
	MPI::Intracomm SolvComm; 
#endif /* USE_MPI */
	int MyRank, SolvCommSize;
  
	integer iPrbmSize;          /* dimensioni totali del problema */
	integer iPrbmBlocks;        /* numero di problemi accoppiati, */
                                    /* ciascuno di dimensione iPrbmSize  */
	integer iBlkSize;	    /* size di ciascun blocco */
	integer* pLocDofs;          /* lista dei dof locali */
	integer* pIntDofs;          /* lista dof interfaccia locale */ 
	int iLocVecDim, iIntVecDim; /* dim. vettori dof locali e interfacce */ 
	int* pRecvDim;              /* vettore dim. interfacce locali; master */
	int* pDispl;                /* vettore disp. x comunicazioni; master */
	integer* pDofsRecvdList;    /* lista dof da ricevere sul master;
				     * i dof del processo i iniziano
				     * da pDofsRecvdList[pDispl[i]] */
	integer* pSchurDofs;        /* lista totale dof interfaccia; master */
	int iSchurIntDim;           /* dimensioni matrice di schur; master */
  
	integer* pGlbToLoc;         /* vett. di trasf. indici Glob->loc */
	integer* pSchGlbToLoc;      /* vett. trasf. idx Glob->Schur; master */
  
	int* pBlockLenght;          /* str. di servizio x datatype; master */ 
#ifdef USE_MPI
	MPI::Aint* pTypeDsp;        /* str. di servizio x datatype; master */ 
	MPI::Datatype** ppNewTypes; /* datatype per la trasmissione dei vettori
				     * soluzione delle interfacce; master */
#endif /* USE_MPI */
	doublereal* pBuffer;        /* buffer di ricezione */

	integer iWorkSpaceSize;


	SchurMatrixHandler* pMH;         /* handler della matrice globale */
	SchurVectorHandler* pRVH;        /* handler del vettore residuo globale */
	/* FIXME: do we need a SchurVectorHandler or is a MyVectorHandler enough? */
	VectorHandler* pSolVH;      /* handler del vettore soluzione globale */

	MatrixHandler* pSchMH;      /*handler matrice delle interfacce */
	VectorHandler* pSchVH;  	   /* vettore delle interfacce  */
	VectorHandler* pSolSchVH;     /* soluzione problema delle interfacce */  
	/* Vettori di lavoro x il Matrix Handler */
	MatrixHandler* pBMH;
	doublereal* pdCM;                 /* puntatore necessario per lo scambio delle schur locali */
	VectorHandler* prVH;
	VectorHandler* pgVH;  
	VectorHandler* pSolrVH;

#ifdef USE_MPI
	MPI::Request* pGSReq;               /* Array di request Send */
	MPI::Request* pGRReq;               /* Array di request Receive */
#endif /* USE_MPI */

	SolutionManager* pLocalSM;           /* Solutore sparso locale */
	SolutionManager* pInterSM;          /* Solutore sparso locale */

  
	bool bNewMatrix;
public:
  
	SchurSolutionManager(integer iSize,
			integer iBlocks,
			integer* pLocalDofs,
			int iDim1,
			integer* pInterfDofs,
			int iDim2,
			SolutionManager* pLSM,
			LinSol& ls);
	
	virtual ~SchurSolutionManager(void);

#ifdef DEBUG
	void IsValid(void) const;
#endif /* DEBUG */

	/* Inizializza il gestore delle  matrici */
	void MatrReset(void);

	/* Inizializzatore "speciale" */
	void MatrInitialize(void);

	/* Risolve i blocchi chiamando il solutore */
	void Solve(void);

	/* sposta il puntatore al vettore del residuo */
	doublereal *pdSetResVec(doublereal* pRes);

	/* sposta il puntatore al vettore del residuo */
	doublereal *pdSetSolVec(doublereal* pSol);

	/* Rende disponibile l'handler per la matrice */
	MatrixHandler* pMatHdl(void) const;

	/* Rende disponibile l'handler per il termine noto */
	VectorHandler* pResHdl(void) const;

	/* Rende disponibile l'handler per la soluzione */
	VectorHandler* pSolHdl(void) const;

	void StartExchIntRes(void);
	void ComplExchIntRes(doublereal& d, const NonlinearSolverTest* t);
	void StartExchIntSol(void);
	void ComplExchIntSol(doublereal& d, const NonlinearSolverTest* t);

private:
	void AssSchur(void);

	void InitializeComm(void);
};

#endif /* SCHSOLMAN_H */

