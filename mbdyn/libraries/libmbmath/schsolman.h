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

/* Solutore Schur */
 
/* 
 * Copyright 1999-2000 Giuseppe Quaranta <quaranta@aero.polimi.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */

#ifndef SCHSOLMAN_H
#define SCHSOLMAN_H

#include <myassert.h>
#include <mynewmem.h>
#include <except.h>
#include <solman.h>
#ifdef USE_UMFPACK3
#include 
#else
#error "Sorry, Parallel MBDyn can not be compiled without UMFPACK3!"
#endif  

#include <mpi++.h>

void InitializeList(int* list, integer dim, integer  value);

/* SchurMatrixHandler - Begin */

const int G_TAG = 100;
const int S_TAG = 200; 

class SchurSolutionManager : public SolutionManager {
 
 public: 
  class ErrGeneric {};
   
 private:

  /* communicator per lo scambio di messaggi */
  MPI::Intracomm SolvComm; 
  int MyRank, SolvCommSize;
  
  integer iPrbmSize;                /*dimensioni totali del problema */
  integer* pLocDofs;                /*lista dei dof locali */
  integer* pIntDofs;                /*lista dof interfaccia locale */ 
  int LocVecDim, IntVecDim;         /* dimensioni vettori dof locali e interfacce */ 
  int* pRecvDim;                    /* vettore dimensioni interfaccie locali; master */
  int* pDispl;                      /* vettore displacement x comunicazioni; master */
  integer* pDofsRecvdList;          /* lista dof da ricevere sul master; i dof del processo i
				       iniziano da pDofsRecvdList[pDispl[i]] */
  integer* pSchurDofs;              /* lista totale dof interfaccia; master */
  int iSchurIntDim;                 /* dimensioni matrice di schur; master */
  
  integer* pGlbToLoc;               /* vettore di trsformazione indici Globale->locale */
  integer* pSchGlbToLoc;            /* vettore di trasformazione indici Globale->Schur; master */
  
  int* pBlockLenght;                 /* struttura di servizio x datatype; master */ 
  MPI::Aint* pTypeDsp;               /* struttura di servizio x datatype; master */ 
  doublereal* pBuffer;               /* buffer di ricezione */
  MPI::Datatype** ppNewTypes;        /* datatype per la trasmissione dei vettori soluzione delle interfacce; master */

  integer iWorkSpaceSize;


  SchurMatrixHandler* pMH;         /* handler della matrice globale */
  SchurVectorHandler* pRVH;        /* handler del vettore residuo  globale */

  
  SpMapMatrixHandler* pSchMH;      /*handler matrice delle interfacce */
  VectorHandler* pSchVH;  	   /* vettore delle interfacce  */

  /* Vettori di lavoro x il Matrix Handler */
  SpMapMatrixHandler* pBMH;
  SpMapMatrixHandler* pEMH;
  SpMapMatrixHandler* pFMH;
  doubleral* pCMH;
  MyVectorHandler* prVH;
  MyVectorHandler* pgVH;  

  /* vettore g (residuo interfaccia) */
  doublereal* pdgVec;

  MPI::Request* pGSReq;               /* Array di request Send */
  MPI::Request* pGRReq;               /* Array di request Receive */
  
#ifdef USE_UMFPACK3
 Umfpack3SparseLUSolutionManager* pLocalSM;           /* Solutore sparso locale */
 Umfpack3SparseLUSolutionManager* pInterSM;          /* Solutore sparso locale */
#else
#error "You need to use UMFPACK3 for parallel solving"
#endif /* USE_UMFPACK3 */

  
  flag fNewMatrix;
 public:
  
  SchurSolutionManager(integer iSize,
		     integer* pLocalDofs,
		     int iDim1,
		     integer* pInterfDofs,
		     int iDim2,
		     integer iWorkSpaceSize = 0,
		     const doublereal& dPivotFactor = 1.0);			 
  
  ~SchurSolutionManager(void);
  
  void IsValid(void) const;
  
  /* Inizializza il gestore delle  matrici */
  void MatrInit(const doublereal& dResetVal);

  /* Risolve i blocchi chiamando il solutore */
  void Solve(void);

  /* Rende disponibile l'handler per la matrice */
  SchurMatrixHandler* pMatHdl(void) const {
    ASSERT(pMH != NULL);	
    return pMH;
  };
   
  /* Rende disponibile l'handler per il termine noto */
  SchurVectorHandler* pResHdl(void) const {
    ASSERT(pRVH != NULL);	
    return pRVH;
  };

  /* Rende disponibile l'handler per la soluzione */
  SchurVectorHandler* pSolHdl(void) const {
    ASSERT(pRVH != NULL);	
    return pRVH;
  };

  void StartExchInt(void);
  
  void ComplExchInt(doublereal& dR, doublereal& dXP);
  

 private:
  
  void AssSchur(void);
  
  void InitializeComm(void);


};

#endif

