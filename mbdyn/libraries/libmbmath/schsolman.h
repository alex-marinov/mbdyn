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

/* Solutore Schur */
 
/* 
 * Copyright 1999-2000 Giuseppe Quaranta <giuquaranta@tiscalinet.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */

#ifndef SCHSOLMAN_H
#define SCHSOLMAN_H


/* per il debugging */
#include "myassert.h"
#include "mynewmem.h"
#include "except.h"
  
#include "solman.h"
#include "harwrap.h"


#include "mschwrap.h"


#include "mpi++.h"




void InitializeList(int* list, integer dim, integer  value);

void InitializeList(doublereal* list, integer dim, doublereal  value);  
/* SchurMatrixHandler - Begin */

const int G_TAG = 100;
const int S_TAG = 200; 

class SchurMatrixHandler : public MatrixHandler {

 public:
  class ErrGeneric{};

 private: 
  integer LSize, ISize; /* dimensioni locali, interfacce */
  MatrixHandler* pB;
  doublereal* pE;
  doublereal* pF;
  doublereal* pC;
  integer* pGTL;
  const doublereal dZero;
  
 public: 
  SchurMatrixHandler(int LocSize, int IntSize,
		     MatrixHandler* pMatB,
		     doublereal* pMatE,
		     doublereal* pMatF,
		     doublereal* pMatC,
		     integer* pGlobToLoc);
  
  ~SchurMatrixHandler(void); 
  
  /* Usata per il debug */
  void IsValid(void) const;
  
  /* Resetta la matrice */
  inline void Init(const doublereal& dResetVal);
  
  /* Inserisce un coefficiente */
  inline flag fPutCoef(integer iRow, integer iCol, const doublereal& dCoef);
  
  /* Incrementa un coefficiente - se non esiste lo crea */
  inline flag fIncCoef(integer iRow, integer iCol, const doublereal& dCoef);
  
  /* Incrementa un coefficiente - se non esiste lo crea */
  inline flag fDecCoef(integer iRow, integer iCol, const doublereal& dCoef);
  
  /* Restituisce un coefficiente - zero se non e' definito */
  inline const doublereal& dGetCoef(integer iRow, integer iCol) const;
  
  /* dimensioni */
  integer iGetNumRows(void) const {
    return LSize+ISize;
  };
  
  integer iGetNumCols(void) const {
    return LSize+ISize;
  };   
  
};

inline void SchurMatrixHandler::Init(const doublereal& dResetVal)
{
#ifdef DEBUG
  IsValid();
#endif
  pB->Init(dResetVal);
  InitializeList(pE,LSize * ISize, dResetVal);
  InitializeList(pF,LSize * ISize, dResetVal);
  InitializeList(pC,ISize * ISize, dResetVal);
}

inline flag SchurMatrixHandler::fPutCoef(integer iRow, 
					 integer iCol, 
					 const doublereal& dCoef) 
{
#ifdef DEBUG
  IsValid();
  if ((pGTL[iRow] == 0) ||(pGTL[iCol] == 0))  {
    cerr << "fPutCoef " << "Process: " << MPI::COMM_WORLD.Get_rank()
      << " Warning, you are trying to operate on a non local value !!!"
      " (Matrix Handler) "<< iRow << " " << iCol  << endl;
    return flag(0);
  }
#endif
  if (pGTL[iRow] > 0) { 
    if (pGTL[iCol] > 0) { 
      pB->fPutCoef(pGTL[iRow], pGTL[iCol], dCoef); 
      // pB[pGTL[iRow]-1 + (pGTL[iCol]-1) *LSize] = dCoef;
    } 
    else {
      pE[pGTL[iRow]-1 + (-pGTL[iCol]-1)*LSize] = dCoef;
    }
  } 
  else {
    if (pGTL[iCol] > 0) {
      pF[-pGTL[iRow]-1 + (pGTL[iCol]-1)*ISize] = dCoef;
    } 
    else {
      pC[-pGTL[iRow]-1 + (-pGTL[iCol]-1)*ISize] = dCoef;
    }
  }
  return flag(0);
}

inline flag SchurMatrixHandler::fIncCoef(integer iRow, 
					 integer iCol, 
					 const doublereal& dCoef)
{
#ifdef DEBUG
  IsValid();
  if ((pGTL[iRow] == 0) ||(pGTL[iCol] == 0))  {
    cerr << "fIncCoef " << "Process: " << MPI::COMM_WORLD.Get_rank()
      << " Warning, you are trying to operate on a non local value !!!"
      " (Matrix Handler) "<< iRow << " " << iCol  << endl;
    return  flag(0);
  }
#endif
  if (pGTL[iRow] > 0) { 
    if (pGTL[iCol] > 0) { 
      pB->fIncCoef(pGTL[iRow], pGTL[iCol], dCoef);
      //      pB[pGTL[iRow]-1 + (pGTL[iCol]-1) *LSize] += dCoef;
    } 
    else {
      pE[pGTL[iRow]-1 + (-pGTL[iCol]-1)*LSize] += dCoef;
    }
  } 
  else {
    if (pGTL[iCol] > 0) {
      pF[-pGTL[iRow]-1 + (pGTL[iCol]-1)*ISize] += dCoef;
    } 
    else {
      pC[-pGTL[iRow]-1 + (-pGTL[iCol]-1)*ISize] += dCoef;
    }
  }
  return flag(0);
}
 

inline flag SchurMatrixHandler::fDecCoef(integer iRow,
					 integer iCol, 
					 const doublereal& dCoef)
{
#ifdef DEBUG
  IsValid();
  if ((pGTL[iRow] == 0) ||(pGTL[iCol] == 0))  {
    cerr << "fDecCoef " << "Process: " << MPI::COMM_WORLD.Get_rank()
      << " Warning, you are trying to operate on a non local value !!!"
      " (Matrix Handler) "<< iRow << " " << iCol << endl;
    return  flag(0);
  }
#endif
  if (pGTL[iRow] > 0) { 
    if (pGTL[iCol] > 0) { 
      pB->fDecCoef(pGTL[iRow], pGTL[iCol], dCoef);
      //pB[pGTL[iRow]-1 + (pGTL[iCol]-1) *LSize] -= dCoef;
    } 
    else {
      pE[pGTL[iRow]-1 + (-pGTL[iCol]-1)*LSize] -= dCoef;
    }
  } 
  else {
    if (iCol > 0) {
      pF[-pGTL[iRow]-1 + (pGTL[iCol]-1)*ISize] -= dCoef;
    } 
    else {
      pC[-pGTL[iRow]-1 + (-pGTL[iCol]-1)*ISize] -= dCoef;
    }
  }
  return flag(0);
}

inline const doublereal& SchurMatrixHandler::dGetCoef(integer iRow, integer iCol) const
{
#ifdef DEBUG
  IsValid();
  if ((pGTL[iRow] == 0) ||(pGTL[iCol] == 0))  {
    cerr << "dGetCoef " << "Process: " << MPI::COMM_WORLD.Get_rank()
      << " Warning, you are trying to operate on a non local value !!!"
      " (Matrix Handler) "<< iRow << " " << iCol << endl;
    return dZero; 
  }
#endif
  if (pGTL[iRow] > 0) { 
    if (pGTL[iCol] > 0) { 
      return pB->dGetCoef(pGTL[iRow], pGTL[iCol]);
      // return pB[pGTL[iRow]-1 + (pGTL[iCol]-1) *LSize];
    } 
    else {
      return pE[pGTL[iRow]-1 + (-pGTL[iCol]-1)*LSize];
    }
  } 
  else {
    if (pGTL[iCol] > 0) {
      return pF[-pGTL[iRow]-1 + (pGTL[iCol]-1)*ISize];
    } 
    else {
      return pC[-pGTL[iRow]-1 + (-pGTL[iCol]-1)*ISize];
    }
  }
}

/* SchurMatrixHandler - End */

/* SchurVectorHandler - Start */

class SchurVectorHandler : public VectorHandler {

 public:

  class ErrGeneric{};

 private:
  integer LSize, ISize;
  doublereal* pLV;
  doublereal* pIV;
  integer* pGTL;
  
  const doublereal dZero;
  
 public:
  SchurVectorHandler(int LocSize, int IntSize,
		     doublereal* pLocVec,
		     doublereal* pIntVec,
		     integer* pGlobToLoc);
  
  ~SchurVectorHandler(void);
  
  /* Usata per il debug */
  void IsValid(void) const;
                                       
  /* restituisce il puntatore al vettore */
  inline doublereal* pdGetVec(void) const {
    cerr << endl << "You shouldn't have asked for the internal pointer to a SchurVector  " << endl; 
    return pLV;
  }; 
  
  /* restituisce le dimensioni del vettore */
  inline integer iGetSize(void) const {
    return LSize+ISize;
  };  
  
  void Resize(integer iNewSize) {
    cerr << endl << " Why ar you trying to resize a SchurVector ????  " << endl;
  };
  
  /* assegna il dResetVal a tutti gli elementi del vettore */
  void Reset(doublereal dResetVal = 0.) {
    InitializeList(pLV, LSize, dResetVal);
    InitializeList(pIV, ISize, dResetVal);
  };
  
  inline flag fPutCoef(integer iRow, const doublereal& dCoef);
  
  inline flag fIncCoef(integer iRow, const doublereal& dCoef);
  
  inline flag fDecCoef(integer iRow, const doublereal& dCoef);
  
  inline const doublereal& dGetCoef(integer iRow) const;
 
};

inline flag SchurVectorHandler::fPutCoef(integer iRow, const doublereal& dCoef)
{
#ifdef DEBUG
  IsValid();
  if (pGTL[iRow] == 0) {
    cerr << "fIncCoef "<< "Process: " << MPI::COMM_WORLD.Get_rank()
      << " Warning, you are trying to operate on a non local value !!!"
      " (Vector Handler) " << iRow << endl;
    return flag(0); 
  }
#endif
  if (pGTL[iRow] > 0) {
    pLV[pGTL[iRow]-1] = dCoef;
  }
  else {
    pIV[-pGTL[iRow]-1] = dCoef;
  }  
}

inline flag SchurVectorHandler::fIncCoef(integer iRow, const doublereal& dCoef)
{
#ifdef DEBUG
  IsValid();
  if (pGTL[iRow] == 0) {
    cerr <<"fDecCoef " << "Process: " << MPI::COMM_WORLD.Get_rank()
      << " Warning, you are trying to operate on a non local value !!!"
      " (Vector Handler) " << iRow << endl;
    return  flag(0);
  }
#endif
  if (pGTL[iRow] > 0) {
    pLV[pGTL[iRow]-1] += dCoef;
  }
  else {
    pIV[-pGTL[iRow]-1] += dCoef;
  }  
  return flag(0); 
}

inline flag SchurVectorHandler::fDecCoef(integer iRow, const doublereal& dCoef)
{
#ifdef DEBUG
  IsValid();
  if (pGTL[iRow] == 0) {
    cerr << " Warning, you are trying to operate on a non local value !!! (Vector Handler) " << iRow <<  endl;
    return  flag(0);
  }
#endif
  if (pGTL[iRow] > 0) {
    pLV[pGTL[iRow]-1] -= dCoef;
  }
  else {
    pIV[-pGTL[iRow]-1] -= dCoef;
  }  
  return flag(0); 
}

inline const doublereal& SchurVectorHandler::dGetCoef(integer iRow) const

{
#ifdef DEBUG
  IsValid();
  if (pGTL[iRow] == 0) {
    cerr <<"dGetCoef "  << "Process: " << MPI::COMM_WORLD.Get_rank()
      << " Warning, you are trying to operate on a non local value !!!"
      " (Vector Handler) " << iRow << endl;
    return dZero;
  } 
#endif
  if (pGTL[iRow] > 0) {
    return pLV[pGTL[iRow]-1];
  }
  else {
    return pIV[-pGTL[iRow]-1];
  }  
}

/* SchurVectorHandler - End */
  


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
  MPI::Datatype** ppNewTypes;        /* datatype per la trasmissione dei vettori soluzioe delle interfacce; master */

  //  integer SchurWorkSpace;
  integer iWorkSpaceSize;

  /* Vettori di lavoro x l'Harwell solver */
  integer*  piBRow;
  integer*  piBCol;
  //  doublereal*  pdBMat;
  HSMatrixHandler* pBMH;
  
  /* Matrice B */
  doublereal*  pdBMat;

  /* vettore r (residuo locale) */
  doublereal* pdrVec;
  
  HarwellLUSolver* pLU;           /* Solutore sparso Harwell  */
  
  /* Matrice E */
  doublereal* pdEMat;
  
  /* Matrice F */
  doublereal* pdFMat;;
  
  /* Matrice C */
  doublereal* pdCMat;

  /* vettore g (residuo interfaccia) */
  doublereal* pdgVec;

  /* Matrice di Schur (master) */
#ifndef USE_MESCHACH
#error "cannot use SchurSolutionManager without meschach"
#endif /* !USE_MESCHACH */
  MeschachSparseLUSolutionManager* pSchSM;
  MatrixHandler* pSchMH;

 /* vettore Schur */
  //doublereal* pdScVec;
  VectorHandler* pSchVH;

  SchurMatrixHandler* pMH;         /* handler della matrice globale */
  SchurVectorHandler* pRVH;        /* handler del vettore residuo  globale */
  
  MPI::Request* pGSReq;               /* Array di request Send */
  MPI::Request* pGRReq;               /* Array di request Receive */
  
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

