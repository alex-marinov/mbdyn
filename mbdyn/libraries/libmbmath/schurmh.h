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
/* 
 * Copyright (C)1996-2001 
 * Giuseppe Quaranta     <quaranta@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 */
/*Schur Matrix Handler */

#ifndef SCHURMH_H
#define SCHURMH_H
#include <myassert.h>
#include <mynewmem.h>
#include <except.h>
#include <spmapmh.h>
#include <solman.h>


class SchurMatrixHandler : public MatrixHandler {

 public:
  class ErrGeneric{};

 protected: 
  integer LSize, ISize; /* dimensioni locali, interfacce */
  MatrixHandler* pB;
  MyVectorHandler* pE;
  doublereal* pdE;
  SpMapMatrixHandler* pF;
  MyVectorHandler* pC;
  doublereal * pdC;                    
  integer* pGTL;       /* Tabella di conversione Global to Local creata da 
                          SchurSolutionManager 
			  i nodi di interfaccia hanno indice negativo per 
			  permetterne la distizione */  
  
   flag extpdE;
 public: 
  inline SchurMatrixHandler(int LocSize, int IntSize,
                     MatrixHandler* pBM, 
		     integer* pGlobToLoc);
  
  inline SchurMatrixHandler(int LocSize, int IntSize,
                     MatrixHandler* pBM, 
		     integer* pGlobToLoc, doublereal* pdEv);

  virtual inline ~SchurMatrixHandler(void); 
  
  /* Usata per il debug */
  virtual inline void IsValid(void) const;
  
  /* Resetta la matrice */
  virtual inline void Init(const doublereal& dResetVal);
  
  /* Resetta la matrice */
  virtual inline void SchurMatrixHandler::MatEFCInit(const doublereal& dResetVal);

  /* Inserisce un coefficiente */
  virtual inline flag fPutCoef(integer iRow, integer iCol, const doublereal& dCoef);
  
  /* Incrementa un coefficiente - se non esiste lo crea */
  virtual inline flag fIncCoef(integer iRow, integer iCol, const doublereal& dCoef);
  
  /* Incrementa un coefficiente - se non esiste lo crea */
  virtual inline flag fDecCoef(integer iRow, integer iCol, const doublereal& dCoef);
  
  /* Restituisce un coefficiente - zero se non e' definito */
  virtual inline const doublereal& dGetCoef(integer iRow, integer iCol) const;
  
  /* dimensioni */
  virtual integer iGetNumRows(void) const {
    return LSize+ISize;
  };
  
  virtual integer iGetNumCols(void) const {
    return LSize+ISize;
  };   
  
   /* Restituisce l'handler alla matrice B */
   virtual MatrixHandler* GetBMat(void) {
      return pB;
   };
   /* Restituisce il puntatore alla alla matrice C 
      che e' un vettore doublereal contenente le righe in successione */
   virtual doublereal* GetCMat(void) {
      return pdC;
   };

  virtual inline doublereal* GetECol(const integer iCol) const {
  	return &pdE[LSize*iCol];
  };
   
  virtual inline doublereal* GetEColSol(const integer iCol) const {
  	return &pdE[LSize*iCol];
  };   
   
   /* calacola g - F*f e lo pone in g */
   virtual inline VectorHandler& CompNewg(VectorHandler& g, const VectorHandler& f) const;  
   
   /* Calcola la matrice di schur locale C-F*E' e la memorizza in C */
   virtual inline void CompLocSchur(void);
   
   /* Calcola  f - E*g e lo pone in f */
   virtual inline VectorHandler& CompNewf(VectorHandler& f, const VectorHandler& g) const; 

   virtual inline void PrintMatrix(void); 
};

SchurMatrixHandler::SchurMatrixHandler(int LocSize, int IntSize,
				       		MatrixHandler* pBM,
				       		integer* pGlobToLoc)
: LSize(LocSize),
ISize(IntSize),
pB(pBM),
pE(NULL),
pdE(NULL),
pF(NULL),
pC(NULL),
pdC(NULL), 
pGTL(pGlobToLoc),
extpdE(1)
{ 
	SAFENEWARR(pdE, doublereal, LSize*ISize); 
	SAFENEWWITHCONSTRUCTOR(pE,
				MyVectorHandler,
				MyVectorHandler(LSize*ISize, pdE));
	SAFENEWARR(pdC, doublereal, ISize*ISize); 
	SAFENEWWITHCONSTRUCTOR(pC,
				MyVectorHandler,
				MyVectorHandler(ISize*ISize, pdC));
	SAFENEWWITHCONSTRUCTOR(pF,
				SpMapMatrixHandler,
				SpMapMatrixHandler(ISize, LSize));
   
#ifdef DEBUG
  IsValid();
#endif /* DEBUG */
}

SchurMatrixHandler::SchurMatrixHandler(int LocSize, int IntSize,
				       		MatrixHandler* pBM,
				       		integer* pGlobToLoc, 
						doublereal* pdEv)
: LSize(LocSize),
ISize(IntSize),
pB(pBM),
pE(NULL),
pdE(NULL),
pF(NULL),
pC(NULL),
pdC(NULL), 
pGTL(pGlobToLoc),
extpdE(0)
{ 
	SAFENEWWITHCONSTRUCTOR(pE,
				MyVectorHandler,
				MyVectorHandler(LSize*ISize, pdEv));
	SAFENEWARR(pdC, doublereal, ISize*ISize); 
	SAFENEWWITHCONSTRUCTOR(pC,
				MyVectorHandler,
				MyVectorHandler(ISize*ISize, pdC));
	SAFENEWWITHCONSTRUCTOR(pF,
				SpMapMatrixHandler,
				SpMapMatrixHandler(ISize, LSize));
   
#ifdef DEBUG
  IsValid();
#endif /* DEBUG */
}

SchurMatrixHandler::~SchurMatrixHandler(void)
{
	if(pE != NULL) {
		SAFEDELETE(pE);
	}
	if(pC != NULL) {	
		SAFEDELETE(pC);
	}
	if(pF != NULL) {
		SAFEDELETE(pF);
	}
	if (extpdE) {
		if(pdE != NULL) {
			SAFEDELETEARR(pdE);
		}
	}
	if(pdC != NULL) {
		SAFEDELETEARR(pdC);
	}
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


inline void SchurMatrixHandler::MatEFCInit(const doublereal& dResetVal)
{
#ifdef DEBUG
  IsValid();
#endif
  pE->Reset(dResetVal);
  pF->Reset(dResetVal);
  pC->Reset(dResetVal);
}

inline void SchurMatrixHandler::Init(const doublereal& dResetVal)
{
#ifdef DEBUG
  IsValid();
#endif
  pB->Init(dResetVal);
  pE->Reset(dResetVal);
  pF->Reset(dResetVal);
  pC->Reset(dResetVal);
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
    } 
    else {
      pE->fPutCoef(pGTL[iRow]-(pGTL[iCol]+1)*LSize, dCoef);
    }
  } 
  else {
    if (pGTL[iCol] > 0) {
      pF->fPutCoef(-pGTL[iRow], pGTL[iCol], dCoef);
    } 
    else {
      pC->fPutCoef(-pGTL[iRow]-(pGTL[iCol]+1)*ISize, dCoef);
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
    } 
    else {
      pE->fIncCoef(pGTL[iRow]-(pGTL[iCol]+1)*LSize, dCoef);
    }
  } 
  else {
    if (pGTL[iCol] > 0) {
      pF->fIncCoef(-pGTL[iRow], pGTL[iCol], dCoef);
    } 
    else {
      pC->fIncCoef(-pGTL[iRow]-(pGTL[iCol]+1)*ISize, dCoef);
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
    } 
    else {
      pE->fDecCoef(pGTL[iRow]-(pGTL[iCol]+1)*LSize, dCoef);
    }
  } 
  else {
    if (iCol > 0) {
      pF->fDecCoef(-pGTL[iRow], pGTL[iCol], dCoef);
    } 
    else {
      pC->fDecCoef(-pGTL[iRow]-(pGTL[iCol]+1)*ISize, dCoef);
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
    return ::dZero; 
  }
#endif
  if (pGTL[iRow] > 0) { 
    if (pGTL[iCol] > 0) { 
      return pB->dGetCoef(pGTL[iRow], pGTL[iCol]);
    } 
    else {
      return pE->dGetCoef(pGTL[iRow]-(pGTL[iCol]+1)*LSize);
    }
  } 
  else {
    if (pGTL[iCol] > 0) {
      return pF->dGetCoef(-pGTL[iRow], pGTL[iCol]);
    } 
    else {
      return pC->dGetCoef(-pGTL[iRow]-(pGTL[iCol]+1)*ISize);
    }
  }
}




inline VectorHandler& SchurMatrixHandler::CompNewg(VectorHandler& g, const VectorHandler& f) const
{
#ifdef DEBUG 
	ASSERT(f.iGetSize() == LSize);
	ASSERT(g.iGetSize() == ISize);
#endif /* DEBUG */
	pF->MatVecDecMul(g, f);
	return g;
}

  /* Calcola le Schur locali */

inline void SchurMatrixHandler::CompLocSchur(void)
{	
  	for(int j=0; j < ISize; j++) {
    		int iColc = j * ISize;
    		int iCole = j * LSize;
    		for (int k=0; k < LSize; k++) {
      			for( int i=0; i < ISize; i++) {
        			pdC[i + iColc] -=  pF->dGetCoef(i+1,k+1) * pdE[k + iCole];
      			}
    		}
  	}
}

inline VectorHandler& SchurMatrixHandler::CompNewf(VectorHandler& f, const VectorHandler& g) const
{
#ifdef DEBUG
	ASSERT(f.iGetSize() == LSize);
	ASSERT(g.iGetSize() == ISize);
#endif
  	for( int j=0; j < ISize; j++) {  
    		int iColx = j * LSize;
    		for (int i=0; i < LSize; i++) {
      			if (pdE[i + iColx] != 0) {
				f.fDecCoef(i+1, pdE[i + iColx]*g.dGetCoef(j+1));
      			}
    		}
  	}
	return f;
} 

inline void SchurMatrixHandler::PrintMatrix(void) {
	std::cout << "Schur Matrix " << std::endl;
	for (int i=0;i < LSize; i++) {
		for (int j=0;j < LSize; j++) {
 			std::cout << pB->dGetCoef(i+1,j+1) << " "; 
		}
		for (int j=0; j < ISize; j++) {
			std::cout << pdE[i+j*LSize] << " ";
		}	
		std::cout << std::endl;
	}
	for (int i=0;i < ISize; i++) {
		for (int j=0;j < LSize; j++) {
 			std::cout << pF->dGetCoef(i+1,j+1) << " "; 
		}
		for (int j=0; j < ISize; j++) {
			std::cout << pdC[i+j*ISize] << " ";
		}	
		std::cout << std::endl;
	}
}	 
/* SchurMatrixHandler - End */

/* SchurVectorHandler - Start */

class SchurVectorHandler : public VectorHandler {

 public:

  class ErrGeneric{};

 private:
  integer LSize, ISize;
  VectorHandler* pLV;
  VectorHandler* pIV;
  doublereal* pIntVec;
  integer* pGTL;
  
 public:
  SchurVectorHandler(int LocSize, int IntSize,
		     VectorHandler* pLocVec,
		     integer* pGlobToLoc)
  :LSize(LocSize),
   ISize(IntSize),
   pLV(pLocVec),
   pIV(NULL),
   pIntVec(NULL),
   pGTL(pGlobToLoc) 
   {
     SAFENEWARR(pIntVec, doublereal, IntSize); 
     SAFENEWWITHCONSTRUCTOR(pIV,
     			MyVectorHandler,
			MyVectorHandler(IntSize, pIntVec));
   };    	     
		     
   SchurVectorHandler(int LocSize, int IntSize,
		     VectorHandler* pLocV,
		     VectorHandler* pIntV,
		     integer* pGlobToLoc)
  :LSize(LocSize),
   ISize(IntSize),
   pLV(pLocV),
   pIV(pIntV),
   pIntVec(NULL),
   pGTL(pGlobToLoc) 
   {};    	     
		     

  ~SchurVectorHandler(void){
  	if (pIntVec != NULL) {
      		SAFEDELETEARR(pIntVec);
	}
	if (pIV != NULL) {
		SAFEDELETE(pIV);
	}	 
   };
   
  /* Usata per il debug */
  void IsValid(void) const {
    NO_OP;
  };  
                                       
  /* restituisce il puntatore al vettore */
  inline doublereal* pdGetVec(void) const {
    std::cerr << std::endl 
    	<< "You shouldn't have asked for the internal pointer to a SchurVector  " << std::endl; 
    return pLV->pdGetVec();
  }; 
  
  /* restituisce le dimensioni del vettore */
  inline integer iGetSize(void) const {
    return LSize+ISize;
  };  
  
  inline void Resize(integer iNewSize) {
    std::cerr << std::endl << " Why are you trying to resize a SchurVector ????  " << std::endl;
    std::cerr << std::endl << "No Operation Performed!!" << std::endl;
  };
  
  /* assegna il dResetVal a tutti gli elementi del vettore */
  inline void Reset(doublereal dResetVal = 0.) {
    pLV->Reset(dResetVal);
    pIV->Reset(dResetVal); 
  };
  
  inline VectorHandler*  GetIVec(void) {
    return pIV; 
  };
  inline VectorHandler* GetLVec(void) {
    return pLV;
  };

  inline flag fPutCoef(integer iRow, const doublereal& dCoef);
  
  inline flag fIncCoef(integer iRow, const doublereal& dCoef);
  
  inline flag fDecCoef(integer iRow, const doublereal& dCoef);
  
  inline const doublereal& dGetCoef(integer iRow) const;
 
  inline void PrintVector(void);  
};

inline flag SchurVectorHandler::fPutCoef(integer iRow, const doublereal& dCoef)
{
#ifdef DEBUG
  IsValid();
  if (pGTL[iRow] == 0) {
    std::cerr << "fIncCoef "<< "Process: " << MPI::COMM_WORLD.Get_rank()
      << " Warning, you are trying to operate on a non local value !!!"
      " (Vector Handler) " << iRow << std::endl;
    return flag(0); 
  }
#endif
  if (pGTL[iRow] > 0) {
    pLV->fPutCoef(pGTL[iRow], dCoef);
  }
  else {
    pIV->fPutCoef(-pGTL[iRow], dCoef);
  }

  return flag(0);
}

inline flag SchurVectorHandler::fIncCoef(integer iRow, const doublereal& dCoef)
{
#ifdef DEBUG
  IsValid();
  if (pGTL[iRow] == 0) {
    std::cerr <<"fDecCoef " << "Process: " << MPI::COMM_WORLD.Get_rank()
      << " Warning, you are trying to operate on a non local value !!!"
      " (Vector Handler) " << iRow << std::endl;
    return  flag(0);
  }
#endif /* DEBUG */

  if (pGTL[iRow] > 0) {
    pLV->fIncCoef(pGTL[iRow], dCoef);
  }
  else {
    pIV->fIncCoef(-pGTL[iRow], dCoef);
  }  
  return flag(0); 
}

inline flag SchurVectorHandler::fDecCoef(integer iRow, const doublereal& dCoef)
{
#ifdef DEBUG
  IsValid();
  if (pGTL[iRow] == 0) {
    std::cerr << 
  	" Warning, you are trying to operate on a non local value !!! (Vector Handler) " 
	<< iRow <<  std::endl;
    return  flag(0);
  }
#endif /* DEBUG */

  if (pGTL[iRow] > 0) {
    pLV->fDecCoef(pGTL[iRow], dCoef);
  }
  else {
    pIV->fDecCoef(-pGTL[iRow], dCoef);
  }  
  return flag(0); 
}

inline const doublereal& SchurVectorHandler::dGetCoef(integer iRow) const

{
#ifdef DEBUG
  IsValid();
  if (pGTL[iRow] == 0) {
    std::cerr <<"dGetCoef "  << "Process: " << MPI::COMM_WORLD.Get_rank()
      << " Warning, you are trying to operate on a non local value !!!"
      " (Vector Handler) " << iRow << std::endl;
    return dZero;
  } 
#endif /* DEBUG */

  if (pGTL[iRow] > 0) {
    return pLV->dGetCoef(pGTL[iRow]);
  }
  else {
    return pIV->dGetCoef(-pGTL[iRow]);
  }  
}

inline void SchurVectorHandler::PrintVector(void) {
	std::cout << "Schur Vector " << std::endl;
	for (int j=0;j < LSize; j++) {
 		std::cout << pLV->dGetCoef(j+1) << " " << std::endl; 
	}
	for (int j=0; j < ISize; j++) {
 		std::cout << pIV->dGetCoef(j+1) << " " << std::endl; 
	}	
}	 
/* SchurVectorHandler - End */
/* SchurMatrixHandlerUm - begin*/

class SchurMatrixHandlerUm : public SchurMatrixHandler {

 public:
  class ErrGeneric{};

 private:
  doublereal* pdEs;
  MyVectorHandler* pEs;
  integer Eflag;  
 
 public: 
  inline SchurMatrixHandlerUm(int LocSize, int IntSize,
                     MatrixHandler* pBM, 
		     integer* pGlobToLoc);
  
  inline ~SchurMatrixHandlerUm(void); 
   
  /* Resetta le matrici E F e C */
  inline void MatEFCInit(const doublereal& dResetVal);

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
  
  inline doublereal* GetECol(const integer iCol) const;
   
  inline doublereal* GetEColSol(const integer iCol) const;   
     

   /* Calcola la matrice di schur locale C-F*E' e la memorizza in C */
   inline void CompLocSchur(void);
   
   /* Calcola  f - E*g e lo pone in f */
   inline VectorHandler& CompNewf(VectorHandler& f, const VectorHandler& g) const; 

   inline void PrintMatrix(void); 
};

SchurMatrixHandlerUm::SchurMatrixHandlerUm(int LocSize, int IntSize,
				       		MatrixHandler* pBM,
				       		integer* pGlobToLoc)
:SchurMatrixHandler(LocSize, IntSize, pBM, pGlobToLoc, NULL), 
pdEs(NULL),
pEs(NULL),
Eflag(1)
{ 
	
	SAFENEWARR(pdEs, doublereal, LSize*(ISize+1));
	pdE = pdEs+LSize;
	pE->Attach(LSize*ISize,pdE,LSize*ISize);
	SAFENEWWITHCONSTRUCTOR(pEs,
				MyVectorHandler,
				MyVectorHandler((LSize*ISize), pdEs));
}

SchurMatrixHandlerUm::~SchurMatrixHandlerUm(void)
{
	if (pE != NULL) {
		SAFEDELETE(pE);
	}
	if (pdEs != NULL) {
		SAFEDELETEARR(pdEs);
	};
}


inline void SchurMatrixHandlerUm::MatEFCInit(const doublereal& dResetVal)
{
#ifdef DEBUG
  IsValid();
#endif
   Eflag= 1;
  for (int i=0; i < LSize*(ISize+1); i++) {
  	pdEs[i] = dResetVal;
  }	
  pF->Reset(dResetVal);
  pC->Reset(dResetVal);
}

inline void SchurMatrixHandlerUm::Init(const doublereal& dResetVal)
{
#ifdef DEBUG
  IsValid();
#endif
  Eflag= 1;
  pB->Init(dResetVal);
  for (int i=0; i < LSize*(ISize+1); i++) {
  	pdEs[i] = dResetVal;
  }	
  pF->Reset(dResetVal);
  pC->Reset(dResetVal);
}

inline flag SchurMatrixHandlerUm::fPutCoef(integer iRow, 
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
    } else {
    	if (Eflag) {
      		pE->fPutCoef(pGTL[iRow]-(pGTL[iCol]+1)*LSize, dCoef);
	} else {
		pEs->fPutCoef(pGTL[iRow]-(pGTL[iCol]+1)*LSize, dCoef);
    	}
    }
  } else {
    if (pGTL[iCol] > 0) {
      pF->fPutCoef(-pGTL[iRow], pGTL[iCol], dCoef);
    } else {
      pC->fPutCoef(-pGTL[iRow]-(pGTL[iCol]+1)*ISize, dCoef);
    }
  }
  return flag(0);
}

inline flag SchurMatrixHandlerUm::fIncCoef(integer iRow, 
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
    } else {
    	if (Eflag) {
      		pE->fIncCoef(pGTL[iRow]-(pGTL[iCol]+1)*LSize, dCoef);
	} else {
		pEs->fIncCoef(pGTL[iRow]-(pGTL[iCol]+1)*LSize, dCoef);
    	}
    }
  } else {
    if (pGTL[iCol] > 0) {
      pF->fIncCoef(-pGTL[iRow], pGTL[iCol], dCoef);
    } else {
      pC->fIncCoef(-pGTL[iRow]-(pGTL[iCol]+1)*ISize, dCoef);
    }
  }
  return flag(0);
}
 

inline flag SchurMatrixHandlerUm::fDecCoef(integer iRow,
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
    } else {
        if (Eflag) {
      		pE->fDecCoef(pGTL[iRow]-(pGTL[iCol]+1)*LSize, dCoef);
	} else {
		pEs->fDecCoef(pGTL[iRow]-(pGTL[iCol]+1)*LSize, dCoef);
    	}
    }
  } else {
    if (iCol > 0) {
      pF->fDecCoef(-pGTL[iRow], pGTL[iCol], dCoef);
    } else {
      pC->fDecCoef(-pGTL[iRow]-(pGTL[iCol]+1)*ISize, dCoef);
    }
  }
  return flag(0);
}

inline const doublereal& SchurMatrixHandlerUm::dGetCoef(integer iRow, integer iCol) const
{
#ifdef DEBUG
  IsValid();
  if ((pGTL[iRow] == 0) ||(pGTL[iCol] == 0))  {
    cerr << "dGetCoef " << "Process: " << MPI::COMM_WORLD.Get_rank()
      << " Warning, you are trying to operate on a non local value !!!"
      " (Matrix Handler) "<< iRow << " " << iCol << endl;
    return ::dZero; 
  }
#endif
  if (pGTL[iRow] > 0) { 
    if (pGTL[iCol] > 0) { 
      return pB->dGetCoef(pGTL[iRow], pGTL[iCol]);
    } else {
    	if (Eflag) {
      		return pE->dGetCoef(pGTL[iRow]-(pGTL[iCol]+1)*LSize);
	} else {
		return pEs->dGetCoef(pGTL[iRow]-(pGTL[iCol]+1)*LSize);
    	}
    }
  } else {
    if (pGTL[iCol] > 0) {
      return pF->dGetCoef(-pGTL[iRow], pGTL[iCol]);
    } 
    else {
      return pC->dGetCoef(-pGTL[iRow]-(pGTL[iCol]+1)*ISize);
    }
  }
}

inline doublereal* SchurMatrixHandlerUm::GetECol(const integer iCol) const  
{	
 	return &pdE[iCol*LSize];
}

inline doublereal* SchurMatrixHandlerUm::GetEColSol(const integer iCol) const  
{	
 	return &pdEs[iCol*LSize];
}



  /* Calcola le Schur locali */

inline void SchurMatrixHandlerUm::CompLocSchur(void)
{	
  	Eflag=0;
	for(int j=0; j < ISize; j++) {
    		int iColc = j * ISize;
    		int iCole = j * LSize;
    		for (int k=0; k < LSize; k++) {
      			for( int i=0; i < ISize; i++) {
        			pdC[i + iColc] -=  pF->dGetCoef(i+1,k+1) * pdEs[k + iCole];
      			}
    		}
  	}
}

inline VectorHandler& SchurMatrixHandlerUm::CompNewf(VectorHandler& f, const VectorHandler& g) const
{
#ifdef DEBUG
	ASSERT(f.iGetSize() == LSize);
	ASSERT(g.iGetSize() == ISize);
#endif
  	for( int j=0; j < ISize; j++) {  
    		int iColx = j * LSize;
    		for (int i=0; i < LSize; i++) {
      			if (pdEs[i + iColx] != 0) {
				f.fDecCoef(i+1, pdEs[i + iColx]*g.dGetCoef(j+1));
      			}
    		}
  	}
	return f;
} 

inline void SchurMatrixHandlerUm::PrintMatrix(void) {
	std::cout << "Schur Matrix " << std::endl;
	for (int i=0;i < LSize; i++) {
		for (int j=0; j < LSize; j++) {
 			std::cout << pB->dGetCoef(i+1,j+1) << " "; 
		}
		for (int j=0; j < ISize; j++) {
			if (Eflag) {
				std::cout << pdE[i+j*LSize] << " ";
			} else {
				std::cout << pdEs[i+j*LSize] << " ";
			}
		}		
		std::cout << std::endl;
	}
	for (int i=0;i < ISize; i++) {
		for (int j=0;j < LSize; j++) {
 			std::cout << pF->dGetCoef(i+1,j+1) << " "; 
		}
		for (int j=0; j < ISize; j++) {
			std::cout << pdC[i+j*ISize] << " ";
		}	
		std::cout << std::endl;
	}
}	 

/* SchurMatrixHandlerUm - End*/
#endif

