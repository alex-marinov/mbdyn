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
#include <mbcomm.h>

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
	doublereal* pdC;                    
	integer* pGTL;       /* Tabella di conversione Global to Local 
			      * creata da SchurSolutionManager 
			      * i nodi di interfaccia hanno indice negativo 
			      * per permetterne la distizione */  
  
	bool extpdE;

public: 
 	SchurMatrixHandler(int LocSize, int IntSize,
   			MatrixHandler* pBM, 
   			integer* pGlobToLoc, doublereal* pdEv = NULL);

	virtual ~SchurMatrixHandler(void); 

	/* Usata per il debug */
	virtual void IsValid(void) const;


	virtual integer iGetNumRows(void) const;
	virtual integer iGetNumCols(void) const;
	MatrixHandler* GetBMat(void);
	doublereal* GetCMat(void);
  
	/* Resetta la matrice */
	virtual inline void Init(const doublereal& dResetVal);
  
	/* Resetta la matrice */
	virtual inline void 
	SchurMatrixHandler::MatEFCInit(const doublereal& dResetVal);

	/* Inserisce un coefficiente */
	virtual inline flag fPutCoef(integer iRow, integer iCol,
			const doublereal& dCoef);
  
	/* Incrementa un coefficiente - se non esiste lo crea */
	virtual inline flag fIncCoef(integer iRow, integer iCol,
			const doublereal& dCoef);
  	
	/* Incrementa un coefficiente - se non esiste lo crea */
	virtual inline flag fDecCoef(integer iRow, integer iCol,
			const doublereal& dCoef);
  
	/* Restituisce un coefficiente - zero se non e' definito */
	virtual inline const doublereal&
	dGetCoef(integer iRow, integer iCol) const;
  
	virtual inline doublereal* GetECol(const integer iCol) const;
	virtual inline doublereal* GetEColSol(const integer iCol) const;
 
	/* calacola g - F*f e lo pone in g */
	virtual inline VectorHandler& CompNewg(VectorHandler& g,
			const VectorHandler& f) const;  
   
	/* Calcola la matrice di schur locale C-F*E' e la memorizza in C */
	virtual inline void CompLocSchur(void);
   
	/* Calcola  f - E*g e lo pone in f */
	virtual inline VectorHandler&
	CompNewf(VectorHandler& f, const VectorHandler& g) const; 

	virtual inline void PrintMatrix(void); 
};

inline doublereal*
SchurMatrixHandler::GetECol(const integer iCol) const
{
	return &pdE[LSize*iCol];
}
   
inline doublereal*
SchurMatrixHandler::GetEColSol(const integer iCol) const
{
  	return &pdE[LSize*iCol];
}

inline void
SchurMatrixHandler::MatEFCInit(const doublereal& dResetVal)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	pE->Reset(dResetVal);
	pF->Reset(dResetVal);
	pC->Reset(dResetVal);
}

inline void
SchurMatrixHandler::Init(const doublereal& dResetVal)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	pB->Init(dResetVal);
	pE->Reset(dResetVal);
	pF->Reset(dResetVal);
	pC->Reset(dResetVal);
}

inline flag
SchurMatrixHandler::fPutCoef(integer iRow, integer iCol, 
		const doublereal& dCoef) 
{
#ifdef DEBUG
	IsValid();
	if ((pGTL[iRow] == 0) ||(pGTL[iCol] == 0))  {
		std::cerr << "fPutCoef - Process(" << MBDynComm.Get_rank()
			<< "): warning, MatrixHandler is trying to operate "
			"on a non local value {" << iRow << "," << iCol << "}"
			<< std::endl;
		return flag(0);
	}
#endif /* DEBUG */

	if (pGTL[iRow] > 0) { 
		if (pGTL[iCol] > 0) { 
			pB->fPutCoef(pGTL[iRow], pGTL[iCol], dCoef); 
		} else {
			pE->fPutCoef(pGTL[iRow]-(pGTL[iCol]+1)*LSize, dCoef);
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

inline flag
SchurMatrixHandler::fIncCoef(integer iRow, integer iCol, 
		const doublereal& dCoef)
{
#ifdef DEBUG
	IsValid();
	if ((pGTL[iRow] == 0) ||(pGTL[iCol] == 0))  {
		std::cerr << "fIncCoef - Process(" << MBDynComm.Get_rank()
			<< "): warning, MatrixHandler is trying to operate "
			"on a non local value {" << iRow << "," << iCol << "}"
			<< std::endl;
		return flag(0);
	}
#endif /* DEBUG */

	if (pGTL[iRow] > 0) { 
		if (pGTL[iCol] > 0) { 
			pB->fIncCoef(pGTL[iRow], pGTL[iCol], dCoef);
		} else {
			pE->fIncCoef(pGTL[iRow]-(pGTL[iCol]+1)*LSize, dCoef);
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
 
inline flag
SchurMatrixHandler::fDecCoef(integer iRow, integer iCol, 
		const doublereal& dCoef)
{
#ifdef DEBUG
	IsValid();
	if ((pGTL[iRow] == 0) ||(pGTL[iCol] == 0))  {
		std::cerr << "fDecCoef - Process(" << MBDynComm.Get_rank()
			<< "): warning, MatrixHandler is trying to operate "
			"on a non local value {" << iRow << "," << iCol << "}"
			<< std::endl;
		return flag(0);
	}
#endif /* DEBUG */

	if (pGTL[iRow] > 0) { 
		if (pGTL[iCol] > 0) { 
			pB->fDecCoef(pGTL[iRow], pGTL[iCol], dCoef);
		} else {
			pE->fDecCoef(pGTL[iRow]-(pGTL[iCol]+1)*LSize, dCoef);
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

inline const doublereal&
SchurMatrixHandler::dGetCoef(integer iRow, integer iCol) const
{
#ifdef DEBUG
	IsValid();
	if ((pGTL[iRow] == 0) ||(pGTL[iCol] == 0))  {
		std::cerr << "dGetCoef - Process(" << MBDynComm.Get_rank()
			<< "): warning, MatrixHandler is trying to operate "
			"on a non local value {"<< iRow << "," << iCol << "}"
			<< std::endl;
		return ::dZero; 
	}
#endif /* DEBUG */

	if (pGTL[iRow] > 0) { 
		if (pGTL[iCol] > 0) { 
			return pB->dGetCoef(pGTL[iRow], pGTL[iCol]);
		} else {
			return pE->dGetCoef(pGTL[iRow]-(pGTL[iCol]+1)*LSize);
		}
	} else {
		if (pGTL[iCol] > 0) {
			return pF->dGetCoef(-pGTL[iRow], pGTL[iCol]);
		} else {
			return pC->dGetCoef(-pGTL[iRow]-(pGTL[iCol]+1)*ISize);
		}
	}
}

inline VectorHandler&
SchurMatrixHandler::CompNewg(VectorHandler& g, const VectorHandler& f) const
{
#ifdef DEBUG 
	ASSERT(f.iGetSize() == LSize);
	ASSERT(g.iGetSize() == ISize);
#endif /* DEBUG */

	pF->MatVecDecMul(g, f);
	return g;
}

/* Calcola le Schur locali */
inline void 
SchurMatrixHandler::CompLocSchur(void)
{	
  	for (int j = 0; j < ISize; j++) {
    		int iColc = j * ISize;
    		int iCole = j * LSize;

    		for (int k = 0; k < LSize; k++) {
      			for (int i = 0; i < ISize; i++) {
        			pdC[i + iColc] -= 
					pF->dGetCoef(i+1,k+1) * pdE[k + iCole];
      			}
    		}
  	}
}

inline VectorHandler&
SchurMatrixHandler::CompNewf(VectorHandler& f, const VectorHandler& g) const
{
#ifdef DEBUG
	ASSERT(f.iGetSize() == LSize);
	ASSERT(g.iGetSize() == ISize);
#endif /* DEBUG */

  	for (int j = 0; j < ISize; j++) {  
    		int iColx = j * LSize;
    		for (int i = 0; i < LSize; i++) {
      			if (pdE[i + iColx] != 0) {
				f.fDecCoef(i+1, pdE[i+iColx]*g.dGetCoef(j+1));
      			}
    		}
  	}
	return f;
} 

inline void
SchurMatrixHandler::PrintMatrix(void)
{
	std::cout << "Schur Matrix " << std::endl;
	for (int i = 0;i < LSize; i++) {
		for (int j = 0;j < LSize; j++) {
 			std::cout << pB->dGetCoef(i+1,j+1) << " "; 
		}
		for (int j = 0; j < ISize; j++) {
			std::cout << pdE[i+j*LSize] << " ";
		}
		std::cout << std::endl;
	}
	for (int i = 0;i < ISize; i++) {
		for (int j = 0;j < LSize; j++) {
 			std::cout << pF->dGetCoef(i+1,j+1) << " "; 
		}
		for (int j = 0; j < ISize; j++) {
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
	integer* pGTL;
  
public:
	SchurVectorHandler(int LocSize, int IntSize,
   			VectorHandler* pLocVec,
   			integer* pGlobToLoc);
	SchurVectorHandler(int LocSize, int IntSize,
   			VectorHandler* pLocV,
   			VectorHandler* pIntV,
   			integer* pGlobToLoc);
	~SchurVectorHandler(void);

	/* Usata per il debug */
	void IsValid(void) const;
 
	/* restituisce il puntatore al vettore */
	inline doublereal* pdGetVec(void) const;

	/* restituisce le dimensioni del vettore */
	inline integer iGetSize(void) const;

	inline void Resize(integer iNewSize);

	/* assegna il dResetVal a tutti gli elementi del vettore */
	inline void Reset(doublereal dResetVal = 0.);
  
	inline VectorHandler* GetIVec(void);
	inline VectorHandler* GetLVec(void);

	inline flag fPutCoef(integer iRow, const doublereal& dCoef);
	inline flag fIncCoef(integer iRow, const doublereal& dCoef);
	inline flag fDecCoef(integer iRow, const doublereal& dCoef);
	inline const doublereal& dGetCoef(integer iRow) const;
	inline void PrintVector(void);  
};

/* restituisce il puntatore al vettore */
inline doublereal*
SchurVectorHandler::pdGetVec(void) const
{
	std::cerr << "You shouldn't have asked for "
		"the internal pointer of a SchurVectorHandler"
		<< std::endl; 
	return pLV->pdGetVec();
}

/* restituisce le dimensioni del vettore */
inline integer
SchurVectorHandler::iGetSize(void) const
{
	return LSize+ISize;
}

inline void
SchurVectorHandler::Resize(integer iNewSize)
{
	std::cerr << "Why are you trying to resize a SchurVector ????"
		<< std::endl
		<< "No Operation Performed!!" << std::endl;
}

/* assegna il dResetVal a tutti gli elementi del vettore */
inline void
SchurVectorHandler::Reset(doublereal dResetVal)
{
	pLV->Reset(dResetVal);
	pIV->Reset(dResetVal); 
}
  
inline VectorHandler*
SchurVectorHandler::GetIVec(void)
{
	return pIV; 
}

inline VectorHandler*
SchurVectorHandler::GetLVec(void)
{
	return pLV;
}

inline flag
SchurVectorHandler::fPutCoef(integer iRow, const doublereal& dCoef)
{
#ifdef DEBUG
	IsValid();
	if (pGTL[iRow] == 0) {
		std::cerr << "fPutCoef - Process(" << MBDynComm.Get_rank()
			<< "): warning, VectorHandler is trying to operate "
			"on a non local value {" << iRow << "}"
			<< std::endl;
		return flag(0); 
	}
#endif /* DEBUG */

	if (pGTL[iRow] > 0) {
		pLV->fPutCoef(pGTL[iRow], dCoef);
	} else {
		pIV->fPutCoef(-pGTL[iRow], dCoef);
	}

	return flag(0);
}

inline flag
SchurVectorHandler::fIncCoef(integer iRow, const doublereal& dCoef)
{
#ifdef DEBUG
	IsValid();
	if (pGTL[iRow] == 0) {
		std::cerr <<"fIncCoef - Process(" << MBDynComm.Get_rank()
			<< "): warning, VectorHandler is trying to operate "
			"on a non local value {" << iRow << "}"
			<< std::endl;
		return flag(0);
	}
#endif /* DEBUG */

	if (pGTL[iRow] > 0) {
		pLV->fIncCoef(pGTL[iRow], dCoef);
	} else {
		pIV->fIncCoef(-pGTL[iRow], dCoef);
	}  
	return flag(0); 
}

inline flag
SchurVectorHandler::fDecCoef(integer iRow, const doublereal& dCoef)
{
#ifdef DEBUG
	IsValid();
	if (pGTL[iRow] == 0) {
		std::cerr <<"fDecCoef - Process(" << MBDynComm.Get_rank()
			<< "): warning, VectorHandler is trying to operate "
			"on a non local value {" << iRow << "}"
			<< std::endl;
		return flag(0);
	}
#endif /* DEBUG */

	if (pGTL[iRow] > 0) {
		pLV->fDecCoef(pGTL[iRow], dCoef);
	} else {
		pIV->fDecCoef(-pGTL[iRow], dCoef);
	}
	return flag(0); 
}

inline const doublereal&
SchurVectorHandler::dGetCoef(integer iRow) const
{
#ifdef DEBUG
	IsValid();
	if (pGTL[iRow] == 0) {
		std::cerr << "dGetCoef - Process:(" << MBDynComm.Get_rank()
			<< "): warning, VectorHandler is trying to operate "
			"on a non local value {" << iRow << "}"
			<< std::endl;
		return dZero;
	} 
#endif /* DEBUG */

	if (pGTL[iRow] > 0) {
		return pLV->dGetCoef(pGTL[iRow]);
	} else {
		return pIV->dGetCoef(-pGTL[iRow]);
	}
}

inline void
SchurVectorHandler::PrintVector(void)
{
	std::cout << "Schur Vector " << std::endl;
	for (int j = 0;j < LSize; j++) {
 		std::cout << pLV->dGetCoef(j+1) << " " << std::endl; 
	}
	for (int j = 0; j < ISize; j++) {
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
	bool Eflag;  

public: 
	SchurMatrixHandlerUm(int LocSize, int IntSize,
   			MatrixHandler* pBM, 
   			integer* pGlobToLoc);
	~SchurMatrixHandlerUm(void); 
	
	/* Resetta le matrici E F e C */
	inline void MatEFCInit(const doublereal& dResetVal);

	/* Resetta la matrice */
	inline void Init(const doublereal& dResetVal);
  
	/* Inserisce un coefficiente */
	inline flag fPutCoef(integer iRow, integer iCol,
			const doublereal& dCoef);

	/* Incrementa un coefficiente - se non esiste lo crea */
	inline flag fIncCoef(integer iRow, integer iCol,
			const doublereal& dCoef);
  
	/* Incrementa un coefficiente - se non esiste lo crea */
	inline flag fDecCoef(integer iRow, integer iCol,
			const doublereal& dCoef);
  
	/* Restituisce un coefficiente - zero se non e' definito */
	inline const doublereal& dGetCoef(integer iRow, integer iCol) const;
  
	inline doublereal* GetECol(const integer iCol) const;
   
	inline doublereal* GetEColSol(const integer iCol) const;   

	/* Calcola la matrice di schur locale C-F*E' e la memorizza in C */
	inline void CompLocSchur(void);
   
	/* Calcola  f - E*g e lo pone in f */
	inline VectorHandler&
	CompNewf(VectorHandler& f, const VectorHandler& g) const; 

	inline void PrintMatrix(void); 
};

inline void
SchurMatrixHandlerUm::MatEFCInit(const doublereal& dResetVal)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	Eflag = true;
	for (int i = 0; i < LSize*(ISize+1); i++) {
		pdEs[i] = dResetVal;
	}
	pF->Reset(dResetVal);
	pC->Reset(dResetVal);
}

inline void SchurMatrixHandlerUm::Init(const doublereal& dResetVal)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	Eflag = true;
	pB->Init(dResetVal);
	for (int i = 0; i < LSize*(ISize+1); i++) {
		pdEs[i] = dResetVal;
	}
	pF->Reset(dResetVal);
	pC->Reset(dResetVal);
}

inline flag
SchurMatrixHandlerUm::fPutCoef(integer iRow, integer iCol, 
		const doublereal& dCoef) 
{
#ifdef DEBUG
	IsValid();
	if ((pGTL[iRow] == 0) ||(pGTL[iCol] == 0))  {
		std::cerr << "fPutCoef - Process(" << MBDynComm.Get_rank()
			<< "): warning, MatrixHandler is trying to operate "
			"on a non local value {"<< iRow << "," << iCol << "}"
			<< std::endl;
		return flag(0);
	}
#endif /* DEBUG */

	if (pGTL[iRow] > 0) { 
		if (pGTL[iCol] > 0) { 
			pB->fPutCoef(pGTL[iRow], pGTL[iCol], dCoef); 
		} else {
			if (Eflag) {
				pE->fPutCoef(pGTL[iRow]-(pGTL[iCol]+1)*LSize,
						dCoef);
			} else {
				pEs->fPutCoef(pGTL[iRow]-(pGTL[iCol]+1)*LSize,
						dCoef);
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

inline flag
SchurMatrixHandlerUm::fIncCoef(integer iRow, integer iCol, 
		const doublereal& dCoef)
{
#ifdef DEBUG
	IsValid();
	if ((pGTL[iRow] == 0) ||(pGTL[iCol] == 0))  {
		std::cerr << "fIncCoef - Process(" << MBDynComm.Get_rank()
			<< "): warning, MatrixHandler is trying to operate "
			"on a non local value {"<< iRow << "," << iCol << "}"
			<< std::endl;
		return flag(0);
	}
#endif /* DEBUG */

	if (pGTL[iRow] > 0) { 
		if (pGTL[iCol] > 0) { 
			pB->fIncCoef(pGTL[iRow], pGTL[iCol], dCoef);
		} else {
			if (Eflag) {
				pE->fIncCoef(pGTL[iRow]-(pGTL[iCol]+1)*LSize,
						dCoef);
			} else {
				pEs->fIncCoef(pGTL[iRow]-(pGTL[iCol]+1)*LSize,
						dCoef);
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
 
inline flag
SchurMatrixHandlerUm::fDecCoef(integer iRow, integer iCol, 
		const doublereal& dCoef)
{
#ifdef DEBUG
	IsValid();
	if ((pGTL[iRow] == 0) ||(pGTL[iCol] == 0))  {
		std::cerr << "fDecCoef - Process(" << MBDynComm.Get_rank()
			<< "): warning, MatrixHandler is trying to operate "
			"on a non local value {"<< iRow << "," << iCol << "}"
			<< std::endl;
		return flag(0);
	}
#endif /* DEBUG */

	if (pGTL[iRow] > 0) { 
		if (pGTL[iCol] > 0) { 
			pB->fDecCoef(pGTL[iRow], pGTL[iCol], dCoef);
		} else {
			if (Eflag) {
				pE->fDecCoef(pGTL[iRow]-(pGTL[iCol]+1)*LSize,
						dCoef);
			} else {
				pEs->fDecCoef(pGTL[iRow]-(pGTL[iCol]+1)*LSize,
						dCoef);
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

inline const doublereal&
SchurMatrixHandlerUm::dGetCoef(integer iRow, integer iCol) const
{
#ifdef DEBUG
	IsValid();
	if ((pGTL[iRow] == 0) ||(pGTL[iCol] == 0))  {
		std::cerr << "dGetCoef - Process(" << MBDynComm.Get_rank()
			<< "): warning, MatrixHandler is trying to operate "
			"on a non local value {" << iRow << "," << iCol << "}"
			<< std::endl;
		return ::dZero; 
	}
#endif /* DEBUG */
	
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
		} else {
			return pC->dGetCoef(-pGTL[iRow]-(pGTL[iCol]+1)*ISize);
		}
	}
}

inline doublereal*
SchurMatrixHandlerUm::GetECol(const integer iCol) const  
{	
 	return &pdE[iCol*LSize];
}

inline doublereal*
SchurMatrixHandlerUm::GetEColSol(const integer iCol) const  
{	
 	return &pdEs[iCol*LSize];
}

/* Calcola le Schur locali */
inline void
SchurMatrixHandlerUm::CompLocSchur(void)
{	
  	Eflag = false;
	for (int j = 0; j < ISize; j++) {
    		int iColc = j * ISize;
    		int iCole = j * LSize;
		
    		for (int k = 0; k < LSize; k++) {
      			for (int i = 0; i < ISize; i++) {
        			pdC[i + iColc] -= 
					pF->dGetCoef(i+1, k+1) * pdEs[k+iCole];
      			}
    		}
  	}
}

inline VectorHandler&
SchurMatrixHandlerUm::CompNewf(VectorHandler& f, const VectorHandler& g) const
{
#ifdef DEBUG
	ASSERT(f.iGetSize() == LSize);
	ASSERT(g.iGetSize() == ISize);
#endif /* DEBUG */

  	for (int j = 0; j < ISize; j++) {  
    		int iColx = j * LSize;
		
    		for (int i = 0; i < LSize; i++) {
      			if (pdEs[i + iColx] != 0) {
				f.fDecCoef(i+1, pdEs[i+iColx]*g.dGetCoef(j+1));
      			}
    		}
  	}
	return f;
} 

inline void
SchurMatrixHandlerUm::PrintMatrix(void)
{
	std::cout << "Schur Matrix " << std::endl;
	for (int i = 0; i < LSize; i++) {
		for (int j = 0; j < LSize; j++) {
 			std::cout << pB->dGetCoef(i+1, j+1) << " "; 
		}
		for (int j = 0; j < ISize; j++) {
			if (Eflag) {
				std::cout << pdE[i+j*LSize] << " ";
			} else {
				std::cout << pdEs[i+j*LSize] << " ";
			}
		}		
		std::cout << std::endl;
	}
	for (int i = 0;i < ISize; i++) {
		for (int j = 0;j < LSize; j++) {
 			std::cout << pF->dGetCoef(i+1, j+1) << " "; 
		}
		for (int j = 0; j < ISize; j++) {
			std::cout << pdC[i+j*ISize] << " ";
		}
		std::cout << std::endl;
	}
}

/* SchurMatrixHandlerUm - End*/

#endif /* SCHURMH_H */

