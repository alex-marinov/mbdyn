/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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
 * Copyright (C) 1999-2014
 * Giuseppe Quaranta     <quaranta@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 */
/* Schur Matrix Handler */

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
	class ErrGeneric : public MBDynErrBase {
	public:
		ErrGeneric(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};

protected:
	integer LSize, ISize; /* dimensioni locali, interfacce */
	MatrixHandler* pB;
	MyVectorHandler* pE;
	doublereal* pdE;
	SpMapMatrixHandler* pF;
	MyVectorHandler* pC;
	doublereal* pdC;

 	/* Tabella di conversione Global to Local
	 * creata da SchurSolutionManager
	 * i nodi di interfaccia hanno indice negativo
	 * per permetterne la distinzione */
	integer* pGTL;

	bool bExtpdE;

public:
 	SchurMatrixHandler(int LocSize, int IntSize,
   			MatrixHandler* pBM,
   			integer* pGlobToLoc, doublereal* pdEv = 0);

	virtual ~SchurMatrixHandler(void);

#ifdef DEBUG
	/* Usata per il debug */
	virtual void IsValid(void) const;
#endif /* DEBUG */

	virtual integer iGetNumRows(void) const;
	virtual integer iGetNumCols(void) const;
	MatrixHandler* GetBMat(void);
	void SetBMat(MatrixHandler* pBM);
	doublereal* GetCMat(void);

	/* Resetta la matrice */
	virtual inline void Reset(void);

	/* Ridimensiona la matrice */
	virtual void Resize(integer, integer);

	/* Resetta la matrice */
	virtual inline void
	MatEFCReset(void);

	/* Inserisce un coefficiente */
	virtual inline void PutCoef(integer iRow, integer iCol,
			const doublereal& dCoef);

	/* Incrementa un coefficiente - se non esiste lo crea */
	virtual inline void IncCoef(integer iRow, integer iCol,
			const doublereal& dCoef);

	/* Incrementa un coefficiente - se non esiste lo crea */
	virtual inline void DecCoef(integer iRow, integer iCol,
			const doublereal& dCoef);

	/* Restituisce un coefficiente - zero se non e' definito */
	virtual inline const doublereal&
	dGetCoef(integer iRow, integer iCol) const;

	virtual const doublereal&
	operator()(integer iRow, integer iCol) const;

	virtual doublereal&
	operator()(integer iRow, integer iCol);

	virtual inline doublereal* GetECol(const integer iCol) const;
	virtual inline doublereal* GetEColSol(const integer iCol) const;

	/* calacola g - F*f e lo pone in g */
	virtual inline VectorHandler&
	CompNewg(VectorHandler& g, const VectorHandler& f) const;

	/* Calcola la matrice di schur locale C - F*E' e la memorizza in C */
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
SchurMatrixHandler::MatEFCReset(void)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	pE->Reset();
	pF->Reset();
	pC->Reset();
}

inline void
SchurMatrixHandler::Reset(void)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	pB->Reset();
	MatEFCReset();
}

inline void
SchurMatrixHandler::PutCoef(integer iRow, integer iCol,
		const doublereal& dCoef)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	if (pGTL[iRow] > 0) {
		if (pGTL[iCol] > 0) {
			pB->PutCoef(pGTL[iRow], pGTL[iCol], dCoef);
			return;

		} else if (pGTL[iCol] < 0) {
			pE->PutCoef(pGTL[iRow] - (pGTL[iCol] + 1)*LSize, dCoef);
			return;
		}

	} else if (pGTL[iRow] < 0) {
		if (pGTL[iCol] > 0) {
			pF->PutCoef(-pGTL[iRow], pGTL[iCol], dCoef);
			return;

		} else if (pGTL[iCol] < 0) {
			pC->PutCoef(-pGTL[iRow] - (pGTL[iCol] + 1)*ISize, dCoef);
			return;
		}
	}

#ifdef USE_MPI
	silent_cerr("SchurMatrixHandler::PutCoef() "
		"Process(" << MBDynComm.Get_rank() << "): "
		"trying to operate on nonlocal indices "
		<< iRow << "," << iCol << std::endl);
#else /* ! USE_MPI */
	silent_cerr("SchurMatrixHandler::PutCoef() "
		"trying to operate on nonlocal indices "
		<< iRow << "," << iCol << std::endl);
#endif /* ! USE_MPI */

	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

inline void
SchurMatrixHandler::IncCoef(integer iRow, integer iCol,
		const doublereal& dCoef)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	if (pGTL[iRow] > 0) {
		if (pGTL[iCol] > 0) {
			pB->IncCoef(pGTL[iRow], pGTL[iCol], dCoef);
			return;

		} else if (pGTL[iCol] < 0) {
			pE->IncCoef(pGTL[iRow] - (pGTL[iCol] + 1)*LSize, dCoef);
			return;
		}

	} else if (pGTL[iRow] < 0) {
		if (pGTL[iCol] > 0) {
			pF->IncCoef(-pGTL[iRow], pGTL[iCol], dCoef);
			return;

		} else if (pGTL[iCol] < 0) {
			pC->IncCoef(-pGTL[iRow] - (pGTL[iCol] + 1)*ISize, dCoef);
			return;
		}
	}

#ifdef USE_MPI
	silent_cerr("SchurMatrixHandler::IncCoef() "
		"Process(" << MBDynComm.Get_rank() << "): "
		"trying to operate on nonlocal indices "
		<< iRow << "," << iCol << std::endl);
#else /* ! USE_MPI */
	silent_cerr("SchurMatrixHandler::IncCoef() "
		"trying to operate on nonlocal indices "
		<< iRow << "," << iCol << std::endl);
#endif /* ! USE_MPI */

	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

inline void
SchurMatrixHandler::DecCoef(integer iRow, integer iCol,
		const doublereal& dCoef)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	if (pGTL[iRow] > 0) {
		if (pGTL[iCol] > 0) {
			pB->DecCoef(pGTL[iRow], pGTL[iCol], dCoef);
			return;

		} else if (pGTL[iCol] < 0) {
			pE->DecCoef(pGTL[iRow] - (pGTL[iCol] + 1)*LSize, dCoef);
			return;
		}

	} else if (pGTL[iRow] < 0) {
		if (pGTL[iCol] > 0) {
			pF->DecCoef(-pGTL[iRow], pGTL[iCol], dCoef);
			return;

		} else if (pGTL[iCol] < 0) {
			pC->DecCoef(-pGTL[iRow] - (pGTL[iCol] + 1)*ISize, dCoef);
			return;
		}
	}

#ifdef USE_MPI
	silent_cerr("SchurMatrixHandler::DecCoef() "
		"Process(" << MBDynComm.Get_rank() << "): "
		"trying to operate on nonlocal indices "
		<< iRow << "," << iCol << std::endl);
#else /* ! USE_MPI */
	silent_cerr("SchurMatrixHandler::DecCoef() "
		"trying to operate on nonlocal indices "
		<< iRow << "," << iCol << std::endl);
#endif /* ! USE_MPI */

	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

inline const doublereal&
SchurMatrixHandler::dGetCoef(integer iRow, integer iCol) const
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	if (pGTL[iRow] > 0) {
		if (pGTL[iCol] > 0) {
			return pB->operator()(pGTL[iRow], pGTL[iCol]);

		} else if(pGTL[iCol] < 0) {
			return pE->operator()(pGTL[iRow] - (pGTL[iCol] + 1)*LSize);
		}

	} else if (pGTL[iRow] < 0) {
		if (pGTL[iCol] > 0) {
			return pF->operator()(-pGTL[iRow], pGTL[iCol]);

		} else if (pGTL[iCol] < 0) {
			return pC->operator()(-pGTL[iRow] - (pGTL[iCol] + 1)*ISize);
		}
	}

#ifdef USE_MPI
	silent_cerr("SchurMatrixHandler::dGetCoef() "
		"Process(" << MBDynComm.Get_rank() << "): "
		"trying to operate on nonlocal indices "
		<< iRow << "," << iCol << std::endl);
#else /* ! USE_MPI */
	silent_cerr("SchurMatrixHandler::dGetCoef() "
		"trying to operate on nonlocal indices "
		<< iRow << "," << iCol << std::endl);
#endif /* ! USE_MPI */

	return ::Zero1;
}

inline doublereal&
SchurMatrixHandler::operator()(integer iRow, integer iCol)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	if (pGTL[iRow] > 0) {
		if (pGTL[iCol] > 0) {
			return (*pB)(pGTL[iRow], pGTL[iCol]);

		} else if (pGTL[iCol] < 0) {
			return (*pE)(pGTL[iRow] - (pGTL[iCol] + 1)*LSize);
		}

	} else {
		if (pGTL[iCol] > 0) {
			return (*pF)(-pGTL[iRow], pGTL[iCol]);

		} else if (pGTL[iCol] < 0) {
			return (*pC)(-pGTL[iRow] - (pGTL[iCol] + 1)*ISize);
		}
	}

#ifdef USE_MPI
	silent_cerr("SchurMatrixHandler::operator() "
		"Process(" << MBDynComm.Get_rank() << "): "
		"trying to operate on nonlocal indices "
		<< iRow << "," << iCol << std::endl);
#else /* ! USE_MPI */
	silent_cerr("SchurMatrixHandler::operator() "
		"trying to operate on nonlocal indices "
		<< iRow << "," << iCol << std::endl);
#endif /* ! USE_MPI */

	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

inline const doublereal&
SchurMatrixHandler::operator () (integer iRow, integer iCol) const
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	if (pGTL[iRow] > 0) {
		if (pGTL[iCol] > 0) {
			return pB->operator()(pGTL[iRow], pGTL[iCol]);

		} else if (pGTL[iCol] < 0) {
			return pE->operator()(pGTL[iRow] - (pGTL[iCol] + 1)*LSize);
		}

	} else {
		if (pGTL[iCol] > 0) {
			return pF->operator()(-pGTL[iRow], pGTL[iCol]);

		} else if (pGTL[iCol] < 0) {
			return pC->operator()(-pGTL[iRow] - (pGTL[iCol] + 1)*ISize);
		}
	}

#ifdef USE_MPI
	silent_cerr("SchurMatrixHandler::operator() "
		"Process(" << MBDynComm.Get_rank() << "): "
		"trying to operate on nonlocal indices "
		<< iRow << "," << iCol << std::endl);
#else /* ! USE_MPI */
	silent_cerr("SchurMatrixHandler::operator() "
		"trying to operate on nonlocal indices "
		<< iRow << "," << iCol << std::endl);
#endif /* ! USE_MPI */

	return ::Zero1;
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
					pF->operator()(i + 1, k + 1) * pdE[k + iCole];
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
				f.DecCoef(i + 1, pdE[i + iColx]*g(j + 1));
      			}
    		}
  	}
	return f;
}

inline void
SchurMatrixHandler::PrintMatrix(void)
{
	silent_cout("Schur Matrix " << std::endl);

	for (int i = 0; i < LSize; i++) {
		for (int j = 0; j < LSize; j++) {
 			silent_cout(pB->operator()(i + 1, j + 1) << " ");
		}

		for (int j = 0; j < ISize; j++) {
			silent_cout(pdE[i + j*LSize] << " ");
		}
		silent_cout(std::endl);
	}

	for (int i = 0; i < ISize; i++) {
		for (int j = 0; j < LSize; j++) {
 			silent_cout(pF->operator()(i + 1, j + 1) << " ");
		}

		for (int j = 0; j < ISize; j++) {
			silent_cout(pdC[i + j*ISize] << " ");
		}
		silent_cout(std::endl);
	}
}

/* SchurMatrixHandler - End */

/* SchurVectorHandler - Start */

class SchurVectorHandler : public VectorHandler {
public:
	class ErrGeneric: public MBDynErrBase {
	public:
		ErrGeneric(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};

private:
	integer LSize, ISize;
	VectorHandler* pLV;
	VectorHandler* pIV;
	bool bExtpIV;
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

#ifdef DEBUG
	/* Usata per il debug */
	void IsValid(void) const;
#endif /* DEBUG */

	/* restituisce il puntatore al vettore */
	inline doublereal* pdGetVec(void) const;

	/* restituisce le dimensioni del vettore */
	inline integer iGetSize(void) const;

	virtual void Reset(void);

	inline void Resize(integer iNewSize);

	inline VectorHandler* GetIVec(void);
	inline VectorHandler* GetLVec(void);

	inline void PutCoef(integer iRow, const doublereal& dCoef);
	inline void IncCoef(integer iRow, const doublereal& dCoef);
	inline void DecCoef(integer iRow, const doublereal& dCoef);
	inline const doublereal& dGetCoef(integer iRow) const;

	inline const doublereal& operator () (integer iRow) const;
	inline doublereal& operator () (integer iRow);

	inline void PrintVector(void);
};

/* restituisce il puntatore al vettore */
inline doublereal*
SchurVectorHandler::pdGetVec(void) const
{
	silent_cerr("You shouldn't have asked for "
		"the internal pointer of a SchurVectorHandler"
		<< std::endl);
	return pLV->pdGetVec();
}

/* restituisce le dimensioni del vettore */
inline integer
SchurVectorHandler::iGetSize(void) const
{
	return LSize + ISize;
}

inline void
SchurVectorHandler::Resize(integer iNewSize)
{
	silent_cerr("Why are you trying to resize a SchurVector ???? "
		<< "No Operation Performed!!" << std::endl);
}

/* assegna 0. a tutti gli elementi del vettore */
inline void
SchurVectorHandler::Reset(void)
{
	pLV->Reset();
	pIV->Reset();
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

inline void
SchurVectorHandler::PutCoef(integer iRow, const doublereal& dCoef)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	if (pGTL[iRow] > 0) {
		pLV->PutCoef(pGTL[iRow], dCoef);

	} else if (pGTL[iRow] < 0 ) {
		pIV->PutCoef(-pGTL[iRow], dCoef);

	} else {
#ifdef USE_MPI
		silent_cerr("SchurVectorHandler::PutCoef() "
			"Process(" << MBDynComm.Get_rank() << "): "
			"trying to operate on nonlocal index "
			<< iRow << std::endl);
#else /* ! USE_MPI */
		silent_cerr("SchurVectorHandler::PutCoef() "
			"trying to operate on nonlocal index "
			<< iRow << std::endl);
#endif /* ! USE_MPI */
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	return;
}

inline void
SchurVectorHandler::IncCoef(integer iRow, const doublereal& dCoef)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	if (pGTL[iRow] > 0) {
		pLV->IncCoef(pGTL[iRow], dCoef);

	} else if (pGTL[iRow] < 0) {
		pIV->IncCoef(-pGTL[iRow], dCoef);

	} else {
#ifdef USE_MPI
		silent_cerr("SchurVectorHandler::IncCoef "
			"Process(" << MBDynComm.Get_rank() << "): "
			"trying to operate on nonlocal index "
			<< iRow << std::endl);
#else /* ! USE_MPI */
		silent_cerr("SchurVectorHandler::IncCoef "
			"trying to operate on nonlocal index "
			<< iRow << std::endl);
#endif /* ! USE_MPI */
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	return;
}

inline void
SchurVectorHandler::DecCoef(integer iRow, const doublereal& dCoef)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	if (pGTL[iRow] > 0) {
		pLV->DecCoef(pGTL[iRow], dCoef);

	} else if (pGTL[iRow] < 0) {
		pIV->DecCoef(-pGTL[iRow], dCoef);

	} else {
#ifdef USE_MPI
		silent_cerr("SchurVectorHandler::DecCoef "
			"Process(" << MBDynComm.Get_rank() << "): "
			"trying to operate on nonlocal index "
			<< iRow << std::endl);
#else /* ! USE_MPI */
		silent_cerr("SchurVectorHandler::DecCoef "
			"trying to operate on nonlocal index "
			<< iRow << std::endl);
#endif /* ! USE_MPI */
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	return;
}

inline const doublereal&
SchurVectorHandler::dGetCoef(integer iRow) const
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	if (pGTL[iRow] > 0) {
		return pLV->operator()(pGTL[iRow]);

	} else if (pGTL[iRow] < 0) {
		return pIV->operator()(-pGTL[iRow]);

	} else {
#ifdef USE_MPI
		silent_cerr("SchurVectorHandler::dGetCoef "
			"Process(" << MBDynComm.Get_rank() << "): "
			"trying to operate on nonlocal index "
			<< iRow << std::endl);
#else /* ! USE_MPI */
		silent_cerr("SchurVectorHandler::dGetCoef "
			"trying to operate on nonlocal index "
			<< iRow << std::endl);
#endif /* ! USE_MPI */
		return Zero1;
	}
}

inline const doublereal&
SchurVectorHandler::operator()(integer iRow) const
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	if (pGTL[iRow] > 0) {
		return pLV->operator()(pGTL[iRow]);

	} else if (pGTL[iRow] < 0) {
		return pIV->operator()(-pGTL[iRow]);

	} else {
#ifdef USE_MPI
		silent_cerr("SchurVectorHandler::operator()"
			"Process(" << MBDynComm.Get_rank() << "): "
			"trying to operate on nonlocal index "
			<< iRow << std::endl);
#else /* ! USE_MPI */
		silent_cerr("SchurVectorHandler::operator()"
			"trying to operate on nonlocal index "
			<< iRow << std::endl);
#endif /* ! USE_MPI */
		return Zero1;
	}
}

inline doublereal&
SchurVectorHandler::operator()(integer iRow)
{
	silent_cerr("SchurVectorHandler::operator() "
		"cannot be used on nonlocal index "
		<< iRow << std::endl);

	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

inline void
SchurVectorHandler::PrintVector(void)
{
	silent_cout("Schur Vector " << std::endl);

	for (int j = 0; j < LSize; j++) {
 		silent_cout(pLV->operator()(j + 1) << " " << std::endl);
	}

	for (int j = 0; j < ISize; j++) {
 		silent_cout(pIV->operator()(j + 1) << " " << std::endl);
	}
}

/* SchurVectorHandler - End */


/* SchurMatrixHandlerUm - begin */

class SchurMatrixHandlerUm : public SchurMatrixHandler {
public:
	class ErrGeneric : public MBDynErrBase {
	public:
		ErrGeneric(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};


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
	inline void MatEFCReset(void);

	/* Resetta la matrice */
	inline void Reset(void);

	/* Inserisce un coefficiente */
	inline void PutCoef(integer iRow, integer iCol,
			const doublereal& dCoef);

	/* Incrementa un coefficiente - se non esiste lo crea */
	inline void IncCoef(integer iRow, integer iCol,
			const doublereal& dCoef);

	/* Incrementa un coefficiente - se non esiste lo crea */
	inline void DecCoef(integer iRow, integer iCol,
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
SchurMatrixHandlerUm::MatEFCReset(void)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	Eflag = true;
	for (int i = 0; i < LSize*(ISize + 1); i++) {
		pdEs[i] = 0.;
	}
	pF->Reset();
	pC->Reset();
}

inline void SchurMatrixHandlerUm::Reset(void)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	Eflag = true;
	pB->Reset();
	MatEFCReset();
}

inline void
SchurMatrixHandlerUm::PutCoef(integer iRow, integer iCol,
		const doublereal& dCoef)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	if (pGTL[iRow] > 0) {
		if (pGTL[iCol] > 0) {
			pB->PutCoef(pGTL[iRow], pGTL[iCol], dCoef);
			return;

		} else if (pGTL[iCol] < 0 ) {
			if (Eflag) {
				pE->PutCoef(pGTL[iRow] - (pGTL[iCol] + 1)*LSize, dCoef);

			} else {
				pEs->PutCoef(pGTL[iRow] - (pGTL[iCol] + 1)*LSize, dCoef);
			}
			return;
		}

	} else if (pGTL[iRow] < 0 ) {
		if (pGTL[iCol] > 0) {
			pF->PutCoef(-pGTL[iRow], pGTL[iCol], dCoef);
			return;

		} else if (pGTL[iCol] < 0) {
			pC->PutCoef(-pGTL[iRow] - (pGTL[iCol] + 1)*ISize, dCoef);
			return;
		}
	}

#ifdef USE_MPI
	silent_cerr("SchurMatrixHandlerUm::PutCoef() "
		"Process(" << MBDynComm.Get_rank() << "): "
		"trying to operate on nonlocal indices "
		<< iRow << "," << iCol << std::endl);
#else /* ! USE_MPI */
	silent_cerr("SchurMatrixHandlerUm::PutCoef() "
		"trying to operate on nonlocal indices "
		<< iRow << "," << iCol << std::endl);
#endif /* ! USE_MPI */

	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

inline void
SchurMatrixHandlerUm::IncCoef(integer iRow, integer iCol,
		const doublereal& dCoef)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	if (pGTL[iRow] > 0) {
		if (pGTL[iCol] > 0) {
			pB->IncCoef(pGTL[iRow], pGTL[iCol], dCoef);
			return;

		} else if (pGTL[iCol] < 0) {
			if (Eflag) {
				pE->IncCoef(pGTL[iRow] - (pGTL[iCol] + 1)*LSize, dCoef);

			} else {
				pEs->IncCoef(pGTL[iRow] - (pGTL[iCol] + 1)*LSize, dCoef);
			}
			return;
		}

	} else if (pGTL[iRow] < 0) {
		if (pGTL[iCol] > 0) {
			pF->IncCoef(-pGTL[iRow], pGTL[iCol], dCoef);
			return;

		} else if (pGTL[iCol] < 0) {
			pC->IncCoef(-pGTL[iRow] - (pGTL[iCol] + 1)*ISize, dCoef);
			return;
		}
	}

#ifdef USE_MPI
	silent_cerr("SchurMatrixHandlerUm::IncCoef() "
		"Process(" << MBDynComm.Get_rank() << "): "
		"trying to operate on nonlocal indices "
		<< iRow << "," << iCol << std::endl);
#else /* ! USE_MPI */
	silent_cerr("SchurMatrixHandlerUm::IncCoef() "
		"trying to operate on nonlocal indices "
		<< iRow << "," << iCol << std::endl);
#endif /* ! USE_MPI */

	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

inline void
SchurMatrixHandlerUm::DecCoef(integer iRow, integer iCol,
		const doublereal& dCoef)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	if (pGTL[iRow] > 0) {
		if (pGTL[iCol] > 0) {
			pB->DecCoef(pGTL[iRow], pGTL[iCol], dCoef);
			return;

		} else if (pGTL[iCol] < 0) {
			if (Eflag) {
				pE->DecCoef(pGTL[iRow] - (pGTL[iCol] + 1)*LSize, dCoef);

			} else {
				pEs->DecCoef(pGTL[iRow] - (pGTL[iCol] + 1)*LSize, dCoef);
			}
			return;
		}

	} else if (pGTL[iRow] < 0) {
		if (iCol > 0) {
			pF->DecCoef(-pGTL[iRow], pGTL[iCol], dCoef);
			return;

		} else if (iCol < 0) {
			pC->DecCoef(-pGTL[iRow] - (pGTL[iCol] + 1)*ISize, dCoef);
			return;
		}
	}

#ifdef USE_MPI
	silent_cerr("SchurMatrixHandlerUm::DecCoef() "
		"Process(" << MBDynComm.Get_rank() << "): "
		"trying to operate on nonlocal indices "
		<< iRow << "," << iCol << std::endl);
#else /* ! USE_MPI */
	silent_cerr("SchurMatrixHandlerUm::DecCoef() "
		"trying to operate on nonlocal indices "
		<< iRow << "," << iCol << std::endl);
#endif /* ! USE_MPI */

	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

inline const doublereal&
SchurMatrixHandlerUm::dGetCoef(integer iRow, integer iCol) const
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	if (pGTL[iRow] > 0) {
		if (pGTL[iCol] > 0) {
			return pB->operator()(pGTL[iRow], pGTL[iCol]);

		} else if (pGTL[iCol] < 0) {
			if (Eflag) {
				return pE->operator()(pGTL[iRow] - (pGTL[iCol] + 1)*LSize);

			} else {
				return pEs->operator()(pGTL[iRow] - (pGTL[iCol] + 1)*LSize);
			}
		}

	} else if (pGTL[iRow] < 0) {
		if (pGTL[iCol] > 0) {
			return pF->operator()(-pGTL[iRow], pGTL[iCol]);

		} else if (pGTL[iCol] < 0) {
			return pC->operator()(-pGTL[iRow] - (pGTL[iCol] + 1)*ISize);
		}
	}

#ifdef USE_MPI
	silent_cerr("SchurMatrixHandlerUm::dGetCoef() "
		"Process(" << MBDynComm.Get_rank() << "): "
		"trying to operate on nonlocal indices "
		<< iRow << "," << iCol << std::endl);
#else /* ! USE_MPI */
	silent_cerr("SchurMatrixHandlerUm::dGetCoef() "
		"trying to operate on nonlocal indices "
		<< iRow << "," << iCol << std::endl);
#endif /* ! USE_MPI */

	return ::Zero1;
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
					pF->operator()(i + 1, k + 1) * pdEs[k + iCole];
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
				f.DecCoef(i + 1, pdEs[i + iColx]*g(j + 1));
      			}
    		}
  	}
	return f;
}

inline void
SchurMatrixHandlerUm::PrintMatrix(void)
{
	silent_cout("Schur Matrix " << std::endl);

	for (int i = 0; i < LSize; i++) {
		for (int j = 0; j < LSize; j++) {
 			silent_cout(pB->operator()(i + 1, j + 1) << " ");
		}

		for (int j = 0; j < ISize; j++) {
			if (Eflag) {
				silent_cout(pdE[i + j*LSize] << " ");
			} else {
				silent_cout(pdEs[i + j*LSize] << " ");
			}
		}
		silent_cout(std::endl);
	}

	for (int i = 0; i < ISize; i++) {
		for (int j = 0;j < LSize; j++) {
 			silent_cout(pF->operator()(i + 1, j + 1) << " ");
		}

		for (int j = 0; j < ISize; j++) {
			silent_cout(pdC[i + j*ISize] << " ");
		}
		silent_cout(std::endl);
	}
}

/* SchurMatrixHandlerUm - End */

#endif /* SCHURMH_H */

