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

#ifndef MSCHWRAP_H
#define MSCHWRAP_H

#ifdef USE_MESCHACH

#include <solman.h>
#include <submat.h>

extern "C" {   
#include <meschach/sparse2.h>
}

/*
 * Meschach defines "catch" as a macro,
 * which conflicts with the reserved C++ word
 */
#undef catch

/* MeschachVectorHandler - begin */

class MeschachVectorHandler : public VectorHandler {
protected:
	VEC* pv;
	double* pdVecm1;
   
public:
	MeschachVectorHandler(int iSize = 0);   
	virtual ~MeschachVectorHandler(void);
   
   	/* Usata per il debug */
   	virtual void IsValid(void) const;
   
   	virtual inline doublereal* pdGetVec(void) const;
   	inline VEC* pGetMeschachVEC(void) const;

   	virtual inline integer iGetSize(void) const;
   
   	virtual void Resize(integer iNewSize);
   	virtual void Reset(doublereal dResetVal = 0.);
   
   	virtual inline flag fPutCoef(integer iRow, const doublereal& dCoef);
   	virtual inline flag fIncCoef(integer iRow, const doublereal& dCoef);
   	virtual inline flag fDecCoef(integer iRow, const doublereal& dCoef);
   	virtual inline const doublereal& dGetCoef(integer iRow) const;
};


inline doublereal*
MeschachVectorHandler::pdGetVec(void) const
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   	return pv->ve;
}

inline VEC*
MeschachVectorHandler::pGetMeschachVEC(void) const
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   	return pv;
}

inline integer
MeschachVectorHandler::iGetSize(void) const
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   	return integer(pv->dim);
}

inline flag
MeschachVectorHandler::fPutCoef(integer iRow, const doublereal& dCoef) 
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   	pdVecm1[iRow] = dCoef;
   	return flag(1);
}

inline flag
MeschachVectorHandler::fIncCoef(integer iRow, const doublereal& dCoef) 
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   	pdVecm1[iRow] += dCoef;
   	return flag(1);
}

inline flag
MeschachVectorHandler::fDecCoef(integer iRow, const doublereal& dCoef) 
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   	pdVecm1[iRow] -= dCoef;
   	return flag(1);
}

inline const doublereal&
MeschachVectorHandler::dGetCoef(integer iRow) const 
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   	return pdVecm1[iRow];
}

/* MeschachVectorHandler - end */

/* MeschachSparseMatrixHandler - begin */

class MeschachSparseMatrixHandler : public MatrixHandler {
private:
   	doublereal dDmy;
   
 	protected:
   	SPMAT* mat;
   
public:
   	MeschachSparseMatrixHandler(int m, int n, int maxlen = 0);
   	~MeschachSparseMatrixHandler(void);

   	/* helpers */
   	inline integer iGetNumRows(void) const;
   	inline integer iGetNumCols(void) const;
   
   	/* costruisce la matrice */
   	void Create(unsigned int m, unsigned int n, unsigned int maxlen = 0);
   	void IsValid(void) const;
   	void Init(const doublereal& d = 0.);
   
   	/* Inserisce un coefficiente */
   	inline flag fPutCoef(integer iRow, integer iCol,
			     const doublereal& dCoef);
			     
   	/* Incrementa un coefficiente - se non esiste lo crea */
   	inline flag fIncCoef(integer iRow, integer iCol,
			     const doublereal& dCoef);
	
   	/* Decrementa un coefficiente - se non esiste lo crea */
   	inline flag fDecCoef(integer iRow, integer iCol,
			     const doublereal& dCoef);
			     
   	/* Restituisce un coefficiente - zero se non e' definito */
   	inline const doublereal& dGetCoef(integer iRow, integer iCol) const;

   	inline SPMAT* pGetMAT(void) const;
};

/* helpers */
inline integer
MeschachSparseMatrixHandler::iGetNumRows(void) const
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   	return integer(mat->m);
}
   
inline integer
MeschachSparseMatrixHandler::iGetNumCols(void) const 
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
  	return integer(mat->n);
}

/* Inserisce un coefficiente */
inline flag
MeschachSparseMatrixHandler::fPutCoef(integer iRow, integer iCol,
				      const doublereal& dCoef)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   	if (dCoef != 0.) {
      		sp_set_val(mat, --iRow, --iCol, dCoef);
      		return flag(0);
   	}
   	return flag(1);
}

/* Incrementa un coefficiente - se non esiste lo crea */
inline flag
MeschachSparseMatrixHandler::fIncCoef(integer iRow, integer iCol,
				      const doublereal& dCoef) 
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   	if (dCoef != 0.) {
		/* FIXME: this is an extension to Meschach */
      		sp_inc_val(mat, --iRow, --iCol, dCoef);
      		return flag(0);
   	}
   	return flag(1);
}

/* Decrementa un coefficiente - se non esiste lo crea */
inline flag
MeschachSparseMatrixHandler::fDecCoef(integer iRow, integer iCol,
				      const doublereal& dCoef) 
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   	if (dCoef != 0.) {
		/* FIXME: this is an extension to Meschach */
      		sp_dec_val(mat, --iRow, --iCol, dCoef);
      		return flag(0);
   	}
   	return flag(1);
}

/* Restituisce un coefficiente - zero se non e' definito */
inline const doublereal&
MeschachSparseMatrixHandler::dGetCoef(integer iRow, integer iCol) const
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   	return ((doublereal&)dDmy = sp_get_val(mat, --iRow, --iCol));
}

inline SPMAT*
MeschachSparseMatrixHandler::pGetMAT(void) const 
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   	return mat;
}

extern ostream&
operator << (ostream& out, const MeschachSparseMatrixHandler& MH);



/* MeschachSparseMatrixHandler -end */

/* MeschachSparseLUSolutionManager - begin */

class MeschachSparseLUSolutionManager : public SolutionManager {
protected:
   	MeschachVectorHandler* prhs;
   	PERM* pivot;
   	MeschachSparseMatrixHandler* pmh;
   
   	enum { RESET, FACTORED } fStatus;
   	double alpha;
   
   	void Create(unsigned int iSize, unsigned int iMaxSize);
   	void Factor(void);
   
public:
   	MeschachSparseLUSolutionManager(int iSize,
					int iMaxSize = 0, 
					double a = 1.);
   	~MeschachSparseLUSolutionManager(void);

   	void IsValid(void) const;
   	void MatrInit(const doublereal& d = 0.);
   
   	void Solve(void);
   
   	/* Rende disponibile l'handler per la matrice */
   	virtual MatrixHandler* pMatHdl(void) const;
	
   	/* Rende disponibile l'handler per il termine noto */
   	virtual VectorHandler* pResHdl(void) const;
	
   	/* Rende disponibile l'handler per la soluzione (e' lo stesso
    	 * del termine noto, ma concettualmente sono separati) */
   	virtual VectorHandler* pSolHdl(void) const;
};

/* MeschachSparseLUSolutionManager - end */

#endif /* USE_MESCHACH */
    
#endif /* MSCHWRAP_H */
