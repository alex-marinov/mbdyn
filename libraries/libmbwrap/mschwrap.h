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

#ifndef MSCHWRAP_H
#define MSCHWRAP_H

#ifdef USE_MESCHACH

#include <iostream>

#include <solman.h>
#include <submat.h>

extern "C" {   
/*
 * Meschach defines "min" and "max" as a macros,
 * screwing up C++ includes.
 * luckily the definitions are guarded by
 * #ifdef max and #ifdef min,
 * thus it should suffice to define them
 * before including anything from meschach
 * and later on undefining them
 */
#ifndef max
#define max max
#define have_defined_max_meschach_workaround
#endif

#ifndef min
#define min min
#define have_defined_max_meschach_workaround
#endif

#ifndef ANSI_C
#define ANSI_C
#define ANSI_C_meshach_defined
#endif

#define m_free int m_free
// we need only sparse2 here; 
// adding also matrix2, that is needed by
// mbdyn/elec/gpc.h to avoid having to dulicate there
// all the dirty tricks made here
#include <meschach/matrix2.h>
#include <meschach/sparse2.h>
#undef m_free
#ifdef ANSI_C_meshach_defined
#undef ANSI_C
#undef ANSI_C_meshach_defined
#endif


#ifdef have_defined_max_meschach_workaround
#undef max
#undef have_defined_max_meschach_workaround
#endif

#ifdef have_defined_min_meschach_workaround
#undef min
#undef have_defined_min_meschach_workaround
#endif
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
	doublereal* pdVecm1;
   
public:
	MeschachVectorHandler(integer iSize = 0);   
	virtual ~MeschachVectorHandler(void);

#ifdef DEBUG
   	/* Usata per il debug */
   	virtual void IsValid(void) const;
#endif /* DEBUG */
   
   	virtual inline doublereal* pdGetVec(void) const;
   	inline VEC* pGetMeschachVEC(void) const;

   	virtual inline integer iGetSize(void) const;
   
   	virtual void Resize(integer iNewSize);
   	virtual void Reset(void);
   
   	virtual inline void PutCoef(integer iRow, const doublereal& dCoef);
   	virtual inline void IncCoef(integer iRow, const doublereal& dCoef);
   	virtual inline void DecCoef(integer iRow, const doublereal& dCoef);
   	virtual inline const doublereal& dGetCoef(integer iRow) const;

	virtual inline const doublereal& operator () (integer iRow) const;
	virtual inline doublereal& operator () (integer iRow);
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

inline void
MeschachVectorHandler::PutCoef(integer iRow, const doublereal& dCoef) 
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   	pdVecm1[iRow] = dCoef;
   	return;
}

inline void
MeschachVectorHandler::IncCoef(integer iRow, const doublereal& dCoef) 
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   	pdVecm1[iRow] += dCoef;
   	return;
}

inline void
MeschachVectorHandler::DecCoef(integer iRow, const doublereal& dCoef) 
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   	pdVecm1[iRow] -= dCoef;
   	return;
}

inline const doublereal&
MeschachVectorHandler::dGetCoef(integer iRow) const 
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   	return pdVecm1[iRow];
}

inline const doublereal&
MeschachVectorHandler::operator () (integer iRow) const 
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   	return pdVecm1[iRow];
}

inline doublereal&
MeschachVectorHandler::operator () (integer iRow)
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
   	MeschachSparseMatrixHandler(integer m, integer n, integer maxlen = 0);
   	~MeschachSparseMatrixHandler(void);

   	/* helpers */
   	inline integer iGetNumRows(void) const;
   	inline integer iGetNumCols(void) const;
   
   	/* costruisce la matrice */
   	void Create(integer m, integer n, integer maxlen = 0);
#ifdef DEBUG
   	void IsValid(void) const;
#endif /* DEBUG */
   	void Reset(void);
   
   	/* Inserisce un coefficiente */
   	inline void PutCoef(integer iRow, integer iCol,
			     const doublereal& dCoef);
			     
   	/* Incrementa un coefficiente - se non esiste lo crea */
   	inline void IncCoef(integer iRow, integer iCol,
			     const doublereal& dCoef);
	
   	/* Decrementa un coefficiente - se non esiste lo crea */
   	inline void DecCoef(integer iRow, integer iCol,
			     const doublereal& dCoef);
			     
   	/* Restituisce un coefficiente - zero se non e' definito */
   	inline const doublereal& dGetCoef(integer iRow, integer iCol) const;
	
	inline const doublereal& operator()(integer iRow, integer iCol) const;
	inline doublereal& operator()(integer iRow, integer iCol);

	void Resize(integer m, integer n);

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
inline void
MeschachSparseMatrixHandler::PutCoef(integer iRow, integer iCol,
				      const doublereal& dCoef)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   	if (dCoef != 0.) {
      		sp_set_val(mat, --iRow, --iCol, dCoef);
      		return;
   	}
   	return;
}

/* Incrementa un coefficiente - se non esiste lo crea */
inline void
MeschachSparseMatrixHandler::IncCoef(integer iRow, integer iCol,
				      const doublereal& dCoef) 
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   	if (dCoef != 0.) {
		/* FIXME: this is an extension to Meschach */
      		/* sp_inc_val(mat, --iRow, --iCol, dCoef); */
		double tmp = sp_get_val(mat, --iRow, --iCol);
      		sp_set_val(mat, iRow, iCol, tmp+dCoef);
      		return;
   	}
   	return;
}

/* Decrementa un coefficiente - se non esiste lo crea */
inline void
MeschachSparseMatrixHandler::DecCoef(integer iRow, integer iCol,
				      const doublereal& dCoef) 
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   	if (dCoef != 0.) {
		/* FIXME: this is an extension to Meschach */
      		/* sp_dec_val(mat, --iRow, --iCol, dCoef); */
		double tmp = sp_get_val(mat, --iRow, --iCol);
      		sp_set_val(mat, iRow, iCol, tmp-dCoef);
      		return;
   	}
   	return;
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

inline const doublereal&
MeschachSparseMatrixHandler::operator()(integer iRow, integer iCol) const
{
	return dGetCoef(iRow, iCol);
}

inline doublereal&
MeschachSparseMatrixHandler::operator()(integer iRow, integer iCol)
{
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	return (doublereal&)dGetCoef(iRow, iCol);
}


inline SPMAT*
MeschachSparseMatrixHandler::pGetMAT(void) const 
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   	return mat;
}

extern std::ostream&
operator << (std::ostream& out, const MeschachSparseMatrixHandler& MH);



/* MeschachSparseMatrixHandler -end */

/* MeschachSparseSolutionManager - begin */

class MeschachSparseSolutionManager : public SolutionManager {
protected:
   	MeschachVectorHandler* prhs;
   	PERM* pivot;
   	MeschachSparseMatrixHandler* pmh;
   
   	enum { RESET, FACTORED } fStatus;
   	doublereal alpha;
   
   	void Create(integer iSize, integer iMaxSize);
   	void Factor(void);
   
public:
   	MeschachSparseSolutionManager(integer iSize,
					integer iMaxSize = 0, 
					const doublereal& a = 1.);
   	~MeschachSparseSolutionManager(void);

#ifdef DEBUG
   	void IsValid(void) const;
#endif /* DEBUG */
   	void MatrReset(void);
   
   	void Solve(void);
   	
	/* sposta il puntatore al vettore del residuo */
   	void pdSetResVec(doublereal* pRes){
		silent_cerr("Sorry Meschach is not available as local parallel solver. "
			<< "Aborting" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	};
   
   	/* sposta il puntatore al vettore del residuo */
   	void pdSetSolVec(doublereal* pSol) {
		silent_cerr("Sorry Meschach is not available as local parallel solver. "
			<< "Aborting" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	};   
   	
	/* Rende disponibile l'handler per la matrice */
   	virtual MatrixHandler* pMatHdl(void) const;
	
   	/* Rende disponibile l'handler per il termine noto */
   	virtual VectorHandler* pResHdl(void) const;
	
   	/* Rende disponibile l'handler per la soluzione (e' lo stesso
    	 * del termine noto, ma concettualmente sono separati) */
   	virtual VectorHandler* pSolHdl(void) const;
};

/* MeschachSparseSolutionManager - end */

#endif /* USE_MESCHACH */
    
#endif /* MSCHWRAP_H */
