/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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

/*****************************************************************************
 *                                                                           *
 *                            SOLUTION MANAGER                               *
 *                                                                           *
 *****************************************************************************/

/* Pierangelo Masarati */


#ifndef SOLMAN_H
#define SOLMAN_H

#include <ac/math.h>
#include <ac/iostream>
#include <ac/f2c.h>

/* per il debugging */
#include <myassert.h>
#include <mynewmem.h>
#include <except.h>

/* Zero for sparse vector and matrix handlers */
extern const doublereal dZero;

/* classi virtuali dichiarate in questo file */
class MatrixHandler;    /* gestore matrice */
class VectorHandler;    /* gestore vettore */
class SolutionManager;  /* gestore della soluzione */

/* classi usate in questo file */
class SubMatrixHandler;
class FullSubMatrixHandler;
class SparseSubMatrixHandler;
class VariableSubMatrixHandler;
class SubVectorHandler;
class Vec3;
class Mat3x3;

/* MatrixHandler - begin */

class MatrixHandler {
public:
	class ErrGeneric {};

public:
	virtual ~MatrixHandler(void);

	/* Usata per il debug */
	virtual void IsValid(void) const = 0;

	/* Resetta la matrice ecc. */
	virtual void Init(const doublereal& dResetVal = 0.) = 0;

	/* Restituisce un puntatore all'array di reali della matrice */
	virtual inline doublereal* pdGetMat(void) const;

	/* Restituisce un puntatore al vettore delle righe */
	virtual inline integer* piGetRows(void) const;

	/* Restituisce un puntatore al vettore delle colonne */
	virtual inline integer* piGetCols(void) const;

	/* Impacchetta la matrice; restituisce il numero di elementi 
	 * diversi da zero */
	virtual integer PacMat(void);

	/* Resetta la matrice ecc. */
	virtual void Reset(const doublereal& dResetVal = 0.);

	/* Inserisce un coefficiente */
	virtual void
	PutCoef(integer iRow, integer iCol, const doublereal& dCoef) = 0;

	/* Incrementa un coefficiente - se non esiste lo crea */
	virtual void
	IncCoef(integer iRow, integer iCol, const doublereal& dCoef) = 0;

	/* Decrementa un coefficiente - se non esiste lo crea */
	virtual void
	DecCoef(integer iRow, integer iCol, const doublereal& dCoef) = 0;

	/* Restituisce un coefficiente - zero se non e' definito */
	virtual const doublereal&
	dGetCoef(integer iRow, integer iCol) const = 0;

	virtual const doublereal&
	operator () (integer iRow, integer iCol) const = 0;

	virtual doublereal&
	operator () (integer iRow, integer iCol) = 0;

	/* dimensioni */
	virtual integer iGetNumRows(void) const = 0;
	virtual integer iGetNumCols(void) const = 0;

	/* Overload di += usato per l'assemblaggio delle matrici */
	virtual MatrixHandler& operator += (const SubMatrixHandler& SubMH);

	/* Overload di -= usato per l'assemblaggio delle matrici */
	virtual MatrixHandler& operator -= (const SubMatrixHandler& SubMH);

	/* Overload di += usato per l'assemblaggio delle matrici
	 * questi li vuole ma non so bene perche'; force per la doppia
	 * derivazione di VariableSubMatrixHandler */
	virtual MatrixHandler&
	operator += (const VariableSubMatrixHandler& SubMH);
	virtual MatrixHandler&
	operator -= (const VariableSubMatrixHandler& SubMH);

	/* */
	virtual MatrixHandler& ScalarMul(const doublereal& d);

	virtual VectorHandler&
	MatVecMul(VectorHandler& out, const VectorHandler& in) const;
	virtual VectorHandler&
	MatTVecMul(VectorHandler& out, const VectorHandler& in) const;
	virtual VectorHandler&
	MatVecIncMul(VectorHandler& out, const VectorHandler& in) const;
	virtual VectorHandler&
	MatTVecIncMul(VectorHandler& out, const VectorHandler& in) const;
};

/* Restituisce un puntatore all'array di reali della matrice */
inline doublereal*
MatrixHandler::pdGetMat(void) const
{
	return NULL;
}

/* Restituisce un puntatore al vettore delle righe */
inline integer*
MatrixHandler::piGetRows(void) const
{
	return NULL;
}

/* Restituisce un puntatore al vettore delle colonne */
inline integer*
MatrixHandler::piGetCols(void) const
{
	return NULL;
}

extern std::ostream&
operator << (std::ostream& out, const MatrixHandler& MH);

/* MatrixHandler - end */


/* VectorHandler - begin */

/* Gestore di vettori. Usa spazio messo a disposizione da altri. */

/* questa classe non e' virtuale pura, ma viene definita nel caso banale di
 * vettore pieno in quanto e' di interesse pratico. Il pacchetto <harwrap>
 * si basa su questa per la gestione dei vettori residuo e soluzione */

class VectorHandler {
public:
	virtual ~VectorHandler(void);

	/* Usata per il debug */
	virtual void IsValid(void) const = 0;

	virtual doublereal* pdGetVec(void) const = 0;

	virtual integer iGetSize(void) const = 0;

	virtual void Resize(integer iNewSize) = 0;

	virtual void Reset(doublereal dResetVal = 0.) = 0;

	virtual void PutCoef(integer iRow, const doublereal& dCoef) = 0;

	virtual void IncCoef(integer iRow, const doublereal& dCoef) = 0;

	virtual void DecCoef(integer iRow, const doublereal& dCoef) = 0;

	virtual const doublereal& dGetCoef(integer iRow) const = 0;

	/* Somma un Vec3 nella posizione desiderata */
	virtual void Add(integer iRow, const Vec3& v);

	/* Sottrae un Vec3 nella posizione desiderata */
	virtual void Sub(integer iRow, const Vec3& v);

	/* Scrive un Vec3 nella posizione desiderata */
	virtual void Put(integer iRow, const Vec3& v);

	/* Somma e moltiplica per uno scalare */
	virtual VectorHandler&
	ScalarAddMul(const VectorHandler& VH, const doublereal& d);

	/* Somma e moltiplica per uno scalare v = VH + d * VH1 */
	virtual VectorHandler&
	ScalarAddMul(const VectorHandler& VH, const VectorHandler& VH1,
			const doublereal& d);

	/* Moltiplica per uno scalare */
	virtual VectorHandler&
	ScalarMul(const VectorHandler& VH, const doublereal& d);

	/* Overload di += usato per la correzione della soluzione */
	virtual VectorHandler& operator += (const VectorHandler& VH);

	/* Overload di += usato per l'assemblaggio del residuo */
	virtual VectorHandler& operator += (const SubVectorHandler& SubVH);

	/* Overload di -= */
	virtual VectorHandler& operator -= (const VectorHandler& VH);

	/* Assegnazione che copia il contenuto della memoria di due handlers */
	virtual VectorHandler& operator = (const VectorHandler& VH);

	/* Norma 2 del vettore */
	virtual doublereal Dot(void) const;

	/* Norma del vettore */
	virtual doublereal Norm(void) const;

	/* Prodotto Scalare di due vettori */
	virtual doublereal InnerProd(const VectorHandler& VH) const;
};

/* VectorHandler - end */


/* MyVectorHandler - begin */

class MyVectorHandler : public VectorHandler {
	friend class SolutionManager;
	friend class Vec3;

protected:
	bool bOwnsMemory;

protected:
	integer iMaxSize;
	integer iCurSize;

	doublereal* pdVec;
	doublereal* pdVecm1;

public:
	MyVectorHandler(integer iSize = 0);

	MyVectorHandler(integer iSize, doublereal* pdTmpVec);

	virtual ~MyVectorHandler(void);

	virtual void Resize(integer iSize);

	void Detach(void);

	void Attach(integer iSize, doublereal* pd, integer iMSize = 0);

	/* Usata per il debug */
	virtual void IsValid(void) const;

	virtual inline doublereal* pdGetVec(void) const;

	virtual inline integer iGetSize(void) const;

	virtual void Reset(doublereal dResetVal = 0.);

	virtual inline void PutCoef(integer iRow, const doublereal& dCoef);

	virtual inline void IncCoef(integer iRow, const doublereal& dCoef);

	virtual inline void DecCoef(integer iRow, const doublereal& dCoef);

	virtual inline const doublereal& dGetCoef(integer iRow) const;

	/* Somma un Vec3 nella posizione desiderata */
	virtual void Add(integer iRow, const Vec3& v);

	/* Sottrae un Vec3 nella posizione desiderata */
	virtual void Sub(integer iRow, const Vec3& v);

	/* Scrive un Vec3 nella posizione desiderata */
	virtual void Put(integer iRow, const Vec3& v);

	/* Somma e moltiplica per uno scalare v = VH + d * VH1 */
	virtual VectorHandler&
	ScalarAddMul(const VectorHandler& VH, const VectorHandler& VH1,
			const doublereal& d);

	/* Somma e moltiplica per uno scalare */
	virtual VectorHandler&
	ScalarAddMul(const VectorHandler& VH, const doublereal& d);

	/* Moltiplica per uno scalare */
	virtual VectorHandler&
	ScalarMul(const VectorHandler& VH, const doublereal& d);

	/* Overload di += usato per la correzione della soluzione */
	virtual VectorHandler& operator += (const VectorHandler& VH);

	/* Overload di += usato per la correzione della soluzione */
	virtual MyVectorHandler& operator += (const MyVectorHandler& VH);

	/* Overload di -= */
	virtual VectorHandler& operator -= (const VectorHandler& VH);

	/* Overload di -= */
	virtual MyVectorHandler& operator -= (const MyVectorHandler& VH);

	/* Assegnazione che copia il contenuto della memoria di due handlers */
	virtual VectorHandler& operator = (const VectorHandler& VH);

	/* Assegnazione che copia il contenuto della memoria di due handlers */
	virtual MyVectorHandler& operator = (const MyVectorHandler& VH);

	/* Norma 2 del vettore */
	doublereal Dot(void) const;

	/* Norma del vettore */
	doublereal Norm(void) const;

	inline doublereal& operator () (integer iRow) const {
#ifdef DEBUG
		IsValid();
		ASSERT((iRow > 0) && (iRow <= iCurSize));
#endif

		return pdVecm1[iRow];
	};
};

inline doublereal*
MyVectorHandler::pdGetVec(void) const
{
#ifdef DEBUG
	IsValid();
#endif

	return pdVec;
}

inline integer
MyVectorHandler::iGetSize(void) const
{
	return iCurSize;
}

inline void
MyVectorHandler::PutCoef(integer iRow, const doublereal& dCoef)
{
	/* Nota: il flag di ritorno e' pleonastico. Lo si e' messo per
	 * analogia con le matrici sparse, in cui l'aggiunta
	 * di un coefficiente puo' risultare in un errore.
	 * Qui, per motivi di efficienza, il controllo sulla validita'
	 * dell'indice e del vettore viene svolto solo in debug */

#ifdef DEBUG
	IsValid();
	ASSERT((iRow > 0) && (iRow <= iCurSize));
#endif /* DEBUG */

	pdVecm1[iRow] = dCoef;

	return;
}

inline void
MyVectorHandler::IncCoef(integer iRow, const doublereal& dCoef)
{
	/* Vedi nota di PutCoef() */

#ifdef DEBUG
	IsValid();
	ASSERT((iRow > 0) && (iRow <= iCurSize));
#endif /* DEBUG */

	pdVecm1[iRow] += dCoef;

	return;
}

inline void
MyVectorHandler::DecCoef(integer iRow, const doublereal& dCoef)
{
	/* Vedi nota di PutCoef() */

#ifdef DEBUG
	IsValid();
	ASSERT((iRow > 0) && (iRow <= iCurSize));
#endif /* DEBUG */

	pdVecm1[iRow] -= dCoef;

	return;
}

inline const doublereal&
MyVectorHandler::dGetCoef(integer iRow) const
{
	/* Vedi nota di PutCoef() */

#ifdef DEBUG
	IsValid();
	ASSERT((iRow > 0) && (iRow <= iCurSize));
#endif

	return pdVecm1[iRow];
}

/* MyVectorHandler - end */


/* SolutionDataManager - begin */

class SolutionDataManager {
public:
	virtual ~SolutionDataManager(void);

	/* Collega il DataManager ed il DriveHandler ai vettori soluzione */
	virtual void
	LinkToSolution(const VectorHandler& XCurr,
			const VectorHandler& XPrimeCurr) = 0;

	/* Assembla il residuo */
	virtual void AssRes(VectorHandler &ResHdl, doublereal dCoef) = 0;
};

/* SolutionDataManager - end */


/* SolutionManager - begin */

class SolutionManager {
public:
	/* Distruttore: dealloca le matrici e distrugge gli oggetti propri */
	virtual ~SolutionManager(void);

	virtual void
	LinkToSolution(const VectorHandler& XCurr,
			const VectorHandler& XPrimeCurr);

	virtual void IsValid(void) const = 0;

	/* Inizializzatore generico */
	virtual void MatrInit(const doublereal& d = 0.) = 0;

	/* Risolve il sistema */
	virtual void Solve(void) = 0;

	/* sposta il puntatore al vettore del residuo */
	virtual void ChangeResPoint(doublereal* pRes) = 0;

	/* sposta il puntatore al vettore del residuo */
	virtual void ChangeSolPoint(doublereal* pSol) = 0;

	/* Rende disponibile l'handler per la matrice */
	virtual MatrixHandler* pMatHdl(void) const = 0;

	/* Rende disponibile l'handler per il termine noto */
	virtual VectorHandler* pResHdl(void) const = 0;

	/* Rende disponibile l'handler per la soluzione (e' lo stesso
	 * del termine noto, ma concettualmente sono separati) */
	virtual VectorHandler* pSolHdl(void) const = 0;
};

/* SolutionManager - end */

#endif /* SOLMAN_H */

