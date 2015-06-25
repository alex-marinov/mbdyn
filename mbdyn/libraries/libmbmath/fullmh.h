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

#ifndef FULLMH_H
#define FULLMH_H

#include "solman.h"
#include "except.h"

#include "spmh.h"
#include "matvec3.h"
#include "matvec3n.h"

/* FullMatrixHandler - begin */

class FullMatrixHandler : public MatrixHandler {
	friend std::ostream&
	operator << (std::ostream& out, const FullMatrixHandler& m);
	friend std::ostream& Write(std::ostream& out,
		const FullMatrixHandler& m,
		const char* s, const char* s2);
	friend class FullSubMatrixHandler;
	friend class SparseSubMatrixHandler;

public:
	// allow copy constructor!
	FullMatrixHandler(const FullMatrixHandler&);
	FullMatrixHandler& operator = (const FullMatrixHandler&);

private:

protected:
	bool bOwnsMemory;

	integer iNumRows;
	integer iNumCols;

	integer iRawSize;
	integer iMaxCols;

	doublereal* pdRaw;
	doublereal* pdRawm1;
	doublereal** ppdCols;
	doublereal** ppdColsm1;

	void CreateColRow(integer iNR, integer iNC);

public:
	class const_iterator {
		friend class FullMatrixHandler;

	private:
		const FullMatrixHandler& m;
		mutable integer i_idx;
		mutable SparseMatrixHandler::SparseMatrixElement elem;

	protected:
		void reset(bool is_end = false);

	public:
		const_iterator(const FullMatrixHandler& m, bool is_end = false);
		~const_iterator(void);
		const FullMatrixHandler::const_iterator& operator ++ (void) const;
		const SparseMatrixHandler::SparseMatrixElement* operator -> (void) const;
		const SparseMatrixHandler::SparseMatrixElement& operator * (void) const;
		bool operator == (const FullMatrixHandler::const_iterator& op) const;
		bool operator != (const FullMatrixHandler::const_iterator& op) const;
	};

private:
	const_iterator m_end;

public:
	FullMatrixHandler::const_iterator begin(void) const {
		return const_iterator(*this);
	};

	const FullMatrixHandler::const_iterator& end(void) const {
		return m_end;
	};

public:
	FullMatrixHandler(doublereal* pd, doublereal** ppd,
			integer iSize, integer iNR, integer iNC,
			integer iMaxCols = 0);

	/* costruttore che si alloca la memoria */
	FullMatrixHandler(integer iNR, integer iNC = 0);

	/* costruttore che non fa nulla */
	FullMatrixHandler(void);

	virtual ~FullMatrixHandler(void);

	/* ridimensiona la matrice (se possiede la memoria) */
	virtual void Resize(integer iNewRows, integer iNewCols);

	/* si stacca dalla memoria a cui e' associato */
	void Detach(void);
	
	/* zero the matrix */
	void Reset(void);

	/* Attacca un nuovo array, con n. righe, n. colonne e dim. massima;
	 * se assente, assunta = nrighe*ncolonne */
	void Attach(integer iNewRows, integer iNewCols,
			doublereal* pd, doublereal** ppd,
			integer iMSize = 0, integer iMaxC = 0);

#ifdef DEBUG
	/* Usata per il debug */
	virtual void IsValid(void) const;
#endif /* DEBUG */

	/* Used to access raw data by c functions */
	const doublereal* pdGetMat(void) const {
		ASSERT(pdRaw != NULL);
		return pdRaw;
	};

	doublereal* pdGetMat(void) {
		ASSERT(pdRaw != NULL);
		return pdRaw;
	};

	virtual inline const doublereal *pdGetVec(integer iCol) const {
#ifdef DEBUG
		ASSERT(pdRaw != NULL);
		ASSERT(ppdColsm1 != NULL);
		ASSERT(iCol > 0 && iCol <= iNumCols);
		ASSERT(ppdColsm1[iCol] != NULL);
#endif /* DEBUG */

		return &ppdColsm1[iCol][1];
	};

	virtual inline doublereal *pdGetVec(integer iCol) {
#ifdef DEBUG
		ASSERT(pdRaw != NULL);
		ASSERT(ppdColsm1 != NULL);
		ASSERT(iCol > 0 && iCol <= iNumCols);
		ASSERT(ppdColsm1[iCol] != NULL);
#endif /* DEBUG */

		return &ppdColsm1[iCol][1];
	};

	/* Inserisce un coefficiente */
	virtual inline void
	PutCoef(integer iRow, integer iCol, const doublereal& dCoef) {
#ifdef DEBUG
		IsValid();
		ASSERT(iRow > 0 && iRow <= iNumRows);
		ASSERT(iCol > 0 && iCol <= iNumCols);
#endif /* DEBUG */

		ppdColsm1[iCol][iRow] = dCoef;
	};

	/* Incrementa un coefficiente - se non esiste lo crea */
	virtual inline void
	IncCoef(integer iRow, integer iCol, const doublereal& dCoef) {
#ifdef DEBUG
		IsValid();
		ASSERT(iRow > 0 && iRow <= iNumRows);
		ASSERT(iCol > 0 && iCol <= iNumCols);
#endif /* DEBUG */

		ppdColsm1[iCol][iRow] += dCoef;
	};

	/* Incrementa un coefficiente - se non esiste lo crea */
	virtual inline void
	DecCoef(integer iRow, integer iCol, const doublereal& dCoef) {
#ifdef DEBUG
		IsValid();
		ASSERT(iRow > 0 && iRow <= iNumRows);
		ASSERT(iCol > 0 && iCol <= iNumCols);
#endif /* DEBUG */

		ppdColsm1[iCol][iRow] -= dCoef;
	};

	/* Restituisce un coefficiente - zero se non e' definito */
	virtual inline const doublereal&
	dGetCoef(integer iRow, integer iCol) const {
#ifdef DEBUG
		IsValid();
		ASSERT(iRow > 0 && iRow <= iNumRows);
		ASSERT(iCol > 0 && iCol <= iNumCols);
#endif /* DEBUG */

		return ppdColsm1[iCol][iRow];
	};

	/* dimensioni */
	virtual integer iGetNumRows(void) const {
		return iNumRows;
	};

	virtual integer iGetNumCols(void) const {
		return iNumCols;
	};

	virtual doublereal&
	operator () (integer iRow, integer iCol) {
#ifdef DEBUG
		IsValid();
		ASSERT(iRow > 0);
		ASSERT(iRow <= iNumRows);
		ASSERT(iCol > 0);
		ASSERT(iCol <= iNumCols);
#endif /* DEBUG */

		return ppdColsm1[iCol][iRow];
	};

	virtual const doublereal&
	operator () (integer iRow, integer iCol) const {
#ifdef DEBUG
		IsValid();
		ASSERT(iRow > 0);
		ASSERT(iRow <= iNumRows);
		ASSERT(iCol > 0);
		ASSERT(iCol <= iNumCols);
#endif /* DEBUG */

		return ppdColsm1[iCol][iRow];
	};

	/* Overload di += usato per l'assemblaggio delle matrici */
	virtual MatrixHandler& operator +=(const SubMatrixHandler& SubMH);

	/* Overload di -= usato per l'assemblaggio delle matrici */
	virtual MatrixHandler& operator -=(const SubMatrixHandler& SubMH);

	/* Overload di += usato per l'assemblaggio delle matrici */
	virtual MatrixHandler& operator +=(const VariableSubMatrixHandler& SubMH);

	/* Overload di -= usato per l'assemblaggio delle matrici */
	virtual MatrixHandler& operator -=(const VariableSubMatrixHandler& SubMH);

	/* Esegue il prodotto tra due matrici e se lo memorizza */
	void MatMul(const FullMatrixHandler& m1, const FullMatrixHandler& m2);

protected:
	MatrixHandler&
	MatMatMul_base(void (MatrixHandler::*op)(integer iRow, integer iCol,
				const doublereal& dCoef),
			MatrixHandler& out, const MatrixHandler& in) const;
	MatrixHandler&
	MatTMatMul_base(void (MatrixHandler::*op)(integer iRow, integer iCol,
				const doublereal& dCoef),
			MatrixHandler& out, const MatrixHandler& in) const;

	virtual VectorHandler&
	MatVecMul_base(void (VectorHandler::*op)(integer iRow,
				const doublereal& dCoef),
			VectorHandler& out, const VectorHandler& in) const;
	virtual VectorHandler&
	MatTVecMul_base(void (VectorHandler::*op)(integer iRow,
				const doublereal& dCoef),
			VectorHandler& out, const VectorHandler& in) const;

public:
	void Add(integer iRow,  integer iCol,
		const FullMatrixHandler & source);
	void Sub(integer iRow,  integer iCol,
		const FullMatrixHandler & source);
	void Put(integer iRow,  integer iCol,
		const FullMatrixHandler & source);
	void Add(integer iRow,  integer iCol,
		const FullMatrixHandler & source, const doublereal dCoef);
	void Sub(integer iRow,  integer iCol,
		const FullMatrixHandler & source, const doublereal dCoef);
	void Put(integer iRow,  integer iCol,
		const FullMatrixHandler & source, const doublereal dCoef);
	void AddT(integer iRow,  integer iCol,
		const FullMatrixHandler & source);
	void SubT(integer iRow,  integer iCol,
		const FullMatrixHandler & source);
	void PutT(integer iRow,  integer iCol,
		const FullMatrixHandler & source);
	void AddT(integer iRow,  integer iCol,
		const FullMatrixHandler & source, const doublereal dCoef);
	void SubT(integer iRow,  integer iCol,
		const FullMatrixHandler & source, const doublereal dCoef);
	void PutT(integer iRow,  integer iCol,
		const FullMatrixHandler & source, const doublereal dCoef);

	/*
	 * Somma un vettore di tipo Vec3 in una data posizione.
	 * Nota: si assume che nella sottomatrice vi sia spazio
	 * per il vettore 3.
	 * Nota: gli indici sono a base 1, in stile FORTRAN.
	 * @param iRow indice di riga della sottomatrice da cui iniziare
	 * @param iCol indice di colonna della sottomatrice da cui iniziare
	 * @param v    Vec3 da sommare
	 */
	void
	Add(integer iRow, integer iCol, const Vec3& v);

	/*
	 * Sottrae un vettore di tipo Vec3 in una data posizione.
	 * Nota: si assume che nella sottomatrice vi sia spazio
	 * per il vettore 3.
	 * Nota: gli indici sono a base 1, in stile FORTRAN.
	 * @param iRow indice di riga della sottomatrice da cui iniziare
	 * @param iCol indice di colonna della sottomatrice da cui iniziare
	 * @param v    Vec3 da sommare
	 */
	void
	Sub(integer iRow, integer iCol, const Vec3& v);

	/*
	 * Scrive un vettore di tipo Vec3 in una data posizione.
	 * Nota: si assume che nella sottomatrice vi sia spazio
	 * per il vettore 3.
	 * Nota: gli indici sono a base 1, in stile FORTRAN.
	 * @param iRow indice di riga della sottomatrice da cui iniziare
	 * @param iCol indice di colonna della sottomatrice da cui iniziare
	 * @param v    Vec3 da sommare
	 */
	void
	Put(integer iRow, integer iCol, const Vec3& v);

	/*
	 * Somma un vettore di tipo Vec3 trasposto in una data posizione.
	 * Nota: si assume che nella sottomatrice vi sia spazio
	 * per il vettore 3 trasposto (riga).
	 * Nota: gli indici sono a base 1, in stile FORTRAN.
	 * @param iRow indice di riga della sottomatrice da cui iniziare
	 * @param iCol indice di colonna della sottomatrice da cui iniziare
	 * @param v    Vec3 da sommare
	 */
	void
	AddT(integer iRow, integer iCol, const Vec3& v);

	/*
	 * Sottrae un vettore di tipo Vec3 trasposto in una data posizione.
	 * Nota: si assume che nella sottomatrice vi sia spazio
	 * per il vettore 3 trasposto (riga).
	 * Nota: gli indici sono a base 1, in stile FORTRAN.
	 * @param iRow indice di riga della sottomatrice da cui iniziare
	 * @param iCol indice di colonna della sottomatrice da cui iniziare
	 * @param v    Vec3 da sommare
	 */
	void
	SubT(integer iRow, integer iCol, const Vec3& v);

	/*
	 * Scrive un vettore di tipo Vec3 trasposto in una data posizione.
	 * Nota: si assume che nella sottomatrice vi sia spazio
	 * per il vettore 3 trasposto (riga).
	 * Nota: gli indici sono a base 1, in stile FORTRAN.
	 * @param iRow indice di riga della sottomatrice da cui iniziare
	 * @param iCol indice di colonna della sottomatrice da cui iniziare
	 * @param v    Vec3 da sommare
	 */
	void
	PutT(integer iRow, integer iCol, const Vec3& v);

#if 0 /* FIXME: replace original? */
	/*
	 * Somma un vettore di tipo Vec3 in una data posizione in diagonale.
	 * Nota: si assume che nella sottomatrice vi sia spazio
	 * per il vettore 3.
	 * Nota: gli indici sono a base 1, in stile FORTRAN.
	 * @param iRow indice di riga della sottomatrice da cui iniziare
	 * @param iCol indice di colonna della sottomatrice da cui iniziare
	 * @param v    Vec3 da sommare
	 */
	void
	AddDiag(integer iRow, integer iCol, const Vec3& v);

	/*
	 * Sottrae un vettore di tipo Vec3 in una data posizione in diagonale.
	 * Nota: si assume che nella sottomatrice vi sia spazio
	 * per il vettore 3.
	 * Nota: gli indici sono a base 1, in stile FORTRAN.
	 * @param iRow indice di riga della sottomatrice da cui iniziare
	 * @param iCol indice di colonna della sottomatrice da cui iniziare
	 * @param v    Vec3 da sommare
	 */
	void
	SubDiag(integer iRow, integer iCol, const Vec3& v);

	/*
	 * Scrive un vettore di tipo Vec3 in una data posizione in diagonale.
	 * Nota: si assume che nella sottomatrice vi sia spazio
	 * per il vettore 3.
	 * Nota: gli indici sono a base 1, in stile FORTRAN.
	 * @param iRow indice di riga della sottomatrice da cui iniziare
	 * @param iCol indice di colonna della sottomatrice da cui iniziare
	 * @param v    Vec3 da sommare
	 */
	void
	PutDiag(integer iRow, integer iCol, const Vec3& v);

	/*
	 * Somma un vettore di tipo Vec3 in una data posizione [ v x ].
	 * Nota: si assume che nella sottomatrice vi sia spazio
	 * per il vettore 3.
	 * Nota: gli indici sono a base 1, in stile FORTRAN.
	 * @param iRow indice di riga della sottomatrice da cui iniziare
	 * @param iCol indice di colonna della sottomatrice da cui iniziare
	 * @param v    Vec3 da sommare
	 */
	void
	AddCross(integer iRow, integer iCol, const Vec3& v);

	/*
	 * Sottrae un vettore di tipo Vec3 in una data posizione [ v x ].
	 * Nota: si assume che nella sottomatrice vi sia spazio
	 * per il vettore 3.
	 * Nota: gli indici sono a base 1, in stile FORTRAN.
	 * @param iRow indice di riga della sottomatrice da cui iniziare
	 * @param iCol indice di colonna della sottomatrice da cui iniziare
	 * @param v    Vec3 da sommare
	 */
	void
	SubCross(integer iRow, integer iCol, const Vec3& v);

	/*
	 * Scrive un vettore di tipo Vec3 in una data posizione [ v x ].
	 * Nota: si assume che nella sottomatrice vi sia spazio
	 * per il vettore 3.
	 * Nota: gli indici sono a base 1, in stile FORTRAN.
	 * @param iRow indice di riga della sottomatrice da cui iniziare
	 * @param iCol indice di colonna della sottomatrice da cui iniziare
	 * @param v    Vec3 da sommare
	 */
	void
	PutCross(integer iRow, integer iCol, const Vec3& v);
#endif

	/*
	 * Somma una matrice di tipo Mat3x3 in una data posizione.
	 * Nota: si assume che nella sottomatrice vi sia spazio
	 * per la matrice 3x3.
	 * Nota: gli indici sono a base 1, in stile FORTRAN.
	 * @param iRow indice di riga della sottomatrice da cui iniziare
	 * @param iCol indice di colonna della sottomatrice da cui iniziare
	 * @param m    Mat3x3 da sommare
	 */
	void
	Add(integer iRow, integer iCol, const Mat3x3& m);
	void
	AddT(integer iRow, integer iCol, const Mat3x3& m);

	/*
	 * Sottrae una matrice di tipo Mat3x3 da una data posizione.
	 * Nota: si assume che nella sottomatrice vi sia spazio
	 * per la matrice 3x3.
	 * Nota: gli indici sono a base 1, in stile FORTRAN.
	 * @param iRow indice di riga della sottomatrice da cui iniziare
	 * @param iCol indice di colonna della sottomatrice da cui iniziare
	 * @param m    Mat3x3 da sottrarre
	 */
	void
	Sub(integer iRow, integer iCol, const Mat3x3& m);
	void
	SubT(integer iRow, integer iCol, const Mat3x3& m);

	/*
	 * Scrive una matrice di tipo Mat3x3 in una data posizione.
	 * Nota: si assume che nella sottomatrice vi sia spazio
	 * per la matrice 3x3.
	 * Nota: gli indici sono a base 1, in stile FORTRAN.
	 * @param iRow indice di riga della sottomatrice da cui iniziare
	 * @param iCol indice di colonna della sottomatrice da cui iniziare
	 * @param m    Mat3x3 da scrivere
	 */
	void
	Put(integer iRow, integer iCol, const Mat3x3& m);
	void
	PutT(integer iRow, integer iCol, const Mat3x3& m);

	/*
	 * Somma una matrice di tipo Mat3xN in una data posizione.
	 * Nota: si assume che nella sottomatrice vi sia spazio
	 * per la matrice 3x3.
	 * Nota: gli indici sono a base 1, in stile FORTRAN.
	 * @param iRow indice di riga della sottomatrice da cui iniziare
	 * @param iCol indice di colonna della sottomatrice da cui iniziare
	 * @param m    Mat3xN
	 */
	void
	Add(integer iRow, integer iCol, const Mat3xN& m);
	void
	AddT(integer iRow, integer iCol, const Mat3xN& m);

	/*
	 * Sottrae una matrice di tipo Mat3xN da una data posizione.
	 * Nota: si assume che nella sottomatrice vi sia spazio
	 * per la matrice 3x3.
	 * Nota: gli indici sono a base 1, in stile FORTRAN.
	 * @param iRow indice di riga della sottomatrice da cui iniziare
	 * @param iCol indice di colonna della sottomatrice da cui iniziare
	 * @param m    Mat3xN
	 */
	void
	Sub(integer iRow, integer iCol, const Mat3xN& m);
	void
	SubT(integer iRow, integer iCol, const Mat3xN& m);

	/*
	 * Scrive una matrice di tipo Mat3xN in una data posizione.
	 * Nota: si assume che nella sottomatrice vi sia spazio
	 * per la matrice 3x3.
	 * Nota: gli indici sono a base 1, in stile FORTRAN.
	 * @param iRow indice di riga della sottomatrice da cui iniziare
	 * @param iCol indice di colonna della sottomatrice da cui iniziare
	 * @param m    Mat3xN
	 */
	void
	Put(integer iRow, integer iCol, const Mat3xN& m);
	void
	PutT(integer iRow, integer iCol, const Mat3xN& m);

	/* come sopra, ma per matrici Nx3 **/
	void Add(integer iRow, integer iCol, const MatNx3& m);
	void Sub(integer iRow, integer iCol, const MatNx3& m);
	void Put(integer iRow, integer iCol, const MatNx3& m);

	void CopyMatrixRow(integer dest_row,
		const FullMatrixHandler & source, integer source_row);
	void CopyMatrixBlock(integer dest_row, integer dest_col,
		const FullMatrixHandler & source, 
		integer source_start_row, integer source_end_row,
		integer source_start_col, integer source_end_col);
};

extern std::ostream&
operator << (std::ostream& out, const FullMatrixHandler& m);

extern std::ostream& Write(std::ostream& out,
	const FullMatrixHandler& m,
	const char* s = " ", 
	const char* s2 = NULL);

/* FullMatrixHandler - end */

#endif /* FULLMH_H */

