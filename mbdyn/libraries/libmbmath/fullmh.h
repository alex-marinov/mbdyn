/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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

/* FullMatrixHandler - begin */

class FullMatrixHandler : public MatrixHandler {
	friend std::ostream&
	operator << (std::ostream& out, const FullMatrixHandler& m);
	friend class FullSubMatrixHandler;
	friend class SparseSubMatrixHandler;

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
	doublereal* pdGetMat(void) const {
		ASSERT(pdRaw != NULL);
		return pdRaw;
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
		ASSERT(iRow > 0 && iRow <= iNumRows);
		ASSERT(iCol > 0 && iCol <= iNumCols);
#endif /* DEBUG */

		return ppdColsm1[iCol][iRow];
	};

	virtual const doublereal&
	operator () (integer iRow, integer iCol) const {
#ifdef DEBUG
		IsValid();
		ASSERT(iRow > 0 && iRow <= iNumRows);
		ASSERT(iCol > 0 && iCol <= iNumCols);
#endif /* DEBUG */

		return ppdColsm1[iCol][iRow];
	};

	/* Overload di += usato per l'assemblaggio delle matrici */
	virtual MatrixHandler& operator +=(const SubMatrixHandler& SubMH);

	/* Overload di -= usato per l'assemblaggio delle matrici */
	virtual MatrixHandler& operator -=(const SubMatrixHandler& SubMH);

	/* Esegue il prodotto tra due matrici e se lo memorizza */
	void MatMul(const FullMatrixHandler& m1, const FullMatrixHandler& m2);
};

extern std::ostream&
operator << (std::ostream& out, const FullMatrixHandler& m);

/* FullMatrixHandler - end */

#endif /* FULLMH_H */

