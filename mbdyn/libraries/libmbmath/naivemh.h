/*
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2003-2004
 * 
 * This code is a partial merge of HmFe and MBDyn.
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 * Paolo Mantegazza     <mantegazza@aero.polimi.it>
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

#ifndef NAIVEMH_H
#define NAIVEMH_H

#include <vector>

#include "myassert.h"
#include "solman.h"

/* Sparse Matrix */
class NaiveMatrixHandler : public MatrixHandler {
protected:
	integer iSize;
	doublereal **ppdRows;
	integer **ppiRows, **ppiCols;
	integer *piNzr, *piNzc;

#ifdef DEBUG
	void IsValid(void) const {
		NO_OP;
	};
#endif /* DEBUG */

public:
	/* FIXME: always square? */
	NaiveMatrixHandler(const integer n);

	virtual ~NaiveMatrixHandler(void);

	integer iGetNumRows(void) const {
		return iSize;
	};

	integer iGetNumCols(void) const {
		return iSize;
	};

	void Init(const doublereal &r = 0.);

	/* Ridimensiona la matrice */
	virtual void Resize(integer, integer);

	/* Inserisce un coefficiente */
	virtual inline void
	PutCoef(integer iRow, integer iCol, const doublereal& dCoef);

	/* Incrementa un coefficiente - se non esiste lo crea */
	virtual inline void
	IncCoef(integer iRow, integer iCol, const doublereal& dCoef);

	/* Decrementa un coefficiente - se non esiste lo crea */
	virtual inline void
	DecCoef(integer iRow, integer iCol, const doublereal& dCoef);

	/* Restituisce un coefficiente - zero se non e' definito */
	virtual inline const doublereal&
	dGetCoef(integer iRow, integer iCol) const;

	virtual inline const doublereal&
	operator () (integer iRow, integer iCol) const;

	virtual inline doublereal&
	operator () (integer iRow, integer iCol);

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

};

/* Inserisce un coefficiente */
void
NaiveMatrixHandler::PutCoef(integer iRow, integer iCol,
		const doublereal& dCoef)
{
	if (dCoef == 0.) {
		return;
	}

	--iRow;
	--iCol;
	if (ppdRows[iRow][iCol] == 0.) {
		ppiRows[iRow][piNzr[iCol]] = iRow;
		ppiCols[iCol][piNzc[iRow]] = iCol;
		piNzr[iCol]++;
		piNzc[iRow]++;
	}

	ppdRows[iRow][iCol] = dCoef;
}

/* Incrementa un coefficiente - se non esiste lo crea */
void
NaiveMatrixHandler::IncCoef(integer iRow, integer iCol,
		const doublereal& dCoef)
{
	if (dCoef == 0.) {
		return;
	}

	--iRow;
	--iCol;
	if (ppdRows[iRow][iCol] == 0.) {
		ppdRows[iRow][iCol] = dCoef;
		ppiRows[iRow][piNzr[iCol]] = iRow;
		ppiCols[iCol][piNzc[iRow]] = iCol;
		piNzr[iCol]++;
		piNzc[iRow]++;
	} else {
		ppdRows[iRow][iCol] += dCoef;
	}
}

/* Decrementa un coefficiente - se non esiste lo crea */
void
NaiveMatrixHandler::DecCoef(integer iRow, integer iCol,
		const doublereal& dCoef)
{
	if (dCoef == 0.) {
		return;
	}

	--iRow;
	--iCol;
	if (ppdRows[iRow][iCol] == 0.) {
		ppdRows[iRow][iCol] = dCoef;
		ppiRows[iRow][piNzr[iCol]] = iRow;
		ppiCols[iCol][piNzc[iRow]] = iCol;
		piNzr[iCol]++;
		piNzc[iRow]++;
	} else {
		ppdRows[iRow][iCol] -= dCoef;
	}
}

/* Restituisce un coefficiente - zero se non e' definito */
const doublereal&
NaiveMatrixHandler::dGetCoef(integer iRow, integer iCol) const
{
	return ppdRows[--iRow][--iCol];
}

const doublereal&
NaiveMatrixHandler::operator () (integer iRow, integer iCol) const
{
	return ppdRows[--iRow][--iCol];
}

doublereal&
NaiveMatrixHandler::operator () (integer iRow, integer iCol)
{
	--iRow;
	--iCol;
	
	if (ppdRows[iRow][iCol] == 0.) {
		ppiRows[iRow][piNzr[iCol]] = iRow;
		ppiCols[iCol][piNzc[iRow]] = iCol;
		piNzr[iCol]++;
		piNzc[iRow]++;
	}

	return ppdRows[iRow][iCol];
}

#endif /* NAIVEMH_H */

