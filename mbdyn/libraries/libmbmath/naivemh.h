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
#include <limits>

#include "myassert.h"
#include "solman.h"

class NaiveSolver;

/* Sparse Matrix */
class NaiveMatrixHandler : public MatrixHandler {
protected:
	friend class NaiveSolver;
	integer iSize;
	doublereal **ppdRows;
	integer **ppiRows, **ppiCols;
	char **ppnonzero;
	integer *piNzr, *piNzc;

#ifdef DEBUG
	void IsValid(void) const {
		NO_OP;
	};
#endif /* DEBUG */

public:
	/* FIXME: always square? yes! */
	NaiveMatrixHandler(const integer n);

	virtual ~NaiveMatrixHandler(void);

	integer iGetNumRows(void) const {
		return iSize;
	};

	integer iGetNumCols(void) const {
		return iSize;
	};

	void Reset(const doublereal r = 0.);
	void Init(const doublereal d = 0.) {
		Reset(d);
	};

	/* Ridimensiona la matrice */
	virtual void Resize(integer, integer) {
		THROW(ErrGeneric());
	};

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


const doublereal&
NaiveMatrixHandler::operator () (integer iRow, integer iCol) const
{
	--iRow;
	--iCol;
	if (ppnonzero[iRow][iCol]) {
		return ppdRows[iRow][iCol];
	}
	return ::dZero;
}

doublereal&
NaiveMatrixHandler::operator () (integer iRow, integer iCol)
{
	--iRow;
	--iCol;
	if (!(ppnonzero[iRow][iCol])) {
		ppnonzero[iRow][iCol] = 1;
		ppiRows[iCol][piNzr[iCol]] = iRow;
		ppiCols[iRow][piNzc[iRow]] = iCol;
		piNzr[iCol]++;
		piNzc[iRow]++;
		ppdRows[iRow][iCol] = 0.;
	}

	return ppdRows[iRow][iCol];
}

#endif /* NAIVEMH_H */

