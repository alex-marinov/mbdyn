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
class MultiThreadDataManager;

/* Sparse Matrix */
class NaiveMatrixHandler : public MatrixHandler {
protected:
	friend void* sum_naive_matrices(void* arg);
	friend class NaiveSolver;
	friend class ParNaiveSolver;
	friend class MultiThreadDataManager;
	integer iSize;
	bool bOwnsMemory;
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
	NaiveMatrixHandler(const integer n, NaiveMatrixHandler *const nmh = 0);

	virtual ~NaiveMatrixHandler(void);

	integer iGetNumRows(void) const {
		return iSize;
	};

	integer iGetNumCols(void) const {
		return iSize;
	};

	void Reset(void);

	/* Ridimensiona la matrice */
	virtual void Resize(integer, integer) {
		throw ErrGeneric();
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
	
	void MakeCCStructure(std::vector<integer>& Ai,
		std::vector<integer>& Ap);

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

/* Sparse Matrix with unknowns permutation*/
class NaivePermMatrixHandler : public NaiveMatrixHandler {
protected:
	const integer* const perm;

#ifdef DEBUG
	void IsValid(void) const {
		NO_OP;
	};
#endif /* DEBUG */

public:
	/* FIXME: always square? yes! */
	NaivePermMatrixHandler(integer iSize,
		const integer *const tperm);

	NaivePermMatrixHandler(NaiveMatrixHandler*const nmh, 
		const integer *const tperm);

	virtual ~NaivePermMatrixHandler(void);

	const integer* const pGetPerm(void) const;

	virtual inline const doublereal&
	operator () (integer iRow, integer iCol) const {
		iCol = perm[iCol - 1] + 1;
		return NaiveMatrixHandler::operator()(iRow, iCol);
	};

	virtual inline doublereal&
	operator () (integer iRow, integer iCol) {
		iCol = perm[iCol - 1] + 1;
		return NaiveMatrixHandler::operator()(iRow, iCol);
	};


};


#endif /* NAIVEMH_H */

