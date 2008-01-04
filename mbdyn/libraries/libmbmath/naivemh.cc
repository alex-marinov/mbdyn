/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2003-2008
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "myassert.h"
#include "mh.h"
#include "submat.h"
#include "naivemh.h"
#include "mthrdslv.h"

#undef CHECK_FOR_ZERO

/* NaiveMatrixHandler begin */

NaiveMatrixHandler::NaiveMatrixHandler(const integer n, 
	NaiveMatrixHandler *const nmh)
:  iSize(n), bOwnsMemory(true),
ppdRows(0), ppiRows(0), ppiCols(0), ppnonzero(0), piNzr(0), piNzc(0)
{
	if (nmh) {
		bOwnsMemory = false;
		iSize = nmh->iSize;
		ppdRows = nmh->ppdRows;
		ppiRows = nmh->ppiRows;
		ppiCols = nmh->ppiCols;
		ppnonzero = nmh->ppnonzero;
		piNzr = nmh->piNzr;
		piNzc = nmh->piNzc;
		
	} else {
		SAFENEWARR(ppiRows, integer *, iSize);
		ppiRows[0] = 0;
		SAFENEWARR(ppiRows[0], integer, iSize*iSize);

		SAFENEWARR(ppiCols, integer *, iSize);
		ppiCols[0] = 0;
		SAFENEWARR(ppiCols[0], integer, iSize*iSize);

		SAFENEWARR(ppdRows, doublereal *, iSize);
		ppdRows[0] = NULL;
		SAFENEWARR(ppdRows[0], doublereal, iSize*iSize);

		SAFENEWARR(ppnonzero, char *, iSize);
		ppnonzero[0] = NULL;
		SAFENEWARR(ppnonzero[0], char, iSize*iSize);

		for (integer i = 1; i < iSize; i++) {
			ppiRows[i] = ppiRows[i - 1] + iSize;
			ppiCols[i] = ppiCols[i - 1] + iSize;
			ppdRows[i] = ppdRows[i - 1] + iSize;
			ppnonzero[i] = ppnonzero[i - 1] + iSize;
		}

		SAFENEWARR(piNzr, integer, iSize);
		SAFENEWARR(piNzc, integer, iSize);

#ifdef HAVE_MEMSET
		memset(ppdRows[0], 0, sizeof(doublereal)*iSize*iSize);
		memset(ppnonzero[0], 0, sizeof(char)*iSize*iSize);
		memset(piNzr, 0, sizeof(integer)*iSize);
		memset(piNzc, 0, sizeof(integer)*iSize);
#else /* ! HAVE_MEMSET */
		for (integer row = 0; row < iSize; row++) {
			for (integer col = 0; col < iSize; col++) {
				ppnonzero[row][col] = 0;
			}
			piNzr[row] = 0;
			piNzc[row] = 0;
		}
#endif /* ! HAVE_MEMSET */		
	}
}

NaiveMatrixHandler::~NaiveMatrixHandler(void)
{
	if (bOwnsMemory) {
		if (ppiRows) {
			if (ppiRows[0]) {
				SAFEDELETEARR(ppiRows[0]);
			}
			SAFEDELETEARR(ppiRows);
		}

		if (ppiCols) {
			if (ppiCols[0]) {
				SAFEDELETEARR(ppiCols[0]);
			}
			SAFEDELETEARR(ppiCols);
		}

		if (ppdRows) {
			if (ppdRows[0]) {
				SAFEDELETEARR(ppdRows[0]);
			}
			SAFEDELETEARR(ppdRows);
		}

		if (ppnonzero) {
			if (ppnonzero[0]) {
				SAFEDELETEARR(ppnonzero[0]);
			}
			SAFEDELETEARR(ppnonzero);
		}

		if (piNzr) {
			SAFEDELETEARR(piNzr);
		}

		if (piNzc) {
			SAFEDELETEARR(piNzc);
		}
	}
}

void
NaiveMatrixHandler::Reset(void)
{
	for (integer row = 0; row < iSize; row++) {
		integer ncols = piNzc[row];
		integer *piCols = ppiCols[row];
		char *pnonzero = ppnonzero[row];
		for (integer col = 0; col < ncols; col++) {
			pnonzero[piCols[col]] = 0;
		}
		
		piNzr[row] = 0;
		piNzc[row] = 0;
	}
}

/* Overload di += usato per l'assemblaggio delle matrici */
MatrixHandler&
NaiveMatrixHandler::operator += (const SubMatrixHandler& SubMH)
{
	integer nr = SubMH.iGetNumRows();
	integer nc = SubMH.iGetNumCols();
	
	for (integer ir = 1; ir <= nr; ir++) {
		integer iRow = SubMH.iGetRowIndex(ir);

		for (integer ic = 1; ic <= nc; ic++) {
			doublereal d = SubMH(ir, ic);

#ifdef CHECK_FOR_ZERO
			if (d != 0.)
#endif /* CHECK_FOR_ZERO */
			{
				integer iCol = SubMH.iGetColIndex(ic);

				operator()(iRow, iCol) += d;
			}
		}
	}

	return *this;
}

/* Overload di -= usato per l'assemblaggio delle matrici */
MatrixHandler&
NaiveMatrixHandler::operator -= (const SubMatrixHandler& SubMH)
{
	integer nr = SubMH.iGetNumRows();
	integer nc = SubMH.iGetNumCols();

	for (integer ir = 1; ir <= nr; ir++) {
		integer iRow = SubMH.iGetRowIndex(ir);

		for (integer ic = 1; ic <= nc; ic++) {
			doublereal d = SubMH(ir, ic);

#ifdef CHECK_FOR_ZERO
			if (d != 0.)
#endif /* CHECK_FOR_ZERO */
			{
				integer iCol = SubMH.iGetColIndex(ic);

				operator()(iRow, iCol) -= d;
			}
		}
	}

	return *this;
}

/* Overload di += usato per l'assemblaggio delle matrici
 * questi li vuole ma non so bene perche'; force per la doppia
 * derivazione di VariableSubMatrixHandler */
MatrixHandler&
NaiveMatrixHandler::operator += (const VariableSubMatrixHandler& SubMH)
{
	switch (SubMH.eStatus) {
	case VariableSubMatrixHandler::FULL:
	{
		const FullSubMatrixHandler& SMH =
			*dynamic_cast<const FullSubMatrixHandler *>(&SubMH);
		/* NOTE: pirm1 is 1-based, for optimization purposes */
		integer *pirm1 = SMH.piRowm1;
		/* NOTE: pic is 0-based, for optimization purposes */
		integer *pic = SMH.piColm1 + 1;

		/* NOTE: ppd is 1-based for rows; access to SMH(iRow, iCol)
		 * results in ppd[iCol - 1][iRow] */
		doublereal **ppd = SMH.ppdCols;

		integer nr = SMH.iGetNumRows();
		integer nc = SMH.iGetNumCols();
		/* NOTE: iR is 1-based, for optimization purposes */
		for (integer iR = 1; iR <= nr; iR++) {
			integer iRow = pirm1[iR];

			/* NOTE: ic is 0-based, for optimization purposes */
			for (integer ic = 0; ic < nc; ic++) {
				doublereal d = ppd[ic][iR];

#ifdef CHECK_FOR_ZERO
				if (d != 0.)
#endif /* CHECK_FOR_ZERO */
				{
					integer iCol = pic[ic];

					operator()(iRow, iCol) += d;
				}
			}
		}
		break;
	}

	case VariableSubMatrixHandler::SPARSE:
	{
		const SparseSubMatrixHandler& SMH =
			*dynamic_cast<const SparseSubMatrixHandler *>(&SubMH);

		for (integer i = 1; i <= SMH.iNumItems; i++) {
			doublereal d = SMH.pdMatm1[i];

#ifdef CHECK_FOR_ZERO
			if (d != 0.)
#endif /* CHECK_FOR_ZERO */
			{
				integer iRow = SMH.piRowm1[i];
				integer iCol = SMH.piColm1[i];

				operator()(iRow, iCol) += d;
			}
		}
		break;
	}

	default:
		break;
	}

	return *this;
}

MatrixHandler&
NaiveMatrixHandler::operator -= (const VariableSubMatrixHandler& SubMH)
{
	switch (SubMH.eStatus) {
	case VariableSubMatrixHandler::FULL:
	{
		const FullSubMatrixHandler& SMH =
			*dynamic_cast<const FullSubMatrixHandler *>(&SubMH);
		integer *pirm1 = SMH.piRowm1;
		integer *pic = SMH.piColm1 + 1;
		doublereal **ppd = SMH.ppdCols;

		integer nr = SMH.iGetNumRows();
		integer nc = SMH.iGetNumCols();
		for (integer iR = 1; iR <= nr; iR++) {
			integer iRow = pirm1[iR];

			for (integer ic = 0; ic < nc; ic++) {
				doublereal d = ppd[ic][iR];

#ifdef CHECK_FOR_ZERO
				if (d != 0.)
#endif /* CHECK_FOR_ZERO */
				{
					integer iCol = pic[ic];

					operator()(iRow, iCol) -= d;
				}
			}
		}
		break;
	}

	case VariableSubMatrixHandler::SPARSE:
	{
		const SparseSubMatrixHandler& SMH =
			*dynamic_cast<const SparseSubMatrixHandler *>(&SubMH);

		for (integer i = 1; i <= SMH.iNumItems; i++) {
			doublereal d = SMH.pdMatm1[i];

#ifdef CHECK_FOR_ZERO
			if (d != 0.)
#endif /* CHECK_FOR_ZERO */
			{
				integer iRow = SMH.piRowm1[i];
				integer iCol = SMH.piColm1[i];

				operator()(iRow, iCol) -= d;
			}
		}
		break;
	}

	default:
		break;
	}

	return *this;
}

void NaiveMatrixHandler::MakeCCStructure(std::vector<integer>& Ai,
		std::vector<integer>& Ap) {
	integer nnz = 0;
	for (integer i = 0; i < iSize; i++) {
		nnz += piNzr[i];
	}
	Ai.resize(nnz);
	Ap.resize(iSize + 1);
	integer x_ptr = 0;
	for (integer col = 0; col < iSize; col++) {
		Ap[col] = x_ptr;
		integer nzr = piNzr[col];
		for (integer row = 0; row < nzr; row++) {
			Ai[x_ptr] = ppiRows[col][row];
			x_ptr++;
		}
	}
	Ap[iSize] = nnz;
};

/* NaiveMatrixHandler end */


/* NaivePermMatrixHandler begin */

NaivePermMatrixHandler::NaivePermMatrixHandler(integer iSize,
		const integer *const tperm)
: NaiveMatrixHandler(iSize), perm(tperm)
{
	NO_OP;
}

NaivePermMatrixHandler::NaivePermMatrixHandler(NaiveMatrixHandler*const nmh, 
		const integer *const tperm)
: NaiveMatrixHandler(0, nmh), perm(tperm)
{
	NO_OP;
}

NaivePermMatrixHandler::~NaivePermMatrixHandler(void)
{
	NO_OP;
}

const integer* const
NaivePermMatrixHandler::pGetPerm(void) const
{
	return perm;
}

/* NaivePermMatrixHandler end */





