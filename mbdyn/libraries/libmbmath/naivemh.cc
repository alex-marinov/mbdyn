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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "mh.h"
#include "submat.h"
#include "naivemh.h"

NaiveMatrixHandler::NaiveMatrixHandler(const integer n)
:  iSize(n),
ppdRows(0), ppiRows(0), ppiCols(0), piNzr(0), piNzc(0)
{
	SAFENEWARR(ppiRows, integer *, iSize);
	ppiRows[0] = NULL;
	SAFENEWARR(ppiRows[0], integer, iSize*iSize);

	SAFENEWARR(ppiCols, integer *, iSize);
	ppiCols[0] = NULL;
	SAFENEWARR(ppiCols[0], integer, iSize*iSize);

	SAFENEWARR(ppdRows, doublereal *, iSize);
	ppdRows[0] = NULL;
	SAFENEWARR(ppdRows[0], doublereal, iSize*iSize);

	for (integer i = 1; i < iSize; i++) {
		ppiRows[i] = ppiRows[i-1] + iSize;
		ppiCols[i] = ppiCols[i-1] + iSize;
		ppdRows[i] = ppdRows[i-1] + iSize;
	}

	SAFENEWARR(piNzr, integer, iSize);
	SAFENEWARR(piNzc, integer, iSize);
}

NaiveMatrixHandler::~NaiveMatrixHandler(void)
{
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

	if (piNzr) {
		SAFEDELETEARR(piNzr);
	}

	if (piNzc) {
		SAFEDELETEARR(piNzc);
	}
}

void
NaiveMatrixHandler::Init(const doublereal& c)
{
#ifdef HAVE_MEMSET_H
	memset(ppiRows[0], 0, sizeof(integer)*iSize*iSize);
	memset(piNzr, 0, sizeof(integer)*iSize);
	memset(piNzc, 0, sizeof(integer)*iSize);
#else /* ! HAVE_MEMSET_H */
	for (integer i = 0; i < iSize*iSize; i++) {
		ppdRows[0][i] = 0.;
	}

	for (integer i = 0; i < iSize; i++) {
		piNzr[i] = 0;
		piNzc[i] = 0;
	}
#endif /* ! HAVE_MEMSET_H */
}

/* Overload di += usato per l'assemblaggio delle matrici */
MatrixHandler&
NaiveMatrixHandler::operator += (const SubMatrixHandler& SubMH)
{
	for (integer ir = 0; ir < SubMH.iGetNumRows(); ir++) {
		integer iRow = SubMH.iGetRowIndex(ir + 1) - 1;

		for (integer ic = 0; ic < SubMH.iGetNumRows(); ic++) {
			doublereal d = SubMH(ir + 1, ic + 1);

			if (d != 0.) {
				integer iCol = SubMH.iGetColIndex(ic + 1) - 1;
				if (ppdRows[iRow][iCol] == 0.) {
					ppdRows[iRow][iCol] = d;

					ppiRows[iCol][piNzr[iCol]] = iRow;
					ppiCols[iRow][piNzc[iRow]] = iCol;
					piNzr[iCol]++;
					piNzc[iRow]++;

				} else {
					ppdRows[iRow][iCol] += d;
				}
			}
		}
	}

	return *this;
}

/* Overload di -= usato per l'assemblaggio delle matrici */
MatrixHandler&
NaiveMatrixHandler::operator -= (const SubMatrixHandler& SubMH)
{
	for (integer ir = 0; ir < SubMH.iGetNumRows(); ir++) {
		integer iRow = SubMH.iGetRowIndex(ir + 1) - 1;

		for (integer ic = 0; ic < SubMH.iGetNumRows(); ic++) {
			doublereal d = SubMH(ir + 1, ic + 1);

			if (d != 0.) {
				integer iCol = SubMH.iGetColIndex(ic + 1) - 1;
				if (ppdRows[iRow][iCol] == 0.) {
					ppdRows[iRow][iCol] = d;

					ppiRows[iCol][piNzr[iCol]] = iRow;
					ppiCols[iRow][piNzc[iRow]] = iCol;
					piNzr[iCol]++;
					piNzc[iRow]++;

				} else {
					ppdRows[iRow][iCol] -= d;
				}
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
		integer *pirm1 = SMH.piRowm1;
		integer *picm1 = SMH.piColm1;
		doublereal **ppd = SMH.ppdCols;

		for (integer ir = 0; ir < SMH.iGetNumRows(); ir++) {
			integer iRow = pirm1[ir + 1] - 1;

			for (integer ic = 0; ic < SMH.iGetNumRows(); ic++) {
				doublereal d = ppd[ic][ir];

				if (d != 0.) {
					integer iCol = picm1[ic + 1] - 1;
					if (ppdRows[iRow][iCol] == 0.) {
						ppdRows[iRow][iCol] = d;

						ppiRows[iRow][piNzr[iCol]] = iRow;
						ppiCols[iRow][piNzc[iRow]] = iCol;
						piNzr[iCol]++;
						piNzc[iRow]++;
					} else {
						ppdRows[iRow][iCol] += d;
					}
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

			if (d != 0.) {
				integer iRow = SMH.piRowm1[i] - 1;
				integer iCol = SMH.piColm1[i] - 1;

				if (ppdRows[iRow][iCol] == 0.) {
					ppdRows[iRow][iCol] = d;

					ppiRows[iRow][piNzr[iCol]] = iRow;
					ppiCols[iRow][piNzc[iRow]] = iCol;
					piNzr[iCol]++;
					piNzc[iRow]++;

				} else {
					ppdRows[iRow][iCol] += d;
				}
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
		integer *picm1 = SMH.piColm1;
		doublereal **ppd = SMH.ppdCols;

		for (integer ir = 0; ir < SMH.iGetNumRows(); ir++) {
			integer iRow = pirm1[ir + 1] - 1;

			for (integer ic = 0; ic < SMH.iGetNumRows(); ic++) {
				doublereal d = ppd[ic][ir];

				if (d != 0.) {
					integer iCol = picm1[ic + 1] - 1;
					if (ppdRows[iRow][iCol] == 0.) {
						ppdRows[iRow][iCol] = d;

						ppiRows[iRow][piNzr[iCol]] = iRow;
						ppiCols[iRow][piNzc[iRow]] = iCol;
						piNzr[iCol]++;
						piNzc[iRow]++;
						
					} else {
						ppdRows[iRow][iCol] -= d;
					}
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

			if (d != 0.) {
				integer iRow = SMH.piRowm1[i] - 1;
				integer iCol = SMH.piColm1[i] - 1;

				if (ppdRows[iRow][iCol] == 0.) {
					ppdRows[iRow][iCol] = d;

					ppiRows[iRow][piNzr[iCol]] = iRow;
					ppiCols[iRow][piNzc[iRow]] = iCol;
					piNzr[iCol]++;
					piNzc[iRow]++;

				} else {
					ppdRows[iRow][iCol] -= d;
				}
			}
		}
		break;
	}

	default:
		break;
	}

	return *this;
}

