/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2003-2015
 * 
 * This code is a partial merge of HmFe and MBDyn.
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 * Paolo Mantegazza     <mantegazza@aero.polimi.it>
 * Marco Morandini  <morandini@aero.polimi.it>
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "spmapmh.h"
#include "ccmh.h"

template <int off>
CColMatrixHandler<off>::CColMatrixHandler(std::vector<doublereal>& x,
		const std::vector<integer>& i,
		const std::vector<integer>& p)
: CompactSparseMatrixHandler_tpl<off>(p.size() - 1, p.size() - 1, x, i, p)
{
	NO_OP;
}

template <int off>
CColMatrixHandler<off>::~CColMatrixHandler()
{
	NO_OP;
}

/* used by MultiThreadDataManager to duplicate the storage array
 * while preserving the CC indices */
template <int off>
CompactSparseMatrixHandler *
CColMatrixHandler<off>::Copy(void) const
{
	std::vector<doublereal> *pax =
		new std::vector<doublereal>(CompactSparseMatrixHandler_tpl<off>::Ax);
	CColMatrixHandler<off> *p =
		new CColMatrixHandler<off>(*pax, CompactSparseMatrixHandler_tpl<off>::Ai,
			CompactSparseMatrixHandler_tpl<off>::Ap);
	p->bMatDuplicate = true;

	return p;
}

template <int off>
void
CColMatrixHandler<off>::Resize(integer n, integer nn)
{
	silent_cerr("CColMatrixHandler<off>::Resize called" << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

/* Estrae una colonna da una matrice */
template <int off>
VectorHandler&
CColMatrixHandler<off>::GetCol(integer icol, VectorHandler& out) const
{
	// NOTE: out must be zeroed by caller
	
        if (icol > SparseMatrixHandler::iGetNumCols()) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	integer idx = CompactSparseMatrixHandler_tpl<off>::Ap[icol - 1] - off;
	integer idxe = CompactSparseMatrixHandler_tpl<off>::Ap[icol] - off;

	for ( ; idx < idxe; idx++) {
		out(CompactSparseMatrixHandler_tpl<off>::Ai[idx] - off + 1) =
			CompactSparseMatrixHandler_tpl<off>::Ax[idx];
	}

	return out;
}

/* Moltiplica per uno scalare e somma a una matrice */
template <int off>
MatrixHandler&
CColMatrixHandler<off>::MulAndSumWithShift(MatrixHandler& out, doublereal s,
		integer drow, integer dcol) const
{
	silent_cerr("CColMatrixHandler<off>::MulAndSumWithShift called"
			<< std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);		
	if ((out.iGetNumCols() < SparseMatrixHandler::iGetNumCols() + dcol)
		|| (out.iGetNumRows() < SparseMatrixHandler::iGetNumRows() + drow)) {
		silent_cerr("Assertion fault "
				"in CColMatrixHandler<off>::MulAndSumWithShift"
				<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	drow = drow + 1;
	for (integer col = 0; col < SparseMatrixHandler::iGetNumCols(); col++) {
		integer idx = CompactSparseMatrixHandler_tpl<off>::Ap[col] - off;
		integer idxe = CompactSparseMatrixHandler_tpl<off>::Ap[col + 1] - off;
		integer newcol = col + dcol + 1;
		for (; idx < idxe; idx++) {
			out.IncCoef(CompactSparseMatrixHandler_tpl<off>::Ai[idx] - off + drow,
				newcol, CompactSparseMatrixHandler_tpl<off>::Ax[idx]*s);
		}
	}
	return out;	
}

template <int off>
MatrixHandler&
CColMatrixHandler<off>::FakeThirdOrderMulAndSumWithShift(MatrixHandler& out, 
		std::vector<bool> b,
		doublereal s,
		integer drow, 
		integer dcol) const
{
	silent_cerr("CColMatrixHandler<off>::FakeThirdOrderMulAndSumWithShift "
			"called" << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);		
	if ((out.iGetNumCols() < SparseMatrixHandler::iGetNumCols() + dcol)
			|| (out.iGetNumRows() < SparseMatrixHandler::iGetNumRows() + drow)) {
		silent_cerr("Assertion fault "
				"in CColMatrixHandler<off>::MulAndSumWithShift"
				<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	drow = drow + 1;
	for (integer col = 0; col < SparseMatrixHandler::iGetNumCols(); col++) {
		integer idx = CompactSparseMatrixHandler_tpl<off>::Ap[col];
		integer idxe = CompactSparseMatrixHandler_tpl<off>::Ap[col + 1] - off;
		integer newcol = col + dcol + 1;
		for (; idx < idxe; idx++) {
			if (b[CompactSparseMatrixHandler_tpl<off>::Ai[idx] - off]) {
				out.IncCoef(CompactSparseMatrixHandler_tpl<off>::Ai[idx] - off + drow,
					newcol, CompactSparseMatrixHandler_tpl<off>::Ax[idx]*s);
			}
		}
	}
	return out;	
}

template class CColMatrixHandler<0>;
template class CColMatrixHandler<1>;
