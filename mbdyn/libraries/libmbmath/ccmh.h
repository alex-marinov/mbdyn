/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2003-2012
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

#ifndef CColMatrixHandler_hh
#define CColMatrixHandler_hh

#include <vector>

#include "myassert.h"
#include "solman.h"
#include "spmh.h"

/* Sparse Matrix in columns form */
template <int off>
class CColMatrixHandler : public CompactSparseMatrixHandler_tpl<off> {
public:
	typedef CompactSparseMatrixHandler_const_iterator<off> const_iterator;

private:
#ifdef DEBUG
	void IsValid(void) const {
		NO_OP;
	};
#endif /* DEBUG */

public:
	CColMatrixHandler(std::vector<doublereal>& x,
			const std::vector<integer>& i,
			const std::vector<integer>& p);

	virtual ~CColMatrixHandler();

	/* used by MultiThreadDataManager to duplicate the storage array
	 * while preserving the CC indices */
	CompactSparseMatrixHandler *Copy(void) const;

public:
	doublereal & operator()(integer i_row, integer i_col) {
		ASSERTMSGBREAK(i_row > 0 && i_row <= SparseMatrixHandler::iGetNumRows(),
				"Error in CColMatrixHandler::operator(), "
				"row index out of range");
		ASSERTMSGBREAK(i_col > 0 && i_col <= SparseMatrixHandler::iGetNumCols(),
				"Error in CColMatrixHandler::operator(), "
				"col index out of range");
		i_row--;
		integer row_begin = CompactSparseMatrixHandler_tpl<off>::Ap[i_col - 1] - off;
		integer row_end = CompactSparseMatrixHandler_tpl<off>::Ap[i_col] - off - 1;
		integer idx;
		integer row;

		if (row_begin == (CompactSparseMatrixHandler_tpl<off>::Ap[i_col] - off)
				|| (CompactSparseMatrixHandler_tpl<off>::Ai[row_begin] - off) > i_row
				|| (CompactSparseMatrixHandler_tpl<off>::Ai[row_end] - off) < i_row)
		{
			/* matrix must be rebuilt */
			throw MatrixHandler::ErrRebuildMatrix(MBDYN_EXCEPT_ARGS);
		}

		while (row_end >= row_begin) {
			idx = (row_begin + row_end)/2;
			row = CompactSparseMatrixHandler_tpl<off>::Ai[idx] - off;
			if (i_row < row) {
				row_end = idx - 1;
			} else if (i_row > row) {
				row_begin = idx + 1;
			} else {
				return CompactSparseMatrixHandler_tpl<off>::Ax[idx];
			}
		}

		/* matrix must be rebuilt */
		throw MatrixHandler::ErrRebuildMatrix(MBDYN_EXCEPT_ARGS);
	};

	const doublereal& operator () (integer i_row, integer i_col) const {
		ASSERTMSGBREAK(i_row > 0 && i_row <= SparseMatrixHandler::iGetNumRows(),
				"Error in CColMatrixHandler::operator(), "
				"row index out of range");
		ASSERTMSGBREAK(i_col > 0 && i_col <= SparseMatrixHandler::iGetNumCols(),
				"Error in CColMatrixHandler::operator(), "
				"col index out of range");
		i_row--;
		integer row_begin = CompactSparseMatrixHandler_tpl<off>::Ap[i_col - 1] - off;
		integer row_end = CompactSparseMatrixHandler_tpl<off>::Ap[i_col] - off - 1;
		integer idx;
		integer row;

		if (row_begin == CompactSparseMatrixHandler_tpl<off>::Ap[i_col] - off
				|| CompactSparseMatrixHandler_tpl<off>::Ai[row_begin] - off > i_row
				|| CompactSparseMatrixHandler_tpl<off>::Ai[row_end] - off < i_row)
		{
			return ::Zero1;
		}

		while (row_end >= row_begin) {
			idx = (row_begin + row_end)/2;
			row = CompactSparseMatrixHandler_tpl<off>::Ai[idx] - off;
			if (i_row < row) {
				row_end = idx - 1;
			} else if (i_row > row) {
				row_begin = idx + 1;
			} else {
				return CompactSparseMatrixHandler_tpl<off>::Ax[idx];
			}
		}

		return ::Zero1;
	};

	void Resize(integer ir, integer ic);

	/* Estrae una colonna da una matrice */
	VectorHandler& GetCol(integer icol, VectorHandler& out) const;

        /* Moltiplica per uno scalare e somma a una matrice */
	MatrixHandler& MulAndSumWithShift(MatrixHandler& out,
			doublereal s = 1.,
			integer drow = 0, integer dcol = 0) const;
	
	MatrixHandler& FakeThirdOrderMulAndSumWithShift(MatrixHandler& out, 
		std::vector<bool> b, doublereal s = 1.,
		integer drow = 0, integer dcol = 0) const;
};

#endif /* CColMatrixHandler_hh */

