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
class CColMatrixHandler : public CompactSparseMatrixHandler {
private:
#ifdef DEBUG
	void IsValid(void) const {
		NO_OP;
	};
#endif /* DEBUG */

public:
	CColMatrixHandler(std::vector<doublereal>& x,
			const std::vector<int>& i,
			const std::vector<int>& p);

	virtual ~CColMatrixHandler();

	/* used by MultiThreadDataManager to duplicate the storage array
	 * while preserving the CC indices */
	CompactSparseMatrixHandler *Copy(void) const;

public:
	doublereal & operator()(integer i_row, integer i_col) {
		ASSERTMSGBREAK(i_row > 0 && i_row <= NRows,
				"Error in CColMatrixHandler::operator(), "
				"row index out of range");
		ASSERTMSGBREAK(i_col > 0 && i_col <= NCols,
				"Error in CColMatrixHandler::operator(), "
				"col index out of range");
		i_row--;
		integer row_begin = Ap[i_col - 1];
		integer row_end = Ap[i_col] - 1;
		integer idx;
		integer row;

		if (row_begin == Ap[i_col]
				|| Ai[row_begin] - off > i_row
				|| Ai[row_end] - off < i_row) {
			/* matrix must be rebuilt */
			THROW(ErrRebuildMatrix());
		}

		while (row_end >= row_begin) {
			idx = (row_begin + row_end)/2;
			row = Ai[idx] - off;
			if (i_row < row) {
				row_end = idx - 1;
			} else if (i_row > row) {
				row_begin = idx + 1;
			} else {
				return Ax[idx];
			}
		}

		/* matrix must be rebuilt */
		THROW(ErrRebuildMatrix());
	};

	const doublereal& operator () (integer i_row, integer i_col) const {
		ASSERTMSGBREAK(i_row > 0 && i_row <= NRows,
				"Error in CColMatrixHandler::operator(), "
				"row index out of range");
		ASSERTMSGBREAK(i_col > 0 && i_col <= NCols,
				"Error in CColMatrixHandler::operator(), "
				"col index out of range");
		i_row--;
		integer row_begin = Ap[i_col - 1];
		integer row_end = Ap[i_col] - 1;
		integer idx;
		integer row;

		if (row_begin == Ap[i_col]
				|| Ai[row_begin] - off > i_row
				|| Ai[row_end] - off < i_row) {
			return ::dZero;
		}

		while (row_end >= row_begin) {
			idx = (row_begin + row_end)/2;
			row = Ai[idx] - off;
			if (i_row < row) {
				row_end = idx - 1;
			} else if (i_row > row) {
				row_begin = idx + 1;
			} else {
				return Ax[idx];
			}
		}

		return ::dZero;
	};

	int MakeCompressedColumnForm(doublereal *const Ax,
			int *const Ai, int *const Ap,
			int offset = 0) const;

        int MakeCompressedColumnForm(std::vector<doublereal>& Ax,
                	std::vector<int>& Ai, std::vector<int>& Ap,
			int offset = 0) const;

	int MakeIndexForm(doublereal *const rAx,
			integer *const Arow, integer *const Acol,
			integer *const AcolSt,
			int offset = 0) const;

        int MakeIndexForm(std::vector<doublereal>& rAx,
                	std::vector<integer>& Arow, std::vector<integer>& Acol,
			std::vector<integer>& AcolSt,
			int offset = 0) const;

	void Resize(integer ir, integer ic);

	/* Estrae una colonna da una matrice */
	VectorHandler& GetCol(integer icol, VectorHandler& out) const;
	
        /* Prodotto Matrice per Matrice */
	SpMapMatrixHandler& MatMatMul(SpMapMatrixHandler& out,
			const SpMapMatrixHandler& in) const;
	
        /* Moltiplica per uno scalare e somma a una matrice */
	MatrixHandler& MulAndSumWithShift(MatrixHandler& out,
			doublereal s = 1.,
			integer drow = 0, integer dcol = 0) const;
	
	MatrixHandler& FakeThirdOrderMulAndSumWithShift(MatrixHandler& out, 
		std::vector<bool> b, doublereal s = 1.,
		integer drow = 0, integer dcol = 0) const;
	
	VectorHandler& MatTVecMul(VectorHandler& out,
			const VectorHandler& in) const;
	
	VectorHandler& MatVecMul(VectorHandler& out,
			const VectorHandler& in) const;

	VectorHandler& MatTVecIncMul(VectorHandler& out,
			const VectorHandler& in) const;
	
	VectorHandler& MatVecIncMul(VectorHandler& out,
			const VectorHandler& in) const;

	VectorHandler& MatVecDecMul(VectorHandler& out,
			const VectorHandler& in) const;
};


#endif /* CColMatrixHandler_hh */

