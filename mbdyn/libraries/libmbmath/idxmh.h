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

#ifndef IDXMH_H
#define IDXMH_H

#include <vector>

#include "myassert.h"
#include "solman.h"
#include "spmh.h"

/* Sparse Matrix in columns form */
class IndexMatrixHandler : public CompactSparseMatrixHandler {
private:
#ifdef DEBUG
	void IsValid(void) const {
		NO_OP;
	};
#endif /* DEBUG */

public:
	IndexMatrixHandler(const int &n,
			std::vector<doublereal>&x,
			const std::vector<int>& i,
			const std::vector<int>& p);

	virtual ~IndexMatrixHandler();

	/* used by MultiThreadDataManager to duplicate the storage array
	 * while preserving the CC indices */
	CompactSparseMatrixHandler *Copy(void) const;

public:
	doublereal & operator()(integer i_row, integer i_col) {
		ASSERTMSGBREAK(i_row > 0 && i_row <= NRows,
				"Error in IndexMatrixHandler::operator(), "
				"row index out of range");
		ASSERTMSGBREAK(i_col > 0 && i_col <= NCols,
				"Error in IndexMatrixHandler::operator(), "
				"col index out of range");

		/* matrix must be rebuilt */
		THROW(ErrRebuildMatrix());
	};

	const doublereal& operator () (integer i_row, integer i_col) const {
		ASSERTMSGBREAK(ix > 0 && ix <= NRows,
				"Error in IndexMatrixHandler::operator(), "
				"row index out of range");
		ASSERTMSGBREAK(iy > 0 && iy <= NCols,
				"Error in IndexMatrixHandler::operator(), "
				"col index out of range");

		return zero;
	};

	int MakeCompressedColumnForm(doublereal *const Ax,
			int *const Ai, int *const Ap,
			integer offset = 0) const;

        int MakeCompressedColumnForm(std::vector<doublereal>& Ax,
                	std::vector<int>& Ai, std::vector<int>& Ap,
			integer offset = 0) const;

	int MakeIndexForm(doublereal *const rAx, integer *const Arow,
			integer *const Acol,
			integer offset = 0) const;

        int MakeIndexForm(std::vector<doublereal>& rAx,
                	std::vector<integer>& Arow, std::vector<integer>& Acol,
			integer offset = 0) const;

	void Reset(const doublereal& c = 0.);

	void Resize(const int &n, const int &nn = 0);

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

#endif /* IDXMH_H */

