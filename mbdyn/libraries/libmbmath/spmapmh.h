/* 
 * HmFe (C) is a FEM analysis code. 
 *
 * Copyright (C) 1996-2004
 *
 * Marco Morandini  <morandini@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */
/* November 2001 
 * Modified to add methods
 * to be used in the parallel MBDyn Solver.
 *
 * Copyright (C) 1996-2004
 *
 * Giuseppe Quaranta  <quaranta@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *      
 */
/*
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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

#ifndef SpMapMatrixHandler_hh
#define SpMapMatrixHandler_hh

#include <map>
#include <vector>
#include "myassert.h"
#include "solman.h"
#include "spmh.h"

/* #define MBDYN_X_KEEP_SPARSITY */

/* Sparse Matrix in columns form */
class SpMapMatrixHandler : public SparseMatrixHandler {
private:
	struct payload { doublereal d; int p; };
	typedef std::map<int, payload> row_cont_type;
	mutable std::vector<row_cont_type> col_indices;

#ifdef DEBUG
	void IsValid(void) const {
		NO_OP;
	};
#endif /* DEBUG */

public:
	SpMapMatrixHandler(const int &n = 0,const int &nn = 0);

	virtual ~SpMapMatrixHandler();

	doublereal & operator()(integer i_row, integer i_col) {
		ASSERTMSGBREAK(i_row > 0 && i_row <= NRows,
				"Error in SpMapMatrixHandler::operator(), "
				"row index out of range");
		ASSERTMSGBREAK(i_col > 0 && i_col <= NCols, 
				"Error in SpMapMatrixHandler::operator(), "
				"col index out of range");
		i_row--;
		i_col--;
		row_cont_type::iterator i;
		row_cont_type& row = col_indices[i_col];
		i = row.find(i_row);
		if (i == row.end()) {
			NZ++;

			row[i_row].d = 0.;
			row[i_row].p = -1;

			return row[i_row].d;

		} else {
			return i->second.d;
		}
	};

	void IncCoef(integer ix, integer iy, const doublereal& inc) {
		ASSERTMSGBREAK(ix > 0 && ix <= NRows,
				"Error in SpMapMatrixHandler::IncCoef(), "
				"row index out of range");
		ASSERTMSGBREAK(iy > 0 && iy <= NCols,
				"Error in SpMapMatrixHandler::IncCoef(), "
				"col index out of range");
#ifdef MBDYN_X_KEEP_SPARSITY
		/* try to keep sparsity */
		if (inc != 0.) {
#endif /* MBDYN_X_KEEP_SPARSITY */
			operator()(ix,iy) += inc;
#ifdef MBDYN_X_KEEP_SPARSITY
		}
#endif /* MBDYN_X_KEEP_SPARSITY */
	};

	void DecCoef(integer ix, integer iy, const doublereal& inc) {
		ASSERTMSGBREAK(ix > 0 && ix <= NRows,
				"Error in SpMapMatrixHandler::DecCoef(), "
				"row index out of range");
		ASSERTMSGBREAK(iy > 0 && iy <= NCols,
				"Error in SpMapMatrixHandler::DecCoef(), "
				"col index out of range");
#ifdef MBDYN_X_KEEP_SPARSITY
		/* try to keep sparsity */
		if (inc != 0.) {
#endif /* MBDYN_X_KEEP_SPARSITY */
			operator()(ix,iy) -= inc;
#ifdef MBDYN_X_KEEP_SPARSITY
		}
#endif /* MBDYN_X_KEEP_SPARSITY */
	};

	void PutCoef(integer ix, integer iy, const doublereal& val) {
		ASSERTMSGBREAK(ix - 1 < NRows,
				"Error in SpMapMatrixHandler::PutCoef(), "
				"row index out of range");
		ASSERTMSGBREAK(iy - 1 < NCols,
				"Error in SpMapMatrixHandler::PutCoef(), "
				"col index out of range");
#ifdef MBDYN_X_KEEP_SPARSITY
		/* try to keep sparsity */
		if (val != 0.) {
#endif /* MBDYN_X_KEEP_SPARSITY */
			operator()(ix,iy) = val;
#ifdef MBDYN_X_KEEP_SPARSITY
		} else {
			row_cont_type::iterator i;
			row_cont_type& row = col_indices[iy-1];
			i = row.find(ix - 1);
			if (i != row.end()) {
				i->second.d = val;
			}
		}
#endif /* MBDYN_X_KEEP_SPARSITY */
	};

	const doublereal& dGetCoef(integer ix, integer iy) const {
		ASSERTMSGBREAK(ix > 0 && ix <= NRows,
				"Error in SpMapMatrixHandler::dGetCoef(), "
				"row index out of range");
		ASSERTMSGBREAK(iy > 0 && iy <= NCols,
				"Error in SpMapMatrixHandler::dGetCoef(), "
				"col index out of range");
		row_cont_type::iterator i;
		row_cont_type& row = col_indices[iy - 1];
		i = row.find(ix - 1);
		if (i == row.end()) {
			return zero;

		} else {
			return i->second.d;
		}
	};

	const doublereal& operator () (integer ix, integer iy) const {
		ASSERTMSGBREAK(ix > 0 && ix <= NRows,
				"Error in SpMapMatrixHandler::operator(), "
				"row index out of range");
		ASSERTMSGBREAK(iy > 0 && iy <= NCols,
				"Error in SpMapMatrixHandler::operator(), "
				"col index out of range");
		row_cont_type::iterator i;
		row_cont_type& row = col_indices[iy - 1];
		i = row.find(ix - 1);
		if (i == row.end()) {
			return zero;

		} else {
			return i->second.d;
		}
	};

	int MakeCompressedColumnForm(doublereal *const Ax,
			int *const Ai, int *const Ap,
			integer offset = 0) const;

        int MakeCompressedColumnForm(std::vector<doublereal>& Ax,
                	std::vector<int>& Ai, std::vector<int>& Ap,
			integer offset = 0) const;

	int MakeIndexForm(doublereal *const Ax,
			integer *const Arow, integer *const Acol,
			integer offset = 0) const;

        int MakeIndexForm(std::vector<doublereal>& Ax,
			std::vector<integer>& Arow, std::vector<integer>& Acol,
			integer offset = 0) const;

	void Reset(const doublereal &r = 0.);

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

#endif /* SpMapMatrixHandler_hh */

