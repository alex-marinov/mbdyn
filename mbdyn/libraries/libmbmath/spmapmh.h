/* 
 * HmFe (C) is a FEM analysis code. 
 *
 * Copyright (C) 1996-2001
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
/*
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
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
#include <myassert.h>
#include <solman.h>

class SpMapMatrixHandler : public MatrixHandler {
private:
//	typedef std::slist<double> x_cont_type;
	typedef std::map<int,double> row_cont_type;
	int N;
	int NZ;
	double zero;
	std::vector<row_cont_type> col_indices;
//	x_cont_type x;	
	
	void IsValid(void) const {
		NO_OP;
	};
public:
	SpMapMatrixHandler(const int &n = 0) : N(n), NZ(0), zero(0.) {
		col_indices.resize(n);
	};
	virtual ~SpMapMatrixHandler() {};
	void Init(const double& = 0.) {
		NO_OP;
	};
	long int iGetNumRows(void) const {
		return N;
	}
	long int iGetNumCols(void) const {
		return N;
	}
	
	double & operator()(const int &i_row, const int &i_col) {
		ASSERTMSGBREAK(i_row < N, "Error in SpMapMatrixHandler::operator()(const int&, const int&), row index out of range");
		ASSERTMSGBREAK(i_col < N, "Error in SpMapMatrixHandler::operator()(const int&, const int&), col index out of range");
		row_cont_type::iterator i;
		row_cont_type & row = col_indices[i_col];
		i = row.find(i_row);
		if (i == row.end()) {
			NZ++;
			return row[i_row] = 0.;
		} else {
			return i->second;
		}
	};
	long int fIncCoef(long int ix, long int iy, const double& inc) {
		ASSERTMSGBREAK(ix-1 < N, "Error in SpMapMatrixHandler::operator()(const int&, const int&), row index out of range");
		ASSERTMSGBREAK(iy-1 < N, "Error in SpMapMatrixHandler::operator()(const int&, const int&), col index out of range");
		//try to keep sparsity
		if (inc != 0.) {
			operator()(ix-1,iy-1) += inc;
		}
		return 1;
	};
	long int fDecCoef(long int ix, long int iy, const double& inc) {
		ASSERTMSGBREAK(ix-1 < N, "Error in SpMapMatrixHandler::operator()(const int&, const int&), row index out of range");
		ASSERTMSGBREAK(iy-1 < N, "Error in SpMapMatrixHandler::operator()(const int&, const int&), col index out of range");
		//try to keep sparsity
		if (inc != 0.) {
			operator()(ix-1,iy-1) -= inc;
		}
		return 1;
	};
	long int fPutCoef(long int ix, long int iy, const double& val) {
		ASSERTMSGBREAK(ix-1 < N, "Error in SpMapMatrixHandler::operator()(const int&, const int&), row index out of range");
		ASSERTMSGBREAK(iy-1 < N, "Error in SpMapMatrixHandler::operator()(const int&, const int&), col index out of range");
		//try to keep sparsity
		if (val != 0.) {
			operator()(ix-1,iy-1) = val;
		} else {
			row_cont_type::iterator i;
			row_cont_type & row = col_indices[iy-1];
			i = row.find(ix-1);
			if (i == row.end()) {
				//do nothing
			} else {
				i->second = val;
			}
		}
		return 1;
	};
	const double& dGetCoef(long int ix, long int iy) const {
		ASSERTMSGBREAK(ix-1 < N, "Error in SpMapMatrixHandler::operator()(const int&, const int&), row index out of range");
		ASSERTMSGBREAK(iy-1 < N, "Error in SpMapMatrixHandler::operator()(const int&, const int&), col index out of range");
		//try to keep sparsity
		row_cont_type::iterator i;
		row_cont_type & row = ((std::vector<row_cont_type>&)col_indices)[iy-1];
		i = row.find(ix-1);
		if (i == row.end()) {
			return zero;
		} else {
			return i->second;
		}
	};
	void MakeCompressedColumnForm(
		double *const Ax,
		int *const Ai,
		int *const Ap) {
		int x_ptr = 0;
		
		row_cont_type::const_iterator ri, re;
		
		for (int col=0; col<N; col++) {
			Ap[col] = x_ptr;
			re = col_indices[col].end();
			for (ri = col_indices[col].begin();ri != re; ri++) {
				Ax[x_ptr] = ri->second;
				Ai[x_ptr] = ri->first;
				x_ptr++;
			}
		}
		ASSERTMSGBREAK(x_ptr == NZ, "Error in SpMapMatrixHandler::MakeCompressedColumnForm");
		Ap[N] = x_ptr;
	};
	void Reset(Real r = 0.) {
		row_cont_type::const_iterator re;
		row_cont_type::iterator ri;
		for (int col=0; col<N; col++) {
			re = col_indices[col].end();
			for (ri = col_indices[col].begin();ri != re; ri++) {
				ri->second = r;
			}
		}
	};
/*
 * 	void Reset() {
 * 		for (int i=0; i<N; i++) {
 * 			col_indices[i].clear();
 * 		};
 * 		NZ = 0;
 * 	};
 */
	void Resize(const int &n) {
		Reset();
		col_indices.resize(n);
		N = n;
	};
	const int Nz() {
		return NZ;
	};
};

/*
 * #include "SubMatrix.hh"
 * void SubMatrix<SpMapMatrixHandler>::addtovalue(
 * 	SpMapMatrixHandler*const m,
 * 	const double &x,
 * 	int ix,
 * 	int iy) const {
 * 	if (x != 0) {
 * 		(*m)(ix,iy) += x;
 * 	}
 * };
 */


#endif //SpMapMatrixHandler_hh
