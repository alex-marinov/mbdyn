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
/* November 2001 
 * Modified to add methods
 * to be used in the parallel MBDyn Solver.
 *
 * Copyright (C) 1996-2001
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
#include <myassert.h>
#include <solman.h>

/* Sparse Matrix in columns form */
class SpMapMatrixHandler : public MatrixHandler {
private:
//	typedef std::slist<doublereal> x_cont_type;
	typedef std::map<int,doublereal> row_cont_type;
	int NRows;
	int NCols;
	int NZ;
	double zero;
	std::vector<row_cont_type> col_indices;
//	x_cont_type x;	
	
	void IsValid(void) const {
		NO_OP;
	};
public:
	SpMapMatrixHandler(const int &n = 0,const int &nn = 0) : NZ(0), zero(0.) {
		int nnn;
		if (nn == 0) {
			nnn = n;
		} else {
			nnn = nn;
		}
		NRows = n; 
		NCols = nnn;
		col_indices.resize(NCols);
	};

	virtual ~SpMapMatrixHandler() {};

	void Init(const double& c = 0.) {
		Reset(c);
	};
	long int iGetNumRows(void) const {
		return NRows;
	}
	long int iGetNumCols(void) const {
		return NCols;
	}

private:
	double & operator()(const int &i_row, const int &i_col) {
		ASSERTMSGBREAK(i_row < NRows, "Error in SpMapMatrixHandler::operator()(const int&, const int&), row index out of range");
		ASSERTMSGBREAK(i_col < NCols, "Error in SpMapMatrixHandler::operator()(const int&, const int&), col index out of range");
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
public:
	long int fIncCoef(long int ix, long int iy, const double& inc) {
		ASSERTMSGBREAK(ix-1 < NRows, "Error in SpMapMatrixHandler::operator()(const int&, const int&), row index out of range");
		ASSERTMSGBREAK(iy-1 < NCols, "Error in SpMapMatrixHandler::operator()(const int&, const int&), col index out of range");
		//try to keep sparsity
		if (inc != 0.) {
			operator()(ix-1,iy-1) += inc;
		}
		return 1;
	};
	long int fDecCoef(long int ix, long int iy, const double& inc) {
		ASSERTMSGBREAK(ix-1 < NRows, "Error in SpMapMatrixHandler::operator()(const int&, const int&), row index out of range");
		ASSERTMSGBREAK(iy-1 < NCols, "Error in SpMapMatrixHandler::operator()(const int&, const int&), col index out of range");
		//try to keep sparsity
		if (inc != 0.) {
			operator()(ix-1,iy-1) -= inc;
		}
		return 1;
	};
	long int fPutCoef(long int ix, long int iy, const double& val) {
		ASSERTMSGBREAK(ix-1 < NRows, "Error in SpMapMatrixHandler::operator()(const int&, const int&), row index out of range");
		ASSERTMSGBREAK(iy-1 < NCols, "Error in SpMapMatrixHandler::operator()(const int&, const int&), col index out of range");
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
		ASSERTMSGBREAK(ix-1 < NRows, "Error in SpMapMatrixHandler::operator()(const int&, const int&), row index out of range");
		ASSERTMSGBREAK(iy-1 < NCols, "Error in SpMapMatrixHandler::operator()(const int&, const int&), col index out of range");
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
		int *const Ap) const {
		int x_ptr = 0;
		
		row_cont_type::const_iterator ri, re;
		
		for (int col=0; col<NCols; col++) {
			Ap[col] = x_ptr;
			re = col_indices[col].end();
			for (ri = col_indices[col].begin();ri != re; ri++) {
				Ax[x_ptr] = ri->second;
				Ai[x_ptr] = ri->first;
				x_ptr++;
			}
		}
		ASSERTMSGBREAK(x_ptr == NZ, "Error in SpMapMatrixHandler::MakeCompressedColumnForm");
		Ap[NCols] = x_ptr;
	};
        void MakeCompressedColumnForm(
                std::vector<double>& Ax,
                std::vector<int>& Ai,
                std::vector<int>& Ap) const {
                Ax.resize(Nz());
                Ai.resize(Nz());
                Ap.resize(iGetNumCols()+1);
                MakeCompressedColumnForm(&(Ax[0]),&(Ai[0]),&(Ap[0]));
        };
	integer MakeIndexForm(
		doublereal *const Ax,
		integer *const Arow,
		integer *const Acol,
		integer offset=0) const {
		
		integer x_ptr = 0;

		row_cont_type::const_iterator ri, re;		
		
		for (int col=0; col<NCols; col++) {
			re = col_indices[col].end();
			for (ri = col_indices[col].begin();ri != re; ri++) {
				Ax[x_ptr] = ri->second;
				Arow[x_ptr] = ri->first+offset;
				Acol[x_ptr] = col+offset;
				x_ptr++;
			}
		}
		ASSERTMSGBREAK(x_ptr == NZ, "Error in SpMapMatrixHandler::MakeIndexForm");
		return Nz();
	};
        integer MakeIndexForm(
                std::vector<doublereal>& Ax,
                std::vector<integer>& Arow,
                std::vector<integer>& Acol,
		integer offset=0) const {
                Ax.resize(Nz());
                Arow.resize(Nz());
                Acol.resize(Nz());
                return MakeIndexForm(&(Ax[0]),&(Arow[0]),&(Acol[0]),offset);
        };
	void Reset(Real r = 0.) {
		row_cont_type::const_iterator re;
		row_cont_type::iterator ri;
		for (int col=0; col<NCols; col++) {
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
	void Resize(const int &n, const int &nn = 0) {
		int nnn;
		if (nn == 0) {
			nnn = n;
		} else {
			nnn = nn;
		}
		col_indices.resize(nn);
		for (int col=0; col<NCols; col++) {
			col_indices[col].clear();
		}
		NCols = n;
		NRows = nnn;
		NZ = 0;
	};
	const int Nz() const {
		return NZ;
	};
	
	/* Estrae una colonna da una matrice */
	VectorHandler& GetCol(long int icol, VectorHandler& out) const {
	        if (icol > iGetNumCols()) {
			THROW(ErrGeneric());
		}
		row_cont_type::const_iterator ri, re;
		re = col_indices[icol].end();
		for (ri = col_indices[icol].begin();ri != re; ri++) {		
			out.fPutCoef(ri->first+1, ri->second);
		}
		return out;
	};		 	
	
        /* Prodotto Matrice per Matrice */
	SpMapMatrixHandler& MatMatMul(SpMapMatrixHandler& out, const SpMapMatrixHandler& in) const {
		if ((in.iGetNumCols() != iGetNumRows())
				|| (in.iGetNumRows() != out.iGetNumRows())
				|| (out.iGetNumCols() != iGetNumCols())) {
			THROW(ErrGeneric());
		}
		out.Reset(0.);
		for (int col=0; col<NCols; col++) {
			row_cont_type::const_iterator ri, re;
			re = col_indices[col].end();
			for (ri = col_indices[col].begin(); ri!=re; ri++) {
				int iend = in.iGetNumCols();
				for (int col2=0; col2<iend;  col2++) {
				out.fIncCoef(ri->first,col2,ri->second*in.dGetCoef(col,col2));
				}
			}
		}
		return out;	
	};
	
	VectorHandler& MatTVecMul(VectorHandler& out, const VectorHandler& in) const {
		if (out.iGetSize() != iGetNumRows()
				|| in.iGetSize() != iGetNumCols()) {
			THROW(ErrGeneric());
		}

		row_cont_type::const_iterator ri, re;
		for (int col=0; col<NCols; col++) {
			double d = 0.;
			re = col_indices[col].end();
			for (ri = col_indices[col].begin();ri != re; ri++) {
				d += ri->second*in.dGetCoef(ri->first+1);
			}
			out.fPutCoef(col+1, d);
		}
		return out;
	};
	
	VectorHandler& MatVecMul(VectorHandler& out, const VectorHandler& in) const {
		if (in.iGetSize() != iGetNumCols()
				|| out.iGetSize() != iGetNumRows()) {
			THROW(ErrGeneric());
  		}

		row_cont_type::const_iterator ri, re;
		out.Reset(0.);
		for (int col=0; col<NCols; col++) {
			re = col_indices[col].end();
			for (ri = col_indices[col].begin();ri != re; ri++) {
				double d = ri->second*in.dGetCoef(col+1);
				out.fIncCoef(ri->first+1, d);
			}
		}
		return out;
	};

	VectorHandler& MatTVecIncMul(VectorHandler& out, const VectorHandler& in) const {
		if (out.iGetSize() != iGetNumRows()
				|| in.iGetSize() != iGetNumCols()) {
			THROW(ErrGeneric());
		}

		row_cont_type::const_iterator ri, re;
		for (int col=0; col<NCols; col++) {
			double d = 0.;
			re = col_indices[col].end();
			for (ri = col_indices[col].begin();ri != re; ri++) {
				d += ri->second*in.dGetCoef(ri->first+1);
			}
			out.fIncCoef(col+1, d);
		}
		return out;
	};
	
	VectorHandler& MatVecIncMul(VectorHandler& out, const VectorHandler& in) const {
		if (in.iGetSize() != iGetNumCols()
				|| out.iGetSize() != iGetNumRows()) {
			THROW(ErrGeneric());
		}

		row_cont_type::const_iterator ri, re;
		for (int col=0; col<NCols; col++) {
			re = col_indices[col].end();
			for (ri = col_indices[col].begin();ri != re; ri++) {
				double d = ri->second*in.dGetCoef(col+1);
				out.fIncCoef(ri->first+1, d);
			}
		}
		return out;
	};

	VectorHandler& MatVecDecMul(VectorHandler& out, const VectorHandler& in) const {
		if (in.iGetSize() != iGetNumCols()
				|| out.iGetSize() != iGetNumRows()) {
			THROW(ErrGeneric());
		}

		row_cont_type::const_iterator ri, re;
		for (int col=0; col<NCols; col++) {
			re = col_indices[col].end();
			for (ri = col_indices[col].begin();ri != re; ri++) {
				double d = ri->second*in.dGetCoef(col+1);
				out.fDecCoef(ri->first+1, d);
			}
		}
		return out;
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
