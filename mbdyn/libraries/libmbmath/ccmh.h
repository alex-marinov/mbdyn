/*
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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
#include <myassert.h>
#include <solman.h>

/* Sparse Matrix in columns form */
class CColMatrixHandler : public MatrixHandler {
private:
	std::vector<doublereal>& Ax;
	bool bMatDuplicate;
	const std::vector<int>& Ai;
	const std::vector<int>& Ap;
	int NRows;
	int NCols;
	int NZ;
	doublereal zero;
	
	void IsValid(void) const {
		NO_OP;
	};
public:
	CColMatrixHandler(std::vector<doublereal>&x,
		const std::vector<int>& i,
		const std::vector<int>& p) : Ax(x), bMatDuplicate(false), 
			Ai(i), Ap(p), 
			NRows(Ap.size()-1), NCols(NRows),
			NZ(Ax.size()), zero(0.) {
	};

	virtual ~CColMatrixHandler() {
		if (bMatDuplicate) {
			delete &Ax;
		}
	};

	/* used by MultiThreadDataManager to duplicate the storage array
	 * while preserving the CC indices */
	CColMatrixHandler * Copy(void) const {
		std::vector<doublereal> *pax = new std::vector<doublereal>(Ax);
		CColMatrixHandler *p = new CColMatrixHandler(*pax, Ai, Ap);
		p->bMatDuplicate = true;
		return p;
	};

	/* used to sum CC matrices with identical indices */
	void AddUnchecked(CColMatrixHandler& m) {
		/* FIXME: put in stl-ish form;
		 * see if we can use something from optimized blas,
		 * e.g. ATLAS, goto or so... */

		/* checks - uncomment to enable */
#if 0
		ASSERT(Ax.size() == m.Ax.size());
		ASSERT(Ai.size() == m.Ai.size());
		ASSERT(Ap.size() == m.Ap.size());
		for (unsigned i = 0; i < Ai.size(); i++) {
			if (Ai[i] != m.Ai[i]) {
				std::cerr << "Ai[" << i << "]" << std::endl;
				THROW(ErrGeneric());
			}
		}
		for (unsigned i = 0; i < Ap.size(); i++) {
			if (Ap[i] != m.Ap[i]) {
				std::cerr << "Ap[" << i << "]" << std::endl;
				THROW(ErrGeneric());
			}
		}
#endif
		
		doublereal *d = &Ax[0], *s = &m.Ax[0];
		for (unsigned long i = 0; i < Ax.size(); i++) {
			d[i] += s[i];
		}
	};

	/* Restituisce un puntatore all'array di reali della matrice */
	virtual inline doublereal* pdGetMat(void) const {
		return &Ax[0];
	};

	void Init(const doublereal& c = 0.) {
		Reset(c);
	};
	integer iGetNumRows(void) const {
		return NRows;
	}
	integer iGetNumCols(void) const {
		return NCols;
	}

public:
	doublereal & operator()(integer i_row, integer i_col) {
		ASSERTMSGBREAK(i_row > 0 && i_row <= NRows, "Error in CColMatrixHandler::operator(), row index out of range");
		ASSERTMSGBREAK(i_col > 0 && i_col <= NCols, "Error in CColMatrixHandler::operator(), col index out of range");
		i_row--;
		integer row_begin = Ap[i_col - 1];
		integer row_end = Ap[i_col] - 1;
		integer idx;
		integer row;

		if (row_begin == Ap[i_col] || Ai[row_begin] > i_row || Ai[row_end] < i_row) {
			THROW(ErrRebuildMatrix());
		}

		while (row_end >= row_begin) {
			idx = (row_begin + row_end)/2;
			row = Ai[idx];
			if (i_row < row) {
				row_end = idx - 1;
			} else if (i_row > row) {
				row_begin = idx + 1;
			} else {
				return Ax[idx];
			}
		}
		/* fixme: handle case */
		THROW(ErrRebuildMatrix());
	};
	const doublereal& operator () (integer i_row, integer i_col) const {
		ASSERTMSGBREAK(ix > 0 && ix <= NRows, "Error in CColMatrixHandler::operator(), row index out of range");
		ASSERTMSGBREAK(iy > 0 && iy <= NCols, "Error in CColMatrixHandler::operator(), col index out of range");
		i_row--;
		integer row_begin = Ap[i_col - 1];
		integer row_end = Ap[i_col] - 1;
		integer idx;
		integer row;

		if (row_begin == Ap[i_col] || Ai[row_begin] > i_row || Ai[row_end] < i_row) {
			return zero;
		}

		while (row_end >= row_begin) {
			idx = (row_begin + row_end)/2;
			row = Ai[idx];
			if (i_row < row) {
				row_end = idx - 1;
			} else if (i_row > row) {
				row_begin = idx + 1;
			} else {
				return Ax[idx];
			}
		}
		/* fixme: handle case */
		return zero;
	};
	void IncCoef(integer ix, integer iy, const doublereal& inc) {
		ASSERTMSGBREAK(ix > 0 && ix <= NRows, "Error in CColMatrixHandler::IncCoef(), row index out of range");
		ASSERTMSGBREAK(iy > 0 && iy <= NCols, "Error in CColMatrixHandler::IncCoef(), col index out of range");
		operator()(ix,iy) += inc;
	};
	void DecCoef(integer ix, integer iy, const doublereal& inc) {
		ASSERTMSGBREAK(ix > 0 && ix <= NRows, "Error in CColMatrixHandler::DecCoef(), row index out of range");
		ASSERTMSGBREAK(iy > 0 && iy <= NCols, "Error in CColMatrixHandler::DecCoef(), col index out of range");
		operator()(ix,iy) -= inc;
	};
	void PutCoef(integer ix, integer iy, const doublereal& val) {
		ASSERTMSGBREAK(ix-1 < NRows, "Error in CColMatrixHandler::PutCoef(), row index out of range");
		ASSERTMSGBREAK(iy-1 < NCols, "Error in CColMatrixHandler::PutCoef(), col index out of range");
		operator()(ix,iy) = val;
	};
	const doublereal& dGetCoef(integer ix, integer iy) const {
		ASSERTMSGBREAK(ix > 0 && ix <= NRows, "Error in CColMatrixHandler::dGetCoef(), row index out of range");
		ASSERTMSGBREAK(iy > 0 && iy <= NCols, "Error in CColMatrixHandler::dGetCoef(), col index out of range");
		return operator()(ix, iy);
	};
	void MakeCompressedColumnForm(
		doublereal *const Ax,
		int *const Ai,
		int *const Ap) const {
		std::cerr << "CColMatrixHandler::MakeCompressedColumnForm called" << std::endl;
		THROW(ErrGeneric());		
	};
        void MakeCompressedColumnForm(
                std::vector<doublereal>& Ax,
                std::vector<int>& Ai,
                std::vector<int>& Ap) const {
		std::cerr << "CColMatrixHandler::MakeCompressedColumnForm called" << std::endl;
		THROW(ErrGeneric());		
        };
	integer MakeIndexForm(
		doublereal *const rAx,
		integer *const Arow,
		integer *const Acol,
		integer offset=0) const {
		std::cerr << "CColMatrixHandler::MakeIndexForm called" << std::endl;
		THROW(ErrGeneric());	
		for (integer col = 0; col < NCols; col++) {
			integer idx = Ap[col];
			integer idxe = Ap[col+1];
			for (; idx < idxe; idx++) {
				Arow[idx] = Ai[idx];
				Acol[idx] = col;
				rAx[idx] = Ax[idx];
			}
		}
		return Nz();
	};
        integer MakeIndexForm(
                std::vector<doublereal>& rAx,
                std::vector<integer>& Arow,
                std::vector<integer>& Acol,
		integer offset=0) const {
		std::cerr << "CColMatrixHandler::MakeIndexForm called" << std::endl;
		THROW(ErrGeneric());		
		rAx.resize(Nz());
		Arow.resize(Nz());
		Acol.resize(Nz());
                return MakeIndexForm(&(rAx[0]),&(Arow[0]),&(Acol[0]),offset);
        };
	void Reset(doublereal r = 0.) {
		std::fill(Ax.begin(), Ax.end(), r);
	};
	void Resize(const int &n, const int &nn = 0) {
		std::cerr << "CColMatrixHandler::Resize called" 
			<< std::endl;
		THROW(ErrGeneric());
	};
	const int Nz() const {
		return NZ;
	};
	
	/* Estrae una colonna da una matrice */
	VectorHandler& GetCol(integer icol, VectorHandler& out) const {
		std::cerr << "CColMatrixHandler::GetCol called" << std::endl;
		THROW(ErrGeneric());		
	        if (icol > iGetNumCols()) {
			THROW(ErrGeneric());
		}
		integer idx = Ap[icol-1];
		integer idxe = Ap[icol];
		for (; idx<idxe; idx++) {
			out.PutCoef(Ai[idx]+1, Ax[idx]);
		}
		return out;
	};
	
        /* Prodotto Matrice per Matrice */
	SpMapMatrixHandler& MatMatMul(SpMapMatrixHandler& out, const SpMapMatrixHandler& in) const {
		std::cerr << "CColMatrixHandler::MatMatMul called" << std::endl;
		THROW(ErrGeneric());		
/*
 * 		if ((in.iGetNumCols() != iGetNumRows())
 * 				|| (in.iGetNumRows() != out.iGetNumRows())
 * 				|| (out.iGetNumCols() != iGetNumCols())) {
 * 			std::cerr << "Assertion fault in SpMapMatrixHandler::MatMatMul" << std::endl;
 * 			THROW(ErrGeneric());
 * 		}
 * 		out.Reset(0.);
 * 		for (int col=0; col<NCols; col++) {
 * 			row_cont_type::const_iterator ri, re;
 * 			re = col_indices[col].end();
 * 			for (ri = col_indices[col].begin(); ri!=re; ri++) {
 * 				int iend = in.iGetNumCols();
 * 				for (int col2=0; col2<iend;  col2++) {
 * 				out.IncCoef(ri->first,col2,ri->second*in.dGetCoef(col,col2));
 * 				}
 * 			}
 * 		}
 */
		return out;	
	};
	
        /* Moltiplica per uno scalare e somma a una matrice */
	MatrixHandler& MulAndSumWithShift(MatrixHandler& out, doublereal s = 1.,
		integer drow = 0, integer dcol = 0) const {
		std::cerr << "CColMatrixHandler::MulAndSumWithShift called" << std::endl;
		THROW(ErrGeneric());		
		if ((out.iGetNumCols() < iGetNumCols()+dcol)
			|| (out.iGetNumRows() < iGetNumRows()+drow)) {
			std::cerr << "Assertion fault in CColMatrixHandler::MulAndSumWithShift" << std::endl;
			THROW(ErrGeneric());
		}
		drow = drow + 1;
		for (integer col = 0; col < NCols; col++) {
			integer idx = Ap[col];
			integer idxe = Ap[col+1];
			integer newcol = col + dcol + 1;
			for (; idx < idxe; idx++) {
				out.IncCoef(Ai[idx]+drow,newcol,Ax[idx]*s);
			}
		}
		return out;	
	};
	
	MatrixHandler& FakeThirdOrderMulAndSumWithShift(
		MatrixHandler& out, 
		std::vector<bool> b,
		doublereal s = 1.,
		integer drow = 0, 
		integer dcol = 0) const {
		std::cerr << "CColMatrixHandler::FakeThirdOrderMulAndSumWithShift called" << std::endl;
		THROW(ErrGeneric());		
		if ((out.iGetNumCols() < iGetNumCols()+dcol)
			|| (out.iGetNumRows() < iGetNumRows()+drow)) {
			std::cerr << "Assertion fault in CColMatrixHandler::MulAndSumWithShift" << std::endl;
			THROW(ErrGeneric());
		}
		drow = drow + 1;
		for (integer col = 0; col < NCols; col++) {
			integer idx = Ap[col];
			integer idxe = Ap[col+1];
			integer newcol = col + dcol + 1;
			for (; idx < idxe; idx++) {
				if (b[Ai[idx]]) {
					out.IncCoef(Ai[idx]+drow,newcol,Ax[idx]*s);
				}
			}
		}
		return out;	
	};
	
	VectorHandler& MatTVecMul(VectorHandler& out, const VectorHandler& in) const {
		std::cerr << "CColMatrixHandler::MatTVecMul called" << std::endl;
		THROW(ErrGeneric());		
		if (out.iGetSize() != iGetNumRows()
				|| in.iGetSize() != iGetNumCols()) {
			THROW(ErrGeneric());
		}

		for (integer col = 0; col < NCols; col++) {
			doublereal d = 0.;
			integer idx = Ap[col];
			integer idxe = Ap[col+1];
			for (; idx < idxe; idx++) {
				d += Ax[idx]*in.dGetCoef(Ai[idx]+1);
			}
			out.PutCoef(col+1, d);
		}
		return out;
	};
	
	VectorHandler& MatVecMul(VectorHandler& out, const VectorHandler& in) const {
		std::cerr << "CColMatrixHandler::MatTVecMul called" << std::endl;
		THROW(ErrGeneric());		
		if (in.iGetSize() != iGetNumCols() 
				|| out.iGetSize() != iGetNumRows()) {
			THROW(ErrGeneric());
  		}

		out.Reset(0.);
		return MatVecIncMul(out, in);
	};

	VectorHandler& MatTVecIncMul(VectorHandler& out, const VectorHandler& in) const {
		std::cerr << "CColMatrixHandler::MatTVecMul called" << std::endl;
		THROW(ErrGeneric());		
		if (out.iGetSize() != iGetNumRows()
				|| in.iGetSize() != iGetNumCols()) {
			THROW(ErrGeneric());
		}

		for (integer col = 0; col < NCols; col++) {
			doublereal d = 0.;
			integer idx = Ap[col];
			integer idxe = Ap[col+1];
			for (; idx < idxe; idx++) {
				d += Ax[idx]*in.dGetCoef(Ai[idx]+1);
			}
			out.IncCoef(col+1, d);
		}
		return out;
	};
	
	VectorHandler& MatVecIncMul(VectorHandler& out, const VectorHandler& in) const {
		std::cerr << "CColMatrixHandler::MatTVecMul called" << std::endl;
		THROW(ErrGeneric());		
		if (in.iGetSize() != iGetNumCols()
				|| out.iGetSize() != iGetNumRows()) {
			THROW(ErrGeneric());
		}

		for (integer col = 0; col < NCols; col++) {
			integer idx = Ap[col];
			integer idxe = Ap[col+1];
			for (; idx < idxe; idx++) {
				doublereal d = Ax[idx]*in.dGetCoef(Ai[idx]+1);
				out.IncCoef(Ai[idx]+1, d);
			}
		}
		return out;
	};

	VectorHandler& MatVecDecMul(VectorHandler& out, const VectorHandler& in) const {
		std::cerr << "CColMatrixHandler::MatTVecMul called" << std::endl;
		THROW(ErrGeneric());		
		if (in.iGetSize() != iGetNumCols()
				|| out.iGetSize() != iGetNumRows()) {
			THROW(ErrGeneric());
		}

		for (integer col = 0; col < NCols; col++) {
			integer idx = Ap[col];
			integer idxe = Ap[col+1];
			for (; idx < idxe; idx++) {
				doublereal d = Ax[idx]*in.dGetCoef(Ai[idx]+1);
				out.DecCoef(Ai[idx]+1, d);
			}
		}
		return out;
	};
};

/*
 * #include "SubMatrix.hh"
 * void SubMatrix<SpMapMatrixHandler>::addtovalue(
 * 	SpMapMatrixHandler*const m,
 * 	const doublereal &x,
 * 	int ix,
 * 	int iy) const {
 * 	if (x != 0) {
 * 		(*m)(ix,iy) += x;
 * 	}
 * };
 */


#endif //CColMatrixHandler_hh
