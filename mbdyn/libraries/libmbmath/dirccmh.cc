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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "spmapmh.h"
#include "dirccmh.h"

template <int off>
DirCColMatrixHandler<off>::DirCColMatrixHandler(std::vector<doublereal>& x,
		const std::vector<integer>& i,
		const std::vector<integer>& p)
: CompactSparseMatrixHandler(p.size() - 1, p.size() - 1, x, i, p),
  pindices(NCols + 1), indices(NRows*NCols, -1)
{
	for (integer col = 1; col <= NCols; col++) {
		pindices[col] = &indices[(col - 1)*NRows] - 1;

		integer row_begin = p[col - 1], row_end = p[col];

		for (integer r = row_begin; r < row_end; r++) {
			pindices[col][i[r] + 1 - off] = r;
		}
	}
}

template <int off>
DirCColMatrixHandler<off>::~DirCColMatrixHandler()
{
	NO_OP;
}

/* used by MultiThreadDataManager to duplicate the storage array
 * while preserving the CC indices */
template <int off>
CompactSparseMatrixHandler *
DirCColMatrixHandler<off>::Copy(void) const
{
	std::vector<doublereal> *pax = new std::vector<doublereal>(Ax);
	DirCColMatrixHandler<off> *p = new DirCColMatrixHandler<off>(*pax, Ai, Ap);
	p->bMatDuplicate = true;

	return p;
}

template <int off>
integer
DirCColMatrixHandler<off>::MakeCompressedColumnForm(doublereal *const Ax,
		integer *const Ai, integer *const Ap,
		int offset) const
{
	return Nz();
}

template <int off>
integer
DirCColMatrixHandler<off>::MakeCompressedColumnForm(std::vector<doublereal>& Ax,
		std::vector<integer>& Ai, std::vector<integer>& Ap,
		int offset) const
{
	return Nz();
}

template <int off>
integer
DirCColMatrixHandler<off>::MakeIndexForm(doublereal *const rAx,
		integer *const Arow, integer *const Acol,
		integer *const AcolSt, int offset) const
{
	silent_cerr("DirCColMatrixHandler<off>::MakeIndexForm called" << std::endl);
	THROW(ErrGeneric());

	/*
	 * copy Ax in rAx;
	 * copy Ai in Arow;
	 * build Acol from Ap;
	 * copy Ap in AcolSt;
	 */
#if 0
	for (integer col = 0; col < NCols; col++) {
		integer idx = Ap[col];
		integer idxe = Ap[col+1];
		for (; idx < idxe; idx++) {
			Arow[idx] = Ai[idx];
			Acol[idx] = col;
			rAx[idx] = Ax[idx];
		}
	}
#endif

	return Nz();
}

template <int off>
integer
DirCColMatrixHandler<off>::MakeIndexForm(std::vector<doublereal>& rAx,
                std::vector<integer>& Arow, std::vector<integer>& Acol,
		std::vector<integer>& AcolSt, int offset) const
{
	silent_cerr("DirCColMatrixHandler<off>::MakeIndexForm called" << std::endl);

	rAx.resize(Nz());
	Arow.resize(Nz());
	Acol.resize(Nz());
	AcolSt.resize(iGetNumCols());

	return MakeIndexForm(&(rAx[0]), &(Arow[0]), &(Acol[0]),
			&(AcolSt[0]), offset);
}

template <int off>
void
DirCColMatrixHandler<off>::Resize(integer n, integer nn)
{
	silent_cerr("DirCColMatrixHandler<off>::Resize called" << std::endl);
	THROW(ErrGeneric());
}

/* Estrae una colonna da una matrice */
template <int off>
VectorHandler&
DirCColMatrixHandler<off>::GetCol(integer icol, VectorHandler& out) const
{
	/*
	 * Note: we assume out has been reset
	 */
	
        if (icol > iGetNumCols()) {
		THROW(ErrGeneric());
	}

	integer idx = Ap[icol - 1];
	integer idxe = Ap[icol];

	for ( ; idx < idxe; idx++) {
		out(Ai[idx] - off + 1) = Ax[idx];
	}

	return out;
}
	
/* Prodotto Matrice per Matrice */
template <int off>
SpMapMatrixHandler& 
DirCColMatrixHandler<off>::MatMatMul(SpMapMatrixHandler& out,
		const SpMapMatrixHandler& in) const
{
	silent_cerr("DirCColMatrixHandler<off>::MatMatMul called" << std::endl);
	THROW(ErrGeneric());		
/*
 * 	if ((in.iGetNumCols() != iGetNumRows())
 * 			|| (in.iGetNumRows() != out.iGetNumRows())
 * 			|| (out.iGetNumCols() != iGetNumCols())) {
 * 		silent_cerr("Assertion fault in SpMapMatrixHandler::MatMatMul"
 * 			<< std::endl);
 * 		THROW(ErrGeneric());
 * 	}
 * 	out.Reset(0.);
 * 	for (integer col=0; col<NCols; col++) {
 * 		row_cont_type::const_iterator ri, re;
 * 		re = col_indices[col].end();
 * 		for (ri = col_indices[col].begin(); ri!=re; ri++) {
 * 			integer iend = in.iGetNumCols();
 * 			for (integer col2=0; col2<iend;  col2++) {
 * 				out.IncCoef(ri->first,col2,ri->second*in.dGetCoef(col,col2));
 * 			}
 * 		}
 * 	}
 */
	return out;	
}
	
/* Moltiplica per uno scalare e somma a una matrice */
template <int off>
MatrixHandler&
DirCColMatrixHandler<off>::MulAndSumWithShift(MatrixHandler& out, doublereal s,
		integer drow, integer dcol) const
{
	silent_cerr("DirCColMatrixHandler<off>::MulAndSumWithShift called"
			<< std::endl);
	THROW(ErrGeneric());		
	if ((out.iGetNumCols() < iGetNumCols()+dcol)
		|| (out.iGetNumRows() < iGetNumRows()+drow)) {
		silent_cerr("Assertion fault "
				"in DirCColMatrixHandler<off>::MulAndSumWithShift"
				<< std::endl);
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
}

template <int off>
MatrixHandler&
DirCColMatrixHandler<off>::FakeThirdOrderMulAndSumWithShift(MatrixHandler& out, 
		std::vector<bool> b,
		doublereal s,
		integer drow, 
		integer dcol) const
{
	silent_cerr("DirCColMatrixHandler<off>::FakeThirdOrderMulAndSumWithShift "
			"called" << std::endl);
	THROW(ErrGeneric());		
	if ((out.iGetNumCols() < iGetNumCols()+dcol)
			|| (out.iGetNumRows() < iGetNumRows()+drow)) {
		silent_cerr("Assertion fault "
				"in DirCColMatrixHandler<off>::MulAndSumWithShift"
				<< std::endl);
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
}
	
template <int off>
VectorHandler&
DirCColMatrixHandler<off>::MatTVecMul(VectorHandler& out,
		const VectorHandler& in) const
{
	silent_cerr("DirCColMatrixHandler<off>::MatTVecMul called" << std::endl);
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
}
	
template <int off>
VectorHandler&
DirCColMatrixHandler<off>::MatVecMul(VectorHandler& out,
		const VectorHandler& in) const
{
	silent_cerr("DirCColMatrixHandler<off>::MatTVecMul called" << std::endl);
	THROW(ErrGeneric());		
	if (in.iGetSize() != iGetNumCols() 
			|| out.iGetSize() != iGetNumRows()) {
		THROW(ErrGeneric());
  	}

	out.Reset(0.);
	return MatVecIncMul(out, in);
}

template <int off>
VectorHandler&
DirCColMatrixHandler<off>::MatTVecIncMul(VectorHandler& out,
		const VectorHandler& in) const
{
	silent_cerr("DirCColMatrixHandler<off>::MatTVecMul called" << std::endl);
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
}
	
template <int off>
VectorHandler&
DirCColMatrixHandler<off>::MatVecIncMul(VectorHandler& out,
		const VectorHandler& in) const {
	silent_cerr("DirCColMatrixHandler<off>::MatTVecMul called" << std::endl);
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
}

template <int off>
VectorHandler&
DirCColMatrixHandler<off>::MatVecDecMul(VectorHandler& out,
		const VectorHandler& in) const
{
	silent_cerr("DirCColMatrixHandler<off>::MatTVecMul called" << std::endl);
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
}

template class DirCColMatrixHandler<0>;
template class DirCColMatrixHandler<1>;
