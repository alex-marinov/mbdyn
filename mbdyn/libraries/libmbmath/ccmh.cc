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
#include "ccmh.h"

template <int off>
CColMatrixHandler<off>::CColMatrixHandler(std::vector<doublereal>& x,
		const std::vector<int>& i,
		const std::vector<int>& p)
: CompactSparseMatrixHandler(p.size() - 1, p.size() - 1, x, i, p)
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
	std::vector<doublereal> *pax = new std::vector<doublereal>(Ax);
	CColMatrixHandler<off> *p = new CColMatrixHandler<off>(*pax, Ai, Ap);
	p->bMatDuplicate = true;

	return p;
}

template <int off>
int
CColMatrixHandler<off>::MakeCompressedColumnForm(doublereal *const Ax,
		int *const Ai, int *const Ap,
		int offset) const
{
	return Nz();
}

template <int off>
int
CColMatrixHandler<off>::MakeCompressedColumnForm(std::vector<doublereal>& Ax,
		std::vector<int>& Ai, std::vector<int>& Ap,
		int offset) const
{
	return Nz();
}

template <int off>
int
CColMatrixHandler<off>::MakeIndexForm(doublereal *const rAx,
		integer *const Arow, integer *const Acol,
		integer *const AcolSt, int offset) const
{
	silent_cerr("CColMatrixHandler<off>::MakeIndexForm called" << std::endl);
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
int
CColMatrixHandler<off>::MakeIndexForm(std::vector<doublereal>& rAx,
                std::vector<integer>& Arow, std::vector<integer>& Acol,
		std::vector<integer>& AcolSt, int offset) const
{
	silent_cerr("CColMatrixHandler<off>::MakeIndexForm called" << std::endl);

	rAx.resize(Nz());
	Arow.resize(Nz());
	Acol.resize(Nz());
	AcolSt.resize(iGetNumCols());

	return MakeIndexForm(&(rAx[0]), &(Arow[0]), &(Acol[0]),
			&(AcolSt[0]), offset);
}

template <int off>
void
CColMatrixHandler<off>::Resize(integer n, integer nn)
{
	silent_cerr("CColMatrixHandler<off>::Resize called" << std::endl);
	THROW(ErrGeneric());
}

/* Estrae una colonna da una matrice */
template <int off>
VectorHandler&
CColMatrixHandler<off>::GetCol(integer icol, VectorHandler& out) const
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
CColMatrixHandler<off>::MatMatMul(SpMapMatrixHandler& out,
		const SpMapMatrixHandler& in) const
{
	silent_cerr("CColMatrixHandler<off>::MatMatMul called" << std::endl);
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
 * 	for (int col=0; col<NCols; col++) {
 * 		row_cont_type::const_iterator ri, re;
 * 		re = col_indices[col].end();
 * 		for (ri = col_indices[col].begin(); ri!=re; ri++) {
 * 			int iend = in.iGetNumCols();
 * 			for (int col2=0; col2<iend;  col2++) {
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
CColMatrixHandler<off>::MulAndSumWithShift(MatrixHandler& out, doublereal s,
		integer drow, integer dcol) const
{
	silent_cerr("CColMatrixHandler<off>::MulAndSumWithShift called"
			<< std::endl);
	THROW(ErrGeneric());		
	if ((out.iGetNumCols() < iGetNumCols()+dcol)
		|| (out.iGetNumRows() < iGetNumRows()+drow)) {
		silent_cerr("Assertion fault "
				"in CColMatrixHandler<off>::MulAndSumWithShift"
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
CColMatrixHandler<off>::FakeThirdOrderMulAndSumWithShift(MatrixHandler& out, 
		std::vector<bool> b,
		doublereal s,
		integer drow, 
		integer dcol) const
{
	silent_cerr("CColMatrixHandler<off>::FakeThirdOrderMulAndSumWithShift "
			"called" << std::endl);
	THROW(ErrGeneric());		
	if ((out.iGetNumCols() < iGetNumCols()+dcol)
			|| (out.iGetNumRows() < iGetNumRows()+drow)) {
		silent_cerr("Assertion fault "
				"in CColMatrixHandler<off>::MulAndSumWithShift"
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
CColMatrixHandler<off>::MatTVecMul(VectorHandler& out,
		const VectorHandler& in) const
{
	silent_cerr("CColMatrixHandler<off>::MatTVecMul called" << std::endl);
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
CColMatrixHandler<off>::MatVecMul(VectorHandler& out,
		const VectorHandler& in) const
{
	silent_cerr("CColMatrixHandler<off>::MatTVecMul called" << std::endl);
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
CColMatrixHandler<off>::MatTVecIncMul(VectorHandler& out,
		const VectorHandler& in) const
{
	silent_cerr("CColMatrixHandler<off>::MatTVecMul called" << std::endl);
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
CColMatrixHandler<off>::MatVecIncMul(VectorHandler& out,
		const VectorHandler& in) const {
	silent_cerr("CColMatrixHandler<off>::MatTVecMul called" << std::endl);
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
CColMatrixHandler<off>::MatVecDecMul(VectorHandler& out,
		const VectorHandler& in) const
{
	silent_cerr("CColMatrixHandler<off>::MatTVecMul called" << std::endl);
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

template class CColMatrixHandler<0>;
template class CColMatrixHandler<1>;
