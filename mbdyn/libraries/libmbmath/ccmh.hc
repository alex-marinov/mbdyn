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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "spmapmh.h"
#include "ccmh.h"

CColMatrixHandler::CColMatrixHandler(const int &n,
		std::vector<doublereal>& x,
		const std::vector<int>& i,
		const std::vector<int>& p)
: CompactSparseMatrixHandler(n, n, x, i, p)
{
	ASSERT(n == Ap.size() - 1);
}

CColMatrixHandler::~CColMatrixHandler()
{
	NO_OP;
}

/* used by MultiThreadDataManager to duplicate the storage array
 * while preserving the CC indices */
CompactSparseMatrixHandler *
CColMatrixHandler::Copy(void) const
{
	std::vector<doublereal> *pax = new std::vector<doublereal>(Ax);
	CColMatrixHandler *p = new CColMatrixHandler(iGetNumRows(),
			*pax, Ai, Ap);
	p->bMatDuplicate = true;

	return p;
}

int
CColMatrixHandler::MakeCompressedColumnForm( doublereal *const Ax,
		int *const Ai,
		int *const Ap,
		integer offset) const
{
	silent_cerr("CColMatrixHandler::MakeCompressedColumnForm called"
			<< std::endl);
	THROW(ErrGeneric());		
	return 0;
}

int
CColMatrixHandler::MakeCompressedColumnForm( std::vector<doublereal>& Ax,
                std::vector<int>& Ai,
                std::vector<int>& Ap,
		integer offset) const
{
	silent_cerr("CColMatrixHandler::MakeCompressedColumnForm called"
			<< std::endl);
	THROW(ErrGeneric());		
	return 0;
}

int
CColMatrixHandler::MakeIndexForm( doublereal *const rAx,
		integer *const Arow,
		integer *const Acol,
		integer offset) const
{
	silent_cerr("CColMatrixHandler::MakeIndexForm called" << std::endl);
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
}

int
CColMatrixHandler::MakeIndexForm( std::vector<doublereal>& rAx,
                std::vector<integer>& Arow,
                std::vector<integer>& Acol,
		integer offset) const
{
	silent_cerr("CColMatrixHandler::MakeIndexForm called" << std::endl);
	THROW(ErrGeneric());		

	rAx.resize(Nz());
	Arow.resize(Nz());
	Acol.resize(Nz());
	return MakeIndexForm(&(rAx[0]), &(Arow[0]), &(Acol[0]), offset);
}

void
CColMatrixHandler::Resize(const int &n, const int &nn)
{
	silent_cerr("CColMatrixHandler::Resize called" << std::endl);
	THROW(ErrGeneric());
}

/* Estrae una colonna da una matrice */
VectorHandler&
CColMatrixHandler::GetCol(integer icol, VectorHandler& out) const
{
	silent_cerr("CColMatrixHandler::GetCol called" << std::endl);
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
}
	
/* Prodotto Matrice per Matrice */
SpMapMatrixHandler& 
CColMatrixHandler::MatMatMul(SpMapMatrixHandler& out,
		const SpMapMatrixHandler& in) const
{
	silent_cerr("CColMatrixHandler::MatMatMul called" << std::endl);
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
MatrixHandler&
CColMatrixHandler::MulAndSumWithShift(MatrixHandler& out, doublereal s,
		integer drow, integer dcol) const
{
	silent_cerr("CColMatrixHandler::MulAndSumWithShift called"
			<< std::endl);
	THROW(ErrGeneric());		
	if ((out.iGetNumCols() < iGetNumCols()+dcol)
		|| (out.iGetNumRows() < iGetNumRows()+drow)) {
		silent_cerr("Assertion fault "
				"in CColMatrixHandler::MulAndSumWithShift"
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

MatrixHandler&
CColMatrixHandler::FakeThirdOrderMulAndSumWithShift(MatrixHandler& out, 
		std::vector<bool> b,
		doublereal s,
		integer drow, 
		integer dcol) const
{
	silent_cerr("CColMatrixHandler::FakeThirdOrderMulAndSumWithShift "
			"called" << std::endl);
	THROW(ErrGeneric());		
	if ((out.iGetNumCols() < iGetNumCols()+dcol)
			|| (out.iGetNumRows() < iGetNumRows()+drow)) {
		silent_cerr("Assertion fault "
				"in CColMatrixHandler::MulAndSumWithShift"
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
	
VectorHandler&
CColMatrixHandler::MatTVecMul(VectorHandler& out,
		const VectorHandler& in) const
{
	silent_cerr("CColMatrixHandler::MatTVecMul called" << std::endl);
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
	
VectorHandler&
CColMatrixHandler::MatVecMul(VectorHandler& out,
		const VectorHandler& in) const
{
	silent_cerr("CColMatrixHandler::MatTVecMul called" << std::endl);
	THROW(ErrGeneric());		
	if (in.iGetSize() != iGetNumCols() 
			|| out.iGetSize() != iGetNumRows()) {
		THROW(ErrGeneric());
  	}

	out.Reset(0.);
	return MatVecIncMul(out, in);
}

VectorHandler&
CColMatrixHandler::MatTVecIncMul(VectorHandler& out,
		const VectorHandler& in) const
{
	silent_cerr("CColMatrixHandler::MatTVecMul called" << std::endl);
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
	
VectorHandler&
CColMatrixHandler::MatVecIncMul(VectorHandler& out,
		const VectorHandler& in) const {
	silent_cerr("CColMatrixHandler::MatTVecMul called" << std::endl);
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

VectorHandler&
CColMatrixHandler::MatVecDecMul(VectorHandler& out,
		const VectorHandler& in) const
{
	silent_cerr("CColMatrixHandler::MatTVecMul called" << std::endl);
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

