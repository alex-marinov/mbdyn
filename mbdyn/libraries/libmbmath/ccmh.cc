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
		const std::vector<integer>& i,
		const std::vector<integer>& p)
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
void
CColMatrixHandler<off>::Resize(integer n, integer nn)
{
	silent_cerr("CColMatrixHandler<off>::Resize called" << std::endl);
	throw ErrGeneric();
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
		throw ErrGeneric();
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
MatrixHandler*
CColMatrixHandler<off>::MatMatMul_base(void (MatrixHandler::*op)(integer iRow,
			integer iCol, const doublereal& dCoef),
			MatrixHandler* out, const MatrixHandler& in) const
{
	silent_cerr("CColMatrixHandler<off>::MatMatMul_base called"
			<< std::endl);
	throw ErrGeneric();		
}

template <int off>
MatrixHandler*
CColMatrixHandler<off>::MatTMatMul_base(void (MatrixHandler::*op)(integer iRow,
			integer iCol, const doublereal& dCoef),
			MatrixHandler* out, const MatrixHandler& in) const
{
	silent_cerr("CColMatrixHandler<off>::MatTMatMul_base called"
			<< std::endl);
	throw ErrGeneric();		
}

/* Moltiplica per uno scalare e somma a una matrice */
template <int off>
MatrixHandler&
CColMatrixHandler<off>::MulAndSumWithShift(MatrixHandler& out, doublereal s,
		integer drow, integer dcol) const
{
	silent_cerr("CColMatrixHandler<off>::MulAndSumWithShift called"
			<< std::endl);
	throw ErrGeneric();		
	if ((out.iGetNumCols() < iGetNumCols()+dcol)
		|| (out.iGetNumRows() < iGetNumRows()+drow)) {
		silent_cerr("Assertion fault "
				"in CColMatrixHandler<off>::MulAndSumWithShift"
				<< std::endl);
		throw ErrGeneric();
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
	throw ErrGeneric();		
	if ((out.iGetNumCols() < iGetNumCols()+dcol)
			|| (out.iGetNumRows() < iGetNumRows()+drow)) {
		silent_cerr("Assertion fault "
				"in CColMatrixHandler<off>::MulAndSumWithShift"
				<< std::endl);
		throw ErrGeneric();
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

/* Prodotto Matrice per Vettore */
template <int off>
VectorHandler&
CColMatrixHandler<off>::MatVecMul_base(void (VectorHandler::*op)(integer iRow,
			const doublereal& dCoef),
			VectorHandler& out, const VectorHandler& in) const
{
	silent_cerr("CColMatrixHandler<off>::MatVecMul_base called"
			<< std::endl);
	throw ErrGeneric();		
}

template <int off>
VectorHandler&
CColMatrixHandler<off>::MatTVecMul_base(void (VectorHandler::*op)(integer iRow,
			const doublereal& dCoef),
			VectorHandler& out, const VectorHandler& in) const
{
	silent_cerr("CColMatrixHandler<off>::MatTVecMul_base called"
			<< std::endl);
	throw ErrGeneric();		
}

template class CColMatrixHandler<0>;
template class CColMatrixHandler<1>;
