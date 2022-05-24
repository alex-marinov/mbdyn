/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2003-2017
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "spmapmh.h"
#include "dirccmh.h"

template <int off, typename idx_type>
DirCColMatrixHandler<off, idx_type>::DirCColMatrixHandler(std::vector<doublereal>& x,
							  const std::vector<idx_type>& i,
							  const std::vector<idx_type>& p)
: CompactSparseMatrixHandler_tpl<off, idx_type>(p.size() - 1, p.size() - 1, x, i, p),
  pindices(this->iGetNumCols() + 1),
  indices(static_cast<size_t>(this->iGetNumRows())*this->iGetNumCols(), -1)
{
	for (integer col = 1; col <= this->iGetNumCols(); col++) {
                pindices[col] = &indices[static_cast<size_t>(col - 1)*this->iGetNumRows()] - 1;

		auto row_begin = p[col - 1] - off, row_end = p[col] - off;

		for (auto r = row_begin; r < row_end; r++) {
			pindices[col][i[r] - off + 1] = r;
		}
	}
}

template <int off, typename idx_type>
DirCColMatrixHandler<off, idx_type>::~DirCColMatrixHandler()
{
	NO_OP;
}

/* used by MultiThreadDataManager to duplicate the storage array
 * while preserving the CC indices */
template <int off, typename idx_type>
CompactSparseMatrixHandler*
DirCColMatrixHandler<off, idx_type>::Copy(void) const
{
	std::vector<doublereal> *pax =
	     new std::vector<doublereal>(this->Ax);
	DirCColMatrixHandler<off, idx_type> *p =
	     new DirCColMatrixHandler<off, idx_type>(*pax, this->Ai,
						     this->Ap);
	p->bMatDuplicate = true;

	return p;
}

template <int off, typename idx_type>
void
DirCColMatrixHandler<off, idx_type>::Resize(integer n, integer nn)
{
	silent_cerr("DirCColMatrixHandler<off>::Resize called" << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

/* Estrae una colonna da una matrice */
template <int off, typename idx_type>
VectorHandler&
DirCColMatrixHandler<off, idx_type>::GetCol(integer icol, VectorHandler& out) const
{
	/*
	 * Note: we assume out has been reset
	 */
	
        if (icol > this->iGetNumCols()) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	auto idx = this->Ap[icol - 1] - off;
	auto idxe = this->Ap[icol] - off;

	for ( ; idx < idxe; idx++) {
	     out(this->Ai[idx] - off + 1) =
		  this->Ax[idx];
	}

	return out;
}
	
/* Moltiplica per uno scalare e somma a una matrice */
template <int off, typename idx_type>
MatrixHandler&
DirCColMatrixHandler<off, idx_type>::MulAndSumWithShift(MatrixHandler& out, doublereal s,
		integer drow, integer dcol) const
{
	silent_cerr("DirCColMatrixHandler<off>::MulAndSumWithShift called"
			<< std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);		
	if ((out.iGetNumCols() < this->iGetNumCols() + dcol)
		|| (out.iGetNumRows() < this->iGetNumRows() + drow))
	{
		silent_cerr("Assertion fault "
				"in DirCColMatrixHandler<off>::MulAndSumWithShift"
				<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	drow = drow + 1;
	for (integer col = 0; col < this->iGetNumCols(); col++) {
	     auto idx = this->Ap[col] - off;
	     auto idxe = this->Ap[col+1] - off;
	     auto newcol = col + dcol + 1;
	     for (; idx < idxe; idx++) {
		  out.IncCoef(this->Ai[idx] - off + drow,
			      newcol, this->Ax[idx]*s);
	     }
	}
	return out;
}

template class DirCColMatrixHandler<0, int32_t>;
template class DirCColMatrixHandler<1, int32_t>;
template class DirCColMatrixHandler<1, int64_t>;
