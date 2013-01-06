/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2003-2013
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "spmh.h"

SparseMatrixHandler::SparseMatrixHandler(const integer &n, const integer &nn)
:  NRows(n), NCols(nn == 0 ? n : nn), NZ(0)
{
	NO_OP;
}

SparseMatrixHandler::~SparseMatrixHandler(void)
{
	NO_OP;
}

CompactSparseMatrixHandler::CompactSparseMatrixHandler(const integer &n,
		const integer &nn,
		std::vector<doublereal>&x,
		const std::vector<integer>& i,
		const std::vector<integer>& p)
: SparseMatrixHandler(n, nn),
bMatDuplicate(false),
Ax(x),
Ai(i),
Ap(p)
{
	NZ = Ax.size();
}

CompactSparseMatrixHandler::~CompactSparseMatrixHandler()
{
	if (bMatDuplicate) {
		delete &Ax;
	}
}

/* used to sum CC matrices with identical indices */
void
CompactSparseMatrixHandler::AddUnchecked(const CompactSparseMatrixHandler& m)
{
	/* FIXME: put in stl-ish form;
	 * see if we can use something from optimized blas,
	 * e.g. ATLAS, goto or so... */

	/* checks - uncomment to enable */
#ifdef DEBUG
	ASSERT(Ax.size() == m.Ax.size());
	ASSERT(Ai.size() == m.Ai.size());
	ASSERT(Ap.size() == m.Ap.size());
	for (std::vector<doublereal>::size_type i = 0; i < Ai.size(); i++) {
		if (Ai[i] != m.Ai[i]) {
			silent_cerr("AddUnchecked: Ai[" << i << "] differs" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	for (std::vector<doublereal>::size_type i = 0; i < Ap.size(); i++) {
		if (Ap[i] != m.Ap[i]) {
			silent_cerr("AddUnchecked: Ap[" << i << "] differs" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
#endif /* DEBUG */
		
	doublereal *d = &Ax[0], *s = &m.Ax[0];
	std::vector<doublereal>::size_type n = Ax.size();
	for (std::vector<doublereal>::size_type i = 0; i < n; i++) {
		d[i] += s[i];
	}
}

void
CompactSparseMatrixHandler::Reset(void)
{
	std::fill(Ax.begin(), Ax.end(), 0.);
}

integer
CompactSparseMatrixHandler::MakeCompressedColumnForm(doublereal *const Ax,
		integer *const Ai, integer *const Ap,
		int offset) const
{
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	return Nz();
}

integer
CompactSparseMatrixHandler::MakeCompressedColumnForm(std::vector<doublereal>& Ax,
		std::vector<integer>& Ai, std::vector<integer>& Ap,
		int offset) const
{
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	return Nz();
}

integer
CompactSparseMatrixHandler::MakeIndexForm(doublereal *const rAx,
		integer *const Arow, integer *const Acol,
		integer *const AcolSt, int offset) const
{
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	return Nz();
}

integer
CompactSparseMatrixHandler::MakeIndexForm(std::vector<doublereal>& rAx,
                std::vector<integer>& Arow, std::vector<integer>& Acol,
		std::vector<integer>& AcolSt, int offset) const
{
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	return Nz();
}

template <int off>
CompactSparseMatrixHandler_tpl<off>::CompactSparseMatrixHandler_tpl(
		const integer &n,
		const integer &nn,
		std::vector<doublereal>&x,
		const std::vector<integer>& i,
		const std::vector<integer>& p)
: CompactSparseMatrixHandler(n, nn, x, i, p),
m_end(*this, true)
{
	NO_OP;
}

template <int off>
CompactSparseMatrixHandler_tpl<off>::~CompactSparseMatrixHandler_tpl(void)
{
	NO_OP;
}

/* Prodotto Matrice per Matrice */
template <int off>
MatrixHandler&
CompactSparseMatrixHandler_tpl<off>::MatMatMul_base(
	void (MatrixHandler::*op)(integer iRow, integer iCol, const doublereal& dCoef),
	MatrixHandler& out, const MatrixHandler& in) const
{
	ASSERT(in.iGetNumRows() == NCols);
	ASSERT(out.iGetNumRows() == NRows);
	ASSERT(in.iGetNumCols() == out.iGetNumCols());

	// NOTE: out must be zeroed by caller

	integer ncols_in = in.iGetNumCols();

	for (integer col_idx = 1; col_idx <= NCols; col_idx++) {
		integer re = Ap[col_idx] - off;
		integer ri = Ap[col_idx - 1] - off;
		for ( ; ri < re; ri++) {
			integer row_idx = Ai[ri] - off + 1;
			const doublereal& d = Ax[ri];
			for (integer col_in = 1; col_in <= ncols_in; col_in++) {
				(out.*op)(row_idx, col_in, d*in(col_idx, col_in));
			}
		}
	}

	return out;
}

/* Prodotto Matrice trasposta per Matrice */
template <int off>
MatrixHandler&
CompactSparseMatrixHandler_tpl<off>::MatTMatMul_base(
	void (MatrixHandler::*op)(integer iRow, integer iCol, const doublereal& dCoef),
	MatrixHandler& out, const MatrixHandler& in) const
{
	ASSERT(in.iGetNumRows() == NRows);
	ASSERT(out.iGetNumRows() == NCols);
	ASSERT(in.iGetNumCols() == out.iGetNumCols());
	// NOTE: out must be zeroed by caller

	integer ncols_in = in.iGetNumCols();

	for (integer col_idx = 1; col_idx <= NCols; col_idx++) {
		integer re = Ap[col_idx] - off;
		integer ri = Ap[col_idx - 1] - off;
		for ( ; ri < re; ri++) {
			integer row_idx = Ai[ri] - off + 1;
			const doublereal& d = Ax[ri];
			for (integer col_in = 1; col_in <= ncols_in; col_in++) {
				(out.*op)(col_idx, col_in, d*in(row_idx, col_in));
			}
		}
	}

	return out;
}

/* Prodotto Matrice per Vettore */
template <int off>
VectorHandler&
CompactSparseMatrixHandler_tpl<off>::MatVecMul_base(
	void (VectorHandler::*op)(integer iRow, const doublereal& dCoef),
	VectorHandler& out, const VectorHandler& in) const
{
	ASSERT(in.iGetSize() == NCols);
	ASSERT(out.iGetSize() == NRows);

	// NOTE: out must be zeroed by caller

	for (integer col_idx = 1; col_idx <= NCols; col_idx++) {
		integer re = Ap[col_idx] - off;
		integer ri = Ap[col_idx - 1] - off;
		for ( ; ri < re; ri++) {
			integer row_idx = Ai[ri] - off + 1;
			const doublereal& d = Ax[ri];
			(out.*op)(row_idx, d*in(col_idx));
		}
	}

	return out;
}

/* Prodotto Matrice trasposta per Vettore */
template <int off>
VectorHandler&
CompactSparseMatrixHandler_tpl<off>::MatTVecMul_base(
	void (VectorHandler::*op)(integer iRow,
	const doublereal& dCoef),
	VectorHandler& out, const VectorHandler& in) const
{
	ASSERT(in.iGetSize() == NRows);
	ASSERT(out.iGetSize() == NCols);

	// NOTE: out must be zeroed by caller

	for (integer col_idx = 1; col_idx <= NCols; col_idx++) {
		integer re = Ap[col_idx] - off;
		integer ri = Ap[col_idx - 1] - off;
		for ( ; ri < re; ri++) {
			integer row_idx = Ai[ri] - off + 1;
			const doublereal& d = Ax[ri];
			(out.*op)(col_idx, d*in(row_idx));
		}
	}

	return out;
}

template class CompactSparseMatrixHandler_tpl<0>;
template class CompactSparseMatrixHandler_tpl<1>;

template <int off>
void
CompactSparseMatrixHandler_const_iterator<off>::reset(bool is_end)
{
	if (is_end) {
		elem.iRow = m.iGetNumRows();
		elem.iCol = m.iGetNumCols();

	} else {
		elem.iCol = 0;
		while (m.Ap[elem.iCol + 1] - off == m.Ap[elem.iCol] - off) {
			if (++elem.iCol == m.iGetNumCols()) {
				elem.iRow = m.iGetNumRows();
				return;
			}
		}
		i_idx = m.Ap[elem.iCol] - off;
		elem.iRow = m.Ai[i_idx] - off;
		elem.dCoef = m.Ax[i_idx];
	}
}

template <int off>
CompactSparseMatrixHandler_const_iterator<off>::CompactSparseMatrixHandler_const_iterator(const CompactSparseMatrixHandler_tpl<off>& m, bool is_end)
: m(m)
{
	reset(is_end);
}

template <int off>
CompactSparseMatrixHandler_const_iterator<off>::~CompactSparseMatrixHandler_const_iterator(void)
{
	NO_OP;
}

template <int off>
const CompactSparseMatrixHandler_const_iterator<off>&
CompactSparseMatrixHandler_const_iterator<off>::operator ++ (void) const
{
	i_idx++;
	while (i_idx == m.Ap[elem.iCol + 1] - off) {
		if (++elem.iCol == m.iGetNumCols()) {
			elem.iRow = m.iGetNumRows();
			return *this;
		}
		i_idx = m.Ap[elem.iCol] - off;
	}
	elem.iRow = m.Ai[i_idx] - off;
	elem.dCoef = m.Ax[i_idx];

	return *this;
}

template <int off>
const SparseMatrixHandler::SparseMatrixElement*
CompactSparseMatrixHandler_const_iterator<off>::operator -> (void)
{
	return &elem;
}

template <int off>
const SparseMatrixHandler::SparseMatrixElement&
CompactSparseMatrixHandler_const_iterator<off>::operator * (void)
{
	return elem;
}

template <int off>
bool
CompactSparseMatrixHandler_const_iterator<off>::operator == (const CompactSparseMatrixHandler_const_iterator<off>& op) const
{
	return elem == op.elem;
}

template <int off>
bool
CompactSparseMatrixHandler_const_iterator<off>::operator != (const CompactSparseMatrixHandler_const_iterator<off>& op) const
{
	return elem != op.elem;
}

template class CompactSparseMatrixHandler_const_iterator<0>;
template class CompactSparseMatrixHandler_const_iterator<1>;
