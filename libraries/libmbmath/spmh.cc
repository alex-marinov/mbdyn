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

#include <algorithm>
#include <iomanip>

#include "spmh.h"

SparseMatrixHandler::SparseMatrixHandler(const integer &n, const integer &nn)
:  NRows(n), NCols(nn == 0 ? n : nn)
{
	NO_OP;
}

SparseMatrixHandler::~SparseMatrixHandler(void)
{
	NO_OP;
}

int32_t SparseMatrixHandler::MakeCompressedColumnForm(doublereal* Ax,
						      int32_t* Ai,
						      int32_t* Ap,
						      int offset) const     
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

int64_t SparseMatrixHandler::MakeCompressedColumnForm(doublereal* Ax,
						      int64_t* Ai,
						      int64_t* Ap,
						      int offset) const     
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}


int32_t SparseMatrixHandler::MakeCompressedColumnForm(std::vector<doublereal>& Ax,
						      std::vector<int32_t>& Ai,
						      std::vector<int32_t>& Ap,
						      int offset) const     
{
     Ax.resize(Nz());
     Ai.resize(Nz());
     Ap.resize(iGetNumCols() + 1);

     return MakeCompressedColumnForm(&Ax.front(), &Ai.front(), &Ap.front(), offset);
}

int64_t SparseMatrixHandler::MakeCompressedColumnForm(std::vector<doublereal>& Ax,
						      std::vector<int64_t>& Ai,
						      std::vector<int64_t>& Ap,
						      int offset) const     
{
     Ax.resize(Nz());
     Ai.resize(Nz());
     Ap.resize(iGetNumCols() + 1);

     return MakeCompressedColumnForm(&Ax.front(), &Ai.front(), &Ap.front(), offset);
}


int32_t SparseMatrixHandler::MakeCompressedRowForm(doublereal *const Ax,
						   int32_t *const Ai,
						   int32_t *const Ap,
						   int offset) const
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

int32_t SparseMatrixHandler::MakeCompressedRowForm(std::vector<doublereal>& Ax,
						   std::vector<int32_t>& Ai,
						   std::vector<int32_t>& Ap,
						   int offset) const
{
     Ax.resize(Nz());
     Ai.resize(Nz());
     Ap.resize(iGetNumRows() + 1);
     
     return MakeCompressedRowForm(&Ax.front(), &Ai.front(), &Ap.front(), offset);
}

int64_t SparseMatrixHandler::MakeCompressedRowForm(doublereal *const Ax,
						   int64_t *const Ai,
						   int64_t *const Ap,
						   int offset) const
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

int64_t SparseMatrixHandler::MakeCompressedRowForm(std::vector<doublereal>& Ax,
						   std::vector<int64_t>& Ai,
						   std::vector<int64_t>& Ap,
						   int offset) const
{
     Ax.resize(Nz());
     Ai.resize(Nz());
     Ap.resize(iGetNumRows() + 1);
     
     return MakeCompressedRowForm(&Ax.front(), &Ai.front(), &Ap.front(), offset);     
}


	
int32_t SparseMatrixHandler::MakeIndexForm(doublereal *const Ax,
					   int32_t *const Arow, int32_t *const Acol,
					   int32_t *const AcolSt,
					   int offset) const
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

	
int32_t SparseMatrixHandler::MakeIndexForm(std::vector<doublereal>& Ax,
					   std::vector<int32_t>& Arow, std::vector<int32_t>& Acol,
					   std::vector<int32_t>& Ap,
					   int offset) const
{
     Ax.resize(Nz());
     Arow.resize(Nz());
     Acol.resize(Nz());
     Ap.resize(iGetNumCols() + 1);

     return MakeIndexForm(&Ax.front(), &Arow.front(), &Acol.front(), &Ap.front(), offset);
}

	
int64_t SparseMatrixHandler::MakeIndexForm(doublereal *const Ax,
					   int64_t *const Arow, int64_t *const Acol,
					   int64_t *const AcolSt,
					   int offset) const
{
     throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

	
int64_t SparseMatrixHandler::MakeIndexForm(std::vector<doublereal>& Ax,
					   std::vector<int64_t>& Arow, std::vector<int64_t>& Acol,
					   std::vector<int64_t>& Ap,
					   int offset) const
{
     Ax.resize(Nz());
     Arow.resize(Nz());
     Acol.resize(Nz());
     Ap.resize(iGetNumCols() + 1);

     return MakeIndexForm(&Ax.front(), &Arow.front(), &Acol.front(), &Ap.front(), offset);
}

CompactSparseMatrixHandler::CompactSparseMatrixHandler(const integer &n,
						       const integer &nn)
: SparseMatrixHandler(n, nn)
{
	
}

CompactSparseMatrixHandler::~CompactSparseMatrixHandler()
{

}

template <int off, typename idx_type>
CompactSparseMatrixHandler_tpl<off, idx_type>::CompactSparseMatrixHandler_tpl(const integer &n,
									      const integer &nn,
									      std::vector<doublereal>&x,
									      const std::vector<idx_type>& i,
									      const std::vector<idx_type>& p)
     : CompactSparseMatrixHandler(n, nn),
       bMatDuplicate(false),
       Ax(x),
       Ai(i),
       Ap(p),
       NZ(x.size()),       
       m_end(*this, true)
{
	NO_OP;
}

template <int off, typename idx_type>
CompactSparseMatrixHandler_tpl<off, idx_type>::~CompactSparseMatrixHandler_tpl(void)
{
     if (bMatDuplicate) {
	  delete &Ax;
     }
}

template <int off, typename idx_type>
const doublereal* CompactSparseMatrixHandler_tpl<off, idx_type>::pdGetMat(void) const {
     return &Ax.front();
};

template <int off, typename idx_type>
void
CompactSparseMatrixHandler_tpl<off, idx_type>::Reset(void)
{
	std::fill(Ax.begin(), Ax.end(), 0.);
}

template <int off, typename idx_type>
integer CompactSparseMatrixHandler_tpl<off, idx_type>::Nz() const
{
     return NZ;
}

/* used to sum CC matrices with identical indices */
template <int off, typename idx_type>
void
CompactSparseMatrixHandler_tpl<off, idx_type>::AddUnchecked(const CompactSparseMatrixHandler& m)
{
	/* FIXME: put in stl-ish form;
	 * see if we can use something from optimized blas,
	 * e.g. ATLAS, goto or so... */

	/* checks - uncomment to enable */
#ifdef DEBUG
        const auto& m2 = dynamic_cast<const CompactSparseMatrixHandler_tpl<off, idx_type>&>(m);
	ASSERT(Ax.size() == m2.Ax.size());
	ASSERT(Ai.size() == m2.Ai.size());
	ASSERT(Ap.size() == m2.Ap.size());
	for (std::vector<doublereal>::size_type i = 0; i < Ai.size(); i++) {
		if (Ai[i] != m2.Ai[i]) {
			silent_cerr("AddUnchecked: Ai[" << i << "] differs" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	for (std::vector<doublereal>::size_type i = 0; i < Ap.size(); i++) {
		if (Ap[i] != m2.Ap[i]) {
			silent_cerr("AddUnchecked: Ap[" << i << "] differs" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
#endif /* DEBUG */
		
	doublereal *d = &Ax.front();
	const doublereal* s = m.pdGetMat();
	std::vector<doublereal>::size_type n = Ax.size();
	for (std::vector<doublereal>::size_type i = 0; i < n; i++) {
		d[i] += s[i];
	}
}

/* Prodotto Matrice per Matrice */
template <int off, typename idx_type>
MatrixHandler&
CompactSparseMatrixHandler_tpl<off, idx_type>::MatMatMul_base(
	void (MatrixHandler::*op)(integer iRow, integer iCol, const doublereal& dCoef),
	MatrixHandler& out, const MatrixHandler& in) const
{
	ASSERT(in.iGetNumRows() == NCols);
	ASSERT(out.iGetNumRows() == NRows);
	ASSERT(in.iGetNumCols() == out.iGetNumCols());

	// NOTE: out must be zeroed by caller

	integer ncols_in = in.iGetNumCols();

	for (integer col_idx = 1; col_idx <= this->NCols; col_idx++) {
		idx_type re = this->Ap[col_idx] - off;
		idx_type ri = this->Ap[col_idx - 1] - off;
		for ( ; ri < re; ri++) {
			idx_type row_idx = this->Ai[ri] - off + 1;
			const doublereal& d = this->Ax[ri];
			for (integer col_in = 1; col_in <= ncols_in; col_in++) {
				(out.*op)(row_idx, col_in, d*in(col_idx, col_in));
			}
		}
	}

	return out;
}

/* Prodotto Matrice trasposta per Matrice */
template <int off, typename idx_type>
MatrixHandler&
CompactSparseMatrixHandler_tpl<off, idx_type>::MatTMatMul_base(
	void (MatrixHandler::*op)(integer iRow, integer iCol, const doublereal& dCoef),
	MatrixHandler& out, const MatrixHandler& in) const
{
	ASSERT(in.iGetNumRows() == NRows);
	ASSERT(out.iGetNumRows() == NCols);
	ASSERT(in.iGetNumCols() == out.iGetNumCols());
	// NOTE: out must be zeroed by caller

	integer ncols_in = in.iGetNumCols();

	for (integer col_idx = 1; col_idx <= this->NCols; col_idx++) {
		idx_type re = this->Ap[col_idx] - off;
		idx_type ri = this->Ap[col_idx - 1] - off;
		for ( ; ri < re; ri++) {
			idx_type row_idx = this->Ai[ri] - off + 1;
			const doublereal& d = this->Ax[ri];
			for (integer col_in = 1; col_in <= ncols_in; col_in++) {
				(out.*op)(col_idx, col_in, d*in(row_idx, col_in));
			}
		}
	}

	return out;
}

/* Prodotto Matrice per Vettore */
template <int off, typename idx_type>
VectorHandler&
CompactSparseMatrixHandler_tpl<off, idx_type>::MatVecMul_base(
	void (VectorHandler::*op)(integer iRow, const doublereal& dCoef),
	VectorHandler& out, const VectorHandler& in) const
{
	ASSERT(in.iGetSize() == NCols);
	ASSERT(out.iGetSize() == NRows);

	// NOTE: out must be zeroed by caller

	for (integer col_idx = 1; col_idx <= this->NCols; col_idx++) {
		idx_type re = this->Ap[col_idx] - off;
		idx_type ri = this->Ap[col_idx - 1] - off;
		for ( ; ri < re; ri++) {
			idx_type row_idx = this->Ai[ri] - off + 1;
			const doublereal& d = this->Ax[ri];
			(out.*op)(row_idx, d*in(col_idx));
		}
	}

	return out;
}

/* Prodotto Matrice trasposta per Vettore */
template <int off, typename idx_type>
VectorHandler&
CompactSparseMatrixHandler_tpl<off, idx_type>::MatTVecMul_base(
	void (VectorHandler::*op)(integer iRow,
	const doublereal& dCoef),
	VectorHandler& out, const VectorHandler& in) const
{
	ASSERT(in.iGetSize() == NRows);
	ASSERT(out.iGetSize() == NCols);

	// NOTE: out must be zeroed by caller

	for (integer col_idx = 1; col_idx <= this->NCols; col_idx++) {
		idx_type re = this->Ap[col_idx] - off;
		idx_type ri = this->Ap[col_idx - 1] - off;
		for ( ; ri < re; ri++) {
			idx_type row_idx = this->Ai[ri] - off + 1;
			const doublereal& d = this->Ax[ri];
			(out.*op)(col_idx, d*in(row_idx));
		}
	}

	return out;
}

template <int off, typename idx_type>
void
CompactSparseMatrixHandler_tpl<off, idx_type>::const_iterator::reset(bool is_end)
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

template <int off, typename idx_type>
CompactSparseMatrixHandler_tpl<off, idx_type>::const_iterator::const_iterator(const CompactSparseMatrixHandler_tpl<off, idx_type>& m, bool is_end)
: m(m)
{
	reset(is_end);
}

template <int off, typename idx_type>
CompactSparseMatrixHandler_tpl<off, idx_type>::const_iterator::~const_iterator(void)
{
	NO_OP;
}

template <int off, typename idx_type>
const typename CompactSparseMatrixHandler_tpl<off, idx_type>::const_iterator&
CompactSparseMatrixHandler_tpl<off, idx_type>::const_iterator::operator ++ (void) const
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

template <int off, typename idx_type>
const SparseMatrixHandler::SparseMatrixElement*
CompactSparseMatrixHandler_tpl<off, idx_type>::const_iterator::operator -> (void)
{
	return &elem;
}

template <int off, typename idx_type>
const SparseMatrixHandler::SparseMatrixElement&
CompactSparseMatrixHandler_tpl<off, idx_type>::const_iterator::operator * (void)
{
	return elem;
}

template <int off, typename idx_type>
bool
CompactSparseMatrixHandler_tpl<off, idx_type>::const_iterator::operator == (const CompactSparseMatrixHandler_tpl<off, idx_type>::const_iterator& op) const
{
	return elem == op.elem;
}

template <int off, typename idx_type>
bool
CompactSparseMatrixHandler_tpl<off, idx_type>::const_iterator::operator != (const CompactSparseMatrixHandler_tpl<off, idx_type>::const_iterator& op) const
{
	return elem != op.elem;
}

template <int off, typename idx_type>
void CompactSparseMatrixHandler_tpl<off, idx_type>::Scale(const std::vector<doublereal>& oRowScale, const std::vector<doublereal>& oColScale)
{
     this->IteratorScale(*this, oRowScale, oColScale);
}

template class CompactSparseMatrixHandler_tpl<0, int32_t>;
template class CompactSparseMatrixHandler_tpl<1, int32_t>;
template class CompactSparseMatrixHandler_tpl<1, int64_t>;

