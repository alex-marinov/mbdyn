/* $Header$ */
/*
 * HmFe (C) is a FEM analysis code.
 *
 * Copyright (C) 1996-2017
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
 * Modified to add methods
 * to be used in the parallel MBDyn Solver.
 *
 * Copyright (C) 2003-2017
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

#include "spmapmh.h"

SpMapMatrixHandler::SpMapMatrixHandler(const integer &n, const integer &nn)
: SparseMatrixHandler(n, nn), m_end(*this, true), NZ(0)
{
	col_indices.resize(NCols);
}

SpMapMatrixHandler::~SpMapMatrixHandler()
{
	NO_OP;
}

template <typename idx_type>

idx_type SpMapMatrixHandler::MakeCompressedColumnFormTpl(doublereal *const Ax,
							 idx_type *const Ai,
							 idx_type *const Ap,
							 int offset) const
{
	idx_type x_ptr = 0;

	row_cont_type::iterator ri;
	row_cont_type::const_iterator re;

	for (integer col = 0; col < NCols; col++) {
		Ap[col] = x_ptr + offset;
		re = col_indices[col].end();
		for (ri = col_indices[col].begin(); ri != re; ++ri) {
			Ax[x_ptr] = ri->second;
			Ai[x_ptr] = ri->first + offset;
			x_ptr++;
		}
	}

	ASSERTMSGBREAK(x_ptr == NZ, "Error in "
			"SpMapMatrixHandler::MakeCompressedColumnForm");

	Ap[NCols] = x_ptr + offset;
	
	return x_ptr;     
}

integer
SpMapMatrixHandler::MakeCompressedColumnForm(doublereal *const Ax,
					     integer *const Ai,
					     integer *const Ap,
					     int offset) const
{
     return MakeCompressedColumnFormTpl(Ax, Ai, Ap, offset);
}

int64_t
SpMapMatrixHandler::MakeCompressedColumnForm(doublereal *const Ax,
					     int64_t *const Ai,
					     int64_t *const Ap,
					     int offset) const
{
     return MakeCompressedColumnFormTpl(Ax, Ai, Ap, offset);
}

integer
SpMapMatrixHandler::MakeIndexForm(doublereal *const Ax,
		integer *const Arow, integer *const Acol,
		integer *const Ap, int offset) const
{
	integer x_ptr = 0;

	row_cont_type::iterator ri;
	row_cont_type::const_iterator re;

	for (integer col = 0; col < NCols; col++) {
		Ap[col] = x_ptr + offset;
		re = col_indices[col].end();
		for (ri = col_indices[col].begin(); ri != re; ++ri) {
			Ax[x_ptr] = ri->second;
			Arow[x_ptr] = ri->first + offset;
			Acol[x_ptr] = col + offset;
			x_ptr++;
		}
	}

	ASSERTMSGBREAK(x_ptr == NZ, "Error in "
			"SpMapMatrixHandler::MakeIndexForm");

	Ap[NCols] = x_ptr + offset;

	return Nz();
}

integer
SpMapMatrixHandler::MakeIndexForm(std::vector<doublereal>& Ax,
                std::vector<integer>& Arow, std::vector<integer>& Acol,
		std::vector<integer>& Ap, int offset) const
{
	Ax.resize(Nz());
	Arow.resize(Nz());
	Acol.resize(Nz());
	Ap.resize(iGetNumCols() + 1);

	return MakeIndexForm(&Ax[0], &Arow[0], &Acol[0], &Ap[0], offset);
}

void
SpMapMatrixHandler::Reset(void)
{
	row_cont_type::const_iterator re;
	row_cont_type::iterator ri;
	for (integer col = 0; col < NCols; col++) {
		re = col_indices[col].end();
		for (ri = col_indices[col].begin(); ri != re; ++ri) {
			ri->second = 0.;
		}
	}
}

void
SpMapMatrixHandler::Resize(integer n, integer nn)
{
	if (nn == 0) {
		nn = n;
	}

	for (integer col = 0; col < NCols; col++) {
		col_indices[col].clear();
	}

	col_indices.resize(nn);
	NRows = n;
	NCols = nn;
	NZ = 0;

	m_end.reset();
}

/* Estrae una colonna da una matrice */
VectorHandler&
SpMapMatrixHandler::GetCol(integer icol, VectorHandler& out) const
{
        if (icol > iGetNumCols()) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	row_cont_type::const_iterator ri, re;
	re = col_indices[icol].end();
	for (ri = col_indices[icol].begin(); ri != re; ++ri) {
		out(ri->first + 1) = ri->second;
	}
	return out;
}

void SpMapMatrixHandler::Scale(const std::vector<doublereal>& oRowScale, const std::vector<doublereal>& oColScale)
{
     IteratorScale(*this, oRowScale, oColScale);
}

integer SpMapMatrixHandler::Nz() const
{
     return NZ;
}

/* Prodotto Matrice per Matrice */
MatrixHandler&
SpMapMatrixHandler::MatMatMul_base(void (MatrixHandler::*op)(integer iRow,
			integer iCol, const doublereal& dCoef),
			MatrixHandler& out, const MatrixHandler& in) const
{
	if ((in.iGetNumRows() != iGetNumCols())
			|| (in.iGetNumCols() != out.iGetNumCols())
			|| (out.iGetNumRows() != iGetNumRows()))
	{
		silent_cerr("Assertion fault "
			"in SpMapMatrixHandler::MatMatIncMul" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	integer ncols_in = in.iGetNumCols();
	for (integer row_in = 0; row_in < NCols; row_in++) {
		row_cont_type::const_iterator ri, re;
		re = col_indices[row_in].end();
		for (ri = col_indices[row_in].begin(); ri != re; ++ri) {
			for (integer col_in = 1; col_in <= ncols_in;  col_in++) {
				(out.*op)(ri->first + 1, col_in,
						ri->second*in(row_in + 1, col_in));
			}
		}
	}

	return out;
}

MatrixHandler&
SpMapMatrixHandler::MatTMatMul_base(void (MatrixHandler::*op)(integer iRow,
			integer iCol, const doublereal& dCoef),
			MatrixHandler& out, const MatrixHandler& in) const
{
	if ((in.iGetNumRows() != iGetNumCols())
			|| (in.iGetNumCols() != out.iGetNumCols())
			|| (out.iGetNumRows() != iGetNumRows()))
	{
		silent_cerr("Assertion fault "
			"in SpMapMatrixHandler::MatTMatMul_base" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	integer ncols_in = in.iGetNumCols();
	for (integer row_out = 0; row_out < NCols; row_out++) {
		row_cont_type::const_iterator ri, re;
		re = col_indices[row_out].end();
		for (ri = col_indices[row_out].begin(); ri != re; ++ri) {
			for (integer col_in = 1; col_in <= ncols_in; col_in++) {
				(out.*op)(row_out + 1, col_in,
						ri->second*in(ri->first + 1, col_in));
			}
		}
	}

	return out;
}

/* Moltiplica per uno scalare e somma a una matrice */
MatrixHandler &
SpMapMatrixHandler::MulAndSumWithShift(MatrixHandler& out, doublereal s ,
		integer drow, integer dcol) const
{
	if ((out.iGetNumCols() < iGetNumCols() + dcol)
			|| (out.iGetNumRows() < iGetNumRows() + drow))
	{
		silent_cerr("Assertion fault "
			"in SpMapMatrixHandler::MulAndSumWithShift"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	drow = drow + 1;

	for (integer col = 0; col < NCols; col++) {
		row_cont_type::const_iterator ri, re;
		re = col_indices[col].end();
		integer newcol = col + dcol + 1;
		for (ri = col_indices[col].begin(); ri != re; ++ri) {
			out.IncCoef(ri->first + drow, newcol, ri->second*s);
		}
	}

	return out;
}

MatrixHandler &
SpMapMatrixHandler::FakeThirdOrderMulAndSumWithShift(MatrixHandler& out,
		std::vector<bool> b, doublereal s, integer drow,
		integer dcol) const
{
	if ((out.iGetNumCols() < iGetNumCols() + dcol)
			|| (out.iGetNumRows() < iGetNumRows() + drow))
	{
		silent_cerr("Assertion fault "
			"in SpMapMatrixHandler::MulAndSumWithShift"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	drow = drow + 1;

	for (integer col = 0; col < NCols; col++) {
		row_cont_type::const_iterator ri, re;
		re = col_indices[col].end();
		integer newcol = col + dcol + 1;
		for (ri = col_indices[col].begin(); ri != re; ++ri) {
			if (b[ri->first]) {
				out.IncCoef(ri->first + drow,
						newcol, ri->second*s);
			}
		}
	}

	return out;
}

VectorHandler &
SpMapMatrixHandler::MatVecMul_base(void (VectorHandler::*op)(integer iRow,
			const doublereal &dCoef),
		VectorHandler& out, const VectorHandler& in) const
{
	if (in.iGetSize() != iGetNumCols()
			|| out.iGetSize() != iGetNumRows())
	{
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
  	}

	row_cont_type::const_iterator ri, re;

	for (integer col = 0; col < NCols; col++) {
		re = col_indices[col].end();
		for (ri = col_indices[col].begin(); ri != re; ++ri) {
			doublereal d = ri->second*in(col + 1);
			(out.*op)(ri->first + 1, d);
		}
	}

	return out;
}

VectorHandler &
SpMapMatrixHandler::MatTVecMul_base(void (VectorHandler::*op)(integer iRow,
			const doublereal &dCoef),
		VectorHandler& out, const VectorHandler& in) const
{
	if (out.iGetSize() != iGetNumCols()
			|| in.iGetSize() != iGetNumRows())
	{
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	row_cont_type::const_iterator ri, re;

	for (integer col = 0; col < NCols; col++) {
		doublereal d = 0.;
		re = col_indices[col].end();
		for (ri = col_indices[col].begin(); ri != re; ++ri) {
			d += ri->second*in(ri->first + 1);
		}
		(out.*op)(col+1, d);
	}

	return out;
}

void
SpMapMatrixHandler::const_iterator::reset(bool is_end)
{
	if (is_end) {
		elem.iRow = m.NRows;
		elem.iCol = m.NCols;

	} else {
		i = m.col_indices[0].begin();
		elem.iRow = i->first;
		elem.iCol = 0;
		elem.dCoef = i->second;
	}
}

SpMapMatrixHandler::const_iterator::const_iterator(const SpMapMatrixHandler& m, bool is_end)
: m(m)
{
	reset(is_end);
}

SpMapMatrixHandler::const_iterator::~const_iterator(void)
{
	NO_OP;
}

const SpMapMatrixHandler::const_iterator&
SpMapMatrixHandler::const_iterator::operator ++ (void) const
{
	++i;
	while (i == m.col_indices[elem.iCol].end()) {
		if (++elem.iCol == m.NCols) {
			elem.iRow = m.NRows;
			return *this;
		}

		i = m.col_indices[elem.iCol].begin();
	}

	elem.iRow = i->first;
	elem.dCoef = i->second;

	return *this;
}

const SparseMatrixHandler::SparseMatrixElement *
SpMapMatrixHandler::const_iterator::operator -> (void) const
{
	return &elem;
}

const SparseMatrixHandler::SparseMatrixElement&
SpMapMatrixHandler::const_iterator::operator * (void) const
{
	return elem;
}

bool
SpMapMatrixHandler::const_iterator::operator == (const SpMapMatrixHandler::const_iterator& op) const
{
	return elem == op.elem;
}

bool
SpMapMatrixHandler::const_iterator::operator != (const SpMapMatrixHandler::const_iterator& op) const
{
	return elem != op.elem;
}

std::ostream& SpMapMatrixHandler::Print(std::ostream& os, MatPrintFormat eFormat) const
{
    if (eFormat == MAT_PRINT_TRIPLET) {
        for (auto i = begin(); i != end(); ++i) {
                if (i->dCoef) { // Output only nonzero entries
                        os << i->iRow + 1 << '\t' << i->iCol + 1 << '\t' << i->dCoef << '\n';
                }
        }
        
        return os;
    } else {
        return SparseMatrixHandler::Print(os, eFormat);
    }
}

void SpMapMatrixHandler::ForceSymmetricGraph()
{
     for (const auto& oItem: *this) {
	  (*this)(oItem.iCol + 1, oItem.iRow + 1);
     }
}
