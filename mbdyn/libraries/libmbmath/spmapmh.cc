/* 
 * HmFe (C) is a FEM analysis code. 
 *
 * Copyright (C) 1996-2004
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
 * Copyright (C) 2003-2004
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
 * Copyright (C) 2003-2004
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "spmapmh.h"

SpMapMatrixHandler::SpMapMatrixHandler(const integer &n, const integer &nn)
: SparseMatrixHandler(n, nn)
{
	col_indices.resize(NCols);
}

SpMapMatrixHandler::~SpMapMatrixHandler()
{
	NO_OP;
}

integer
SpMapMatrixHandler::MakeCompressedColumnForm(doublereal *const Ax,
		integer *const Ai, integer *const Ap, int offset) const
{
	integer x_ptr = 0;
		
	row_cont_type::iterator ri;
	row_cont_type::const_iterator re;
		
	for (integer col = 0; col < NCols; col++) {
		Ap[col] = x_ptr;
		re = col_indices[col].end();
		for (ri = col_indices[col].begin(); ri != re; ri++) {
			Ax[x_ptr] = ri->second;
			Ai[x_ptr] = ri->first + offset;
			x_ptr++;
		}
	}

	ASSERTMSGBREAK(x_ptr == NZ, "Error in "
			"SpMapMatrixHandler::MakeCompressedColumnForm");

	Ap[NCols] = x_ptr;

	return Nz();
}

integer
SpMapMatrixHandler::MakeCompressedColumnForm(std::vector<doublereal>& Ax,
                std::vector<integer>& Ai, std::vector<integer>& Ap,
		int offset) const
{
	Ax.resize(Nz());
	Ai.resize(Nz());
	Ap.resize(iGetNumCols() + 1);

	return MakeCompressedColumnForm(&Ax[0], &Ai[0], &Ap[0], offset);
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
		Ap[col] = x_ptr;
		re = col_indices[col].end();
		for (ri = col_indices[col].begin(); ri != re; ri++) {
			Ax[x_ptr] = ri->second;
			Arow[x_ptr] = ri->first + offset;
			Acol[x_ptr] = col + offset;
			x_ptr++;
		}
	}

	ASSERTMSGBREAK(x_ptr == NZ, "Error in "
			"SpMapMatrixHandler::MakeIndexForm");

	Ap[NCols] = x_ptr;

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
		for (ri = col_indices[col].begin();ri != re; ri++) {
			ri->second = 0.;
		}
	}
}

void
SpMapMatrixHandler::Resize(integer n, integer nn)
{
	integer nnn;

	if (nn == 0) {
		nnn = n;
	} else {
		nnn = nn;
	}

	for (integer col = 0; col < NCols; col++) {
		col_indices[col].clear();
	}

	col_indices.resize(nnn);
	NRows = n;
	NCols = nnn;
	NZ = 0;
}

/* Estrae una colonna da una matrice */
VectorHandler&
SpMapMatrixHandler::GetCol(integer icol, VectorHandler& out) const
{
        if (icol > iGetNumCols()) {
		throw ErrGeneric();
	}
	row_cont_type::const_iterator ri, re;
	re = col_indices[icol].end();
	for (ri = col_indices[icol].begin();ri != re; ri++) {		
		out.PutCoef(ri->first + 1, ri->second);
	}
	return out;
}
	
/* Prodotto Matrice per Matrice */
SpMapMatrixHandler &
SpMapMatrixHandler::MatMatMul(SpMapMatrixHandler& out,
		const SpMapMatrixHandler& in) const
{
	if ((in.iGetNumCols() != iGetNumRows())
			|| (in.iGetNumRows() != out.iGetNumRows())
			|| (out.iGetNumCols() != iGetNumCols())) {
		silent_cerr("Assertion fault "
			"in SpMapMatrixHandler::MatMatMul" << std::endl);
		throw ErrGeneric();
	}

	out.Reset();

	for (integer col = 0; col < NCols; col++) {
		row_cont_type::const_iterator ri, re;
		re = col_indices[col].end();
		for (ri = col_indices[col].begin(); ri != re; ri++) {
			integer iend = in.iGetNumCols();
			for (integer col2 = 0; col2 < iend;  col2++) {
				out.IncCoef(ri->first,col2,
						ri->second*in(col, col2));
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
	if ((out.iGetNumCols() < iGetNumCols()+dcol)
			|| (out.iGetNumRows() < iGetNumRows()+drow)) {
		silent_cerr("Assertion fault "
			"in SpMapMatrixHandler::MulAndSumWithShift"
			<< std::endl);
		throw ErrGeneric();
	}

	drow = drow + 1;

	for (integer col = 0; col < NCols; col++) {
		row_cont_type::const_iterator ri, re;
		re = col_indices[col].end();
		integer newcol = col + dcol + 1;
		for (ri = col_indices[col].begin(); ri != re; ri++) {
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
	if ((out.iGetNumCols() < iGetNumCols()+dcol)
			|| (out.iGetNumRows() < iGetNumRows()+drow)) {
		silent_cerr("Assertion fault "
			"in SpMapMatrixHandler::MulAndSumWithShift"
			<< std::endl);
		throw ErrGeneric();
	}

	drow = drow + 1;

	for (integer col = 0; col < NCols; col++) {
		row_cont_type::const_iterator ri, re;
		re = col_indices[col].end();
		integer newcol = col + dcol + 1;
		for (ri = col_indices[col].begin(); ri != re; ri++) {
			if (b[ri->first]) {
				out.IncCoef(ri->first + drow,
						newcol, ri->second*s);
			}
		}
	}

	return out;	
}

VectorHandler &
SpMapMatrixHandler::MatTVecMul(VectorHandler& out,
		const VectorHandler& in) const
{
	if (out.iGetSize() != iGetNumRows()
			|| in.iGetSize() != iGetNumCols()) {
		throw ErrGeneric();
	}

	row_cont_type::const_iterator ri, re;

	for (integer col = 0; col < NCols; col++) {
		doublereal d = 0.;
		re = col_indices[col].end();
		for (ri = col_indices[col].begin();ri != re; ri++) {
			d += ri->second*in(ri->first + 1);
		}
		out.PutCoef(col+1, d);
	}

	return out;
}
	
VectorHandler &
SpMapMatrixHandler::MatVecMul(VectorHandler& out,
		const VectorHandler& in) const
{
	if (in.iGetSize() != iGetNumCols()
			|| out.iGetSize() != iGetNumRows()) {
		throw ErrGeneric();
  	}

	row_cont_type::const_iterator ri, re;
	out.Reset();

	for (integer col = 0; col < NCols; col++) {
		re = col_indices[col].end();
		for (ri = col_indices[col].begin();ri != re; ri++) {
			doublereal d = ri->second*in(col + 1);
			out.IncCoef(ri->first + 1, d);
		}
	}

	return out;
}

VectorHandler &
SpMapMatrixHandler::MatTVecIncMul(VectorHandler& out,
		const VectorHandler& in) const
{
	if (out.iGetSize() != iGetNumRows()
			|| in.iGetSize() != iGetNumCols()) {
		throw ErrGeneric();
	}

	row_cont_type::const_iterator ri, re;

	for (integer col = 0; col < NCols; col++) {
		doublereal d = 0.;
		re = col_indices[col].end();
		for (ri = col_indices[col].begin();ri != re; ri++) {
			d += ri->second*in(ri->first + 1);
		}
		out.IncCoef(col+1, d);
	}

	return out;
}
	
VectorHandler &
SpMapMatrixHandler::MatVecIncMul(VectorHandler& out,
		const VectorHandler& in) const
{
	if (in.iGetSize() != iGetNumCols()
			|| out.iGetSize() != iGetNumRows()) {
		throw ErrGeneric();
	}

	row_cont_type::const_iterator ri, re;

	for (integer col = 0; col < NCols; col++) {
		re = col_indices[col].end();
		for (ri = col_indices[col].begin();ri != re; ri++) {
			doublereal d = ri->second*in(col + 1);
			out.IncCoef(ri->first + 1, d);
		}
	}

	return out;
}

VectorHandler &
SpMapMatrixHandler::MatVecDecMul(VectorHandler& out,
		const VectorHandler& in) const
{
	if (in.iGetSize() != iGetNumCols()
			|| out.iGetSize() != iGetNumRows()) {
		throw ErrGeneric();
	}

	row_cont_type::const_iterator ri, re;

	for (integer col = 0; col < NCols; col++) {
		re = col_indices[col].end();
		for (ri = col_indices[col].begin();ri != re; ri++) {
			doublereal d = ri->second*in(col + 1);
			out.DecCoef(ri->first + 1, d);
		}
	}

	return out;
}

