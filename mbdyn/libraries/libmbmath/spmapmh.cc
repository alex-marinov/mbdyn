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

#if 0	/* column-oriented */
integer
SpMapMatrixHandler::MakeNaiveForm(doublereal *const Ax,
		integer *const Arow, integer *const Acol,
		integer *const curCol, int offset) const
{
	row_cont_type::iterator ri;
	row_cont_type::const_iterator re;

	integer *picol = Acol;
	integer size = iGetNumCols();
	doublereal *pd = Ax;

	for (integer col = 0; col < NCols; col++) {
		integer row_ptr = 0;

		re = col_indices[col].end();
		for (ri = col_indices[col].begin(); ri != re; ri++) {
			/* the location in column <col>, row <ri->first>
			 * is set to <ri->second>, the matrix coefficient */
			pd[ri->first + offset] = ri->second;

			/* the location in column <col>, row <row_ptr>
			 * is set to <ri->first>, the row index */
			picol[row_ptr++] = ri->first + offset;

			/* the location in column <col_ptr>, row <ri->first>
			 * is set to <col>, the column index
			 *
			 * FIXME: this is the weak part; I'm using
			 * a linear search throughout the row, maybe
			 * a binary search between columns 0 and <col>
			 * is better.
			 */
#if 1
			Arow[ri->first + curCol[ri->first]*size] = col;
			curCol[ri->first]++;
#else
			integer *pdrow = &Arow[ri->first];
			for (integer ic = 0; ic < size; ic++) {
				if (pdrow[0] == -offset) {
					break;
				}
				pdrow += size;
			}
			ASSERT(pdrow[0] == -offset);
			pdrow[0] = col;
#endif
		}

		picol += size;
		pd += size;
	}

	return Nz();
}
#endif /* 0 */

integer
SpMapMatrixHandler::MakeNaiveForm(doublereal *const Ax,
		integer *const Arow, integer *const Acol,
		integer *const nzr, integer *const nzc,
		int offset) const
{
	integer	size = iGetNumRows();
	
	for (integer col = 0; col < NCols; col++) {
		row_cont_type::iterator row_it = col_indices[col].begin();
		row_cont_type::const_iterator row_end = col_indices[col].end();

		for (; row_it != row_end; row_it++) {
			integer irow = row_it->first;
			integer nirow = size*irow;
// 			integer icol = col;
			
			if (row_it->second) {
				Ax[nirow + col] = row_it->second;

				Arow[nirow + col] |= 0xF0000000;
				Arow[size*col + nzr[col]] |= irow;
				Acol[nirow + nzc[irow]] = col;
// 				std::cerr << "inserisco in " <<
// 					nirow + nzc[irow] <<
// 					" " << col << std::endl;

				nzc[irow]++;
				nzr[col]++;
			}
		}
	}

#if 0
	std::cerr << "/////////////////////////////////////" << std::endl;
	std::cerr << "int neq = " << size << ";" << std::endl;
	std::cerr << "int nzr[] = {" << std::endl;
	for (integer ir = 0; ir < size - 1; ir++) {
		std::cerr << "\t" << nzr[ir] << "," << std::endl;
	}
	std::cerr << "\t" << nzr[size - 1] << std::endl
		<< "};" << std::endl;

	std::cerr << "int nzc[] = {" << std::endl;
	for (integer ir = 0; ir < size - 1; ir++) {
		std::cerr << "\t" << nzc[ir] << "," << std::endl;
	}
	std::cerr << "\t" << nzc[size - 1] << std::endl
		<< "};" << std::endl;

	std::cerr << "int ri[]["<< size << "] = {" << std::endl;
	for (integer ir = 0; ir < size; ir++) {
		std::cerr << "\t{";
		for (integer iz = 0; iz < nzr[ir] - 1; iz++) {
			std::cerr << Arow[size*ir + iz] << ", ";
		}
		std::cerr << Arow[size*ir + nzr[ir] - 1] << "}," << std::endl;
	}
	std::cerr << "\t{0} // last + 1" << std::endl
		<< "};" << std::endl;

	std::cerr << "int ci[]["<< size << "] = {" << std::endl;
	for (integer ic = 0; ic < size; ic++) {
		std::cerr << "\t{";
		for (integer iz = 0; iz < nzc[ic] - 1; iz++) {
			std::cerr << Acol[size*ic + iz] << ", ";
		}
		std::cerr << Arow[size*ic + nzc[ic] - 1] << "}," << std::endl;
	}
	std::cerr << "\t{0} // last + 1" << std::endl
		<< "};" << std::endl;
	std::cerr << "double A[][" << size << "] = {" << std::endl;
	for (integer ir = 0; ir < size; ir++) {
		std::cerr << "\t{";
		for (integer ic = 0; ic < size-1; ic++) {
			std::cerr << Ax[size*ir + ic] << ", ";
		}
		std::cerr << Ax[size*ir + size - 1] << "}," << std::endl;
	}
	std::cerr << "\t{0} // last + 1" << std::endl
		<< "};" << std::endl;
#endif

	return Nz();
}

integer
SpMapMatrixHandler::MakeNaiveForm(std::vector<doublereal>& Ax,
                std::vector<integer>& Arow, std::vector<integer>& Acol,
		std::vector<integer>& Nzr, std::vector<integer>& Nzc,
		int offset) const
{
	integer s = iGetNumCols();
	integer s2 = s*s;

	Ax.resize(s2);
	//std::fill(Ax.begin(), Ax.end(), 0.);
	Arow.resize(s2);
	std::fill(Arow.begin(), Arow.end(), 0);
	Acol.resize(s2);

	Nzr.resize(s, 0);
	std::fill(Nzr.begin(), Nzr.end(), 0);
	Nzc.resize(s, 0);
	std::fill(Nzc.begin(), Nzc.end(), 0);
	

	return MakeNaiveForm(&Ax[0], &Arow[0], &Acol[0],
			&Nzr[0], &Nzc[0], offset);
}

void
SpMapMatrixHandler::Reset(const doublereal &r)
{
	row_cont_type::const_iterator re;
	row_cont_type::iterator ri;
	for (integer col = 0; col < NCols; col++) {
		re = col_indices[col].end();
		for (ri = col_indices[col].begin();ri != re; ri++) {
			ri->second = r;
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
		THROW(ErrGeneric());
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
		std::cerr << "Assertion fault "
			"in SpMapMatrixHandler::MatMatMul" << std::endl;
		THROW(ErrGeneric());
	}

	out.Reset(0.);

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
		std::cerr << "Assertion fault "
			"in SpMapMatrixHandler::MulAndSumWithShift"
			<< std::endl;
		THROW(ErrGeneric());
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
		std::cerr << "Assertion fault "
			"in SpMapMatrixHandler::MulAndSumWithShift"
			<< std::endl;
		THROW(ErrGeneric());
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
		THROW(ErrGeneric());
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
		THROW(ErrGeneric());
  	}

	row_cont_type::const_iterator ri, re;
	out.Reset(0.);

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
		THROW(ErrGeneric());
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
		THROW(ErrGeneric());
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
		THROW(ErrGeneric());
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

