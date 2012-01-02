/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
 *
 * Pierangelo Masarati	<masarati@aero.polimi.it>
 * Paolo Mantegazza	<mantegazza@aero.polimi.it>
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

/* solution manager */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cstring>	/* for memset() */

#include <iostream>
#include <iomanip>

#include "solman.h"
#include "submat.h"
#include "matvec3.h"

/* MatrixHandler - begin */

MatrixHandler::~MatrixHandler(void)
{
	NO_OP;
}

/* Ridimensiona ed inizializza. */
void
MatrixHandler::ResizeReset(integer iNewRow, integer iNewCol)
{
	Resize(iNewRow, iNewCol);
	Reset();
}

/* Impacchetta la matrice; restituisce il numero di elementi 
 * diversi da zero */
integer
MatrixHandler::PacMat(void)
{
	return 0L;
}

/* Overload di = */
MatrixHandler& 
MatrixHandler::operator = (const MatrixHandler& MH)
{
	integer nr = MH.iGetNumRows();
	integer nc = MH.iGetNumCols();

	Resize(nr, nc);

	for (integer i = 1; i <= nr; i++) { 
		for (integer ii = 1; ii <= nc; ii++) { 
			this->operator()(i, ii) = MH(i, ii);
		}
	}

	return *this;
}



/* Overload di += usato per l'assemblaggio delle matrici */
MatrixHandler&
MatrixHandler::operator +=(const SubMatrixHandler& SubMH)
{
	return SubMH.AddTo(*this);
}

/* Overload di -= usato per l'assemblaggio delle matrici */
MatrixHandler&
MatrixHandler::operator -=(const SubMatrixHandler& SubMH)
{
	return SubMH.SubFrom(*this);
}

/* Overload di += usato per l'assemblaggio delle matrici */
MatrixHandler&
MatrixHandler::operator +=(const VariableSubMatrixHandler& SubMH)
{
	return SubMH.AddTo(*this);
}

/* Overload di -= usato per l'assemblaggio delle matrici */
MatrixHandler&
MatrixHandler::operator -=(const VariableSubMatrixHandler& SubMH)
{
	return SubMH.SubFrom(*this);
}

MatrixHandler&
MatrixHandler::ScalarMul(const doublereal& d)
{
	integer nr = iGetNumRows();
	integer nc = iGetNumCols();

	for (integer i = 1; i <= nr; i++) {
		for (integer j = 1; j <= nc; j++) {
			this->operator()(i, j) *= d;
		}
	}

	return *this;
}

/* Matrix Matrix product */
MatrixHandler&
MatrixHandler::MatMatMul_base(void (MatrixHandler::*op)(integer iRow, 
			integer iCol, const doublereal& dCoef),
		MatrixHandler& out, const MatrixHandler& in) const
{
	integer out_nc = out.iGetNumCols();
	integer out_nr = out.iGetNumRows();
	integer in_nr = in.iGetNumRows();

	if (out_nr != iGetNumRows()
		|| out_nc != in.iGetNumCols()
		|| in_nr != iGetNumCols())
	{
		const char *strop;

		if (op == &MatrixHandler::IncCoef) {
			strop = "+=";
		} else if (op == &MatrixHandler::DecCoef) {
			strop = "-=";
		} else if (op == &MatrixHandler::PutCoef) {
			strop = "=";
		} else {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		silent_cerr("MatrixHandler::MatMatMul_base: size mismatch "
			"out(" << out_nr << ", " << out_nc << ") "
			<< strop << " this(" << iGetNumRows() << ", " << iGetNumCols() << ") "
			"* in(" << in_nr << ", " << in.iGetNumCols() << ")"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	for (integer c = 1; c <= out_nc; c++) {
		for (integer r = 1; r <= out_nr; r++) {
			doublereal d = 0.;

			for (integer k = 1; k <= in_nr; k++) {
				d += this->operator()(r, k)*in(k, c);
			}

			(out.*op)(r, c, d);
		}
	}

	return out;
}

MatrixHandler&
MatrixHandler::MatTMatMul_base(void (MatrixHandler::*op)(integer iRow, 
			integer iCol, const doublereal& dCoef),
		MatrixHandler& out, const MatrixHandler& in) const
{
	integer out_nc = out.iGetNumCols();
	integer out_nr = out.iGetNumRows();
	integer in_nr = in.iGetNumRows();

	if (out_nr != iGetNumCols()
		|| out_nc != in.iGetNumCols()
		|| in_nr != iGetNumRows())
	{
		const char *strop;

		if (op == &MatrixHandler::IncCoef) {
			strop = "+=";
		} else if (op == &MatrixHandler::DecCoef) {
			strop = "-=";
		} else if (op == &MatrixHandler::PutCoef) {
			strop = "=";
		} else {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		silent_cerr("MatrixHandler::MatTMatMul_base: size mismatch "
			"out(" << out_nr << ", " << out_nc << ") "
			<< strop << " this(" << iGetNumRows() << ", " << iGetNumCols() << ")^T "
			"* in(" << in_nr << ", " << in.iGetNumCols() << ")"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	for (integer c = 1; c <= out_nc; c++) {
		for (integer r = 1; r <= out_nr; r++) {
			doublereal d = 0.;

			for (integer k = 1; k <= in_nr; k++) {
				d += this->operator()(k, r)*in(k, c);
			}

			(out.*op)(r, c, d);
		}
	}

	return out;
}

MatrixHandler&
MatrixHandler::MatMatMul(MatrixHandler& out, const MatrixHandler& in) const
{
	/* Put is implemented resetting out first, then passing IncCoef()
	 * so that out-of-order assignments work */
	out.Reset();
	return MatMatMul_base(&MatrixHandler::IncCoef, out, in);
}

MatrixHandler&
MatrixHandler::MatTMatMul(MatrixHandler& out, const MatrixHandler& in) const
{
	/* Put is implemented resetting out first, then passing IncCoef()
	 * so that out-of-order assignments work */
	out.Reset();
	return MatTMatMul_base(&MatrixHandler::IncCoef, out, in);
}

MatrixHandler&
MatrixHandler::MatMatIncMul(MatrixHandler& out, const MatrixHandler& in) const
{
	return MatMatMul_base(&MatrixHandler::IncCoef, out, in);
}

MatrixHandler&
MatrixHandler::MatTMatIncMul(MatrixHandler& out, const MatrixHandler& in) const
{
	return MatTMatMul_base(&MatrixHandler::IncCoef, out, in);
}

MatrixHandler&
MatrixHandler::MatMatDecMul(MatrixHandler& out, const MatrixHandler& in) const
{
	return MatMatMul_base(&MatrixHandler::DecCoef, out, in);
}

MatrixHandler&
MatrixHandler::MatTMatDecMul(MatrixHandler& out, const MatrixHandler& in) const
{
	return MatTMatMul_base(&MatrixHandler::DecCoef, out, in);
}

/* Matrix Vector product */
VectorHandler&
MatrixHandler::MatVecMul_base(
	void (VectorHandler::*op)(integer iRow, const doublereal& dCoef),
	VectorHandler& out, const VectorHandler& in) const
{
	integer nr = iGetNumRows();
	integer nc = iGetNumCols();

	if (out.iGetSize() != nr || in.iGetSize() != nc) {
		silent_cerr("MatrixHandler::MatVecMul_base(): size mismatch "
			"out(" << out.iGetSize() << ", 1) "
			"= this(" << nr << ", " << nc << ") "
			"* in(" << in.iGetSize() << ", 1)" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	for (integer r = 1; r <= nr; r++) {
		doublereal d = 0.;

		for (integer c = 1; c <= nc; c++) {
			d += this->operator()(r, c)*in(c);
		}
		(out.*op)(r, d);
	}

	return out;

}

VectorHandler&
MatrixHandler::MatTVecMul_base(
	void (VectorHandler::*op)(integer iRow, const doublereal& dCoef),
	VectorHandler& out, const VectorHandler& in) const
{
	integer nr = iGetNumRows();
	integer nc = iGetNumCols();

	if (out.iGetSize() != nc || in.iGetSize() != nr) {
		silent_cerr("MatrixHandler::MatVecMul_base(): size mismatch "
			"out(" << out.iGetSize() << ", 1) "
			"= this(" << nr << ", " << nc << ")^T "
			"* in(" << in.iGetSize() << ", 1)" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	for (integer r = 1; r <= nc; r++) {
		doublereal d = 0.;

		for (integer c = 1; c <= nr; c++) {
			d += this->operator()(c, r)*in(c);
		}
		(out.*op)(r, d);
	}

	return out;
}

VectorHandler&
MatrixHandler::MatVecMul(VectorHandler& out, const VectorHandler& in) const
{
	/* Put is implemented resetting out first, then passing IncCoef()
	 * so that out-of-order assignments work */
	out.Reset();
	return MatVecMul_base(&VectorHandler::IncCoef, out, in);
}

VectorHandler&
MatrixHandler::MatTVecMul(VectorHandler& out, const VectorHandler& in) const
{
	/* Put is implemented resetting out first, then passing IncCoef()
	 * so that out-of-order assignments work */
	out.Reset();
	return MatTVecMul_base(&VectorHandler::IncCoef, out, in);
}

VectorHandler&
MatrixHandler::MatVecIncMul(VectorHandler& out, const VectorHandler& in) const
{
	return MatVecMul_base(&VectorHandler::IncCoef, out, in);
}

VectorHandler&
MatrixHandler::MatTVecIncMul(VectorHandler& out, const VectorHandler& in) const
{
	return MatTVecMul_base(&VectorHandler::IncCoef, out, in);
}

VectorHandler&
MatrixHandler::MatVecDecMul(VectorHandler& out, const VectorHandler& in) const
{
	return MatVecMul_base(&VectorHandler::DecCoef, out, in);
}

VectorHandler&
MatrixHandler::MatTVecDecMul(VectorHandler& out, const VectorHandler& in) const
{
	return MatTVecMul_base(&VectorHandler::DecCoef, out, in);
}

void 
MatrixHandler::IncCoef(integer ix, integer iy, const doublereal& inc) {
	operator()(ix, iy) += inc;
}

void 
MatrixHandler::DecCoef(integer ix, integer iy, const doublereal& inc) {
	operator()(ix, iy) -= inc;
}

void 
MatrixHandler::PutCoef(integer ix, integer iy, const doublereal& val) {
	operator()(ix, iy) = val;
}

const doublereal& 
MatrixHandler::dGetCoef(integer ix, integer iy) const {
	return operator()(ix, iy);
}

std::ostream&
operator << (std::ostream& out, const MatrixHandler& MH)
{
	integer nr = MH.iGetNumRows();
	integer nc = MH.iGetNumCols();

	for (integer i = 1; i <= nr; i++) {
		for (integer j = 1; j <= nc; j++) {
			out << std::setw(16) << MH(i, j);
		}
		out << std::endl;
	}

	return out;
}

/* MatrixHandler - end */

