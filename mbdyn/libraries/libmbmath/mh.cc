/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <string.h>	/* for memset() */

#include "ac/iostream"
#include "ac/iomanip"

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
	for (integer i = 1; i <= iGetNumRows(); i++) {
		for (integer j = 0; j <= iGetNumCols(); j++) {
			this->operator()(i, j) *= d;
		}
	}

	return *this;
}

VectorHandler&
MatrixHandler::MatVecMul(VectorHandler& out, const VectorHandler& in) const
{
	if (out.iGetSize() != iGetNumRows()
			|| in.iGetSize() != iGetNumCols()) {
		THROW(ErrGeneric());
	}

	for (integer r = 1; r <= iGetNumRows(); r++) {
		doublereal d = 0.;

		for (integer c = 1; c <= in.iGetSize(); c++) {
			d += this->operator()(r, c)*in(c);
		}
		out(r) = d;
	}

	return out;
}

VectorHandler&
MatrixHandler::MatTVecMul(VectorHandler& out, const VectorHandler& in) const
{
	if (out.iGetSize() != iGetNumCols()
			|| in.iGetSize() != iGetNumRows()) {
		THROW(ErrGeneric());
	}

	for (integer r = 1; r <= iGetNumCols(); r++) {
		doublereal d = 0.;

		for (integer c = 1; c <= in.iGetSize(); c++) {
			d += this->operator()(c, r)*in(c);
		}
		out(r) = d;
	}
	return out;
}

VectorHandler&
MatrixHandler::MatVecIncMul(VectorHandler& out, const VectorHandler& in) const
{
	if (out.iGetSize() != iGetNumRows()
			|| in.iGetSize() != iGetNumCols()) {
		THROW(ErrGeneric());
	}

	for (integer r = 1; r <= iGetNumRows(); r++) {
		doublereal d = 0.;

		for (integer c = 1; c <= in.iGetSize(); c++) {
			d += this->operator()(r, c)*in(c);
		}
		out(r) += d;
	}
	return out;
}

VectorHandler&
MatrixHandler::MatTVecIncMul(VectorHandler& out, const VectorHandler& in) const
{
	if (out.iGetSize() != iGetNumCols()
			|| in.iGetSize() != iGetNumRows()) {
		THROW(ErrGeneric());
	}

	for (integer r = 1; r <= iGetNumCols(); r++) {
		doublereal d = 0.;

		for (integer c = 1; c <= in.iGetSize(); c++) {
			d += this->operator()(c, r)*in(c);
		}
		out(r) += d;
	}
	return out;
}

void 
MatrixHandler::IncCoef(integer ix, integer iy, const doublereal& inc) {
	operator()(ix,iy) += inc;
};

void 
MatrixHandler::DecCoef(integer ix, integer iy, const doublereal& inc) {
	operator()(ix,iy) -= inc;
};

void 
MatrixHandler::PutCoef(integer ix, integer iy, const doublereal& val) {
	operator()(ix,iy) = val;
};

const doublereal& 
MatrixHandler::dGetCoef(integer ix, integer iy) const {
	return operator()(ix, iy);
};


std::ostream&
operator << (std::ostream& out, const MatrixHandler& MH)
{
	for (integer i = 1; i <= MH.iGetNumRows(); i++) {
		for (integer j = 1; j <= MH.iGetNumCols(); j++) {
			out << std::setw(16) << MH(i, j);
		}
		out << std::endl;
	}
	return out;
}

/* MatrixHandler - end */

