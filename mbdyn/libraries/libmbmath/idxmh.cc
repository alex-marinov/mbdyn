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
#include "idxmh.h"

IndexMatrixHandler::IndexMatrixHandler(const int &n,
		std::vector<doublereal>& x,
		const std::vector<int>& i,
		const std::vector<int>& p)
: CompactSparseMatrixHandler(n, n, x, i, p)
{
	NO_OP;
}

IndexMatrixHandler::~IndexMatrixHandler()
{
	NO_OP;
}

/* used by MultiThreadDataManager to duplicate the storage array
 * while preserving the CC indices */
CompactSparseMatrixHandler *
IndexMatrixHandler::Copy(void) const
{
	std::vector<doublereal> *pax = new std::vector<doublereal>(Ax);
	IndexMatrixHandler *p = new IndexMatrixHandler(NRows, *pax, Ai, Ap);
	p->bMatDuplicate = true;

	return p;
}

int
IndexMatrixHandler::MakeCompressedColumnForm( doublereal *const Ax,
		int *const Ai,
		int *const Ap,
		integer offset) const
{
	silent_cerr("IndexMatrixHandler::MakeCompressedColumnForm called"
			<< std::endl);
	THROW(ErrGeneric());		
	return 0;
}

int
IndexMatrixHandler::MakeCompressedColumnForm( std::vector<doublereal>& Ax,
                std::vector<int>& Ai,
                std::vector<int>& Ap,
		integer offset) const
{
	silent_cerr("IndexMatrixHandler::MakeCompressedColumnForm called"
			<< std::endl);
	THROW(ErrGeneric());		
	return 0;
}

int
IndexMatrixHandler::MakeIndexForm( doublereal *const rAx,
		integer *const Arow,
		integer *const Acol,
		integer offset) const
{
	return Nz();
}

int
IndexMatrixHandler::MakeIndexForm( std::vector<doublereal>& rAx,
                std::vector<integer>& Arow,
                std::vector<integer>& Acol,
		integer offset) const
{
	return MakeIndexForm(&(rAx[0]), &(Arow[0]), &(Acol[0]), offset);
}

void
IndexMatrixHandler::Resize(const int &n, const int &nn)
{
	silent_cerr("IndexMatrixHandler::Resize called" << std::endl);
	THROW(ErrGeneric());
}

/* Estrae una colonna da una matrice */
VectorHandler&
IndexMatrixHandler::GetCol(integer icol, VectorHandler& out) const
{
	THROW(ErrGeneric());		
	return out;
}
	
/* Prodotto Matrice per Matrice */
SpMapMatrixHandler& 
IndexMatrixHandler::MatMatMul(SpMapMatrixHandler& out,
		const SpMapMatrixHandler& in) const
{
	THROW(ErrGeneric());		
	return out;	
}
	
/* Moltiplica per uno scalare e somma a una matrice */
MatrixHandler&
IndexMatrixHandler::MulAndSumWithShift(MatrixHandler& out, doublereal s,
		integer drow, integer dcol) const
{
	THROW(ErrGeneric());		
	return out;	
}

MatrixHandler&
IndexMatrixHandler::FakeThirdOrderMulAndSumWithShift(MatrixHandler& out, 
		std::vector<bool> b,
		doublereal s,
		integer drow, 
		integer dcol) const
{
	THROW(ErrGeneric());		
	return out;	
}
	
VectorHandler&
IndexMatrixHandler::MatTVecMul(VectorHandler& out,
		const VectorHandler& in) const
{
	THROW(ErrGeneric());		
	return out;
}
	
VectorHandler&
IndexMatrixHandler::MatVecMul(VectorHandler& out,
		const VectorHandler& in) const
{
	THROW(ErrGeneric());		
	return MatVecIncMul(out, in);
}

VectorHandler&
IndexMatrixHandler::MatTVecIncMul(VectorHandler& out,
		const VectorHandler& in) const
{
	THROW(ErrGeneric());		
	return out;
}
	
VectorHandler&
IndexMatrixHandler::MatVecIncMul(VectorHandler& out,
		const VectorHandler& in) const {
	THROW(ErrGeneric());		
	return out;
}

VectorHandler&
IndexMatrixHandler::MatVecDecMul(VectorHandler& out,
		const VectorHandler& in) const
{
	THROW(ErrGeneric());		
	return out;
}

