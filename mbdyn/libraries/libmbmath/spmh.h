/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2003-2008
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

#ifndef SPMH_H
#define SPMH_H

#include <vector>

#include "myassert.h"
#include "solman.h"

/* Sparse Matrix */
class SparseMatrixHandler : public MatrixHandler {
protected:
	integer NRows;
	integer NCols;
	integer NZ;

#ifdef DEBUG
	void IsValid(void) const {
		NO_OP;
	};
#endif /* DEBUG */

public:
	struct SparseMatrixElement {
		integer iRow;
		integer iCol;
		doublereal dCoef;

		SparseMatrixElement(void)
		: iRow(0), iCol(0), dCoef(0.) {};
		SparseMatrixElement(integer iRow, integer iCol, const doublereal& dCoef)
		: iRow(iRow), iCol(iCol), dCoef(dCoef) {};
	};

public:
	const integer Nz() const {
		return NZ;
	};
	
	/* FIXME: always square? */
	SparseMatrixHandler(const integer &n, const integer &nn = 0);

	virtual ~SparseMatrixHandler(void);

	integer iGetNumRows(void) const {
		return NRows;
	};

	integer iGetNumCols(void) const {
		return NCols;
	};

	virtual
	integer MakeCompressedColumnForm(doublereal *const Ax,
			integer *const Ai, integer *const Ap,
			int offset = 0) const = 0;

	virtual
        integer MakeCompressedColumnForm(std::vector<doublereal>& Ax,
                	std::vector<integer>& Ai, std::vector<integer>& Ap,
			int offset = 0) const = 0;

	virtual
	integer MakeIndexForm(doublereal *const Ax,
			integer *const Arow, integer *const Acol,
			integer *const AcolSt,
			int offset = 0) const = 0;

	virtual
        integer MakeIndexForm(std::vector<doublereal>& Ax,
			std::vector<integer>& Arow, std::vector<integer>& Acol,
			std::vector<integer>& AcolSt,
			int offset = 0) const = 0;

	/* Estrae una colonna da una matrice */
	virtual VectorHandler& GetCol(integer icol,
			VectorHandler& out) const = 0;
};

/* Sparse Matrix in compact form */
class CompactSparseMatrixHandler : public SparseMatrixHandler {
protected:
	bool bMatDuplicate;
	std::vector<doublereal>& Ax;
	const std::vector<integer>& Ai;
	const std::vector<integer>& Ap;

#ifdef DEBUG
	void IsValid(void) const {
		NO_OP;
	};
#endif /* DEBUG */

public:
	CompactSparseMatrixHandler(const integer &n, const integer &nn,
			std::vector<doublereal>& x,
			const std::vector<integer>& i,
			const std::vector<integer>& p);

	virtual ~CompactSparseMatrixHandler();

	/* used by MultiThreadDataManager to duplicate the storage array
	 * while preserving the CC indices */
	virtual CompactSparseMatrixHandler *Copy(void) const = 0;

	/* used to sum CC matrices with identical indices */
	void AddUnchecked(const CompactSparseMatrixHandler& m);

	/* Restituisce un puntatore all'array di reali della matrice */
	virtual inline doublereal* pdGetMat(void) const {
		return &Ax[0];
	};

public:
	void Reset(void);
	
	virtual
	integer MakeCompressedColumnForm(doublereal *const Ax,
			integer *const Ai, integer *const Ap,
			int offset = 0) const;

	virtual
        integer MakeCompressedColumnForm(std::vector<doublereal>& Ax,
                	std::vector<integer>& Ai, std::vector<integer>& Ap,
			int offset = 0) const;

	virtual
	integer MakeIndexForm(doublereal *const Ax,
			integer *const Arow, integer *const Acol,
			integer *const AcolSt,
			int offset = 0) const;

	virtual
        integer MakeIndexForm(std::vector<doublereal>& Ax,
			std::vector<integer>& Arow, std::vector<integer>& Acol,
			std::vector<integer>& AcolSt,
			int offset = 0) const;
};

/* Sparse Matrix in compact form */
template <int off>
class CompactSparseMatrixHandler_tpl : public CompactSparseMatrixHandler {
public:
	CompactSparseMatrixHandler_tpl(const integer &n, const integer &nn,
			std::vector<doublereal>& x,
			const std::vector<integer>& i,
			const std::vector<integer>& p);
	virtual ~CompactSparseMatrixHandler_tpl(void);

protected:
	/* Matrix Matrix product */
	MatrixHandler*
	MatMatMul_base(void (MatrixHandler::*op)(integer iRow, integer iCol,
				const doublereal& dCoef),
			MatrixHandler* out, const MatrixHandler& in) const;
	MatrixHandler*
	MatTMatMul_base(void (MatrixHandler::*op)(integer iRow, integer iCol,
				const doublereal& dCoef),
			MatrixHandler* out, const MatrixHandler& in) const;

	/* Matrix Vector product */
	virtual VectorHandler&
	MatVecMul_base(void (VectorHandler::*op)(integer iRow,
				const doublereal& dCoef),
			VectorHandler& out, const VectorHandler& in) const;
	virtual VectorHandler&
	MatTVecMul_base(void (VectorHandler::*op)(integer iRow,
				const doublereal& dCoef),
			VectorHandler& out, const VectorHandler& in) const;
};

#endif /* SPMH_H */
