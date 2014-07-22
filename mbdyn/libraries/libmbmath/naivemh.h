/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2003-2014
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

#ifndef NAIVEMH_H
#define NAIVEMH_H

#include <vector>

#include "myassert.h"
#include "solman.h"
#include "spmh.h"

class NaiveSolver;
class MultiThreadDataManager;

/* Sparse Matrix */
class NaiveMatrixHandler : public MatrixHandler {
	friend void* sum_naive_matrices(void* arg);
	friend class NaiveSolver;
	friend class ParNaiveSolver;
	friend class MultiThreadDataManager;

private:
	// don't allow copy constructor!
	NaiveMatrixHandler(const NaiveMatrixHandler&);

protected:
	integer iSize;
	bool bOwnsMemory;
	doublereal **ppdRows;
	integer **ppiRows, **ppiCols;
	char **ppnonzero;
	integer *piNzr, *piNzc;

#ifdef DEBUG
	void IsValid(void) const {
		NO_OP;
	};
#endif /* DEBUG */

public:
	class const_iterator {
		friend class NaiveMatrixHandler;

	private:
		const NaiveMatrixHandler& m;
		mutable integer i_row;
		mutable SparseMatrixHandler::SparseMatrixElement elem;

	protected:
		void reset(bool is_end = false);

	public:
		const_iterator(const NaiveMatrixHandler& m, bool is_end = false);
		~const_iterator(void);
		const NaiveMatrixHandler::const_iterator& operator ++ (void) const;
		const SparseMatrixHandler::SparseMatrixElement* operator -> (void) const;
		const SparseMatrixHandler::SparseMatrixElement& operator * (void) const;
		bool operator == (const NaiveMatrixHandler::const_iterator& op) const;
		bool operator != (const NaiveMatrixHandler::const_iterator& op) const;
	};

protected:
	const_iterator m_end;

public:
	NaiveMatrixHandler::const_iterator begin(void) const {
		return const_iterator(*this);
	};

	const NaiveMatrixHandler::const_iterator& end(void) const {
		return m_end;
	};

public:
	/* FIXME: always square? yes! */
	NaiveMatrixHandler(const integer n, NaiveMatrixHandler *const nmh = 0);

	virtual ~NaiveMatrixHandler(void);

	integer iGetNumRows(void) const {
		return iSize;
	};

	integer iGetNumCols(void) const {
		return iSize;
	};

	void Reset(void);

	/* Ridimensiona la matrice */
	virtual void Resize(integer, integer) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	};

	virtual inline const doublereal&
	operator () (integer iRow, integer iCol) const;

	virtual inline doublereal&
	operator () (integer iRow, integer iCol);

	/* Overload di += usato per l'assemblaggio delle matrici */
	virtual MatrixHandler& operator += (const SubMatrixHandler& SubMH);

	/* Overload di -= usato per l'assemblaggio delle matrici */
	virtual MatrixHandler& operator -= (const SubMatrixHandler& SubMH);

	/* Overload di += usato per l'assemblaggio delle matrici
	 * questi li vuole ma non so bene perche'; force per la doppia
	 * derivazione di VariableSubMatrixHandler */
	virtual MatrixHandler&
	operator += (const VariableSubMatrixHandler& SubMH);
	virtual MatrixHandler&
	operator -= (const VariableSubMatrixHandler& SubMH);
	
	void MakeCCStructure(std::vector<integer>& Ai,
		std::vector<integer>& Ap);

protected:
        /* Matrix Matrix product */
	virtual MatrixHandler&
	MatMatMul_base(void (MatrixHandler::*op)(integer iRow, integer iCol,
				const doublereal& dCoef),
			MatrixHandler& out, const MatrixHandler& in) const;
	virtual MatrixHandler&
	MatTMatMul_base(void (MatrixHandler::*op)(integer iRow, integer iCol,
				const doublereal& dCoef),
			MatrixHandler& out, const MatrixHandler& in) const;

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


const doublereal&
NaiveMatrixHandler::operator () (integer iRow, integer iCol) const
{
	ASSERT(iRow > 0);
	ASSERT(iRow <= iGetNumRows());
	ASSERT(iCol > 0);
	ASSERT(iCol <= iGetNumCols());

	--iRow;
	--iCol;
	if (ppnonzero[iRow][iCol]) {
		return ppdRows[iRow][iCol];
	}
	return ::Zero1;
}

doublereal&
NaiveMatrixHandler::operator () (integer iRow, integer iCol)
{
	ASSERT(iRow > 0);
	ASSERT(iRow <= iGetNumRows());
	ASSERT(iCol > 0);
	ASSERT(iCol <= iGetNumCols());

	--iRow;
	--iCol;
	if (!(ppnonzero[iRow][iCol])) {
		ppnonzero[iRow][iCol] = 1;
		ppiRows[iCol][piNzr[iCol]] = iRow;
		ppiCols[iRow][piNzc[iRow]] = iCol;
		piNzr[iCol]++;
		piNzc[iRow]++;
		ppdRows[iRow][iCol] = 0.;
	}

	return ppdRows[iRow][iCol];
}

/* Sparse Matrix with unknowns permutation*/
class NaivePermMatrixHandler : public NaiveMatrixHandler {
protected:
	const std::vector<integer>& perm;
	const std::vector<integer>& invperm;

#ifdef DEBUG
	void IsValid(void) const {
		NO_OP;
	};
#endif /* DEBUG */

public:
	class const_iterator {
		friend class NaivePermMatrixHandler;

	private:
		const NaivePermMatrixHandler& m;
		mutable integer i_row;
		mutable SparseMatrixHandler::SparseMatrixElement elem;

	protected:
		void reset(bool is_end = false);

	public:
		const_iterator(const NaivePermMatrixHandler& m);
		const_iterator(const NaivePermMatrixHandler& m, bool);
		~const_iterator(void);
		const NaivePermMatrixHandler::const_iterator& operator ++ (void) const;
		const SparseMatrixHandler::SparseMatrixElement* operator -> (void) const;
		const SparseMatrixHandler::SparseMatrixElement& operator * (void) const;
		bool operator == (const NaivePermMatrixHandler::const_iterator& op) const;
		bool operator != (const NaivePermMatrixHandler::const_iterator& op) const;
	};

protected:
	const_iterator m_end;

public:
	NaivePermMatrixHandler::const_iterator begin(void) const {
		return const_iterator(*this);
	};

	const NaivePermMatrixHandler::const_iterator& end(void) const {
		return m_end;
	};

public:
	/* FIXME: always square? yes! */
	NaivePermMatrixHandler(integer iSize,
		const std::vector<integer>& tperm,
		const std::vector<integer>& invperm);

	NaivePermMatrixHandler(NaiveMatrixHandler*const nmh, 
		const std::vector<integer>& tperm,
		const std::vector<integer>& invperm);

	virtual ~NaivePermMatrixHandler(void);

	const std::vector<integer>& GetPerm(void) const;

	const std::vector<integer>& GetInvPerm(void) const;

	virtual inline const doublereal&
	operator () (integer iRow, integer iCol) const {
		ASSERT(iRow > 0);
		ASSERT(iRow <= iGetNumRows());
		ASSERT(iCol > 0);
		ASSERT(iCol <= iGetNumCols());

		ASSERT(perm.size() == (size_t)iGetNumRows());
		ASSERT(perm.size() == (size_t)iGetNumCols());

		/* FIXME: stupid 0/1 based arrays... */
		iCol = perm[iCol - 1] + 1;
		return NaiveMatrixHandler::operator()(iRow, iCol);
	};

	virtual inline doublereal&
	operator () (integer iRow, integer iCol) {
		ASSERT(iRow > 0);
		ASSERT(iRow <= iGetNumRows());
		ASSERT(iCol > 0);
		ASSERT(iCol <= iGetNumCols());

		ASSERT(perm.size() == (size_t)iGetNumRows());
		ASSERT(perm.size() == (size_t)iGetNumCols());

		/* FIXME: stupid 0/1 based arrays... */
		iCol = perm[iCol - 1] + 1;
		return NaiveMatrixHandler::operator()(iRow, iCol);
	};

protected:
        /* Matrix Matrix product */
	virtual MatrixHandler&
	MatMatMul_base(void (MatrixHandler::*op)(integer iRow, integer iCol,
				const doublereal& dCoef),
			MatrixHandler& out, const MatrixHandler& in) const;
	virtual MatrixHandler&
	MatTMatMul_base(void (MatrixHandler::*op)(integer iRow, integer iCol,
				const doublereal& dCoef),
			MatrixHandler& out, const MatrixHandler& in) const;

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


#endif /* NAIVEMH_H */

