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
/* November 2001 
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
 * Copyright (C) 1996-2017
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

#ifndef SpMapMatrixHandler_hh
#define SpMapMatrixHandler_hh

#include <map>
#include <vector>
#include "myassert.h"
#include "solman.h"
#include "spmh.h"

/* #define MBDYN_X_KEEP_SPARSITY */

/* Sparse Matrix in columns form */
class SpMapMatrixHandler : public SparseMatrixHandler {
private:
        integer NZ;
	typedef std::map<integer, doublereal> row_cont_type;
	mutable std::vector<row_cont_type> col_indices;

	// don't allow copy constructor!
	SpMapMatrixHandler(const SpMapMatrixHandler&);

	template <typename idx_type>
	idx_type MakeCompressedColumnFormTpl(doublereal *const Ax,
					     idx_type *const Ai,
					     idx_type *const Ap,
					     int offset) const;

	template <typename idx_type>
	idx_type MakeCompressedRowFormTpl(doublereal *const Ax,
					  idx_type *const Ai,
					  idx_type *const Ap,
					  int offset) const;
     
        template <typename idx_type>
	idx_type MakeIndexFormTpl(doublereal *const Ax,
				  idx_type *const Arow, idx_type *const Acol,
				  idx_type *const AcolSt,
				  int offset = 0) const;
#ifdef DEBUG
	void IsValid(void) const {
		NO_OP;
	};
#endif /* DEBUG */

public:
	class const_iterator {
		friend class SpMapMatrixHandler;

	private:
		const SpMapMatrixHandler& m;
		mutable row_cont_type::const_iterator i;
		mutable SparseMatrixHandler::SparseMatrixElement elem;

	protected:
		void reset(bool is_end = false);

	public:
		const_iterator(const SpMapMatrixHandler& m, bool is_end = false);
		~const_iterator(void);
		const SpMapMatrixHandler::const_iterator& operator ++ (void) const;
		const SparseMatrixHandler::SparseMatrixElement* operator -> (void) const;
		const SparseMatrixHandler::SparseMatrixElement& operator * (void) const;
		bool operator == (const SpMapMatrixHandler::const_iterator& op) const;
		bool operator != (const SpMapMatrixHandler::const_iterator& op) const;
	};

private:
	const_iterator m_end;

public:
	SpMapMatrixHandler::const_iterator begin(void) const {
		return const_iterator(*this);
	};

	const SpMapMatrixHandler::const_iterator& end(void) const {
		return m_end;
	};

public:
	SpMapMatrixHandler(const integer &n = 0,const integer &nn = 0);

	virtual ~SpMapMatrixHandler(void);

	doublereal & operator()(integer i_row, integer i_col) {
		ASSERTMSGBREAK(i_row > 0 && i_row <= NRows,
				"Error in SpMapMatrixHandler::operator(), "
				"row index out of range");
		ASSERTMSGBREAK(i_col > 0 && i_col <= NCols, 
				"Error in SpMapMatrixHandler::operator(), "
				"col index out of range");
		i_row--;
		i_col--;
		row_cont_type::iterator i;
		row_cont_type& row = col_indices[i_col];
		i = row.find(i_row);
		if (i == row.end()) {
			NZ++;

			row[i_row] = 0.;

			return row[i_row];

		} else {
			return i->second;
		}
	};

	void IncCoef(integer ix, integer iy, const doublereal& inc) {
		ASSERTMSGBREAK(ix > 0 && ix <= NRows,
				"Error in SpMapMatrixHandler::IncCoef(), "
				"row index out of range");
		ASSERTMSGBREAK(iy > 0 && iy <= NCols,
				"Error in SpMapMatrixHandler::IncCoef(), "
				"col index out of range");
#ifdef MBDYN_X_KEEP_SPARSITY
		/* try to keep sparsity */
		if (inc != 0.) {
#endif /* MBDYN_X_KEEP_SPARSITY */
			operator()(ix,iy) += inc;
#ifdef MBDYN_X_KEEP_SPARSITY
		}
#endif /* MBDYN_X_KEEP_SPARSITY */
	};

	void DecCoef(integer ix, integer iy, const doublereal& inc) {
		ASSERTMSGBREAK(ix > 0 && ix <= NRows,
				"Error in SpMapMatrixHandler::DecCoef(), "
				"row index out of range");
		ASSERTMSGBREAK(iy > 0 && iy <= NCols,
				"Error in SpMapMatrixHandler::DecCoef(), "
				"col index out of range");
#ifdef MBDYN_X_KEEP_SPARSITY
		/* try to keep sparsity */
		if (inc != 0.) {
#endif /* MBDYN_X_KEEP_SPARSITY */
			operator()(ix,iy) -= inc;
#ifdef MBDYN_X_KEEP_SPARSITY
		}
#endif /* MBDYN_X_KEEP_SPARSITY */
	};

	void PutCoef(integer ix, integer iy, const doublereal& val) {
		ASSERTMSGBREAK(ix - 1 < NRows,
				"Error in SpMapMatrixHandler::PutCoef(), "
				"row index out of range");
		ASSERTMSGBREAK(iy - 1 < NCols,
				"Error in SpMapMatrixHandler::PutCoef(), "
				"col index out of range");
#ifdef MBDYN_X_KEEP_SPARSITY
		/* try to keep sparsity */
		if (val != 0.) {
#endif /* MBDYN_X_KEEP_SPARSITY */
			operator()(ix,iy) = val;
#ifdef MBDYN_X_KEEP_SPARSITY
		} else {
			row_cont_type::iterator i;
			row_cont_type& row = col_indices[iy-1];
			i = row.find(ix - 1);
			if (i != row.end()) {
				i->second = val;
			}
		}
#endif /* MBDYN_X_KEEP_SPARSITY */
	};

	const doublereal& dGetCoef(integer ix, integer iy) const {
		ASSERTMSGBREAK(ix > 0 && ix <= NRows,
				"Error in SpMapMatrixHandler::dGetCoef(), "
				"row index out of range");
		ASSERTMSGBREAK(iy > 0 && iy <= NCols,
				"Error in SpMapMatrixHandler::dGetCoef(), "
				"col index out of range");
		row_cont_type::iterator i;
		row_cont_type& row = col_indices[iy - 1];
		i = row.find(ix - 1);
		if (i == row.end()) {
			return ::Zero1;

		} else {
			return i->second;
		}
	};

	const doublereal& operator () (integer ix, integer iy) const {
		ASSERTMSGBREAK(ix > 0 && ix <= NRows,
				"Error in SpMapMatrixHandler::operator(), "
				"row index out of range");
		ASSERTMSGBREAK(iy > 0 && iy <= NCols,
				"Error in SpMapMatrixHandler::operator(), "
				"col index out of range");
		row_cont_type::iterator i;
		row_cont_type& row = col_indices[iy - 1];
		i = row.find(ix - 1);
		if (i == row.end()) {
			return ::Zero1;

		} else {
			return i->second;
		}
	};

	using SparseMatrixHandler::MakeCompressedColumnForm;
		
	int32_t MakeCompressedColumnForm(doublereal *const Ax,
					 int32_t *const Ai,
					 int32_t *const Ap,
					 int offset = 0) const override;

	int64_t MakeCompressedColumnForm(doublereal *const Ax,
					 int64_t *const Ai,
					 int64_t *const Ap,
					 int offset = 0) const override;

        int32_t MakeCompressedRowForm(doublereal *const Ax,
				      int32_t *const Ai,
				      int32_t *const Ap,
				      int offset = 0) const override;

       int64_t MakeCompressedRowForm(doublereal *const Ax,
				     int64_t *const Ai,
				     int64_t *const Ap,
				     int offset = 0) const override;
     
       int32_t MakeIndexForm(doublereal *const Ax,
			     int32_t *const Arow, int32_t *const Acol,
			     int32_t *const AcolSt,
			     int offset = 0) const override;

       int64_t MakeIndexForm(doublereal *const Ax,
			     int64_t *const Arow, int64_t *const Acol,
			     int64_t *const AcolSt,
			     int offset = 0) const override;


	void Reset(void);

	void Resize(integer ir, integer ic);

	/* Estrae una colonna da una matrice */
	VectorHandler& GetCol(integer icol, VectorHandler& out) const;

        virtual void Scale(const std::vector<doublereal>& oRowScale, const std::vector<doublereal>& oColScale) override;

        virtual integer Nz() const override;
	/* Matrix Matrix product */
protected:
	MatrixHandler&
	MatMatMul_base(void (MatrixHandler::*op)(integer iRow, integer iCol,
				const doublereal& dCoef),
			MatrixHandler& out, const MatrixHandler& in) const;
	MatrixHandler&
	MatTMatMul_base(void (MatrixHandler::*op)(integer iRow, integer iCol,
				const doublereal& dCoef),
			MatrixHandler& out, const MatrixHandler& in) const;
public:

        /* Moltiplica per uno scalare e somma a una matrice */
	MatrixHandler& MulAndSumWithShift(MatrixHandler& out,
			doublereal s = 1.,
			integer drow = 0, integer dcol = 0) const;
	
	MatrixHandler& FakeThirdOrderMulAndSumWithShift(MatrixHandler& out, 
			std::vector<bool> b, doublereal s = 1.,
			integer drow = 0, integer dcol = 0) const;
	
	/* Matrix Vector product */
protected:
	virtual VectorHandler&
	MatVecMul_base(void (VectorHandler::*op)(integer iRow,
				const doublereal& dCoef),
			VectorHandler& out, const VectorHandler& in) const;
	virtual VectorHandler&
	MatTVecMul_base(void (VectorHandler::*op)(integer iRow,
				const doublereal& dCoef),
			VectorHandler& out, const VectorHandler& in) const;

public:
        virtual void EnumerateNz(const std::function<EnumerateNzCallback>& func) const override;
        virtual SpMapMatrixHandler* Copy() const override;
};

#endif /* SpMapMatrixHandler_hh */

