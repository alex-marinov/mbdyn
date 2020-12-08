/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2020
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

/*
 AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
        Copyright (C) 2020(-2020) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#ifndef __SP_GRADIENT_SPARSE_MATRIX_HANDLER_H__INCLUDED__
#define __SP_GRADIENT_SPARSE_MATRIX_HANDLER_H__INCLUDED__

#ifdef USE_SPARSE_AUTODIFF
#include <vector>

#include "myassert.h"
#include "spmh.h"

#ifdef USE_MULTITHREAD
#include <atomic>
#endif

#include "sp_gradient.h"

class SpGradientSparseMatrixHandler: public SparseMatrixHandler {
     struct SparseRow;
public:
     SpGradientSparseMatrixHandler(const integer& iNumRows, const integer& iNumCols);
     virtual ~SpGradientSparseMatrixHandler();

     class const_iterator {
     public:
          inline explicit const_iterator(const std::vector<SparseRow>& oRows,
                                         std::vector<SparseRow>::const_iterator pCurrRow,
                                         const sp_grad::SpDerivRec* pCurrCol);
          inline const const_iterator& operator++ (void);
          inline const SparseMatrixHandler::SparseMatrixElement* operator->() const;
          inline const SparseMatrixHandler::SparseMatrixElement& operator*() const;
          inline bool operator==(const const_iterator& op) const;
          inline bool operator!=(const const_iterator& op) const;

     private:
          inline void Update();
          const std::vector<SparseRow>& oRows;
          std::vector<SparseRow>::const_iterator pCurrRow;
          const sp_grad::SpDerivRec* pCurrCol;
          mutable SparseMatrixHandler::SparseMatrixElement elem;
     };

     inline const_iterator begin() const;
     inline const_iterator end() const;

#ifdef DEBUG
     virtual void IsValid() const override;
#endif

     virtual void Resize(integer, integer) override;

     virtual void ResizeReset(integer, integer) override;

     virtual void Reset() override;

     virtual const doublereal&
     operator()(integer iRow, integer iCol) const override;

     virtual doublereal&
     operator()(integer iRow, integer iCol) override;

     using SparseMatrixHandler::MakeCompressedColumnForm;
     using SparseMatrixHandler::MakeCompressedRowForm;

     virtual
     int32_t MakeCompressedColumnForm(doublereal *const Ax,
				      int32_t *const Ai,
				      int32_t *const Ap,
				      int offset = 0) const override;

     virtual
     int64_t MakeCompressedColumnForm(doublereal *const Ax,
				      int64_t *const Ai,
				      int64_t *const Ap,
				      int offset = 0) const override;

     virtual
     int32_t MakeCompressedRowForm(doublereal *const Ax,
				   int32_t *const Ai,
				   int32_t *const Ap,
				   int offset = 0) const override;

     virtual
     int64_t MakeCompressedRowForm(doublereal *const Ax,
				   int64_t *const Ai,
				   int64_t *const Ap,
				   int offset = 0) const override;

     virtual
     int32_t MakeIndexForm(doublereal *const Ax,
			   int32_t *const Arow, int32_t *const Acol,
			   int32_t *const AcolSt,
			   int offset = 0) const override;

     virtual
     int64_t MakeIndexForm(doublereal *const Ax,
			   int64_t *const Arow, int64_t *const Acol,
			   int64_t *const AcolSt,
			   int offset = 0) const override;

     virtual
     VectorHandler& GetCol(integer icol,
			   VectorHandler& out) const override;

     virtual void Scale(const std::vector<doublereal>& oRowScale, const std::vector<doublereal>& oColScale) override;

     virtual bool AddItem(integer iRow, const sp_grad::SpGradient& oItem) override;

     std::ostream& Print(std::ostream& os, MatPrintFormat eFormat) const override;

     virtual doublereal Norm(Norm_t eNorm = NORM_1) const override;

     virtual integer Nz() const override;

protected:
     virtual MatrixHandler&
     MatMatMul_base(void (MatrixHandler::*op)(integer iRow, integer iCol,
					      const doublereal& dCoef),
		    MatrixHandler& out, const MatrixHandler& in) const override;
     virtual MatrixHandler&
     MatTMatMul_base(void (MatrixHandler::*op)(integer iRow, integer iCol,
					       const doublereal& dCoef),
		     MatrixHandler& out, const MatrixHandler& in) const override;

     VectorHandler&
     MatVecMul_base(
	  void (VectorHandler::*op)(integer iRow, const doublereal& dCoef),
	  VectorHandler& out, const VectorHandler& in) const override;

     VectorHandler&
     MatTVecMul_base(
	  void (VectorHandler::*op)(integer iRow, const doublereal& dCoef),
	  VectorHandler& out, const VectorHandler& in) const override;

private:
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
			       idx_type *const Arow,
			       idx_type* const Acol,
			       idx_type *const Ap,
			       int offset) const;

     struct SparseRow: sp_grad::SpGradient {
#ifdef USE_MULTITHREAD
	  SparseRow() {
	       bLocked = false;
	  }

	  mutable std::atomic<bool> bLocked;
#endif
     };

#ifdef USE_MULTITHREAD
     std::atomic<integer> NZ;
#else
     integer NZ;
#endif
     std::vector<SparseRow> oRows;
};

SpGradientSparseMatrixHandler::const_iterator SpGradientSparseMatrixHandler::begin() const {
     return const_iterator(oRows, oRows.begin(), oRows.begin()->begin());
}

SpGradientSparseMatrixHandler::const_iterator SpGradientSparseMatrixHandler::end() const {
     return const_iterator(oRows, oRows.end(), nullptr);
}

SpGradientSparseMatrixHandler::const_iterator::const_iterator(const std::vector<SparseRow>& oRows,
							      std::vector<SparseRow>::const_iterator pCurrRow,
							      const sp_grad::SpDerivRec* pCurrCol)
     :oRows(oRows), pCurrRow(pCurrRow), pCurrCol(pCurrCol) {

     SP_GRAD_ASSERT(pCurrRow <= oRows.end());
     SP_GRAD_ASSERT(pCurrRow != oRows.end() ? pCurrCol < pCurrRow->end() : true);

     Update();
}

void SpGradientSparseMatrixHandler::const_iterator::Update()
{
     if (pCurrRow < oRows.end() && pCurrCol < pCurrRow->end()) {
	  elem.iRow = pCurrRow - oRows.begin();
	  elem.iCol = pCurrCol->iDof - 1;
	  elem.dCoef = pCurrCol->dDer;

	  SP_GRAD_ASSERT(elem.iRow >= 0);
	  SP_GRAD_ASSERT(static_cast<size_t>(elem.iRow) < oRows.size());
	  SP_GRAD_ASSERT(elem.iCol >= 0);
	  SP_GRAD_ASSERT(static_cast<size_t>(elem.iCol) < oRows.size());
     } else {
#ifdef SP_GRAD_DEBUG
	  elem.iRow = -std::numeric_limits<integer>::max();
	  elem.iCol = -std::numeric_limits<integer>::max();
	  elem.dCoef = std::numeric_limits<doublereal>::infinity();
#endif
     }
}

const SpGradientSparseMatrixHandler::const_iterator& SpGradientSparseMatrixHandler::const_iterator::operator++()
{
     SP_GRAD_ASSERT(pCurrRow < oRows.end());
     SP_GRAD_ASSERT(pCurrCol < pCurrRow->end());

     if (++pCurrCol >= pCurrRow->end()) {
	  if (++pCurrRow < oRows.end()) {
	       pCurrCol = pCurrRow->begin();
	  }
     }

     Update();

     return *this;
}
const SparseMatrixHandler::SparseMatrixElement* SpGradientSparseMatrixHandler::const_iterator::operator->() const
{
     SP_GRAD_ASSERT(pCurrRow >= oRows.begin());
     SP_GRAD_ASSERT(pCurrRow < oRows.end());
     SP_GRAD_ASSERT(pCurrCol >= pCurrRow->begin());
     SP_GRAD_ASSERT(pCurrCol < pCurrRow->end());

     return &elem;
}

const SparseMatrixHandler::SparseMatrixElement& SpGradientSparseMatrixHandler::const_iterator::operator*() const
{
     SP_GRAD_ASSERT(pCurrRow >= oRows.begin());
     SP_GRAD_ASSERT(pCurrRow < oRows.end());
     SP_GRAD_ASSERT(pCurrCol >= pCurrRow->begin());
     SP_GRAD_ASSERT(pCurrCol < pCurrRow->end());

     return elem;
}

bool SpGradientSparseMatrixHandler::const_iterator::operator==(const const_iterator& op) const
{
     SP_GRAD_ASSERT(pCurrRow >= oRows.begin());
     SP_GRAD_ASSERT(pCurrRow <= oRows.end());
     SP_GRAD_ASSERT(op.pCurrRow >= oRows.begin());
     SP_GRAD_ASSERT(op.pCurrRow <= oRows.end());

     return pCurrRow < oRows.end() && op.pCurrRow < oRows.end() ? pCurrRow == op.pCurrRow && pCurrCol == op.pCurrCol : pCurrRow == op.pCurrRow;
}

bool SpGradientSparseMatrixHandler::const_iterator::operator!=(const const_iterator& op) const
{
     return !(*this == op);
}

#ifdef USE_MULTITHREAD
class SpGradientSparseMatrixWrapper: public MatrixHandler
{
public:
     explicit SpGradientSparseMatrixWrapper(SpGradientSparseMatrixHandler* pMH = nullptr);
     ~SpGradientSparseMatrixWrapper();

#ifdef DEBUG
     virtual void IsValid() const override;
#endif
     void SetMatrixHandler(SpGradientSparseMatrixHandler* pMH);

     virtual integer iGetNumRows(void) const override;
     virtual integer iGetNumCols(void) const override;

     virtual void Resize(integer, integer) override;

     virtual void ResizeReset(integer, integer) override;

     virtual void Reset() override;

     virtual const doublereal&
     operator()(integer iRow, integer iCol) const override;

     virtual doublereal&
     operator()(integer iRow, integer iCol) override;

     virtual void Scale(const std::vector<doublereal>& oRowScale, const std::vector<doublereal>& oColScale) override;

     virtual bool AddItem(integer iRow, const sp_grad::SpGradient& oItem) override;

     std::ostream& Print(std::ostream& os, MatPrintFormat eFormat) const override;

     virtual doublereal Norm(Norm_t eNorm = NORM_1) const override;

protected:
     virtual MatrixHandler&
     MatMatMul_base(void (MatrixHandler::*op)(integer iRow, integer iCol,
					      const doublereal& dCoef),
		    MatrixHandler& out, const MatrixHandler& in) const override;
     virtual MatrixHandler&
     MatTMatMul_base(void (MatrixHandler::*op)(integer iRow, integer iCol,
					       const doublereal& dCoef),
		     MatrixHandler& out, const MatrixHandler& in) const override;

     VectorHandler&
     MatVecMul_base(
	  void (VectorHandler::*op)(integer iRow, const doublereal& dCoef),
	  VectorHandler& out, const VectorHandler& in) const override;

     VectorHandler&
     MatTVecMul_base(
	  void (VectorHandler::*op)(integer iRow, const doublereal& dCoef),
	  VectorHandler& out, const VectorHandler& in) const override;

private:
     SpGradientSparseMatrixHandler* pMH;
};
#endif

#endif
#endif
