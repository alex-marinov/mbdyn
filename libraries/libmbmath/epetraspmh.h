/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2022
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
        Copyright (C) 2022(-2022) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#ifndef __EPETRA_SPARSE_MATRIX_HANDLER_H__INCLUDED__
#define __EPETRA_SPARSE_MATRIX_HANDLER_H__INCLUDED__

#ifdef USE_TRILINOS
#include "myassert.h"
#include "spmh.h"
#include "cscmhtpl.h"

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>

class EpetraSparseMatrixHandler: public SparseMatrixHandler {
public:
     EpetraSparseMatrixHandler(const integer& iNumRows,
                               const integer& iNumCols,
                               integer iNumColsAlloc,
                               const Epetra_Comm& oComm);
     virtual ~EpetraSparseMatrixHandler();

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
     int32_t MakeCompressedRowForm(std::vector<doublereal>& Ax,
                                   std::vector<int32_t>& Ai,
                                   std::vector<int32_t>& Ap,
                                   int offset = 0) const override;

     virtual
     int64_t MakeCompressedRowForm(std::vector<doublereal>& Ax,
                                   std::vector<int64_t>& Ai,
                                   std::vector<int64_t>& Ap,
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
     int32_t MakeCompressedColumnForm(std::vector<doublereal>& Ax,
                                      std::vector<int32_t>& Ai,
                                      std::vector<int32_t>& Ap,
                                      int offset = 0) const override;

     virtual
     int64_t MakeCompressedColumnForm(std::vector<doublereal>& Ax,
                                      std::vector<int64_t>& Ai,
                                      std::vector<int64_t>& Ap,
                                      int offset = 0) const override;

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
     VectorHandler& GetCol(integer icol,
                           VectorHandler& out) const override;

     virtual void Scale(const std::vector<doublereal>& oRowScale, const std::vector<doublereal>& oColScale) override;

     virtual bool AddItem(integer iRow, const sp_grad::SpGradient& oItem) override;

     virtual void EnumerateNz(const std::function<EnumerateNzCallback>& func) const override;

     virtual doublereal Norm(Norm_t eNorm = NORM_1) const override;

     virtual integer Nz() const override;

     virtual EpetraSparseMatrixHandler* Copy() const override;

     virtual integer PacMat() override;

     class const_iterator {
          friend class EpetraSparseMatrixHandler;
          
          const integer* const rowptr;
          const integer* const colind;
          const doublereal* const values;
          const integer NRows;
          const integer NZ;
	  integer iIdx;
	  SparseMatrixHandler::SparseMatrixElement elem;
          
#ifdef DEBUG
	  bool bInvariant() const {
	       ASSERT(iIdx >= 0);
	       ASSERT(iIdx <= NZ);
	       ASSERT(iIdx == NZ || (elem.iRow >= 0 && elem.iRow < NRows));
	       ASSERT(iIdx == NZ || rowptr[elem.iRow] <= iIdx);
	       ASSERT(iIdx == NZ || rowptr[elem.iRow + 1] > iIdx);
               
	       return true;
	  }
#endif
	  void UpdateElem() {
	       if (iIdx < NZ) {
		    elem.iCol = colind[iIdx];
		    elem.dCoef = values[iIdx];

		    if (iIdx >= rowptr[elem.iRow + 1]) {
			 ++elem.iRow;
		    }
	       } else {
		    ASSERT(iIdx == NZ);
#ifdef DEBUG
		    elem.iRow = std::numeric_limits<decltype(elem.iRow)>::min();
		    elem.iCol = std::numeric_limits<decltype(elem.iCol)>::min();
		    elem.dCoef = -std::numeric_limits<decltype(elem.dCoef)>::max();
#endif
	       }
	  }

	  const_iterator(const integer* rowptr,
                         const integer* colind,
                         const doublereal* values,
                         integer NRows,                         
                         integer iIdx,
                         integer iRow)
	       :rowptr(rowptr),
                colind(colind),
                values(values),
                NRows(NRows),
                NZ(rowptr[NRows] - rowptr[0]),
                iIdx(iIdx) {

	       elem.iRow = iRow;

	       UpdateElem();

	       ASSERT(bInvariant());
	  }

     public:
	  ~const_iterator() {
	       ASSERT(bInvariant());
	  }

	  const const_iterator& operator++() {
	       ASSERT(bInvariant());

	       ++iIdx;

	       UpdateElem();

	       ASSERT(bInvariant());

	       return *this;
	  }

	  const SparseMatrixHandler::SparseMatrixElement* operator->() const {
	       ASSERT(bInvariant());

	       return &elem;
	  }

	  const SparseMatrixHandler::SparseMatrixElement& operator*() const {
	       ASSERT(bInvariant());

	       return elem;
	  }

	  bool operator == (const const_iterator& op) const {
	       ASSERT(bInvariant());
	       ASSERT(rowptr == op.rowptr);
               ASSERT(colind == op.colind);
               ASSERT(values == op.values);

	       return iIdx == op.iIdx;
	  }

	  bool operator != (const const_iterator& op) const {
	       ASSERT(bInvariant());
	       ASSERT(rowptr == op.rowptr);
               ASSERT(colind == op.colind);
               ASSERT(values == op.values);
               
	       return iIdx != op.iIdx;
	  }
     };

     const_iterator begin() const;
     const_iterator end() const;

     const Epetra_CrsMatrix* pGetEpetraCrsMatrix() const {
          return &oEPM;
     }

     Epetra_CrsMatrix* pGetEpetraCrsMatrix() {
          return &oEPM;
     }
protected:
     virtual VectorHandler&
     MatVecMul_base(void (VectorHandler::*op)(integer iRow,
                                              const doublereal& dCoef),
                    VectorHandler& out, const VectorHandler& in) const override;
     virtual VectorHandler&
     MatTVecMul_base(void (VectorHandler::*op)(integer iRow,
                                               const doublereal& dCoef),
                     VectorHandler& out, const VectorHandler& in) const override;
private:
     void ExtractCrsDataPointers(integer*& rowptr, integer*& colind, doublereal*& values) const;
     CSCMatrixHandlerTpl<doublereal, integer, 0>& GetTransposedCSC() const;
     
     template <typename idx_type>
     idx_type MakeCompressedRowFormTpl(doublereal *const Ax,
                                          idx_type *const Ai,
                                          idx_type *const Ap,
                                          int offset) const;

     template <typename idx_type>
     idx_type MakeCompressedColumnFormTpl(doublereal *const Ax,
                                          idx_type *const Ai,
                                          idx_type *const Ap,
                                          int offset) const;
     
     inline integer PacMat() const {
          return const_cast<EpetraSparseMatrixHandler*>(this)->PacMat();
     }
     
     const Epetra_Comm& oComm;
     mutable Epetra_CrsMatrix oEPM;
     const integer iNumColsAlloc;
     mutable CSCMatrixHandlerTpl<doublereal, integer, 0> oCscT;
};

#endif
#endif
