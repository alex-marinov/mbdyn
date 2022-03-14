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

#include "mbconfig.h"

#ifdef USE_TRILINOS
#include <algorithm>
#include <iomanip>
#include <new>

#ifdef USE_SPARSE_AUTODIFF
#include "sp_gradient.h"
#endif
#undef HAVE_BLAS
#include "epetraspmh.h"
#include "epetravh.h"
#include "cscmhtpl.h"

EpetraSparseMatrixHandler::EpetraSparseMatrixHandler(const integer& iNumRows, const integer& iNumCols, integer iNumColsAlloc, const Epetra_Comm& oComm)
     :SparseMatrixHandler(iNumRows, iNumCols),
      oComm(oComm),
      oEPM(::Copy, Epetra_Map(NRows, 1, oComm), Epetra_Map(NCols, 1, oComm), iNumColsAlloc),
      iNumColsAlloc(iNumColsAlloc)
{
     if (iNumRows != iNumCols) {
          silent_cerr("EpetraSparseMatrixHandler: matrix must be square!\n");
          throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
     }

#ifdef DEBUG
     IsValid();
#endif     
}

EpetraSparseMatrixHandler::~EpetraSparseMatrixHandler()
{
#ifdef DEBUG
     IsValid();
#endif
}

#ifdef DEBUG
void EpetraSparseMatrixHandler::IsValid() const
{
     ASSERT(NRows == oEPM.NumGlobalRows());
     ASSERT(NCols == oEPM.NumGlobalCols());
}
#endif

void EpetraSparseMatrixHandler::Resize(integer, integer)
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

void EpetraSparseMatrixHandler::ResizeReset(integer iNumRows, integer iNumCols)
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

void EpetraSparseMatrixHandler::Reset()
{
#ifdef DEBUG
     IsValid();
#endif
     
     if (oEPM.Filled()) {
          integer* rowptr;
          integer* colind;
          doublereal* values;

          ExtractCrsDataPointers(rowptr, colind, values);
          
          std::vector<integer> rowSize;

          rowSize.reserve(NRows);

          for (integer i = 0; i < NRows; ++i) {
               rowSize.push_back(rowptr[i + 1] - rowptr[i]);
          }

          Epetra_Map oMapRows(NRows, 1, oComm);
          Epetra_Map oMapCols(NCols, 1, oComm);
          
          oCscT = CSCMatrixHandlerTpl<doublereal, integer, 0>();

          oEPM.~Epetra_CrsMatrix();
          new (&oEPM) Epetra_CrsMatrix(::Copy, oMapRows, oMapCols, &rowSize.front());
     } else {
          oEPM.PutScalar(0.);
     }
     
     ASSERT(!oEPM.Filled());

#ifdef DEBUG
     IsValid();
#endif     
}

VectorHandler&
EpetraSparseMatrixHandler::MatVecMul_base(void (VectorHandler::*op)(integer iRow,
                                                                    const doublereal& dCoef),
                                          VectorHandler& out, const VectorHandler& in) const
{
     return GetTransposedCSC().MatTVecMul_base(op, out, in);
}

VectorHandler&
EpetraSparseMatrixHandler::MatTVecMul_base(void (VectorHandler::*op)(integer iRow,
                                          const doublereal& dCoef),
                                           VectorHandler& out, const VectorHandler& in) const
{
     return GetTransposedCSC().MatVecMul_base(op, out, in);
}

void EpetraSparseMatrixHandler::FillComplete() const
{
#ifdef DEBUG
     IsValid();
#endif

     if (!oEPM.Filled()) {
          oEPM.FillComplete(true);

          integer* rowptr;
          integer* colind;
          doublereal* values;

          ExtractCrsDataPointers(rowptr, colind, values);

          oCscT = CSCMatrixHandlerTpl<doublereal, integer, 0>(values, colind, rowptr, NRows, Nz());
     }
#ifdef DEBUG
     IsValid();
#endif
}

EpetraSparseMatrixHandler::const_iterator EpetraSparseMatrixHandler::begin() const {
     integer* rowptr;
     integer* colind;
     doublereal* values;

     ExtractCrsDataPointers(rowptr, colind, values);

     ASSERT(NRows == oEPM.NumGlobalRows());
     
     return const_iterator(rowptr, colind, values, oEPM.NumGlobalRows(), 0, 0);
}

EpetraSparseMatrixHandler::const_iterator EpetraSparseMatrixHandler::end() const {
     integer* rowptr;
     integer* colind;
     doublereal* values;

     ExtractCrsDataPointers(rowptr, colind, values);

     ASSERT(NRows == oEPM.NumGlobalRows());
     
     return const_iterator(rowptr, colind, values, NRows, Nz(), NRows - 1);
}    

void EpetraSparseMatrixHandler::ExtractCrsDataPointers(integer*& rowptr, integer*& colind, doublereal*& values) const
{
     ASSERT(oEPM.Filled());
     ASSERT(oEPM.StorageOptimized());
     
     oEPM.ExtractCrsDataPointers(rowptr, colind, values);

     ASSERT(rowptr != nullptr);
     ASSERT(colind != nullptr);
     ASSERT(values != nullptr);
}

CSCMatrixHandlerTpl<doublereal, integer, 0>& EpetraSparseMatrixHandler::GetTransposedCSC() const
{
     ASSERT(oEPM.Filled());
     
#ifdef DEBUG
     IsValid();
     oCscT.IsValid();
#endif
     
     return oCscT;
}

const doublereal&
EpetraSparseMatrixHandler::operator()(integer iRow, integer iCol) const
{
     return GetTransposedCSC()(iCol, iRow);
}

doublereal&
EpetraSparseMatrixHandler::operator()(integer iRow, integer iCol)
{
     return GetTransposedCSC()(iCol, iRow);
}

template <typename idx_type>
idx_type EpetraSparseMatrixHandler::MakeCompressedColumnFormTpl(doublereal *const Ax,
                                                                idx_type *const Ai,
                                                                idx_type *const Ap,
                                                                const int offset) const
{
#ifdef DEBUG
     IsValid();
#endif
     integer* rowptr;
     integer* colind;
     doublereal* values;

     ASSERT(oEPM.Filled());
     ASSERT(oEPM.StorageOptimized());
     
     ExtractCrsDataPointers(rowptr, colind, values);

     const integer iNumRows = iGetNumRows();
     const integer iNumCols = iGetNumCols();
     const integer iNumNz = oEPM.NumGlobalNonzeros();
     std::vector<idx_type> oColSize(iNumCols, 0);

     for (integer i = 0; i < iNumNz; ++i) {
          ASSERT(colind[i] >= 0);
          ASSERT(colind[i] < iNumCols);
          
          ++oColSize[colind[i]];
     }

     idx_type iPtr = offset;

     for (integer iCol = 0; iCol < iNumCols; ++iCol) {
	  Ap[iCol] = iPtr;

          ASSERT(oColSize[iCol] >= 0);
          ASSERT(oColSize[iCol] <= iNumRows);
          
	  iPtr += oColSize[iCol];
     }

     Ap[iNumCols] = iPtr;

     ASSERT(iPtr - offset == iNumNz);

     std::fill(oColSize.begin(), oColSize.end(), 0);

#ifdef DEBUG
     constexpr idx_type iInvalidIndex = -1;
     constexpr doublereal dInvalidValue = std::numeric_limits<doublereal>::infinity();

     std::fill(Ai, Ai + iNumNz, iInvalidIndex);
     std::fill(Ax, Ax + iNumNz, dInvalidValue);
#endif

     for (integer iRow = 0; iRow < iNumRows; ++iRow) {
	  for (integer iItem = rowptr[iRow]; iItem < rowptr[iRow + 1]; ++iItem) {
	       const integer iColIdx = colind[iItem];
	       const idx_type iSlot = Ap[iColIdx] - offset + oColSize[iColIdx]++;

	       ASSERT(iSlot >= 0);
	       ASSERT(iSlot < iNumNz);
	       ASSERT(Ai[iSlot] == iInvalidIndex);
	       ASSERT(Ax[iSlot] == dInvalidValue);

	       Ai[iSlot] = iRow + offset;
	       Ax[iSlot] = values[iItem];
	  }
     }

#ifdef DEBUG
     idx_type iSizeTot = 0;

     for (idx_type iColSize: oColSize) {
	  iSizeTot += iColSize;
     }

     ASSERT(iSizeTot == iNumNz);

     for (integer iCol = 0; iCol < iNumCols; ++iCol) {
	  idx_type iRowPrev = -1;
	  for (idx_type iSlot = Ap[iCol] - offset; iSlot < Ap[iCol + 1] - offset; ++iSlot) {
	       idx_type iRow = Ai[iSlot] - offset;

	       ASSERT(iRow > iRowPrev);
	       ASSERT(iRow >= 0);
	       ASSERT(iRow < NRows);
	       ASSERT(std::isfinite(Ax[iSlot]));

	       for (integer iItem = rowptr[iRow]; iItem < rowptr[iRow + 1]; ++iItem) {
		    if (colind[iItem] == iCol) {
			 ASSERT(Ax[iSlot] == values[iItem]);
		    }
	       }

	       iRowPrev = iRow;
	  }
     }
#endif

     return iPtr - offset;     
}


int32_t EpetraSparseMatrixHandler::MakeCompressedColumnForm(std::vector<doublereal>& Ax,
                                                            std::vector<int32_t>& Ai,
                                                            std::vector<int32_t>& Ap,
                                                            int offset) const
{
     FillComplete(); // FillComplete must be called before Nz()

     Ai.resize(Nz());
     Ax.resize(Nz());
     Ap.resize(iGetNumCols() + 1);

     return MakeCompressedColumnForm(&Ax.front(), &Ai.front(), &Ap.front(), offset);
}


int64_t EpetraSparseMatrixHandler::MakeCompressedColumnForm(std::vector<doublereal>& Ax,
                                                            std::vector<int64_t>& Ai,
                                                            std::vector<int64_t>& Ap,
                                                            int offset) const
{
     FillComplete(); // FillComplete must be called before Nz()

     Ai.resize(Nz());
     Ax.resize(Nz());
     Ap.resize(iGetNumCols() + 1);

     return MakeCompressedColumnForm(&Ax.front(), &Ai.front(), &Ap.front(), offset);     
}


int32_t EpetraSparseMatrixHandler::MakeCompressedColumnForm(doublereal *const Ax,
								int32_t *const Ai,
								int32_t *const Ap,
								int offset) const
{
     return MakeCompressedColumnFormTpl(Ax, Ai, Ap, offset);
}

int64_t EpetraSparseMatrixHandler::MakeCompressedColumnForm(doublereal *const Ax,
								int64_t *const Ai,
								int64_t *const Ap,
								int offset) const
{
     return MakeCompressedColumnFormTpl(Ax, Ai, Ap, offset);
}

template <typename idx_type>
idx_type EpetraSparseMatrixHandler::MakeCompressedRowFormTpl(doublereal *const Ax,
                                                             idx_type *const Ai,
                                                             idx_type *const Ap,
                                                             int offset) const
{
#ifdef DEBUG
     IsValid();
#endif
     
     integer* rowptr;
     integer* colind;
     doublereal* values;

     ASSERT(oEPM.Filled());
     ASSERT(oEPM.StorageOptimized());
     
     ExtractCrsDataPointers(rowptr, colind, values);

     const integer iNumRows = iGetNumRows();
     const integer iNumNz = oEPM.NumGlobalNonzeros();

     for (integer iRow = 0; iRow <= iNumRows; ++iRow) {
          Ap[iRow] = rowptr[iRow] + offset;
     }
     
     for (integer iItem = 0; iItem < iNumNz; ++iItem) {
          Ax[iItem] = values[iItem];
          Ai[iItem] = colind[iItem] + offset;
     }

     return iNumNz;
}


int32_t EpetraSparseMatrixHandler::MakeCompressedRowForm(std::vector<doublereal>& Ax,
                                                         std::vector<int32_t>& Ai,
                                                         std::vector<int32_t>& Ap,
                                                         int offset) const
{
     FillComplete(); // Must be called before Nz()

     Ai.resize(Nz());
     Ax.resize(Nz());
     Ap.resize(iGetNumCols() + 1);

     return MakeCompressedRowForm(&Ax.front(), &Ai.front(), &Ap.front(), offset);
}

int64_t EpetraSparseMatrixHandler::MakeCompressedRowForm(std::vector<doublereal>& Ax,
                                                         std::vector<int64_t>& Ai,
                                                         std::vector<int64_t>& Ap,
                                                         int offset) const
{
     FillComplete(); // Must be called before Nz()

     Ai.resize(Nz());
     Ax.resize(Nz());
     Ap.resize(iGetNumCols() + 1);

     return MakeCompressedRowForm(&Ax.front(), &Ai.front(), &Ap.front(), offset);
}

int32_t
EpetraSparseMatrixHandler::MakeCompressedRowForm(doublereal *const Ax,
                                                 int32_t *const Ai,
                                                 int32_t *const Ap,
                                                 int offset) const
{
     return MakeCompressedRowFormTpl(Ax, Ai, Ap, offset);
}

int64_t
EpetraSparseMatrixHandler::MakeCompressedRowForm(doublereal *const Ax,
                                                 int64_t *const Ai,
                                                 int64_t *const Ap,
                                                 int offset) const
{
     return MakeCompressedRowFormTpl(Ax, Ai, Ap, offset);
}

#ifdef USE_SPARSE_AUTODIFF
bool EpetraSparseMatrixHandler::AddItem(integer iRow, const sp_grad::SpGradient& oItem)
{
#ifdef DEBUG
     IsValid();
#endif     
     ASSERT(iRow >= 1);
     ASSERT(iRow <= iGetNumRows());

     std::vector<doublereal> rgValues;
     std::vector<sp_grad::index_type> rgColIndex;

     rgValues.reserve(oItem.iGetSize());
     rgColIndex.reserve(oItem.iGetSize());
     
     for (const auto& oDer: oItem) {
          rgValues.push_back(oDer.dDer);
          rgColIndex.push_back(oDer.iDof);
     }

     integer ierr;
          
     if (oEPM.Filled()) {
          ierr = oEPM.SumIntoGlobalValues(iRow, oItem.iGetSize(), &rgValues.front(), &rgColIndex.front());
     } else {
          ierr = oEPM.InsertGlobalValues(iRow, oItem.iGetSize(), &rgValues.front(), &rgColIndex.front());
     }

     if (ierr < 0) {
          ASSERT(0);
          throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }
     
     return true;
}
#endif

VectorHandler& EpetraSparseMatrixHandler::GetCol(integer icol, VectorHandler& out) const
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

void EpetraSparseMatrixHandler::Scale(const std::vector<doublereal>& oRowScale, const std::vector<doublereal>& oColScale)
{
#ifdef DEBUG
     IsValid();
#endif

     integer* rowptr;
     integer* colind;
     doublereal* values;

     ExtractCrsDataPointers(rowptr, colind, values);
     
     CSCMatrixHandlerTpl<doublereal, integer, 0> A_T(values, colind, rowptr, NRows, Nz());
     
     A_T.Scale(oColScale, oRowScale); // Exchange rows by cols
     
#ifdef DEBUG
     IsValid();
#endif     
}

doublereal EpetraSparseMatrixHandler::Norm(Norm_t eNorm) const
{
#ifdef DEBUG
     IsValid();
#endif
     
     doublereal dNorm = 0.;
     switch (eNorm) {
     case NORM_1: {
          return oEPM.NormOne();
     } break;
     case NORM_INF: {
          return oEPM.NormInf();
     } break;
     default:
	  throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
     }

     return dNorm;
}

integer EpetraSparseMatrixHandler::Nz() const {
#ifdef DEBUG
     IsValid();
#endif
     
     ASSERT(oEPM.Filled());
     
     return oEPM.NumGlobalNonzeros();
}

void EpetraSparseMatrixHandler::EnumerateNz(const std::function<EnumerateNzCallback>& func) const
{
#ifdef DEBUG
     IsValid();
#endif
     FillComplete();
     
     integer* rowptr;
     integer* colind;
     doublereal* values;

     ASSERT(oEPM.Filled());
     ASSERT(oEPM.StorageOptimized());
     
     ExtractCrsDataPointers(rowptr, colind, values);

     const integer iNumRows = iGetNumRows();

     for (integer iRow = 0; iRow < iNumRows; ++iRow) {
          for (integer iItem = rowptr[iRow]; iItem < rowptr[iRow + 1]; ++iItem) {
               func(iRow + 1, colind[iItem] + 1, values[iItem]);
          }
     }
}

EpetraSparseMatrixHandler* EpetraSparseMatrixHandler::Copy() const
{
#ifdef DEBUG
     IsValid();
#endif
     
     EpetraSparseMatrixHandler* pMH = nullptr;

     SAFENEWWITHCONSTRUCTOR(pMH,
                            EpetraSparseMatrixHandler,
                            EpetraSparseMatrixHandler(iGetNumRows(),
                                                      iGetNumCols(),
                                                      iNumColsAlloc,
                                                      oEPM.Comm()));

     return pMH;
}
#endif
