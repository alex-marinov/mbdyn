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

#include "mbconfig.h"

#include <iomanip>

#include "sp_gradient_spmh.h"

SpGradientSparseMatrixHandler::SpGradientSparseMatrixHandler(const integer& iNumRows, const integer& iNumCols)
     :SparseMatrixHandler(iNumRows, iNumCols),
      NZ(0), oRows(iNumRows) {
}

SpGradientSparseMatrixHandler::~SpGradientSparseMatrixHandler()
{
}

#ifdef DEBUG
void SpGradientSparseMatrixHandler::IsValid() const
{
#ifdef USE_MULTITHREAD
     std::vector<bool> bValidated(oRows.size(), false);
     bool bNotValidated;
#endif

#ifdef USE_MULTITHREAD
     do {
	  bNotValidated = false;
#endif
	  for (size_t i = 0; i < oRows.size(); ++i) {
#ifdef USE_MULTITHREAD
	       if (bValidated[i]) {
		    continue;
	       }

	       if (std::atomic_exchange(&oRows[i].bLocked, true)) {
		    bNotValidated = true;
		    continue;
	       }
#endif

	       SP_GRAD_ASSERT(oRows[i].bValid());

	       for (const auto& oDer: oRows[i]) {
		    SP_GRAD_ASSERT(oDer.iDof >= 1);
		    SP_GRAD_ASSERT(oDer.iDof <= NCols);
	       }

#ifdef USE_MULTITHREAD
	       bValidated[i] = true;
	       std::atomic_exchange(&oRows[i].bLocked, false);
#endif
	  }
#ifdef USE_MULTITHREAD
     } while (bNotValidated);
#endif
}
#endif

void SpGradientSparseMatrixHandler::Resize(integer, integer)
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

void SpGradientSparseMatrixHandler::ResizeReset(integer iNumRows, integer iNumCols)
{
     oRows.resize(iNumRows);
     Reset();
     
     NRows = iNumRows;
     NCols = iNumCols;
}

void SpGradientSparseMatrixHandler::Reset()
{
     // Not thread safe
     for (auto& oRow: oRows) {
	  oRow.ResetNumeric();
     }
}

void
SpGradientSparseMatrixHandler::IncCoef(integer iRow, integer iCol, const doublereal& dCoef)
{
     // Not recommended for general use
     // Needed for DataManager::InitialJointAssembly to make it independent from the bUseAutoDiff() flag
     sp_grad::SpGradient g;

     g.Reset(0., iCol, dCoef);

     AddItem(iRow, g);
}

MatrixHandler&
SpGradientSparseMatrixHandler::MatMatMul_base(void (MatrixHandler::*op)(integer iRow, integer iCol,
					 const doublereal& dCoef),
	       MatrixHandler& out, const MatrixHandler& in) const
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

MatrixHandler&
SpGradientSparseMatrixHandler::MatTMatMul_base(void (MatrixHandler::*op)(integer iRow, integer iCol,
					  const doublereal& dCoef),
		MatrixHandler& out, const MatrixHandler& in) const
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

VectorHandler&
SpGradientSparseMatrixHandler::MatVecMul_base(
     void (VectorHandler::*op)(integer iRow, const doublereal& dCoef),
     VectorHandler& out, const VectorHandler& in) const
{
     integer nr = iGetNumRows();
     integer nc = iGetNumCols();

     if (out.iGetSize() != nr || in.iGetSize() != nc) {
	  silent_cerr("SpGradientSparseMatrixHandler::MatVecMul_base(): size mismatch "
		      "out(" << out.iGetSize() << ", 1) "
		      "= this(" << nr << ", " << nc << ") "
		      "* in(" << in.iGetSize() << ", 1)" << std::endl);
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     for (integer r = 1; r <= nr; ++r) {
	  doublereal d = 0.;

	  for (const auto& oDer: oRows[r - 1]) {
	       d += oDer.dDer * in(oDer.iDof);
	  }

	  (out.*op)(r, d);
     }

     return out;

}

VectorHandler&
SpGradientSparseMatrixHandler::MatTVecMul_base(
     void (VectorHandler::*op)(integer iRow, const doublereal& dCoef),
     VectorHandler& out, const VectorHandler& in) const
{
     integer nr = iGetNumRows();
     integer nc = iGetNumCols();

     if (out.iGetSize() != nc || in.iGetSize() != nr) {
	  silent_cerr("SpGradientSparseMatrixHandler::MatVecMul_base(): size mismatch "
		      "out(" << out.iGetSize() << ", 1) "
		      "= this(" << nr << ", " << nc << ")^T "
		      "* in(" << in.iGetSize() << ", 1)" << std::endl);
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     for (integer r = 1; r <= nr; ++r) {
	  for (const auto& oDer: oRows[r - 1]) {
	       (out.*op)(oDer.iDof, oDer.dDer * in(r));
	  }
     }

     return out;
}

const doublereal&
SpGradientSparseMatrixHandler::operator()(integer iRow, integer iCol) const
{
     for (const auto& oDer: oRows[iRow - 1]) {
          if (oDer.iDof == iCol) {
               return oDer.dDer;
          }
     }

     static constexpr doublereal dZero = 0.;
     
     return dZero;
}

doublereal&
SpGradientSparseMatrixHandler::operator()(integer iRow, integer iCol)
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

template <typename Operation>
bool SpGradientSparseMatrixHandler::ItemOperation(integer iRow, const Operation& oper, const sp_grad::SpGradient& oItem)
{
     SP_GRAD_ASSERT(iRow >= 1);
     SP_GRAD_ASSERT(iRow <= iGetNumRows());
     SP_GRAD_ASSERT(oRows.size() == static_cast<size_t>(iGetNumRows()));

#ifdef SP_GRAD_DEBUG
     for (const auto& r: oItem) {
	  SP_GRAD_ASSERT(r.iDof >= 1);
	  SP_GRAD_ASSERT(r.iDof <= iGetNumCols());
     }
#endif

     SparseRow& oRow = oRows[iRow - 1];

#ifdef USE_MULTITHREAD
     if (std::atomic_exchange(&oRow.bLocked, true)) {
	  return false;
     }
#endif
     const integer iSizeRowPrev = oRow.iGetSize();

     oper(oRow, oItem);

     SP_GRAD_ASSERT(oRow.bIsUnique());

     NZ += oRow.iGetSize() - iSizeRowPrev;

#ifdef SP_GRAD_DEBUG
     for (const auto& r: oRows[iRow - 1]) {
	  SP_GRAD_ASSERT(r.iDof >= 1);
	  SP_GRAD_ASSERT(r.iDof <= iGetNumCols());
     }
#endif

#ifdef USE_MULTITHREAD
     std::atomic_exchange(&oRow.bLocked, false);
#endif

     return true;     
}

void SpGradientSparseMatrixHandler::AddCoef(sp_grad::SpGradient& a, const sp_grad::SpGradient& b)
{
     a += sp_grad::EvalUnique(b);
}

void SpGradientSparseMatrixHandler::SubCoef(sp_grad::SpGradient& a, const sp_grad::SpGradient& b)
{
     a -= sp_grad::EvalUnique(b);
}

template <typename idx_type>
idx_type SpGradientSparseMatrixHandler::MakeCompressedColumnFormTpl(doublereal *const Ax,
								    idx_type *const Ai,
								    idx_type *const Ap,
								    int offset) const
{
    std::vector<idx_type> oColSize(NCols, 0);

     for (const auto& oRow: oRows) {
	  SP_GRAD_ASSERT(oRow.bValid());
	  SP_GRAD_ASSERT(oRow.iGetSize() > 0);
	  SP_GRAD_ASSERT(oRow.bIsUnique());

	  for (const auto& oDer: oRow) {
	       SP_GRAD_ASSERT(oDer.iDof >= 1);
	       SP_GRAD_ASSERT(oDer.iDof <= NCols);

	       ++oColSize[oDer.iDof - 1];
	  }
     }

     idx_type iPtr = offset;

     for (integer iCol = 0; iCol < NCols; ++iCol) {
	  Ap[iCol] = iPtr;

	  SP_GRAD_ASSERT(oColSize[iCol] >= 0);
	  SP_GRAD_ASSERT(oColSize[iCol] <= NRows);

	  iPtr += oColSize[iCol];
     }

     Ap[NCols] = iPtr;

     SP_GRAD_ASSERT(iPtr - offset == Nz());

     std::fill(oColSize.begin(), oColSize.end(), 0);

#ifdef SP_GRAD_DEBUG
     constexpr idx_type iInvalidIndex = -1;
     constexpr doublereal dInvalidValue = std::numeric_limits<doublereal>::infinity();

     std::fill(Ai, Ai + Nz(), iInvalidIndex);
     std::fill(Ax, Ax + Nz(), dInvalidValue);
#endif

     for (integer iRow = 0; iRow < NRows; ++iRow) {
	  for (const auto& oDer: oRows[iRow]) {
	       const integer iColIdx = oDer.iDof - 1;

	       SP_GRAD_ASSERT(iColIdx >= 0);
	       SP_GRAD_ASSERT(iColIdx < NCols);

	       const idx_type iSlot = Ap[iColIdx] - offset + oColSize[iColIdx]++;

	       SP_GRAD_ASSERT(iSlot >= 0);
	       SP_GRAD_ASSERT(iSlot < Nz());
	       SP_GRAD_ASSERT(Ai[iSlot] == iInvalidIndex);
	       SP_GRAD_ASSERT(Ax[iSlot] == dInvalidValue);

	       Ai[iSlot] = iRow + offset;
	       Ax[iSlot] = oDer.dDer;
	  }
     }

#ifdef SP_GRAD_DEBUG
     idx_type iSizeTot = 0;

     for (idx_type iColSize: oColSize) {
	  iSizeTot += iColSize;
     }

     SP_GRAD_ASSERT(iSizeTot == Nz());

     for (integer iCol = 0; iCol < NCols; ++iCol) {
	  idx_type iRowPrev = -1;
	  for (idx_type iSlot = Ap[iCol] - offset; iSlot < Ap[iCol + 1] - offset; ++iSlot) {
	       idx_type iRow = Ai[iSlot] - offset;

	       SP_GRAD_ASSERT(iRow > iRowPrev);
	       SP_GRAD_ASSERT(iRow >= 0);
	       SP_GRAD_ASSERT(iRow < NRows);
	       SP_GRAD_ASSERT(std::isfinite(Ax[iSlot]));

	       for (const auto& oDer: oRows[iRow]) {
		    if (oDer.iDof - 1 == iCol) {
			 SP_GRAD_ASSERT(Ax[iSlot] == oDer.dDer);
		    }
	       }

	       iRowPrev = iRow;
	  }
     }
#endif

     return iPtr - offset;
}

int32_t SpGradientSparseMatrixHandler::MakeCompressedColumnForm(doublereal *const Ax,
								int32_t *const Ai,
								int32_t *const Ap,
								int offset) const
{
     return MakeCompressedColumnFormTpl(Ax, Ai, Ap, offset);
}

int64_t SpGradientSparseMatrixHandler::MakeCompressedColumnForm(doublereal *const Ax,
								int64_t *const Ai,
								int64_t *const Ap,
								int offset) const
{
     return MakeCompressedColumnFormTpl(Ax, Ai, Ap, offset);
}

bool SpGradientSparseMatrixHandler::AddItem(integer iRow, const sp_grad::SpGradient& oItem)
{
     return ItemOperation(iRow, AddCoef, oItem);
}

bool SpGradientSparseMatrixHandler::SubItem(integer iRow, const sp_grad::SpGradient& oItem)
{
     return ItemOperation(iRow, SubCoef, oItem);
}

int32_t SpGradientSparseMatrixHandler::MakeIndexForm(doublereal *const Ax,
						     int32_t *const Arow, int32_t *const Acol,
						     int32_t *const Ap,
						     int offset) const
{
     return MakeIndexFormTpl(Ax, Arow, Acol, Ap, offset);
}


int64_t SpGradientSparseMatrixHandler::MakeIndexForm(doublereal *const Ax,
						     int64_t *const Arow, int64_t *const Acol,
						     int64_t *const Ap,
						     int offset) const
{
     return MakeIndexFormTpl(Ax, Arow, Acol, Ap, offset);
}

template <typename idx_type>
idx_type SpGradientSparseMatrixHandler::MakeIndexFormTpl(doublereal *const Ax,
							 idx_type *const Arow, idx_type *const Acol,
							 idx_type *const Ap,
							 int offset) const
{
    std::vector<idx_type> oColSize(NCols, 0);

     for (const auto& oRow: oRows) {
	  SP_GRAD_ASSERT(oRow.bValid());
	  SP_GRAD_ASSERT(oRow.iGetSize() > 0);
	  SP_GRAD_ASSERT(oRow.bIsUnique());

	  for (const auto& oDer: oRow) {
	       SP_GRAD_ASSERT(oDer.iDof >= 1);
	       SP_GRAD_ASSERT(oDer.iDof <= NCols);

	       ++oColSize[oDer.iDof - 1];
	  }
     }

     idx_type iPtr = offset;

     for (integer iCol = 0; iCol < NCols; ++iCol) {
	  Ap[iCol] = iPtr;

	  SP_GRAD_ASSERT(oColSize[iCol] >= 0);
	  SP_GRAD_ASSERT(oColSize[iCol] <= NRows);

	  iPtr += oColSize[iCol];
     }

     Ap[NCols] = iPtr;

     SP_GRAD_ASSERT(iPtr - offset == Nz());

     std::fill(oColSize.begin(), oColSize.end(), 0);

#ifdef SP_GRAD_DEBUG
     constexpr idx_type iInvalidIndex = -1;
     constexpr doublereal dInvalidValue = std::numeric_limits<doublereal>::infinity();

     std::fill(Arow, Arow + Nz(), iInvalidIndex);
     std::fill(Acol, Acol + Nz(), iInvalidIndex);
     std::fill(Ax, Ax + Nz(), dInvalidValue);
#endif

     for (integer iRow = 0; iRow < NRows; ++iRow) {
	  for (const auto& oDer: oRows[iRow]) {
	       const integer iColIdx = oDer.iDof - 1;

	       SP_GRAD_ASSERT(iColIdx >= 0);
	       SP_GRAD_ASSERT(iColIdx < NCols);

	       const idx_type iSlot = Ap[iColIdx] - offset + oColSize[iColIdx]++;

	       SP_GRAD_ASSERT(iSlot >= 0);
	       SP_GRAD_ASSERT(iSlot < Nz());
	       SP_GRAD_ASSERT(Arow[iSlot] == iInvalidIndex);
	       SP_GRAD_ASSERT(Acol[iSlot] == iInvalidIndex);
	       SP_GRAD_ASSERT(Ax[iSlot] == dInvalidValue);

	       Arow[iSlot] = iRow + offset;
	       Acol[iSlot] = iColIdx + offset;
	       Ax[iSlot] = oDer.dDer;
	  }
     }

     return iPtr - offset;
}


VectorHandler& SpGradientSparseMatrixHandler::GetCol(integer icol,
						     VectorHandler& out) const
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

void SpGradientSparseMatrixHandler::Scale(const std::vector<doublereal>& oRowScale, const std::vector<doublereal>& oColScale)
{
     const bool bScaleRows = !oRowScale.empty();
     const bool bScaleCols = !oColScale.empty();

     ASSERT(!bScaleRows || oRowScale.size() == static_cast<size_t>(iGetNumRows()));
     ASSERT(!bScaleCols || oColScale.size() == static_cast<size_t>(iGetNumCols()));
     ASSERT(bScaleRows || bScaleCols);

     if (bScaleRows && bScaleCols) {
	  for (integer iRow = 0; iRow < NRows; ++iRow) {
	       oRows[iRow].Scale(oRowScale[iRow], oColScale);
	  }
     } else if (bScaleCols) {
	  for (integer iRow = 0; iRow < NRows; ++iRow) {
	       oRows[iRow].Scale(1., oColScale);
	  }
     } else if (bScaleRows) {
	  for (integer iRow = 0; iRow < NRows; ++iRow) {
	       oRows[iRow] *= oRowScale[iRow];
	  }
     }
}

doublereal SpGradientSparseMatrixHandler::Norm(Norm_t eNorm) const
{
     doublereal dNorm = 0.;
     switch (eNorm) {
     case NORM_1: {
	  std::vector<doublereal> csum(iGetNumCols(), 0.);

	  for (const auto& oItem: *this) {
	       ASSERT(oItem.iCol >= 0);
	       ASSERT(oItem.iCol < iGetNumCols());
	       csum[oItem.iCol] += std::abs(oItem.dCoef);
	  }

	  for (doublereal d: csum) {
	       dNorm = std::max(d, dNorm);
	  }
     } break;
     case NORM_INF: {
	  std::vector<doublereal> rsum(iGetNumRows(), 0.);

	  for (const auto& oItem: *this) {
	       ASSERT(oItem.iRow >= 0);
	       ASSERT(oItem.iRow < iGetNumRows());
	       rsum[oItem.iRow] += std::abs(oItem.dCoef);
	  }

	  for (doublereal d: rsum) {
	       dNorm = std::max(d, dNorm);
	  }
     } break;
     default:
	  throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
     }

     return dNorm;
}

integer SpGradientSparseMatrixHandler::Nz() const {
     return NZ;
}

void SpGradientSparseMatrixHandler::EnumerateNz(const std::function<EnumerateNzCallback>& func) const
{
     for (integer iRow = 1; iRow <= iGetNumRows(); ++iRow) {
	  for (const auto& oDer: oRows[iRow - 1]) {
	       func(iRow, oDer.iDof, oDer.dDer);
	  }
     }
}

SpGradientSparseMatrixHandler* SpGradientSparseMatrixHandler::Copy() const
{
     SpGradientSparseMatrixHandler* pMH = nullptr;

     SAFENEWWITHCONSTRUCTOR(pMH, SpGradientSparseMatrixHandler, SpGradientSparseMatrixHandler(iGetNumRows(), iGetNumCols()));

     return pMH;
}

#ifdef USE_MULTITHREAD
SpGradientSparseMatrixWrapper::SpGradientSparseMatrixWrapper(SpGradientSparseMatrixHandler* pMH)
     :pMH(pMH)
{
}

SpGradientSparseMatrixWrapper::~SpGradientSparseMatrixWrapper()
{
}

#ifdef DEBUG
void SpGradientSparseMatrixWrapper::IsValid() const
{
     pMH->IsValid();
}
#endif

void SpGradientSparseMatrixWrapper::SetMatrixHandler(SpGradientSparseMatrixHandler* pMatHdl)
{
     pMH = pMatHdl;
}

integer SpGradientSparseMatrixWrapper::iGetNumRows() const
{
     return pMH->SpGradientSparseMatrixHandler::iGetNumRows();
}

integer SpGradientSparseMatrixWrapper::iGetNumCols() const
{
     return pMH->SpGradientSparseMatrixHandler::iGetNumCols();
}

void SpGradientSparseMatrixWrapper::Resize(integer, integer)
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

void SpGradientSparseMatrixWrapper::ResizeReset(integer, integer)
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

void SpGradientSparseMatrixWrapper::Reset()
{
     // FIXME: pMH->Reset() must be called in advance by the main thread before starting multithreaded assembly.
}

const doublereal&
SpGradientSparseMatrixWrapper::operator()(integer iRow, integer iCol) const
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

doublereal&
SpGradientSparseMatrixWrapper::operator()(integer iRow, integer iCol)
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

void SpGradientSparseMatrixWrapper::Scale(const std::vector<doublereal>& oRowScale, const std::vector<doublereal>& oColScale)
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

bool SpGradientSparseMatrixWrapper::AddItem(integer iRow, const sp_grad::SpGradient& oItem)
{
     return pMH->SpGradientSparseMatrixHandler::AddItem(iRow, oItem);
}

void SpGradientSparseMatrixWrapper::EnumerateNz(const std::function<EnumerateNzCallback>& func) const
{
     pMH->EnumerateNz(func);
}

doublereal SpGradientSparseMatrixWrapper::Norm(Norm_t eNorm) const
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

SpGradientSparseMatrixWrapper* SpGradientSparseMatrixWrapper::Copy() const
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

MatrixHandler&
SpGradientSparseMatrixWrapper::MatMatMul_base(void (MatrixHandler::*op)(integer iRow, integer iCol,
									const doublereal& dCoef),
					      MatrixHandler& out, const MatrixHandler& in) const
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

MatrixHandler&
SpGradientSparseMatrixWrapper::MatTMatMul_base(void (MatrixHandler::*op)(integer iRow, integer iCol,
									 const doublereal& dCoef),
					       MatrixHandler& out, const MatrixHandler& in) const
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

VectorHandler&
SpGradientSparseMatrixWrapper::MatVecMul_base(
     void (VectorHandler::*op)(integer iRow, const doublereal& dCoef),
     VectorHandler& out, const VectorHandler& in) const
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

VectorHandler&
SpGradientSparseMatrixWrapper::MatTVecMul_base(
     void (VectorHandler::*op)(integer iRow, const doublereal& dCoef),
     VectorHandler& out, const VectorHandler& in) const
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}
#endif
