/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2020
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

/*
 AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
        Copyright (C) 2020(-2020) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#ifndef __SP_MATRIX_BASE_H__INCLUDED__
#define __SP_MATRIX_BASE_H__INCLUDED__

#include <algorithm>
#include <cstdlib>
#include <new>

#include "sp_matrix_base_fwd.h"
#include "sp_gradient_func.h"
#include "sp_gradient.h"

namespace sp_grad {
     index_type SpMatrixBaseData::iGetNumRows() const noexcept {
	  SP_GRAD_ASSERT(iNumRows >= 0);
	  SP_GRAD_ASSERT(iNumCols >= 0);
	  return iNumRows;
     }

     index_type SpMatrixBaseData::iGetNumCols() const noexcept {
	  SP_GRAD_ASSERT(iNumRows >= 0);
	  SP_GRAD_ASSERT(iNumCols >= 0);

	  return iNumCols;
     }

     index_type SpMatrixBaseData::iGetNumElem() const noexcept {
	  SP_GRAD_ASSERT(iNumRows >= 0);
	  SP_GRAD_ASSERT(iNumCols >= 0);

	  return iNumRows * iNumCols;
     }

     index_type SpMatrixBaseData::iGetRefCnt() const noexcept {
	  SP_GRAD_ASSERT(iRefCnt >= 0);

	  return iRefCnt;
     }

     index_type SpMatrixBaseData::iGetMaxDeriv() const noexcept {
	  SP_GRAD_ASSERT(iNumDeriv >= 0);

	  return iNumDeriv;
     }

     SpMatrixBaseData* SpMatrixBaseData::pGetNullData() noexcept {
	  SP_GRAD_ASSERT(oNullData.iRefCnt >= 0);
	  SP_GRAD_ASSERT(oNullData.iNumRows == 0);
	  SP_GRAD_ASSERT(oNullData.iNumCols == 0);

	  return &oNullData;
     }

     SpMatrixBaseData::SpMatrixBaseData(index_type iNumRows, index_type iNumCols, index_type iRefCnt, index_type iNumDeriv) noexcept
	  :iNumRows(iNumRows),
	   iNumCols(iNumCols),
	   iRefCnt(iRefCnt),
	   iNumDeriv(iNumDeriv) {
     }

     namespace util {
	  template <typename ValueType>
	  void SpMatrixDataTraits<ValueType>::Construct(SpMatrixData<ValueType>& oData, index_type iNumDeriv, void* pExtraMem) {
	       const index_type iNumElem = oData.iGetNumElem();

	       for (index_type i = 1; i <= iNumElem; ++i) {
		    new(oData.pGetData(i)) ValueType{}; // Enforce default constructor also for doublereal!
	       }
	  }

	  template <typename ValueType>
	  size_t constexpr SpMatrixDataTraits<ValueType>::uOffsetDeriv(size_t uSizeGrad) noexcept {
	       return 0;
	  }

	  template <typename ValueType>
	  constexpr size_t SpMatrixDataTraits<ValueType>::uSizeDeriv(index_type iNumDeriv, index_type iNumItems) noexcept {
	       return 0;
	  }

	  template <typename ValueType>
	  constexpr size_t SpMatrixDataTraits<ValueType>::uAlign() noexcept {
	       return 0;
	  }

	  void SpMatrixDataTraits<SpGradient>::Construct(SpMatrixData<SpGradient>& oData, index_type iNumDeriv, void* pExtraMem) {
	       static_assert(alignof(SpDerivData) == alignof(SpDerivRec));

	       SP_GRAD_ASSERT(reinterpret_cast<size_t>(pExtraMem) % alignof(SpDerivData) == 0);

	       const size_t uBytesDeriv = sizeof(SpDerivData) + iNumDeriv * sizeof(SpDerivRec);
	       const index_type iNumElem = oData.iGetNumElem();

	       for (index_type i = 1; i <= iNumElem; ++i) {
		    size_t uBytesOffset = (i - 1) * uBytesDeriv;
		    auto pDerivData = reinterpret_cast<SpDerivData*>(reinterpret_cast<char*>(pExtraMem) + uBytesOffset);

		    SP_GRAD_ASSERT(reinterpret_cast<size_t>(pDerivData) % alignof(SpDerivData) == 0);

		    new(pDerivData) SpDerivData(0., iNumDeriv, 0, true, 1, &oData);
		    new(oData.pGetData(i)) SpGradient(pDerivData);
	       }
	  }

	  constexpr size_t SpMatrixDataTraits<SpGradient>::uOffsetDeriv(size_t uSizeGrad) noexcept {
	       return uSizeGrad % alignof(SpDerivData)
		    ? alignof(SpDerivData) - uSizeGrad % alignof(SpDerivData)
		    : 0;
	  }

	  constexpr size_t SpMatrixDataTraits<SpGradient>::uSizeDeriv(index_type iNumDeriv, index_type iNumItems) noexcept {
	       return (sizeof(SpDerivData) + sizeof(SpDerivRec) * iNumDeriv) * iNumItems;
	  }

	  constexpr size_t SpMatrixDataTraits<SpGradient>::uAlign() noexcept {
	       return alignof(SpDerivData);
	  }
     }

     template <typename ValueType>
     SpMatrixData<ValueType>::SpMatrixData(index_type iNumRows,
					   index_type iNumCols,
					   index_type iRefCnt,
					   index_type iNumDeriv,
					   void* pExtraMem)
	  :SpMatrixBaseData(iNumRows, iNumCols, iRefCnt, iNumDeriv) {

	  util::SpMatrixDataTraits<ValueType>::Construct(*this, iNumDeriv, pExtraMem);
     }

     template <typename ValueType>
     SpMatrixData<ValueType>::~SpMatrixData() {
	  for (auto& g: *this) {
	       g.~ValueType();
	  }
     }

     template <typename ValueType>
     const ValueType* SpMatrixData<ValueType>::pGetData(index_type iRow, index_type iCol) const noexcept {
	  SP_GRAD_ASSERT(iRow >= 1);
	  SP_GRAD_ASSERT(iRow <= iGetNumRows());
	  SP_GRAD_ASSERT(iCol >= 1);
	  SP_GRAD_ASSERT(iCol <= iGetNumCols());

	  return pGetData() + (iCol - 1) * iGetNumRows() + iRow - 1;
     }

     template <typename ValueType>
     ValueType* SpMatrixData<ValueType>::pGetData(index_type iRow, index_type iCol) noexcept {
	  SP_GRAD_ASSERT(iRow >= 1);
	  SP_GRAD_ASSERT(iRow <= iGetNumRows());
	  SP_GRAD_ASSERT(iCol >= 1);
	  SP_GRAD_ASSERT(iCol <= iGetNumCols());

	  return pGetData() + (iCol - 1) * iGetNumRows() + iRow - 1;
     }

     template <typename ValueType>
     ValueType* SpMatrixData<ValueType>::pGetData() noexcept {
	  auto pVal = reinterpret_cast<ValueType*>(this + 1);

	  SP_GRAD_ASSERT(reinterpret_cast<size_t>(pVal) % alignof(ValueType) == 0);

	  return pVal;
     }

     template <typename ValueType>
     ValueType* SpMatrixData<ValueType>::pGetData(index_type iRow) noexcept {
	  SP_GRAD_ASSERT(iRow >= 1);
	  SP_GRAD_ASSERT(iRow <= iNumRows * iNumCols);

	  return pGetData() + iRow - 1;
     }

     template <typename ValueType>
     const ValueType* SpMatrixData<ValueType>::pGetData() const noexcept {
	  auto pVal = reinterpret_cast<const ValueType*>(this + 1);

	  SP_GRAD_ASSERT(reinterpret_cast<size_t>(pVal) % alignof(ValueType) == 0);

	  return pVal;
     }

     template <typename ValueType>
     const ValueType* SpMatrixData<ValueType>::pGetData(index_type iRow) const noexcept {
	  SP_GRAD_ASSERT(iRow >= 1);
	  SP_GRAD_ASSERT(iRow <= iNumRows * iNumCols);

	  return pGetData() + iRow - 1;
     }

     template <typename ValueType>
     void SpMatrixData<ValueType>::IncRef() noexcept {
	  SP_GRAD_ASSERT(iRefCnt >= 0);
	  ++iRefCnt;
     }

     template <typename ValueType>
     void SpMatrixData<ValueType>::DecRef() {
	  SP_GRAD_ASSERT(iRefCnt >= 1);

	  if (--iRefCnt == 0) {
	       if (this != &oNullData) {	       
		    this->~SpMatrixData();

#if defined(HAVE_ALIGNED_MALLOC)
		    _aligned_free(this);
#else
		    free(this);
#endif
	       }
	  }
     }

     template <typename ValueType>
     void SpMatrixData<ValueType>::Attach(ValueType* pData) noexcept {
	  if (!bIsOwnerOf(pData)) {
	       IncRef();
	  }
     }

     template <typename ValueType>
     void SpMatrixData<ValueType>::Detach(ValueType* pData) {
	  if (!bIsOwnerOf(pData)) {
	       DecRef();
	  }
     }

     template <typename ValueType>
     ValueType* SpMatrixData<ValueType>::begin() noexcept {
	  return pGetData();
     }

     template <typename ValueType>
     ValueType* SpMatrixData<ValueType>::end() noexcept {
	  return pGetData() + iGetNumElem();
     }

     template <typename ValueType>
     const ValueType* SpMatrixData<ValueType>::begin() const noexcept {
	  return pGetData();
     }

     template <typename ValueType>
     const ValueType* SpMatrixData<ValueType>::end() const noexcept {
	  return pGetData() + iGetNumElem();
     }

     template <typename ValueType>
     constexpr bool SpMatrixData<ValueType>::bIsOwnerOf(const ValueType* pData) const noexcept {
	  return pData >= begin() && pData < end();
     }

     namespace util {
	  template <index_type iSizeStatic>
	  struct MatrixDataSizeHelper {
	       static_assert(iSizeStatic > 0);

	       static constexpr index_type iGetSizeStatic(index_type iSize) {
		    return iSizeStatic;
	       }
	  };

	  template <>
	  struct MatrixDataSizeHelper<SpMatrixSize::DYNAMIC> {
	       static_assert(SpMatrixSize::DYNAMIC < 0);

	       static constexpr index_type iGetSizeStatic(index_type iSize) {
		    return iSize;
	       }
	  };
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     SpMatrixDataCTD<ValueType, NumRows, NumCols>::SpMatrixDataCTD(index_type iNumRows,
								   index_type iNumCols,
								   index_type iRefCnt,
								   index_type iNumDeriv,
								   void* pExtraMem)
	  :SpMatrixData<ValueType>(iNumRows,
				   iNumCols,
				   iRefCnt,
				   iNumDeriv,
				   pExtraMem) {
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     SpMatrixDataCTD<ValueType, NumRows, NumCols>*
     SpMatrixDataCTD<ValueType, NumRows, NumCols>::pAllocate(index_type iNumRows, index_type iNumCols, index_type iNumDeriv) {
	  SP_GRAD_ASSERT(NumRows == SpMatrixSize::DYNAMIC || NumRows == iNumRows);
	  SP_GRAD_ASSERT(NumCols == SpMatrixSize::DYNAMIC || NumCols == iNumCols);

	  iNumRows = util::MatrixDataSizeHelper<iNumRowsStatic>::iGetSizeStatic(iNumRows);
	  iNumCols = util::MatrixDataSizeHelper<iNumColsStatic>::iGetSizeStatic(iNumCols);
	  const index_type iNumItems = iNumRows * iNumCols;
	  const size_t uSizeGrad = sizeof(SpMatrixDataCTD) + sizeof(ValueType) * iNumItems;
	  const size_t uOffsetDeriv = util::SpMatrixDataTraits<ValueType>::uOffsetDeriv(uSizeGrad);
	  const size_t uSizeDeriv = util::SpMatrixDataTraits<ValueType>::uSizeDeriv(iNumDeriv, iNumItems);
	  const size_t uSize = uSizeGrad + uOffsetDeriv + uSizeDeriv;
	  constexpr size_t uAlignVal = alignof(ValueType) > alignof(SpMatrixDataCTD)
	       ? alignof(ValueType) : alignof(SpMatrixDataCTD);
	  constexpr size_t uAlignDer = util::SpMatrixDataTraits<ValueType>::uAlign();
	  constexpr size_t uAlign = uAlignVal > uAlignDer ? uAlignVal : uAlignDer;

#if defined(HAVE_POSIX_MEMALIGN)
	  char* pMem;
	  if (0 != posix_memalign(reinterpret_cast<void**>(&pMem), uAlign, uSize)) {
	       pMem = nullptr;
	  }
#elif defined(HAVE_MEMALIGN)
	  char* pMem = reinterpret_cast<char*>(memalign(uAlign, uSize));	  
#elif defined(HAVE_ALIGNED_MALLOC)
	  char* pMem = reinterpret_cast<char*>(_aligned_malloc(uSize, uAlign));
#elif defined(HAVE_ALIGNED_ALLOC)
	  char* pMem = reinterpret_cast<char*>(aligned_alloc(uAlign, uSize));
#else
	  char* pMem = malloc(uSize);
#endif
	  
	  if (!pMem) {
	       throw std::bad_alloc();
	  }

	  void* pExtraMem = pMem + uSizeGrad + uOffsetDeriv;

	  return new(pMem) SpMatrixDataCTD(iNumRows, iNumCols, 0, iNumDeriv, pExtraMem);
     }

     namespace util {
	  template <typename MatEvalType, typename Expr, typename ValueType, index_type NumRows, index_type NumCols>
	  void MatEvalHelperCompr<SpGradCommon::ExprEvalUncompressed>::ElemEval(const Expr& oExpr, SpMatrixBase<ValueType, NumRows, NumCols>& A) {
	       static_assert(eComprFlag == SpGradCommon::ExprEvalUncompressed);

	       oExpr.template ElemEvalUncompr<MatEvalType>(A);
	  }

	  template <typename MatEvalType, typename Expr, typename ValueType, index_type NumRows, index_type NumCols>
	  void MatEvalHelperCompr<SpGradCommon::ExprEvalCompressed>::ElemEval(const Expr& oExpr, SpMatrixBase<ValueType, NumRows, NumCols>& A) {
	       static_assert(eComprFlag == SpGradCommon::ExprEvalCompressed);

	       oExpr.template ElemEvalCompr<MatEvalType>(A);
	  }

	  template <typename ValueType, index_type NumRows, index_type NumCols>
	  void MatEvalHelperTransp<MatTranspEvalFlag::DIRECT>::ResizeReset(SpMatrixBase<ValueType, NumRows, NumCols>& A,
									   index_type iNumRows,
									   index_type iNumCols,
									   index_type iMaxSize) {
	       A.ResizeReset(iNumRows, iNumCols, iMaxSize);
	  }

	  template <typename ValueType, index_type NumRows, index_type NumCols>
	  ValueType& MatEvalHelperTransp<MatTranspEvalFlag::DIRECT>::GetElem(SpMatrixBase<ValueType, NumRows, NumCols>& A, index_type i, index_type j) {
	       return A.GetElem(i, j);
	  }

	  template <typename ValueType, index_type NumRows, index_type NumCols>
	  void MatEvalHelperTransp<MatTranspEvalFlag::TRANSPOSED>::ResizeReset(SpMatrixBase<ValueType, NumRows, NumCols>& A,
									       index_type iNumRows,
									       index_type iNumCols,
									       index_type iMaxSize) {
	       A.ResizeReset(iNumCols, iNumRows, iMaxSize);
	  }

	  template <typename ValueType, index_type NumRows, index_type NumCols>
	  ValueType& MatEvalHelperTransp<MatTranspEvalFlag::TRANSPOSED>::GetElem(SpMatrixBase<ValueType, NumRows, NumCols>& A, index_type i, index_type j) {
	       return A.GetElem(j, i);
	  }

	  template <typename ValueType, SpGradCommon::ExprEvalFlags eCompr>
	  struct ComprEvalHelper {
	       static constexpr SpGradCommon::ExprEvalFlags eComprFlag = eCompr;
	  };

	  template <SpGradCommon::ExprEvalFlags eCompr>
	  struct ComprEvalHelper<doublereal, eCompr> {
	       // It makes no sense to allocate a class SpGradExprDofMap for a doublereal
	       static const SpGradCommon::ExprEvalFlags eComprFlag = SpGradCommon::ExprEvalUncompressed;
	  };
     };

     template <typename ValueType, typename DERIVED>
     template <typename MatEvalType, index_type NumRows, index_type NumCols>
     void SpMatElemExprBase<ValueType, DERIVED>::ElemEval(SpMatrixBase<ValueType, NumRows, NumCols>& A) const {
	  util::MatEvalHelperCompr<util::ComprEvalHelper<ValueType, eComprFlag>::eComprFlag>::template ElemEval<MatEvalType>(*this, A);
     }

     template <typename ValueType, typename DERIVED>
     template <typename MatEvalType, index_type NumRows, index_type NumCols>
     void SpMatElemExprBase<ValueType, DERIVED>::ElemEvalUncompr(SpMatrixBase<ValueType, NumRows, NumCols>& A) const {
	  const index_type iNumRows = iGetNumRows();
	  const index_type iNumCols = iGetNumCols();
	  const index_type iMaxSize = iGetMaxSize();

	  constexpr bool bTransposed = MatEvalType::eTranspFlag == util::MatTranspEvalFlag::TRANSPOSED;

	  static_assert(NumRows == (!bTransposed ? iNumRowsStatic : iNumColsStatic), "Number of rows does not match");
	  static_assert(NumCols == (!bTransposed ? iNumColsStatic : iNumRowsStatic), "Number of columns does not match");

	  SpMatrixBase<ValueType, NumRows, NumCols> Atmp;

	  if (!bHaveRefTo(A)) {
	       Atmp = std::move(A);
	  }

	  MatEvalType::ResizeReset(Atmp, iNumRows, iNumCols, iMaxSize);

	  SP_GRAD_ASSERT(Atmp.iGetNumRows() == (!bTransposed ? iNumRows : iNumCols));
	  SP_GRAD_ASSERT(Atmp.iGetNumCols() == (!bTransposed ? iNumCols : iNumRows));

	  for (index_type j = 1; j <= iNumCols; ++j) {
	       for (index_type i = 1; i <= iNumRows; ++i) {
		    ValueType& Ai = MatEvalType::GetElem(Atmp, i, j);

		    SP_GRAD_ASSERT(SpGradient::iGetSize(Ai) == 0);
		    SP_GRAD_ASSERT(SpGradient::dGetValue(Ai) == 0.);

		    SpGradient::SetValuePreserve(Ai, dGetValue(i, j));
		    InsertDeriv(Ai, 1., i, j);
	       }
	  }

	  A = std::move(Atmp);
     }

     template <typename ValueType, typename DERIVED>
     template <typename MatEvalType, index_type NumRows, index_type NumCols>
     void SpMatElemExprBase<ValueType, DERIVED>::ElemEvalCompr(SpMatrixBase<ValueType, NumRows, NumCols>& A) const {
	  const index_type iNumRows = iGetNumRows();
	  const index_type iNumCols = iGetNumCols();

	  constexpr bool bTransposed = MatEvalType::eTranspFlag == util::MatTranspEvalFlag::TRANSPOSED;

	  static_assert(NumRows == (!bTransposed ? iNumRowsStatic : iNumColsStatic), "Number of rows does not match");
	  static_assert(NumCols == (!bTransposed ? iNumColsStatic : iNumRowsStatic), "Number of columns does not match");

	  SpMatrixBase<ValueType, NumRows, NumCols> Atmp;

	  if (!bHaveRefTo(A)) {
	       Atmp = std::move(A);
	  }

	  SpGradDofStat oDofStat;

	  for (index_type j = 1; j <= iNumCols; ++j) {
	       for (index_type i = 1; i <= iNumRows; ++i) {
		    GetDofStat(oDofStat, i, j);
	       }
	  }

	  SpGradExpDofMap oDofMap(oDofStat);

	  for (index_type j = 1; j <= iNumCols; ++j) {
	       for (index_type i = 1; i <= iNumRows; ++i) {
		    InsertDof(oDofMap, i, j);
	       }
	  }

	  oDofMap.InsertDone();

	  MatEvalType::ResizeReset(Atmp, iNumRows, iNumCols, oDofStat.iNumNz);

	  SP_GRAD_ASSERT(Atmp.iGetNumRows() == (!bTransposed ? iNumRows : iNumCols));
	  SP_GRAD_ASSERT(Atmp.iGetNumCols() == (!bTransposed ? iNumCols : iNumRows));

	  for (index_type j = 1; j <= iNumCols; ++j) {
	       for (index_type i = 1; i <= iNumRows; ++i) {
		    ValueType& Ai = MatEvalType::GetElem(Atmp, i, j);

		    SP_GRAD_ASSERT(SpGradient::iGetSize(Ai) == 0);
		    SP_GRAD_ASSERT(SpGradient::dGetValue(Ai) == 0.);

		    Ai.SetValuePreserve(dGetValue(i, j));
		    Ai.InitDeriv(oDofMap);
		    AddDeriv(Ai, 1., oDofMap, i, j);
	       }
	  }

	  A = std::move(Atmp);
     }

     namespace util {
	  struct MatMulExprHelper {
	       SpGradDofStat oDofStat;
	       SpGradExpDofMap oDofMap;

	       template <typename TA, typename TB>
	       inline void InnerProduct(typename ResultType<TA, TB>::Type& Aij,
					const TA* pAFirst,
					const TA* pALast,
					index_type iAOffset,
					const TB* pBFirst,
					const TB* pBLast,
					index_type iBOffset) const {
		    util::MapInnerProduct(Aij,
					  pAFirst,
					  pALast,
					  iAOffset,
					  pBFirst,
					  pBLast,
					  iBOffset,
					  oDofMap);
	       }

	       inline void ResetDofStat() {
		    oDofStat = SpGradDofStat{};
	       }

	       inline void ResetDofMap() {
		    oDofMap.Reset(oDofStat);
	       }

	       inline void InsertDone() {
		    oDofMap.InsertDone();
	       }

	       template <typename T>
	       inline void GetDofStat(const T* pFirst, const T* const pLast, const index_type iOffset) noexcept {
		    while (pFirst < pLast) {
			 SpGradient::GetDofStat(*pFirst, oDofStat);
			 pFirst += iOffset;
		    }
	       }

	       template <typename T>
	       inline void InsertDof(const T* pFirst, const T* const pLast, const index_type iOffset) noexcept {
		    while (pFirst < pLast) {
			 SpGradient::InsertDof(*pFirst, oDofMap);
			 pFirst += iOffset;
		    }
	       }
	  };


	  template <typename MatEvalType, bool bIsGradientLhs, bool bIsGradientRhs>
	  struct MatMulExprLoop;

	  template <typename MatEvalType>
	  struct MatMulExprLoop<MatEvalType, true, true> {
	       template <typename MatA, typename MatU, typename MatV>
	       static inline void InnerProduct(MatA& Atmp,
					       const MatU& utmp,
					       const MatV& vtmp,
					       const index_type iNumRows,
					       const index_type iNumCols,
					       const index_type iRowSizeU,
					       const index_type iRowOffsetU,
					       const index_type iColOffsetU,
					       const index_type iRowOffsetV,
					       const index_type iColOffsetV,
					       const index_type iColSizeV) {
		    MatMulExprHelper oMatMulHelper;

		    for (index_type j = 1; j <= iNumCols; ++j) {
			 const SpGradient* const pvtmp_j = vtmp.begin() + (j - 1) * iColOffsetV;

			 for (index_type i = 1; i <= iNumRows; ++i) {
			      const SpGradient* const putmp_i = utmp.begin() + (i - 1) * iRowOffsetU;

			      // Generic method: rebuild dof map for every row and every column
			      oMatMulHelper.ResetDofStat();
			      oMatMulHelper.GetDofStat(putmp_i, putmp_i + iRowSizeU, iColOffsetU);
			      oMatMulHelper.GetDofStat(pvtmp_j, pvtmp_j + iColSizeV, iRowOffsetV);
			      oMatMulHelper.ResetDofMap();
			      oMatMulHelper.InsertDof(putmp_i, putmp_i + iRowSizeU, iColOffsetU);
			      oMatMulHelper.InsertDof(pvtmp_j, pvtmp_j + iColSizeV, iRowOffsetV);
			      oMatMulHelper.InsertDone();

			      oMatMulHelper.InnerProduct(MatEvalType::GetElem(Atmp, i, j),
							 putmp_i,
							 putmp_i + iRowSizeU,
							 iColOffsetU,
							 pvtmp_j,
							 pvtmp_j + iColSizeV,
							 iRowOffsetV);
			 }
		    }
	       }
	  };

	  template <typename MatEvalType>
	  struct MatMulExprLoop<MatEvalType, true, false> {
	       template <typename MatA, typename MatU, typename MatV>
	       static inline void InnerProduct(MatA& Atmp,
					       const MatU& utmp,
					       const MatV& vtmp,
					       const index_type iNumRows,
					       const index_type iNumCols,
					       const index_type iRowSizeU,
					       const index_type iRowOffsetU,
					       const index_type iColOffsetU,
					       const index_type iRowOffsetV,
					       const index_type iColOffsetV,
					       const index_type iColSizeV) {
		    MatMulExprHelper oMatMulHelper;

		    for (index_type i = 1; i <= iNumRows; ++i) {
			 const SpGradient* const putmp_i = utmp.begin() + (i - 1) * iRowOffsetU;

			 // Build dof map once for each row and reuse it for all columns
			 oMatMulHelper.ResetDofStat();
			 oMatMulHelper.GetDofStat(putmp_i, putmp_i + iRowSizeU, iColOffsetU);
			 oMatMulHelper.ResetDofMap();
			 oMatMulHelper.InsertDof(putmp_i, putmp_i + iRowSizeU, iColOffsetU);
			 oMatMulHelper.InsertDone();

			 for (index_type j = 1; j <= iNumCols; ++j) {
			      const doublereal* const pvtmp_j = vtmp.begin() + (j - 1) * iColOffsetV;

			      oMatMulHelper.InnerProduct(MatEvalType::GetElem(Atmp, i, j),
							 putmp_i,
							 putmp_i + iRowSizeU,
							 iColOffsetU,
							 pvtmp_j,
							 pvtmp_j + iColSizeV,
							 iRowOffsetV);
			 }
		    }

	       }
	  };

	  template <typename MatEvalType>
	  struct MatMulExprLoop<MatEvalType, false, true> {
	       template <typename MatA, typename MatU, typename MatV>
	       static inline void InnerProduct(MatA& Atmp,
					       const MatU& utmp,
					       const MatV& vtmp,
					       const index_type iNumRows,
					       const index_type iNumCols,
					       const index_type iRowSizeU,
					       const index_type iRowOffsetU,
					       const index_type iColOffsetU,
					       const index_type iRowOffsetV,
					       const index_type iColOffsetV,
					       const index_type iColSizeV) {
		    MatMulExprHelper oMatMulHelper;

		    for (index_type j = 1; j <= iNumCols; ++j) {
			 const SpGradient* const pvtmp_j = vtmp.begin() + (j - 1) * iColOffsetV;

			 // Build dof map for each column and reuse it for all rows
			 oMatMulHelper.ResetDofStat();
			 oMatMulHelper.GetDofStat(pvtmp_j, pvtmp_j + iColSizeV, iRowOffsetV);
			 oMatMulHelper.ResetDofMap();
			 oMatMulHelper.InsertDof(pvtmp_j, pvtmp_j + iColSizeV, iRowOffsetV);
			 oMatMulHelper.InsertDone();

			 for (index_type i = 1; i <= iNumRows; ++i) {
			      const doublereal* const putmp_i = utmp.begin() + (i - 1) * iRowOffsetU;

			      oMatMulHelper.InnerProduct(MatEvalType::GetElem(Atmp, i, j),
							 putmp_i,
							 putmp_i + iRowSizeU,
							 iColOffsetU,
							 pvtmp_j,
							 pvtmp_j + iColSizeV,
							 iRowOffsetV);
			 }
		    }
	       }
	  };

	  template <typename MatEvalType>
	  struct MatMulExprLoop<MatEvalType, false, false> {
	       template <typename MatA, typename MatU, typename MatV>
	       static inline void InnerProduct(MatA& Atmp,
					       const MatU& utmp,
					       const MatV& vtmp,
					       const index_type iNumRows,
					       const index_type iNumCols,
					       const index_type iRowSizeU,
					       const index_type iRowOffsetU,
					       const index_type iColOffsetU,
					       const index_type iRowOffsetV,
					       const index_type iColOffsetV,
					       const index_type iColSizeV) noexcept {

		    for (index_type j = 1; j <= iNumCols; ++j) {
			 const doublereal* const pvtmp_j = vtmp.begin() + (j - 1) * iColOffsetV;

			 for (index_type i = 1; i <= iNumRows; ++i) {
			      const doublereal* const putmp_i = utmp.begin() + (i - 1) * iRowOffsetU;

			      util::InnerProduct(MatEvalType::GetElem(Atmp, i, j),
						 putmp_i,
						 putmp_i + iRowSizeU,
						 iColOffsetU,
						 pvtmp_j,
						 pvtmp_j + iColSizeV,
						 iRowOffsetV);
			 }
		    }
	       }
	  };
     }

     template <typename LhsValue, typename RhsValue, typename LhsExpr, typename RhsExpr>
     template <typename MatEvalType, index_type NumRowsA, index_type NumColsA>
     void SpMatMulExpr<LhsValue, RhsValue, LhsExpr, RhsExpr>::Eval(SpMatrixBase<ValueType, NumRowsA, NumColsA>& A) const {
	  constexpr bool bTransposedEval = MatEvalType::eTranspFlag == util::MatTranspEvalFlag::TRANSPOSED;

	  static_assert(NumRowsA == (!bTransposedEval ? iNumRowsStatic : iNumColsStatic), "Number of rows does not match");
	  static_assert(NumColsA == (!bTransposedEval ? iNumColsStatic : iNumRowsStatic), "Number of columns does not match");

	  SpMatrixBase<ValueType, NumRowsA, NumColsA> Atmp;

	  if (!bHaveRefTo(A)) {
	       Atmp = std::move(A);
	  }

	  MatEvalType::ResizeReset(Atmp, iGetNumRows(), iGetNumCols(), 0);

	  typedef typename util::remove_all<LhsExpr>::type LhsExprType;
	  typedef typename util::remove_all<RhsExpr>::type RhsExprType;

	  constexpr bool bLhsUsesIterators = (LhsExprType::uMatAccess & util::MatAccessFlag::ITERATORS) != 0;
	  constexpr bool bRhsUsesIterators = (RhsExprType::uMatAccess & util::MatAccessFlag::ITERATORS) != 0;

	  constexpr bool bUseTmpExprLhs = !(LhsExprType::iNumElemOps == 0 && bLhsUsesIterators);
	  constexpr bool bUseTmpExprRhs = !(RhsExprType::iNumElemOps == 0 && bRhsUsesIterators);

	  typename util::TempExprHelper<LhsExpr, bUseTmpExprLhs>::Type utmp = u;
	  typename util::TempExprHelper<RhsExpr, bUseTmpExprRhs>::Type vtmp = v;

	  static_assert(utmp.iNumElemOps == 0);
	  static_assert(vtmp.iNumElemOps == 0);
	  static_assert((utmp.uMatAccess & util::MatAccessFlag::ITERATORS) != 0);
	  static_assert((vtmp.uMatAccess & util::MatAccessFlag::ITERATORS) != 0);

	  const index_type iRowOffsetU = utmp.iGetRowOffset();
	  const index_type iColOffsetU = utmp.iGetColOffset();
	  const index_type iRowSizeU = utmp.iGetColOffset() * utmp.iGetNumCols();
	  const index_type iRowOffsetV = vtmp.iGetRowOffset();
	  const index_type iColOffsetV = vtmp.iGetColOffset();
	  const index_type iColSizeV = vtmp.iGetRowOffset() * vtmp.iGetNumRows();

	  constexpr bool bIsGradientLhs = std::is_same<SpGradient, LhsValue>::value;
	  constexpr bool bIsGradientRhs = std::is_same<SpGradient, RhsValue>::value;

	  typedef util::MatMulExprLoop<MatEvalType, bIsGradientLhs, bIsGradientRhs> MatMulExprLoop;

	  MatMulExprLoop::InnerProduct(Atmp,
				       utmp,
				       vtmp,
				       iGetNumRows(),
				       iGetNumCols(),
				       iRowSizeU,
				       iRowOffsetU,
				       iColOffsetU,
				       iRowOffsetV,
				       iColOffsetV,
				       iColSizeV);

	  A = std::move(Atmp);
     }

     template <typename ValueType, typename ScalarExpr>
     void SpMatElemScalarExpr<ValueType, ScalarExpr>::InsertDof(SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
	  SpGradient::InsertDof(u, oExpDofMap);
     }

     template <typename ValueType, typename ScalarExpr>
     template <typename ValueTypeB>
     void SpMatElemScalarExpr<ValueType, ScalarExpr>::AddDeriv(ValueTypeB& g, doublereal dCoef, const SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
	  SpGradient::AddDeriv(u, g, dCoef, oExpDofMap);
     }

#ifdef SP_GRAD_FORCE_DEFAULT_MATRIX_SIZE
     namespace util {
	  template <bool bStaticSize> struct MatrixDefaultSizeHelper;

	  template <>
	  struct MatrixDefaultSizeHelper<true> {
	       template <typename ValueType, index_type NumRows, index_type NumCols>
	       static void ResizeToDefaultSize(SpMatrixBase<ValueType, NumRows, NumCols>& A) {
		    static_assert(NumRows != SpMatrixSize::DYNAMIC, "No default size for dynamic number of rows");
		    static_assert(NumCols != SpMatrixSize::DYNAMIC, "No default size for dynamic number of columns");

		    A.ResizeReset(NumRows, NumCols, 0);
	       }
	  };

	  template <>
	  struct MatrixDefaultSizeHelper<false> {
	       template <typename ValueType, index_type NumRows, index_type NumCols>
	       static void ResizeToDefaultSize(SpMatrixBase<ValueType, NumRows, NumCols>&) {
		    static_assert(NumRows == SpMatrixSize::DYNAMIC || NumCols == SpMatrixSize::DYNAMIC);
	       }
	  };
     }
#endif

     template <typename ValueType, index_type NumRows, index_type NumCols>
     SpMatrixBase<ValueType, NumRows, NumCols>::SpMatrixBase()
	  :SpMatrixBase(pGetNullData()) {

#ifdef SP_GRAD_FORCE_DEFAULT_MATRIX_SIZE
	  util::MatrixDefaultSizeHelper<bStaticSize>::ResizeToDefaultSize(*this);
#endif
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     SpMatrixBase<ValueType, NumRows, NumCols>::SpMatrixBase(index_type iNumRows, index_type iNumCols, index_type iNumDeriv)
	  :SpMatrixBase(pGetNullData()) {
	  ResizeReset(iNumRows, iNumCols, iNumDeriv);
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     SpMatrixBase<ValueType, NumRows, NumCols>::SpMatrixBase(const SpMatrixBase& oMat)
	  :SpMatrixBase(pGetNullData()) {
	  oMat.Eval(*this);
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     template <typename ValueTypeExpr, typename Expr>
     SpMatrixBase<ValueType, NumRows, NumCols>::SpMatrixBase(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr)
	  :SpMatrixBase(pGetNullData()) {
	  oExpr.Eval(*this);
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     void SpMatrixBase<ValueType, NumRows, NumCols>::ResizeReset(index_type iNumRows, index_type iNumCols, index_type iNumDeriv) {
	  SP_GRAD_ASSERT(iNumRows == NumRows || NumRows == SpMatrixSize::DYNAMIC);
	  SP_GRAD_ASSERT(iNumCols == NumCols || NumCols == SpMatrixSize::DYNAMIC);

	  if (iNumRows == pData->iGetNumRows() &&
	      iNumCols == pData->iGetNumCols() &&
	      iNumDeriv <= pData->iGetMaxDeriv() &&
	      pData->iGetRefCnt() <= 1) {

	       for (auto& g: *pData) {
		    SpGradient::ResizeReset(g, 0., iNumDeriv);
	       }
	  } else {
	       auto pNewData = SpMatrixDataType::pAllocate(iNumRows, iNumCols, iNumDeriv);

	       pData->DecRef();
	       pData = pNewData;
	       pData->IncRef();
	  }
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     SpMatrixBase<ValueType, NumRows, NumCols>::SpMatrixBase(SpMatrixBase&& oMat)
	  :SpMatrixBase(pGetNullData()) {

	  *this = std::move(oMat);

	  SP_GRAD_ASSERT(!bStaticSize);
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     SpMatrixBase<ValueType, NumRows, NumCols>& SpMatrixBase<ValueType, NumRows, NumCols>::operator=(const SpMatrixBase& oMat) {
	  oMat.Eval(*this);
	  return *this;
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     SpMatrixBase<ValueType, NumRows, NumCols>& SpMatrixBase<ValueType, NumRows, NumCols>::operator=(SpMatrixBase&& oMat) {
	  std::swap(pData, oMat.pData);

	  return *this;
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     template <typename ValueTypeExpr, typename Expr>
     SpMatrixBase<ValueType, NumRows, NumCols>& SpMatrixBase<ValueType, NumRows, NumCols>::operator=(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr) {
	  oExpr.Eval(*this);
	  return *this;
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     SpMatrixBase<ValueType, NumRows, NumCols>::~SpMatrixBase() {
	  pData->DecRef();
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     index_type SpMatrixBase<ValueType, NumRows, NumCols>::iGetNumRows() const {
	  return pData->iGetNumRows();
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     index_type SpMatrixBase<ValueType, NumRows, NumCols>::iGetNumCols() const {
	  return pData->iGetNumCols();
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     index_type SpMatrixBase<ValueType, NumRows, NumCols>::iGetSize(index_type i, index_type j) const {
	  return SpGradient::iGetSize(GetElem(i, j));
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     doublereal SpMatrixBase<ValueType, NumRows, NumCols>::dGetValue(index_type i, index_type j) const {
	  return SpGradient::dGetValue(GetElem(i, j));
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     template <typename ValueTypeB>
     void SpMatrixBase<ValueType, NumRows, NumCols>::InsertDeriv(ValueTypeB& g, doublereal dCoef, index_type i, index_type j) const {
	  SpGradient::InsertDeriv(GetElem(i, j), g, dCoef);
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     void SpMatrixBase<ValueType, NumRows, NumCols>::GetDofStat(SpGradDofStat& s, index_type i, index_type j) const {
	  SpGradient::GetDofStat(GetElem(i, j), s);
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     void SpMatrixBase<ValueType, NumRows, NumCols>::InsertDof(SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
	  SpGradient::InsertDof(GetElem(i, j), oExpDofMap);
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     template <typename ValueTypeB>
     void SpMatrixBase<ValueType, NumRows, NumCols>::AddDeriv(ValueTypeB& g, doublereal dCoef, const SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
	  SpGradient::AddDeriv(GetElem(i, j), g, dCoef, oExpDofMap);
     }

     namespace util {
	  template <typename ValueType, typename ValueTypeExpr, typename Expr>
	  struct MatrixBaseRefHelper {
	       template <index_type NumRows, index_type NumCols>
	       static constexpr bool bHaveRefTo(const SpMatrixBase<ValueType, NumRows, NumCols>& A, const SpMatElemExprBase<ValueTypeExpr, Expr>& B) {
		    return false;
	       }
	  };

	  template <typename ValueType>
	  struct MatrixBaseRefHelper<ValueType, ValueType, SpMatrixBase<ValueType> > {
	       template <index_type NumRows, index_type NumCols>
	       static constexpr bool bHaveRefTo(const SpMatrixBase<ValueType, NumRows, NumCols>& A, const SpMatElemExprBase<ValueType, SpMatrixBase<ValueType, NumRows, NumCols> >& B) {
		    return A.bHaveSameRep(*B.pGetRep());
	       }
	  };
     };

     template <typename ValueType, typename ScalarExpr>
     constexpr SpMatElemScalarExpr<ValueType, ScalarExpr>::SpMatElemScalarExpr(const ScalarExpr& u) noexcept
	  :u(u) {

     }

     template <typename ValueType, typename ScalarExpr>
     constexpr doublereal SpMatElemScalarExpr<ValueType, ScalarExpr>::dGetValue(index_type, index_type) const noexcept {
	  return SpGradient::dGetValue(u);
     }

     template <typename ValueType, typename ScalarExpr>
     constexpr index_type SpMatElemScalarExpr<ValueType, ScalarExpr>::iGetSize(index_type, index_type) const noexcept {
	  return SpGradient::iGetSize(u);
     }

     template <typename ValueType, typename ScalarExpr>
     constexpr index_type SpMatElemScalarExpr<ValueType, ScalarExpr>::iGetNumRows() noexcept {
	  return 1;
     }

     template <typename ValueType, typename ScalarExpr>
     constexpr index_type SpMatElemScalarExpr<ValueType, ScalarExpr>::iGetNumCols() noexcept {
	  return 1;
     }

     template <typename ValueType, typename ScalarExpr>
     constexpr index_type SpMatElemScalarExpr<ValueType, ScalarExpr>::iGetMaxSize() const noexcept {
	  return SpGradient::iGetSize(u);
     }

     template <typename ValueType, typename ScalarExpr>
     template <typename ValueTypeB>
     void SpMatElemScalarExpr<ValueType, ScalarExpr>::InsertDeriv(ValueTypeB& g, doublereal dCoef, index_type, index_type) const {
	  SpGradient::InsertDeriv(u, g, dCoef);
     }

     template <typename ValueType, typename ScalarExpr>
     void SpMatElemScalarExpr<ValueType, ScalarExpr>::GetDofStat(SpGradDofStat& s, index_type, index_type) const noexcept {
	  SpGradient::GetDofStat(u, s);
     }

     template <typename ValueType, typename ScalarExpr>
     template <typename ExprTypeB, typename ExprB>
     constexpr bool SpMatElemScalarExpr<ValueType, ScalarExpr>::bHaveRefTo(const SpMatElemExprBase<ExprTypeB, ExprB>& A) const noexcept {
	  return SpGradient::bHaveRefTo(u, *A.pGetRep());
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     template <typename ExprType, typename Expr>
     constexpr bool SpMatrixBase<ValueType, NumRows, NumCols>::bHaveRefTo(const SpMatElemExprBase<ExprType, Expr>& A) const noexcept {
	  return util::MatrixBaseRefHelper<ValueType, ExprType, Expr>::bHaveRefTo(*this, A);
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     ValueType* SpMatrixBase<ValueType, NumRows, NumCols>::begin() {
	  return pData->begin();
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     ValueType* SpMatrixBase<ValueType, NumRows, NumCols>::end() {
	  return pData->end();
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     const ValueType* SpMatrixBase<ValueType, NumRows, NumCols>::begin() const {
	  return pData->begin();
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     const ValueType* SpMatrixBase<ValueType, NumRows, NumCols>::end() const {
	  return pData->end();
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     bool SpMatrixBase<ValueType, NumRows, NumCols>::bHaveSameRep(const SpMatrixBase<ValueType, NumRows, NumCols>& A) const noexcept {
	  return pData == A.pData && pData != pGetNullData();
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     bool SpMatrixBase<ValueType, NumRows, NumCols>::bIsOwnerOf(const SpMatrixData<ValueType>* pDataB) const noexcept {
	  return pData == pDataB;
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     const ValueType& SpMatrixBase<ValueType, NumRows, NumCols>::GetElem(index_type i) const {
	  SP_GRAD_ASSERT(pData != pGetNullData());

	  return *pData->pGetData(i);
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     ValueType& SpMatrixBase<ValueType, NumRows, NumCols>::GetElem(index_type i) {
	  SP_GRAD_ASSERT(pData != pGetNullData());

	  return *pData->pGetData(i);
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     const ValueType& SpMatrixBase<ValueType, NumRows, NumCols>::GetElem(index_type i, index_type j) const {
	  SP_GRAD_ASSERT(pData != pGetNullData());

	  return *pData->pGetData(i, j);
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     ValueType& SpMatrixBase<ValueType, NumRows, NumCols>::GetElem(index_type i, index_type j) {
	  SP_GRAD_ASSERT(pData != pGetNullData());

	  return *pData->pGetData(i, j);
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     typename SpMatrixBase<ValueType, NumRows, NumCols>::SpMatrixDataType*
     SpMatrixBase<ValueType, NumRows, NumCols>::pGetNullData() {
	  return static_cast<SpMatrixDataType*>(SpMatrixBaseData::pGetNullData());
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     template <typename MatEvalType, index_type NumRowsA, index_type NumColsA>
     void SpMatrixBase<ValueType, NumRows, NumCols>::Eval(SpMatrixBase<ValueType, NumRowsA, NumColsA>& A) const {
	  this->template ElemEval<MatEvalType>(A);
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     index_type SpMatrixBase<ValueType, NumRows, NumCols>::iGetMaxSize() const {
	  return this->iGetMaxSizeElem();
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     SpMatrixBase<ValueType, NumRows, NumCols>::SpMatrixBase(SpMatrixDataType* pData)
	  :pData(pData) {
	  pData->IncRef();
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     SpMatrix<ValueType, NumRows, NumCols>::SpMatrix(index_type iNumRows, index_type iNumCols, index_type iNumDeriv)
	  :SpMatrixBase<ValueType, NumRows, NumCols>(iNumRows, iNumCols, iNumDeriv) {
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     ValueType& SpMatrix<ValueType, NumRows, NumCols>::operator()(index_type iRow, index_type iCol) {
	  return this->GetElem(iRow, iCol);
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     inline const ValueType& SpMatrix<ValueType, NumRows, NumCols>::operator() (index_type iRow, index_type iCol) const {
	  return this->GetElem(iRow, iCol);
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     template <typename ValueTypeExpr, typename Expr>
     SpMatrix<ValueType, NumRows, NumCols>::SpMatrix(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr)
	  :SpMatrixBase<ValueType, NumRows, NumCols>(oExpr) {
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     template <typename ValueTypeExpr, typename Expr>
     SpMatrix<ValueType, NumRows, NumCols>& SpMatrix<ValueType, NumRows, NumCols>::operator=(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr) {
	  SpMatrixBase<ValueType, NumRows, NumCols>::operator=(oExpr);
	  return *this;
     }

     template <typename ValueType, index_type NumRows>
     SpColVector<ValueType, NumRows>::SpColVector(index_type iNumRows, index_type iNumDeriv)
	  :SpMatrixBase<ValueType, NumRows, 1>(iNumRows, 1, iNumDeriv) {
     }

     template <typename ValueType, index_type NumRows>
     template <typename ValueTypeExpr, typename Expr>
     SpColVector<ValueType, NumRows>::SpColVector(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr)
	  :SpMatrixBase<ValueType, NumRows, 1>(oExpr) {
     }

     template <typename ValueType, index_type NumRows>
     ValueType& SpColVector<ValueType, NumRows>::operator()(index_type iRow) {
	  return this->GetElem(iRow);
     }

     template <typename ValueType, index_type NumRows>
     const ValueType& SpColVector<ValueType, NumRows>::operator() (index_type iRow) const {
	  return this->GetElem(iRow);
     }

     template <typename ValueType, index_type NumRows>
     template <typename ValueTypeExpr, typename Expr>
     SpColVector<ValueType, NumRows>& SpColVector<ValueType, NumRows>::operator=(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr) {
	  SpMatrixBase<ValueType, NumRows, 1>::operator=(oExpr);
	  return *this;
     }

     template <typename ValueType, index_type NumCols>
     SpRowVector<ValueType, NumCols>::SpRowVector(index_type iNumCols, index_type iNumDeriv)
	  :SpMatrixBase<ValueType, 1, NumCols>(1, iNumCols, iNumDeriv) {
     }

     template <typename ValueType, index_type NumCols>
     template <typename ValueTypeExpr, typename Expr>
     SpRowVector<ValueType, NumCols>::SpRowVector(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr)
	  :SpMatrixBase<ValueType, 1, NumCols>(oExpr) {
     }

     template <typename ValueType, index_type NumCols>
     ValueType& SpRowVector<ValueType, NumCols>::operator()(index_type iCol) {
	  return this->GetElem(iCol);
     }

     template <typename ValueType, index_type NumCols>
     const ValueType& SpRowVector<ValueType, NumCols>::operator() (index_type iCol) const {
	  return this->GetElem(iCol);
     }

     template <typename ValueType, index_type NumCols>
     template <typename ValueTypeExpr, typename Expr>
     SpRowVector<ValueType, NumCols>& SpRowVector<ValueType, NumCols>::operator=(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr) {
	  SpMatrixBase<ValueType, 1, NumCols>::operator=(oExpr);
	  return *this;
     }

     template <typename LhsValue, typename RhsValue, typename LhsExpr, typename RhsExpr>
     inline constexpr
     SpMatElemBinExpr<typename util::ResultType<LhsValue, RhsValue>::Type,
		      SpGradBinPlus,
		      const SpMatElemExprBase<LhsValue, LhsExpr>&,
		      const SpMatElemExprBase<RhsValue, RhsExpr>&>
     operator+(const SpMatElemExprBase<LhsValue, LhsExpr>& A,
	       const SpMatElemExprBase<RhsValue, RhsExpr>& B) noexcept {
	  return decltype(operator+(A, B))(A, B);
     }

     template <typename LhsValue, typename RhsValue, typename LhsExpr, typename RhsExpr>
     inline constexpr
     SpMatElemBinExpr<typename util::ResultType<LhsValue, RhsValue>::Type,
		      SpGradBinMinus,
		      const SpMatElemExprBase<LhsValue, LhsExpr>&,
		      const SpMatElemExprBase<RhsValue, RhsExpr>&>
     operator-(const SpMatElemExprBase<LhsValue, LhsExpr>& A,
	       const SpMatElemExprBase<RhsValue, RhsExpr>& B) noexcept {
	  return decltype(operator-(A, B))(A, B);
     }

     template <typename Value, typename Expr>
     inline constexpr
     SpMatElemUnaryExpr<Value,
			SpGradUnaryMinus,
			const SpMatElemExprBase<Value, Expr>& >
     operator-(const SpMatElemExprBase<Value, Expr>& A) noexcept {
	  return decltype(operator-(A))(A);
     }

     template <typename Value, typename Expr>
     inline constexpr
     SpMatElemTranspExpr<Value, const SpMatElemExprBase<Value, Expr>&>
     Transpose(const SpMatElemExprBase<Value, Expr>& A) noexcept {
	  return decltype(Transpose(A))(A);
     }

     template <typename Value, typename Expr>
     inline constexpr
     const SpMatElemExprBase<Value, Expr>&
     Transpose(const SpMatElemTranspExpr<Value, const SpMatElemExprBase<Value, Expr>&>& A) noexcept {
	  return decltype(Transpose(A))(A.Transpose());
     }

     template <typename Value, typename Expr>
     inline constexpr
     SpMatElemComprExpr<Value, const SpMatElemExprBase<Value, Expr>&>
     EvalCompressed(const SpMatElemExprBase<Value, Expr>& A) noexcept {
	  return decltype(EvalCompressed(A))(A);
     }

     template <typename LhsValue, typename LhsExpr, typename RhsExpr>
     inline constexpr
     SpMatElemBinExpr<typename util::ResultType<LhsValue, SpGradient>::Type,
		      SpGradBinMult,
		      const SpMatElemExprBase<LhsValue, LhsExpr>&,
		      const SpMatElemScalarExpr<SpGradient, SpGradient> >
     operator*(const SpMatElemExprBase<LhsValue, LhsExpr>& A,
	       const SpGradBase<RhsExpr>& b) noexcept {
	  return decltype(operator*(A, b))(A, SpGradient{b}); // Avoid multiple evaluations of b!
     }

     template <typename LhsValue, typename LhsExpr>
     inline constexpr
     SpMatElemBinExpr<typename util::ResultType<LhsValue, doublereal>::Type,
		      SpGradBinMult,
		      const SpMatElemExprBase<LhsValue, LhsExpr>&,
		      const SpMatElemScalarExpr<doublereal, doublereal> >
     operator*(const SpMatElemExprBase<LhsValue, LhsExpr>& A,
	       const doublereal b) noexcept {
	  return decltype(operator*(A, b))(A, b);
     }

     template <typename LhsValue, typename LhsExpr, typename RhsExpr>
     inline constexpr
     SpMatElemBinExpr<typename util::ResultType<LhsValue, SpGradient>::Type,
		      SpGradBinDiv,
		      const SpMatElemExprBase<LhsValue, LhsExpr>&,
		      const SpMatElemScalarExpr<SpGradient, SpGradient> >
     operator/(const SpMatElemExprBase<LhsValue, LhsExpr>& A,
	       const SpGradBase<RhsExpr>& b) noexcept {
	  return decltype(operator/(A, b))(A, SpGradient{b}); // Avoid multiple evaluations of b!
     }

     template <typename LhsValue, typename LhsExpr>
     inline constexpr
     SpMatElemBinExpr<typename util::ResultType<LhsValue, doublereal>::Type,
		      SpGradBinDiv,
		      const SpMatElemExprBase<LhsValue, LhsExpr>&,
		      const SpMatElemScalarExpr<doublereal, doublereal> >
     operator/(const SpMatElemExprBase<LhsValue, LhsExpr>& A,
	       const doublereal b) noexcept {
	  return decltype(operator/(A, b))(A, b);
     }

     template <typename LhsValue, typename RhsValue, typename LhsExpr, typename RhsExpr>
     inline constexpr
     SpMatMulExpr<LhsValue,
		  RhsValue,
		  const SpMatElemExprBase<LhsValue, LhsExpr>&,
		  const SpMatElemExprBase<RhsValue, RhsExpr>&>
     operator*(const SpMatElemExprBase<LhsValue, LhsExpr>& A,
	       const SpMatElemExprBase<RhsValue, RhsExpr>& B) noexcept {
	  return decltype(operator*(A, B))(A, B);
     }

     template <typename LhsValue, typename RhsValue, typename LhsExpr, typename RhsExpr>
     inline constexpr
     typename util::ResultType<LhsValue, RhsValue>::Type
     Dot(const SpMatElemExprBase<LhsValue, LhsExpr>& u, const SpMatElemExprBase<RhsValue, RhsExpr>& v) {
	  return *SpMatrixBase<typename util::ResultType<LhsValue, RhsValue>::Type, 1, 1>(Transpose(u) * v).begin();
     }
}

#endif
