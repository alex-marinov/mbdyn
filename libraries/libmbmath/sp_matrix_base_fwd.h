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

#ifndef __SP_MATRIX_BASE_FWD_H__INCLUDED__
#define __SP_MATRIX_BASE_FWD_H__INCLUDED__

#include "sp_gradient_base.h"

namespace sp_grad {
     struct SpMatrixSize {
	  enum SpMatrixSizeType: index_type {
	       DYNAMIC = -1
	  };
     };

     enum class SpMatOpType: index_type {
	  SCALAR,
	  MATRIX
     };

     template <typename ValueType, index_type NumRows = SpMatrixSize::DYNAMIC, index_type NumCols = SpMatrixSize::DYNAMIC>
     class SpMatrixBase;

     template <typename ValueType>
     class SpMatrixData;

     namespace util {
	  template <typename T>
	  struct remove_all {
	       typedef typename std::remove_cv<typename std::remove_reference<T>::type>::type type;
	  };

	  template <typename ValueType>
	  struct SpMatrixDataTraits {
	       static inline void Construct(SpMatrixData<ValueType>& oData, index_type iNumDeriv, void* pExtraMem);
	       static constexpr inline size_t uOffsetDeriv(size_t uSizeGrad) noexcept;
	       static constexpr inline size_t uSizeDeriv(index_type iNumDeriv, index_type iNumItems) noexcept;
	       static constexpr inline size_t uAlign() noexcept;
	  };

	  template <>
	  struct SpMatrixDataTraits<SpGradient> {
	       static inline void Construct(SpMatrixData<SpGradient>& oData, index_type iNumDeriv, void* pExtraMem);
	       static constexpr inline size_t uOffsetDeriv(size_t uSizeGrad) noexcept;
	       static constexpr inline size_t uSizeDeriv(index_type iNumDeriv, index_type iNumItems) noexcept;
	       static constexpr inline size_t uAlign() noexcept;
	  };
     }

     class SpMatrixBaseData {
     protected:
	  inline SpMatrixBaseData(index_type iNumRows,
				  index_type iNumCols,
				  index_type iRefCnt,
				  index_type iNumDeriv) noexcept;

     public:
	  inline index_type iGetNumRows() const noexcept;
	  inline index_type iGetNumCols() const noexcept;
	  inline index_type iGetNumElem() const noexcept;
	  inline index_type iGetRefCnt() const noexcept;
	  inline index_type iGetMaxDeriv() const noexcept;
	  inline static SpMatrixBaseData* pGetNullData() noexcept;

     protected:
	  index_type iNumRows;
	  index_type iNumCols;
	  index_type iRefCnt;
	  index_type iNumDeriv;

	  static SP_GRAD_THREAD_LOCAL SpMatrixBaseData oNullData;
     };

     template <typename ValueType>
     class SpMatrixData: public SpMatrixBaseData {
     protected:
	  inline SpMatrixData(index_type iNumRows,
			      index_type iNumCols,
			      index_type iRefCnt,
			      index_type iNumDeriv,
			      void* pExtraMem);
	  inline ~SpMatrixData();

     public:
	  inline ValueType* pGetData() noexcept;
	  inline ValueType* pGetData(index_type iRow) noexcept;
	  inline const ValueType* pGetData() const noexcept;
	  inline const ValueType* pGetData(index_type iRow) const noexcept;
	  inline ValueType* pGetData(index_type iRow, index_type iCol) noexcept;
	  inline const ValueType* pGetData(index_type iRow, index_type iCol) const noexcept;
	  inline void Attach(ValueType* pData) noexcept;
	  inline void Detach(ValueType* pData);
	  inline void IncRef() noexcept;
	  inline void DecRef();
	  inline ValueType* begin() noexcept;
	  inline ValueType* end() noexcept;
	  inline const ValueType* begin() const noexcept;
	  inline const ValueType* end() const noexcept;
	  constexpr inline bool bIsOwnerOf(const ValueType* pData) const noexcept;
     };

     template <typename ValueType, index_type NumRows, index_type NumCols>
     class SpMatrixDataCTD: public SpMatrixData<ValueType> {
     private:
	  static constexpr index_type iNumRowsStatic = NumRows;
	  static constexpr index_type iNumColsStatic = NumCols;

	  inline SpMatrixDataCTD(index_type iNumRows,
				 index_type iNumCols,
				 index_type iRefCnt,
				 index_type iNumDeriv,
				 void* pExtraMem);
     public:
	  static inline SpMatrixDataCTD* pAllocate(index_type iNumRows, index_type iNumCols, index_type iNumDeriv);
     };

     namespace util {
	  enum class MatTranspEvalFlag: index_type  {
	       DIRECT,
	       TRANSPOSED
	  };

	  template <MatTranspEvalFlag eTransp, SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalUncompressed>
	  struct MatEvalHelper;

	  template <SpGradCommon::ExprEvalFlags eCompr>
	  struct MatEvalHelperCompr;

	  template <>
	  struct MatEvalHelperCompr<SpGradCommon::ExprEvalUncompressed> {
	       static constexpr SpGradCommon::ExprEvalFlags eComprFlag = SpGradCommon::ExprEvalUncompressed;

	       template <typename MatEvalType, typename Expr, typename ValueType, index_type NumRows, index_type NumCols>
	       static inline void ElemEval(const Expr& oExpr, SpMatrixBase<ValueType, NumRows, NumCols>& A);
	  };

	  template <>
	  struct MatEvalHelperCompr<SpGradCommon::ExprEvalCompressed> {
	       static constexpr SpGradCommon::ExprEvalFlags eComprFlag = SpGradCommon::ExprEvalCompressed;

	       template <typename MatEvalType, typename Expr, typename ValueType, index_type NumRows, index_type NumCols>
	       static inline void ElemEval(const Expr& oExpr, SpMatrixBase<ValueType, NumRows, NumCols>& A);
	  };

	  template <MatTranspEvalFlag eTransp>
	  struct MatEvalHelperTransp;

	  template <>
	  struct MatEvalHelperTransp<MatTranspEvalFlag::DIRECT> {
	       static constexpr MatTranspEvalFlag eTranspFlag = MatTranspEvalFlag::DIRECT;

	       template <typename ValueType, index_type NumRows, index_type NumCols>
	       static inline void ResizeReset(SpMatrixBase<ValueType, NumRows, NumCols>& A,
					      index_type iNumRows,
					      index_type iNumCols,
					      index_type iMaxSize);

	       template <typename ValueType, index_type NumRows, index_type NumCols>
	       static inline ValueType& GetElem(SpMatrixBase<ValueType, NumRows, NumCols>& A,
						index_type i,
						index_type j);
	  };

	  template <>
	  struct MatEvalHelperTransp<MatTranspEvalFlag::TRANSPOSED> {
	       static constexpr MatTranspEvalFlag eTranspFlag = MatTranspEvalFlag::TRANSPOSED;

	       template <typename ValueType, index_type NumRows, index_type NumCols>
	       static inline void ResizeReset(SpMatrixBase<ValueType, NumRows, NumCols>& A,
					      index_type iNumRows,
					      index_type iNumCols,
					      index_type iMaxSize);

	       template <typename ValueType, index_type NumRows, index_type NumCols>
	       static inline ValueType& GetElem(SpMatrixBase<ValueType, NumRows, NumCols>& A,
						index_type i,
						index_type j);
	  };

	  template <typename MatEvalType> struct TranspMatEvalHelper;

	  template <>
	  struct TranspMatEvalHelper<MatEvalHelperTransp<MatTranspEvalFlag::DIRECT> > {
	       typedef MatEvalHelperTransp<MatTranspEvalFlag::TRANSPOSED> Type;
	  };

	  template <>
	  struct TranspMatEvalHelper<MatEvalHelperTransp<MatTranspEvalFlag::TRANSPOSED> > {
	       typedef MatEvalHelperTransp<MatTranspEvalFlag::DIRECT> Type;
	  };

	  enum MatAccessFlag: unsigned {
	       ELEMENT_WISE = 0x1,
	       ITERATORS = 0x2,
	       MATRIX_WISE = 0x4
	  };

	  template <typename Expr, bool bUseTempExpr>
	  struct TempExprHelper;

	  template <typename Expr>
	  struct TempExprHelper<Expr, true> {
	       typedef typename remove_all<Expr>::type ExprType;
	       typedef typename ExprType::ValueType ValueType;

	       // Don't use const here in order to enable move constructors
	       typedef SpMatrixBase<ValueType, ExprType::iNumRowsStatic, ExprType::iNumColsStatic> Type;
	  };

	  template <typename Expr>
	  struct TempExprHelper<Expr, false> {
	       typedef const Expr Type;
	  };
     };

     template <typename VALUE, typename DERIVED>
     class SpMatElemExprBase
     {
     protected:
	  constexpr SpMatElemExprBase() noexcept {}
	  ~SpMatElemExprBase() noexcept {}

     public:
	  typedef VALUE ValueType;
	  typedef DERIVED DerivedType;
	  static constexpr index_type iNumElemOps = DERIVED::iNumElemOps;
	  static constexpr unsigned uMatAccess = DERIVED::uMatAccess;
	  static constexpr index_type iNumRowsStatic = DERIVED::iNumRowsStatic;
	  static constexpr index_type iNumColsStatic = DERIVED::iNumColsStatic;
	  static constexpr SpMatOpType eMatOpType = DERIVED::eMatOpType;
	  static constexpr SpGradCommon::ExprEvalFlags eComprFlag = DERIVED::eComprFlag;

	  template <typename MatEvalType = util::MatEvalHelper<util::MatTranspEvalFlag::DIRECT>, index_type NumRows, index_type NumCols>
	  inline void ElemEval(SpMatrixBase<ValueType, NumRows, NumCols>& A) const;

	  template <typename MatEvalType, index_type NumRows, index_type NumCols>
	  inline void ElemEvalUncompr(SpMatrixBase<ValueType, NumRows, NumCols>& A) const;

	  template <typename MatEvalType, index_type NumRows, index_type NumCols>
	  inline void ElemEvalCompr(SpMatrixBase<ValueType, NumRows, NumCols>& A) const;

	  template <typename MatEvalType = util::MatEvalHelperTransp<util::MatTranspEvalFlag::DIRECT>, index_type NumRows, index_type NumCols>
	  inline void Eval(SpMatrixBase<ValueType, NumRows, NumCols>& A) const {
	       pGetRep()->template Eval<MatEvalType>(A);
	  }

	  index_type iGetNumRows() const {
	       return pGetRep()->iGetNumRows();
	  }

	  index_type iGetNumCols() const {
	       return pGetRep()->iGetNumCols();
	  }

	  constexpr doublereal dGetValue(index_type i, index_type j) const {
	       SP_GRAD_ASSERT(&SpMatElemExprBase::dGetValue != &DERIVED::dGetValue);
	       return pGetRep()->dGetValue(i, j);
	  }

	  template <typename ExprType, typename Expr>
	  constexpr bool bHaveRefTo(const SpMatElemExprBase<ExprType, Expr>& A) const {
	       return pGetRep()->bHaveRefTo(A);
	  }

	  constexpr index_type iGetSize(index_type i, index_type j) const {
	       SP_GRAD_ASSERT(&SpMatElemExprBase::iGetSize != &DERIVED::iGetSize);
	       return pGetRep()->iGetSize(i, j);
	  }

	  index_type iGetMaxSizeElem() const {
	       index_type iMaxSize = 0;

	       for (index_type j = 1; j <= iGetNumCols(); ++j) {
		    for (index_type i = 1; i <= iGetNumRows(); ++i) {
			 iMaxSize = std::max(iMaxSize, iGetSize(i, j));
		    }
	       }

	       return iMaxSize;
	  }

	  index_type iGetMaxSize() const {
	       return pGetRep()->iGetMaxSize();
	  }

	  template <typename ValueTypeB>
	  void InsertDeriv(ValueTypeB& g, doublereal dCoef, index_type i, index_type j) const {
	       pGetRep()->InsertDeriv(g, dCoef, i, j);
	  }

	  void GetDofStat(SpGradDofStat& s, index_type i, index_type j) const {
	       pGetRep()->GetDofStat(s, i, j);
	  }

	  void InsertDof(SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
	       pGetRep()->InsertDof(oExpDofMap, i, j);
	  }

	  template <typename ValueTypeB>
	  void AddDeriv(ValueTypeB& g, doublereal dCoef, const SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
	       pGetRep()->AddDeriv(g, dCoef, oExpDofMap, i, j);
	  }

	  const ValueType* begin() const {
	       return pGetRep()->begin();
	  }

	  const ValueType* end() const {
	       return pGetRep()->end();
	  }

	  inline constexpr index_type iGetRowOffset() const noexcept { return pGetRep()->iGetRowOffset(); }
	  inline constexpr index_type iGetColOffset() const noexcept { return pGetRep()->iGetColOffset(); }

	  constexpr const DERIVED* pGetRep() const noexcept {
	       return static_cast<const DERIVED*>(this);
	  }
     };

     template <typename ValueType, typename ScalarExpr>
     class SpMatElemScalarExpr;

     namespace util {
	  template <typename UValue, typename UExpr, typename VValue, typename VExpr>
	  struct MatrixSizeHelper {
	       static index_type iGetNumRows(const SpMatElemExprBase<UValue, UExpr>& u,
					     const SpMatElemExprBase<VValue, VExpr>& v) noexcept {
		    SP_GRAD_ASSERT(u.iGetNumRows() == v.iGetNumRows());
		    return u.iGetNumRows();
	       }

	       static index_type iGetNumCols(const SpMatElemExprBase<UValue, UExpr>& u,
					     const SpMatElemExprBase<VValue, VExpr>& v) noexcept {
		    SP_GRAD_ASSERT(u.iGetNumCols() == v.iGetNumCols());
		    return u.iGetNumCols();
	       }
	  };

	  template <typename UValue, typename UExpr, typename VValue, typename VExpr>
	  struct MatrixSizeHelper<UValue, SpMatElemScalarExpr<UValue, UExpr>, VValue, VExpr> {
	       static index_type iGetNumRows(const SpMatElemExprBase<UValue, SpMatElemScalarExpr<UValue, UExpr> >& u,
					     const SpMatElemExprBase<VValue, VExpr>& v) noexcept {
		    SP_GRAD_ASSERT(u.iGetNumRows() == 1);
		    return v.iGetNumRows();
	       }

	       static index_type iGetNumCols(const SpMatElemExprBase<UValue, SpMatElemScalarExpr<UValue, UExpr> >& u,
					     const SpMatElemExprBase<VValue, VExpr>& v) noexcept {
		    SP_GRAD_ASSERT(u.iGetNumCols() == 1);
		    return v.iGetNumCols();
	       }
	  };

	  template <typename UValue, typename UExpr, typename VValue, typename VExpr>
	  struct MatrixSizeHelper<UValue, UExpr, VValue, SpMatElemScalarExpr<VValue, VExpr> > {
	       static index_type iGetNumRows(const SpMatElemExprBase<UValue, UExpr>& u,
					     const SpMatElemExprBase<VValue, SpMatElemScalarExpr<VValue, VExpr> >& v) {
		    SP_GRAD_ASSERT(v.iGetNumRows() == 1);
		    return u.iGetNumRows();
	       }

	       static index_type iGetNumCols(const SpMatElemExprBase<UValue, UExpr>& u,
					     const SpMatElemExprBase<VValue, SpMatElemScalarExpr<VValue, VExpr> >& v) {
		    SP_GRAD_ASSERT(v.iGetNumCols() == 1);
		    return u.iGetNumCols();
	       }
	  };
     }

     template <typename ValueType, typename BinaryFunc, typename LhsExpr, typename RhsExpr>
     class SpMatElemBinExpr
	  :public SpMatElemExprBase<ValueType, SpMatElemBinExpr<ValueType, BinaryFunc, LhsExpr, RhsExpr> > {
	  typedef typename util::remove_all<LhsExpr>::type LhsExprType;
	  typedef typename util::remove_all<RhsExpr>::type RhsExprType;
	  static constexpr bool bUseTempExprLhs = !(LhsExprType::uMatAccess & util::MatAccessFlag::ELEMENT_WISE);
	  static constexpr bool bUseTempExprRhs = !(RhsExprType::uMatAccess & util::MatAccessFlag::ELEMENT_WISE);
	  typedef typename util::TempExprHelper<LhsExpr, bUseTempExprLhs>::Type LhsTempExpr;
	  typedef typename util::TempExprHelper<RhsExpr, bUseTempExprRhs>::Type RhsTempExpr;
	  typedef typename util::remove_all<LhsTempExpr>::type LhsTempExprType;
	  typedef typename util::remove_all<RhsTempExpr>::type RhsTempExprType;
	  typedef typename LhsTempExprType::ValueType LhsValueType;
	  typedef typename RhsTempExprType::ValueType RhsValueType;
	  typedef typename LhsTempExprType::DerivedType LhsDerivedType;
	  typedef typename RhsTempExprType::DerivedType RhsDerivedType;
	  typedef util::MatrixSizeHelper<LhsValueType, LhsDerivedType, RhsValueType, RhsDerivedType> MatSizeHelp;
	  static constexpr bool bMinOneScalarOp = LhsExprType::eMatOpType == SpMatOpType::SCALAR ||
	       RhsExprType::eMatOpType == SpMatOpType::SCALAR;
     public:
	  static constexpr index_type iNumElemOps = LhsTempExprType::iNumElemOps + RhsTempExprType::iNumElemOps + 1;
	  static constexpr unsigned uMatAccess = util::MatAccessFlag::ELEMENT_WISE;
	  static constexpr index_type iNumRowsStatic = LhsExprType::iNumRowsStatic;
	  static constexpr index_type iNumColsStatic = LhsExprType::iNumColsStatic;
	  static constexpr SpMatOpType eMatOpType = SpMatOpType::MATRIX;
	  static constexpr SpGradCommon::ExprEvalFlags eComprFlag = util::ExprEvalFlagsHelper<LhsExprType::eComprFlag,
											      RhsExprType::eComprFlag>::eExprEvalFlags;
	  static_assert(LhsExprType::eMatOpType == SpMatOpType::MATRIX || RhsExprType::eMatOpType == SpMatOpType::MATRIX,
			"At least one operand must be a matrix! Use SpGradient instead if both operands are scalar!");
	  static_assert(LhsExprType::iNumRowsStatic == RhsExprType::iNumRowsStatic || bMinOneScalarOp,
			"Number of rows of two matrix operands do not match!");
	  static_assert(LhsExprType::iNumColsStatic == RhsExprType::iNumColsStatic || bMinOneScalarOp,
			"Number of columns of two matrix operands do not match!");

	  constexpr SpMatElemBinExpr(const LhsExprType& u, const RhsExprType& v) noexcept
	       :u(u), v(v) {
	  }

	  SpMatElemBinExpr(SpMatElemBinExpr&& oExpr)
	       :u(std::move(oExpr.u)), v(std::move(oExpr.v)) {
	  }

	  template <typename MatEvalType = util::MatEvalHelperTransp<util::MatTranspEvalFlag::DIRECT>, index_type iNumRowsA, index_type iNumColsA >
	  void Eval(SpMatrixBase<ValueType, iNumRowsA, iNumColsA>& A) const {
	       this->template ElemEval<MatEvalType>(A);
	  }

	  constexpr doublereal dGetValue(index_type i, index_type j) const {
	       return BinaryFunc::f(u.dGetValue(i, j), v.dGetValue(i, j));
	  }

	  constexpr index_type iGetSize(index_type i, index_type j) const {
	       return u.iGetSize(i, j) + v.iGetSize(i, j);
	  }

	  constexpr index_type iGetNumRows() const {
	       return MatSizeHelp::iGetNumRows(u, v);
	  }

	  constexpr index_type iGetNumCols() const {
	       return MatSizeHelp::iGetNumCols(u, v);
	  }

	  constexpr index_type iGetMaxSize() const {
	       return this->iGetMaxSizeElem();
	  }

	  template <typename ValueTypeB>
	  void InsertDeriv(ValueTypeB& g, doublereal dCoef, index_type i, index_type j) const {
	       const doublereal ui = u.dGetValue(i, j);
	       const doublereal vi = v.dGetValue(i, j);

	       u.InsertDeriv(g, BinaryFunc::df_du(ui, vi) * dCoef, i, j);
	       v.InsertDeriv(g, BinaryFunc::df_dv(ui, vi) * dCoef, i, j);
	  }

	  void GetDofStat(SpGradDofStat& s, index_type i, index_type j) const {
	       u.GetDofStat(s, i, j);
	       v.GetDofStat(s, i, j);
	  }

	  void InsertDof(SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
	       u.InsertDof(oExpDofMap, i, j);
	       v.InsertDof(oExpDofMap, i, j);
	  }

	  template <typename ValueTypeB>
	  void AddDeriv(ValueTypeB& g, doublereal dCoef, const SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
	       const doublereal ui = u.dGetValue(i, j);
	       const doublereal vi = v.dGetValue(i, j);

	       u.AddDeriv(g, BinaryFunc::df_du(ui, vi) * dCoef, oExpDofMap, i, j);
	       v.AddDeriv(g, BinaryFunc::df_dv(ui, vi) * dCoef, oExpDofMap, i, j);
	  }


	  template <typename ExprType, typename Expr>
	  constexpr bool bHaveRefTo(const SpMatElemExprBase<ExprType, Expr>& A) const {
	       return u.bHaveRefTo(A) || v.bHaveRefTo(A);
	  }

	  const ValueType* begin() const = delete;
	  const ValueType* end() const = delete;
	  inline constexpr index_type iGetRowOffset() const = delete;
	  inline constexpr index_type iGetColOffset() const = delete;
     private:
	  LhsTempExpr u;
	  RhsTempExpr v;
     };

     template <typename ValueType, typename UnaryFunc, typename Expr>
     class SpMatElemUnaryExpr
	  :public SpMatElemExprBase<ValueType, SpMatElemUnaryExpr<ValueType, UnaryFunc, Expr> > {
	  typedef typename util::remove_all<Expr>::type ExprType;
	  static constexpr bool bUseTempExpr = !(ExprType::uMatAccess & util::MatAccessFlag::ELEMENT_WISE);
	  typedef typename util::TempExprHelper<Expr, bUseTempExpr>::Type TempExpr;
	  typedef typename util::remove_all<TempExpr>::type TempExprType;
     public:
	  static constexpr index_type iNumElemOps = TempExprType::iNumElemOps + 1;
	  static constexpr unsigned uMatAccess = util::MatAccessFlag::ELEMENT_WISE;
	  static constexpr index_type iNumRowsStatic = ExprType::iNumRowsStatic;
	  static constexpr index_type iNumColsStatic = ExprType::iNumColsStatic;
	  static constexpr SpMatOpType eMatOpType = SpMatOpType::MATRIX;
	  static constexpr SpGradCommon::ExprEvalFlags eComprFlag = ExprType::eComprFlag;

	  static_assert(ExprType::eMatOpType == SpMatOpType::MATRIX,
			"Operand must be a matrix! Use SpGradient for scalar expressions!");

	  constexpr SpMatElemUnaryExpr(const Expr& u) noexcept
	       :u(u) {
	  }

	  SpMatElemUnaryExpr(SpMatElemUnaryExpr&& oExpr) noexcept
	       :u(std::move(oExpr.u)) {
	  }

	  template <typename MatEvalType = util::MatEvalHelperTransp<util::MatTranspEvalFlag::DIRECT> >
	  void Eval(SpMatrixBase<ValueType, iNumRowsStatic, iNumColsStatic>& A) const {
	       this->template ElemEval<MatEvalType>(A);
	  }

	  constexpr doublereal dGetValue(index_type i, index_type j) const {
	       return UnaryFunc::f(u.dGetValue(i, j));
	  }

	  constexpr index_type iGetSize(index_type i, index_type j) const {
	       return u.iGetSize(i, j);
	  }

	  constexpr index_type iGetNumRows() const {
	       return u.iGetNumRows();
	  }

	  constexpr index_type iGetNumCols() const {
	       return u.iGetNumCols();
	  }

	  constexpr index_type iGetMaxSize() const {
	       return u.iGetMaxSize();
	  }

	  template <typename ValueTypeB>
	  void InsertDeriv(ValueTypeB& g, doublereal dCoef, index_type i, index_type j) const {
	       u.InsertDeriv(g, UnaryFunc::df_du(u.dGetValue(i, j)) * dCoef, i, j);
	  }

	  void GetDofStat(SpGradDofStat& s, index_type i, index_type j) const {
	       u.GetDofStat(s, i, j);
	  }

	  void InsertDof(SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
	       u.InsertDof(oExpDofMap, i, j);
	  }

	  template <typename ValueTypeB>
	  void AddDeriv(ValueTypeB& g, doublereal dCoef, const SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
	       u.AddDeriv(g, UnaryFunc::df_du(u.dGetValue(i, j)) * dCoef, oExpDofMap, i, j);
	  }

	  template <typename ExprTypeB, typename ExprB>
	  constexpr bool bHaveRefTo(const SpMatElemExprBase<ExprTypeB, ExprB>& A) const {
	       return u.bHaveRefTo(A);
	  }

	  const ValueType* begin() const = delete;
	  const ValueType* end() const = delete;
	  inline constexpr index_type iGetRowOffset() const = delete;
	  inline constexpr index_type iGetColOffset() const = delete;
     private:
	  TempExpr u;
     };

     template <typename ValueType, typename Expr>
     class SpMatElemComprExpr
	  :public SpMatElemExprBase<ValueType, SpMatElemComprExpr<ValueType, Expr> > {
	  typedef typename util::remove_all<Expr>::type ExprType;
     public:
	  static constexpr index_type iNumElemOps = ExprType::iNumElemOps;
	  static constexpr unsigned uMatAccess = ExprType::uMatAccess;
	  static constexpr index_type iNumRowsStatic = ExprType::iNumRowsStatic;
	  static constexpr index_type iNumColsStatic = ExprType::iNumColsStatic;
	  static constexpr SpMatOpType eMatOpType = SpMatOpType::MATRIX;
	  static constexpr SpGradCommon::ExprEvalFlags eComprFlag = SpGradCommon::ExprEvalCompressed;

	  static_assert(ExprType::eMatOpType == SpMatOpType::MATRIX,
			"Operand must be a matrix! Use SpGradient for scalar expressions!");

	  constexpr SpMatElemComprExpr(const Expr& u) noexcept
	       :u(u) {
	  }

	  template <typename MatEvalType = util::MatEvalHelperTransp<util::MatTranspEvalFlag::DIRECT> >
	  void Eval(SpMatrixBase<ValueType, iNumRowsStatic, iNumColsStatic>& A) const {
	       this->template ElemEval<MatEvalType>(A);
	  }

	  constexpr doublereal dGetValue(index_type i, index_type j) const {
	       return u.dGetValue(i, j);
	  }

	  constexpr index_type iGetSize(index_type i, index_type j) const {
	       return u.iGetSize(i, j);
	  }

	  constexpr index_type iGetNumRows() const {
	       return u.iGetNumRows();
	  }

	  constexpr index_type iGetNumCols() const {
	       return u.iGetNumCols();
	  }

	  constexpr index_type iGetMaxSize() const {
	       return u.iGetMaxSize();
	  }

	  template <typename ValueTypeB>
	  void InsertDeriv(ValueTypeB& g, doublereal dCoef, index_type i, index_type j) const {
	       u.InsertDeriv(g, dCoef, i, j);
	  }

	  void GetDofStat(SpGradDofStat& s, index_type i, index_type j) const {
	       u.GetDofStat(s, i, j);
	  }

	  void InsertDof(SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
	       u.InsertDof(oExpDofMap, i, j);
	  }

	  template <typename ValueTypeB>
	  void AddDeriv(ValueTypeB& g, doublereal dCoef, const SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
	       u.AddDeriv(g, dCoef, oExpDofMap, i, j);
	  }

	  template <typename ExprTypeB, typename ExprB>
	  constexpr bool bHaveRefTo(const SpMatElemExprBase<ExprTypeB, ExprB>& A) const {
	       return u.bHaveRefTo(A);
	  }

	  const ValueType* begin() const { return u.begin(); }
	  const ValueType* end() const { return u.end(); }
	  inline constexpr index_type iGetRowOffset() const { return u.iGetRowOffset(); }
	  inline constexpr index_type iGetColOffset() const { return u.iGetColOffset(); }
     private:
	  const Expr u;
     };

     template <typename ValueType, typename Expr>
     class SpMatElemTranspExpr
	  :public SpMatElemExprBase<ValueType, SpMatElemTranspExpr<ValueType, Expr> > {
	  typedef typename util::remove_all<Expr>::type ExprType;
     public:
	  static constexpr index_type iNumElemOps = ExprType::iNumElemOps;
	  static constexpr unsigned uMatAccess = ExprType::uMatAccess;
	  static constexpr index_type iNumRowsStatic = ExprType::iNumColsStatic;
	  static constexpr index_type iNumColsStatic = ExprType::iNumRowsStatic;
	  static constexpr SpMatOpType eMatOpType = SpMatOpType::MATRIX;
	  static constexpr SpGradCommon::ExprEvalFlags eComprFlag = ExprType::eComprFlag;

	  static_assert(ExprType::eMatOpType == SpMatOpType::MATRIX, "Operand must be a matrix! A scalar cannot be transposed!");

	  constexpr SpMatElemTranspExpr(const Expr& u) noexcept
	       :u(u) {
	  }

	  template <typename MatEvalType = util::MatEvalHelperTransp<util::MatTranspEvalFlag::DIRECT>, index_type NumRows, index_type NumCols>
	  void Eval(SpMatrixBase<ValueType, NumRows, NumCols>& A) const {
	       u.template Eval<typename util::TranspMatEvalHelper<MatEvalType>::Type>(A);
	  }

	  constexpr doublereal dGetValue(index_type i, index_type j) const {
	       return u.dGetValue(j, i);
	  }

	  constexpr index_type iGetSize(index_type i, index_type j) const {
	       return u.iGetSize(j, i);
	  }

	  constexpr index_type iGetNumRows() const noexcept {
	       return u.iGetNumCols();
	  }

	  constexpr index_type iGetNumCols() const noexcept {
	       return u.iGetNumRows();
	  }

	  constexpr index_type iGetMaxSize() const noexcept {
	       return u.iGetMaxSize();
	  }

	  template <typename ValueTypeB>
	  void InsertDeriv(ValueTypeB& g, doublereal dCoef, index_type i, index_type j) const {
	       u.InsertDeriv(g, dCoef, j, i);
	  }

	  void GetDofStat(SpGradDofStat& s, index_type i, index_type j) const noexcept {
	       u.GetDofStat(s, j, i);
	  }

	  void InsertDof(SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
	       u.InsertDof(oExpDofMap, j, i);
	  }

	  template <typename ValueTypeB>
	  void AddDeriv(ValueTypeB& g, doublereal dCoef, const SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
	       u.AddDeriv(g, dCoef, oExpDofMap, j, i);
	  }

	  template <typename ExprTypeB, typename ExprB>
	  constexpr bool bHaveRefTo(const SpMatElemExprBase<ExprTypeB, ExprB>& A) const noexcept {
	       return u.bHaveRefTo(A);
	  }

	  const ValueType* begin() const noexcept { return u.begin(); }
	  const ValueType* end() const noexcept { return u.end(); }
	  inline constexpr index_type iGetRowOffset() const noexcept { return u.iGetColOffset(); }
	  inline constexpr index_type iGetColOffset() const noexcept { return u.iGetRowOffset(); }
     private:
	  const Expr u;
     };

     template <typename ValueType, typename Expr>
     class SpMatColVecExpr
	  :public SpMatElemExprBase<ValueType, SpMatColVecExpr<ValueType, Expr> > {
	  typedef typename util::remove_all<Expr>::type ExprType;
     public:
	  static constexpr index_type iNumElemOps = ExprType::iNumElemOps;
	  static constexpr unsigned uMatAccess = ExprType::uMatAccess;
	  static constexpr index_type iNumRowsStatic = ExprType::iNumRowsStatic;
	  static constexpr index_type iNumColsStatic = 1;
	  static constexpr SpMatOpType eMatOpType = SpMatOpType::MATRIX;
	  static constexpr SpGradCommon::ExprEvalFlags eComprFlag = ExprType::eComprFlag;

	  static_assert(ExprType::eMatOpType == SpMatOpType::MATRIX, "Operand must be a matrix!");

	  constexpr SpMatColVecExpr(const Expr& u, index_type iCol) noexcept
	       :u(u), iCol(iCol) {
	  }

	  template <typename MatEvalType = util::MatEvalHelperTransp<util::MatTranspEvalFlag::DIRECT>, index_type NumRows, index_type NumCols>
	  void Eval(SpMatrixBase<ValueType, NumRows, NumCols>& A) const {
	       this->template ElemEval<MatEvalType>(A);
	  }

	  constexpr doublereal dGetValue(index_type i, index_type j) const {
	       SP_GRAD_ASSERT(j == 1);

	       return u.dGetValue(i, iCol);
	  }

	  constexpr index_type iGetSize(index_type i, index_type j) const {
	       SP_GRAD_ASSERT(j == 1);

	       return u.iGetSize(i, iCol);
	  }

	  constexpr index_type iGetNumRows() const noexcept {
	       return u.iGetNumRows();
	  }

	  static constexpr index_type iGetNumCols() noexcept {
	       return 1;
	  }

	  constexpr index_type iGetMaxSize() const noexcept {
	       index_type iMaxSize = 0;

	       for (index_type i = 1; i <= u.iGetNumRows(); ++i) {
		    iMaxSize = std::max(iMaxSize, u.iGetSize(i, iCol));
	       }

	       return iMaxSize;
	  }

	  template <typename ValueTypeB>
	  void InsertDeriv(ValueTypeB& g, doublereal dCoef, index_type i, index_type j) const {
	       SP_GRAD_ASSERT(j == 1);

	       u.InsertDeriv(g, dCoef, i, iCol);
	  }

	  void GetDofStat(SpGradDofStat& s, index_type i, index_type j) const noexcept {
	       SP_GRAD_ASSERT(j == 1);

	       u.GetDofStat(s, i, iCol);
	  }

	  void InsertDof(SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
	       SP_GRAD_ASSERT(j == 1);

	       u.InsertDof(oExpDofMap, i, iCol);
	  }

	  template <typename ValueTypeB>
	  void AddDeriv(ValueTypeB& g, doublereal dCoef, const SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
	       SP_GRAD_ASSERT(j == 1);

	       u.AddDeriv(g, dCoef, oExpDofMap, i, iCol);
	  }

	  template <typename ExprTypeB, typename ExprB>
	  constexpr bool bHaveRefTo(const SpMatElemExprBase<ExprTypeB, ExprB>& A) const noexcept {
	       return u.bHaveRefTo(A);
	  }

	  const ValueType* begin() const { return u.begin() + (iCol - 1) * u.iGetColOffset(); }
	  const ValueType* end() const { return u.begin() + iCol * u.iGetColOffset(); }
	  inline constexpr index_type iGetRowOffset() const { return u.iGetRowOffset(); }
	  inline constexpr index_type iGetColOffset() const { return u.iGetColOffset(); }
     private:
	  const Expr u;
	  const index_type iCol;
     };

     template <typename ValueType, typename Expr>
     class SpMatRowVecExpr
	  :public SpMatElemExprBase<ValueType, SpMatRowVecExpr<ValueType, Expr> > {
	  typedef typename util::remove_all<Expr>::type ExprType;
     public:
	  static constexpr index_type iNumElemOps = ExprType::iNumElemOps;
	  static constexpr unsigned uMatAccess = ExprType::uMatAccess;
	  static constexpr index_type iNumRowsStatic = 1;
	  static constexpr index_type iNumColsStatic = ExprType::iNumColsStatic;
	  static constexpr SpMatOpType eMatOpType = SpMatOpType::MATRIX;
	  static constexpr SpGradCommon::ExprEvalFlags eComprFlag = ExprType::eComprFlag;

	  static_assert(ExprType::eMatOpType == SpMatOpType::MATRIX, "Operand must be a matrix!");

	  constexpr SpMatRowVecExpr(const Expr& u, index_type iRow) noexcept
	       :u(u), iRow(iRow) {
	  }

	  template <typename MatEvalType = util::MatEvalHelperTransp<util::MatTranspEvalFlag::DIRECT>, index_type NumRows, index_type NumCols>
	  void Eval(SpMatrixBase<ValueType, NumRows, NumCols>& A) const {
	       this->template ElemEval<MatEvalType>(A);
	  }

	  constexpr doublereal dGetValue(index_type i, index_type j) const {
	       SP_GRAD_ASSERT(i == 1);

	       return u.dGetValue(iRow, j);
	  }

	  constexpr index_type iGetSize(index_type i, index_type j) const {
	       SP_GRAD_ASSERT(i == 1);

	       return u.iGetSize(iRow, j);
	  }

	  static constexpr index_type iGetNumRows() noexcept {
	       return 1;
	  }

	  constexpr index_type iGetNumCols() const noexcept {
	       return u.iGetNumCols();
	  }

	  constexpr index_type iGetMaxSize() const noexcept {
	       index_type iMaxSize = 0;

	       for (index_type j = 1; j <= u.iGetNumCols(); ++j) {
		    iMaxSize = std::max(iMaxSize, u.iGetSize(iRow, j));
	       }

	       return iMaxSize;
	  }

	  template <typename ValueTypeB>
	  void InsertDeriv(ValueTypeB& g, doublereal dCoef, index_type i, index_type j) const {
	       SP_GRAD_ASSERT(i == 1);

	       u.InsertDeriv(g, dCoef, iRow, j);
	  }

	  void GetDofStat(SpGradDofStat& s, index_type i, index_type j) const noexcept {
	       SP_GRAD_ASSERT(i == 1);

	       u.GetDofStat(s, iRow, j);
	  }

	  void InsertDof(SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
	       SP_GRAD_ASSERT(i == 1);

	       u.InsertDof(oExpDofMap, iRow, j);
	  }

	  template <typename ValueTypeB>
	  void AddDeriv(ValueTypeB& g, doublereal dCoef, const SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
	       SP_GRAD_ASSERT(i == 1);

	       u.AddDeriv(g, dCoef, oExpDofMap, iRow, j);
	  }

	  template <typename ExprTypeB, typename ExprB>
	  constexpr bool bHaveRefTo(const SpMatElemExprBase<ExprTypeB, ExprB>& A) const noexcept {
	       return u.bHaveRefTo(A);
	  }

	  const ValueType* begin() const { return u.begin() + (iRow - 1) * u.iGetRowOffset(); }
	  const ValueType* end() const { return u.begin() + iRow * u.iGetRowOffset(); }
	  inline constexpr index_type iGetRowOffset() const { return u.iGetRowOffset(); }
	  inline constexpr index_type iGetColOffset() const { return u.iGetColOffset(); }

     private:
	  const Expr u;
	  const index_type iRow;
     };

     template <typename LhsValue, typename RhsValue, typename LhsExpr, typename RhsExpr>
     class SpMatMulExpr
	  :public SpMatElemExprBase<typename util::ResultType<LhsValue, RhsValue>::Type,
				    SpMatMulExpr<LhsValue, RhsValue, LhsExpr, RhsExpr> > {
	  typedef typename util::remove_all<LhsExpr>::type LhsExprType;
	  typedef typename util::remove_all<RhsExpr>::type RhsExprType;
	  enum OpCntEstimate: index_type {
	       ASSUMED_OPERATION_COUNT = LhsExprType::iNumColsStatic == SpMatrixSize::DYNAMIC ? 10 : LhsExprType::iNumColsStatic
	  };
     public:
	  typedef typename util::ResultType<LhsValue, RhsValue>::Type ValueType;

	  static constexpr index_type iNumElemOps = (util::remove_all<LhsExpr>::type::iNumElemOps +
						     util::remove_all<RhsExpr>::type::iNumElemOps + 2) * ASSUMED_OPERATION_COUNT;

	  static constexpr unsigned uMatAccess = util::MatAccessFlag::MATRIX_WISE;
	  static constexpr index_type iNumRowsStatic = LhsExprType::iNumRowsStatic;
	  static constexpr index_type iNumColsStatic = RhsExprType::iNumColsStatic;
	  static constexpr SpMatOpType eMatOpType = SpMatOpType::MATRIX;
	  static constexpr SpGradCommon::ExprEvalFlags eComprFlag = util::ExprEvalFlagsHelper<LhsExprType::eComprFlag,
											      RhsExprType::eComprFlag>::eExprEvalFlags;

	  static_assert(LhsExprType::eMatOpType == SpMatOpType::MATRIX, "Left hand side of matrix product must be a matrix");
	  static_assert(RhsExprType::eMatOpType == SpMatOpType::MATRIX, "Right hand side of matrix product must be a matrix");
	  static_assert(LhsExprType::iNumColsStatic == RhsExprType::iNumRowsStatic, "Incompatible matrix sizes in matrix product");

	  SpMatMulExpr(const LhsExpr& u, const RhsExpr& v) noexcept
	       :u(u), v(v) {
	       SP_GRAD_ASSERT(u.iGetNumCols() == v.iGetNumRows());
	  }

	  template <typename MatEvalType = util::MatEvalHelperTransp<util::MatTranspEvalFlag::DIRECT>, index_type NumRowsA, index_type NumColsA>
	  inline void Eval(SpMatrixBase<ValueType, NumRowsA, NumColsA>& A) const;

	  doublereal dGetValue(index_type, index_type) const = delete;

	  index_type iGetSize(index_type i, index_type j) const = delete;

	  constexpr index_type iGetNumRows() const {
	       return u.iGetNumRows();
	  }

	  constexpr index_type iGetNumCols() const {
	       return v.iGetNumCols();
	  }

	  constexpr index_type iGetMaxSize() const = delete;

	  template <typename ValueTypeB>
	  void InsertDeriv(ValueTypeB& g, doublereal dCoef, index_type i, index_type j) const = delete;

	  void GetDofStat(SpGradDofStat& s, index_type i, index_type j) const = delete;

	  void InsertDof(SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const = delete;

	  template <typename ValueTypeB>
	  void AddDeriv(ValueTypeB& g, doublereal dCoef, const SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const = delete;

	  template <typename ExprType, typename Expr>
	  constexpr bool bHaveRefTo(const SpMatElemExprBase<ExprType, Expr>& A) const {
	       return u.bHaveRefTo(A) || v.bHaveRefTo(A);
	  }

	  const ValueType* begin() const = delete;
	  const ValueType* end() const = delete;
	  inline constexpr index_type iGetRowOffset() const = delete;
	  inline constexpr index_type iGetColOffset() const = delete;
     private:
	  const LhsExpr u;
	  const RhsExpr v;
     };

     namespace util {
	  template <typename Expr>
	  struct ScalarExprEvalFlagsHelper {
	       static constexpr SpGradCommon::ExprEvalFlags eComprFlag = Expr::eExprEvalFlags;
	  };

	  template <>
	  struct ScalarExprEvalFlagsHelper<doublereal> {
	       static constexpr SpGradCommon::ExprEvalFlags eComprFlag = SpGradCommon::ExprEvalUncompressed;
	  };

     }

     template <typename ValueType, typename ScalarExpr>
     class SpMatElemScalarExpr: public SpMatElemExprBase<ValueType, SpMatElemScalarExpr<ValueType, ScalarExpr> > {
	  enum OpCntEstimate: index_type {
	       EST_NUMBER_OF_SCALAR_OPS = 1
	  };
     public:
	  static constexpr index_type iNumElemOps = EST_NUMBER_OF_SCALAR_OPS;
	  static constexpr unsigned uMatAccess = util::MatAccessFlag::ELEMENT_WISE;
	  static constexpr index_type iNumRowsStatic = 1;
	  static constexpr index_type iNumColsStatic = 1;
	  static constexpr SpMatOpType eMatOpType = SpMatOpType::SCALAR;
	  static constexpr SpGradCommon::ExprEvalFlags eComprFlag = util::ScalarExprEvalFlagsHelper<ScalarExpr>::eComprFlag;

	  inline constexpr SpMatElemScalarExpr(const ScalarExpr& u) noexcept;

	  template <typename MatEvalType = util::MatEvalHelperTransp<util::MatTranspEvalFlag::DIRECT> >
	  void Eval(SpMatrixBase<ValueType, iNumRowsStatic, iNumColsStatic>& A) const {
	       this->template ElemEval<MatEvalType>(A);
	  }

	  inline constexpr doublereal dGetValue(index_type, index_type) const noexcept;

	  inline constexpr index_type iGetSize(index_type, index_type) const noexcept;

	  inline static constexpr index_type iGetNumRows() noexcept;

	  inline static constexpr index_type iGetNumCols() noexcept;

	  inline constexpr index_type iGetMaxSize() const noexcept;

	  template <typename ValueTypeB>
	  inline void InsertDeriv(ValueTypeB& g, doublereal dCoef, index_type, index_type) const;

	  inline void GetDofStat(SpGradDofStat& s, index_type, index_type) const noexcept;

	  inline void InsertDof(SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const;

	  template <typename ValueTypeB>
	  inline void AddDeriv(ValueTypeB& g, doublereal dCoef, const SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const;

	  template <typename ExprTypeB, typename ExprB>
	  inline constexpr bool bHaveRefTo(const SpMatElemExprBase<ExprTypeB, ExprB>& A) const noexcept;

	  const ValueType* begin() const = delete;
	  const ValueType* end() const = delete;
	  inline constexpr index_type iGetRowOffset() const = delete;
	  inline constexpr index_type iGetColOffset() const = delete;
     private:
	  const ScalarExpr u;
     };

     template <typename ValueType, index_type NumRows, index_type NumCols>
     class SpMatrixBase: public SpMatElemExprBase<ValueType, SpMatrixBase<ValueType, NumRows, NumCols> > {
	  static_assert(NumRows >= 1 || NumRows == SpMatrixSize::DYNAMIC);
	  static_assert(NumCols >= 1 || NumCols == SpMatrixSize::DYNAMIC);
	  friend class SpMatrixBase<const ValueType, NumRows, NumCols>;
     public:
	  static constexpr index_type iNumElemOps = 0;
	  static constexpr unsigned uMatAccess = util::MatAccessFlag::ELEMENT_WISE |
	       util::MatAccessFlag::ITERATORS |
	       util::MatAccessFlag::MATRIX_WISE;
	  static constexpr index_type iNumRowsStatic = NumRows;
	  static constexpr index_type iNumColsStatic = NumCols;
	  static constexpr SpMatOpType eMatOpType = SpMatOpType::MATRIX;
	  static constexpr SpGradCommon::ExprEvalFlags eComprFlag = SpGradCommon::ExprEvalUncompressed;

	  inline SpMatrixBase();
	  inline SpMatrixBase(const SpMatrixBase& oMat);
	  inline SpMatrixBase(SpMatrixBase&& oMat);
	  inline SpMatrixBase(index_type iNumRows, index_type iNumCols, index_type iNumDeriv);
	  template <typename ValueTypeExpr, typename Expr>
	  inline SpMatrixBase(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr);
	  inline ~SpMatrixBase();
	  inline void ResizeReset(index_type iNumRows, index_type iNumCols, index_type iNumDeriv);
	  inline SpMatrixBase& operator=(SpMatrixBase&& oMat);
	  template <typename ValueTypeExpr, typename Expr>
	  inline SpMatrixBase& operator=(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr);
	  inline SpMatrixBase& operator=(const SpMatrixBase& oMat);
	  inline index_type iGetNumRows() const;
	  inline index_type iGetNumCols() const;
	  inline index_type iGetSize(index_type i, index_type j) const;
	  inline doublereal dGetValue(index_type i, index_type j) const;
	  template <typename ValueType_B>
	  inline void InsertDeriv(ValueType_B& g, doublereal dCoef, index_type i, index_type j) const;
	  inline void GetDofStat(SpGradDofStat& s, index_type i, index_type j) const;
	  inline void InsertDof(SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const;
	  template <typename ValueTypeB>
	  inline void AddDeriv(ValueTypeB& g, doublereal dCoef, const SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const;
	  template <typename ExprType, typename Expr>
	  inline constexpr bool bHaveRefTo(const SpMatElemExprBase<ExprType, Expr>& A) const noexcept;
	  inline const ValueType& GetElem(index_type i) const;
	  inline ValueType& GetElem(index_type i);
	  inline const ValueType& GetElem(index_type i, index_type j) const;
	  inline ValueType& GetElem(index_type i, index_type j);
	  inline ValueType* begin();
	  inline ValueType* end();
	  inline const ValueType* begin() const;
	  inline const ValueType* end() const;
	  static inline constexpr index_type iGetRowOffset() noexcept { return 1; }
	  inline constexpr index_type iGetColOffset() const noexcept { return iGetNumRows(); }
	  inline bool bHaveSameRep(const SpMatrixBase& A) const noexcept;
	  template <typename ValueTypeB>
	  constexpr static inline bool bIsOwnerOf(const SpMatrixData<ValueTypeB>*) noexcept { return false; }
	  inline bool bIsOwnerOf(const SpMatrixData<ValueType>* pDataB) const noexcept;
	  template <typename MatEvalType = util::MatEvalHelperTransp<util::MatTranspEvalFlag::DIRECT>,
		    index_type NumRowsA, index_type NumColsA>
	  inline void Eval(SpMatrixBase<ValueType, NumRowsA, NumColsA>& A) const;
	  inline index_type iGetMaxSize() const;
	  inline constexpr SpMatColVecExpr<ValueType, const SpMatrixBase&> GetCol(index_type iCol) const {
	       return decltype(GetCol(iCol))(*this, iCol);
	  }
	  inline constexpr SpMatRowVecExpr<ValueType, const SpMatrixBase&> GetRow(index_type iRow) const {
	       return decltype(GetRow(iRow))(*this, iRow);
	  }
     private:
	  typedef SpMatrixDataCTD<typename std::remove_cv<ValueType>::type, iNumRowsStatic, iNumColsStatic> SpMatrixDataType;

	  explicit inline SpMatrixBase(SpMatrixDataType* pData);

	  static inline SpMatrixDataType* pGetNullData();

	  static constexpr bool bStaticSize = iNumRowsStatic != SpMatrixSize::DYNAMIC && iNumColsStatic != SpMatrixSize::DYNAMIC;
	  SpMatrixDataType* pData;
     };

     template <typename ValueType, index_type NumRows = SpMatrixSize::DYNAMIC, index_type NumCols = SpMatrixSize::DYNAMIC>
     class SpMatrix: public SpMatrixBase<ValueType, NumRows, NumCols> {
     public:
	  inline SpMatrix()=default;
	  inline SpMatrix(index_type iNumRows, index_type iNumCols, index_type iNumDeriv);
	  inline SpMatrix(const SpMatrix& oMat)=default;
	  inline SpMatrix(SpMatrix&& oMat)=default;
	  template <typename ValueTypeExpr, typename Expr>
	  inline SpMatrix(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr);
	  inline SpMatrix& operator=(SpMatrix&& oMat)=default;
	  inline SpMatrix& operator=(const SpMatrix& oMat)=default;
	  template <typename ValueTypeExpr, typename Expr>
	  inline SpMatrix& operator=(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr);
	  inline ValueType& operator()(index_type iRow, index_type iCol);
	  inline const ValueType& operator() (index_type iRow, index_type iCol) const;
     };

     template <typename ValueType, index_type NumRows = SpMatrixSize::DYNAMIC>
     class SpColVector: public SpMatrixBase<ValueType, NumRows, 1> {
     public:
	  inline SpColVector()=default;
	  inline SpColVector(index_type iNumRows, index_type iNumDeriv);
	  inline SpColVector(const SpColVector& oVec)=default;
	  inline SpColVector(SpColVector&& oVec)=default;
	  template <typename ValueTypeExpr, typename Expr>
	  inline SpColVector(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr);
	  inline SpColVector& operator=(const SpColVector& oVec)=default;
	  inline SpColVector& operator=(SpColVector&& oVec)=default;
	  template <typename ValueTypeExpr, typename Expr>
	  inline SpColVector& operator=(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr);
	  inline ValueType& operator()(index_type iRow);
	  inline const ValueType& operator() (index_type iRow) const;
     };

     template <typename ValueType, index_type NumCols = SpMatrixSize::DYNAMIC>
     class SpRowVector: public SpMatrixBase<ValueType, 1, NumCols> {
     public:
	  inline SpRowVector()=default;
	  inline SpRowVector(index_type iNumCols, index_type iNumDeriv);
	  inline SpRowVector(const SpRowVector& oVec)=default;
	  inline SpRowVector(SpRowVector&& oVec)=default;
	  template <typename ValueTypeExpr, typename Expr>
	  inline SpRowVector(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr);
	  inline SpRowVector& operator=(const SpRowVector& oVec)=default;
	  inline SpRowVector& operator=(SpRowVector&& oVec)=default;
	  template <typename ValueTypeExpr, typename Expr>
	  inline SpRowVector& operator=(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr);
	  inline ValueType& operator()(index_type iCol);
	  inline const ValueType& operator() (index_type iCol) const;
     };
}
#endif
