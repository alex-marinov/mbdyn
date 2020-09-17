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

#include "matvec3.h"
#include "matvec3n.h"
#include "RotCoeff.hh"
#include "Rot.hh"

#include "sp_matrix_base_fwd.h"
#include "sp_gradient_func.h"
#include "sp_gradient.h"

namespace sp_grad {
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

	  template <typename ValueType, SpGradCommon::ExprEvalFlags eCompr>
	  struct ComprEvalHelper {
	       static constexpr SpGradCommon::ExprEvalFlags eComprFlag = eCompr;
	  };

	  template <SpGradCommon::ExprEvalFlags eCompr>
	  struct ComprEvalHelper<doublereal, eCompr> {
	       // It makes no sense to allocate a class SpGradExprDofMap for a doublereal
	       static const SpGradCommon::ExprEvalFlags eComprFlag = SpGradCommon::ExprEvalUncompressed;
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
	  static constexpr index_type iNumRowsStatic = LhsExprType::eMatOpType == SpMatOpType::MATRIX ? LhsExprType::iNumRowsStatic : RhsExprType::iNumRowsStatic;
	  static constexpr index_type iNumColsStatic = LhsExprType::eMatOpType == SpMatOpType::MATRIX ? LhsExprType::iNumColsStatic : RhsExprType::iNumColsStatic;
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

	  template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
		    SpGradCommon::ExprEvalFlags eCompr = eComprFlag,
		    index_type iNumRowsA, index_type iNumColsA>
	  void Eval(SpMatrixBase<ValueType, iNumRowsA, iNumColsA>& A) const {
	       this->template ElemEval<eTransp, eCompr>(A);
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

     template <typename ValueType, typename LhsExpr, typename RhsExpr>
     class SpMatCrossExpr
	  :public SpMatElemExprBase<ValueType, SpMatCrossExpr<ValueType, LhsExpr, RhsExpr> > {
	  typedef typename util::remove_all<LhsExpr>::type LhsExprType;
	  typedef typename util::remove_all<RhsExpr>::type RhsExprType;
	  static constexpr bool bUseTempExprLhs = !(LhsExprType::uMatAccess & util::MatAccessFlag::ELEMENT_WISE) || LhsExprType::iNumElemOps > 0;
	  static constexpr bool bUseTempExprRhs = !(RhsExprType::uMatAccess & util::MatAccessFlag::ELEMENT_WISE) || RhsExprType::iNumElemOps > 0;
	  typedef typename util::TempExprHelper<LhsExpr, bUseTempExprLhs>::Type LhsTempExpr;
	  typedef typename util::TempExprHelper<RhsExpr, bUseTempExprRhs>::Type RhsTempExpr;
	  typedef typename util::remove_all<LhsTempExpr>::type LhsTempExprType;
	  typedef typename util::remove_all<RhsTempExpr>::type RhsTempExprType;
	  typedef typename LhsTempExprType::ValueType LhsValueType;
	  typedef typename RhsTempExprType::ValueType RhsValueType;
	  typedef typename LhsTempExprType::DerivedType LhsDerivedType;
	  typedef typename RhsTempExprType::DerivedType RhsDerivedType;
	  typedef util::MatrixSizeHelper<LhsValueType, LhsDerivedType, RhsValueType, RhsDerivedType> MatSizeHelp;

	  static constexpr index_type x = 1, y = 2, z = 3;
	  static constexpr index_type e1[3] = {y, z, x};
	  static constexpr index_type e2[3] = {z, x, y};
	  static constexpr index_type e3[3] = {z, x, y};
	  static constexpr index_type e4[3] = {y, z, x};
     public:
	  static constexpr index_type iNumElemOps = LhsTempExprType::iNumElemOps + RhsTempExprType::iNumElemOps + 1;
	  static constexpr unsigned uMatAccess = util::MatAccessFlag::ELEMENT_WISE;
	  static constexpr index_type iNumRowsStatic = 3;
	  static constexpr index_type iNumColsStatic = 1;
	  static constexpr SpMatOpType eMatOpType = SpMatOpType::MATRIX;
	  static constexpr SpGradCommon::ExprEvalFlags eComprFlag = util::ExprEvalFlagsHelper<LhsExprType::eComprFlag,
											      RhsExprType::eComprFlag>::eExprEvalFlags;
	  static_assert(LhsExprType::eMatOpType == SpMatOpType::MATRIX && RhsExprType::eMatOpType == SpMatOpType::MATRIX,
			"Both operands must be vectors!");
	  static_assert(LhsExprType::iNumRowsStatic == 3 && RhsExprType::iNumRowsStatic == 3,
			"Both operands must be 3x1 vectors");
	  static_assert(LhsExprType::iNumColsStatic == 1 && RhsExprType::iNumColsStatic == 1,
			"Both operands must be 3x1 vectors");

	  constexpr SpMatCrossExpr(const LhsExprType& u, const RhsExprType& v) noexcept
	       :u(u), v(v) {
	  }

	  SpMatCrossExpr(SpMatCrossExpr&& oExpr)
	       :u(std::move(oExpr.u)), v(std::move(oExpr.v)) {
	  }

	  template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
		    SpGradCommon::ExprEvalFlags eCompr = eComprFlag,
		    index_type iNumRowsA, index_type iNumColsA>
	  void Eval(SpMatrixBase<ValueType, iNumRowsA, iNumColsA>& A) const {
	       this->template ElemEval<eTransp, eCompr>(A);
	  }

	  constexpr doublereal dGetValue(index_type i, index_type j) const {
	       SP_GRAD_ASSERT(j == 1);
	       SP_GRAD_ASSERT(i >= 1);
	       SP_GRAD_ASSERT(i <= 3);

	       static_assert(sizeof(e1) / sizeof(e1[0]) == 3);
	       static_assert(sizeof(e2) / sizeof(e2[0]) == 3);
	       static_assert(sizeof(e3) / sizeof(e3[0]) == 3);
	       static_assert(sizeof(e4) / sizeof(e4[0]) == 3);

	       --i;
	       
	       return u.dGetValue(e1[i], j) * v.dGetValue(e2[i], j) - u.dGetValue(e3[i], j) * v.dGetValue(e4[i], j);
	  }

	  constexpr index_type iGetSize(index_type i, index_type j) const {
	       SP_GRAD_ASSERT(j == 1);
	       SP_GRAD_ASSERT(i >= 1);
	       SP_GRAD_ASSERT(i <= 3);

	       static_assert(sizeof(e1) / sizeof(e1[0]) == 3);
	       static_assert(sizeof(e2) / sizeof(e2[0]) == 3);
	       static_assert(sizeof(e3) / sizeof(e3[0]) == 3);
	       static_assert(sizeof(e4) / sizeof(e4[0]) == 3);

	       --i;
	       
	       return u.iGetSize(e1[i], j)
		    + v.iGetSize(e2[i], j)
		    + u.iGetSize(e3[i], j)
		    + v.iGetSize(e4[i], j);
	  }

	  constexpr index_type iGetNumRows() const {
	       return iNumRowsStatic;
	  }

	  constexpr index_type iGetNumCols() const {
	       return iNumColsStatic;
	  }

	  constexpr index_type iGetMaxSize() const {
	       return this->iGetMaxSizeElem();
	  }

	  template <typename ValueTypeB>
	  void InsertDeriv(ValueTypeB& g, doublereal dCoef, index_type i, index_type j) const {
	       SP_GRAD_ASSERT(j == 1);
	       SP_GRAD_ASSERT(i >= 1);
	       SP_GRAD_ASSERT(i <= 3);

	       static_assert(sizeof(e1) / sizeof(e1[0]) == 3);
	       static_assert(sizeof(e2) / sizeof(e2[0]) == 3);
	       static_assert(sizeof(e3) / sizeof(e3[0]) == 3);
	       static_assert(sizeof(e4) / sizeof(e4[0]) == 3);
	       
	       --i;
	       
	       u.InsertDeriv(g, dCoef * v.dGetValue(e2[i], j), e1[i], j);
	       v.InsertDeriv(g, dCoef * u.dGetValue(e1[i], j), e2[i], j);
	       u.InsertDeriv(g, -dCoef * v.dGetValue(e4[i], j), e3[i], j);
	       v.InsertDeriv(g, -dCoef * u.dGetValue(e3[i], j), e4[i], j);
	  }

	  void GetDofStat(SpGradDofStat& s, index_type i, index_type j) const {
	       SP_GRAD_ASSERT(j == 1);
	       SP_GRAD_ASSERT(i >= 1);
	       SP_GRAD_ASSERT(i <= 3);

	       static_assert(sizeof(e1) / sizeof(e1[0]) == 3);
	       static_assert(sizeof(e2) / sizeof(e2[0]) == 3);
	       static_assert(sizeof(e3) / sizeof(e3[0]) == 3);
	       static_assert(sizeof(e4) / sizeof(e4[0]) == 3);
	       
	       --i;
	       
	       u.GetDofStat(s, e1[i], j);
	       v.GetDofStat(s, e2[i], j);
	       u.GetDofStat(s, e3[i], j);
	       v.GetDofStat(s, e4[i], j);
	  }

	  void InsertDof(SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
	       SP_GRAD_ASSERT(j == 1);
	       SP_GRAD_ASSERT(i >= 1);
	       SP_GRAD_ASSERT(i <= 3);

	       static_assert(sizeof(e1) / sizeof(e1[0]) == 3);
	       static_assert(sizeof(e2) / sizeof(e2[0]) == 3);
	       static_assert(sizeof(e3) / sizeof(e3[0]) == 3);
	       static_assert(sizeof(e4) / sizeof(e4[0]) == 3);
	       
	       --i;
	       
	       u.InsertDof(oExpDofMap, e1[i], j);
	       v.InsertDof(oExpDofMap, e2[i], j);
	       u.InsertDof(oExpDofMap, e3[i], j);
	       v.InsertDof(oExpDofMap, e4[i], j);
	  }

	  template <typename ValueTypeB>
	  void AddDeriv(ValueTypeB& g, doublereal dCoef, const SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
	       SP_GRAD_ASSERT(j == 1);
	       SP_GRAD_ASSERT(i >= 1);
	       SP_GRAD_ASSERT(i <= 3);

	       static_assert(sizeof(e1) / sizeof(e1[0]) == 3);
	       static_assert(sizeof(e2) / sizeof(e2[0]) == 3);
	       static_assert(sizeof(e3) / sizeof(e3[0]) == 3);
	       static_assert(sizeof(e4) / sizeof(e4[0]) == 3);
	       
	       --i;
	       
	       u.AddDeriv(g, dCoef * v.dGetValue(e2[i], j), oExpDofMap, e1[i], j);
	       v.AddDeriv(g, dCoef * u.dGetValue(e1[i], j), oExpDofMap, e2[i], j);
	       u.AddDeriv(g, -dCoef * v.dGetValue(e4[i], j), oExpDofMap, e3[i], j);
	       v.AddDeriv(g, -dCoef * u.dGetValue(e3[i], j), oExpDofMap, e4[i], j);
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

     template <typename ValueType, typename LhsExpr, typename RhsExpr>
     constexpr index_type SpMatCrossExpr<ValueType, LhsExpr, RhsExpr>::e1[3];
     template <typename ValueType, typename LhsExpr, typename RhsExpr>
     constexpr index_type SpMatCrossExpr<ValueType, LhsExpr, RhsExpr>::e2[3];
     template <typename ValueType, typename LhsExpr, typename RhsExpr>
     constexpr index_type SpMatCrossExpr<ValueType, LhsExpr, RhsExpr>::e3[3];
     template <typename ValueType, typename LhsExpr, typename RhsExpr>
     constexpr index_type SpMatCrossExpr<ValueType, LhsExpr, RhsExpr>::e4[3];

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

	  constexpr explicit SpMatElemUnaryExpr(const Expr& u) noexcept
	       :u(u) {
	  }

	  SpMatElemUnaryExpr(SpMatElemUnaryExpr&& oExpr) noexcept
	       :u(std::move(oExpr.u)) {
	  }

	  template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
		    SpGradCommon::ExprEvalFlags eCompr = eComprFlag>
	  void Eval(SpMatrixBase<ValueType, iNumRowsStatic, iNumColsStatic>& A) const {
	       this->template ElemEval<eTransp, eCompr>(A);
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
	  static constexpr SpMatOpType eMatOpType = ExprType::eMatOpType;
	  static constexpr SpGradCommon::ExprEvalFlags eComprFlag = util::ComprEvalHelper<ValueType, SpGradCommon::ExprEvalCompressed>::eComprFlag;

	  static_assert(ExprType::eMatOpType == SpMatOpType::MATRIX,
			"Operand must be a matrix! Use SpGradient for scalar expressions!");

	  constexpr explicit SpMatElemComprExpr(const Expr& u) noexcept
	       :u(u) {
	  }

	  template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
		    SpGradCommon::ExprEvalFlags eCompr = eComprFlag>
	  void Eval(SpMatrixBase<ValueType, iNumRowsStatic, iNumColsStatic>& A) const {
	       u.template Eval<eTransp, eComprFlag>(A);
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

	  constexpr explicit SpMatElemTranspExpr(const Expr& u) noexcept
	       :u(u) {
	  }

	  template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
		    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalUncompressed,
		    index_type NumRows, index_type NumCols>
	  void Eval(SpMatrixBase<ValueType, NumRows, NumCols>& A) const {
	       u.template Eval<util::TranspMatEvalHelper<eTransp>::Value, eCompr>(A);
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

	  template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
		    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalUncompressed,
		    index_type NumRows, index_type NumCols>
	  void Eval(SpMatrixBase<ValueType, NumRows, NumCols>& A) const {
	       this->template ElemEval<eTransp, eCompr>(A);
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

	  template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
		    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalUncompressed,
		    index_type NumRows, index_type NumCols>
	  void Eval(SpMatrixBase<ValueType, NumRows, NumCols>& A) const {
	       this->template ElemEval<eTransp, eCompr>(A);
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

	  template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
		    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalUncompressed,
		    index_type NumRowsA, index_type NumColsA>
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

	  inline constexpr explicit SpMatElemScalarExpr(const ScalarExpr& u) noexcept;

	  template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
		    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalUncompressed>
	  void Eval(SpMatrixBase<ValueType, iNumRowsStatic, iNumColsStatic>& A) const {
	       this->template ElemEval<eTransp, eCompr>(A);
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
	  inline SpMatrixBase& operator*=(const ValueType& b);
	  inline SpMatrixBase& operator/=(const ValueType& b);
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
	  template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
		    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalUncompressed,
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
	  inline SpMatrix& operator=(const Mat3x3& oMat);
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
	  inline SpColVector& operator=(const Vec3& oVec);
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

     template <typename ValueType, index_type NumRows, index_type NumCols>
     inline std::ostream& operator<<(std::ostream& os, const SpMatrixBase<ValueType, NumRows, NumCols>& A);
     
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
	  template <MatTranspEvalFlag eTransp, typename Expr, typename ValueType, index_type NumRows, index_type NumCols>
	  void MatEvalHelperCompr<SpGradCommon::ExprEvalUncompressed>::ElemEval(const Expr& oExpr, SpMatrixBase<ValueType, NumRows, NumCols>& A) {
	       static_assert(eComprFlag == SpGradCommon::ExprEvalUncompressed);

	       oExpr.template ElemEvalUncompr<eTransp>(A);
	  }

	  template <MatTranspEvalFlag eTransp, typename Expr, typename ValueType, index_type NumRows, index_type NumCols>
	  void MatEvalHelperCompr<SpGradCommon::ExprEvalCompressed>::ElemEval(const Expr& oExpr, SpMatrixBase<ValueType, NumRows, NumCols>& A) {
	       static_assert(eComprFlag == SpGradCommon::ExprEvalCompressed);

	       oExpr.template ElemEvalCompr<eTransp>(A);
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
     };

     template <typename ValueType, typename DERIVED>
     template <util::MatTranspEvalFlag eTransp, SpGradCommon::ExprEvalFlags eCompr, index_type NumRows, index_type NumCols>
     void SpMatElemExprBase<ValueType, DERIVED>::ElemEval(SpMatrixBase<ValueType, NumRows, NumCols>& A) const {
	  util::MatEvalHelperCompr<eCompr>::template ElemEval<eTransp>(*this, A);
     }

     template <typename ValueType, typename DERIVED>
     template <util::MatTranspEvalFlag eTransp, index_type NumRows, index_type NumCols>
     void SpMatElemExprBase<ValueType, DERIVED>::ElemEvalUncompr(SpMatrixBase<ValueType, NumRows, NumCols>& A) const {
	  const index_type iNumRows = iGetNumRows();
	  const index_type iNumCols = iGetNumCols();
	  const index_type iMaxSize = iGetMaxSize();

	  typedef util::MatEvalHelperTransp<eTransp> MatEvalType;
	  
	  constexpr bool bTransposed = eTransp == util::MatTranspEvalFlag::TRANSPOSED;

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
     template <util::MatTranspEvalFlag eTransp, index_type NumRows, index_type NumCols>
     void SpMatElemExprBase<ValueType, DERIVED>::ElemEvalCompr(SpMatrixBase<ValueType, NumRows, NumCols>& A) const {
	  const index_type iNumRows = iGetNumRows();
	  const index_type iNumCols = iGetNumCols();

	  typedef util::MatEvalHelperTransp<eTransp> MatEvalType;
	  constexpr bool bTransposed = eTransp == util::MatTranspEvalFlag::TRANSPOSED;

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
		    
		    SP_GRAD_ASSERT(SpGradient::bIsUnique(Ai));
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


	  template <util::MatTranspEvalFlag eTransp, bool bIsGradientLhs, bool bIsGradientRhs>
	  struct MatMulExprLoop;

	  template <util::MatTranspEvalFlag eTransp>
	  struct MatMulExprLoop<eTransp, true, true> {
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
		    typedef MatEvalHelperTransp<eTransp> MatEvalType;

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

	  template <util::MatTranspEvalFlag eTransp>
	  struct MatMulExprLoop<eTransp, true, false> {
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
		    typedef MatEvalHelperTransp<eTransp> MatEvalType;
		    
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

	  template <util::MatTranspEvalFlag eTransp>
	  struct MatMulExprLoop<eTransp, false, true> {
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
		    typedef MatEvalHelperTransp<eTransp> MatEvalType;
		    
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

	  template <util::MatTranspEvalFlag eTransp>
	  struct MatMulExprLoop<eTransp, false, false> {
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

		    typedef util::MatEvalHelperTransp<eTransp> MatEvalType;
		    
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
     template <util::MatTranspEvalFlag eTransp, SpGradCommon::ExprEvalFlags eCompr, index_type NumRowsA, index_type NumColsA>
     void SpMatMulExpr<LhsValue, RhsValue, LhsExpr, RhsExpr>::Eval(SpMatrixBase<ValueType, NumRowsA, NumColsA>& A) const {
	  typedef util::MatEvalHelperTransp<eTransp> MatEvalType;
	  
	  constexpr bool bTransposedEval = eTransp == util::MatTranspEvalFlag::TRANSPOSED;

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

	  typedef util::MatMulExprLoop<eTransp, bIsGradientLhs, bIsGradientRhs> MatMulExprLoop;

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
	  oExpr.template Eval<util::MatTranspEvalFlag::DIRECT, Expr::eComprFlag>(*this);

#ifdef SP_GRAD_DEBUG
	  if (Expr::eComprFlag == SpGradCommon::ExprEvalCompressed) {
	       for (const auto& a: *this) {
		    SP_GRAD_ASSERT(SpGradient::bIsUnique(a));
	       }
	  }
#endif
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
	  oExpr.template Eval<util::MatTranspEvalFlag::DIRECT, Expr::eComprFlag>(*this);

#ifdef SP_GRAD_DEBUG
	  if (Expr::eComprFlag == SpGradCommon::ExprEvalCompressed) {
	       for (const auto& a: *this) {
		    SP_GRAD_ASSERT(SpGradient::bIsUnique(a));
	       }
	  }
#endif
	  
	  return *this;
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     SpMatrixBase<ValueType, NumRows, NumCols>& SpMatrixBase<ValueType, NumRows, NumCols>::operator*=(const ValueType& b) {
	  for (auto& a:*this) {
	       a *= b;
	  }

	  return *this;
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     SpMatrixBase<ValueType, NumRows, NumCols>& SpMatrixBase<ValueType, NumRows, NumCols>::operator/=(const ValueType& b) {
	  for (auto& a:*this) {
	       a /= b;
	  }

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
	  template <typename ValueType, typename ValueTypeExpr, typename Expr, index_type NumRows, index_type NumCols>
	  struct MatrixBaseRefHelper {
	       static constexpr bool bHaveRefTo(const SpMatrixBase<ValueType, NumRows, NumCols>& A, const SpMatElemExprBase<ValueTypeExpr, Expr>& B) {
		    return false;
	       }
	  };

	  template <typename ValueType, index_type NumRows, index_type NumCols>
	  struct MatrixBaseRefHelper<ValueType, ValueType, SpMatrixBase<ValueType, NumRows, NumCols>, NumRows, NumCols> {
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
	  return util::MatrixBaseRefHelper<ValueType, ExprType, Expr, NumRows, NumCols>::bHaveRefTo(*this, A);
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
     template <util::MatTranspEvalFlag eTransp, SpGradCommon::ExprEvalFlags eCompr, index_type NumRowsA, index_type NumColsA>
     void SpMatrixBase<ValueType, NumRows, NumCols>::Eval(SpMatrixBase<ValueType, NumRowsA, NumColsA>& A) const {
	  this->template ElemEval<eTransp, eCompr>(A);
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

     template <typename ValueType, index_type NumRows, index_type NumCols>
     SpMatrix<ValueType, NumRows, NumCols>& SpMatrix<ValueType, NumRows, NumCols>::operator=(const Mat3x3& oMat) {
	  SpMatrixBase<ValueType, NumRows, NumCols>::operator=(oMat);
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

     template <typename ValueType, index_type NumRows>
     SpColVector<ValueType, NumRows>& SpColVector<ValueType, NumRows>::operator=(const Vec3& oVec) {
	  SpMatrixBase<ValueType, NumRows, 1>::operator=(oVec);
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
     
     template <typename LhsValue, typename RhsValue, typename LhsExpr, typename RhsExpr>
     inline constexpr
     SpMatCrossExpr<typename util::ResultType<LhsValue, RhsValue>::Type,
		    const SpMatElemExprBase<LhsValue, LhsExpr>&,
		    const SpMatElemExprBase<RhsValue, RhsExpr>&>
     Cross(const SpMatElemExprBase<LhsValue, LhsExpr>& A,
	   const SpMatElemExprBase<RhsValue, RhsExpr>& B) noexcept {
	  return decltype(Cross(A, B))(A, B);
     }

     template <typename Value, typename Expr>
     inline constexpr
     SpMatElemTranspExpr<Value, const SpMatElemExprBase<Value, Expr>&>
     Transpose(const SpMatElemExprBase<Value, Expr>& A) noexcept {
	  return decltype(Transpose(A)){A};
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
	  return decltype(EvalCompressed(A)){A};
     }

     template <typename LhsValue, typename LhsExpr, typename RhsExpr>
     inline constexpr
     SpMatElemBinExpr<typename util::ResultType<LhsValue, SpGradient>::Type,
		      SpGradBinMult,
		      const SpMatElemExprBase<LhsValue, LhsExpr>&,
		      const SpMatElemScalarExpr<SpGradient, SpGradient> >
     operator*(const SpMatElemExprBase<LhsValue, LhsExpr>& A,
	       const SpGradBase<RhsExpr>& b) noexcept {
	  return decltype(operator*(A, b)){A, SpMatElemScalarExpr<SpGradient, SpGradient>{SpGradient{b}}}; // Avoid multiple evaluations of b!
     }
     
     template <typename LhsValue, typename LhsExpr>
     inline constexpr
     SpMatElemBinExpr<typename util::ResultType<LhsValue, doublereal>::Type,
		      SpGradBinMult,
		      const SpMatElemExprBase<LhsValue, LhsExpr>&,
		      const SpMatElemScalarExpr<doublereal, doublereal> >
     operator*(const SpMatElemExprBase<LhsValue, LhsExpr>& A,
	       const doublereal b) noexcept {
	  return decltype(operator*(A, b)){A, SpMatElemScalarExpr<doublereal, doublereal>{b}};
     }

     template <typename LhsExpr, typename RhsValue, typename RhsExpr>
     inline constexpr
     SpMatElemBinExpr<typename util::ResultType<SpGradient, RhsValue>::Type,
		      SpGradBinMult,			 
		      const SpMatElemScalarExpr<SpGradient, SpGradient>,
		      const SpMatElemExprBase<RhsValue, RhsExpr>&>
     operator*(const SpGradBase<LhsExpr>& a,
	       const SpMatElemExprBase<RhsValue, RhsExpr>& B) noexcept {
	  return decltype(operator*(a, B)){SpMatElemScalarExpr<SpGradient, SpGradient>{SpGradient{a}}, B}; // Avoid multiple evaluations of a!
     }

     template <typename RhsValue, typename RhsExpr>
     inline constexpr
     SpMatElemBinExpr<typename util::ResultType<doublereal, RhsValue>::Type,
		      SpGradBinMult,			 
		      const SpMatElemScalarExpr<doublereal, doublereal>,
		      const SpMatElemExprBase<RhsValue, RhsExpr>&>
     operator*(doublereal a,
	       const SpMatElemExprBase<RhsValue, RhsExpr>& B) noexcept {
	  return decltype(operator*(a, B)){SpMatElemScalarExpr<doublereal, doublereal>{a}, B};
     }
     
     template <typename LhsValue, typename LhsExpr, typename RhsExpr>
     inline constexpr
     SpMatElemBinExpr<typename util::ResultType<LhsValue, SpGradient>::Type,
		      SpGradBinDiv,
		      const SpMatElemExprBase<LhsValue, LhsExpr>&,
		      const SpMatElemScalarExpr<SpGradient, SpGradient> >
     operator/(const SpMatElemExprBase<LhsValue, LhsExpr>& A,
	       const SpGradBase<RhsExpr>& b) noexcept {
	  return decltype(operator/(A, b)){A, SpMatElemScalarExpr<SpGradient, SpGradient>{SpGradient{b}}}; // Avoid multiple evaluations of b!
     }

     template <typename LhsValue, typename LhsExpr>
     inline constexpr
     SpMatElemBinExpr<typename util::ResultType<LhsValue, doublereal>::Type,
		      SpGradBinDiv,
		      const SpMatElemExprBase<LhsValue, LhsExpr>&,
		      const SpMatElemScalarExpr<doublereal, doublereal> >
     operator/(const SpMatElemExprBase<LhsValue, LhsExpr>& A,
	       const doublereal b) noexcept {
	  return decltype(operator/(A, b)){A, SpMatElemScalarExpr<doublereal, doublereal>{b}};
     }

     template <typename LhsValue, typename RhsValue, typename LhsExpr, typename RhsExpr>
     inline constexpr
     SpMatMulExpr<LhsValue,
		  RhsValue,
		  const SpMatElemExprBase<LhsValue, LhsExpr>&,
		  const SpMatElemExprBase<RhsValue, RhsExpr>&>
     operator*(const SpMatElemExprBase<LhsValue, LhsExpr>& A,
	       const SpMatElemExprBase<RhsValue, RhsExpr>& B) noexcept {
	  return decltype(operator*(A, B)){A, B};
     }
     
     template <typename LhsValue, typename RhsValue, typename LhsExpr, typename RhsExpr>
     inline constexpr
     typename util::ResultType<LhsValue, RhsValue>::Type
     Dot(const SpMatElemExprBase<LhsValue, LhsExpr>& u, const SpMatElemExprBase<RhsValue, RhsExpr>& v) {
	  return *SpMatrixBase<typename util::ResultType<LhsValue, RhsValue>::Type, 1, 1>{Transpose(u) * v}.begin();
     }


     template <typename ValueType>
     inline SpMatrix<ValueType, 3, 3>
     MatRVec(const SpColVector<ValueType, 3>& g) {
	  SpMatrix<ValueType, 3, 3> RDelta(3, 3, g.iGetMaxSize());
	  
	  const ValueType d = 4. / (4. + Dot(g, g));
	  const ValueType tmp1 = -g(3) * g(3);
	  const ValueType tmp2 = -g(2) * g(2);
	  const ValueType tmp3 = -g(1) * g(1);
	  const ValueType tmp4 = g(1) * g(2) * 0.5;
	  const ValueType tmp5 = g(2) * g(3) * 0.5;
	  const ValueType tmp6 = g(1) * g(3) * 0.5;

	  RDelta(1,1) = (tmp1 + tmp2) * d * 0.5 + 1;
	  RDelta(1,2) = (tmp4 - g(3)) * d;
	  RDelta(1,3) = (tmp6 + g(2)) * d;
	  RDelta(2,1) = (g(3) + tmp4) * d;
	  RDelta(2,2) = (tmp1 + tmp3) * d * 0.5 + 1.;
	  RDelta(2,3) = (tmp5 - g(1)) * d;
	  RDelta(3,1) = (tmp6 - g(2)) * d;
	  RDelta(3,2) = (tmp5 + g(1)) * d;
	  RDelta(3,3) = (tmp2 + tmp3) * d * 0.5 + 1.;

	  return RDelta;
     }

     template <typename ValueType>
     inline SpMatrix<ValueType, 3, 3> MatRotVec(const SpColVector<ValueType, 3>& p) {
	  constexpr index_type cid = RotCoeff::COEFF_B;
	  using std::sqrt;
	  using std::sin;
	  using std::cos;
            
	  ValueType phip[10];
	  ValueType phi2(Dot(p, p));
	  ValueType pd(sqrt(phi2));
	  ValueType cf[RotCoeff::COEFF_B];
	  index_type k, j;

	  if (pd < RotCoeff::SerThrsh[cid-1]) {
	       SpGradient::ResizeReset(phip[0], 1., SpGradient::iGetSize(phi2));
	       for (j = 1; j <= 9; j++) {
                    phip[j] = phip[j-1]*phi2;
	       }
	       for (k = 0; k < cid; k++) {
		    SpGradient::ResizeReset(cf[k], 0., 0);
                    for (j = 0; j < RotCoeff::SerTrunc[k]; j++) {
			 cf[k] += phip[j]/RotCoeff::SerCoeff[k][j];
                    }
	       }
	  } else {
	       cf[0] = sin(pd) / pd;                 // a = sin(phi)/phi
	       cf[1]=(1. - cos(pd)) / phi2;           // b = (1.-cos(phi))/phi2
	  }

	  SpMatrix<ValueType, 3, 3> R(3, 3, p.iGetMaxSize());
	  
	  R(1,1) = cf[1]*((-p(3)*p(3))-p(2)*p(2))+1.;
	  R(1,2) = cf[1]*p(1)*p(2)-cf[0]*p(3);
	  R(1,3) = cf[1]*p(1)*p(3)+cf[0]*p(2);
	  R(2,1) = cf[0]*p(3)+cf[1]*p(1)*p(2);
	  R(2,2) = cf[1]*((-p(3)*p(3))-p(1)*p(1))+1.;
	  R(2,3) = cf[1]*p(2)*p(3)-cf[0]*p(1);
	  R(3,1) = cf[1]*p(1)*p(3)-cf[0]*p(2);
	  R(3,2) = cf[1]*p(2)*p(3)+cf[0]*p(1);
	  R(3,3) = cf[1]*((-p(2)*p(2))-p(1)*p(1))+1.;

	  return R;
     }

     template <typename ValueType>
     inline void RotCo(const ValueType& phi, ValueType& cf) {
	  // This algorithm is a simplified version of RotCo in RotCoeff.hc
	  // from Marco Morandini  <morandini@aero.polimi.it>
	  // and Teodoro Merlini  <merlini@aero.polimi.it>
	  using std::sin;
	  using std::cos;
	  using std::sqrt;
	  using std::fabs;
	
	  ValueType phip[10];
	  ValueType phi2(phi * phi);

	  if (fabs(phi) < RotCoeff::SerThrsh[0]) {
	       SpGradient::ResizeReset(phip[0], 1., SpGradient::iGetSize(phi));
            
	       for (index_type j = 1; j <= 9; j++) {
                    phip[j] = phip[j - 1] * phi2;
	       }

	       SpGradient::ResizeReset(cf, 0., SpGradient::iGetSize(phi));
            
	       for (index_type j = 0; j < RotCoeff::SerTrunc[0]; j++) {
                    cf += phip[j] / RotCoeff::SerCoeff[0][j];
	       }

	       return;
	  } 
	
	  const ValueType pd(sqrt(phi2));
	  cf = sin(pd) / pd;                 // a = sin(phi)/phi
     }
     
     template <typename ValueType>
     inline SpColVector<ValueType, 3>
     VecRot(const SpMatrix<ValueType, 3, 3>& R) {
	  // Modified from Appendix 2.4 of
	  //
	  // author = {Marco Borri and Lorenzo Trainelli and Carlo L. Bottasso},
	  // title = {On Representations and Parameterizations of Motion},
	  // journal = {Multibody System Dynamics},
	  // volume = {4},
	  // pages = {129--193},
	  // year = {2000}

	  SpColVector<ValueType, 3> unit(3, R.iGetMaxSize());
	  
	  using std::atan2;
	  using std::sqrt;
            
	  const ValueType cosphi = 0.5 * (R(1, 1) + R(2, 2) + R(3, 3) - 1.);
            
	  if (cosphi > 0.) {
	       unit(1) = 0.5*(R(3, 2) - R(2, 3));
	       unit(2) = 0.5*(R(1, 3) - R(3, 1));
	       unit(3) = 0.5*(R(2, 1) - R(1, 2));
                
	       const ValueType sinphi2 = Dot(unit, unit);
	       ValueType sinphi;
                
	       if (sinphi2 != 0) {
                    sinphi = sqrt(sinphi2);
	       } else {
                    sinphi = unit(1);
	       }
                
	       const ValueType phi = atan2(sinphi, cosphi);
	       ValueType a;
	       RotCo(phi, a);
	       unit /= a;
	  } else {
	       // -1 <= cosphi <= 0
	       SpMatrix<ValueType, 3, 3> eet = (R + Transpose(R)) * 0.5;
	       eet(1, 1) -= cosphi;
	       eet(2, 2) -= cosphi;
	       eet(3, 3) -= cosphi;
	       // largest (abs) component of unit vector phi/|phi|
	       index_type maxcol = 1;
	       if (eet(2, 2) > eet(1, 1)) {
                    maxcol = 2;
	       }
	       if (eet(3, 3) > eet(maxcol, maxcol)) {
                    maxcol = 3;
	       }
	       unit = (eet.GetCol(maxcol)/sqrt(eet(maxcol, maxcol)*(1. - cosphi)));
	       
	       ValueType sinphi{};

	       for (index_type i = 1; i <= 3; ++i) {
		    SpColVector<ValueType, 3> x = Cross(unit, R.GetCol(i));
                    sinphi -= x(i) * 0.5;
	       }

	       unit *= atan2(sinphi, cosphi);
	  }

	  return unit;
     }
     
     template <typename ValueType, index_type NumRows, index_type NumCols>
     inline std::ostream& operator<<(std::ostream& os, const SpMatrixBase<ValueType, NumRows, NumCols>& A) {
	  for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
	       for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
		    os << A.dGetValue(i, j) << ' ';
	       }
	  }

	  return os;
     }
}

#endif
