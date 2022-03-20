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
#include "sp_gradient.h"

namespace sp_grad {
     template <typename ValueType, typename ScalarExpr, index_type NumRows = 1, index_type NumCols = 1>
     class SpMatElemScalarExpr;

     template <typename ValueType, typename Expr>
     class SpMatElemUniqueExpr;

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
               static constexpr SpGradCommon::ExprEvalFlags eExprEvalFlags = eCompr;
          };

          template <SpGradCommon::ExprEvalFlags eCompr>
          struct ComprEvalHelper<doublereal, eCompr> {
               // It makes no sense to allocate a class SpGradExprDofMap for a doublereal
               static const SpGradCommon::ExprEvalFlags eExprEvalFlags = SpGradCommon::ExprEvalDuplicate;
          };

          template <typename Expr, bool bUseTempExpr>
          struct TempExprHelper;

          template <typename Expr>
          struct TempExprHelper<Expr, true> {
               typedef typename remove_all<Expr>::type ExprType;
               typedef typename ExprType::ValueType ValueType;

               template <typename Value, typename ExprType>
               constexpr static SpMatElemUniqueExpr<Value, const SpMatElemExprBase<Value, ExprType>&>
               EvalUnique(const SpMatElemExprBase<Value, ExprType>& oExpr) {
                    static_assert(std::is_same<decltype(oExpr), Expr>::value);
                    return decltype(EvalUnique(oExpr)){oExpr};
               }

               // Don't use const here in order to enable move constructors
               typedef SpMatrixBase<ValueType, ExprType::iNumRowsStatic, ExprType::iNumColsStatic> Type;
          };

          template <typename Expr>
          struct TempExprHelper<Expr, false> {
               template <typename Value, typename ExprType>
               constexpr static const SpMatElemExprBase<Value, ExprType>&
               EvalUnique(const SpMatElemExprBase<Value, ExprType>& oExpr) {
                    static_assert(std::is_same<decltype(oExpr), Expr>::value);
                    return oExpr;
               }

               typedef const Expr Type;
          };

          template <index_type iSizeStatic>
          struct MatrixDataSizeHelper {
               static_assert(iSizeStatic > 0);

               static constexpr index_type iGetSizeStatic(index_type) {
                    return iSizeStatic;
               }
          };

          template <>
          struct MatrixDataSizeHelper<SpMatrixSize::DYNAMIC> {
               static_assert(SpMatrixSize::DYNAMIC < 0);

               static index_type iGetSizeStatic(index_type iSize) {
                    SP_GRAD_ASSERT(iSize >= 0);
                    return iSize;
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
          static constexpr index_type iNumRowsStatic = LhsExprType::eMatOpType == SpMatOpType::MATRIX ? LhsExprType::iNumRowsStatic : RhsExprType::iNumRowsStatic;
          static constexpr index_type iNumColsStatic = LhsExprType::eMatOpType == SpMatOpType::MATRIX ? LhsExprType::iNumColsStatic : RhsExprType::iNumColsStatic;
          static constexpr SpMatOpType eMatOpType = SpMatOpType::MATRIX;
          static constexpr SpGradCommon::ExprEvalFlags eExprEvalFlags = util::ExprEvalFlagsHelper<LhsExprType::eExprEvalFlags,
                                                                                              RhsExprType::eExprEvalFlags>::eExprEvalFlags;
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
                    SpGradCommon::ExprEvalFlags eCompr = eExprEvalFlags,
                    typename ValueTypeA, index_type iNumRowsA, index_type iNumColsA>
          void Eval(SpMatrixBase<ValueTypeA, iNumRowsA, iNumColsA>& A) const {
               this->template ElemEval<eTransp, eCompr>(A);
          }

          template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
                    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalDuplicate,
                    typename Func,
                    typename ValueTypeA,
                    index_type NumRowsA, index_type NumColsA>
          inline void AssignEval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
               this->template ElemAssign<eTransp, eCompr, Func>(A);
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

          typedef util::TempExprHelper<LhsExpr, bUseTempExprLhs> LhsTempHelper;
          typedef util::TempExprHelper<RhsExpr, bUseTempExprRhs> RhsTempHelper;
          typedef typename LhsTempHelper::Type LhsTempExpr;
          typedef typename RhsTempHelper::Type RhsTempExpr;
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
          static constexpr SpGradCommon::ExprEvalFlags eExprEvalFlags = util::ExprEvalFlagsHelper<LhsExprType::eExprEvalFlags,
                                                                                              RhsExprType::eExprEvalFlags>::eExprEvalFlags;
          static_assert(LhsExprType::eMatOpType == SpMatOpType::MATRIX && RhsExprType::eMatOpType == SpMatOpType::MATRIX,
                        "Both operands must be vectors!");
          static_assert(LhsExprType::iNumRowsStatic == 3 && RhsExprType::iNumRowsStatic == 3,
                        "Both operands must be 3x1 vectors");
          static_assert(LhsExprType::iNumColsStatic == 1 && RhsExprType::iNumColsStatic == 1,
                        "Both operands must be 3x1 vectors");

          constexpr SpMatCrossExpr(const LhsExprType& u, const RhsExprType& v) noexcept
               :u(LhsTempHelper::EvalUnique(u)), v(RhsTempHelper::EvalUnique(v)) {
          }

          SpMatCrossExpr(SpMatCrossExpr&& oExpr)
               :u(std::move(oExpr.u)), v(std::move(oExpr.v)) {
          }

          template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
                    SpGradCommon::ExprEvalFlags eCompr = eExprEvalFlags,
                    typename ValueTypeA, index_type iNumRowsA, index_type iNumColsA>
          void Eval(SpMatrixBase<ValueTypeA, iNumRowsA, iNumColsA>& A) const {
               this->template ElemEval<eTransp, eCompr>(A);
          }


          template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
                    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalDuplicate,
                    typename Func,
                    typename ValueTypeA,
                    index_type NumRowsA, index_type NumColsA>
          inline void AssignEval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
               this->template ElemAssign<eTransp, eCompr, Func>(A);
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
          static constexpr SpGradCommon::ExprEvalFlags eExprEvalFlags = ExprType::eExprEvalFlags;

          static_assert(ExprType::eMatOpType == SpMatOpType::MATRIX,
                        "Operand must be a matrix! Use SpGradient for scalar expressions!");

          constexpr explicit SpMatElemUnaryExpr(const Expr& u) noexcept
               :u(u) {
          }

          SpMatElemUnaryExpr(SpMatElemUnaryExpr&& oExpr) noexcept
               :u(std::move(oExpr.u)) {
          }

          template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
                    SpGradCommon::ExprEvalFlags eCompr = eExprEvalFlags, typename ValueTypeA, index_type NumRowsA, index_type NumColsA>
          void Eval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
               this->template ElemEval<eTransp, eCompr>(A);
          }

          template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
                    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalDuplicate,
                    typename Func,
                    typename ValueTypeA,
                    index_type NumRowsA, index_type NumColsA>
          void AssignEval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
               this->template ElemAssign<eTransp, eCompr, Func>(A);
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
     class SpMatElemUniqueExpr
          :public SpMatElemExprBase<ValueType, SpMatElemUniqueExpr<ValueType, Expr> > {
          typedef typename util::remove_all<Expr>::type ExprType;
     public:
          static constexpr index_type iNumElemOps = ExprType::iNumElemOps;
          static constexpr unsigned uMatAccess = ExprType::uMatAccess;
          static constexpr index_type iNumRowsStatic = ExprType::iNumRowsStatic;
          static constexpr index_type iNumColsStatic = ExprType::iNumColsStatic;
          static constexpr SpMatOpType eMatOpType = ExprType::eMatOpType;
          static constexpr SpGradCommon::ExprEvalFlags eExprEvalFlags = util::ComprEvalHelper<ValueType, SpGradCommon::ExprEvalUnique>::eExprEvalFlags;

          static_assert(ExprType::eMatOpType == SpMatOpType::MATRIX,
                        "Operand must be a matrix! Use SpGradient for scalar expressions!");

          constexpr explicit SpMatElemUniqueExpr(const Expr& u) noexcept
               :u(u) {
          }

          template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
                    SpGradCommon::ExprEvalFlags eCompr = eExprEvalFlags, typename ValueTypeA, index_type NumRowsA, index_type NumColsA>
          void Eval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
               u.template Eval<eTransp, eExprEvalFlags>(A);
          }

          template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
                    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalDuplicate,
                    typename Func,
                    typename ValueTypeA,
                    index_type NumRowsA, index_type NumColsA>
          inline void AssignEval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
               this->template ElemAssign<eTransp, eCompr, Func>(A);
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
          static constexpr SpGradCommon::ExprEvalFlags eExprEvalFlags = ExprType::eExprEvalFlags;

          static_assert(ExprType::eMatOpType == SpMatOpType::MATRIX, "Operand must be a matrix! A scalar cannot be transposed!");

          constexpr explicit SpMatElemTranspExpr(const Expr& u) noexcept
               :u(u) {
          }

          template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
                    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalDuplicate,
                    typename ValueTypeA, index_type NumRowsA, index_type NumColsA>
          void Eval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
               u.template Eval<util::TranspMatEvalHelper<eTransp>::Value, eCompr>(A);
          }

          template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
                    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalDuplicate,
                    typename Func,
                    typename ValueTypeA,
                    index_type NumRowsA, index_type NumColsA>
          inline void AssignEval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
               this->template ElemAssign<util::TranspMatEvalHelper<eTransp>::Value, eCompr, Func>(A);
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

          const Expr& Transpose() const { return u; }

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
          static constexpr SpGradCommon::ExprEvalFlags eExprEvalFlags = ExprType::eExprEvalFlags;

          static_assert(ExprType::eMatOpType == SpMatOpType::MATRIX, "Operand must be a matrix!");

          constexpr SpMatColVecExpr(const Expr& u, index_type iCol) noexcept
               :u(u), iCol(iCol) {
          }

          template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
                    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalDuplicate,
                    typename ValueTypeA, index_type NumRowsA, index_type NumColsA>
          void Eval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
               this->template ElemEval<eTransp, eCompr>(A);
          }

          template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
                    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalDuplicate,
                    typename Func,
                    typename ValueTypeA,
                    index_type NumRowsA, index_type NumColsA>
          inline void AssignEval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
               this->template ElemAssign<eTransp, eCompr, Func>(A);
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
          static constexpr SpGradCommon::ExprEvalFlags eExprEvalFlags = ExprType::eExprEvalFlags;

          static_assert(ExprType::eMatOpType == SpMatOpType::MATRIX, "Operand must be a matrix!");

          constexpr SpMatRowVecExpr(const Expr& u, index_type iRow) noexcept
               :u(u), iRow(iRow) {
          }

          template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
                    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalDuplicate,
                    typename ValueTypeA, index_type NumRowsA, index_type NumColsA>
          void Eval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
               this->template ElemEval<eTransp, eCompr>(A);
          }

          template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
                    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalDuplicate,
                    typename Func,
                    typename ValueTypeA,
                    index_type NumRowsA, index_type NumColsA>
          inline void AssignEval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
               this->template ElemAssign<eTransp, eCompr, Func>(A);
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

     template <typename ValueType, typename Expr>
     class SpSubMatDynExpr
          :public SpMatElemExprBase<ValueType, SpSubMatDynExpr<ValueType, Expr> > {
          typedef typename util::remove_all<Expr>::type ExprType;
     public:
          static constexpr index_type iNumElemOps = ExprType::iNumElemOps;
          static constexpr unsigned uMatAccess = ExprType::uMatAccess;
          static constexpr index_type iNumRowsStatic = SpMatrixSize::DYNAMIC;
          static constexpr index_type iNumColsStatic = SpMatrixSize::DYNAMIC;
          static constexpr SpMatOpType eMatOpType = SpMatOpType::MATRIX;
          static constexpr SpGradCommon::ExprEvalFlags eExprEvalFlags = ExprType::eExprEvalFlags;

          static_assert(ExprType::eMatOpType == SpMatOpType::MATRIX, "Operand must be a matrix! A scalar cannot be transposed!");

          explicit SpSubMatDynExpr(const Expr& u, index_type iRowStart, index_type iRowStep, index_type iNumRows, index_type iColStart, index_type iColStep, index_type iNumCols) noexcept
               :u(u),
                iRowStart(iRowStart),
                iRowStep(iRowStep),
                iNumRows(iNumRows),
                iColStart(iColStart),
                iColStep(iColStep),
                iNumCols(iNumCols)
          {
               SP_GRAD_ASSERT(iRowStart >= 1);
               SP_GRAD_ASSERT(iColStart >= 1);
               SP_GRAD_ASSERT(iRowStart <= u.iGetNumRows());
               SP_GRAD_ASSERT(iColStart <= u.iGetNumCols());
               SP_GRAD_ASSERT(iRowStep >= 1);
               SP_GRAD_ASSERT(iColStep >= 1);
               SP_GRAD_ASSERT(iNumRows >= 1);
               SP_GRAD_ASSERT(iNumCols >= 1);
               SP_GRAD_ASSERT(iNumRows <= u.iGetNumRows());
               SP_GRAD_ASSERT(iNumCols <= u.iGetNumCols());

#ifdef SP_GRAD_DEBUG
               for (index_type i = 1; i <= iGetNumRows(); ++i) {
                    SP_GRAD_ASSERT(iGetRowIndex(i) >= 1);
                    SP_GRAD_ASSERT(iGetRowIndex(i) <= u.iGetNumRows());
               }

               for (index_type j = 1; j <= iGetNumCols(); ++j) {
                    SP_GRAD_ASSERT(iGetColIndex(j) >= 1);
                    SP_GRAD_ASSERT(iGetColIndex(j) <= u.iGetNumCols());
               }
#endif
          }

          template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
                    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalDuplicate,
                    typename ValueTypeA, index_type NumRowsA, index_type NumColsA>
          void Eval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
               this->template ElemEval<eTransp, eCompr>(A);
          }

          template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
                    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalDuplicate,
                    typename Func,
                    typename ValueTypeA,
                    index_type NumRowsA, index_type NumColsA>
          inline void AssignEval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
               this->template ElemAssign<eTransp, eCompr, Func>(A);
          }

          constexpr doublereal dGetValue(index_type i, index_type j) const {
               return u.dGetValue(iGetRowIndex(i), iGetColIndex(j));
          }

          constexpr index_type iGetSize(index_type i, index_type j) const {
               return u.iGetSize(iGetRowIndex(i), iGetColIndex(j));
          }

          constexpr index_type iGetNumRows() const noexcept {
               return iNumRows;
          }

          constexpr index_type iGetNumCols() const noexcept {
               return iNumCols;
          }

          constexpr index_type iGetMaxSize() const noexcept {
               return u.iGetMaxSize();
          }

          template <typename ValueTypeB>
          void InsertDeriv(ValueTypeB& g, doublereal dCoef, index_type i, index_type j) const {
               u.InsertDeriv(g, dCoef, iGetRowIndex(i), iGetColIndex(j));
          }

          void GetDofStat(SpGradDofStat& s, index_type i, index_type j) const noexcept {
               u.GetDofStat(s, iGetRowIndex(i), iGetColIndex(j));
          }

          void InsertDof(SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
               u.InsertDof(oExpDofMap, iGetRowIndex(i), iGetColIndex(j));
          }

          template <typename ValueTypeB>
          void AddDeriv(ValueTypeB& g, doublereal dCoef, const SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
               u.AddDeriv(g, dCoef, oExpDofMap, iGetRowIndex(i), iGetColIndex(j));
          }

          template <typename ExprTypeB, typename ExprB>
          constexpr bool bHaveRefTo(const SpMatElemExprBase<ExprTypeB, ExprB>& A) const noexcept {
               return u.bHaveRefTo(A);
          }

          const ValueType* begin() const noexcept { return u.begin() + (iRowStart - 1) * u.iGetRowOffset() + (iColStart - 1) * u.iGetColOffset(); }
          const ValueType* end() const noexcept { return begin() + iGetRowOffset() * iNumRows * iGetColOffset() * iNumCols; }
          inline constexpr index_type iGetRowOffset() const noexcept { return u.iGetRowOffset() * iRowStep; }
          inline constexpr index_type iGetColOffset() const noexcept { return u.iGetColOffset() * iColStep; }

     private:
          index_type iGetRowIndex(index_type i) const noexcept {
               SP_GRAD_ASSERT(i >= 1);
               SP_GRAD_ASSERT(i <= iNumRows);

               const index_type iRowIndex = iRowStart + (i - 1) * iRowStep;

               SP_GRAD_ASSERT(iRowIndex >= 1);
               SP_GRAD_ASSERT(iRowIndex <= u.iGetNumRows());

               return iRowIndex;
          }

          index_type iGetColIndex(index_type j) const noexcept {
               SP_GRAD_ASSERT(j >= 1);
               SP_GRAD_ASSERT(j <= iNumCols);

               const index_type iColIndex = iColStart + (j - 1) * iColStep;

               SP_GRAD_ASSERT(iColIndex >= 1);
               SP_GRAD_ASSERT(iColIndex <= u.iGetNumCols());

               return iColIndex;
          }
          const Expr u;
          const index_type iRowStart, iRowStep, iNumRows, iColStart, iColStep, iNumCols;
     };

     template <typename ValueType, typename Expr, index_type iRowStart, index_type iRowStep, index_type iNumRows, index_type iColStart, index_type iColStep, index_type iNumCols>
     class SpSubMatStatExpr
          :public SpMatElemExprBase<ValueType, SpSubMatStatExpr<ValueType, Expr, iRowStart, iRowStep, iNumRows, iColStart, iColStep, iNumCols> > {
          typedef typename util::remove_all<Expr>::type ExprType;
     public:
          static constexpr index_type iNumElemOps = ExprType::iNumElemOps;
          static constexpr unsigned uMatAccess = ExprType::uMatAccess;
          static constexpr index_type iNumRowsStatic = iNumRows;
          static constexpr index_type iNumColsStatic = iNumCols;
          static constexpr SpMatOpType eMatOpType = SpMatOpType::MATRIX;
          static constexpr SpGradCommon::ExprEvalFlags eExprEvalFlags = ExprType::eExprEvalFlags;

          static_assert(ExprType::iNumRowsStatic != SpMatrixSize::DYNAMIC, "Operand row size must dynamic");
          static_assert(ExprType::iNumColsStatic != SpMatrixSize::DYNAMIC, "Operand column size must be dynamic");
          static_assert(ExprType::eMatOpType == SpMatOpType::MATRIX, "Operand must be a matrix! A scalar cannot be transposed!");
          static_assert(iRowStart >= 1);
          static_assert(iColStart >= 1);
          static_assert(iRowStart <= ExprType::iNumRowsStatic);
          static_assert(iColStart <= ExprType::iNumColsStatic);
          static_assert(iRowStep >= 1);
          static_assert(iColStep >= 1);
          static_assert(iNumRows >= 1);
          static_assert(iNumCols >= 1);
          static_assert(iRowStart + (iNumRows - 1) * iRowStep <= ExprType::iNumRowsStatic);
          static_assert(iColStart + (iNumCols - 1) * iColStep <= ExprType::iNumColsStatic);

          constexpr explicit SpSubMatStatExpr(const Expr& u) noexcept
               :u(u)
               {
#ifdef SP_GRAD_DEBUG
                    for (index_type i = 1; i <= iGetNumRows(); ++i) {
                         SP_GRAD_ASSERT(iGetRowIndex(i) >= 1);
                         SP_GRAD_ASSERT(iGetRowIndex(i) <= u.iGetNumRows());
                    }

                    for (index_type j = 1; j <= iGetNumCols(); ++j) {
                         SP_GRAD_ASSERT(iGetColIndex(j) >= 1);
                         SP_GRAD_ASSERT(iGetColIndex(j) <= u.iGetNumCols());
                    }
#endif
               }

          template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
                    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalDuplicate,
                    typename ValueTypeA, index_type NumRowsA, index_type NumColsA>
          void Eval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
               this->template ElemEval<eTransp, eCompr>(A);
          }

          template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
                    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalDuplicate,
                    typename Func,
                    typename ValueTypeA,
                    index_type NumRowsA, index_type NumColsA>
          inline void AssignEval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
               this->template ElemAssign<eTransp, eCompr, Func>(A);
          }

          constexpr doublereal dGetValue(index_type i, index_type j) const {
               return u.dGetValue(iGetRowIndex(i), iGetColIndex(j));
          }

          constexpr index_type iGetSize(index_type i, index_type j) const {
               return u.iGetSize(iGetRowIndex(i), iGetColIndex(j));
          }

          static constexpr index_type iGetNumRows() noexcept {
               return iNumRowsStatic;
          }

          static constexpr index_type iGetNumCols() noexcept {
               return iNumColsStatic;
          }

          constexpr index_type iGetMaxSize() const noexcept {
               return u.iGetMaxSize();
          }

          template <typename ValueTypeB>
          void InsertDeriv(ValueTypeB& g, doublereal dCoef, index_type i, index_type j) const {
               u.InsertDeriv(g, dCoef, iGetRowIndex(i), iGetColIndex(j));
          }

          void GetDofStat(SpGradDofStat& s, index_type i, index_type j) const noexcept {
               u.GetDofStat(s, iGetRowIndex(i), iGetColIndex(j));
          }

          void InsertDof(SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
               u.InsertDof(oExpDofMap, iGetRowIndex(i), iGetColIndex(j));
          }

          template <typename ValueTypeB>
          void AddDeriv(ValueTypeB& g, doublereal dCoef, const SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
               u.AddDeriv(g, dCoef, oExpDofMap, iGetRowIndex(i), iGetColIndex(j));
          }

          template <typename ExprTypeB, typename ExprB>
          constexpr bool bHaveRefTo(const SpMatElemExprBase<ExprTypeB, ExprB>& A) const noexcept {
               return u.bHaveRefTo(A);
          }

          const ValueType* begin() const noexcept { return u.begin() + (iRowStart - 1) * u.iGetRowOffset() + (iColStart - 1) * u.iGetColOffset(); }
          const ValueType* end() const noexcept { return begin() + iGetRowOffset() * iNumRows * iGetColOffset() * iNumCols; }
          inline constexpr index_type iGetRowOffset() const noexcept { return u.iGetRowOffset() * iRowStep; }
          inline constexpr index_type iGetColOffset() const noexcept { return u.iGetColOffset() * iColStep; }

     private:
          static index_type iGetRowIndex(index_type i) noexcept {
               SP_GRAD_ASSERT(i >= 1);
               SP_GRAD_ASSERT(i <= iNumRows);

               const index_type iRowIndex = iRowStart + (i - 1) * iRowStep;

               SP_GRAD_ASSERT(iRowIndex >= 1);
               SP_GRAD_ASSERT(iRowIndex <= ExprType::iNumRowsStatic);

               return iRowIndex;
          }

          static index_type iGetColIndex(index_type j) noexcept {
               SP_GRAD_ASSERT(j >= 1);
               SP_GRAD_ASSERT(j <= iNumCols);

               const index_type iColIndex = iColStart + (j - 1) * iColStep;

               SP_GRAD_ASSERT(iColIndex >= 1);
               SP_GRAD_ASSERT(iColIndex <= ExprType::iNumColsStatic);

               return iColIndex;
          }

          const Expr u;
     };

     template <typename ValueType, typename Expr, index_type iNumRows, index_type iNumCols>
     class SpSubMatStatResExpr: public SpMatElemExprBase<ValueType, SpSubMatStatResExpr<ValueType, Expr, iNumRows, iNumCols> > {
          typedef typename util::remove_all<Expr>::type ExprType;
     public:
          static constexpr index_type iNumElemOps = ExprType::iNumElemOps;
          static constexpr unsigned uMatAccess = ExprType::uMatAccess;
          static constexpr index_type iNumRowsStatic = iNumRows;
          static constexpr index_type iNumColsStatic = iNumCols;
          static constexpr SpMatOpType eMatOpType = SpMatOpType::MATRIX;
          static constexpr SpGradCommon::ExprEvalFlags eExprEvalFlags = ExprType::eExprEvalFlags;

          static_assert(ExprType::eMatOpType == SpMatOpType::MATRIX, "Operand must be a matrix! A scalar cannot be transposed!");
          static_assert(iNumRows >= 1);
          static_assert(iNumCols >= 1);
          static_assert(iNumRows != SpMatrixSize::DYNAMIC);
          static_assert(iNumCols != SpMatrixSize::DYNAMIC);

          constexpr explicit SpSubMatStatResExpr(const Expr& u, index_type iRowStart, index_type iRowStep, index_type iColStart, index_type iColStep) noexcept
               :u(u),
                iRowStart(iRowStart),
                iRowStep(iRowStep),
                iColStart(iColStart),
                iColStep(iColStep)
               {
#ifdef SP_GRAD_DEBUG
                    SP_GRAD_ASSERT(iRowStart >= 1);
                    SP_GRAD_ASSERT(iColStart >= 1);
                    SP_GRAD_ASSERT(iRowStart <= u.iGetNumRows());
                    SP_GRAD_ASSERT(iColStart <= u.iGetNumCols());
                    SP_GRAD_ASSERT(iRowStep >= 1);
                    SP_GRAD_ASSERT(iColStep >= 1);
                    SP_GRAD_ASSERT(iNumRows >= 1);
                    SP_GRAD_ASSERT(iNumCols >= 1);
                    SP_GRAD_ASSERT(iNumRows <= u.iGetNumRows());
                    SP_GRAD_ASSERT(iNumCols <= u.iGetNumCols());

                    for (index_type i = 1; i <= iGetNumRows(); ++i) {
                         SP_GRAD_ASSERT(iGetRowIndex(i) >= 1);
                         SP_GRAD_ASSERT(iGetRowIndex(i) <= u.iGetNumRows());
                    }

                    for (index_type j = 1; j <= iGetNumCols(); ++j) {
                         SP_GRAD_ASSERT(iGetColIndex(j) >= 1);
                         SP_GRAD_ASSERT(iGetColIndex(j) <= u.iGetNumCols());
                    }
#endif
               }

          template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
                    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalDuplicate,
                    typename ValueTypeA, index_type NumRowsA, index_type NumColsA>
          void Eval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
               this->template ElemEval<eTransp, eCompr>(A);
          }

          template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
                    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalDuplicate,
                    typename Func,
                    typename ValueTypeA,
                    index_type NumRowsA, index_type NumColsA>
          inline void AssignEval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
               this->template ElemAssign<eTransp, eCompr, Func>(A);
          }

          constexpr doublereal dGetValue(index_type i, index_type j) const {
               return u.dGetValue(iGetRowIndex(i), iGetColIndex(j));
          }

          constexpr index_type iGetSize(index_type i, index_type j) const {
               return u.iGetSize(iGetRowIndex(i), iGetColIndex(j));
          }

          static constexpr index_type iGetNumRows() noexcept {
               return iNumRowsStatic;
          }

          static constexpr index_type iGetNumCols() noexcept {
               return iNumColsStatic;
          }

          constexpr index_type iGetMaxSize() const noexcept {
               return u.iGetMaxSize();
          }

          template <typename ValueTypeB>
          void InsertDeriv(ValueTypeB& g, doublereal dCoef, index_type i, index_type j) const {
               u.InsertDeriv(g, dCoef, iGetRowIndex(i), iGetColIndex(j));
          }

          void GetDofStat(SpGradDofStat& s, index_type i, index_type j) const noexcept {
               u.GetDofStat(s, iGetRowIndex(i), iGetColIndex(j));
          }

          void InsertDof(SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
               u.InsertDof(oExpDofMap, iGetRowIndex(i), iGetColIndex(j));
          }

          template <typename ValueTypeB>
          void AddDeriv(ValueTypeB& g, doublereal dCoef, const SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
               u.AddDeriv(g, dCoef, oExpDofMap, iGetRowIndex(i), iGetColIndex(j));
          }

          template <typename ExprTypeB, typename ExprB>
          constexpr bool bHaveRefTo(const SpMatElemExprBase<ExprTypeB, ExprB>& A) const noexcept {
               return u.bHaveRefTo(A);
          }

          const ValueType* begin() const noexcept { return u.begin() + (iRowStart - 1) * u.iGetRowOffset() + (iColStart - 1) * u.iGetColOffset(); }
          const ValueType* end() const noexcept { return begin() + iGetRowOffset() * iNumRows * iGetColOffset() * iNumCols; }
          inline constexpr index_type iGetRowOffset() const noexcept { return u.iGetRowOffset() * iRowStep; }
          inline constexpr index_type iGetColOffset() const noexcept { return u.iGetColOffset() * iColStep; }

     private:
          index_type iGetRowIndex(index_type i) const noexcept {
               SP_GRAD_ASSERT(i >= 1);
               SP_GRAD_ASSERT(i <= iNumRows);

               const index_type iRowIndex = iRowStart + (i - 1) * iRowStep;

               SP_GRAD_ASSERT(iRowIndex >= 1);
               SP_GRAD_ASSERT(iRowIndex <= u.iGetNumRows());

               return iRowIndex;
          }

          index_type iGetColIndex(index_type j) const noexcept {
               SP_GRAD_ASSERT(j >= 1);
               SP_GRAD_ASSERT(j <= iNumCols);

               const index_type iColIndex = iColStart + (j - 1) * iColStep;

               SP_GRAD_ASSERT(iColIndex >= 1);
               SP_GRAD_ASSERT(iColIndex <= u.iGetNumCols());

               return iColIndex;
          }

          const Expr u;
          const index_type iRowStart, iRowStep, iColStart, iColStep;
     };

     template <typename ValueType, typename Expr, index_type iRowStart, index_type iRowStep, index_type iNumRows>
     class SpSubMatStatRowExpr
     :public SpMatElemExprBase<ValueType, SpSubMatStatRowExpr<ValueType, Expr, iRowStart, iRowStep, iNumRows> > {
          typedef typename util::remove_all<Expr>::type ExprType;
     public:
          static constexpr index_type iNumElemOps = ExprType::iNumElemOps;
          static constexpr unsigned uMatAccess = ExprType::uMatAccess;
          static constexpr index_type iNumRowsStatic = iNumRows;
          static constexpr index_type iNumColsStatic = SpMatrixSize::DYNAMIC;
          static constexpr SpMatOpType eMatOpType = SpMatOpType::MATRIX;
          static constexpr SpGradCommon::ExprEvalFlags eExprEvalFlags = ExprType::eExprEvalFlags;

          static_assert(ExprType::iNumRowsStatic >= 1, "Number of rows of operand must be static");
          static_assert(ExprType::iNumColsStatic == iNumColsStatic, "Number of columns of operand must be dynamic");
          static_assert(iNumRows >= 1, "Number of rows must be static");
          static_assert(iRowStart >= 1);
          static_assert(iRowStep >= 1);
          static_assert(iRowStart + (iNumRows - 1) * iRowStep <= ExprType::iNumRowsStatic);
          static_assert(ExprType::eMatOpType == SpMatOpType::MATRIX, "Operand must be a matrix! A scalar cannot be transposed!");

          explicit SpSubMatStatRowExpr(const Expr& u, index_type iColStart, index_type iColStep, index_type iNumCols) noexcept
               :u(u),
                iColStart(iColStart),
                iColStep(iColStep),
                iNumCols(iNumCols)
          {
               SP_GRAD_ASSERT(iRowStart >= 1);
               SP_GRAD_ASSERT(iColStart >= 1);
               SP_GRAD_ASSERT(iRowStart <= u.iGetNumRows());
               SP_GRAD_ASSERT(iColStart <= u.iGetNumCols());
               SP_GRAD_ASSERT(iRowStep >= 1);
               SP_GRAD_ASSERT(iColStep >= 1);
               SP_GRAD_ASSERT(iNumRows >= 1);
               SP_GRAD_ASSERT(iNumCols >= 1);
               SP_GRAD_ASSERT(iNumRows <= u.iGetNumRows());
               SP_GRAD_ASSERT(iNumCols <= u.iGetNumCols());

#ifdef SP_GRAD_DEBUG
               for (index_type i = 1; i <= iGetNumRows(); ++i) {
                    SP_GRAD_ASSERT(iGetRowIndex(i) >= 1);
                    SP_GRAD_ASSERT(iGetRowIndex(i) <= u.iGetNumRows());
               }

               for (index_type j = 1; j <= iGetNumCols(); ++j) {
                    SP_GRAD_ASSERT(iGetColIndex(j) >= 1);
                    SP_GRAD_ASSERT(iGetColIndex(j) <= u.iGetNumCols());
               }
#endif
          }

          template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
                    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalDuplicate,
                    typename ValueTypeA, index_type NumRowsA, index_type NumColsA>
          void Eval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
               this->template ElemEval<eTransp, eCompr>(A);
          }

          template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
                    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalDuplicate,
                    typename Func,
                    typename ValueTypeA,
                    index_type NumRowsA, index_type NumColsA>
          inline void AssignEval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
               this->template ElemAssign<eTransp, eCompr, Func>(A);
          }

          constexpr doublereal dGetValue(index_type i, index_type j) const {
               return u.dGetValue(iGetRowIndex(i), iGetColIndex(j));
          }

          constexpr index_type iGetSize(index_type i, index_type j) const {
               return u.iGetSize(iGetRowIndex(i), iGetColIndex(j));
          }

          static constexpr index_type iGetNumRows() noexcept {
               return iNumRows;
          }

          constexpr index_type iGetNumCols() const noexcept {
               return iNumCols;
          }

          constexpr index_type iGetMaxSize() const noexcept {
               return u.iGetMaxSize();
          }

          template <typename ValueTypeB>
          void InsertDeriv(ValueTypeB& g, doublereal dCoef, index_type i, index_type j) const {
               u.InsertDeriv(g, dCoef, iGetRowIndex(i), iGetColIndex(j));
          }

          void GetDofStat(SpGradDofStat& s, index_type i, index_type j) const noexcept {
               u.GetDofStat(s, iGetRowIndex(i), iGetColIndex(j));
          }

          void InsertDof(SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
               u.InsertDof(oExpDofMap, iGetRowIndex(i), iGetColIndex(j));
          }

          template <typename ValueTypeB>
          void AddDeriv(ValueTypeB& g, doublereal dCoef, const SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
               u.AddDeriv(g, dCoef, oExpDofMap, iGetRowIndex(i), iGetColIndex(j));
          }

          template <typename ExprTypeB, typename ExprB>
          constexpr bool bHaveRefTo(const SpMatElemExprBase<ExprTypeB, ExprB>& A) const noexcept {
               return u.bHaveRefTo(A);
          }

          const ValueType* begin() const noexcept { return u.begin() + (iRowStart - 1) * u.iGetRowOffset() + (iColStart - 1) * u.iGetColOffset(); }
          const ValueType* end() const noexcept { return begin() + iGetRowOffset() * iNumRows * iGetColOffset() * iNumCols; }
          inline constexpr index_type iGetRowOffset() const noexcept { return u.iGetRowOffset() * iRowStep; }
          inline constexpr index_type iGetColOffset() const noexcept { return u.iGetColOffset() * iColStep; }

     private:
          static index_type iGetRowIndex(index_type i) noexcept {
               SP_GRAD_ASSERT(i >= 1);
               SP_GRAD_ASSERT(i <= iNumRows);

               const index_type iRowIndex = iRowStart + (i - 1) * iRowStep;

               return iRowIndex;
          }

          index_type iGetColIndex(index_type j) const noexcept {
               SP_GRAD_ASSERT(j >= 1);
               SP_GRAD_ASSERT(j <= iNumCols);

               const index_type iColIndex = iColStart + (j - 1) * iColStep;

               SP_GRAD_ASSERT(iColIndex >= 1);
               SP_GRAD_ASSERT(iColIndex <= u.iGetNumCols());

               return iColIndex;
          }

          const Expr u;
          const index_type iColStart, iColStep, iNumCols;
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
          static constexpr SpGradCommon::ExprEvalFlags eExprEvalFlags = util::ExprEvalFlagsHelper<LhsExprType::eExprEvalFlags,
                                                                                              RhsExprType::eExprEvalFlags>::eExprEvalFlags;

          static_assert(LhsExprType::eMatOpType == SpMatOpType::MATRIX, "Left hand side of matrix product must be a matrix");
          static_assert(RhsExprType::eMatOpType == SpMatOpType::MATRIX, "Right hand side of matrix product must be a matrix");
          static_assert(LhsExprType::iNumColsStatic == RhsExprType::iNumRowsStatic, "Incompatible matrix sizes in matrix product");

          SpMatMulExpr(const LhsExpr& u, const RhsExpr& v) noexcept
               :u(u), v(v) {
               SP_GRAD_ASSERT(u.iGetNumCols() == v.iGetNumRows());
          }

          template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
                    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalDuplicate,
                    typename ValueTypeA, index_type NumRowsA, index_type NumColsA>
          inline void Eval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const;

          template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
                    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalDuplicate,
                    typename Func,
                    typename ValueTypeA,
                    index_type NumRowsA, index_type NumColsA>
          void AssignEval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const=delete;

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
               static constexpr SpGradCommon::ExprEvalFlags eExprEvalFlags = Expr::eExprEvalFlags;
          };

          template <>
          struct ScalarExprEvalFlagsHelper<doublereal> {
               static constexpr SpGradCommon::ExprEvalFlags eExprEvalFlags = SpGradCommon::ExprEvalDuplicate;
          };

     }

     template <typename ValueType, typename ScalarExpr, index_type NumRows, index_type NumCols>
     class SpMatElemScalarExpr: public SpMatElemExprBase<ValueType, SpMatElemScalarExpr<ValueType, ScalarExpr, NumRows, NumCols> > {
          enum OpCntEstimate: index_type {
               EST_NUMBER_OF_SCALAR_OPS = 1
          };
     public:
          static constexpr index_type iNumElemOps = EST_NUMBER_OF_SCALAR_OPS;
          static constexpr unsigned uMatAccess = util::MatAccessFlag::ELEMENT_WISE;
          static constexpr index_type iNumRowsStatic = NumRows;
          static constexpr index_type iNumColsStatic = NumCols;
          static constexpr SpMatOpType eMatOpType = SpMatOpType::SCALAR;
          static constexpr SpGradCommon::ExprEvalFlags eExprEvalFlags = util::ScalarExprEvalFlagsHelper<typename util::remove_all<ScalarExpr>::type>::eExprEvalFlags;

          inline constexpr explicit SpMatElemScalarExpr(const ScalarExpr& u) noexcept;

          template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
                    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalDuplicate,
                    typename ValueTypeA, index_type NumRowsA, index_type NumColsA>
          void Eval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
               this->template ElemEval<eTransp, eCompr>(A);
          }

          template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
                    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalDuplicate,
                    typename Func,
                    typename ValueTypeA,
                    index_type NumRowsA, index_type NumColsA>
          inline void AssignEval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
               this->template ElemAssign<eTransp, eCompr, Func>(A);
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
          static constexpr SpGradCommon::ExprEvalFlags eExprEvalFlags = SpGradCommon::ExprEvalDuplicate;

          inline SpMatrixBase();
          inline SpMatrixBase(const SpMatrixBase& oMat);
          inline SpMatrixBase(SpMatrixBase&& oMat);
          inline SpMatrixBase(index_type iNumRows, index_type iNumCols, index_type iNumDeriv);
          template <typename ValueTypeExpr, typename Expr>
          inline SpMatrixBase(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr);
          inline SpMatrixBase(const std::initializer_list<ValueType>& rgValues);
          inline ~SpMatrixBase();
          inline void ResizeReset(index_type iNumRows, index_type iNumCols, index_type iNumDeriv);
          inline SpMatrixBase& operator=(SpMatrixBase&& oMat);
          template <typename ValueTypeExpr, typename Expr>
          inline SpMatrixBase& operator=(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr);
          inline SpMatrixBase& operator=(const SpMatrixBase& oMat);
          template <typename Expr>
          inline SpMatrixBase& operator*=(const SpGradBase<Expr>& b);
          inline SpMatrixBase& operator*=(const SpGradient& b);
          inline SpMatrixBase& operator*=(doublereal b);
          template <typename Expr>
          inline SpMatrixBase& operator/=(const SpGradBase<Expr>& b);
          inline SpMatrixBase& operator/=(const SpGradient& b);
          inline SpMatrixBase& operator/=(const doublereal b);
          template <typename ValueTypeExpr, typename Expr>
          inline SpMatrixBase& operator+=(const SpMatElemExprBase<ValueTypeExpr, Expr>& b);
          template <typename ValueTypeExpr, typename Expr>
          inline SpMatrixBase& operator-=(const SpMatElemExprBase<ValueTypeExpr, Expr>& b);
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
          inline constexpr index_type iGetColOffset() const noexcept { return util::MatrixDataSizeHelper<NumRows>::iGetSizeStatic(iGetNumRows()); }
          inline bool bHaveSameRep(const SpMatrixBase& A) const noexcept;
          template <typename ValueTypeB>
          constexpr static inline bool bIsOwnerOf(const SpMatrixData<ValueTypeB>*) noexcept { return false; }
          inline bool bIsOwnerOf(const SpMatrixData<ValueType>* pDataB) const noexcept;
          template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
                    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalDuplicate,
                    typename ValueTypeA, index_type NumRowsA, index_type NumColsA>
          inline void Eval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const;
          template <util::MatTranspEvalFlag eTransp = util::MatTranspEvalFlag::DIRECT,
                    SpGradCommon::ExprEvalFlags eCompr = SpGradCommon::ExprEvalDuplicate,
                    typename Func,
                    typename ValueTypeA,
                    index_type NumRowsA, index_type NumColsA>
          void AssignEval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const;
          inline index_type iGetMaxSize() const;
          inline constexpr SpMatColVecExpr<ValueType, const SpMatrixBase&> GetCol(index_type iCol) const {
               return decltype(GetCol(iCol))(*this, iCol);
          }
          inline constexpr SpMatRowVecExpr<ValueType, const SpMatrixBase&> GetRow(index_type iRow) const {
               return decltype(GetRow(iRow))(*this, iRow);
          }
          SpMatrixBase<doublereal, NumRows, NumCols> GetValue() const {
               SpMatrixBase<doublereal, NumRows, NumCols> A(iGetNumRows(), iGetNumCols(), 0);

               for (index_type i = 1; i <= pData->iGetNumElem(); ++i) {
                    A.GetElem(i) = SpGradient::dGetValue(GetElem(i));
               }

               return A;
          }

          void MakeUnique();

#ifdef SP_GRAD_DEBUG
          bool bValid() const;
          static bool bValid(const SpGradient& g);
          static bool bValid(doublereal d);
#endif
     private:
          typedef SpMatrixDataCTD<ValueType, iNumRowsStatic, iNumColsStatic> SpMatrixDataType;

          explicit inline SpMatrixBase(SpMatrixDataType* pData);

          static inline SpMatrixDataType* pGetNullData();

          static constexpr bool bStaticSize = iNumRowsStatic != SpMatrixSize::DYNAMIC && iNumColsStatic != SpMatrixSize::DYNAMIC;
          SpMatrixDataType* pData;
     };

     template <typename ValueType, index_type NumRows, index_type NumCols>
     inline std::ostream& operator<<(std::ostream& os, const SpMatrixBase<ValueType, NumRows, NumCols>& A);

     template <typename ValueType,
               index_type NumRows = SpMatrixSize::DYNAMIC,
               index_type NumCols = SpMatrixSize::DYNAMIC>
     class SpMatrix: public SpMatrixBase<ValueType, NumRows, NumCols> {
     public:
          inline SpMatrix() {
               static_assert(NumRows == SpMatrixSize::DYNAMIC || NumCols == SpMatrixSize::DYNAMIC,
                             "Do not forget to allocate memory for this matrix!\n"
                             "It will not be allocated by default.\n"
                             "If you need a preallocated matrix, use SpMatrixA instead!");
          }
          inline SpMatrix(index_type iNumRows, index_type iNumCols, index_type iNumDeriv);
          inline SpMatrix(const SpMatrix& oMat)=default;
          inline SpMatrix(SpMatrix&& oMat)=default;
          inline SpMatrix(SpMatrixBase<ValueType, NumRows, NumCols>&& oMat)
               :SpMatrixBase<ValueType, NumRows, NumCols>(std::move(oMat)) {
          }
          inline SpMatrix(const std::initializer_list<ValueType>& rgValues);
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

     template <typename ValueType,
               index_type NumRows = SpMatrixSize::DYNAMIC>
     class SpColVector: public SpMatrixBase<ValueType, NumRows, 1> {
     public:
          inline SpColVector() {
               static_assert(NumRows == SpMatrixSize::DYNAMIC,
                             "Do not forget to allocate memory for this vector!\n"
                             "It will not be allocated by default.\n"
                             "If you need a preallocated vector, use SpColVectorA instead!");
          }
          inline SpColVector(index_type iNumRows, index_type iNumDeriv);
          inline SpColVector(const SpColVector& oVec)=default;
          inline SpColVector(SpColVector&& oVec)=default;
          inline SpColVector(SpMatrixBase<ValueType, NumRows, 1>&& oMat)
               :SpMatrixBase<ValueType, NumRows, 1>(std::move(oMat)) {
          }
          template <typename ValueTypeExpr, typename Expr>
          inline SpColVector(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr);
          inline SpColVector(const std::initializer_list<ValueType>& rgValues);
          inline void ResizeReset(index_type iNumRows, index_type iNumDeriv);
          inline SpColVector& operator=(const SpColVector& oVec)=default;
          inline SpColVector& operator=(SpColVector&& oVec)=default;
          inline SpColVector& operator=(const Vec3& oVec);
          template <typename ValueTypeExpr, typename Expr>
          inline SpColVector& operator=(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr);
          inline ValueType& operator()(index_type iRow);
          inline const ValueType& operator() (index_type iRow) const;
     };

     template <typename ValueType,
               index_type NumCols = SpMatrixSize::DYNAMIC>
     class SpRowVector: public SpMatrixBase<ValueType, 1, NumCols> {
     public:
          inline SpRowVector() {
               static_assert(NumCols == SpMatrixSize::DYNAMIC,
                             "Do not forget to allocate memory for this vector!\n"
                             "It will not be allocated by default.\n"
                             "If you need a preallocated vector, use SpRowVectorA instead!");
          }
          inline SpRowVector(index_type iNumCols, index_type iNumDeriv);
          inline SpRowVector(const SpRowVector& oVec)=default;
          inline SpRowVector(SpRowVector&& oVec)=default;
          inline SpRowVector(SpMatrixBase<ValueType, 1, NumCols>&& oMat)
               :SpMatrixBase<ValueType, 1, NumCols>(std::move(oMat)) {
          }
          template <typename ValueTypeExpr, typename Expr>
          inline SpRowVector(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr);
          inline SpRowVector(const std::initializer_list<ValueType>& rgValues);
          inline void ResizeReset(index_type iNumCols, index_type iNumDeriv);
          inline SpRowVector& operator=(const SpRowVector& oVec)=default;
          inline SpRowVector& operator=(SpRowVector&& oVec)=default;
          template <typename ValueTypeExpr, typename Expr>
          inline SpRowVector& operator=(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr);
          inline ValueType& operator()(index_type iCol);
          inline const ValueType& operator() (index_type iCol) const;
     };

     template <typename ValueType,
               index_type NumRows = SpMatrixSize::DYNAMIC,
               index_type NumCols = SpMatrixSize::DYNAMIC,
               index_type NumDer = 0>
     class SpMatrixA: public SpMatrix<ValueType, NumRows, NumCols> {
     public:
          inline SpMatrixA()
               :SpMatrix<ValueType, NumRows, NumCols>(NumRows, NumCols, NumDer) {
               static_assert(NumRows != SpMatrixSize::DYNAMIC, "matrix size cannot be dynamic for SpMatrixA; use SpMatrix instead");
               static_assert(NumCols != SpMatrixSize::DYNAMIC, "matrix size cannot be dynamic for SpMatrixA; use SpMatrix instead");
               static_assert(NumDer != SpMatrixSize::DYNAMIC, "matrix size cannot be dynamic for SpMatrixA; use SpMatrix instead");
          }

          inline explicit SpMatrixA(index_type iNumDeriv)
               :SpMatrix<ValueType, NumRows, NumCols>(NumRows, NumCols, iNumDeriv) {
          }

          inline SpMatrixA(const SpMatrixA& oMat)=default;
          inline SpMatrixA(SpMatrixA&& oMat)=default;
          inline SpMatrixA(const std::initializer_list<ValueType>& rgValues)
               :SpMatrix<ValueType, NumRows, NumCols>(rgValues) {
          }

          template <typename ValueTypeExpr, typename Expr>
          inline SpMatrixA(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr)
               :SpMatrix<ValueType, NumRows, NumCols>(oExpr) {
          }

          inline SpMatrixA(SpMatrixBase<ValueType, NumRows, NumCols>&& oMat)
               :SpMatrix<ValueType, NumRows, NumCols>(std::move(oMat)) {
          }

          inline void ResizeReset(index_type iNumDer) {
               SpMatrix<ValueType, NumRows, NumCols>::ResizeReset(NumRows, NumCols, iNumDer);
          }

          inline SpMatrixA& operator=(SpMatrixA&& oMat)=default;
          inline SpMatrixA& operator=(const SpMatrixA& oMat)=default;
          inline SpMatrixA& operator=(const Mat3x3& oMat) {
               SpMatrix<ValueType, NumRows, NumCols>::operator=(oMat);
               return *this;
          }
          template <typename ValueTypeExpr, typename Expr>
          inline SpMatrixA& operator=(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr) {
               SpMatrix<ValueType, NumRows, NumCols>::operator=(oExpr);
               return *this;
          }
     };

     template <typename ValueType,
               index_type NumRows = SpMatrixSize::DYNAMIC,
               index_type NumDer = 0>
     class SpColVectorA: public SpColVector<ValueType, NumRows> {
     public:
          inline SpColVectorA()
               :SpColVector<ValueType, NumRows>(NumRows, NumDer) {
               static_assert(NumRows != SpMatrixSize::DYNAMIC, "matrix size cannot be dynamic for SpColVectorA; use SpColVector instead");
               static_assert(NumDer != SpMatrixSize::DYNAMIC, "matrix size cannot be dynamic for SpColVectorA; use SpColVector instead");
          }

          inline explicit SpColVectorA(index_type iNumDeriv)
               :SpColVector<ValueType, NumRows>(NumRows, iNumDeriv) {
          }

          inline SpColVectorA(const SpColVectorA&)=default;
          inline SpColVectorA(SpColVectorA&&)=default;
          inline SpColVectorA(const std::initializer_list<ValueType>& rgValues)
               :SpColVector<ValueType, NumRows>(rgValues) {
          }
          template <typename ValueTypeExpr, typename Expr>
          inline SpColVectorA(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr)
               :SpColVector<ValueType, NumRows>(oExpr) {
          }

          inline SpColVectorA(SpMatrixBase<ValueType, NumRows, 1>&& oMat)
               :SpColVector<ValueType, NumRows>(std::move(oMat)) {
          }

          inline void ResizeReset(index_type iNumDeriv) {
               SpColVector<ValueType, NumRows>::ResizeReset(NumRows, iNumDeriv);
          }

          inline SpColVectorA& operator=(SpColVectorA&& oMat)=default;
          inline SpColVectorA& operator=(const SpColVectorA& oMat)=default;
          inline SpColVectorA& operator=(const Mat3x3& oMat) {
               SpColVector<ValueType, NumRows>::operator=(oMat);
               return *this;
          }
          template <typename ValueTypeExpr, typename Expr>
          inline SpColVectorA& operator=(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr) {
               SpColVector<ValueType, NumRows>::operator=(oExpr);
               return *this;
          }
     };

     template <typename ValueType,
               index_type NumCols = SpMatrixSize::DYNAMIC,
               index_type NumDer = 0>
     class SpRowVectorA: public SpRowVector<ValueType, NumCols> {
     public:
          inline SpRowVectorA()
               :SpRowVector<ValueType, NumCols>(NumCols, NumDer) {
               static_assert(NumCols != SpMatrixSize::DYNAMIC, "matrix size cannot be dynamic for SpRowVectorA; use SpRowVector instead");
               static_assert(NumDer != SpMatrixSize::DYNAMIC, "matrix size cannot be dynamic for SpRowVectorA; use SpRowVector instead");
          }

          inline explicit SpRowVectorA(index_type iNumDeriv)
               :SpRowVector<ValueType, NumCols>(NumCols, iNumDeriv) {
          }

          inline SpRowVectorA(const SpRowVectorA&)=default;
          inline SpRowVectorA(SpRowVectorA&&)=default;

          inline SpRowVectorA(const std::initializer_list<ValueType>& rgValues)
               :SpRowVector<ValueType, NumCols>(rgValues) {
          }

          template <typename ValueTypeExpr, typename Expr>
          inline SpRowVectorA(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr)
               :SpRowVector<ValueType, NumCols>(oExpr) {
          }

          inline SpRowVectorA(SpMatrixBase<ValueType, 1, NumCols>&& oMat)
               :SpRowVector<ValueType, NumCols>(std::move(oMat)) {
          }

          inline void ResizeReset(index_type iNumDeriv) {
               SpRowVector<ValueType, NumCols>::ResizeReset(NumCols, iNumDeriv);
          }

          inline SpRowVectorA& operator=(SpRowVectorA&& oMat)=default;
          inline SpRowVectorA& operator=(const SpRowVectorA& oMat)=default;

          inline SpRowVectorA& operator=(const Mat3x3& oMat) {
               SpRowVector<ValueType, NumCols>::operator=(oMat);
               return *this;
          }

          template <typename ValueTypeExpr, typename Expr>
          inline SpRowVectorA& operator=(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr) {
               SpRowVector<ValueType, NumCols>::operator=(oExpr);
               return *this;
          }
     };

     template <typename ValueType>
     index_type SpMatrixBaseData<ValueType>::iGetNumRows() const noexcept {
          SP_GRAD_ASSERT(iNumRows >= 0);
          SP_GRAD_ASSERT(iNumCols >= 0);
          return iNumRows;
     }

     template <typename ValueType>
     index_type SpMatrixBaseData<ValueType>::iGetNumCols() const noexcept {
          SP_GRAD_ASSERT(iNumRows >= 0);
          SP_GRAD_ASSERT(iNumCols >= 0);

          return iNumCols;
     }

     template <typename ValueType>
     index_type SpMatrixBaseData<ValueType>::iGetNumElem() const noexcept {
          SP_GRAD_ASSERT(iNumRows >= 0);
          SP_GRAD_ASSERT(iNumCols >= 0);

          return iNumRows * iNumCols;
     }

     template <typename ValueType>
     index_type SpMatrixBaseData<ValueType>::iGetRefCnt() const noexcept {
          SP_GRAD_ASSERT(iRefCnt >= 0);

          return iRefCnt;
     }

     template <typename ValueType>
     index_type SpMatrixBaseData<ValueType>::iGetMaxDeriv() const noexcept {
          SP_GRAD_ASSERT(iNumDeriv >= 0);

          return iNumDeriv;
     }

     template <typename ValueType>
     bool SpMatrixBaseData<ValueType>::bCheckSize(index_type iNumRowsReq, index_type iNumColsReq, index_type iNumDerivReq) const noexcept {
          return iNumRowsReq == iNumRows && iNumColsReq == iNumCols && iNumDerivReq <= iNumDeriv && iRefCnt <= 1;
     }

     template <typename ValueType>
     SpMatrixBaseData<ValueType>::SpMatrixBaseData(index_type iNumRows, index_type iNumCols, index_type iRefCnt, index_type iNumDeriv) noexcept
          :iNumRows(iNumRows),
           iNumCols(iNumCols),
           iRefCnt(iRefCnt),
           iNumDeriv(iNumDeriv) {
     }

     template <typename ValueType>
     SpMatrixBaseData<ValueType>* SpMatrixBaseData<ValueType>::pGetNullData() noexcept {
          SP_GRAD_ASSERT(oNullData.pNullData != nullptr);
          SP_GRAD_ASSERT(oNullData.pNullData->iRefCnt >= 0);
          SP_GRAD_ASSERT(oNullData.pNullData->iNumRows == 0);
          SP_GRAD_ASSERT(oNullData.pNullData->iNumCols == 0);

          return oNullData.pNullData;
     }

     template <typename ValueType>
     SpMatrixBaseData<ValueType>::NullData::NullData()
          :pNullData{new SpMatrixBaseData<ValueType>{0, 0, 1, 0}} {
     }

     template <typename ValueType>
     SpMatrixBaseData<ValueType>::NullData::~NullData()
     {
          SP_GRAD_ASSERT(pNullData != nullptr);
          SP_GRAD_ASSERT(pNullData->iRefCnt == 1);
          SP_GRAD_ASSERT(pNullData->iNumRows == 0);
          SP_GRAD_ASSERT(pNullData->iNumCols == 0);
          SP_GRAD_ASSERT(pNullData->iNumDeriv == 0);

          delete pNullData;

          pNullData = nullptr;
     }

     template <typename ValueType>
     typename SpMatrixBaseData<ValueType>::NullData SpMatrixBaseData<ValueType>::oNullData;

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

                    new(pDerivData) SpDerivData(0., iNumDeriv, 0, true, 0, &oData);
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
          :SpMatrixBaseData<ValueType>(iNumRows, iNumCols, iRefCnt, iNumDeriv) {

          util::SpMatrixDataTraits<ValueType>::Construct(*this, iNumDeriv, pExtraMem);
     }

     template <typename ValueType>
     SpMatrixData<ValueType>::~SpMatrixData() {
          for (auto& g: *this) {
               g.~ValueType();
          }
     }

     template <typename ValueType>
     ValueType* SpMatrixData<ValueType>::pGetData() noexcept {
          return this->rgData;
     }

     template <typename ValueType>
     ValueType* SpMatrixData<ValueType>::pGetData(index_type iRow) noexcept {
          SP_GRAD_ASSERT(iRow >= 1);
          SP_GRAD_ASSERT(iRow <= this->iNumRows * this->iNumCols);

          return pGetData() + iRow - 1;
     }

     template <typename ValueType>
     const ValueType* SpMatrixData<ValueType>::pGetData() const noexcept {
          return this->rgData;
     }

     template <typename ValueType>
     const ValueType* SpMatrixData<ValueType>::pGetData(index_type iRow) const noexcept {
          SP_GRAD_ASSERT(iRow >= 1);
          SP_GRAD_ASSERT(iRow <= this->iNumRows * this->iNumCols);

          return pGetData() + iRow - 1;
     }

     template <typename ValueType>
     void SpMatrixData<ValueType>::IncRef() noexcept {
          SP_GRAD_ASSERT(this->iRefCnt >= 0);
          ++this->iRefCnt;
     }

     template <typename ValueType>
     void SpMatrixData<ValueType>::DecRef() {
          SP_GRAD_ASSERT(this->iRefCnt >= 1);

          if (--this->iRefCnt == 0) {
               SP_GRAD_ASSERT(this != this->pGetNullData());
               this->~SpMatrixData();

#if defined(HAVE_ALIGNED_MALLOC)
               _aligned_free(this);
#else
               free(this);
#endif
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
          return pGetData() + this->iGetNumElem();
     }

     template <typename ValueType>
     const ValueType* SpMatrixData<ValueType>::begin() const noexcept {
          return pGetData();
     }

     template <typename ValueType>
     const ValueType* SpMatrixData<ValueType>::end() const noexcept {
          return pGetData() + this->iGetNumElem();
     }

     template <typename ValueType>
     constexpr bool SpMatrixData<ValueType>::bIsOwnerOf(const ValueType* pData) const noexcept {
          return pData >= begin() && pData < end();
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
     constexpr index_type SpMatrixDataCTD<ValueType, NumRows, NumCols>::iGetNumElem() const noexcept {
          SP_GRAD_ASSERT(iGetNumRows() >= 0);
          SP_GRAD_ASSERT(iGetNumCols() >= 0);

          constexpr bool bDynamicSize = NumRows == SpMatrixSize::DYNAMIC || NumCols == SpMatrixSize::DYNAMIC;
          constexpr index_type iNumElemStatic = bDynamicSize ? SpMatrixSize::DYNAMIC : NumRows * NumCols;
          return util::MatrixDataSizeHelper<iNumElemStatic>::iGetSizeStatic(iGetNumRows() * iGetNumCols());
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     const ValueType* SpMatrixDataCTD<ValueType, NumRows, NumCols>::pGetData(index_type iRow, index_type iCol) const noexcept {
          SP_GRAD_ASSERT(iRow >= 1);
          SP_GRAD_ASSERT(iRow <= iGetNumRows());
          SP_GRAD_ASSERT(iCol >= 1);
          SP_GRAD_ASSERT(iCol <= iGetNumCols());

          return pGetData() + (iCol - 1) * iGetNumRows() + iRow - 1;
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     ValueType* SpMatrixDataCTD<ValueType, NumRows, NumCols>::pGetData(index_type iRow, index_type iCol) noexcept {
          SP_GRAD_ASSERT(iRow >= 1);
          SP_GRAD_ASSERT(iRow <= iGetNumRows());
          SP_GRAD_ASSERT(iCol >= 1);
          SP_GRAD_ASSERT(iCol <= iGetNumCols());

          return pGetData() + (iCol - 1) * iGetNumRows() + iRow - 1;
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     ValueType* SpMatrixDataCTD<ValueType, NumRows, NumCols>::end() noexcept {
          return pGetData() + iGetNumElem();
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     const ValueType* SpMatrixDataCTD<ValueType, NumRows, NumCols>::end() const noexcept {
          return pGetData() + iGetNumElem();
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     inline constexpr index_type SpMatrixDataCTD<ValueType, NumRows, NumCols>::iGetNumRows() const noexcept {
          return util::MatrixDataSizeHelper<NumRows>::iGetSizeStatic(SpMatrixData<ValueType>::iGetNumRows());
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     inline constexpr index_type SpMatrixDataCTD<ValueType, NumRows, NumCols>::iGetNumCols() const noexcept {
          return util::MatrixDataSizeHelper<NumCols>::iGetSizeStatic(SpMatrixData<ValueType>::iGetNumCols());
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
          char* pMem = reinterpret_cast<char*>(malloc(uSize));
#endif

          if (!pMem) {
               throw std::bad_alloc();
          }

          void* pExtraMem = pMem + uSizeGrad + uOffsetDeriv;

          return new(pMem) SpMatrixDataCTD(iNumRows, iNumCols, 0, iNumDeriv, pExtraMem);
     }

     namespace util {
          template <MatTranspEvalFlag eTransp, typename Expr, typename ValueType, index_type NumRows, index_type NumCols>
          void MatEvalHelperCompr<SpGradCommon::ExprEvalDuplicate>::ElemEval(const Expr& oExpr, SpMatrixBase<ValueType, NumRows, NumCols>& A) {
               static_assert(eExprEvalFlags == SpGradCommon::ExprEvalDuplicate);

               oExpr.template ElemEvalUncompr<eTransp>(A);
          }

          template <MatTranspEvalFlag eTransp, typename Expr, typename ValueType, index_type NumRows, index_type NumCols>
          void MatEvalHelperCompr<SpGradCommon::ExprEvalUnique>::ElemEval(const Expr& oExpr, SpMatrixBase<ValueType, NumRows, NumCols>& A) {
               static_assert(eExprEvalFlags == SpGradCommon::ExprEvalUnique);

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

          template <typename ValueTypeA, typename ValueTypeB>
          struct ElemAssignHelper {
               template <MatTranspEvalFlag eTransp, typename Func, typename ExprB, index_type NumRowsA, index_type NumColsA>
               static inline void ElemAssignCompr(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A, const SpMatElemExprBase<ValueTypeB, ExprB>& B) = delete;

               template <MatTranspEvalFlag eTransp, typename Func, typename ExprB, index_type NumRowsA, index_type NumColsA>
               static inline void ElemAssignUncompr(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A, const SpMatElemExprBase<ValueTypeB, ExprB>& B) = delete;
          };

          template <>
          struct ElemAssignHelper<doublereal, doublereal> {
               template <MatTranspEvalFlag eTransp, typename Func, typename ExprB, index_type NumRowsA, index_type NumColsA>
               static inline void ElemAssignCompr(SpMatrixBase<doublereal, NumRowsA, NumColsA>& A, const SpMatElemExprBase<doublereal, ExprB>& B) {
                    typedef typename remove_all<ExprB>::type ExprTypeB;
                    constexpr index_type iNumRowsStatic = ExprTypeB::iNumRowsStatic;
                    constexpr index_type iNumColsStatic = ExprTypeB::iNumColsStatic;

                    const index_type iNumRows = A.iGetNumRows();
                    const index_type iNumCols = A.iGetNumCols();

                    static_assert(NumRowsA == iNumRowsStatic, "Number of rows does not match");
                    static_assert(NumColsA == iNumColsStatic, "Number of columns does not match");

                    SP_GRAD_ASSERT(B.iGetNumRows() == iNumRows || B.iGetNumRows() == SpMatrixSize::DYNAMIC);
                    SP_GRAD_ASSERT(B.iGetNumCols() == iNumCols || B.iGetNumCols() == SpMatrixSize::DYNAMIC);

                    SpMatrixBase<doublereal, NumRowsA, NumColsA> Atmp;

                    if (!B.bHaveRefTo(A)) {
                         Atmp = std::move(A);
                    } else {
                         Atmp = A;
                    }

                    for (index_type j = 1; j <= iNumCols; ++j) {
                         for (index_type i = 1; i <= iNumRows; ++i) {
                              doublereal& Atmpij = Atmp.GetElem(i, j);
                              Atmpij = Func::f(Atmpij, B.dGetValue(i, j));
                         }
                    }

                    A = std::move(Atmp);
               }

               template <MatTranspEvalFlag eTransp, typename Func, typename ExprB, index_type NumRowsA, index_type NumColsA>
               static inline void ElemAssignUncompr(SpMatrixBase<doublereal, NumRowsA, NumColsA>& A, const SpMatElemExprBase<doublereal, ExprB>& B) {
                    ElemAssignCompr<eTransp, Func>(A, B);
               }
          };


          template <>
          struct ElemAssignHelper<SpGradient, doublereal> {
               template <MatTranspEvalFlag eTransp, typename Func, typename ExprB, index_type NumRowsA, index_type NumColsA>
               static inline void ElemAssignCompr(SpMatrixBase<SpGradient, NumRowsA, NumColsA>& A, const SpMatElemExprBase<doublereal, ExprB>& B) {
                    typedef typename remove_all<ExprB>::type ExprTypeB;
                    constexpr index_type iNumRowsStatic = ExprTypeB::iNumRowsStatic;
                    constexpr index_type iNumColsStatic = ExprTypeB::iNumColsStatic;

                    const index_type iNumRows = A.iGetNumRows();
                    const index_type iNumCols = A.iGetNumCols();

                    static_assert(NumRowsA == iNumRowsStatic, "Number of rows does not match");
                    static_assert(NumColsA == iNumColsStatic, "Number of columns does not match");

                    SP_GRAD_ASSERT(B.iGetNumRows() == iNumRows || B.iGetNumRows() == SpMatrixSize::DYNAMIC);
                    SP_GRAD_ASSERT(B.iGetNumCols() == iNumCols || B.iGetNumCols() == SpMatrixSize::DYNAMIC);

                    for (index_type j = 1; j <= iNumCols; ++j) {
                         for (index_type i = 1; i <= iNumRows; ++i) {
                              SpGradient& Aij = A.GetElem(i, j);
                              const doublereal uij = Aij.dGetValue();
                              const doublereal vij = B.dGetValue(i, j);
                              const doublereal fij = Func::f(uij, vij);
                              const doublereal dfij_du = Func::df_du(uij, vij);
                              Aij.template InitDerivAssign<Func>(fij, dfij_du, Aij.iGetSize());
                         }
                    }
               }

               template <MatTranspEvalFlag eTransp, typename Func, typename ExprB, index_type NumRowsA, index_type NumColsA>
               static inline void ElemAssignUncompr(SpMatrixBase<SpGradient, NumRowsA, NumColsA>& A, const SpMatElemExprBase<doublereal, ExprB>& B) {
                    ElemAssignCompr<eTransp, Func>(A, B);
               }
          };

          template <>
          struct ElemAssignHelper<SpGradient, SpGradient> {
               template <MatTranspEvalFlag eTransp, typename Func, typename ExprB, index_type NumRowsA, index_type NumColsA>
               static inline void ElemAssignCompr(SpMatrixBase<SpGradient, NumRowsA, NumColsA>& A, const SpMatElemExprBase<SpGradient, ExprB>& B) {
                    typedef typename remove_all<ExprB>::type ExprTypeB;
                    constexpr index_type iNumRowsStatic = ExprTypeB::iNumRowsStatic;
                    constexpr index_type iNumColsStatic = ExprTypeB::iNumColsStatic;

                    const index_type iNumRows = A.iGetNumRows();
                    const index_type iNumCols = A.iGetNumCols();

                    static_assert(NumRowsA == iNumRowsStatic, "Number of rows does not match");
                    static_assert(NumColsA == iNumColsStatic, "Number of columns does not match");

                    SP_GRAD_ASSERT(B.iGetNumRows() == iNumRows || B.iGetNumRows() == SpMatrixSize::DYNAMIC);
                    SP_GRAD_ASSERT(B.iGetNumCols() == iNumCols || B.iGetNumCols() == SpMatrixSize::DYNAMIC);

                    SpMatrixBase<SpGradient, NumRowsA, NumColsA> Atmp;

                    if (!B.bHaveRefTo(A)) {
                         Atmp = std::move(A);
                    } else {
                         Atmp = A;
                    }

                    SP_GRAD_ASSERT(Atmp.bValid());

                    SpGradExpDofMap oDofMap;

                    for (index_type j = 1; j <= iNumCols; ++j) {
                         for (index_type i = 1; i <= iNumRows; ++i) {
                              SpGradient& Atmpij = Atmp.GetElem(i, j);

                              if (!Atmpij.bIsUnique()) {
                                   Atmpij.MakeUnique();
                              }

                              SpGradDofStat oDofStat;

                              Atmpij.GetDofStat(oDofStat);
                              B.GetDofStat(oDofStat, i, j);

                              oDofMap.Reset(oDofStat);

                              Atmpij.InsertDof(oDofMap);
                              B.InsertDof(oDofMap, i, j);

                              oDofMap.InsertDone();

                              const doublereal uij = Atmpij.dGetValue();
                              const doublereal vij = B.dGetValue(i, j);
                              const doublereal fij = Func::f(uij, vij);
                              const doublereal dfij_du = Func::df_du(uij, vij);
                              const doublereal dfij_dv = Func::df_dv(uij, vij);

                              Atmpij.template InitDerivAssign<Func>(fij, dfij_du, oDofMap);

                              B.AddDeriv(Atmpij, dfij_dv, oDofMap, i, j);

                              SP_GRAD_ASSERT(Atmpij.bIsUnique());
                         }
                    }

                    A = std::move(Atmp);

                    SP_GRAD_ASSERT(A.bValid());
               }

               template <MatTranspEvalFlag eTransp, typename Func, typename ExprB, index_type NumRowsA, index_type NumColsA>
               static inline void ElemAssignUncompr(SpMatrixBase<SpGradient, NumRowsA, NumColsA>& A, const SpMatElemExprBase<SpGradient, ExprB>& B) {
                    typedef typename remove_all<ExprB>::type ExprTypeB;
                    constexpr index_type iNumRowsStatic = ExprTypeB::iNumRowsStatic;
                    constexpr index_type iNumColsStatic = ExprTypeB::iNumColsStatic;

                    const index_type iNumRows = A.iGetNumRows();
                    const index_type iNumCols = A.iGetNumCols();

                    static_assert(NumRowsA == iNumRowsStatic, "Number of rows does not match");
                    static_assert(NumColsA == iNumColsStatic, "Number of columns does not match");

                    SP_GRAD_ASSERT(B.iGetNumRows() == iNumRows || B.iGetNumRows() == SpMatrixSize::DYNAMIC);
                    SP_GRAD_ASSERT(B.iGetNumCols() == iNumCols || B.iGetNumCols() == SpMatrixSize::DYNAMIC);

                    SpMatrixBase<SpGradient, NumRowsA, NumColsA> Atmp;

                    if (!B.bHaveRefTo(A)) {
                         Atmp = std::move(A);
                    } else {
                         Atmp = A;
                    }

                    SP_GRAD_ASSERT(Atmp.bValid());

                    for (index_type j = 1; j <= iNumCols; ++j) {
                         for (index_type i = 1; i <= iNumRows; ++i) {
                              SpGradient& Atmpij = Atmp.GetElem(i, j);

                              const doublereal uij = Atmpij.dGetValue();
                              const doublereal vij = B.dGetValue(i, j);
                              const doublereal fij = Func::f(uij, vij);
                              const doublereal dfij_du = Func::df_du(uij, vij);
                              const doublereal dfij_dv = Func::df_dv(uij, vij);

                              Atmpij.template InitDerivAssign<Func>(fij, dfij_du, Atmpij.iGetSize() + B.iGetSize(i, j));

                              B.InsertDeriv(Atmpij, dfij_dv, i, j);

                              SP_GRAD_ASSERT(!Atmpij.bIsUnique());
                         }
                    }

                    A = std::move(Atmp);

                    SP_GRAD_ASSERT(A.bValid());
               }
          };

          template <SpGradCommon::ExprEvalFlags eFlags>
          struct ElemAssignHelperCompr {
               template <MatTranspEvalFlag eTransp, typename Func, typename ValueA, typename ValueB, typename ExprB, index_type NumRowsA, index_type NumColsA>
               static inline void ElemAssign(SpMatrixBase<ValueA, NumRowsA, NumColsA>& A, const SpMatElemExprBase<ValueB, ExprB>& B) = delete;
          };

          template <>
          struct ElemAssignHelperCompr<SpGradCommon::ExprEvalUnique> {
               template <MatTranspEvalFlag eTransp, typename Func, typename ValueA, typename ValueB, typename ExprB, index_type NumRowsA, index_type NumColsA>
               static inline void ElemAssign(SpMatrixBase<ValueA, NumRowsA, NumColsA>& A, const SpMatElemExprBase<ValueB, ExprB>& B) {
                    ElemAssignHelper<ValueA, ValueB>::template ElemAssignCompr<eTransp, Func>(A, B);
               }
          };

          template <>
          struct ElemAssignHelperCompr<SpGradCommon::ExprEvalDuplicate> {
               template <MatTranspEvalFlag eTransp, typename Func, typename ValueA, typename ValueB, typename ExprB, index_type NumRowsA, index_type NumColsA>
               static inline void ElemAssign(SpMatrixBase<ValueA, NumRowsA, NumColsA>& A, const SpMatElemExprBase<ValueB, ExprB>& B) {
                    ElemAssignHelper<ValueA, ValueB>::template ElemAssignUncompr<eTransp, Func>(A, B);
               }
          };
     };

     template <typename ValueType, typename DERIVED>
     template <util::MatTranspEvalFlag eTransp,
               SpGradCommon::ExprEvalFlags eCompr,
               typename ValueTypeA,
               index_type NumRowsA, index_type NumColsA>
     void SpMatElemExprBase<ValueType, DERIVED>::ElemEval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
          util::MatEvalHelperCompr<eCompr>::template ElemEval<eTransp>(*this, A);
     }

     template <typename ValueType, typename DERIVED>
     template <util::MatTranspEvalFlag eTransp,
               SpGradCommon::ExprEvalFlags eCompr,
               typename Func,
               typename ValueTypeA,
               index_type NumRowsA, index_type NumColsA>
     void SpMatElemExprBase<ValueType, DERIVED>::ElemAssign(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
          util::ElemAssignHelperCompr<eCompr>::template ElemAssign<eTransp, Func>(A, *this);
     }

     template <typename ValueType, typename DERIVED>
     template <util::MatTranspEvalFlag eTransp, typename ValueTypeA, index_type NumRowsA, index_type NumColsA>
     void SpMatElemExprBase<ValueType, DERIVED>::ElemEvalUncompr(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
          const index_type iNumRows = iGetNumRows();
          const index_type iNumCols = iGetNumCols();
          const index_type iMaxSize = iGetMaxSize();

          typedef util::MatEvalHelperTransp<eTransp> MatEvalType;

          constexpr bool bTransposed = eTransp == util::MatTranspEvalFlag::TRANSPOSED;

          static_assert(NumRowsA == (!bTransposed ? iNumRowsStatic : iNumColsStatic), "Number of rows does not match");
          static_assert(NumColsA == (!bTransposed ? iNumColsStatic : iNumRowsStatic), "Number of columns does not match");

          SpMatrixBase<ValueTypeA, NumRowsA, NumColsA> Atmp;

          if (!bHaveRefTo(A)) {
               Atmp = std::move(A);
          }

          MatEvalType::ResizeReset(Atmp, iNumRows, iNumCols, iMaxSize);

          SP_GRAD_ASSERT(Atmp.iGetNumRows() == (!bTransposed ? iNumRows : iNumCols));
          SP_GRAD_ASSERT(Atmp.iGetNumCols() == (!bTransposed ? iNumCols : iNumRows));

          for (index_type j = 1; j <= iNumCols; ++j) {
               for (index_type i = 1; i <= iNumRows; ++i) {
                    ValueTypeA& Ai = MatEvalType::GetElem(Atmp, i, j);

                    SP_GRAD_ASSERT(SpGradient::iGetSize(Ai) == 0);
                    SP_GRAD_ASSERT(SpGradient::dGetValue(Ai) == 0.);

                    SpGradient::ResizeReset(Ai, dGetValue(i, j), iGetSize(i, j));
                    InsertDeriv(Ai, 1., i, j);
               }
          }

          A = std::move(Atmp);
     }

     template <typename ValueType, typename DERIVED>
     template <util::MatTranspEvalFlag eTransp, typename ValueTypeA, index_type NumRowsA, index_type NumColsA>
     void SpMatElemExprBase<ValueType, DERIVED>::ElemEvalCompr(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
          const index_type iNumRows = iGetNumRows();
          const index_type iNumCols = iGetNumCols();

          typedef util::MatEvalHelperTransp<eTransp> MatEvalType;
          constexpr bool bTransposed = eTransp == util::MatTranspEvalFlag::TRANSPOSED;

          static_assert(NumRowsA == (!bTransposed ? iNumRowsStatic : iNumColsStatic), "Number of rows does not match");
          static_assert(NumColsA == (!bTransposed ? iNumColsStatic : iNumRowsStatic), "Number of columns does not match");

          SpMatrixBase<ValueTypeA, NumRowsA, NumColsA> Atmp;

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
                    ValueTypeA& Ai = MatEvalType::GetElem(Atmp, i, j);

                    SP_GRAD_ASSERT(SpGradient::iGetSize(Ai) == 0);
                    SP_GRAD_ASSERT(SpGradient::dGetValue(Ai) == 0.);

                    Ai.ResizeReset(dGetValue(i, j), oDofStat.iNumNz);
                    Ai.InitDeriv(oDofMap);
                    AddDeriv(Ai, 1., oDofMap, i, j);

                    SP_GRAD_ASSERT(SpGradient::bIsUnique(Ai));
               }
          }

          A = std::move(Atmp);
     }

     template <typename ValueType, typename DERIVED>
     template <util::MatTranspEvalFlag eTransp, typename Func, typename ValueTypeA, index_type NumRowsA, index_type NumColsA>
     void SpMatElemExprBase<ValueType, DERIVED>::ElemAssignCompr(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
          util::ElemAssignHelper<ValueTypeA, ValueType>::template ElemAssignCompr<eTransp, Func, DERIVED, NumRowsA, NumColsA>(A, *this);
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
     template <util::MatTranspEvalFlag eTransp, SpGradCommon::ExprEvalFlags eCompr, typename ValueTypeA, index_type NumRowsA, index_type NumColsA>
     void SpMatMulExpr<LhsValue, RhsValue, LhsExpr, RhsExpr>::Eval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
          typedef util::MatEvalHelperTransp<eTransp> MatEvalType;

          constexpr bool bTransposedEval = eTransp == util::MatTranspEvalFlag::TRANSPOSED;

          static_assert(NumRowsA == (!bTransposedEval ? iNumRowsStatic : iNumColsStatic), "Number of rows does not match");
          static_assert(NumColsA == (!bTransposedEval ? iNumColsStatic : iNumRowsStatic), "Number of columns does not match");

          SpMatrixBase<ValueTypeA, NumRowsA, NumColsA> Atmp;

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

          typedef util::TempExprHelper<LhsExpr, bUseTmpExprLhs> UTmpType;
          typedef util::TempExprHelper<RhsExpr, bUseTmpExprRhs> VTmpType;

          typename UTmpType::Type utmp{UTmpType::EvalUnique(u)};
          typename VTmpType::Type vtmp{VTmpType::EvalUnique(v)};

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

     template <typename ValueType, typename ScalarExpr, index_type NumRows, index_type NumCols>
     void SpMatElemScalarExpr<ValueType, ScalarExpr, NumRows, NumCols>::InsertDof(SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
          SpGradient::InsertDof(u, oExpDofMap);
     }

     template <typename ValueType, typename ScalarExpr, index_type NumRows, index_type NumCols>
     template <typename ValueTypeB>
     void SpMatElemScalarExpr<ValueType, ScalarExpr, NumRows, NumCols>::AddDeriv(ValueTypeB& g, doublereal dCoef, const SpGradExpDofMap& oExpDofMap, index_type i, index_type j) const {
          SpGradient::AddDeriv(u, g, dCoef, oExpDofMap);
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     SpMatrixBase<ValueType, NumRows, NumCols>::SpMatrixBase()
          :SpMatrixBase(pGetNullData()) {
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     SpMatrixBase<ValueType, NumRows, NumCols>::SpMatrixBase(index_type iNumRows, index_type iNumCols, index_type iNumDeriv)
          :SpMatrixBase(pGetNullData()) {

          ResizeReset(iNumRows, iNumCols, iNumDeriv);
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     SpMatrixBase<ValueType, NumRows, NumCols>::SpMatrixBase(const SpMatrixBase& oMat)
          :SpMatrixBase(pGetNullData()) {

          *this = oMat;	// Enable inexpensive shadow copies of SpGradient
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     template <typename ValueTypeExpr, typename Expr>
     SpMatrixBase<ValueType, NumRows, NumCols>::SpMatrixBase(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr)
          :SpMatrixBase(pGetNullData()) {

          constexpr bool bThisIsGradient = std::is_same<ValueType, SpGradient>::value;
          constexpr bool bThisIsDouble = std::is_same<ValueType, doublereal>::value;
          constexpr bool bExprIsGradient = std::is_same<ValueTypeExpr, SpGradient>::value;
          constexpr bool bExprIsDouble = std::is_same<ValueTypeExpr, doublereal>::value;
          static_assert(bThisIsGradient || bThisIsDouble);
          static_assert(bExprIsGradient || bExprIsDouble);
          static_assert(bThisIsGradient || !bExprIsGradient, "Cannot convert SpGradient to doublereal");

          oExpr.template Eval<util::MatTranspEvalFlag::DIRECT, Expr::eExprEvalFlags>(*this);

#ifdef SP_GRAD_DEBUG
          if (Expr::eExprEvalFlags == SpGradCommon::ExprEvalUnique) {
               for (const auto& a: *this) {
                    SP_GRAD_ASSERT(SpGradient::bIsUnique(a));
               }
          }
#endif
          SP_GRAD_ASSERT(bValid());
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     SpMatrixBase<ValueType, NumRows, NumCols>::SpMatrixBase(const std::initializer_list<ValueType>& rgValues)
          :SpMatrixBase(NumRows, NumCols, 0)
     {
          static_assert(NumRows != SpMatrixSize::DYNAMIC);
          static_assert(NumCols != SpMatrixSize::DYNAMIC);

          SP_GRAD_ASSERT(rgValues.size() == static_cast<size_t>(NumRows * NumCols));

          index_type iRow = 1, iCol = 1;

          for (const ValueType& aij: rgValues) {
               GetElem(iRow, iCol) = aij;

               // Insert values column wise like Mat3x3(a11, a21, a31, a12, a22, a32, ...)
               if (++iRow > NumRows) {
                    iRow = 1;
                    ++iCol;
               }
          }

          SP_GRAD_ASSERT(bValid());
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     void SpMatrixBase<ValueType, NumRows, NumCols>::ResizeReset(index_type iNumRows, index_type iNumCols, index_type iNumDeriv) {
          SP_GRAD_ASSERT(bValid());
          SP_GRAD_ASSERT(iNumRows == NumRows || NumRows == SpMatrixSize::DYNAMIC);
          SP_GRAD_ASSERT(iNumCols == NumCols || NumCols == SpMatrixSize::DYNAMIC);

          if (pData->bCheckSize(iNumRows, iNumCols, iNumDeriv)) {
               for (ValueType& g: *pData) {
                    SpGradient::ResizeReset(g, 0, iNumDeriv);
               }
          } else {
               auto pNewData = SpMatrixDataType::pAllocate(iNumRows, iNumCols, iNumDeriv);

               pData->DecRef();
               pData = pNewData;
               pData->IncRef();
          }

          SP_GRAD_ASSERT(bValid());
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     SpMatrixBase<ValueType, NumRows, NumCols>::SpMatrixBase(SpMatrixBase&& oMat)
          :SpMatrixBase(pGetNullData()) {

          *this = std::move(oMat);
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     SpMatrixBase<ValueType, NumRows, NumCols>& SpMatrixBase<ValueType, NumRows, NumCols>::operator=(const SpMatrixBase& oMat) {
          SP_GRAD_ASSERT(bValid());

          if (&oMat == this) {
               return *this;
          }

          ResizeReset(oMat.iGetNumRows(), oMat.iGetNumCols(), 0);

          const index_type iNumElem = pData->iGetNumElem();

          for (index_type i = 1; i <= iNumElem; ++i) {
               GetElem(i) = oMat.GetElem(i); // Enable inexpensive shadow copies of SpGradient
          }

          SP_GRAD_ASSERT(bValid());

          return *this;
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     SpMatrixBase<ValueType, NumRows, NumCols>& SpMatrixBase<ValueType, NumRows, NumCols>::operator=(SpMatrixBase&& oMat) {
          SP_GRAD_ASSERT(bValid());
          SP_GRAD_ASSERT(oMat.bValid());

          std::swap(pData, oMat.pData);

          SP_GRAD_ASSERT(oMat.bValid());
          SP_GRAD_ASSERT(bValid());

          return *this;
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     template <typename ValueTypeExpr, typename Expr>
     SpMatrixBase<ValueType, NumRows, NumCols>& SpMatrixBase<ValueType, NumRows, NumCols>::operator=(const SpMatElemExprBase<ValueTypeExpr, Expr>& oExpr) {
          SP_GRAD_ASSERT(bValid());

          constexpr bool bThisIsGradient = std::is_same<ValueType, SpGradient>::value;
          constexpr bool bThisIsDouble = std::is_same<ValueType, doublereal>::value;
          constexpr bool bExprIsGradient = std::is_same<ValueTypeExpr, SpGradient>::value;
          constexpr bool bExprIsDouble = std::is_same<ValueTypeExpr, doublereal>::value;
          static_assert(bThisIsGradient || bThisIsDouble);
          static_assert(bExprIsGradient || bExprIsDouble);
          static_assert(bThisIsGradient || !bExprIsGradient, "Cannot convert SpGradient to doublereal in assignment");

          oExpr.template Eval<util::MatTranspEvalFlag::DIRECT, Expr::eExprEvalFlags>(*this);

#ifdef SP_GRAD_DEBUG
          if (Expr::eExprEvalFlags == SpGradCommon::ExprEvalUnique) {
               for (const auto& a: *this) {
                    SP_GRAD_ASSERT(SpGradient::bIsUnique(a));
               }
          }
#endif
          SP_GRAD_ASSERT(bValid());

          return *this;
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     template <typename Expr>
     SpMatrixBase<ValueType, NumRows, NumCols>& SpMatrixBase<ValueType, NumRows, NumCols>::operator*=(const SpGradBase<Expr>& b) {
          SP_GRAD_ASSERT(bValid());

          constexpr SpGradCommon::ExprEvalFlags eEvalFlags = util::remove_all<Expr>::type::eExprEvalFlags;
          constexpr bool bThisIsGradient = std::is_same<ValueType, SpGradient>::value;

          static_assert(bThisIsGradient, "Cannot convert SpGradient to doublereal");
          static_assert(!std::is_same<typename util::remove_all<Expr>::type, SpGradient>::value);

          const SpMatElemScalarExpr<SpGradient, const SpGradient, NumRows, NumCols> btmp{b};

          btmp.template AssignEval<util::MatTranspEvalFlag::DIRECT, eEvalFlags, SpGradBinMult>(*this);

          SP_GRAD_ASSERT(bValid());

          return *this;
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     SpMatrixBase<ValueType, NumRows, NumCols>& SpMatrixBase<ValueType, NumRows, NumCols>::operator*=(const SpGradient& b) {
          SP_GRAD_ASSERT(bValid());

          constexpr SpGradCommon::ExprEvalFlags eEvalFlags = SpGradient::eExprEvalFlags;
          constexpr bool bThisIsGradient = std::is_same<ValueType, SpGradient>::value;

          static_assert(bThisIsGradient, "Cannot convert SpGradient to doublereal");

          const SpMatElemScalarExpr<SpGradient, const SpGradient&, NumRows, NumCols> btmp{b};

          btmp.template AssignEval<util::MatTranspEvalFlag::DIRECT, eEvalFlags, SpGradBinMult>(*this);

          SP_GRAD_ASSERT(bValid());

          return *this;
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     SpMatrixBase<ValueType, NumRows, NumCols>& SpMatrixBase<ValueType, NumRows, NumCols>::operator*=(doublereal b) {
          SP_GRAD_ASSERT(bValid());

          for (auto& ai: *this) {
               ai *= b;
          }

          SP_GRAD_ASSERT(bValid());

          return *this;
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     template <typename Expr>
     SpMatrixBase<ValueType, NumRows, NumCols>& SpMatrixBase<ValueType, NumRows, NumCols>::operator/=(const SpGradBase<Expr>& b) {
          SP_GRAD_ASSERT(bValid());

          constexpr SpGradCommon::ExprEvalFlags eEvalFlags = util::remove_all<Expr>::type::eExprEvalFlags;
          constexpr bool bThisIsGradient = std::is_same<ValueType, SpGradient>::value;

          static_assert(bThisIsGradient, "Cannot convert SpGradient to doublereal");
          static_assert(!std::is_same<typename util::remove_all<Expr>::type, SpGradient>::value);

          const SpMatElemScalarExpr<SpGradient, const SpGradient, NumRows, NumCols> btmp{b};

          btmp.template AssignEval<util::MatTranspEvalFlag::DIRECT, eEvalFlags, SpGradBinDiv>(*this);

          SP_GRAD_ASSERT(bValid());

          return *this;
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     SpMatrixBase<ValueType, NumRows, NumCols>& SpMatrixBase<ValueType, NumRows, NumCols>::operator/=(const SpGradient& b) {
          SP_GRAD_ASSERT(bValid());

          constexpr SpGradCommon::ExprEvalFlags eEvalFlags = SpGradient::eExprEvalFlags;
          constexpr bool bThisIsGradient = std::is_same<ValueType, SpGradient>::value;
          static_assert(bThisIsGradient, "Cannot convert SpGradient to doublereal");

          const SpMatElemScalarExpr<SpGradient, const SpGradient&, NumRows, NumCols> btmp{b};

          btmp.template AssignEval<util::MatTranspEvalFlag::DIRECT, eEvalFlags, SpGradBinDiv>(*this);

          SP_GRAD_ASSERT(bValid());

          return *this;
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     SpMatrixBase<ValueType, NumRows, NumCols>& SpMatrixBase<ValueType, NumRows, NumCols>::operator/=(doublereal b) {
          SP_GRAD_ASSERT(bValid());

          for (auto& ai: *this) {
               ai /= b;
          }

          SP_GRAD_ASSERT(bValid());

          return *this;
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     template <typename ValueTypeExpr, typename Expr>
     SpMatrixBase<ValueType, NumRows, NumCols>& SpMatrixBase<ValueType, NumRows, NumCols>::operator+=(const SpMatElemExprBase<ValueTypeExpr, Expr>& b)
     {
          SP_GRAD_ASSERT(bValid());

          constexpr SpGradCommon::ExprEvalFlags eEvalFlags = util::remove_all<Expr>::type::eExprEvalFlags;
          constexpr bool bThisIsGradient = std::is_same<ValueType, SpGradient>::value;
          constexpr bool bThisIsDouble = std::is_same<ValueType, doublereal>::value;
          constexpr bool bExprIsGradient = std::is_same<ValueTypeExpr, SpGradient>::value;
          constexpr bool bExprIsDouble = std::is_same<ValueTypeExpr, doublereal>::value;

          static_assert(bThisIsGradient || bThisIsDouble);
          static_assert(bExprIsGradient || bExprIsDouble);
          static_assert(bThisIsGradient || !bExprIsGradient, "Cannot convert SpGradient to doublereal");

          typedef typename util::remove_all<Expr>::type ExprType;
          static constexpr bool bUseTempExpr = !(ExprType::uMatAccess & util::MatAccessFlag::ELEMENT_WISE);
          typedef typename util::TempExprHelper<decltype(b), bUseTempExpr>::Type TempExpr;

          const TempExpr btmp{b};

          btmp.template AssignEval<util::MatTranspEvalFlag::DIRECT, eEvalFlags, SpGradBinPlus>(*this);

          SP_GRAD_ASSERT(bValid());

          return *this;
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     template <typename ValueTypeExpr, typename Expr>
     SpMatrixBase<ValueType, NumRows, NumCols>& SpMatrixBase<ValueType, NumRows, NumCols>::operator-=(const SpMatElemExprBase<ValueTypeExpr, Expr>& b)
     {
          SP_GRAD_ASSERT(bValid());

          constexpr SpGradCommon::ExprEvalFlags eEvalFlags = util::remove_all<Expr>::type::eExprEvalFlags;
          constexpr bool bThisIsGradient = std::is_same<ValueType, SpGradient>::value;
          constexpr bool bThisIsDouble = std::is_same<ValueType, doublereal>::value;
          constexpr bool bExprIsGradient = std::is_same<ValueTypeExpr, SpGradient>::value;
          constexpr bool bExprIsDouble = std::is_same<ValueTypeExpr, doublereal>::value;

          static_assert(bThisIsGradient || bThisIsDouble);
          static_assert(bExprIsGradient || bExprIsDouble);
          static_assert(bThisIsGradient || !bExprIsGradient, "Cannot convert SpGradient to doublereal");

          typedef typename util::remove_all<Expr>::type ExprType;
          static constexpr bool bUseTempExpr = !(ExprType::uMatAccess & util::MatAccessFlag::ELEMENT_WISE);
          typedef typename util::TempExprHelper<decltype(b), bUseTempExpr>::Type TempExpr;

          const TempExpr btmp{b};

          btmp.template AssignEval<util::MatTranspEvalFlag::DIRECT, eEvalFlags, SpGradBinMinus>(*this);

          SP_GRAD_ASSERT(bValid());

          return *this;
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     SpMatrixBase<ValueType, NumRows, NumCols>::~SpMatrixBase() {
          SP_GRAD_ASSERT(bValid());

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

     template <typename ValueType, index_type NumRows, index_type NumCols>
     void SpMatrixBase<ValueType, NumRows, NumCols>::MakeUnique() {
          *this = EvalUnique(*this);
     }

#ifdef SP_GRAD_DEBUG
     template <typename ValueType, index_type NumRows, index_type NumCols>
     bool SpMatrixBase<ValueType, NumRows, NumCols>::bValid() const {
          const SpMatrixData<ValueType>& oData = *pData; // Do not account for fixed size matrices

          SP_GRAD_ASSERT(oData.iGetRefCnt() >= 1);
          SP_GRAD_ASSERT(oData.iGetNumRows() >= 0);
          SP_GRAD_ASSERT(oData.iGetNumCols() >= 0);
          SP_GRAD_ASSERT(oData.iGetMaxDeriv() >= 0);

          if (&oData != pGetNullData()) {
               SP_GRAD_ASSERT(NumRows == SpMatrixSize::DYNAMIC || NumRows == oData.iGetNumRows());
               SP_GRAD_ASSERT(NumCols == SpMatrixSize::DYNAMIC || NumCols == oData.iGetNumCols());
          } else {
               SP_GRAD_ASSERT(oData.iGetNumRows() == 0);
               SP_GRAD_ASSERT(oData.iGetNumCols() == 0);
               SP_GRAD_ASSERT(oData.iGetMaxDeriv() == 0);
          }

          for (const auto& ai: oData) {
               SP_GRAD_ASSERT(bValid(ai));
          }

          return true;
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     bool SpMatrixBase<ValueType, NumRows, NumCols>::bValid(const SpGradient& g) {
          return g.bValid();
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     bool SpMatrixBase<ValueType, NumRows, NumCols>::bValid(doublereal d) {
          return std::isfinite(d);
     }
#endif

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

     template <typename ValueType, typename ScalarExpr, index_type NumRows, index_type NumCols>
     constexpr SpMatElemScalarExpr<ValueType, ScalarExpr, NumRows, NumCols>::SpMatElemScalarExpr(const ScalarExpr& u) noexcept
          :u(u) {

     }

     template <typename ValueType, typename ScalarExpr, index_type NumRows, index_type NumCols>
     constexpr doublereal SpMatElemScalarExpr<ValueType, ScalarExpr, NumRows, NumCols>::dGetValue(index_type, index_type) const noexcept {
          return SpGradient::dGetValue(u);
     }

     template <typename ValueType, typename ScalarExpr, index_type NumRows, index_type NumCols>
     constexpr index_type SpMatElemScalarExpr<ValueType, ScalarExpr, NumRows, NumCols>::iGetSize(index_type, index_type) const noexcept {
          return SpGradient::iGetSize(u);
     }

     template <typename ValueType, typename ScalarExpr, index_type NumRows, index_type NumCols>
     constexpr index_type SpMatElemScalarExpr<ValueType, ScalarExpr, NumRows, NumCols>::iGetNumRows() noexcept {
          return iNumRowsStatic;
     }

     template <typename ValueType, typename ScalarExpr, index_type NumRows, index_type NumCols>
     constexpr index_type SpMatElemScalarExpr<ValueType, ScalarExpr, NumRows, NumCols>::iGetNumCols() noexcept {
          return iNumColsStatic;
     }

     template <typename ValueType, typename ScalarExpr, index_type NumRows, index_type NumCols>
     constexpr index_type SpMatElemScalarExpr<ValueType, ScalarExpr, NumRows, NumCols>::iGetMaxSize() const noexcept {
          return SpGradient::iGetSize(u);
     }

     template <typename ValueType, typename ScalarExpr, index_type NumRows, index_type NumCols>
     template <typename ValueTypeB>
     void SpMatElemScalarExpr<ValueType, ScalarExpr, NumRows, NumCols>::InsertDeriv(ValueTypeB& g, doublereal dCoef, index_type, index_type) const {
          SpGradient::InsertDeriv(u, g, dCoef);
     }

     template <typename ValueType, typename ScalarExpr, index_type NumRows, index_type NumCols>
     void SpMatElemScalarExpr<ValueType, ScalarExpr, NumRows, NumCols>::GetDofStat(SpGradDofStat& s, index_type, index_type) const noexcept {
          SpGradient::GetDofStat(u, s);
     }

     template <typename ValueType, typename ScalarExpr, index_type NumRows, index_type NumCols>
     template <typename ExprTypeB, typename ExprB>
     constexpr bool SpMatElemScalarExpr<ValueType, ScalarExpr, NumRows, NumCols>::bHaveRefTo(const SpMatElemExprBase<ExprTypeB, ExprB>& A) const noexcept {
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
          return static_cast<SpMatrixDataType*>(SpMatrixData<ValueType>::pGetNullData());
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     template <util::MatTranspEvalFlag eTransp, SpGradCommon::ExprEvalFlags eCompr, typename ValueTypeA, index_type NumRowsA, index_type NumColsA>
     void SpMatrixBase<ValueType, NumRows, NumCols>::Eval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
          this->template ElemEval<eTransp, eCompr>(A);
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     template <util::MatTranspEvalFlag eTransp,
               SpGradCommon::ExprEvalFlags eCompr,
               typename Func,
               typename ValueTypeA,
               index_type NumRowsA, index_type NumColsA>
     void SpMatrixBase<ValueType, NumRows, NumCols>::AssignEval(SpMatrixBase<ValueTypeA, NumRowsA, NumColsA>& A) const {
          SP_GRAD_ASSERT(bValid());
          SP_GRAD_ASSERT(A.bValid());

          this->template ElemAssign<eTransp, eCompr, Func>(A);

          SP_GRAD_ASSERT(bValid());
          SP_GRAD_ASSERT(A.bValid());
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     index_type SpMatrixBase<ValueType, NumRows, NumCols>::iGetMaxSize() const {
          return this->iGetMaxSizeElem();
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     SpMatrixBase<ValueType, NumRows, NumCols>::SpMatrixBase(SpMatrixDataType* pData)
          :pData(pData) {
          pData->IncRef();

          SP_GRAD_ASSERT(bValid());
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
     SpMatrix<ValueType, NumRows, NumCols>::SpMatrix(const std::initializer_list<ValueType>& rgValues)
          :SpMatrixBase<ValueType, NumRows, NumCols>(rgValues) {
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
     SpColVector<ValueType, NumRows>::SpColVector(const std::initializer_list<ValueType>& rgValues)
          :SpMatrixBase<ValueType, NumRows, 1>(rgValues) {
     }

     template <typename ValueType, index_type NumRows>
     void SpColVector<ValueType, NumRows>::ResizeReset(index_type iNumRows, index_type iNumDeriv) {
          SpMatrixBase<ValueType, NumRows, 1>::ResizeReset(iNumRows, 1, iNumDeriv);
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
     SpRowVector<ValueType, NumCols>::SpRowVector(const std::initializer_list<ValueType>& rgValues)
          :SpMatrixBase<ValueType, 1, NumCols>(rgValues) {
     }

     template <typename ValueType, index_type NumCols>
     void SpRowVector<ValueType, NumCols>::ResizeReset(index_type iNumCols, index_type iNumDeriv) {
          SpMatrixBase<ValueType, 1, NumCols>::ResizeReset(1, iNumCols, iNumDeriv);
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
     SpMatElemUniqueExpr<Value, const SpMatElemExprBase<Value, Expr>&>
     EvalUnique(const SpMatElemExprBase<Value, Expr>& A) noexcept {
          return decltype(EvalUnique(A)){A};
     }

     template <typename LhsValue, typename LhsExpr>
     inline constexpr
     SpMatElemBinExpr<typename util::ResultType<LhsValue, SpGradient>::Type,
                      SpGradBinMult,
                      const SpMatElemExprBase<LhsValue, LhsExpr>&,
                      const SpMatElemScalarExpr<SpGradient, const SpGradient&> >
     operator*(const SpMatElemExprBase<LhsValue, LhsExpr>& A,
               const SpGradient& b) noexcept {
          return decltype(operator*(A, b)){A, SpMatElemScalarExpr<SpGradient, const SpGradient&>{b}};
     }

     template <typename LhsValue, typename LhsExpr, typename RhsExpr>
     inline constexpr
     SpMatElemBinExpr<typename util::ResultType<LhsValue, SpGradient>::Type,
                      SpGradBinMult,
                      const SpMatElemExprBase<LhsValue, LhsExpr>&,
                      const SpMatElemScalarExpr<SpGradient, SpGradient> >
     operator*(const SpMatElemExprBase<LhsValue, LhsExpr>& A,
               const SpGradBase<RhsExpr>& b) noexcept {
          static_assert(!std::is_same<typename util::remove_all<RhsExpr>::type, SpGradient>::value);
          return decltype(operator*(A, b)){A, SpMatElemScalarExpr<SpGradient, SpGradient>{SpGradient{EvalUnique(b)}}}; // Avoid multiple evaluations of b!
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

     template <typename RhsValue, typename RhsExpr>
     inline constexpr
     SpMatElemBinExpr<typename util::ResultType<SpGradient, RhsValue>::Type,
                      SpGradBinMult,
                      const SpMatElemScalarExpr<SpGradient, const SpGradient&>,
                      const SpMatElemExprBase<RhsValue, RhsExpr>&>
     operator*(const SpGradient& a,
               const SpMatElemExprBase<RhsValue, RhsExpr>& B) noexcept {
          return decltype(operator*(a, B)){SpMatElemScalarExpr<SpGradient, const SpGradient&>{a}, B};
     }

     template <typename LhsExpr, typename RhsValue, typename RhsExpr>
     inline constexpr
     SpMatElemBinExpr<typename util::ResultType<SpGradient, RhsValue>::Type,
                      SpGradBinMult,
                      const SpMatElemScalarExpr<SpGradient, SpGradient>,
                      const SpMatElemExprBase<RhsValue, RhsExpr>&>
     operator*(const SpGradBase<LhsExpr>& a,
               const SpMatElemExprBase<RhsValue, RhsExpr>& B) noexcept {
          static_assert(!std::is_same<typename util::remove_all<LhsExpr>::type, SpGradient>::value);
          return decltype(operator*(a, B)){SpMatElemScalarExpr<SpGradient, SpGradient>{SpGradient{EvalUnique(a)}}, B}; // Avoid multiple evaluations of a!
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

     template <typename LhsValue, typename LhsExpr>
     inline constexpr
     SpMatElemBinExpr<typename util::ResultType<LhsValue, SpGradient>::Type,
                      SpGradBinDiv,
                      const SpMatElemExprBase<LhsValue, LhsExpr>&,
                      const SpMatElemScalarExpr<SpGradient, const SpGradient&> >
     operator/(const SpMatElemExprBase<LhsValue, LhsExpr>& A,
               const SpGradient& b) noexcept {
          return decltype(operator/(A, b)){A, SpMatElemScalarExpr<SpGradient, const SpGradient&>{b}};
     }

     template <typename LhsValue, typename LhsExpr, typename RhsExpr>
     inline constexpr
     SpMatElemBinExpr<typename util::ResultType<LhsValue, SpGradient>::Type,
                      SpGradBinDiv,
                      const SpMatElemExprBase<LhsValue, LhsExpr>&,
                      const SpMatElemScalarExpr<SpGradient, SpGradient> >
     operator/(const SpMatElemExprBase<LhsValue, LhsExpr>& A,
               const SpGradBase<RhsExpr>& b) noexcept {
          return decltype(operator/(A, b)){A, SpMatElemScalarExpr<SpGradient, SpGradient>{SpGradient{EvalUnique(b)}}}; // Avoid multiple evaluations of b!
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
          static_assert(u.iNumRowsStatic == v.iNumRowsStatic);
          static_assert(u.iNumColsStatic == 1);
          static_assert(v.iNumColsStatic == 1);
          SP_GRAD_ASSERT(u.iGetNumRows() == v.iGetNumRows());

          typedef const SpMatElemExprBase<LhsValue, LhsExpr>& LhsTmpExpr;
          typedef const SpMatElemExprBase<RhsValue, RhsExpr>& RhsTmpExpr;
          typedef typename util::remove_all<LhsTmpExpr>::type LhsExprType;
          typedef typename util::remove_all<RhsTmpExpr>::type RhsExprType;

          constexpr bool bLhsUsesIterators = (LhsExprType::uMatAccess & util::MatAccessFlag::ITERATORS) != 0;
          constexpr bool bRhsUsesIterators = (RhsExprType::uMatAccess & util::MatAccessFlag::ITERATORS) != 0;

          constexpr bool bUseTmpExprLhs = !(LhsExprType::iNumElemOps == 0 && bLhsUsesIterators);
          constexpr bool bUseTmpExprRhs = !(RhsExprType::iNumElemOps == 0 && bRhsUsesIterators);

          typedef util::TempExprHelper<LhsTmpExpr, bUseTmpExprLhs> UTmpType;
          typedef util::TempExprHelper<RhsTmpExpr, bUseTmpExprRhs> VTmpType;

          typename UTmpType::Type utmp{UTmpType::EvalUnique(u)};
          typename VTmpType::Type vtmp{VTmpType::EvalUnique(v)};

          static_assert(utmp.iNumElemOps == 0);
          static_assert(vtmp.iNumElemOps == 0);
          static_assert((utmp.uMatAccess & util::MatAccessFlag::ITERATORS) != 0);
          static_assert((vtmp.uMatAccess & util::MatAccessFlag::ITERATORS) != 0);

          const index_type iRowOffsetU = utmp.iGetRowOffset();
          const index_type iColSizeU = utmp.iGetRowOffset() * utmp.iGetNumRows();
          const index_type iRowOffsetV = vtmp.iGetRowOffset();
          const index_type iColSizeV = vtmp.iGetNumRows() * iRowOffsetV;

          typename util::ResultType<LhsValue, RhsValue>::Type a;

          util::MapInnerProduct(a,
                                utmp.begin(),
                                utmp.begin() + iColSizeU,
                                iRowOffsetU,
                                vtmp.begin(),
                                vtmp.begin() + iColSizeV,
                                iRowOffsetV);

          return a;
     }

     template <typename Value, typename Expr>
     inline constexpr Value Norm(const SpMatElemExprBase<Value, Expr>& u) {
          using std::sqrt;
          return sqrt(Dot(u, u));
     }

     template <typename ValueType, typename Expr>
     inline constexpr SpSubMatDynExpr<ValueType, const SpMatElemExprBase<ValueType, Expr>&>
     SubMatrix(const SpMatElemExprBase<ValueType, Expr>& A, index_type iRowStart, index_type iRowStep, index_type iNumRows, index_type iColStart, index_type iColStep, index_type iNumCols) {
          return decltype(SubMatrix(A, iRowStart, iRowStep, iNumRows, iColStart, iColStep, iNumCols))(A, iRowStart, iRowStep, iNumRows, iColStart, iColStep, iNumCols);
     }

     template <index_type iRowStart, index_type iRowStep, index_type iNumRows, index_type iColStart, index_type iColStep, index_type iNumCols, typename ValueType, typename Expr>
     inline constexpr SpSubMatStatExpr<ValueType, const SpMatElemExprBase<ValueType, Expr>&, iRowStart, iRowStep, iNumRows, iColStart, iColStep, iNumCols>
     SubMatrix(const SpMatElemExprBase<ValueType, Expr>& A) {
          return decltype(SubMatrix<iRowStart, iRowStep, iNumRows, iColStart, iColStep, iNumCols>(A))(A);
     }

     template <index_type iNumRows, index_type iNumCols, typename ValueType, typename Expr>
     inline constexpr SpSubMatStatResExpr<ValueType, const SpMatElemExprBase<ValueType, Expr>&, iNumRows, iNumCols>
     SubMatrix(const SpMatElemExprBase<ValueType, Expr>& A, index_type iRowStart, index_type iRowStep, index_type iColStart, index_type iColStep) {
          return decltype(SubMatrix<iNumRows, iNumCols>(A, iRowStart, iRowStep, iColStart, iColStep))(A, iRowStart, iRowStep, iColStart, iColStep);
     }

     template <index_type iRowStart, index_type iRowStep, index_type iNumRows, typename ValueType, typename Expr>
     inline constexpr SpSubMatStatRowExpr<ValueType, const SpMatElemExprBase<ValueType, Expr>&, iRowStart, iRowStep, iNumRows>
     SubMatrix(const SpMatElemExprBase<ValueType, Expr>& A, index_type iColStart, index_type iColStep, index_type iNumCols) {
          return decltype(SubMatrix<iRowStart, iRowStep, iNumRows>(A, iColStart, iColStep, iNumCols))(A, iColStart, iColStep, iNumCols);
     }

     template <index_type iRowStart, index_type iRowStep, index_type iNumRows, typename ValueType, typename Expr>
     inline constexpr SpSubMatStatExpr<ValueType, const SpMatElemExprBase<ValueType, Expr>&, iRowStart, iRowStep, iNumRows, 1, 1, 1>
     SubColVector(const SpMatElemExprBase<ValueType, Expr>& A) {
          static_assert(A.iNumColsStatic == 1);
          return decltype(SubColVector<iRowStart, iRowStep, iNumRows>(A))(A);
     }

     template <index_type iColStart, index_type iColStep, index_type iNumCols, typename ValueType, typename Expr>
     inline constexpr SpSubMatStatExpr<ValueType, const SpMatElemExprBase<ValueType, Expr>&, 1, 1, 1, iColStart, iColStep, iNumCols>
     SubRowVector(const SpMatElemExprBase<ValueType, Expr>& A) {
          static_assert(A.iNumRowsStatic == 1);
          return decltype(SubRowVector<iColStart, iColStep, iNumCols>(A))(A);
     }

     template <typename ValueType>
     SpMatrix<ValueType, 3, 3> MatGVec(const SpColVector<ValueType, 3>& g) {
          const ValueType d = (4./(4.+Dot(g, g)));

          SpMatrix<ValueType, 3, 3> G(3, 3, 0);

          G(1, 1) = d;
          G(2, 1) = g(3) * d / 2.;
          G(3, 1) = -g(2) * d / 2.;
          G(1, 2) = -g(3) * d / 2.;
          G(2, 2) = d;
          G(3, 2) = g(1) * d / 2.;
          G(1, 3) = g(2) * d / 2.;
          G(2, 3) = -g(1) * d / 2.;
          G(3, 3) = d;

          return G;
     }

     template <typename ValueType>
     inline SpMatrix<ValueType, 3, 3>
     MatRVec(const SpColVector<ValueType, 3>& g) {
          SpMatrix<ValueType, 3, 3> RDelta(3, 3, 0);

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
     inline ValueType RotCo(const ValueType& phi) {
          // This algorithm is a simplified version of RotCo in RotCoeff.hc
          // from Marco Morandini  <morandini@aero.polimi.it>
          // and Teodoro Merlini  <merlini@aero.polimi.it>
          using std::sin;
          using std::cos;
          using std::sqrt;
          using std::fabs;

          constexpr index_type N = 10;
          ValueType phip[N];
          ValueType phi2{EvalUnique(phi * phi)};

          if (fabs(phi) < RotCoeff::SerThrsh[0]) {
               SpGradient::ResizeReset(phip[0], 1., 0);

               for (index_type j = 1; j <= N - 1; j++) {
                    phip[j] = EvalUnique(phip[j - 1] * phi2);
               }

               ValueType cf;
               SpGradient::ResizeReset(cf, 0., SpGradient::iGetSize(phip[N - 1]));

               for (index_type j = 0; j < RotCoeff::SerTrunc[0]; j++) {
                    cf += EvalUnique(phip[j] / RotCoeff::SerCoeff[0][j]);
               }

               return cf;
          }

          const ValueType pd{sqrt(phi2)};
          return sin(pd) / pd;                 // a = sin(phi)/phi
     }

     template <typename ValueType>
     inline SpColVector<ValueType, 3>
     VecRotMat(const SpMatrix<ValueType, 3, 3>& R) {
          // Modified from Appendix 2.4 of
          //
          // author = {Marco Borri and Lorenzo Trainelli and Carlo L. Bottasso},
          // title = {On Representations and Parameterizations of Motion},
          // journal = {Multibody System Dynamics},
          // volume = {4},
          // pages = {129--193},
          // year = {2000}

          SpColVector<ValueType, 3> unit(3, 0);

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
               unit /= RotCo(phi);
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
                    SpColVector<ValueType, 1> xi = SubMatrix<1, 1>(Cross(unit, R.GetCol(i)), i, 1, 1, 1);
                    sinphi -= xi(1) * 0.5;
               }

               unit *= atan2(sinphi, cosphi);
          }

          return unit;
     }

     template <typename ValueType>
     inline SpMatrix<ValueType, 3, 3>
     MatCrossVec(const SpColVector<ValueType, 3>& v, doublereal d) {
          SpMatrix<ValueType, 3, 3> A(3, 3, v.iGetMaxSize());

          SpGradient::ResizeReset(A(1, 1), d, 0);
          SpGradient::ResizeReset(A(2, 2), d, 0);
          SpGradient::ResizeReset(A(3, 3), d, 0);

          A(1, 2) = -v(3);
          A(1, 3) =  v(2);
          A(2, 1) =  v(3);
          A(2, 3) = -v(1);
          A(3, 1) = -v(2);
          A(3, 2) =  v(1);

          return A;
     }

     template <typename T>
     inline T Det(const SpMatrix<T, 2, 2>& A) {
          return A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1);
     }

     template <typename T>
     inline SpMatrix<T, 2, 2> Inv(const SpMatrix<T, 2, 2>& A) {
          const T detA = Det(A);

          return SpMatrix<T, 2, 2>{A(2, 2) / detA,
                    -A(2, 1) / detA,
                    -A(1, 2) / detA,
                    A(1, 1) / detA};
     }

     template <typename ValueType, index_type NumRows, index_type NumCols>
     inline std::ostream& operator<<(std::ostream& os, const SpMatrixBase<ValueType, NumRows, NumCols>& A) {
          for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
               for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
                    os << A.GetElem(i, j) << ' ';
               }
          }

          return os;
     }
}

#endif
