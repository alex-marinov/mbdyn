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

#ifndef __SP_GRADIENT_EXPR__H_INCLUDED__
#define __SP_GRADIENT_EXPR__H_INCLUDED__

#include "sp_gradient_base.h"
#include "sp_gradient_util.h"

namespace sp_grad {
     template <typename BinaryFunc, typename LhsExpr, typename RhsExpr>
     class SpGradBinExpr: public SpGradBase<SpGradBinExpr<BinaryFunc, LhsExpr, RhsExpr> > {
     public:
          static constexpr SpGradCommon::ExprEvalFlags
          eExprEvalFlags = util::ExprEvalFlagsHelper<std::remove_reference<LhsExpr>::type::eExprEvalFlags,
                                                     std::remove_reference<RhsExpr>::type::eExprEvalFlags>::eExprEvalFlags;

          constexpr SpGradBinExpr(const LhsExpr& u, const RhsExpr&  v)
               :u(u),
                v(v),
                f(BinaryFunc::f(u.dGetValue(), v.dGetValue())),
                df_du(BinaryFunc::df_du(u.dGetValue(), v.dGetValue())),
                df_dv(BinaryFunc::df_dv(u.dGetValue(), v.dGetValue())) {
          }

          constexpr doublereal dGetValue() const {
               return f;
          }

          constexpr index_type iGetSize() const {
               return u.iGetSize() + v.iGetSize();
          }

          void GetDofStat(SpGradDofStat& s) const {
               u.GetDofStat(s);
               v.GetDofStat(s);
          }

          void InsertDeriv(SpGradient& g, doublereal dCoef) const {
               u.InsertDeriv(g, df_du * dCoef);
               v.InsertDeriv(g, df_dv * dCoef);
          }

          template <typename ExprOther>
          constexpr bool bHaveRefTo(const SpGradBase<ExprOther>& g) const {
               return u.bHaveRefTo(g) || v.bHaveRefTo(g);
          }

          void InsertDof(SpGradExpDofMap& oExpDofMap) const {
               u.InsertDof(oExpDofMap);
               v.InsertDof(oExpDofMap);
          }

          void AddDeriv(SpGradient& g, doublereal dCoef, const SpGradExpDofMap& oExpDofMap) const {
               u.AddDeriv(g, df_du * dCoef, oExpDofMap);
               v.AddDeriv(g, df_dv * dCoef, oExpDofMap);
          }

#ifdef SP_GRAD_DEBUG
          void PrintValue(std::ostream& os) const {
               BinaryFunc::Print(os);
               os << '(';
               u.PrintValue(os);
               os << ',';
               v.PrintValue(os);
               os << ')';
          }

          void PrintDeriv(std::ostream& os, doublereal dCoef) const {
               BinaryFunc::Print(os);
               os << '(';
               u.PrintDeriv(os, df_du * dCoef);
               os << ',';
               v.PrintDeriv(os, df_dv * dCoef);
               os << ')';
          }
#endif

     private:
          const LhsExpr u;
          const RhsExpr v;
          const doublereal f;
          const doublereal df_du, df_dv;
     };

     template <typename UnaryFunc, typename Expr>
     class SpGradUnExpr: public SpGradBase<SpGradUnExpr<UnaryFunc, Expr> > {
     public:
          static constexpr SpGradCommon::ExprEvalFlags eExprEvalFlags = std::remove_reference<Expr>::type::eExprEvalFlags;

          explicit constexpr SpGradUnExpr(const Expr&  u)
               :u(u),
                f(UnaryFunc::f(u.dGetValue())),
                df_du(UnaryFunc::df_du(u.dGetValue())) {
          }

          constexpr doublereal dGetValue() const {
               return f;
          }

          constexpr index_type iGetSize() const {
               return u.iGetSize();
          }

          void GetDofStat(SpGradDofStat& s) const {
               u.GetDofStat(s);
          }

          void InsertDeriv(SpGradient& g, doublereal dCoef) const {
               u.InsertDeriv(g, df_du * dCoef);
          }

          template <typename ExprOther>
          constexpr bool bHaveRefTo(const SpGradBase<ExprOther>& g) const {
               return u.bHaveRefTo(g);
          }

          void InsertDof(SpGradExpDofMap& oExpDofMap) const {
               u.InsertDof(oExpDofMap);
          }

          void AddDeriv(SpGradient& g, doublereal dCoef, const SpGradExpDofMap& oExpDofMap) const {
               u.AddDeriv(g, df_du * dCoef, oExpDofMap);
          }

#ifdef SP_GRAD_DEBUG
          void PrintValue(std::ostream& os) const {
               UnaryFunc::Print(os);
               os << '(';
               u.PrintValue(os);
               os << ')';
          }

          void PrintDeriv(std::ostream& os, doublereal dCoef) const {
               UnaryFunc::Print(os);
               os << '(';
               u.PrintDeriv(os, df_du * dCoef);
               os << ')';
          }
#endif
     private:
          const Expr u;
          const doublereal f;
          const doublereal df_du;
     };

     template <typename ScalarType>
     class SpGradConstExpr: public SpGradBase<SpGradConstExpr<ScalarType> > {
     public:
          static constexpr SpGradCommon::ExprEvalFlags eExprEvalFlags = SpGradCommon::ExprEvalDuplicate;

          explicit constexpr SpGradConstExpr(ScalarType  u)
               :u(u) {
          }

          constexpr ScalarType dGetValue() const {
               return u;
          }

          static constexpr index_type iGetSize() {
               return 0;
          }

          static void GetDofStat(SpGradDofStat&) {
          }

          static void InsertDeriv(SpGradient&, doublereal) {
          }

          template <typename Expr>
          static constexpr bool bHaveRefTo(const SpGradBase<Expr>&) {
               return false;
          }

          static void InsertDof(SpGradExpDofMap&) {
          }

          static void AddDeriv(SpGradient&, doublereal, const SpGradExpDofMap&) {
          }

#ifdef SP_GRAD_DEBUG
          void PrintValue(std::ostream& os) const {
               os << u << ' ';
          }

          void PrintDeriv(std::ostream& os, doublereal dCoef) const {
               os << 0 << ' ';
          }
#endif

     private:
          const ScalarType u;
     };

     template <typename Expr>
     class SpGradComprExpr: public SpGradBase<SpGradComprExpr<Expr> > {
     public:
          static constexpr SpGradCommon::ExprEvalFlags eExprEvalFlags = SpGradCommon::ExprEvalUnique;

          constexpr SpGradComprExpr(const Expr&  u)
               :u(u) {
          }

          constexpr doublereal dGetValue() const {
               return u.dGetValue();
          }

          constexpr index_type iGetSize() const {
               return u.iGetSize();
          }

          void GetDofStat(SpGradDofStat& s) const {
               u.GetDofStat(s);
          }

          void InsertDeriv(SpGradient& g, doublereal dCoef) const {
               u.InsertDeriv(g, dCoef);
          }

          template <typename ExprOther>
          constexpr bool bHaveRefTo(const SpGradBase<ExprOther>& g) const {
               return u.bHaveRefTo(g);
          }

          void InsertDof(SpGradExpDofMap& oExpDofMap) const {
               u.InsertDof(oExpDofMap);
          }

          void AddDeriv(SpGradient& g, doublereal dCoef, const SpGradExpDofMap& oExpDofMap) const {
               u.AddDeriv(g, dCoef, oExpDofMap);
          }

     private:
          const Expr u;
     };

     template <typename BinaryFunc, typename LhsExpr, typename RhsExpr>
     class SpGradBoolExpr: public SpGradBase<SpGradBoolExpr<BinaryFunc, LhsExpr, RhsExpr> > {
     public:
          static constexpr SpGradCommon::ExprEvalFlags
          eExprEvalFlags = util::ExprEvalFlagsHelper<std::remove_reference<LhsExpr>::type::eExprEvalFlags,
                                                     std::remove_reference<RhsExpr>::type::eExprEvalFlags>::eExprEvalFlags;

          constexpr SpGradBoolExpr(const LhsExpr& u, const RhsExpr&  v)
               :f(BinaryFunc::f(u.dGetValue(), v.dGetValue())) {
          }

          operator bool() const {
               return f;
          }

          constexpr bool dGetValue() const {
               return f;
          }

          static constexpr index_type iGetSize() {
               return 0;
          }

          static void GetDofStat(SpGradDofStat& s) {
          }

          static void InsertDeriv(SpGradient& g, doublereal dCoef) {
          }

          template <typename ExprOther>
          static constexpr bool bHaveRefTo(const SpGradBase<ExprOther>& g) {
               return false;
          }

          static void InsertDof(SpGradExpDofMap& oExpDofMap) {
          }

          static void AddDeriv(SpGradient& g, doublereal dCoef, const SpGradExpDofMap& oExpDofMap) {
          }

#ifdef SP_GRAD_DEBUG
          static void PrintValue(std::ostream& os) {
          }

          static void PrintDeriv(std::ostream& os, doublereal dCoef) {
          }
#endif

     private:
          const bool f;
     };

     namespace util {
          template <typename BinaryFunc, bool bIsScalarConstLhs, bool bIsScalarConstRhs>
          struct GpGradProdBinExprEval {
               static void dEval() = delete;
          };

          template <typename BinaryFunc>
          struct GpGradProdBinExprEval<BinaryFunc, false, false> {
               template <typename LhsExpr, typename RhsExpr>
               static doublereal dEval(const LhsExpr& u, const RhsExpr& v, doublereal& dDer) {
                    const doublereal uv = u.dGetValue();
                    const doublereal vv = v.dGetValue();
                    const doublereal ud = u.dGetDeriv();
                    const doublereal vd = v.dGetDeriv();

                    doublereal dVal = BinaryFunc::f(uv, vv);

                    dDer = BinaryFunc::df_du(uv, vv) * ud + BinaryFunc::df_dv(uv, vv) * vd;

                    return dVal;
               }
          };

          template <typename BinaryFunc>
          struct GpGradProdBinExprEval<BinaryFunc, true, false> {
               template <typename LhsExpr, typename RhsExpr>
               static doublereal dEval(const LhsExpr& u, const RhsExpr& v, doublereal& dDer) {
                    const doublereal uv = u.dGetValue();
                    const doublereal vv = v.dGetValue();
                    const doublereal vd = v.dGetDeriv();

                    static_assert(u.dGetDeriv() == 0., "u must be a scalar constant");

                    doublereal dVal = BinaryFunc::f(uv, vv);

                    dDer = BinaryFunc::df_dv(uv, vv) * vd;

                    return dVal;
               }
          };

          template <typename BinaryFunc>
          struct GpGradProdBinExprEval<BinaryFunc, false, true> {
               template <typename LhsExpr, typename RhsExpr>
               static doublereal dEval(const LhsExpr& u, const RhsExpr& v, doublereal& dDer) {
                    const doublereal uv = u.dGetValue();
                    const doublereal vv = v.dGetValue();
                    const doublereal ud = u.dGetDeriv();

                    static_assert(v.dGetDeriv() == 0., "v must be a scalar constant");

                    doublereal dVal = BinaryFunc::f(uv, vv);

                    dDer = BinaryFunc::df_du(uv, vv) * ud;

                    return dVal;
               }
          };
     }

     template <typename BinaryFunc, typename LhsExpr, typename RhsExpr>
     class GpGradProdBinExpr: public GpGradProdBase<GpGradProdBinExpr<BinaryFunc, LhsExpr, RhsExpr> > {
     public:
          GpGradProdBinExpr(const LhsExpr& u, const RhsExpr&  v) {
               typedef typename util::remove_all<LhsExpr>::type LhsExprType;
               typedef typename util::remove_all<RhsExpr>::type RhsExprType;
               constexpr bool bIsScalarConstLhs = LhsExprType::bIsScalarConst;
               constexpr bool bIsScalarConstRhs = RhsExprType::bIsScalarConst;

               static_assert(!(bIsScalarConstLhs && bIsScalarConstRhs), "At least one argument must be no constant");

               typedef util::GpGradProdBinExprEval<BinaryFunc, bIsScalarConstLhs, bIsScalarConstRhs> ExprEvalHelper;

               dVal = ExprEvalHelper::dEval(u, v, dDer);
          }

          constexpr doublereal dGetValue() const {
               return dVal;
          }

          constexpr doublereal dGetDeriv() const {
               return dDer;
          }

          static constexpr bool bIsScalarConst = false;

     private:
          doublereal dVal;
          doublereal dDer;
     };

     template <typename UnaryFunc, typename Expr>
     class GpGradProdUnExpr: public GpGradProdBase<GpGradProdUnExpr<UnaryFunc, Expr> > {
     public:
          explicit GpGradProdUnExpr(const Expr& u) {
               const doublereal uv = u.dGetValue();
               const doublereal ud = u.dGetDeriv();

               dVal = UnaryFunc::f(uv);
               dDer = UnaryFunc::df_du(uv) * ud;
          }

          constexpr doublereal dGetValue() const {
               return dVal;
          }

          constexpr doublereal dGetDeriv() const {
               return dDer;
          }

          static constexpr bool bIsScalarConst = false;
     private:
          doublereal dVal;
          doublereal dDer;
     };

     template <typename ScalarType>
     class GpGradProdConstExpr: public GpGradProdBase<GpGradProdConstExpr<ScalarType>> {
     public:
          constexpr explicit GpGradProdConstExpr(doublereal dVal)
               :dVal(dVal) {
          }

          constexpr ScalarType dGetValue() const {
               return dVal;
          }

          static constexpr doublereal dGetDeriv() {
               return 0.;
          }

          static constexpr bool bIsScalarConst = true;
     private:
          const ScalarType dVal;
     };

     template <typename BinaryFunc, typename LhsExpr, typename RhsExpr>
     class GpGradProdBoolExpr: public GpGradProdBase<GpGradProdBoolExpr<BinaryFunc, LhsExpr, RhsExpr>> {
     public:
          constexpr GpGradProdBoolExpr(const LhsExpr& u, const RhsExpr& v)
               :bVal(BinaryFunc::f(u.dGetValue(), v.dGetValue())) {
          }

          constexpr operator bool() const {
               return bVal;
          }

          constexpr bool dGetValue() const {
               return bVal;
          }

          static constexpr doublereal dGetDeriv() {
               return 0.;
          }

          static constexpr bool bIsScalarConst = true;
     private:
          const bool bVal;
     };
}
#endif
