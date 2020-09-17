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
	  static constexpr SpGradCommon::ExprEvalFlags eExprEvalFlags = SpGradCommon::ExprEvalUncompressed;
                
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
	  static constexpr SpGradCommon::ExprEvalFlags eExprEvalFlags = SpGradCommon::ExprEvalCompressed;
                
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

}
#endif
