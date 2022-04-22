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
        Copyright (C) 2020(-2022) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#ifndef __SP_GRADIENT_OP__H_INCLUDED__
#define __SP_GRADIENT_OP__H_INCLUDED__

#include "sp_gradient_func.h"
#include "sp_gradient.h"

namespace sp_grad {
#define SP_GRADIENT_DEFINE_BINARY_OPERATOR_CONST_ARG_RHS(SP_GRAD_EXPR, SP_GRAD_OP_FUNC, SP_GRAD_OP_CLASS, SP_CONST_ARG_TYPE) \
     template <typename LhsExpr>					\
     constexpr inline SP_GRAD_EXPR<SP_GRAD_OP_CLASS, const SpGradBase<LhsExpr>&, const SpGradConstExpr<SP_CONST_ARG_TYPE> > \
     SP_GRAD_OP_FUNC(const SpGradBase<LhsExpr>& u, SP_CONST_ARG_TYPE v) noexcept { \
	  return decltype(SP_GRAD_OP_FUNC(u, v))(u, SpGradConstExpr<SP_CONST_ARG_TYPE>(v)); \
     }

#define SP_GRADIENT_DEFINE_BINARY_OPERATOR_CONST_ARG_LHS(SP_GRAD_EXPR, SP_GRAD_OP_FUNC, SP_GRAD_OP_CLASS, SP_CONST_ARG_TYPE) \
     template <typename RhsExpr>					\
     constexpr inline SP_GRAD_EXPR<SP_GRAD_OP_CLASS, const SpGradConstExpr<SP_CONST_ARG_TYPE>, const SpGradBase<RhsExpr>& > \
     SP_GRAD_OP_FUNC(SP_CONST_ARG_TYPE u, const SpGradBase<RhsExpr>& v) noexcept { \
	  return decltype(SP_GRAD_OP_FUNC(u, v))(SpGradConstExpr<SP_CONST_ARG_TYPE>(u), v); \
     }

#define SP_GRADIENT_DEFINE_BINARY_OPERATOR(SP_GRAD_EXPR, SP_GRAD_OP_FUNC, SP_GRAD_OP_CLASS) \
     SP_GRADIENT_DEFINE_BINARY_OPERATOR_CONST_ARG_LHS(SP_GRAD_EXPR, SP_GRAD_OP_FUNC, SP_GRAD_OP_CLASS, doublereal) \
     SP_GRADIENT_DEFINE_BINARY_OPERATOR_CONST_ARG_RHS(SP_GRAD_EXPR, SP_GRAD_OP_FUNC, SP_GRAD_OP_CLASS, doublereal) \
									\
     template <typename LhsExpr, typename RhsExpr>			\
     constexpr inline SP_GRAD_EXPR<SP_GRAD_OP_CLASS, const SpGradBase<LhsExpr>&, const SpGradBase<RhsExpr>& > \
     SP_GRAD_OP_FUNC(const SpGradBase<LhsExpr>& u, const SpGradBase<RhsExpr>& v) noexcept { \
	  return decltype(SP_GRAD_OP_FUNC(u, v))(u, v);			\
     }

#define SP_GRADIENT_DEFINE_UNARY_OPERATOR(SP_GRAD_OP_FUNC, SP_GRAD_OP_CLASS) \
     template <typename Expr>						\
     constexpr inline SpGradUnExpr<SP_GRAD_OP_CLASS, const SpGradBase<Expr>& > \
     SP_GRAD_OP_FUNC(const SpGradBase<Expr>& u) noexcept {		\
	  return decltype(SP_GRAD_OP_FUNC(u))(u);			\
     }

     SP_GRADIENT_DEFINE_BINARY_OPERATOR(SpGradBinExpr, operator +, SpGradBinPlus)
     SP_GRADIENT_DEFINE_BINARY_OPERATOR(SpGradBinExpr, operator -, SpGradBinMinus)
     SP_GRADIENT_DEFINE_BINARY_OPERATOR(SpGradBinExpr, operator *, SpGradBinMult)
     SP_GRADIENT_DEFINE_BINARY_OPERATOR(SpGradBinExpr, operator /, SpGradBinDiv)
     SP_GRADIENT_DEFINE_BINARY_OPERATOR(SpGradBinExpr, pow, SpGradBinPow)
     SP_GRADIENT_DEFINE_BINARY_OPERATOR(SpGradBinExpr, atan2, SpGradBinAtan2)
     SP_GRADIENT_DEFINE_BINARY_OPERATOR(SpGradBinExpr, copysign, SpGradBinCopysign)
     SP_GRADIENT_DEFINE_BINARY_OPERATOR(SpGradBinExpr, fmod, SpGradBinFmod)

     SP_GRADIENT_DEFINE_BINARY_OPERATOR_CONST_ARG_RHS(SpGradBinExpr, pow, SpGradBinPowInt, integer)

     SP_GRADIENT_DEFINE_BINARY_OPERATOR(SpGradBoolExpr, operator <, SpGradBoolLessThan)
     SP_GRADIENT_DEFINE_BINARY_OPERATOR(SpGradBoolExpr, operator <=, SpGradBoolLessEqual)
     SP_GRADIENT_DEFINE_BINARY_OPERATOR(SpGradBoolExpr, operator >, SpGradBoolGreaterThan)
     SP_GRADIENT_DEFINE_BINARY_OPERATOR(SpGradBoolExpr, operator >=, SpGradBoolGreaterEqual)
     SP_GRADIENT_DEFINE_BINARY_OPERATOR(SpGradBoolExpr, operator ==, SpGradBoolEqualTo)
     SP_GRADIENT_DEFINE_BINARY_OPERATOR(SpGradBoolExpr, operator !=, SpGradBoolNotEqualTo)

     SP_GRADIENT_DEFINE_UNARY_OPERATOR(operator-, SpGradUnaryMinus)
     SP_GRADIENT_DEFINE_UNARY_OPERATOR(fabs, SpGradFabs)
     SP_GRADIENT_DEFINE_UNARY_OPERATOR(sqrt, SpGradSqrt)
     SP_GRADIENT_DEFINE_UNARY_OPERATOR(exp, SpGradExp)
     SP_GRADIENT_DEFINE_UNARY_OPERATOR(log, SpGradLog)
     SP_GRADIENT_DEFINE_UNARY_OPERATOR(sin, SpGradSin)
     SP_GRADIENT_DEFINE_UNARY_OPERATOR(cos, SpGradCos)
     SP_GRADIENT_DEFINE_UNARY_OPERATOR(tan, SpGradTan)
     SP_GRADIENT_DEFINE_UNARY_OPERATOR(sinh, SpGradSinh)
     SP_GRADIENT_DEFINE_UNARY_OPERATOR(cosh, SpGradCosh)
     SP_GRADIENT_DEFINE_UNARY_OPERATOR(tanh, SpGradTanh)
     SP_GRADIENT_DEFINE_UNARY_OPERATOR(asin, SpGradAsin)
     SP_GRADIENT_DEFINE_UNARY_OPERATOR(acos, SpGradAcos)
     SP_GRADIENT_DEFINE_UNARY_OPERATOR(atan, SpGradAtan)
     SP_GRADIENT_DEFINE_UNARY_OPERATOR(asinh, SpGradAsinh)
     SP_GRADIENT_DEFINE_UNARY_OPERATOR(acosh, SpGradAcosh)
     SP_GRADIENT_DEFINE_UNARY_OPERATOR(atanh, SpGradAtanh)

     template <typename Expr>
     constexpr inline SpGradComprExpr<const SpGradBase<Expr>&>
     EvalUnique(const SpGradBase<Expr>& g) noexcept {
	  return decltype(EvalUnique(g))(g);
     }

     inline constexpr doublereal EvalUnique(doublereal d) noexcept { return d; }

     template <typename Expr>
     constexpr inline const GpGradProdBase<Expr>&
     EvalUnique(const GpGradProdBase<Expr>& g) noexcept {
	  return g;
     }

     constexpr inline const GpGradProd&
     EvalUnique(const GpGradProd& g) noexcept {
          return g;
     }
     
#ifdef SP_GRAD_DEBUG
     template <typename Expr>
     std::ostream& operator<<(std::ostream& os, const SpGradBase<Expr>& g) {
	  os << "f=(";
	  g.PrintValue(os);
	  os << ")\n";
	  os << "df=(";
	  g.PrintDeriv(os, 1.);
	  os << ")";

	  return os;
     }

     std::ostream& operator<<(std::ostream& os, const SpGradient& g);
#endif

#define GRAD_PROD_DEFINE_BINARY_OPERATOR_CONST_ARG_RHS(SP_GRAD_EXPR, SP_GRAD_OP_FUNC, SP_GRAD_OP_CLASS, SP_CONST_ARG_TYPE) \
     template <typename LhsExpr>					\
     constexpr inline SP_GRAD_EXPR<SP_GRAD_OP_CLASS, const GpGradProdBase<LhsExpr>&, const GpGradProdConstExpr<SP_CONST_ARG_TYPE> > \
     SP_GRAD_OP_FUNC(const GpGradProdBase<LhsExpr>& u, SP_CONST_ARG_TYPE v) noexcept { \
	  return decltype(SP_GRAD_OP_FUNC(u, v))(u, GpGradProdConstExpr<SP_CONST_ARG_TYPE>(v)); \
     }

#define GRAD_PROD_DEFINE_BINARY_OPERATOR_CONST_ARG_LHS(SP_GRAD_EXPR, SP_GRAD_OP_FUNC, SP_GRAD_OP_CLASS, SP_CONST_ARG_TYPE) \
     template <typename RhsExpr>					\
     constexpr inline SP_GRAD_EXPR<SP_GRAD_OP_CLASS, const GpGradProdConstExpr<SP_CONST_ARG_TYPE>, const GpGradProdBase<RhsExpr>& > \
     SP_GRAD_OP_FUNC(SP_CONST_ARG_TYPE u, const GpGradProdBase<RhsExpr>& v) noexcept { \
	  return decltype(SP_GRAD_OP_FUNC(u, v))(GpGradProdConstExpr<SP_CONST_ARG_TYPE>(u), v); \
     }

#define GRAD_PROD_DEFINE_BINARY_OPERATOR(SP_GRAD_EXPR, SP_GRAD_OP_FUNC, SP_GRAD_OP_CLASS) \
     GRAD_PROD_DEFINE_BINARY_OPERATOR_CONST_ARG_LHS(SP_GRAD_EXPR, SP_GRAD_OP_FUNC, SP_GRAD_OP_CLASS, doublereal) \
     GRAD_PROD_DEFINE_BINARY_OPERATOR_CONST_ARG_RHS(SP_GRAD_EXPR, SP_GRAD_OP_FUNC, SP_GRAD_OP_CLASS, doublereal) \
									\
     template <typename LhsExpr, typename RhsExpr>			\
     constexpr inline SP_GRAD_EXPR<SP_GRAD_OP_CLASS, const GpGradProdBase<LhsExpr>&, const GpGradProdBase<RhsExpr>& > \
     SP_GRAD_OP_FUNC(const GpGradProdBase<LhsExpr>& u, const GpGradProdBase<RhsExpr>& v) noexcept { \
	  return decltype(SP_GRAD_OP_FUNC(u, v))(u, v);			\
     }

#define GRAD_PROD_DEFINE_UNARY_OPERATOR(SP_GRAD_OP_FUNC, SP_GRAD_OP_CLASS) \
     template <typename Expr>						\
     constexpr inline GpGradProdUnExpr<SP_GRAD_OP_CLASS, const GpGradProdBase<Expr>& > \
     SP_GRAD_OP_FUNC(const GpGradProdBase<Expr>& u) noexcept {		\
	  return decltype(SP_GRAD_OP_FUNC(u))(u);			\
     }

     GRAD_PROD_DEFINE_BINARY_OPERATOR(GpGradProdBinExpr, operator +, SpGradBinPlus)
     GRAD_PROD_DEFINE_BINARY_OPERATOR(GpGradProdBinExpr, operator -, SpGradBinMinus)
     GRAD_PROD_DEFINE_BINARY_OPERATOR(GpGradProdBinExpr, operator *, SpGradBinMult)
     GRAD_PROD_DEFINE_BINARY_OPERATOR(GpGradProdBinExpr, operator /, SpGradBinDiv)
     GRAD_PROD_DEFINE_BINARY_OPERATOR(GpGradProdBinExpr, pow, SpGradBinPow)
     GRAD_PROD_DEFINE_BINARY_OPERATOR(GpGradProdBinExpr, atan2, SpGradBinAtan2)
     GRAD_PROD_DEFINE_BINARY_OPERATOR(GpGradProdBinExpr, copysign, SpGradBinCopysign)
     GRAD_PROD_DEFINE_BINARY_OPERATOR(GpGradProdBinExpr, fmod, SpGradBinFmod)

     GRAD_PROD_DEFINE_BINARY_OPERATOR_CONST_ARG_RHS(GpGradProdBinExpr, pow, SpGradBinPowInt, integer)

     GRAD_PROD_DEFINE_BINARY_OPERATOR(GpGradProdBoolExpr, operator <, SpGradBoolLessThan)
     GRAD_PROD_DEFINE_BINARY_OPERATOR(GpGradProdBoolExpr, operator <=, SpGradBoolLessEqual)
     GRAD_PROD_DEFINE_BINARY_OPERATOR(GpGradProdBoolExpr, operator >, SpGradBoolGreaterThan)
     GRAD_PROD_DEFINE_BINARY_OPERATOR(GpGradProdBoolExpr, operator >=, SpGradBoolGreaterEqual)
     GRAD_PROD_DEFINE_BINARY_OPERATOR(GpGradProdBoolExpr, operator ==, SpGradBoolEqualTo)
     GRAD_PROD_DEFINE_BINARY_OPERATOR(GpGradProdBoolExpr, operator !=, SpGradBoolNotEqualTo)

     GRAD_PROD_DEFINE_UNARY_OPERATOR(operator-, SpGradUnaryMinus)
     GRAD_PROD_DEFINE_UNARY_OPERATOR(fabs, SpGradFabs)
     GRAD_PROD_DEFINE_UNARY_OPERATOR(sqrt, SpGradSqrt)
     GRAD_PROD_DEFINE_UNARY_OPERATOR(exp, SpGradExp)
     GRAD_PROD_DEFINE_UNARY_OPERATOR(log, SpGradLog)
     GRAD_PROD_DEFINE_UNARY_OPERATOR(sin, SpGradSin)
     GRAD_PROD_DEFINE_UNARY_OPERATOR(cos, SpGradCos)
     GRAD_PROD_DEFINE_UNARY_OPERATOR(tan, SpGradTan)
     GRAD_PROD_DEFINE_UNARY_OPERATOR(sinh, SpGradSinh)
     GRAD_PROD_DEFINE_UNARY_OPERATOR(cosh, SpGradCosh)
     GRAD_PROD_DEFINE_UNARY_OPERATOR(tanh, SpGradTanh)
     GRAD_PROD_DEFINE_UNARY_OPERATOR(asin, SpGradAsin)
     GRAD_PROD_DEFINE_UNARY_OPERATOR(acos, SpGradAcos)
     GRAD_PROD_DEFINE_UNARY_OPERATOR(atan, SpGradAtan)
     GRAD_PROD_DEFINE_UNARY_OPERATOR(asinh, SpGradAsinh)
     GRAD_PROD_DEFINE_UNARY_OPERATOR(acosh, SpGradAcosh)
     GRAD_PROD_DEFINE_UNARY_OPERATOR(atanh, SpGradAtanh)

#ifdef SP_GRAD_DEBUG
     std::ostream& operator<<(std::ostream& os, const GpGradProd& g);
#endif
}
#endif
