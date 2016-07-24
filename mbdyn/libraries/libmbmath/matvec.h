/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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
 AUTHOR: Reinhard Resch <r.resch@secop.com>
        Copyright (C) 2013(-2015) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#ifndef ___MAT_VEC_H__INCLUDED___
#define ___MAT_VEC_H__INCLUDED___

#include <cassert>
#include <iostream>
#include "myassert.h"
#include "gradient.h"
#include "matvec3.h" // FIXME: We need that for compatibility reasons
#include "matvec6.h"

#ifndef MATVEC_DEBUG
    #ifdef DEBUG
	    #define MATVEC_DEBUG 1
    #else
	    #define MATVEC_DEBUG 0
    #endif
#endif

#if MATVEC_DEBUG == 0 || defined(DEBUG)
	#define MATVEC_ASSERT(expr) ASSERT(expr)
#elif MATVEC_DEBUG > 0
	#define MATVEC_ASSERT(expr) assert(expr)
#endif

namespace grad {

template <typename T>
struct VectorSize
{
	static const int N = -1; // Note must not be zero because that would mean variable size!
};

template <>
struct VectorSize<scalar_func_type>
{
	static const int N = 1;
};

template <>
struct VectorSize<Vec3>
{
	static const int N = 3;
};

template <>
struct VectorSize<Vec6>
{
	static const int N = 6;
};

template <typename T, index_type N_rows, index_type N_cols>
class Matrix;

template <typename T, index_type N_rows>
class Vector;

/**
 * This is a default implementation that handles built in types like float doublereal and long doublereal
 * and also classes which do not need expression templates
 */
template <typename ScalarBinaryFunction, typename T, typename ScalarLhsExpr, typename ScalarRhsExpr>
class GenericBinaryExpression {
public:
	static const bool bAlias = false;
	static const index_type iMaxDerivatives = 0;
	static const bool bVectorize = true;
	static const index_type iDimension = -1;
	typedef char vector_deriv_type;
	typedef T ScalarType;
	typedef T ExpressionType;

	GenericBinaryExpression(const ScalarLhsExpr& lhs, const ScalarRhsExpr& rhs)
		:a(ScalarBinaryFunction::f(lhs, rhs)) {
	};

	operator ExpressionType() const {
		return a;
	}

    scalar_func_type dGetValue() const {
        return a;
    }

    scalar_deriv_type dGetDerivativeLocal(index_type iLocalDof) const {
        return 0;
    }

    vector_deriv_type dGetDerivativeLocalVector(index_type iLocalVecDof) const {
        return vector_deriv_type();
    }

    index_type iGetStartIndexLocal() const {
        return std::numeric_limits<index_type>::max();
    }

    index_type iGetEndIndexLocal() const {
        return 0;
    }

    index_type iGetStartIndexLocalVector() const {
        return std::numeric_limits<index_type>::max();
    }

    index_type iGetEndIndexLocalVector() const {
        return 0;
    }

    LocalDofMap* pGetDofMap() const {
        return 0;
    }

    bool bHaveReferenceTo(const void* p) const {
        return false;
    }

    static index_type iGetMaxDerivatives() {
    	return iMaxDerivatives;
    }

    void Compute() const {}

private:
	const ExpressionType a;
};

template <typename ScalarUnaryFunction, typename T, typename ScalarExpr>
class GenericUnaryExpression {
public:
	static const bool bAlias = false;
	static const index_type iMaxDerivatives = 0;
	static const bool bVectorize = true;
	static const index_type iDimension = -1;
	typedef char vector_deriv_type;
	typedef T ScalarType;
	typedef T ExpressionType;

	GenericUnaryExpression(const ScalarExpr& expr)
		:a(ScalarUnaryFunction::f(expr)) {
	};

	operator ExpressionType() const {
		return a;
	}

    scalar_func_type dGetValue() const {
        return a;
    }

    scalar_deriv_type dGetDerivativeLocal(index_type iLocalDof) const {
        return 0;
    }

    vector_deriv_type dGetDerivativeLocalVector(index_type iLocalVecDof) const {
        return vector_deriv_type();
    }

    index_type iGetStartIndexLocal() const {
        return std::numeric_limits<index_type>::max();
    }

    index_type iGetEndIndexLocal() const {
        return 0;
    }

    index_type iGetStartIndexLocalVector() const {
        return std::numeric_limits<index_type>::max();
    }

    index_type iGetEndIndexLocalVector() const {
        return 0;
    }

    LocalDofMap* pGetDofMap() const {
        return 0;
    }

    bool bHaveReferenceTo(const void* p) const {
        return false;
    }

    static index_type iGetMaxDerivatives() {
    	return iMaxDerivatives;
    }

    void Compute() const {}

private:
	const ExpressionType a;
};

template <typename ScalarTypeLhs, typename ScalarTypeRhs>
struct CommonScalarType {
private:
	typedef void ScalarType;
};

template <index_type N_SIZE>
struct CommonScalarType<Gradient<N_SIZE>, Gradient<N_SIZE> > {
	typedef Gradient<N_SIZE> ScalarType;
};

template <index_type N_SIZE>
struct CommonScalarType<Gradient<N_SIZE>, scalar_func_type> {
	typedef Gradient<N_SIZE> ScalarType;
};

template <index_type N_SIZE>
struct CommonScalarType<scalar_func_type, Gradient<N_SIZE> > {
	typedef Gradient<N_SIZE> ScalarType;
};

template <>
struct CommonScalarType<scalar_func_type, scalar_func_type> {
	typedef scalar_func_type ScalarType;
};

template <typename T>
struct BasicScalarType {
	typedef T ScalarType;
};

template <typename BinFunc, typename LhsExpr, typename RhsExpr>
struct BasicScalarType<BinaryExpr<BinFunc, LhsExpr, RhsExpr> > {
	typedef typename CommonScalarType<typename BasicScalarType<LhsExpr>::ScalarType, typename BasicScalarType<RhsExpr>::ScalarType>::ScalarType ScalarType;
};

template <typename UnFunc, typename Expr>
struct BasicScalarType<UnaryExpr<UnFunc, Expr> > {
	typedef typename BasicScalarType<Expr>::ScalarType ScalarType;
};

template <typename Expression>
struct BasicScalarType<GradientExpression<Expression> > {
	typedef typename BasicScalarType<Expression>::ScalarType ScalarType;
};

template <index_type N_SIZE, bool ALIAS>
struct BasicScalarType<DirectExpr<Gradient<N_SIZE>, ALIAS> > {
	typedef Gradient<N_SIZE> ScalarType;
};

template <index_type N_SIZE>
struct BasicScalarType<ConstExpr<Gradient<N_SIZE> > > {
	typedef scalar_func_type ScalarType;
};

template <typename ScalarBinaryFunction, typename T, typename ScalarLhsExpr, typename ScalarRhsExpr>
struct BasicScalarType<GenericBinaryExpression<ScalarBinaryFunction, T, ScalarLhsExpr, ScalarRhsExpr> > {
	typedef typename GenericBinaryExpression<ScalarBinaryFunction, T, ScalarLhsExpr, ScalarRhsExpr>::ScalarType ScalarType;
};

template <typename ScalarUnaryFunction, typename T, typename ScalarExpr>
struct BasicScalarType<GenericUnaryExpression<ScalarUnaryFunction, T, ScalarExpr> > {
	typedef typename GenericUnaryExpression<ScalarUnaryFunction, T, ScalarExpr>::ScalarType ScalarType;
};

template <typename T>
struct ScalarTypeTraits {
	typedef T ScalarType;
	typedef T DirectExpressionType;

	template <typename Expression, typename PointerType>
	static bool bHaveReferenceTo(const Expression&, const PointerType*, const PointerType*) {
		return false;
	}
};

template <index_type N_SIZE>
struct ScalarTypeTraits<Gradient<N_SIZE> > {
	typedef Gradient<N_SIZE> ScalarType;
	typedef GradientExpression<DirectExpr<Gradient<N_SIZE> > > DirectExpressionType;

	template <typename Expression, typename PointerType>
	static bool bHaveReferenceTo(const GradientExpression<Expression>& g, const PointerType* pFirst, const PointerType* pLast) {
		for (const PointerType* p = pFirst; p <= pLast; ++p) {
			if (g.bHaveReferenceTo(p)) {
				return true;
			}
		}

		return false;
	}
};

template <typename ScalarBinaryFunction, typename T, typename ScalarLhsExpr, typename ScalarRhsExpr>
struct ScalarBinaryExpressionTraits: ScalarTypeTraits<T> {
	typedef GenericBinaryExpression<ScalarBinaryFunction, T, ScalarLhsExpr, ScalarRhsExpr> ExpressionType;
};

template <typename ScalarBinaryFunction, index_type N_SIZE, typename ScalarLhsExpr, typename ScalarRhsExpr>
struct ScalarBinaryExpressionTraits<ScalarBinaryFunction, Gradient<N_SIZE>, ScalarLhsExpr, ScalarRhsExpr>: ScalarTypeTraits<Gradient<N_SIZE> > {
	typedef GradientExpression<BinaryExpr<ScalarBinaryFunction, ScalarLhsExpr, ScalarRhsExpr> > ExpressionType;
};

template <typename ScalarBinaryFunction, index_type N_SIZE, typename ScalarLhsExpr>
struct ScalarBinaryExpressionTraits<ScalarBinaryFunction, Gradient<N_SIZE>, ScalarLhsExpr, Gradient<N_SIZE> >: ScalarTypeTraits<Gradient<N_SIZE> > {
	typedef GradientExpression<BinaryExpr<ScalarBinaryFunction, ScalarLhsExpr, DirectExpr<Gradient<N_SIZE> > > > ExpressionType;
};

template <typename ScalarBinaryFunction, index_type N_SIZE, typename ScalarLhsExpr>
struct ScalarBinaryExpressionTraits<ScalarBinaryFunction, Gradient<N_SIZE>, ScalarLhsExpr, scalar_func_type>: ScalarTypeTraits<Gradient<N_SIZE> > {
	typedef GradientExpression<BinaryExpr<ScalarBinaryFunction, ScalarLhsExpr, ConstExpr<Gradient<N_SIZE> > > > ExpressionType;
};

template <typename ScalarBinaryFunction, index_type N_SIZE, typename ScalarRhsExpr>
struct ScalarBinaryExpressionTraits<ScalarBinaryFunction, Gradient<N_SIZE>, Gradient<N_SIZE>, ScalarRhsExpr>: ScalarTypeTraits<Gradient<N_SIZE> > {
	typedef GradientExpression<BinaryExpr<ScalarBinaryFunction, DirectExpr<Gradient<N_SIZE> >, ScalarRhsExpr> > ExpressionType;
};

template <typename ScalarBinaryFunction, index_type N_SIZE, typename ScalarRhsExpr>
struct ScalarBinaryExpressionTraits<ScalarBinaryFunction, Gradient<N_SIZE>, scalar_func_type, ScalarRhsExpr>: ScalarTypeTraits<Gradient<N_SIZE> > {
	typedef GradientExpression<BinaryExpr<ScalarBinaryFunction, ConstExpr<Gradient<N_SIZE> >, ScalarRhsExpr> > ExpressionType;
};

template <typename ScalarBinaryFunction, index_type N_SIZE>
struct ScalarBinaryExpressionTraits<ScalarBinaryFunction, Gradient<N_SIZE>, Gradient<N_SIZE>, Gradient<N_SIZE> >: ScalarTypeTraits<Gradient<N_SIZE> > {
	typedef GradientExpression<BinaryExpr<ScalarBinaryFunction, DirectExpr<Gradient<N_SIZE> >, DirectExpr<Gradient<N_SIZE> > > > ExpressionType;
};

template <typename ScalarBinaryFunction, index_type N_SIZE>
struct ScalarBinaryExpressionTraits<ScalarBinaryFunction, Gradient<N_SIZE>, Gradient<N_SIZE>, scalar_func_type>: ScalarTypeTraits<Gradient<N_SIZE> > {
	typedef GradientExpression<BinaryExpr<ScalarBinaryFunction, DirectExpr<Gradient<N_SIZE> >, ConstExpr<Gradient<N_SIZE> > > > ExpressionType;
};

template <typename ScalarBinaryFunction, index_type N_SIZE>
struct ScalarBinaryExpressionTraits<ScalarBinaryFunction, Gradient<N_SIZE>, scalar_func_type, Gradient<N_SIZE> >: ScalarTypeTraits<Gradient<N_SIZE> > {
	typedef GradientExpression<BinaryExpr<ScalarBinaryFunction, ConstExpr<Gradient<N_SIZE> >, DirectExpr<Gradient<N_SIZE> > > > ExpressionType;
};

template <typename ScalarUnaryFunction, typename T, typename ScalarExpr>
struct ScalarUnaryExpressionTraits: ScalarTypeTraits<T> {
	typedef GenericUnaryExpression<ScalarUnaryFunction, T, ScalarExpr> ExpressionType;
};

template <typename ScalarUnaryFunction, index_type N_SIZE, typename ScalarExpr>
struct ScalarUnaryExpressionTraits<ScalarUnaryFunction, Gradient<N_SIZE>, ScalarExpr>: ScalarTypeTraits<Gradient<N_SIZE> > {
	typedef GradientExpression<UnaryExpr<ScalarUnaryFunction, ScalarExpr> > ExpressionType;
};

template <typename ScalarUnaryFunction, index_type N_SIZE>
struct ScalarUnaryExpressionTraits<ScalarUnaryFunction, Gradient<N_SIZE>, Gradient<N_SIZE> >: ScalarTypeTraits<Gradient<N_SIZE> > {
	typedef GradientExpression<UnaryExpr<ScalarUnaryFunction, DirectExpr<Gradient<N_SIZE> > > > ExpressionType;
};

template <typename ScalarBinaryFunction, typename ScalarLhsExpr, typename ScalarRhsExpr>
class ScalarBinaryOperation {
private:
	typedef typename BasicScalarType<ScalarLhsExpr>::ScalarType ScalarTypeLhs;
	typedef typename BasicScalarType<ScalarRhsExpr>::ScalarType ScalarTypeRhs;

public:
	typedef typename CommonScalarType<ScalarTypeLhs, ScalarTypeRhs>::ScalarType ScalarType;

    typedef typename ScalarBinaryExpressionTraits<ScalarBinaryFunction,
    					ScalarType,
    					ScalarLhsExpr,
    					ScalarRhsExpr>::ExpressionType ExpressionType;

    static ExpressionType f(const ScalarLhsExpr& u, const ScalarRhsExpr& v) {
        return ExpressionType(u, v);
    }
};

template <typename ScalarUnaryFunction, typename ScalarExpr>
class ScalarUnaryOperation {
public:
	typedef typename BasicScalarType<ScalarExpr>::ScalarType ScalarType;

    typedef typename ScalarUnaryExpressionTraits<ScalarUnaryFunction,
    					ScalarType,
    					ScalarExpr>::ExpressionType ExpressionType;

    static ExpressionType f(const ScalarExpr& u) {
        return ExpressionType(u);
    }
};

/**
 * Check if two arrays are overlapping
 *
 * pFirstArray points to the first element in the array.
 * pLastArray points to the last element in the array, not beyond the last element like in STL.
 *
 *           ^         ^
 * [ 0,  1, (2,  3, 4, 5), 6 ]
 * [(0,  1,  2), 3, 4, 5,  6 ]
 *   ^       ^
 */
template <typename ScalarType>
inline bool bArrayOverlap(const ScalarType* pFirstArray1,
				   const ScalarType* pLastArray1,
				   const ScalarType* pFirstArray2,
				   const ScalarType* pLastArray2) {
	MATVEC_ASSERT(pLastArray1 >= pFirstArray1);
	MATVEC_ASSERT(pLastArray2 >= pFirstArray2);

	return (pFirstArray1 >= pFirstArray2 && pFirstArray1 <= pLastArray2)
		|| (pLastArray1 >= pFirstArray2 && pLastArray1 <= pLastArray2)
		|| (pFirstArray2 >= pFirstArray1 && pFirstArray2 <= pLastArray1)
		|| (pLastArray2 >= pFirstArray1 && pLastArray2 <= pLastArray1);
}

template <typename ScalarType1, typename ScalarType2>
inline bool bArrayOverlap(const ScalarType1*, const ScalarType1*, const ScalarType2*, const ScalarType2*) {
	// Self reference with different data types is not possible since unions are not supported
	return false;
}

template <typename T>
inline void ZeroInit(T* first, T* last) {
	NO_OP;
}

template <>
inline void ZeroInit<float>(float* first, float* last) {
	array_fill(first, last, 0.0F);
}
 
template <>
inline void ZeroInit<double>(double* first, double* last) {
	array_fill(first, last, 0.0);
}

template <>
inline void ZeroInit<long double>(long double* first, long double* last) {
	array_fill(first, last, 0.0L);
}

/**
 * FIXME: In order to reduce the number of matching function calls
 * for template operators, we have to provide the number of rows and columns
 * twice in some situations.
 *
 * In order to check for inconsistent definitions
 * of the number of rows and columns at compile time,
 * the following struct is used.
 */
template <long DIFF>
struct IndexCheck {
private:
	/* Whenever we get an compilation error here,
	 * the value of iNumRows or iNumCols in VectorExpression or MatrixExpression
	 * is not consistent with their template argument Expression
	 */
	enum CheckType {INDEX_CHECK};
};

template <>
struct IndexCheck<0L> {
	enum CheckType {INDEX_CHECK};
};

namespace MatVecHelp
{
	template <typename T>
	struct AliasTypeHelper
	{

	};

	template <>
	struct AliasTypeHelper<scalar_func_type>
	{
		static const bool bAlias = false;
	};

	template <typename Expression>
	struct AliasTypeHelper<GradientExpression<Expression> >
	{
		static const bool bAlias = GradientExpression<Expression>::bAlias;
	};

	template <index_type N_SIZE, bool ALIAS>
	struct AliasTypeHelper<DirectExpr<Gradient<N_SIZE>, ALIAS> >
	{
		static const bool bAlias = DirectExpr<Gradient<N_SIZE>, ALIAS>::bAlias;
	};

	template <index_type N_SIZE>
	struct AliasTypeHelper<Gradient<N_SIZE> >
	{
		static const bool bAlias = false;
	};

	template <typename BinFunc, typename LhsExpr, typename RhsExpr>
	struct AliasTypeHelper<BinaryExpr<BinFunc, LhsExpr, RhsExpr> >
	{
		static const bool bAlias = BinaryExpr<BinFunc, LhsExpr, RhsExpr>::bAlias;
	};

	template <typename UnFunc, typename Expr>
	struct AliasTypeHelper<UnaryExpr<UnFunc, Expr> >
	{
		static const bool bAlias = UnaryExpr<UnFunc, Expr>::bAlias;
	};

	template <bool bAlias>
	struct ApplyAliasHelperMatrix;

	template <>
	struct ApplyAliasHelperMatrix<false>
	{
	    template <typename MatrixType, typename Func, typename Expression>
	    static inline void ApplyMatrixFunc(MatrixType& A, const Expression& B, const Func& f) {
	    	A.ApplyMatrixFuncNoAlias(B, f);
	    }
	};

	template <>
	struct ApplyAliasHelperMatrix<true>
	{
	    template <typename MatrixType, typename Func, typename Expression>
		static inline void ApplyMatrixFunc(MatrixType& A, const Expression& B, const Func& f) {
			A.ApplyMatrixFuncAlias(B, f);
		}
	};

	struct Assign
	{
		template <typename T, typename U>
		static inline void Eval(T& b, const U& a) {
			b = a;
		}
	};

	struct Add
	{
		template <typename T, typename U>
		static inline void Eval(T& b, const U& a) {
			b += a;
		}
	};

	struct Sub
	{
		template <typename T, typename U>
		static inline void Eval(T& b, const U& a) {
			b -= a;
		}
	};

	struct Mul
	{
		template <typename T, typename U>
		static inline void Eval(T& b, const U& a) {
			b *= a;
		}
	};

	struct Div
	{
		template <typename T, typename U>
		static inline void Eval(T& b, const U& a) {
			b /= a;
		}
	};
}

template <typename Expression, index_type N_rows>
class VectorExpression: public Expression {
public:
	static const index_type iNumRows = N_rows;
	typedef typename Expression::ScalarType ScalarType;
	typedef typename Expression::ExpressionType ExpressionType;

	explicit VectorExpression(const Expression& e)
		:Expression(e) {
#if MATVEC_DEBUG > 0
		AssertValid();
#endif
	}

    template <typename Expr>
    explicit VectorExpression(const Expr& u)
        :Expression(u) {
#if MATVEC_DEBUG > 0
		AssertValid();
#endif
    }
    
    template <typename LhsExpr, typename RhsExpr>
    explicit VectorExpression(const LhsExpr& u, const RhsExpr& v)
        :Expression(u, v) {
#if MATVEC_DEBUG > 0
		AssertValid();
#endif
    }

    static index_type iGetNumRows() { return iNumRows; }

private:
#if MATVEC_DEBUG > 0
    void AssertValid() const {
    	MATVEC_ASSERT(Expression::iGetNumRows() == iNumRows);
    }
#endif
    typedef typename IndexCheck<iNumRows - Expression::iNumRows>::CheckType check_iNumRows;
};

template <typename Expression, index_type N_rows, index_type N_cols, bool CLEAR_ALIAS=false>
class MatrixExpression: public Expression {
public:
	static const bool bAlias = CLEAR_ALIAS ? false : Expression::bAlias;
    static const index_type iNumRows = N_rows;
    static const index_type iNumCols = N_cols;
	typedef typename Expression::ScalarType ScalarType;
	typedef typename Expression::ExpressionType ExpressionType;

	explicit MatrixExpression(const Expression& e)
		:Expression(e) {
#if MATVEC_DEBUG > 0
		AssertValid();
#endif
	}

    template <typename Expr>
    explicit MatrixExpression(const Expr& u)
        :Expression(u) {
#if MATVEC_DEBUG > 0
		AssertValid();
#endif
    }

    template <typename LhsExpr, typename RhsExpr>
    explicit MatrixExpression(const LhsExpr& u, const RhsExpr& v)
        :Expression(u, v) {
#if MATVEC_DEBUG > 0
		AssertValid();
#endif
    }

    static index_type iGetNumRows() {
    	return iNumRows;
    }
    static index_type iGetNumCols() {
    	return iNumCols;
    }

private:

#if MATVEC_DEBUG > 0
    void AssertValid() const {
    	MATVEC_ASSERT(Expression::iGetNumRows() == iNumRows);
        MATVEC_ASSERT(Expression::iGetNumCols() == iNumCols);
    }
#endif

    typedef typename IndexCheck<iNumRows - Expression::iNumRows>::CheckType check_iNumRows;
    typedef typename IndexCheck<iNumCols - Expression::iNumCols>::CheckType check_iNumCols;
};

/**
 * This class handles binary expressions of the form
 * f(vector1, vector2) = vector3
 */
template <typename ScalarBinFunc, typename VectorLhsExpr, typename VectorRhsExpr>
class VectorVectorVectorBinaryExpr {
public:
	static const bool bAlias = VectorLhsExpr::bAlias || VectorRhsExpr::bAlias;
	static const index_type iNumRows = VectorLhsExpr::iNumRows;
	typedef typename ScalarBinFunc::ScalarType ScalarType;
    typedef typename ScalarBinFunc::ExpressionType ExpressionType;

    VectorVectorVectorBinaryExpr(const VectorLhsExpr& u, const VectorRhsExpr& v)
        :oU(u), oV(v) {
            
    }
    
    ExpressionType operator()(index_type i) const {
    	MATVEC_ASSERT(i >= 1);
    	MATVEC_ASSERT(i <= iGetNumRows());
        return ScalarBinFunc::f(oU(i), oV(i));
    }
    
    index_type iGetNumRows() const {
        MATVEC_ASSERT(oU.iGetNumRows() == oV.iGetNumRows());
        MATVEC_ASSERT(oU.iGetNumRows() == iNumRows);
        return oU.iGetNumRows();
    }
    
    template <typename ScalarType2>
    bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
    	return oU.bHaveReferenceTo(pFirst, pLast) || oV.bHaveReferenceTo(pFirst, pLast);
    }

private:
    const VectorLhsExpr oU;
    const VectorRhsExpr oV;
};

/**
 * This class handles expressions of the form
 * f(vector1, scalar1) = vector2
 */
template <typename ScalarBinFunc, typename VectorLhsExpr, typename ScalarRhsExpr>
class VectorScalarVectorBinaryExpr {
public:
	static const bool bAlias = VectorLhsExpr::bAlias || MatVecHelp::AliasTypeHelper<ScalarRhsExpr>::bAlias;
	static const index_type iNumRows = VectorLhsExpr::iNumRows;
	typedef typename ScalarBinFunc::ScalarType ScalarType;
    typedef typename ScalarBinFunc::ExpressionType ExpressionType;

    VectorScalarVectorBinaryExpr(const VectorLhsExpr& u, const ScalarRhsExpr& v)
        :oU(u), oV(v) {

    }

    ExpressionType operator()(index_type i) const {
    	MATVEC_ASSERT(i >= 1);
    	MATVEC_ASSERT(i <= iGetNumRows());
        return ScalarBinFunc::f(oU(i), oV);
    }

    index_type iGetNumRows() const {
    	MATVEC_ASSERT(oU.iGetNumRows() == iNumRows);
        return oU.iGetNumRows();
    }

    template <typename ScalarType2>
    bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
    	typedef ScalarTypeTraits<BasicScalarType<ScalarRhsExpr> > ScalarTraits;
    	return oU.bHaveReferenceTo(pFirst, pLast)
    			|| ScalarTraits::bHaveReferenceTo(oV, pFirst, pLast);
    }

private:
    const VectorLhsExpr oU;
    const ScalarRhsExpr oV;
};

/**
 * This class handles expressions of the form
 * f(vector1) = vector2
 */
template <typename ScalarUnaryFunc, typename VectorExpr>
class VectorVectorUnaryExpr {
public:
	static const bool bAlias = VectorExpr::bAlias;
	static const index_type iNumRows = VectorExpr::iNumRows;
	typedef typename ScalarUnaryFunc::ScalarType ScalarType;
    typedef typename ScalarUnaryFunc::ExpressionType ExpressionType;

    VectorVectorUnaryExpr(const VectorExpr& u)
        :oU(u) {

    }

    ExpressionType operator()(index_type i) const {
    	MATVEC_ASSERT(i >= 1);
    	MATVEC_ASSERT(i <= iNumRows);
        return ScalarUnaryFunc::f(oU(i));
    }

    index_type iGetNumRows() const {
    	MATVEC_ASSERT(oU.iGetNumRows() == iNumRows);
        return oU.iGetNumRows();
    }

    template <typename ScalarType2>
    bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
    	return oU.bHaveReferenceTo(pFirst, pLast);
    }

private:
    const VectorExpr oU;
};

template <index_type iStartIndex, index_type iEndIndex, typename VectorExpr>
class SubVectorExpr {
public:
	static const bool bAlias = VectorExpr::bAlias;
	static const index_type iNumRows = iEndIndex - iStartIndex + 1;
	typedef typename VectorExpr::ScalarType ScalarType;
    typedef typename VectorExpr::ExpressionType ExpressionType;

    SubVectorExpr(const VectorExpr& u)
        :oU(u) {

    }

    ExpressionType operator()(index_type i) const {
    	MATVEC_ASSERT(i >= 1);
    	MATVEC_ASSERT(i <= iNumRows);
        return oU(i + iStartIndex - 1);
    }

    index_type iGetNumRows() const {
    	MATVEC_ASSERT(oU.iGetNumRows() == VectorExpr::iNumRows);
        return iNumRows;
    }

    template <typename ScalarType2>
    bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
    	return oU.bHaveReferenceTo(pFirst, pLast);
    }

private:
    typedef typename MaxSizeCheck<iEndIndex <= VectorExpr::iNumRows>::CheckType check_iEndIndex;
    typedef typename MaxSizeCheck<iStartIndex >= 1>::CheckType check_iStartIndex;
    const VectorExpr oU;
};

template <typename VectorType, bool ALIAS=false>
class VectorDirectExpr {
public:
	static const bool bAlias = ALIAS;
	static const index_type iNumRows = VectorType::iNumRows;
    typedef typename VectorType::ScalarType ScalarType;
    typedef typename VectorType::ExpressionType ExpressionType;

    VectorDirectExpr(const VectorType& u)
        :oU(u) {
    }
    
    const ScalarType& operator()(index_type i) const {
    	MATVEC_ASSERT(i >= 1);
    	MATVEC_ASSERT(i <= iNumRows);
        return oU(i);
    }
    
    index_type iGetNumRows() const {
    	MATVEC_ASSERT(oU.iGetNumRows() == iNumRows);
        return oU.iGetNumRows();
    }
    
    template <typename ScalarType2>
    bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
    	return oU.bHaveReferenceTo(pFirst, pLast);
    }

private:
    const VectorType& oU;
};

class Vec3DirectExpr {
public:
	static const bool bAlias = false;
	static const index_type iNumRows = 3;
    typedef ScalarTypeTraits<doublereal>::ScalarType ScalarType;
    typedef ScalarTypeTraits<doublereal>::DirectExpressionType ExpressionType;

    Vec3DirectExpr(const Vec3& u)
    	:oU(u) {

    }

    const ScalarType& operator()(index_type i) const {
    	MATVEC_ASSERT(i >= 1);
    	MATVEC_ASSERT(i <= iNumRows);
    	return oU(i);
    }

    index_type iGetNumRows() const {
    	return iNumRows;
    }

    template <typename ScalarType2>
    bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
    	return bArrayOverlap(oU.pGetVec(), oU.pGetVec() + 2, pFirst, pLast);
    }

private:
    const Vec3& oU;
};

template <typename T, index_type N_rows, index_type N_offset>
class SliceVector {
public:
	static const bool bAlias = MatVecHelp::AliasTypeHelper<T>::bAlias;
	static const index_type iNumRows = N_rows;
    typedef typename ScalarTypeTraits<T>::ScalarType ScalarType;
    typedef typename ScalarTypeTraits<T>::DirectExpressionType ExpressionType;

    SliceVector(const ScalarType* p)
    	:pVec(p) {

    }

    const ScalarType& operator()(index_type iRow) const {
    	--iRow;	// row index is 1-based
    	MATVEC_ASSERT(iRow >= 0);
    	MATVEC_ASSERT(iRow < iGetNumRows());
        return *(pVec + iRow * N_offset);
    }

    static index_type iGetNumRows() { return iNumRows; }

    template <typename ScalarType2>
    bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
    	return bArrayOverlap(pVec, pVec + iNumRows - 1, pFirst, pLast);
    }

private:
    const ScalarType* const pVec;
};

template <typename MatrixExpr>
class ColumnVectorExpr {
public:
	static const bool bAlias = MatrixExpr::bAlias;
	static const index_type iNumRows = MatrixExpr::iNumRows;
	typedef typename MatrixExpr::ScalarType ScalarType;
	typedef typename MatrixExpr::ExpressionType ExpressionType;

	ColumnVectorExpr(const MatrixExpr& A, index_type iCol)
		:A(A), iCol(iCol) {
		MATVEC_ASSERT(iNumRows == A.iGetNumRows());
	}

    ExpressionType operator()(index_type iRow) const {
    	MATVEC_ASSERT(iRow >= 1);
    	MATVEC_ASSERT(iRow <= iNumRows);
        return A(iRow, iCol);
    }

    static index_type iGetNumRows() { return iNumRows; }

    template <typename ScalarType2>
    bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
    	return A.bHaveReferenceTo(pFirst, pLast);
    }

private:
	const MatrixExpr A;
	const index_type iCol;
};

template <typename MatrixExpr>
class RowVectorExpr {
public:
	static const bool bAlias = MatrixExpr::bAlias;
	static const index_type iNumRows = MatrixExpr::iNumCols;
	typedef typename MatrixExpr::ScalarType ScalarType;
	typedef typename MatrixExpr::ExpressionType ExpressionType;

	RowVectorExpr(const MatrixExpr& A, index_type iRow)
		:A(A), iRow(iRow) {
		MATVEC_ASSERT(iNumRows == A.iGetNumCols());
	}

    ExpressionType operator()(index_type iCol) const {
    	MATVEC_ASSERT(iCol >= 1);
    	MATVEC_ASSERT(iCol <= iNumRows);
        return A(iRow, iCol);
    }

    static index_type iGetNumRows() { return iNumRows; }

    template <typename ScalarType2>
    bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
    	return A.bHaveReferenceTo(pFirst, pLast);
    }

private:
	const MatrixExpr A;
	const index_type iRow;
};

template <typename MatrixExpr>
class TransposedMatrix {
public:
	static const bool bAlias = MatrixExpr::bAlias;
    static const index_type iNumRows = MatrixExpr::iNumCols;
    static const index_type iNumCols = MatrixExpr::iNumRows;
	typedef typename MatrixExpr::RowVectorType ColumnVectorType;
	typedef typename MatrixExpr::ColumnVectorType RowVectorType;
    typedef typename MatrixExpr::ScalarType ScalarType;
    typedef typename MatrixExpr::ExpressionType ExpressionType;

	TransposedMatrix(const MatrixExpr& A)
		:A(A) {
		MATVEC_ASSERT(iNumRows == A.iGetNumCols());
		MATVEC_ASSERT(iNumCols == A.iGetNumRows());
	}

    const ScalarType& operator()(index_type i, index_type j) const {
    	MATVEC_ASSERT(i >= 1);
    	MATVEC_ASSERT(i <= iGetNumRows());
    	MATVEC_ASSERT(j >= 1);
    	MATVEC_ASSERT(j <= iGetNumCols());
        return A(j, i);
    }

    index_type iGetNumRows() const {
    	MATVEC_ASSERT(iNumRows == A.iGetNumCols());
        return A.iGetNumCols();
    }

    index_type iGetNumCols() const {
    	MATVEC_ASSERT(iNumCols == A.iGetNumRows());
    	return A.iGetNumRows();
    }

    RowVectorType GetRow(index_type iRow) const {
    	MATVEC_ASSERT(iRow >= 1);
    	MATVEC_ASSERT(iRow <= iGetNumRows());
    	return A.GetCol(iRow);
    }

    ColumnVectorType GetCol(index_type iCol) const {
    	MATVEC_ASSERT(iCol >= 1);
    	MATVEC_ASSERT(iCol <= iGetNumCols());
    	return A.GetRow(iCol);
    }

    template <typename ScalarType2>
    bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
    	return A.bHaveReferenceTo(pFirst, pLast);
    }

private:
	const MatrixExpr A;
};

template <index_type iRowStart,
		  index_type iRowEnd,
		  index_type iColStart,
		  index_type iColEnd,
		  typename MatrixExpr>
class SubMatrixExpr {
public:
	static const bool bAlias = MatrixExpr::bAlias;
    static const index_type iNumRows = iRowEnd - iRowStart + 1;
    static const index_type iNumCols = iColEnd - iColStart + 1;

	typedef VectorExpression<SubVectorExpr<iColStart, iColEnd, typename MatrixExpr::RowVectorType>, iNumCols> RowVectorType;
	typedef VectorExpression<SubVectorExpr<iRowStart, iRowEnd, typename MatrixExpr::ColumnVectorType>, iNumRows> ColumnVectorType;
    typedef typename MatrixExpr::ScalarType ScalarType;
    typedef typename MatrixExpr::ExpressionType ExpressionType;

    SubMatrixExpr(const MatrixExpr& A)
		:A(A) {
		MATVEC_ASSERT(MatrixExpr::iNumCols == A.iGetNumCols());
		MATVEC_ASSERT(MatrixExpr::iNumRows == A.iGetNumRows());
	}

    const ScalarType& operator()(index_type i, index_type j) const {
    	MATVEC_ASSERT(i >= 1);
    	MATVEC_ASSERT(i <= iGetNumRows());
    	MATVEC_ASSERT(j >= 1);
    	MATVEC_ASSERT(j <= iGetNumCols());
        return A(i + iRowStart - 1, j + iColStart - 1);
    }

    index_type iGetNumRows() const {
        return iNumRows;
    }

    index_type iGetNumCols() const {
    	return iNumCols;
    }

    RowVectorType GetRow(index_type iRow) const {
    	MATVEC_ASSERT(iRow >= 1);
    	MATVEC_ASSERT(iRow <= iGetNumRows());
    	return RowVectorType(A.GetRow(iRow + iRowStart - 1));
    }

    ColumnVectorType GetCol(index_type iCol) const {
    	MATVEC_ASSERT(iCol >= 1);
    	MATVEC_ASSERT(iCol <= iGetNumCols());
    	return ColumnVectorType(A.GetCol(iCol + iColStart - 1));
    }

    template <typename ScalarType2>
    bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
    	return A.bHaveReferenceTo(pFirst, pLast);
    }

private:
    typedef typename MaxSizeCheck<iRowStart >= 1>::CheckType check_iRowStart;
    typedef typename MaxSizeCheck<iRowEnd <= MatrixExpr::iNumRows>::CheckType check_iRowEnd;
    typedef typename MaxSizeCheck<iColStart >= 1>::CheckType check_iColStart;
    typedef typename MaxSizeCheck<iColEnd <= MatrixExpr::iNumCols>::CheckType check_iColEnd;
	const MatrixExpr A;
};

template <typename ScalarBinFunc, typename MatrixLhsExpr, typename MatrixRhsExpr>
class MatrixMatrixMatrixBinaryExpr {
public:
	static const bool bAlias = MatrixLhsExpr::bAlias || MatrixRhsExpr::bAlias;
	static const index_type iNumRows = MatrixLhsExpr::iNumRows;
	static const index_type iNumCols = MatrixLhsExpr::iNumCols;
	typedef typename ScalarBinFunc::ScalarType ScalarType;
    typedef typename ScalarBinFunc::ExpressionType ExpressionType;
    typedef VectorExpression<RowVectorExpr<MatrixMatrixMatrixBinaryExpr>, iNumCols> RowVectorType;
    typedef VectorExpression<ColumnVectorExpr<MatrixMatrixMatrixBinaryExpr>, iNumRows> ColumnVectorType;

    MatrixMatrixMatrixBinaryExpr(const MatrixLhsExpr& u, const MatrixRhsExpr& v)
        :oU(u), oV(v) {

    	MATVEC_ASSERT(iNumRows == oU.iGetNumRows());
    	MATVEC_ASSERT(iNumCols == oU.iGetNumCols());
    	MATVEC_ASSERT(iNumRows == oV.iGetNumRows());
    	MATVEC_ASSERT(iNumCols == oV.iGetNumCols());
    }

    ExpressionType operator()(index_type i, index_type j) const {
    	MATVEC_ASSERT(i >= 1);
    	MATVEC_ASSERT(i <= iGetNumRows());
    	MATVEC_ASSERT(j >= 1);
    	MATVEC_ASSERT(j <= iGetNumCols());
        return ScalarBinFunc::f(oU(i, j), oV(i, j));
    }

    index_type iGetNumRows() const {
        MATVEC_ASSERT(oU.iGetNumRows() == oV.iGetNumRows());
        MATVEC_ASSERT(iNumRows == oU.iGetNumRows());
        return oU.iGetNumRows();
    }

    index_type iGetNumCols() const {
    	MATVEC_ASSERT(oU.iGetNumCols() == oV.iGetNumCols());
    	MATVEC_ASSERT(iNumCols == oU.iGetNumCols());
    	return oU.iGetNumCols();
    }

    RowVectorType GetRow(index_type iRow) const {
    	return RowVectorType(*this, iRow);
    }

    ColumnVectorType GetCol(index_type iCol) const {
    	return ColumnVectorType(*this, iCol);
    }

    template <typename ScalarType2>
    bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
    	return oU.bHaveReferenceTo(pFirst, pLast) || oV.bHaveReferenceTo(pFirst, pLast);
    }

private:
    const MatrixLhsExpr oU;
    const MatrixRhsExpr oV;

    // check if the dimensions of both matrices are the same
    typedef typename IndexCheck<MatrixLhsExpr::iNumRows - MatrixRhsExpr::iNumRows>::CheckType check_iNumRows;
    typedef typename IndexCheck<MatrixLhsExpr::iNumCols - MatrixRhsExpr::iNumCols>::CheckType check_iNumCols;
};

/**
 * This class handles expressions of the form
 * f(matrix1, scalar1) = matrix2
 */
template <typename ScalarBinFunc, typename MatrixLhsExpr, typename ScalarRhsExpr>
class MatrixScalarMatrixBinaryExpr {
public:
	static const bool bAlias = MatrixLhsExpr::bAlias || MatVecHelp::AliasTypeHelper<ScalarRhsExpr>::bAlias;
	static const index_type iNumRows = MatrixLhsExpr::iNumRows;
	static const index_type iNumCols = MatrixLhsExpr::iNumCols;
	typedef typename ScalarBinFunc::ScalarType ScalarType;
    typedef typename ScalarBinFunc::ExpressionType ExpressionType;
    typedef VectorExpression<RowVectorExpr<MatrixScalarMatrixBinaryExpr>, iNumCols> RowVectorType;
    typedef VectorExpression<ColumnVectorExpr<MatrixScalarMatrixBinaryExpr>, iNumRows> ColumnVectorType;

    MatrixScalarMatrixBinaryExpr(const MatrixLhsExpr& u, const ScalarRhsExpr& v)
        :oU(u), oV(v) {

    	MATVEC_ASSERT(iNumRows == oU.iGetNumRows());
    	MATVEC_ASSERT(iNumCols == oU.iGetNumCols());
    }

    ExpressionType operator()(index_type i, index_type j) const {
    	MATVEC_ASSERT(i >= 1);
    	MATVEC_ASSERT(i <= iGetNumRows());
    	MATVEC_ASSERT(j >= 1);
    	MATVEC_ASSERT(j <= iGetNumCols());
        return ScalarBinFunc::f(oU(i, j), oV);
    }

    index_type iGetNumRows() const {
    	MATVEC_ASSERT(iNumRows == oU.iGetNumRows());
        return oU.iGetNumRows();
    }

    index_type iGetNumCols() const {
    	MATVEC_ASSERT(iNumCols == oU.iGetNumCols());
    	return oU.iGetNumCols();
    }

    RowVectorType GetRow(index_type iRow) const {
    	MATVEC_ASSERT(iRow >= 1);
    	MATVEC_ASSERT(iRow <= iGetNumRows());
    	return RowVectorType(*this, iRow);
    }

    ColumnVectorType GetCol(index_type iCol) const {
    	MATVEC_ASSERT(iCol >= 1);
    	MATVEC_ASSERT(iCol <= iGetNumCols());
    	return ColumnVectorType(*this, iCol);
    }

    template <typename ScalarType2>
    bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
    	typedef ScalarTypeTraits<BasicScalarType<ScalarRhsExpr> > ScalarTraits;
    	return oU.bHaveReferenceTo(pFirst, pLast)
    			|| ScalarTraits::bHaveReferenceTo(oV, pFirst, pLast);
    }

private:
    const MatrixLhsExpr oU;
    const ScalarRhsExpr oV;
};

template <typename ScalarUnaryFunc, typename MatrixExpr>
class MatrixMatrixUnaryExpr {
public:
	static const bool bAlias = MatrixExpr::bAlias;
	static const index_type iNumRows = MatrixExpr::iNumRows;
	static const index_type iNumCols = MatrixExpr::iNumCols;
	typedef typename ScalarUnaryFunc::ScalarType ScalarType;
    typedef typename ScalarUnaryFunc::ExpressionType ExpressionType;
    typedef VectorExpression<RowVectorExpr<MatrixMatrixUnaryExpr>, iNumCols> RowVectorType;
    typedef VectorExpression<ColumnVectorExpr<MatrixMatrixUnaryExpr>, iNumRows> ColumnVectorType;

    MatrixMatrixUnaryExpr(const MatrixExpr& u)
        :oU(u) {
    	MATVEC_ASSERT(iNumRows == oU.iGetNumRows());
    	MATVEC_ASSERT(iNumCols == oU.iGetNumCols());
    }

    ExpressionType operator()(index_type i, index_type j) const {
    	MATVEC_ASSERT(i >= 1);
    	MATVEC_ASSERT(i <= iGetNumRows());
    	MATVEC_ASSERT(j >= 1);
    	MATVEC_ASSERT(j <= iGetNumCols());
        return ScalarUnaryFunc::f(oU(i, j));
    }

    index_type iGetNumRows() const {
    	MATVEC_ASSERT(iNumRows == oU.iGetNumRows());
        return oU.iGetNumRows();
    }

    index_type iGetNumCols() const {
    	MATVEC_ASSERT(iNumCols == oU.iGetNumCols());
    	return oU.iGetNumCols();
    }

    RowVectorType GetRow(index_type iRow) const {
    	MATVEC_ASSERT(iRow >= 1);
    	MATVEC_ASSERT(iRow <= iGetNumRows());
    	return RowVectorType(*this, iRow);
    }

    ColumnVectorType GetCol(index_type iCol) const {
    	MATVEC_ASSERT(iCol >= 1);
    	MATVEC_ASSERT(iCol <= iGetNumCols());
    	return ColumnVectorType(*this, iCol);
    }

    template <typename ScalarType2>
    bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
    	return oU.bHaveReferenceTo(pFirst, pLast);
    }

private:
    const MatrixExpr oU;
};

template <typename MatrixType, bool ALIAS=false>
class MatrixDirectExpr {
public:
	static const bool bAlias = ALIAS;
	static const index_type iNumRows = MatrixType::iNumRows;
	static const index_type iNumCols = MatrixType::iNumCols;
    typedef typename MatrixType::ScalarType ScalarType;
    typedef typename MatrixType::ExpressionType ExpressionType;
    typedef typename MatrixType::RowVectorType RowVectorType;
    typedef typename MatrixType::ColumnVectorType ColumnVectorType;

    MatrixDirectExpr(const MatrixType& A)
        :A(A) {
    	MATVEC_ASSERT(iNumRows == A.iGetNumRows());
    	MATVEC_ASSERT(iNumCols == A.iGetNumCols());
    }

    const ScalarType& operator()(index_type i, index_type j) const {
    	MATVEC_ASSERT(i >= 1);
    	MATVEC_ASSERT(i <= iGetNumRows());
    	MATVEC_ASSERT(j >= 1);
    	MATVEC_ASSERT(j <= iGetNumCols());
        return A(i, j);
    }

    index_type iGetNumRows() const {
    	MATVEC_ASSERT(iNumRows == A.iGetNumRows());
        return A.iGetNumRows();
    }

    index_type iGetNumCols() const {
    	MATVEC_ASSERT(iNumCols == A.iGetNumCols());
    	return A.iGetNumCols();
    }

    RowVectorType GetRow(index_type iRow) const {
    	MATVEC_ASSERT(iRow >= 1);
    	MATVEC_ASSERT(iRow <= iGetNumRows());
    	return A.GetRow(iRow);
    }

    ColumnVectorType GetCol(index_type iCol) const {
    	MATVEC_ASSERT(iCol >= 1);
    	MATVEC_ASSERT(iCol <= iGetNumCols());
    	return A.GetCol(iCol);
    }

    template <typename ScalarType2>
    bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
    	return A.bHaveReferenceTo(pFirst, pLast);
    }

private:
    const MatrixType& A;
};

class Mat3x3DirectExpr {
public:
	static const bool bAlias = false;
	static const index_type iNumRows = 3;
	static const index_type iNumCols = 3;

    typedef ScalarTypeTraits<doublereal>::ScalarType ScalarType;
    typedef ScalarTypeTraits<doublereal>::DirectExpressionType ExpressionType;

    typedef VectorExpression<SliceVector<doublereal, iNumCols, iNumRows>, iNumCols> RowVectorType;
    typedef VectorExpression<SliceVector<doublereal, iNumRows, 1>, iNumRows> ColumnVectorType;

    Mat3x3DirectExpr(const Mat3x3& A)
        :A(A) {

    }

    const ScalarType& operator()(index_type i, index_type j) const {
    	MATVEC_ASSERT(i >= 1);
    	MATVEC_ASSERT(i <= iGetNumRows());
    	MATVEC_ASSERT(j >= 1);
    	MATVEC_ASSERT(j <= iGetNumCols());
        return A(i, j);
    }

    index_type iGetNumRows() const {
        return iNumRows;
    }

    index_type iGetNumCols() const {
    	return iNumCols;
    }

    RowVectorType GetRow(index_type iRow) const {
    	MATVEC_ASSERT(iRow >= 1);
    	MATVEC_ASSERT(iRow <= iGetNumRows());
    	--iRow;
    	return RowVectorType(A.pGetMat() + iRow);
    }

    ColumnVectorType GetCol(index_type iCol) const {
    	MATVEC_ASSERT(iCol >= 1);
    	MATVEC_ASSERT(iCol <= iGetNumCols());
    	--iCol;
    	return ColumnVectorType(A.pGetMat() + iCol * iNumRows);
    }

    template <typename ScalarType2>
    bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
    	return bArrayOverlap(A.pGetMat(), A.pGetMat() + (iNumRows - 1) * (iNumCols - 1), pFirst, pLast);
    }

private:
    const Mat3x3& A;
};

template <typename InitClass, typename T, index_type N_rows, index_type N_cols>
class MatrixInit: public InitClass {
public:
	template <typename InitArg>
	explicit MatrixInit(const InitArg& v)
		:InitClass(v) {

	}
};

template <typename T, index_type N_rows, index_type N_cols>
class Matrix {    
public:
    static const index_type iNumRows = N_rows;
    static const index_type iNumCols = N_cols;

    typedef typename ScalarTypeTraits<T>::ScalarType ScalarType;
    typedef typename ScalarTypeTraits<T>::DirectExpressionType ExpressionType;

    typedef VectorExpression<SliceVector<T, N_cols, 1>, N_cols> RowVectorType;
    typedef VectorExpression<SliceVector<T, N_rows, N_cols>, N_rows> ColumnVectorType;

    Matrix() {
    	for (index_type i = 0; i < iNumRows; ++i) {
    		ZeroInit(&rgMat[i][0], &rgMat[i][iNumCols - 1] + 1);
    	}
    }

    Matrix(const T& A11, const T& A21, const T& A12, const T& A22) {
    	typedef typename IndexCheck<iNumRows - 2>::CheckType check_iNumRows;
    	typedef typename IndexCheck<iNumCols - 2>::CheckType check_iNumCols;

    	(*this)(1, 1) = A11;
    	(*this)(2, 1) = A21;
    	(*this)(1, 2) = A12;
    	(*this)(2, 2) = A22;
    }

    explicit inline Matrix(const Mat3x3& A);

    template <typename InitClass>
    explicit Matrix(const MatrixInit<InitClass, T, N_rows, N_cols>& func) {
    	func.Initialize(*this);
    }

    template <typename Expression>
    Matrix(const MatrixExpression<Expression, N_rows, N_cols>& A) {
    	// No aliases are possible because the object did not exist before
    	using namespace MatVecHelp;
    	ApplyMatrixFuncNoAlias(A, Assign());
    }

    template <typename T2>
    Matrix(const Matrix<T2, N_rows, N_cols>& A, LocalDofMap* pDofMap) {
    	Copy(A, pDofMap);
    }

    template <typename T2>
    void Copy(const Matrix<T2, N_rows, N_cols>& A, LocalDofMap* pDofMap) {
    	MATVEC_ASSERT(A.iGetNumRows() == iGetNumRows());
    	MATVEC_ASSERT(A.iGetNumCols() == iGetNumCols());

    	for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
    		for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
    			grad::Copy((*this)(i, j), A(i, j), pDofMap);
    		}
    	}
    }

    void Reset() {
       	for (index_type i = 1; i <= iGetNumRows(); ++i) {
			for (index_type j = 1; j <= iGetNumCols(); ++j) {
				grad::Reset((*this)(i, j));
			}
		}
    }

    template <typename Expression>
    Matrix& operator=(const MatrixExpression<Expression, N_rows, N_cols>& A) {
    	using namespace MatVecHelp;

    	ApplyMatrixFunc<Assign>(A);

    	return *this;
    }

    template <typename T_Rhs>
    Matrix& operator+=(const Matrix<T_Rhs, N_rows, N_cols>& A) {
		using namespace MatVecHelp;

		ApplyMatrixFunc<Add>(MatrixExpression<MatrixDirectExpr<Matrix<T_Rhs, N_rows, N_cols> >, N_rows, N_cols>(A));

    	return *this;
    }

    template <typename T_Rhs>
    Matrix& operator-=(const Matrix<T_Rhs, N_rows, N_cols>& A) {
		using namespace MatVecHelp;

		ApplyMatrixFunc<Sub>(MatrixExpression<MatrixDirectExpr<Matrix<T_Rhs, N_rows, N_cols> >, N_rows, N_cols>(A));

    	return *this;
    }

    template <typename Expression>
    Matrix& operator+=(const MatrixExpression<Expression, N_rows, N_cols>& A) {
		using namespace MatVecHelp;

		ApplyMatrixFunc<Add>(A);

    	return *this;
    }

    template <typename Expression>
    Matrix& operator-=(const MatrixExpression<Expression, N_rows, N_cols>& A) {
		using namespace MatVecHelp;

		ApplyMatrixFunc<Sub>(A);

    	return *this;
    }

    template <typename T_Rhs>
    Matrix& operator*=(const T_Rhs& a) {
    	using namespace MatVecHelp;

    	ApplyScalarFunc<Mul>(a);

    	return *this;
    }

    template <typename T_Rhs>
    Matrix& operator/=(const T_Rhs& a) {
    	using namespace MatVecHelp;

    	ApplyScalarFunc<Div>(a);

    	return *this;
    }

    template <typename ScalarExpression>
    Matrix& operator*=(const GradientExpression<ScalarExpression>& a) {
    	using namespace MatVecHelp;

    	ApplyScalarFunc<Mul>(ScalarType(a));

    	return *this;
    }

    template <typename ScalarExpression>
    Matrix& operator/=(const GradientExpression<ScalarExpression>& a) {
    	using namespace MatVecHelp;

    	ApplyScalarFunc<Div>(ScalarType(a));

    	return *this;
    }

    inline Matrix& operator=(const Mat3x3& A);

    const ScalarType& operator()(index_type iRow, index_type iCol) const {
    	--iRow;	// row and column indices are 1-based for compatibility reasons
    	--iCol;
        MATVEC_ASSERT(iRow >= 0);
        MATVEC_ASSERT(iRow < iNumRows);
        MATVEC_ASSERT(iCol >= 0);
        MATVEC_ASSERT(iCol < iNumCols);
        return rgMat[iRow][iCol];
    }
    
    ScalarType& operator()(index_type iRow, index_type iCol) {
    	--iRow;	// row and column indices are 1-based for compatibility reasons
    	--iCol;
        MATVEC_ASSERT(iRow >= 0);
        MATVEC_ASSERT(iRow < iNumRows);
        MATVEC_ASSERT(iCol >= 0);
        MATVEC_ASSERT(iCol < iNumCols);
        return rgMat[iRow][iCol];
    }    
    
    RowVectorType
    GetRow(index_type iRow) const {
    	MATVEC_ASSERT(iRow >= 1);
    	MATVEC_ASSERT(iRow <= iNumRows);
    	--iRow;
    	return RowVectorType(&rgMat[iRow][0]);
    }

    ColumnVectorType
    GetCol(index_type iCol) const {
    	MATVEC_ASSERT(iCol >= 1);
    	MATVEC_ASSERT(iCol <= iNumCols);
    	--iCol;
    	return ColumnVectorType(&rgMat[0][iCol]);
    }

    index_type iGetNumRows() const { return iNumRows; }
    index_type iGetNumCols() const { return iNumCols; }
    const ScalarType* pGetMat() const { return &rgMat[0][0]; }

    template <typename ScalarType2>
    bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
    	return bArrayOverlap(&rgMat[0][0],
    						 &rgMat[iNumRows - 1][iNumCols -1],
    						 pFirst,
    						 pLast);
    }

private:
    const ScalarType* pGetFirstElem() const {
    	return &rgMat[0][0];
    }

    const ScalarType* pGetLastElem() const {
    	return &rgMat[iNumRows - 1][iNumCols - 1];
    }

    friend struct MatVecHelp::ApplyAliasHelperMatrix<false>;
    friend struct MatVecHelp::ApplyAliasHelperMatrix<true>;

    template <typename Func, typename U>
    void ApplyScalarFunc(const U& a) {
		for (index_type i = 1; i <= iGetNumRows(); ++i) {
			for (index_type j = 1; j <= iGetNumCols(); ++j) {
				Func::Eval((*this)(i, j), a);
			}
		}
    }

    template <typename Func, typename Expression>
    void ApplyMatrixFunc(const MatrixExpression<Expression, N_rows, N_cols>& A) {
    	using namespace MatVecHelp;

    	ApplyAliasHelperMatrix<MatrixExpression<Expression, N_rows, N_cols>::bAlias>::ApplyMatrixFunc(*this, A, Func());
    }

    template <typename Func, typename Expression>
    void ApplyMatrixFuncNoAlias(const MatrixExpression<Expression, N_rows, N_cols>& A, const Func&) {
			MATVEC_ASSERT(N_rows == A.iGetNumRows());
			MATVEC_ASSERT(N_cols == A.iGetNumCols());
		MATVEC_ASSERT(A.iGetNumRows() == iGetNumRows());
		MATVEC_ASSERT(A.iGetNumCols() == iGetNumCols());
    	MATVEC_ASSERT(!A.bHaveReferenceTo(pGetFirstElem(), pGetLastElem()));

		for (index_type i = 1; i <= iGetNumRows(); ++i) {
			for (index_type j = 1; j <= iGetNumCols(); ++j) {
				Func::Eval((*this)(i, j), A(i, j));
			}
    	}
    }

    template <typename Func, typename Expression>
    void ApplyMatrixFuncAlias(const MatrixExpression<Expression, N_rows, N_cols>& A, const Func& f) {
    	ApplyMatrixFuncNoAlias(MatrixExpression<MatrixDirectExpr<Matrix>, N_rows, N_cols>(Matrix(A)), f);
    }
private:
    ScalarType rgMat[iNumRows][iNumCols];
};

template <typename T, index_type N_rows, index_type N_cols>
inline MatrixExpression<TransposedMatrix<MatrixDirectExpr<Matrix<T, N_rows, N_cols> > >, N_cols, N_rows>
Transpose(const Matrix<T, N_rows, N_cols>& A) {
	return MatrixExpression<TransposedMatrix<MatrixDirectExpr<Matrix<T, N_rows, N_cols> > >, N_cols, N_rows>(Direct(A));
}

template <typename MatrixExpr, index_type N_rows, index_type N_cols>
inline MatrixExpression<TransposedMatrix<MatrixExpr>, N_cols, N_rows>
Transpose(const MatrixExpression<MatrixExpr, N_rows, N_cols>& A) {
	return MatrixExpression<TransposedMatrix<MatrixExpr>, N_cols, N_rows>(A);
}

template <typename T, index_type N_rows, index_type N_cols>
inline MatrixExpression<MatrixDirectExpr<Matrix<T, N_rows, N_cols> >, N_rows, N_cols>
Direct(const Matrix<T, N_rows, N_cols>& A) {
	return MatrixExpression<MatrixDirectExpr<Matrix<T, N_rows, N_cols> >, N_rows, N_cols>(A);
}

inline MatrixExpression<Mat3x3DirectExpr, 3, 3>
Direct(const Mat3x3& A) {
	return MatrixExpression<Mat3x3DirectExpr, 3, 3>(A);
}

template <typename T, index_type N_rows, index_type N_cols>
inline MatrixExpression<MatrixDirectExpr<Matrix<T, N_rows, N_cols>, true>, N_rows, N_cols>
Alias(const Matrix<T, N_rows, N_cols>& A) {
	return MatrixExpression<MatrixDirectExpr<Matrix<T, N_rows, N_cols>, true>, N_rows, N_cols>(A);
}

template <typename T, index_type N_rows, index_type N_cols>
inline Matrix<T, N_rows, N_cols>::Matrix(const Mat3x3& A) {
	using namespace MatVecHelp;

	ApplyMatrixFunc<Assign>(Direct(A));
}

template <typename T, index_type N_rows, index_type N_cols>
inline Matrix<T, N_rows, N_cols>&
Matrix<T, N_rows, N_cols>::operator=(const Mat3x3& A) {
	using namespace MatVecHelp;

	ApplyMatrixFunc<Assign>(Direct(A));

	return *this;
}

template <typename T, index_type N_rows>
class Vector {
public:
    static const index_type iNumRows = N_rows;
    typedef typename ScalarTypeTraits<T>::ScalarType ScalarType;
    typedef typename ScalarTypeTraits<T>::DirectExpressionType ExpressionType;

    Vector() {
    	ZeroInit(rgVec, rgVec + iGetNumRows());
    }

    Vector(const T& v1, const T& v2) {
    	typedef typename IndexCheck<iNumRows - 2>::CheckType check_iNumRows;

    	(*this)(1) = v1;
    	(*this)(2) = v2;
    }

    template <typename Expr1, typename Expr2>
    Vector(const GradientExpression<Expr1>& v1, const GradientExpression<Expr2>& v2)
    {
    	typedef typename IndexCheck<iNumRows - 2>::CheckType check_iNumRows;

    	(*this)(1) = v1;
    	(*this)(2) = v2;
    }

    Vector(const T& v1, const T& v2, const T& v3) {
    	typedef typename IndexCheck<iNumRows - 3>::CheckType check_iNumRows;

    	(*this)(1) = v1;
    	(*this)(2) = v2;
    	(*this)(3) = v3;
    }

    template <typename Expr1, typename Expr2, typename Expr3>
    Vector(const GradientExpression<Expr1>& v1, const GradientExpression<Expr2>& v2, const GradientExpression<Expr3>& v3)
    {
    	typedef typename IndexCheck<iNumRows - 3>::CheckType check_iNumRows;

    	(*this)(1) = v1;
    	(*this)(2) = v2;
    	(*this)(3) = v3;
    }

    explicit inline Vector(const Vec3& v);

    template <typename Expression>
    Vector(const VectorExpression<Expression, N_rows>& f) {
    	using namespace MatVecHelp;

    	ApplyMatrixFuncNoAlias(f, Assign());
    }
    
    template <typename T2>
    explicit Vector(const Vector<T2, N_rows>& v) {
    	MATVEC_ASSERT(v.iGetNumRows() == iGetNumRows());

    	for (index_type i = 1; i <= v.iGetNumRows(); ++i) {
    		grad::Convert((*this)(i), v(i));
    	}
    }

    template <typename T2>
    Vector(const Vector<T2, N_rows>& v, LocalDofMap* pDofMap) {
    	Copy(v, pDofMap);
    }

    template <typename T2>
    void Copy(const Vector<T2, N_rows>& v, LocalDofMap* pDofMap) {
    	MATVEC_ASSERT(v.iGetNumRows() == iGetNumRows());

    	for (index_type i = 1; i <= iGetNumRows(); ++i) {
    		grad::Copy((*this)(i), v(i), pDofMap);
    	}
    }

    void Reset() {
    	for (index_type i = 1; i <= iGetNumRows(); ++i) {
    		grad::Reset((*this)(i));
    	}
    }

    template <typename Expression>
    Vector& operator=(const VectorExpression<Expression, N_rows>& f) {
    	using namespace MatVecHelp;

    	ApplyMatrixFunc<Assign>(f);

        return *this;
    }

    template <typename T_Rhs>
    Vector& operator+=(const Vector<T_Rhs, N_rows>& v) {
    	using namespace MatVecHelp;

    	ApplyMatrixFunc<Add>(VectorExpression<VectorDirectExpr<Vector>, N_rows>(v));

    	return *this;
    }

    template <typename T_Rhs>
    Vector& operator-=(const Vector<T_Rhs, N_rows>& v) {
    	using namespace MatVecHelp;

    	ApplyMatrixFunc<Sub>(VectorExpression<VectorDirectExpr<Vector>, N_rows>(v));

    	return *this;
    }

    template <typename Expression>
    Vector& operator+=(const VectorExpression<Expression, N_rows>& f) {
    	using namespace MatVecHelp;

    	ApplyMatrixFunc<Add>(f);

    	return *this;
    }

    template <typename Expression>
    Vector& operator-=(const VectorExpression<Expression, N_rows>& f) {
    	using namespace MatVecHelp;

    	ApplyMatrixFunc<Sub>(f);

    	return *this;
    }

    template <typename T_Rhs>
    Vector& operator*=(const T_Rhs& g) {
    	using namespace MatVecHelp;

    	ApplyScalarFunc<Mul>(g);

    	return *this;
    }

    template <typename T_Rhs>
    Vector& operator/=(const T_Rhs& g) {
    	using namespace MatVecHelp;

    	ApplyScalarFunc<Div>(g);

    	return *this;
    }

    template <typename ScalarExpression>
    Vector& operator*=(const GradientExpression<ScalarExpression>& f) {
    	using namespace MatVecHelp;

    	ApplyScalarFunc<Mul>(ScalarType(f));

    	return *this;
    }

    template <typename ScalarExpression>
    Vector& operator/=(const GradientExpression<ScalarExpression>& f) {
    	using namespace MatVecHelp;

    	ApplyScalarFunc<Div>(ScalarType(f));

    	return *this;
    }

    inline Vector& operator=(const Vec3& v);

    const ScalarType& operator()(index_type iRow) const {
    	--iRow; // Row index is 1-based for compatibility reasons
        MATVEC_ASSERT(iRow >= 0);
        MATVEC_ASSERT(iRow < iNumRows);
        return rgVec[iRow];
    }
    
    ScalarType& operator()(index_type iRow) {
    	--iRow; // Row index is 1-based for compatibility reasons
        MATVEC_ASSERT(iRow >= 0);
        MATVEC_ASSERT(iRow < iNumRows);
        return rgVec[iRow];
    }
    
    index_type iGetNumRows() const { return iNumRows; }
    
    ScalarType* pGetVec() { return rgVec; }
    const ScalarType* pGetVec() const { return rgVec; }

    template <typename ScalarType2>
    bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
    	return bArrayOverlap(&rgVec[0],
    						 &rgVec[iNumRows - 1],
    						 pFirst,
    						 pLast);
    }

private:
    const ScalarType* pGetFirstElem() const {
    	return &rgVec[0];
    }

    const ScalarType* pGetLastElem() const {
    	return &rgVec[iNumRows - 1];
    }

    friend struct MatVecHelp::ApplyAliasHelperMatrix<false>;
	friend struct MatVecHelp::ApplyAliasHelperMatrix<true>;

	template <typename Func, typename U>
	void ApplyScalarFunc(const U& a) {
		for (index_type i = 1; i <= iGetNumRows(); ++i) {
			Func::Eval((*this)(i), a);
			}
        }

	template <typename Func, typename Expression>
	void ApplyMatrixFunc(const VectorExpression<Expression, N_rows>& A) {
		using namespace MatVecHelp;

		ApplyAliasHelperMatrix<VectorExpression<Expression, N_rows>::bAlias>::ApplyMatrixFunc(*this, A, Func());
	}

	template <typename Func, typename Expression>
	void ApplyMatrixFuncNoAlias(const VectorExpression<Expression, N_rows>& A, const Func&) {
		MATVEC_ASSERT(N_rows == A.iGetNumRows());
		MATVEC_ASSERT(A.iGetNumRows() == iGetNumRows());
		MATVEC_ASSERT(!A.bHaveReferenceTo(pGetFirstElem(), pGetLastElem()));

		for (index_type i = 1; i <= iGetNumRows(); ++i) {
			Func::Eval((*this)(i), A(i));
		}
	}

	template <typename Func, typename Expression>
	void ApplyMatrixFuncAlias(const VectorExpression<Expression, N_rows>& A, const Func& f) {
		ApplyMatrixFuncNoAlias(VectorExpression<VectorDirectExpr<Vector>, N_rows>(Vector(A)), f);
    }

private:
    ScalarType rgVec[iNumRows];
};

template <typename T, index_type N_rows>
inline VectorExpression<VectorDirectExpr<Vector<T, N_rows> >, N_rows>
Direct(const Vector<T, N_rows>& v) {
	return VectorExpression<VectorDirectExpr<Vector<T, N_rows> >, N_rows>(v);
}

inline VectorExpression<Vec3DirectExpr, 3>
Direct(const Vec3& v) {
	return VectorExpression<Vec3DirectExpr, 3>(v);
}

template <typename T, index_type N_rows>
inline VectorExpression<VectorDirectExpr<Vector<T, N_rows>, true>, N_rows>
Alias(const Vector<T, N_rows>& v) {
	return VectorExpression<VectorDirectExpr<Vector<T, N_rows>, true>, N_rows>(v);
}

template <typename T, index_type N_rows>
inline
Vector<T, N_rows>::Vector(const Vec3& v) {
	using namespace MatVecHelp;

	ApplyMatrixFunc<Assign>(Direct(v));
}

template <typename T, index_type N_rows>
inline Vector<T, N_rows>&
Vector<T, N_rows>::operator=(const Vec3& v) {
	using namespace MatVecHelp;

	ApplyMatrixFunc<Assign>(Direct(v));

	return *this;
}

template <index_type iStartIndex, index_type iEndIndex, typename T, index_type N_rows>
VectorExpression<SubVectorExpr<iStartIndex, iEndIndex, VectorDirectExpr<Vector<T, N_rows> > >, iEndIndex - iStartIndex + 1>
SubVector(const Vector<T, N_rows>& v) {
	MATVEC_ASSERT(iStartIndex >= 1);
	MATVEC_ASSERT(iEndIndex <= v.iGetNumRows());
	return VectorExpression<SubVectorExpr<iStartIndex, iEndIndex, VectorDirectExpr<Vector<T, N_rows> > >, iEndIndex - iStartIndex + 1>(Direct(v));
}

template <index_type iStartIndex, index_type iEndIndex, typename VectorExpr, index_type N_rows>
VectorExpression<SubVectorExpr<iStartIndex, iEndIndex, VectorExpr> , iEndIndex - iStartIndex + 1>
SubVector(const VectorExpression<VectorExpr, N_rows>& v) {
	MATVEC_ASSERT(iStartIndex >= 1);
	MATVEC_ASSERT(iEndIndex <= v.iGetNumRows());
	return VectorExpression<SubVectorExpr<iStartIndex, iEndIndex, VectorExpr> , iEndIndex - iStartIndex + 1>(v);
}

template <index_type iRowStart, index_type iRowEnd, index_type iColStart, index_type iColEnd, typename T, index_type N_rows, index_type N_cols>
MatrixExpression<SubMatrixExpr<iRowStart, iRowEnd, iColStart, iColEnd, MatrixExpression<MatrixDirectExpr<Matrix<T, N_rows, N_cols> >, N_rows, N_cols> >, iRowEnd - iRowStart + 1, iColEnd - iColStart + 1>
SubMatrix(const Matrix<T, N_rows, N_cols>& A) {
	return MatrixExpression<SubMatrixExpr<iRowStart, iRowEnd, iColStart, iColEnd, MatrixExpression<MatrixDirectExpr<Matrix<T, N_rows, N_cols> >, N_rows, N_cols> >, iRowEnd - iRowStart + 1, iColEnd - iColStart + 1>(Direct(A));
}

template <index_type iRowStart, index_type iRowEnd, index_type iColStart, index_type iColEnd, typename MatrixExpr, index_type N_rows, index_type N_cols>
MatrixExpression<SubMatrixExpr<iRowStart, iRowEnd, iColStart, iColEnd, MatrixExpression<MatrixExpr, N_rows, N_cols> >, iRowEnd - iRowStart + 1, iColEnd - iColStart + 1>
SubMatrix(const MatrixExpression<MatrixExpr, N_rows, N_cols>& A) {
	return MatrixExpression<SubMatrixExpr<iRowStart, iRowEnd, iColStart, iColEnd, MatrixExpression<MatrixExpr, N_rows, N_cols> >, iRowEnd - iRowStart + 1, iColEnd - iColStart + 1>(A);
}

template <typename T, typename VectorExpr>
class MatCrossInit {
public:
	MatCrossInit(const VectorExpression<VectorExpr, 3>& v)
		:v(v) {

	}

	void Initialize(Matrix<T, 3, 3>& A) const {
		const index_type x = 1, y = 2, z = 3;
		/*
		 * A = [ 0, -z,  y;
		 *       z,  0, -x;
		 *      -y,  x,  0];
		 */
		A(x, x) = A(y, y) = A(z, z) = T(0.);
		A(x, y) = -v(z);
		A(x, z) =  v(y);
		A(y, x) =  v(z);
		A(y, z) = -v(x);
		A(z, x) = -v(y);
		A(z, y) =  v(x);
	}

private:
	const VectorExpression<VectorExpr, 3> v;
};

template <typename T, typename VectorExpr>
class MatCrossCrossInit {
public:
	MatCrossCrossInit(const VectorExpression<VectorExpr, 3>& v)
		:v(v) {

	}

	void Initialize(Matrix<T, 3, 3>& A) const {
		const index_type x = 1, y = 2, z = 3;

		const T vxvx = v(x) * v(x);
		const T vyvy = v(y) * v(y);
		const T vzvz = v(z) * v(z);
		const T vxvy = v(x) * v(y);
		const T vxvz = v(x) * v(z);
		const T vyvz = v(y) * v(z);

		A(x, x) = -vzvz - vyvy;
		A(x, y) = vxvy;
		A(x, z) = vxvz;
		A(y, x) = vxvy;
		A(y, y) = -vzvz - vxvx;
		A(y, z) = vyvz;
		A(z, x) = vxvz;
		A(z, y) = vyvz;
		A(z, z) = -vyvy - vxvx;
	}

private:
	const VectorExpression<VectorExpr, 3> v;
};

template <typename T, typename VectorExpr>
class MatGInit {
public:
	MatGInit(const VectorExpression<VectorExpr, 3>& g)
		:g(g) {

	}

	void Initialize(Matrix<T, 3, 3>& G) const {
/*
 * G = 4/(4 + g' * g) * (eye(3) + 1/2*skew(g))
 *
 * skew(g) = [     0, -g(3),   g(2);
 *              g(3),     0,  -g(1);
 *             -g(2),  g(1),     0]
 */
	    const T d = (4./(4.+Dot(g, g)));

		G(1, 1) = d;
		G(2, 1) = g(3) * d / 2.;
		G(3, 1) = -g(2) * d / 2.;
		G(1, 2) = -g(3) * d / 2.;
		G(2, 2) = d;
		G(3, 2) = g(1) * d / 2.;
		G(1, 3) = g(2) * d / 2.;
		G(2, 3) = -g(1) * d / 2.;
		G(3, 3) = d;
	}

private:
	const VectorExpression<VectorExpr, 3> g;
};

template <typename T>
inline MatrixInit<MatCrossInit<T, VectorDirectExpr<Vector<T, 3> > >, T, 3, 3>
MatCrossVec(const Vector<T, 3>& v) {
	return MatrixInit<MatCrossInit<T, VectorDirectExpr<Vector<T, 3> > >, T, 3, 3>(Direct(v));
}

template <typename VectorExpr>
inline MatrixInit<MatCrossInit<typename VectorExpr::ScalarType, VectorExpr>, typename VectorExpr::ScalarType, 3, 3>
MatCrossVec(const VectorExpression<VectorExpr, 3>& v) {
	return MatrixInit<MatCrossInit<typename VectorExpr::ScalarType, VectorExpr>, typename VectorExpr::ScalarType, 3, 3>(v);
}

template <typename T>
inline MatrixInit<MatCrossCrossInit<T, VectorDirectExpr<Vector<T, 3> > >, T, 3, 3>
MatCrossCrossVec(const Vector<T, 3>& v) {
	return MatrixInit<MatCrossCrossInit<T, VectorDirectExpr<Vector<T, 3> > >, T, 3, 3>(Direct(v));
}

template <typename VectorExpr>
inline MatrixInit<MatCrossCrossInit<typename VectorExpr::ScalarType, VectorExpr>, typename VectorExpr::ScalarType, 3, 3>
MatCrossCrossVec(const VectorExpression<VectorExpr, 3>& v) {
	return MatrixInit<MatCrossCrossInit<typename VectorExpr::ScalarType, VectorExpr>, typename VectorExpr::ScalarType, 3, 3>(v);
}

template <typename T>
inline MatrixInit<MatGInit<T, VectorDirectExpr<Vector<T, 3> > >, T, 3, 3>
MatGVec(const Vector<T, 3>& g) {
	return MatrixInit<MatGInit<T, VectorDirectExpr<Vector<T, 3> > >, T, 3, 3>(Direct(g));
}

template <typename VectorExpr>
inline MatrixInit<MatGInit<typename VectorExpr::ScalarType, VectorExpr>, typename VectorExpr::ScalarType, 3, 3>
MatGVec(const VectorExpression<VectorExpr, 3>& g) {
	return MatrixInit<MatGInit<typename VectorExpr::ScalarType, VectorExpr>, typename VectorExpr::ScalarType, 3, 3>(g);
}

template <typename T, index_type N_rows, index_type N_cols>
class TabularMatrixView {
public:
	TabularMatrixView(const Matrix<T, N_rows, N_cols>& A, int iColWidth)
		:A(A), iColWidth(iColWidth) {

	}

	void Print(std::ostream& os) const {
		for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
			for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
				os << std::setw(iColWidth) << A(i, j) << ' ';
			}
			os << std::endl;
		}
	}

private:
	const Matrix<T, N_rows, N_cols>& A;
	const int iColWidth;
};

template <typename T, index_type N_rows, index_type N_cols>
inline TabularMatrixView<T, N_rows, N_cols>
Tabular(const Matrix<T, N_rows, N_cols>& A, int iColWidth=10) {
	return TabularMatrixView<T, N_rows, N_cols>(A, iColWidth);
}

template <typename T, index_type N_rows, index_type N_cols>
inline std::ostream& operator<<(std::ostream& os, const Matrix<T, N_rows, N_cols>& A) {
    for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
        for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
            os << A(i, j);

            if (i < A.iGetNumRows() || j < A.iGetNumCols())
            	os << ' ';
        }
    }
    
    return os;
}

template <typename T, index_type N_rows, index_type N_cols>
inline std::ostream& operator<<(std::ostream& os, const TabularMatrixView<T, N_rows, N_cols>& tabA) {
	tabA.Print(os);

    return os;
}

template <typename T, index_type N_rows>
inline std::ostream& operator<<(std::ostream& os, const Vector<T, N_rows>& x) {
    for (index_type i = 1; i <= x.iGetNumRows(); ++i) {
        os << x(i);

        if (i < x.iGetNumRows())
        	os << ' ';
    }
    
    return os;
}

template <typename VectorExpressionType, index_type N_rows, index_type N_index>
struct SumTraits {
	typedef typename VectorExpression<VectorExpressionType, N_rows>::ScalarType ScalarType;
	typedef typename VectorExpression<VectorExpressionType, N_rows>::ExpressionType ScalarExpressionType;
	typedef SumTraits<VectorExpressionType, N_rows, N_index - 1> SumTraitsN_minus_1;
	typedef typename SumTraitsN_minus_1::ExpressionType ExpressionTypeN_minus_1;
	typedef typename ScalarBinaryExpressionTraits<FuncPlus, ScalarType, ExpressionTypeN_minus_1, ScalarExpressionType>::ExpressionType ExpressionType;

	static ExpressionType
	Sum(const VectorExpression<VectorExpressionType, N_rows>& v) {
		return ExpressionType(SumTraitsN_minus_1::Sum(v), v(N_index));
	}
};

template <typename VectorExpressionType, index_type N_rows>
struct SumTraits<VectorExpressionType, N_rows, 1> {
	typedef typename VectorExpression<VectorExpressionType, N_rows>::ScalarType ScalarType;
	typedef typename VectorExpression<VectorExpressionType, N_rows>::ExpressionType ExpressionType;

	static ExpressionType
	Sum(const VectorExpression<VectorExpressionType, N_rows>& v) {
		return v(1);
	}
};

template <typename VectorExpressionType, index_type N_rows>
inline typename SumTraits<VectorExpressionType, N_rows, N_rows>::ExpressionType
Sum(const VectorExpression<VectorExpressionType, N_rows>& v) {
	return SumTraits<VectorExpressionType, N_rows, N_rows>::Sum(v);
}

template <typename T, index_type N_rows>
inline typename SumTraits<VectorDirectExpr<Vector<T, N_rows> >, N_rows, N_rows>::ExpressionType
Sum(const Vector<T, N_rows>& v) {
	return SumTraits<VectorDirectExpr<Vector<T, N_rows> >, N_rows, N_rows>::Sum(Direct(v));
}

template <typename VectorExprLhs, typename VectorExprRhs, index_type N_rows, index_type N_index>
struct DotTraits {
	typedef typename VectorExpression<VectorExprLhs, N_rows>::ScalarType ScalarTypeLhs;
	typedef typename VectorExpression<VectorExprLhs, N_rows>::ExpressionType ScalarExprLhs;
	typedef typename VectorExpression<VectorExprRhs, N_rows>::ScalarType ScalarTypeRhs;
	typedef typename VectorExpression<VectorExprRhs, N_rows>::ExpressionType ScalarExprRhs;
	typedef typename CommonScalarType<ScalarTypeLhs, ScalarTypeRhs>::ScalarType ScalarType;

	typedef typename ScalarBinaryExpressionTraits<FuncMult,
												  ScalarType,
	  	  	  	  	  							  ScalarExprLhs,
	  	  	  	  	  							  ScalarExprRhs
	  	  	  	  	  							  >::ExpressionType MultExpressionType;

	typedef typename ScalarBinaryExpressionTraits<FuncPlus,
												  ScalarType,
												  typename DotTraits<VectorExprLhs, VectorExprRhs, N_rows, N_index - 1>::ExpressionType,
												  MultExpressionType
												  >::ExpressionType ExpressionType;

	static ExpressionType
	Dot(const VectorExpression<VectorExprLhs, N_rows>& u, const VectorExpression<VectorExprRhs, N_rows>& v) {
		return ExpressionType(DotTraits<VectorExprLhs, VectorExprRhs, N_rows, N_index - 1>::Dot(u, v), MultExpressionType(u(N_index), v(N_index)));
	}
};

template <typename VectorExprLhs, typename VectorExprRhs, index_type N_rows>
struct DotTraits<VectorExprLhs, VectorExprRhs, N_rows, 1> {
	typedef typename VectorExpression<VectorExprLhs, N_rows>::ScalarType ScalarTypeLhs;
	typedef typename VectorExpression<VectorExprLhs, N_rows>::ExpressionType ScalarExprLhs;
	typedef typename VectorExpression<VectorExprRhs, N_rows>::ScalarType ScalarTypeRhs;
	typedef typename VectorExpression<VectorExprRhs, N_rows>::ExpressionType ScalarExprRhs;
	typedef typename CommonScalarType<ScalarTypeLhs, ScalarTypeRhs>::ScalarType ScalarType;
	typedef typename ScalarBinaryExpressionTraits<FuncMult,
												  ScalarType,
	  	  	  	  	  							  ScalarExprLhs,
	  	  	  	  	  							  ScalarExprRhs
	  	  	  	  	  							  >::ExpressionType ExpressionType;

	static ExpressionType
	Dot(const VectorExpression<VectorExprLhs, N_rows>& u, const VectorExpression<VectorExprRhs, N_rows>& v) {
		return ExpressionType(u(1), v(1));
	}
};

template <typename VectorExprLhs, typename VectorExprRhs, index_type N_rows>
inline typename DotTraits<VectorExprLhs, VectorExprRhs, N_rows, N_rows>::ExpressionType
Dot(const VectorExpression<VectorExprLhs, N_rows>& u, const VectorExpression<VectorExprRhs, N_rows>& v) {
	return DotTraits<VectorExprLhs, VectorExprRhs, N_rows, N_rows>::Dot(u, v);
}

template <typename VectorExprLhs, typename T, index_type N_rows>
inline typename DotTraits<VectorExprLhs, VectorDirectExpr<Vector<T, N_rows> >, N_rows, N_rows>::ExpressionType
Dot(const VectorExpression<VectorExprLhs, N_rows>& u, const Vector<T, N_rows>& v) {
	return DotTraits<VectorExprLhs, VectorDirectExpr<Vector<T, N_rows> >, N_rows, N_rows>::Dot(u, Direct(v));
}

template <typename T, typename VectorExprRhs, index_type N_rows>
inline typename DotTraits<VectorDirectExpr<Vector<T, N_rows> >, VectorExprRhs, N_rows, N_rows>::ExpressionType
Dot(const Vector<T, N_rows>& u, const VectorExpression<VectorExprRhs, N_rows>& v) {
	return DotTraits<VectorDirectExpr<Vector<T, N_rows> >, VectorExprRhs, N_rows, N_rows>::Dot(Direct(u), v);
}

template <typename T_Lhs, typename T_Rhs, index_type N_rows>
inline typename DotTraits<VectorDirectExpr<Vector<T_Lhs, N_rows> >, VectorDirectExpr<Vector<T_Rhs, N_rows> >, N_rows, N_rows>::ExpressionType
Dot(const Vector<T_Lhs, N_rows>& u, const Vector<T_Rhs, N_rows>& v) {
	return DotTraits<VectorDirectExpr<Vector<T_Lhs, N_rows> >, VectorDirectExpr<Vector<T_Rhs, N_rows> >, N_rows, N_rows>::Dot(Direct(u), Direct(v));
}

template <typename VectorExpr, index_type N_rows>
inline typename VectorExpression<VectorExpr, N_rows>::ScalarType
Norm(const VectorExpression<VectorExpr, N_rows>& u) {
	// avoid double evaluation
	const Vector<typename VectorExpression<VectorExpr, N_rows>::ScalarType, N_rows> u1(u);
	using std::sqrt; // needed for g++-4.8
	return sqrt(Dot(u1, u1));
}

template <typename T, index_type N_rows>
inline T
Norm(const Vector<T, N_rows>& u) {
	return Norm(Direct(u));
}

template <typename VectorLhsExpr, typename VectorRhsExpr>
struct CrossTraits {
	typedef typename VectorLhsExpr::ScalarType ScalarTypeLhs;
	typedef typename VectorLhsExpr::ExpressionType ExpressionTypeLhs;
	typedef typename VectorRhsExpr::ScalarType ScalarTypeRhs;
	typedef typename VectorRhsExpr::ExpressionType ExpressionTypeRhs;
	typedef typename CommonScalarType<ScalarTypeLhs, ScalarTypeRhs>::ScalarType ScalarType;
	typedef typename ScalarBinaryExpressionTraits<FuncMult,
												  ScalarType,
												  ExpressionTypeLhs,
												  ExpressionTypeRhs
												  >::ExpressionType ExprMult;
	typedef typename ScalarBinaryExpressionTraits<FuncMinus,
												  ScalarType,
												  ExprMult,
												  ExprMult
												  >::ExpressionType ExpressionType;
private:
	typedef typename IndexCheck<VectorLhsExpr::iNumRows - 3>::CheckType check_iNumRowsLhs;
	typedef typename IndexCheck<VectorRhsExpr::iNumRows - 3>::CheckType check_iNumRowsRhs;
};

template <typename VectorLhsExpr, typename VectorRhsExpr>
class VectorCrossExpr {
public:
	static const bool bAlias = VectorLhsExpr::bAlias || VectorRhsExpr::bAlias;
	static const index_type iNumRows = 3;
	typedef typename CrossTraits<VectorLhsExpr, VectorRhsExpr>::ExprMult ExprMult;
	typedef typename CrossTraits<VectorLhsExpr, VectorRhsExpr>::ScalarType ScalarType;
	typedef typename CrossTraits<VectorLhsExpr, VectorRhsExpr>::ExpressionType ExpressionType;

	VectorCrossExpr(const VectorLhsExpr& u, const VectorRhsExpr& v)
		:oU(u), oV(v) {

		MATVEC_ASSERT(iNumRows == oU.iGetNumRows());
		MATVEC_ASSERT(iNumRows == oV.iGetNumRows());
	}

	ExpressionType operator()(index_type i) const {
		--i;
		MATVEC_ASSERT(i >= 0);
		MATVEC_ASSERT(i < iNumRows);

		static const index_type x = 1, y = 2, z = 3;
		static const index_type e1[3] = {y, z, x},
								e2[3] = {z, x, y},
								e3[3] = {z, x, y},
								e4[3] = {y, z, x};

		return ExpressionType(ExprMult(oU(e1[i]), oV(e2[i])),
							  ExprMult(oU(e3[i]), oV(e4[i])));
	}

	index_type iGetNumRows() const {
		return iNumRows;
	}

	template <typename ScalarType2>
	bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
		return oU.bHaveReferenceTo(pFirst, pLast) || oV.bHaveReferenceTo(pFirst, pLast);
	}

private:
	const VectorLhsExpr oU;
	const VectorRhsExpr oV;

	typedef typename IndexCheck<iNumRows - VectorLhsExpr::iNumRows>::CheckType check_VectorLhsExpr;
	typedef typename IndexCheck<iNumRows - VectorRhsExpr::iNumRows>::CheckType check_VectorRhsExpr;
};

template <typename VectorLhsExpr, typename VectorRhsExpr>
VectorExpression<VectorCrossExpr<VectorLhsExpr, VectorRhsExpr>, 3>
inline Cross(const VectorExpression<VectorLhsExpr, 3>& u, const VectorExpression<VectorRhsExpr, 3>& v) {
	return VectorExpression<VectorCrossExpr<VectorLhsExpr, VectorRhsExpr>, 3>(u, v);
}

template <typename VectorLhsExpr, typename T_Rhs>
VectorExpression<VectorCrossExpr<VectorLhsExpr, VectorDirectExpr<Vector<T_Rhs, 3> > >, 3>
inline Cross(const VectorExpression<VectorLhsExpr, 3>& u, const Vector<T_Rhs, 3>& v) {
	return VectorExpression<VectorCrossExpr<VectorLhsExpr, VectorDirectExpr<Vector<T_Rhs, 3> > >, 3>(u, Direct(v));
}

template <typename T_Lhs, typename VectorRhsExpr>
VectorExpression<VectorCrossExpr<VectorDirectExpr<Vector<T_Lhs, 3> >, VectorRhsExpr>, 3>
inline Cross(const Vector<T_Lhs, 3>& u, const VectorExpression<VectorRhsExpr, 3>& v) {
	return VectorExpression<VectorCrossExpr<VectorDirectExpr<Vector<T_Lhs, 3> >, VectorRhsExpr>, 3>(Direct(u), v);
}

template <typename T_Lhs, typename T_Rhs>
VectorExpression<VectorCrossExpr<VectorDirectExpr<Vector<T_Lhs, 3> >, VectorDirectExpr<Vector<T_Rhs, 3> > >, 3>
inline Cross(const Vector<T_Lhs, 3>& u, const Vector<T_Rhs, 3>& v) {
	return VectorExpression<VectorCrossExpr<VectorDirectExpr<Vector<T_Lhs, 3> >, VectorDirectExpr<Vector<T_Rhs, 3> > >, 3>(Direct(u), Direct(v));
}

template <typename T>
inline T Det(const Matrix<T, 2, 2>& A) {
	return A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1);
}

template <typename T>
inline Matrix<T, 2, 2> Inv(const Matrix<T, 2, 2>& A) {
	const T detA = Det(A);

	return Matrix<T, 2, 2>(A(2, 2) / detA,
						  -A(2, 1) / detA,
						  -A(1, 2) / detA,
						   A(1, 1) / detA);
}

template <typename MatrixLhsExpr, typename VectorRhsExpr>
class MatrixVectorProduct {
public:
	static const bool bAlias = MatrixLhsExpr::bAlias || VectorRhsExpr::bAlias;
	static const index_type iNumRows = MatrixLhsExpr::iNumRows;
	typedef typename MatrixLhsExpr::ScalarType MatrixLhsScalarExpr;
	typedef typename VectorRhsExpr::ScalarType VectorRhsScalarExpr;
	typedef typename MatrixLhsExpr::RowVectorType MatrixLhsRowVector;
	typedef typename CommonScalarType<typename BasicScalarType<MatrixLhsScalarExpr>::ScalarType,
									  typename BasicScalarType<VectorRhsScalarExpr>::ScalarType>::ScalarType ScalarType;
    typedef typename DotTraits<MatrixLhsRowVector, VectorRhsExpr, MatrixLhsExpr::iNumCols, MatrixLhsExpr::iNumCols>::ExpressionType ExpressionType;

    MatrixVectorProduct(const MatrixLhsExpr& A, const VectorRhsExpr& x)
        :A(A), x(x) {
    	MATVEC_ASSERT(iNumRows == A.iGetNumRows());
    }

    ExpressionType operator()(index_type i) const {
    	MATVEC_ASSERT(i >= 1);
    	MATVEC_ASSERT(i <= iGetNumRows());
        return Dot(A.GetRow(i), x);
    }

    index_type iGetNumRows() const {
    	MATVEC_ASSERT(iNumRows == A.iGetNumRows());
        return A.iGetNumRows();
    }

    template <typename ScalarType2>
    bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
    	return A.bHaveReferenceTo(pFirst, pLast) || x.bHaveReferenceTo(pFirst, pLast);
    }

private:
    const MatrixLhsExpr A;
    const VectorRhsExpr x;
};

template <typename MatrixLhsExpr, typename MatrixRhsExpr>
class MatrixMatrixProduct {
public:
	static const bool bAlias = MatrixLhsExpr::bAlias || MatrixRhsExpr::bAlias;
	static const index_type iNumRows = MatrixLhsExpr::iNumRows;
	static const index_type iNumCols = MatrixRhsExpr::iNumCols;
	typedef typename MatrixLhsExpr::ScalarType MatrixLhsScalarExpr;
	typedef typename MatrixRhsExpr::ScalarType MatrixRhsScalarExpr;
	typedef typename MatrixLhsExpr::RowVectorType MatrixLhsRowVector;
	typedef typename MatrixRhsExpr::ColumnVectorType MatrixRhsColumnVector;
	typedef VectorExpression<RowVectorExpr<MatrixMatrixProduct>, iNumCols> RowVectorType;
	typedef VectorExpression<ColumnVectorExpr<MatrixMatrixProduct>, iNumRows> ColumnVectorType;
	typedef typename CommonScalarType<typename BasicScalarType<MatrixLhsScalarExpr>::ScalarType,
									  typename BasicScalarType<MatrixRhsScalarExpr>::ScalarType>::ScalarType ScalarType;
    typedef typename DotTraits<MatrixLhsRowVector, MatrixRhsColumnVector, MatrixLhsExpr::iNumCols, MatrixLhsExpr::iNumCols>::ExpressionType ExpressionType;

    MatrixMatrixProduct(const MatrixLhsExpr& A, const MatrixRhsExpr& B)
        :A(A), B(B) {
    	MATVEC_ASSERT(A.iGetNumRows() == iNumRows);
    	MATVEC_ASSERT(B.iGetNumCols() == iNumCols);
    	MATVEC_ASSERT(A.iGetNumCols() == B.iGetNumRows());
    }

    ExpressionType operator()(index_type i, index_type j) const {
    	MATVEC_ASSERT(i >= 1);
    	MATVEC_ASSERT(i <= iGetNumRows());
    	MATVEC_ASSERT(j >= 1);
    	MATVEC_ASSERT(j <= iGetNumCols());
        return Dot(A.GetRow(i), B.GetCol(j));
    }

    index_type iGetNumRows() const {
    	MATVEC_ASSERT(iNumRows == A.iGetNumRows());
        return A.iGetNumRows();
    }

    index_type iGetNumCols() const {
    	MATVEC_ASSERT(iNumCols == B.iGetNumCols());
    	return B.iGetNumCols();
    }

    RowVectorType GetRow(index_type iRow) const {
    	MATVEC_ASSERT(iRow >= 1);
    	MATVEC_ASSERT(iRow <= iNumRows);
    	return RowVectorType(*this, iRow);
    }

    ColumnVectorType GetCol(index_type iCol) const {
    	MATVEC_ASSERT(iCol >= 1);
    	MATVEC_ASSERT(iCol <= iNumCols);
    	return ColumnVectorType(*this, iCol);
    }

    template <typename ScalarType2>
    bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
    	return A.bHaveReferenceTo(pFirst, pLast) || B.bHaveReferenceTo(pFirst, pLast);
    }

private:
    const MatrixLhsExpr A;
    const MatrixRhsExpr B;
};

/****************************************************************************************************************
 * matrix vector product
 ****************************************************************************************************************/

template <typename MatrixLhsExpr, index_type N_rows, index_type N_cols, typename VectorRhsExpr>
inline VectorExpression<MatrixVectorProduct<MatrixExpression<MatrixLhsExpr, N_rows, N_cols>, VectorExpression<VectorRhsExpr, N_cols> >, N_rows>
operator* (const MatrixExpression<MatrixLhsExpr, N_rows, N_cols>& A, const VectorExpression<VectorRhsExpr, N_cols>& x) {
	return VectorExpression<MatrixVectorProduct<MatrixExpression<MatrixLhsExpr, N_rows, N_cols>, VectorExpression<VectorRhsExpr, N_cols> >, N_rows>(A, x);
}

template <typename T, index_type N_rows, index_type N_cols, typename VectorRhsExpr>
inline VectorExpression<MatrixVectorProduct<MatrixExpression<MatrixDirectExpr<Matrix<T, N_rows, N_cols> >, N_rows, N_cols>, VectorExpression<VectorRhsExpr, N_cols> >, N_rows>
operator* (const Matrix<T, N_rows, N_cols>& A, const VectorExpression<VectorRhsExpr, N_cols>& x) {
	return VectorExpression<MatrixVectorProduct<MatrixExpression<MatrixDirectExpr<Matrix<T, N_rows, N_cols> >, N_rows, N_cols>, VectorExpression<VectorRhsExpr, N_cols> >, N_rows>(Direct(A), x);
}

template <typename MatrixLhsExpr, index_type N_rows, index_type N_cols, typename T>
inline VectorExpression<MatrixVectorProduct<MatrixExpression<MatrixLhsExpr, N_rows, N_cols>, VectorExpression<VectorDirectExpr<Vector<T, N_cols> >, N_cols> >, N_rows>
operator* (const MatrixExpression<MatrixLhsExpr, N_rows, N_cols>& A, const Vector<T, N_cols>& x) {
	return VectorExpression<MatrixVectorProduct<MatrixExpression<MatrixLhsExpr, N_rows, N_cols>, VectorExpression<VectorDirectExpr<Vector<T, N_cols> >, N_cols> >, N_rows>(A, Direct(x));
}

template <typename T1, typename T2, index_type N_rows, index_type N_cols>
inline VectorExpression<MatrixVectorProduct<MatrixExpression<MatrixDirectExpr<Matrix<T1, N_rows, N_cols> >, N_rows, N_cols>, VectorExpression<VectorDirectExpr<Vector<T2, N_cols> >, N_cols> >, N_rows>
operator* (const Matrix<T1, N_rows, N_cols>& A, const Vector<T2, N_cols>& x) {
	return VectorExpression<MatrixVectorProduct<MatrixExpression<MatrixDirectExpr<Matrix<T1, N_rows, N_cols> >, N_rows, N_cols>, VectorExpression<VectorDirectExpr<Vector<T2, N_cols> >, N_cols> >, N_rows>(Direct(A), Direct(x));
}

template <typename T, index_type N_rows>
inline VectorExpression<MatrixVectorProduct<MatrixExpression<MatrixDirectExpr<Matrix<T, N_rows, 3> >, N_rows, 3>, VectorExpression<Vec3DirectExpr, 3> >, N_rows>
operator* (const Matrix<T, N_rows, 3>& A, const Vec3& x) {
	return VectorExpression<MatrixVectorProduct<MatrixExpression<MatrixDirectExpr<Matrix<T, N_rows, 3> >, N_rows, 3>, VectorExpression<Vec3DirectExpr, 3> >, N_rows>(Direct(A), Direct(x));
}

template <typename T>
inline VectorExpression<MatrixVectorProduct<MatrixExpression<Mat3x3DirectExpr, 3, 3>, VectorExpression<VectorDirectExpr<Vector<T, 3> >, 3> >, 3>
operator* (const Mat3x3& A, const Vector<T, 3>& x) {
	return VectorExpression<MatrixVectorProduct<MatrixExpression<Mat3x3DirectExpr, 3, 3>, VectorExpression<VectorDirectExpr<Vector<T, 3> >, 3> >, 3>(Direct(A), Direct(x));
}

/****************************************************************************************************************
 * matrix matrix product
 ****************************************************************************************************************/

template <typename MatrixLhsExpr, index_type N_rows_Lhs, index_type N_cols_Lhs, typename MatrixRhsExpr, index_type N_cols_Rhs>
inline MatrixExpression<MatrixMatrixProduct<MatrixExpression<MatrixLhsExpr, N_rows_Lhs, N_cols_Lhs>, MatrixExpression<MatrixRhsExpr, N_cols_Lhs, N_cols_Rhs> >, N_rows_Lhs, N_cols_Rhs>
operator* (const MatrixExpression<MatrixLhsExpr, N_rows_Lhs, N_cols_Lhs>& A, const MatrixExpression<MatrixRhsExpr, N_cols_Lhs, N_cols_Rhs>& B) {
	return MatrixExpression<MatrixMatrixProduct<MatrixExpression<MatrixLhsExpr, N_rows_Lhs, N_cols_Lhs>, MatrixExpression<MatrixRhsExpr, N_cols_Lhs, N_cols_Rhs> >, N_rows_Lhs, N_cols_Rhs>(A, B);
}

template <typename MatrixLhsExpr, index_type N_rows_Lhs, index_type N_cols_Lhs, typename T, index_type N_cols_Rhs>
inline MatrixExpression<MatrixMatrixProduct<MatrixExpression<MatrixLhsExpr, N_rows_Lhs, N_cols_Lhs>, MatrixDirectExpr<Matrix<T, N_cols_Lhs, N_cols_Rhs> > >, N_rows_Lhs, N_cols_Rhs>
operator* (const MatrixExpression<MatrixLhsExpr, N_rows_Lhs, N_cols_Lhs>& A, const Matrix<T, N_cols_Lhs, N_cols_Rhs>& B) {
	return MatrixExpression<MatrixMatrixProduct<MatrixExpression<MatrixLhsExpr, N_rows_Lhs, N_cols_Lhs>, MatrixDirectExpr<Matrix<T, N_cols_Lhs, N_cols_Rhs> > >, N_rows_Lhs, N_cols_Rhs>(A, Direct(B));
}

template <typename T, index_type N_rows_Lhs, index_type N_cols_Lhs, typename MatrixRhsExpr, index_type N_cols_Rhs>
inline MatrixExpression<MatrixMatrixProduct<MatrixDirectExpr<Matrix<T, N_rows_Lhs, N_cols_Lhs> >, MatrixExpression<MatrixRhsExpr, N_cols_Lhs, N_cols_Rhs> >, N_rows_Lhs, N_cols_Rhs>
operator* (const Matrix<T, N_rows_Lhs, N_cols_Lhs>& A, const MatrixExpression<MatrixRhsExpr, N_cols_Lhs, N_cols_Rhs>& B) {
	return MatrixExpression<MatrixMatrixProduct<MatrixDirectExpr<Matrix<T, N_rows_Lhs, N_cols_Lhs> >, MatrixExpression<MatrixRhsExpr, N_cols_Lhs, N_cols_Rhs> >, N_rows_Lhs, N_cols_Rhs>(Direct(A), B);
}

template <typename T_Lhs, index_type N_rows_Lhs, index_type N_cols_Lhs, typename T_Rhs, index_type N_cols_Rhs>
inline MatrixExpression<MatrixMatrixProduct<MatrixDirectExpr<Matrix<T_Lhs, N_rows_Lhs, N_cols_Lhs> >, MatrixDirectExpr<Matrix<T_Rhs, N_cols_Lhs, N_cols_Rhs> > >, N_rows_Lhs, N_cols_Rhs>
operator* (const Matrix<T_Lhs, N_rows_Lhs, N_cols_Lhs>& A, const Matrix<T_Rhs, N_cols_Lhs, N_cols_Rhs>& B) {
	return MatrixExpression<MatrixMatrixProduct<MatrixDirectExpr<Matrix<T_Lhs, N_rows_Lhs, N_cols_Lhs> >, MatrixDirectExpr<Matrix<T_Rhs, N_cols_Lhs, N_cols_Rhs> > >, N_rows_Lhs, N_cols_Rhs>(Direct(A), Direct(B));
}

template <typename T, index_type N_rows_Lhs>
inline MatrixExpression<MatrixMatrixProduct<MatrixDirectExpr<Matrix<T, N_rows_Lhs, 3> >, Mat3x3DirectExpr>, N_rows_Lhs, 3>
operator* (const Matrix<T, N_rows_Lhs, 3>& A, const Mat3x3& B) {
	return MatrixExpression<MatrixMatrixProduct<MatrixDirectExpr<Matrix<T, N_rows_Lhs, 3> >, Mat3x3DirectExpr>, N_rows_Lhs, 3>(Direct(A), Direct(B));
}

template <typename T, index_type N_cols_Rhs>
inline MatrixExpression<MatrixMatrixProduct<Mat3x3DirectExpr, MatrixDirectExpr<Matrix<T, 3, N_cols_Rhs> > >, 3, N_cols_Rhs>
operator* (const Mat3x3& A, const Matrix<T, 3, N_cols_Rhs>& B) {
	return MatrixExpression<MatrixMatrixProduct<Mat3x3DirectExpr, MatrixDirectExpr<Matrix<T, 3, N_cols_Rhs> > >, 3, N_cols_Rhs>(Direct(A), Direct(B));
}

#define VECTOR_VECTOR_DEFINE_BINARY_FUNCTION(ExpressionName, FunctionName, FunctionClass) \
	template <typename VectorLhsExpr, typename VectorRhsExpr, index_type N_rows> \
	inline VectorExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename VectorLhsExpr::ExpressionType, typename VectorRhsExpr::ExpressionType>, VectorLhsExpr, VectorRhsExpr>, N_rows> \
	FunctionName(const VectorExpression<VectorLhsExpr, N_rows>& u, const VectorExpression<VectorRhsExpr, N_rows>& v) { \
		return VectorExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename VectorLhsExpr::ExpressionType, typename VectorRhsExpr::ExpressionType>, VectorLhsExpr, VectorRhsExpr>, N_rows>(u, v); \
	} \
	\
	template <typename T, index_type N_rows> \
	inline VectorExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename ScalarTypeTraits<T>::DirectExpressionType, typename ScalarTypeTraits<T>::DirectExpressionType>, VectorDirectExpr<Vector<T, N_rows> >, VectorDirectExpr<Vector<T, N_rows> > >, N_rows> \
	FunctionName(const Vector<T, N_rows>& u, const Vector<T, N_rows>& v) { \
		return VectorExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename ScalarTypeTraits<T>::DirectExpressionType, typename ScalarTypeTraits<T>::DirectExpressionType>, VectorDirectExpr<Vector<T, N_rows> >, VectorDirectExpr<Vector<T, N_rows> > >, N_rows>(u, v); \
	} \
	\
	template <typename T1, typename T2, index_type N_rows> \
	inline VectorExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename ScalarTypeTraits<T1>::DirectExpressionType, typename ScalarTypeTraits<T2>::DirectExpressionType>, VectorDirectExpr<Vector<T1, N_rows> >, VectorDirectExpr<Vector<T2, N_rows> > >, N_rows> \
	FunctionName(const Vector<T1, N_rows>& u, const Vector<T2, N_rows>& v) { \
		return VectorExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename ScalarTypeTraits<T1>::DirectExpressionType, typename ScalarTypeTraits<T2>::DirectExpressionType>, VectorDirectExpr<Vector<T1, N_rows> >, VectorDirectExpr<Vector<T2, N_rows> > >, N_rows>(u, v); \
	} \
	\
	template <index_type N_rows> \
	inline VectorExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, scalar_func_type, scalar_func_type>, VectorDirectExpr<Vector<scalar_func_type, N_rows> >, VectorDirectExpr<Vector<scalar_func_type, N_rows> > >, N_rows> \
	FunctionName(const Vector<scalar_func_type, N_rows>& u, const Vector<scalar_func_type, N_rows>& v) { \
		return VectorExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, scalar_func_type, scalar_func_type>, VectorDirectExpr<Vector<scalar_func_type, N_rows> >, VectorDirectExpr<Vector<scalar_func_type, N_rows> > >, N_rows>(u, v); \
	} \
	\
	template <typename VectorLhsExpr, typename T, index_type N_rows> \
	inline VectorExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename VectorLhsExpr::ExpressionType, typename ScalarTypeTraits<T>::DirectExpressionType >, VectorLhsExpr, VectorDirectExpr<Vector<T, N_rows> > >, N_rows> \
	FunctionName(const VectorExpression<VectorLhsExpr, N_rows>& u, const Vector<T, N_rows>& v) { \
		return VectorExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename VectorLhsExpr::ExpressionType, typename ScalarTypeTraits<T>::DirectExpressionType >, VectorLhsExpr, VectorDirectExpr<Vector<T, N_rows> > >, N_rows>(u, v); \
	} \
	\
	template <typename T, index_type N_rows, typename VectorRhsExpr> \
	inline VectorExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename ScalarTypeTraits<T>::DirectExpressionType, typename VectorRhsExpr::ExpressionType>, VectorDirectExpr<Vector<T, N_rows> >, VectorRhsExpr>, N_rows> \
	FunctionName(const Vector<T, N_rows>& u, const VectorExpression<VectorRhsExpr, N_rows>& v) { \
		return VectorExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename ScalarTypeTraits<T>::DirectExpressionType, typename VectorRhsExpr::ExpressionType>, VectorDirectExpr<Vector<T, N_rows> >, VectorRhsExpr>, N_rows>(u, v); \
	}

#define VECTOR_SCALAR_DEFINE_BINARY_FUNCTION(ExpressionName, FunctionName, FunctionClass) \
	template <typename VectorLhsExpr, typename ScalarRhsExpr, index_type N_rows> \
	inline VectorExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename VectorLhsExpr::ExpressionType, ScalarRhsExpr>, VectorLhsExpr, ScalarRhsExpr>, N_rows> \
	FunctionName(const VectorExpression<VectorLhsExpr, N_rows>& u, const GradientExpression<ScalarRhsExpr>& v) { \
		return VectorExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename VectorLhsExpr::ExpressionType, ScalarRhsExpr>, VectorLhsExpr, ScalarRhsExpr>, N_rows>(u, v); \
	} \
	\
	template <typename VectorLhsExpr, index_type N_rows, index_type N_SIZE> \
	inline VectorExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename VectorLhsExpr::ExpressionType, DirectExpr<Gradient<N_SIZE> > >, VectorLhsExpr, DirectExpr<Gradient<N_SIZE> > >, N_rows> \
	FunctionName(const VectorExpression<VectorLhsExpr, N_rows>& u, const Gradient<N_SIZE>& v) { \
		return VectorExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename VectorLhsExpr::ExpressionType, DirectExpr<Gradient<N_SIZE> > >, VectorLhsExpr, DirectExpr<Gradient<N_SIZE> > >, N_rows>(u, v); \
	} \
	\
	template <typename VectorLhsExpr, index_type N_rows> \
	inline VectorExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename VectorLhsExpr::ExpressionType, scalar_func_type>, VectorLhsExpr, scalar_func_type>, N_rows> \
	FunctionName(const VectorExpression<VectorLhsExpr, N_rows>& u, scalar_func_type v) { \
		return VectorExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename VectorLhsExpr::ExpressionType, scalar_func_type >, VectorLhsExpr, scalar_func_type>, N_rows>(u, v); \
	} \
	\
	template <typename T, index_type N_rows> \
	inline VectorExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename ScalarTypeTraits<T>::DirectExpressionType, scalar_func_type >, VectorDirectExpr<Vector<T, N_rows> >, scalar_func_type>, N_rows> \
	FunctionName(const Vector<T, N_rows>& u, scalar_func_type v) { \
		return VectorExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename ScalarTypeTraits<T>::DirectExpressionType, scalar_func_type >, VectorDirectExpr<Vector<T, N_rows> >, scalar_func_type>, N_rows>(u, v); \
	} \
	\
	template <typename T, index_type N_rows, index_type N_SIZE> \
	inline VectorExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename ScalarTypeTraits<T>::DirectExpressionType, DirectExpr<Gradient<N_SIZE> > >, VectorDirectExpr<Vector<T, N_rows> >, DirectExpr<Gradient<N_SIZE> > >, N_rows> \
	FunctionName(const Vector<T, N_rows>& u, const Gradient<N_SIZE>& v) { \
		return VectorExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename ScalarTypeTraits<T>::DirectExpressionType, DirectExpr<Gradient<N_SIZE> > >, VectorDirectExpr<Vector<T, N_rows> >, DirectExpr<Gradient<N_SIZE> > >, N_rows>(u, v); \
	} \
	\
	template <typename T, index_type N_rows, typename RhsExpr> \
	inline VectorExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename ScalarTypeTraits<T>::DirectExpressionType, GradientExpression<RhsExpr> >, VectorDirectExpr<Vector<T, N_rows> >, GradientExpression<RhsExpr> >, N_rows> \
	FunctionName(const Vector<T, N_rows>& u, const GradientExpression<RhsExpr>& v) { \
		return VectorExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename ScalarTypeTraits<T>::DirectExpressionType, GradientExpression<RhsExpr> >, VectorDirectExpr<Vector<T, N_rows> >, GradientExpression<RhsExpr> >, N_rows>(u, v); \
	}

#define VECTOR_DEFINE_UNARY_FUNCTION(ExpressionName, FunctionName, FunctionClass) \
	template <typename VectorExpr, index_type N_rows> \
	inline VectorExpression<ExpressionName<ScalarUnaryOperation<FunctionClass, typename VectorExpr::ExpressionType>, VectorExpr>, N_rows> \
	FunctionName(const VectorExpression<VectorExpr, N_rows>& u) { \
		return VectorExpression<ExpressionName<ScalarUnaryOperation<FunctionClass, typename VectorExpr::ExpressionType>, VectorExpr>, N_rows>(u); \
	} \
	\
	template <typename T, index_type N_rows> \
	inline VectorExpression<ExpressionName<ScalarUnaryOperation<FunctionClass, typename ScalarTypeTraits<T>::DirectExpressionType>, VectorDirectExpr<Vector<T, N_rows> > >, N_rows> \
	FunctionName(const Vector<T, N_rows>& u) { \
		return VectorExpression<ExpressionName<ScalarUnaryOperation<FunctionClass, typename ScalarTypeTraits<T>::DirectExpressionType>, VectorDirectExpr<Vector<T, N_rows> > >, N_rows>(u); \
	}

#define MATRIX_DEFINE_UNARY_FUNCTION(ExpressionName, FunctionName, FunctionClass) \
	template <typename MatrixExpr, index_type N_rows, index_type N_cols> \
	inline MatrixExpression<ExpressionName<ScalarUnaryOperation<FunctionClass, typename MatrixExpr::ExpressionType>, MatrixExpr>, N_rows, N_cols> \
	FunctionName(const MatrixExpression<MatrixExpr, N_rows, N_cols>& u) { \
		return MatrixExpression<ExpressionName<ScalarUnaryOperation<FunctionClass, typename MatrixExpr::ExpressionType>, MatrixExpr>, N_rows, N_cols>(u); \
	} \
	\
	template <typename T, index_type N_rows, index_type N_cols> \
	inline MatrixExpression<ExpressionName<ScalarUnaryOperation<FunctionClass, typename ScalarTypeTraits<T>::DirectExpressionType>, MatrixDirectExpr<Matrix<T, N_rows, N_cols> > >, N_rows, N_cols> \
	FunctionName(const Matrix<T, N_rows, N_cols>& u) { \
		return MatrixExpression<ExpressionName<ScalarUnaryOperation<FunctionClass, typename ScalarTypeTraits<T>::DirectExpressionType>, MatrixDirectExpr<Matrix<T, N_rows, N_cols> > >, N_rows, N_cols>(u); \
	}

#define MATRIX_MATRIX_DEFINE_BINARY_FUNCTION(ExpressionName, FunctionName, FunctionClass) \
	template <typename MatrixLhsExpr, typename MatrixRhsExpr, index_type N_rows, index_type N_cols> \
	inline MatrixExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename MatrixLhsExpr::ExpressionType, typename MatrixRhsExpr::ExpressionType>, MatrixLhsExpr, MatrixRhsExpr>, N_rows, N_cols> \
	FunctionName(const MatrixExpression<MatrixLhsExpr, N_rows, N_cols>& u, const MatrixExpression<MatrixRhsExpr, N_rows, N_cols>& v) { \
		return MatrixExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename MatrixLhsExpr::ExpressionType, typename MatrixRhsExpr::ExpressionType>, MatrixLhsExpr, MatrixRhsExpr>, N_rows, N_cols>(u, v); \
	} \
	\
	template <typename T, index_type N_rows, index_type N_cols> \
	inline MatrixExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename ScalarTypeTraits<T>::DirectExpressionType, typename ScalarTypeTraits<T>::DirectExpressionType>, MatrixDirectExpr<Matrix<T, N_rows, N_cols> >, MatrixDirectExpr<Matrix<T, N_rows, N_cols> > >, N_rows, N_cols> \
	FunctionName(const Matrix<T, N_rows, N_cols>& u, const Matrix<T, N_rows, N_cols>& v) { \
		return MatrixExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename ScalarTypeTraits<T>::DirectExpressionType, typename ScalarTypeTraits<T>::DirectExpressionType>, MatrixDirectExpr<Matrix<T, N_rows, N_cols> >, MatrixDirectExpr<Matrix<T, N_rows, N_cols> > >, N_rows, N_cols>(u, v); \
	} \
	\
	template <typename MatrixLhsExpr, typename T, index_type N_rows, index_type N_cols> \
	inline MatrixExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename MatrixLhsExpr::ExpressionType, typename ScalarTypeTraits<T>::DirectExpressionType >, MatrixLhsExpr, MatrixDirectExpr<Matrix<T, N_rows, N_cols> > >, N_rows, N_cols> \
	FunctionName(const MatrixExpression<MatrixLhsExpr, N_rows, N_cols>& u, const Matrix<T, N_rows, N_cols>& v) { \
		return MatrixExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename MatrixLhsExpr::ExpressionType, typename ScalarTypeTraits<T>::DirectExpressionType >, MatrixLhsExpr, MatrixDirectExpr<Matrix<T, N_rows, N_cols> > >, N_rows, N_cols>(u, v); \
	} \
	\
	template <typename T, index_type N_rows, index_type N_cols, typename MatrixRhsExpr> \
	inline MatrixExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename ScalarTypeTraits<T>::DirectExpressionType, typename MatrixRhsExpr::ExpressionType>, MatrixDirectExpr<Matrix<T, N_rows, N_cols> >, MatrixRhsExpr>, N_rows, N_cols> \
	FunctionName(const Matrix<T, N_rows, N_cols>& u, const MatrixExpression<MatrixRhsExpr, N_rows, N_cols>& v) { \
		return MatrixExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename ScalarTypeTraits<T>::DirectExpressionType, typename MatrixRhsExpr::ExpressionType>, MatrixDirectExpr<Matrix<T, N_rows, N_cols> >, MatrixRhsExpr>, N_rows, N_cols>(u, v); \
	}

#define MATRIX_SCALAR_DEFINE_BINARY_FUNCTION(ExpressionName, FunctionName, FunctionClass) \
	template <typename MatrixLhsExpr, typename ScalarRhsExpr, index_type N_rows, index_type N_cols> \
	inline MatrixExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename MatrixLhsExpr::ExpressionType, ScalarRhsExpr>, MatrixLhsExpr, ScalarRhsExpr>, N_rows, N_cols> \
	FunctionName(const MatrixExpression<MatrixLhsExpr, N_rows, N_cols>& u, const GradientExpression<ScalarRhsExpr>& v) { \
		return MatrixExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename MatrixLhsExpr::ExpressionType, ScalarRhsExpr>, MatrixLhsExpr, ScalarRhsExpr>, N_rows, N_cols>(u, v); \
	} \
	\
	template <typename MatrixLhsExpr, index_type N_rows, index_type N_cols, index_type N_SIZE> \
	inline MatrixExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename MatrixLhsExpr::ExpressionType, DirectExpr<Gradient<N_SIZE> > >, MatrixLhsExpr, DirectExpr<Gradient<N_SIZE> > >, N_rows, N_cols> \
	FunctionName(const MatrixExpression<MatrixLhsExpr, N_rows, N_cols>& u, const Gradient<N_SIZE>& v) { \
		return MatrixExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename MatrixLhsExpr::ExpressionType, DirectExpr<Gradient<N_SIZE> > >, MatrixLhsExpr, DirectExpr<Gradient<N_SIZE> > >, N_rows, N_cols>(u, v); \
	} \
	\
	template <typename MatrixLhsExpr, index_type N_rows, index_type N_cols> \
	inline MatrixExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename MatrixLhsExpr::ExpressionType, scalar_func_type>, MatrixLhsExpr, scalar_func_type>, N_rows, N_cols> \
	FunctionName(const MatrixExpression<MatrixLhsExpr, N_rows, N_cols>& u, scalar_func_type v) { \
		return MatrixExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename MatrixLhsExpr::ExpressionType, scalar_func_type >, MatrixLhsExpr, scalar_func_type>, N_rows, N_cols>(u, v); \
	} \
	\
	template <typename T, index_type N_rows, index_type N_cols> \
	inline MatrixExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename ScalarTypeTraits<T>::DirectExpressionType, scalar_func_type >, MatrixDirectExpr<Matrix<T, N_rows, N_cols> >, scalar_func_type>, N_rows, N_cols> \
	FunctionName(const Matrix<T, N_rows, N_cols>& u, scalar_func_type v) { \
		return MatrixExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename ScalarTypeTraits<T>::DirectExpressionType, scalar_func_type >, MatrixDirectExpr<Matrix<T, N_rows, N_cols> >, scalar_func_type>, N_rows, N_cols>(u, v); \
	} \
	\
	template <typename T, index_type N_rows, index_type N_cols, index_type N_SIZE> \
	inline MatrixExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename ScalarTypeTraits<T>::DirectExpressionType, DirectExpr<Gradient<N_SIZE> > >, MatrixDirectExpr<Matrix<T, N_rows, N_cols> >, DirectExpr<Gradient<N_SIZE> > >, N_rows, N_cols> \
	FunctionName(const Matrix<T, N_rows, N_cols>& u, const Gradient<N_SIZE>& v) { \
		return MatrixExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename ScalarTypeTraits<T>::DirectExpressionType, DirectExpr<Gradient<N_SIZE> > >, MatrixDirectExpr<Matrix<T, N_rows, N_cols> >, DirectExpr<Gradient<N_SIZE> > >, N_rows, N_cols>(u, v); \
	} \
	\
	template <typename T, index_type N_rows, index_type N_cols, typename RhsExpr> \
	inline MatrixExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename ScalarTypeTraits<T>::DirectExpressionType, GradientExpression<RhsExpr> >, MatrixDirectExpr<Matrix<T, N_rows, N_cols> >, GradientExpression<RhsExpr> >, N_rows, N_cols> \
	FunctionName(const Matrix<T, N_rows, N_cols>& u, const GradientExpression<RhsExpr>& v) { \
		return MatrixExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename ScalarTypeTraits<T>::DirectExpressionType, GradientExpression<RhsExpr> >, MatrixDirectExpr<Matrix<T, N_rows, N_cols> >, GradientExpression<RhsExpr> >, N_rows, N_cols>(u, v); \
	}

VECTOR_VECTOR_DEFINE_BINARY_FUNCTION(VectorVectorVectorBinaryExpr, operator+, FuncPlus)
VECTOR_VECTOR_DEFINE_BINARY_FUNCTION(VectorVectorVectorBinaryExpr, operator-, FuncMinus)
VECTOR_SCALAR_DEFINE_BINARY_FUNCTION(VectorScalarVectorBinaryExpr, operator*, FuncMult)
VECTOR_SCALAR_DEFINE_BINARY_FUNCTION(VectorScalarVectorBinaryExpr, operator/, FuncDiv)
VECTOR_DEFINE_UNARY_FUNCTION(VectorVectorUnaryExpr, operator-, FuncUnaryMinus)

MATRIX_MATRIX_DEFINE_BINARY_FUNCTION(MatrixMatrixMatrixBinaryExpr, operator+, FuncPlus)
MATRIX_MATRIX_DEFINE_BINARY_FUNCTION(MatrixMatrixMatrixBinaryExpr, operator-, FuncMinus)
MATRIX_SCALAR_DEFINE_BINARY_FUNCTION(MatrixScalarMatrixBinaryExpr, operator*, FuncMult)
MATRIX_SCALAR_DEFINE_BINARY_FUNCTION(MatrixScalarMatrixBinaryExpr, operator/, FuncDiv)
MATRIX_DEFINE_UNARY_FUNCTION(MatrixMatrixUnaryExpr, operator-, FuncUnaryMinus)

#undef VECTOR_VECTOR_DEFINE_BINARY_FUNCTION
#undef VECTOR_SCALAR_DEFINE_BINARY_FUNCTION
#undef VECTOR_DEFINE_UNARY_FUNCTION

#undef MATRIX_MATRIX_DEFINE_BINARY_FUNCTION
#undef MATRIX_SCALAR_DEFINE_BINARY_FUNCTION
#undef MATRIX_DEFINE_UNARY_FUNCTION
}

#endif
