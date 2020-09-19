/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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
        Copyright (C) 2013(-2020) all rights reserved.

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
#include "matvec3n.h"
#include "RotCoeff.hh"
#include "Rot.hh"

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
    struct VectorSize {
	static const int N = -1; // Note must not be zero because that would mean variable size!
};

template <>
    struct VectorSize<scalar_func_type> {
	static const int N = 1;
};

template <>
    struct VectorSize<Vec3> {
	static const int N = 3;
};

template <>
    struct VectorSize<Vec6> {
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

namespace MatVecHelp
{
	template <typename T>
	struct AliasTypeHelper {

	};

	template <>
	struct AliasTypeHelper<scalar_func_type> {
		static const bool bAlias = false;
	};

	template <typename Expression>
	struct AliasTypeHelper<GradientExpression<Expression> > {
		static const bool bAlias = GradientExpression<Expression>::bAlias;
	};

	template <index_type N_SIZE, bool ALIAS>
	struct AliasTypeHelper<DirectExpr<Gradient<N_SIZE>, ALIAS> > {
		static const bool bAlias = DirectExpr<Gradient<N_SIZE>, ALIAS>::bAlias;
	};

	template <index_type N_SIZE>
	struct AliasTypeHelper<Gradient<N_SIZE> > {
		static const bool bAlias = false;
	};

	template <typename BinFunc, typename LhsExpr, typename RhsExpr>
	struct AliasTypeHelper<BinaryExpr<BinFunc, LhsExpr, RhsExpr> > {
		static const bool bAlias = BinaryExpr<BinFunc, LhsExpr, RhsExpr>::bAlias;
	};

	template <typename UnFunc, typename Expr>
	struct AliasTypeHelper<UnaryExpr<UnFunc, Expr> > {
		static const bool bAlias = UnaryExpr<UnFunc, Expr>::bAlias;
	};

	template <bool bAlias>
	struct ApplyAliasHelperMatrix;

	template <>
	struct ApplyAliasHelperMatrix<false> {
	    template <typename MatrixType, typename Func, typename Expression>
	    static inline void ApplyMatrixFunc(MatrixType& A, const Expression& B, const Func& f) {
	    	A.ApplyMatrixFuncNoAlias(B, f);
	    }
	};

	template <>
	struct ApplyAliasHelperMatrix<true> {
	    template <typename MatrixType, typename Func, typename Expression>
		static inline void ApplyMatrixFunc(MatrixType& A, const Expression& B, const Func& f) {
			A.ApplyMatrixFuncAlias(B, f);
		}
	};

	struct Assign {
		template <typename T, typename U>
		static inline void Eval(T& b, const U& a) {
			b = a;
		}
	};

	struct Add {
		template <typename T, typename U>
		static inline void Eval(T& b, const U& a) {
			b += a;
		}
	};

	struct Sub {
		template <typename T, typename U>
		static inline void Eval(T& b, const U& a) {
			b -= a;
		}
	};

	struct Mul {
		template <typename T, typename U>
		static inline void Eval(T& b, const U& a) {
			b *= a;
		}
	};

	struct Div {
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

        template <typename ScalarType>
        explicit VectorExpression(const ScalarType* pArray, index_type iRows, index_type iOffset)
            :Expression(pArray, iRows, iOffset) {
#if MATVEC_DEBUG > 0
            AssertValid();
#endif            
        }

        index_type iGetNumRows() const { return Expression::iGetNumRows(); }

private:
#if MATVEC_DEBUG > 0
    void AssertValid() const {
            MATVEC_ASSERT((Expression::iGetNumRows() == iNumRows) || (iNumRows == DYNAMIC_SIZE && Expression::iGetNumRows() >= 0));
    }
#endif
     static_assert(iNumRows == Expression::iNumRows);
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

        index_type iGetNumRows() const {
            return Expression::iGetNumRows();
    }
        index_type iGetNumCols() const {
            return Expression::iGetNumCols();
    }

private:

#if MATVEC_DEBUG > 0
    void AssertValid() const {
            MATVEC_ASSERT((Expression::iGetNumRows() == iNumRows) || (iNumRows == DYNAMIC_SIZE && Expression::iGetNumRows() >= 0));
            MATVEC_ASSERT((Expression::iGetNumCols() == iNumCols) || (iNumCols == DYNAMIC_SIZE && Expression::iGetNumCols() >= 0));
    }
#endif

     static_assert(iNumRows == DYNAMIC_SIZE || Expression::iNumRows == DYNAMIC_SIZE ? true : iNumRows == Expression::iNumRows);
     static_assert(iNumCols == DYNAMIC_SIZE || Expression::iNumCols == DYNAMIC_SIZE ? true : iNumCols == Expression::iNumCols);
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
            MATVEC_ASSERT((oU.iGetNumRows() == iNumRows) || (iNumRows == DYNAMIC_SIZE && oU.iGetNumRows() >= 0));
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
            MATVEC_ASSERT((oU.iGetNumRows() == iNumRows) || (iNumRows == DYNAMIC_SIZE && oU.iGetNumRows() >= 0));
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
            MATVEC_ASSERT(i <= iGetNumRows());
        return ScalarUnaryFunc::f(oU(i));
    }

    index_type iGetNumRows() const {
            MATVEC_ASSERT((oU.iGetNumRows() == iNumRows) || (iNumRows == DYNAMIC_SIZE && oU.iGetNumRows() >= 0));
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
            MATVEC_ASSERT((oU.iGetNumRows() == VectorExpr::iNumRows) || (VectorExpr::iNumRows == DYNAMIC_SIZE && oU.iGetNumRows() >= 0));
        return iNumRows;
    }

    template <typename ScalarType2>
    bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
    	return oU.bHaveReferenceTo(pFirst, pLast);
    }

private:
     static_assert(iEndIndex <= VectorExpr::iNumRows);
     static_assert(iStartIndex >= 1);
     static_assert(VectorExpr::iNumRows != DYNAMIC_SIZE);
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
            MATVEC_ASSERT(i <= iGetNumRows());
        return oU(i);
    }
    
    index_type iGetNumRows() const {
            MATVEC_ASSERT((oU.iGetNumRows() == iNumRows) || (iNumRows == DYNAMIC_SIZE && oU.iGetNumRows() >= 0));
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

        SliceVector(const ScalarType* p, index_type iRows, index_type iOffset)
    	:pVec(p) {
            MATVEC_ASSERT((iRows == iNumRows) || (iNumRows == DYNAMIC_SIZE && iRows >= 0));
            MATVEC_ASSERT(iOffset == N_offset);
        }

        const ScalarType& operator()(index_type iRow) const {
            --iRow;	// row index is 1-based
            MATVEC_ASSERT(iRow >= 0);
            MATVEC_ASSERT(iRow < iGetNumRows());
            return *(pVec + iRow * N_offset);
        }

        index_type iGetNumRows() const { return iNumRows; }

        template <typename ScalarType2>
        bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
            return bArrayOverlap(pVec, pVec + (iNumRows - 1) * N_offset, pFirst, pLast);
        }

    private:
        const ScalarType* const pVec;
    
        static_assert(iNumRows != DYNAMIC_SIZE);
        static_assert(N_offset != DYNAMIC_SIZE);
    };

    template <typename T, index_type N_rows>
    class SliceVector<T, N_rows, DYNAMIC_SIZE> {
    public:
        static const bool bAlias = MatVecHelp::AliasTypeHelper<T>::bAlias;
        static const index_type iNumRows = N_rows;
        typedef typename ScalarTypeTraits<T>::ScalarType ScalarType;
        typedef typename ScalarTypeTraits<T>::DirectExpressionType ExpressionType;

        SliceVector(const ScalarType* p, index_type iRows, index_type iOffset)
            :pVec(p), iOffset(iOffset) {
            
            MATVEC_ASSERT((iRows == iNumRows) || (iNumRows == DYNAMIC_SIZE && iRows >= 0));
        }

        const ScalarType& operator()(index_type iRow) const {
            --iRow;	// row index is 1-based
            MATVEC_ASSERT(iRow >= 0);
            MATVEC_ASSERT(iRow < iGetNumRows());
            return *(pVec + iRow * iOffset);
        }

        index_type iGetNumRows() const { return iNumRows; }

        template <typename ScalarType2>
        bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
            return bArrayOverlap(pVec, pVec + (iNumRows - 1) * iOffset, pFirst, pLast);
        }

    private:
        const ScalarType* const pVec;
        const index_type iOffset;
    
	 static_assert(iNumRows != DYNAMIC_SIZE);
    };

    template <typename T, index_type N_offset>
    class SliceVector<T, DYNAMIC_SIZE, N_offset> {
    public:
        static const bool bAlias = MatVecHelp::AliasTypeHelper<T>::bAlias;
        static const index_type iNumRows = DYNAMIC_SIZE;
        typedef typename ScalarTypeTraits<T>::ScalarType ScalarType;
        typedef typename ScalarTypeTraits<T>::DirectExpressionType ExpressionType;

        SliceVector(const ScalarType* p, index_type iRows, index_type iOffset)
            :pVec(p), iCurrRows(iRows) {
            
            MATVEC_ASSERT(iOffset == N_offset);
    }

    const ScalarType& operator()(index_type iRow) const {
    	--iRow;	// row index is 1-based
    	MATVEC_ASSERT(iRow >= 0);
    	MATVEC_ASSERT(iRow < iGetNumRows());
        return *(pVec + iRow * N_offset);
    }

        index_type iGetNumRows() const { return iCurrRows; }

        template <typename ScalarType2>
        bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
            return bArrayOverlap(pVec, pVec + (iCurrRows - 1) * N_offset, pFirst, pLast);
        }

    private:
        const ScalarType* const pVec;
        const index_type iCurrRows;
    
	 static_assert(N_offset != DYNAMIC_SIZE);
    };

    
    template <typename T>
    class SliceVector<T, DYNAMIC_SIZE, DYNAMIC_SIZE> {
    public:
        static const bool bAlias = MatVecHelp::AliasTypeHelper<T>::bAlias;
        static const index_type iNumRows = DYNAMIC_SIZE;
        typedef typename ScalarTypeTraits<T>::ScalarType ScalarType;
        typedef typename ScalarTypeTraits<T>::DirectExpressionType ExpressionType;

        SliceVector(const ScalarType* p, index_type iRows, index_type iOffset)
            :pVec(p), iCurrRows(iRows), iOffset(iOffset) {

        }

        const ScalarType& operator()(index_type iRow) const {
            --iRow;	// row index is 1-based
            MATVEC_ASSERT(iRow >= 0);
            MATVEC_ASSERT(iRow < iGetNumRows());
            return *(pVec + iRow * iOffset);
        }

        index_type iGetNumRows() const { return iCurrRows; }

    template <typename ScalarType2>
    bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
            return bArrayOverlap(pVec, pVec + (iCurrRows - 1) * iOffset, pFirst, pLast);
    }

private:
    const ScalarType* const pVec;
        const index_type iCurrRows;
        const index_type iOffset;
    
	 static_assert(iNumRows == DYNAMIC_SIZE);
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
            MATVEC_ASSERT((iNumRows == A.iGetNumRows()) || (iNumRows == DYNAMIC_SIZE && A.iGetNumRows() >= 0));
	}

    ExpressionType operator()(index_type iRow) const {
    	MATVEC_ASSERT(iRow >= 1);
            MATVEC_ASSERT(iRow <= iGetNumRows());
        return A(iRow, iCol);
    }

        index_type iGetNumRows() const { return A.iGetNumRows(); }

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
            MATVEC_ASSERT((iNumRows == A.iGetNumCols()) || (iNumRows == DYNAMIC_SIZE && A.iGetNumCols() >= 0));
	}

    ExpressionType operator()(index_type iCol) const {
    	MATVEC_ASSERT(iCol >= 1);
            MATVEC_ASSERT(iCol <= iGetNumRows());
        return A(iRow, iCol);
    }

        index_type iGetNumRows() const { return A.iGetNumCols(); }

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
            MATVEC_ASSERT((iNumRows == A.iGetNumCols()) || (iNumRows == DYNAMIC_SIZE && A.iGetNumCols() >= 0));
            MATVEC_ASSERT((iNumCols == A.iGetNumRows()) || (iNumCols == DYNAMIC_SIZE && A.iGetNumRows() >= 0));
	}

    const ScalarType& operator()(index_type i, index_type j) const {
    	MATVEC_ASSERT(i >= 1);
    	MATVEC_ASSERT(i <= iGetNumRows());
    	MATVEC_ASSERT(j >= 1);
    	MATVEC_ASSERT(j <= iGetNumCols());
        return A(j, i);
    }

    index_type iGetNumRows() const {
            MATVEC_ASSERT((iNumRows == A.iGetNumCols()) || (iNumRows == DYNAMIC_SIZE && A.iGetNumCols() >= 0));
        return A.iGetNumCols();
    }

    index_type iGetNumCols() const {
            MATVEC_ASSERT((iNumCols == A.iGetNumRows()) || (iNumCols == DYNAMIC_SIZE && A.iGetNumRows() >= 0));
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
            MATVEC_ASSERT((MatrixExpr::iNumCols == A.iGetNumCols()) || (MatrixExpr::iNumCols == DYNAMIC_SIZE && A.iGetNumCols() >= 0));
            MATVEC_ASSERT((MatrixExpr::iNumRows == A.iGetNumRows()) || (MatrixExpr::iNumRows == DYNAMIC_SIZE && A.iGetNumRows() >= 0));
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
     static_assert(iRowStart >= 1);
     static_assert(iRowEnd <= MatrixExpr::iNumRows);
     static_assert(iColStart >= 1);
     static_assert(iColEnd <= MatrixExpr::iNumCols);
     static_assert(iNumRows != DYNAMIC_SIZE);
     static_assert(iNumCols != DYNAMIC_SIZE);
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

            MATVEC_ASSERT(oU.iGetNumRows() == oV.iGetNumRows());
            MATVEC_ASSERT(oU.iGetNumCols() == oV.iGetNumCols());
            
            MATVEC_ASSERT((iNumRows == DYNAMIC_SIZE && oU.iGetNumRows() >= 0) || (iNumRows == oU.iGetNumRows()));
            MATVEC_ASSERT((iNumCols == DYNAMIC_SIZE && oU.iGetNumCols() >= 0) || (iNumCols == oU.iGetNumCols()));
            MATVEC_ASSERT((iNumRows == DYNAMIC_SIZE && oV.iGetNumRows() >= 0) || (iNumRows == oV.iGetNumRows()));
            MATVEC_ASSERT((iNumCols == DYNAMIC_SIZE && oV.iGetNumCols() >= 0) || (iNumCols == oV.iGetNumCols()));
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
            MATVEC_ASSERT((iNumRows == DYNAMIC_SIZE && oU.iGetNumRows() >= 0) || (iNumRows == oU.iGetNumRows()));
        return oU.iGetNumRows();
    }

    index_type iGetNumCols() const {
    	MATVEC_ASSERT(oU.iGetNumCols() == oV.iGetNumCols());
            MATVEC_ASSERT((iNumCols == DYNAMIC_SIZE && oU.iGetNumCols() >= 0) || (iNumCols == oU.iGetNumCols()));
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
     static_assert((MatrixLhsExpr::iNumRows == DYNAMIC_SIZE || MatrixRhsExpr::iNumRows == DYNAMIC_SIZE) ? true : MatrixLhsExpr::iNumRows == MatrixRhsExpr::iNumRows);
     static_assert((MatrixLhsExpr::iNumCols == DYNAMIC_SIZE || MatrixRhsExpr::iNumCols == DYNAMIC_SIZE) ? true : MatrixLhsExpr::iNumCols == MatrixRhsExpr::iNumCols);
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

            MATVEC_ASSERT((iNumRows == oU.iGetNumRows()) || (iNumRows == DYNAMIC_SIZE && oU.iGetNumRows() >= 0));
            MATVEC_ASSERT((iNumCols == oU.iGetNumCols()) || (iNumCols == DYNAMIC_SIZE && oU.iGetNumCols() >= 0)) ;
    }

    ExpressionType operator()(index_type i, index_type j) const {
    	MATVEC_ASSERT(i >= 1);
    	MATVEC_ASSERT(i <= iGetNumRows());
    	MATVEC_ASSERT(j >= 1);
    	MATVEC_ASSERT(j <= iGetNumCols());
        return ScalarBinFunc::f(oU(i, j), oV);
    }

    index_type iGetNumRows() const {
            MATVEC_ASSERT((iNumRows == oU.iGetNumRows()) || (iNumRows == DYNAMIC_SIZE && oU.iGetNumRows() >= 0));
        return oU.iGetNumRows();
    }

    index_type iGetNumCols() const {
            MATVEC_ASSERT((iNumCols == oU.iGetNumCols()) || (iNumCols == DYNAMIC_SIZE && oU.iGetNumCols() >= 0));
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
            MATVEC_ASSERT((iNumRows == oU.iGetNumRows()) || (iNumRows == DYNAMIC_SIZE && oU.iGetNumRows() >= 0));
            MATVEC_ASSERT((iNumCols == oU.iGetNumCols()) || (iNumCols == DYNAMIC_SIZE && oU.iGetNumCols() >= 0));
    }

    ExpressionType operator()(index_type i, index_type j) const {
    	MATVEC_ASSERT(i >= 1);
    	MATVEC_ASSERT(i <= iGetNumRows());
    	MATVEC_ASSERT(j >= 1);
    	MATVEC_ASSERT(j <= iGetNumCols());
        return ScalarUnaryFunc::f(oU(i, j));
    }

    index_type iGetNumRows() const {
            MATVEC_ASSERT((iNumRows == oU.iGetNumRows()) || (iNumRows == DYNAMIC_SIZE && oU.iGetNumRows() >= 0));
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
            MATVEC_ASSERT((iNumRows == A.iGetNumRows()) || (iNumRows == DYNAMIC_SIZE && A.iGetNumRows() >= 0));
            MATVEC_ASSERT((iNumCols == A.iGetNumCols()) || (iNumCols == DYNAMIC_SIZE && A.iGetNumCols() >= 0));
    }

    const ScalarType& operator()(index_type i, index_type j) const {
    	MATVEC_ASSERT(i >= 1);
    	MATVEC_ASSERT(i <= iGetNumRows());
    	MATVEC_ASSERT(j >= 1);
    	MATVEC_ASSERT(j <= iGetNumCols());
        return A(i, j);
    }

    index_type iGetNumRows() const {
            MATVEC_ASSERT((iNumRows == A.iGetNumRows()) || (iNumRows == DYNAMIC_SIZE && A.iGetNumRows() >= 0));
        return A.iGetNumRows();
    }

    index_type iGetNumCols() const {
            MATVEC_ASSERT((iNumCols == A.iGetNumCols()) || (iNumCols == DYNAMIC_SIZE && A.iGetNumCols() >= 0));
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
            return RowVectorType(A.pGetMat() + iRow, iGetNumCols(), iGetNumRows());
    }

    ColumnVectorType GetCol(index_type iCol) const {
    	MATVEC_ASSERT(iCol >= 1);
    	MATVEC_ASSERT(iCol <= iGetNumCols());
    	--iCol;
            return ColumnVectorType(A.pGetMat() + iCol * iNumRows, iGetNumRows(), 1);
    }

    template <typename ScalarType2>
    bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
    	return bArrayOverlap(A.pGetMat(), A.pGetMat() + (iNumRows - 1) * (iNumCols - 1), pFirst, pLast);
    }

private:
    const Mat3x3& A;
};

    class Mat3xNDirectExpr {
public:
        static const bool bAlias = false;
        static const index_type iNumRows = 3;
        static const index_type iNumCols = DYNAMIC_SIZE;

        typedef void RowVectorType;
        typedef void ColumnVectorType;
        typedef ScalarTypeTraits<doublereal>::ScalarType ScalarType;
        typedef ScalarTypeTraits<doublereal>::DirectExpressionType ExpressionType;

        Mat3xNDirectExpr(const Mat3xN& A)
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
            return A.iGetNumCols();
        }

        template <typename ScalarType2>
        bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
            return bArrayOverlap(&A(1, 1), &A(A.iGetNumRows(), A.iGetNumCols()), pFirst, pLast);
	}

    private:
        const Mat3xN& A;
};

    class MatNxNDirectExpr {
public:
        static const bool bAlias = false;
        static const index_type iNumRows = DYNAMIC_SIZE;
        static const index_type iNumCols = DYNAMIC_SIZE;

        typedef void RowVectorType;
        typedef void ColumnVectorType;
        typedef ScalarTypeTraits<doublereal>::ScalarType ScalarType;
        typedef ScalarTypeTraits<doublereal>::DirectExpressionType ExpressionType;

        MatNxNDirectExpr(const MatNxN& A)
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
            return A.iGetNumRows();
        }

        index_type iGetNumCols() const {
            return A.iGetNumCols();
    }

        template <typename ScalarType2>
        bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
            return bArrayOverlap(&A(1, 1), &A(A.iGetNumRows(), A.iGetNumCols()), pFirst, pLast);
        }

    private:
        const MatNxN& A;
    };
    
    template <typename InitClass, typename T, index_type N_rows, index_type N_cols>
    class MatrixInit: public InitClass {
    public:
	template <typename InitArg>
	explicit MatrixInit(const InitArg& v)
            :InitClass(v) {

    }

	template <typename InitArg1, typename InitArg2>
	explicit MatrixInit(const InitArg1& v1, const InitArg2& a2)
            :InitClass(v1, a2) {

    }
    };

    template <typename InitClass, typename T, index_type N_rows>
    class VectorInit: public InitClass {
    public:
	template <typename InitArg>
	explicit VectorInit(const InitArg& R)
            :InitClass(R) {

    }
    };    

    template <typename T, index_type N_rows>
    class VectorData {
    public:
        explicit VectorData(index_type N, bool bZeroInit) {
            MATVEC_ASSERT(N == N_rows);

            if (bZeroInit) {
                ZeroInit(&rgData[0], &rgData[N_rows]);
    		}
    	}

        VectorData(const VectorData& oData) {
            std::copy(&oData.rgData[0], &oData.rgData[N_rows], &rgData[0]);
    }

        template <typename U>
        VectorData(const VectorData<U, N_rows>& oData) {
            std::copy(&oData.rgData[0], &oData.rgData[N_rows], &rgData[0]);
			}

        T& operator[](index_type i) {
            MATVEC_ASSERT(i >= 0);
            MATVEC_ASSERT(i < N_rows);
            return rgData[i];
		}

        const T& operator[](integer i) const {
            MATVEC_ASSERT(i >= 0);
            MATVEC_ASSERT(i < N_rows);
            return rgData[i];
    }

        index_type iGetNumRows() const {
            return N_rows;
        }

        void Resize(index_type iNumRows) {
            MATVEC_ASSERT(iNumRows == N_rows);
        }

        const T* pGetFirstElem() const {
            return &rgData[0];
    }

        const T* pGetLastElem() const {
            return &rgData[N_rows - 1];
        }

    private:
        T rgData[N_rows];
    };

    template <typename T>
    class VectorData<T, DYNAMIC_SIZE> {
    public:
        explicit VectorData(index_type N, bool)
            :rgData(N) {
    }

        VectorData(const VectorData& oVec)
            :rgData(oVec.rgData) {
        }

        T& operator[](index_type i) {
            MATVEC_ASSERT(i >= 0);
            MATVEC_ASSERT(i < index_type(rgData.size()));
            return rgData[i];
        }

        const T& operator[](integer i) const {
            MATVEC_ASSERT(i >= 0);
            MATVEC_ASSERT(i < index_type(rgData.size()));
            return rgData[i];
        }

        index_type iGetNumRows() const {
            return rgData.size();
        }

        void Resize(index_type iNumRows) {
            if (iNumRows != iGetNumRows()) {
                rgData.resize(iNumRows);
            }
        }

        const T* pGetFirstElem() const {
            return &rgData.front();
        }

        const T* pGetLastElem() const {
            return &rgData.back();
        }
        
    private:
        std::vector<T, GradientAllocator<T> > rgData;
    };

    
    template <typename T, index_type N_rows, index_type N_cols>
    class MatrixData
    {
    public:
        MatrixData(index_type iNumRows, index_type iNumCols, bool bZeroInit) {
            MATVEC_ASSERT(iNumRows == N_rows);
            MATVEC_ASSERT(iNumCols == N_cols);

            if (bZeroInit) {
                ZeroInit(&rgData[0][0], &rgData[iNumRows - 1][iNumCols - 1] + 1);
            }
        }

        index_type iGetNumRows() const { return N_rows; }
        index_type iGetNumCols() const { return N_cols; }
    
        T& operator() (index_type i, index_type j) {
            MATVEC_ASSERT(i >= 1 && i <= N_rows);
            MATVEC_ASSERT(j >= 1 && j <= N_cols);
            return rgData[i - 1][j - 1];
        }

        const T& operator() (index_type i, index_type j) const {
            MATVEC_ASSERT(i >= 1 && i <= N_rows);
            MATVEC_ASSERT(j >= 1 && j <= N_cols);        
            return rgData[i - 1][j - 1];
        }

        void Resize(index_type iRows, index_type iCols) {
            MATVEC_ASSERT(iRows == N_rows);
            MATVEC_ASSERT(iCols == N_cols);
        }

        const T* pGetFirstElem() const {
            return &rgData[0][0];
        }

        const T* pGetLastElem() const {
            return &rgData[N_rows - 1][N_cols - 1];
        }
        
    private:
        T rgData[N_rows][N_cols];
    };

    template <typename T>
    class MatrixData<T, DYNAMIC_SIZE, DYNAMIC_SIZE>
    {
    public:
        MatrixData(index_type iNumRows, index_type iNumCols, bool bZeroInit)
            :iNumRows(iNumRows),
             iNumCols(iNumCols),
             rgData(iNumRows * iNumCols, bZeroInit) {
        }

        index_type iGetNumRows() const { return iNumRows; }
        index_type iGetNumCols() const { return iNumCols; }       
        
        T& operator() (index_type i, index_type j) {
            MATVEC_ASSERT(i >= 1 && i <= iNumRows);
            MATVEC_ASSERT(j >= 1 && j <= iNumCols);
            return rgData[(i - 1) * iNumCols + (j - 1)];
        }

        const T& operator() (index_type i, index_type j) const {
            MATVEC_ASSERT(i >= 1 && i <= iNumRows);
            MATVEC_ASSERT(j >= 1 && j <= iNumCols);
            return rgData[(i - 1) * iNumCols + (j - 1)];
        }

        void Resize(index_type iRows, index_type iCols) {
            if (iRows != iNumRows || iCols != iNumCols) {
                rgData.Resize(iRows * iCols);
                iNumRows = iRows;
                iNumCols = iCols;
            }
        }

        const T* pGetFirstElem() const {
            return rgData.pGetFirstElem();
        }

        const T* pGetLastElem() const {
            return rgData.pGetLastElem();
        }        
    private:
        integer iNumRows, iNumCols;
        VectorData<T, DYNAMIC_SIZE> rgData;
    };

    template <typename T, index_type N_rows>
    class MatrixData<T, N_rows, DYNAMIC_SIZE>
    {
    public:
        MatrixData(index_type iNumRows, index_type iNumCols, bool bZeroInit)
            :iNumCols(iNumCols),
             rgData(N_rows * iNumCols, bZeroInit) {
            MATVEC_ASSERT(iNumRows == N_rows);
        }

        index_type iGetNumRows() const { return N_rows; }
        index_type iGetNumCols() const { return iNumCols; }       
        
        T& operator() (index_type i, index_type j) {
            MATVEC_ASSERT(i >= 1 && i <= N_rows);
            MATVEC_ASSERT(j >= 1 && j <= iNumCols);
            return rgData[(i - 1) * iNumCols + (j - 1)];
        }

        const T& operator() (index_type i, index_type j) const {
            MATVEC_ASSERT(i >= 1 && i <= N_rows);
            MATVEC_ASSERT(j >= 1 && j <= iNumCols);
            return rgData[(i - 1) * iNumCols + (j - 1)];
        }

        void Resize(index_type iRows, index_type iCols) {
            MATVEC_ASSERT(iRows == N_rows);
            
            if (iCols != iNumCols) {
                rgData.Resize(iRows * iCols);
                iNumCols = iCols;
            }
        }

        const T* pGetFirstElem() const {
            return rgData.pGetFirstElem();
        }

        const T* pGetLastElem() const {
            return rgData.pGetLastElem();
        }        
    private:
        integer iNumCols;
        VectorData<T, DYNAMIC_SIZE> rgData;
    };

    template <typename T, index_type N_cols>
    class MatrixData<T, DYNAMIC_SIZE, N_cols>
    {
    public:
        MatrixData(index_type iNumRows, index_type iNumCols, bool bZeroInit)
            :iNumRows(iNumRows),
             rgData(iNumRows * N_cols, bZeroInit) {
            MATVEC_ASSERT(iNumCols == N_cols);
        }

        index_type iGetNumRows() const { return iNumRows; }
        index_type iGetNumCols() const { return N_cols; }       
        
        T& operator() (index_type i, index_type j) {
            MATVEC_ASSERT(i >= 1 && i <= iNumRows);
            MATVEC_ASSERT(j >= 1 && j <= N_cols);
            return rgData[(i - 1) * N_cols + (j - 1)];
        }

        const T& operator() (index_type i, index_type j) const {
            MATVEC_ASSERT(i >= 1 && i <= iNumRows);
            MATVEC_ASSERT(j >= 1 && j <= N_cols);
            return rgData[(i - 1) * N_cols + (j - 1)];
        }

        void Resize(index_type iRows, index_type iCols) {
            MATVEC_ASSERT(iCols == N_cols);
            
            if (iRows != iNumRows) {
                rgData.Resize(iRows * iCols);
                iNumRows = iRows;
            }
        }

        const T* pGetFirstElem() const {
            return rgData.pGetFirstElem();
        }

        const T* pGetLastElem() const {
            return rgData.pGetLastElem();
        }        
    private:
        integer iNumRows;
        VectorData<T, DYNAMIC_SIZE> rgData;
    };
    
    template <typename T, index_type N_rows, index_type N_cols>
    class Matrix {    
    public:
        static const index_type iNumRows = N_rows;
        static const index_type iNumCols = N_cols;
        static const index_type iInitNumRows = iNumRows == DYNAMIC_SIZE ? 0 : iNumRows;
        static const index_type iInitNumCols = iNumCols == DYNAMIC_SIZE ? 0 : iNumCols;
        typedef typename ScalarTypeTraits<T>::ScalarType ScalarType;
        typedef typename ScalarTypeTraits<T>::DirectExpressionType ExpressionType;

        typedef VectorExpression<SliceVector<T, N_cols, 1>, N_cols> RowVectorType;
        typedef VectorExpression<SliceVector<T, N_rows, N_cols>, N_rows> ColumnVectorType;

        Matrix()
            :rgMat(iInitNumRows, iInitNumCols, true) {
        }

        Matrix(index_type iRows, index_type iCols)
            :rgMat(iRows, iCols, true) {
        }
        
        Matrix(const T& A11, const T& A21,
               const T& A12, const T& A22)
            :rgMat(iNumRows, iNumCols, false) {
	     static_assert(iNumRows == 2);
	     static_assert(iNumCols == 2);

            (*this)(1, 1) = A11;
            (*this)(2, 1) = A21;
            (*this)(1, 2) = A12;
            (*this)(2, 2) = A22;
        }

        Matrix(const T& A11, const T& A21, const T& A31,
               const T& A12, const T& A22, const T& A32,
               const T& A13, const T& A23, const T& A33)
            :rgMat(iNumRows, iNumCols, false) {
	     static_assert(iNumRows == 3);
	     static_assert(iNumCols == 3);

            (*this)(1, 1) = A11;
            (*this)(2, 1) = A21;
            (*this)(3, 1) = A31;
            (*this)(1, 2) = A12;
            (*this)(2, 2) = A22;
            (*this)(3, 2) = A32;
            (*this)(1, 3) = A13;
            (*this)(2, 3) = A23;
            (*this)(3, 3) = A33;
        }
        
        Matrix(const Matrix& A)
            :rgMat(A.rgMat) {
        }
        
        explicit inline Matrix(const Mat3x3& A);

        template <typename InitClass>
        explicit Matrix(const MatrixInit<InitClass, T, N_rows, N_cols>& func)
            :rgMat(func.iGetNumRows(), func.iGetNumCols(), false) {
            func.Initialize(*this);
        }

        template <typename Expression>
        Matrix(const MatrixExpression<Expression, N_rows, N_cols>& A)
            :rgMat(A.iGetNumRows(), A.iGetNumCols(), false) {
            // No aliases are possible because the object did not exist before
            using namespace MatVecHelp;
            ApplyMatrixFuncNoAlias(A, Assign());
        }

        template <typename T2>
        explicit Matrix(const Matrix<T2, N_rows, N_cols>& A)
            :rgMat(A.iGetNumRows(), A.iGetNumCols(), false) {
        
            MATVEC_ASSERT(A.iGetNumRows() == iGetNumRows());
            MATVEC_ASSERT(A.iGetNumCols() == iGetNumCols());
            
            for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
                for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
                    grad::Convert((*this)(i, j), A(i, j));
                }
            }
        }
        
        template <typename T2>
        Matrix(const Matrix<T2, N_rows, N_cols>& A, LocalDofMap* pDofMap)
            :rgMat(A.iGetNumRows(), A.iGetNumCols(), false) {
            Copy(A, pDofMap);
        }

        template <typename T2>
        void Copy(const Matrix<T2, N_rows, N_cols>& A, LocalDofMap* pDofMap) {
            rgMat.Resize(A.iGetNumRows(), A.iGetNumCols());
    
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

        void Resize(index_type iRows, index_type iCols) {
            rgMat.Resize(iRows, iCols);
        }

        template <typename Expression>
        Matrix& operator=(const MatrixExpression<Expression, N_rows, N_cols>& A) {

            rgMat.Resize(A.iGetNumRows(), A.iGetNumCols());
    
            using namespace MatVecHelp;

            ApplyMatrixFunc<Assign>(A);

            return *this;
        }

        template <typename T_Rhs>
        Matrix& operator+=(const Matrix<T_Rhs, N_rows, N_cols>& A) {
            MATVEC_ASSERT(A.iGetNumRows() == iGetNumRows());
            MATVEC_ASSERT(A.iGetNumCols() == iGetNumCols());
            using namespace MatVecHelp;

            ApplyMatrixFunc<Add>(MatrixExpression<MatrixDirectExpr<Matrix<T_Rhs, N_rows, N_cols> >, N_rows, N_cols>(A));

            return *this;
        }

        template <typename T_Rhs>
        Matrix& operator-=(const Matrix<T_Rhs, N_rows, N_cols>& A) {
            MATVEC_ASSERT(A.iGetNumRows() == iGetNumRows());
            MATVEC_ASSERT(A.iGetNumCols() == iGetNumCols());
    
		using namespace MatVecHelp;

		ApplyMatrixFunc<Sub>(MatrixExpression<MatrixDirectExpr<Matrix<T_Rhs, N_rows, N_cols> >, N_rows, N_cols>(A));

    	return *this;
    }

    template <typename Expression>
    Matrix& operator+=(const MatrixExpression<Expression, N_rows, N_cols>& A) {
            MATVEC_ASSERT(A.iGetNumRows() == iGetNumRows());
            MATVEC_ASSERT(A.iGetNumCols() == iGetNumCols());
		using namespace MatVecHelp;

		ApplyMatrixFunc<Add>(A);

    	return *this;
    }

    template <typename Expression>
    Matrix& operator-=(const MatrixExpression<Expression, N_rows, N_cols>& A) {
            MATVEC_ASSERT(A.iGetNumRows() == iGetNumRows());
            MATVEC_ASSERT(A.iGetNumCols() == iGetNumCols());
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
            MATVEC_ASSERT(iRow >= 1);
            MATVEC_ASSERT(iRow <= iGetNumRows());
            MATVEC_ASSERT(iCol >= 1);
            MATVEC_ASSERT(iCol <= iGetNumCols());
            return rgMat(iRow, iCol);
    }
    
    ScalarType& operator()(index_type iRow, index_type iCol) {
            MATVEC_ASSERT(iRow >= 1);
            MATVEC_ASSERT(iRow <= iGetNumRows());
            MATVEC_ASSERT(iCol >= 1);
            MATVEC_ASSERT(iCol <= iGetNumCols());
            return rgMat(iRow, iCol);
    }    
    
    RowVectorType
    GetRow(index_type iRow) const {
    	MATVEC_ASSERT(iRow >= 1);
            MATVEC_ASSERT(iRow <= iGetNumRows());
            return RowVectorType(&rgMat(iRow, 1), iGetNumCols(), 1);
    }

    ColumnVectorType
    GetCol(index_type iCol) const {
    	MATVEC_ASSERT(iCol >= 1);
            MATVEC_ASSERT(iCol <= iGetNumCols());
            return ColumnVectorType(&rgMat(1, iCol), iGetNumRows(), iGetNumCols());
    }

        index_type iGetNumRows() const { return rgMat.iGetNumRows(); }
        index_type iGetNumCols() const { return rgMat.iGetNumCols(); }
        const ScalarType* pGetMat() const { return pGetFirstElem(); }

    template <typename ScalarType2>
    bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
            return bArrayOverlap(pGetFirstElem(),
                                 pGetLastElem(),
    						 pFirst,
    						 pLast);
    }

private:
    const ScalarType* pGetFirstElem() const {
            return rgMat.pGetFirstElem();
    }

    const ScalarType* pGetLastElem() const {
            return rgMat.pGetLastElem();
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
            MATVEC_ASSERT(!A.bHaveReferenceTo(pGetFirstElem(), pGetLastElem()));
            
            rgMat.Resize(A.iGetNumRows(), A.iGetNumCols());
            
		MATVEC_ASSERT(A.iGetNumRows() == iGetNumRows());
		MATVEC_ASSERT(A.iGetNumCols() == iGetNumCols());

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
        MatrixData<ScalarType, iNumRows, iNumCols> rgMat;
};

template <typename T, index_type N_rows, index_type N_cols>
    inline MatrixExpression<MatrixDirectExpr<Matrix<T, N_rows, N_cols> >, N_rows, N_cols>
    Direct(const Matrix<T, N_rows, N_cols>& A) {
	return MatrixExpression<MatrixDirectExpr<Matrix<T, N_rows, N_cols> >, N_rows, N_cols>(A);
    }

    inline MatrixExpression<Mat3x3DirectExpr, 3, 3>
    Direct(const Mat3x3& A) {
	return MatrixExpression<Mat3x3DirectExpr, 3, 3>(A);
    }

    inline MatrixExpression<Mat3xNDirectExpr, 3, DYNAMIC_SIZE>
    Direct(const Mat3xN& A) {
        return MatrixExpression<Mat3xNDirectExpr, 3, DYNAMIC_SIZE>(A);
    }

    inline MatrixExpression<MatNxNDirectExpr, DYNAMIC_SIZE, DYNAMIC_SIZE>
    Direct(const MatNxN& A) {
        return MatrixExpression<MatNxNDirectExpr, DYNAMIC_SIZE, DYNAMIC_SIZE>(A);
    }
    
    template <typename T, index_type N_rows, index_type N_cols>
inline MatrixExpression<TransposedMatrix<MatrixDirectExpr<Matrix<T, N_rows, N_cols> > >, N_cols, N_rows>
Transpose(const Matrix<T, N_rows, N_cols>& A) {
	return MatrixExpression<TransposedMatrix<MatrixDirectExpr<Matrix<T, N_rows, N_cols> > >, N_cols, N_rows>(Direct(A));
}

    inline MatrixExpression<TransposedMatrix<Mat3xNDirectExpr>, DYNAMIC_SIZE, 3>
    Transpose(const Mat3xN& A) {
	return MatrixExpression<TransposedMatrix<Mat3xNDirectExpr>, DYNAMIC_SIZE, 3>(Direct(A));
}

    inline MatrixExpression<TransposedMatrix<MatNxNDirectExpr>, DYNAMIC_SIZE, DYNAMIC_SIZE>
    Transpose(const MatNxN& A) {
	return MatrixExpression<TransposedMatrix<MatNxNDirectExpr>, DYNAMIC_SIZE, DYNAMIC_SIZE>(Direct(A));
}

    template <typename MatrixExpr, index_type N_rows, index_type N_cols>
    inline MatrixExpression<TransposedMatrix<MatrixExpr>, N_cols, N_rows>
    Transpose(const MatrixExpression<MatrixExpr, N_rows, N_cols>& A) {
	return MatrixExpression<TransposedMatrix<MatrixExpr>, N_cols, N_rows>(A);
}

template <typename T, index_type N_rows, index_type N_cols>
inline MatrixExpression<MatrixDirectExpr<Matrix<T, N_rows, N_cols>, true>, N_rows, N_cols>
Alias(const Matrix<T, N_rows, N_cols>& A) {
	return MatrixExpression<MatrixDirectExpr<Matrix<T, N_rows, N_cols>, true>, N_rows, N_cols>(A);
}

template <typename T, index_type N_rows, index_type N_cols>
    inline Matrix<T, N_rows, N_cols>::Matrix(const Mat3x3& A)
        :rgMat(N_rows, N_cols, false) {
	using namespace MatVecHelp;

	ApplyMatrixFunc<Assign>(Direct(A));
}

template <typename T, index_type N_rows, index_type N_cols>
inline Matrix<T, N_rows, N_cols>&
Matrix<T, N_rows, N_cols>::operator=(const Mat3x3& A) {
     static_assert(N_rows == 3 || N_rows == DYNAMIC_SIZE);
     static_assert(N_cols == 3 || N_cols == DYNAMIC_SIZE);
        
        rgMat.Resize(3, 3);

	using namespace MatVecHelp;

	ApplyMatrixFunc<Assign>(Direct(A));

	return *this;
}

template <typename T, index_type N_rows>
class Vector {
public:
    static const index_type iNumRows = N_rows;
        static const index_type iInitNumRows = iNumRows == DYNAMIC_SIZE ? 0 : iNumRows;
    typedef typename ScalarTypeTraits<T>::ScalarType ScalarType;
    typedef typename ScalarTypeTraits<T>::DirectExpressionType ExpressionType;

        Vector()
            :rgVec(iInitNumRows, true) {
    }

        explicit Vector(index_type N)
            :rgVec(N, true) {
        }
    
        Vector(const Vector& oVec)
            :rgVec(oVec.rgVec) {
        }
    
        Vector(const T& v1, const T& v2)
            :rgVec(iInitNumRows, false) {
	     static_assert(iNumRows == 2);

    	(*this)(1) = v1;
    	(*this)(2) = v2;
    }

    template <typename Expr1, typename Expr2>
    Vector(const GradientExpression<Expr1>& v1, const GradientExpression<Expr2>& v2)
            :rgVec(iNumRows, false) {
	 static_assert(iNumRows == 2);

    	(*this)(1) = v1;
    	(*this)(2) = v2;
    }

        Vector(const T& v1, const T& v2, const T& v3)
            :rgVec(iNumRows, false) {
	     static_assert(iNumRows == 3);

    	(*this)(1) = v1;
    	(*this)(2) = v2;
    	(*this)(3) = v3;
    }

    template <typename Expr1, typename Expr2, typename Expr3>
    Vector(const GradientExpression<Expr1>& v1, const GradientExpression<Expr2>& v2, const GradientExpression<Expr3>& v3)
            :rgVec(iNumRows, false) {
	 static_assert(iNumRows == 3);

    	(*this)(1) = v1;
    	(*this)(2) = v2;
    	(*this)(3) = v3;
    }

    explicit inline Vector(const Vec3& v);

    template <typename Expression>
        Vector(const VectorExpression<Expression, N_rows>& f)
            :rgVec(f.iGetNumRows(), false) {
    	using namespace MatVecHelp;

    	ApplyMatrixFuncNoAlias(f, Assign());
    }
    
        template <typename InitClass>
        explicit Vector(const VectorInit<InitClass, T, N_rows>& func)
            :rgVec(func.iGetNumRows(), false) {
            func.Initialize(*this);
        }
        
    template <typename T2>
        explicit Vector(const Vector<T2, N_rows>& v)
            :rgVec(v.iGetNumRows(), false) {
        
    	MATVEC_ASSERT(v.iGetNumRows() == iGetNumRows());

    	for (index_type i = 1; i <= v.iGetNumRows(); ++i) {
    		grad::Convert((*this)(i), v(i));
    	}
    }

    template <typename T2>
        Vector(const Vector<T2, N_rows>& v, LocalDofMap* pDofMap)
            :rgVec(v.iGetNumRows(), false) {
        
    	Copy(v, pDofMap);
    }

    template <typename T2>
    void Copy(const Vector<T2, N_rows>& v, LocalDofMap* pDofMap) {
            rgVec.Resize(v.iGetNumRows());
    
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

        void Resize(index_type iNumRows) {
            rgVec.Resize(iNumRows);
        }

    template <typename Expression>
    Vector& operator=(const VectorExpression<Expression, N_rows>& f) {
            rgVec.Resize(f.iGetNumRows());
    
    	using namespace MatVecHelp;

    	ApplyMatrixFunc<Assign>(f);

        return *this;
    }

    template <typename T_Rhs>
    Vector& operator+=(const Vector<T_Rhs, N_rows>& v) {
            MATVEC_ASSERT(v.iGetNumRows() == iGetNumRows());
            using namespace MatVecHelp;

            ApplyMatrixFunc<Add>(VectorExpression<VectorDirectExpr<Vector<T_Rhs, N_rows> >, N_rows>(v));

            return *this;
        }

        Vector& operator+=(const Vec3& v) {
            MATVEC_ASSERT(iGetNumRows() == 3);
            static_assert(iNumRows == 3 || iNumRows == DYNAMIC_SIZE);
            
    	using namespace MatVecHelp;

            ApplyMatrixFunc<Add>(VectorExpression<Vec3DirectExpr, N_rows>(v));

    	return *this;
    }

    template <typename T_Rhs>
    Vector& operator-=(const Vector<T_Rhs, N_rows>& v) {
            MATVEC_ASSERT(v.iGetNumRows() == iGetNumRows());
    	using namespace MatVecHelp;

            ApplyMatrixFunc<Sub>(VectorExpression<VectorDirectExpr<Vector<T_Rhs, N_rows> >, N_rows>(v));

    	return *this;
    }

    template <typename Expression>
    Vector& operator+=(const VectorExpression<Expression, N_rows>& f) {
            MATVEC_ASSERT(f.iGetNumRows() == iGetNumRows());
    	using namespace MatVecHelp;

    	ApplyMatrixFunc<Add>(f);

    	return *this;
    }

    template <typename Expression>
    Vector& operator-=(const VectorExpression<Expression, N_rows>& f) {
            MATVEC_ASSERT(f.iGetNumRows() == iGetNumRows());
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
            MATVEC_ASSERT(iRow < iGetNumRows());
        return rgVec[iRow];
    }
    
    ScalarType& operator()(index_type iRow) {
    	--iRow; // Row index is 1-based for compatibility reasons
        MATVEC_ASSERT(iRow >= 0);
            MATVEC_ASSERT(iRow < iGetNumRows());
        return rgVec[iRow];
    }
    
        index_type iGetNumRows() const {
            const index_type iRows = rgVec.iGetNumRows();
    
            MATVEC_ASSERT((N_rows == DYNAMIC_SIZE) || (iRows == N_rows));
        
            return iRows;
        }
    
        ScalarType* pGetVec() { return &rgVec[0]; }
        const ScalarType* pGetVec() const { return pGetFirstElem(); }

    template <typename ScalarType2>
    bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
            return bArrayOverlap(pGetFirstElem(),
                                 pGetLastElem(),
    						 pFirst,
    						 pLast);
    }

private:
    const ScalarType* pGetFirstElem() const {
            return rgVec.pGetFirstElem();
    }

    const ScalarType* pGetLastElem() const {
            return rgVec.pGetLastElem();
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
            MATVEC_ASSERT((N_rows == A.iGetNumRows()) || (N_rows == DYNAMIC_SIZE && A.iGetNumRows() >= 0));
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
        VectorData<ScalarType, N_rows> rgVec;
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
    Vector<T, N_rows>::Vector(const Vec3& v)
        :rgVec(iNumRows, false) {
	using namespace MatVecHelp;

        static_assert(iNumRows == 3);
        
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
        static const index_type iNumRows = 3;
        static const index_type iNumCols = 3;
        
	MatCrossInit(const VectorExpression<VectorExpr, 3>& v, doublereal d)
            :v(v), d(d) {

	}

	void Initialize(Matrix<T, 3, 3>& A) const {
		const index_type x = 1, y = 2, z = 3;
		/*
		 * A = [ 0, -z,  y;
		 *       z,  0, -x;
		 *      -y,  x,  0];
		 */
            A(x, x) = A(y, y) = A(z, z) = d;
		A(x, y) = -v(z);
		A(x, z) =  v(y);
		A(y, x) =  v(z);
		A(y, z) = -v(x);
		A(z, x) = -v(y);
		A(z, y) =  v(x);
	}

        index_type iGetNumRows() const { return iNumRows; }
        index_type iGetNumCols() const { return iNumCols; }
private:
	const VectorExpression<VectorExpr, 3> v;
        const doublereal d;
    };

    template <typename T, typename VectorExpr>
    class MatRInit {
    public:
        static const index_type iNumRows = 3;
        static const index_type iNumCols = 3;

        MatRInit(const VectorExpression<VectorExpr, 3>& g)
            :g(g) {
        }

        void Initialize(Matrix<T, 3, 3>& RDelta) const {
            const T d = 4. / (4. + Dot(g, g));

            const T tmp1 = -g(3) * g(3);
            const T tmp2 = -g(2) * g(2);
            const T tmp3 = -g(1) * g(1);
            const T tmp4 = g(1) * g(2) * 0.5;
            const T tmp5 = g(2) * g(3) * 0.5;
            const T tmp6 = g(1) * g(3) * 0.5;

            RDelta(1,1) = (tmp1 + tmp2) * d * 0.5 + 1;
            RDelta(1,2) = (tmp4 - g(3)) * d;
            RDelta(1,3) = (tmp6 + g(2)) * d;
            RDelta(2,1) = (g(3) + tmp4) * d;
            RDelta(2,2) = (tmp1 + tmp3) * d * 0.5 + 1.;
            RDelta(2,3) = (tmp5 - g(1)) * d;
            RDelta(3,1) = (tmp6 - g(2)) * d;
            RDelta(3,2) = (tmp5 + g(1)) * d;
            RDelta(3,3) = (tmp2 + tmp3) * d * 0.5 + 1.;
        }

        index_type iGetNumRows() const { return iNumRows; }
        index_type iGetNumCols() const { return iNumCols; }
        
    private:
	const VectorExpression<VectorExpr, 3> g;
};


    template <typename T, typename VectorExpr>
    class MatRotVecInit {
    public:
        static const index_type iNumRows = 3;
        static const index_type iNumCols = 3;

        MatRotVecInit(const VectorExpression<VectorExpr, 3>& Phi)
            :p(Phi) {
        }

        void Initialize(Matrix<T, 3, 3>& R) const {
            const index_type cid = RotCoeff::COEFF_B;
            using std::sqrt;
            using std::sin;
            using std::cos;
            
            T phip[10];
            T phi2(Dot(p, p));
            T pd(sqrt(phi2));
            T cf[RotCoeff::COEFF_B];
            index_type k, j;

            if (pd < RotCoeff::SerThrsh[cid-1]) {
                phip[0] = 1.;
                for (j = 1; j <= 9; j++) {
                    phip[j] = phip[j-1]*phi2;
                }
                for (k = 0; k < cid; k++) {
                    cf[k] = 0.;
                    for (j = 0; j < RotCoeff::SerTrunc[k]; j++) {
                        cf[k] += phip[j]/RotCoeff::SerCoeff[k][j];
                    }
                }
            } else {
                cf[0] = sin(pd) / pd;                 // a = sin(phi)/phi
                cf[1]=(1. - cos(pd)) / phi2;           // b = (1.-cos(phi))/phi2
            }

            R(1,1) = cf[1]*((-p(3)*p(3))-p(2)*p(2))+1.;
            R(1,2) = cf[1]*p(1)*p(2)-cf[0]*p(3);
            R(1,3) = cf[1]*p(1)*p(3)+cf[0]*p(2);
            R(2,1) = cf[0]*p(3)+cf[1]*p(1)*p(2);
            R(2,2) = cf[1]*((-p(3)*p(3))-p(1)*p(1))+1.;
            R(2,3) = cf[1]*p(2)*p(3)-cf[0]*p(1);
            R(3,1) = cf[1]*p(1)*p(3)-cf[0]*p(2);
            R(3,2) = cf[1]*p(2)*p(3)+cf[0]*p(1);
            R(3,3) = cf[1]*((-p(2)*p(2))-p(1)*p(1))+1.;
        }

        index_type iGetNumRows() const { return iNumRows; }
        index_type iGetNumCols() const { return iNumCols; }
        
    private:
	const VectorExpression<VectorExpr, 3> p;
    };
        
template <typename T, typename VectorExpr>
class MatCrossCrossInit {
public:
        static const index_type iNumRows = 3;
        static const index_type iNumCols = 3;
        
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

        index_type iGetNumRows() const { return iNumRows; }
        index_type iGetNumCols() const { return iNumCols; }
private:
	const VectorExpression<VectorExpr, 3> v;
};

template <typename T, typename VectorExpr>
class MatGInit {
public:
        static const index_type iNumRows = 3;
        static const index_type iNumCols = 3;
        
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

        index_type iGetNumRows() const { return iNumRows; }
        index_type iGetNumCols() const { return iNumCols; }
private:
	const VectorExpression<VectorExpr, 3> g;
};

template <typename T>
inline MatrixInit<MatCrossInit<T, VectorDirectExpr<Vector<T, 3> > >, T, 3, 3>
    MatCrossVec(const Vector<T, 3>& v, doublereal d=0.) {
	return MatrixInit<MatCrossInit<T, VectorDirectExpr<Vector<T, 3> > >, T, 3, 3>(Direct(v), d);
    }

    template <typename T>
    inline MatrixInit<MatCrossInit<T, Vec3DirectExpr>, T, 3, 3>
    MatCrossVec(const Vec3& v, doublereal d=0.) {
	return MatrixInit<MatCrossInit<T, Vec3DirectExpr>, T, 3, 3>(Direct(v), d);
}

template <typename VectorExpr>
inline MatrixInit<MatCrossInit<typename VectorExpr::ScalarType, VectorExpr>, typename VectorExpr::ScalarType, 3, 3>
    MatCrossVec(const VectorExpression<VectorExpr, 3>& v, doublereal d = 0.) {
	return MatrixInit<MatCrossInit<typename VectorExpr::ScalarType, VectorExpr>, typename VectorExpr::ScalarType, 3, 3>(v, d);
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

    template <typename T>
    inline MatrixInit<MatRInit<T, VectorDirectExpr<Vector<T, 3> > >, T, 3, 3>
    MatRVec(const Vector<T, 3>& g) {
	return MatrixInit<MatRInit<T, VectorDirectExpr<Vector<T, 3> > >, T, 3, 3>(Direct(g));
    }

    template <typename VectorExpr>
    inline MatrixInit<MatRInit<typename VectorExpr::ScalarType, VectorExpr>, typename VectorExpr::ScalarType, 3, 3>
    MatRVec(const VectorExpression<VectorExpr, 3>& g) {
	return MatrixInit<MatRInit<typename VectorExpr::ScalarType, VectorExpr>, typename VectorExpr::ScalarType, 3, 3>(g);
    }    
    
    template <typename T>
    inline MatrixInit<MatRotVecInit<T, VectorDirectExpr<Vector<T, 3> > >, T, 3, 3>
    MatRotVec(const Vector<T, 3>& Phi) {
	return MatrixInit<MatRotVecInit<T, VectorDirectExpr<Vector<T, 3> > >, T, 3, 3>(Direct(Phi));
    }

    template <typename VectorExpr>
    inline MatrixInit<MatRotVecInit<typename VectorExpr::ScalarType, VectorExpr>, typename VectorExpr::ScalarType, 3, 3>
    MatRotVec(const VectorExpression<VectorExpr, 3>& Phi) {
	return MatrixInit<MatRotVecInit<typename VectorExpr::ScalarType, VectorExpr>, typename VectorExpr::ScalarType, 3, 3>(Phi);
    }
            
    template <typename T, typename MatrixExpr>
    class VecRotInit {
    public:
        static const index_type iNumRows = 3;

        VecRotInit(const MatrixExpression<MatrixExpr, 3, 3>& R)
            :R(R) {
        }

        index_type iGetNumRows() const { return iNumRows; }
#if 0
        void Initialize(Vector<T, 3>& Phi) const {
            using std::sqrt;
            using std::fabs;
            using std::acos;
            /*
              This algorithm was taken from
              Rebecca M. Brannon and coauthors
              Computational Physics and Mechanics
              Sandia National Laboratories
              Albuquerque, NM 87185-0820
              http://www.me.unm.edu/~rmbrann/gobag.html
            */
            T cth = (R(1, 1) + R(2, 2) + R(3, 3) - 1.) * 0.5;
            T sth(0.);

            if (cth < 1.) {
                sth = sqrt(1. - cth * cth);
            }
            
            const scalar_func_type puny = 1e-12;

            if (fabs(sth) > puny) {
                const T angle = acos(cth);
                const T a1 = 0.5 * angle / sth;
                Phi(1) = (R(3, 2) - R(2, 3)) * a1;
                Phi(2) = (R(1, 3) - R(3, 1)) * a1;
                Phi(3) = (R(2, 1) - R(1, 2)) * a1;
            } else if (cth > 0.) {
                Phi(1) = R(3, 2);
                Phi(2) = R(1, 3);
                Phi(3) = R(2, 1);
            } else {
                const T angle = acos(cth);
                Matrix<T, 3, 3> scr = R;

                cth = 0.;

                index_type j = -1;

                for (index_type k = 1; k <= 3; ++k) {
                    scr(k, k) += 1.;
                    sth = Dot(scr.GetCol(k), scr.GetCol(k));

                    if (sth > cth) {
                        cth = sqrt(sth);
                        j = k;

                        for (index_type i = 1; i <= 3; ++i) {
                            scr(i, j) /= cth;
                        }
                    }
                }

                if (j < 1) {
                    MATVEC_ASSERT(false);
                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                for (index_type i = 1; i <= 3; ++i) {
                    Phi(i) = scr(i, j) * angle;
                }
            }
        }
#else
        void Initialize(Vector<T, 3>& unit) const {
            // Modified from Appendix 2.4 of
            //
            // author = {Marco Borri and Lorenzo Trainelli and Carlo L. Bottasso},
            // title = {On Representations and Parameterizations of Motion},
            // journal = {Multibody System Dynamics},
            // volume = {4},
            // pages = {129--193},
            // year = {2000}
            using std::atan2;
            using std::sqrt;
            
            const T cosphi = 0.5 * (R(1, 1) + R(2, 2) + R(3, 3) - 1.);
            
            if (cosphi > 0.) {
                unit(1) = 0.5*(R(3, 2) - R(2, 3));
		unit(2) = 0.5*(R(1, 3) - R(3, 1));
                unit(3) = 0.5*(R(2, 1) - R(1, 2));
                
		const T sinphi2 = Dot(unit, unit);
                T sinphi;
                
                if (sinphi2 != 0) {
                    sinphi = sqrt(sinphi2);
                } else {
                    sinphi = unit(1);
                }
                
		const T phi = atan2(sinphi, cosphi);
                T a;
                RotCo(phi, a);
		unit /= a;
            } else {
		// -1 <= cosphi <= 0
		Matrix<T, 3, 3> eet = (R + Transpose(R)) * 0.5;
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
                T sinphi(0.);
                for (index_type i = 1; i <= 3; ++i) {
                    sinphi -= Cross(unit, R.GetCol(i))(i) * 0.5;
                }

		unit *= atan2(sinphi, cosphi);
            }
        }
        
    private:        
        static void RotCo(const T& phi, T& cf) {
            // This algorithm is a simplified version of RotCo in RotCoeff.hc
            // from Marco Morandini  <morandini@aero.polimi.it>
            // and Teodoro Merlini  <merlini@aero.polimi.it>
            using std::sin;
            using std::cos;
            using std::sqrt;
            using std::fabs;
	
            T phip[10];
            T phi2(phi * phi);

            if (fabs(phi) < RotCoeff::SerThrsh[0]) {
                phip[0] = 1.;
            
                for (index_type j = 1; j <= 9; j++) {
                    phip[j] = phip[j - 1] * phi2;
                }

                cf = 0.;
            
                for (index_type j = 0; j < RotCoeff::SerTrunc[0]; j++) {
                    cf += phip[j] / RotCoeff::SerCoeff[0][j];
                }

                return;
            } 
	
            const T pd(sqrt(phi2));
            cf = sin(pd) / pd;                 // a = sin(phi)/phi
        };
#endif
        
    private:
        const MatrixExpression<MatrixExpr, 3, 3> R;
    };

    template <typename T>
    inline VectorInit<VecRotInit<T, MatrixDirectExpr<Matrix<T, 3, 3> > >, T, 3>
    VecRotMat(const Matrix<T, 3, 3>& R) {
	return VectorInit<VecRotInit<T, MatrixDirectExpr<Matrix<T, 3, 3> > >, T, 3>(Direct(R));
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

    template <typename T_Lhs, index_type N_rows>
    inline typename DotTraits<VectorDirectExpr<Vector<T_Lhs, N_rows> >, Vec3DirectExpr, N_rows, N_rows>::ExpressionType
    Dot(const Vector<T_Lhs, N_rows>& u, const Vec3& v) {
	return DotTraits<VectorDirectExpr<Vector<T_Lhs, N_rows> >, Vec3DirectExpr, N_rows, N_rows>::Dot(Direct(u), Direct(v));
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
        using std::sqrt;
	return sqrt(Dot(u, u));
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
     static_assert(VectorLhsExpr::iNumRows == 3);
     static_assert(VectorRhsExpr::iNumRows == 3);
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

     static_assert(iNumRows == VectorLhsExpr::iNumRows);
     static_assert(iNumRows == VectorRhsExpr::iNumRows);
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

    template <typename T_Lhs>
    VectorExpression<VectorCrossExpr<VectorDirectExpr<Vector<T_Lhs, 3> >, Vec3DirectExpr>, 3>
    inline Cross(const Vector<T_Lhs, 3>& u, const Vec3& v) {
	return VectorExpression<VectorCrossExpr<VectorDirectExpr<Vector<T_Lhs, 3> >, Vec3DirectExpr>, 3>(Direct(u), Direct(v));
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

    template <index_type N_rows, index_type N_cols, typename MatrixLhsExpr, typename VectorRhsExpr>
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
	 static_assert(N_rows == MatrixLhsExpr::iNumRows);
	 static_assert(N_cols == MatrixLhsExpr::iNumCols);
	 static_assert(N_cols == VectorRhsExpr::iNumRows);
	 static_assert(N_cols != DYNAMIC_SIZE);
        
        const MatrixLhsExpr A;
        const VectorRhsExpr x;
    };

    template <typename MatrixLhsExpr, typename VectorRhsExpr>
    class MatrixVectorProduct<DYNAMIC_SIZE, DYNAMIC_SIZE, MatrixLhsExpr, VectorRhsExpr> {
    public:
        static const bool bAlias = MatrixLhsExpr::bAlias || VectorRhsExpr::bAlias;
        static const index_type iNumRows = MatrixLhsExpr::iNumRows;
        typedef typename MatrixLhsExpr::ScalarType MatrixLhsScalarExpr;
        typedef typename VectorRhsExpr::ScalarType VectorRhsScalarExpr;

        typedef typename CommonScalarType<typename BasicScalarType<MatrixLhsScalarExpr>::ScalarType,
                                          typename BasicScalarType<VectorRhsScalarExpr>::ScalarType>::ScalarType ScalarType;
        typedef ScalarType ExpressionType;

        MatrixVectorProduct(const MatrixLhsExpr& A, const VectorRhsExpr& x)
            :A(A), x(x) {

        }
    
        ScalarType operator()(index_type i) const {
            MATVEC_ASSERT(i >= 1);
            MATVEC_ASSERT(i <= iGetNumRows());
            MATVEC_ASSERT(A.iGetNumCols() == x.iGetNumRows());
        
            ScalarType b_i(0);
        
            for (integer j = 1; j <= A.iGetNumCols(); ++j) {
                b_i += A(i, j) * x(j);
            }

            return b_i;
        }

        index_type iGetNumRows() const {
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

    template <index_type N_rows, typename MatrixLhsExpr, typename VectorRhsExpr>
    class MatrixVectorProduct<N_rows, DYNAMIC_SIZE, MatrixLhsExpr, VectorRhsExpr>: public MatrixVectorProduct<DYNAMIC_SIZE, DYNAMIC_SIZE, MatrixLhsExpr, VectorRhsExpr> {
    public:
        MatrixVectorProduct(const MatrixLhsExpr& A, const VectorRhsExpr& x)
            :MatrixVectorProduct<DYNAMIC_SIZE, DYNAMIC_SIZE, MatrixLhsExpr, VectorRhsExpr>(A, x) {

        }
    };

    template <index_type N_cols, typename MatrixLhsExpr, typename VectorRhsExpr>
    class MatrixVectorProduct<DYNAMIC_SIZE, N_cols, MatrixLhsExpr, VectorRhsExpr>: public MatrixVectorProduct<DYNAMIC_SIZE, DYNAMIC_SIZE, MatrixLhsExpr, VectorRhsExpr> {
    public:
        MatrixVectorProduct(const MatrixLhsExpr& A, const VectorRhsExpr& x)
            :MatrixVectorProduct<DYNAMIC_SIZE, DYNAMIC_SIZE, MatrixLhsExpr, VectorRhsExpr>(A, x) {

        }
    };
    
    template <index_type N_rows_Lhs, index_type N_cols_Lhs, index_type N_cols_Rhs, typename MatrixLhsExpr, typename MatrixRhsExpr>
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
	 static_assert(N_rows_Lhs == MatrixLhsExpr::iNumRows);
	 static_assert(N_cols_Lhs == MatrixLhsExpr::iNumCols);
	 static_assert(N_cols_Rhs == MatrixRhsExpr::iNumCols);
	 static_assert(MatrixLhsExpr::iNumCols == MatrixRhsExpr::iNumRows);
        
        const MatrixLhsExpr A;
        const MatrixRhsExpr B;
    };

    template <typename MatrixLhsExpr, typename MatrixRhsExpr>
    class MatrixMatrixProduct<DYNAMIC_SIZE, DYNAMIC_SIZE, DYNAMIC_SIZE, MatrixLhsExpr, MatrixRhsExpr> {
    public:
	static const bool bAlias = MatrixLhsExpr::bAlias || MatrixRhsExpr::bAlias;
	static const index_type iNumRows = DYNAMIC_SIZE;
	static const index_type iNumCols = DYNAMIC_SIZE;
	typedef typename MatrixLhsExpr::ScalarType MatrixLhsScalarExpr;
	typedef typename MatrixRhsExpr::ScalarType MatrixRhsScalarExpr;
	typedef typename CommonScalarType<typename BasicScalarType<MatrixLhsScalarExpr>::ScalarType,
                                          typename BasicScalarType<MatrixRhsScalarExpr>::ScalarType>::ScalarType ScalarType;
        typedef ScalarType ExpressionType;

        MatrixMatrixProduct(const MatrixLhsExpr& A, const MatrixRhsExpr& B)
            :A(A), B(B) {
            MATVEC_ASSERT(A.iGetNumCols() == B.iGetNumRows());
        }

        ExpressionType operator()(index_type i, index_type j) const {
            MATVEC_ASSERT(i >= 1);
            MATVEC_ASSERT(i <= iGetNumRows());
            MATVEC_ASSERT(j >= 1);
            MATVEC_ASSERT(j <= iGetNumCols());
            MATVEC_ASSERT(A.iGetNumCols() == B.iGetNumRows());
            
            ScalarType Cij(0.);

            for (index_type k = 1; k <= A.iGetNumCols(); ++k) {
                Cij += A(i, k) * B(k, j);
            }

            return Cij;
        }

        index_type iGetNumRows() const {
            return A.iGetNumRows();
        }

        index_type iGetNumCols() const {
            return B.iGetNumCols();
        }

        template <typename ScalarType2>
        bool bHaveReferenceTo(const ScalarType2* pFirst, const ScalarType2* pLast) const {
            return A.bHaveReferenceTo(pFirst, pLast) || B.bHaveReferenceTo(pFirst, pLast);
        }

    private:        
    const MatrixLhsExpr A;
    const MatrixRhsExpr B;
};

    template <index_type N_rows_lhs, typename MatrixLhsExpr, typename MatrixRhsExpr>
    class MatrixMatrixProduct<N_rows_lhs, DYNAMIC_SIZE, DYNAMIC_SIZE, MatrixLhsExpr, MatrixRhsExpr>
        :public MatrixMatrixProduct<DYNAMIC_SIZE, DYNAMIC_SIZE, DYNAMIC_SIZE, MatrixLhsExpr, MatrixRhsExpr> {
    public:
        MatrixMatrixProduct(const MatrixLhsExpr& A, const MatrixRhsExpr& B)
            :MatrixMatrixProduct<DYNAMIC_SIZE, DYNAMIC_SIZE, DYNAMIC_SIZE, MatrixLhsExpr, MatrixRhsExpr>(A, B) {
        }
    };

    template <index_type N_cols_lhs, typename MatrixLhsExpr, typename MatrixRhsExpr>
    class MatrixMatrixProduct<DYNAMIC_SIZE, N_cols_lhs, DYNAMIC_SIZE, MatrixLhsExpr, MatrixRhsExpr>
        :public MatrixMatrixProduct<DYNAMIC_SIZE, DYNAMIC_SIZE, DYNAMIC_SIZE, MatrixLhsExpr, MatrixRhsExpr> {
    public:
        MatrixMatrixProduct(const MatrixLhsExpr& A, const MatrixRhsExpr& B)
            :MatrixMatrixProduct<DYNAMIC_SIZE, DYNAMIC_SIZE, DYNAMIC_SIZE, MatrixLhsExpr, MatrixRhsExpr>(A, B) {
        }
    };

    template <index_type N_cols_rhs, typename MatrixLhsExpr, typename MatrixRhsExpr>
    class MatrixMatrixProduct<DYNAMIC_SIZE, DYNAMIC_SIZE, N_cols_rhs, MatrixLhsExpr, MatrixRhsExpr>
        :public MatrixMatrixProduct<DYNAMIC_SIZE, DYNAMIC_SIZE, DYNAMIC_SIZE, MatrixLhsExpr, MatrixRhsExpr> {
    public:
        MatrixMatrixProduct(const MatrixLhsExpr& A, const MatrixRhsExpr& B)
            :MatrixMatrixProduct<DYNAMIC_SIZE, DYNAMIC_SIZE, DYNAMIC_SIZE, MatrixLhsExpr, MatrixRhsExpr>(A, B) {
        }
    };
    
    template <index_type N_rows_lhs, index_type N_cols_lhs, typename MatrixLhsExpr, typename MatrixRhsExpr>
    class MatrixMatrixProduct<N_rows_lhs, N_cols_lhs, DYNAMIC_SIZE, MatrixLhsExpr, MatrixRhsExpr>
        :public MatrixMatrixProduct<DYNAMIC_SIZE, DYNAMIC_SIZE, DYNAMIC_SIZE, MatrixLhsExpr, MatrixRhsExpr> {
    public:
        MatrixMatrixProduct(const MatrixLhsExpr& A, const MatrixRhsExpr& B)
            :MatrixMatrixProduct<DYNAMIC_SIZE, DYNAMIC_SIZE, DYNAMIC_SIZE, MatrixLhsExpr, MatrixRhsExpr>(A, B) {
        }
    };

    template <index_type N_cols_lhs, index_type N_cols_rhs, typename MatrixLhsExpr, typename MatrixRhsExpr>
    class MatrixMatrixProduct<DYNAMIC_SIZE, N_cols_lhs, N_cols_rhs, MatrixLhsExpr, MatrixRhsExpr>
        :public MatrixMatrixProduct<DYNAMIC_SIZE, DYNAMIC_SIZE, DYNAMIC_SIZE, MatrixLhsExpr, MatrixRhsExpr> {
    public:
        MatrixMatrixProduct(const MatrixLhsExpr& A, const MatrixRhsExpr& B)
            :MatrixMatrixProduct<DYNAMIC_SIZE, DYNAMIC_SIZE, DYNAMIC_SIZE, MatrixLhsExpr, MatrixRhsExpr>(A, B) {
        }
    };
    
    template <index_type N_rows_lhs, index_type N_cols_rhs, typename MatrixLhsExpr, typename MatrixRhsExpr>
    class MatrixMatrixProduct<N_rows_lhs, DYNAMIC_SIZE, N_cols_rhs, MatrixLhsExpr, MatrixRhsExpr>
        :public MatrixMatrixProduct<DYNAMIC_SIZE, DYNAMIC_SIZE, DYNAMIC_SIZE, MatrixLhsExpr, MatrixRhsExpr> {
    public:
        MatrixMatrixProduct(const MatrixLhsExpr& A, const MatrixRhsExpr& B)
            :MatrixMatrixProduct<DYNAMIC_SIZE, DYNAMIC_SIZE, DYNAMIC_SIZE, MatrixLhsExpr, MatrixRhsExpr>(A, B) {
        }
    };
/****************************************************************************************************************
 * matrix vector product
 ****************************************************************************************************************/

template <typename MatrixLhsExpr, index_type N_rows, index_type N_cols, typename VectorRhsExpr>
    inline VectorExpression<MatrixVectorProduct<N_rows, N_cols, MatrixExpression<MatrixLhsExpr, N_rows, N_cols>, VectorExpression<VectorRhsExpr, N_cols> >, N_rows>
operator* (const MatrixExpression<MatrixLhsExpr, N_rows, N_cols>& A, const VectorExpression<VectorRhsExpr, N_cols>& x) {
	return VectorExpression<MatrixVectorProduct<N_rows, N_cols, MatrixExpression<MatrixLhsExpr, N_rows, N_cols>, VectorExpression<VectorRhsExpr, N_cols> >, N_rows>(A, x);
}

template <typename T, index_type N_rows, index_type N_cols, typename VectorRhsExpr>
    inline VectorExpression<MatrixVectorProduct<N_rows, N_cols, MatrixExpression<MatrixDirectExpr<Matrix<T, N_rows, N_cols> >, N_rows, N_cols>, VectorExpression<VectorRhsExpr, N_cols> >, N_rows>
operator* (const Matrix<T, N_rows, N_cols>& A, const VectorExpression<VectorRhsExpr, N_cols>& x) {
	return VectorExpression<MatrixVectorProduct<N_rows, N_cols, MatrixExpression<MatrixDirectExpr<Matrix<T, N_rows, N_cols> >, N_rows, N_cols>, VectorExpression<VectorRhsExpr, N_cols> >, N_rows>(Direct(A), x);
}

template <typename MatrixLhsExpr, index_type N_rows, index_type N_cols, typename T>
    inline VectorExpression<MatrixVectorProduct<N_rows, N_cols, MatrixExpression<MatrixLhsExpr, N_rows, N_cols>, VectorExpression<VectorDirectExpr<Vector<T, N_cols> >, N_cols> >, N_rows>
operator* (const MatrixExpression<MatrixLhsExpr, N_rows, N_cols>& A, const Vector<T, N_cols>& x) {
	return VectorExpression<MatrixVectorProduct<N_rows, N_cols, MatrixExpression<MatrixLhsExpr, N_rows, N_cols>, VectorExpression<VectorDirectExpr<Vector<T, N_cols> >, N_cols> >, N_rows>(A, Direct(x));
}

template <typename T1, typename T2, index_type N_rows, index_type N_cols>
    inline VectorExpression<MatrixVectorProduct<N_rows, N_cols, MatrixExpression<MatrixDirectExpr<Matrix<T1, N_rows, N_cols> >, N_rows, N_cols>, VectorExpression<VectorDirectExpr<Vector<T2, N_cols> >, N_cols> >, N_rows>
operator* (const Matrix<T1, N_rows, N_cols>& A, const Vector<T2, N_cols>& x) {
	return VectorExpression<MatrixVectorProduct<N_rows, N_cols, MatrixExpression<MatrixDirectExpr<Matrix<T1, N_rows, N_cols> >, N_rows, N_cols>, VectorExpression<VectorDirectExpr<Vector<T2, N_cols> >, N_cols> >, N_rows>(Direct(A), Direct(x));
}

template <typename T, index_type N_rows>
    inline VectorExpression<MatrixVectorProduct<N_rows, 3, MatrixExpression<MatrixDirectExpr<Matrix<T, N_rows, 3> >, N_rows, 3>, VectorExpression<Vec3DirectExpr, 3> >, N_rows>
operator* (const Matrix<T, N_rows, 3>& A, const Vec3& x) {
	return VectorExpression<MatrixVectorProduct<N_rows, 3, MatrixExpression<MatrixDirectExpr<Matrix<T, N_rows, 3> >, N_rows, 3>, VectorExpression<Vec3DirectExpr, 3> >, N_rows>(Direct(A), Direct(x));
}

template <typename T>
    inline VectorExpression<MatrixVectorProduct<3, 3, MatrixExpression<Mat3x3DirectExpr, 3, 3>, VectorExpression<VectorDirectExpr<Vector<T, 3> >, 3> >, 3>
operator* (const Mat3x3& A, const Vector<T, 3>& x) {
	return VectorExpression<MatrixVectorProduct<3, 3, MatrixExpression<Mat3x3DirectExpr, 3, 3>, VectorExpression<VectorDirectExpr<Vector<T, 3> >, 3> >, 3>(Direct(A), Direct(x));
    }

    template <typename T>
    inline VectorExpression<MatrixVectorProduct<DYNAMIC_SIZE, DYNAMIC_SIZE, MatrixExpression<MatNxNDirectExpr, DYNAMIC_SIZE, DYNAMIC_SIZE>, VectorExpression<VectorDirectExpr<Vector<T, DYNAMIC_SIZE> >, DYNAMIC_SIZE> >, DYNAMIC_SIZE> operator*(const MatNxN& A, const Vector<T, DYNAMIC_SIZE>& x) {
        MATVEC_ASSERT(A.iGetNumCols() == x.iGetNumRows());
        return VectorExpression<MatrixVectorProduct<DYNAMIC_SIZE, DYNAMIC_SIZE, MatrixExpression<MatNxNDirectExpr, DYNAMIC_SIZE, DYNAMIC_SIZE>, VectorExpression<VectorDirectExpr<Vector<T, DYNAMIC_SIZE> >, DYNAMIC_SIZE> >, DYNAMIC_SIZE>(Direct(A), Direct(x));
    }

    template <typename T>
    inline VectorExpression<MatrixVectorProduct<3, DYNAMIC_SIZE, MatrixExpression<Mat3xNDirectExpr, 3, DYNAMIC_SIZE>, VectorExpression<VectorDirectExpr<Vector<T, DYNAMIC_SIZE> >, DYNAMIC_SIZE> >, 3> operator*(const Mat3xN& A, const Vector<T, DYNAMIC_SIZE>& x) {
        MATVEC_ASSERT(A.iGetNumCols() == x.iGetNumRows());
        MATVEC_ASSERT(A.iGetNumRows() == 3);
        return VectorExpression<MatrixVectorProduct<3, DYNAMIC_SIZE, MatrixExpression<Mat3xNDirectExpr, 3, DYNAMIC_SIZE>, VectorExpression<VectorDirectExpr<Vector<T, DYNAMIC_SIZE> >, DYNAMIC_SIZE> >, 3> (Direct(A), Direct(x));
}

/****************************************************************************************************************
 * matrix matrix product
 ****************************************************************************************************************/

template <typename MatrixLhsExpr, index_type N_rows_Lhs, index_type N_cols_Lhs, typename MatrixRhsExpr, index_type N_cols_Rhs>
    inline MatrixExpression<MatrixMatrixProduct<N_rows_Lhs, N_cols_Lhs, N_cols_Rhs, MatrixExpression<MatrixLhsExpr, N_rows_Lhs, N_cols_Lhs>, MatrixExpression<MatrixRhsExpr, N_cols_Lhs, N_cols_Rhs> >, N_rows_Lhs, N_cols_Rhs>
operator* (const MatrixExpression<MatrixLhsExpr, N_rows_Lhs, N_cols_Lhs>& A, const MatrixExpression<MatrixRhsExpr, N_cols_Lhs, N_cols_Rhs>& B) {
	return MatrixExpression<MatrixMatrixProduct<N_rows_Lhs, N_cols_Lhs, N_cols_Rhs, MatrixExpression<MatrixLhsExpr, N_rows_Lhs, N_cols_Lhs>, MatrixExpression<MatrixRhsExpr, N_cols_Lhs, N_cols_Rhs> >, N_rows_Lhs, N_cols_Rhs>(A, B);
}

template <typename MatrixLhsExpr, index_type N_rows_Lhs, index_type N_cols_Lhs, typename T, index_type N_cols_Rhs>
    inline MatrixExpression<MatrixMatrixProduct<N_rows_Lhs, N_cols_Lhs, N_cols_Rhs, MatrixExpression<MatrixLhsExpr, N_rows_Lhs, N_cols_Lhs>, MatrixDirectExpr<Matrix<T, N_cols_Lhs, N_cols_Rhs> > >, N_rows_Lhs, N_cols_Rhs>
operator* (const MatrixExpression<MatrixLhsExpr, N_rows_Lhs, N_cols_Lhs>& A, const Matrix<T, N_cols_Lhs, N_cols_Rhs>& B) {
	return MatrixExpression<MatrixMatrixProduct<N_rows_Lhs, N_cols_Lhs, N_cols_Rhs, MatrixExpression<MatrixLhsExpr, N_rows_Lhs, N_cols_Lhs>, MatrixDirectExpr<Matrix<T, N_cols_Lhs, N_cols_Rhs> > >, N_rows_Lhs, N_cols_Rhs>(A, Direct(B));
}

template <typename T, index_type N_rows_Lhs, index_type N_cols_Lhs, typename MatrixRhsExpr, index_type N_cols_Rhs>
    inline MatrixExpression<MatrixMatrixProduct<N_rows_Lhs, N_cols_Lhs, N_cols_Rhs, MatrixDirectExpr<Matrix<T, N_rows_Lhs, N_cols_Lhs> >, MatrixExpression<MatrixRhsExpr, N_cols_Lhs, N_cols_Rhs> >, N_rows_Lhs, N_cols_Rhs>
operator* (const Matrix<T, N_rows_Lhs, N_cols_Lhs>& A, const MatrixExpression<MatrixRhsExpr, N_cols_Lhs, N_cols_Rhs>& B) {
	return MatrixExpression<MatrixMatrixProduct<N_rows_Lhs, N_cols_Lhs, N_cols_Rhs, MatrixDirectExpr<Matrix<T, N_rows_Lhs, N_cols_Lhs> >, MatrixExpression<MatrixRhsExpr, N_cols_Lhs, N_cols_Rhs> >, N_rows_Lhs, N_cols_Rhs>(Direct(A), B);
}

template <typename T_Lhs, index_type N_rows_Lhs, index_type N_cols_Lhs, typename T_Rhs, index_type N_cols_Rhs>
    inline MatrixExpression<MatrixMatrixProduct<N_rows_Lhs, N_cols_Lhs, N_cols_Rhs, MatrixDirectExpr<Matrix<T_Lhs, N_rows_Lhs, N_cols_Lhs> >, MatrixDirectExpr<Matrix<T_Rhs, N_cols_Lhs, N_cols_Rhs> > >, N_rows_Lhs, N_cols_Rhs>
operator* (const Matrix<T_Lhs, N_rows_Lhs, N_cols_Lhs>& A, const Matrix<T_Rhs, N_cols_Lhs, N_cols_Rhs>& B) {
	return MatrixExpression<MatrixMatrixProduct<N_rows_Lhs, N_cols_Lhs, N_cols_Rhs, MatrixDirectExpr<Matrix<T_Lhs, N_rows_Lhs, N_cols_Lhs> >, MatrixDirectExpr<Matrix<T_Rhs, N_cols_Lhs, N_cols_Rhs> > >, N_rows_Lhs, N_cols_Rhs>(Direct(A), Direct(B));
}

template <typename T, index_type N_rows_Lhs>
    inline MatrixExpression<MatrixMatrixProduct<N_rows_Lhs, 3, 3, MatrixDirectExpr<Matrix<T, N_rows_Lhs, 3> >, Mat3x3DirectExpr>, N_rows_Lhs, 3>
operator* (const Matrix<T, N_rows_Lhs, 3>& A, const Mat3x3& B) {
	return MatrixExpression<MatrixMatrixProduct<N_rows_Lhs, 3, 3, MatrixDirectExpr<Matrix<T, N_rows_Lhs, 3> >, Mat3x3DirectExpr>, N_rows_Lhs, 3>(Direct(A), Direct(B));
}

template <typename T, index_type N_cols_Rhs>
    inline MatrixExpression<MatrixMatrixProduct<3, 3, N_cols_Rhs, Mat3x3DirectExpr, MatrixDirectExpr<Matrix<T, 3, N_cols_Rhs> > >, 3, N_cols_Rhs>
operator* (const Mat3x3& A, const Matrix<T, 3, N_cols_Rhs>& B) {
	return MatrixExpression<MatrixMatrixProduct<3, 3, N_cols_Rhs, Mat3x3DirectExpr, MatrixDirectExpr<Matrix<T, 3, N_cols_Rhs> > >, 3, N_cols_Rhs>(Direct(A), Direct(B));
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
    template <typename T1, index_type N_rows>              \
    inline VectorExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename ScalarTypeTraits<T1>::DirectExpressionType, typename ScalarTypeTraits<doublereal>::DirectExpressionType>, VectorDirectExpr<Vector<T1, N_rows> >, Vec3DirectExpr>, N_rows> \
    FunctionName(const Vector<T1, N_rows>& u, const Vec3& v) { \
        return VectorExpression<ExpressionName<ScalarBinaryOperation<FunctionClass, typename ScalarTypeTraits<T1>::DirectExpressionType, typename ScalarTypeTraits<doublereal>::DirectExpressionType>, VectorDirectExpr<Vector<T1, N_rows> >, Vec3DirectExpr>, N_rows>(u, v); \
    }                                                                   \
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
