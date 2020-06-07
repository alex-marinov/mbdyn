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
        Copyright (C) 2013(-2017) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#include <algorithm>
#include <cassert>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <set>
#include <typeinfo>

#ifdef HAVE_BLITZ
#include <blitz/blitz.h>
#include <blitz/array.h>
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>
#include <blitz/matrix.h>
#endif

#ifdef HAVE_FEENABLEEXCEPT
#define _GNU_SOURCE 1
#include <fenv.h>
#endif

#ifndef MATVEC_DEBUG
	#define MATVEC_DEBUG 1
#endif

#ifndef GRADIENT_DEBUG
	#define GRADIENT_DEBUG 1
#endif

#include "ac/f2c.h"
#include "clock_time.h"
#include "matvec.h"
#include "matvec3.h"
#include "Rot.hh"
#include "matvecass.h"
#include "spmapmh.h"

using namespace grad;

int NLoops = 1;
int NLoopsAss = 1;
void tic();
void tic(doublereal& dTime);
doublereal toc();

doublereal random1() {
	return 2 * (doublereal(rand()) / RAND_MAX) - 1;
}

template <typename T>
bool bCompare(const T& a, const T& b, doublereal dTolRel = 0.) {
    assert(!std::isnan(a));
    assert(!std::isnan(b));
    doublereal dTolAbs = std::max<T>(1., std::max<T>(std::abs(a), std::abs(b))) * dTolRel;
    return std::abs(a - b) <= dTolAbs;
}

template <index_type N_SIZE>
bool bCompare(const Gradient<N_SIZE>& a, const Gradient<N_SIZE>& b, doublereal dTolRel = 0.) {
    assert(!std::isnan(a.dGetValue()));
    assert(!std::isnan(b.dGetValue()));
    doublereal dTolAbs = std::max(1., std::max(std::abs(a.dGetValue()), std::abs(b.dGetValue()))) * dTolRel;

    if (std::abs(a.dGetValue() - b.dGetValue()) > dTolAbs) {
        std::cerr << "a = " << a.dGetValue() << " != b = " << b.dGetValue() << std::endl;
        return false;
    }

    index_type iStart = std::min(a.iGetStartIndexLocal(), b.iGetStartIndexLocal());
    index_type iEnd = std::max(a.iGetEndIndexLocal(), b.iGetEndIndexLocal());

    for (index_type i = iStart; i < iEnd; ++i) {
        doublereal dTolAbs = std::max<scalar_deriv_type>(1., std::max<scalar_deriv_type>(std::abs(a.dGetDerivativeLocal(i)), std::abs(b.dGetDerivativeLocal(i)))) * dTolRel;
        assert(!std::isnan(a.dGetDerivativeLocal(i)));
        assert(!std::isnan(b.dGetDerivativeLocal(i)));
        if (std::abs(a.dGetDerivativeLocal(i) - b.dGetDerivativeLocal(i)) > dTolAbs) {
            std::cerr << "ad(" << i << ") = " << a.dGetDerivativeLocal(i) << " != bd(" << i << ") = " << b.dGetDerivativeLocal(i) << std::endl;
            return false;
        }
    }

    return true;
}

#define TEST_GROUP1
#define TEST_GROUP2

#ifdef TEST_GROUP1
void testScalarTypeTraits() {
	typedef ScalarBinaryExpressionTraits<FuncPlus, doublereal, doublereal, doublereal>::ExpressionType Expr001;
	typedef ScalarBinaryExpressionTraits<FuncPlus, Gradient<0>, Gradient<0>, Gradient<0> >::ExpressionType Expr002;
	typedef ScalarBinaryExpressionTraits<FuncPlus, Gradient<0>, doublereal, Gradient<0> >::ExpressionType Expr003;
	typedef ScalarBinaryExpressionTraits<FuncPlus, Gradient<0>, Gradient<0>, doublereal>::ExpressionType Expr004;
	typedef ScalarBinaryExpressionTraits<FuncPlus, Gradient<0>, Expr003, Expr004>::ExpressionType Expr005;
	typedef ScalarBinaryExpressionTraits<FuncPlus, Gradient<0>, Expr005, Expr001>::ExpressionType Expr006;
	typedef ScalarBinaryExpressionTraits<FuncPlus, Gradient<0>, Expr005, Expr006>::ExpressionType Expr007;
	typedef ScalarBinaryExpressionTraits<FuncPlus, doublereal, Expr001, Expr001>::ExpressionType Expr008;
	typedef ScalarBinaryExpressionTraits<FuncPlus, doublereal, Expr001, Expr008>::ExpressionType Expr009;
	typedef ScalarBinaryExpressionTraits<FuncPlus, doublereal, Expr008, Expr009>::ExpressionType Expr010;
	typedef ScalarBinaryExpressionTraits<FuncPlus, Gradient<0>, Expr010, Expr007>::ExpressionType Expr011;
	typedef BasicScalarType<Expr001>::ScalarType T001;
	typedef BasicScalarType<Expr002>::ScalarType T002;
	typedef BasicScalarType<Expr003>::ScalarType T003;
	typedef BasicScalarType<Expr004>::ScalarType T004;
	typedef BasicScalarType<Expr005>::ScalarType T005;
	typedef BasicScalarType<Expr006>::ScalarType T006;
	typedef BasicScalarType<Expr007>::ScalarType T007;
	typedef BasicScalarType<Expr008>::ScalarType T008;
	typedef BasicScalarType<Expr009>::ScalarType T009;
	typedef BasicScalarType<Expr010>::ScalarType T010;
	typedef BasicScalarType<Expr011>::ScalarType T011;
	typedef CrossTraits<VectorDirectExpr<Vector<Gradient<0>, 3> >, VectorDirectExpr<Vector<Gradient<0>, 3> > >::ExpressionType C001;
	typedef CrossTraits<VectorDirectExpr<Vector<doublereal, 3> >, VectorDirectExpr<Vector<doublereal, 3> > >::ExpressionType C002;
	std::cout << "Expr001=" << typeid(Expr001).name() << std::endl;
	std::cout << "Expr002=" << typeid(Expr002).name() << std::endl;
	std::cout << "Expr003=" << typeid(Expr003).name() << std::endl;
	std::cout << "Expr004=" << typeid(Expr004).name() << std::endl;
	std::cout << "Expr005=" << typeid(Expr005).name() << std::endl;
	std::cout << "Expr006=" << typeid(Expr006).name() << std::endl;
	std::cout << "Expr007=" << typeid(Expr007).name() << std::endl;
	std::cout << "T001=" << typeid(T001).name() << std::endl;
	std::cout << "T002=" << typeid(T002).name() << std::endl;
	std::cout << "T003=" << typeid(T003).name() << std::endl;
	std::cout << "T004=" << typeid(T004).name() << std::endl;
	std::cout << "T005=" << typeid(T005).name() << std::endl;
	std::cout << "T006=" << typeid(T006).name() << std::endl;
	std::cout << "T007=" << typeid(T007).name() << std::endl;
	std::cout << "T008=" << typeid(T008).name() << std::endl;
	std::cout << "T009=" << typeid(T009).name() << std::endl;
	std::cout << "T010=" << typeid(T010).name() << std::endl;
	std::cout << "T011=" << typeid(T011).name() << std::endl;
	assert(typeid(T001) == typeid(doublereal));
	assert(typeid(T002) == typeid(Gradient<0>));
	assert(typeid(T003) == typeid(Gradient<0>));
	assert(typeid(T004) == typeid(Gradient<0>));
	assert(typeid(T005) == typeid(Gradient<0>));
	assert(typeid(T006) == typeid(Gradient<0>));
	assert(typeid(T007) == typeid(Gradient<0>));
	assert(typeid(T008) == typeid(doublereal));
	assert(typeid(T009) == typeid(doublereal));
	assert(typeid(T010) == typeid(doublereal));
	assert(typeid(T011) == typeid(Gradient<0>));

	const index_type N_rows = 3;
	typedef VectorExpression<VectorVectorVectorBinaryExpr<ScalarBinaryOperation<FuncPlus, ScalarTypeTraits<Gradient<3> >::DirectExpressionType, ScalarTypeTraits<Gradient<3> >::DirectExpressionType>, VectorDirectExpr<Vector<Gradient<3>, N_rows> >, VectorDirectExpr<Vector<Gradient<3>, N_rows> > >, N_rows> V001;
	typedef SumTraits<V001, N_rows, 2> Sum003;
	Sum003 s003;
	std::cout << typeid(s003).name() << std::endl;

	std::cout << "sizeof(Matrix<Gradient<12>, 3, 3>())="
			<< sizeof(Matrix<Gradient<12>, 3, 3>) << std::endl;

	std::cout << "sizeof(Vector<Gradient<12> >())="
			<< sizeof(Vector<Gradient<12>, 3>) << std::endl;

	std::cout << "sizeof(Matrix<Gradient<12>, 3, 3>() * Vector<Gradient<12>, 3>())="
			<< sizeof(Matrix<Gradient<12>, 3, 3>() * Vector<Gradient<12>, 3>()) << std::endl;

	RangeVector<doublereal, 3> v3;
	std::cout << "v3.iGetMaxSize()=" << v3.iGetMaxSize() << std::endl;
	RangeVector<doublereal, 0> v0;
	std::cout << "v0.iGetMaxSize()=" << v0.iGetMaxSize() << std::endl;

	std::cout << typeid(C001).name() << std::endl;
	std::cout << typeid(C002).name() << std::endl;

	Gradient<0> g0;
	std::cout << "g0 max derivatives: " << g0.iGetMaxDerivatives() << std::endl;
}

#ifdef HAVE_BLITZ
template <typename T, index_type N>
void func(LocalDofMap* pDofMap, const blitz::TinyVector<T, N>& a, const blitz::TinyVector<T, N>& b, blitz::TinyVector<T, N>& c) {
	doublereal r1, r2, r3;
	r1 = random1();
	r2 = random1();
	r3 = random1();
#if 1
	const T d1 = blitz::dot(a, b) * (2. / 1.5);

	blitz::TinyVector<T, N> b_d1;

	for (int i = 0; i < N; ++i) {
		b_d1(i) = b(i) * d1;
	}

	c = (a * r1 + b * r2) / 2. - (a * r3 - b) * r1 + b_d1;
#else
	T d1 = blitz::dot(a, b);
	d1 *=  (2. / 1.5);
	c = (a * r1 + b * r2) / 2. - (a * r3 - b) * r1 + b * d1;
#endif
}
#endif

template <typename T, index_type N>
void func(const Vector<T, N>& a, const Vector<T, N>& b, Vector<T, N>& c) {
	doublereal r1, r2, r3;
	r1 = random1();
	r2 = random1();
	r3 = random1();
	const T d1 = Dot(a, b) * (2. / 1.5);
	c = (a * (r1 / 2) + b * (r2 / 2.)) - (a * (r3 * r1) - b * r1) + b * d1;
}

template <typename T>
void func(const T a[], const T b[], T c[], index_type N) {
	doublereal r1, r2, r3;
	r1 = random1();
	r2 = random1();
	r3 = random1();

	T d1 = T(0.);

	for (int i = 0; i < N; ++i) {
		d1 += a[i] * b[i] * (2. / 1.5);
	}

	for (int i = 0; i < N; ++i) {
		c[i] = (a[i] * r1 + b[i] * r2) / 2. - (a[i] * r3 - b[i]) * r1 + b[i] * d1;
	}
}

template <typename T, index_type N>
void func2(const Matrix<T, N, N>& A, const Vector<T, N>& b, const Vector<T, N>& c, Vector<T, N>& d, doublereal e, doublereal& dt) {
	doublereal start, stop;
	tic(start);
	d = (Transpose(A) * Vector<T, N>(A * Vector<T, N>((b - c) * e)));
	tic(stop);
	dt += stop - start;
}

template <typename T>
void func2(const T *A, const T *b, const T *c__, T *d__, const doublereal& e, doublereal& dt) {
	doublereal start, stop;
	tic(start);
    /* Parameter adjustments */
    --d__;
    --c__;
    --b;
    A -= 4;
    const T tmp1 = b[3] - c__[3];
    const T tmp2 = b[2] - c__[2];
    const T tmp3 = b[1] - c__[1];
    const T tmp4 = tmp3 * A[7];
    const T tmp5 = tmp1 * A[9];
    const T tmp6 = tmp1 * A[12];
    const T tmp7 = tmp1 * A[6];
    const T tmp8 = tmp2 * A[11];
    const T tmp10 = tmp3 * A[10];
    const T tmp11 = tmp6 + tmp8 + tmp10;
    const T tmp13 = tmp5 + tmp2 * A[8] + tmp4;
    const T tmp14 = tmp7 + tmp2 * A[5] + tmp3 * A[4];
    /* Function Body */
    d__[1] = (A[10] * tmp11 + A[7] * tmp13 + A[4] * tmp14) * e;
    d__[2] = (A[11] * tmp11 + A[8] * tmp13 + A[5] * tmp14) * e;
    d__[3] = (A[12] * tmp11 + A[9] * tmp13 + A[6] * tmp14) * e;

    tic(stop);
    dt += stop - start;
} /* func2_ */


#ifdef HAVE_FC_F77
/*
	  SUBROUTINE FUNC2AD(A, b, c, d, e)
      IMPLICIT NONE
      DOUBLE PRECISION A(3, 3), b(3), c(3), d(3), e
 */
extern "C" void __FC_DECL__(func2ad)(const doublereal A[3][3],
						 const doublereal b[3],
						 const doublereal c[3],
						 doublereal d[3],
						 const doublereal& e);

const integer nbdirsmax = 12;
/*
     SUBROUTINE FUNC2AD_DV(a, ad, b, bd, c, cd, d, dd, e, nbdirs)
      IMPLICIT INTEGER (n-n)
      PARAMETER (nbdirsmax = 12)
      DOUBLE PRECISION a(3, 3), b(3), c(3), d(3), e
      DOUBLE PRECISION ad(nbdirsmax, 3, 3), bd(nbdirsmax, 3), cd(
     +                 nbdirsmax, 3), dd(nbdirsmax, 3)
*/
extern "C" void __FC_DECL__(func2ad_dv)(const doublereal A[3][3],
							const doublereal Ad[3][3][nbdirsmax],
							const doublereal b[3],
							const doublereal bd[3][nbdirsmax],
							const doublereal c[3],
							const doublereal cd[3][nbdirsmax],
							doublereal d[3],
							doublereal dd[3][nbdirsmax],
							const doublereal& e,
							const integer& nbdirs);

inline void func2ad(const Matrix<doublereal, 3, 3>& A, const Vector<doublereal, 3>& b, const Vector<doublereal, 3>& c, Vector<doublereal, 3>& d, const doublereal& e, LocalDofMap*, doublereal& dt) {

	doublereal A_F[3][3];

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			A_F[j][i] = A(i + 1, j + 1);
		}
	}

	doublereal start, stop;

	tic(start);
	__FC_DECL__(func2ad)(A_F, b.pGetVec(), c.pGetVec(), d.pGetVec(), e);
	tic(stop);

	dt += stop - start;
}

template <index_type N_SIZE>
inline void func2ad(const Matrix<Gradient<N_SIZE>, 3, 3>& A, const Vector<Gradient<N_SIZE>, 3>& b, const Vector<Gradient<N_SIZE>, 3>& c, Vector<Gradient<N_SIZE>, 3>& d, const doublereal& e, LocalDofMap* pDofMap, doublereal& dt) {
     static_assert(N_SIZE <= index_type(nbdirsmax));

	doublereal A_F[3][3], Ad_F[3][3][nbdirsmax], b_F[3], bd_F[3][nbdirsmax], c_F[3], cd_F[3][nbdirsmax], d_F[3], dd_F[3][nbdirsmax];

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			const Gradient<N_SIZE>& A_ij = A(i + 1, j + 1);
			A_F[j][i] = A_ij.dGetValue();
			for (index_type k = 0; k < N_SIZE; ++k) {
				Ad_F[j][i][k] = A_ij.dGetDerivativeGlobal(k);
			}
		}

		b_F[i] = b(i + 1).dGetValue();
		c_F[i] = c(i + 1).dGetValue();

		for (index_type k = 0; k < N_SIZE; ++k) {
			bd_F[i][k] = b(i + 1).dGetDerivativeGlobal(k);
			cd_F[i][k] = c(i + 1).dGetDerivativeGlobal(k);
		}
	}

	doublereal start, stop;
	tic(start);
	__FC_DECL__(func2ad_dv)(A_F, Ad_F, b_F, bd_F, c_F, cd_F, d_F, dd_F, e, N_SIZE);
	tic(stop);
	dt += stop - start;

	for (int i = 0; i < 3; ++i) {
		d(i + 1).SetValuePreserve(d_F[i]);
		d(i + 1).DerivativeResizeReset(pDofMap, 0L, N_SIZE, MapVectorBase::GLOBAL, 0.);
		for (index_type k = 0; k < N_SIZE; ++k) {
			d(i + 1).SetDerivativeGlobal(k, dd_F[i][k]);
		}
	}
}
#endif // HAVE_FC_F77

template <typename T>
void func3(const Matrix<T, 3, 3>& R1, const Matrix<T, 3, 3>& R2, const Vector<T, 3>& a, const Vector<T, 3>& b, Vector<T, 3>& c, doublereal e) {
	c = (Cross(R1 * a, R2 * b) + Cross(R2 * a, R1 * b)) * e;
}

template
void func3<doublereal>(const Matrix<doublereal, 3, 3>& R1, const Matrix<doublereal, 3, 3>& R2, const Vector<doublereal, 3>& a, const Vector<doublereal, 3>& b, Vector<doublereal, 3>& c, doublereal e);

template
void func3<Gradient<0> >(const Matrix<Gradient<0> , 3, 3>& R1, const Matrix<Gradient<0> , 3, 3>& R2, const Vector<Gradient<0> , 3>& a, const Vector<Gradient<0> , 3>& b, Vector<Gradient<0> , 3>& c, doublereal e);


template <typename T, index_type N>
void callFunc(LocalDofMap* pDofMap, const Vector<T, N>& a, const Vector<T, N>& b, Vector<T, N>& c, Vector<T, N>& c1) {
	srand(0);
    tic();
    for (int i = 0; i < NLoops; ++i) {
    	func(a, b, c);
    }

    std::cerr << "Vector (" << typeid(T).name() << "): " << toc() << "s" << std::endl;

    srand(0);
    tic();
    for (int i = 0; i < NLoops; ++i) {
    	func(a.pGetVec(), b.pGetVec(), c1.pGetVec(), N);
    }

    std::cerr << "Array (" << typeid(T).name() << "): " << toc() << "s" << std::endl;

    std::cout << "a=" << std::endl << a << std::endl;
    std::cout << "b=" << std::endl << b << std::endl;
    std::cout << "c=" << std::endl << c << std::endl;
    std::cout << "c1=" << std::endl << c1 << std::endl;

    const doublereal dTol = sqrt(std::numeric_limits<doublereal>::epsilon());

    for (int i = 1; i <= N; ++i) {
    	assert(bCompare(c(i), c1(i), dTol));
    }
}

template <typename T>
void callFunc2(LocalDofMap* pDofMap, const Matrix<T, 3, 3>& A, const Vector<T, 3>& b, const Vector<T, 3>& c, Vector<T, 3>& d, Vector<T, 3>& d_C, Vector<T, 3>& d_F) {
	srand(0);
    doublereal dtMatVec = 0., dtC = 0.;
    for (int i = 0; i < NLoops; ++i) {
    	func2(A, b, c, d, random1(), dtMatVec);
    }

    std::cerr << "matvec (" << typeid(T).name() << "): " << dtMatVec << "s" << std::endl;

    srand(0);
    for (int i = 0; i < NLoops; ++i) {
    	func2(A.pGetMat(), b.pGetVec(), c.pGetVec(), d_C.pGetVec(), random1(), dtC);
    }

    std::cerr << "C (" << typeid(T).name() << "): " << dtC << "s" << std::endl;

#ifdef HAVE_FC_F77
    doublereal dtF = 0.;
    srand(0);
    for (int i = 0; i < NLoops; ++i) {
    	func2ad(A, b, c, d_F, random1(), pDofMap, dtF);
    }

    std::cerr << "F77 (" << typeid(T).name() << "): " << dtF << "s" << std::endl;

    std::cerr << "overhead matvec:" << dtMatVec / std::max(std::numeric_limits<doublereal>::epsilon(), dtF) << std::endl;
    std::cerr << "overhead C:" << dtC / std::max(std::numeric_limits<doublereal>::epsilon(), dtF) << std::endl;
#endif // HAVE_FC_F77

    std::cout << "A=" << std::endl << A << std::endl;
    std::cout << "b=" << std::endl << b << std::endl;
    std::cout << "c=" << std::endl << c << std::endl;
    std::cout << "d=" << std::endl << d << std::endl;
    std::cout << "d_C=" << std::endl << d_C << std::endl;
    std::cout << "d_F=" << std::endl << d_F << std::endl;

    const doublereal dTol = 10 * std::numeric_limits<scalar_deriv_type>::epsilon();

    for (int i = 1; i <= 3; ++i) {
    	assert(bCompare(d(i), d_C(i), dTol));
#ifdef HAVE_FC_F77
    	assert(bCompare(d(i), d_F(i), dTol));
#endif // HAVE_FC_F77
    }
}

template <index_type N>
void testMatVecGradient(doublereal c_C[N], doublereal cd_C[N][N]) {
    LocalDofMap dof;
    Vector<Gradient<N>, N> a, b;
    
    for (index_type i = 0; i < N; ++i) {
        a(i + 1).SetValuePreserve(100*(i + 1));
        b(i + 1).SetValuePreserve(1000*(i + 10));

        a(i + 1).DerivativeResizeReset(&dof, i, MapVectorBase::GLOBAL, 0.);
        a(i + 1).SetDerivativeGlobal(i, -1. - 10.);
        b(i + 1).DerivativeResizeReset(&dof, i, MapVectorBase::GLOBAL, 0.);
        b(i + 1).SetDerivativeGlobal(i, -2. - 20.);
    }
    
    Vector<Gradient<N>, N> c, c1;

    callFunc(&dof, a, b, c, c1);

    Vector<Gradient<N>, N> d = -c;

    std::cout << "c=" << std::endl << c << std::endl;
    std::cout << "d=" << std::endl << d << std::endl;

    const doublereal dTol = sqrt(std::numeric_limits<scalar_deriv_type>::epsilon());

    for (int i = 1; i <= N; ++i) {
    	assert(bCompare(c(i), c1(i), dTol));
    	assert(bCompare(d(i), Gradient<N>(-c(i))));
    }

    for (int i = 0; i < N; ++i) {
    	c_C[i] = c(i + 1).dGetValue();

    	for (int j = 0; j < N; ++j) {
    		cd_C[i][j] = c(i + 1).dGetDerivativeGlobal(j);
    	}
    }

    Gradient<N> s1;

    for (int i = 1; i <= N; ++i) {
    	s1 += c(i);
    }

    Gradient<N> s2 = Sum(c);

    Gradient<N> s4;

    for (int i = 1; i <= N; ++i) {
    	s4 += (a(i) * 2 - b(i) * 3) / 1.5 + c(i);
    }

    Gradient<N> s3 = Sum((a * 2 - b * 3) / 1.5 + c);

    Gradient<N> d1 = Dot(a, b);
    Gradient<N> d2;

    for (int i = 1; i <= N; ++i) {
    	d2 += a(i) * b(i);
    }

    Gradient<N> d3 = Dot(a + b, a - b);
    Gradient<N> d4;

    for (int i = 1; i <= N; ++i) {
    	d4 += (a(i) + b(i)) * (a(i) - b(i));
    }

    assert(d3.bIsEqual(d4));

    d3 = Dot(a + b, b);

    d4.SetValue(0.);

    for (int i = 1; i <= N; ++i) {
    	d4 += (a(i) + b(i)) * b(i);
    }

    assert(d3.bIsEqual(d4));

    d3 = Dot(b, a + b);

    std::cout << "s1=" << s1 << std::endl;
    std::cout << "s2=" << s2 << std::endl;
    std::cout << "s3=" << s3 << std::endl;
    std::cout << "s4=" << s4 << std::endl;
    std::cout << "d1=" << d1 << std::endl;
    std::cout << "d2=" << d2 << std::endl;
    std::cout << "d3=" << d3 << std::endl;
    std::cout << "d4=" << d4 << std::endl;

    assert(bCompare(d3, d4, dTol));
    assert(bCompare(s2, s1, dTol));
    assert(bCompare(s4, s3, dTol));
    assert(bCompare(d2, d1, dTol));

    const Vector<Gradient<N>, N> d5 = (d + c + c1) * (3.5 / (c(1) + c(2)) * c(1) / 2. * (c(1) - c(3)) / c(2));

    d += c + c1;

    d *= 3.5;

    d /= c(1) + c(2);

    d *= c(1);

    d /= 2.;

    d *= c(1) - c(3);

    d /= c(2);

    std::cout << "d=" << d << std::endl;
    std::cout << "d5=" << d5 << std::endl;

    for (int i = 1; i <= N; ++i) {
    	assert(bCompare(d(i), d5(i), dTol));
    }
}

#ifdef HAVE_BLITZ
template <index_type N>
void testMatVecGradientBlitz(doublereal c_C[N], doublereal cd_C[N][N]) {
    LocalDofMap dof;
    blitz::TinyVector<Gradient<N>, N> a, b;

    for (index_type i = 0; i < N; ++i) {
        a(i).SetValuePreserve(100*(i + 1));
        b(i).SetValuePreserve(1000*(i + 10));

        a(i).DerivativeResizeReset(&dof, i, MapVectorBase::GLOBAL, 0.);
        a(i).SetDerivativeGlobal(i, -1. - 10.);
        b(i).DerivativeResizeReset(&dof, i, MapVectorBase::GLOBAL, 0.);
        b(i).SetDerivativeGlobal(i, -2. - 20.);
    }

    blitz::TinyVector<Gradient<N>, N> c;

	srand(0);
    tic();
    for (int i = 0; i < NLoops; ++i) {
    	func(&dof, a, b, c);
    }

    for (int i = 0; i < N; ++i) {
    	c_C[i] = c(i).dGetValue();

    	for (int j = 0; j < N; ++j) {
    		cd_C[i][j] = c(i).dGetDerivativeGlobal(j);
    	}
    }

    std::cerr << "blitz vector (Gradient): " << toc() << "s" << std::endl;

    std::cout << "a=" << std::endl << a << std::endl;
    std::cout << "b=" << std::endl << b << std::endl;
    std::cout << "c=" << std::endl << c << std::endl;
}
#endif

template <index_type N>
void testMatVecDouble(doublereal c_C[N]) {
    Vector<doublereal, N> a, b;
    
    for (index_type i = 0; i < N; ++i) {
        a(i + 1) = 100*(i + 1);
        b(i + 1) = 1000*(i + 10);
    }
    
    Vector<doublereal, N> c, c1;

    callFunc(0, a, b, c, c1);

    Vector<doublereal, N> d = -c;

    const doublereal dTol = sqrt(std::numeric_limits<scalar_deriv_type>::epsilon());

    for (int i = 1; i <= N; ++i) {
    	assert(bCompare(c(i), c1(i), dTol));
    	assert(d(i) == -c(i));
    }

    for (int i = 0; i < N; ++i) {
    	c_C[i] = c(i + 1);
    }

    doublereal s1 = 0.;

    for (int i = 0; i < N; ++i) {
    	s1 += c(i + 1);
    }

    doublereal s2 = Sum(c);

    assert(bCompare(s2, s1, dTol));

    doublereal s3 = Sum((a * 2 - b * 3) / 1.5 + c);

    doublereal s4 = 0.;

    for (int i = 1; i <= N; ++i) {
    	s4 += (a(i) * 2 - b(i) * 3) / 1.5 + c(i);
    }

    assert(bCompare(s2, s1, dTol));
    assert(bCompare(s4, s3, dTol));

    doublereal d1 = Dot(a, b);
    doublereal d2 = 0.;

    for (int i = 1; i <= N; ++i) {
    	d2 += a(i) * b(i);
    }

    assert(bCompare(d2, d1, std::numeric_limits<scalar_deriv_type>::epsilon()));

    std::cout << "s1=" << s1 << std::endl;
    std::cout << "s2=" << s2 << std::endl;
    std::cout << "s3=" << s3 << std::endl;
    std::cout << "s4=" << s4 << std::endl;
    std::cout << "d1=" << d1 << std::endl;
    std::cout << "d2=" << d2 << std::endl;

    const Vector<doublereal, N> d5 = (d + c + c1) * (3.5 / (c(1) + c(2)) * c(1) / 2. * (c(1) - c(3)) / c(2));

    d += c + c1;

    d *= 3.5;

    d /= c(1) + c(2);

    d *= c(1);

    d /= 2.;

    d *= c(1) - c(3);

    d /= c(2);

    for (int i = 1; i <= N; ++i) {
    	assert(bCompare(d(i), d5(i), dTol));
    }

    std::cout << "d=" << d << std::endl;
    std::cout << "d5=" << d5 << std::endl;
}

template <index_type N_SIZE>
void testMatVecGradient2() {
    LocalDofMap dof;
    Matrix<Gradient<N_SIZE>, 3, 3> A;
    Vector<Gradient<N_SIZE>, 3> b, c;

    for (index_type i = 0; i < 3; ++i) {
        b(i + 1).SetValuePreserve(100*(i + 1));
        b(i + 1).DerivativeResizeReset(&dof, 0L, N_SIZE, MapVectorBase::GLOBAL, 0.);
        c(i + 1).SetValuePreserve(1000*(i + 10));
        c(i + 1).DerivativeResizeReset(&dof, 0L, N_SIZE, MapVectorBase::GLOBAL, 0.);

        for (index_type k = 0; k < N_SIZE; ++k) {
        	c(i + 1).SetDerivativeGlobal(k, 3.5 * (k + 1) + 3.2 * (i + 1));
        	b(i + 1).SetDerivativeGlobal(k, -17.4 * (k + 1) + 7.5 *(i + 1));
        }

        for (index_type j = 0; j < 3; ++j) {
            A(i + 1, j + 1).SetValuePreserve(100*(i + 1) + j + 1);
            A(i + 1, j + 1).DerivativeResizeReset(&dof, 0L, N_SIZE, MapVectorBase::GLOBAL, 0.);
            for (index_type k = 0; k < N_SIZE; ++k) {
            	A(i + 1, j + 1).SetDerivativeGlobal(k, -9.2 * (k + 1.) + 10.5 * (i + 1) + 3.2 * (j + 1) + 7.5);
            }
        }
    }

    Vector<Gradient<N_SIZE>, 3> d, d_C, d_F;

    callFunc2(&dof, A, b, c, d, d_C, d_F);

    Matrix<Gradient<N_SIZE>, 3, 3> A_2 = A * 2.;
    Vector<Gradient<N_SIZE>, 3> b_2 = b * 2.;

    A = Alias(A) * 2.; // Test self assignment
    b = Alias(b) * 2.;

    for (index_type i = 1; i < 3; ++i) {
    	assert(bCompare(b(i), b_2(i)));

    	for (index_type j = 1; j < 3; ++j) {
    		assert(bCompare(A(i, j), A_2(i, j)));
    	}
    }

    const doublereal dTol = sqrt(std::numeric_limits<scalar_deriv_type>::epsilon());

    Matrix<Gradient<N_SIZE>, 3, 3> A_3 = (A + Transpose(A) * 0.75) * (2. / b(1));

    A += Transpose(Alias(A)) * 0.75;
    A *= 2.;
    A /= b(1);

    for (index_type i = 1; i < 3; ++i) {
    	for (index_type j = 1; j < 3; ++j) {
    		assert(bCompare(A(i, j), A_3(i, j), dTol));
    	}
    }

    Matrix<Gradient<N_SIZE>, 3, 3> A_4 = A * A(2, 3);
    A *= Alias(A(2, 3));

    for (index_type i = 1; i < 3; ++i) {
    	for (index_type j = 1; j < 3; ++j) {
    		assert(bCompare(A(i, j), A_4(i, j), dTol));
    	}
    }

    Vector<Gradient<N_SIZE>, 3> b_3 = b * b(2);

    b *= Alias(b(2));

    for (index_type i = 1; i <= 3; ++i) {
    	assert(bCompare(b(i), b_3(i), dTol));
    }

    Vector<Gradient<N_SIZE>, 3> b_4 = b / sqrt(Dot(b, b));

    b /= sqrt(Dot(Alias(b), b));

    for (index_type i = 1; i <= 3; ++i) {
    	assert(bCompare(b(i), b_4(i), dTol));
    }

    Vector<Gradient<N_SIZE>, 3> b_5 = ((b + b * 0.75 / b(1)) * (2. / (A(1, 3) + A(2, 3))) * A(2, 1) + Transpose(A).GetCol(3)) / A(1, 1) / 3.8;

    b += b * 0.75 / Alias(b(1));
    b *= 2.;
    b /= (A(1, 3) + A(2, 3));
    b *= A(2, 1);
    b += Transpose(A).GetCol(3);
    b /= A(1, 1);
    b /= 3.8;

    for (index_type i = 1; i <= 3; ++i) {
    	assert(bCompare(b(i), b_5(i), dTol));
    }

    Vector<Gradient<N_SIZE>, 2> b23 = SubVector<2, 3>(b);
    Vector<Gradient<N_SIZE>, 2> b12 = SubVector<1, 2>(b);
    assert(bCompare(b23(1), b(2)));
    assert(bCompare(b23(2), b(3)));
    assert(bCompare(b12(1), b(1)));
    assert(bCompare(b12(2), b(2)));
    Vector<Gradient<N_SIZE>, 2> e = SubVector<2, 3>(b / 2.) + SubVector<1, 2>(b * 3.);
    assert(bCompare(e(1), Gradient<N_SIZE>(b(2) / 2. + b(1) * 3), dTol));
    assert(bCompare(e(2), Gradient<N_SIZE>(b(3) / 2. + b(2) * 3), dTol));

    A = Alias(A) * 0.5 + Transpose(A) * 3.5;
    b = Alias(b) * 2.3 + c * 5.0;
}

void testMatVecDouble2() {
    Matrix<doublereal, 3, 3> A;
    Vector<doublereal, 3> b, c;

    for (index_type i = 0; i < 3; ++i) {
        b(i + 1)=(100*(i + 1));
        c(i + 1)=(1000*(i + 10));

        for (index_type j = 0; j < 3; ++j) {
            A(i + 1, j + 1) = (100*(i + 1) + j + 1);
        }
    }

    Vector<doublereal, 3> d, d_C, d_F;

    callFunc2(0, A, b, c, d, d_C, d_F);
}

void testMatVecProduct() {
    Matrix<doublereal, 3, 4> A;
    Matrix<doublereal, 3, 3> R;
    Matrix<doublereal, 4, 2> B;

    for (index_type i = 0; i < A.iGetNumRows(); ++i) {
    	for (index_type j = 0; j < A.iGetNumCols(); ++j) {
    		A(i + 1, j + 1) = 10 * (i + 1) + j + 1;
    	}
    }

    for (index_type i = 0; i < R.iGetNumRows(); ++i) {
    	for (index_type j = 0; j < R.iGetNumCols(); ++j) {
    		R(i + 1, j + 1) = 10 * (i + 1) + j + 1;
    	}
    }

    for (index_type i = 0; i < B.iGetNumRows(); ++i) {
    	for (index_type j = 0; j < B.iGetNumCols(); ++j) {
    		B(i + 1, j + 1) = (10 * (i + 1) + j + 1);
    	}
    }

    std::cout << "A=" << std::endl << A << std::endl;

    Matrix<doublereal, 4, 3> A_T = Transpose(A);

    std::cout << "A^T=" << std::endl << A_T << std::endl;

    Matrix<doublereal, 4, 3> A_T2 = Transpose(Direct(A));

    for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
    	for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
    		assert(bCompare(A(i, j), A_T(j, i)));
    		assert(bCompare(A(i, j), A_T2(j, i)));
    	}
    }

    Matrix<doublereal, 3, 4> A2 = Transpose(Transpose(A));

    for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
    	for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
    		assert(bCompare(A2(i, j), A(i, j)));
    	}
    }

    std::cout << "\nrows of A:" << std::endl;

    for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
    	Matrix<doublereal, 3, 4>::RowVectorType r = A.GetRow(i);
    	std::cout << i - 1 << ": ";

    	for (index_type j = 1; j <= r.iGetNumRows(); ++j) {
    		assert(r(j) == A(i, j));
    		std::cout << r(j) << " ";
    	}

    	std::cout << std::endl;
    }

    std::cout << "\ncolumns of A:" << std::endl;

    for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
    	Matrix<doublereal, 3, 4>::ColumnVectorType c = A.GetCol(j);
    	std::cout << j - 1 << ": ";

    	for (index_type i = 1; i <= c.iGetNumRows(); ++i) {
    		assert(c(i) == A(i, j));
    		std::cout << c(i) << " ";
    	}

    	std::cout << std::endl;
    }

    Matrix<doublereal, 3, 4> D = -A;
    Matrix<doublereal, 3, 3> E = -(A * Transpose(A));


    for (int i = 1; i <= 3; ++i) {
    	for (int j = 1; j <= 3; ++j) {
    		assert(bCompare(D(i, j), -A(i, j)));
    	}
    }

    assert(E(1, 1) == -630);
    assert(E(1, 2) == -1130);
    assert(E(1, 3) == -1630);
    assert(E(2, 1) == -1130);
    assert(E(2, 2) == -2030);
    assert(E(2, 3) == -2930);
    assert(E(3, 1) == -1630);
    assert(E(3, 2) == -2930);
    assert(E(3, 3) == -4230);

    Vector<doublereal, 4> x;
    Vector<doublereal, 3> b;

    for (index_type i = 0; i < x.iGetNumRows(); ++i) {
    	x(i + 1) = 100 * (i + 1);
    }

    for (index_type i = 1; i <= b.iGetNumRows(); ++i) {
    	b(i) = Dot(A.GetRow(i), x);
    }

    Vector<Gradient<0>, 3> grad_b = Direct(b);
    Vector<doublereal, 3> b2 = A * x;
    Vector<doublereal, 3> b3 = b + A * x + b;
    Vector<doublereal, 3> b4 = Direct(A) * Direct(x);
    Vector<doublereal, 3> b5 = Direct(A) * x;
    Vector<doublereal, 3> b6 = A * Direct(x);
    Vector<doublereal, 3> c = R * (A * x);
    Vector<doublereal, 4> d = Transpose(A) * b;
    Vector<doublereal, 4> e = Transpose(A) * (A * x);
    Vector<doublereal, 3> f = A * (Transpose(A) * (A * (x + x)));

    for (index_type i = 1; i <= d.iGetNumRows(); ++i) {
    	assert(d(i) == e(i));
    }

    Matrix<doublereal, 3, 2> C = Direct(A) * Direct(B);
    Vector<doublereal, 2> g;
    g(1) = 523;
    g(2) = -786;
    Vector<doublereal, 3> h = A * (B * g);
    Vector<doublereal, 3> h2 = (Direct(A) * Direct(B)) * g;
    Vector<doublereal, 3> h3 = (Direct(A) * B) * g;
    Vector<doublereal, 3> h4 = (A * Direct(B)) * g;
    Vector<doublereal, 3> h5 = (A * B) * g;
    Vector<doublereal, 3> h6 = (A * B) * Direct(g);

    for (index_type i = 1; i <= h.iGetNumRows(); ++i) {
    	assert(h(i) == h2(i));
    	assert(h(i) == h3(i));
    	assert(h(i) == h4(i));
    	assert(h(i) == h5(i));
    	assert(h(i) == h6(i));
    }

    std::cout << "B=" << std::endl << B << std::endl;
    std::cout << "A * B  = C" << std::endl << C << std::endl;

    std::cout << "R=" << std::endl << R << std::endl;
    std::cout << "x = " << std::endl << x << std::endl;
    std::cout << "A * x = b" << std::endl << b << std::endl;
    std::cout << "R * A * x = c" << std::endl << c << std::endl;
    std::cout << "A^T * b = d" << std::endl << d << std::endl;
    std::cout << "A^T * A * x = e" << std::endl << e << std::endl;
    std::cout << "A * A^T * A * 2 * x = f" << std::endl << f << std::endl;
    std::cout << "A * B * g = h" << std::endl << h << std::endl;

    for (index_type i = 1; i <= b2.iGetNumRows(); ++i) {
    	assert(b2(i) == b(i));
    	assert(b3(i) == 3 * b(i));
    	assert(b4(i) == b(i));
    	assert(b5(i) == b(i));
    	assert(b6(i) == b(i));
    }

    assert(b(1) == 13000);
    assert(b(2) == 23000);
    assert(b(3) == 33000);

    assert(c(1) == 848000);
    assert(c(2) == 1538000);
    assert(c(3) == 2228000);

    assert(d(1) == 1649000.00);
    assert(d(2) == 1718000.00);
    assert(d(3) == 1787000.00);
    assert(d(4) == 1856000.00);

    assert(C(1, 1) == 1.35e3);
    assert(C(2, 1) == 2.39e3);
    assert(C(3, 1) == 3.43e3);
    assert(C(1, 2) == 1.4e3);
    assert(C(2, 2) == 2.48e3);
    assert(C(3, 2) == 3.56e3);

    assert(h(1) == -394350);
    assert(h(2) == -699310);
    assert(h(3) == -1004270);

    Vector<doublereal, 2> b23 = SubVector<2, 3>(b);
    Vector<doublereal, 2> b12 = SubVector<1, 2>(b);
    assert(bCompare(b23(1), b(2)));
    assert(bCompare(b23(2), b(3)));
    assert(bCompare(b12(1), b(1)));
    assert(bCompare(b12(2), b(2)));
    Vector<doublereal, 2> e123 = SubVector<2, 3>(b / 2.) + SubVector<1, 2>(b * 3.);
    assert(bCompare(e123(1), (b(2) / 2. + b(1) * 3)));
    assert(bCompare(e123(2), (b(3) / 2. + b(2) * 3)));

    Matrix<doublereal, 5, 7> F;

    for (index_type i = 1; i <= 5; ++i) {
    	for (index_type j = 1; j <= 7; ++j) {
    		F(i, j) = i * 10 + j;
    	}
    }

    Matrix<doublereal, 2, 3> G = SubMatrix<3, 4, 5, 7>(F);

    for (index_type i = 1; i <= 2; ++i) {
    	for (index_type j = 1; j <= 3; ++j) {
    		assert(bCompare(G(i, j), F(i + 2, j + 4)));
    	}
    }

    Vector<doublereal, 3> G1 = SubMatrix<4, 4, 2, 4>(F).GetRow(1);

    assert(bCompare(G1(1), F(4, 2)));
    assert(bCompare(G1(2), F(4, 3)));
    assert(bCompare(G1(3), F(4, 4)));

    Vector<doublereal, 2> G2 = SubMatrix<1, 2, 5, 7>(F).GetCol(1);

    assert(bCompare(G2(1), F(1, 5)));
    assert(bCompare(G2(2), F(2, 5)));

    Vector<doublereal, 2> G3 = SubMatrix<1, 2, 5, 7>(F).GetCol(2);

    assert(bCompare(G3(1), F(1, 6)));
    assert(bCompare(G3(2), F(2, 6)));


    Vector<doublereal, 2> G4 = SubMatrix<1, 2, 5, 7>(F).GetCol(3);

    assert(bCompare(G4(1), F(1, 7)));
    assert(bCompare(G4(2), F(2, 7)));

    Matrix<doublereal, 3, 4> M = SubMatrix<2, 4, 2, 5>(F) * SubMatrix<2, 5, 4, 7>(F);

    assert(bCompare(M(1, 1), 3716.));
    assert(bCompare(M(1, 2), 3810.));
    assert(bCompare(M(1, 3), 3904.));
    assert(bCompare(M(1, 4), 3998.));
    assert(bCompare(M(2, 1), 5276.));
    assert(bCompare(M(2, 2), 5410.));
    assert(bCompare(M(2, 3), 5544.));
    assert(bCompare(M(2, 4), 5678.));
    assert(bCompare(M(3, 1), 6836.));
    assert(bCompare(M(3, 2), 7010.));
    assert(bCompare(M(3, 3), 7184.));
    assert(bCompare(M(3, 4), 7358.));

    std::cout << "F=\n" << Tabular(F, 5) << std::endl;
    std::cout << "G=\n" << Tabular(G, 5) << std::endl;
    std::cout << "G1=\n" << G1 << std::endl;
    std::cout << "G2=\n" << G2 << std::endl;
    std::cout << "G3=\n" << G3 << std::endl;
    std::cout << "G4=\n" << G4 << std::endl;
    std::cout << "M=\n" << Tabular(M, 5) << std::endl;
}

namespace testMatVecProductGradient_testData {

template <typename S, index_type N_rows, index_type N_SIZE>
void testGradient(const S& ref, const Vector<Gradient<N_SIZE>, N_rows>& v, doublereal dTol) {
	for (index_type i = 0; i < index_type(sizeof(ref.val)/sizeof(ref.val[0])); ++i) {
		assert(bCompare(v(i + 1).dGetValue(), ref.val[i], dTol));
		for (index_type j = 0; j < index_type(sizeof(ref.der[0])/sizeof(ref.der[0][0])); ++j) {
			std::cout << "v(" << i + 1 << ")=" << v(i + 1).dGetDerivativeGlobal(j + 1) << std::endl;
			std::cout << "ref.der[" << i << "][" << j << "]=" << ref.der[i][j] << std::endl;

			const bool bOK = bCompare<scalar_deriv_type>(v(i + 1).dGetDerivativeGlobal(j + 1), ref.der[i][j], dTol);
			std::cout << "dTol=" << dTol << " :[" << (bOK ? "OK" : "NOK") << "]" << std::endl;
			assert(bOK);
		}
	}
}

template <typename S, index_type N_SIZE>
void testGradient(const S& ref, const Gradient<N_SIZE>& g, doublereal dTol) {
	for (index_type i = 0; i < index_type(sizeof(ref.val)/sizeof(ref.val[0])); ++i) {
		assert(bCompare(g.dGetValue(), ref.val[i], dTol));
		for (index_type j = 0; j < index_type(sizeof(ref.der[0])/sizeof(ref.der[0][0])); ++j) {
			assert(bCompare<scalar_deriv_type>(g.dGetDerivativeGlobal(j + 1), ref.der[i][j], dTol));
		}
	}
}

static const struct test_b {
doublereal val[3];
doublereal der[3][12];
} oct_b = {
{6300,
11300,
16300},
{{220, 0, 0, 240, 0, 0, 260, 0, 0, 280, 0, 0},
{210, 110, 0, 220, 120, 0, 230, 130, 0, 240, 140, 0},
{310, 0, 110, 320, 0, 120, 330, 0, 130, 340, 0, 140}}};

static const struct test_c {
doublereal val[3];
doublereal der[3][12];
} oct_c = {
{-6300,
-11300,
-16300},
{{-220, 0, 0, -240, 0, 0, -260, 0, 0, -280, 0, 0},
{-210, -110, 0, -220, -120, 0, -230, -130, 0, -240, -140, 0},
{-310, 0, -110, -320, 0, -120, -330, 0, -130, -340, 0, -140}}};

static const struct test_d {
doublereal val[3];
doublereal der[3][12];
} oct_d = {
{0,
0,
0},
{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}};

static const struct test_e {
doublereal val[3];
doublereal der[3][12];
} oct_e = {
{0,
0,
0},
{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}};

static const struct test_f {
doublereal val[3];
doublereal der[3][12];
} oct_f = {
{0,
0,
0},
{{0, 0, -0, 0, 0, -0, 0, 0, -0, 0, 0, -0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, -0, 0, 0, -0, 0, 0, -0, 0, 0, -0, 0}}};

static const struct test_g {
doublereal val[3];
doublereal der[3][12];
} oct_g = {
{6300,
11300,
16300},
{{220, 0, 0, 240, 0, 0, 260, 0, 0, 280, 0, 0},
{210, 110, 0, 220, 120, 0, 230, 130, 0, 240, 140, 0},
{310, 0, 110, 320, 0, 120, 330, 0, 130, 340, 0, 140}}};

static const struct test_i {
doublereal val[3];
doublereal der[3][12];
} oct_i = {
{0,
0,
-0},
{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{-0, -0, -0, -0, -0, -0, -0, -0, -0, -0, -0, -0}}};

static const struct test_j {
doublereal val[3];
doublereal der[3][12];
} oct_j = {
{-74031.42857142857,
-50797.14285714286,
63828.57142857143},
{{-1390.571428571429, -389.7142857142857, -229.4285714285714, -1446.857142857143, -425.1428571428572, -250.2857142857143, -1503.142857142857, -460.5714285714286, -271.1428571428572, -1559.428571428571, -496, -292},
{-611.1428571428571, 0, -493.4285714285714, -585.1428571428571, 0, -538.2857142857143, -559.1428571428571, 0, -583.1428571428571, -533.1428571428571, 0, -628},
{1400.857142857143, 493.4285714285714, 0, 1487.428571428571, 538.2857142857143, 0, 1574, 583.1428571428571, 0, 1660.571428571429, 628, 0}}};

static const struct test_l {
doublereal val[3];
doublereal der[3][12];
} oct_l = {
{163936.0000000001,
-1920098.285714285,
-1337944.571428571},
{{-2648.085714285714, -3602.028571428571, 6118.514285714286, -3602.457142857142, -3929.485714285714, 6674.742857142858, -4556.82857142857, -4256.942857142857, 7230.971428571428, -5511.200000000002, -4584.4, 7787.2},
{-39236.54285714286, -12579.28571428571, -2844.914285714286, -41293.65714285714, -13722.85714285714, -3103.542857142857, -43350.77142857143, -14866.42857142857, -3362.171428571429, -45407.88571428572, -16010, -3620.8},
{-19746.11428571428, -2844.914285714286, -9421.657142857142, -19748.8, -3103.542857142857, -10278.17142857143, -19751.48571428571, -3362.171428571428, -11134.68571428571, -19754.17142857143, -3620.8, -11991.2}}};

static const struct test_m {
doublereal val[3];
doublereal der[3][12];
} oct_m = {
{819680.0000000002,
-9600491.428571427,
-6689722.857142856},
{{-13240.42857142857, -18010.14285714286, 30592.57142857143, -18012.28571428571, -19647.42857142857, 33373.71428571429, -22784.14285714285, -21284.71428571429, 36154.85714285714, -27556.00000000001, -22922, 38936},
{-196182.7142857143, -62896.42857142857, -14224.57142857143, -206468.2857142857, -68614.28571428572, -15517.71428571428, -216753.8571428572, -74332.14285714284, -16810.85714285714, -227039.4285714286, -80050, -18104},
{-98730.57142857142, -14224.57142857143, -47108.28571428571, -98743.99999999997, -15517.71428571429, -51390.85714285714, -98757.42857142858, -16810.85714285714, -55673.42857142857, -98770.85714285713, -18104, -59956}}};

static const struct test_r {
doublereal val[1];
doublereal der[1][12];
} oct_r = {
{218540},
{{5765, -803, 1364, 6130, -876, 1488, 6495, -949, 1612, 6860, -1022, 1736}}};

static const struct test_s {
doublereal val[1];
doublereal der[1][12];
} oct_s = {
{218540},
{{5765, -803, 1364, 6130, -876, 1488, 6495, -949, 1612, 6860, -1022, 1736}}};

static const struct test_t {
doublereal val[1];
doublereal der[1][12];
} oct_t = {
{433070000},
{{17624000, 2486000, 3586000, 18428000, 2712000, 3912000, 19232000, 2938000, 4238000, 20036000, 3164000, 4564000}}};

static const struct test_norm_g {
        doublereal val[1];
doublereal der[1][12];
} oct_norm_g = {
{20810.33397137105},
{{423.4434686210582, 59.72994002450923, 86.15911702650448, 442.760794357062, 65.15993457219189, 93.99176402891398, 462.0781200930658, 70.58992911987455, 101.8244110313235, 481.3954458290696, 76.01992366755721, 109.657058033733}}};

}

template <typename T, index_type N_rows>
inline T
Norm_1(const Vector<T, N_rows>& u) {
	return sqrt(Dot(u, u));
}

#ifdef HAVE_FC_F77
extern "C" void __FC_DECL__(func2addad_dv)(const doublereal x[],
										   const doublereal xd[],
										   const doublereal y[],
										   const doublereal yd[],
										   doublereal z[],
										   doublereal zd[],
										   const integer& n,
										   const integer& nbdirs);

extern "C" void __FC_DECL__(func2mulad_dv)(const doublereal x[],
										   const doublereal xd[],
										   const doublereal y[],
										   const doublereal yd[],
										   doublereal z[],
										   doublereal zd[],
										   const integer& n,
										   const integer& nbdirs);
#endif // HAVE_FC_F77

template <typename T, index_type iRowCount>
void doVecAdd(const Vector<T, iRowCount>& x, const Vector<T, iRowCount>& y, Vector<T, iRowCount>& z)
{
	z = x + y;
}

template <typename T, index_type iRowCount>
void doVecMul(const Vector<T, iRowCount>& x, const Vector<T, iRowCount>& y, Vector<T, iRowCount>& z)
{
	for (index_type i = 1; i <= z.iGetNumRows(); ++i)
	{
		z(i) = x(i) * y(i);
	}
}

template <index_type iRowCount, index_type iMaxDeriv, typename Function, typename FunctionF77>
void testVecOp(const int M, const int N, Function f, FunctionF77 f77, const char* function) {
    Vector<Gradient<iMaxDeriv>, iRowCount> x, y, z;
    LocalDofMap dof;

    for (index_type i = 1; i <= x.iGetNumRows(); ++i) {
		x(i).SetValuePreserve(i * 100);
		x(i).DerivativeResizeReset(&dof, 0, N, MapVectorBase::GLOBAL, 0.);
		y(i).SetValuePreserve(i);
		y(i).DerivativeResizeReset(&dof, 0, N, MapVectorBase::GLOBAL, 0.);

        for (int j = 0; j < N; ++j) {
			x(i).SetDerivativeGlobal(j, (j + 1));
			y(i).SetDerivativeGlobal(j, (j + 1));
		}
    }

    double start = mbdyn_clock_time();

    for (int i = 0; i < M; ++i) {
    	f(x, y, z);
    }

    const double dt = mbdyn_clock_time() - start;

    std::cerr << "C++: testVecAdd<" << iRowCount << "," << iMaxDeriv << ">(" << M << ", " << N << ", \"" << function << "\"):" << dt << std::endl;

    doublereal xF[iRowCount], yF[iRowCount], zF[iRowCount];
    doublereal *xdF = new doublereal[iRowCount*N];
    doublereal *ydF = new doublereal[iRowCount*N];
    doublereal *zdF = new doublereal[iRowCount*N];

    for (int i = 0; i < iRowCount; ++i) {
    	xF[i] = x(i + 1).dGetValue();
    	yF[i] = y(i + 1).dGetValue();
    	zF[i] = 0.;

    	for (int j = 0; j < N; ++j) {
    		xdF[i * N + j] = x(i + 1).dGetDerivativeLocal(j);
    		ydF[i * N + j] = y(i + 1).dGetDerivativeLocal(j);
    		zdF[i * N + j] = 0.;
    	}
    }

    start = mbdyn_clock_time();

    for (int i = 0; i < M; ++i) {
    	f77(xF, xdF, yF, ydF, zF, zdF, iRowCount, N);
    }

    const double dtF77 = mbdyn_clock_time() - start;

    std::cerr << "F77: testVecAdd<" << iRowCount << "," << iMaxDeriv << ">(" << M << ", " << N << ", \"" << function << "\"):" << dtF77 << std::endl;
    std::cerr << "overhead=" << dt/std::max(std::numeric_limits<doublereal>::epsilon(), dtF77) << std::endl;

    for (int i = 0; i < iRowCount; ++i) {
    	assert(xF[i] == x(i + 1).dGetValue());
    	assert(yF[i] == y(i + 1).dGetValue());
    	assert(zF[i] == z(i + 1).dGetValue());

    	for (int j = 0; j < N; ++j) {
    		assert(xdF[i * N + j] == x(i + 1).dGetDerivativeLocal(j));
    		assert(ydF[i * N + j] == y(i + 1).dGetDerivativeLocal(j));
    		assert(zdF[i * N + j] == z(i + 1).dGetDerivativeLocal(j));
    	}
    }

    delete [] xdF;
    delete [] ydF;
    delete [] zdF;
}

void testMatVecProductGradient() {
    Matrix<Gradient<0>, 3, 4> A;
    LocalDofMap dof;

    for (index_type i = 0; i < A.iGetNumRows(); ++i) {
    	for (index_type j = 0; j < A.iGetNumCols(); ++j) {
    		A(i + 1, j + 1).SetValue(10 * (i + 1) + j + 1);
    		A(i + 1, j + 1).DerivativeResizeReset(&dof, j * A.iGetNumRows() + i + 1, MapVectorBase::GLOBAL, 0.);
    		A(i + 1, j + 1).SetDerivativeGlobal(j * A.iGetNumRows() + i + 1, 1.);
    	}
    }

    Matrix<Gradient<0>, 3, 4> B1;
    Matrix<doublereal, 3, 4> Ad;

    for (int i = 1; i <= 3; ++i) {
    	for (int j = 1; j <= 4; ++j) {
    		B1(i, j) = 3.5 * A(i, j);
    		Ad(i, j) = A(i, j).dGetValue();
    	}
    }

    Matrix<Gradient<0>, 3, 4> C1 = A + B1;
    Matrix<Gradient<0>, 3, 4> D1 = (A + B1) - (C1 + B1);
    Matrix<Gradient<0>, 3, 4> E1 = A + (C1 - B1);
    Matrix<Gradient<0>, 3, 4> F1 = (C1 - B1) - A;
    Matrix<Gradient<0>, 3, 4> G1 = A - B1;
    Matrix<Gradient<0>, 3, 4> H1 = A * 3.5 + B1 / 4.;
    Matrix<Gradient<0>, 3, 4> I1 = (A * 2. - B1 / 5.) * 3.5;
    Matrix<Gradient<0>, 3, 4> J1 = A * B1(1, 1) + I1 * H1(2, 2);
    Matrix<Gradient<0>, 3, 4> K1 = (I1 + J1) * H1(3, 4);

    const doublereal dTol = 10 * std::numeric_limits<scalar_deriv_type>::epsilon();

    for (int i = 1; i <= 3; ++i) {
    	for (int j = 1; j <= 4; ++j) {
    		assert(bCompare(C1(i, j), Gradient<0>(A(i, j) + B1(i, j)), dTol));
    		assert(bCompare(D1(i, j), Gradient<0>(A(i, j) + B1(i, j) - C1(i, j) - B1(i, j)), dTol));
    		assert(bCompare(E1(i, j), Gradient<0>(A(i, j) + C1(i, j) - B1(i, j)), dTol));
    		assert(bCompare(F1(i, j), Gradient<0>(C1(i, j) - B1(i, j) - A(i, j)), dTol));
    		assert(bCompare(G1(i, j), Gradient<0>(A(i, j) - B1(i, j)), dTol));
    		assert(bCompare(H1(i, j), Gradient<0>(A(i, j) * 3.5 + B1(i, j) / 4.), dTol));
    		assert(bCompare(I1(i, j), Gradient<0>((A(i, j) * 2. - B1(i, j) / 5.) * 3.5), dTol));
    		assert(bCompare(J1(i, j), Gradient<0>(A(i, j) * B1(1, 1) + I1(i, j) * H1(2, 2)), dTol));
    		assert(bCompare(K1(i, j), Gradient<0>((I1(i, j) + J1(i, j)) * H1(3, 4)), dTol));
    	}
    }

    std::cout << "A=" << std::endl << A << std::endl;

    std::cout << "A=" << std::endl << A << std::endl;

    Matrix<Gradient<0>, 4, 3> A_T = Transpose(A);

    std::cout << "A^T=" << std::endl << A_T << std::endl;

    Matrix<Gradient<0>, 4, 3> A_T2 = Transpose(Direct(A));

    Matrix<Gradient<0>, 3, 4> B = -A;

    Matrix<Gradient<0>, 3, 3> D = -(A * Transpose(A));

    for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
    	for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
    		assert(bCompare(A(i, j), A_T(j, i)));
    		assert(bCompare(A(i, j), A_T2(j, i)));
    		assert(bCompare(B(i, j), Gradient<0>(-A(i, j))));
    	}
    }

    Matrix<Gradient<0>, 3, 4> A2 = Transpose(Transpose(A));

    for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
    	for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
    		assert(bCompare(A2(i, j), A(i, j)));
    	}
    }

    std::cout << "\nrows of A:" << std::endl;

    for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
    	Matrix<Gradient<0>, 3, 4>::RowVectorType r = A.GetRow(i);
    	std::cout << i - 1 << ": ";

    	for (index_type j = 1; j <= r.iGetNumRows(); ++j) {
    		assert(bCompare(r(j), A(i, j), std::numeric_limits<scalar_deriv_type>::epsilon()));
    		std::cout << r(j) << " ";
    	}

    	std::cout << std::endl;
    }

    std::cout << "\ncolumns of A:" << std::endl;

    for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
    	Matrix<Gradient<0>, 3, 4>::ColumnVectorType c = A.GetCol(j);
    	std::cout << j - 1 << ": ";

    	for (index_type i = 1; i <= c.iGetNumRows(); ++i) {
    		assert(bCompare(c(i), A(i, j), std::numeric_limits<scalar_deriv_type>::epsilon()));
    		std::cout << c(i) << " ";
    	}

    	std::cout << std::endl;
    }

    Vector<Gradient<0>, 4> x;
    Vector<Gradient<0>, 3> b;
    Vector<doublereal, 4> v;

    for (index_type i = 1; i <= x.iGetNumRows(); ++i) {
    	x(i) = A(1, i) * 10.;
    	v(i) = x(i).dGetValue();
    }

    Vector<Gradient<0>, 3> b_1;
    Gradient<0> b1_2;
    b1_2 = A(1, 1) * x(1);
    b1_2 += A(1, 2) * x(2);
    b1_2 += A(1, 3) * x(3);
    b1_2 += A(1, 4) * x(4);

    b_1(1) = A(1, 1) * x(1) + A(1, 2) * x(2) + A(1, 3) * x(3) + A(1, 4) * x(4);

    assert(b1_2.bIsEqual(b_1(1)));

    b_1(2) = A(2, 1) * x(1) + A(2, 2) * x(2) + A(2, 3) * x(3) + A(2, 4) * x(4);
    b_1(3) = A(3, 1) * x(1) + A(3, 2) * x(2) + A(3, 3) * x(3) + A(3, 4) * x(4);


    using namespace testMatVecProductGradient_testData;
    testGradient(oct_b, b_1, dTol);

    for (index_type i = 1; i <= b.iGetNumRows(); ++i) {
    	b(i) = Dot(A.GetRow(i), x);
    }

    testGradient(oct_b, b, dTol);

    Vector<Gradient<0>, 3> b2 = A * x;
    Vector<Gradient<0>, 3> b2d = Ad * x;
    Vector<Gradient<0>, 3> b2d1 = -Ad * x - Ad * x;
    Vector<Gradient<0>, 3> b3 = b + A * x + b2;
    Vector<Gradient<0>, 3> b4 = A * (x + x);
    Vector<Gradient<0>, 3> b5 = Direct(A) * Direct(x);
    Vector<Gradient<0>, 3> b6 = Direct(A) * x;
    Vector<Gradient<0>, 3> b7 = A * Direct(x);
    Vector<Gradient<0>, 3> b8 = A * v;
    Vector<Gradient<0>, 3> c = -(A * x);
    Vector<Gradient<0>, 3> d = Cross(c, b) * 5.;
    Vector<Gradient<0>, 3> e = Cross(c + d, b) / 2.;
    Vector<Gradient<0>, 3> f = Cross(e, d + c) - e;
    Vector<Gradient<0>, 3> g = Cross(d + e, f - c) + b;
    Matrix<Gradient<0>, 3, 3> skew_g(MatCrossVec(g));
    Matrix<Gradient<0>, 3, 3> skew_gmf(MatCrossVec(g - f));
    Vector<Gradient<0>, 3> gmf = g - f;
    Gradient<0> norm_Ax = Norm(A * x);
    Gradient<0> z = Norm(A * x) / Norm(Ad * x) + (Norm(A * x) / Norm(Ad * x) - Norm(A * x) / Norm(Ad * x)) * exp(-pow(Norm(x) / 5., 3./2.));
    doublereal zd = Norm(Ad * v) / Norm(Ad * v) + (Norm(Ad * v) / Norm(Ad * v) - Norm(Ad * v) / Norm(Ad * v)) * exp(-pow(Norm(v) / 5., 3./2.));

    for (int i = 1; i <= 3; ++i) {
    	for (int j = 1; j <= 3; ++j) {
    		if (i == j) {
    			assert(skew_g(i, j).bIsEqual(Gradient<0>()));
    			assert(skew_gmf(i, j).bIsEqual(Gradient<0>()));
    		} else {
    			assert(skew_g(i, j).bIsEqual(-skew_g(j, i)));
    			assert(skew_gmf(i, j).bIsEqual(-skew_gmf(j, i)));
    		}
    	}
    }

    assert(skew_g(1, 2).bIsEqual(-g(3)));
    assert(skew_g(2, 1).bIsEqual(g(3)));
    assert(skew_g(1, 3).bIsEqual(g(2)));
    assert(skew_g(3, 1).bIsEqual(-g(2)));
    assert(skew_g(2, 3).bIsEqual(-g(1)));
    assert(skew_g(3, 2).bIsEqual(g(1)));

    assert(skew_gmf(1, 2).bIsEqual(-gmf(3)));
    assert(skew_gmf(2, 1).bIsEqual(gmf(3)));
    assert(skew_gmf(1, 3).bIsEqual(gmf(2)));
    assert(skew_gmf(3, 1).bIsEqual(-gmf(2)));
    assert(skew_gmf(2, 3).bIsEqual(-gmf(1)));
    assert(skew_gmf(3, 2).bIsEqual(gmf(1)));

    Vector<doublereal, 3> h;
    h(1) = 15.7;
    h(2) = -7.3;
    h(3) = 12.4;

    Vector<Gradient<0>, 3> i = Cross(d, h) * 5.7;
    Vector<Gradient<0>, 3> j = Cross(h, g) / 3.5;
    Vector<doublereal, 3> k = Cross(h, h) + h * 3.;
    Vector<Gradient<0>, 3> l = Cross(h + h, j) / 2.;
    Vector<Gradient<0>, 3> m = Cross(h, j - i) * 5.;
    Vector<doublereal, 3> n = Cross(h + h, h) + h;
    Vector<doublereal, 3> o = Cross(h, h + h) * 2.;
    Vector<doublereal, 3> p = Cross(h * 2.5, h/3.) - h;
    Gradient<0> r = Dot(h, g);
    Gradient<0> s = Dot(g, h);
    Gradient<0> t = Dot(g, g);
    doublereal u = Dot(h, h);
    Gradient<0> r1 = Dot(Direct(h), g);
    Gradient<0> r2 = Dot(h, Direct(g));
    Gradient<0> r3 = Dot(Direct(h), Direct(g));
    Gradient<0> norm_g1 = Norm(g);
    Gradient<0> norm_g2 = Norm(Cross(d + e, f - c) + b);

    assert(r.bIsEqual(r1));
    assert(r.bIsEqual(r2));
    assert(r.bIsEqual(r3));

    std::cout << "x = " << std::endl << x << std::endl;
    std::cout << "b = A * x" << std::endl << b << std::endl;
    std::cout << "c = -A * x" << std::endl << c << std::endl;
    std::cout << "d = cross(c, b) * 5" << std::endl << d << std::endl;
    std::cout << "e = cross(c + d, b) / 2" << std::endl << e << std::endl;
    std::cout << "f = cross(e, d + c) - e" << std::endl << f << std::endl;
    std::cout << "g = cross(d + e, f - c) + b" << std::endl << g << std::endl;
    std::cout << "h = " << std::endl << h << std::endl;
    std::cout << "i = cross(d, h) * 5.7" << std::endl << i << std::endl;
    std::cout << "j = cross(h, g) / 3.5" << std::endl << j << std::endl;
    std::cout << "k = cross(h, h) + h * 3" << std::endl << k << std::endl;
    std::cout << "l = cross(h + h, j) / 2" << std::endl << l << std::endl;
    std::cout << "m = cross(h, j - i) * 5" << std::endl << m << std::endl;
    std::cout << "n = cross(h + h, h) + h" << std::endl << n << std::endl;
    std::cout << "o = cross(h, h + h) * 2" << std::endl << o << std::endl;
    std::cout << "p = cross(h * 2.5, h / 3) - h" << std::endl << p << std::endl;
    std::cout << "r = Dot(h, g) = " << r << std::endl;
    std::cout << "s = Dot(g, h) = " << s << std::endl;
    std::cout << "t = Dot(g, g) = " << t << std::endl;
    std::cout << "u = Dot(h, h) = " << u << std::endl;
    std::cout << "norm_g1 = Norm(g) = " << norm_g1 << std::endl;
    std::cout << "norm_g2 = Norm(Cross(d + e, f - c) + b) = " << norm_g2 << std::endl;
    std::cout << "z=" << z << std::endl;
    std::cout << "zd=" << zd << std::endl;

    for (index_type i = 1; i <= b2.iGetNumRows(); ++i) {
    	assert(bCompare(b2d(i).dGetValue(), b2(i).dGetValue(), dTol));
    	assert(bCompare(b2d1(i).dGetValue(), -2 * b2(i).dGetValue(), dTol));
    	assert(bCompare(b2(i), b(i), dTol));
    	assert(bCompare(b3(i), Gradient<0>(3 * b(i)), dTol));
    	assert(bCompare(b4(i), Gradient<0>(2 * b(i)), dTol));
    	assert(bCompare(b5(i), b(i), dTol));
    	assert(bCompare(b6(i), b(i), dTol));
    	assert(bCompare(b7(i), b(i), dTol));
    	assert(bCompare<scalar_deriv_type>(b8(i).dGetValue(), b(i).dGetValue(), dTol));
    	assert(bCompare(c(i), Gradient<0>(-b(i)), dTol));
    }


    assert(bCompare(b(1).dGetValue(), 6300., dTol));
    assert(bCompare(b(2).dGetValue(), 11300., dTol));
    assert(bCompare(b(3).dGetValue(), 16300., dTol));

    assert(bCompare(D(1, 1).dGetValue(), -630., dTol));
    assert(bCompare(D(1, 2).dGetValue(), -1130., dTol));
    assert(bCompare(D(1, 3).dGetValue(), -1630., dTol));
    assert(bCompare(D(2, 1).dGetValue(), -1130., dTol));
    assert(bCompare(D(2, 2).dGetValue(), -2030., dTol));
    assert(bCompare(D(2, 3).dGetValue(), -2930., dTol));
    assert(bCompare(D(3, 1).dGetValue(), -1630., dTol));
    assert(bCompare(D(3, 2).dGetValue(), -2930., dTol));
    assert(bCompare(D(3, 3).dGetValue(), -4230., dTol));

    testGradient(oct_b, b, dTol);
    testGradient(oct_c, c, dTol);
    testGradient(oct_d, d, dTol);
    testGradient(oct_e, e, dTol);
    testGradient(oct_f, f, dTol);
    testGradient(oct_g, g, dTol);
    testGradient(oct_i, i, dTol);
    testGradient(oct_j, j, dTol);
    testGradient(oct_l, l, dTol);
    testGradient(oct_m, m, dTol);
    testGradient(oct_r, r, dTol);
    testGradient(oct_s, s, dTol);
    testGradient(oct_t, t, dTol);
    testGradient(oct_norm_g, norm_g1, dTol);
    testGradient(oct_norm_g, norm_g2, dTol);

    A.Reset();
    x.Reset();
    b.Reset();
    std::cout << "after reset ...\n";
    std::cout << "A=" << A << std::endl;
    std::cout << "x=" << x << std::endl;
    std::cout << "b=" << b << std::endl;
}

template <index_type N_SIZE>
void testMatVecProductGradient2(index_type iNumDeriv, int N) {
    Matrix<Gradient<N_SIZE>, 3, 5> A;
    Matrix<Gradient<N_SIZE>, 5, 7> B;
    LocalDofMap dof;

    srand(0);
    
    for (index_type i = 0; i < A.iGetNumRows(); ++i) {
    	for (index_type j = 0; j < A.iGetNumCols(); ++j) {
    		A(i + 1, j + 1).SetValue(10. * (i + 1) + j + 1 + rand());
    		A(i + 1, j + 1).DerivativeResizeReset(&dof, 0, iNumDeriv, MapVectorBase::GLOBAL, 0.);
            for (int k = 0; k < iNumDeriv; ++k) {
                A(i + 1, j + 1).SetDerivativeGlobal(k, 1000. * (i + 1) + 100. * (j + 1) + k + 1 + rand());
            }
    	}
    }
    
    for (index_type i = 0; i < B.iGetNumRows(); ++i) {
    	for (index_type j = 0; j < B.iGetNumCols(); ++j) {
    		B(i + 1, j + 1).SetValue(1000. * (i + 1) + 100. * (j + 1) + rand());
    		B(i + 1, j + 1).DerivativeResizeReset(&dof, 0, iNumDeriv, MapVectorBase::GLOBAL, 0.);
            for (int k = 0; k < iNumDeriv; ++k) {
                B(i + 1, j + 1).SetDerivativeGlobal(k, 10000. * (i + 1) + 1000. * (j + 1) + 10. * (k + 1) + rand());
            }
    	}
    }    

    Matrix<Gradient<N_SIZE>, 3, 7> C;
    
    double start = mbdyn_clock_time();
    
    for (int i = 0; i < N; ++i) {
        C = A * B;
    }

    double dt = (mbdyn_clock_time() - start) / N;
    
    std::cerr << "iMaxDerivatives=" << iNumDeriv << std::endl;
    std::cerr << "dt=" << dt << std::endl;
    std::cout << "A=\n" << A << std::endl;
    std::cout << "B=\n" << B << std::endl;
    std::cout << "C=\n" << C << std::endl;
}

#ifdef HAVE_BLITZ
template <index_type N>
void testMatVecDoubleBlitz(doublereal c_C[N]) {
    blitz::TinyVector<doublereal, N> a, b;

    for (index_type i = 0; i < N; ++i) {
        a(i) = 100*(i + 1);
        b(i) = 1000*(i + 10);
    }

    blitz::TinyVector<doublereal, N> c, c1;

	srand(0);
    tic();
    for (int i = 0; i < NLoops; ++i) {
    	func(0, a, b, c);
    }

    for (int i = 0; i < N; ++i) {
    	c_C[i] = c(i);
    }

    std::cerr << "blitz vector (doublereal): " << toc() << "s" << std::endl;

    std::cout << "a=" << std::endl << a << std::endl;
    std::cout << "b=" << std::endl << b << std::endl;
    std::cout << "c=" << std::endl << c << std::endl;
}
#endif

template <index_type N_rows>
void testMatVecCopy() {
	std::cerr << "---------------------------\ntestMatVecCopy<" << N_rows << ">()\n";

	LocalDofMap dofMap1;
	Vector<Gradient<0>, N_rows> v;

	for (int i = 1; i <= v.iGetNumRows(); ++i) {
		v(i).SetValuePreserve(i);
		v(i).DerivativeResizeReset(&dofMap1, i, MapVectorBase::GLOBAL, i);
	}

	std::cout << "v=" << std::endl << v << std::endl;

	for (int i = 2; i <= v.iGetNumRows() - 1; ++i) {
		LocalDofMap dofMap2;
		Gradient<4> g1(v(1), &dofMap2), gi(v(i), &dofMap2), gim1(v(i - 1), &dofMap2), gip1(v(i + 1), &dofMap2);
		Gradient<4> x = 1000 * g1 + 100 * gi + 10 * gip1 + gim1;
		std::cout << "x(" << i << ")=" << x << std::endl;

		assert(x.dGetValue() == 1000 + 100 * i + 10 * (i + 1) + i - 1);

		if (i == 2) {
			assert(x.dGetDerivativeGlobal(1) == 1001.);
		} else {
			assert(x.dGetDerivativeGlobal(1) == 1000.);
			assert(x.dGetDerivativeGlobal(i - 1) == i - 1);
		}

		assert(x.dGetDerivativeGlobal(i) == 100. * i);
		assert(x.dGetDerivativeGlobal(i + 1) == 10. * (i + 1));
	}

	v.Reset();
	std::cout << "v=" << v << std::endl;
}

namespace gradVecAssTest {
const int I1 = 1, I2 = 2, I3 = 3;

template <typename T>
Matrix<T, 3, 3>& Euler123ToMatR(const Vector<T, 3>& v, Matrix<T, 3, 3>& R) {
	T d = v(1);
	T dCosAlpha(cos(d));
	T dSinAlpha(sin(d));
	d = v(2);
	T dCosBeta(cos(d));
	T dSinBeta(sin(d));
	d = v(3);
	T dCosGamma(cos(d));
	T dSinGamma(sin(d));

	R(1, 1) = dCosBeta*dCosGamma;
	R(2, 1) = dCosAlpha*dSinGamma + dSinAlpha*dSinBeta*dCosGamma;
	R(3, 1) = dSinAlpha*dSinGamma - dCosAlpha*dSinBeta*dCosGamma;
	R(1, 2) = -dCosBeta*dSinGamma;
	R(2, 2) = dCosAlpha*dCosGamma - dSinAlpha*dSinBeta*dSinGamma;
	R(3, 2) = dSinAlpha*dCosGamma + dCosAlpha*dSinBeta*dSinGamma;
	R(1, 3) = dSinBeta;
	R(2, 3) = -dSinAlpha*dCosBeta;
	R(3, 3) = dCosAlpha*dCosBeta;

	return R;
}

#ifdef HAVE_BLITZ

class Node {
public:
    Node(const Vector<doublereal, 3>& X_0,
         const Vector<doublereal, 3>& XP_0,
         const Matrix<doublereal, 3, 3>& R_0,
         const Vector<doublereal, 3>& W_0)
        :iFirstDofIndex(-1), R0(R_0), W0(W_0) {

        for (int i = 0; i < 3; ++i) {
            X(i + 1).SetValuePreserve(X_0(i + 1));
            X(i + 1).DerivativeResizeReset(&dof, i, MapVectorBase::LOCAL, 0.);

            XP(i + 1).SetValuePreserve(XP_0(i + 1));
            XP(i + 1).DerivativeResizeReset(&dof, i, MapVectorBase::LOCAL, 0.);
            XP(i + 1).SetDerivativeLocal(i, -1.); // derivative will be always -1

            g(i + 1).SetValuePreserve(0.);
            g(i + 1).DerivativeResizeReset(&dof, i, MapVectorBase::LOCAL, 0.);

            gP(i + 1).SetValuePreserve(0.);
            gP(i + 1).DerivativeResizeReset(&dof, i, MapVectorBase::LOCAL, 0.);
            gP(i + 1).SetDerivativeLocal(i, -1.); // derivative will be always -1
        }
    }

    void SetValue(blitz::Array<doublereal, 1>& XCurr, blitz::Array<doublereal, 1>& XPrimeCurr) {
        assert(iFirstDofIndex != -1);

        for (int i = 0; i < 3; ++i) {
            XCurr(iFirstDofIndex + i) = X(i + 1).dGetValue();
            XPrimeCurr(iFirstDofIndex + i) = XP(i + 1).dGetValue();
            XCurr(iFirstDofIndex + i + 3) = g(i + 1).dGetValue();
            XPrimeCurr(iFirstDofIndex + i + 3) = gP(i + 1).dGetValue();
        }
    }

    void Update(const blitz::Array<doublereal, 1>& XCurr, const blitz::Array<doublereal, 1>& XPrimeCurr, doublereal dCoef) {
        assert(iFirstDofIndex != -1);

        for (int i = 0; i < 3; ++i) {
            X(i + 1).SetValuePreserve(XCurr(iFirstDofIndex + i));
            X(i + 1).SetDerivativeLocal(i, -dCoef);
            XP(i + 1).SetValuePreserve(XPrimeCurr(iFirstDofIndex + i));
            g(i + 1).SetValuePreserve(XCurr(iFirstDofIndex + i + 3));
            g(i + 1).SetDerivativeLocal(i, -dCoef);
            gP(i + 1).SetValuePreserve(XPrimeCurr(iFirstDofIndex + i + 3));
        }

        UpdateRotation();
    }

    void UpdateRotation() {
        const Matrix<Gradient<NADVars>, 3, 3> skew_g(MatCrossVec(g));
        const Matrix<Gradient<NADVars>, 3, 3> skew_skew_g(MatCrossCrossVec(g));

        const Gradient<NADVars> d = 4. / (4. + Dot(g, g));
#if 0
        Matrix<Gradient<NADVars>, 3, 3> RDelta, G;
        for (int i = 1; i <= 3; ++i) {
            RDelta(i, i).SetValuePreserve(1.);
            G(i, i) = d;
        }

        for (int i = 1; i <= 3; ++i) {
            for (int j = 1; j <= 3; ++j) {
                RDelta(i, j) += d * (skew_g(i, j) + 0.5 * skew_skew_g(i, j));
                G(i, j) += 0.5 * d * skew_g(i, j);
            }
        }
#else
        Matrix<Gradient<NADVars>, 3, 3> RDelta = (skew_g + skew_skew_g * 0.5) * d;
        Matrix<Gradient<NADVars>, 3, 3> G = skew_g * Gradient<NADVars>(0.5 * d);

        for (int i = 1; i <= 3; ++i) {
            RDelta(i, i) += 1.;
            G(i, i) += d;
        }
#endif
        R = RDelta * R0;
        W = G * gP + RDelta * W0;
    }

    void AfterConvergence(const blitz::Array<doublereal, 1>& XCurr, const blitz::Array<doublereal, 1>& XPrimeCurr) {
        for (int i = 1; i <= 3; ++i) {
            W0(i) = W(i).dGetValue();
            g(i).SetValuePreserve(0.);
            gP(i).SetValuePreserve(0.);

            for (int j = 1; j <= 3; ++j) {
                R0(i, j) = R(i, j).dGetValue();
            }
        }
    }

    void SetFirstDofIndex(int iDofIndex) {
        iFirstDofIndex = iDofIndex;
    }

    int iGetFirstIndex() const {
        return iFirstDofIndex;
    }

    int iGetNumDof() const {
        return 6;
    }

    void GetXCurr(Vector<doublereal, 3>& XCurr, LocalDofMap*) const {
        for (int i = 1; i <= 3; ++i) {
            XCurr(i) = X(i).dGetValue();
        }
    }

    template <index_type N_SIZE>
    void GetXCurr(Vector<Gradient<N_SIZE>, 3>& XCurr, LocalDofMap* pDofMap) const {
        assert(iFirstDofIndex != -1);
        assert(pDofMap != 0);

        for (int i = 1; i <= 3; ++i) {
            XCurr(i).SetValuePreserve(X(i).dGetValue());
            XCurr(i).DerivativeResizeReset(pDofMap, iFirstDofIndex + X(i).iGetStartIndexLocal(), iFirstDofIndex + X(i).iGetEndIndexLocal(), MapVectorBase::GLOBAL, 0.);

            for (index_type j = X(i).iGetStartIndexLocal(); j < X(i).iGetEndIndexLocal(); ++j) {
            	XCurr(i).SetDerivativeGlobal(iFirstDofIndex + j, X(i).dGetDerivativeLocal(j));
            }
        }
    }

    void GetVCurr(Vector<doublereal, 3>& VCurr, LocalDofMap*) const {
        for (int i = 1; i <= 3; ++i) {
            VCurr(i) = XP(i).dGetValue();
        }
    }

    template <index_type N_SIZE>
    void GetVCurr(Vector<Gradient<N_SIZE>, 3>& VCurr, LocalDofMap* pDofMap) const {
        assert(iFirstDofIndex != -1);
        assert(pDofMap != 0);

        for (int i = 1; i <= 3; ++i) {
            VCurr(i).SetValuePreserve(XP(i).dGetValue());
            VCurr(i).DerivativeResizeReset(pDofMap, iFirstDofIndex + XP(i).iGetStartIndexLocal(), iFirstDofIndex + XP(i).iGetEndIndexLocal(), MapVectorBase::GLOBAL, 0.);

            for (index_type j = XP(i).iGetStartIndexLocal(); j < XP(i).iGetEndIndexLocal(); ++j) {
            	VCurr(i).SetDerivativeGlobal(iFirstDofIndex + j, XP(i).dGetDerivativeLocal(j));
            }
        }
    }

    void GetRCurr(Matrix<doublereal, 3, 3>& RCurr, LocalDofMap*) const {
        for (int i = 1; i <= 3; ++i) {
            for (int j = 1; j <= 3; ++j) {
                RCurr(i, j) = R(i, j).dGetValue();
            }
        }
    }

    template <index_type N_SIZE>
    void GetRCurr(Matrix<Gradient<N_SIZE>, 3, 3>& RCurr, LocalDofMap* pDofMap) const {
        assert(iFirstDofIndex != -1);
        assert(pDofMap != 0);

        for (int i = 1; i <= 3; ++i) {
            for (int j = 1; j <= 3; ++j) {
                RCurr(i, j).SetValuePreserve(R(i, j).dGetValue());
                RCurr(i, j).DerivativeResizeReset(pDofMap, iFirstDofIndex + R(i, j).iGetStartIndexLocal() + 3, iFirstDofIndex + R(i, j).iGetEndIndexLocal() + 3, MapVectorBase::GLOBAL, 0.);

                for (index_type k = R(i, j).iGetStartIndexLocal(); k < R(i, j).iGetEndIndexLocal(); ++k) {
                	RCurr(i, j).SetDerivativeGlobal(iFirstDofIndex + k + 3, R(i, j).dGetDerivativeLocal(k));
                }
            }
        }
    }

    const Matrix<doublereal, 3, 3>& GetRRef() const {
        return R0;
    }

    void GetWCurr(Vector<doublereal, 3>& WCurr, LocalDofMap*) const {
        for (int i = 1; i <= 3; ++i) {
            WCurr(i) = W(i).dGetValue();
        }
    }

    template <index_type N_SIZE>
    void GetWCurr(Vector<Gradient<N_SIZE>, 3>& WCurr, LocalDofMap* pDofMap) const {
        assert(iFirstDofIndex != -1);
        assert(pDofMap != 0);

        for (int i = 1; i <= 3; ++i) {
            WCurr(i).SetValuePreserve(W(i).dGetValue());
            WCurr(i).DerivativeResizeReset(pDofMap, iFirstDofIndex + W(i).iGetStartIndexLocal() + 3, iFirstDofIndex + W(i).iGetEndIndexLocal() + 3, MapVectorBase::GLOBAL, 0.);

            for (index_type j = W(i).iGetStartIndexLocal(); j < W(i).iGetEndIndexLocal(); ++j) {
            	WCurr(i).SetDerivativeGlobal(iFirstDofIndex + j + 3, W(i).dGetDerivativeLocal(j));
            }
        }
    }

    const Vector<doublereal, 3>& GetWRef() const {
        return W0;
    }

private:
    int iFirstDofIndex;
    Matrix<doublereal, 3, 3> R0;
    Vector<doublereal, 3> W0;
    static const int NADVars = 3;
    Vector<Gradient<NADVars>, 3> X, XP, g, gP, W;
    Matrix<Gradient<NADVars>, 3, 3> R;
    LocalDofMap dof;
};

template <typename T>
struct ResItem {
    int iEquIndex;
    T dCoef;

    ResItem(int iEquIndex_=-1, T dCoef_=T(0.))
        :iEquIndex(iEquIndex_), dCoef(dCoef_) {
    }
};

class FullSubMatrixHandler {
public:
    FullSubMatrixHandler(index_type iNumRows=0, index_type iNumCols=0) {
        ResizeReset(iNumRows, iNumCols);
    }

    void ResizeReset(index_type iNumRows, index_type iNumCols) {
        oWorkMat.resize(iNumRows, iNumCols);
        oRowIndex.resize(iNumRows);
        oColIndex.resize(iNumCols);
        oWorkMat.initialize(0.);
        oRowIndex.initialize(-1);
        oColIndex.initialize(-1);
    }

    void PutRowIndex(int iSubRow, int iRow) {
        assert(iSubRow < oRowIndex.rows());
        oRowIndex(iSubRow) = iRow;
    }

    void PutColIndex(int iSubCol, int iCol) {
        assert(iSubCol < oColIndex.rows());
        oColIndex(iSubCol) = iCol;
    }

    void PutCoef(int iSubRow, int iSubCol, doublereal dCoef) {
        assert(iSubRow >= 0 && iSubRow < oWorkMat.rows());
        assert(iSubCol >= 0 && iSubCol < oWorkMat.cols());

        oWorkMat(iSubRow, iSubCol) = dCoef;
    }

    void IncCoef(int iSubRow, int iSubCol, doublereal dCoef) {
        assert(iSubRow >= 0 && iSubRow < oWorkMat.rows());
        assert(iSubCol >= 0 && iSubCol < oWorkMat.cols());

        oWorkMat(iSubRow, iSubCol) += dCoef;
    }

    void AddTo(blitz::Array<doublereal, 2>& JacMat) const {
        for (int i = 0; i < oWorkMat.rows(); ++i) {
            for (int j = 0; j < oWorkMat.cols(); ++j) {
                assert(oRowIndex(i) >= 0);
                assert(oRowIndex(i) < JacMat.rows());
                assert(oColIndex(j) >= 0);
                assert(oColIndex(j) < JacMat.cols());
                JacMat(oRowIndex(i), oColIndex(j)) += oWorkMat(i, j);
            }
        }
    }

private:
    blitz::Array<doublereal, 2> oWorkMat;
    blitz::Array<int, 1> oRowIndex, oColIndex;
};

class SparseSubMatrixHandler {
public:
    struct JacItem {
        int iEquIndex;
        int iDofIndex;
        doublereal dCoef;

        JacItem(int iEquIndex=-1, int iDofIndex=-1, doublereal dCoef=0.)
            :iEquIndex(iEquIndex), iDofIndex(iDofIndex), dCoef(dCoef) {
        }
    };

    typedef std::vector<JacItem> VectorType;
    typedef VectorType::const_iterator const_iterator;

    SparseSubMatrixHandler(int iNumItems=0) {
        if (iNumItems > 0) {
            WorkMat.reserve(iNumItems);
        }
    }

    template <index_type N_SIZE, typename T>
    SparseSubMatrixHandler& AssJac(T* pElem, LocalDofMap* pDofMap, blitz::Array<ResItem<Gradient<N_SIZE> >, 1>& WorkVec, doublereal dCoef) {
        pElem->AssRes(WorkVec, dCoef, pDofMap);
        ResizeReset(0);

        for (int i = 0; i < WorkVec.rows(); ++i) {
            const ResItem<Gradient<N_SIZE> >& resItem = WorkVec(i);

            for (index_type j = resItem.dCoef.iGetStartIndexLocal(); j < resItem.dCoef.iGetEndIndexLocal(); ++j) {
                index_type iDofIndex = resItem.dCoef.iGetGlobalDof(j);
                InsertItem(JacItem(resItem.iEquIndex, iDofIndex, resItem.dCoef.dGetDerivativeLocal(j)));
            }
        }

        return *this;
    }

    void ResizeReset(int iNumItems) {
        WorkMat.resize(iNumItems);
    }

    int iGetSize() const { return WorkMat.size(); }

    void InsertItem(const JacItem& item) {
        WorkMat.push_back(item);
    }

    void AddTo(blitz::Array<doublereal, 2>& JacMat) const {
        for (const_iterator j = begin(); j != end(); ++j) {
            JacMat(j->iEquIndex, j->iDofIndex) += j->dCoef;
        }
    }
    const_iterator begin() const { return WorkMat.begin(); }
    const_iterator end() const { return WorkMat.end(); }

private:
    VectorType WorkMat;
};

class Element {
public:
    virtual blitz::Array<ResItem<doublereal>, 1>& AssRes(blitz::Array<ResItem<doublereal>, 1>& WorkVec, doublereal dCoef)=0;
    virtual SparseSubMatrixHandler& AssJac(SparseSubMatrixHandler& WorkMat, doublereal dCoef)=0;
    virtual FullSubMatrixHandler& AssJac(FullSubMatrixHandler& WorkMat, doublereal dCoef)=0;
    virtual index_type iGetNumRows() const=0;
    virtual index_type iGetNumCols() const=0;
    virtual ~Element(){ }
};

class Element1: public Element {
private:
    Node* node1;
    Node* node2;
    Vector<doublereal, 3> o1, o2;
    Matrix<doublereal, 3, 3> S, D;
    static const int NADVars = 12;
    LocalDofMap dofMap;

public:
    Element1(Node* node1_,
             const Vector<doublereal, 3>& o1_,
             Node* node2_,
             const Vector<doublereal, 3>& o2_,
             doublereal s,
             doublereal d)
        :node1(node1_),
         node2(node2_),
         o1(o1_),
         o2(o2_),
         dofMap(iGetNumCols()) {

        /*
            S=[ s,  -s,    0;
               -s, 2*s,   -s;
                0,  -s,  2*s];
        */

        S(1, 1) = s;
        S(2, 1) = -s;
        S(1, 2) = -s;
        S(2, 2) = 2*s;
        S(3, 2) = -s;
        S(2, 3) = -s;
        S(3, 3) = 2*s;

        /*
            D=[ d,     -d,      0;
               -d,  2 * d,     -d;
                0,     -d,  2 * d];
        */

        D(1, 1) = d;
        D(2, 1) = -d;
        D(1, 2) = -d;
        D(2, 2) = 2*d;
        D(3, 2) = -d;
        D(2, 3) = -d;
        D(3, 3) = 2*d;
    }

    template <typename T>
    blitz::Array<ResItem<T>, 1>& AssRes(blitz::Array<ResItem<T>, 1>& WorkVec, doublereal dCoef, LocalDofMap *pDofMap) {
        WorkVec.resize(iGetNumRows());
        typedef Vector<T, 3> Vec3;
        typedef Matrix<T, 3, 3> Mat3x3;
        Vec3 X1, X2, V1, V2, W1, W2;
        Mat3x3 R1, R2;

        node1->GetXCurr(X1, pDofMap);
        node1->GetVCurr(V1, pDofMap);
        node1->GetRCurr(R1, pDofMap);
        node1->GetWCurr(W1, pDofMap);
        node2->GetXCurr(X2, pDofMap);
        node2->GetVCurr(V2, pDofMap);
        node2->GetRCurr(R2, pDofMap);
        node2->GetWCurr(W2, pDofMap);

        const Vec3 R1o1 = R1 * o1;
        const Vec3 R2o2 = R2 * o2;
        const Vec3 dX = Transpose(R1) * Vec3(X1 + R1o1 - X2 - R2o2);
        const Vec3 dV = Transpose(R1) * Vec3(V1 + Cross(W1, R1o1) - V2 - Cross(W2, R2o2));
        const Vec3 F1 = R1 * Vec3(-S * dX - D * dV);
        const Vec3 M1 = Cross(R1o1, F1), M2 = Cross(R2o2, -F1);

        for (int i = 0; i < 6; ++i) {
            WorkVec(i).iEquIndex = node1->iGetFirstIndex() + i;
            WorkVec(i + 6).iEquIndex = node2->iGetFirstIndex() + i;
        }

        for (int i = 0; i < 3; ++i) {
            WorkVec(i).dCoef = F1(i + 1);
            WorkVec(i + 3).dCoef = M1(i + 1);
            WorkVec(i + 6).dCoef = -F1(i + 1);
            WorkVec(i + 9).dCoef = M2(i + 1);
        }

        return WorkVec;
    }

    virtual blitz::Array<ResItem<doublereal>, 1>& AssRes(blitz::Array<ResItem<doublereal>, 1>& WorkVec, doublereal dCoef) {
        return AssRes(WorkVec, dCoef, 0);
    }

    virtual SparseSubMatrixHandler& AssJac(SparseSubMatrixHandler& WorkMat, doublereal dCoef) {
        blitz::Array<ResItem<Gradient<NADVars> >, 1> WorkVec;
        return WorkMat.AssJac(this, &dofMap, WorkVec, dCoef);
    }

    virtual FullSubMatrixHandler& AssJac(FullSubMatrixHandler& WorkMat, doublereal dCoef) {
        WorkMat.ResizeReset(iGetNumRows(), iGetNumCols());
        typedef Matrix<doublereal, 3, 3> Mat3x3;
        typedef Vector<doublereal, 3> Vec3;

        for (int i = 0; i < 6; ++i) {
            WorkMat.PutColIndex(i, node1->iGetFirstIndex() + i);
            WorkMat.PutColIndex(i + 6, node2->iGetFirstIndex() + i);
            WorkMat.PutRowIndex(i, node1->iGetFirstIndex() + i);
            WorkMat.PutRowIndex(i + 6, node2->iGetFirstIndex() + i);
        }

        const Vec3& W1_0 = node1->GetWRef();
        const Vec3& W2_0 = node2->GetWRef();
        const Mat3x3& R1_0 = node1->GetRRef();
        const Mat3x3& R2_0 = node2->GetRRef();

        Vec3 X1, X2, V1, V2, W1, W2;
        Mat3x3 R1, R2;

        node1->GetXCurr(X1, 0);
        node1->GetVCurr(V1, 0);
        node1->GetRCurr(R1, 0);
        node1->GetWCurr(W1, 0);

        node2->GetXCurr(X2, 0);
        node2->GetVCurr(V2, 0);
        node2->GetRCurr(R2, 0);
        node2->GetWCurr(W2, 0);

#ifdef ASS_JAC_USE_TEMP_EXPR
        const Mat3x3 skew_W2_0(MatCrossVec(W2_0));
        const Vec3 R1o1 = Vec3(R1 * o1);
        const Mat3x3 skew_R1o1(MatCrossVec(R1o1));
        const Vec3 R1_0o1 = Vec3(R1_0 * o1);
        const Mat3x3 skew_R1_0o1(MatCrossVec(R1_0o1));
        const Vec3 R2o2 = Vec3(R2 * o2);
        const Mat3x3 skew_R2o2(MatCrossVec(R2o2));
        const Vec3 R2_0o2 = Vec3(R2_0 * o2);
        const Mat3x3 skew_R2_0o2(MatCrossVec(R2_0o2));
        const Vec3 dX = Vec3(Mat3x3(Transpose(R1)) * Vec3(Vec3(Vec3(X1 + R1o1) - X2) - R2o2));
        const Vec3 dV = Vec3(Mat3x3(Transpose(R1)) * Vec3(Vec3(Vec3(V1 + Vec3(Cross(W1, R1o1))) - V2) - Vec3(Cross(W2, R2o2))));
        const Vec3 F1_R1 = Vec3(-Vec3(Vec3(S * dX) + Vec3(D * dV)));
        const Vec3 F1 = Vec3(R1 * F1_R1);
        const Vec3 F2 = Vec3(-F1);

        const Mat3x3 dF1_dX1 = Mat3x3(-Mat3x3(R1 * Mat3x3(S * Transpose(R1))));

        const Mat3x3 ddX_dg1 = Mat3x3(Mat3x3(Mat3x3(Transpose(R1_0)) * Mat3x3(MatCrossVec(Vec3(Vec3(Vec3(X1 + R1o1) - X2) - R2o2)))) - Mat3x3(Mat3x3(Transpose(R1)) * skew_R1_0o1));
        const Mat3x3 ddV_dg1 = Mat3x3(Mat3x3(Transpose(R1_0)) * Mat3x3(MatCrossVec(Vec3(Vec3(Vec3(V1 + Vec3(Cross(W1, R1o1))) - V2) - Vec3(Cross(W2, R2o2))))))
                                               + Mat3x3(Mat3x3(Transpose(R1)) * Mat3x3(Mat3x3(skew_R1o1 * Mat3x3(MatCrossVec(W1_0))) - Mat3x3(Mat3x3(MatCrossVec(W1)) * skew_R1_0o1)));
        const Mat3x3 dF1_dg1 = Mat3x3(Mat3x3(MatCrossVec(Vec3(R1_0 * Vec3(-F1_R1)))) - Mat3x3(R1 * Mat3x3(Mat3x3(S * ddX_dg1) + Mat3x3(D * ddV_dg1))));

        const Mat3x3 dF1_dX2 = Mat3x3(R1 * Mat3x3(S * Mat3x3(Transpose(R1))));
        const Mat3x3 ddX_dg2 = Mat3x3(Mat3x3(Transpose(R1)) * skew_R2_0o2);
        const Mat3x3 ddV_dg2 = Mat3x3(Mat3x3(Transpose(R1)) * Mat3x3(Mat3x3(skew_R2o2 * Mat3x3(-skew_W2_0)) + Mat3x3(skew_W2_0 * skew_R2_0o2)));
        const Mat3x3 dF1_dg2 = Mat3x3(Mat3x3(-R1) * Mat3x3(Mat3x3(S * ddX_dg2) + Mat3x3(D * ddV_dg2)));

        const Mat3x3 dF2_dX1 = Mat3x3(-dF1_dX1);
        const Mat3x3 dF2_dg1 = Mat3x3(-dF1_dg1);
        const Mat3x3 dF2_dX2 = Mat3x3(-dF1_dX2);
        const Mat3x3 dF2_dg2 = Mat3x3(-dF1_dg2);

        const Mat3x3 dM1_dX1 = Mat3x3(skew_R1o1 * dF1_dX1);
        const Mat3x3 dM1_dg1 = Mat3x3(Mat3x3(Mat3x3(MatCrossVec(F1)) * skew_R1_0o1) + Mat3x3(skew_R1o1 * dF1_dg1));
        const Mat3x3 dM1_dX2 = Mat3x3(skew_R1o1 * dF1_dX2);
        const Mat3x3 dM1_dg2 = Mat3x3(skew_R1o1 * dF1_dg2);

        const Mat3x3 dM2_dX1 = Mat3x3(skew_R2o2 * dF2_dX1);
        const Mat3x3 dM2_dg1 = Mat3x3(skew_R2o2 * dF2_dg1);
        const Mat3x3 dM2_dX2 = Mat3x3(skew_R2o2 * dF2_dX2);
        const Mat3x3 dM2_dg2 = Mat3x3(Mat3x3(Mat3x3(MatCrossVec(F2)) * skew_R2_0o2) + Mat3x3(skew_R2o2 * dF2_dg2));

        const Mat3x3 dF1_dV1 = Mat3x3(R1 * Mat3x3(Mat3x3(-D) * Mat3x3(Transpose(R1))));
        const Mat3x3 ddV_dgP1 = Mat3x3(Mat3x3(-Mat3x3(Transpose(R1))) * skew_R1o1);
        const Mat3x3 dF1_dgP1 = Mat3x3(R1 * Mat3x3(Mat3x3(-D) * ddV_dgP1));
        const Mat3x3 dF1_dV2 = Mat3x3(R1 * Mat3x3(D * Mat3x3(Transpose(R1))));
        const Mat3x3 ddV_dgP2 = Mat3x3(Mat3x3(Transpose(R1)) * skew_R2o2);
        const Mat3x3 dF1_dgP2 = Mat3x3(R1 * Mat3x3(Mat3x3(-D) * ddV_dgP2));

        const Mat3x3 dM1_dV1 = Mat3x3(skew_R1o1 * dF1_dV1);
        const Mat3x3 dM1_dgP1 = Mat3x3(skew_R1o1 * dF1_dgP1);
        const Mat3x3 dM1_dV2 = Mat3x3(skew_R1o1 * dF1_dV2);
        const Mat3x3 dM1_dgP2 = Mat3x3(skew_R1o1 * dF1_dgP2);

        const Mat3x3 dF2_dV1 = Mat3x3(-dF1_dV1);
        const Mat3x3 dF2_dgP1 = Mat3x3(-dF1_dgP1);
        const Mat3x3 dF2_dV2 = Mat3x3(-dF1_dV2);
        const Mat3x3 dF2_dgP2 = Mat3x3(-dF1_dgP2);

        const Mat3x3 dM2_dV1 = Mat3x3(skew_R2o2 * dF2_dV1);
        const Mat3x3 dM2_dgP1 = Mat3x3(skew_R2o2 * dF2_dgP1);
        const Mat3x3 dM2_dV2 = Mat3x3(skew_R2o2 * dF2_dV2);
        const Mat3x3 dM2_dgP2 = Mat3x3(skew_R2o2 * dF2_dgP2);
#else
        const Mat3x3 skew_W2_0(MatCrossVec(W2_0));
        const Vec3 R1o1 = R1 * o1;
        const Mat3x3 skew_R1o1(MatCrossVec(R1o1));
        const Vec3 R1_0o1 = R1_0 * o1;
        const Mat3x3 skew_R1_0o1(MatCrossVec(R1_0o1));
        const Vec3 R2o2 = R2 * o2;
        const Mat3x3 skew_R2o2(MatCrossVec(R2o2));
        const Vec3 R2_0o2 = R2_0 * o2;
        const Mat3x3 skew_R2_0o2(MatCrossVec(R2_0o2));
        const Vec3 dX = Transpose(R1) * Vec3(X1 + R1o1 - X2 - R2o2);
        const Vec3 dV = Transpose(R1) * Vec3(V1 + Cross(W1, R1o1) - V2 - Cross(W2, R2o2));
        const Vec3 F1_R1 = -(S * dX + D * dV);
        const Vec3 F1 = R1 * F1_R1;
        const Vec3 F2 = -F1;

        const Mat3x3 dF1_dX1 = -R1 * Mat3x3(S * Transpose(R1));

        const Mat3x3 ddX_dg1 = Transpose(R1_0) * Mat3x3(MatCrossVec(X1 + R1o1 - X2 - R2o2)) - Transpose(R1) * skew_R1_0o1;
        const Mat3x3 ddV_dg1 = Transpose(R1_0) * Mat3x3(MatCrossVec(V1 + Cross(W1, R1o1) - V2 - Cross(W2, R2o2)))
                                               + Transpose(R1) * Mat3x3(skew_R1o1 * Mat3x3(MatCrossVec(W1_0)) - Mat3x3(MatCrossVec(W1)) * skew_R1_0o1);
        const Mat3x3 dF1_dg1 = Mat3x3(MatCrossVec(R1_0 * (-F1_R1))) - R1 * Mat3x3(S * ddX_dg1 + D * ddV_dg1);

        const Mat3x3 dF1_dX2 = R1 * Mat3x3(S * Transpose(R1));
        const Mat3x3 ddX_dg2 = Transpose(R1) * skew_R2_0o2;
        const Mat3x3 ddV_dg2 = Transpose(R1) * Mat3x3(skew_R2o2 * (-skew_W2_0) + skew_W2_0 * skew_R2_0o2);
        const Mat3x3 dF1_dg2 = -R1 * Mat3x3(S * ddX_dg2 + D * ddV_dg2);

        const Mat3x3 dF2_dX1 = -dF1_dX1;
        const Mat3x3 dF2_dg1 = -dF1_dg1;
        const Mat3x3 dF2_dX2 = -dF1_dX2;
        const Mat3x3 dF2_dg2 = -dF1_dg2;

        const Mat3x3 dM1_dX1 = skew_R1o1 * dF1_dX1;
        const Mat3x3 dM1_dg1 = Mat3x3(MatCrossVec(F1)) * skew_R1_0o1 + skew_R1o1 * dF1_dg1;
        const Mat3x3 dM1_dX2 = skew_R1o1 * dF1_dX2;
        const Mat3x3 dM1_dg2 = skew_R1o1 * dF1_dg2;

        const Mat3x3 dM2_dX1 = skew_R2o2 * dF2_dX1;
        const Mat3x3 dM2_dg1 = skew_R2o2 * dF2_dg1;
        const Mat3x3 dM2_dX2 = skew_R2o2 * dF2_dX2;
        const Mat3x3 dM2_dg2 = Mat3x3(MatCrossVec(F2)) * skew_R2_0o2 + skew_R2o2 * dF2_dg2;

        const Mat3x3 dF1_dV1 = R1 * Mat3x3((-D) * Transpose(R1));
        const Mat3x3 ddV_dgP1 = -Transpose(R1) * skew_R1o1;
        const Mat3x3 dF1_dgP1 = R1 * Mat3x3((-D) * ddV_dgP1);
        const Mat3x3 dF1_dV2 = R1 * Mat3x3(D * Transpose(R1));
        const Mat3x3 ddV_dgP2 = Transpose(R1) * skew_R2o2;
        const Mat3x3 dF1_dgP2 = R1 * Mat3x3((-D) * ddV_dgP2);

        const Mat3x3 dM1_dV1 = skew_R1o1 * dF1_dV1;
        const Mat3x3 dM1_dgP1 = skew_R1o1 * dF1_dgP1;
        const Mat3x3 dM1_dV2 = skew_R1o1 * dF1_dV2;
        const Mat3x3 dM1_dgP2 = skew_R1o1 * dF1_dgP2;

        const Mat3x3 dF2_dV1 = -dF1_dV1;
        const Mat3x3 dF2_dgP1 = -dF1_dgP1;
        const Mat3x3 dF2_dV2 = -dF1_dV2;
        const Mat3x3 dF2_dgP2 = -dF1_dgP2;

        const Mat3x3 dM2_dV1 = skew_R2o2 * dF2_dV1;
        const Mat3x3 dM2_dgP1 = skew_R2o2 * dF2_dgP1;
        const Mat3x3 dM2_dV2 = skew_R2o2 * dF2_dV2;
        const Mat3x3 dM2_dgP2 = skew_R2o2 * dF2_dgP2;
#endif

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                WorkMat.PutCoef(i, j,     -dF1_dV1(i + 1, j + 1)  - dCoef * dF1_dX1(i + 1, j + 1));
                WorkMat.PutCoef(i, j + 3, -dF1_dgP1(i + 1, j + 1) - dCoef * dF1_dg1(i + 1, j + 1));
                WorkMat.PutCoef(i, j + 6, -dF1_dV2(i + 1, j + 1)  - dCoef * dF1_dX2(i + 1, j + 1));
                WorkMat.PutCoef(i, j + 9, -dF1_dgP2(i + 1, j + 1) - dCoef * dF1_dg2(i + 1, j + 1));

                WorkMat.PutCoef(i + 3, j,     -dM1_dV1(i + 1, j + 1)  - dCoef * dM1_dX1(i + 1, j + 1));
                WorkMat.PutCoef(i + 3, j + 3, -dM1_dgP1(i + 1, j + 1) - dCoef * dM1_dg1(i + 1, j + 1));
                WorkMat.PutCoef(i + 3, j + 6, -dM1_dV2(i + 1, j + 1)  - dCoef * dM1_dX2(i + 1, j + 1));
                WorkMat.PutCoef(i + 3, j + 9, -dM1_dgP2(i + 1, j + 1) - dCoef * dM1_dg2(i + 1, j + 1));

                WorkMat.PutCoef(i + 6, j,     -dF2_dV1(i + 1, j + 1)  - dCoef * dF2_dX1(i + 1, j + 1));
                WorkMat.PutCoef(i + 6, j + 3, -dF2_dgP1(i + 1, j + 1) - dCoef * dF2_dg1(i + 1, j + 1));
                WorkMat.PutCoef(i + 6, j + 6, -dF2_dV2(i + 1, j + 1)  - dCoef * dF2_dX2(i + 1, j + 1));
                WorkMat.PutCoef(i + 6, j + 9, -dF2_dgP2(i + 1, j + 1)  - dCoef * dF2_dg2(i + 1, j + 1));

                WorkMat.PutCoef(i + 9, j,     -dM2_dV1(i + 1, j + 1)  - dCoef * dM2_dX1(i + 1, j + 1));
                WorkMat.PutCoef(i + 9, j + 3, -dM2_dgP1(i + 1, j + 1) - dCoef * dM2_dg1(i + 1, j + 1));
                WorkMat.PutCoef(i + 9, j + 6, -dM2_dV2(i + 1, j + 1)  - dCoef * dM2_dX2(i + 1, j + 1));
                WorkMat.PutCoef(i + 9, j + 9, -dM2_dgP2(i + 1, j + 1) - dCoef * dM2_dg2(i + 1, j + 1));
            }
        }

        return WorkMat;
    }

    index_type iGetNumRows() const { return 12; }
    index_type iGetNumCols() const { return NADVars; }
};

void testAssembly() {
    doublereal tckRes = 0;
    doublereal tckJacAD = 0;
    doublereal tckJac = 0;
    doublereal tckStart;
    long iIterCnt = 0;
    tic(tckStart);

    for (int loop = 0; loop < NLoopsAss; ++loop) {
        const int iNumNodes = 3;

        Node* nodes[iNumNodes] = { 0 };

        for (int i = 0; i < iNumNodes; ++i) {
            Vector<doublereal, 3> X0, XP0, Phi0, W0;

            for (int j = 0; j < 3; ++j) {
                X0(j + 1) = ((i + 1) * 10 + j + 1);
                XP0(j + 1) = ((i + 1) * 1000 + (j + 1) * 100);
                Phi0(j + 1) = ((i + 1) * 0.1 + (j + 1) * 0.01);
                W0(j + 1) = ((i + 1) * 0.1 + (j + 1) * 0.01);
            }

            Matrix<doublereal, 3, 3> R0;

            Euler123ToMatR(Phi0, R0);

            nodes[i] = new Node(X0, XP0, R0, W0);
        }

        int iNumDof = 0;

        for (int i = 0; i < iNumNodes; ++i) {
            nodes[i]->SetFirstDofIndex(iNumDof);
            iNumDof += nodes[i]->iGetNumDof();
        }

        const int iNumElem = iNumNodes - 1;

        Element* elements[iNumElem] = {0};

        for (int i = 0; i < iNumElem; ++i) {
            Vector<doublereal, 3> o1, o2;

            o1(1) = 1.;
            o1(2) = 2.;
            o1(3) = 3.;
            o2(1) = 4.;
            o2(2) = 5.;
            o2(3) = 6.;

            doublereal s = 100 * (i + 1);
            doublereal d = 10 * (i + 1);
            elements[i] = new Element1(nodes[i], o1, nodes[i + 1], o2, s, d);
        }

        blitz::Array<doublereal, 1> XCurr(iNumDof), XPrimeCurr(iNumDof);

        XCurr.initialize(0.);
        XPrimeCurr.initialize(0.);

        for (int i = 0; i < iNumNodes; ++i) {
            nodes[i]->SetValue(XCurr, XPrimeCurr);
        }

        blitz::Array<ResItem<doublereal>, 1> WorkVec;
        SparseSubMatrixHandler WorkMatSp;
        FullSubMatrixHandler WorkMatFull;
        blitz::Array<doublereal, 1> ResVec;
        blitz::Array<doublereal, 2> JacMatAD, JacMat;

        ResVec.resize(iNumDof);
        JacMatAD.resize(iNumDof, iNumDof);
        JacMat.resize(iNumDof, iNumDof);

        ResVec.initialize(0.);
        JacMatAD.initialize(0.);
        JacMat.initialize(0.);

        const doublereal dCoef = 0.001;

        for (int iTime = 0; iTime < 100; ++iTime) {
        	iIterCnt++;

            for (int i = 0; i < iNumNodes; ++i) {
                nodes[i]->Update(XCurr, XPrimeCurr, dCoef);
            }

            for (int i = 0; i < iNumElem; ++i) {
                tic();

                elements[i]->AssRes(WorkVec, dCoef);

                tckRes += toc();

                for (int j = 0; j < WorkVec.rows(); ++j) {
                    ResVec(WorkVec(j).iEquIndex) += WorkVec(j).dCoef;
                }

                tic();

                elements[i]->AssJac(WorkMatSp, dCoef);

                tckJacAD += toc();

                WorkMatSp.AddTo(JacMatAD);

                tic();
                elements[i]->AssJac(WorkMatFull, dCoef);
                tckJac += toc();

                WorkMatFull.AddTo(JacMat);
            }

            for (int i = 0; i < JacMat.rows(); ++i) {
                for (int j = 0; j < JacMat.cols(); ++j) {
                    const doublereal dTol = sqrt(std::numeric_limits<scalar_deriv_type>::epsilon()) * std::max(1., std::max(std::abs(JacMat(i, j)), std::abs(JacMatAD(i, j))));
                    if(std::abs(JacMat(i, j) - JacMatAD(i, j)) > dTol) {
                    	throw std::runtime_error("testAssembly(): incorrect result");
                    }
                }
            }
        }

        if (loop == 0) {
            std::cout << "ResVec = [" << std::endl;

            for (int i = 0; i < iNumDof; ++i) {
                std::cout << std::setw(10) << ResVec(i) << std::endl;
            }

            std::cout << "]" << std::endl;

            std::cout << "JacMatAD = [" << std::endl;

            for (int i = 0; i < iNumDof; ++i) {
                for (int j = 0; j < iNumDof; ++j) {
                    std::cout << std::setw(10) << std::setprecision(16) << JacMatAD(i, j) << " ";
                }
                std::cout << std::endl;
            }

            std::cout << "]" << std::endl;

            std::cout << "JacMat = [" << std::endl;

            for (int i = 0; i < iNumDof; ++i) {
                for (int j = 0; j < iNumDof; ++j) {
                    std::cout << std::setw(10) << std::setprecision(16) << JacMat(i, j) << " ";
                }
                std::cout << std::endl;
            }

            std::cout << "]" << std::endl;
        }

        for (int i = 0; i < iNumElem; ++i) {
            delete elements[i];
        }

        for (int i = 0; i < iNumNodes; ++i) {
            delete nodes[i];
        }
    }

    doublereal tckEnd;
    tic(tckEnd);
    std::cerr << "number of iterations:" << iIterCnt << std::endl;
    std::cerr << "AssRes:" << tckRes / iIterCnt << std::endl;
    std::cerr << "AssJacAD:" << tckJacAD / iIterCnt << std::endl;
    std::cerr << "AssJac:" << tckJac / iIterCnt << std::endl;
    std::cerr << "overhead:" << tckJacAD / tckJac << std::endl;
    std::cerr << "total time:" << (tckEnd - tckStart) / iIterCnt << std::endl;
}

#endif
}

void testMatVec3() {
	srand(0);

	for (int n = 0; n < NLoops; ++n) {
		Vector<doublereal, 3> g1, v1;
		Vec3 g2, v2;
		for (index_type i = 1; i <= 3; ++i) {
			g2(i) = g1(i) = random1() * 1e-1;
			v2(i) = v1(i) = random1() * 1e1;
		}

		Matrix<doublereal, 3, 3> R1(MatGVec(g1));
		Mat3x3 R2(CGR_Rot::MatG, g2);
		Matrix<doublereal, 3, 3> C1(MatCrossVec(g1));
		Mat3x3 C2(MatCross, g2);
		Matrix<doublereal, 3, 3> CC1(MatCrossCrossVec(g1));
		Mat3x3 CC2(MatCrossCross, g2, g2);
		Vector<doublereal, 3> CCv1 = CC1 * v1;
		Vec3 CCv2 = CC2 * v2;
		Vector<doublereal, 3> CCv3 = CC1 * v2;
		Vector<doublereal, 3> CCv4 = CC2 * v1;
		Matrix<doublereal, 3, 3> A1 = R1 * C1;
		Mat3x3 A2 = R2 * C2;
		Matrix<doublereal, 3, 3> A3 = R1 * C2;
		Matrix<doublereal, 3, 3> A4 = R2 * C1;
		Vector<doublereal, 3> v3(v1(1), v1(2), v1(3));
		Vector<doublereal, 2> v4(v1(1), v1(2));
        Matrix<doublereal, 3, 3> X1(MatCrossVec(g1, 1.));
        Mat3x3 X2(1., g2);

		if (n == 0) {
			std::cout << "g1=" << std::endl << g1 << std::endl;
			std::cout << "g2=" << std::endl << g2 << std::endl;
			std::cout << "R1=" << std::endl << R1 << std::endl;
			std::cout << "R2=" << std::endl << R2 << std::endl;
			std::cout << "C1=" << std::endl << C1 << std::endl;
			std::cout << "C2=" << std::endl << C2 << std::endl;
			std::cout << "CC1=" << std::endl << CC1 << std::endl;
			std::cout << "CC2=" << std::endl << CC2 << std::endl;
			std::cout << "CCv1=" << std::endl << CCv1 << std::endl;
			std::cout << "CCv2=" << std::endl << CCv2 << std::endl;
			std::cout << "CCv3=" << std::endl << CCv3 << std::endl;
			std::cout << "CCv4=" << std::endl << CCv4 << std::endl;
			std::cout << "A1=" << std::endl << A1 << std::endl;
			std::cout << "A2=" << std::endl << A2 << std::endl;
			std::cout << "A3=" << std::endl << A3 << std::endl;
			std::cout << "A4=" << std::endl << A4 << std::endl;
		}

        const doublereal dTol = 10*std::numeric_limits<scalar_deriv_type>::epsilon();

		for (index_type i = 1; i <= 3; ++i) {
			for (index_type j = 1; j <= 3; ++j) {
				assert(std::abs(R1(i, j) - R2(i, j)) < dTol);
				assert(std::abs(C1(i, j) - C2(i, j)) < dTol);
				assert(std::abs(CC1(i, j) - CC2(i, j)) < dTol);
				assert(std::abs(A1(i, j) - A2(i, j)) < dTol);
				assert(std::abs(A1(i, j) - A3(i, j)) < dTol);
				assert(std::abs(A1(i, j) - A4(i, j)) < dTol);
                assert(std::abs(X1(i, j) - X2(i, j)) < dTol);
			}
			assert(std::abs(CCv1(i) - CCv2(i)) < dTol);
			assert(std::abs(CCv3(i) - CCv1(i)) < dTol);
			assert(std::abs(CCv4(i) - CCv1(i)) < dTol);
			assert(v3(i) == v1(i));

			if (i < 3) {
				assert(v4(i) == v1(i));
			}
		}
	}
}

void testSubVecAss() {
	LocalDofMap dof;
	const integer N = 3;
	MySubVectorHandler vh(N);

	GradientAssVec<doublereal> WorkVec(vh);

	for (integer i = 1; i <= N; ++i) {
		WorkVec.AddItem(i, i * 10.);
	}

	SparseSubMatrixHandler mh(4*N);
	GradientAssVec<Gradient<4> > WorkMat(mh);

	for (integer i = 1; i <= N; ++i) {
		Gradient<4> g;
		g.SetValue(i * 10.);
		g.DerivativeResizeReset(&dof, 1, 5, MapVectorBase::GLOBAL, 0.);
		for (integer k = 1; k <= 4; ++k) {
			g.SetDerivativeGlobal(k, i + k * 0.1);
		}

		WorkMat.AddItem(i, g);
	}

	std::cout << "WorkVec=" << std::endl;
	for (integer i = 1; i <= vh.iGetSize(); ++i) {
		std::cout << " " << vh.dGetCoef(i) << std::endl;
	}

	std::cout << std::endl;

	std::cout << "WorkMat=" << std::endl;

	for (integer i = 1; i <= mh.iGetNumRows(); ++i) {
		std::cout << " " << mh.iGetRowIndex(i)
				  << " " << mh.iGetColIndex(i)
				  << " " << mh.dGetCoef(i, 0) << std::endl;
	}

	std::cout << std::endl;
}

void testSubVecAssMatVec() {
	LocalDofMap dof;
	const integer N = 3;
	MySubVectorHandler vh(2*N);

	Vector<doublereal, 3> v;

	for (integer i = 1; i <= 3; ++i) {
		v(i) = i * 10;
	}

	GradientAssVec<doublereal> WorkVec(vh);

	WorkVec.AddItem(1, v);

	SparseSubMatrixHandler mh(2*4*N);
	GradientAssVec<Gradient<4> > WorkMat(mh);

	Vector<Gradient<4>, 3> g;

	for (integer i = 1; i <= 3; ++i) {
		g(i).SetValue(i * 10.);
		g(i).DerivativeResizeReset(&dof, 1, 5, MapVectorBase::GLOBAL, 0.);
		for (integer k = 1; k <= 4; ++k) {
			g(i).SetDerivativeGlobal(k, i + k * 0.1);
		}
	}

	WorkMat.AddItem(1, g);

	GradientAssVec<doublereal> WorkVec2(vh, GradientAssVecBase::APPEND);
	Vector<doublereal, 3> v2 = v * 2.;
	WorkVec2.AddItem(4, v2);

	GradientAssVec<Gradient<4> > WorkMat2(mh, GradientAssVecBase::APPEND);
	Vector<Gradient<4>, 3> g2 = g * 2.;
	WorkMat2.AddItem(4, g2);

	std::cout << "v=" << std::endl << v << std::endl;
	std::cout << "WorkVec=" << std::endl;
	for (integer i = 1; i <= vh.iGetSize(); ++i) {
		std::cout << " " << vh.iGetRowIndex(i) << " " << vh.dGetCoef(i) << std::endl;
	}

	std::cout << std::endl;

	std::cout << "g=" << std::endl << g << std::endl;

	std::cout << "WorkMat=" << std::endl;

	for (integer i = 1; i <= mh.iGetNumRows(); ++i) {
		std::cout << " " << mh.iGetRowIndex(i)
				  << " " << mh.iGetColIndex(i)
				  << " " << mh.dGetCoef(i, 0) << std::endl;
	}

	std::cout << std::endl;
}

void testInv() {
	Matrix<doublereal, 2, 2> A;

	A(1, 1) =  0.658371336838182;
	A(1, 2) = 0.733036075010795;
	A(2, 1) = 0.483830962404444;
	A(2, 2) = 0.395950513263802;

	const doublereal detA = Det(A);
	const Matrix<doublereal, 2, 2> invA = Inv(A);
	const Matrix<doublereal, 2, 2> B1 = invA * A;
	const Matrix<doublereal, 2, 2> B2 = A * invA;
	std::cout << "A=" << Tabular(A) << std::endl;
	std::cout << "Inv(A)=" << Tabular(invA) << std::endl;
	std::cout << "Det(A)=" << detA << std::endl;
	std::cout << "Inv(A)*A=" << Tabular(B1) << std::endl;
	std::cout << "A * Inv(A)=" << Tabular(B2) << std::endl;

	const doublereal dTol = sqrt(std::numeric_limits<doublereal>::epsilon());

	assert(bCompare(B1(1, 1), 1., dTol));
	assert(bCompare(B2(1, 1), 1., dTol));
	assert(bCompare(B1(2, 2), 1., dTol));
	assert(bCompare(B2(2, 2), 1., dTol));
	assert(bCompare(B1(1, 2), 0., dTol));
	assert(bCompare(B2(1, 2), 0., dTol));
	assert(bCompare(B1(2, 1), 0., dTol));
	assert(bCompare(B2(2, 1), 0., dTol));

	assert(bCompare(invA(1, 1), -4.21299780160758, dTol));
	assert(bCompare(invA(2, 2), -7.00521126207687, dTol));
	assert(bCompare(invA(1, 2), 7.79965998039245, dTol));
	assert(bCompare(invA(2, 1), 5.14806449967026, dTol));
	assert(bCompare(detA, -0.0939830809103952, dTol));
}

void testSolve() {
	const Mat3x3 A(1.32972137393521,      0.61849905020148,     0.709385530146435,
				   0.61849905020148,     0.290559435134808,     0.344471283014357,
				   0.709385530146435,     0.344471283014357,     0.752776020323268);

	const Vec3 b(0.815664130323409,
		     	 0.255816061333836,
		     	 0.416955203644826);

	const doublereal dTol = sqrt(std::numeric_limits<doublereal>::epsilon());

	const Vec3 x1 = A.Solve(b);
	const Vec3 x2 = A.LDLSolve(b);
	const doublereal f1 = (A * x1 - b).Norm();
	const doublereal f2 = (A * x2 - b).Norm();
	std::cout << "A=" << Tabular(Matrix<doublereal, 3, 3>(A)) << std::endl;
	std::cout << "b=" << b << std::endl;
	std::cout << "x1=" << x1 << std::endl;
	std::cout << "x2=" << x2 << std::endl;
	std::cout << "f1=" << f1 << std::endl;
	std::cout << "f2=" << f2 << std::endl;
	assert(f1 < dTol);
	assert(f2 < dTol);
}

doublereal dStartTime;

void tic(doublereal& dTime) {
	dTime = mbdyn_clock_time();
}

void tic() {
    dStartTime = mbdyn_clock_time();
}

doublereal toc() {
	return mbdyn_clock_time() - dStartTime;
}

template <index_type N>
void testMatVec() {
	doublereal c_MatVecGradient[N], cd_MatVecGradient[N][N];
	doublereal c_MatVecDouble[N];

	std::cerr << "---------------------------\ntestMatVecGradient<" << N << ">()\n";
    testMatVecGradient<N>(c_MatVecGradient, cd_MatVecGradient);

    std::cerr << "---------------------------\ntestMatVecDouble<" << N << ">()\n";
    testMatVecDouble<N>(c_MatVecDouble);

#ifdef HAVE_BLITZ
	doublereal c_MatVecDoubleBlitz[N];
	doublereal c_MatVecGradientBlitz[N], cd_MatVecGradientBlitz[N][N];
    std::cerr << "---------------------------\ntestMatVecDoubleBlitz<" << N << ">()\n";
    testMatVecDoubleBlitz<N>(c_MatVecDoubleBlitz);

	std::cerr << "---------------------------\ntestMatVecGradientBlitz<" << N << ">()\n";
    testMatVecGradientBlitz<N>(c_MatVecGradientBlitz, cd_MatVecGradientBlitz);
#endif // HAVE_BLITZ

    const doublereal dTol = 10 * std::numeric_limits<scalar_deriv_type>::epsilon();

    for (int i = 0; i < N; ++i) {
    	assert(bCompare(c_MatVecGradient[i], c_MatVecDouble[i], dTol));
#ifdef HAVE_BLITZ
    	assert(bCompare(c_MatVecGradient[i], c_MatVecDoubleBlitz[i], dTol));
    	assert(bCompare(c_MatVecGradient[i], c_MatVecGradientBlitz[i], dTol));

    	for (int j = 0; j < N; ++j) {
    		assert(bCompare(cd_MatVecGradient[i][j], cd_MatVecGradientBlitz[i][j]));
    	}
#endif
    }
}

template <int iNumDofMax>
void cppad_benchmark1(const int N) {
    const int iNumDof = 6;
    
    assert(iNumDofMax == 0 || iNumDofMax >= iNumDof);
    
    Matrix<Gradient<iNumDofMax>, 3, 3> R;
    Vector<Gradient<iNumDofMax>, 3> X, Y, Phi;
    Vector<doublereal, 3> o;
    LocalDofMap dof;
        
    o(1) = 1.;
    o(2) = 2.;
    o(3) = 3.;
    
    Matrix<doublereal, 3, iNumDof> jac;
    
    const double start = mbdyn_clock_time();
    double calc = 0.;
    
    for (int loop = 0; loop < N; ++loop) {            
        for (int i = 0; i < 3; ++i) {
            X(i + 1).SetValuePreserve(0.);
            X(i + 1).DerivativeResizeReset(&dof, i, MapVectorBase::GLOBAL, 1.);
            Phi(i + 1).SetValuePreserve(0.);
            Phi(i + 1).DerivativeResizeReset(&dof, i + 3, MapVectorBase::GLOBAL, 1.);
        }
        
        const double start_calc = mbdyn_clock_time();
        
        gradVecAssTest::Euler123ToMatR(Phi, R);
        
        Y = X + R * o;
        
        calc += mbdyn_clock_time() - start_calc;
        
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < iNumDof; ++j) {
                jac(i + 1, j + 1) = Y(i + 1).dGetDerivativeGlobal(j);
            }
        }
    }
    
    std::cout << "calculation time: " << calc/N << "s\n";
    std::cout << "elapsed time: " << (mbdyn_clock_time() - start)/N << "s\n";
    
    std::cout << "x=" << X << std::endl;
    std::cout << "y=" << Y << std::endl;
    std::cout << "jac=\n";

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < iNumDof; ++j) {
            std::cout << jac(i + 1, j + 1) << '\t';
        }
        std::cout << std::endl;
    }
}

template <int iNumDofMax>
void cppad_benchmark2(const int N) {
    const int iNumDof = 12;
    
    assert(iNumDofMax == 0 || iNumDofMax >= iNumDof);
    
    Matrix<Gradient<iNumDofMax>, 3, 3> R1, R2;
    Vector<Gradient<iNumDofMax>, 3> X1, X2, Y, Phi1, Phi2;
    Vector<doublereal, 3> o1, o2;
    LocalDofMap dof;
        
    o1(1) = 1.;
    o1(2) = 2.;
    o1(3) = 3.;

    o2(1) = 1.;
    o2(2) = 2.;
    o2(3) = 3.;
    
    Matrix<doublereal, 3, iNumDof> jac;
    
    const double start = mbdyn_clock_time();
    double calc = 0.;
    
    for (int loop = 0; loop < N; ++loop) {            
        for (int i = 0; i < 3; ++i) {
            X1(i + 1).SetValuePreserve(0.);
            X1(i + 1).DerivativeResizeReset(&dof, i, MapVectorBase::GLOBAL, 1.);
            Phi1(i + 1).SetValuePreserve(0.);
            Phi1(i + 1).DerivativeResizeReset(&dof, i + 3, MapVectorBase::GLOBAL, 1.);
            X2(i + 1).SetValuePreserve(0.);
            X2(i + 1).DerivativeResizeReset(&dof, i + 6, MapVectorBase::GLOBAL, 1.);
            Phi2(i + 1).SetValuePreserve(0.);
            Phi2(i + 1).DerivativeResizeReset(&dof, i + 9, MapVectorBase::GLOBAL, 1.);            
        }
        
        const double start_calc = mbdyn_clock_time();
        
        gradVecAssTest::Euler123ToMatR(Phi1, R1);
        gradVecAssTest::Euler123ToMatR(Phi2, R2);
        
        Y = X1 + R1 * o1 - X2 - R2 * o2;
        
        calc += mbdyn_clock_time() - start_calc;
        
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < iNumDof; ++j) {
                jac(i + 1, j + 1) = Y(i + 1).dGetDerivativeGlobal(j);
            }
        }
    }
    
    std::cout << "calculation time: " << calc/N << "s\n";
    std::cout << "elapsed time: " << (mbdyn_clock_time() - start)/N << "s\n";
    
    std::cout << "X1=" << X1 << std::endl;
    std::cout << "X2=" << X2 << std::endl;
    std::cout << "Y=" << Y << std::endl;
    std::cout << "jac=\n";

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < iNumDof; ++j) {
            std::cout << jac(i + 1, j + 1) << '\t';
        }
        std::cout << std::endl;
    }
}

template <int iNumDofMax>
void cppad_benchmark3(const int N) {
    const int iNumDof = 12;
    
    assert(iNumDofMax == 0 || iNumDofMax >= iNumDof);
    
    typedef Matrix<Gradient<iNumDofMax>, 3, 3> Mat3x3;
    typedef Vector<Gradient<iNumDofMax>, 3> Vec3;
    typedef Vector<doublereal, 3> CVec3;
    
    Mat3x3 R1, R2;
    Vec3 X1, X2, Y, Phi1, Phi2;
    CVec3 o1, o2;
    LocalDofMap dof;
        
    o1(1) = 1.;
    o1(2) = 2.;
    o1(3) = 3.;

    o2(1) = 1.;
    o2(2) = 2.;
    o2(3) = 3.;
    
    Matrix<doublereal, 3, iNumDof> jac;
    
    const double start = mbdyn_clock_time();
    double calc = 0.;
    
    for (int loop = 0; loop < N; ++loop) {            
        for (int i = 0; i < 3; ++i) {
            X1(i + 1).SetValuePreserve(0.);
            X1(i + 1).DerivativeResizeReset(&dof, i, MapVectorBase::GLOBAL, 1.);
            Phi1(i + 1).SetValuePreserve(0.);
            Phi1(i + 1).DerivativeResizeReset(&dof, i + 3, MapVectorBase::GLOBAL, 1.);
            X2(i + 1).SetValuePreserve(0.);
            X2(i + 1).DerivativeResizeReset(&dof, i + 6, MapVectorBase::GLOBAL, 1.);
            Phi2(i + 1).SetValuePreserve(0.);
            Phi2(i + 1).DerivativeResizeReset(&dof, i + 9, MapVectorBase::GLOBAL, 1.);            
        }
        
        const double start_calc = mbdyn_clock_time();
        
        gradVecAssTest::Euler123ToMatR(Phi1, R1);
        gradVecAssTest::Euler123ToMatR(Phi2, R2);
        
        Y = Transpose(R2) * Vec3(X1 + R1 * o1 - X2) - o2;
        
        calc += mbdyn_clock_time() - start_calc;
        
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < iNumDof; ++j) {
                jac(i + 1, j + 1) = Y(i + 1).dGetDerivativeGlobal(j);
            }
        }
    }
    
    std::cout << "calculation time: " << calc/N << "s\n";
    std::cout << "elapsed time: " << (mbdyn_clock_time() - start)/N << "s\n";
    
    std::cout << "X1=" << X1 << std::endl;
    std::cout << "X2=" << X2 << std::endl;
    std::cout << "Y=" << Y << std::endl;
    std::cout << "jac=\n";

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < iNumDof; ++j) {
            std::cout << jac(i + 1, j + 1) << '\t';
        }
        std::cout << std::endl;
    }
}

void Mat3xN_test(int N, int M) {
    Mat3xN A(N, 0.);

    for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
        for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
            A(i, j) = 100 * i + j;
        }
    }

    Vector<doublereal, DYNAMIC_SIZE> x1(A.iGetNumCols()), x2(A.iGetNumCols()), x(A.iGetNumCols());

    for (index_type i = 1; i <= x.iGetNumRows(); ++i) {
        x1(i) = i;
        x2(i) = 10 * i;
    }
    
    Vector<doublereal, 3> b = A * x;
    
    const double start = mbdyn_clock_time();
    
    for (int i = 0; i < M; ++i) {
        x = x1 * 3. + x2 * 5.;
        b = A * x;
    }

    std::cerr << "Mat3xN: " << (mbdyn_clock_time() - start) / M << "s\n";

    const doublereal tol = sqrt(std::numeric_limits<doublereal>::epsilon());
    
    for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
        doublereal b_i = 0.;
        
        for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
            b_i += A(i, j) * (3. * x1(j) + 5. * x2(j));
        }

        assert(bCompare(b_i, b(i), tol));
    }
}

void MatNxN_test(int N, int M) {
    MatNxN A(N, 0.);

    for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
        for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
            A(i, j) = 100 * i + j;
        }
    }

    Vector<doublereal, DYNAMIC_SIZE> x1(A.iGetNumCols()), x2(A.iGetNumCols()), x(A.iGetNumCols());

    for (index_type i = 1; i <= x.iGetNumRows(); ++i) {
        x1(i) = i;
        x2(i) = 10 * i;
    }
    
    Vector<doublereal, DYNAMIC_SIZE> b = A * x;
    
    const double start = mbdyn_clock_time();
    
    for (int i = 0; i < M; ++i) {
        x = x1 * 3. + x2 * 5.;
        b = A * x;
    }

    std::cerr << "MatNxN: " << (mbdyn_clock_time() - start) / M << "s\n";

    const doublereal tol = sqrt(std::numeric_limits<doublereal>::epsilon());
    
    for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
        doublereal b_i = 0.;
        
        for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
            b_i += A(i, j) * (3. * x1(j) + 5. * x2(j));
        }

        assert(bCompare(b_i, b(i), tol));
    }
}

void MatDynamic_test(index_type iNumRows, index_type iNumCols, index_type iNumLoops) {
    Matrix<doublereal, DYNAMIC_SIZE, DYNAMIC_SIZE> A(iNumRows, iNumCols);

    for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
        for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
            A(i, j) = 100 * i + j;
        }
    }

    Vector<doublereal, DYNAMIC_SIZE> x1(A.iGetNumCols()), x2(A.iGetNumCols()), x(A.iGetNumCols());

    for (index_type i = 1; i <= x.iGetNumRows(); ++i) {
        x1(i) = i;
        x2(i) = 10 * i;
    }
    
    Vector<doublereal, DYNAMIC_SIZE> b = A * x;
    
    const double start = mbdyn_clock_time();
    
    for (int i = 0; i < iNumLoops; ++i) {
        x = x1 * 3. + x2 * 5.;
        b = A * x;
    }

    std::cerr << "Matrix<DYNAMIC_SIZE, DYNAMIC_SIZE>: " << (mbdyn_clock_time() - start) / iNumLoops << "s\n";

    const doublereal tol = sqrt(std::numeric_limits<doublereal>::epsilon());
    
    for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
        doublereal b_i = 0.;
        
        for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
            b_i += A(i, j) * (3. * x1(j) + 5. * x2(j));
        }

        assert(bCompare(b_i, b(i), tol));
    }
}


template <index_type N_SIZE>
void Mat3xN_test_grad(int iNumDeriv, int N, int M) {
    assert((N_SIZE == 0) || (N_SIZE >= iNumDeriv));
    
    LocalDofMap dof;
    
    Mat3xN A(N, 0.);

    for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
        for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
            A(i, j) = 100 * i + j;
        }
    }

    Vector<Gradient<N_SIZE>, DYNAMIC_SIZE> x1(A.iGetNumCols()), x2(A.iGetNumCols()), x(A.iGetNumCols());

    for (index_type i = 1; i <= x.iGetNumRows(); ++i) {
        x1(i).SetValuePreserve(i);
        x1(i).DerivativeResizeReset(&dof, 0, iNumDeriv, MapVectorBase::LOCAL, 0.);
        
        for (index_type k = 0; k < iNumDeriv; ++k) {
            x1(i).SetDerivativeLocal(k, -1. * k);
        }
        
        x2(i).SetValuePreserve(10 * i);
        x2(i).DerivativeResizeReset(&dof, 0, iNumDeriv, MapVectorBase::LOCAL, 0.);
        
        for (index_type k = 0; k < iNumDeriv; ++k) {
            x2(i).SetDerivativeLocal(k, 10. * k);
        }
    }
    
    Vector<Gradient<N_SIZE>, 3> b = A * x;
    
    const double start = mbdyn_clock_time();
    
    for (int i = 0; i < M; ++i) {
        x = x1 * 3. + x2 * 5.;
        b = A * x;
    }

    std::cerr << "Mat3xN * Gradient: " << (mbdyn_clock_time() - start) / M << "s\n";

    const doublereal tol = sqrt(std::numeric_limits<doublereal>::epsilon());
    
    for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
        Gradient<N_SIZE> b_i;
        
        for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
            b_i += A(i, j) * (3. * x1(j) + 5. * x2(j));
        }

        assert(bCompare(b_i, b(i), tol));
    }
}

template <index_type N_SIZE>
void Mat3xNT_test_grad(int iNumDeriv, int N, int M) {
    assert((N_SIZE == 0) || (N_SIZE >= iNumDeriv));
    
    LocalDofMap dof;
    
    Mat3xN A(N, 0.);

    for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
        for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
            A(i, j) = 100 * i + j;
        }
    }

    Vector<Gradient<N_SIZE>, 3> x1(A.iGetNumRows()), x2(A.iGetNumRows());

    for (index_type i = 1; i <= x1.iGetNumRows(); ++i) {
        x1(i).SetValuePreserve(i);
        x1(i).DerivativeResizeReset(&dof, 0, iNumDeriv, MapVectorBase::LOCAL, 0.);
        
        for (index_type k = 0; k < iNumDeriv; ++k) {
            x1(i).SetDerivativeLocal(k, -1. * k);
        }
        
        x2(i).SetValuePreserve(10 * i);
        x2(i).DerivativeResizeReset(&dof, 0, iNumDeriv, MapVectorBase::LOCAL, 0.);
        
        for (index_type k = 0; k < iNumDeriv; ++k) {
            x2(i).SetDerivativeLocal(k, 10. * k);
        }
    }
    
    Vector<Gradient<N_SIZE>, DYNAMIC_SIZE> b = Transpose(A) * x1;
    
    const double start = mbdyn_clock_time();
    
    for (int i = 0; i < M; ++i) {
        b = Transpose(A) * (x1 * 3. + x2 * 5.);
    }

    std::cerr << "Transpose(Mat3xN) * Gradient: " << (mbdyn_clock_time() - start) / M << "s\n";

    const doublereal tol = sqrt(std::numeric_limits<doublereal>::epsilon());
    
    for (index_type i = 1; i <= A.iGetNumCols(); ++i) {
        Gradient<N_SIZE> b_i;
        
        for (index_type j = 1; j <= A.iGetNumRows(); ++j) {
            b_i += A(j, i) * (3. * x1(j) + 5. * x2(j));
        }

        assert(bCompare(b_i, b(i), tol));
    }
}

template <index_type N_SIZE>
void MatNxNT_test_grad(int iNumDeriv, int N, int M) {
    assert((N_SIZE == 0) || (N_SIZE >= iNumDeriv));
    
    LocalDofMap dof;
    
    MatNxN A(N, 0.);

    for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
        for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
            A(i, j) = 100 * i + j;
        }
    }

    Vector<Gradient<N_SIZE>, DYNAMIC_SIZE> x1(A.iGetNumRows()), x2(A.iGetNumRows());

    for (index_type i = 1; i <= x1.iGetNumRows(); ++i) {
        x1(i).SetValuePreserve(i);
        x1(i).DerivativeResizeReset(&dof, 0, iNumDeriv, MapVectorBase::LOCAL, 0.);
        
        for (index_type k = 0; k < iNumDeriv; ++k) {
            x1(i).SetDerivativeLocal(k, -1. * k);
        }
        
        x2(i).SetValuePreserve(10 * i);
        x2(i).DerivativeResizeReset(&dof, 0, iNumDeriv, MapVectorBase::LOCAL, 0.);
        
        for (index_type k = 0; k < iNumDeriv; ++k) {
            x2(i).SetDerivativeLocal(k, 10. * k);
        }
    }
    
    Vector<Gradient<N_SIZE>, DYNAMIC_SIZE> b = Transpose(A) * x1;
    
    const double start = mbdyn_clock_time();
    
    for (int i = 0; i < M; ++i) {
        b = Transpose(A) * (x1 * 3. + x2 * 5.);
    }

    std::cerr << "Transpose(MatNxN) * Gradient: " << (mbdyn_clock_time() - start) / M << "s\n";

    const doublereal tol = sqrt(std::numeric_limits<doublereal>::epsilon());
    
    for (index_type i = 1; i <= A.iGetNumCols(); ++i) {
        Gradient<N_SIZE> b_i;
        
        for (index_type j = 1; j <= A.iGetNumRows(); ++j) {
            b_i += A(j, i) * (3. * x1(j) + 5. * x2(j));
        }

        assert(bCompare(b_i, b(i), tol));
    }
}

template <index_type N_SIZE>
void MatNxN_test_grad(int iNumDeriv, int N, int M)
{
    assert((N_SIZE == 0) || (N_SIZE >= iNumDeriv));
    
    LocalDofMap dof;
    
    MatNxN A(N, 0.);

    for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
        for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
            A(i, j) = 100 * i + j;
        }
    }

    Vector<Gradient<N_SIZE>, DYNAMIC_SIZE> x1(A.iGetNumCols()), x2(A.iGetNumCols()), x(A.iGetNumCols());

    for (index_type i = 1; i <= x.iGetNumRows(); ++i) {
        x1(i).SetValuePreserve(i);
        x1(i).DerivativeResizeReset(&dof, 0, iNumDeriv, MapVectorBase::LOCAL, 0.);
        
        for (index_type k = 0; k < iNumDeriv; ++k) {
            x1(i).SetDerivativeLocal(k, -1. * k);
        }
        
        x2(i).SetValuePreserve(10 * i);
        x2(i).DerivativeResizeReset(&dof, 0, iNumDeriv, MapVectorBase::LOCAL, 0.);
        
        for (index_type k = 0; k < iNumDeriv; ++k) {
            x2(i).SetDerivativeLocal(k, 10. * k);
        }
    }
    
    Vector<Gradient<N_SIZE>, DYNAMIC_SIZE> b = A * x;
    
    const double start = mbdyn_clock_time();
    
    for (int i = 0; i < M; ++i) {
        x = x1 * 3. + x2 * 5.;
        b = A * x;
    }

    std::cerr << "MatNxN * Gradient: " << (mbdyn_clock_time() - start) / M << "s\n";

    const doublereal tol = sqrt(std::numeric_limits<doublereal>::epsilon());
    
    for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
        Gradient<N_SIZE> b_i;
        
        for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
            b_i += A(i, j) * (3. * x1(j) + 5. * x2(j));
        }

        assert(bCompare(b_i, b(i), tol));
    }
}

template <index_type N_SIZE>
void MatDynamic_test_grad(index_type iNumDeriv, index_type iNumRows, index_type iNumCols, int iNumLoops)
{
    assert((N_SIZE == 0) || (N_SIZE >= iNumDeriv));
    
    LocalDofMap dof;
    
    Matrix<Gradient<N_SIZE>, DYNAMIC_SIZE, DYNAMIC_SIZE> A(iNumRows, iNumCols);

    for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
        for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
            A(i, j).SetValuePreserve(100 * i + j);
            A(i, j).DerivativeResizeReset(&dof, 0, iNumDeriv, MapVectorBase::LOCAL, 0.);
            
            for (index_type k = 0; k < iNumDeriv; ++k) {
                A(i, j).SetDerivativeLocal(k, 1000. * i + 100. * j + k + 1);
            }
        }
    }

    Vector<Gradient<N_SIZE>, DYNAMIC_SIZE> x1(A.iGetNumCols()), x2(A.iGetNumCols()), x(A.iGetNumCols());

    for (index_type i = 1; i <= x.iGetNumRows(); ++i) {
        x1(i).SetValuePreserve(i);
        x1(i).DerivativeResizeReset(&dof, 0, iNumDeriv, MapVectorBase::LOCAL, 0.);
        
        for (index_type k = 0; k < iNumDeriv; ++k) {
            x1(i).SetDerivativeLocal(k, -1. * k);
        }
        
        x2(i).SetValuePreserve(10 * i);
        x2(i).DerivativeResizeReset(&dof, 0, iNumDeriv, MapVectorBase::LOCAL, 0.);
        
        for (index_type k = 0; k < iNumDeriv; ++k) {
            x2(i).SetDerivativeLocal(k, 10. * k);
        }
    }
    
    Vector<Gradient<N_SIZE>, DYNAMIC_SIZE> b = A * x;
    
    const double start = mbdyn_clock_time();
    
    for (int i = 0; i < iNumLoops; ++i) {
        x = x1 * 3. + x2 * 5.;
        b = A * x;
    }

    std::cerr << "Matrix<Gradient,DYNAMIC_SIZE,DYNAMIC_SIZE> * Gradient: " << (mbdyn_clock_time() - start) / iNumLoops << "s\n";

    const doublereal tol = sqrt(std::numeric_limits<doublereal>::epsilon());
    
    for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
        Gradient<N_SIZE> b_i;
        
        for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
            b_i += A(i, j) * (3. * x1(j) + 5. * x2(j));
        }

        assert(bCompare(b_i, b(i), tol));
    }
}

template <index_type N_ROWS, index_type N_COLS>
void MatDynamicT_test(index_type iNumRows, index_type iNumCols, int iNumLoops)
{
    assert((N_ROWS == iNumRows) || (N_ROWS == DYNAMIC_SIZE));
    assert((N_COLS == iNumCols) || (N_COLS == DYNAMIC_SIZE));
     
    Matrix<doublereal, N_ROWS, N_COLS> A(iNumRows, iNumCols);

    for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
        for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
            A(i, j) = 100 * i + j;
        }
    }

    Vector<doublereal, N_ROWS> x1(A.iGetNumRows()), x2(A.iGetNumRows());

    for (index_type i = 1; i <= x1.iGetNumRows(); ++i) {
        x1(i) = i;
        x2(i) = 10 * i;
    }
    
    Vector<doublereal, N_COLS> b = Transpose(A) * x1;
    
    const double start = mbdyn_clock_time();
    
    for (int i = 0; i < iNumLoops; ++i) {
        b = Transpose(A) * (x1 * 3. + x2 * 5.);
    }

    std::cerr << "Transpose(Matrix<doublereal>,"
              << iNumRows << "(" << N_ROWS << ")," << iNumCols << "(" << N_COLS << ")"
              << ">) * Vector<doublereal>,"
              << iNumRows << "(" << N_ROWS << ")>: " << (mbdyn_clock_time() - start) / iNumLoops << "s\n";

    const doublereal tol = sqrt(std::numeric_limits<doublereal>::epsilon());
    
    for (index_type i = 1; i <= A.iGetNumCols(); ++i) {
        doublereal b_i = 0.;
        
        for (index_type j = 1; j <= A.iGetNumRows(); ++j) {
            b_i += A(j, i) * (3. * x1(j) + 5. * x2(j));
        }

        assert(bCompare(b_i, b(i), tol));
    }
}

template <index_type N_DERIV, index_type N_ROWS, index_type N_COLS>
void MatDynamicT_test_grad(index_type iNumDeriv, index_type iNumRows, index_type iNumCols, int iNumLoops)
{
    assert((N_DERIV == 0) || (N_DERIV >= iNumDeriv));
    assert((N_ROWS == iNumRows) || (N_ROWS == DYNAMIC_SIZE));
    assert((N_COLS == iNumCols) || (N_COLS == DYNAMIC_SIZE));
    
    LocalDofMap dof;
    
    Matrix<Gradient<N_DERIV>, N_ROWS, N_COLS> A(iNumRows, iNumCols);

    for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
        for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
            A(i, j).SetValuePreserve(100 * i + j);
            A(i, j).DerivativeResizeReset(&dof, 0, iNumDeriv, MapVectorBase::LOCAL, 0.);
            
            for (index_type k = 0; k < iNumDeriv; ++k) {
                A(i, j).SetDerivativeLocal(k, 1000. * i + 100. * j + k + 1);
            }
        }
    }

    Vector<Gradient<N_DERIV>, N_ROWS> x1(A.iGetNumRows()), x2(A.iGetNumRows());

    for (index_type i = 1; i <= x1.iGetNumRows(); ++i) {
        x1(i).SetValuePreserve(i);
        x1(i).DerivativeResizeReset(&dof, 0, iNumDeriv, MapVectorBase::LOCAL, 0.);
        
        for (index_type k = 0; k < iNumDeriv; ++k) {
            x1(i).SetDerivativeLocal(k, -1. * k);
        }
        
        x2(i).SetValuePreserve(10 * i);
        x2(i).DerivativeResizeReset(&dof, 0, iNumDeriv, MapVectorBase::LOCAL, 0.);
        
        for (index_type k = 0; k < iNumDeriv; ++k) {
            x2(i).SetDerivativeLocal(k, 10. * k);
        }
    }
    
    Vector<Gradient<N_DERIV>, N_COLS> b = Transpose(A) * x1;
    
    const double start = mbdyn_clock_time();
    
    for (int i = 0; i < iNumLoops; ++i) {
        b = Transpose(A) * (x1 * 3. + x2 * 5.);
    }

    std::cerr << "Transpose(Matrix<Gradient<" << iNumDeriv << "(" << N_DERIV << ")"
              << ">," << iNumRows << "(" << N_ROWS << ")," << iNumCols << "(" << N_COLS << ")"
              << ">) * Vector<Gradient<" << iNumDeriv << "(" << N_DERIV << ")"
              << ">," << iNumRows << "(" << N_ROWS << ")>: " << (mbdyn_clock_time() - start) / iNumLoops << "s\n";

    const doublereal tol = sqrt(std::numeric_limits<doublereal>::epsilon());
    
    for (index_type i = 1; i <= A.iGetNumCols(); ++i) {
        Gradient<N_DERIV> b_i;
        
        for (index_type j = 1; j <= A.iGetNumRows(); ++j) {
            b_i += A(j, i) * (3. * x1(j) + 5. * x2(j));
        }

        assert(bCompare(b_i, b(i), tol));
    }

    for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
        for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
            assert(bCompare(A(i, j), A.GetRow(i)(j), 0.));
            assert(bCompare(A(i, j), A.GetCol(j)(i), 0.));
            assert(bCompare(A(i, j), Transpose(A)(j, i), 0.));
            assert(bCompare(A(i, j), Transpose(A).GetCol(i)(j), 0.));
            assert(bCompare(A(i, j), Transpose(A).GetRow(j)(i), 0.));
        }
    }
}

void MatManip_test(int NLoops) {
    const doublereal tol1 = 10. * std::numeric_limits<doublereal>::epsilon();
    const doublereal tol2 = 1e-4 * sqrt(tol1);
    const doublereal alpha = 1. - sqrt(tol2);
    
    srand(0);
    
    for (int k = 0; k < 3 * NLoops; ++k) {
        Vector<doublereal, 3> g0;
        Vec3 g0_;
        
        for (int i = 1; i <= 3; ++i) {
            g0_(i) = (2. * rand() / RAND_MAX - 1.);
        }

        g0_ *= (2. * rand() / RAND_MAX - 1.) * alpha * M_PI / g0_.Norm();

        if (k >= NLoops && k < 2 * NLoops) {
            g0_ *= sqrt(std::numeric_limits<doublereal>::epsilon());
        } else if (k >= 2 * NLoops) {
            g0_ *= std::numeric_limits<doublereal>::epsilon();
        }
        
        g0 = g0_;
        
        const Mat3x3 G1(CGR_Rot::MatG, g0_);
        const Mat3x3 R1(CGR_Rot::MatR, g0_);
        const Mat3x3 X1(1., g0_);
        const Matrix<doublereal, 3, 3> G2(MatGVec(g0));
        Matrix<doublereal, 3, 3> R2(MatRVec(g0));
        const Matrix<doublereal, 3, 3> X2(MatCrossVec(g0, 1.));
        
        for (index_type i = 1; i <= 3; ++i) {
            for (index_type j = 1; j <= 3; ++j) {
                assert(bCompare(G1(i, j), G2(i, j), tol1));
                assert(bCompare(R1(i, j), R2(i, j), tol1));
                assert(bCompare(X1(i, j), X2(i, j), tol1));
            }
        }

        if (k == 0) {
            R2 = ::Eye3;
        } else if (k == 1) {
            R2(1, 1) = -1;
            R2(1, 2) = 0;
            R2(1, 3) = 0;
            R2(2, 1) = 0;
            R2(2, 2) = -1;
            R2(2, 3) = 0;
            R2(3, 1) = 0;
            R2(3, 2) = 0;
            R2(3, 3) = 1;
        } else if (k == 2) {
            R2(1, 1) = -1;
            R2(1, 2) = 0;
            R2(1, 3) = 0;
            R2(2, 1) = 0;
            R2(2, 2) = 1;
            R2(2, 3) = 0;
            R2(3, 1) = 0;
            R2(3, 2) = 0;
            R2(3, 3) = -1;            
        } else if (k == 3) {
            R2(1, 1) = 1;
            R2(1, 2) = 0;
            R2(1, 3) = 0;
            R2(2, 1) = 0;
            R2(2, 2) = -1;
            R2(2, 3) = 0;
            R2(3, 1) = 0;
            R2(3, 2) = 0;
            R2(3, 3) = -1;
        } else if (k == 4) {            
            R2(1,1)=-2.2841125213377644e-01;
            R2(1,2)=9.5997603033895429e-01;
            R2(1,3)=1.6209355654480440e-01;
            R2(2,1)=9.5997603033895418e-01;
            R2(2,2)=1.9435901751267337e-01;
            R2(2,3)=2.0166951551033063e-01;
            R2(3,1)=1.6209355654480423e-01;
            R2(3,2)=2.0166951551033083e-01;
            R2(3,3)=-9.6594776537889704e-01;         
        }
        
        Mat3x3 R2_;

        for (index_type i = 1; i <= 3; ++i) {
            for (index_type j = 1; j <= 3; ++j) {
                R2_(i, j) = R2(i, j);
            }
        }
        
        const Vec3 g1(RotManip::VecRot(R2_));
        const Vector<doublereal, 3> g2(VecRotMat(R2));
        const Matrix<doublereal, 3, 3> R2__(MatRotVec(g2));
        for (index_type i = 1; i <= 3; ++i) {
            assert(bCompare(g1(i),
                            g2(i),
                            std::numeric_limits<doublereal>::epsilon())
                   || bCompare(fabs(g1(i) - g2(i)),
                               2 * M_PI,
                               2 * M_PI * std::numeric_limits<doublereal>::epsilon()));
        }

        for (index_type i = 1; i <= 3; ++i) {
                for (index_type j = 1; j <= 3; ++j) {
                        assert(bCompare(R2__(i, j), R2(i, j), std::pow(std::numeric_limits<doublereal>::epsilon(), 0.9)));
                }
        }     
    }
}

doublereal MatIndexValGen(index_type i, index_type j, index_type k) {
    return 100. * i + 10. * j + k;
}

doublereal MatRandValGen(index_type, index_type, index_type) {
    return 2. * (static_cast<doublereal>(rand()) / RAND_MAX) - 1.;
}

template<index_type N_rows_lhs, index_type N_cols_lhs, index_type N_cols_rhs, doublereal (*f)(index_type, index_type, index_type)>
void MatrixMatrixProduct_test(index_type iNumRowsLhs, index_type iNumColsLhs, index_type iNumColsRhs)
{
    assert(N_rows_lhs == DYNAMIC_SIZE || N_rows_lhs == iNumRowsLhs);
    assert(N_cols_lhs == DYNAMIC_SIZE || N_cols_lhs == iNumColsLhs);
    assert(N_cols_rhs == DYNAMIC_SIZE || N_cols_rhs == iNumColsRhs);
    
    Matrix<doublereal, N_rows_lhs, N_cols_lhs> A(iNumRowsLhs, iNumColsLhs);
    Matrix<doublereal, N_cols_lhs, N_cols_rhs> B(iNumColsLhs, iNumColsRhs);
    Matrix<doublereal, N_rows_lhs, N_cols_rhs> D(iNumRowsLhs, iNumColsRhs);
    
    for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
        for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
            A(i, j) = f(i, j, 1);
        }
    }

    for (index_type i = 1; i <= B.iGetNumRows(); ++i) {
        for (index_type j = 1; j <= B.iGetNumCols(); ++j) {
            B(i, j) = f(i, j, 2);
        }
    }

    for (index_type i = 1; i <= D.iGetNumRows(); ++i) {
        for (index_type j = 1; j <= D.iGetNumCols(); ++j) {
            D(i, j) = f(i, j, 2);
        }
    }
    
    Matrix<doublereal, N_rows_lhs, N_cols_rhs> C = A * B + D;

    Matrix<doublereal, N_rows_lhs, N_cols_rhs> C2(iNumRowsLhs, iNumColsRhs);

    for (index_type i = 1; i <= C2.iGetNumRows(); ++i) {
        for (index_type j = 1; j <= C2.iGetNumCols(); ++j) {
            C2(i, j) = D(i, j);
            
            for (index_type k = 1; k <= A.iGetNumCols(); ++k) {
                C2(i, j) += A(i, k) * B(k, j);
            }
        }
    }

    assert(C2.iGetNumRows() == C.iGetNumRows());
    assert(C2.iGetNumCols() == C.iGetNumCols());

    bool bOk = true;
    
    for (index_type i = 1; i <= C.iGetNumRows(); ++i) {
        for (index_type j = 1; j <= C.iGetNumCols(); ++j) {
            bOk = bOk && bCompare(C2(i, j), C(i, j), sqrt(std::numeric_limits<doublereal>::epsilon()));
        }
    }
    
    std::cerr << "A=\n" << Tabular(A) << std::endl;
    std::cerr << "B=\n" << Tabular(B) << std::endl;
    std::cerr << "C=\n" << Tabular(C) << std::endl;
    std::cerr << "D=\n" << Tabular(D) << std::endl;    
    std::cerr << "C2=\n" << Tabular(C2) << std::endl;

    if (!bOk) {
        std::cerr << "test failed: " << __FILE__ << ":" << __LINE__ << std::endl;        
        assert(bOk);
    }
}
#endif

#ifdef TEST_GROUP2
void ComplianceModelTest(const index_type iNumNodes, const index_type iNumModes)
{
    std::cerr << "ComplianceModelTest(" << iNumNodes << ")\n";
    
    LocalDofMap oDofPress, oDofPressAsp, oDofDef;

    Vector<Gradient<0>, DYNAMIC_SIZE> p(iNumNodes), paspn(iNumNodes), pasp(iNumNodes), w(iNumNodes), a(iNumModes);

    const index_type iFirstDofX = 3;
    const index_type iFirstDofPress = iFirstDofX + 12 + 5;
    const index_type iFirstDofDef = iFirstDofPress + iNumNodes + 3;
    const index_type iFirstDofMod = iFirstDofDef + iNumNodes + 3;
    const index_type iNumDofTotal = iFirstDofMod + iNumModes - 1;
    
    Vector<Gradient<0>, 12> X;
    std::vector<bool> bPasp(iNumNodes);

    for (index_type i = 0; i < iNumNodes; ++i) {
        bPasp[i] = random1() >= 0;
    }

    std::vector<bool> bCav(iNumNodes);

    for (index_type i = 0; i < iNumNodes; ++i) {
        bCav[i] = random1() >= -0.2;
    }
    
    doublereal pEnv = 0.1;
    doublereal dScale = 10.;
    doublereal dScaleP = 5;

    std::vector<index_type> perm(iNumNodes), iperm(iNumNodes);

    for (index_type i = 1; i <= iNumNodes; ++i) {
        perm[i - 1] = i;
    }

    std::random_shuffle(perm.begin(), perm.end());

    for (index_type i = 1; i <= iNumNodes; ++i) {
        iperm[perm[i - 1] - 1] = i;
    }

    for (index_type i = 1; i <= iNumNodes; ++i) {
        assert(iperm[perm[i - 1] - 1] == i);
    }

    for (index_type i = 1; i <= iNumNodes; ++i) {
        std::cerr << i << "->" << perm[i - 1] << std::endl;
    }
    
    for (index_type i = 1; i <= X.iGetNumRows(); ++i) {
        X(i).SetValuePreserve(i);
        X(i).DerivativeResizeReset(&oDofDef, iFirstDofX + i - 1, MapVectorBase::GLOBAL, -1.);
    }
    
    for (index_type i = 1; i <= iNumNodes; ++i) {
        p(i).SetValuePreserve(perm[i - 1] * 100);
        p(i).DerivativeResizeReset(&oDofPress, iFirstDofPress + perm[i - 1] - 1, MapVectorBase::GLOBAL, -1.);

        if (bCav[i - 1]) {
            p(i) = 0.;
        }
        
        p(i) -= pEnv;
        p(i) *= dScaleP;
    }

    for (index_type i = 1; i <= iNumNodes; ++i) {
        w(i).SetValuePreserve(perm[i - 1] * 10);
        w(i).DerivativeResizeReset(&oDofDef, iFirstDofDef + perm[i - 1] - 1, MapVectorBase::GLOBAL, -1.);
    }

    for (index_type i = 1; i <= iNumModes; ++i) {
        a(i).SetValuePreserve(i * 100);
        a(i).DerivativeResizeReset(&oDofPress, iFirstDofMod + i - 1, MapVectorBase::GLOBAL, -1.);
    }
    
    for (index_type i = 1; i <= iNumNodes; ++i) {
        if (bPasp[i - 1]) {
            for (index_type j = 1; j <= 12; ++j) {
                paspn(i) += 1e-6 * pow(w(i) + X(j), 2.0);
            }
        }
    }
    
    for (index_type i = 1; i <= iNumNodes; ++i) {
        pasp(i).Copy(paspn(i), &oDofPressAsp);
        pasp(i) -= pEnv;
        pasp(i) *= dScaleP;
    }

    Matrix<doublereal, DYNAMIC_SIZE, DYNAMIC_SIZE> C(iNumNodes, iNumNodes), D(iNumNodes, iNumModes);

    for (index_type i = 1; i <= iNumNodes; ++i) {
        for (index_type j = 1; j <= iNumNodes; ++j) {
            C(i, j) = random1();
        }

        for (index_type j = 1; j <= iNumModes; ++j) {
            D(i, j) = random1();
        }
    }

    SparseSubMatrixHandler mh(4 * iNumNodes * (iNumNodes + iNumModes));
    GradientAssVec<Gradient<0> > WorkVec(mh);

    Gradient<0> f;

    index_type iEqIndex = 0;
    
    for (index_type i = 1; i <= iNumNodes; ++i) {
        f.SetValuePreserve(perm[iEqIndex] * 10);
        f.DerivativeResizeReset(&oDofPress, iFirstDofDef + perm[iEqIndex] - 1, MapVectorBase::GLOBAL, -1.);

        f *= dScale;

        for (index_type j = 1; j <= iNumNodes; ++j) {
            const doublereal Cij = dScale * C(i, j);
            f -= Cij * p(j);
            f -= Cij * pasp(j);
        }

        for (index_type j = 1; j <= iNumModes; ++j) {
            f -= dScale * D(i, j) * a(j);
        }
        
        WorkVec.AddItem(perm[iEqIndex++], f);
    }

    SpMapMatrixHandler A(iNumNodes, iNumDofTotal);

    A += mh;

    Matrix<doublereal, DYNAMIC_SIZE, 12> dpasp_dX(iNumNodes, 12);

    for (index_type i = 1; i <= iNumNodes; ++i) {
        if (bPasp[i - 1]) {
            for (index_type j = 1; j <= 12; ++j) {
                dpasp_dX(i, j) = 2e-6 * (w(i).dGetValue() + X(j).dGetValue());
            }
        }
    }

    Matrix<doublereal, DYNAMIC_SIZE, DYNAMIC_SIZE> dpasp_dw(iNumNodes, iNumNodes);

    for (index_type i = 1; i <= iNumNodes; ++i) {
        if (bPasp[i - 1]) {
            for (index_type j = 1; j <= 12; ++j) {
                dpasp_dw(i, i) += 2e-6 * (w(i).dGetValue() + X(j).dGetValue());
            }
        }
    }

    Matrix<doublereal, DYNAMIC_SIZE, 12> df_dX = C * dpasp_dX * (-dScale * dScaleP);
    Matrix<doublereal, DYNAMIC_SIZE, DYNAMIC_SIZE> df_dp = C * (-dScale * dScaleP);
    Matrix<doublereal, DYNAMIC_SIZE, DYNAMIC_SIZE> df_dw = C * dpasp_dw * (-dScale * dScaleP);
    Matrix<doublereal, DYNAMIC_SIZE, DYNAMIC_SIZE> df_da = D * (-dScale);
    
    for (index_type i = 1; i <= iNumNodes; ++i) {
        df_dw(i, i) += dScale;
    }

    for (index_type j = 1; j <= iNumNodes; ++j) {
        if (bCav[j - 1]) {
            for (index_type i = 1; i <= iNumNodes; ++i) {
                df_dp(i, j) = 0;
            }
        }
    }

    Matrix<doublereal, DYNAMIC_SIZE, DYNAMIC_SIZE> Aref(iNumNodes, iNumDofTotal);

    for (index_type i = 1; i <= iNumNodes; ++i) {
        for (index_type j = 1; j <= 12; ++j) {
            Aref(perm[i - 1], j - 1 + iFirstDofX) = -df_dX(i, j);
        }

        for (index_type j = 1; j <= iNumNodes; ++j) {
            Aref(perm[i - 1], perm[j - 1] - 1 + iFirstDofPress) = -df_dp(i, j);
        }

        for (index_type j = 1; j <= iNumNodes; ++j) {
            Aref(perm[i - 1], perm[j - 1] - 1 + iFirstDofDef) = -df_dw(i, j);
        }

        for (index_type j = 1; j <= iNumModes; ++j) {
            Aref(perm[i - 1], j - 1 + iFirstDofMod) = -df_da(i, j);
        }
    }
    
    const SpMapMatrixHandler& Aro = A;
    
    std::cerr << "A=[" << std::endl;
    
    for (index_type i = 1; i <= iNumNodes; ++i) {
        for (index_type j = 1; j <= iNumDofTotal; ++j) {
            std::cerr << Aro(i, j) << " ";
        }

        std::cerr << std::endl;
    }

    std::cerr << "];" << std::endl;

    std::cerr << "Aref=[" << std::endl;
    
    for (index_type i = 1; i <= iNumNodes; ++i) {
        for (index_type j = 1; j <= iNumDofTotal; ++j) {
            std::cerr << Aref(i, j) << " ";
        }

        std::cerr << std::endl;
    }
    
    std::cerr << "];" << std::endl;

    const doublereal rtol = sqrt(std::numeric_limits<doublereal>::epsilon());
    const doublereal atol = sqrt(std::numeric_limits<doublereal>::epsilon());
    
    for (index_type i = 1; i <= iNumNodes; ++i) {
        for (index_type j = 1; j <= iNumDofTotal; ++j) {
            assert(fabs(Aref(i, j) - Aro(i, j)) <= rtol * fabs(Aref(i, j)) + atol);
        }
    }
}

void ComplianceModelTest2(index_type iNumDeriv)
{
    doublereal dVal = 0.;
    std::map<index_type, doublereal> oValMap;

    LocalDofMap oDofMap;

    Vector<Gradient<0>, DYNAMIC_SIZE> v(iNumDeriv);
    Vector<doublereal, DYNAMIC_SIZE> s(iNumDeriv);
    Gradient<0> g;

    for (index_type i = 1; i <= iNumDeriv; ++i) {
        const doublereal vi = random1();
        const doublereal si = random1();
        const doublereal dvi_dx = random1();
        const index_type idx = rand() + 1;
        s(i) = si;
        dVal += si * vi;
        oValMap[idx] += si * dvi_dx;
        
        v(i).SetValuePreserve(vi);
        v(i).DerivativeResizeReset(&oDofMap, idx, MapVectorBase::GLOBAL, dvi_dx);        
    }

    std::vector<index_type> perm(iNumDeriv);

    for (index_type i = 0; i < iNumDeriv; ++i) {
        perm[i] = i + 1;
    }

    std::random_shuffle(perm.begin(), perm.end());
    
    for (index_type i = 1; i <= iNumDeriv; ++i) {
        g += s(perm[i - 1]) * v(perm[i - 1]);
    }

    const doublereal dTol = sqrt(std::numeric_limits<doublereal>::epsilon());

    std::cerr << "scalar value = " << dVal << std::endl;

    std::cerr << "map derivatives:\n";
    
    for (auto i = oValMap.begin(); i != oValMap.end(); ++i) {
        std::cerr << i->first << "->" << i->second << std::endl;
    }

    std::cerr << "gradient:\n";
    std::cerr << "value = " << g.dGetValue() << std::endl;
    std::cerr << "derivatives:\n";    
    for (index_type i = g.iGetStartIndexLocal(); i < g.iGetEndIndexLocal(); ++i) {
        const doublereal dDeriv = g.dGetDerivativeLocal(i);
        const index_type iDofIndex = g.iGetGlobalDof(i);
        
        if (dDeriv) {
            std::cerr << iDofIndex << "->" << dDeriv << std::endl;
        }
    }
    
    assert(fabs(g.dGetValue() - dVal) < dTol);

    const std::map<index_type, doublereal>& oValMapRef = oValMap;
    
    for (auto i = oValMapRef.begin(); i != oValMapRef.end(); ++i) {
        const index_type iDof = i->first;
        const doublereal dDeriv = i->second;
        const doublereal dg_dxi = g.dGetDerivativeGlobal(iDof);
        assert(fabs(dg_dxi - dDeriv) < dTol);
    }
    
    for (index_type i = g.iGetStartIndexLocal(); i < g.iGetEndIndexLocal(); ++i) {
        const doublereal dDeriv = g.dGetDerivativeLocal(i);
        const index_type iDofIndex = g.iGetGlobalDof(i);
        const auto j = oValMapRef.find(iDofIndex);
        doublereal dDerivRef = 0;
        
        if (j != oValMapRef.end()) {
            dDerivRef = j->second;
        }
        assert(fabs(dDeriv - dDerivRef) < dTol);
    }
}

void ComplianceModelTest3(index_type iNumDeriv)
{
    doublereal dVal = 0.;
    std::map<index_type, doublereal> oValMap;

    LocalDofMap oDofMap, oDofMap2;

    Vector<Gradient<0>, DYNAMIC_SIZE> v(iNumDeriv);
    Vector<doublereal, DYNAMIC_SIZE> s(iNumDeriv);

    for (index_type i = 1; i <= iNumDeriv; ++i) {
        const doublereal vi = random1();
        const doublereal si = random1();
        const doublereal dvi_dx = random1();
        const index_type idx = rand() + 1;
        s(i) = si;
        dVal += si * vi;
        oValMap[idx] += si * dvi_dx;
        
        v(i).SetValuePreserve(vi);
        v(i).DerivativeResizeReset(&oDofMap2, idx, MapVectorBase::GLOBAL, dvi_dx);        
    }

    std::vector<index_type> perm(iNumDeriv);

    for (index_type i = 0; i < iNumDeriv; ++i) {
        perm[i] = i + 1;
    }

    std::random_shuffle(perm.begin(), perm.end());

    Gradient<0> g(0., &oDofMap);
    
    for (index_type i = 1; i <= iNumDeriv; ++i) {
        g += s(perm[i - 1]) * v(perm[i - 1]);
    }

    const doublereal dTol = sqrt(std::numeric_limits<doublereal>::epsilon());

    std::cerr << "scalar value = " << dVal << std::endl;

    std::cerr << "map derivatives:\n";
    
    for (auto i = oValMap.begin(); i != oValMap.end(); ++i) {
        std::cerr << i->first << "->" << i->second << std::endl;
    }

    std::cerr << "gradient:\n";
    std::cerr << "value = " << g.dGetValue() << std::endl;
    std::cerr << "derivatives:\n";    
    for (index_type i = g.iGetStartIndexLocal(); i < g.iGetEndIndexLocal(); ++i) {
        const doublereal dDeriv = g.dGetDerivativeLocal(i);
        const index_type iDofIndex = g.iGetGlobalDof(i);
        
        if (dDeriv) {
            std::cerr << iDofIndex << "->" << dDeriv << std::endl;
        }
    }
    
    assert(fabs(g.dGetValue() - dVal) < dTol);

    const std::map<index_type, doublereal>& oValMapRef = oValMap;
    
    for (auto i = oValMapRef.begin(); i != oValMapRef.end(); ++i) {
        const index_type iDof = i->first;
        const doublereal dDeriv = i->second;
        const doublereal dg_dxi = g.dGetDerivativeGlobal(iDof);
        assert(fabs(dg_dxi - dDeriv) < dTol);
    }
    
    for (index_type i = g.iGetStartIndexLocal(); i < g.iGetEndIndexLocal(); ++i) {
        const doublereal dDeriv = g.dGetDerivativeLocal(i);
        const index_type iDofIndex = g.iGetGlobalDof(i);
        const auto j = oValMapRef.find(iDofIndex);
        doublereal dDerivRef = 0;
        
        if (j != oValMapRef.end()) {
            dDerivRef = j->second;
        }
        assert(fabs(dDeriv - dDerivRef) < dTol);
    }
}

void ComplianceModelTest4(index_type iNumDeriv)
{
    doublereal dVal = 1.;
    std::map<index_type, doublereal> oValMap;

    LocalDofMap oDofMap, oDofMap2;

    Vector<Gradient<0>, DYNAMIC_SIZE> v(iNumDeriv);
    Vector<doublereal, DYNAMIC_SIZE> s(iNumDeriv);

    for (index_type i = 1; i <= iNumDeriv; ++i) {
        const doublereal vi = random1();
        const doublereal si = random1();
        const doublereal dvi_dx = random1();
        const index_type idx = rand() + 1;
        s(i) = si;
        dVal += si * vi;
        oValMap[idx] += si * dvi_dx;

        LocalDofMap* pDofMap = random1() > 0. ? &oDofMap : &oDofMap2;
        v(i).SetValuePreserve(vi);
        v(i).DerivativeResizeReset(pDofMap, idx, MapVectorBase::GLOBAL, dvi_dx);        
    }

    std::vector<index_type> perm(iNumDeriv);

    for (index_type i = 0; i < iNumDeriv; ++i) {
        perm[i] = i + 1;
    }

    std::random_shuffle(perm.begin(), perm.end());
    
    Gradient<0> g(1., &oDofMap);
    
    for (index_type i = 1; i <= iNumDeriv; ++i) {
        g += s(perm[i - 1]) * v(perm[i - 1]);
    }

    const doublereal dTol = sqrt(std::numeric_limits<doublereal>::epsilon());

    std::cerr << "scalar value = " << dVal << std::endl;

    std::cerr << "map derivatives:\n";
    
    for (auto i = oValMap.begin(); i != oValMap.end(); ++i) {
        std::cerr << i->first << "->" << i->second << std::endl;
    }

    std::cerr << "gradient:\n";
    std::cerr << "value = " << g.dGetValue() << std::endl;
    std::cerr << "derivatives:\n";    
    for (index_type i = g.iGetStartIndexLocal(); i < g.iGetEndIndexLocal(); ++i) {
        const doublereal dDeriv = g.dGetDerivativeLocal(i);
        const index_type iDofIndex = g.iGetGlobalDof(i);
        
        if (dDeriv) {
            std::cerr << iDofIndex << "->" << dDeriv << std::endl;
        }
    }
    
    assert(fabs(g.dGetValue() - dVal) < dTol);

    const std::map<index_type, doublereal>& oValMapRef = oValMap;
    
    for (auto i = oValMapRef.begin(); i != oValMapRef.end(); ++i) {
        const index_type iDof = i->first;
        const doublereal dDeriv = i->second;
        const doublereal dg_dxi = g.dGetDerivativeGlobal(iDof);
        assert(fabs(dg_dxi - dDeriv) < dTol);
    }
    
    for (index_type i = g.iGetStartIndexLocal(); i < g.iGetEndIndexLocal(); ++i) {
        const doublereal dDeriv = g.dGetDerivativeLocal(i);
        const index_type iDofIndex = g.iGetGlobalDof(i);
        const auto j = oValMapRef.find(iDofIndex);
        doublereal dDerivRef = 0;
        
        if (j != oValMapRef.end()) {
            dDerivRef = j->second;
        }
        assert(fabs(dDeriv - dDerivRef) < dTol);
    }
}

void ComplianceModelTest5(const index_type iNumDeriv, const index_type iSliceLen)
{
    doublereal dVal = 1.;
    std::map<index_type, doublereal> oValMap;

    LocalDofMap oDofMap;

    Vector<Gradient<0>, DYNAMIC_SIZE> v(iNumDeriv);
    Vector<doublereal, DYNAMIC_SIZE> s(iNumDeriv);

    std::vector<index_type> idx(iSliceLen);
    std::vector<doublereal> dvi_dx(iSliceLen);
    
    for (index_type i = 1; i <= iNumDeriv; ++i) {
        const doublereal vi = random1();
        const index_type idx0 = rand() + 1;
        
        for (index_type j = 0; j < iSliceLen; ++j) {            
            dvi_dx[j] = random1();
            idx[j] = idx0 + j;
        }
        
        const doublereal si = random1();
                
        s(i) = si;
        dVal += si * vi;

        for (index_type j = 0; j < iSliceLen; ++j) {
            oValMap[idx[j]] += si * dvi_dx[j];
        }
        
        v(i).SetValuePreserve(vi);
        v(i).DerivativeResizeReset(&oDofMap, idx[0], idx[iSliceLen - 1] + 1, MapVectorBase::GLOBAL, 0.);

        for (index_type j = 0; j < iSliceLen; ++j) {
            v(i).SetDerivativeGlobal(idx[j], dvi_dx[j]);
        }
    }
    
    std::vector<index_type> perm(iNumDeriv);

    for (index_type i = 0; i < iNumDeriv; ++i) {
        perm[i] = i + 1;
    }

    std::random_shuffle(perm.begin(), perm.end());
    
    Gradient<0> g(1., &oDofMap);
    
    for (index_type i = 1; i <= iNumDeriv; ++i) {
        g += s(perm[i - 1]) * v(perm[i - 1]);
    }

    const doublereal dTol = sqrt(std::numeric_limits<doublereal>::epsilon());

    std::cerr << "scalar value = " << dVal << std::endl;

    std::cerr << "map derivatives:\n";
    
    for (auto i = oValMap.begin(); i != oValMap.end(); ++i) {
        std::cerr << i->first << "->" << i->second << std::endl;
    }

    std::cerr << "gradient:\n";
    std::cerr << "value = " << g.dGetValue() << std::endl;
    std::cerr << "derivatives:\n";    
    for (index_type i = g.iGetStartIndexLocal(); i < g.iGetEndIndexLocal(); ++i) {
        const doublereal dDeriv = g.dGetDerivativeLocal(i);
        const index_type iDofIndex = g.iGetGlobalDof(i);
        
        if (dDeriv) {
            std::cerr << iDofIndex << "->" << dDeriv << std::endl;
        }
    }
    
    assert(fabs(g.dGetValue() - dVal) < dTol);

    const std::map<index_type, doublereal>& oValMapRef = oValMap;
    
    for (auto i = oValMapRef.begin(); i != oValMapRef.end(); ++i) {
        const index_type iDof = i->first;
        const doublereal dDeriv = i->second;
        const doublereal dg_dxi = g.dGetDerivativeGlobal(iDof);
        assert(fabs(dg_dxi - dDeriv) < dTol);
    }
    
    for (index_type i = g.iGetStartIndexLocal(); i < g.iGetEndIndexLocal(); ++i) {
        const doublereal dDeriv = g.dGetDerivativeLocal(i);
        const index_type iDofIndex = g.iGetGlobalDof(i);
        const auto j = oValMapRef.find(iDofIndex);
        doublereal dDerivRef = 0;
        
        if (j != oValMapRef.end()) {
            dDerivRef = j->second;
        }
        assert(fabs(dDeriv - dDerivRef) < dTol);
    }
}

void ComplianceModelTest6(const index_type iNumDeriv, const index_type iSliceLen)
{
    doublereal dVal = 1.;
    std::map<index_type, doublereal> oValMap;

    LocalDofMap oDofMap, oDofMap2;

    Vector<Gradient<0>, DYNAMIC_SIZE> u(iNumDeriv), v(iNumDeriv);
    Vector<doublereal, DYNAMIC_SIZE> s(iNumDeriv);

    std::vector<index_type> vidx(iSliceLen);
    std::vector<doublereal> dvi_dx(iSliceLen);
    
    for (index_type i = 1; i <= iNumDeriv; ++i) {
        const doublereal vi = random1();
        const index_type vidx0 = rand() + 1;
        
        for (index_type j = 0; j < iSliceLen; ++j) {            
            dvi_dx[j] = random1();
            vidx[j] = vidx0 + j;
        }
        
        const doublereal si = random1();
                
        s(i) = si;
        dVal += si * vi;

        for (index_type j = 0; j < iSliceLen; ++j) {
            oValMap[vidx[j]] += si * dvi_dx[j];
        }
        
        v(i).SetValuePreserve(vi);
        v(i).DerivativeResizeReset(&oDofMap, vidx[0], vidx[iSliceLen - 1] + 1, MapVectorBase::GLOBAL, 0.);

        for (index_type j = 0; j < iSliceLen; ++j) {
            v(i).SetDerivativeGlobal(vidx[j], dvi_dx[j]);
        }
    }

    std::vector<doublereal> dui_dx(iSliceLen);
    std::vector<index_type> uidx(iSliceLen);
    
    for (index_type i = 1; i <= iNumDeriv; ++i) {
        const doublereal ui = random1();
        const index_type uidx0 = rand() + 1;
        
        for (index_type j = 0; j < iSliceLen; ++j) {            
            dui_dx[j] = random1();
            uidx[j] = uidx0 + j;
        }
                        
        dVal += s(i) * ui;

        for (index_type j = 0; j < iSliceLen; ++j) {
            oValMap[uidx[j]] += s(i) * dui_dx[j];
        }
        
        u(i).SetValuePreserve(ui);
        u(i).DerivativeResizeReset(&oDofMap2, uidx[0], uidx[iSliceLen - 1] + 1, MapVectorBase::GLOBAL, 0.);

        for (index_type j = 0; j < iSliceLen; ++j) {
            u(i).SetDerivativeGlobal(uidx[j], dui_dx[j]);
        }
        }
    
    std::vector<index_type> perm(iNumDeriv);

    for (index_type i = 0; i < iNumDeriv; ++i) {
        perm[i] = i + 1;
    }

    std::random_shuffle(perm.begin(), perm.end());
    
    Gradient<0> g(1., &oDofMap);
    
    for (index_type i = 1; i <= iNumDeriv; ++i) {
        g += s(perm[i - 1]) * v(perm[i - 1]);
        g += s(perm[i - 1]) * u(perm[i - 1]);
}

    const doublereal dTol = sqrt(std::numeric_limits<doublereal>::epsilon());

    std::cerr << "scalar value = " << dVal << std::endl;

    std::cerr << "map derivatives:\n";
    
    for (auto i = oValMap.begin(); i != oValMap.end(); ++i) {
        std::cerr << i->first << "->" << i->second << std::endl;
    }

    std::cerr << "gradient:\n";
    std::cerr << "value = " << g.dGetValue() << std::endl;
    std::cerr << "derivatives:\n";    
    for (index_type i = g.iGetStartIndexLocal(); i < g.iGetEndIndexLocal(); ++i) {
        const doublereal dDeriv = g.dGetDerivativeLocal(i);
        const index_type iDofIndex = g.iGetGlobalDof(i);
        
        if (dDeriv) {
            std::cerr << iDofIndex << "->" << dDeriv << std::endl;
        }
    }
    
    assert(fabs(g.dGetValue() - dVal) < dTol);

    const std::map<index_type, doublereal>& oValMapRef = oValMap;
    
    for (auto i = oValMapRef.begin(); i != oValMapRef.end(); ++i) {
        const index_type iDof = i->first;
        const doublereal dDeriv = i->second;
        const doublereal dg_dxi = g.dGetDerivativeGlobal(iDof);
        assert(fabs(dg_dxi - dDeriv) < dTol);
    }
    
    for (index_type i = g.iGetStartIndexLocal(); i < g.iGetEndIndexLocal(); ++i) {
        const doublereal dDeriv = g.dGetDerivativeLocal(i);
        const index_type iDofIndex = g.iGetGlobalDof(i);
        const auto j = oValMapRef.find(iDofIndex);
        doublereal dDerivRef = 0;
        
        if (j != oValMapRef.end()) {
            dDerivRef = j->second;
        }
        assert(fabs(dDeriv - dDerivRef) < dTol);
    }
}

void ComplianceModelTest7(const index_type iNumDeriv, const index_type iSliceLen)
{
    doublereal dVal = 1.;
    std::map<index_type, doublereal> oValMap;

    LocalDofMap oDofMap, oDofMap2;

    Vector<Gradient<0>, DYNAMIC_SIZE> u(iNumDeriv), v(iNumDeriv);
    Vector<doublereal, DYNAMIC_SIZE> s(iNumDeriv);

    std::vector<index_type> vidx(iSliceLen);
    std::vector<doublereal> dvi_dx(iSliceLen);
    
    for (index_type i = 1; i <= iNumDeriv; ++i) {
        const doublereal vi = random1();
        const index_type vidx0 = rand() + 1;
        
        for (index_type j = 0; j < iSliceLen; ++j) {            
            dvi_dx[j] = random1();
            vidx[j] = vidx0 + j;
        }
        
        const doublereal si = random1();
                
        s(i) = si;
        dVal += si * vi;

        for (index_type j = 0; j < iSliceLen; ++j) {
            oValMap[vidx[j]] += si * dvi_dx[j];
        }
        
        v(i).SetValuePreserve(vi);
        v(i).DerivativeResizeReset(&oDofMap, vidx[0], vidx[iSliceLen - 1] + 1, MapVectorBase::GLOBAL, 0.);

        for (index_type j = 0; j < iSliceLen; ++j) {
            v(i).SetDerivativeGlobal(vidx[j], dvi_dx[j]);
        }
    }

    std::vector<doublereal> dui_dx(iSliceLen);
    std::vector<index_type> uidx(iSliceLen);
    
    for (index_type i = 1; i <= iNumDeriv; ++i) {
        const doublereal ui = random1();
        const index_type uidx0 = rand() + 1;
        
        for (index_type j = 0; j < iSliceLen; ++j) {            
            dui_dx[j] = random1();
            uidx[j] = uidx0 + j;
        }
                        
        dVal += s(i) * ui;

        for (index_type j = 0; j < iSliceLen; ++j) {
            oValMap[uidx[j]] += s(i) * dui_dx[j];
        }
        
        u(i).SetValuePreserve(ui);
        u(i).DerivativeResizeReset(&oDofMap2, uidx[0], uidx[iSliceLen - 1] + 1, MapVectorBase::GLOBAL, 0.);

        for (index_type j = 0; j < iSliceLen; ++j) {
            u(i).SetDerivativeGlobal(uidx[j], dui_dx[j]);
        }

        for (index_type j = 0; j < iSliceLen; ++j) {
            doublereal wi = random1();
            doublereal dwi_dx = random1();
            index_type widx = rand() + 1;
        
            dVal += s(i) * wi;
            oValMap[widx] += s(i) * dwi_dx;

            Gradient<0> w;

            w.SetValuePreserve(wi);
            w.DerivativeResizeReset(&oDofMap2, widx, MapVectorBase::GLOBAL, dwi_dx);
            u(i) += w;
        }
    }
    
    std::vector<index_type> perm1(iNumDeriv), perm2(iNumDeriv);

    for (index_type i = 0; i < iNumDeriv; ++i) {
        perm1[i] = perm2[i] = i + 1;
    }

    std::random_shuffle(perm1.begin(), perm1.end());
    std::random_shuffle(perm2.begin(), perm2.end());
    
    Gradient<0> g(1., &oDofMap);
    
    for (index_type i = 1; i <= iNumDeriv; ++i) {
        g += s(perm1[i - 1]) * v(perm1[i - 1]);
        g += s(perm2[i - 1]) * u(perm2[i - 1]);
    }
    
    const doublereal dTol = sqrt(std::numeric_limits<doublereal>::epsilon());

    std::cerr << "scalar value = " << dVal << std::endl;

    std::cerr << "map derivatives:\n";
    
    for (auto i = oValMap.begin(); i != oValMap.end(); ++i) {
        std::cerr << i->first << "->" << i->second << std::endl;
    }

    std::cerr << "gradient:\n";
    std::cerr << "value = " << g.dGetValue() << std::endl;
    std::cerr << "derivatives:\n";    
    for (index_type i = g.iGetStartIndexLocal(); i < g.iGetEndIndexLocal(); ++i) {
        const doublereal dDeriv = g.dGetDerivativeLocal(i);
        const index_type iDofIndex = g.iGetGlobalDof(i);
        
        if (dDeriv) {
            std::cerr << iDofIndex << "->" << dDeriv << std::endl;
        }
    }

    const std::map<index_type, doublereal>& oValMapRef = oValMap;
    
    for (auto i = oValMapRef.begin(); i != oValMapRef.end(); ++i) {
        const index_type iDof = i->first;
        const doublereal dDeriv = i->second;
        const doublereal dg_dxi = g.dGetDerivativeGlobal(iDof);
        assert(fabs(dg_dxi - dDeriv) < dTol);
    }
    
    for (index_type i = g.iGetStartIndexLocal(); i < g.iGetEndIndexLocal(); ++i) {
        const doublereal dDeriv = g.dGetDerivativeLocal(i);
        const index_type iDofIndex = g.iGetGlobalDof(i);
        const auto j = oValMapRef.find(iDofIndex);
        doublereal dDerivRef = 0;
        
        if (j != oValMapRef.end()) {
            dDerivRef = j->second;
        }
        assert(fabs(dDeriv - dDerivRef) < dTol);
    }

    assert(fabs(g.dGetValue() - dVal) < dTol);
}

void ComplianceModelTest8(const index_type iNumDeriv, const index_type iSliceLen)
{
    doublereal dVal = 1.;
    std::map<index_type, doublereal> oValMap;

    LocalDofMap oDofMap, oDofMap2, oDofMap3;

    Vector<Gradient<0>, DYNAMIC_SIZE> u(iNumDeriv), v(iNumDeriv), u2(iNumDeriv);
    Vector<doublereal, DYNAMIC_SIZE> s(iNumDeriv);

    std::vector<index_type> vidx(iSliceLen);
    std::vector<doublereal> dvi_dx(iSliceLen);
    
    for (index_type i = 1; i <= iNumDeriv; ++i) {
        const doublereal vi = random1();
        const index_type vidx0 = rand() + 1;
        
        for (index_type j = 0; j < iSliceLen; ++j) {            
            dvi_dx[j] = random1();
            vidx[j] = vidx0 + j;
        }
        
        const doublereal si = random1();
                
        s(i) = si;
        dVal += si * vi;

        for (index_type j = 0; j < iSliceLen; ++j) {
            oValMap[vidx[j]] += si * dvi_dx[j];
        }
        
        v(i).SetValuePreserve(vi);
        v(i).DerivativeResizeReset(&oDofMap, vidx[0], vidx[iSliceLen - 1] + 1, MapVectorBase::GLOBAL, 0.);

        for (index_type j = 0; j < iSliceLen; ++j) {
            v(i).SetDerivativeGlobal(vidx[j], dvi_dx[j]);
        }
    }

    std::vector<doublereal> dui_dx(iSliceLen);
    std::vector<index_type> uidx(iSliceLen);
    
    for (index_type i = 1; i <= iNumDeriv; ++i) {
        const doublereal ui = random1();
        const index_type uidx0 = rand() + 1;
        
        for (index_type j = 0; j < iSliceLen; ++j) {            
            dui_dx[j] = random1();
            uidx[j] = uidx0 + j;
        }
                        
        dVal += s(i) * ui;

        for (index_type j = 0; j < iSliceLen; ++j) {
            oValMap[uidx[j]] += s(i) * dui_dx[j];
        }
        
        u(i).SetValuePreserve(ui);
        u(i).DerivativeResizeReset(&oDofMap2, uidx[0], uidx[iSliceLen - 1] + 1, MapVectorBase::GLOBAL, 0.);

        for (index_type j = 0; j < iSliceLen; ++j) {
            u(i).SetDerivativeGlobal(uidx[j], dui_dx[j]);
        }

        for (index_type j = 0; j < iSliceLen; ++j) {
            doublereal wi = random1();
            doublereal dwi_dx = random1();
            index_type widx = rand() + 1;
        
            dVal += s(i) * wi;
            oValMap[widx] += s(i) * dwi_dx;

            Gradient<0> w;

            w.SetValuePreserve(wi);
            w.DerivativeResizeReset(&oDofMap2, widx, MapVectorBase::GLOBAL, dwi_dx);
            u(i) += w;
        }
    }

    std::vector<index_type> perm1(iNumDeriv), perm2(iNumDeriv), perm3(iNumDeriv);

    for (index_type i = 0; i < iNumDeriv; ++i) {
        perm1[i] = perm2[i] = perm3[i] = i + 1;
    }

    std::random_shuffle(perm1.begin(), perm1.end());
    std::random_shuffle(perm2.begin(), perm2.end());
    std::random_shuffle(perm3.begin(), perm3.end());
    
    for (index_type i = 1; i <= iNumDeriv; ++i) {
        u2(perm3[i - 1]).Copy(u(perm3[i - 1]), &oDofMap3);
    }
    
    Gradient<0> g(1., &oDofMap);
    
    for (index_type i = 1; i <= iNumDeriv; ++i) {
        g += s(perm1[i - 1]) * v(perm1[i - 1]);
        g += s(perm2[i - 1]) * u2(perm2[i - 1]);
    }
    
    const doublereal dTol = sqrt(std::numeric_limits<doublereal>::epsilon());

    std::cerr << "scalar value = " << dVal << std::endl;

    std::cerr << "map derivatives:\n";
    
    for (auto i = oValMap.begin(); i != oValMap.end(); ++i) {
        std::cerr << i->first << "->" << i->second << std::endl;
    }

    std::cerr << "gradient:\n";
    std::cerr << "value = " << g.dGetValue() << std::endl;
    std::cerr << "derivatives:\n";    
    for (index_type i = g.iGetStartIndexLocal(); i < g.iGetEndIndexLocal(); ++i) {
        const doublereal dDeriv = g.dGetDerivativeLocal(i);
        const index_type iDofIndex = g.iGetGlobalDof(i);
        
        if (dDeriv) {
            std::cerr << iDofIndex << "->" << dDeriv << std::endl;
        }
    }

    const std::map<index_type, doublereal>& oValMapRef = oValMap;
    
    for (auto i = oValMapRef.begin(); i != oValMapRef.end(); ++i) {
        const index_type iDof = i->first;
        const doublereal dDeriv = i->second;
        const doublereal dg_dxi = g.dGetDerivativeGlobal(iDof);
        assert(fabs(dg_dxi - dDeriv) < dTol);
    }
    
    for (index_type i = g.iGetStartIndexLocal(); i < g.iGetEndIndexLocal(); ++i) {
        const doublereal dDeriv = g.dGetDerivativeLocal(i);
        const index_type iDofIndex = g.iGetGlobalDof(i);
        const auto j = oValMapRef.find(iDofIndex);
        doublereal dDerivRef = 0;
        
        if (j != oValMapRef.end()) {
            dDerivRef = j->second;
        }
        assert(fabs(dDeriv - dDerivRef) < dTol);
    }

    assert(fabs(g.dGetValue() - dVal) < dTol);
}

void ComplianceModelTest9(const index_type iNumDeriv, const index_type iSliceLen)
{
    doublereal dVal = 1.;
    std::map<index_type, doublereal> oValMap;

    LocalDofMap oDofMap, oDofMap2, oDofMap3, oDofMap4;

    Vector<Gradient<0>, DYNAMIC_SIZE> u(iNumDeriv), u2(iNumDeriv), v(iNumDeriv), v2(iNumDeriv);
    Vector<doublereal, DYNAMIC_SIZE> s(iNumDeriv);

    std::vector<index_type> vidx(iSliceLen);
    std::vector<doublereal> dvi_dx(iSliceLen);
    
    for (index_type i = 1; i <= iNumDeriv; ++i) {
        const doublereal vi = random1();
        const index_type vidx0 = rand() + 1;
        
        for (index_type j = 0; j < iSliceLen; ++j) {            
            dvi_dx[j] = random1();
            vidx[j] = vidx0 + j;
        }
        
        const doublereal si = random1();
                
        s(i) = si;
        dVal += si * vi;

        for (index_type j = 0; j < iSliceLen; ++j) {
            oValMap[vidx[j]] += si * dvi_dx[j];
        }
        
        v(i).SetValuePreserve(vi);
        v(i).DerivativeResizeReset(&oDofMap, vidx[0], vidx[iSliceLen - 1] + 1, MapVectorBase::GLOBAL, 0.);

        for (index_type j = 0; j < iSliceLen; ++j) {
            v(i).SetDerivativeGlobal(vidx[j], dvi_dx[j]);
        }
    }

    std::vector<doublereal> dui_dx(iSliceLen);
    std::vector<index_type> uidx(iSliceLen);
    
    for (index_type i = 1; i <= iNumDeriv; ++i) {
        const doublereal ui = random1();
        const index_type uidx0 = rand() + 1;
        
        for (index_type j = 0; j < iSliceLen; ++j) {            
            dui_dx[j] = random1();
            uidx[j] = uidx0 + j;
        }
                        
        dVal += s(i) * ui;

        for (index_type j = 0; j < iSliceLen; ++j) {
            oValMap[uidx[j]] += s(i) * dui_dx[j];
        }
        
        u(i).SetValuePreserve(ui);
        u(i).DerivativeResizeReset(&oDofMap2, uidx[0], uidx[iSliceLen - 1] + 1, MapVectorBase::GLOBAL, 0.);

        for (index_type j = 0; j < iSliceLen; ++j) {
            u(i).SetDerivativeGlobal(uidx[j], dui_dx[j]);
        }

        for (index_type j = 0; j < iSliceLen; ++j) {
            doublereal wi = random1();
            doublereal dwi_dx = random1();
            index_type widx = rand() + 1;
        
            dVal += s(i) * wi;
            oValMap[widx] += s(i) * dwi_dx;

            Gradient<0> w;

            w.SetValuePreserve(wi);
            w.DerivativeResizeReset(&oDofMap2, widx, MapVectorBase::GLOBAL, dwi_dx);
            u(i) += w;
        }
    }

    std::vector<index_type> perm1(iNumDeriv), perm2(iNumDeriv), perm3(iNumDeriv);

    for (index_type i = 0; i < iNumDeriv; ++i) {
        perm1[i] = perm2[i] = perm3[i] = i + 1;
    }

    std::random_shuffle(perm1.begin(), perm1.end());
    std::random_shuffle(perm2.begin(), perm2.end());
    std::random_shuffle(perm3.begin(), perm3.end());
    
    for (index_type i = 1; i <= iNumDeriv; ++i) {
        u2(perm3[i - 1]).Copy(u(perm3[i - 1]), &oDofMap3);
        v2(perm1[i - 1]).Copy(v(perm1[i - 1]), &oDofMap4);
    }
    
    Gradient<0> g(1., &oDofMap);
    
    for (index_type i = 1; i <= iNumDeriv; ++i) {
        g += s(perm1[i - 1]) * v2(perm1[i - 1]);
        g += s(perm2[i - 1]) * u2(perm2[i - 1]);
    }
    
    const doublereal dTol = sqrt(std::numeric_limits<doublereal>::epsilon());

    std::cerr << "scalar value = " << dVal << std::endl;

    std::cerr << "map derivatives:\n";
    
    for (auto i = oValMap.begin(); i != oValMap.end(); ++i) {
        std::cerr << i->first << "->" << i->second << std::endl;
    }

    std::cerr << "gradient:\n";
    std::cerr << "value = " << g.dGetValue() << std::endl;
    std::cerr << "derivatives:\n";    
    for (index_type i = g.iGetStartIndexLocal(); i < g.iGetEndIndexLocal(); ++i) {
        const doublereal dDeriv = g.dGetDerivativeLocal(i);
        const index_type iDofIndex = g.iGetGlobalDof(i);
        
        if (dDeriv) {
            std::cerr << iDofIndex << "->" << dDeriv << std::endl;
        }
    }

    const std::map<index_type, doublereal>& oValMapRef = oValMap;
    
    for (auto i = oValMapRef.begin(); i != oValMapRef.end(); ++i) {
        const index_type iDof = i->first;
        const doublereal dDeriv = i->second;
        const doublereal dg_dxi = g.dGetDerivativeGlobal(iDof);
        assert(fabs(dg_dxi - dDeriv) < dTol);
    }
    
    for (index_type i = g.iGetStartIndexLocal(); i < g.iGetEndIndexLocal(); ++i) {
        const doublereal dDeriv = g.dGetDerivativeLocal(i);
        const index_type iDofIndex = g.iGetGlobalDof(i);
        const auto j = oValMapRef.find(iDofIndex);
        doublereal dDerivRef = 0;
        
        if (j != oValMapRef.end()) {
            dDerivRef = j->second;
        }
        assert(fabs(dDeriv - dDerivRef) < dTol);
    }

    assert(fabs(g.dGetValue() - dVal) < dTol);
}

void RangeVectorTest1(index_type N)
{
    RangeVector<doublereal, 0> oRange;

    for (index_type k = 0; k < N; ++k) {
        index_type j = rand() % 1000;
        index_type m = std::min(j, rand() % 100);
        index_type n = rand() % 100;
        
        oRange.ResizeReset(j, j + 1, 1);

        std::cerr << "range=";
        for (index_type i = oRange.iGetStartIndex(); i < oRange.iGetEndIndex(); ++i) {
            std::cerr << i << ":" << oRange.GetValue(i) << std::endl;
        }
        std::cerr << std::endl;
    
        oRange.ResizePreserve(j - m, j + n + 1);

        std::cerr << "range=";
        for (index_type i = oRange.iGetStartIndex(); i < oRange.iGetEndIndex(); ++i) {
            std::cerr << i << ":" << oRange.GetValue(i) << std::endl;
        }
        std::cerr << std::endl;

        for (index_type i = oRange.iGetStartIndex(); i < oRange.iGetEndIndex(); ++i) {
            if (i == j) {
                assert(oRange.GetValue(i) == 1);
            } else {
                assert(oRange.GetValue(i) == 0);
            }
        }
    }
}

void RangeVectorTest2(index_type N)
{
    RangeVector<doublereal, 0> oRange;

    for (index_type k = 0; k < N; ++k) {
        index_type o = std::max(1, rand() % 100);

        std::set<index_type> sidx;
        
        for (index_type l = 0; l < o; ++l) {
            sidx.insert(rand() % 1000);
        }

        o = sidx.size();
        
        std::vector<index_type> idx;
        std::vector<doublereal> v;

        idx.reserve(o);
        v.reserve(o);
        
        for (auto l = sidx.begin(); l != sidx.end(); ++l) {
            idx.push_back(*l);
            v.push_back(random1());
        }        
        
        index_type jmin = *std::min_element(idx.begin(), idx.end());
        index_type jmax = *std::max_element(idx.begin(), idx.end());
        index_type m = std::min(jmin, rand() % 100);
        index_type n = rand() % 100;
        
        oRange.ResizeReset(jmin, jmax + 1, 0.);

        for (index_type l = 0; l < o; ++l) {
            oRange.SetValue(idx[l], v[l]);
        }

        std::cerr << "range=";
        for (index_type i = oRange.iGetStartIndex(); i < oRange.iGetEndIndex(); ++i) {
            std::cerr << i << ":" << oRange.GetValue(i) << std::endl;
        }
        std::cerr << std::endl;
    
        oRange.ResizePreserve(jmin - m, jmax + n + 1);

        std::cerr << "range=";
        for (index_type i = oRange.iGetStartIndex(); i < oRange.iGetEndIndex(); ++i) {
            std::cerr << i << ":" << oRange.GetValue(i) << std::endl;
        }
        std::cerr << std::endl;

        for (index_type i = oRange.iGetStartIndex(); i < oRange.iGetEndIndex(); ++i) {
            bool bFound = false;
            
            for (index_type l = 0; l < o; ++l) {
                if (idx[l] == i) {
                    bFound = true;
                    assert(oRange.GetValue(i) == v[l]);
                }
            }

            if (!bFound) {
                assert(oRange.GetValue(i) == 0.);
            }
        }
    }
}
#endif

int main(int argc, char* argv[]) {
#ifdef HAVE_FEENABLEEXCEPT
	feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
#endif
	if (argc > 1) {
		NLoops = atol(argv[1]);
	}

    if (argc > 2) {
        NLoopsAss = atol(argv[2]);
    }

    if (NLoops < 1) {
        NLoops = 1;
    }

    if (NLoopsAss < 1) {
        NLoopsAss = 1;
    }

#ifdef TEST_GROUP1
    std::cerr << "MatManip_test()\n";
    
    MatManip_test(NLoops);
    
	std::cerr << "---------------------------\ntestScalarTypeTraits()\n";
	testScalarTypeTraits();

	testMatVec<3>();
	testMatVec<6>();
	testMatVec<8>();
	testMatVec<12>();
	testMatVec<24>();

    std::cerr << "---------------------------\ntestMatVecGradient2<1>()\n";
	testMatVecGradient2<1>();

    std::cerr << "---------------------------\ntestMatVecGradient2<2>()\n";
	testMatVecGradient2<2>();

    std::cerr << "---------------------------\ntestMatVecGradient2<3>()\n";
	testMatVecGradient2<3>();

    std::cerr << "---------------------------\ntestMatVecGradient2<4>()\n";
	testMatVecGradient2<4>();

    std::cerr << "---------------------------\ntestMatVecGradient2<5>()\n";
	testMatVecGradient2<5>();

    std::cerr << "---------------------------\ntestMatVecGradient2<6>()\n";
	testMatVecGradient2<6>();

    std::cerr << "---------------------------\ntestMatVecGradient2<8>()\n";
	testMatVecGradient2<8>();

    std::cerr << "---------------------------\ntestMatVecGradient2<10>()\n";
	testMatVecGradient2<10>();

    std::cerr << "---------------------------\ntestMatVecGradient2<12>()\n";
	testMatVecGradient2<12>();

	std::cerr << "---------------------------\ntestMatVecDouble2()\n";
	testMatVecDouble2();

    std::cerr << "---------------------------\ntestMatVecProduct()\n";
    testMatVecProduct();
    std::cerr << "---------------------------\ntestMatVecProductGradient()\n";
    testMatVecProductGradient();

    std::cerr << "---------------------------\ntestMatVecProductGradient2()\n";
    testMatVecProductGradient2<1>(1, NLoops);
    testMatVecProductGradient2<4>(4, NLoops);
    testMatVecProductGradient2<8>(8, NLoops);
    testMatVecProductGradient2<12>(12, NLoops);
    testMatVecProductGradient2<32>(32, NLoops);
    testMatVecProductGradient2<64>(64, NLoops);
    testMatVecProductGradient2<0>(256, NLoops);
    
    std::cerr << "---------------------------\ntestMatVecCopy()\n";
    testMatVecCopy<1>();
    testMatVecCopy<3>();
    testMatVecCopy<5>();
    testMatVecCopy<9>();

#ifdef HAVE_BLITZ
    std::cerr << "---------------------------\ntestAssembly()\n";
    gradVecAssTest::testAssembly();
#endif

    std::cerr << "---------------------------\ntestMatVec3()\n";
    testMatVec3();

    std::cerr << "---------------------------\ntestSubVecAss()\n";
    testSubVecAss();

    std::cerr << "---------------------------\ntestSubVecAssMatVec()\n";
    testSubVecAssMatVec();

    std::cerr << "---------------------------\ntestInv()\n";
    testInv();

    std::cerr << "----------------------------\ntestSolve()\n";
    testSolve();

      {  
            const int N = argc > 3 ? atoi(argv[3]) : 1;
            const int nr = 4, nc = 4;
            Matrix<doublereal, nr, nc> A, B, C;

            for (int i = 1; i <= nr; ++i)
            {
                for (int j = 1; j <= nc; ++j)
                {
                    B(i, j) = rand();
                    C(i, j) = rand();
                }
            }

            clock_t start = clock();

            for (int i = 0; i < N; ++i)
            {
                A = B * C;
            }

            std::cerr << "time: " << double(clock()-start)/N/CLOCKS_PER_SEC << std::endl;
      }
      
    std::cerr << "----------------------------\ncppad_benchmark1<0>()\n";
    cppad_benchmark1<0>(NLoops);

    std::cerr << "----------------------------\ncppad_benchmark1<6>()\n";
    cppad_benchmark1<6>(NLoops);

    std::cerr << "----------------------------\ncppad_benchmark2<0>()\n";
    cppad_benchmark2<0>(NLoops);

    std::cerr << "----------------------------\ncppad_benchmark2<12>()\n";
    cppad_benchmark2<12>(NLoops);

    std::cerr << "----------------------------\ncppad_benchmark3<0>()\n";
    cppad_benchmark3<0>(NLoops);

    std::cerr << "----------------------------\ncppad_benchmark3<12>()\n";
    cppad_benchmark3<12>(NLoops);
      
#ifdef HAVE_FC_F77
    testVecOp<3, 2>(NLoops, 2, doVecAdd<Gradient<2>, 3>, __FC_DECL__(func2addad_dv), "add");
	testVecOp<3, 4>(NLoops, 4, doVecAdd<Gradient<4>, 3>, __FC_DECL__(func2addad_dv), "add");
	testVecOp<3, 8>(NLoops, 8, doVecAdd<Gradient<8>, 3>, __FC_DECL__(func2addad_dv), "add");
	testVecOp<3, 16>(NLoops, 16, doVecAdd<Gradient<16>, 3>, __FC_DECL__(func2addad_dv), "add");

    testVecOp<12, 2>(NLoops, 2, doVecAdd<Gradient<2>, 12>, __FC_DECL__(func2addad_dv), "add");
	testVecOp<12, 4>(NLoops, 4, doVecAdd<Gradient<4>, 12>, __FC_DECL__(func2addad_dv), "add");
	testVecOp<12, 8>(NLoops, 8, doVecAdd<Gradient<8>, 12>, __FC_DECL__(func2addad_dv), "add");
	testVecOp<12, 16>(NLoops, 16, doVecAdd<Gradient<16>, 12>, __FC_DECL__(func2addad_dv), "add");

    testVecOp<3, 2>(NLoops, 2, doVecMul<Gradient<2>, 3>, __FC_DECL__(func2mulad_dv), "mul");
	testVecOp<3, 4>(NLoops, 4, doVecMul<Gradient<4>, 3>, __FC_DECL__(func2mulad_dv), "mul");
	testVecOp<3, 8>(NLoops, 8, doVecMul<Gradient<8>, 3>, __FC_DECL__(func2mulad_dv), "mul");
	testVecOp<3, 16>(NLoops, 16, doVecMul<Gradient<16>, 3>, __FC_DECL__(func2mulad_dv), "mul");

    testVecOp<12, 2>(NLoops, 2, doVecAdd<Gradient<2>, 12>, __FC_DECL__(func2addad_dv), "add");
	testVecOp<12, 4>(NLoops, 4, doVecAdd<Gradient<4>, 12>, __FC_DECL__(func2addad_dv), "add");
	testVecOp<12, 8>(NLoops, 8, doVecAdd<Gradient<8>, 12>, __FC_DECL__(func2addad_dv), "add");
	testVecOp<12, 16>(NLoops, 16, doVecAdd<Gradient<16>, 12>, __FC_DECL__(func2addad_dv), "add");

    testVecOp<3, 0>(NLoops, 2, doVecAdd<Gradient<0>, 3>, __FC_DECL__(func2addad_dv), "add");
    testVecOp<3, 0>(NLoops, 4, doVecAdd<Gradient<0>, 3>, __FC_DECL__(func2addad_dv), "add");
    testVecOp<3, 0>(NLoops, 8, doVecAdd<Gradient<0>, 3>, __FC_DECL__(func2addad_dv), "add");
    testVecOp<3, 0>(NLoops, 16, doVecAdd<Gradient<0>, 3>, __FC_DECL__(func2addad_dv), "add");
    testVecOp<3, 0>(NLoops, 256, doVecAdd<Gradient<0>, 3>, __FC_DECL__(func2addad_dv), "add");
    testVecOp<3, 0>(NLoops, 512, doVecAdd<Gradient<0>, 3>, __FC_DECL__(func2addad_dv), "add");
    testVecOp<3, 0>(NLoops, 1024, doVecAdd<Gradient<0>, 3>, __FC_DECL__(func2addad_dv), "add");

    testVecOp<12, 0>(NLoops, 2, doVecAdd<Gradient<0>, 12>, __FC_DECL__(func2addad_dv), "add");
    testVecOp<12, 0>(NLoops, 4, doVecAdd<Gradient<0>, 12>, __FC_DECL__(func2addad_dv), "add");
    testVecOp<12, 0>(NLoops, 8, doVecAdd<Gradient<0>, 12>, __FC_DECL__(func2addad_dv), "add");
    testVecOp<12, 0>(NLoops, 16, doVecAdd<Gradient<0>, 12>, __FC_DECL__(func2addad_dv), "add");
    testVecOp<12, 0>(NLoops, 256, doVecAdd<Gradient<0>, 12>, __FC_DECL__(func2addad_dv), "add");
    testVecOp<12, 0>(NLoops, 512, doVecAdd<Gradient<0>, 12>, __FC_DECL__(func2addad_dv), "add");
    testVecOp<12, 0>(NLoops, 1024, doVecAdd<Gradient<0>, 12>, __FC_DECL__(func2addad_dv), "add");

    testVecOp<3, 0>(NLoops, 2, doVecMul<Gradient<0>, 3>, __FC_DECL__(func2mulad_dv), "mul");
    testVecOp<3, 0>(NLoops, 4, doVecMul<Gradient<0>, 3>, __FC_DECL__(func2mulad_dv), "mul");
    testVecOp<3, 0>(NLoops, 8, doVecMul<Gradient<0>, 3>, __FC_DECL__(func2mulad_dv), "mul");
    testVecOp<3, 0>(NLoops, 16, doVecMul<Gradient<0>, 3>, __FC_DECL__(func2mulad_dv), "mul");
    testVecOp<3, 0>(NLoops, 256, doVecMul<Gradient<0>, 3>, __FC_DECL__(func2mulad_dv), "mul");
    testVecOp<3, 0>(NLoops, 512, doVecMul<Gradient<0>, 3>, __FC_DECL__(func2mulad_dv), "mul");
    testVecOp<3, 0>(NLoops, 1024, doVecMul<Gradient<0>, 3>, __FC_DECL__(func2mulad_dv), "mul");

    testVecOp<12, 0>(NLoops, 2, doVecMul<Gradient<0>, 12>, __FC_DECL__(func2mulad_dv), "mul");
    testVecOp<12, 0>(NLoops, 4, doVecMul<Gradient<0>, 12>, __FC_DECL__(func2mulad_dv), "mul");
    testVecOp<12, 0>(NLoops, 8, doVecMul<Gradient<0>, 12>, __FC_DECL__(func2mulad_dv), "mul");
    testVecOp<12, 0>(NLoops, 16, doVecMul<Gradient<0>, 12>, __FC_DECL__(func2mulad_dv), "mul");
    testVecOp<12, 0>(NLoops, 256, doVecMul<Gradient<0>, 12>, __FC_DECL__(func2mulad_dv), "mul");
    testVecOp<12, 0>(NLoops, 512, doVecMul<Gradient<0>, 12>, __FC_DECL__(func2mulad_dv), "mul");
    testVecOp<12, 0>(NLoops, 1024, doVecMul<Gradient<0>, 12>, __FC_DECL__(func2mulad_dv), "mul");
#endif // HAVE_FC_F77

    Mat3xN_test(3, NLoops);
    Mat3xN_test(10, NLoops);
    Mat3xN_test(100, NLoops);

    MatNxN_test(3, NLoops);
    MatNxN_test(10, NLoops);
    MatNxN_test(100, NLoops);

    Mat3xN_test_grad<4>(4, 3, NLoops);
    Mat3xN_test_grad<5>(5, 10, NLoops);
    Mat3xN_test_grad<0>(20, 100, NLoops);

    Mat3xNT_test_grad<4>(4, 3, NLoops);
    Mat3xNT_test_grad<5>(5, 10, NLoops);
    Mat3xNT_test_grad<0>(20, 100, NLoops);
    
    MatNxN_test_grad<4>(4, 3, NLoops);
    MatNxN_test_grad<5>(5, 10, NLoops);
    MatNxN_test_grad<0>(20, 100, NLoops);

    MatNxNT_test_grad<4>(4, 3, NLoops);
    MatNxNT_test_grad<5>(5, 10, NLoops);
    MatNxNT_test_grad<0>(20, 100, NLoops);
    
    MatDynamic_test(3, 4, NLoops);
    MatDynamic_test(10, 20, NLoops);
    
    MatDynamic_test_grad<3>(3, 4, 5, NLoops);
    MatDynamic_test_grad<5>(5, 10, 20, NLoops);
    MatDynamic_test_grad<0>(15, 20, 30, NLoops);

    MatDynamicT_test<10, 20>(10, 20, NLoops);
    MatDynamicT_test<DYNAMIC_SIZE, 20>(10, 20, NLoops);
    MatDynamicT_test<10, DYNAMIC_SIZE>(10, 20, NLoops);
    MatDynamicT_test<DYNAMIC_SIZE, DYNAMIC_SIZE>(10, 20, NLoops);

    MatDynamicT_test_grad<5, 10, 20>(5, 10, 20, NLoops);
    MatDynamicT_test_grad<5, DYNAMIC_SIZE, 20>(5, 10, 20, NLoops);
    MatDynamicT_test_grad<5, 10, DYNAMIC_SIZE>(5, 10, 20, NLoops);
    MatDynamicT_test_grad<5, DYNAMIC_SIZE, DYNAMIC_SIZE>(5, 10, 20, NLoops);
    MatDynamicT_test_grad<0, 10, 20>(5, 10, 20, NLoops);
    MatDynamicT_test_grad<0, DYNAMIC_SIZE, 20>(5, 10, 20, NLoops);
    MatDynamicT_test_grad<0, 10, DYNAMIC_SIZE>(5, 10, 20, NLoops);
    MatDynamicT_test_grad<0, DYNAMIC_SIZE, DYNAMIC_SIZE>(5, 10, 20, NLoops);
    
    
    MatrixMatrixProduct_test<2, 3, 4, MatIndexValGen>(2, 3, 4);
    MatrixMatrixProduct_test<4, 2, 3, MatIndexValGen>(4, 2, 3);
    MatrixMatrixProduct_test<3, 4, 2, MatIndexValGen>(3, 4, 2);
    MatrixMatrixProduct_test<10, 8, 7, MatIndexValGen>(10, 8, 7);

    
    MatrixMatrixProduct_test<DYNAMIC_SIZE, DYNAMIC_SIZE, DYNAMIC_SIZE, MatIndexValGen>(2, 3, 4);
    MatrixMatrixProduct_test<DYNAMIC_SIZE, DYNAMIC_SIZE, DYNAMIC_SIZE, MatIndexValGen>(4, 2, 3);
    MatrixMatrixProduct_test<DYNAMIC_SIZE, DYNAMIC_SIZE, DYNAMIC_SIZE, MatIndexValGen>(3, 4, 2);
    MatrixMatrixProduct_test<DYNAMIC_SIZE, DYNAMIC_SIZE, DYNAMIC_SIZE, MatIndexValGen>(10, 8, 7);

    
    MatrixMatrixProduct_test<2, DYNAMIC_SIZE, DYNAMIC_SIZE, MatIndexValGen>(2, 3, 4);
    MatrixMatrixProduct_test<4, DYNAMIC_SIZE, DYNAMIC_SIZE, MatIndexValGen>(4, 2, 3);
    MatrixMatrixProduct_test<3, DYNAMIC_SIZE, DYNAMIC_SIZE, MatIndexValGen>(3, 4, 2);
    MatrixMatrixProduct_test<10, DYNAMIC_SIZE, DYNAMIC_SIZE, MatIndexValGen>(10, 8, 7);

    
    MatrixMatrixProduct_test<DYNAMIC_SIZE, 3, DYNAMIC_SIZE, MatIndexValGen>(2, 3, 4);
    MatrixMatrixProduct_test<DYNAMIC_SIZE, 2, DYNAMIC_SIZE, MatIndexValGen>(4, 2, 3);
    MatrixMatrixProduct_test<DYNAMIC_SIZE, 4, DYNAMIC_SIZE, MatIndexValGen>(3, 4, 2);
    MatrixMatrixProduct_test<DYNAMIC_SIZE, 8, DYNAMIC_SIZE, MatIndexValGen>(10, 8, 7);

    
    MatrixMatrixProduct_test<DYNAMIC_SIZE, DYNAMIC_SIZE, 4, MatIndexValGen>(2, 3, 4);
    MatrixMatrixProduct_test<DYNAMIC_SIZE, DYNAMIC_SIZE, 3, MatIndexValGen>(4, 2, 3);
    MatrixMatrixProduct_test<DYNAMIC_SIZE, DYNAMIC_SIZE, 2, MatIndexValGen>(3, 4, 2);
    MatrixMatrixProduct_test<DYNAMIC_SIZE, DYNAMIC_SIZE, 7, MatIndexValGen>(10, 8, 7);

    
    MatrixMatrixProduct_test<2, 3, DYNAMIC_SIZE, MatIndexValGen>(2, 3, 4);
    MatrixMatrixProduct_test<4, 2, DYNAMIC_SIZE, MatIndexValGen>(4, 2, 3);
    MatrixMatrixProduct_test<3, 4, DYNAMIC_SIZE, MatIndexValGen>(3, 4, 2);
    MatrixMatrixProduct_test<10, 8, DYNAMIC_SIZE, MatIndexValGen>(10, 8, 7);

    
    MatrixMatrixProduct_test<DYNAMIC_SIZE, 3, 4, MatIndexValGen>(2, 3, 4);
    MatrixMatrixProduct_test<DYNAMIC_SIZE, 2, 3, MatIndexValGen>(4, 2, 3);
    MatrixMatrixProduct_test<DYNAMIC_SIZE, 4, 2, MatIndexValGen>(3, 4, 2);
    MatrixMatrixProduct_test<DYNAMIC_SIZE, 8, 7, MatIndexValGen>(10, 8, 7);

    
    MatrixMatrixProduct_test<2, DYNAMIC_SIZE, 4, MatIndexValGen>(2, 3, 4);
    MatrixMatrixProduct_test<4, DYNAMIC_SIZE, 3, MatIndexValGen>(4, 2, 3);
    MatrixMatrixProduct_test<3, DYNAMIC_SIZE, 2, MatIndexValGen>(3, 4, 2);
    MatrixMatrixProduct_test<10, DYNAMIC_SIZE, 7, MatIndexValGen>(10, 8, 7);

    srand(0);
    
    MatrixMatrixProduct_test<2, 3, 4, MatRandValGen>(2, 3, 4);
    MatrixMatrixProduct_test<4, 2, 3, MatRandValGen>(4, 2, 3);
    MatrixMatrixProduct_test<3, 4, 2, MatRandValGen>(3, 4, 2);
    MatrixMatrixProduct_test<10, 8, 7, MatRandValGen>(10, 8, 7);

    srand(0);
    
    MatrixMatrixProduct_test<DYNAMIC_SIZE, DYNAMIC_SIZE, DYNAMIC_SIZE, MatRandValGen>(2, 3, 4);
    MatrixMatrixProduct_test<DYNAMIC_SIZE, DYNAMIC_SIZE, DYNAMIC_SIZE, MatRandValGen>(4, 2, 3);
    MatrixMatrixProduct_test<DYNAMIC_SIZE, DYNAMIC_SIZE, DYNAMIC_SIZE, MatRandValGen>(3, 4, 2);
    MatrixMatrixProduct_test<DYNAMIC_SIZE, DYNAMIC_SIZE, DYNAMIC_SIZE, MatRandValGen>(10, 8, 7);

    srand(0);
    
    MatrixMatrixProduct_test<2, DYNAMIC_SIZE, DYNAMIC_SIZE, MatRandValGen>(2, 3, 4);
    MatrixMatrixProduct_test<4, DYNAMIC_SIZE, DYNAMIC_SIZE, MatRandValGen>(4, 2, 3);
    MatrixMatrixProduct_test<3, DYNAMIC_SIZE, DYNAMIC_SIZE, MatRandValGen>(3, 4, 2);
    MatrixMatrixProduct_test<10, DYNAMIC_SIZE, DYNAMIC_SIZE, MatRandValGen>(10, 8, 7);

    srand(0);
    
    MatrixMatrixProduct_test<DYNAMIC_SIZE, 3, DYNAMIC_SIZE, MatRandValGen>(2, 3, 4);
    MatrixMatrixProduct_test<DYNAMIC_SIZE, 2, DYNAMIC_SIZE, MatRandValGen>(4, 2, 3);
    MatrixMatrixProduct_test<DYNAMIC_SIZE, 4, DYNAMIC_SIZE, MatRandValGen>(3, 4, 2);
    MatrixMatrixProduct_test<DYNAMIC_SIZE, 8, DYNAMIC_SIZE, MatRandValGen>(10, 8, 7);

    srand(0);
    
    MatrixMatrixProduct_test<DYNAMIC_SIZE, DYNAMIC_SIZE, 4, MatRandValGen>(2, 3, 4);
    MatrixMatrixProduct_test<DYNAMIC_SIZE, DYNAMIC_SIZE, 3, MatRandValGen>(4, 2, 3);
    MatrixMatrixProduct_test<DYNAMIC_SIZE, DYNAMIC_SIZE, 2, MatRandValGen>(3, 4, 2);
    MatrixMatrixProduct_test<DYNAMIC_SIZE, DYNAMIC_SIZE, 7, MatRandValGen>(10, 8, 7);

    srand(0);
    
    MatrixMatrixProduct_test<2, 3, DYNAMIC_SIZE, MatRandValGen>(2, 3, 4);
    MatrixMatrixProduct_test<4, 2, DYNAMIC_SIZE, MatRandValGen>(4, 2, 3);
    MatrixMatrixProduct_test<3, 4, DYNAMIC_SIZE, MatRandValGen>(3, 4, 2);
    MatrixMatrixProduct_test<10, 8, DYNAMIC_SIZE, MatRandValGen>(10, 8, 7);

    srand(0);
    
    MatrixMatrixProduct_test<DYNAMIC_SIZE, 3, 4, MatRandValGen>(2, 3, 4);
    MatrixMatrixProduct_test<DYNAMIC_SIZE, 2, 3, MatRandValGen>(4, 2, 3);
    MatrixMatrixProduct_test<DYNAMIC_SIZE, 4, 2, MatRandValGen>(3, 4, 2);
    MatrixMatrixProduct_test<DYNAMIC_SIZE, 8, 7, MatRandValGen>(10, 8, 7);

    srand(0);
    
    MatrixMatrixProduct_test<2, DYNAMIC_SIZE, 4, MatRandValGen>(2, 3, 4);
    MatrixMatrixProduct_test<4, DYNAMIC_SIZE, 3, MatRandValGen>(4, 2, 3);
    MatrixMatrixProduct_test<3, DYNAMIC_SIZE, 2, MatRandValGen>(3, 4, 2);
    MatrixMatrixProduct_test<10, DYNAMIC_SIZE, 7, MatRandValGen>(10, 8, 7);
#endif
    
#ifdef TEST_GROUP2
    std::cerr << "ComplianceModelTest:\n";
    srand(0);
    ComplianceModelTest(11, 20);
    srand(0);
    ComplianceModelTest(117, 25);
    srand(0);
    ComplianceModelTest(239, 30);

    std::cerr << "ComplianceModelTest2:\n";
    srand(0);
    ComplianceModelTest2(5);
    ComplianceModelTest2(17);
    ComplianceModelTest2(128);
    ComplianceModelTest2(17384);

    std::cerr << "ComplianceModelTest3:\n";
    srand(0);
    ComplianceModelTest3(5);
    ComplianceModelTest3(17);
    ComplianceModelTest3(128);
    ComplianceModelTest3(17384);

    std::cerr << "ComplianceModelTest4:\n";
    srand(0);
    ComplianceModelTest4(5);
    ComplianceModelTest4(17);
    ComplianceModelTest4(128);
    ComplianceModelTest4(17384);

    std::cerr << "ComplianceModelTest5:\n";
    srand(0);
    ComplianceModelTest5(5, 2);
    ComplianceModelTest5(17, 5);
    ComplianceModelTest5(128, 6);
    ComplianceModelTest5(17384, 3);

    std::cerr << "ComplianceModelTest6:\n";
    srand(0);
    ComplianceModelTest6(5, 2);
    ComplianceModelTest6(17, 5);
    ComplianceModelTest6(128, 6);
    ComplianceModelTest6(17384, 3);

    std::cerr << "ComplianceModelTest7:\n";
    srand(0);
    ComplianceModelTest7(5, 2);
    ComplianceModelTest7(17, 5);
    ComplianceModelTest7(35, 5);
    ComplianceModelTest7(1200, 12);

    std::cerr << "ComplianceModelTest8:\n";
    srand(0);
    ComplianceModelTest8(5, 2);
    ComplianceModelTest8(17, 5);
    ComplianceModelTest8(35, 5);
    ComplianceModelTest8(1200, 12);

    std::cerr << "ComplianceModelTest9:\n";
    srand(0);
    ComplianceModelTest9(5, 2);
    ComplianceModelTest9(17, 5);
    ComplianceModelTest9(35, 5);
    ComplianceModelTest9(1200, 12);

    std::cerr << "RangeVectorTest1():\n";
    srand(0);
    RangeVectorTest1(1000);

    std::cerr << "RangeVectorTest2():\n";
    srand(0);
    RangeVectorTest2(100);
#endif
    
#ifdef TEST_GROUP1
    std::cerr << "TEST_GROUP1 executed" << std::endl;
#endif

#ifdef TEST_GROUP2
    std::cerr << "TEST_GROUP2 executed" << std::endl;
#endif
    
#ifdef NDEBUG
    std::cerr << "\nNo tests have been done" << std::endl;
#else
    std::cerr << "\nAll tests passed" << std::endl;
#endif

    std::cerr << "MATVEC_DEBUG=" << MATVEC_DEBUG << std::endl;
    std::cerr << "GRADIENT_DEBUG=" << GRADIENT_DEBUG << std::endl;

    return 0;
}

/*
%!test
%! printf("\n\ntestMatVecProduct()\n");
%!
%! A=[11 12 13 14
%! 21 22 23 24
%! 31 32 33 34]
%!
%! B= [0.011e3 0.012e3
%!    0.021e3 0.022e3
%!	  0.031e3 0.032e3
%!	  0.041e3 0.042e3];
%!
%! R=[11 12 13
%! 21 22 23
%! 31 32 33]
%!
%! x=[100
%! 200
%! 300
%! 400]
%!
%! g = [523;
%!      -786 ]
%!
%! b=A*x
%!
%! c = R * A * x
%!
%! C = A * B
%!
%! D = -A * A.'
%!
%! h = A * B * g
%!
%! d=A.'*b
%!
%! e=A.' * A * x
%! f=A*A.'*A*(x+x)

%!function print_matrix(x)
%!  printf("{");
%!  for i=1:rows(x)
%! 		if (columns(x) > 1)
%!			printf("{");
%!		endif
%!		for j=1:columns(x)
%!			printf("%.16g", x(i, j));
%!			if (j < columns(x))
%!				printf(", ");
%!			endif
%!		endfor
%!		if (columns(x) > 1)
%!			printf("}");
%!		endif
%!		if (i < rows(x))
%!			printf(",\n");
%!		endif
%!  endfor
%!  printf("}");

%!function print_gradient(g)
%!  printf("\nstatic const struct test_%s {\ndouble val[%d];\n", inputname(1), length(g.x));
%!  printf("doublereal der[%d][%d];\n}", rows(g.J), columns(g.J));
%!  printf(" oct_%s = {\n", inputname(1));
%!  print_matrix(g.x);
%!  printf(",\n");
%!  print_matrix(g.J);
%!  printf("};\n");
%!  % printf("\ntestGradient(oct_%s, %s, dTol);\n", inputname(1), inputname(1));
%!

%!test
%! A = [11, 12, 13, 14;
%!      21, 22, 23, 24;
%! 		31, 32, 33, 34]
%!
%! h = [15.7; -7.3; 12.4]
%!
%! A = gradinit(A);
%!
%! x = (A(1, :) * 10).'
%!
%! b = A * x
%! c = -A * x
%! d = cross(c, b) * 5
%! e = cross(c + d, b) / 2
%! f = cross(e, d + c) - e
%! g = cross(d + e, f - c) + b
%! i = cross(d, h) * 5.7
%! j = cross(h, g) / 3.5
%! k = cross(h, h) + h * 3
%! l = cross(h + h, j) / 2
%! m = cross(h, j - i) * 5
%! n = cross(h + h, h) + h
%! o = cross(h, h + h) * 2
%! p = cross(h * 2.5, h / 2) - h
%! r = dot(h, g)
%! s = dot(g, h)
%! s1 = dot(g.x, h)
%! t = dot(g, g)
%! u = dot(h, h);
%! norm_g = norm(g)
%! printf("\n\ntestMatVecProductGradient()\n");
%! print_gradient(b);
%! print_gradient(c);
%! print_gradient(d);
%! print_gradient(e);
%! print_gradient(f);
%! print_gradient(g);
%! print_gradient(i);
%! print_gradient(j);
%! print_gradient(l);
%! print_gradient(m);
%! print_gradient(r);
%! print_gradient(s);
%! print_gradient(t);
%! print_gradient(norm_g);
*/

