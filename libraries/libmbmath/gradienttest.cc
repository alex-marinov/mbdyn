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
 AUTHOR: Reinhard Resch <r.resch@a1.net>
        Copyright (C) 2013(-2017) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <functional>
#include <typeinfo>
#include <cmath>

#ifdef HAVE_FEENABLEEXCEPT
#define _GNU_SOURCE 1
#include <fenv.h>
#endif

#ifdef HAVE_BLITZ
#include <blitz/blitz.h>
#include <blitz/array.h>
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>
#include <blitz/matrix.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979
#endif

#ifndef GRADIENT_DEBUG
	#define GRADIENT_DEBUG 1
#endif

#include "clock_time.h"
#include "gradient.h"

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
	doublereal dTolAbs = std::max<T>(1., std::max<T>(std::abs(a), std::abs(b))) * dTolRel;

	if (!(std::abs(a - b) <= dTolAbs)) {
		std::cerr << "a = " << a << " != b = " << b << std::endl;
		return false;
	}

	return true;
}

template <index_type N_SIZE>
bool bCompare(const Gradient<N_SIZE>& a, const Gradient<N_SIZE>& b, doublereal dTolRel = 0.) {
	doublereal dTolAbs = std::max(1., std::max(std::abs(a.dGetValue()), std::abs(b.dGetValue()))) * dTolRel;

	if (std::abs(a.dGetValue() - b.dGetValue()) > dTolAbs) {
		std::cerr << "a = " << a.dGetValue() << " != b = " << b.dGetValue() << std::endl;
		return false;
	}

	index_type iStart = std::min(a.iGetStartIndexLocal(), b.iGetStartIndexLocal());
	index_type iEnd = std::max(a.iGetEndIndexLocal(), b.iGetEndIndexLocal());

	for (index_type i = iStart; i < iEnd; ++i) {
		doublereal dTolAbs = std::max(1., std::max(std::abs(a.dGetDerivativeLocal(i)), std::abs(b.dGetDerivativeLocal(i)))) * dTolRel;

		if (std::abs(a.dGetDerivativeLocal(i) - b.dGetDerivativeLocal(i)) > dTolAbs) {
			std::cerr << "ad(" << i << ") = " << a.dGetDerivativeLocal(i) << " != bd(" << i << ") = " << b.dGetDerivativeLocal(i) << std::endl;
			return false;
		}
	}

	return true;
}

template <int N_SIZE>
void testRangeVector() {
	std::cout << "testRangeVector<" << N_SIZE << ">();\n";
    assert(N_SIZE == 0 || N_SIZE >= 6);
    index_type start = 3;
    index_type end = 6;
    RangeVector<doublereal, N_SIZE> v(start, end, 0.);
    RangeVector<float, 2 * N_SIZE> vf(v), vf2;
    RangeVector<doublereal, N_SIZE> v2(v);

    int k = 0;
    
    for (index_type i = start; i < end; ++i) {
        v.SetValue(i, ++k);
    }
    
    v2 = v;
    vf = v;
    vf.Reset();
    vf.Reserve(N_SIZE == 0 ? vf.iGetEndIndex() * 2 : 2 * N_SIZE);
    vf = v2;
    const RangeVector<float, N_SIZE> v3(v2);
    v2.Reset();
    v2 = v3;
    v2.Reserve(N_SIZE == 0 ? v2.iGetEndIndex() * 2 : N_SIZE);

    std::cout << "v(" << v.iGetStartIndex() << ":" << v.iGetEndIndex() - 1 << ")\n";

    for (index_type i = v.iGetStartIndex(); i < v.iGetEndIndex(); ++i) {
        std::cout << "v[" << i << "]=" << v.GetValue(i) << std::endl;
    }

    for (index_type i = vf.iGetStartIndex(); i < vf.iGetEndIndex(); ++i) {
        std::cout << "vf[" << i << "]=" << vf.GetValue(i) << std::endl;
    }

    for (index_type i = v2.iGetStartIndex(); i < v2.iGetEndIndex(); ++i) {
        std::cout << "v2[" << i << "]=" << v2.GetValue(i) << std::endl;
    }

    for (index_type i = v3.iGetStartIndex(); i < v3.iGetEndIndex(); ++i) {
        std::cout << "v3[" << i << "]=" << v3.GetValue(i) << std::endl;
    }

    std::cout << std::endl;

    for (index_type i = 0; i < v.iGetEndIndex(); ++i)
    {
    	assert(v.GetValue(i) == vf.GetValue(i));
    	assert(v2.GetValue(i) == v.GetValue(i));
    }

    v2.ResizeReset(v.iGetStartIndex(), v.iGetEndIndex(), 0.);

    for (index_type i = v.iGetStartIndexVector(); i < v.iGetEndIndexVector(); ++i)
    {
    	v2.SetVectorValue(i, v.GetVectorValue(i));
    }

    for (index_type i = v.iGetStartIndex(); i < v.iGetEndIndex(); ++i)
    {
    	assert(v2.GetValue(i) == v.GetValue(i));
    }

    vf.ResizeReset(0, 12, 0.0f);

    for (index_type i = vf.iGetStartIndex(); i < vf.iGetEndIndex(); ++i)
    {
    	vf.SetValue(i, float(i + 1));
    }

    vf2.ResizeReset(vf.iGetStartIndex(), vf.iGetEndIndex(), 0.0f);

    for (index_type i = vf.iGetStartIndexVector(); i < vf.iGetEndIndexVector(); ++i)
    {
    	vf2.SetVectorValue(i, vf.GetVectorValue(i));
    }

    for (index_type i = vf.iGetStartIndex(); i < vf.iGetEndIndex(); ++i)
    {
    	assert(vf.GetValue(i) == vf2.GetValue(i));
    }

    for (index_type i = vf2.iGetStartIndex(); i < vf2.iGetEndIndex(); ++i) {
        std::cout << "vf2[" << i << "]=" << vf2.GetValue(i) << std::endl;
    }

    std::cout << std::endl;

    v.ResizePreserve(1, 6);

    std::cout << "v.ResizePreserve(1, 6);\n";
    std::cout << "v(" << v.iGetStartIndex() << ":" << v.iGetEndIndex() - 1 << ")\n";

    for (index_type i = v.iGetStartIndex(); i < v.iGetEndIndex(); ++i) {
        std::cout << "v[" << i << "]=" << v.GetValue(i) << std::endl;
    }

    std::cout << std::endl;

    v.ResizePreserve(1, 5);

    std::cout << "v.ResizePreserve(1, 5);\n";
    std::cout << "v(" << v.iGetStartIndex() << ":" << v.iGetEndIndex() - 1 << ")\n";

    for (index_type i = v.iGetStartIndex(); i < v.iGetEndIndex(); ++i) {
        std::cout << "v[" << i << "]=" << v.GetValue(i) << std::endl;
    }

    std::cout << std::endl;

    v.ResizePreserve(0, 6);
    
    std::cout << "v.ResizePreserve(0, 6);\n";
    std::cout << "v(" << v.iGetStartIndex() << ":" << v.iGetEndIndex() - 1 << ")\n";

    for (index_type i = v.iGetStartIndex(); i < v.iGetEndIndex(); ++i) {
        std::cout << "v[" << i << "]=" << v.GetValue(i) << std::endl;
    }

    std::cout << std::endl;

    v.ResizePreserve(4, 5);

    std::cout << " v.ResizePreserve(4, 5);\n";
    std::cout << "v(" << v.iGetStartIndex() << ":" << v.iGetEndIndex() - 1 << ")\n";

    for (index_type i = v.iGetStartIndex(); i < v.iGetEndIndex(); ++i) {
        std::cout << "v[" << i << "]=" << v.GetValue(i) << std::endl;
    }

    std::cout << std::endl;


    v.ResizePreserve(0, 6);

    std::cout << "v.ResizePreserve(0, 6);\n";
    std::cout << "v(" << v.iGetStartIndex() << ":" << v.iGetEndIndex() - 1 << ")\n";

    for (index_type i = v.iGetStartIndex(); i < v.iGetEndIndex(); ++i) {
        std::cout << "v[" << i << "]=" << v.GetValue(i) << std::endl;
    }

    std::cout << std::endl;
}

template <int N_SIZE>
void testMapVector() {
    assert(N_SIZE == 0 || N_SIZE == 19);
    LocalDofMap dof;
    
    const int iCount = 4;
    
    MapVector<N_SIZE> v[iCount] = {   MapVector<N_SIZE>(&dof, 501, 506+1, MapVector<N_SIZE>::GLOBAL, 1.),
                                      MapVector<N_SIZE>(&dof, 601, 606+1, MapVector<N_SIZE>::GLOBAL, 2.),
                                      MapVector<N_SIZE>(&dof, 701, 706+1, MapVector<N_SIZE>::GLOBAL, 3.),
                                      MapVector<N_SIZE>(&dof, 801, 801+1, MapVector<N_SIZE>::GLOBAL, 4.)   };

    MapVector<N_SIZE> v1(&dof, 0, 1+1, MapVector<N_SIZE>::LOCAL, 5.);
    MapVector<N_SIZE> v2(&dof, 1, 2+1, MapVector<N_SIZE>::LOCAL, 6.);                           

    for (index_type j = 0; j < dof.Size(); ++j) {
        for (int i = 0; i < iCount; ++i) {
            if (j >= v[i].iGetStartIndexLocal() && j < v[i].iGetEndIndexLocal()) {
                doublereal dVal = v[i].dGetLocalVector(j);

                std::cout << "v[" << i << "][" << j << "->" << dof.iGetGlobalDof(j) << "]=" << dVal << std::endl;
            }
        }
    }

    for (index_type j = 0; j < dof.Size(); ++j) {
        if (j >= v1.iGetStartIndexLocal() && j < v1.iGetEndIndexLocal()) {
            doublereal dVal = v1.dGetLocalVector(j);

            std::cout << "v1[" << j << "->" << dof.iGetGlobalDof(j) << "]=" << dVal << std::endl;
        }
    }

    for (index_type j = 0; j < dof.Size(); ++j) {
        if (j >= v2.iGetStartIndexLocal() && j < v2.iGetEndIndexLocal()) {
            doublereal dVal = v2.dGetLocalVector(j);

            std::cout << "v2[" << j << "->" << dof.iGetGlobalDof(j) << "]=" << dVal << std::endl;
        }
    }
}

template <typename T>
void func(const T& u, const T& v, const T& w, doublereal e, T& f) {
    f = ((((3 * u + 2 * v) * (u - v) / (1 - w) * pow(u/w, v - 1) * sin(v) * cos(w) * (1 - tan(-w + v))) * e + 1. - 11. + 4.5 - 1.) * 3.5 / 2.8 + u - v) * w / u;
}

template <typename T>
void func_tmp(const T& u, const T& v, const T& w, doublereal e, T& f) {
	f = u;
	f *= 3;
	f += 2 * v;
	f *= (u - v);
	f /= (1 - w);
	f *= pow(u/w, v - 1);
	f *= sin(v);
	f *= cos(w);
	f *= (1 - tan(-w + v));
	f *= e;

#if GRADIENT_DEBUG > 0
	const T tmp1 =
#endif
	f++;

	GRADIENT_ASSERT(tmp1 == f - 1);

	f -= 11.;
	f += 4.5;

#if GRADIENT_DEBUG > 0
	const T tmp2 =
#endif
	--f;

	GRADIENT_ASSERT(tmp2 == f);

#if GRADIENT_DEBUG > 0
	const T tmp3 = f;
	f += f;
	GRADIENT_ASSERT(bCompare(f, 2 * tmp3));
	f -= f;
	GRADIENT_ASSERT(bCompare(f, T()));
	f = tmp3;
#endif

	f *= 3.5;
	f /= 2.8;
	f += u;
	f -= v;
	f *= w;
	f /= u;
}

doublereal sec(doublereal x) {
    return 1./cos(x);
}

void funcDeriv(index_type N, doublereal u, const doublereal ud[], doublereal v,  const doublereal vd[], doublereal w, const doublereal wd[], doublereal e, doublereal& f, doublereal fd[]) {
	const doublereal tmp1 = sin(v);
	const doublereal tmp2 = cos(w);
	const doublereal tmp3 = u/w;
	const doublereal tmp6 = pow(tmp3, v-1);
	const doublereal tmp4 = tmp1*tmp6*tmp2;
	const doublereal tmp5 = (tan(w-v)+1);
	const doublereal tmp7 = ((2*v+3*u)*tmp4*tmp5);
	const doublereal tmp8 = pow(sec(w-v),2);
	const doublereal tmp9 = (u-v);
	const doublereal tmp10 = (1-w);
	const doublereal tmp11 = (v-1);
	const doublereal tmp12 = (2*v+3*u);
	const doublereal tmp13 = tmp9*tmp4*tmp5;
	const doublereal tmp14 = tmp12*tmp4*tmp5;
	const doublereal tmp15 = (tmp9*tmp12*tmp4*tmp8);
	const doublereal tmp16 = (tmp9*tmp11*tmp14);
	const doublereal tmp17 = tmp15/tmp10;
	const doublereal tmp18 = tmp6*tmp2*tmp5;
	const doublereal tmp19 = tmp9*tmp12*tmp1;
	const doublereal tmp20 = tmp7/tmp10;
	const doublereal tmp21 = 3.5 / 2.8;
	const doublereal tmp22 = e * tmp21;
	const doublereal tmp23 = ((((tmp9*tmp14)/tmp10) * e + 1. - 11. + 4.5 - 1.) * tmp21 + u - v);

    f = tmp23 * w / u;
    
    doublereal df_du = ((tmp16/(u*tmp10)+tmp20+(3*tmp13)/tmp10)*tmp22 + 1.) * w / u - tmp23 * w / (u * u);
    
    doublereal df_dv = (((tmp19*log(tmp3)*tmp18)/tmp10-tmp20+(2*tmp13)/tmp10+(tmp9*tmp12*cos(v)*tmp18)/tmp10-tmp17)*tmp22 - 1.) * w / u;
    
    doublereal df_dw = (-(tmp19*tmp6*sin(w)*tmp5)/tmp10-tmp16/(tmp10*w)+(tmp9*tmp14)/pow(1-w,2)+tmp17)*tmp22 * w / u + tmp23 / u;

    for (index_type i = 0; i < N; ++i) {
    	fd[i] = df_du * ud[i] + df_dv * vd[i] + df_dw * wd[i];
    }
}

template <typename T>
void func2(const T& a, const T& b, const T& c, const T& d, integer e, T& f)
{
    f = d - pow(a / b + c * d, e) / d + a * b;
}

double func2Deriv(double a, double b, double c, double d, integer e, double ad, double bd, double cd, double dd)
{
    return dd*((-(c*pow(c*d+a/b,e-1)*e)/d)+pow(c*d+a/b,e)/pow(d,2)+1)+ad*(b-(pow(c*d+a/b,e-1)*e)/(b*d))+bd*((a*pow(c*d+a/b,e-1)*e)/(pow(b,2)*d)+a)-cd*pow(c*d+a/b,e-1)*e;
}

template <index_type N_SIZE>
void test_func2(index_type N) {
    Gradient<N_SIZE> a, b, c, d, f;
    doublereal f_tmp;
    LocalDofMap dof;

    srand(0);

    a.SetValuePreserve(fabs(random1() * 3.5));
    b.SetValuePreserve(fabs(random1() * 1.5));
    c.SetValuePreserve(fabs(random1() * 5.5));
    d.SetValuePreserve(fabs(random1() * 3.5));

    a.DerivativeResizeReset(&dof, 0, N, MapVectorBase::GLOBAL, 0.);
    b.DerivativeResizeReset(&dof, 0, N, MapVectorBase::GLOBAL, 0.);
    c.DerivativeResizeReset(&dof, 0, N, MapVectorBase::GLOBAL, 0.);
    d.DerivativeResizeReset(&dof, 0, N, MapVectorBase::GLOBAL, 0.);
    
    for (index_type i = 0; i < N; ++i) {
    	a.SetDerivativeGlobal(i, random1() * 0.1);
        b.SetDerivativeGlobal(i, random1() * 3);
        c.SetDerivativeGlobal(i, random1() * 5);
        d.SetDerivativeGlobal(i, random1() * 9);
    }

    for (index_type e = -3; e <= 10; ++e) {
        func2(a, b, c, d, e, f);
        func2(a.dGetValue(),
              b.dGetValue(),
              c.dGetValue(),
              d.dGetValue(),
              e,
              f_tmp);

        assert(bCompare(f.dGetValue(), f_tmp, sqrt(std::numeric_limits<doublereal>::epsilon())));

        std::cout << "f=" << f.dGetValue() << std::endl;
        std::cout << "f_tmp=" << f_tmp << std::endl;
        
        for (index_type i = 0; i < N; ++i) {
            double fd_tmp = func2Deriv(a.dGetValue(),
                                       b.dGetValue(),
                                       c.dGetValue(),
                                       d.dGetValue(),
                                       e,
                                       a.dGetDerivativeGlobal(i),
                                       b.dGetDerivativeGlobal(i),
                                       c.dGetDerivativeGlobal(i),
                                       d.dGetDerivativeGlobal(i));
                                       
            assert(bCompare(f.dGetDerivativeGlobal(i), fd_tmp, sqrt(std::numeric_limits<doublereal>::epsilon())));

            std::cout << "df/dX" << i << "=" << f.dGetDerivativeGlobal(i) << std::endl;
            std::cout << "fd_tmp" << i << "=" << fd_tmp << std::endl;
        }
    }
}

template <int N_SIZE>
void testGradient(index_type N) {
    assert(N_SIZE == 0 || N_SIZE == N);

    srand(0);

    LocalDofMap dof;
    Gradient<N_SIZE> u, v, w, f, f_tmp;
    doublereal f_d;
    
    u.SetValuePreserve(fabs(random1()) * 10);
    v.SetValuePreserve(random1() * 3);
    w.SetValuePreserve(random1() * 0.01);
    
    u.DerivativeResizeReset(&dof, 0, N, MapVectorBase::GLOBAL, 0.);
    v.DerivativeResizeReset(&dof, 0, N, MapVectorBase::GLOBAL, 0.);
    w.DerivativeResizeReset(&dof, 0, N, MapVectorBase::GLOBAL, 0.);
    
    for (index_type i = 0; i < N; ++i) {
    	u.SetDerivativeGlobal(i, random1() * 0.01);
    	v.SetDerivativeGlobal(i, random1() * 3);
    	w.SetDerivativeGlobal(i, random1() * 10);
    }
    
    doublereal* ud_C = new doublereal[N];
    doublereal* vd_C = new doublereal[N];
    doublereal* wd_C = new doublereal[N];
    doublereal f_C;
    doublereal* fd_C = new doublereal[N];

    for (index_type i = 0; i < N; ++i) {
    	ud_C[i] = u.dGetDerivativeGlobal(i);
    	vd_C[i] = v.dGetDerivativeGlobal(i);
    	wd_C[i] = w.dGetDerivativeGlobal(i);
    }

    doublereal t_AD = 0., t_AD_tmp = 0., t_C = 0., t_d = 0.;

    for (int i = 0; i < NLoops; ++i) {
    	const doublereal e = random1();
    	doublereal start, stop;

    	tic(start);
    	func(u, v, w, e, f);
    	tic(stop);

    	t_AD += stop - start;

    	tic(start);
    	func(u.dGetValue(), v.dGetValue(), w.dGetValue(), e, f_d);
    	tic(stop);

    	t_d += stop - start;

    	tic(start);
    	func(u, v, w, e, f_tmp);
    	tic(stop);

    	t_AD_tmp += stop - start;

		tic(start);
		funcDeriv(N,
				  u.dGetValue(),
				  ud_C,
				  v.dGetValue(),
				  vd_C,
				  w.dGetValue(),
				  wd_C,
				  e,
				  f_C,
				  fd_C);
		tic(stop);

		t_C += stop - start;

		assert(bCompare(f, f_tmp, std::numeric_limits<doublereal>::epsilon()));
	    assert(bCompare(f.dGetValue(), f_C, sqrt(std::numeric_limits<doublereal>::epsilon())));
	    assert(bCompare(f.dGetValue(), f_d, std::numeric_limits<doublereal>::epsilon()));

	    for (index_type j = 0; j < N; ++j) {
	    	doublereal err = std::abs(f.dGetDerivativeGlobal(j) / fd_C[j] - 1.);
	    	assert(err < sqrt(std::numeric_limits<scalar_deriv_type>::epsilon()));
	    }
    }
    
    std::cerr << "testGradient<" << N_SIZE << ">(" << N << ")" << std::endl;

    std::cout << "dof map:" << std::endl << dof << std::endl;
    std::cout << "u=" << u << std::endl;
    std::cout << "v=" << v << std::endl;
    std::cout << "w=" << w << std::endl;
    std::cout << "f=" << f << std::endl;
    std::cout << "f_tmp=" << f_tmp << std::endl;

    std::cout << "f_C=" << f_C << std::endl;
    std::cout << "fd_C=" << std::endl;

    for (int i = 0; i < N; ++i) {
    	std::cout << " " << fd_C[i] << std::endl;
    }

    std::cout << std::endl;

    std::cerr << "t_AD=" << t_AD << "s" << std::endl;
    std::cerr << "t_d=" << t_d << "s" << std::endl;
    std::cerr << "t_AD_tmp=" << t_AD_tmp << "s" << std::endl;
    std::cerr << "t_C=" << t_C << "s" << std::endl;
    std::cerr << "overhead:" << t_AD / t_C << std::endl;
    std::cerr << "overhead tmp:" << t_AD_tmp / t_C << std::endl;

    delete [] ud_C;
    delete [] vd_C;
    delete [] wd_C;
    delete [] fd_C;
}

template <index_type N_SIZE>
void testGradient2() {
	assert(N_SIZE == 0 || N_SIZE == 3);

	const doublereal dTol = sqrt(std::numeric_limits<doublereal>::epsilon());

    srand(0);

    LocalDofMap dof;
    Gradient<N_SIZE> u, v, w, f;

    u.SetValuePreserve(fabs(random1()) * 10);
    v.SetValuePreserve(fabs(random1()) * 3);
    w.SetValuePreserve(fabs(random1()) * 0.01);

    u.DerivativeResizeReset(&dof, 0, MapVectorBase::GLOBAL, random1() * 0.01);
    v.DerivativeResizeReset(&dof, 1, MapVectorBase::GLOBAL, random1() * 3);
    w.DerivativeResizeReset(&dof, 2, MapVectorBase::GLOBAL, random1() * 10);

    doublereal ud_C[3];
    doublereal vd_C[3];
    doublereal wd_C[3];
    doublereal f_C;
    doublereal fd_C[3];

    for (index_type i = 0; i < 3; ++i) {
    	ud_C[i] = u.dGetDerivativeGlobal(i);
    	vd_C[i] = v.dGetDerivativeGlobal(i);
    	wd_C[i] = w.dGetDerivativeGlobal(i);
    }

    doublereal t_AD = 0., t_C = 0.;

    for (int i = 0; i < NLoops; ++i) {
    	const doublereal e = random1();
    	doublereal start, stop;

    	tic(start);
    	func(u, v, w, e, f);
    	tic(stop);

    	t_AD += stop - start;

		tic(start);
		funcDeriv(3,
				  u.dGetValue(),
				  ud_C,
				  v.dGetValue(),
				  vd_C,
				  w.dGetValue(),
				  wd_C,
				  e,
				  f_C,
				  fd_C);
		tic(stop);

		t_C += stop - start;

	    assert(bCompare(f.dGetValue(), f_C, dTol));
	    //assert(bCompare(Y[0], f_C, dTol));

	    for (index_type j = 0; j < 3; ++j) {
	    	assert(bCompare<scalar_deriv_type>(f.dGetDerivativeGlobal(j), fd_C[j], sqrt(std::numeric_limits<scalar_deriv_type>::epsilon())));
	    //	assert(bCompare(jac[j], fd_C[j], dTol));
	    }
    }

    std::cerr << "testGradient2<" << N_SIZE << ">(3)" << std::endl;

    std::cout << "dof map:" << std::endl << dof << std::endl;
    std::cout << "u=" << u << std::endl;
    std::cout << "v=" << v << std::endl;
    std::cout << "w=" << w << std::endl;
    std::cout << "f=" << f << std::endl;

    std::cout << "f_C=" << f_C << std::endl;
    std::cout << "fd_C=" << std::endl;

    for (int i = 0; i < 3; ++i) {
    	std::cout << " " << fd_C[i] << std::endl;
    }

    std::cout << std::endl;

    std::cerr << "t_AD=" << t_AD << "s" << std::endl;
    std::cerr << "t_C=" << t_C << "s" << std::endl;
    std::cerr << "overhead:" << t_AD / t_C << std::endl;

    Gradient<N_SIZE> w2 = 3 * u + 2 * v + 0.5 * w;

    w = 3 * u + 2 * v + 0.5 * Alias(w);

    GRADIENT_ASSERT(bCompare(w, w2));

    std::cout << "w=" << w << std::endl;
    std::cout << "w2=" << w2 << std::endl;
}
#if 0
void testGradient3(index_type N)
{
    srand(0);

    LocalDofMap dof;
    Gradient<N_SIZE> u, v, w, f;

    u.SetValuePreserve(random1() * 10);
    v.SetValuePreserve(random1() * 3);
    w.SetValuePreserve(random1() * 0.01);

    u.DerivativeResizeReset(&dof, N, MapVectorBase::GLOBAL, 0.);
    v.DerivativeResizeReset(&dof, N, MapVectorBase::GLOBAL, 0.);
    w.DerivativeResizeReset(&dof, N, MapVectorBase::GLOBAL, 0.);


    doublereal t_AD = 0., t_C = 0.;

    for (int i = 0; i < NLoops; ++i) {
    	const doublereal e = random1();
    	doublereal start, stop;

    	tic(start);
    	func3(u, v, w, e, f);
    	tic(stop);

    	t_AD += stop - start;
    }
}
#endif
void testGradientCopy(int N) {
	LocalDofMap map1, map2;
	bool bFirst = true;

	std::cerr << "testGradientCopy(" << N << ")\n";

	tic();

	for (int i = 0; i < NLoops; ++i) {
		MapVector<0> ud(&map1, 100, 102, MapVectorBase::GLOBAL, 1.);

		for (int i = 0; i < N; ++i) {
			MapVector<0> dummy(&map1, 102 + i, MapVectorBase::GLOBAL, 1.5);
		}

		MapVector<0> vd(&map1, 200, 203, MapVectorBase::GLOBAL, 2.);

		for (int i = 0; i < N; ++i) {
			MapVector<0> dummy(&map1, 203 + i, MapVectorBase::GLOBAL, 2.5);
		}

		MapVector<0> wd(&map1, 300, 302, MapVectorBase::GLOBAL, 3.);

		for (int i = 0; i < N; ++i) {
			MapVector<0> dummy(&map1, 302 + i, MapVectorBase::GLOBAL, 3.5);
		}

		Gradient<0> u(100., ud);
		Gradient<0> v(200., vd);
		Gradient<0> w(300., wd);

		assert(u != 100.001);
		assert(u >= 100.);
		assert(u <= 100.);
		assert(u == 100.);
		assert(u > 99.999);
		assert(u < 100.001);
		assert(100. == u);
		assert(100.001 != u);
		assert(100. <= u);
		assert(100. >= u);
		assert(99.999 < u);
		assert(100.001 > u);

		assert(v != 200.001);
		assert(v >= 200.);
		assert(v <= 200.);
		assert(v == 200.);
		assert(v > 199.999);
		assert(v < 200.001);
		assert(200. == v);
		assert(200.001 != v);
		assert(200. <= v);
		assert(200. >= v);
		assert(199.999 < v);
		assert(200.001 > v);

		Gradient<0> x = v + w;

		if (bFirst) {
			std::cout << std::endl;
			std::cout << "map1=" << std::endl << map1 << std::endl;
			std::cout << "u=" << u << std::endl;
			std::cout << "v=" << v << std::endl;
			std::cout << "w=" << w << std::endl;
			std::cout << "x=" << x << std::endl;
			std::cout << std::endl;
		}

		Gradient<5> y(x, &map2);

		if (bFirst) {
			std::cout << "map2=" << std::endl << map2 << std::endl;
			std::cout << "y=" << y << std::endl;
			std::cout << std::endl;
		}

		bFirst = false;
	}

	std::cerr << "t=" << toc() << "s" << std::endl;
}

template <index_type N_SIZE>
void testDifferentDofMaps(index_type N) {
	{
		LocalDofMap dofMap1, dofMap2;

		Gradient<N_SIZE> g1(1000., MapVector<N_SIZE>(&dofMap1, 100, 100 + N, MapVectorBase::GLOBAL, 1.));
		Gradient<N_SIZE> g2(2000., MapVector<N_SIZE>(&dofMap2, 102, 102 + N, MapVectorBase::GLOBAL, 2.));

		std::cerr << "testDifferentDofMaps<" << N_SIZE << ">(" << N << "): operator+=()\n";
		std::cerr << "g1=" << g1 << std::endl;
		std::cerr << "g2=" << g2 << std::endl;

		g1 += g2;

		std::cerr << "g1=" << g1 << std::endl;

		assert(g1.dGetValue() == 3000.);

		for (index_type i = 100; i < 102 + N; ++i) {
			if (i < 100 + N || i >= 102) {
				const doublereal d = g1.dGetDerivativeGlobal(i);
				if (i >= 102 && i < 100 + N) {
					assert(d == 3.);
				} else if (i < 100 + N) {
					assert(d == 1.);
				} else {
					assert(d == 2.);
				}
			}
		}
	}
	{
		LocalDofMap dofMap1, dofMap2;

		const doublereal u = 3, v = 4, du = 5., dv = 6.;
		Gradient<N_SIZE> g1(u, MapVector<N_SIZE>(&dofMap1, 100, 100 + N, MapVectorBase::GLOBAL, du));
		Gradient<N_SIZE> g2(v, MapVector<N_SIZE>(&dofMap2, 102, 102 + N, MapVectorBase::GLOBAL, dv));

		std::cerr << "testDifferentDofMaps<" << N_SIZE << ">(" << N << "): operator*=()\n";
		std::cerr << "g1=" << g1 << std::endl;
		std::cerr << "g2=" << g2 << std::endl;

		g1 *= g2;

		std::cerr << "g1=" << g1 << std::endl;

		assert(g1.dGetValue() == u * v);

		for (index_type i = 100; i < 102 + N; ++i) {
			if (i < 100 + N || i >= 102) {
				const doublereal d = g1.dGetDerivativeGlobal(i);
				if (i >= 102 && i < 100 + N) {
					assert(d == du * v + u * dv);
				} else if (i < 100 + N) {
					assert(d == du * v);
				} else {
					assert(d == u * dv);
				}
			}
		}
	}
}

template <index_type N_SIZE>
void testGradientLin(int N) {
	Gradient<N_SIZE> A11, A12, A13, A14;
	Gradient<N_SIZE> x1, x2, x3, x4;
	LocalDofMap oDof;

	A11.SetValuePreserve(1);
	A12.SetValuePreserve(2);
	A13.SetValuePreserve(3);
	A14.SetValuePreserve(4);

	A11.DerivativeResizeReset(&oDof, 0, MapVectorBase::GLOBAL, 1);
	A12.DerivativeResizeReset(&oDof, 1, MapVectorBase::GLOBAL, 1);
	A13.DerivativeResizeReset(&oDof, 2, MapVectorBase::GLOBAL, 1);
	A14.DerivativeResizeReset(&oDof, 3, MapVectorBase::GLOBAL, 1);

	std::cout << "A11=" << A11 << std::endl;
	std::cout << "A12=" << A12 << std::endl;
	std::cout << "A13=" << A13 << std::endl;
	std::cout << "A14=" << A14 << std::endl;

	x1.SetValuePreserve(1);
	x2.SetValuePreserve(1);
	x3.SetValuePreserve(1);
	x4.SetValuePreserve(1);

	x1.DerivativeResizeReset(&oDof, 0, MapVectorBase::GLOBAL, 1);
	x2.DerivativeResizeReset(&oDof, 1, MapVectorBase::GLOBAL, 1);
	x3.DerivativeResizeReset(&oDof, 2, MapVectorBase::GLOBAL, 1);
	x4.DerivativeResizeReset(&oDof, 3, MapVectorBase::GLOBAL, 1);

	std::cout << "x1=" << x1 << std::endl;
	std::cout << "x2=" << x2 << std::endl;
	std::cout << "x3=" << x3 << std::endl;
	std::cout << "x4=" << x4 << std::endl;

	Gradient<N_SIZE> y1 = A11 * x1;
	Gradient<N_SIZE> y2 = A12 * x2;
	Gradient<N_SIZE> y3 = A13 * x3;
	Gradient<N_SIZE> y4 = A14 * x4;

	std::cout << "y1=" << y1 << std::endl;
	std::cout << "y2=" << y2 << std::endl;
	std::cout << "y3=" << y3 << std::endl;
	std::cout << "y4=" << y4 << std::endl;

	Gradient<N_SIZE> y_2 = y1 + y2;
	Gradient<N_SIZE> y_3 = y_2 + y3;
	Gradient<N_SIZE> y_4 = y_3 + y4;

	Gradient<N_SIZE> y_1 = y1 + y2 + y3 + y4;
	Gradient<N_SIZE> y = A11 * x1 + A12 * x2 + A13 * x3 + A14 * x4;

	std::cout << "y=" << y << std::endl;
	std::cout << "y_1=" << y_1 << std::endl;
	std::cout << "y_2=" << y_2 << std::endl;
	std::cout << "y_3=" << y_3 << std::endl;
	std::cout << "y_4=" << y_4 << std::endl;

	assert(y.bIsEqual(y_1));
	assert(y.bIsEqual(y_4));
	assert(y.dGetValue() == 10);
	assert(y.dGetDerivativeGlobal(0) == 2);
	assert(y.dGetDerivativeGlobal(1) == 3);
	assert(y.dGetDerivativeGlobal(2) == 4);
	assert(y.dGetDerivativeGlobal(3) == 5);
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

#ifdef HAVE_BLITZ
namespace gradAssTest {

using namespace blitz;

const int I1 = 0, I2 = 1, I3 = 2;

template <typename T>
void Euler123ToMatR(const TinyVector<T, 3>& theta, TinyMatrix<T, 3, 3>& R1) {
    const T cos_theta2 = cos(theta(I2));
    const T cos_theta3 = cos(theta(I3));
    const T sin_theta2 = sin(theta(I2));
    const T sin_theta3 = sin(theta(I3));
    const T sin_theta1 = sin(theta(I1));
    const T cos_theta1 = cos(theta(I1));
    const T sin_theta1_sin_theta2 = sin_theta1 * sin_theta2;
    const T cos_theta1_sin_theta3 = cos_theta1 * sin_theta3;
    const T cos_theta1_cos_theta3 = cos_theta1 * cos_theta3;
    R1(I1, I1) = cos_theta2 * cos_theta3;
    R1(I1, I2) = cos_theta2 * sin_theta3;
    R1(I1, I3) = -sin_theta2;
    R1(I2, I1) = sin_theta1_sin_theta2 * cos_theta3 - cos_theta1_sin_theta3;
    R1(I2, I2) = sin_theta1_sin_theta2 * sin_theta3 + cos_theta1_cos_theta3;
    R1(I2, I3) = sin_theta1 * cos_theta2;
    R1(I3, I1) = sin_theta1 * sin_theta3 + cos_theta1_cos_theta3 * sin_theta2;
    R1(I3, I2) = cos_theta1_sin_theta3 * sin_theta2 - sin_theta1 * cos_theta3;
    R1(I3, I3) = cos_theta1 * cos_theta2;
}

template <typename T, int N_rows, int N_cols>
TinyMatrix<T, N_cols, N_rows>& transpose(const TinyMatrix<T, N_rows, N_cols>& A, TinyMatrix<T, N_cols, N_rows>& A_T) {
    for (int i = 0; i < N_rows; ++i) {
        for (int j = 0; j < N_cols; ++j) {
            A_T(j, i) = A(i, j);
        }
    }

    return A_T;
}

template <typename T, int N_rows, int N_cols>
TinyVector<T, N_cols>& MatTDotVec(const TinyMatrix<T, N_rows, N_cols>& A, const TinyVector<T, N_rows>& x, TinyVector<T, N_cols>& v) {
    for (int i = 0; i < N_cols; ++i) {
        v(i) = T(0.);

        for (int j = 0; j < N_rows; ++j) {
            v(i) += A(j, i) * x(j);
        }
    }

    return v;
}

template <typename T, int R1, int C1, int C2>
TinyMatrix<T, R1, C2> product(const TinyMatrix<T, R1, C1>& A, const TinyMatrix<T, C1, C2>& B) {
    TinyMatrix<T, R1, C2> C;
    C.initialize(T(0.));

    for (int i = 0; i < R1; ++i) {
        for (int j = 0; j < C2; ++j) {
            for (int k = 0; k < C1; ++k) {
                C(i, j) += A(i, k) * B(k, j);
            }
        }
    }

    return C;
}

template <typename T>
TinyMatrix<T, 3, 3>& skew(TinyMatrix<T, 3, 3>& skew_g, const TinyVector<T, 3>& g) {
    skew_g(0, 0) = T(0.);
    skew_g(1, 1) = T(0.);
    skew_g(2, 2) = T(0.);

    skew_g(0, 1) = -g(2);
    skew_g(1, 0) = g(2);

    skew_g(0, 2) = g(1);
    skew_g(2, 0) = -g(1);

    skew_g(1, 2) = -g(0);
    skew_g(2, 1) = g(0);

    return skew_g;
}

template <typename T>
TinyMatrix<T, 3, 3> skew(const TinyVector<T, 3>& g) {
    TinyMatrix<T, 3, 3> skew_g;
    /*
    skew(g) = [ 0, -z, y;
                z, 0, -x;
                -y, x, 0 ];
    */

    skew(skew_g, g);

    return skew_g;
}

template <typename T>
TinyMatrix<T, 3, 3>& skew_skew(TinyMatrix<T, 3, 3>& skew_g_skew_g, const TinyVector<T, 3>& g) {
    skew_g_skew_g(0, 0) = -g(2)*g(2)-g(1)*g(1);
    skew_g_skew_g(0, 1) = g(0)*g(1);
    skew_g_skew_g(0, 2) = g(0)*g(2);
    skew_g_skew_g(1, 0) = g(0)*g(1);
    skew_g_skew_g(1, 1) = -g(2)*g(2)-g(0)*g(0);
    skew_g_skew_g(1, 2) = g(1)*g(2);
    skew_g_skew_g(2, 0) = g(0)*g(2);
    skew_g_skew_g(2, 1) = g(1)*g(2);
    skew_g_skew_g(2, 2) = -g(1)*g(1)-g(0)*g(0);

    return skew_g_skew_g;
}

template <typename T>
TinyMatrix<T, 3, 3> transpose(const TinyMatrix<T, 3, 3>& A) {
    TinyMatrix<T, 3, 3> At;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            At(i, j) = A(j, i);
        }
    }

    return At;
}

template <typename T, int N_rows, int N_cols>
TinyMatrix<T, N_rows, N_cols> operator-(const TinyMatrix<T, N_rows, N_cols>& A) {
    TinyMatrix<T, N_rows, N_cols> B(A);

    for (int i = 0; i < N_rows; ++i) {
        for (int j = 0; j < N_cols; ++j) {
            B(i, j) *= -1.;
        }
    }

    return B;
}

template <typename T, int N_rows, int N_cols>
TinyMatrix<T, N_rows, N_cols> operator-(const TinyMatrix<T, N_rows, N_cols>& A, const TinyMatrix<T, N_rows, N_cols>& B) {
    TinyMatrix<T, N_rows, N_cols> C(A);

    for (int i = 0; i < N_rows; ++i) {
        for (int j = 0; j < N_cols; ++j) {
            C(i, j) -= B(i, j);
        }
    }

    return C;
}

template <typename T, int N_rows, int N_cols>
TinyMatrix<T, N_rows, N_cols> operator+(const TinyMatrix<T, N_rows, N_cols>& A, const TinyMatrix<T, N_rows, N_cols>& B) {
    TinyMatrix<T, N_rows, N_cols> C(A);

    for (int i = 0; i < N_rows; ++i) {
        for (int j = 0; j < N_cols; ++j) {
            C(i, j) += B(i, j);
        }
    }

    return C;
}

class Node {
public:
    Node(const TinyVector<doublereal, 3>& X_0,
         const TinyVector<doublereal, 3>& XP_0,
         const TinyMatrix<doublereal, 3, 3>& R_0,
         const TinyVector<doublereal, 3>& W_0)
        :iFirstDofIndex(-1), R0(R_0), W0(W_0) {

        for (int i = 0; i < 3; ++i) {
            X(i).SetValuePreserve(X_0(i));
            X(i).DerivativeResizeReset(&dof, i, MapVectorBase::LOCAL, 0.);

            XP(i).SetValuePreserve(XP_0(i));
            XP(i).DerivativeResizeReset(&dof, i, MapVectorBase::LOCAL, 0.);
            XP(i).SetDerivativeLocal(i, -1.); // derivative will be always -1

            g(i).SetValuePreserve(0.);
            g(i).DerivativeResizeReset(&dof, i, MapVectorBase::LOCAL, 0.);

            gP(i).SetValuePreserve(0.);
            gP(i).DerivativeResizeReset(&dof, i, MapVectorBase::LOCAL, 0.);
            gP(i).SetDerivativeLocal(i, -1.); // derivative will be always -1
        }
    }

    void SetValue(Array<doublereal, 1>& XCurr, Array<doublereal, 1>& XPrimeCurr) {
        assert(iFirstDofIndex != -1);

        for (int i = 0; i < 3; ++i) {
            XCurr(iFirstDofIndex + i) = X(i).dGetValue();
            XPrimeCurr(iFirstDofIndex + i) = XP(i).dGetValue();
            XCurr(iFirstDofIndex + i + 3) = g(i).dGetValue();
            XPrimeCurr(iFirstDofIndex + i + 3) = gP(i).dGetValue();
        }
    }

    void Update(const Array<doublereal, 1>& XCurr, const Array<doublereal, 1>& XPrimeCurr, doublereal dCoef) {
        assert(iFirstDofIndex != -1);

        for (int i = 0; i < 3; ++i) {
            X(i).SetValuePreserve(XCurr(iFirstDofIndex + i));
            X(i).SetDerivativeLocal(i, -dCoef);
            XP(i).SetValuePreserve(XPrimeCurr(iFirstDofIndex + i));
            g(i).SetValuePreserve(XCurr(iFirstDofIndex + i + 3));
            g(i).SetDerivativeLocal(i, -dCoef);
            gP(i).SetValuePreserve(XPrimeCurr(iFirstDofIndex + i + 3));
        }

        UpdateRotation();
    }

    void UpdateRotation() {
        TinyMatrix<Gradient<NADVars>, 3, 3> RDelta, G;
        TinyMatrix<Gradient<NADVars>, 3, 3> skew_g;
        TinyMatrix<Gradient<NADVars>, 3, 3> skew_skew_g;

        skew(skew_g, g);
        skew_skew(skew_skew_g, g);

        const Gradient<NADVars> d = 4. / (4. + dot(g, g));

        for (int i = 0; i < 3; ++i) {
            RDelta(i, i).SetValuePreserve(1.);
            G(i, i) = d;
        }

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                RDelta(i, j) += d * (skew_g(i, j) + 0.5 * skew_skew_g(i, j));
                G(i, j) += 0.5 * d * skew_g(i, j);
            }
        }

        R = product(RDelta, R0);
        W = product(G, gP) + product(RDelta, W0);
    }

    void AfterConvergence(const Array<doublereal, 1>& XCurr, const Array<doublereal, 1>& XPrimeCurr) {
        for (int i = 0; i < 3; ++i) {
            W0(i) = W(i).dGetValue();
            g(i).SetValuePreserve(0.);
            gP(i).SetValuePreserve(0.);

            for (int j = 0; j < 3; ++j) {
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

    void GetXCurr(TinyVector<doublereal, 3>& XCurr, LocalDofMap*) const {
        for (int i = 0; i < 3; ++i) {
            XCurr(i) = X(i).dGetValue();
        }
    }

    template <index_type N_SIZE>
    void GetXCurr(TinyVector<Gradient<N_SIZE>, 3>& XCurr, LocalDofMap* pDofMap) const {
        assert(iFirstDofIndex != -1);
        assert(pDofMap != 0);

        for (int i = 0; i < 3; ++i) {
            XCurr(i).SetValuePreserve(X(i).dGetValue());
            XCurr(i).DerivativeResizeReset(pDofMap, iFirstDofIndex + X(i).iGetStartIndexLocal(), iFirstDofIndex + X(i).iGetEndIndexLocal(), MapVectorBase::GLOBAL, 0.);

            for (index_type j = X(i).iGetStartIndexLocal(); j < X(i).iGetEndIndexLocal(); ++j) {
            	XCurr(i).SetDerivativeGlobal(iFirstDofIndex + j, X(i).dGetDerivativeLocal(j));
            }
        }
    }

    void GetVCurr(TinyVector<doublereal, 3>& VCurr, LocalDofMap*) const {
        for (int i = 0; i < 3; ++i) {
            VCurr(i) = XP(i).dGetValue();
        }
    }

    template <index_type N_SIZE>
    void GetVCurr(TinyVector<Gradient<N_SIZE>, 3>& VCurr, LocalDofMap* pDofMap) const {
        assert(iFirstDofIndex != -1);
        assert(pDofMap != 0);

        for (int i = 0; i < 3; ++i) {
            VCurr(i).SetValuePreserve(XP(i).dGetValue());
            VCurr(i).DerivativeResizeReset(pDofMap, iFirstDofIndex + XP(i).iGetStartIndexLocal(), iFirstDofIndex + XP(i).iGetEndIndexLocal(), MapVectorBase::GLOBAL, 0.);

            for (index_type j = XP(i).iGetStartIndexLocal(); j < XP(i).iGetEndIndexLocal(); ++j) {
            	VCurr(i).SetDerivativeGlobal(iFirstDofIndex + j, XP(i).dGetDerivativeLocal(j));
            }
        }
    }

    void GetRCurr(TinyMatrix<doublereal, 3, 3>& RCurr, LocalDofMap*) const {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                RCurr(i, j) = R(i, j).dGetValue();
            }
        }
    }

    template <index_type N_SIZE>
    void GetRCurr(TinyMatrix<Gradient<N_SIZE>, 3, 3>& RCurr, LocalDofMap* pDofMap) const {
        assert(iFirstDofIndex != -1);
        assert(pDofMap != 0);

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                RCurr(i, j).SetValuePreserve(R(i, j).dGetValue());
                RCurr(i, j).DerivativeResizeReset(pDofMap, iFirstDofIndex + R(i, j).iGetStartIndexLocal() + 3, iFirstDofIndex + R(i, j).iGetEndIndexLocal() + 3, MapVectorBase::GLOBAL, 0.);

                for (index_type k = R(i, j).iGetStartIndexLocal(); k < R(i, j).iGetEndIndexLocal(); ++k) {
                	RCurr(i, j).SetDerivativeGlobal(iFirstDofIndex + k + 3, R(i, j).dGetDerivativeLocal(k));
                }
            }
        }
    }

    const TinyMatrix<doublereal, 3, 3>& GetRRef() const {
        return R0;
    }

    void GetWCurr(TinyVector<doublereal, 3>& WCurr, LocalDofMap*) const {
        for (int i = 0; i < 3; ++i) {
            WCurr(i) = W(i).dGetValue();
        }
    }

    template <index_type N_SIZE>
    void GetWCurr(TinyVector<Gradient<N_SIZE>, 3>& WCurr, LocalDofMap* pDofMap) const {
        assert(iFirstDofIndex != -1);
        assert(pDofMap != 0);

        for (int i = 0; i < 3; ++i) {
            WCurr(i).SetValuePreserve(W(i).dGetValue());
            WCurr(i).DerivativeResizeReset(pDofMap, iFirstDofIndex + W(i).iGetStartIndexLocal() + 3, iFirstDofIndex + W(i).iGetEndIndexLocal() + 3, MapVectorBase::GLOBAL, 0.);

            for (index_type j = W(i).iGetStartIndexLocal(); j < W(i).iGetEndIndexLocal(); ++j) {
            	WCurr(i).SetDerivativeGlobal(iFirstDofIndex + j + 3, W(i).dGetDerivativeLocal(j));
            }
        }
    }

    const TinyVector<doublereal, 3>& GetWRef() const {
        return W0;
    }

private:
    int iFirstDofIndex;
    TinyMatrix<doublereal, 3, 3> R0;
    TinyVector<doublereal, 3> W0;
    static const int NADVars = 3;
    TinyVector<Gradient<NADVars>, 3> X, XP, g, gP, W;
    TinyMatrix<Gradient<NADVars>, 3, 3> R;
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

    void AddTo(Array<doublereal, 2>& JacMat) const {
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
    Array<doublereal, 2> oWorkMat;
    Array<int, 1> oRowIndex, oColIndex;
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
    SparseSubMatrixHandler& AssJac(T* pElem, LocalDofMap* pDofMap, Array<ResItem<Gradient<N_SIZE> >, 1>& WorkVec, doublereal dCoef) {
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

    void AddTo(Array<doublereal, 2>& JacMat) const {
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
    virtual ~Element() { }
};

class Element1: public Element {
private:
    Node* node1;
    Node* node2;
    TinyVector<doublereal, 3> o1, o2;
    TinyMatrix<doublereal, 3, 3> S, D;
    static const int NADVars = 12;
    LocalDofMap dofMap;

public:
    Element1(Node* node1_,
             const TinyVector<doublereal, 3>& o1_,
             Node* node2_,
             const TinyVector<doublereal, 3>& o2_,
             doublereal s,
             doublereal d)
        :node1(node1_),
         node2(node2_),
         o1(o1_),
         o2(o2_),
         dofMap(iGetNumCols()) {

        S.initialize(0.);
        D.initialize(0.);

        /*
            S=[ s,  -s,    0;
               -s, 2*s,   -s;
                0,  -s,  2*s];
        */

        S(0, 0) = s;
        S(1, 0) = -s;
        S(0, 1) = -s;
        S(1, 1) = 2*s;
        S(2, 1) = -s;
        S(1, 2) = -s;
        S(2, 2) = 2*s;

        /*
            D=[ d,     -d,      0;
               -d,  2 * d,     -d;
                0,     -d,  2 * d];
        */

        D(0, 0) = d;
        D(1, 0) = -d;
        D(0, 1) = -d;
        D(1, 1) = 2*d;
        D(2, 1) = -d;
        D(1, 2) = -d;
        D(2, 2) = 2*d;
    }

    template <typename T>
    Array<ResItem<T>, 1>& AssRes(Array<ResItem<T>, 1>& WorkVec, doublereal dCoef, LocalDofMap *pDofMap) {
        WorkVec.resize(iGetNumRows());

        TinyVector<T, 3> X1, X2, V1, V2, W1, W2;
        TinyMatrix<T, 3, 3> R1, R2;

        node1->GetXCurr(X1, pDofMap);
        node1->GetVCurr(V1, pDofMap);
        node1->GetRCurr(R1, pDofMap);
        node1->GetWCurr(W1, pDofMap);
        node2->GetXCurr(X2, pDofMap);
        node2->GetVCurr(V2, pDofMap);
        node2->GetRCurr(R2, pDofMap);
        node2->GetWCurr(W2, pDofMap);

        const TinyVector<T, 3> R1o1 = product(R1, o1);
        const TinyVector<T, 3> R2o2 = product(R2, o2);
        const TinyMatrix<T, 3, 3> R1_T = transpose(R1);
        const TinyVector<T, 3> dX = product(R1_T, TinyVector<T, 3>(X1 + R1o1 - X2 - R2o2));
        const TinyVector<T, 3> dV = product(R1_T, TinyVector<T, 3>(V1 + cross(W1, R1o1) - V2 - cross(W2, R2o2)));
        const TinyVector<T, 3> F1 = product(R1, TinyVector<T, 3>(product(-S, dX) - product(D, dV)));
        const TinyVector<T, 3> M1 = cross(R1o1, F1), M2 = -cross(R2o2, F1);

        for (int i = 0; i < 6; ++i) {
            WorkVec(i).iEquIndex = node1->iGetFirstIndex() + i;
            WorkVec(i + 6).iEquIndex = node2->iGetFirstIndex() + i;
        }

        for (int i = 0; i < 3; ++i) {
            WorkVec(i).dCoef = F1(i);
            WorkVec(i + 3).dCoef = M1(i);
            WorkVec(i + 6).dCoef = -F1(i);
            WorkVec(i + 9).dCoef = M2(i);
        }

        return WorkVec;
    }

    virtual Array<ResItem<doublereal>, 1>& AssRes(Array<ResItem<doublereal>, 1>& WorkVec, doublereal dCoef) {
        return AssRes(WorkVec, dCoef, 0);
    }

    virtual SparseSubMatrixHandler& AssJac(SparseSubMatrixHandler& WorkMat, doublereal dCoef) {
        Array<ResItem<Gradient<NADVars> >, 1> WorkVec;
        return WorkMat.AssJac(this, &dofMap, WorkVec, dCoef);
    }

    virtual FullSubMatrixHandler& AssJac(FullSubMatrixHandler& WorkMat, doublereal dCoef) {
        WorkMat.ResizeReset(iGetNumRows(), iGetNumCols());

        for (int i = 0; i < 6; ++i) {
            WorkMat.PutColIndex(i, node1->iGetFirstIndex() + i);
            WorkMat.PutColIndex(i + 6, node2->iGetFirstIndex() + i);
            WorkMat.PutRowIndex(i, node1->iGetFirstIndex() + i);
            WorkMat.PutRowIndex(i + 6, node2->iGetFirstIndex() + i);
        }

        const TinyVector<doublereal, 3>& W1_0 = node1->GetWRef();
        const TinyVector<doublereal, 3>& W2_0 = node2->GetWRef();
        const TinyMatrix<doublereal, 3, 3>& R1_0 = node1->GetRRef();
        const TinyMatrix<doublereal, 3, 3>& R2_0 = node2->GetRRef();

        TinyVector<doublereal, 3> X1, X2, V1, V2, W1, W2;
        TinyMatrix<doublereal, 3, 3> R1, R2;

        node1->GetXCurr(X1, 0);
        node1->GetVCurr(V1, 0);
        node1->GetRCurr(R1, 0);
        node1->GetWCurr(W1, 0);

        node2->GetXCurr(X2, 0);
        node2->GetVCurr(V2, 0);
        node2->GetRCurr(R2, 0);
        node2->GetWCurr(W2, 0);

        const TinyMatrix<doublereal, 3, 3> skew_W2_0 = skew(W2_0);
        const TinyVector<doublereal, 3> R1o1 = product(R1, o1);
        const TinyMatrix<doublereal, 3, 3> skew_R1o1 = skew(R1o1);
        const TinyVector<doublereal, 3> R1_0o1 = product(R1_0, o1);
        const TinyMatrix<doublereal, 3, 3> skew_R1_0o1 = skew(R1_0o1);
        const TinyVector<doublereal, 3> R2o2 = product(R2, o2);
        const TinyMatrix<doublereal, 3, 3> skew_R2o2 = skew(R2o2);
        const TinyVector<doublereal, 3> R2_0o2 = product(R2_0, o2);
        const TinyMatrix<doublereal, 3, 3> skew_R2_0o2 = skew(R2_0o2);
        const TinyMatrix<doublereal, 3, 3> R1_T = transpose(R1);
        const TinyVector<doublereal, 3> dX = product(R1_T, TinyVector<doublereal, 3>(X1 + R1o1 - X2 - R2o2));
        const TinyVector<doublereal, 3> dV = product(R1_T, TinyVector<doublereal, 3>(V1 + cross(W1, R1o1) - V2 - cross(W2, R2o2)));
        const TinyVector<doublereal, 3> F1_R1 = -(product(S, dX) + product(D, dV));
        const TinyVector<doublereal, 3> F1 = product(R1, F1_R1);
        const TinyVector<doublereal, 3> F2 = -F1;
        const TinyVector<doublereal, 3> M1 = cross(R1o1, F1), M2 = -cross(R2o2, F1);
        const TinyMatrix<doublereal, 3, 3> R1_0_T = transpose(R1_0);
        const TinyMatrix<doublereal, 3, 3> dF1_dX1 = product(R1, product(-S, R1_T));

        const TinyMatrix<doublereal, 3, 3> ddX_dg1 = product(R1_0_T, skew(TinyVector<doublereal, 3>(X1 + R1o1 - X2 - R2o2))) - product(R1_T, skew_R1_0o1);
        const TinyMatrix<doublereal, 3, 3> ddV_dg1 = product(R1_0_T, skew(TinyVector<doublereal, 3>(V1 + cross(W1, R1o1) - V2 - cross(W2, R2o2))))
                                               + product(R1_T, TinyMatrix<doublereal, 3, 3>(product(skew(R1o1), skew(W1_0)) - product(skew(W1), skew_R1_0o1)));
        const TinyMatrix<doublereal, 3, 3> dF1_dg1 = skew(TinyVector<doublereal, 3>(product(R1_0, TinyVector<doublereal, 3>(-F1_R1))))
                                               - TinyMatrix<doublereal, 3, 3>(product(R1, TinyMatrix<doublereal, 3, 3>(product(S, ddX_dg1) + product(D, ddV_dg1))));

        const TinyMatrix<doublereal, 3, 3> dF1_dX2 = product(R1, product(S, R1_T));
        const TinyMatrix<doublereal, 3, 3> ddX_dg2 = product(R1_T, skew_R2_0o2);
        const TinyMatrix<doublereal, 3, 3> ddV_dg2 = product(R1_T, TinyMatrix<doublereal, 3, 3>(product(skew_R2o2, -skew_W2_0)) + TinyMatrix<doublereal, 3, 3>(product(skew_W2_0, skew_R2_0o2)));
        const TinyMatrix<doublereal, 3, 3> dF1_dg2 = -product(R1, product(S, ddX_dg2) + product(D, ddV_dg2));

        const TinyMatrix<doublereal, 3, 3> dF2_dX1 = -dF1_dX1;
        const TinyMatrix<doublereal, 3, 3> dF2_dg1 = -dF1_dg1;
        const TinyMatrix<doublereal, 3, 3> dF2_dX2 = -dF1_dX2;
        const TinyMatrix<doublereal, 3, 3> dF2_dg2 = -dF1_dg2;

        const TinyMatrix<doublereal, 3, 3> dM1_dX1 = product(skew_R1o1, dF1_dX1);
        const TinyMatrix<doublereal, 3, 3> dM1_dg1 = product(skew(F1), skew_R1_0o1) + product(skew_R1o1, dF1_dg1);
        const TinyMatrix<doublereal, 3, 3> dM1_dX2 = product(skew_R1o1, dF1_dX2);
        const TinyMatrix<doublereal, 3, 3> dM1_dg2 = product(skew_R1o1, dF1_dg2);

        const TinyMatrix<doublereal, 3, 3> dM2_dX1 = product(skew_R2o2, dF2_dX1);
        const TinyMatrix<doublereal, 3, 3> dM2_dg1 = product(skew_R2o2, dF2_dg1);
        const TinyMatrix<doublereal, 3, 3> dM2_dX2 = product(skew_R2o2, dF2_dX2);
        const TinyMatrix<doublereal, 3, 3> dM2_dg2 = product(skew(F2), skew_R2_0o2) + product(skew_R2o2, dF2_dg2);

        const TinyMatrix<doublereal, 3, 3> dF1_dV1 = product(R1, product(-D, R1_T));
        const TinyMatrix<doublereal, 3, 3> ddV_dgP1 = product(-R1_T, skew_R1o1);
        const TinyMatrix<doublereal, 3, 3> dF1_dgP1 = product(R1, product(-D, ddV_dgP1));
        const TinyMatrix<doublereal, 3, 3> dF1_dV2 = product(R1, product(D, R1_T));
        const TinyMatrix<doublereal, 3, 3> ddV_dgP2 = product(R1_T, skew_R2o2);
        const TinyMatrix<doublereal, 3, 3> dF1_dgP2 = product(R1, product(-D, ddV_dgP2));

        const TinyMatrix<doublereal, 3, 3> dM1_dV1 = product(skew_R1o1, dF1_dV1);
        const TinyMatrix<doublereal, 3, 3> dM1_dgP1 = product(skew_R1o1, dF1_dgP1);
        const TinyMatrix<doublereal, 3, 3> dM1_dV2 = product(skew_R1o1, dF1_dV2);
        const TinyMatrix<doublereal, 3, 3> dM1_dgP2 = product(skew_R1o1, dF1_dgP2);

        const TinyMatrix<doublereal, 3, 3> dF2_dV1 = -dF1_dV1;
        const TinyMatrix<doublereal, 3, 3> dF2_dgP1 = -dF1_dgP1;
        const TinyMatrix<doublereal, 3, 3> dF2_dV2 = -dF1_dV2;
        const TinyMatrix<doublereal, 3, 3> dF2_dgP2 = -dF1_dgP2;

        const TinyMatrix<doublereal, 3, 3> dM2_dV1 = product(skew_R2o2, dF2_dV1);
        const TinyMatrix<doublereal, 3, 3> dM2_dgP1 = product(skew_R2o2, dF2_dgP1);
        const TinyMatrix<doublereal, 3, 3> dM2_dV2 = product(skew_R2o2, dF2_dV2);
        const TinyMatrix<doublereal, 3, 3> dM2_dgP2 = product(skew_R2o2, dF2_dgP2);

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                WorkMat.PutCoef(i, j,     -dF1_dV1(i, j)  - dCoef * dF1_dX1(i, j));
                WorkMat.PutCoef(i, j + 3, -dF1_dgP1(i, j) - dCoef * dF1_dg1(i, j));
                WorkMat.PutCoef(i, j + 6, -dF1_dV2(i, j)  - dCoef * dF1_dX2(i, j));
                WorkMat.PutCoef(i, j + 9, -dF1_dgP2(i, j) - dCoef * dF1_dg2(i, j));

                WorkMat.PutCoef(i + 3, j,     -dM1_dV1(i, j)  - dCoef * dM1_dX1(i, j));
                WorkMat.PutCoef(i + 3, j + 3, -dM1_dgP1(i, j) - dCoef * dM1_dg1(i, j));
                WorkMat.PutCoef(i + 3, j + 6, -dM1_dV2(i, j)  - dCoef * dM1_dX2(i, j));
                WorkMat.PutCoef(i + 3, j + 9, -dM1_dgP2(i, j) - dCoef * dM1_dg2(i, j));

                WorkMat.PutCoef(i + 6, j,     -dF2_dV1(i, j)  - dCoef * dF2_dX1(i, j));
                WorkMat.PutCoef(i + 6, j + 3, -dF2_dgP1(i, j) - dCoef * dF2_dg1(i, j));
                WorkMat.PutCoef(i + 6, j + 6, -dF2_dV2(i, j)  - dCoef * dF2_dX2(i, j));
                WorkMat.PutCoef(i + 6, j + 9, -dF2_dgP2(i, j)  - dCoef * dF2_dg2(i, j));

                WorkMat.PutCoef(i + 9, j,     -dM2_dV1(i, j)  - dCoef * dM2_dX1(i, j));
                WorkMat.PutCoef(i + 9, j + 3, -dM2_dgP1(i, j) - dCoef * dM2_dg1(i, j));
                WorkMat.PutCoef(i + 9, j + 6, -dM2_dV2(i, j)  - dCoef * dM2_dX2(i, j));
                WorkMat.PutCoef(i + 9, j + 9, -dM2_dgP2(i, j) - dCoef * dM2_dg2(i, j));
            }
        }
        return WorkMat;
    }

    index_type iGetNumRows() const { return 12; }
    index_type iGetNumCols() const { return NADVars; }
};

TinyMatrix<doublereal, 3, 3> Euler123ToMatR(const TinyVector<doublereal, 3>& v) {
	doublereal d = v(0);
	doublereal dCosAlpha(cos(d));
	doublereal dSinAlpha(sin(d));
	d = v(1);
	doublereal dCosBeta(cos(d));
	doublereal dSinBeta(sin(d));
	d = v(2);
	doublereal dCosGamma(cos(d));
	doublereal dSinGamma(sin(d));

    TinyMatrix<doublereal, 3, 3> R;


    R(0, 0) = dCosBeta*dCosGamma;
    R(1, 0) = dCosAlpha*dSinGamma + dSinAlpha*dSinBeta*dCosGamma;
	R(2, 0) = dSinAlpha*dSinGamma - dCosAlpha*dSinBeta*dCosGamma;
    R(0, 1) = -dCosBeta*dSinGamma;
    R(1, 1) = dCosAlpha*dCosGamma - dSinAlpha*dSinBeta*dSinGamma;
    R(2, 1) = dSinAlpha*dCosGamma + dCosAlpha*dSinBeta*dSinGamma;
    R(0, 2) = dSinBeta;
    R(1, 2) = -dSinAlpha*dCosBeta;
    R(2, 2) = dCosAlpha*dCosBeta;

    return R;
}

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
            TinyVector<doublereal, 3> X0, XP0, Phi0, W0;

            for (int j = 0; j < 3; ++j) {
                X0(j) = ((i + 1) * 10 + j + 1);
                XP0(j) = ((i + 1) * 1000 + (j + 1) * 100);
                Phi0(j) = ((i + 1) * 0.1 + (j + 1) * 0.01);
                W0(j) = ((i + 1) * 0.1 + (j + 1) * 0.01);
            }

            TinyMatrix<doublereal, 3, 3> R0 = Euler123ToMatR(Phi0);

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
            TinyVector<doublereal, 3> o1, o2;

            o1(0) = 1.;
            o1(1) = 2.;
            o1(2) = 3.;
            o2(0) = 4.;
            o2(1) = 5.;
            o2(2) = 6.;

            doublereal s = 100 * (i + 1);
            doublereal d = 10 * (i + 1);
            elements[i] = new Element1(nodes[i], o1, nodes[i + 1], o2, s, d);
        }

        Array<doublereal, 1> XCurr(iNumDof), XPrimeCurr(iNumDof);

        XCurr.initialize(0.);
        XPrimeCurr.initialize(0.);

        for (int i = 0; i < iNumNodes; ++i) {
            nodes[i]->SetValue(XCurr, XPrimeCurr);
        }

        Array<ResItem<doublereal>, 1> WorkVec;
        SparseSubMatrixHandler WorkMatSp;
        FullSubMatrixHandler WorkMatFull;
        Array<doublereal, 1> ResVec;
        Array<doublereal, 2> JacMatAD, JacMat;

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
                std::cout << setw(10) << ResVec(i) << std::endl;
            }

            std::cout << "]" << std::endl;

            std::cout << "JacMatAD = [" << std::endl;

            for (int i = 0; i < iNumDof; ++i) {
                for (int j = 0; j < iNumDof; ++j) {
                    std::cout << setw(10) << std::setprecision(16) << JacMatAD(i, j) << " ";
                }
                std::cout << std::endl;
            }

            std::cout << "]" << std::endl;

            std::cout << "JacMat = [" << std::endl;

            for (int i = 0; i < iNumDof; ++i) {
                for (int j = 0; j < iNumDof; ++j) {
                    std::cout << setw(10) << std::setprecision(16) << JacMat(i, j) << " ";
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

} // namespace
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

    testRangeVector<0>();
    testRangeVector<6>();
    testMapVector<0>();
    testMapVector<19>();

    for (index_type i = 1; i <= 12; ++i) {
    	testGradient<0>(i);
    }

    testGradient<1>(1);
    testGradient<2>(2);
    testGradient<4>(4);
    testGradient<6>(6);
    testGradient<8>(8);
    testGradient<10>(10);
    testGradient<12>(12);
    testGradient2<0>();
    testGradient2<3>();

    std::cerr << "test_func2:" << std::endl;
    
    test_func2<1>(1);
    test_func2<2>(2);
    test_func2<4>(4);
    test_func2<6>(6);
    test_func2<8>(8);
    test_func2<10>(10);
    test_func2<12>(12);
    test_func2<0>(12);
    
    testGradientCopy(0);
    testGradientCopy(1);
    testGradientCopy(2);
    testGradientCopy(4);
    testGradientCopy(6);
    testGradientCopy(8);
    testGradientCopy(10);
    testGradientCopy(15);
    testGradientLin<4>(4);
    testGradientLin<0>(4);

#ifdef HAVE_BLITZ
    gradAssTest::testAssembly();
#endif

    testDifferentDofMaps<3>(1);
    testDifferentDofMaps<4>(2);
    testDifferentDofMaps<5>(3);
    testDifferentDofMaps<6>(4);
    testDifferentDofMaps<7>(5);
    testDifferentDofMaps<8>(6);

    for (integer i = 1; i <= 10; ++i) {
    	testDifferentDofMaps<0>(i);
    	testDifferentDofMaps<12>(i);
    	testDifferentDofMaps<14>(i);
    	testDifferentDofMaps<16>(i);
    }

#ifndef NDEBUG
    std::cerr << "All tests passed" << std::endl;
#else
    std::cerr << "No tests have been done" << std::endl;
#endif

    std::cerr << "GRADIENT_DEBUG=" << GRADIENT_DEBUG << std::endl;

    return 0;
}
