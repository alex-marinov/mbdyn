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

#include "mbconfig.h"
#ifdef USE_AUTODIFF
#include "gradient.h"
#include "matvec.h"
#endif
#include "sp_gradient.h"
#include "sp_gradient_op.h"
#include "sp_gradient_test_func.h"

namespace sp_grad_test {
     template <typename T>
     void func_scalar1(const T& u, const T& v, const T& w, doublereal e, T& f) {
          f = ((((3 * u + 2 * v) * (u - v) / (1 - w) * pow(fabs(u/w), v - 1) * sin(v) * cos(w) * (1 - tan(-w + v))) * e + 1. - 11. + 4.5 - 1.) * 3.5 / 2.8 + u - v) * w / u;
     }

     template
     void func_scalar1<doublereal>(const doublereal& u, const doublereal& v, const doublereal& w, doublereal e, doublereal& f);

     template
     void func_scalar1<SpGradient>(const SpGradient& u, const SpGradient& v, const SpGradient& w, doublereal e, SpGradient& f);
#ifdef USE_AUTODIFF
     template
     void func_scalar1<grad::Gradient<0> >(const grad::Gradient<0>& u, const grad::Gradient<0>& v, const grad::Gradient<0>& w, doublereal e, grad::Gradient<0>& f);
#endif
     template <typename T>
     void func_scalar1_compressed(const T& u, const T& v, const T& w, doublereal e, T& f) {
          f = EvalUnique(((((3 * u + 2 * v) * (u - v) / (1 - w) * pow(fabs(u/w), v - 1) * sin(v) * cos(w) * (1 - tan(-w + v))) * e + 1. - 11. + 4.5 - 1.) * 3.5 / 2.8 + u - v) * w / u);
     }

     template
     void func_scalar1_compressed<doublereal>(const doublereal& u, const doublereal& v, const doublereal& w, doublereal e, doublereal& f);

     template
     void func_scalar1_compressed<SpGradient>(const SpGradient& u, const SpGradient& v, const SpGradient& w, doublereal e, SpGradient& f);

     template <typename T>
     void func_scalar2(const T& u, const T& v, const T& w, doublereal e, T& f) {
          f = u;
          f *= 3.;
          f += 2 * v;
          f *= (u - v);
          f /= (1 - w);
          f *= pow(fabs(u/w), v - 1);
          f *= sin(v);
          f *= cos(w);
          f *= (1 - tan(-w + v));
          f *= e;
          f += 1;
          f -= 11;
          f += 4.5 - 1;
          f *= 3.5;
          f /= 2.8;
          f += u;
          f -= v;
          f *= w;
          f /= u;
     }

     template
     void func_scalar2<doublereal>(const doublereal& u, const doublereal& v, const doublereal& w, doublereal e, doublereal& f);

     template
     void func_scalar2<SpGradient>(const SpGradient& u, const SpGradient& v, const SpGradient& w, doublereal e, SpGradient& f);

     template <typename U, typename V, typename W>
     bool func_bool1(const U& u, const V& v, const W& w, doublereal e) {
          return u + v >= e * (v - w);
     }

     template
     bool func_bool1(const SpGradient& u, const SpGradient& v, const SpGradient& w, doublereal e);

     template
     bool func_bool1(const doublereal& u, const doublereal& v, const doublereal& w, doublereal e);

     template
     bool func_bool1(const SpGradient& u, const doublereal& v, const doublereal& w, doublereal e);

     template
     bool func_bool1(const doublereal& u, const SpGradient& v, const doublereal& w, doublereal e);

     template
     bool func_bool1(const doublereal& u, const doublereal& v, const SpGradient& w, doublereal e);

     template
     bool func_bool1(const SpGradient& u, const SpGradient& v, const doublereal& w, doublereal e);

     template
     bool func_bool1(const doublereal& u, const SpGradient& v, const SpGradient& w, doublereal e);

     template
     bool func_bool1(const SpGradient& u, const doublereal& v, const SpGradient& w, doublereal e);

     template <typename U, typename V, typename W>
     bool func_bool2(const U& u, const V& v, const W& w, doublereal e) {
          return u + v <= e * (v - w);
     }

     template
     bool func_bool2(const SpGradient& u, const SpGradient& v, const SpGradient& w, doublereal e);

     template
     bool func_bool2(const doublereal& u, const doublereal& v, const doublereal& w, doublereal e);

     template
     bool func_bool2(const SpGradient& u, const doublereal& v, const doublereal& w, doublereal e);

     template
     bool func_bool2(const doublereal& u, const SpGradient& v, const doublereal& w, doublereal e);

     template
     bool func_bool2(const doublereal& u, const doublereal& v, const SpGradient& w, doublereal e);

     template
     bool func_bool2(const SpGradient& u, const SpGradient& v, const doublereal& w, doublereal e);

     template
     bool func_bool2(const doublereal& u, const SpGradient& v, const SpGradient& w, doublereal e);

     template
     bool func_bool2(const SpGradient& u, const doublereal& v, const SpGradient& w, doublereal e);

     template <typename U, typename V, typename W>
     bool func_bool3(const U& u, const V& v, const W& w, doublereal e) {
          return u + v > e * (v - w);
     }

     template
     bool func_bool3(const SpGradient& u, const SpGradient& v, const SpGradient& w, doublereal e);

     template
     bool func_bool3(const doublereal& u, const doublereal& v, const doublereal& w, doublereal e);

     template
     bool func_bool3(const SpGradient& u, const doublereal& v, const doublereal& w, doublereal e);

     template
     bool func_bool3(const doublereal& u, const SpGradient& v, const doublereal& w, doublereal e);

     template
     bool func_bool3(const doublereal& u, const doublereal& v, const SpGradient& w, doublereal e);

     template
     bool func_bool3(const SpGradient& u, const SpGradient& v, const doublereal& w, doublereal e);

     template
     bool func_bool3(const doublereal& u, const SpGradient& v, const SpGradient& w, doublereal e);

     template
     bool func_bool3(const SpGradient& u, const doublereal& v, const SpGradient& w, doublereal e);

     template <typename U, typename V, typename W>
     bool func_bool4(const U& u, const V& v, const W& w, doublereal e) {
          return u + v < e * (v - w);
     }

     template
     bool func_bool4(const SpGradient& u, const SpGradient& v, const SpGradient& w, doublereal e);

     template
     bool func_bool4(const doublereal& u, const doublereal& v, const doublereal& w, doublereal e);

     template
     bool func_bool4(const SpGradient& u, const doublereal& v, const doublereal& w, doublereal e);

     template
     bool func_bool4(const doublereal& u, const SpGradient& v, const doublereal& w, doublereal e);

     template
     bool func_bool4(const doublereal& u, const doublereal& v, const SpGradient& w, doublereal e);

     template
     bool func_bool4(const SpGradient& u, const SpGradient& v, const doublereal& w, doublereal e);

     template
     bool func_bool4(const doublereal& u, const SpGradient& v, const SpGradient& w, doublereal e);

     template
     bool func_bool4(const SpGradient& u, const doublereal& v, const SpGradient& w, doublereal e);

     template <typename U, typename V, typename W>
     bool func_bool5(const U& u, const V& v, const W& w, doublereal e) {
          return (u + v - w) * e == e * (v + u - w);
     }

     template
     bool func_bool5(const SpGradient& u, const SpGradient& v, const SpGradient& w, doublereal e);

     template
     bool func_bool5(const doublereal& u, const doublereal& v, const doublereal& w, doublereal e);

     template
     bool func_bool5(const SpGradient& u, const doublereal& v, const doublereal& w, doublereal e);

     template
     bool func_bool5(const doublereal& u, const SpGradient& v, const doublereal& w, doublereal e);

     template
     bool func_bool5(const doublereal& u, const doublereal& v, const SpGradient& w, doublereal e);

     template
     bool func_bool5(const SpGradient& u, const SpGradient& v, const doublereal& w, doublereal e);

     template
     bool func_bool5(const doublereal& u, const SpGradient& v, const SpGradient& w, doublereal e);

     template
     bool func_bool5(const SpGradient& u, const doublereal& v, const SpGradient& w, doublereal e);

     template <typename U, typename V, typename W>
     bool func_bool6(const U& u, const V& v, const W& w, doublereal e) {
          return (u + v - w) * e != e * (v + u - w);
     }

     template
     bool func_bool6(const SpGradient& u, const SpGradient& v, const SpGradient& w, doublereal e);

     template
     bool func_bool6(const doublereal& u, const doublereal& v, const doublereal& w, doublereal e);

     template
     bool func_bool6(const SpGradient& u, const doublereal& v, const doublereal& w, doublereal e);

     template
     bool func_bool6(const doublereal& u, const SpGradient& v, const doublereal& w, doublereal e);

     template
     bool func_bool6(const doublereal& u, const doublereal& v, const SpGradient& w, doublereal e);

     template
     bool func_bool6(const SpGradient& u, const SpGradient& v, const doublereal& w, doublereal e);

     template
     bool func_bool6(const doublereal& u, const SpGradient& v, const SpGradient& w, doublereal e);

     template
     bool func_bool6(const SpGradient& u, const doublereal& v, const SpGradient& w, doublereal e);

     doublereal sec(doublereal x) {
          return 1./cos(x);
     }

     void func_scalar1_dv(const index_type nbdirs,
                          const doublereal u,
                          const doublereal ud[],
                          const doublereal v,
                          const doublereal vd[],
                          const doublereal w,
                          const doublereal wd[],
                          doublereal e,
                          doublereal& f,
                          doublereal fd[],
                          doublereal work[]) {
          auto fabs0d = work;
          auto pwy1d = work + nbdirs;
          auto pwr1d = work + 2 * nbdirs;
          auto arg1d = work + 3 * nbdirs;
          doublereal fabs0;
          doublereal pwy1;
          doublereal pwr1;
          doublereal arg1;
          index_type nd;
          if (u/w >= 0.0) {
               for (nd = 0; nd < nbdirs; ++nd)
                    fabs0d[nd] = (ud[nd]*w-u*wd[nd])/(w*w);
               fabs0 = u/w;
          } else {
               for (nd = 0; nd < nbdirs; ++nd)
                    fabs0d[nd] = -((ud[nd]*w-u*wd[nd])/(w*w));
               fabs0 = -(u/w);
          }
          pwy1 = v - 1;
          pwr1 = pow(fabs0, pwy1);
          arg1 = -w + v;
          for (nd = 0; nd < nbdirs; ++nd) {
               pwy1d[nd] = vd[nd];
               if (fabs0 > 0.0)
                    pwr1d[nd] = pow(fabs0, pwy1)*(log(fabs0)*pwy1d[nd]+pwy1*fabs0d[nd]
                                                  /fabs0);
               else if (fabs0 == 0.0)
                    if (pwy1 == 1.0)
                         pwr1d[nd] = fabs0d[nd];
                    else
                         pwr1d[nd] = 0.0;
               else if (pwy1 == (int)pwy1)
                    pwr1d[nd] = pwy1*pow(fabs0, (pwy1-1))*fabs0d[nd];
               else
                    pwr1d[nd] = 0.0;
               arg1d[nd] = vd[nd] - wd[nd];
               fd[nd] = (((3.5*e*(((((3*ud[nd]+2*vd[nd])*(u-v)+(3*u+2*v)*(ud[nd]-
                                                                          vd[nd]))*(1-w)+(3*u+2*v)*(u-v)*wd[nd])*pwr1*sin(v)/((1-w)*(1-w))+(
                                                                               3*u+2*v)*(u-v)*(pwr1d[nd]*sin(v)+pwr1*vd[nd]*cos(v))/(1-w))*cos(w)
                                  *(1-tan(arg1))+(3*u+2*v)*(u-v)*pwr1*sin(v)*(-(wd[nd]*sin(w)*(1-tan
                                                                                               (arg1)))-cos(w)*arg1d[nd]*(1.0+tan(arg1)*tan(arg1)))/(1-w))/2.8+ud
                           [nd]-vd[nd])*w+(((3*u+2*v)*(u-v)/(1-w)*pwr1*sin(v)*cos(w)*(1-tan(
                                                                                           arg1))*e+1.-11.+4.5-1.)*3.5/2.8+u-v)*wd[nd])*u-(((3*u+2*v)*(u-v)/(
                                                                                                                                                 1-w)*pwr1*sin(v)*cos(w)*(1-tan(arg1))*e+1.-11.+4.5-1.)*3.5/2.8+u-v
                                                                                                )*w*ud[nd])/(u*u);
          }
          f = (((3*u+2*v)*(u-v)/(1-w)*pwr1*sin(v)*cos(w)*(1-tan(arg1))*e+1.-11.+4.5
                -1.)*3.5/2.8+u-v)*w/u;
     }

     template <typename TA, typename TB>
     void func_mat_mul1(typename util::ResultType<TA, TB>::Type& g, const std::vector<TA>& A, const std::vector<TB>& B)
     {
          util::InnerProduct(g, A.begin(), A.end(), 1, B.begin(), B.end(), 1);
     }

     template
     void func_mat_mul1<SpGradient, SpGradient>(SpGradient& g, const std::vector<SpGradient>& A, const std::vector<SpGradient>& B);

     template
     void func_mat_mul1<doublereal, SpGradient>(SpGradient& g, const std::vector<doublereal>& A, const std::vector<SpGradient>& B);

     template
     void func_mat_mul1<SpGradient, doublereal>(SpGradient& g, const std::vector<SpGradient>& A, const std::vector<doublereal>& B);

     template
     void func_mat_mul1<doublereal, doublereal>(doublereal& g, const std::vector<doublereal>& A, const std::vector<doublereal>& B);

     template <typename TA, typename TB>
     void func_mat_mul2(typename util::ResultType<TA, TB>::Type& s, const std::vector<TA>& a, const std::vector<TB>& b) {
          SP_GRAD_ASSERT(a.size() == b.size());

          const index_type n = static_cast<index_type>(a.size());

          typedef typename util::ResultType<TA, TB>::Type T;

          s = T();

          for (index_type i = 0; i < n; ++i) {
               s += a[i] * b[i];
          }
     }

     template
     void func_mat_mul2<doublereal, doublereal>(doublereal& s, const std::vector<doublereal>& a, const std::vector<doublereal>& b);

     template
     void func_mat_mul2<SpGradient, doublereal>(SpGradient& s, const std::vector<SpGradient>& a, const std::vector<doublereal>& b);

     template
     void func_mat_mul2<doublereal, SpGradient>(SpGradient& s, const std::vector<doublereal>& a, const std::vector<SpGradient>& b);

     template
     void func_mat_mul2<SpGradient, SpGradient>(SpGradient& s, const std::vector<SpGradient>& a, const std::vector<SpGradient>& b);


     template <typename TA, typename TB>
     void
     func_mat_mul3(typename util::ResultType<TA, TB>::Type& g, const std::vector<TA>& A, const std::vector<TB>& B)
     {
          util::MapInnerProduct(g, A.begin(), A.end(), 1, B.begin(), B.end(), 1);
     }

     template
     void func_mat_mul3<SpGradient, SpGradient>(SpGradient& s, const std::vector<SpGradient>& A, const std::vector<SpGradient>& B);

     template
     void func_mat_mul3<doublereal, SpGradient>(SpGradient& s, const std::vector<doublereal>& A, const std::vector<SpGradient>& B);

     template
     void func_mat_mul3<SpGradient, doublereal>(SpGradient& s, const std::vector<SpGradient>& A, const std::vector<doublereal>& B);

     template
     void func_mat_mul3<doublereal, doublereal>(doublereal& s, const std::vector<doublereal>& A, const std::vector<doublereal>& B);

     void func_mat_mul1(index_type n, const doublereal a[], const doublereal b[], doublereal& s) {
          index_type i;

          s = 0.;

          for (i = 0; i < n; ++i) {
               s += a[i] * b[i];
          }
     }

     void func_mat_mul1_dv(index_type n, const doublereal a[], const doublereal ad[], const doublereal b[], const doublereal bd[], doublereal& s, doublereal sd[], index_type nbdirs) {
          index_type i, j;

          s = 0.;

          for (i = 0; i < n; ++i) {
               s += a[i] * b[i];
          }

          for (j = 0; j < nbdirs; ++j) {
               sd[j] = 0.;

               for (i = 0; i < n; ++i) {
                    sd[j] += ad[j * n + i] * b[i] + a[i] * bd[j * n + i];
               }
          }
     }

     template <typename TA, typename TX>
     void func_mat_mul4(const std::vector<TA>& A, const std::vector<TX>& x, std::vector<typename util::ResultType<TA, TX>::Type>& b) {
          SP_GRAD_ASSERT(A.size() == x.size() * b.size());

          const index_type iRowsA = static_cast<index_type>(b.size());
          const index_type iColsA = static_cast<index_type>(x.size());

          for (index_type iRow = 0; iRow < iRowsA; ++iRow) {
               util::MapInnerProduct(b[iRow],
                                     A.begin() + iRow,
                                     A.begin() + iRow + iRowsA * iColsA,
                                     iRowsA,
                                     x.begin(),
                                     x.end(),
                                     1);
          }
     }

     template
     void func_mat_mul4<doublereal, doublereal>(const std::vector<doublereal>& A, const std::vector<doublereal>& x, std::vector<doublereal>& b);

     template
     void func_mat_mul4<SpGradient, SpGradient>(const std::vector<SpGradient>& A, const std::vector<SpGradient>& x, std::vector<SpGradient>& b);

     template
     void func_mat_mul4<SpGradient, doublereal>(const std::vector<SpGradient>& A, const std::vector<doublereal>& x, std::vector<SpGradient>& b);

     template <>
     void func_mat_mul4<doublereal, SpGradient>(const std::vector<doublereal>& A, const std::vector<SpGradient>& x, std::vector<SpGradient>& b) {
          SP_GRAD_ASSERT(A.size() == x.size() * b.size());

          const index_type iRowsA = static_cast<index_type>(b.size());
          const index_type iColsA = static_cast<index_type>(x.size());

          SpGradDofStat s;

          for (const auto& xi:x) {
               SpGradient::GetDofStat(xi, s);
          }

          SpGradExpDofMap oDofMap(s);

          for (const auto& xi:x) {
               SpGradient::InsertDof(xi, oDofMap);
          }

          oDofMap.InsertDone();

          for (index_type iRow = 0; iRow < iRowsA; ++iRow) {
               util::MapInnerProduct(b[iRow],
                                     A.begin() + iRow,
                                     A.begin() + iRow + iRowsA * iColsA,
                                     iRowsA,
                                     x.begin(),
                                     x.end(),
                                     1,
                                     oDofMap);
          }
     }

     template <>
     void func_mat_mul4<SpGradient, SpGradient>(const index_type irows,
                                                const index_type icols,
                                                const index_type nbdirs,
                                                const doublereal* const A,
                                                const doublereal* const Ad,
                                                const doublereal* const x,
                                                const doublereal* const xd,
                                                doublereal* const b,
                                                doublereal* const bd)
     {
          func_mat_mul4_gg(irows, icols, nbdirs, A, Ad, x, xd, b, bd);
     }

     template <>
     void func_mat_mul4<doublereal, doublereal>(const index_type irows,
                                                const index_type icols,
                                                const index_type nbdirs,
                                                const doublereal* const A,
                                                const doublereal* const Ad,
                                                const doublereal* const x,
                                                const doublereal* const xd,
                                                doublereal* const b,
                                                doublereal* const bd)
     {
          func_mat_mul4_dd(irows, icols, A, x, b);
     }


     template <>
     void func_mat_mul4<doublereal, SpGradient>(const index_type irows,
                                                const index_type icols,
                                                const index_type nbdirs,
                                                const doublereal* const A,
                                                const doublereal* const Ad,
                                                const doublereal* const x,
                                                const doublereal* const xd,
                                                doublereal* const b,
                                                doublereal* const bd)
     {
          func_mat_mul4_dg(irows, icols, nbdirs, A, x, xd, b, bd);
     }


     template <>
     void func_mat_mul4<SpGradient, doublereal>(const index_type irows,
                                                const index_type icols,
                                                const index_type nbdirs,
                                                const doublereal* const A,
                                                const doublereal* const Ad,
                                                const doublereal* const x,
                                                const doublereal* const xd,
                                                doublereal* const b,
                                                doublereal* const bd)
     {
          func_mat_mul4_gd(irows, icols, nbdirs, A, Ad, x, b, bd);
     }

     void func_mat_mul4_dd(const index_type irows,
                           const index_type icols,
                           const doublereal* const A,
                           const doublereal* const x,
                           doublereal* const b) {
          for (index_type i = 0; i < irows; ++i) {
               b[i] = 0.;
          }

          for (index_type j = 0; j < icols; ++j) {
               for (index_type i = 0; i < irows; ++i) {
                    b[i] += A[j * irows + i] * x[j];
               }
          }
     }

     void func_mat_mul4_dg(const index_type irows,
                           const index_type icols,
                           const index_type nbdirs,
                           const doublereal* const A,
                           const doublereal* const x,
                           const doublereal* const xd,
                           doublereal* const b,
                           doublereal* const bd) {

          for (index_type i = 0; i < irows; ++i) {
               b[i] = 0.;
          }

          for (index_type k = 0; k < nbdirs; ++k) {
               for (index_type i = 0; i < irows; ++i) {
                    bd[k * irows + i] = 0.;
               }
          }

          for (index_type j = 0; j < icols; ++j) {
               for (index_type i = 0; i < irows; ++i) {
                    b[i] += A[j * irows + i] * x[j];
               }
          }

          for (index_type k = 0; k < nbdirs; ++k) {
               for (index_type j = 0; j < icols; ++j) {
                    for (index_type i = 0; i < irows; ++i) {
                         bd[k * irows + i] += A[j * irows + i] * xd[k * icols + j];
                    }
               }
          }
     }

     void func_mat_mul4_gg(const index_type irows,
                           const index_type icols,
                           const index_type nbdirs,
                           const doublereal* const A,
                           const doublereal* const Ad,
                           const doublereal* const x,
                           const doublereal* const xd,
                           doublereal* const b,
                           doublereal* const bd) {

          for (index_type i = 0; i < irows; ++i) {
               b[i] = 0.;
          }

          for (index_type k = 0; k < nbdirs; ++k) {
               for (index_type i = 0; i < irows; ++i) {
                    bd[k * irows + i] = 0.;
               }
          }

          for (index_type j = 0; j < icols; ++j) {
               for (index_type i = 0; i < irows; ++i) {
                    b[i] += A[j * irows + i] * x[j];
               }
          }

          for (index_type k = 0; k < nbdirs; ++k) {
               for (index_type j = 0; j < icols; ++j) {
                    for (index_type i = 0; i < irows; ++i) {
                         bd[k * irows + i] += A[j * irows + i] * xd[k * icols + j] + Ad[(k * icols + j) * irows + i] * x[j];
                    }
               }
          }
     }

     void func_mat_mul4_gd(const index_type irows,
                           const index_type icols,
                           const index_type nbdirs,
                           const doublereal* const A,
                           const doublereal* const Ad,
                           const doublereal* const x,
                           doublereal* const b,
                           doublereal* const bd) {

          for (index_type i = 0; i < irows; ++i) {
               b[i] = 0.;
          }

          for (index_type k = 0; k < nbdirs; ++k) {
               for (index_type i = 0; i < irows; ++i) {
                    bd[k * irows + i] = 0.;
               }
          }

          for (index_type j = 0; j < icols; ++j) {
               for (index_type i = 0; i < irows; ++i) {
                    b[i] += A[j * irows + i] * x[j];
               }
          }

          for (index_type k = 0; k < nbdirs; ++k) {
               for (index_type j = 0; j < icols; ++j) {
                    for (index_type i = 0; i < irows; ++i) {
                         bd[k * irows + i] += Ad[(k * icols + j) * irows + i] * x[j];
                    }
               }
          }
     }

     template <typename TA, typename TX>
     void func_mat_mul4m(const std::vector<TA>& A, const std::vector<TX>& x, std::vector<typename util::ResultType<TA,TX>::Type>& b) {
          SP_GRAD_ASSERT(A.size() == x.size() * b.size());

          const index_type iRowsA = static_cast<index_type>(b.size());
          const index_type iColsA = static_cast<index_type>(x.size());

          for (index_type iRow = 0; iRow < iRowsA; ++iRow) {
               util::InnerProduct(b[iRow],
                                  A.begin() + iRow,
                                  A.begin() + iRow + iRowsA * iColsA,
                                  iRowsA,
                                  x.begin(),
                                  x.end(),
                                  1);
          }
     }

     template
     void func_mat_mul4m<doublereal, doublereal>(const std::vector<doublereal>& A, const std::vector<doublereal>& x, std::vector<doublereal>& b);

     template
     void func_mat_mul4m<doublereal, SpGradient>(const std::vector<doublereal>& A, const std::vector<SpGradient>& x, std::vector<SpGradient>& b);

     template
     void func_mat_mul4m<SpGradient, SpGradient>(const std::vector<SpGradient>& A, const std::vector<SpGradient>& x, std::vector<SpGradient>& b);

     template
     void func_mat_mul4m<SpGradient, doublereal>(const std::vector<SpGradient>& A, const std::vector<doublereal>& x, std::vector<SpGradient>& b);

     template <typename TA, typename TX>
     void func_mat_mul4s(const std::vector<TA>& A, const std::vector<TX>& x, std::vector<typename util::ResultType<TA, TX>::Type>& b) {
          SP_GRAD_ASSERT(A.size() == x.size() * b.size());

          typedef typename util::ResultType<TA, TX>::Type TB;

          const sp_grad::index_type iNumRows = static_cast<index_type>(b.size());
          const sp_grad::index_type iNumCols = static_cast<index_type>(x.size());

          for (sp_grad::index_type i = 0; i < iNumRows; ++i) {
               b[i] = TB();
          }

          for (sp_grad::index_type j = 0; j < iNumCols; ++j) {
               for (sp_grad::index_type i = 0; i < iNumRows; ++i) {
                    b[i] += A[j * iNumRows + i] * x[j];
               }
          }
     }
     
#ifdef USE_AUTODIFF
     template
     void func_mat_mul4s<grad::Gradient<0>, grad::Gradient<0> >(const std::vector<grad::Gradient<0> >&, const std::vector<grad::Gradient<0> >&, std::vector<grad::Gradient<0> >&);

     template
     void func_mat_mul4s<grad::Gradient<0>, doublereal>(const std::vector<grad::Gradient<0> >&, const std::vector<doublereal>&, std::vector<grad::Gradient<0> >&);

     template
     void func_mat_mul4s<doublereal, grad::Gradient<0> >(const std::vector<doublereal>&, const std::vector<grad::Gradient<0> >&, std::vector<grad::Gradient<0> >&);
#endif
     template
     void func_mat_mul4s<SpGradient, SpGradient >(const std::vector<SpGradient >&, const std::vector<SpGradient >&, std::vector<SpGradient >&);

     template
     void func_mat_mul4s<SpGradient, doublereal>(const std::vector<SpGradient >&, const std::vector<doublereal>&, std::vector<SpGradient >&);

     template
     void func_mat_mul4s<doublereal, SpGradient >(const std::vector<doublereal>&, const std::vector<SpGradient >&, std::vector<SpGradient >&);

     template
     void func_mat_mul4s<doublereal, doublereal>(const std::vector<doublereal>&, const std::vector<doublereal>&, std::vector<doublereal>&);

     template <typename TA, typename TX>
     void func_mat_mul4sm(const std::vector<TA>& A, const std::vector<TX>& x, std::vector<typename util::ResultType<TA, TX>::Type>& b) {
          SP_GRAD_ASSERT(A.size() == x.size() * b.size());

          typedef typename util::ResultType<TA, TX>::Type TB;

          const index_type iNumRows = static_cast<index_type>(b.size());
          const index_type iNumCols = static_cast<index_type>(x.size());

          for (index_type i = 0; i < iNumRows; ++i) {
               b[i] = TB();
          }

          for (index_type j = 0; j < iNumCols; ++j) {
               for (index_type i = 0; i < iNumRows; ++i) {
                    b[i] += EvalUnique(A[j * iNumRows + i] * x[j]);
               }
          }
     }

     template
     void func_mat_mul4sm<SpGradient, SpGradient >(const std::vector<SpGradient >&, const std::vector<SpGradient >&, std::vector<SpGradient >&);

     template
     void func_mat_mul4sm<SpGradient, doublereal>(const std::vector<SpGradient >&, const std::vector<doublereal>&, std::vector<SpGradient >&);

     template
     void func_mat_mul4sm<doublereal, SpGradient >(const std::vector<doublereal>&, const std::vector<SpGradient >&, std::vector<SpGradient >&);

     template
     void func_mat_mul4sm<doublereal, doublereal>(const std::vector<doublereal>&, const std::vector<doublereal>&, std::vector<doublereal>&);

     void func_mat_add7_ggg(const index_type n,
                            const index_type nbdirs,
                            const doublereal A[],
                            const doublereal Ad[],
                            const doublereal B[],
                            const doublereal Bd[],
                            const doublereal C[],
                            const doublereal Cd[],
                            doublereal D[],
                            doublereal Dd[]) {

          for (index_type i = 0; i < n; ++i) {
               D[i] = -A[i] / 3. + (B[i] * 5. - C[i] / 4.) * 1.5;
          }

          for (index_type k = 0; k < nbdirs; ++k) {
               for (index_type i = 0; i < n; ++i) {
                    Dd[k * n + i] = -Ad[k * n + i] / 3. + (Bd[k * n + i] * 5. - Cd[k * n + i] / 4.) * 1.5;
               }
          }
     }

     void func_mat_add7_ddd(const index_type n, const doublereal A[], const doublereal B[], const doublereal C[], doublereal D[]) {
          for (index_type i = 0; i < n; ++i) {
               D[i] = -A[i] / 3. + (B[i] * 5. - C[i] / 4.) * 1.5;
          }
     }

     void func_mat_add7_gdd(const index_type n,
                            const index_type nbdirs,
                            const doublereal A[],
                            const doublereal Ad[],
                            const doublereal B[],
                            const doublereal C[],
                            doublereal D[],
                            doublereal Dd[]) {

          for (index_type i = 0; i < n; ++i) {
               D[i] = -A[i] / 3. + (B[i] * 5.- C[i] / 4.) * 1.5;
          }

          for (index_type k = 0; k < nbdirs; ++k) {
               for (index_type i = 0; i < n; ++i) {
                    Dd[k * n + i] = -Ad[k * n + i] / 3.;
               }
          }
     }

     void func_mat_add7_dgd(const index_type n,
                            const index_type nbdirs,
                            const doublereal A[],
                            const doublereal B[],
                            const doublereal Bd[],
                            const doublereal C[],
                            doublereal D[],
                            doublereal Dd[]) {

          for (index_type i = 0; i < n; ++i) {
               D[i] = -A[i] / 3. + (B[i] * 5. - C[i] / 4.) * 1.5;
          }

          for (index_type k = 0; k < nbdirs; ++k) {
               for (index_type i = 0; i < n; ++i) {
                    Dd[k * n + i] = Bd[k * n + i] * 5 * 1.5;
               }
          }
     }

     void func_mat_add7_ddg(const index_type n,
                            const index_type nbdirs,
                            const doublereal A[],
                            const doublereal B[],
                            const doublereal C[],
                            const doublereal Cd[],
                            doublereal D[],
                            doublereal Dd[]) {

          for (index_type i = 0; i < n; ++i) {
               D[i] = -A[i] / 3. + (B[i] * 5.- C[i] / 4.) * 1.5;
          }

          for (index_type k = 0; k < nbdirs; ++k) {
               for (index_type i = 0; i < n; ++i) {
                    Dd[k * n + i] = -Cd[k * n + i] / 4. * 1.5;
               }
          }
     }

     void func_mat_add7_gdg(const index_type n,
                            const index_type nbdirs,
                            const doublereal A[],
                            const doublereal Ad[],
                            const doublereal B[],
                            const doublereal C[],
                            const doublereal Cd[],
                            doublereal D[],
                            doublereal Dd[]) {

          for (index_type i = 0; i < n; ++i) {
               D[i] = -A[i] / 3. + (B[i] * 5.- C[i] / 4.) * 1.5;
          }

          for (index_type k = 0; k < nbdirs; ++k) {
               for (index_type i = 0; i < n; ++i) {
                    Dd[k * n + i] = -Ad[k * n + i] / 3. - Cd[k * n + i] / 4. * 1.5;
               }
          }
     }

     void func_mat_add7_ggd(const index_type n,
                            const index_type nbdirs,
                            const doublereal A[],
                            const doublereal Ad[],
                            const doublereal B[],
                            const doublereal Bd[],
                            const doublereal C[],
                            doublereal D[],
                            doublereal Dd[]) {

          for (index_type i = 0; i < n; ++i) {
               D[i] = -A[i] / 3. + (B[i] * 5.- C[i] / 4.) * 1.5;
          }

          for (index_type k = 0; k < nbdirs; ++k) {
               for (index_type i = 0; i < n; ++i) {
                    Dd[k * n + i] = -Ad[k * n + i] / 3. + Bd[k * n + i] * 5. * 1.5;
               }
          }
     }

     void func_mat_add7_dgg(const index_type n,
                            const index_type nbdirs,
                            const doublereal A[],
                            const doublereal B[],
                            const doublereal Bd[],
                            const doublereal C[],
                            const doublereal Cd[],
                            doublereal D[],
                            doublereal Dd[]) {

          for (index_type i = 0; i < n; ++i) {
               D[i] = -A[i] / 3. + (B[i] * 5.- C[i] / 4.) * 1.5;
          }

          for (index_type k = 0; k < nbdirs; ++k) {
               for (index_type i = 0; i < n; ++i) {
                    Dd[k * n + i] = (Bd[k * n + i] * 5. - Cd[k * n + i] / 4.) * 1.5;
               }
          }
     }

     template <>
     void func_mat_add7<SpGradient, SpGradient, SpGradient>(index_type n,
                                                            index_type nbdirs,
                                                            const doublereal A[],
                                                            const doublereal Ad[],
                                                            const doublereal B[],
                                                            const doublereal Bd[],
                                                            const doublereal C[],
                                                            const doublereal Cd[],
                                                            doublereal D[],
                                                            doublereal Dd[]) {
          func_mat_add7_ggg(n, nbdirs, A, Ad, B, Bd, C, Cd, D, Dd);
     }

     template <>
     void func_mat_add7<doublereal, doublereal, doublereal>(index_type n,
                                                            index_type nbdirs,
                                                            const doublereal A[],
                                                            const doublereal Ad[],
                                                            const doublereal B[],
                                                            const doublereal Bd[],
                                                            const doublereal C[],
                                                            const doublereal Cd[],
                                                            doublereal D[],
                                                            doublereal Dd[]) {
          func_mat_add7_ddd(n, A, B, C, D);
     }

     template <>
     void func_mat_add7<SpGradient, doublereal, doublereal>(index_type n,
                                                            index_type nbdirs,
                                                            const doublereal A[],
                                                            const doublereal Ad[],
                                                            const doublereal B[],
                                                            const doublereal Bd[],
                                                            const doublereal C[],
                                                            const doublereal Cd[],
                                                            doublereal D[],
                                                            doublereal Dd[]) {
          func_mat_add7_gdd(n, nbdirs, A, Ad, B, C, D, Dd);
     }

     template <>
     void func_mat_add7<doublereal, SpGradient, doublereal>(index_type n,
                                                            index_type nbdirs,
                                                            const doublereal A[],
                                                            const doublereal Ad[],
                                                            const doublereal B[],
                                                            const doublereal Bd[],
                                                            const doublereal C[],
                                                            const doublereal Cd[],
                                                            doublereal D[],
                                                            doublereal Dd[]) {
          func_mat_add7_dgd(n, nbdirs, A, B, Bd, C, D, Dd);
     }

     template <>
     void func_mat_add7<doublereal, doublereal, SpGradient>(index_type n,
                                                            index_type nbdirs,
                                                            const doublereal A[],
                                                            const doublereal Ad[],
                                                            const doublereal B[],
                                                            const doublereal Bd[],
                                                            const doublereal C[],
                                                            const doublereal Cd[],
                                                            doublereal D[],
                                                            doublereal Dd[]) {
          func_mat_add7_ddg(n, nbdirs, A, B, C, Cd, D, Dd);
     }

     template <>
     void func_mat_add7<SpGradient, SpGradient, doublereal>(index_type n,
                                                            index_type nbdirs,
                                                            const doublereal A[],
                                                            const doublereal Ad[],
                                                            const doublereal B[],
                                                            const doublereal Bd[],
                                                            const doublereal C[],
                                                            const doublereal Cd[],
                                                            doublereal D[],
                                                            doublereal Dd[]) {
          func_mat_add7_ggd(n, nbdirs, A, Ad, B, Bd, C, D, Dd);
     }

     template <>
     void func_mat_add7<doublereal, SpGradient, SpGradient>(index_type n,
                                                            index_type nbdirs,
                                                            const doublereal A[],
                                                            const doublereal Ad[],
                                                            const doublereal B[],
                                                            const doublereal Bd[],
                                                            const doublereal C[],
                                                            const doublereal Cd[],
                                                            doublereal D[],
                                                            doublereal Dd[]) {
          func_mat_add7_dgg(n, nbdirs, A, B, Bd, C, Cd, D, Dd);
     }

     template <>
     void func_mat_add7<SpGradient, doublereal, SpGradient>(index_type n,
                                                            index_type nbdirs,
                                                            const doublereal A[],
                                                            const doublereal Ad[],
                                                            const doublereal B[],
                                                            const doublereal Bd[],
                                                            const doublereal C[],
                                                            const doublereal Cd[],
                                                            doublereal D[],
                                                            doublereal Dd[]) {
          func_mat_add7_gdg(n, nbdirs, A, Ad, B, C, Cd, D, Dd);
     }


     template <typename TA, typename TB, typename TC, index_type NumRows, index_type NumCols>
     void func_mat_add7(const SpMatrixBase<TA, NumRows, NumCols>& A,
                        const SpMatrixBase<TB, NumRows, NumCols>& B,
                        const SpMatrixBase<TC, NumRows, NumCols>& C,
                        SpMatrixBase<typename util::ResultType<typename util::ResultType<TA, TB>::Type, TC>::Type, NumRows, NumCols>& D) {

          SP_GRAD_ASSERT(A.iGetNumRows() == B.iGetNumRows());
          SP_GRAD_ASSERT(A.iGetNumCols() == B.iGetNumCols());
          SP_GRAD_ASSERT(A.iGetNumRows() == C.iGetNumRows());
          SP_GRAD_ASSERT(A.iGetNumCols() == C.iGetNumCols());

          D = -A / 3. + (B * 5. - C / 4.) * 1.5;
     }

     template
     void func_mat_add7(const SpMatrixBase<SpGradient>&,
                        const SpMatrixBase<SpGradient>&,
                        const SpMatrixBase<SpGradient>&,
                        SpMatrixBase<SpGradient>&);

     template
     void func_mat_add7(const SpMatrixBase<doublereal>&,
                        const SpMatrixBase<doublereal>&,
                        const SpMatrixBase<doublereal>&,
                        SpMatrixBase<doublereal>&);

     template
     void func_mat_add7(const SpMatrixBase<SpGradient>&,
                        const SpMatrixBase<doublereal>&,
                        const SpMatrixBase<doublereal>&,
                        SpMatrixBase<SpGradient>&);

     template
     void func_mat_add7(const SpMatrixBase<doublereal>&,
                        const SpMatrixBase<SpGradient>&,
                        const SpMatrixBase<doublereal>&,
                        SpMatrixBase<SpGradient>&);

     template
     void func_mat_add7(const SpMatrixBase<doublereal>&,
                        const SpMatrixBase<doublereal>&,
                        const SpMatrixBase<SpGradient>&,
                        SpMatrixBase<SpGradient>&);

     template
     void func_mat_add7(const SpMatrixBase<SpGradient>&,
                        const SpMatrixBase<SpGradient>&,
                        const SpMatrixBase<doublereal>&,
                        SpMatrixBase<SpGradient>&);

     template
     void func_mat_add7(const SpMatrixBase<doublereal>&,
                        const SpMatrixBase<SpGradient>&,
                        const SpMatrixBase<SpGradient>&,
                        SpMatrixBase<SpGradient>&);

     template
     void func_mat_add7(const SpMatrixBase<SpGradient>&,
                        const SpMatrixBase<doublereal>&,
                        const SpMatrixBase<SpGradient>&,
                        SpMatrixBase<SpGradient>&);

     template
     void func_mat_add7(const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add7(const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add7(const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add7(const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add7(const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add7(const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add7(const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add7(const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);



     template
     void func_mat_add7(const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add7(const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add7(const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add7(const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add7(const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add7(const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add7(const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add7(const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template <typename TA, typename TB, typename TC, index_type NumRows, index_type NumCols>
     void func_mat_add7a(const SpMatrixBase<TA, NumRows, NumCols>& A,
                         const SpMatrixBase<TB, NumRows, NumCols>& B,
                         const SpMatrixBase<TC, NumRows, NumCols>& C,
                         SpMatrixBase<typename util::ResultType<typename util::ResultType<TA, TB>::Type, TC>::Type, NumRows, NumCols>& D) {

          SP_GRAD_ASSERT(A.iGetNumRows() == B.iGetNumRows());
          SP_GRAD_ASSERT(A.iGetNumCols() == B.iGetNumCols());
          SP_GRAD_ASSERT(A.iGetNumRows() == C.iGetNumRows());
          SP_GRAD_ASSERT(A.iGetNumCols() == C.iGetNumCols());

          D = EvalUnique(-A / 3. + (B * 5. - C / 4.) * 1.5);
     }

     template
     void func_mat_add7a(const SpMatrixBase<SpGradient>&,
                         const SpMatrixBase<SpGradient>&,
                         const SpMatrixBase<SpGradient>&,
                         SpMatrixBase<SpGradient>&);

     template
     void func_mat_add7a(const SpMatrixBase<doublereal>&,
                         const SpMatrixBase<doublereal>&,
                         const SpMatrixBase<doublereal>&,
                         SpMatrixBase<doublereal>&);

     template
     void func_mat_add7a(const SpMatrixBase<SpGradient>&,
                         const SpMatrixBase<doublereal>&,
                         const SpMatrixBase<doublereal>&,
                         SpMatrixBase<SpGradient>&);

     template
     void func_mat_add7a(const SpMatrixBase<doublereal>&,
                         const SpMatrixBase<SpGradient>&,
                         const SpMatrixBase<doublereal>&,
                         SpMatrixBase<SpGradient>&);

     template
     void func_mat_add7a(const SpMatrixBase<doublereal>&,
                         const SpMatrixBase<doublereal>&,
                         const SpMatrixBase<SpGradient>&,
                         SpMatrixBase<SpGradient>&);

     template
     void func_mat_add7a(const SpMatrixBase<SpGradient>&,
                         const SpMatrixBase<SpGradient>&,
                         const SpMatrixBase<doublereal>&,
                         SpMatrixBase<SpGradient>&);

     template
     void func_mat_add7a(const SpMatrixBase<doublereal>&,
                         const SpMatrixBase<SpGradient>&,
                         const SpMatrixBase<SpGradient>&,
                         SpMatrixBase<SpGradient>&);

     template
     void func_mat_add7a(const SpMatrixBase<SpGradient>&,
                         const SpMatrixBase<doublereal>&,
                         const SpMatrixBase<SpGradient>&,
                         SpMatrixBase<SpGradient>&);

     template
     void func_mat_add7a(const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                         const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                         const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                         SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add7a(const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                         const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                         const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                         SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add7a(const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                         const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                         const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                         SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add7a(const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                         const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                         const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                         SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add7a(const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                         const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                         const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                         SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add7a(const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                         const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                         const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                         SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add7a(const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                         const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                         const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                         SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add7a(const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                         const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                         const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                         SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);



     template
     void func_mat_add7a(const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                         const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                         const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                         SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add7a(const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                         const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                         const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                         SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add7a(const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                         const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                         const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                         SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add7a(const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                         const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                         const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                         SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add7a(const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                         const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                         const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                         SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add7a(const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                         const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                         const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                         SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add7a(const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                         const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                         const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                         SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add7a(const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                         const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                         const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                         SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template <typename TA, typename TB, typename TC, index_type NumRows, index_type NumCols>
     void func_mat_add8(const SpMatrixBase<TA, NumRows, NumCols>& A,
                        const SpMatrixBase<TB, NumRows, NumCols>& B,
                        const TC& c,
                        SpMatrixBase<typename util::ResultType<typename util::ResultType<TA, TB>::Type, TC>::Type, NumRows, NumCols>& D) {

          SP_GRAD_ASSERT(A.iGetNumRows() == B.iGetNumRows());
          SP_GRAD_ASSERT(A.iGetNumCols() == B.iGetNumCols());

          D = (A / (2. * exp(c)) + B * 3) * sin(c - 1) + (B / 4. - A * 5.) * (cos(c + 1) / pow(c - 1., 2));
     }

     void func_mat_add8(index_type n,
                        index_type nbdirs,
                        const doublereal A[],
                        const doublereal Ad[],
                        const doublereal B[],
                        const doublereal Bd[],
                        doublereal c,
                        const doublereal cd[],
                        doublereal D[],
                        doublereal Dd[]) {

          for (index_type i = 0; i < n; ++i) {
               D[i] = (A[i] / (2. * exp(c)) + B[i] * 3.) * sin(c - 1.) + (B[i] / 4. - A[i] * 5.) * cos(c + 1.) / pow(c - 1., 2);
          }

          for (index_type k = 0; k < nbdirs; ++k) {
               for (index_type i = 0; i < n; ++i) {
                    Dd[k * n + i] = (Ad[k * n + i] / (2. * exp(c)) - A[i] / (2. * exp(c)) * cd[k]
                                     + Bd[k * n + i] * 3.) * sin(c - 1)
                         + (A[i] / (2. * exp(c)) + B[i] * 3.) * cos(c - 1.) * cd[k]
                         + (Bd[k * n + i] / 4. - Ad[k * n + i] * 5.) * cos(c + 1) / pow(c - 1., 2)
                         + (B[i] / 4. - A[i] * 5.) * ((-sin(c + 1) / pow(c - 1., 2) - 2. * cos(c + 1.) / pow(c - 1., 3)) * cd[k]);
               }
          }
     }

     template
     void func_mat_add8(const SpMatrixBase<SpGradient>&,
                        const SpMatrixBase<SpGradient>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient>&);

     template
     void func_mat_add8(const SpMatrixBase<doublereal>&,
                        const SpMatrixBase<doublereal>&,
                        const doublereal&,
                        SpMatrixBase<doublereal>&);

     template
     void func_mat_add8(const SpMatrixBase<SpGradient>&,
                        const SpMatrixBase<doublereal>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient>&);

     template
     void func_mat_add8(const SpMatrixBase<doublereal>&,
                        const SpMatrixBase<SpGradient>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient>&);

     template
     void func_mat_add8(const SpMatrixBase<doublereal>&,
                        const SpMatrixBase<doublereal>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient>&);

     template
     void func_mat_add8(const SpMatrixBase<SpGradient>&,
                        const SpMatrixBase<SpGradient>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient>&);

     template
     void func_mat_add8(const SpMatrixBase<doublereal>&,
                        const SpMatrixBase<SpGradient>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient>&);

     template
     void func_mat_add8(const SpMatrixBase<SpGradient>&,
                        const SpMatrixBase<doublereal>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient>&);

     template
     void func_mat_add8(const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add8(const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const doublereal&,
                        SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add8(const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add8(const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add8(const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add8(const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add8(const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add8(const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add8(const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add8(const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const doublereal&,
                        SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add8(const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add8(const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add8(const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add8(const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add8(const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add8(const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);


     template <typename TA, typename TB, typename TC, index_type NumRows, index_type NumCols>
     void func_mat_add8a(const SpMatrixBase<TA, NumRows, NumCols>& A,
                        const SpMatrixBase<TB, NumRows, NumCols>& B,
                        const TC& c,
                        SpMatrixBase<typename util::ResultType<typename util::ResultType<TA, TB>::Type, TC>::Type, NumRows, NumCols>& D) {

	  typedef typename util::ResultType<typename util::ResultType<TA, TB>::Type, TC>::Type TE;
	  
          SP_GRAD_ASSERT(A.iGetNumRows() == B.iGetNumRows());
          SP_GRAD_ASSERT(A.iGetNumCols() == B.iGetNumCols());

	  D = A;
	  D /= 2.;
	  D /= exp(c);
	  D += B * 3.;
	  D *= sin(c - 1.);

	  SpMatrixBase<TE, NumRows, NumCols> E = B;
	  
	  E *= 0.25;
	  E -= A * 5.;
	  E *= cos(c + 1.);
	  E /= pow(c - 1., 2);
	  D += E;
     }

     template
     void func_mat_add8a(const SpMatrixBase<SpGradient>&,
                        const SpMatrixBase<SpGradient>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient>&);

     template
     void func_mat_add8a(const SpMatrixBase<doublereal>&,
                        const SpMatrixBase<doublereal>&,
                        const doublereal&,
                        SpMatrixBase<doublereal>&);

     template
     void func_mat_add8a(const SpMatrixBase<SpGradient>&,
                        const SpMatrixBase<doublereal>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient>&);

     template
     void func_mat_add8a(const SpMatrixBase<doublereal>&,
                        const SpMatrixBase<SpGradient>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient>&);

     template
     void func_mat_add8a(const SpMatrixBase<doublereal>&,
                        const SpMatrixBase<doublereal>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient>&);

     template
     void func_mat_add8a(const SpMatrixBase<SpGradient>&,
                        const SpMatrixBase<SpGradient>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient>&);

     template
     void func_mat_add8a(const SpMatrixBase<doublereal>&,
                        const SpMatrixBase<SpGradient>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient>&);

     template
     void func_mat_add8a(const SpMatrixBase<SpGradient>&,
                        const SpMatrixBase<doublereal>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient>&);

     template
     void func_mat_add8a(const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add8a(const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const doublereal&,
                        SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add8a(const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add8a(const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add8a(const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add8a(const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add8a(const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add8a(const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add8a(const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add8a(const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const doublereal&,
                        SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add8a(const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add8a(const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add8a(const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add8a(const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add8a(const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add8a(const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);
     
     template <typename TA, typename TB, typename TC, index_type NumRows, index_type NumCols>
     void func_mat_add8b(const SpMatrixBase<TA, NumRows, NumCols>& A,
                        const SpMatrixBase<TB, NumRows, NumCols>& B,
                        const TC& c,
                        SpMatrixBase<typename util::ResultType<typename util::ResultType<TA, TB>::Type, TC>::Type, NumRows, NumCols>& D) {

	  typedef typename util::ResultType<typename util::ResultType<TA, TB>::Type, TC>::Type TE;
	  
          SP_GRAD_ASSERT(A.iGetNumRows() == B.iGetNumRows());
          SP_GRAD_ASSERT(A.iGetNumCols() == B.iGetNumCols());

	  D = A;
	  D /= 2.;
	  D /= EvalUnique(exp(c));
	  D += EvalUnique(B * 3.);
	  D *= EvalUnique(sin(c - 1.));

	  SpMatrixBase<TE, NumRows, NumCols> E = B;
	  
	  E *= 0.25;
	  E -= EvalUnique(A * 5.);
	  E *= EvalUnique(cos(c + 1.));
	  E /= EvalUnique(pow(c - 1., 2));
	  D += E;
     }

     template
     void func_mat_add8b(const SpMatrixBase<SpGradient>&,
                        const SpMatrixBase<SpGradient>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient>&);

     template
     void func_mat_add8b(const SpMatrixBase<doublereal>&,
                        const SpMatrixBase<doublereal>&,
                        const doublereal&,
                        SpMatrixBase<doublereal>&);

     template
     void func_mat_add8b(const SpMatrixBase<SpGradient>&,
                        const SpMatrixBase<doublereal>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient>&);

     template
     void func_mat_add8b(const SpMatrixBase<doublereal>&,
                        const SpMatrixBase<SpGradient>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient>&);

     template
     void func_mat_add8b(const SpMatrixBase<doublereal>&,
                        const SpMatrixBase<doublereal>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient>&);

     template
     void func_mat_add8b(const SpMatrixBase<SpGradient>&,
                        const SpMatrixBase<SpGradient>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient>&);

     template
     void func_mat_add8b(const SpMatrixBase<doublereal>&,
                        const SpMatrixBase<SpGradient>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient>&);

     template
     void func_mat_add8b(const SpMatrixBase<SpGradient>&,
                        const SpMatrixBase<doublereal>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient>&);

     template
     void func_mat_add8b(const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add8b(const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const doublereal&,
                        SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add8b(const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add8b(const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add8b(const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add8b(const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add8b(const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add8b(const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&);

     template
     void func_mat_add8b(const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add8b(const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const doublereal&,
                        SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add8b(const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add8b(const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add8b(const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add8b(const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add8b(const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);

     template
     void func_mat_add8b(const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&);
     
     template <typename TA, typename TB, typename TC, index_type NumRows, index_type NumCols>
     void func_mat_add9(const SpMatrixBase<TA, NumRows, NumCols>& A,
                        const SpMatrixBase<TB, NumCols, NumRows>& B,
                        const TC& c,
                        SpMatrixBase<typename util::ResultType<typename util::ResultType<TA, TB>::Type, TC>::Type, NumCols, NumRows>& D) {

          SP_GRAD_ASSERT(A.iGetNumRows() == B.iGetNumCols());
          SP_GRAD_ASSERT(A.iGetNumCols() == B.iGetNumRows());

          D = Transpose((Transpose(Transpose(A) / (2. * exp(c))) + Transpose(B * 3)) * sin(c - 1) + (Transpose(B) / 4. - A * 5.) * (cos(c + 1) / pow(c - 1., 2)));
     }
     
     template
     void func_mat_add9(const SpMatrixBase<SpGradient, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC>&,
                        const SpMatrixBase<SpGradient, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC>&);

     template
     void func_mat_add9(const SpMatrixBase<doublereal, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC>&,
                        const SpMatrixBase<doublereal, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC>&,
                        const doublereal&,
                        SpMatrixBase<doublereal, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC>&);

     template
     void func_mat_add9(const SpMatrixBase<SpGradient, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC>&,
                        const SpMatrixBase<doublereal, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC>&);

     template
     void func_mat_add9(const SpMatrixBase<doublereal, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC>&,
                        const SpMatrixBase<SpGradient, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC>&);

     template
     void func_mat_add9(const SpMatrixBase<doublereal, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC>&,
                        const SpMatrixBase<doublereal, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC>&);

     template
     void func_mat_add9(const SpMatrixBase<SpGradient, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC>&,
                        const SpMatrixBase<SpGradient, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC>&);

     template
     void func_mat_add9(const SpMatrixBase<doublereal, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC>&,
                        const SpMatrixBase<SpGradient, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC>&);

     template
     void func_mat_add9(const SpMatrixBase<SpGradient, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC>&,
                        const SpMatrixBase<doublereal, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC>&);

     template
     void func_mat_add9(const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<SpGradient, iNumColsStatic1, iNumRowsStatic1>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumColsStatic1, iNumRowsStatic1>&);

     template
     void func_mat_add9(const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<doublereal, iNumColsStatic1, iNumRowsStatic1>&,
                        const doublereal&,
                        SpMatrixBase<doublereal, iNumColsStatic1, iNumRowsStatic1>&);

     template
     void func_mat_add9(const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<doublereal, iNumColsStatic1, iNumRowsStatic1>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient, iNumColsStatic1, iNumRowsStatic1>&);

     template
     void func_mat_add9(const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<SpGradient, iNumColsStatic1, iNumRowsStatic1>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient, iNumColsStatic1, iNumRowsStatic1>&);

     template
     void func_mat_add9(const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<doublereal, iNumColsStatic1, iNumRowsStatic1>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumColsStatic1, iNumRowsStatic1>&);

     template
     void func_mat_add9(const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<SpGradient, iNumColsStatic1, iNumRowsStatic1>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient, iNumColsStatic1, iNumRowsStatic1>&);

     template
     void func_mat_add9(const SpMatrixBase<doublereal, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<SpGradient, iNumColsStatic1, iNumRowsStatic1>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumColsStatic1, iNumRowsStatic1>&);

     template
     void func_mat_add9(const SpMatrixBase<SpGradient, iNumRowsStatic1, iNumColsStatic1>&,
                        const SpMatrixBase<doublereal, iNumColsStatic1, iNumRowsStatic1>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumColsStatic1, iNumRowsStatic1>&);

     template
     void func_mat_add9(const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<SpGradient, iNumColsStatic2, iNumRowsStatic2>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumColsStatic2, iNumRowsStatic2>&);

     template
     void func_mat_add9(const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<doublereal, iNumColsStatic2, iNumRowsStatic2>&,
                        const doublereal&,
                        SpMatrixBase<doublereal, iNumColsStatic2, iNumRowsStatic2>&);

     template
     void func_mat_add9(const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<doublereal, iNumColsStatic2, iNumRowsStatic2>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient, iNumColsStatic2, iNumRowsStatic2>&);

     template
     void func_mat_add9(const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<SpGradient, iNumColsStatic2, iNumRowsStatic2>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient, iNumColsStatic2, iNumRowsStatic2>&);

     template
     void func_mat_add9(const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<doublereal, iNumColsStatic2, iNumRowsStatic2>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumColsStatic2, iNumRowsStatic2>&);

     template
     void func_mat_add9(const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<SpGradient, iNumColsStatic2, iNumRowsStatic2>&,
                        const doublereal&,
                        SpMatrixBase<SpGradient, iNumColsStatic2, iNumRowsStatic2>&);

     template
     void func_mat_add9(const SpMatrixBase<doublereal, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<SpGradient, iNumColsStatic2, iNumRowsStatic2>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumColsStatic2, iNumRowsStatic2>&);

     template
     void func_mat_add9(const SpMatrixBase<SpGradient, iNumRowsStatic2, iNumColsStatic2>&,
                        const SpMatrixBase<doublereal, iNumColsStatic2, iNumRowsStatic2>&,
                        const SpGradient&,
                        SpMatrixBase<SpGradient, iNumColsStatic2, iNumRowsStatic2>&);

     
     template <typename TA, typename TB>
     void func_mat_mul10(const SpMatrixBase<TA>& A,
                         const SpMatrixBase<TB>& B,
                         SpMatrixBase<typename util::ResultType<TA, TB>::Type>& C) {
          SP_GRAD_ASSERT(A.iGetNumCols() == B.iGetNumRows());

          C = A * B;
     }

     template
     void func_mat_mul10(const SpMatrixBase<SpGradient>&,
                         const SpMatrixBase<SpGradient>&,
                         SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul10(const SpMatrixBase<doublereal>&,
                         const SpMatrixBase<doublereal>&,
                         SpMatrixBase<doublereal>&);

     template
     void func_mat_mul10(const SpMatrixBase<SpGradient>&,
                         const SpMatrixBase<doublereal>&,
                         SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul10(const SpMatrixBase<doublereal>&,
                         const SpMatrixBase<SpGradient>&,
                         SpMatrixBase<SpGradient>&);


     template <typename TA, typename TB>
     void func_mat_mul10_trans(const SpMatrixBase<TA>& A,
                               const SpMatrixBase<TB>& B,
                               SpMatrixBase<typename util::ResultType<TA, TB>::Type>& C) {
          SP_GRAD_ASSERT(A.iGetNumCols() == B.iGetNumRows());

          C = Transpose(Transpose(B) * Transpose(A));
     }

     template
     void func_mat_mul10_trans(const SpMatrixBase<SpGradient>&,
                               const SpMatrixBase<SpGradient>&,
                               SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul10_trans(const SpMatrixBase<doublereal>&,
                               const SpMatrixBase<doublereal>&,
                               SpMatrixBase<doublereal>&);

     template
     void func_mat_mul10_trans(const SpMatrixBase<SpGradient>&,
                               const SpMatrixBase<doublereal>&,
                               SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul10_trans(const SpMatrixBase<doublereal>&,
                               const SpMatrixBase<SpGradient>&,
                               SpMatrixBase<SpGradient>&);

     template <typename TA, typename TB>
     void func_mat_mul10_trans_add(const SpMatrixBase<TA>& A,
                                   const SpMatrixBase<TB>& B,
                                   SpMatrixBase<typename util::ResultType<TA, TB>::Type>& C) {
          SP_GRAD_ASSERT(A.iGetNumCols() == B.iGetNumRows());

          C = Transpose((Transpose(A * B) + Transpose(B) * Transpose(A)) * 0.5);
     }

     template
     void func_mat_mul10_trans_add(const SpMatrixBase<SpGradient>&,
                                   const SpMatrixBase<SpGradient>&,
                                   SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul10_trans_add(const SpMatrixBase<doublereal>&,
                                   const SpMatrixBase<doublereal>&,
                                   SpMatrixBase<doublereal>&);

     template
     void func_mat_mul10_trans_add(const SpMatrixBase<SpGradient>&,
                                   const SpMatrixBase<doublereal>&,
                                   SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul10_trans_add(const SpMatrixBase<doublereal>&,
                                   const SpMatrixBase<SpGradient>&,
                                   SpMatrixBase<SpGradient>&);

     void func_mat_mul10(index_type nra,
                         index_type nca,
                         index_type ncb,
                         index_type nbdirs,
                         const doublereal A[],
                         const doublereal Ad[],
                         const doublereal B[],
                         const doublereal Bd[],
                         doublereal C[],
                         doublereal Cd[]) {
          for (index_type j = 0; j < ncb; ++j) {
               for (index_type i = 0; i < nra; ++i) {
                    doublereal cij = 0.;

                    for (index_type k = 0; k < nca; ++k) {
                         cij += A[i + nra * k] * B[k + nca * j];
                    }

                    C[i + nra * j] = cij;
               }
          }

          for (index_type l = 0; l < nbdirs; ++l) {
               for (index_type j = 0; j < ncb; ++j) {
                    for (index_type i = 0; i < nra; ++i) {
                         doublereal cijd = 0.;

                         for (index_type k = 0; k < nca; ++k) {
                              cijd += Ad[l * nra * nca + i + nra * k] * B[k + nca * j]
                                   + A[i + nra * k] * Bd[l * nca * ncb + k + nca * j];
                         }

                         Cd[l * nra * ncb + i + nra * j] = cijd;
                    }
               }
          }
     }
#ifdef USE_AUTODIFF
     template <typename TA, typename TB>
     void func_mat_mul10(const grad::Matrix<TA, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>& A,
                         const grad::Matrix<TB, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>& B,
                         grad::Matrix<typename util::ResultType<TA, TB>::Type, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>& C) {
          SP_GRAD_ASSERT(A.iGetNumCols() == B.iGetNumRows());

          C = A * B;
     }

     template
     void func_mat_mul10(const grad::Matrix<grad::Gradient<0>, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&,
                         const grad::Matrix<grad::Gradient<0>, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&,
                         grad::Matrix<grad::Gradient<0>, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&);

     template
     void func_mat_mul10(const grad::Matrix<doublereal, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&,
                         const grad::Matrix<doublereal, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&,
                         grad::Matrix<doublereal, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&);

     template
     void func_mat_mul10(const grad::Matrix<grad::Gradient<0>, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&,
                         const grad::Matrix<doublereal, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&,
                         grad::Matrix<grad::Gradient<0>, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&);

     template
     void func_mat_mul10(const grad::Matrix<doublereal, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&,
                         const grad::Matrix<grad::Gradient<0>, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&,
                         grad::Matrix<grad::Gradient<0>, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&);


     template <typename TA, typename TB>
     void func_mat_mul10_dof_map(const grad::Matrix<TA, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>& A,
                                 const grad::Matrix<TB, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>& B,
                                 grad::Matrix<typename util::ResultType<TA, TB>::Type, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>& C) {
          SP_GRAD_ASSERT(A.iGetNumCols() == B.iGetNumRows());

          for (grad::index_type j = 1; j <= C.iGetNumCols(); ++j) {
               for (grad::index_type i = 1; i <= C.iGetNumRows(); ++i) {
                    grad::Reset(C(i, j));
               }

               for (grad::index_type k = 1; k <= A.iGetNumCols(); ++k) {
                    for (grad::index_type i = 1; i <= C.iGetNumRows(); ++i) {
                         C(i, j) += A(i, k) * B(k, j);
                    }
               }
          }
     }

     template <>
     void func_mat_mul10_dof_map<doublereal, grad::Gradient<0>>(const grad::Matrix<doublereal, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>& A,
                                                                const grad::Matrix<grad::Gradient<0>, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>& B,
                                                                grad::Matrix<grad::Gradient<0>, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>& C) {
          SP_GRAD_ASSERT(A.iGetNumCols() == B.iGetNumRows());

          grad::LocalDofMap oDofMap;
          grad::Gradient<0> Bkj;

          for (grad::index_type j = 1; j <= C.iGetNumCols(); ++j) {
               for (grad::index_type i = 1; i <= C.iGetNumRows(); ++i) {
                    grad::Reset(C(i, j));
               }

               for (grad::index_type k = 1; k <= A.iGetNumCols(); ++k) {
                    oDofMap.Reset();
                    Copy(Bkj, B(k, j), &oDofMap);

                    for (grad::index_type i = 1; i <= C.iGetNumRows(); ++i) {
                         C(i, j) += A(i, k) * Bkj;
                    }
               }
          }
     }

     template
     void func_mat_mul10_dof_map(const grad::Matrix<grad::Gradient<0>, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&,
                                 const grad::Matrix<grad::Gradient<0>, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&,
                                 grad::Matrix<grad::Gradient<0>, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&);

     template
     void func_mat_mul10_dof_map(const grad::Matrix<doublereal, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&,
                                 const grad::Matrix<doublereal, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&,
                                 grad::Matrix<doublereal, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&);

     template
     void func_mat_mul10_dof_map(const grad::Matrix<grad::Gradient<0>, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&,
                                 const grad::Matrix<doublereal, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&,
                                 grad::Matrix<grad::Gradient<0>, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&);

     template
     void func_mat_mul10_dof_map(const grad::Matrix<doublereal, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&,
                                 const grad::Matrix<grad::Gradient<0>, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&,
                                 grad::Matrix<grad::Gradient<0>, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&);

     template <typename TA, typename TB>
     void func_mat_mul10_trans_add(const grad::Matrix<TA, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>& A,
                                   const grad::Matrix<TB, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>& B,
                                   grad::Matrix<typename util::ResultType<TA, TB>::Type, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>& C) {
          SP_GRAD_ASSERT(A.iGetNumCols() == B.iGetNumRows());

          typedef typename util::ResultType<TA, TB>::Type TC;
          typedef grad::Matrix<TC, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE> MatProdType;

          const MatProdType AB = A * B;
          const MatProdType B_TA_T = Transpose(B) * Transpose(A);
          const MatProdType AB_T = (Transpose(AB) + B_TA_T) * 0.5;
          C = Transpose(AB_T);
     }

     template
     void func_mat_mul10_trans_add(const grad::Matrix<grad::Gradient<0>, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&,
                                   const grad::Matrix<grad::Gradient<0>, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&,
                                   grad::Matrix<grad::Gradient<0>, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&);

     template
     void func_mat_mul10_trans_add(const grad::Matrix<doublereal, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&,
                                   const grad::Matrix<doublereal, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&,
                                   grad::Matrix<doublereal, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&);

     template
     void func_mat_mul10_trans_add(const grad::Matrix<grad::Gradient<0>, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&,
                                   const grad::Matrix<doublereal, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&,
                                   grad::Matrix<grad::Gradient<0>, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&);

     template
     void func_mat_mul10_trans_add(const grad::Matrix<doublereal, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&,
                                   const grad::Matrix<grad::Gradient<0>, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&,
                                   grad::Matrix<grad::Gradient<0>, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>&);
#endif
     template <typename TA, typename TB, typename TC>
     void func_mat_mul11(const SpMatrixBase<TA>& A,
                         const SpMatrixBase<TB>& B,
                         const SpMatrixBase<TC>& C,
                         SpMatrixBase<typename util::ResultType<typename util::ResultType<TA, TB>::Type, TC>::Type>& D) {
          SP_GRAD_ASSERT(A.iGetNumRows() == B.iGetNumRows());
          SP_GRAD_ASSERT(A.iGetNumCols() == B.iGetNumCols());
          SP_GRAD_ASSERT(A.iGetNumRows() == C.iGetNumRows());
          SP_GRAD_ASSERT(A.iGetNumCols() == C.iGetNumCols());

          D = Transpose(A * 3. - B / 2.) * (B / 4. + C * 5.);
     }

     void func_mat_mul11(index_type n,
                         index_type nbdirs,
                         const doublereal A[],
                         const doublereal Ad[],
                         const doublereal B[],
                         const doublereal Bd[],
                         const doublereal C[],
                         const doublereal Cd[],
                         doublereal D[],
                         doublereal Dd[],
                         doublereal Tmp1[],
                         doublereal Tmp1d[],
                         doublereal Tmp2[],
                         doublereal Tmp2d[]) {
          for (index_type j = 0; j < n; ++j) {
               for (index_type i = 0; i < n; ++i) {
                    Tmp1[i * n + j] = A[j * n + i] * 3. - B[j * n + i] / 2.;
               }
          }

          for (index_type k = 0; k < nbdirs; ++k) {
               for (index_type j = 0; j < n; ++j) {
                    for (index_type i = 0; i < n; ++i) {
                         Tmp1d[k * n * n + i * n + j] = Ad[k * n * n + j * n + i] * 3. - Bd[k * n * n + j * n + i] / 2.;
                    }
               }
          }

          for (index_type j = 0; j < n; ++j) {
               for (index_type i = 0; i < n; ++i) {
                    Tmp2[j * n + i] = B[j * n + i] / 4. + C[j * n + i] * 5.;
               }
          }

          for (index_type k = 0; k < nbdirs; ++k) {
               for (index_type j = 0; j < n; ++j) {
                    for (index_type i = 0; i < n; ++i) {
                         Tmp2d[k * n * n + j * n + i] = Bd[k * n * n + j * n + i] / 4. + Cd[k * n * n + j * n + i] * 5.;
                    }
               }
          }

          func_mat_mul10(n, n, n, nbdirs, Tmp1, Tmp1d, Tmp2, Tmp2d, D, Dd);
     }

     template
     void func_mat_mul11(const SpMatrixBase<SpGradient>&,
                         const SpMatrixBase<SpGradient>&,
                         const SpMatrixBase<SpGradient>&,
                         SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul11(const SpMatrixBase<doublereal>&,
                         const SpMatrixBase<doublereal>&,
                         const SpMatrixBase<doublereal>&,
                         SpMatrixBase<doublereal>&);


     template
     void func_mat_mul11(const SpMatrixBase<SpGradient>&,
                         const SpMatrixBase<doublereal>&,
                         const SpMatrixBase<doublereal>&,
                         SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul11(const SpMatrixBase<doublereal>&,
                         const SpMatrixBase<SpGradient>&,
                         const SpMatrixBase<doublereal>&,
                         SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul11(const SpMatrixBase<doublereal>&,
                         const SpMatrixBase<doublereal>&,
                         const SpMatrixBase<SpGradient>&,
                         SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul11(const SpMatrixBase<SpGradient>&,
                         const SpMatrixBase<SpGradient>&,
                         const SpMatrixBase<doublereal>&,
                         SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul11(const SpMatrixBase<doublereal>&,
                         const SpMatrixBase<SpGradient>&,
                         const SpMatrixBase<SpGradient>&,
                         SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul11(const SpMatrixBase<SpGradient>&,
                         const SpMatrixBase<doublereal>&,
                         const SpMatrixBase<SpGradient>&,
                         SpMatrixBase<SpGradient>&);

     template <typename TA, typename TB, typename TC>
     void func_mat_mul11_trans(const SpMatrixBase<TA>& A,
                               const SpMatrixBase<TB>& B,
                               const SpMatrixBase<TC>& C,
                               SpMatrixBase<typename util::ResultType<typename util::ResultType<TA, TB>::Type, TC>::Type>& D) {
          SP_GRAD_ASSERT(A.iGetNumRows() == B.iGetNumRows());
          SP_GRAD_ASSERT(A.iGetNumCols() == B.iGetNumCols());
          SP_GRAD_ASSERT(A.iGetNumRows() == C.iGetNumRows());
          SP_GRAD_ASSERT(A.iGetNumCols() == C.iGetNumCols());

          D = Transpose(Transpose(B / 4. + C * 5.) * (A * 3. - B / 2.));
     }

     template
     void func_mat_mul11_trans(const SpMatrixBase<SpGradient>&,
                               const SpMatrixBase<SpGradient>&,
                               const SpMatrixBase<SpGradient>&,
                               SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul11_trans(const SpMatrixBase<doublereal>&,
                               const SpMatrixBase<doublereal>&,
                               const SpMatrixBase<doublereal>&,
                               SpMatrixBase<doublereal>&);


     template
     void func_mat_mul11_trans(const SpMatrixBase<SpGradient>&,
                               const SpMatrixBase<doublereal>&,
                               const SpMatrixBase<doublereal>&,
                               SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul11_trans(const SpMatrixBase<doublereal>&,
                               const SpMatrixBase<SpGradient>&,
                               const SpMatrixBase<doublereal>&,
                               SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul11_trans(const SpMatrixBase<doublereal>&,
                               const SpMatrixBase<doublereal>&,
                               const SpMatrixBase<SpGradient>&,
                               SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul11_trans(const SpMatrixBase<SpGradient>&,
                               const SpMatrixBase<SpGradient>&,
                               const SpMatrixBase<doublereal>&,
                               SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul11_trans(const SpMatrixBase<doublereal>&,
                               const SpMatrixBase<SpGradient>&,
                               const SpMatrixBase<SpGradient>&,
                               SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul11_trans(const SpMatrixBase<SpGradient>&,
                               const SpMatrixBase<doublereal>&,
                               const SpMatrixBase<SpGradient>&,
                               SpMatrixBase<SpGradient>&);


     template <typename TA, typename TB, typename TC>
     void func_mat_mul12a(const SpMatrixBase<TA>& A,
                          const SpMatrixBase<TB>& B,
                          const SpMatrixBase<TC>& C,
                          SpMatrixBase<typename util::ResultType<typename util::ResultType<TA, TB>::Type, TC>::Type>& D) {
          D = Transpose(A * B * C);
     }

     void func_mat_mul12(index_type nra,
                         index_type nca,
                         index_type ncb,
                         index_type ncc,
                         index_type nbdirs,
                         const doublereal A[],
                         const doublereal Ad[],
                         const doublereal B[],
                         const doublereal Bd[],
                         const doublereal C[],
                         const doublereal Cd[],
                         doublereal D[],
                         doublereal Dd[],
                         doublereal Tmp1[],
                         doublereal Tmp1d[]) {
          func_mat_mul10(nra, nca, ncb, nbdirs, A, Ad, B, Bd, Tmp1, Tmp1d);
          func_mat_mul10(nra, ncb, ncc, nbdirs, Tmp1, Tmp1d, C, Cd, D, Dd);
     }

     template
     void func_mat_mul12a(const SpMatrixBase<SpGradient>&,
                          const SpMatrixBase<SpGradient>&,
                          const SpMatrixBase<SpGradient>&,
                          SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul12a(const SpMatrixBase<doublereal>&,
                          const SpMatrixBase<doublereal>&,
                          const SpMatrixBase<doublereal>&,
                          SpMatrixBase<doublereal>&);


     template
     void func_mat_mul12a(const SpMatrixBase<SpGradient>&,
                          const SpMatrixBase<doublereal>&,
                          const SpMatrixBase<doublereal>&,
                          SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul12a(const SpMatrixBase<doublereal>&,
                          const SpMatrixBase<SpGradient>&,
                          const SpMatrixBase<doublereal>&,
                          SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul12a(const SpMatrixBase<doublereal>&,
                          const SpMatrixBase<doublereal>&,
                          const SpMatrixBase<SpGradient>&,
                          SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul12a(const SpMatrixBase<SpGradient>&,
                          const SpMatrixBase<SpGradient>&,
                          const SpMatrixBase<doublereal>&,
                          SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul12a(const SpMatrixBase<doublereal>&,
                          const SpMatrixBase<SpGradient>&,
                          const SpMatrixBase<SpGradient>&,
                          SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul12a(const SpMatrixBase<SpGradient>&,
                          const SpMatrixBase<doublereal>&,
                          const SpMatrixBase<SpGradient>&,
                          SpMatrixBase<SpGradient>&);

     template <typename TA, typename TB, typename TC>
     void func_mat_mul12b(const SpMatrixBase<TA>& A,
                          const SpMatrixBase<TB>& B,
                          const SpMatrixBase<TC>& C,
                          SpMatrixBase<typename util::ResultType<typename util::ResultType<TA, TB>::Type, TC>::Type>& D) {
          D = Transpose(C) * Transpose(A * B);
     }

     template
     void func_mat_mul12b(const SpMatrixBase<SpGradient>&,
                          const SpMatrixBase<SpGradient>&,
                          const SpMatrixBase<SpGradient>&,
                          SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul12b(const SpMatrixBase<doublereal>&,
                          const SpMatrixBase<doublereal>&,
                          const SpMatrixBase<doublereal>&,
                          SpMatrixBase<doublereal>&);


     template
     void func_mat_mul12b(const SpMatrixBase<SpGradient>&,
                          const SpMatrixBase<doublereal>&,
                          const SpMatrixBase<doublereal>&,
                          SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul12b(const SpMatrixBase<doublereal>&,
                          const SpMatrixBase<SpGradient>&,
                          const SpMatrixBase<doublereal>&,
                          SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul12b(const SpMatrixBase<doublereal>&,
                          const SpMatrixBase<doublereal>&,
                          const SpMatrixBase<SpGradient>&,
                          SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul12b(const SpMatrixBase<SpGradient>&,
                          const SpMatrixBase<SpGradient>&,
                          const SpMatrixBase<doublereal>&,
                          SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul12b(const SpMatrixBase<doublereal>&,
                          const SpMatrixBase<SpGradient>&,
                          const SpMatrixBase<SpGradient>&,
                          SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul12b(const SpMatrixBase<SpGradient>&,
                          const SpMatrixBase<doublereal>&,
                          const SpMatrixBase<SpGradient>&,
                          SpMatrixBase<SpGradient>&);

     template <typename TA, typename TB, typename TC>
     void func_mat_mul12c(const SpMatrixBase<TA>& A,
                          const SpMatrixBase<TB>& B,
                          const SpMatrixBase<TC>& C,
                          SpMatrixBase<typename util::ResultType<typename util::ResultType<TA, TB>::Type, TC>::Type>& D) {
          D = Transpose(C) * Transpose(B) * Transpose(A);
     }

     template
     void func_mat_mul12c(const SpMatrixBase<SpGradient>&,
                          const SpMatrixBase<SpGradient>&,
                          const SpMatrixBase<SpGradient>&,
                          SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul12c(const SpMatrixBase<doublereal>&,
                          const SpMatrixBase<doublereal>&,
                          const SpMatrixBase<doublereal>&,
                          SpMatrixBase<doublereal>&);


     template
     void func_mat_mul12c(const SpMatrixBase<SpGradient>&,
                          const SpMatrixBase<doublereal>&,
                          const SpMatrixBase<doublereal>&,
                          SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul12c(const SpMatrixBase<doublereal>&,
                          const SpMatrixBase<SpGradient>&,
                          const SpMatrixBase<doublereal>&,
                          SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul12c(const SpMatrixBase<doublereal>&,
                          const SpMatrixBase<doublereal>&,
                          const SpMatrixBase<SpGradient>&,
                          SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul12c(const SpMatrixBase<SpGradient>&,
                          const SpMatrixBase<SpGradient>&,
                          const SpMatrixBase<doublereal>&,
                          SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul12c(const SpMatrixBase<doublereal>&,
                          const SpMatrixBase<SpGradient>&,
                          const SpMatrixBase<SpGradient>&,
                          SpMatrixBase<SpGradient>&);

     template
     void func_mat_mul12c(const SpMatrixBase<SpGradient>&,
                          const SpMatrixBase<doublereal>&,
                          const SpMatrixBase<SpGradient>&,
                          SpMatrixBase<SpGradient>&);
}
