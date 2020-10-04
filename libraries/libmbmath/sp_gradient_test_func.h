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

#ifndef __SP_GRADIENT_TEST_FUNC__H_INCLUDED__
#define __SP_GRADIENT_TEST_FUNC__H_INCLUDED__

#include "sp_gradient.h"

namespace sp_grad {
     namespace util {
#ifdef USE_AUTODIFF
	  template <>
	  struct ResultType<grad::Gradient<0>, grad::Gradient<0> > {
	       typedef grad::Gradient<0> Type;
	  };

	  template <>
	  struct ResultType<grad::Gradient<0>, doublereal> {
	       typedef grad::Gradient<0> Type;
	  };

	  template <>
	  struct ResultType<doublereal, grad::Gradient<0> > {
	       typedef grad::Gradient<0> Type;
	  };
#endif
     }
}

namespace sp_grad_test {
     using namespace sp_grad;

     template <typename T>
     void func_scalar1(const T& u, const T& v, const T& w, doublereal e, T& f);

     template <typename T>
     void func_scalar1_compressed(const T& u, const T& v, const T& w, doublereal e, T& f);

     template <typename T>
     void func_scalar2(const T& u, const T& v, const T& w, doublereal e, T& f);

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
			  doublereal work[]);

     template <typename U, typename V, typename W>
     bool func_bool1(const U& u, const V& v, const W& w, doublereal e);

     template <typename U, typename V, typename W>
     bool func_bool2(const U& u, const V& v, const W& w, doublereal e);

     template <typename U, typename V, typename W>
     bool func_bool3(const U& u, const V& v, const W& w, doublereal e);

     template <typename U, typename V, typename W>
     bool func_bool4(const U& u, const V& v, const W& w, doublereal e);

     template <typename U, typename V, typename W>
     bool func_bool5(const U& u, const V& v, const W& w, doublereal e);

     template <typename U, typename V, typename W>
     bool func_bool6(const U& u, const V& v, const W& w, doublereal e);

     template <typename TA, typename TB>
     void
     func_mat_mul1(typename util::ResultType<TA, TB>::Type& g, const std::vector<TA>& A, const std::vector<TB>& B);

     void func_mat_mul1(index_type n, const doublereal a[], const doublereal b[], doublereal& s);

     void func_mat_mul1_dv(index_type n,
			   const doublereal a[],
			   const doublereal ad[],
			   const doublereal b[],
			   const doublereal bd[],
			   doublereal& s,
			   doublereal sd[],
			   index_type nbdirs);

     template <typename TA, typename TB>
     void
     func_mat_mul3(typename util::ResultType<TA, TB>::Type& g, const std::vector<TA>& A, const std::vector<TB>& B);

     template <typename TA, typename TB>
     void
     func_mat_mul2(typename util::ResultType<TA, TB>::Type& g, const std::vector<TA>& a, const std::vector<TB>& b);

     template <typename TA, typename TX>
     void func_mat_mul4(const std::vector<TA>& A, const std::vector<TX>& x, std::vector<typename util::ResultType<TA, TX>::Type>& b);

     void func_mat_mul4(index_type irows, index_type icols, const doublereal* A, const doublereal* x, doublereal* b);

     void func_mat_mul4_dd(index_type irows,
			   index_type icols,
			   const doublereal* A,
			   const doublereal* x,
			   doublereal* b);


     void func_mat_mul4_dg(index_type irows,
			   index_type icols,
			   index_type nbdirs,
			   const doublereal* A,
			   const doublereal* x,
			   const doublereal* xd,
			   doublereal* b,
			   doublereal* bd);

     void func_mat_mul4_gg(index_type irows,
			   index_type icols,
			   index_type nbdirs,
			   const doublereal* A,
			   const doublereal* Ad,
			   const doublereal* x,
			   const doublereal* xd,
			   doublereal* b,
			   doublereal* bd);

     void func_mat_mul4_gd(index_type irows,
			   index_type icols,
			   index_type nbdirs,
			   const doublereal* A,
			   const doublereal* Ad,
			   const doublereal* x,
			   doublereal* b,
			   doublereal* bd);

     template <typename TA, typename TB>
     void func_mat_mul4(const index_type irows,
			const index_type icols,
			const index_type nbdirs,
			const doublereal* A,
			const doublereal* Ad,
			const doublereal* x,
			const doublereal* xd,
			doublereal* b,
			doublereal* bd);

     template <typename TA, typename TX>
     void func_mat_mul4m(const std::vector<TA>& A, const std::vector<TX>& x, std::vector<typename util::ResultType<TA,TX>::Type>& b);

     template <typename TA, typename TX>
     void func_mat_mul4s(const std::vector<TA>& A, const std::vector<TX>& x, std::vector<typename util::ResultType<TA, TX>::Type>& b);

     template <typename TA, typename TX>
     void func_mat_mul4sm(const std::vector<TA>& A, const std::vector<TX>& x, std::vector<typename util::ResultType<TA, TX>::Type>& b);

     static constexpr index_type iNumRowsStatic1 = 10, iNumColsStatic1 = 8, iNumRowsStatic2 = 7, iNumColsStatic2 = 9;
     
     template <typename TA, typename TB, typename TC, index_type NumRows, index_type NumCols>
     void func_mat_add7(const SpMatrixBase<TA, NumRows, NumCols>& A,
			const SpMatrixBase<TB, NumRows, NumCols>& B,
			const SpMatrixBase<TC, NumRows, NumCols>& C,
			SpMatrixBase<typename util::ResultType<typename util::ResultType<TA, TB>::Type, TC>::Type, NumRows, NumCols>& D);

     template <typename TA, typename TB, typename TC, index_type NumRows, index_type NumCols>
     void func_mat_add7a(const SpMatrixBase<TA, NumRows, NumCols>& A,
			 const SpMatrixBase<TB, NumRows, NumCols>& B,
			 const SpMatrixBase<TC, NumRows, NumCols>& C,
			 SpMatrixBase<typename util::ResultType<typename util::ResultType<TA, TB>::Type, TC>::Type, NumRows, NumCols>& D);

     template <typename TA, typename TB, typename TC>
     void func_mat_add7(index_type n,
			index_type nbdirs,
			const doublereal A[],
			const doublereal Ad[],
			const doublereal B[],
			const doublereal Bd[],
			const doublereal C[],
			const doublereal Cd[],
			doublereal D[],
			doublereal Dd[]);

     template <typename TA, typename TB, typename TC, index_type NumRows, index_type NumCols>
     void func_mat_add8(const SpMatrixBase<TA, NumRows, NumCols>& A,
			const SpMatrixBase<TB, NumRows, NumCols>& B,
			const TC& c,
			SpMatrixBase<typename util::ResultType<typename util::ResultType<TA, TB>::Type, TC>::Type, NumRows, NumCols>& D);

     void func_mat_add8(index_type n,
			index_type nbdirs,
			const doublereal A[],
			const doublereal Ad[],
			const doublereal B[],
			const doublereal Bd[],
			doublereal c,
			const doublereal cd[],
			doublereal D[],
			doublereal Dd[]);

     template <typename TA, typename TB, typename TC, index_type NumRows, index_type NumCols>
     void func_mat_add8a(const SpMatrixBase<TA, NumRows, NumCols>& A,
			 const SpMatrixBase<TB, NumRows, NumCols>& B,
			 const TC& c,
			 SpMatrixBase<typename util::ResultType<typename util::ResultType<TA, TB>::Type, TC>::Type, NumRows, NumCols>& D);

     template <typename TA, typename TB, typename TC, index_type NumRows, index_type NumCols>
     void func_mat_add8b(const SpMatrixBase<TA, NumRows, NumCols>& A,
			 const SpMatrixBase<TB, NumRows, NumCols>& B,
			 const TC& c,
			 SpMatrixBase<typename util::ResultType<typename util::ResultType<TA, TB>::Type, TC>::Type, NumRows, NumCols>& D);
     
     template <typename TA, typename TB, typename TC, index_type NumRows, index_type NumCols>
     void func_mat_add9(const SpMatrixBase<TA, NumRows, NumCols>& A,
			const SpMatrixBase<TB, NumCols, NumRows>& B,
			const TC& c,
			SpMatrixBase<typename util::ResultType<typename util::ResultType<TA, TB>::Type, TC>::Type, NumCols, NumRows>& D);

     template <typename TA, typename TB>
     void func_mat_mul10(const SpMatrixBase<TA>& A,
			 const SpMatrixBase<TB>& B,
			 SpMatrixBase<typename util::ResultType<TA, TB>::Type>& C);
#ifdef USE_AUTODIFF
     template <typename TA, typename TB>
     void func_mat_mul10(const grad::Matrix<TA, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>& A,
			 const grad::Matrix<TB, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>& B,
			 grad::Matrix<typename util::ResultType<TA, TB>::Type, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>& C);
#endif
     void func_mat_mul10(index_type nra,
			 index_type nca,
			 index_type ncb,
			 index_type nbdirs,
			 const doublereal A[],
			 const doublereal Ad[],
			 const doublereal B[],
			 const doublereal Bd[],
			 doublereal C[],
			 doublereal Cd[]);

     template <typename TA, typename TB>
     void func_mat_mul10_trans(const SpMatrixBase<TA>& A,
			       const SpMatrixBase<TB>& B,
			       SpMatrixBase<typename util::ResultType<TA, TB>::Type>& C_T);

     template <typename TA, typename TB>
     void func_mat_mul10_trans_add(const SpMatrixBase<TA>& A,
				   const SpMatrixBase<TB>& B,
				   SpMatrixBase<typename util::ResultType<TA, TB>::Type>& C);
#ifdef USE_AUTODIFF
     template <typename TA, typename TB>
     void func_mat_mul10_trans_add(const grad::Matrix<TA, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>& A,
				   const grad::Matrix<TB, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>& B,
				   grad::Matrix<typename util::ResultType<TA, TB>::Type, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>& C);

     template <typename TA, typename TB>
     void func_mat_mul10_dof_map(const grad::Matrix<TA, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>& A,
				 const grad::Matrix<TB, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>& B,
				 grad::Matrix<typename util::ResultType<TA, TB>::Type, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE>& C);
#endif
     template <typename TA, typename TB, typename TC>
     void func_mat_mul11(const SpMatrixBase<TA>& A,
			 const SpMatrixBase<TB>& B,
			 const SpMatrixBase<TC>& C,
			 SpMatrixBase<typename util::ResultType<typename util::ResultType<TA, TB>::Type, TC>::Type>& D);

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
			 doublereal Tmp2d[]);

     template <typename TA, typename TB, typename TC>
     void func_mat_mul11_trans(const SpMatrixBase<TA>& A,
			       const SpMatrixBase<TB>& B,
			       const SpMatrixBase<TC>& C,
			       SpMatrixBase<typename util::ResultType<typename util::ResultType<TA, TB>::Type, TC>::Type>& D);

     template <typename TA, typename TB, typename TC>
     void func_mat_mul12a(const SpMatrixBase<TA>& A,
			  const SpMatrixBase<TB>& B,
			  const SpMatrixBase<TC>& C,
			  SpMatrixBase<typename util::ResultType<typename util::ResultType<TA, TB>::Type, TC>::Type>& D);

     template <typename TA, typename TB, typename TC>
     void func_mat_mul12b(const SpMatrixBase<TA>& A,
			  const SpMatrixBase<TB>& B,
			  const SpMatrixBase<TC>& C,
			  SpMatrixBase<typename util::ResultType<typename util::ResultType<TA, TB>::Type, TC>::Type>& D);

     template <typename TA, typename TB, typename TC>
     void func_mat_mul12c(const SpMatrixBase<TA>& A,
			  const SpMatrixBase<TB>& B,
			  const SpMatrixBase<TC>& C,
			  SpMatrixBase<typename util::ResultType<typename util::ResultType<TA, TB>::Type, TC>::Type>& D);

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
			 doublereal Tmp1d[]);
}
#endif
