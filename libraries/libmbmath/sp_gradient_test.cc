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

#ifdef HAVE_FEENABLEEXCEPT
#define _GNU_SOURCE 1
#include <fenv.h>
#endif // HAVE_FENV_H

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "matvec3.h"
#include "sp_gradient.h"
#include "sp_matrix_base.h"
#include "sp_gradient_op.h"
#include "sp_matvecass.h"

#ifdef USE_AUTODIFF
#include "gradient.h"
#include "matvec.h"
#endif

#include "sp_gradient_test_func.h"

namespace sp_grad_test {
     inline doublereal* front(std::vector<doublereal>& v) {
	  return v.size() ? &v.front() : nullptr;
     }

     inline const doublereal* front(const std::vector<doublereal>& v) {
	  return v.size() ? &v.front() : nullptr;
     }
     
#ifdef USE_AUTODIFF
	inline void SetDofMap(doublereal& g, grad::LocalDofMap*) {
		g = 0;
	}

	template <grad::index_type NSIZE>
	inline void SetDofMap(grad::Gradient<NSIZE>& g, grad::LocalDofMap* pDofMap) {
		g = grad::Gradient<NSIZE>(0., pDofMap);
	}
#endif
	template <typename RAND_NZ_T, typename RAND_DOF_T, typename RAND_VAL_T, typename GEN_T>
	index_type sp_grad_rand_gen(SpGradient& u, RAND_NZ_T& randnz, RAND_DOF_T& randdof, RAND_VAL_T& randval, GEN_T& gen) {
		std::vector<SpDerivRec> ud;

		index_type N = randnz(gen);

		ud.reserve(N > 0 ? N : 0);

		for (index_type i = 0; i < N; ++i) {
			auto d = randval(gen);

			ud.emplace_back(randdof(gen), d);
		}

		doublereal v = randval(gen);

#ifdef SP_GRAD_RANDGEN_NO_ZEROS
		if (v == 0.) {
			v = sqrt(std::numeric_limits<doublereal>::epsilon());
		}
#endif

		u.Reset(v, ud);

		u.Sort();

		return u.iGetSize();
	}

	template <typename RAND_NZ_T, typename RAND_DOF_T, typename RAND_VAL_T, typename GEN_T>
	index_type sp_grad_rand_gen(doublereal& u, RAND_NZ_T& randnz, RAND_DOF_T& randdof, RAND_VAL_T& randval, GEN_T& gen) {
		u = randval(gen);
		return 0;
	}

	void sp_grad_assert_equal(doublereal u, doublereal v, doublereal dTol) {
		assert(fabs(u - v) / std::max(1., fabs(u) + fabs(v)) < dTol);
	}

	void sp_grad_assert_equal(const SpGradient& u, const SpGradient& v, doublereal dTol) {
		SpGradDofStat s;

		u.GetDofStat(s);
		v.GetDofStat(s);

		sp_grad_assert_equal(u.dGetValue(), v.dGetValue(), dTol);

		for (index_type i = s.iMinDof; i <= s.iMaxDof; ++i) {
			sp_grad_assert_equal(u.dGetDeriv(i), v.dGetDeriv(i), dTol);
		}
	}

	void test0(index_type inumloops, index_type inumnz, index_type inumdof) {
	     using namespace std;
	     using namespace std::chrono;

	     cerr << __PRETTY_FUNCTION__ << ":\n";

	     random_device rd;
	     mt19937 gen(rd());
	     uniform_real_distribution<doublereal> randval(-1., 1.);
	     uniform_int_distribution<index_type> randdof(1, inumdof);
	     uniform_int_distribution<index_type> randnz(0, inumnz - 1);

	     gen.seed(0);

	     for (index_type iloop = 0; iloop < inumloops; ++iloop) {
		  Mat3x3 A, C;
		  Mat3xN A3n(3), C3n(3);
		  MatNx3 An3(3), Cn3(3);
		  MatNxN Ann(3, 3), Cnn(3, 3);
		  Vec3 g;
		  Vec3 b;
		  doublereal d;
		  SpMatrix<doublereal, 2, 2> E2x2(2, 2, 0);
		  SpMatrix<SpGradient, 3, 3> X1(3, 3, 0);
		  SpMatrix<doublereal, 3, 3> X2(3, 3, 0);
		  SpMatrix<SpGradient, 3, 3> X3a(3, 3, 0), X4a(3, 3, 0), X5a(3, 3, 0);
		  SpMatrix<SpGradient, 3, 3> X3b(3, 3, 0), X4b(3, 3, 0), X5b(3, 3, 0);
		  SpMatrixA<SpGradient, 3, 3, 10> X6a, X6b;
		  SpMatrixA<SpGradient, 3, 3> X7a(10), X7b(0);
		  SpColVectorA<SpGradient, 3, 10> X8a, X8b;
		  SpRowVectorA<SpGradient, 3, 10> X9a, X9b;
		  SpMatrix<SpGradient> X10a, X10b, X11a(3, 3, 10), X11b(3, 3, 10);
		  SpMatrixA<SpGradient, 2, 3> X12a, X12b, X12c;
		  SpMatrixA<SpGradient, 3, 2> X13a, X12d;
		  SpRowVectorA<SpGradient, 3> X14a, X14b;
		  SpMatrixA<SpGradient, 3, 3> X15a, X15b;
		  
		  X10a.ResizeReset(3, 3, 10);
		  X10b.ResizeReset(3, 3, 10);
		  
		  sp_grad_rand_gen(d, randnz, randdof, randval, gen);

		  for (index_type i = 1; i <= 2; ++i) {
		       for (index_type j = 1; j <= 2; ++j) {
			    sp_grad_rand_gen(E2x2(i, j), randnz, randdof, randval, gen);
		       }
		  }

		  for (index_type i = 1; i <= 3; ++i) {
		       for (index_type j = 1; j <= 3; ++j) {
			    sp_grad_rand_gen(A(i, j), randnz, randdof, randval, gen);
			    sp_grad_rand_gen(C(i, j), randnz, randdof, randval, gen);
			    sp_grad_rand_gen(X1(i, j), randnz, randdof, randval, gen);
			    sp_grad_rand_gen(X2(i, j), randnz, randdof, randval, gen);
			    Ann(i, j) = An3(i, j) = A3n(i, j) = A(i, j);
			    Cnn(i, j) = Cn3(i, j) = C3n(i, j) = C(i, j);
		       }

		       sp_grad_rand_gen(b(i), randnz, randdof, randval, gen);
		       sp_grad_rand_gen(g(i), randnz, randdof, randval, gen);
		  }

		  for (index_type i = 1; i <= 3; ++i) {
		       for (index_type j = 1; j <= 2; ++j) {
			    sp_grad_rand_gen(X12a(j, i), randnz, randdof, randval, gen);
			    sp_grad_rand_gen(X13a(i, j), randnz, randdof, randval, gen);
		       }
		  }

		  SpMatrix<doublereal, 3, 3> Asp(A), Csp(C);
		  SpMatrix<doublereal, 2, 2> InvE2x2 = Inv(E2x2);
		  SpMatrix<doublereal, 2, 2> InvE2x2_E2x2 = InvE2x2 * E2x2;
		  SpMatrix<doublereal, 3, 3> Asp4{A(1, 1), A(2, 1), A(3, 1),
						  A(1, 2), A(2, 2), A(3, 2),
						  A(1, 3), A(2, 3), A(3, 3)};
		  Mat3x3 A4(A(1, 1), A(2, 1), A(3, 1),
			    A(1, 2), A(2, 2), A(3, 2),
			    A(1, 3), A(2, 3), A(3, 3));
		  SpColVector<doublereal, 3> gsp(g);
		  SpColVector<doublereal, 3> gsp4{g(1), g(2), g(3)};
		  SpRowVector<doublereal, 3> gsp5{g(1), g(2), g(3)};
		  Mat3x3 RDelta(CGR_Rot::MatR, g);
		  Mat3x3 RVec = RotManip::Rot(g);
		  Vec3 g2 = RotManip::VecRot(RVec);
		  Vec3 g3 = -g2;
		  SpMatrix<SpGradient, 3, 3> Asp2{Asp}, Asp3{A};
		  SpMatrix<doublereal, 3, 3> Asp5{Asp2.GetValue()};
		  SpColVector<SpGradient, 3> gsp3{g};
		  Asp2 = Transpose(A);
		  gsp3 = 2 * g / Dot(g, g);
		  SpMatrix<doublereal, 3, 3> RDeltasp = MatRVec(gsp);
		  SpMatrix<doublereal, 3, 3> RVecsp = MatRotVec(gsp);
		  SpColVector<doublereal, 3> gsp2 = VecRotMat(RVecsp);
		  SpMatrix<doublereal, SpMatrixSize::DYNAMIC, 3> Aspn(3, 3, 0), Cspn(3, 3, 0);
		  SpMatrix<doublereal, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC> Aspnn(3, 3, 0), Cspnn(3, 3, 0);
		  SpColVector<SpGradient, 3> xsp{g};
		  xsp /= Norm(xsp);
		  SpGradient nx = Norm(xsp);


		  for (index_type i = 1; i <= 3; ++i) {
		       for (index_type j = 1; j <= 3; ++j) {
			    Aspnn(i, j) = Aspn(i, j) = A(i, j);
			    Cspnn(i, j) = Cspn(i, j) = C(i, j);
		       }
		  }

		  SpColVector<doublereal, 3> bsp(b);
		  SpColVector<doublereal, SpMatrixSize::DYNAMIC> bspn(3, 0);

		  Asp = A;
		  bsp = b;

		  for (index_type i = 1; i <= 3; ++i) {
		       bspn(i) = bsp(i);
		  }

		  SpColVector<doublereal, 3> csp1 = Asp * b;
		  SpColVector<doublereal, 3> csp2 = A * bsp;
		  SpColVector<doublereal, 3> csp3 = Asp * bsp;
		  SpColVector<doublereal, 3> csp4{A * b};

		  SpColVector<doublereal, 3> q = csp1;
		  q = 3. * q + 1.5 * q;

		  Vec3 c = A * b;

		  Mat3x3 D = A * C;

		  SpMatrix<doublereal, 3, 3> Dsp1{A * C};
		  SpMatrix<doublereal, 3, 3> Dsp2 = Asp * C;
		  SpMatrix<doublereal, 3, 3> Dsp3 = A * Csp;
		  SpMatrix<doublereal, 3, 3> Dsp4 = Asp * Csp;
		  SpMatrix<doublereal, 3, SpMatrixSize::DYNAMIC> Dsp5 = Asp * C3n;
		  SpMatrix<doublereal, 3, 3> Dsp6 = A3n * Cspn;
		  SpMatrix<doublereal, 3, 3> Dsp7 = Transpose(Transpose(Csp) * Transpose(Asp));
		  SpMatrix<doublereal, 3, 3> Dsp8 = Transpose(C.Transpose() * Transpose(Asp));
		  SpMatrix<doublereal, 3, 3> Dsp9 = Transpose(Transpose(Csp) * A.Transpose());
		  SpMatrix<doublereal, 3, SpMatrixSize::DYNAMIC> Dsp10 = Transpose(Transpose(C3n) * Transpose(Asp));
		  SpMatrix<doublereal, 3, 3> Dsp11 = Transpose(Transpose(Cspn) * Transpose(A3n));
		  SpMatrix<doublereal, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC> Dsp12 = An3 * C3n;
		  SpMatrix<doublereal, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC> Dsp13 = Transpose(Transpose(C3n) * Transpose(An3));
		  SpMatrix<doublereal, 3, 3> Dsp14 = A3n * Cn3;
		  SpMatrix<doublereal, 3, 3> Dsp15 = Transpose(Transpose(Cn3) * Transpose(A3n));
		  SpMatrix<doublereal, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC> Dsp16 = Ann * Cnn;
		  SpMatrix<doublereal, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC> Dsp17 = Aspnn * Cspnn;
		  SpMatrix<doublereal, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC> Dsp18 = Aspnn * Cnn;
		  SpMatrix<doublereal, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC> Dsp19 = Ann * Cspnn;
		  SpMatrix<doublereal, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC> Dsp20 = Transpose(Transpose(Cspnn) * Transpose(Aspnn));
		  SpMatrix<doublereal, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC> Dsp21 = Transpose(Transpose(Cnn) * Transpose(Ann));
		  SpMatrix<doublereal, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC> Dsp22 = Transpose(Transpose(Cspnn) * Transpose(Ann));
		  SpMatrix<doublereal, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC> Dsp23 = Transpose(Transpose(Cnn) * Transpose(Aspnn));

		  Mat3x3 E = A.Transpose() * C;
		  SpMatrix<doublereal, 3, 3> Esp1 = Transpose(Asp) * Csp;
		  SpMatrix<doublereal, 3, 3> Esp2 = A.Transpose() * Csp;
		  SpMatrix<doublereal, 3, 3> Esp3 = Transpose(Asp) * C;
		  SpMatrix<doublereal, 3, 3> Esp4{A.Transpose() * C};
		  SpMatrix<doublereal, 3, 3> Esp5 = Transpose(A) * Csp;
		  SpMatrix<doublereal, 3, 3> Esp6 = Transpose(A) * C;
		  SpMatrix<doublereal, 3, 3> Esp7 = Transpose(Transpose(C) * A);
		  SpMatrix<doublereal, 3, 3> Esp8 = Transpose(Transpose(Csp) * Asp);
		  SpMatrix<doublereal, 3, 3> Esp9 = Transpose(Transpose(C) * Asp);
		  SpMatrix<doublereal, 3, 3> Esp10 = Transpose(Transpose(Csp) * A);

		  Mat3x3 F = A.Transpose() + C;
		  SpMatrix<doublereal, 3, 3> Fsp1 = Transpose(Asp) + Csp;
		  SpMatrix<doublereal, 3, 3> Fsp2 = A.Transpose() + Csp;
		  SpMatrix<doublereal, 3, 3> Fsp3 = Transpose(Asp) + C;
		  SpMatrix<doublereal, 3, 3> Fsp4{A.Transpose() + C};

		  Mat3x3 G = A.Transpose() - C;
		  SpMatrix<doublereal, 3, 3> Gsp1 = Transpose(Asp) - Csp;
		  SpMatrix<doublereal, 3, 3> Gsp2 = A.Transpose() - Csp;
		  SpMatrix<doublereal, 3, 3> Gsp3 = Transpose(Asp) - C;
		  SpMatrix<doublereal, 3, 3> Gsp4{A.Transpose() - C};

		  Mat3x3 H = A.Transpose() + C;
		  SpMatrix<doublereal, 3, 3> Hsp1 = Transpose(Asp) + Csp;
		  SpMatrix<doublereal, 3, 3> Hsp2 = A.Transpose() + Csp;
		  SpMatrix<doublereal, 3, 3> Hsp3 = Transpose(Asp) + C;
		  SpMatrix<doublereal, 3, 3> Hsp4{A.Transpose() + C};

		  Mat3x3 I = A.Transpose() * d;
		  SpMatrix<doublereal, 3, 3> Isp1 = Transpose(Asp) * d;
		  SpMatrix<doublereal, 3, 3> Isp2{A.Transpose() * d};
		  SpMatrix<doublereal, 3, 3> Isp3 = d * Transpose(Asp);

		  Mat3x3 J = C.Transpose() * (A.Transpose() * 5. * (C * 0.5 + b.Tens(b) / b.Dot() * 2.5)).Transpose() * (3. / 2.) * C;
		  SpMatrix<doublereal, 3, 3> Jsp1 = Transpose(Csp) * Transpose(5. * Transpose(Asp) * (Csp * 0.5 + bsp * Transpose(bsp) / Dot(bsp, bsp) * 2.5)) * (3. / 2.) * Csp;
		  SpMatrix<doublereal, 3, 3> Jsp2 = C.Transpose() * Transpose(5. * Transpose(Asp) * (Csp * 0.5 + bsp * Transpose(bsp) / Dot(bsp, bsp) * 2.5)) * (3. / 2.) * Csp;
		  SpMatrix<doublereal, 3, 3> Jsp3 = Transpose(Csp) * Transpose(5. * Transpose(Asp) * (C * 0.5 + bsp * Transpose(bsp) / Dot(bsp, bsp) * 2.5)) * (3. / 2.) * Csp;
		  SpMatrix<doublereal, 3, 3> Jsp4 = Transpose(Csp) * Transpose(5. * Transpose(Asp) * (Csp * 0.5 + bsp * Transpose(bsp) / Dot(bsp, bsp) * 2.5)) * (3. / 2.) * C;
		  SpMatrix<doublereal, 3, 3> Jsp5 = Transpose(Csp) * Transpose(A.Transpose() * 5. * (Csp * 0.5 + bsp * Transpose(bsp) / Dot(bsp, bsp) * 2.5)) * (3. / 2.) * Csp;
		  SpMatrix<doublereal, 3, 3> Jsp6 = Transpose(Csp) * Transpose(5. * Transpose(Asp) * (Csp * 0.5 + b * Transpose(bsp) / Dot(bsp, bsp) * 2.5)) * (3. / 2.) * Csp;
		  SpMatrix<doublereal, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC> Jsp7 = Transpose(Cnn) * Transpose(5. * Transpose(Ann) * (Cnn * 0.5 + bspn * Transpose(bspn) / Dot(bspn, bspn) * 2.5)) * (3. / 2.) * Cnn;

		  Mat3x3 K = A.Transpose() / d;
		  SpMatrix<doublereal, 3, 3> Ksp1 = Transpose(Asp) / d;
		  SpMatrix<doublereal, 3, 3> Ksp2{A.Transpose() / d};

		  SpMatrix<doublereal, 3, 3> L = b.Tens();
		  SpMatrix<doublereal, 3, 3> Lsp1 = bsp * Transpose(bsp);
		  SpMatrix<doublereal, SpMatrixSize::DYNAMIC, SpMatrixSize::DYNAMIC> Lsp2 = bspn * Transpose(bspn);

		  constexpr doublereal dTol = sqrt(std::numeric_limits<doublereal>::epsilon());

		  sp_grad_assert_equal(nx.dGetValue(), 1.0, dTol);

		  for (index_type i = 1; i <= 2; ++i) {
		       for (index_type j = 1; j <= 2; ++j) {
			    sp_grad_assert_equal(InvE2x2_E2x2(i, j), ::Eye3(i, j), dTol);
		       }
		  }

		  X3a = X1 + X2;
		  X3b = X1;
		  X3b += X2;
		  X4a = X1 * 2. + X2 / 3.;
		  X4b = X1;
		  X4b *= 2.;
		  X4b += X2 / 3.;
		  X5a = (X1 + X2) * 3. + (X3a - X4a) / 4. + (X1 - X2) / 2.;
		  X5b = X1;
		  X5b += X2;
		  X5b *= 3.;
		  X5b += X3b / 4.;
		  X5b -= X4b * 0.25;
		  X5b += X1 * 0.5;
		  X5b -= X2 / 2.;
		  X6a = (X1 + X2) * 3. + (X3a - X4a) / 4. + (X1 - X2) / 2.;
		  X6b = X1;
		  X6b += X2;
		  X6b *= 3.;
		  X6b += X3b / 4.;
		  X6b -= X4b * 0.25;
		  X6b += X1 * 0.5;
		  X6b -= X2 / 2.;
		  X7a = (X1 + X2) * 3. + (X3a - X4a) / 4. + (X1 - X2) / 2.;
		  X7b = X1;
		  X7b += X2;
		  X7b *= 3.;
		  X7b += X3b / 4.;
		  X7b -= X4b * 0.25;
		  X7b += X1 * 0.5;
		  X7b -= X2 / 2.;

		  X8a = (X1.GetCol(3) + X2.GetCol(3)) * 3. + (X3a.GetCol(3) - X4a.GetCol(3)) / 4. + (X1.GetCol(3) - X2.GetCol(3)) / 2.;

		  X8b = X1.GetCol(3);
		  X8b += X2.GetCol(3);
		  X8b *= 3.;
		  X8b += X3b.GetCol(3) / 4.;
		  X8b -= X4b.GetCol(3) * 0.25;
		  X8b += X1.GetCol(3) * 0.5;
		  X8b -= X2.GetCol(3) / 2.;

		  X9a = (Transpose(X1.GetCol(3)) + Transpose(X2.GetCol(3))) * 3. + (Transpose(X3a.GetCol(3)) - Transpose(X4a.GetCol(3))) / 4. + Transpose(X1.GetCol(3) - X2.GetCol(3)) / 2.;

		  X9b = Transpose(X1.GetCol(3));
		  X9b += Transpose(X2.GetCol(3));
		  X9b *= 3.;
		  X9b += Transpose(X3b.GetCol(3) / 4.);
		  X9b -= Transpose(X4b.GetCol(3) * 0.25);
		  X9b += Transpose(X1.GetCol(3)) * 0.5;
		  X9b -= Transpose(X2.GetCol(3)) / 2.;

		  X12a = Transpose(X13a) * 3.;
		  X12b = Transpose(X13a);
		  X12b += Transpose(X13a) * 2.;
		  X12c = Transpose(X13a) / Norm(X13a.GetCol(1));
		  X12c *= 3. * Norm(X13a.GetCol(1));
		  X12d += Transpose(Transpose(X13a) / Norm(X13a.GetCol(1)));
		  X12d *= 3. * Norm(X13a.GetCol(1));
		  
		  X14a = Transpose(X13a.GetCol(1) + X13a.GetCol(2));

		  for (index_type i = 1; i <= X13a.iGetNumCols(); ++i) {
		       X14b += Transpose(X13a.GetCol(i));
		  }

		  X15a = Transpose(X1) * X2 + X3a * Transpose(X4a) - X5a / Norm(X5a.GetCol(1)) + X2 * Transpose(X2);
		  X15b += Transpose(X1) * X2;
		  X15b += X3a * Transpose(X4a);
		  X15b -= X5a / Norm(X5a.GetCol(1));
		  X15b += X2 * Transpose(X2);
		  
		  for (index_type i = 0; i < 10; ++i) {
		       SpMatrix<SpGradient, 3, 3> X3c{X3a};
		       SpGradient X3a11 = X3c(1, 1), X3a0;

		       X3a11 *= 3.;
		       X3a11 += X3a0;
		       X3a0 += X3a11;

		       {
			    SpGradient X3a01 = std::move(X3a0);
			    SpGradient X3a02 = std::move(X3a01);
			    X3a0 = std::move(X3a02);		     
			    X3a11 = X3a0;
			    X3a11 = X3c(1, 1);
			    X3a01 = std::move(X3a11);
			    X3a02 = std::move(X3a01);
			    X3a11 = std::move(X3a02);
		       }
		       SpMatrix<SpGradient, 3, 3> X3d = std::move(X3c);
		       SpMatrix<SpGradient, 3, 3> X3e = std::move(X3d);
		       SpGradient X3a03 = X3e(1, 1);
		       SpGradient X3a04 = X3a03 + X3a11;
		       SpGradient X3a05 = std::move(X3a04);
		       SpGradient X3a07 = std::move(X3e(1,1)) + std::move(X3e(2,2));
		       SpGradient X3a08 = X3e(1,1) + X3e(2,2) + X3e(3,3);
		       SpGradient X3a09 = std::move(X3e(1, 1));
		       {
			    SpGradient X3e10 = std::move(X3a09);
		       }
		       X3e(1,1) = std::move(X3a09);
		       X3e(2,2) = std::move(X3a05);
		       X3e(3,3) += X3e(1, 1) + X3e(2, 2) + 2. * X3e(3, 3);
		       X3a07 = std::move(X3a08);
		       SpMatrix<SpGradient, 3, 3> X3f = X3e;
		       X3e += X3f;
		       X3f = std::move(X3e);
		       {
			    SpGradient X3e11;
			    X3a07 += X3e11;
			    X3e11 *= 10.;
			    X3e11 += 10.;
			    X3a07 = std::move(X3e11);
		       }
		       X3e(1,1) = std::move(X3a07);

		       {
			    std::vector<SpGradient> rgVec;
		       
			    rgVec.reserve(10);

			    for (index_type i = 1; i <= 3; ++i) {
				 for (index_type j = 1; j <= 3; ++j) {
				      rgVec.emplace_back(std::move(X3e(i,j)));
				 }
			    }

			    rgVec.push_back(rgVec.front());
			    rgVec.push_back(rgVec.back());
			    rgVec.push_back(X3e(1, 1));
			    rgVec.emplace_back(X3e(3, 2));
			    rgVec.push_back(X3e(1, 1) + X3e(2,2));
			    rgVec.emplace_back(X3e(1, 1) + X3e(2,2));
			    X3e(3, 3) = std::move(rgVec.front());
			    X3e(1, 2) = std::move(rgVec[3]);
		       }
		  }
		  for (index_type i = 1; i <= 3; ++i) {
		       for (index_type j = 1; j <= 3; ++j) {
			    X10a(i, j) = (X1(i, j) + X2(i, j)) * 3. + (X3a(i, j) - X4a(i, j)) / 4. + (X1(i, j) - X2(i, j)) / 2.;
			    X10b(i, j) = X1(i, j);
			    X10b(i, j) += X2(i, j);
			    X10b(i, j) *= 3.;
			    X10b(i, j) += X3b(i, j) / 4.;
			    X10b(i, j) -= X4b(i, j) * 0.25;
			    X10b(i, j) += X1(i, j) * 0.5;
			    X10b(i, j) -= X2(i, j) / 2.;		   
		  
			    X11a(i, j) = (X1(i, j) + X2(i, j)) * 3. + (X3a(i, j) - X4a(i, j)) / 4. + (X1(i, j) - X2(i, j)) / 2.;
			    X11b(i, j) = X1(i, j);
			    X11b(i, j) += X2(i, j);
			    X11b(i, j) *= 3.;
			    X11b(i, j) += X3b(i, j) / 4.;
			    X11b(i, j) -= X4b(i, j) * 0.25;
			    X11b(i, j) += X1(i, j) * 0.5;
			    X11b(i, j) -= X2(i, j) / 2.;
		       }
		  }
		  
		  for (index_type i = 1; i <= 3; ++i) {
		       sp_grad_assert_equal(csp1(i), c(i), dTol);
		       sp_grad_assert_equal(csp2(i), c(i), dTol);
		       sp_grad_assert_equal(csp3(i), c(i), dTol);
		       sp_grad_assert_equal(csp4(i), c(i), dTol);
		       sp_grad_assert_equal(q(i), 4.5 * c(i), dTol);
		       sp_grad_assert_equal(gsp2(i), gsp(i), dTol);
		       sp_grad_assert_equal(gsp2(i), g(i), dTol);
		       sp_grad_assert_equal(gsp4(i), g(i), dTol);
		       sp_grad_assert_equal(gsp5(i), g(i), dTol);

		       for (index_type j = 1; j <= 3; ++j) {
			    sp_grad_assert_equal(Dsp1(i, j), D(i, j), dTol);
			    sp_grad_assert_equal(Dsp2(i, j), D(i, j), dTol);
			    sp_grad_assert_equal(Dsp3(i, j), D(i, j), dTol);
			    sp_grad_assert_equal(Dsp4(i, j), D(i, j), dTol);
			    sp_grad_assert_equal(Dsp5(i, j), D(i, j), dTol);
			    sp_grad_assert_equal(Dsp6(i, j), D(i, j), dTol);
			    sp_grad_assert_equal(Dsp7(i, j), D(i, j), dTol);
			    sp_grad_assert_equal(Dsp8(i, j), D(i, j), dTol);
			    sp_grad_assert_equal(Dsp9(i, j), D(i, j), dTol);
			    sp_grad_assert_equal(Dsp10(i, j), D(i, j), dTol);
			    sp_grad_assert_equal(Dsp11(i, j), D(i, j), dTol);
			    sp_grad_assert_equal(Dsp12(i, j), D(i, j), dTol);
			    sp_grad_assert_equal(Dsp13(i, j), D(i, j), dTol);
			    sp_grad_assert_equal(Dsp14(i, j), D(i, j), dTol);
			    sp_grad_assert_equal(Dsp15(i, j), D(i, j), dTol);
			    sp_grad_assert_equal(Dsp16(i, j), D(i, j), dTol);
			    sp_grad_assert_equal(Dsp17(i, j), D(i, j), dTol);
			    sp_grad_assert_equal(Dsp18(i, j), D(i, j), dTol);
			    sp_grad_assert_equal(Dsp19(i, j), D(i, j), dTol);
			    sp_grad_assert_equal(Dsp20(i, j), D(i, j), dTol);
			    sp_grad_assert_equal(Dsp21(i, j), D(i, j), dTol);
			    sp_grad_assert_equal(Dsp22(i, j), D(i, j), dTol);
			    sp_grad_assert_equal(Dsp23(i, j), D(i, j), dTol);

			    sp_grad_assert_equal(Esp1(i, j), E(i, j), dTol);
			    sp_grad_assert_equal(Esp2(i, j), E(i, j), dTol);
			    sp_grad_assert_equal(Esp3(i, j), E(i, j), dTol);
			    sp_grad_assert_equal(Esp4(i, j), E(i, j), dTol);
			    sp_grad_assert_equal(Esp5(i, j), E(i, j), dTol);
			    sp_grad_assert_equal(Esp6(i, j), E(i, j), dTol);
			    sp_grad_assert_equal(Esp7(i, j), E(i, j), dTol);
			    sp_grad_assert_equal(Esp8(i, j), E(i, j), dTol);
			    sp_grad_assert_equal(Esp9(i, j), E(i, j), dTol);
			    sp_grad_assert_equal(Esp10(i, j), E(i, j), dTol);

			    sp_grad_assert_equal(Fsp1(i, j), F(i, j), dTol);
			    sp_grad_assert_equal(Fsp2(i, j), F(i, j), dTol);
			    sp_grad_assert_equal(Fsp3(i, j), F(i, j), dTol);
			    sp_grad_assert_equal(Fsp4(i, j), F(i, j), dTol);

			    sp_grad_assert_equal(Gsp1(i, j), G(i, j), dTol);
			    sp_grad_assert_equal(Gsp2(i, j), G(i, j), dTol);
			    sp_grad_assert_equal(Gsp3(i, j), G(i, j), dTol);
			    sp_grad_assert_equal(Gsp4(i, j), G(i, j), dTol);

			    sp_grad_assert_equal(Hsp1(i, j), H(i, j), dTol);
			    sp_grad_assert_equal(Hsp2(i, j), H(i, j), dTol);
			    sp_grad_assert_equal(Hsp3(i, j), H(i, j), dTol);
			    sp_grad_assert_equal(Hsp4(i, j), H(i, j), dTol);

			    sp_grad_assert_equal(Isp1(i, j), I(i, j), dTol);
			    sp_grad_assert_equal(Isp2(i, j), I(i, j), dTol);
			    sp_grad_assert_equal(Isp3(i, j), I(i, j), dTol);

			    sp_grad_assert_equal(Jsp1(i, j), J(i, j), dTol);
			    sp_grad_assert_equal(Jsp2(i, j), J(i, j), dTol);
			    sp_grad_assert_equal(Jsp3(i, j), J(i, j), dTol);
			    sp_grad_assert_equal(Jsp4(i, j), J(i, j), dTol);
			    sp_grad_assert_equal(Jsp5(i, j), J(i, j), dTol);
			    sp_grad_assert_equal(Jsp6(i, j), J(i, j), dTol);
			    sp_grad_assert_equal(Jsp7(i, j), J(i, j), dTol);

			    sp_grad_assert_equal(Ksp1(i, j), K(i, j), dTol);
			    sp_grad_assert_equal(Ksp2(i, j), K(i, j), dTol);

			    sp_grad_assert_equal(Lsp1(i, j), L(i, j), dTol);

			    sp_grad_assert_equal(RDeltasp(i, j), RDelta(i, j), dTol);
			    sp_grad_assert_equal(RVecsp(i, j), RVec(i, j), dTol);
			    sp_grad_assert_equal(Asp4(i, j), A(i, j), dTol);
			    sp_grad_assert_equal(Asp4(i, j), A4(i, j), dTol);
			    sp_grad_assert_equal(X3a(i, j), X3b(i, j), dTol);
			    sp_grad_assert_equal(X4a(i, j), X4b(i, j), dTol);
			    sp_grad_assert_equal(X5b(i, j), X5a(i, j), dTol);
			    sp_grad_assert_equal(X6a(i, j), X5a(i, j), dTol);
			    sp_grad_assert_equal(X6b(i, j), X5a(i, j), dTol);
			    sp_grad_assert_equal(X7a(i, j), X5a(i, j), dTol);
			    sp_grad_assert_equal(X7b(i, j), X5a(i, j), dTol);
			    sp_grad_assert_equal(X10a(i, j), X5a(i, j), dTol);
			    sp_grad_assert_equal(X10b(i, j), X5a(i, j), dTol);
			    sp_grad_assert_equal(X11a(i, j), X5a(i, j), dTol);
			    sp_grad_assert_equal(X11b(i, j), X5a(i, j), dTol);
			    sp_grad_assert_equal(Asp5(i, j), A(i, j), dTol);
			    sp_grad_assert_equal(X15b(i, j), X15a(i, j), dTol);			    
		       }
		       sp_grad_assert_equal(X8a(i), X5a(i, 3), dTol);
		       sp_grad_assert_equal(X8b(i), X5a(i, 3), dTol);
                       sp_grad_assert_equal(X9a(i), X5a(i, 3), dTol);
                       sp_grad_assert_equal(X9b(i), X5a(i, 3), dTol);
		       sp_grad_assert_equal(X14b(i), X14a(i), dTol);
		       sp_grad_assert_equal(X14a(i), X13a(i, 1) + X13a(i, 2), dTol);
		  }

		  for (index_type i = 1; i <= X12a.iGetNumRows(); ++i) {
		       for (index_type j = 1; j <= X12a.iGetNumCols(); ++j) {
			    sp_grad_assert_equal(X12b(i, j), X12a(i, j), dTol);
			    sp_grad_assert_equal(X12c(i, j), X12a(i, j), dTol);
			    sp_grad_assert_equal(X12d(j, i), X12a(i, j), dTol);
		       }
		  }
	     }


	}

	void test1()
	{
		using namespace std;
		using namespace std::chrono;

		cerr << __PRETTY_FUNCTION__ << ":\n";
		vector<SpDerivRec> v{{3, 3.2}, {4, 4.2}, {5, 5.2}, {6, 6.2}, {12, 12.2}, {13, 13.2}, {15, 15.2}, {14, 14.2}};
		SpGradient g1(10., {{13, 13.1}, {1, 1.1}, {2, 2.1}, {3, 3.1}, {5, 5.1}, {6, 6.1}, {10, 10.1}, {12, 12.1}});
		SpGradient g2, g3(20., v), g4, g5;

		g2.Reset(20., {{3, 3.2}, {4, 4.2}, {5, 5.2}, {6, 6.2}, {12, 12.2}, {13, 13.2}, {15, 15.2}, {14, 14.2}});
		g4.Reset(g2.dGetValue(), 3, g2.dGetDeriv(3));
#ifdef SP_GRAD_DEBUG
		cout << "g1=" << g1 << endl;
		cout << "g2=" << g2 << endl;
		cout << "g3=" << g3 << endl;
		cout << "g4=" << g4 << endl;
#endif
		SpGradient r1 = pow(2. * (g1 + g2), 1e-2 * g1);
		SpGradient r2 = -cos(r1) * sin(g2 - g1 / (1. - sqrt(pow(g1, 2) + pow(g2, 2))));
		SpGradient r3 = r2, r4, r5, r6;

		SpGradient a{1.0, {{1, 10.}}};
		SpGradient b{2.0, {{1, 20.}}};
		auto f = sqrt(a + b) / (a - b);

#ifdef SP_GRAD_DEBUG
		cout << f << endl;
#endif
		r4 = r1;
		r5 = r4;
		r6 = move(r4);

		SpGradient r7(move(r5));
		SpGradient r8(r7);

		r1 = r8;

		cout << fixed << setprecision(3);
#ifdef SP_GRAD_DEBUG
		cout << "r=" << r1 << endl;

		cout << "g1=" << g1 << endl;
		cout << "g2=" << g2 << endl;
#endif
		constexpr int w = 8;

		cout << "|" << setw(w) << left << "idof"
		     << "|" << setw(w) << left << "g1"
		     << "|" << setw(w) << left << "g2"
		     << "|" << setw(w) << left << "r1"
		     << "|" << setw(w) << left << "r2"
		     << "|" << endl;

		cout << "+" << setw(w) << setfill('-') << '-'
		     << "+" << setw(w) << setfill('-') << '-'
		     << "+" << setw(w) << setfill('-') << '-'
		     << "+" << setw(w) << setfill('-') << '-'
		     << "+" << setw(w) << setfill('-') << '-'
		     << "+" << endl << setfill(' ');

		cout << "|" << setw(w) << right << "#"
		     << "|" << setw(w) << right << g1.dGetValue()
		     << "|" << setw(w) << right << g2.dGetValue()
		     << "|" << setw(w) << right << r1.dGetValue()
		     << "|" << setw(w) << right << r2.dGetValue()
		     << "|" << endl;

		cout << "+" << setw(w) << setfill('-') << '-'
		     << "+" << setw(w) << setfill('-') << '-'
		     << "+" << setw(w) << setfill('-') << '-'
		     << "+" << setw(w) << setfill('-') << '-'
		     << "+" << setw(w) << setfill('-') << '-'
		     << "+" << endl << setfill(' ');

		for (index_type idof = 1; idof <= 15; ++idof) {
			cout << "|" << setw(w) << right << idof
			     << "|" << setw(w) << right << g1.dGetDeriv(idof)
			     << "|" << setw(w) << right << g2.dGetDeriv(idof)
			     << "|" << setw(w) << right << r1.dGetDeriv(idof)
			     << "|" << setw(w) << right << r2.dGetDeriv(idof)
			     << "|" << endl;
		}
	}

	void test2(index_type inumloops, index_type inumnz, index_type inumdof)
	{
		using namespace std;
		using namespace std::chrono;

		cerr << __PRETTY_FUNCTION__ << ":\n";

		random_device rd;
		mt19937 gen(rd());
		uniform_real_distribution<doublereal> randval(-1., 1.);
		uniform_int_distribution<index_type> randdof(1, inumdof);
		uniform_int_distribution<index_type> randnz(0, inumnz - 1);

		gen.seed(0);

		SpGradient u, v, w, f, fc, f2;

#ifdef USE_AUTODIFF
		grad::Gradient<0> us, vs, ws, fs;
		grad::LocalDofMap oDofMap;
		index_type fsnz = 0;
#endif

		doublereal e = randval(gen), fVal;
		index_type unz = 0, vnz = 0, wnz = 0;
		index_type fnz = 0, f2nz = 0, funz = 0, f2unz = 0, fcnz = 0;
		constexpr doublereal dTol = pow(numeric_limits<doublereal>::epsilon(), 0.5);
		vector<doublereal> ud, vd, wd, fd, work;
		duration<long long, ratio<1L, 1000000000L> > sp_grad_time(0), sp_grad_compr_time(0), sp_grad2_time(0), c_full_time(0);

		decltype(sp_grad_time) sp_grad_s_time(0);

		for (index_type iloop = 0; iloop < inumloops; ++iloop) {
			sp_grad_rand_gen(u, randnz, randdof, randval, gen);
			sp_grad_rand_gen(v, randnz, randdof, randval, gen);
			sp_grad_rand_gen(w, randnz, randdof, randval, gen);

			unz += u.iGetSize();
			vnz += v.iGetSize();
			wnz += w.iGetSize();

			SpGradDofStat s;

			u.GetDofStat(s);
			v.GetDofStat(s);
			w.GetDofStat(s);

			const index_type nbdirs = s.iNumNz ? s.iMaxDof : 0;
#ifdef USE_AUTODIFF
			oDofMap.Reset();
			us.SetValue(u.dGetValue());
			vs.SetValue(v.dGetValue());
			ws.SetValue(w.dGetValue());

			if (s.iNumNz) {
				us.DerivativeResizeReset(&oDofMap, s.iMinDof, s.iMaxDof + 1, grad::MapVectorBase::GLOBAL, 0.);
				vs.DerivativeResizeReset(&oDofMap, s.iMinDof, s.iMaxDof + 1, grad::MapVectorBase::GLOBAL, 0.);
				ws.DerivativeResizeReset(&oDofMap, s.iMinDof, s.iMaxDof + 1, grad::MapVectorBase::GLOBAL, 0.);

				for (index_type i = s.iMinDof; i <= s.iMaxDof; ++i) {
					us.SetDerivativeGlobal(i, u.dGetDeriv(i));
					vs.SetDerivativeGlobal(i, v.dGetDeriv(i));
					ws.SetDerivativeGlobal(i, w.dGetDeriv(i));
				}
			}
#endif

			func_scalar1(u.dGetValue(), v.dGetValue(), w.dGetValue(), e, fVal);

			auto sp_grad_start = high_resolution_clock::now();

			func_scalar1(u, v, w, e, f);

			sp_grad_time += high_resolution_clock::now() - sp_grad_start;

			funz += f.iGetSize();

			f.Sort();

			fnz += f.iGetSize();
#ifdef USE_AUTODIFF
			auto sp_grad_s_start = high_resolution_clock::now();

			func_scalar1(us, vs, ws, e, fs);

			sp_grad_s_time += high_resolution_clock::now() - sp_grad_s_start;

			fsnz += fs.iGetLocalSize();
#endif
			auto sp_grad_compr_start = high_resolution_clock::now();

			func_scalar1_compressed(u, v, w, e, fc);

			sp_grad_compr_time += high_resolution_clock::now() - sp_grad_compr_start;

			fcnz += fc.iGetSize();

			auto sp_grad2_start = high_resolution_clock::now();

			func_scalar2(u, v, w, e, f2);

			sp_grad2_time += high_resolution_clock::now() - sp_grad2_start;

			f2unz += f2.iGetSize();

			f2.Sort();

			f2nz += f2.iGetSize();

			assert(fabs(f.dGetValue() / fVal - 1.) < dTol);

			assert(fabs(fc.dGetValue() / fVal - 1.) < dTol);

			assert(fabs(f2.dGetValue() / fVal - 1.) < dTol);

			ud.clear();
			vd.clear();
			wd.clear();
			fd.clear();
			work.clear();
			ud.resize(nbdirs);
			vd.resize(nbdirs);
			wd.resize(nbdirs);
			fd.resize(nbdirs);
			work.resize(4 * nbdirs);

			for (index_type i = 1; i <= s.iMaxDof; ++i) {
				ud[i - 1] = u.dGetDeriv(i);
				vd[i - 1] = v.dGetDeriv(i);
				wd[i - 1] = w.dGetDeriv(i);
			}

			auto c_full_start = high_resolution_clock::now();

			func_scalar1_dv(nbdirs,
					u.dGetValue(),
					front(ud),
					v.dGetValue(),
					front(vd),
					w.dGetValue(),
					front(wd),
					e,
					fVal,
					front(fd),
					front(work));

			c_full_time += high_resolution_clock::now() - c_full_start;

			assert(fabs(f.dGetValue() - fVal) < dTol * max(1., fabs(fVal)));

			SP_GRAD_TRACE("fref f\n");
			SP_GRAD_TRACE(fVal << " " << f.dGetValue() << endl);

			for (index_type i = 1; i <= s.iMaxDof; ++i) {
				assert(fabs(fd[i - 1] - f.dGetDeriv(i)) < dTol * max(1.0, fabs(fd[i - 1])));
				assert(fabs(fd[i - 1] - fc.dGetDeriv(i)) < dTol * max(1.0, fabs(fd[i - 1])));
				assert(fabs(fd[i - 1] - f2.dGetDeriv(i)) < dTol * max(1.0, fabs(fd[i - 1])));
				SP_GRAD_TRACE(fd[i - 1] << " " << f.dGetDeriv(i) << endl);
			}
#ifdef USE_AUTODIFF
			for (index_type i = fs.iGetStartIndexLocal(); i < fs.iGetEndIndexLocal(); ++i) {
				index_type iDof = fs.pGetDofMap()->iGetGlobalDof(i);
				assert(fabs(fd[iDof - 1] - fs.dGetDerivativeLocal(i)) < dTol * max(1.0, fabs(fd[iDof - 1])));
			}
#endif
		}

#ifdef USE_AUTODIFF
		auto sp_grad_s_time_ns = duration_cast<nanoseconds>(sp_grad_s_time).count();
#endif
		auto sp_grad_time_ns = duration_cast<nanoseconds>(sp_grad_time).count();
		auto sp_grad_compr_time_ns = duration_cast<nanoseconds>(sp_grad_compr_time).count();
		auto sp_grad2_time_ns = duration_cast<nanoseconds>(sp_grad2_time).count();
		auto c_full_time_ns = duration_cast<nanoseconds>(c_full_time).count();

		cerr << "test2: test passed with tolerance "
		     << scientific << setprecision(6)
		     << dTol
		     << " and nz=" << inumnz
		     << " unz=" << fixed << setprecision(0) << ceil(static_cast<doublereal>(unz) / inumloops)
		     << " vnz=" << fixed << setprecision(0) << ceil(static_cast<doublereal>(vnz) / inumloops)
		     << " wnz=" << fixed << setprecision(0) << ceil(static_cast<doublereal>(wnz) / inumloops)
		     << endl;
#ifdef USE_AUTODIFF
		cerr << "test2: sp_grad_s_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_s_time_ns) / 1e9
		     << " nz=" << fixed << setprecision(0) << ceil(static_cast<doublereal>(fsnz) / inumloops)
		     << endl;
#endif
		cerr << "test2: sp_grad_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_time_ns) / 1e9
		     << " nz=" << fixed << setprecision(0) << ceil(static_cast<doublereal>(fnz) / inumloops)
		     << " unz=" << fixed << setprecision(0) << ceil(static_cast<doublereal>(funz) / inumloops)
		     << endl;
		cerr << "test2: sp_grad_compressed_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_compr_time_ns) / 1e9
		     << " nz=" << fixed << setprecision(0) << ceil(static_cast<doublereal>(fcnz) / inumloops)
		     << endl;

		cerr << "test2: sp_grad2_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad2_time_ns) / 1e9
		     << " nz=" << fixed << setprecision(0) << ceil(static_cast<doublereal>(f2nz) / inumloops)
		     << " unz=" << fixed << setprecision(0) << ceil(static_cast<doublereal>(f2unz) / inumloops)
		     << endl;
		cerr << "test2: c_full_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(c_full_time_ns) / 1e9 << endl;
#ifdef USE_AUTODIFF
		cerr << "test2: sp_grad_s_time / c_full_time = "
		     << static_cast<doublereal>(sp_grad_s_time_ns) / max<int64_t>(1L, c_full_time_ns) // Avoid division by zero
		     << endl;
#endif
		cerr << "test2: sp_grad_time / c_full_time = "
		     << static_cast<doublereal>(sp_grad_time_ns) / max<int64_t>(1L, c_full_time_ns) // Avoid division by zero
		     << endl;
		cerr << "test2: sp_grad_compr_time / c_full_time = "
		     << static_cast<doublereal>(sp_grad_compr_time_ns) / max<int64_t>(1L, c_full_time_ns) // Avoid division by zero
		     << endl;
		cerr << "test2: sp_grad2_time / c_full_time = "
		     << static_cast<doublereal>(sp_grad2_time_ns) / max<int64_t>(1L, c_full_time_ns) // Avoid division by zero
		     << endl;
	}

	void test_bool1(index_type inumloops, index_type inumnz, index_type inumdof)
	{
		using namespace std;
		using namespace std::chrono;

		cerr << __PRETTY_FUNCTION__ << ":\n";

		random_device rd;
		mt19937 gen(rd());
		uniform_real_distribution<doublereal> randval(-1., 1.);
		uniform_int_distribution<index_type> randdof(1, inumdof);
		uniform_int_distribution<index_type> randnz(0, inumnz - 1);

		gen.seed(0);

		SpGradient u, v, w;
		doublereal e = randval(gen);

		for (index_type iloop = 0; iloop < inumloops; ++iloop) {
			sp_grad_rand_gen(u, randnz, randdof, randval, gen);
			sp_grad_rand_gen(v, randnz, randdof, randval, gen);
			sp_grad_rand_gen(w, randnz, randdof, randval, gen);

#define SP_GRAD_BOOL_TEST_FUNC(func)                                    \
			{                                               \
				bool f1 = func(u, v, w, e);             \
				bool f2 = func(u.dGetValue(), v.dGetValue(), w.dGetValue(), e); \
				bool f3 = func(u, v.dGetValue(), w.dGetValue(), e); \
				bool f4 = func(u.dGetValue(), v, w.dGetValue(), e); \
				bool f5 = func(u.dGetValue(), v.dGetValue(), w, e); \
				bool f6 = func(u, v, w.dGetValue(), e); \
				bool f7 = func(u.dGetValue(), v, w, e); \
				bool f8 = func(u, v.dGetValue(), w, e); \
									\
				assert(f1 == f2);                       \
				assert(f1 == f3);                       \
				assert(f1 == f4);                       \
				assert(f1 == f5);                       \
				assert(f1 == f6);                       \
				assert(f1 == f7);                       \
				assert(f1 == f8);                       \
			}

			SP_GRAD_BOOL_TEST_FUNC(func_bool1);
			SP_GRAD_BOOL_TEST_FUNC(func_bool2);
			SP_GRAD_BOOL_TEST_FUNC(func_bool3);
			SP_GRAD_BOOL_TEST_FUNC(func_bool4);
			SP_GRAD_BOOL_TEST_FUNC(func_bool5);
			SP_GRAD_BOOL_TEST_FUNC(func_bool6);
		}
	}

	template <typename TA, typename TB>
	void test3(const index_type inumloops, const index_type inumnz, const index_type inumdof, const index_type imatsize)
	{
		using namespace std;
		using namespace std::chrono;

		cerr << __PRETTY_FUNCTION__ << ":\n";

		random_device rd;
		mt19937 gen(rd());
		uniform_real_distribution<doublereal> randval(-1., 1.);
		uniform_int_distribution<index_type> randdof(1, inumdof);
		uniform_int_distribution<index_type> randnz(0, inumnz - 1);

		gen.seed(0);

		constexpr doublereal dTol = pow(numeric_limits<doublereal>::epsilon(), 0.5);

		vector<TA> A;
		vector<TB> B;
		vector<doublereal> Av(imatsize), Bv(imatsize), Ad, Bd, fd;
		index_type fnz = 0, f2nz = 0, funz = 0, f2unz = 0, f3nz = 0, f3unz = 0;
		index_type anz = 0, bnz = 0;

		A.resize(imatsize);
		B.resize(imatsize);

		duration<long long, ratio<1L, 1000000000L> >
			sp_grad_time(0),
			sp_grad2_time(0),
			sp_grad3_time(0),
			c_full_time(0);

		typename util::ResultType<TA, TB>::Type f, f2, f3;

		for (index_type iloop = 0; iloop < inumloops; ++iloop) {
			for (auto& a:A) {
				anz += sp_grad_rand_gen(a, randnz, randdof, randval, gen);
			}

			for (auto& b:B) {
				bnz += sp_grad_rand_gen(b, randnz, randdof, randval, gen);
			}

			auto sp_grad_start = high_resolution_clock::now();

			func_mat_mul1(f, A, B);

			sp_grad_time += high_resolution_clock::now() - sp_grad_start;

			funz += SpGradient::iGetSize(f);

			SpGradient::Sort(f);

			fnz += SpGradient::iGetSize(f);

			auto sp_grad2_start = high_resolution_clock::now();

			func_mat_mul2(f2, A, B);

			sp_grad2_time += high_resolution_clock::now() - sp_grad2_start;

			f2unz += SpGradient::iGetSize(f2);

			SpGradient::Sort(f2);

			f2nz += SpGradient::iGetSize(f2);

			auto sp_grad3_start = high_resolution_clock::now();

			func_mat_mul3(f3, A, B);

			sp_grad3_time += high_resolution_clock::now() - sp_grad3_start;

			f3unz += SpGradient::iGetSize(f3);

			SpGradient::Sort(f3);

			f3nz += SpGradient::iGetSize(f3);

			SpGradDofStat s;

			for (const auto& a: A) {
				SpGradient::GetDofStat(a, s);
			}

			for (const auto& b: B) {
				SpGradient::GetDofStat(b, s);
			}

			const index_type nbdirs = s.iNumNz ? s.iMaxDof : 0;
			doublereal fVal;

			Ad.resize(nbdirs * imatsize);
			Bd.resize(nbdirs * imatsize);
			fd.resize(nbdirs);

			for (index_type i = 0; i < imatsize; ++i) {
				Av[i] = SpGradient::dGetValue(A[i]);
				Bv[i] = SpGradient::dGetValue(B[i]);

				for (index_type j = 1; j <= s.iMaxDof; ++j) {
					Ad[i + (j - 1) * imatsize] = SpGradient::dGetDeriv(A[i], j);
					Bd[i + (j - 1) * imatsize] = SpGradient::dGetDeriv(B[i], j);
				}
			}

			auto c_full_start = high_resolution_clock::now();

			func_mat_mul1_dv(imatsize,
					 front(Av),
					 front(Ad),
					 front(Bv),
					 front(Bd),
					 fVal,
					 front(fd),
					 s.iMaxDof);

			c_full_time += high_resolution_clock::now() - c_full_start;

			assert(fabs(SpGradient::dGetValue(f) - fVal) < dTol * max(1., fabs(fVal)));
			assert(fabs(SpGradient::dGetValue(f2) - fVal) < dTol * max(1., fabs(fVal)));
			assert(fabs(SpGradient::dGetValue(f3) - fVal) < dTol * max(1., fabs(fVal)));

			for (index_type i = 1; i <= s.iMaxDof; ++i) {
				assert(fabs(fd[i - 1] - SpGradient::dGetDeriv(f, i)) < dTol * max(1.0, fabs(fd[i - 1])));
			}

			assert(fabs(SpGradient::dGetValue(f2) - fVal) < dTol * max(1., fabs(fVal)));

			for (index_type i = 1; i <= s.iMaxDof; ++i) {
				assert(fabs(fd[i - 1] - SpGradient::dGetDeriv(f2, i)) < dTol * max(1.0, fabs(fd[i - 1])));
			}

			assert(fabs(SpGradient::dGetValue(f3) - fVal) < dTol * max(1., fabs(fVal)));

			for (index_type i = 1; i <= s.iMaxDof; ++i) {
				assert(fabs(fd[i - 1] - SpGradient::dGetDeriv(f3, i)) < dTol * max(1.0, fabs(fd[i - 1])));
			}
		}

		cerr << "test3: test passed with tolerance "
		     << scientific << dTol
		     << " and nz=" << inumnz
		     << " anz=" << fixed << setprecision(0) << ceil(static_cast<doublereal>(anz) / (inumloops * imatsize))
		     << " bnz=" << fixed << setprecision(0) << ceil(static_cast<doublereal>(bnz) / (inumloops * imatsize))
		     << endl;

		auto sp_grad_time_ns = duration_cast<nanoseconds>(sp_grad_time).count();
		auto sp_grad2_time_ns = duration_cast<nanoseconds>(sp_grad2_time).count();
		auto sp_grad3_time_ns = duration_cast<nanoseconds>(sp_grad3_time).count();
		auto c_full_time_ns = duration_cast<nanoseconds>(c_full_time).count();

		cerr << "test3: sp_grad_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_time_ns) / 1e9
		     << " nz=" << fixed << setprecision(0) << ceil(static_cast<doublereal>(fnz) / inumloops)
		     << " unz=" << fixed << setprecision(0) << ceil(static_cast<doublereal>(funz) / inumloops) << endl;
		cerr << "test3: sp_grad2_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad2_time_ns) / 1e9
		     << " nz=" << fixed << setprecision(0) << ceil(static_cast<doublereal>(f2nz) / inumloops)
		     << " unz=" << fixed << setprecision(0) << ceil(static_cast<doublereal>(f2unz) / inumloops) << endl;
		cerr << "test3: sp_grad3_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad3_time_ns) / 1e9
		     << " nz=" << fixed << setprecision(0) << ceil(static_cast<doublereal>(f3nz) / inumloops)
		     << " unz=" << fixed << setprecision(0) << ceil(static_cast<doublereal>(f3unz) / inumloops) << endl;
		cerr << "test3: c_full_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(c_full_time_ns) / 1e9 << endl;
		cerr << "test3: sp_grad_time / c_full_time = "
		     << static_cast<doublereal>(sp_grad_time_ns) / max<int64_t>(1L, c_full_time_ns) // Avoid division by zero
		     << endl;
		cerr << "test3: sp_grad2_time / c_full_time = "
		     << static_cast<doublereal>(sp_grad2_time_ns) / max<int64_t>(1L, c_full_time_ns) // Avoid division by zero
		     << endl;
		cerr << "test3: sp_grad3_time / c_full_time = "
		     << static_cast<doublereal>(sp_grad3_time_ns) / max<int64_t>(1L, c_full_time_ns) // Avoid division by zero
		     << endl;
	}
#ifdef USE_AUTODIFF
	template <typename T>
	struct SpGradient2SliceGradient;

	template <>
	struct SpGradient2SliceGradient<sp_grad::SpGradient> {
		typedef grad::Gradient<0> Type;

		static void Copy(const sp_grad::SpGradient& g1, Type& g2, grad::LocalDofMap& oDofMap) {
			sp_grad::SpGradDofStat s;
			g1.GetDofStat(s);

			g2.SetValue(g1.dGetValue());

			if (s.iNumNz) {
				g2.DerivativeResizeReset(&oDofMap, s.iMinDof, s.iMaxDof + 1, grad::MapVectorBase::GLOBAL, 0.);

				for (const auto& d: g1) {
					g2.SetDerivativeGlobal(d.iDof, d.dDer);
				}
			}
		}

		static doublereal dGetDerivativeGlobal(const grad::Gradient<0>& g, index_type iDof) {
			if (g.pGetDofMap()->iGetLocalIndex(iDof) >= 0) {
				return g.dGetDerivativeGlobal(iDof);
			} else {
				return 0.;
			}
		}
	};

	template <>
	struct SpGradient2SliceGradient<doublereal> {
		typedef doublereal Type;

		static void Copy(doublereal g1, doublereal& g2, grad::LocalDofMap&) {
			g2 = g1;
		}

		static doublereal dGetDerivativeGlobal(doublereal, index_type) {
			return 0.;
		}
	};
#endif
	template <typename TA, typename TX>
	void test4(const index_type inumloops,
		   const index_type inumnz,
		   const index_type inumdof,
		   const index_type imatrows,
		   const index_type imatcols)
	{
		using namespace std;
		using namespace std::chrono;

		cerr << __PRETTY_FUNCTION__ << ":\n";
		random_device rd;
		mt19937 gen(rd());
		uniform_real_distribution<doublereal> randval(-1., 1.);
		uniform_int_distribution<index_type> randdof(1, inumdof);
		uniform_int_distribution<index_type> randnz(0, inumnz - 1);

		gen.seed(0);

		constexpr doublereal dTol = pow(numeric_limits<doublereal>::epsilon(), 0.5);

		typedef typename util::ResultType<TA,TX>::Type TB;
		vector<TA> A;
		vector<TX> x;
		vector<TB> b4, b5, b6, b7;
#ifdef USE_AUTODIFF
		typedef typename SpGradient2SliceGradient<TA>::Type TAS;
		typedef typename SpGradient2SliceGradient<TX>::Type TXS;
		typedef typename util::ResultType<TAS, TXS>::Type TBS;
		vector<TAS> As;
		vector<TXS> xs;
		vector<TBS> bs;
		grad::LocalDofMap oDofMap;
		As.resize(imatrows * imatcols);
		xs.resize(imatcols);
		bs.resize(imatrows);
#endif
		vector<doublereal> Av, Ad, xv, bv, xd, bd;
		index_type xnz = 0, anz = 0, b4nz = 0, b5unz = 0, b5nz = 0;

		A.resize(imatrows * imatcols);
		Av.resize(imatrows * imatcols);
		x.resize(imatcols);
		b4.resize(imatrows);
		b5.resize(imatrows);
		b6.resize(imatrows);
		b7.resize(imatrows);
		xv.resize(imatcols);
		bv.resize(imatrows);

		duration<long long, ratio<1L, 1000000000L> >
			sp_grad4_time(0),
			sp_grad4m_time(0),
			sp_grad4sp_time(0),
			sp_grad4sm_time(0),
#ifdef USE_AUTODIFF
			sp_grad4s_time(0),
#endif
			c_full_time(0);

		for (index_type iloop = 0; iloop < inumloops; ++iloop) {
			for (auto& ai: A) {
				sp_grad_rand_gen(ai, randnz, randdof, randval, gen);
				anz += SpGradient::iGetSize(ai);
			}

			for (auto& xi: x) {
				sp_grad_rand_gen(xi, randnz, randdof, randval, gen);
				xnz += SpGradient::iGetSize(xi);
			}
#ifdef USE_AUTODIFF
			for (index_type i = 0; i < imatrows * imatcols; ++i) {
				SpGradient2SliceGradient<TA>::Copy(A[i], As[i], oDofMap);
			}

			for (index_type i = 0; i < imatcols; ++i) {
				SpGradient2SliceGradient<TX>::Copy(x[i], xs[i], oDofMap);
			}
#endif
			auto sp_grad4_start = high_resolution_clock::now();

			func_mat_mul4(A, x, b4);

			sp_grad4_time += high_resolution_clock::now() - sp_grad4_start;

			for (const auto& b4i: b4) {
				b4nz += SpGradient::iGetSize(b4i);
			}

			auto sp_grad4m_start = high_resolution_clock::now();

			func_mat_mul4m(A, x, b5);

			sp_grad4m_time += high_resolution_clock::now() - sp_grad4m_start;

			auto sp_grad4sp_start = high_resolution_clock::now();

			func_mat_mul4s(A, x, b6);

			sp_grad4sp_time += high_resolution_clock::now() - sp_grad4sp_start;

			auto sp_grad4sm_start = high_resolution_clock::now();

			func_mat_mul4sm(A, x, b7);

			sp_grad4sm_time += high_resolution_clock::now() - sp_grad4sm_start;
#ifdef USE_AUTODIFF
			auto sp_grad4s_start = high_resolution_clock::now();

			func_mat_mul4s(As, xs, bs);

			sp_grad4s_time += high_resolution_clock::now() - sp_grad4s_start;
#endif
			for (const auto& b5i: b5) {
				b5unz += SpGradient::iGetSize(b5i);
			}

			for (auto& b5i: b5) {
				SpGradient::Sort(b5i);
			}

			for (const auto& b5i: b5) {
				b5nz += SpGradient::iGetSize(b5i);
			}

			SpGradDofStat s;

			for (const auto& xi: x) {
				SpGradient::GetDofStat(xi, s);
			}

			const index_type nbdirs = s.iNumNz ? s.iMaxDof : 0;

			Ad.resize(imatrows * imatcols * nbdirs);
			xd.resize(nbdirs * imatcols);
			bd.resize(nbdirs * imatrows);

			for (index_type i = 0; i < imatrows * imatcols; ++i) {
				Av[i] = SpGradient::dGetValue(A[i]);

				for (index_type j = 1; j <= nbdirs; ++j) {
					Ad[i + (j - 1) * imatcols * imatrows] = SpGradient::dGetDeriv(A[i], j);
				}
			}

			for (index_type i = 0; i < imatcols; ++i) {
				xv[i] = SpGradient::dGetValue(x[i]);

				for (index_type j = 1; j <= nbdirs; ++j) {
					xd[i + (j - 1) * imatcols] = SpGradient::dGetDeriv(x[i], j);
				}
			}

			auto c_full_start = high_resolution_clock::now();

			func_mat_mul4<TA, TX>(imatrows,
					      imatcols,
					      nbdirs,
					      front(Av),
					      front(Ad),
					      front(xv),
					      front(xd),
					      front(bv),
					      front(bd));

			c_full_time += high_resolution_clock::now() - c_full_start;

			for (index_type i = 0; i < imatrows; ++i) {
				for (index_type j = 1; j <= nbdirs; ++j) {
					assert(fabs(bd[i + (j - 1) * imatrows] - SpGradient::dGetDeriv(b4[i], j)) / max(1., fabs(bd[i + (j - 1) * imatrows])) < dTol);
					assert(fabs(bd[i + (j - 1) * imatrows] - SpGradient::dGetDeriv(b5[i], j)) / max(1., fabs(bd[i + (j - 1) * imatrows])) < dTol);
					assert(fabs(bd[i + (j - 1) * imatrows] - SpGradient::dGetDeriv(b6[i], j)) / max(1., fabs(bd[i + (j - 1) * imatrows])) < dTol);
					assert(fabs(bd[i + (j - 1) * imatrows] - SpGradient::dGetDeriv(b7[i], j)) / max(1., fabs(bd[i + (j - 1) * imatrows])) < dTol);
#ifdef USE_AUTODIFF
					assert(fabs(bd[i + (j - 1) * imatrows] - SpGradient2SliceGradient<TB>::dGetDerivativeGlobal(bs[i], j)) / max(1., fabs(bd[i + (j - 1) * imatrows])) < dTol);
#endif
				}
			}

			for (index_type i = 0; i < imatrows; ++i) {
				assert(fabs(bv[i] - SpGradient::dGetValue(b4[i])) / max(1., fabs(bv[i])) < dTol);
				assert(fabs(bv[i] - SpGradient::dGetValue(b5[i])) / max(1., fabs(bv[i])) < dTol);
				assert(fabs(bv[i] - SpGradient::dGetValue(b6[i])) / max(1., fabs(bv[i])) < dTol);
				assert(fabs(bv[i] - SpGradient::dGetValue(b7[i])) / max(1., fabs(bv[i])) < dTol);
#ifdef USE_AUTODIFF
				assert(fabs(bv[i] - grad::dGetValue(bs[i])) / max(1., fabs(bv[i])) < dTol);
#endif
			}
		}

		auto sp_grad4_time_ns = duration_cast<nanoseconds>(sp_grad4_time).count();
		auto sp_grad4m_time_ns = duration_cast<nanoseconds>(sp_grad4m_time).count();
		auto sp_grad4sp_time_ns = duration_cast<nanoseconds>(sp_grad4sp_time).count();
		auto sp_grad4sm_time_ns = duration_cast<nanoseconds>(sp_grad4sm_time).count();
#ifdef USE_AUTODIFF
		auto sp_grad4s_time_ns = duration_cast<nanoseconds>(sp_grad4s_time).count();
#endif
		auto c_full_time_ns = duration_cast<nanoseconds>(c_full_time).count();

		cerr << "test4: sp_grad4_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad4_time_ns) / 1e9
		     << " bnz=" << fixed << setprecision(0) << ceil(static_cast<doublereal>(b4nz) / (imatrows * inumloops))
		     << " anz=" << fixed << setprecision(0) << ceil(static_cast<doublereal>(anz) / (imatcols * inumloops))
		     << " xnz=" << fixed << setprecision(0) << ceil(static_cast<doublereal>(xnz) / (imatcols * inumloops))
		     << endl;
		cerr << "test4: sp_grad4m_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad4m_time_ns) / 1e9
		     << " bnz=" << fixed << setprecision(0) << ceil(static_cast<doublereal>(b5nz) / (imatrows * inumloops))
		     << " bunz=" << fixed << setprecision(0) << ceil(static_cast<doublereal>(b5unz) / (imatrows * inumloops))
		     << endl;
		cerr << "test4: sp_grad4sp_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad4sp_time_ns) / 1e9
		     << endl;
		cerr << "test4: sp_grad4sm_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad4sm_time_ns) / 1e9
		     << endl;
#ifdef USE_AUTODIFF
		cerr << "test4: sp_grad4s_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad4s_time_ns) / 1e9
		     << endl;
#endif
		cerr << "test4: c_full_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(c_full_time_ns) / 1e9
		     << endl;
		cerr << "test4: sp_grad4_time / c_full_time = "
		     << static_cast<doublereal>(sp_grad4_time_ns) / max<int64_t>(1L, c_full_time_ns) // Avoid division by zero
		     << endl;
		cerr << "test4: sp_grad4m_time / c_full_time = "
		     << static_cast<doublereal>(sp_grad4m_time_ns) / max<int64_t>(1L, c_full_time_ns) // Avoid division by zero
		     << endl;
		cerr << "test4: sp_grad4sp_time / c_full_time = "
		     << static_cast<doublereal>(sp_grad4sp_time_ns) / max<int64_t>(1L, c_full_time_ns) // Avoid division by zero
		     << endl;
		cerr << "test4: sp_grad4sm_time / c_full_time = "
		     << static_cast<doublereal>(sp_grad4sm_time_ns) / max<int64_t>(1L, c_full_time_ns) // Avoid division by zero
		     << endl;
#ifdef USE_AUTODIFF
		cerr << "test4: sp_grad4s_time / c_full_time = "
		     << static_cast<doublereal>(sp_grad4s_time_ns) / max<int64_t>(1L, c_full_time_ns) // Avoid division by zero
		     << endl;
#endif
	}

	void test6(index_type inumloops, index_type inumnz, index_type inumdof, index_type imatrows, index_type imatcols)
	{
		using namespace std;
		using namespace std::chrono;

		cerr << __PRETTY_FUNCTION__ << ":\n";

		random_device rd;
		mt19937 gen(rd());
		uniform_real_distribution<doublereal> randval(-1., 1.);
		uniform_int_distribution<index_type> randdof(1, inumdof);
		uniform_int_distribution<index_type> randnz(0, inumnz - 1);
		uniform_int_distribution<index_type> randresizecnt(0, imatrows * imatcols - 1);
		uniform_int_distribution<index_type> randresizeidx(0, imatrows * imatcols - 1);

		gen.seed(0);

		std::vector<SpGradient> rgVec;

		for (index_type iloop = 0; iloop < std::min<index_type>(10, inumloops); ++iloop) {
			SpMatrixBase<SpGradient> A(imatrows, imatcols, inumnz);

			for (auto& g: A) {
				sp_grad_rand_gen(g, randnz, randdof, randval, gen);
			}


			SpMatrixBase<SpGradient> A2(A.iGetNumRows(), A.iGetNumCols(), inumnz);

			std::copy(A.begin(), A.end(), A2.begin());

			rgVec.clear();

			rgVec.reserve(A.iGetNumRows() * A.iGetNumCols() * 5);

			std::copy(A.begin(), A.end(), std::back_inserter(rgVec));

			SpMatrixBase<SpGradient> B(std::move(A));

			std::copy(B.begin(), B.end(), std::back_inserter(rgVec));

			SpMatrixBase<SpGradient> C, D, E;

			C = std::move(B);

			std::copy(C.begin(), C.end(), std::back_inserter(rgVec));

			D = std::move(C);

			std::copy(D.begin(), D.end(), std::back_inserter(rgVec));

			E = std::move(D);

			std::copy(E.begin(), E.end(), std::back_inserter(rgVec));

			index_type irescnt = randresizecnt(gen);

			for (index_type i = 0; i < irescnt; ++i) {
				rgVec[randresizeidx(gen)] += EvalUnique(rgVec[randresizeidx(gen)] * randval(gen) * rgVec[randresizeidx(gen)]);
				rgVec[randresizeidx(gen)] *= EvalUnique(rgVec[randresizeidx(gen)] * randval(gen) - rgVec[randresizeidx(gen)] * randval(gen));
				rgVec[randresizeidx(gen)] -= EvalUnique(rgVec[randresizeidx(gen)] * randval(gen) / (rgVec[randresizeidx(gen)] * randval(gen) + randval(gen)));
				rgVec[randresizeidx(gen)] /= EvalUnique((rgVec[randresizeidx(gen)] * randval(gen) + randval(gen)));

				rgVec[randresizeidx(gen)] *= randval(gen);
				rgVec[randresizeidx(gen)] /= randval(gen);
				rgVec[randresizeidx(gen)] += randval(gen);
				rgVec[randresizeidx(gen)] -= randval(gen);
				rgVec[randresizeidx(gen)] = rgVec[randresizeidx(gen)];
				rgVec[randresizeidx(gen)] = EvalUnique(rgVec[randresizeidx(gen)] * (1 + randval(gen)) / (2. + randval(gen) + pow(rgVec[randresizeidx(gen)], 2)));
			}

			for (auto& vi: rgVec) {
				vi.Sort();
			}

			A = std::move(E);

			assert(A2.iGetNumRows() == A.iGetNumRows());
			assert(A2.iGetNumCols() == A.iGetNumCols());

			for (index_type i = 0; i < A.iGetNumRows() * A.iGetNumCols(); ++i) {
				const SpGradient& ai = *(A.begin() + i);
				const SpGradient& a2i = *(A2.begin() + i);

				assert(ai.dGetValue() == a2i.dGetValue());

				for (index_type j = 0; j < inumdof; ++j) {
					assert(ai.dGetDeriv(j) == a2i.dGetDeriv(j));
				}
			}

			for (index_type i = 0; i < irescnt; ++i) {
				*(A.begin() + randresizeidx(gen)) += *(A.begin() + randresizeidx(gen)) * randval(gen) * *(A.begin() + randresizeidx(gen));
				*(A.begin() + randresizeidx(gen)) *= *(A.begin() + randresizeidx(gen)) * randval(gen) - *(A.begin() + randresizeidx(gen)) * randval(gen);
				*(A.begin() + randresizeidx(gen)) -= *(A.begin() + randresizeidx(gen)) * randval(gen) / (*(A.begin() + randresizeidx(gen)) * randval(gen) + randval(gen));
				*(A.begin() + randresizeidx(gen)) /= (*(A.begin() + randresizeidx(gen)) * randval(gen) + randval(gen));

				*(A.begin() + randresizeidx(gen)) *= randval(gen);
				*(A.begin() + randresizeidx(gen)) /= randval(gen);
				*(A.begin() + randresizeidx(gen)) += randval(gen);
				*(A.begin() + randresizeidx(gen)) -= randval(gen);
				*(A.begin() + randresizeidx(gen)) = *(A.begin() + randresizeidx(gen));
				*(A.begin() + randresizeidx(gen)) = EvalUnique(*(A.begin() + randresizeidx(gen)) * (1 + randval(gen)) / (2. + randval(gen) + pow(*(A.begin() + randresizeidx(gen)), 2)));

				*(A.begin() + randresizeidx(gen)) += EvalUnique(rgVec[randresizeidx(gen)] * randval(gen) * *(A.begin() + randresizeidx(gen)));
				rgVec[randresizeidx(gen)] *= *(A.begin() + randresizeidx(gen)) * randval(gen) - rgVec[randresizeidx(gen)] * randval(gen);
				*(A.begin() + randresizeidx(gen)) -= rgVec[randresizeidx(gen)] * randval(gen) / (*(A.begin() + randresizeidx(gen)) * randval(gen) + randval(gen));
				*(A.begin() + randresizeidx(gen)) /= (rgVec[randresizeidx(gen)] * randval(gen) + randval(gen));

				*(A.begin() + randresizeidx(gen)) *= randval(gen);
				*(A.begin() + randresizeidx(gen)) = rgVec[randresizeidx(gen)];
				*(A.begin() + randresizeidx(gen)) = rgVec[randresizeidx(gen)] * (1 + randval(gen)) / (2. + randval(gen) + pow(rgVec[randresizeidx(gen)], 2));
			}
		}
	}

	template <typename TA, typename TB, typename TC, index_type NumRows = SpMatrixSize::DYNAMIC, index_type NumCols = SpMatrixSize::DYNAMIC>
	void test7(const index_type inumloops,
		   const index_type inumnz,
		   const index_type inumdof,
		   const index_type imatrows,
		   const index_type imatcols)
	{
		using namespace std;
		using namespace std::chrono;

		cerr << __PRETTY_FUNCTION__ << ":\n";
		random_device rd;
		mt19937 gen(rd());
		uniform_real_distribution<doublereal> randval(-1., 1.);
		uniform_int_distribution<index_type> randdof(1, inumdof);
		uniform_int_distribution<index_type> randnz(0, inumnz - 1);

		gen.seed(0);

		constexpr doublereal dTol = pow(numeric_limits<doublereal>::epsilon(), 0.5);

		typedef typename util::ResultType<typename util::ResultType<TA, TB>::Type, TC>::Type TD;

		SpMatrixBase<TA, NumRows, NumCols> A(imatrows, imatcols, inumnz);
		SpMatrixBase<TB, NumRows, NumCols> B(imatrows, imatcols, inumnz);
		SpMatrixBase<TC, NumRows, NumCols> C(imatrows, imatcols, inumnz);
		SpMatrixBase<TD, NumRows, NumCols> D(imatrows, imatcols, 3 * inumnz);
		SpMatrixBase<TD, NumRows, NumCols> Da(imatrows, imatcols, 3 * inumnz);
		vector<doublereal> Av(imatrows * imatcols), Ad,
			Bv(imatrows * imatcols), Bd,
			Cv(imatrows * imatcols), Cd,
			Dv(imatrows * imatcols), Dd;

		Ad.reserve(imatrows * imatcols * inumnz);
		Bd.reserve(imatrows * imatcols * inumnz);
		Cd.reserve(imatrows * imatcols * inumnz);
		Dd.reserve(imatrows * imatcols * inumnz);

		duration<long long, ratio<1L, 1000000000L> > sp_grad_time(0), sp_grad_a_time(0), c_full_time(0);

		for (index_type iloop = 0; iloop < inumloops; ++iloop) {
			SpGradDofStat s;

			for (auto& ai: A) {
				sp_grad_rand_gen(ai, randnz, randdof, randval, gen);
				SpGradient::GetDofStat(ai, s);
			}

			for (auto& bi: B) {
				sp_grad_rand_gen(bi, randnz, randdof, randval, gen);
				SpGradient::GetDofStat(bi, s);
			}

			for (auto& ci: C) {
				sp_grad_rand_gen(ci, randnz, randdof, randval, gen);
				SpGradient::GetDofStat(ci, s);
			}

			const index_type nbdirs = s.iNumNz ? s.iMaxDof : 0;

			Ad.resize(imatrows * imatcols * nbdirs);
			Bd.resize(imatrows * imatcols * nbdirs);
			Cd.resize(imatrows * imatcols * nbdirs);
			Dd.resize(imatrows * imatcols * nbdirs);

			for (index_type i = 0; i < imatrows * imatcols; ++i) {
				Av[i] = SpGradient::dGetValue(A.GetElem(i + 1));

				for (index_type j = 1; j <= s.iMaxDof; ++j) {
					Ad[i + (j - 1) * imatcols * imatrows] = SpGradient::dGetDeriv(A.GetElem(i + 1), j);
				}
			}

			for (index_type i = 0; i < imatrows * imatcols; ++i) {
				Bv[i] = SpGradient::dGetValue(B.GetElem(i + 1));

				for (index_type j = 1; j <= s.iMaxDof; ++j) {
					Bd[i + (j - 1) * imatcols * imatrows] = SpGradient::dGetDeriv(B.GetElem(i + 1), j);
				}
			}

			for (index_type i = 0; i < imatrows * imatcols; ++i) {
				Cv[i] = SpGradient::dGetValue(C.GetElem(i + 1));

				for (index_type j = 1; j <= s.iMaxDof; ++j) {
					Cd[i + (j - 1) * imatcols * imatrows] = SpGradient::dGetDeriv(C.GetElem(i + 1), j);
				}
			}

			auto start = high_resolution_clock::now();

			func_mat_add7<TA, TB, TC>(imatrows * imatcols,
						  nbdirs,
						  front(Av),
						  front(Ad),
						  front(Bv),
						  front(Bd),
						  front(Cv),
						  front(Cd),
						  front(Dv),
						  front(Dd));

			c_full_time += high_resolution_clock::now() - start;

			start = high_resolution_clock::now();

			func_mat_add7(A, B, C, D);

			sp_grad_time += high_resolution_clock::now() - start;

			start = high_resolution_clock::now();

			func_mat_add7a(A, B, C, D);

			sp_grad_a_time += high_resolution_clock::now() - start;

			for (index_type j = 1; j <= imatcols; ++j) {
				for (index_type i = 1; i <= imatrows; ++i) {
					const doublereal aij = SpGradient::dGetValue(A.GetElem(i, j));
					const doublereal bij = SpGradient::dGetValue(B.GetElem(i, j));
					const doublereal cij = SpGradient::dGetValue(C.GetElem(i, j));
					const doublereal dijr1 = SpGradient::dGetValue(D.GetElem(i, j));
					const doublereal dijr2 = Dv[(j - 1) * imatrows + i - 1];
					const doublereal dij = -aij / 3. + (bij * 5. - cij / 4.) * 1.5;

					assert(fabs(dij - dijr1) / std::max(1., fabs(dij)) < dTol);
					assert(fabs(dij - dijr2) / std::max(1., fabs(dij)) < dTol);

					for (index_type k = s.iMinDof; k <= s.iMaxDof; ++k) {
						const doublereal daij = SpGradient::dGetDeriv(A.GetElem(i, j), k);
						const doublereal dbij = SpGradient::dGetDeriv(B.GetElem(i, j), k);
						const doublereal dcij = SpGradient::dGetDeriv(C.GetElem(i, j), k);
						const doublereal ddijr1 = SpGradient::dGetDeriv(D.GetElem(i, j), k);
						const doublereal ddijr2 = Dd[((j - 1) * imatrows + (i - 1) + (k - 1) * imatrows * imatcols)];
						const doublereal ddij = -daij / 3. + (dbij * 5. - dcij / 4.) * 1.5;

						assert(fabs(ddij - ddijr1) / std::max(1., fabs(ddij)) < dTol);
						assert(fabs(ddij - ddijr2) / std::max(1., fabs(ddij)) < dTol);
					}
				}
			}
		}

		auto sp_grad_time_ns = duration_cast<nanoseconds>(sp_grad_time).count();
		auto sp_grad_a_time_ns = duration_cast<nanoseconds>(sp_grad_a_time).count();
		auto c_full_time_ns = duration_cast<nanoseconds>(c_full_time).count();

		cerr << "test7: sp_grad_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_time_ns) / 1e9
		     << endl;
		cerr << "test7: sp_grad_a_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_a_time_ns) / 1e9
		     << endl;
		cerr << "test7: c_full_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(c_full_time_ns) / 1e9
		     << endl;
		cerr << "test7: sp_grad_time / c_full_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_time_ns) / static_cast<doublereal>(c_full_time_ns)
		     << endl;
		cerr << "test7: sp_grad_a_time / c_full_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_a_time_ns) / static_cast<doublereal>(c_full_time_ns)
		     << endl;
	}

	template <typename TA, typename TB, typename TC, index_type NumRows = SpMatrixSize::DYNAMIC, index_type NumCols = SpMatrixSize::DYNAMIC>
	void test8(const index_type inumloops,
		   const index_type inumnz,
		   const index_type inumdof,
		   const index_type imatrows,
		   const index_type imatcols)
	{
		using namespace std;
		using namespace std::chrono;

		cerr << __PRETTY_FUNCTION__ << ":\n";
		random_device rd;
		mt19937 gen(rd());
		uniform_real_distribution<doublereal> randval(-1., 1.);
		uniform_int_distribution<index_type> randdof(1, inumdof);
		uniform_int_distribution<index_type> randnz(0, inumnz - 1);

		gen.seed(0);

		constexpr doublereal dTol = pow(numeric_limits<doublereal>::epsilon(), 0.5);

		typedef typename util::ResultType<typename util::ResultType<TA, TB>::Type, TC>::Type TD;

		SpMatrixBase<TA, NumRows, NumCols> A(imatrows, imatcols, inumnz);
		SpMatrixBase<TB, NumRows, NumCols> B(imatrows, imatcols, inumnz);
		SpMatrixBase<TB, NumCols, NumRows> B_T(imatcols, imatrows, inumnz);
		TC C;
		SpMatrixBase<TD, NumRows, NumCols> D(imatrows, imatcols, 3 * inumnz);
		SpMatrixBase<TD, NumRows, NumCols> Da(imatrows, imatcols, 3 * inumnz);
		SpMatrixBase<TD, NumRows, NumCols> Db(imatrows, imatcols, 3 * inumnz);
		SpMatrixBase<TD, NumCols, NumRows> D_T(imatcols, imatrows, 3 * inumnz);
		vector<doublereal> Av(imatrows * imatcols), Ad,
			Bv(imatrows * imatcols), Bd,
			Cd,
			Dv(imatrows * imatcols), Dd;
		doublereal Cv;

		Ad.reserve(imatrows * imatcols * inumnz);
		Bd.reserve(imatrows * imatcols * inumnz);
		Cd.reserve(inumnz);
		Dd.reserve(imatrows * imatcols * inumnz);

		duration<long long, ratio<1L, 1000000000L> > sp_grad_time(0), sp_grad_a_time(0), sp_grad_b_time(0), sp_grad_time9(0), c_full_time(0);

		for (index_type iloop = 0; iloop < inumloops; ++iloop) {
			SpGradDofStat s;

			for (auto& ai: A) {
				sp_grad_rand_gen(ai, randnz, randdof, randval, gen);
				SpGradient::GetDofStat(ai, s);
			}

			for (auto& bi: B) {
				sp_grad_rand_gen(bi, randnz, randdof, randval, gen);
				SpGradient::GetDofStat(bi, s);
			}

			B_T = Transpose(B);

			sp_grad_rand_gen(C, randnz, randdof, randval, gen);

			SpGradient::GetDofStat(C, s);

			const index_type nbdirs = s.iNumNz ? s.iMaxDof : 0;

			Ad.resize(imatrows * imatcols * nbdirs);
			Bd.resize(imatrows * imatcols * nbdirs);
			Cd.resize(nbdirs);
			Dd.resize(imatrows * imatcols * nbdirs);

			for (index_type i = 0; i < imatrows * imatcols; ++i) {
				Av[i] = SpGradient::dGetValue(A.GetElem(i + 1));

				for (index_type j = 1; j <= nbdirs; ++j) {
					Ad[i + (j - 1) * imatcols * imatrows] = SpGradient::dGetDeriv(A.GetElem(i + 1), j);
				}
			}

			for (index_type i = 0; i < imatrows * imatcols; ++i) {
				Bv[i] = SpGradient::dGetValue(B.GetElem(i + 1));

				for (index_type j = 1; j <= nbdirs; ++j) {
					Bd[i + (j - 1) * imatcols * imatrows] = SpGradient::dGetDeriv(B.GetElem(i + 1), j);
				}
			}


			Cv = SpGradient::dGetValue(C);

			for (index_type j = 1; j <= nbdirs; ++j) {
				Cd[j - 1] = SpGradient::dGetDeriv(C, j);
			}

			auto start = high_resolution_clock::now();

			func_mat_add8(imatrows * imatcols,
				      nbdirs,
				      front(Av),
				      front(Ad),
				      front(Bv),
				      front(Bd),
				      Cv,
				      front(Cd),
				      front(Dv),
				      front(Dd));

			c_full_time += high_resolution_clock::now() - start;

			start = high_resolution_clock::now();

			func_mat_add8(A, B, C, D);

			sp_grad_time += high_resolution_clock::now() - start;

			start = high_resolution_clock::now();

			func_mat_add8a(A, B, C, Da);

			sp_grad_a_time += high_resolution_clock::now() - start;

			start = high_resolution_clock::now();

			func_mat_add8b(A, B, C, Db);

			sp_grad_b_time += high_resolution_clock::now() - start;

			start = high_resolution_clock::now();

			func_mat_add9(A, B_T, C, D_T);

			sp_grad_time9 += high_resolution_clock::now() - start;

			for (index_type j = 1; j <= imatcols; ++j) {
				for (index_type i = 1; i <= imatrows; ++i) {
					const doublereal dij1 = SpGradient::dGetValue(D.GetElem(i, j));
					const doublereal dij2 = Dv[(j - 1) * imatrows + i - 1];
					const doublereal dij3 = SpGradient::dGetValue(D_T.GetElem(j, i));
					const doublereal dij4 = SpGradient::dGetValue(Da.GetElem(i, j));
					const doublereal dij5 = SpGradient::dGetValue(Db.GetElem(i, j));

					assert(fabs(dij1 - dij2) / std::max(1., fabs(dij2)) < dTol);
					assert(fabs(dij3 - dij2) / std::max(1., fabs(dij2)) < dTol);
					assert(fabs(dij4 - dij2) / std::max(1., fabs(dij2)) < dTol);
					assert(fabs(dij5 - dij2) / std::max(1., fabs(dij2)) < dTol);

					for (index_type k = 1; k <= nbdirs; ++k) {
						const doublereal ddij1 = SpGradient::dGetDeriv(D.GetElem(i, j), k);
						const doublereal ddij2 = Dd[((j - 1) * imatrows + (i - 1) + (k - 1) * imatrows * imatcols)];
						const doublereal ddij3 = SpGradient::dGetDeriv(D_T.GetElem(j, i), k);
						const doublereal ddij4 = SpGradient::dGetDeriv(Da.GetElem(i, j), k);
						const doublereal ddij5 = SpGradient::dGetDeriv(Db.GetElem(i, j), k);
						assert(fabs(ddij1 - ddij2) / std::max(1., fabs(ddij2)) < dTol);
						assert(fabs(ddij3 - ddij2) / std::max(1., fabs(ddij2)) < dTol);
						assert(fabs(ddij4 - ddij2) / std::max(1., fabs(ddij2)) < dTol);
						assert(fabs(ddij5 - ddij2) / std::max(1., fabs(ddij2)) < dTol);
					}
				}
			}
		}

		auto sp_grad_time_ns = duration_cast<nanoseconds>(sp_grad_time).count();
		auto sp_grad_a_time_ns = duration_cast<nanoseconds>(sp_grad_a_time).count();
		auto sp_grad_b_time_ns = duration_cast<nanoseconds>(sp_grad_b_time).count();
		auto sp_grad_time9_ns = duration_cast<nanoseconds>(sp_grad_time9).count();
		auto c_full_time_ns = duration_cast<nanoseconds>(c_full_time).count();

		cerr << "test8: sp_grad_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_time_ns) / 1e9
		     << endl;
		cerr << "test8: sp_grad_a_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_a_time_ns) / 1e9
		     << endl;
		cerr << "test8: sp_grad_b_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_b_time_ns) / 1e9
		     << endl;
		cerr << "test8: sp_grad_time9 = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_time9_ns) / 1e9
		     << endl;
		cerr << "test8: c_full_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(c_full_time_ns) / 1e9
		     << endl;
		cerr << "test8: sp_grad_time / c_full_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_time_ns) / static_cast<doublereal>(c_full_time_ns)
		     << endl;
		cerr << "test8: sp_grad_time9 / c_full_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_time9_ns) / static_cast<doublereal>(c_full_time_ns)
		     << endl;
	}

	template <typename TA, typename TB>
	void test10(const index_type inumloops,
		    const index_type inumnz,
		    const index_type inumdof,
		    const index_type iamatrows,
		    const index_type iamatcols,
		    const index_type ibmatcols)
	{
		using namespace std;
		using namespace std::chrono;

		cerr << __PRETTY_FUNCTION__ << ":\n";
		random_device rd;
		mt19937 gen(rd());
		uniform_real_distribution<doublereal> randval(-1., 1.);
		uniform_int_distribution<index_type> randdof(1, inumdof);
		uniform_int_distribution<index_type> randnz(0, inumnz - 1);

		gen.seed(0);

		constexpr doublereal dTol = pow(numeric_limits<doublereal>::epsilon(), 0.5);

		typedef typename util::ResultType<TA, TB>::Type TC;
#ifdef USE_AUTODIFF
		typedef typename SpGradient2SliceGradient<TA>::Type TAs;
		typedef typename SpGradient2SliceGradient<TB>::Type TBs;
		typedef typename util::ResultType<TAs, TBs>::Type TCs;
#endif
		SpMatrixBase<TA> A(iamatrows, iamatcols, inumnz);
		SpMatrixBase<TB> B(iamatcols, ibmatcols, inumnz);
		SpMatrixBase<TC> C(iamatrows, ibmatcols, 2 * inumnz);
		SpMatrixBase<TC> C_T(iamatrows, ibmatcols, 2 * inumnz);
		SpMatrixBase<TC> C_TA(iamatrows, ibmatcols, 2 * inumnz);
#ifdef USE_AUTODIFF
		grad::Matrix<TAs, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE> As(iamatrows, iamatcols);
		grad::Matrix<TBs, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE> Bs(iamatcols, ibmatcols);
		grad::Matrix<TCs, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE> Cs(iamatrows, ibmatcols);
		grad::Matrix<TCs, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE> Cs_TAs(iamatrows, ibmatcols);
		grad::Matrix<TCs, grad::DYNAMIC_SIZE, grad::DYNAMIC_SIZE> Cs_dof_map(iamatrows, ibmatcols);
		grad::LocalDofMap oDofMap;

		for (index_type j = 1; j <= Cs_dof_map.iGetNumCols(); ++j) {
			for (index_type i = 1; i <= Cs_dof_map.iGetNumRows(); ++i) {
				TCs& Csij = Cs_dof_map(i, j);
				sp_grad_test::SetDofMap(Csij, &oDofMap);
			}
		}
#endif
		vector<doublereal> Av(iamatrows * iamatcols), Ad,
			Bv(iamatcols * ibmatcols), Bd,
			Cv(iamatrows * ibmatcols), Cd;

		Ad.reserve(iamatrows * iamatcols * inumnz);
		Bd.reserve(iamatcols * ibmatcols * inumnz);
		Cd.reserve(iamatrows * ibmatcols * inumnz);

		duration<long long, ratio<1L, 1000000000L> > sp_grad_time(0),
			sp_grad_trans_time(0),
			sp_grad_trans_add_time(0),
#ifdef USE_AUTODIFF
			sp_grad_s_time(0),
			sp_grad_s_trans_add_time(0),
			sp_grad_s_dof_map_time(0),
#endif
			c_full_time(0);

		for (index_type iloop = 0; iloop < inumloops; ++iloop) {
			SpGradDofStat s;

			for (auto& ai: A) {
				sp_grad_rand_gen(ai, randnz, randdof, randval, gen);
				SpGradient::GetDofStat(ai, s);
			}

			for (auto& bi: B) {
				sp_grad_rand_gen(bi, randnz, randdof, randval, gen);
				SpGradient::GetDofStat(bi, s);
			}
#ifdef USE_AUTODIFF
			for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
				for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
					SpGradient2SliceGradient<TA>::Copy(A.GetElem(i, j), As(i, j), oDofMap);
				}
			}

			for (index_type j = 1; j <= B.iGetNumCols(); ++j) {
				for (index_type i = 1; i <= B.iGetNumRows(); ++i) {
					SpGradient2SliceGradient<TB>::Copy(B.GetElem(i, j), Bs(i, j), oDofMap);
				}
			}
#endif
			const index_type nbdirs = s.iNumNz ? s.iMaxDof : 0;

			Ad.resize(iamatrows * iamatcols * nbdirs);
			Bd.resize(iamatcols * ibmatcols * nbdirs);
			Cd.resize(iamatrows * ibmatcols * nbdirs);

			for (index_type i = 0; i < iamatrows * iamatcols; ++i) {
				Av[i] = SpGradient::dGetValue(A.GetElem(i + 1));

				for (index_type j = 1; j <= nbdirs; ++j) {
					Ad[i + (j - 1) * iamatcols * iamatrows] = SpGradient::dGetDeriv(A.GetElem(i + 1), j);
				}
			}

			for (index_type i = 0; i < iamatcols * ibmatcols; ++i) {
				Bv[i] = SpGradient::dGetValue(B.GetElem(i + 1));

				for (index_type j = 1; j <= nbdirs; ++j) {
					Bd[i + (j - 1) * iamatcols * ibmatcols] = SpGradient::dGetDeriv(B.GetElem(i + 1), j);
				}
			}

			auto start = high_resolution_clock::now();

			func_mat_mul10(iamatrows,
				       iamatcols,
				       ibmatcols,
				       nbdirs,
				       front(Av),
				       front(Ad),
				       front(Bv),
				       front(Bd),
				       front(Cv),
				       front(Cd));

			c_full_time += high_resolution_clock::now() - start;

			start = high_resolution_clock::now();

			func_mat_mul10(A, B, C);

			sp_grad_time += high_resolution_clock::now() - start;

			start = high_resolution_clock::now();

			func_mat_mul10_trans(A, B, C_T);

			sp_grad_trans_time += high_resolution_clock::now() - start;

			start = high_resolution_clock::now();

			func_mat_mul10_trans_add(A, B, C_TA);

			sp_grad_trans_add_time += high_resolution_clock::now() - start;
#ifdef USE_AUTODIFF
			start = high_resolution_clock::now();

			func_mat_mul10(As, Bs, Cs);

			sp_grad_s_time += high_resolution_clock::now() - start;

			start = high_resolution_clock::now();

			func_mat_mul10_trans_add(As, Bs, Cs_TAs);

			sp_grad_s_trans_add_time += high_resolution_clock::now() - start;

			start = high_resolution_clock::now();

			func_mat_mul10_dof_map(As, Bs, Cs_dof_map);

			sp_grad_s_dof_map_time += high_resolution_clock::now() - start;
#endif
			for (index_type j = 1; j <= ibmatcols; ++j) {
				for (index_type i = 1; i <= iamatrows; ++i) {
					const doublereal cij1 = SpGradient::dGetValue(C.GetElem(i, j));
					const doublereal cij2 = Cv[(j - 1) * iamatrows + i - 1];
					const doublereal cij3 = SpGradient::dGetValue(C_T.GetElem(i, j));
					const doublereal cij4 = SpGradient::dGetValue(C_TA.GetElem(i, j));

					assert(fabs(cij1 - cij2) / std::max(1., fabs(cij2)) < dTol);
					assert(fabs(cij3 - cij2) / std::max(1., fabs(cij2)) < dTol);
					assert(fabs(cij4 - cij2) / std::max(1., fabs(cij2)) < dTol);

					for (index_type k = 1; k <= nbdirs; ++k) {
						const doublereal dcij1 = SpGradient::dGetDeriv(C.GetElem(i, j), k);
						const doublereal dcij2 = Cd[((j - 1) * iamatrows + (i - 1) + (k - 1) * iamatrows * ibmatcols)];
						const doublereal dcij3 = SpGradient::dGetDeriv(C_T.GetElem(i, j), k);
						const doublereal dcij4 = SpGradient::dGetDeriv(C_TA.GetElem(i, j), k);

						assert(fabs(dcij1 - dcij2) / std::max(1., fabs(dcij2)) < dTol);
						assert(fabs(dcij3 - dcij2) / std::max(1., fabs(dcij2)) < dTol);
						assert(fabs(dcij4 - dcij2) / std::max(1., fabs(dcij2)) < dTol);
					}
				}
			}
#ifdef USE_AUTODIFF
			for (index_type j = 1; j <= Cs.iGetNumCols(); ++j) {
				for (index_type i = 1; i <= Cs.iGetNumRows(); ++i) {
					const doublereal cij1 = Cv[(j - 1) * iamatrows + i - 1];
					const doublereal cij2 = grad::dGetValue(Cs(i, j));
					const doublereal cij4 = grad::dGetValue(Cs_TAs(i, j));
					const doublereal cij5 = grad::dGetValue(Cs_dof_map(i, j));

					assert(fabs(cij2 - cij1) / std::max(1., fabs(cij1)) < dTol);
					assert(fabs(cij4 - cij1) / std::max(1., fabs(cij1)) < dTol);
					assert(fabs(cij5 - cij1) / std::max(1., fabs(cij1)) < dTol);

					for (index_type k = 1; k <= nbdirs; ++k) {
						const doublereal dcij1 = Cd[((j - 1) * iamatrows + (i - 1) + (k - 1) * iamatrows * ibmatcols)];
						const doublereal dcij2 = SpGradient2SliceGradient<TC>::dGetDerivativeGlobal(Cs(i, j), k);
						const doublereal dcij4 = SpGradient2SliceGradient<TC>::dGetDerivativeGlobal(Cs_TAs(i, j), k);
						const doublereal dcij5 = SpGradient2SliceGradient<TC>::dGetDerivativeGlobal(Cs_dof_map(i, j), k);

						assert(fabs(dcij2 - dcij1) / std::max(1., fabs(dcij1)) < dTol);
						assert(fabs(dcij4 - dcij1) / std::max(1., fabs(dcij1)) < dTol);
						assert(fabs(dcij5 - dcij1) / std::max(1., fabs(dcij1)) < dTol);
					}
				}
			}
#endif
		}

		auto sp_grad_time_ns = duration_cast<nanoseconds>(sp_grad_time).count();
		auto sp_grad_trans_time_ns = duration_cast<nanoseconds>(sp_grad_trans_time).count();
		auto sp_grad_trans_add_time_ns = duration_cast<nanoseconds>(sp_grad_trans_add_time).count();
#ifdef USE_AUTODIFF
		auto sp_grad_s_time_ns = duration_cast<nanoseconds>(sp_grad_s_time).count();
		auto sp_grad_s_trans_add_time_ns = duration_cast<nanoseconds>(sp_grad_s_trans_add_time).count();
		auto sp_grad_s_dof_map_time_ns = duration_cast<nanoseconds>(sp_grad_s_dof_map_time).count();
#endif
		auto c_full_time_ns = duration_cast<nanoseconds>(c_full_time).count();

		cerr << "test10: sp_grad_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_time_ns) / 1e9
		     << endl;

		cerr << "test10: sp_grad_trans_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_trans_time_ns) / 1e9
		     << endl;

		cerr << "test10: sp_grad_trans_add_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_trans_add_time_ns) / 1e9
		     << endl;
#ifdef USE_AUTODIFF
		cerr << "test10: sp_grad_s_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_s_time_ns) / 1e9
		     << endl;

		cerr << "test10: sp_grad_s_trans_add_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_s_trans_add_time_ns) / 1e9
		     << endl;

		cerr << "test10: sp_grad_s_dof_map_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_s_dof_map_time_ns) / 1e9
		     << endl;
#endif
		cerr << "test10: c_full_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(c_full_time_ns) / 1e9
		     << endl;

		cerr << "test10: sp_grad_time / c_full_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_time_ns) / static_cast<doublereal>(c_full_time_ns)
		     << endl;

		cerr << "test10: sp_grad_trans_time / c_full_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_trans_time_ns) / static_cast<doublereal>(c_full_time_ns)
		     << endl;

		cerr << "test10: sp_grad_trans_add_time / c_full_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_trans_add_time_ns) / static_cast<doublereal>(c_full_time_ns)
		     << endl;
#ifdef USE_AUTODIFF
		cerr << "test10: sp_grad_s_time / c_full_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_s_time_ns) / static_cast<doublereal>(c_full_time_ns)
		     << endl;
#endif
	}


	template <typename TA, typename TB, typename TC>
	void test11(const index_type inumloops,
		    const index_type inumnz,
		    const index_type inumdof,
		    const index_type imatrows)
	{
		using namespace std;
		using namespace std::chrono;

		cerr << __PRETTY_FUNCTION__ << ":\n";
		random_device rd;
		mt19937 gen(rd());
		uniform_real_distribution<doublereal> randval(-1., 1.);
		uniform_int_distribution<index_type> randdof(1, inumdof);
		uniform_int_distribution<index_type> randnz(0, inumnz - 1);

		gen.seed(0);

		constexpr doublereal dTol = pow(numeric_limits<doublereal>::epsilon(), 0.5);

		typedef typename util::ResultType<typename util::ResultType<TA, TB>::Type, TC>::Type TD;

		const index_type imatcols = imatrows; // Needed for func_mat_mul11

		SpMatrixBase<TA> A(imatrows, imatcols, inumnz);
		SpMatrixBase<TB> B(imatcols, imatcols, inumnz);
		SpMatrixBase<TC> C(imatrows, imatcols, inumnz);
		SpMatrixBase<TD> D(imatrows, imatcols, inumnz * 4);
		SpMatrixBase<TD> D_T(imatrows, imatcols, inumnz * 4);

		vector<doublereal>
			Av(imatrows * imatcols), Ad,
			Bv(imatrows * imatcols), Bd,
			Cv(imatrows * imatcols), Cd,
			Dv(imatrows * imatcols), Dd,
			Tmp1(imatrows * imatcols), Tmp1d,
			Tmp2(imatrows * imatcols), Tmp2d;

		duration<long long, ratio<1L, 1000000000L> > sp_grad_time(0), sp_grad_trans_time(0), c_full_time(0);

		for (index_type iloop = 0; iloop < inumloops; ++iloop) {
			SpGradDofStat s;

			for (auto& ai: A) {
				sp_grad_rand_gen(ai, randnz, randdof, randval, gen);
				SpGradient::GetDofStat(ai, s);
			}

			for (auto& bi: B) {
				sp_grad_rand_gen(bi, randnz, randdof, randval, gen);
				SpGradient::GetDofStat(bi, s);
			}

			for (auto& ci: C) {
				sp_grad_rand_gen(ci, randnz, randdof, randval, gen);
				SpGradient::GetDofStat(ci, s);
			}

			const index_type nbdirs = s.iNumNz ? s.iMaxDof : 0;

			Ad.resize(imatrows * imatcols * nbdirs);
			Bd.resize(imatrows * imatcols * nbdirs);
			Cd.resize(imatrows * imatcols * nbdirs);
			Dd.resize(imatrows * imatcols * nbdirs);
			Tmp1d.resize(imatrows * imatcols * nbdirs);
			Tmp2d.resize(imatrows * imatcols * nbdirs);

			for (index_type i = 0; i < imatrows * imatcols; ++i) {
				Av[i] = SpGradient::dGetValue(A.GetElem(i + 1));

				for (index_type j = 1; j <= nbdirs; ++j) {
					Ad[i + (j - 1) * imatcols * imatrows] = SpGradient::dGetDeriv(A.GetElem(i + 1), j);
				}
			}

			for (index_type i = 0; i < imatcols * imatcols; ++i) {
				Bv[i] = SpGradient::dGetValue(B.GetElem(i + 1));

				for (index_type j = 1; j <= nbdirs; ++j) {
					Bd[i + (j - 1) * imatcols * imatcols] = SpGradient::dGetDeriv(B.GetElem(i + 1), j);
				}
			}

			for (index_type i = 0; i < imatcols * imatcols; ++i) {
				Cv[i] = SpGradient::dGetValue(C.GetElem(i + 1));

				for (index_type j = 1; j <= nbdirs; ++j) {
					Cd[i + (j - 1) * imatcols * imatcols] = SpGradient::dGetDeriv(C.GetElem(i + 1), j);
				}
			}

			auto start = high_resolution_clock::now();

			func_mat_mul11(imatrows,
				       nbdirs,
				       front(Av),
				       front(Ad),
				       front(Bv),
				       front(Bd),
				       front(Cv),
				       front(Cd),
				       front(Dv),
				       front(Dd),
				       front(Tmp1),
				       front(Tmp1d),
				       front(Tmp2),
				       front(Tmp2d));

			c_full_time += high_resolution_clock::now() - start;

			start = high_resolution_clock::now();

			func_mat_mul11(A, B, C, D);

			sp_grad_time += high_resolution_clock::now() - start;

			start = high_resolution_clock::now();

			func_mat_mul11_trans(A, B, C, D_T);

			sp_grad_trans_time += high_resolution_clock::now() - start;

			for (index_type j = 1; j <= imatcols; ++j) {
				for (index_type i = 1; i <= imatrows; ++i) {
					const doublereal dij1 = SpGradient::dGetValue(D.GetElem(i, j));
					const doublereal dij2 = Dv[(j - 1) * imatrows + i - 1];
					const doublereal dij3 = SpGradient::dGetValue(D_T.GetElem(i, j));

					assert(fabs(dij1 - dij2) / std::max(1., fabs(dij2)) < dTol);
					assert(fabs(dij3 - dij2) / std::max(1., fabs(dij2)) < dTol);

					for (index_type k = 1; k <= nbdirs; ++k) {
						const doublereal ddij1 = SpGradient::dGetDeriv(D.GetElem(i, j), k);
						const doublereal ddij2 = Dd[((j - 1) * imatrows + (i - 1) + (k - 1) * imatrows * imatcols)];
						const doublereal ddij3 = SpGradient::dGetDeriv(D_T.GetElem(i, j), k);
						assert(fabs(ddij1 - ddij2) / std::max(1., fabs(ddij2)) < dTol);
						assert(fabs(ddij3 - ddij2) / std::max(1., fabs(ddij2)) < dTol);
					}
				}
			}
		}

		auto sp_grad_time_ns = duration_cast<nanoseconds>(sp_grad_time).count();
		auto sp_grad_trans_time_ns = duration_cast<nanoseconds>(sp_grad_trans_time).count();
		auto c_full_time_ns = duration_cast<nanoseconds>(c_full_time).count();

		cerr << "test11: sp_grad_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_time_ns) / 1e9
		     << endl;

		cerr << "test11: sp_grad_trans_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_trans_time_ns) / 1e9
		     << endl;

		cerr << "test11: c_full_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(c_full_time_ns) / 1e9
		     << endl;

		cerr << "test11: sp_grad_time / c_full_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_time_ns) / static_cast<doublereal>(c_full_time_ns)
		     << endl;

		cerr << "test11: sp_grad_trans_time / c_full_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_trans_time_ns) / static_cast<doublereal>(c_full_time_ns)
		     << endl;
	}

	template <typename TA, typename TB, typename TC>
	void test12(const index_type inumloops,
		    const index_type inumnz,
		    const index_type inumdof,
		    const index_type imatrowsa,
		    const index_type imatcolsa,
		    const index_type imatcolsb,
		    const index_type imatcolsc)
	{
		using namespace std;
		using namespace std::chrono;

		cerr << __PRETTY_FUNCTION__ << ":\n";
		random_device rd;
		mt19937 gen(rd());
		uniform_real_distribution<doublereal> randval(-1., 1.);
		uniform_int_distribution<index_type> randdof(1, inumdof);
		uniform_int_distribution<index_type> randnz(0, inumnz - 1);

		gen.seed(0);

		constexpr doublereal dTol = pow(numeric_limits<doublereal>::epsilon(), 0.5);

		typedef typename util::ResultType<typename util::ResultType<TA, TB>::Type, TC>::Type TD;

		SpMatrixBase<TA> A(imatrowsa, imatcolsa, inumnz);
		SpMatrixBase<TB> B(imatcolsa, imatcolsb, inumnz);
		SpMatrixBase<TC> C(imatcolsb, imatcolsc, inumnz);
		SpMatrixBase<TD> Da(imatcolsc, imatrowsa, 3 * inumnz),
			Db(imatcolsc, imatrowsa, 3 * inumnz),
			Dc(imatcolsc, imatrowsa, 3 * inumnz);

		vector<doublereal>
			Av(imatrowsa * imatcolsa), Ad,
			Bv(imatcolsa * imatcolsb), Bd,
			Cv(imatcolsb * imatcolsc), Cd,
			Dv(imatrowsa * imatcolsc), Dd,
			Tmp1(imatrowsa * imatcolsb), Tmp1d;

		duration<long long, ratio<1L, 1000000000L> > sp_grad_time_a(0),
			sp_grad_time_b(0),
			sp_grad_time_c(0),
			c_full_time(0);

		for (index_type iloop = 0; iloop < inumloops; ++iloop) {
			SpGradDofStat s;

			for (auto& ai: A) {
				sp_grad_rand_gen(ai, randnz, randdof, randval, gen);
				SpGradient::GetDofStat(ai, s);
			}

			for (auto& bi: B) {
				sp_grad_rand_gen(bi, randnz, randdof, randval, gen);
				SpGradient::GetDofStat(bi, s);
			}

			for (auto& ci: C) {
				sp_grad_rand_gen(ci, randnz, randdof, randval, gen);
				SpGradient::GetDofStat(ci, s);
			}

			const index_type nbdirs = s.iNumNz ? s.iMaxDof : 0;

			Ad.resize(imatrowsa * imatcolsa * nbdirs);
			Bd.resize(imatcolsa * imatcolsb * nbdirs);
			Cd.resize(imatcolsb * imatcolsc * nbdirs);
			Dd.resize(imatrowsa * imatcolsc * nbdirs);
			Tmp1d.resize(imatrowsa * imatcolsb * nbdirs);

			for (index_type i = 0; i < imatrowsa * imatcolsa; ++i) {
				Av[i] = SpGradient::dGetValue(A.GetElem(i + 1));

				for (index_type j = 1; j <= nbdirs; ++j) {
					Ad[i + (j - 1) * imatcolsa * imatrowsa] = SpGradient::dGetDeriv(A.GetElem(i + 1), j);
				}
			}

			for (index_type i = 0; i < imatcolsa * imatcolsb; ++i) {
				Bv[i] = SpGradient::dGetValue(B.GetElem(i + 1));

				for (index_type j = 1; j <= nbdirs; ++j) {
					Bd[i + (j - 1) * imatcolsa * imatcolsb] = SpGradient::dGetDeriv(B.GetElem(i + 1), j);
				}
			}

			for (index_type i = 0; i < imatcolsb * imatcolsc; ++i) {
				Cv[i] = SpGradient::dGetValue(C.GetElem(i + 1));

				for (index_type j = 1; j <= nbdirs; ++j) {
					Cd[i + (j - 1) * imatcolsb * imatcolsc] = SpGradient::dGetDeriv(C.GetElem(i + 1), j);
				}
			}

			auto start = high_resolution_clock::now();

			func_mat_mul12(imatrowsa,
				       imatcolsa,
				       imatcolsb,
				       imatcolsc,
				       nbdirs,
				       front(Av),
				       front(Ad),
				       front(Bv),
				       front(Bd),
				       front(Cv),
				       front(Cd),
				       front(Dv),
				       front(Dd),
				       front(Tmp1),
				       front(Tmp1d));

			c_full_time += high_resolution_clock::now() - start;

			start = high_resolution_clock::now();

			func_mat_mul12a(A, B, C, Da);

			sp_grad_time_a += high_resolution_clock::now() - start;

			start = high_resolution_clock::now();

			func_mat_mul12b(A, B, C, Db);

			sp_grad_time_b += high_resolution_clock::now() - start;

			start = high_resolution_clock::now();

			func_mat_mul12c(A, B, C, Dc);

			sp_grad_time_c += high_resolution_clock::now() - start;

			for (index_type j = 1; j <= imatcolsc; ++j) {
				for (index_type i = 1; i <= imatrowsa; ++i) {
					const doublereal dija = SpGradient::dGetValue(Da.GetElem(j, i));
					const doublereal dijb = SpGradient::dGetValue(Db.GetElem(j, i));
					const doublereal dijc = SpGradient::dGetValue(Dc.GetElem(j, i));
					const doublereal dij2 = Dv[(j - 1) * imatrowsa + i - 1];

					assert(fabs(dija - dij2) / std::max(1., fabs(dij2)) < dTol);
					assert(fabs(dijb - dij2) / std::max(1., fabs(dij2)) < dTol);
					assert(fabs(dijc - dij2) / std::max(1., fabs(dij2)) < dTol);

					for (index_type k = 1; k <= nbdirs; ++k) {
						const doublereal ddija = SpGradient::dGetDeriv(Da.GetElem(j, i), k);
						const doublereal ddijb = SpGradient::dGetDeriv(Db.GetElem(j, i), k);
						const doublereal ddijc = SpGradient::dGetDeriv(Dc.GetElem(j, i), k);
						const doublereal ddij2 = Dd[((j - 1) * imatrowsa + (i - 1) + (k - 1) * imatrowsa * imatcolsc)];
						assert(fabs(ddija - ddij2) / std::max(1., fabs(ddij2)) < dTol);
						assert(fabs(ddijb - ddij2) / std::max(1., fabs(ddij2)) < dTol);
						assert(fabs(ddijc - ddij2) / std::max(1., fabs(ddij2)) < dTol);
					}
				}
			}
		}

		auto sp_grad_time_a_ns = duration_cast<nanoseconds>(sp_grad_time_a).count();
		auto sp_grad_time_b_ns = duration_cast<nanoseconds>(sp_grad_time_b).count();
		auto sp_grad_time_c_ns = duration_cast<nanoseconds>(sp_grad_time_c).count();
		auto c_full_time_ns = duration_cast<nanoseconds>(c_full_time).count();

		cerr << "test12: sp_grad_time_a = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_time_a_ns) / 1e9
		     << endl;

		cerr << "test12: sp_grad_time_b = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_time_b_ns) / 1e9
		     << endl;

		cerr << "test12: sp_grad_time_c = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_time_c_ns) / 1e9
		     << endl;

		cerr << "test12: c_full_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(c_full_time_ns) / 1e9
		     << endl;

		cerr << "test12: sp_grad_time_a / c_full_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_time_a_ns) / static_cast<doublereal>(c_full_time_ns)
		     << endl;

		cerr << "test12: sp_grad_time_b / c_full_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_time_b_ns) / static_cast<doublereal>(c_full_time_ns)
		     << endl;

		cerr << "test12: sp_grad_time_c / c_full_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_time_c_ns) / static_cast<doublereal>(c_full_time_ns)
		     << endl;
	}

	template <typename T>
	void test13(const index_type inumloops,
		    const index_type inumnz,
		    const index_type inumdof)
	{
		using namespace std;
		using namespace std::chrono;

		cerr << __PRETTY_FUNCTION__ << ":\n";

		random_device rd;
		mt19937 gen(rd());
		uniform_real_distribution<doublereal> randval(-1., 1.);
		uniform_int_distribution<index_type> randdof(1, inumdof);
		uniform_int_distribution<index_type> randnz(0, inumnz - 1);

		gen.seed(0);

		duration<long long, ratio<1L, 1000000000L> > sp_grad_time(0);

		for (index_type iloop = 0; iloop < inumloops; ++iloop) {
			SpMatrix<T, 3, 4> A(3, 4, inumnz);

			for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
				for (index_type j = 1; j <= A.iGetNumCols(); ++j) {
					sp_grad_rand_gen(A(i, j), randnz, randdof, randval, gen);
				}
			}

			SpMatrix<T, 4, 2> B(4, 2, inumnz);

			for (index_type i = 1; i <= B.iGetNumRows(); ++i) {
				for (index_type j = 1; j <= B.iGetNumCols(); ++j) {
					sp_grad_rand_gen(B(i, j), randnz, randdof, randval, gen);
				}
			}

			SpColVector<T, 4> v1(4, inumnz);

			for (index_type i = 1; i <= v1.iGetNumRows(); ++i) {
				sp_grad_rand_gen(v1(i), randnz, randdof, randval, gen);
			}

			SpColVector<T, 3> v2(3, inumnz);

			for (index_type i = 1; i <= v2.iGetNumRows(); ++i) {
				sp_grad_rand_gen(v2(i), randnz, randdof, randval, gen);
			}

			auto start = high_resolution_clock::now();

			SpMatrix<T, 3, 2> C = A * B;
			SpMatrix<T, 2, 3> D = Transpose(A * B);
			SpMatrix<T, 2, 3> E = Transpose(B) * Transpose(A);
			SpMatrix<T, 3, 3> F = A * B * Transpose(A * B);
			SpMatrix<T, 3, 3> G = A * B * Transpose(B) * Transpose(A);
			SpMatrix<T, 3, 3> H = C * Transpose(C);
			SpMatrix<T, 2, 2> I = Transpose(C) * C;
			SpMatrix<T, 3, 3> J = Transpose(A * B * Transpose(A * B));
			SpColVector<T, 3> w = F * v2;
			SpRowVector<T, 1> u = Transpose(v1) * v1;
			SpMatrix<T, 4, 4> t = v1 * Transpose(v1);
			SpRowVector<T, 3> s = Transpose(v1) * Transpose(A);
			SpRowVector<T, 3> r = Transpose(A * v1);
			SpColVector<T, 3> x = C.GetCol(1) + C.GetCol(2);
			SpColVector<T, 3> vm1 = Cross(Transpose(r), Transpose(s));
			SpColVector<T, 3> vm2 = EvalUnique(Cross(Transpose(r), Transpose(s)));
			SpColVector<T, 3> vr1(3, 0), vr2(3, 0);

			const T x1 = r(1), y1 = r(2), z1 = r(3);
			const T x2 = s(1), y2 = s(2), z2 = s(3);

			vr1(1) = y1 * z2 - y2 * z1;
			vr1(2) = -(x1 * z2 - x2 * z1);
			vr1(3) = x1 * y2 - x2 * y1;

			vr2(1) = EvalUnique(y1 * z2 - y2 * z1);
			vr2(2) = -EvalUnique(x1 * z2 - x2 * z1);
			vr2(3) = EvalUnique(x1 * y2 - x2 * y1);

			SpMatrix<T, 3, 3> y(3, 3, 0);

			for (index_type i = 1; i <= F.iGetNumRows(); ++i) {
				for (index_type j = 1; j <= F.iGetNumCols(); ++j) {
					y = F.GetCol(j) * F.GetRow(i);
				}
			}

			for (index_type i = 1; i <= F.iGetNumRows(); ++i) {
				for (index_type j = 1; j <= F.iGetNumCols(); ++j) {
					y = Transpose(F.GetRow(i)) * Transpose(F.GetCol(j));
				}
			}

			SpMatrix<T, A.iNumRowsStatic, B.iNumColsStatic> K(A.iNumRowsStatic, B.iNumColsStatic, 0);
			SpMatrix<T, A.iNumRowsStatic, B.iNumColsStatic> L(A.iNumRowsStatic, B.iNumColsStatic, 0);

			for (index_type j = 1; j <= B.iGetNumCols(); ++j) {
				for (index_type i = 1; i <= A.iGetNumRows(); ++i) {
					K(i, j) = Dot(Transpose(A.GetRow(i)), B.GetCol(j));
					L(i, j) = Dot(B.GetCol(j), Transpose(A.GetRow(i)));
				}
			}

			T q1 = Dot(v2, F * v2);
			T q2 = Dot(v2, w);
			T q3 = Dot(F * v2, v2);

			sp_grad_time += high_resolution_clock::now() - start;

			constexpr doublereal dTol = sqrt(std::numeric_limits<doublereal>::epsilon());

			sp_grad_assert_equal(q1, q2, dTol);
			sp_grad_assert_equal(q1, q3, dTol);

			assert(SpGradient::dGetValue(q1) >= 0.);

			for (index_type i = 1; i <= 3; ++i) {
				for (index_type j = 1; j <= 2; ++j) {
					sp_grad_assert_equal(C(i, j), D(j, i), dTol);
					sp_grad_assert_equal(C(i, j), E(j, i), dTol);
					sp_grad_assert_equal(C(i, j), K(i, j), dTol);
					sp_grad_assert_equal(C(i, j), L(i, j), dTol);
				}
			}

			for (index_type i = 1; i <= 3; ++i) {
				for (index_type j = 1; j <= 3; ++j) {
					sp_grad_assert_equal(F(i, j), G(i, j), dTol);
					sp_grad_assert_equal(F(i, j), H(i, j), dTol);
					sp_grad_assert_equal(F(i, j), J(i, j), dTol);
				}
			}

			for (index_type i = 1; i <= H.iGetNumRows(); ++i) {
				for (index_type j = i; j <= H.iGetNumCols(); ++j) {
					sp_grad_assert_equal(H(i, j), H(j, i), dTol);
				}
			}

			for (index_type i = 1; i <= I.iGetNumRows(); ++i) {
				for (index_type j = i; j <= I.iGetNumCols(); ++j) {
					sp_grad_assert_equal(I(i, j), I(j, i), dTol);
				}
			}

			for (index_type i = 1; i <= t.iGetNumRows(); ++i) {
				for (index_type j = i; j <= t.iGetNumCols(); ++j) {
					sp_grad_assert_equal(t(i, j), t(j, i), dTol);
				}
			}

			for (index_type i = 1; i <= 3; ++i) {
				sp_grad_assert_equal(vm1(i), vr1(i), dTol);
				sp_grad_assert_equal(vm2(i), vr1(i), dTol);
				sp_grad_assert_equal(vr2(i), vr1(i), dTol);
			}

			start = high_resolution_clock::now();

			I = Transpose(C) * C;
			H = C * Transpose(C);
			F = A * B * Transpose(A * B);
			G = A * B * Transpose(B) * Transpose(A);
			C = A * B;
			D = Transpose(A * B);
			E = Transpose(B) * Transpose(A);
			w = F * v2;
			r = Transpose(w);

			sp_grad_time += high_resolution_clock::now() - start;

			for (index_type i = 1; i <= 3; ++i) {
				for (index_type j = 1; j <= 2; ++j) {
					sp_grad_assert_equal(C(i, j), D(j, i), dTol);
					sp_grad_assert_equal(C(i, j), E(j, i), dTol);
				}
			}

			for (index_type i = 1; i <= 3; ++i) {
				for (index_type j = 1; j <= 3; ++j) {
					sp_grad_assert_equal(F(i, j), G(i, j), dTol);
					sp_grad_assert_equal(F(i, j), H(i, j), dTol);
				}
			}
		}

		auto sp_grad_time_ns = duration_cast<nanoseconds>(sp_grad_time).count();

		cerr << "test13: sp_grad_time = " << fixed << setprecision(6)
		     << static_cast<doublereal>(sp_grad_time_ns) / 1e9
		     << endl;
	}

	template <typename T>
	void test15(const index_type inumloops,
		    const index_type inumnz,
		    const index_type inumdof)
	{
		using namespace std;
		using namespace std::chrono;

		cerr << __PRETTY_FUNCTION__ << ":\n";

		random_device rd;
		mt19937 gen(rd());
		uniform_real_distribution<doublereal> randval(-M_PI, M_PI);
		uniform_int_distribution<index_type> randdof(1, inumdof);
		uniform_int_distribution<index_type> randnz(0, inumnz - 1);

		gen.seed(0);

		constexpr doublereal dTol = pow(numeric_limits<doublereal>::epsilon(), 0.5);

		SpMatrix<T, 3, 3> R(3, 3, inumnz), R2(3, 3, inumnz);
		SpColVector<T, 3> Phi(3, inumnz), Phi2(3, inumnz);

		for (index_type iloop = 0; iloop < inumloops; ++iloop) {
			for (auto& phii: Phi) {
			     sp_grad_rand_gen(phii, randnz, randdof, randval, gen);
			}

			R = MatRotVec(Phi);
			Phi2 = VecRotMat(R);
			R2 = MatRotVec(Phi2);

			for (index_type i = 1; i <= 3; ++i) {
			     if (sqrt(Dot(Phi, Phi)) < M_PI) {
				  sp_grad_assert_equal(Phi2(i), Phi(i), dTol);
			     }

			     for (index_type j = 1; j <= 3; ++j) {
				  sp_grad_assert_equal(R2(i, j), R(i, j), dTol);
			     }
			}

		}
	}
}

int main(int argc, char* argv[])
{
	using namespace sp_grad_test;
	using namespace std;

	try {
#ifdef HAVE_FEENABLEEXCEPT
		feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
#endif
		const index_type inumloops = argc > 1 ? atoi(argv[1]) : 1;
		const index_type inumnz = argc > 2 ? atoi(argv[2]) : 100;
		const index_type inumdof = argc > 3 ? atoi(argv[3]) : 200;
		const index_type imatrows = argc > 4 ? atoi(argv[4]) : 10;
		const index_type imatcols = argc > 5 ? atoi(argv[5]) : 8;
		const index_type imatcolsb = argc > 6 ? atoi(argv[6]) : 5;
		const index_type imatcolsc = argc > 7 ? atoi(argv[7]) : 7;

		constexpr doublereal dall_tests = -1.0;
		const doublereal dtest = argc > 8 ? atof(argv[8]) : dall_tests;

#define SP_GRAD_RUN_TEST(number)                                \
		(dtest == dall_tests || dtest == (number))

		if (SP_GRAD_RUN_TEST(0.1)) test0(inumloops, inumnz, inumdof);
		if (SP_GRAD_RUN_TEST(1.1)) test1();
		if (SP_GRAD_RUN_TEST(2.1)) test2(inumloops, inumnz, inumdof);
		if (SP_GRAD_RUN_TEST(3.1)) test3<SpGradient, SpGradient>(inumloops, inumnz, inumdof, imatrows);
		if (SP_GRAD_RUN_TEST(3.2)) test3<doublereal, SpGradient>(inumloops, inumnz, inumdof, imatrows);
		if (SP_GRAD_RUN_TEST(3.3)) test3<SpGradient, doublereal>(inumloops, inumnz, inumdof, imatrows);
		if (SP_GRAD_RUN_TEST(3.4)) test3<doublereal, doublereal>(inumloops, inumnz, inumdof, imatrows);
		if (SP_GRAD_RUN_TEST(4.1)) test4<doublereal, doublereal>(inumloops, inumnz, inumdof, imatrows, imatcols);
		if (SP_GRAD_RUN_TEST(4.2)) test4<doublereal, SpGradient>(inumloops, inumnz, inumdof, imatrows, imatcols);
		if (SP_GRAD_RUN_TEST(4.3)) test4<SpGradient, doublereal>(inumloops, inumnz, inumdof, imatrows, imatcols);
		if (SP_GRAD_RUN_TEST(4.4)) test4<SpGradient, SpGradient>(inumloops, inumnz, inumdof, imatrows, imatcols);
		if (SP_GRAD_RUN_TEST(6.1)) test6(inumloops, inumnz, inumdof, imatrows, imatcols);
		if (SP_GRAD_RUN_TEST(7.1)) test7<SpGradient, SpGradient, SpGradient>(inumloops, inumnz, inumdof, imatrows, imatcols);
		if (SP_GRAD_RUN_TEST(7.2)) test7<doublereal, doublereal, doublereal>(inumloops, inumnz, inumdof, imatrows, imatcols);
		if (SP_GRAD_RUN_TEST(7.3)) test7<SpGradient, doublereal, doublereal>(inumloops, inumnz, inumdof, imatrows, imatcols);
		if (SP_GRAD_RUN_TEST(7.4)) test7<doublereal, SpGradient, doublereal>(inumloops, inumnz, inumdof, imatrows, imatcols);
		if (SP_GRAD_RUN_TEST(7.5)) test7<doublereal, doublereal, SpGradient>(inumloops, inumnz, inumdof, imatrows, imatcols);
		if (SP_GRAD_RUN_TEST(7.6)) test7<SpGradient, SpGradient, doublereal>(inumloops, inumnz, inumdof, imatrows, imatcols);
		if (SP_GRAD_RUN_TEST(7.7)) test7<doublereal, SpGradient, SpGradient>(inumloops, inumnz, inumdof, imatrows, imatcols);
		if (SP_GRAD_RUN_TEST(7.8)) test7<SpGradient, doublereal, SpGradient>(inumloops, inumnz, inumdof, imatrows, imatcols);

		if (SP_GRAD_RUN_TEST(7.1)) test7<SpGradient, SpGradient, SpGradient, iNumRowsStatic1, iNumColsStatic1>(inumloops, inumnz, inumdof, iNumRowsStatic1, iNumColsStatic1);
		if (SP_GRAD_RUN_TEST(7.2)) test7<doublereal, doublereal, doublereal, iNumRowsStatic1, iNumColsStatic1>(inumloops, inumnz, inumdof, iNumRowsStatic1, iNumColsStatic1);
		if (SP_GRAD_RUN_TEST(7.3)) test7<SpGradient, doublereal, doublereal, iNumRowsStatic1, iNumColsStatic1>(inumloops, inumnz, inumdof, iNumRowsStatic1, iNumColsStatic1);
		if (SP_GRAD_RUN_TEST(7.4)) test7<doublereal, SpGradient, doublereal, iNumRowsStatic1, iNumColsStatic1>(inumloops, inumnz, inumdof, iNumRowsStatic1, iNumColsStatic1);
		if (SP_GRAD_RUN_TEST(7.5)) test7<doublereal, doublereal, SpGradient, iNumRowsStatic1, iNumColsStatic1>(inumloops, inumnz, inumdof, iNumRowsStatic1, iNumColsStatic1);
		if (SP_GRAD_RUN_TEST(7.6)) test7<SpGradient, SpGradient, doublereal, iNumRowsStatic1, iNumColsStatic1>(inumloops, inumnz, inumdof, iNumRowsStatic1, iNumColsStatic1);
		if (SP_GRAD_RUN_TEST(7.7)) test7<doublereal, SpGradient, SpGradient, iNumRowsStatic1, iNumColsStatic1>(inumloops, inumnz, inumdof, iNumRowsStatic1, iNumColsStatic1);
		if (SP_GRAD_RUN_TEST(7.8)) test7<SpGradient, doublereal, SpGradient, iNumRowsStatic1, iNumColsStatic1>(inumloops, inumnz, inumdof, iNumRowsStatic1, iNumColsStatic1);

		if (SP_GRAD_RUN_TEST(7.1)) test7<SpGradient, SpGradient, SpGradient, 7, 9>(inumloops, inumnz, inumdof, 7, 9);
		if (SP_GRAD_RUN_TEST(7.2)) test7<doublereal, doublereal, doublereal, 7, 9>(inumloops, inumnz, inumdof, 7, 9);
		if (SP_GRAD_RUN_TEST(7.3)) test7<SpGradient, doublereal, doublereal, 7, 9>(inumloops, inumnz, inumdof, 7, 9);
		if (SP_GRAD_RUN_TEST(7.4)) test7<doublereal, SpGradient, doublereal, 7, 9>(inumloops, inumnz, inumdof, 7, 9);
		if (SP_GRAD_RUN_TEST(7.5)) test7<doublereal, doublereal, SpGradient, 7, 9>(inumloops, inumnz, inumdof, 7, 9);
		if (SP_GRAD_RUN_TEST(7.6)) test7<SpGradient, SpGradient, doublereal, 7, 9>(inumloops, inumnz, inumdof, 7, 9);
		if (SP_GRAD_RUN_TEST(7.7)) test7<doublereal, SpGradient, SpGradient, 7, 9>(inumloops, inumnz, inumdof, 7, 9);
		if (SP_GRAD_RUN_TEST(7.8)) test7<SpGradient, doublereal, SpGradient, 7, 9>(inumloops, inumnz, inumdof, 7, 9);

		if (SP_GRAD_RUN_TEST(8.1)) test8<SpGradient, SpGradient, SpGradient>(inumloops, inumnz, inumdof, imatrows, imatcols);
		if (SP_GRAD_RUN_TEST(8.2)) test8<doublereal, doublereal, doublereal>(inumloops, inumnz, inumdof, imatrows, imatcols);
		if (SP_GRAD_RUN_TEST(8.3)) test8<SpGradient, doublereal, doublereal>(inumloops, inumnz, inumdof, imatrows, imatcols);
		if (SP_GRAD_RUN_TEST(8.4)) test8<doublereal, SpGradient, doublereal>(inumloops, inumnz, inumdof, imatrows, imatcols);
		if (SP_GRAD_RUN_TEST(8.5)) test8<doublereal, doublereal, SpGradient>(inumloops, inumnz, inumdof, imatrows, imatcols);
		if (SP_GRAD_RUN_TEST(8.6)) test8<SpGradient, SpGradient, doublereal>(inumloops, inumnz, inumdof, imatrows, imatcols);
		if (SP_GRAD_RUN_TEST(8.8)) test8<doublereal, SpGradient, SpGradient>(inumloops, inumnz, inumdof, imatrows, imatcols);
		if (SP_GRAD_RUN_TEST(8.9)) test8<SpGradient, doublereal, SpGradient>(inumloops, inumnz, inumdof, imatrows, imatcols);

		if (SP_GRAD_RUN_TEST(8.1)) test8<SpGradient, SpGradient, SpGradient, iNumRowsStatic1, iNumColsStatic1>(inumloops, inumnz, inumdof, iNumRowsStatic1, iNumColsStatic1);
		if (SP_GRAD_RUN_TEST(8.2)) test8<doublereal, doublereal, doublereal, iNumRowsStatic1, iNumColsStatic1>(inumloops, inumnz, inumdof, iNumRowsStatic1, iNumColsStatic1);
		if (SP_GRAD_RUN_TEST(8.3)) test8<SpGradient, doublereal, doublereal, iNumRowsStatic1, iNumColsStatic1>(inumloops, inumnz, inumdof, iNumRowsStatic1, iNumColsStatic1);
		if (SP_GRAD_RUN_TEST(8.4)) test8<doublereal, SpGradient, doublereal, iNumRowsStatic1, iNumColsStatic1>(inumloops, inumnz, inumdof, iNumRowsStatic1, iNumColsStatic1);
		if (SP_GRAD_RUN_TEST(8.5)) test8<doublereal, doublereal, SpGradient, iNumRowsStatic1, iNumColsStatic1>(inumloops, inumnz, inumdof, iNumRowsStatic1, iNumColsStatic1);
		if (SP_GRAD_RUN_TEST(8.6)) test8<SpGradient, SpGradient, doublereal, iNumRowsStatic1, iNumColsStatic1>(inumloops, inumnz, inumdof, iNumRowsStatic1, iNumColsStatic1);
		if (SP_GRAD_RUN_TEST(8.8)) test8<doublereal, SpGradient, SpGradient, iNumRowsStatic1, iNumColsStatic1>(inumloops, inumnz, inumdof, iNumRowsStatic1, iNumColsStatic1);
		if (SP_GRAD_RUN_TEST(8.9)) test8<SpGradient, doublereal, SpGradient, iNumRowsStatic1, iNumColsStatic1>(inumloops, inumnz, inumdof, iNumRowsStatic1, iNumColsStatic1);

		if (SP_GRAD_RUN_TEST(8.1)) test8<SpGradient, SpGradient, SpGradient, iNumRowsStatic2, iNumColsStatic2>(inumloops, inumnz, inumdof, iNumRowsStatic2, iNumColsStatic2);
		if (SP_GRAD_RUN_TEST(8.2)) test8<doublereal, doublereal, doublereal, iNumRowsStatic2, iNumColsStatic2>(inumloops, inumnz, inumdof, iNumRowsStatic2, iNumColsStatic2);
		if (SP_GRAD_RUN_TEST(8.3)) test8<SpGradient, doublereal, doublereal, iNumRowsStatic2, iNumColsStatic2>(inumloops, inumnz, inumdof, iNumRowsStatic2, iNumColsStatic2);
		if (SP_GRAD_RUN_TEST(8.4)) test8<doublereal, SpGradient, doublereal, iNumRowsStatic2, iNumColsStatic2>(inumloops, inumnz, inumdof, iNumRowsStatic2, iNumColsStatic2);
		if (SP_GRAD_RUN_TEST(8.5)) test8<doublereal, doublereal, SpGradient, iNumRowsStatic2, iNumColsStatic2>(inumloops, inumnz, inumdof, iNumRowsStatic2, iNumColsStatic2);
		if (SP_GRAD_RUN_TEST(8.6)) test8<SpGradient, SpGradient, doublereal, iNumRowsStatic2, iNumColsStatic2>(inumloops, inumnz, inumdof, iNumRowsStatic2, iNumColsStatic2);
		if (SP_GRAD_RUN_TEST(8.8)) test8<doublereal, SpGradient, SpGradient, iNumRowsStatic2, iNumColsStatic2>(inumloops, inumnz, inumdof, iNumRowsStatic2, iNumColsStatic2);
		if (SP_GRAD_RUN_TEST(8.9)) test8<SpGradient, doublereal, SpGradient, iNumRowsStatic2, iNumColsStatic2>(inumloops, inumnz, inumdof, iNumRowsStatic2, iNumColsStatic2);

		if (SP_GRAD_RUN_TEST(10.1)) test10<SpGradient, SpGradient>(inumloops, inumnz, inumdof, imatrows, imatcols, imatcolsb);
		if (SP_GRAD_RUN_TEST(10.2)) test10<doublereal, doublereal>(inumloops, inumnz, inumdof, imatrows, imatcols, imatcolsb);
		if (SP_GRAD_RUN_TEST(10.3)) test10<SpGradient, doublereal>(inumloops, inumnz, inumdof, imatrows, imatcols, imatcolsb);
		if (SP_GRAD_RUN_TEST(10.4)) test10<doublereal, SpGradient>(inumloops, inumnz, inumdof, imatrows, imatcols, imatcolsb);

		if (SP_GRAD_RUN_TEST(11.1)) test11<SpGradient, SpGradient, SpGradient>(inumloops, inumnz, inumdof, imatrows);
		if (SP_GRAD_RUN_TEST(11.2)) test11<doublereal, doublereal, doublereal>(inumloops, inumnz, inumdof, imatrows);
		if (SP_GRAD_RUN_TEST(11.3)) test11<SpGradient, doublereal, doublereal>(inumloops, inumnz, inumdof, imatrows);
		if (SP_GRAD_RUN_TEST(11.4)) test11<doublereal, SpGradient, doublereal>(inumloops, inumnz, inumdof, imatrows);
		if (SP_GRAD_RUN_TEST(11.5)) test11<doublereal, doublereal, SpGradient>(inumloops, inumnz, inumdof, imatrows);
		if (SP_GRAD_RUN_TEST(11.6)) test11<SpGradient, SpGradient, doublereal>(inumloops, inumnz, inumdof, imatrows);
		if (SP_GRAD_RUN_TEST(11.7)) test11<doublereal, SpGradient, SpGradient>(inumloops, inumnz, inumdof, imatrows);
		if (SP_GRAD_RUN_TEST(11.8)) test11<SpGradient, doublereal, SpGradient>(inumloops, inumnz, inumdof, imatrows);

		if (SP_GRAD_RUN_TEST(12.1)) test12<SpGradient, SpGradient, SpGradient>(inumloops, inumnz, inumdof, imatrows, imatcols, imatcolsb, imatcolsc);
		if (SP_GRAD_RUN_TEST(12.2)) test12<doublereal, doublereal, doublereal>(inumloops, inumnz, inumdof, imatrows, imatcols, imatcolsb, imatcolsc);
		if (SP_GRAD_RUN_TEST(12.3)) test12<SpGradient, doublereal, doublereal>(inumloops, inumnz, inumdof, imatrows, imatcols, imatcolsb, imatcolsc);
		if (SP_GRAD_RUN_TEST(12.4)) test12<doublereal, SpGradient, doublereal>(inumloops, inumnz, inumdof, imatrows, imatcols, imatcolsb, imatcolsc);
		if (SP_GRAD_RUN_TEST(12.5)) test12<doublereal, doublereal, SpGradient>(inumloops, inumnz, inumdof, imatrows, imatcols, imatcolsb, imatcolsc);
		if (SP_GRAD_RUN_TEST(12.6)) test12<SpGradient, SpGradient, doublereal>(inumloops, inumnz, inumdof, imatrows, imatcols, imatcolsb, imatcolsc);
		if (SP_GRAD_RUN_TEST(12.7)) test12<doublereal, SpGradient, SpGradient>(inumloops, inumnz, inumdof, imatrows, imatcols, imatcolsb, imatcolsc);
		if (SP_GRAD_RUN_TEST(12.8)) test12<SpGradient, doublereal, SpGradient>(inumloops, inumnz, inumdof, imatrows, imatcols, imatcolsb, imatcolsc);

		if (SP_GRAD_RUN_TEST(13.1)) test13<doublereal>(inumloops, inumnz, inumdof);
		if (SP_GRAD_RUN_TEST(13.1)) test13<SpGradient>(inumloops, inumnz, inumdof);

		if (SP_GRAD_RUN_TEST(14.1)) test_bool1(inumloops, inumnz, inumdof);
		if (SP_GRAD_RUN_TEST(15.1)) test15<doublereal>(inumloops, inumnz, inumdof);
		if (SP_GRAD_RUN_TEST(15.2)) test15<SpGradient>(inumloops, inumnz, inumdof);

#ifdef SP_GRAD_DEBUG
		cerr << "remaining references for SpGradient::oNullData: " << SpGradient::iGetRefCntNullData() << endl;
		cerr << "remaining references for SpMatrixBase::oNullData: " << SpMatrixBaseData::iGetRefCntNullData() << endl;
#endif
		SP_GRAD_ASSERT(SpGradient::iGetRefCntNullData() == 1);
		SP_GRAD_ASSERT(SpMatrixBaseData::iGetRefCntNullData() == 1);
		
		cerr << "All tests passed\n"
		     << "\n\tloops performed: " << inumloops
		     << "\n\tmax nonzeros: " << inumnz
		     << "\n\tmax dof: " << inumdof
		     << "\n\tmatrix size: " << imatrows << " x " << imatcols << " x " << imatcolsb << " x " << imatcolsc
		     << endl;

		return 0;
	} catch (const std::exception& err) {
		cerr << "an exception occured: " << err.what() << endl;
		return 1;
	}
}
