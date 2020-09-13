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

#ifndef __SP_GRADIENT_FUNCTIONS__
#define __SP_GRADIENT_FUNCTIONS__

#include <cmath>

#include "sp_gradient_base.h"

namespace sp_grad {
     class SpGradBinPlus {
     public:
	  static constexpr doublereal f(doublereal u, doublereal v) {
	       return u + v;
	  }

	  static constexpr doublereal df_du(doublereal u, doublereal v) {
	       return 1.;
	  }

	  static constexpr doublereal df_dv(doublereal u, doublereal v) {
	       return 1.;
	  }

#ifdef DEBUG
	  static void Print(std::ostream& os) {
	       os << "+";
	  }
#endif
     };

     class SpGradBinMinus {
     public:
	  static constexpr doublereal f(doublereal u, doublereal v) {
	       return u - v;
	  }

	  static constexpr doublereal df_du(doublereal u, doublereal v) {
	       return 1.;
	  }

	  static constexpr doublereal df_dv(doublereal u, doublereal v) {
	       return -1.;
	  }

#ifdef DEBUG
	  static void Print(std::ostream& os) {
	       os << "-";
	  }
#endif
     };

     class SpGradBinMult {
     public:
	  static constexpr doublereal f(doublereal u, doublereal v) {
	       return u * v;
	  }

	  static constexpr doublereal df_du(doublereal u, doublereal v) {
	       return v;
	  }

	  static constexpr doublereal df_dv(doublereal u, doublereal v) {
	       return u;
	  }

#ifdef DEBUG
	  static void Print(std::ostream& os) {
	       os << "*";
	  }
#endif
     };

     class SpGradBinDiv {
     public:
	  static constexpr doublereal f(doublereal u, doublereal v) {
	       return u / v;
	  }

	  static constexpr doublereal df_du(doublereal u, doublereal v) {
	       return 1. / v;
	  }

	  static constexpr doublereal df_dv(doublereal u, doublereal v) {
	       return -u / (v * v);
	  }
#ifdef DEBUG
	  static void Print(std::ostream& os) {
	       os << "/";
	  }
#endif                
     };

     class SpGradBinPow {
     public:
	  static constexpr doublereal f(doublereal u, doublereal v) {
	       return pow(u, v);
	  }

	  static constexpr doublereal df_du(doublereal u, doublereal v) {
	       return v * pow(u, v - 1.);
	  }

	  static constexpr doublereal df_dv(doublereal u, doublereal v) {
	       return pow(u, v) * log(u);
	  }
#ifdef DEBUG
	  static void Print(std::ostream& os) {
	       os << "pow";
	  }
#endif                
     };

     class SpGradBinPowInt {
     public:
	  static constexpr doublereal f(doublereal u, integer v) {
	       return pow(u, v);
	  }

	  static constexpr doublereal df_du(doublereal u, integer v) {
	       return v * pow(u, v - 1);
	  }

	  static constexpr doublereal df_dv(doublereal u, integer v) {
	       return 0.;
	  }
#ifdef DEBUG
	  static void Print(std::ostream& os) {
	       os << "powi";
	  }
#endif                
     };

     class SpGradBinAtan2 {
     public:
	  static constexpr doublereal f(doublereal u, doublereal v) {
	       return atan2(u, v);
	  }

	  static constexpr doublereal df_du(doublereal u, doublereal v) {
	       return v / (v * v + u * u);
	  }

	  static constexpr doublereal df_dv(doublereal u, doublereal v) {
	       return -u / (v * v + u * u);
	  }
                
#ifdef DEBUG
	  static void Print(std::ostream& os) {
	       os << "atan2";
	  }
#endif                
     };

     class SpGradBinCopysign {
     public:
	  static doublereal f(doublereal u, doublereal v) {
	       return std::copysign(u, v);
	  }

	  static doublereal df_du(doublereal u, doublereal v) {
	       return std::copysign(1., u) * std::copysign(1., v);
	  }

	  static constexpr doublereal df_dv(doublereal u, doublereal v) {
	       return 0.;
	  }
                
#ifdef DEBUG
	  static void Print(std::ostream& os) {
	       os << "copysign";
	  }
#endif                
     };

     class SpGradBinFmod {
     public:
	  static constexpr doublereal f(doublereal u, doublereal v) {
	       return fmod(u, v);
	  }

	  static constexpr doublereal df_du(doublereal u, doublereal v) {
	       return 1.;
	  }

	  static constexpr doublereal df_dv(doublereal u, doublereal v) {
	       return -int(u / v);
	  }

#ifdef DEBUG
	  static void Print(std::ostream& os) {
	       os << "fmod";
	  }
#endif                
     };

     struct SpGradUnaryMinus {
	  static constexpr doublereal f(doublereal u) {
	       return -u;
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return -1;
	  }
#ifdef DEBUG
	  static void Print(std::ostream& os) {
	       os << "-";
	  }
#endif                
     };


     class SpGradFabs {
     public:
	  static doublereal f(doublereal u) {
	       return fabs(u);
	  }

	  static doublereal df_du(doublereal u) {
	       return std::copysign(1., u);
	  }
                
#ifdef DEBUG
	  static void Print(std::ostream& os) {
	       os << "fabs";
	  }
#endif                
     };

     class SpGradSqrt {
     public:
	  static constexpr doublereal f(doublereal u) {
	       return sqrt(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return 1. / (2. * sqrt(u));
	  }
#ifdef DEBUG
	  static void Print(std::ostream& os) {
	       os << "sqrt";
	  }
#endif                
     };

     class SpGradExp {
     public:
	  static constexpr doublereal f(doublereal u) {
	       return exp(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return exp(u);
	  }
#ifdef DEBUG
	  static void Print(std::ostream& os) {
	       os << "exp";
	  }
#endif                
     };

     class SpGradLog {
     public:
	  static constexpr doublereal f(doublereal u) {
	       return log(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return 1. / u;
	  }
#ifdef DEBUG
	  static void Print(std::ostream& os) {
	       os << "log";
	  }
#endif                
     };

     class SpGradSin {
     public:
	  static constexpr doublereal f(doublereal u) {
	       return sin(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return cos(u);
	  }
#ifdef DEBUG
	  static void Print(std::ostream& os) {
	       os << "sin";
	  }
#endif                
     };

     class SpGradCos {
     public:
	  static constexpr doublereal f(doublereal u) {
	       return cos(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return -sin(u);
	  }
#ifdef DEBUG
	  static void Print(std::ostream& os) {
	       os << "cos";
	  }
#endif                
     };

     class SpGradTan {
     public:
	  static constexpr doublereal f(doublereal u) {
	       return tan(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return 1. + tan(u) * tan(u);
	  }
#ifdef DEBUG
	  static void Print(std::ostream& os) {
	       os << "tan";
	  }
#endif                
     };

     class SpGradSinh {
     public:
	  static constexpr doublereal f(doublereal u) {
	       return sinh(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return cosh(u);
	  }
#ifdef DEBUG
	  static void Print(std::ostream& os) {
	       os << "sinh";
	  }
#endif                
     };

     class SpGradCosh {
     public:
	  static constexpr doublereal f(doublereal u) {
	       return cosh(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return sinh(u);
	  }
#ifdef DEBUG
	  static void Print(std::ostream& os) {
	       os << "cosh";
	  }
#endif                
     };

     class SpGradTanh {
     public:
	  static constexpr doublereal f(doublereal u) {
	       return tanh(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return 1. - tanh(u) * tanh(u);
	  }
#ifdef DEBUG
	  static void Print(std::ostream& os) {
	       os << "tanh";
	  }
#endif                
     };

     class SpGradAsin {
     public:
	  static constexpr doublereal f(doublereal u) {
	       return asin(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return 1. / sqrt(1 - u * u);
	  }
#ifdef DEBUG
	  static void Print(std::ostream& os) {
	       os << "asin";
	  }
#endif                
     };

     class SpGradAcos {
     public:
	  static constexpr doublereal f(doublereal u) {
	       return acos(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return -1. / sqrt(1 - u * u);
	  }
#ifdef DEBUG
	  static void Print(std::ostream& os) {
	       os << "acos";
	  }
#endif                
     };

     class SpGradAtan {
     public:
	  static constexpr doublereal f(doublereal u) {
	       return atan(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return 1. / (1. + u * u);
	  }
#ifdef DEBUG
	  static void Print(std::ostream& os) {
	       os << "atan";
	  }
#endif                
     };

     class SpGradAsinh {
     public:
	  static constexpr doublereal f(doublereal u) {
	       return asinh(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return 1. / sqrt(1. + u * u);
	  }
#ifdef DEBUG
	  static void Print(std::ostream& os) {
	       os << "asinh";
	  }
#endif                
     };

     class SpGradAcosh {
     public:
	  static constexpr doublereal f(doublereal u) {
	       return acosh(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return 1. / sqrt(u * u - 1.);
	  }
#ifdef DEBUG
	  static void Print(std::ostream& os) {
	       os << "acosh";
	  }
#endif                
     };

     class SpGradAtanh {
     public:
	  static constexpr doublereal f(doublereal u) {
	       return atanh(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return 1. / (1. - u * u);
	  }
#ifdef DEBUG
	  static void Print(std::ostream& os) {
	       os << "atanh";
	  }
#endif                
     };

     class SpGradBoolLessThan {
     public:
	  static constexpr bool f(doublereal u, doublereal v) {
	       return u < v;
	  }
     };

     class SpGradBoolLessEqual {
     public:
	  static constexpr bool f(doublereal u, doublereal v) {
	       return u <= v;
	  }
     };

     class SpGradBoolGreaterThan {
     public:
	  static constexpr bool f(doublereal u, doublereal v) {
	       return u > v;
	  }
     };

     class SpGradBoolGreaterEqual {
     public:
	  static constexpr bool f(doublereal u, doublereal v) {
	       return u >= v;
	  }
     };

     class SpGradBoolEqualTo {
     public:
	  static constexpr bool f(doublereal u, doublereal v) {
	       return u == v;
	  }
     };

     class SpGradBoolNotEqualTo {
     public:
	  static constexpr bool f(doublereal u, doublereal v) {
	       return u != v;
	  }
     };
}

#endif
