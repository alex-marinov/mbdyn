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
     struct SpGradAssignNoUpdateU {
     	  static constexpr void update_u(doublereal df_du, SpDerivRec* pFirstU, SpDerivRec* pLastU) {
	       // u* = f(u, v)
	       // f(u, v) = {u + v, u - v}
	       // df/du == 1
	       // u*' = u' + df/dv * v'	       
	  }
     };

     struct SpGradAssignUpdateU {
	  static constexpr void update_u(doublereal df_du, SpDerivRec* pFirstU, SpDerivRec* pLastU) {
	       // u* = f(u, v)
	       // f(u, v) = {u * v, u / v}
	       // u*' = df/du * u' + df/dv * v'
	       while (pFirstU < pLastU) {
		    pFirstU->dDer *= df_du;
		    ++pFirstU;
	       }
	  }
     };
     
     struct SpGradBinPlus: SpGradAssignNoUpdateU {
	  static constexpr doublereal f(doublereal u, doublereal v) {
	       return u + v;
	  }

	  static constexpr doublereal df_du(doublereal u, doublereal v) {
	       return 1.;
	  }

	  static constexpr doublereal df_dv(doublereal u, doublereal v) {
	       return 1.;
	  }
	  
#ifdef SP_GRAD_DEBUG
	  static void Print(std::ostream& os) {
	       os << "+";
	  }
#endif
     };

     struct SpGradBinMinus: SpGradAssignNoUpdateU {
	  static constexpr doublereal f(doublereal u, doublereal v) {
	       return u - v;
	  }

	  static constexpr doublereal df_du(doublereal u, doublereal v) {
	       return 1.;
	  }

	  static constexpr doublereal df_dv(doublereal u, doublereal v) {
	       return -1.;
	  }

#ifdef SP_GRAD_DEBUG
	  static void Print(std::ostream& os) {
	       os << "-";
	  }
#endif
     };

     struct SpGradBinMult: SpGradAssignUpdateU {
	  static constexpr doublereal f(doublereal u, doublereal v) {
	       return u * v;
	  }

	  static constexpr doublereal df_du(doublereal u, doublereal v) {
	       return v;
	  }

	  static constexpr doublereal df_dv(doublereal u, doublereal v) {
	       return u;
	  }
#ifdef SP_GRAD_DEBUG
	  static void Print(std::ostream& os) {
	       os << "*";
	  }
#endif
     };

     struct SpGradBinDiv: SpGradAssignUpdateU {
	  static constexpr doublereal f(doublereal u, doublereal v) {
	       return u / v;
	  }

	  static constexpr doublereal df_du(doublereal u, doublereal v) {
	       return 1. / v;
	  }

	  static constexpr doublereal df_dv(doublereal u, doublereal v) {
	       return -u / (v * v);
	  }
	  
#ifdef SP_GRAD_DEBUG
	  static void Print(std::ostream& os) {
	       os << "/";
	  }
#endif                
     };

     struct SpGradBinPow {
	  static constexpr doublereal f(doublereal u, doublereal v) {
	       return pow(u, v);
	  }

	  static constexpr doublereal df_du(doublereal u, doublereal v) {
	       return v * pow(u, v - 1.);
	  }

	  static constexpr doublereal df_dv(doublereal u, doublereal v) {
	       return pow(u, v) * log(u);
	  }
#ifdef SP_GRAD_DEBUG
	  static void Print(std::ostream& os) {
	       os << "pow";
	  }
#endif                
     };

     struct SpGradBinPowInt {
	  static constexpr doublereal f(doublereal u, integer v) {
	       return pow(u, v);
	  }

	  static constexpr doublereal df_du(doublereal u, integer v) {
	       return v * pow(u, v - 1);
	  }

	  static constexpr doublereal df_dv(doublereal u, integer v) {
	       return 0.;
	  }
#ifdef SP_GRAD_DEBUG
	  static void Print(std::ostream& os) {
	       os << "powi";
	  }
#endif                
     };

     struct SpGradBinAtan2 {
	  static constexpr doublereal f(doublereal u, doublereal v) {
	       return atan2(u, v);
	  }

	  static constexpr doublereal df_du(doublereal u, doublereal v) {
	       return v / (v * v + u * u);
	  }

	  static constexpr doublereal df_dv(doublereal u, doublereal v) {
	       return -u / (v * v + u * u);
	  }
                
#ifdef SP_GRAD_DEBUG
	  static void Print(std::ostream& os) {
	       os << "atan2";
	  }
#endif                
     };

     struct SpGradBinCopysign {
	  static doublereal f(doublereal u, doublereal v) {
	       return std::copysign(u, v);
	  }

	  static doublereal df_du(doublereal u, doublereal v) {
	       return std::copysign(1., u) * std::copysign(1., v);
	  }

	  static constexpr doublereal df_dv(doublereal u, doublereal v) {
	       return 0.;
	  }
                
#ifdef SP_GRAD_DEBUG
	  static void Print(std::ostream& os) {
	       os << "copysign";
	  }
#endif                
     };

     struct SpGradBinFmod {
	  static constexpr doublereal f(doublereal u, doublereal v) {
	       return fmod(u, v);
	  }

	  static constexpr doublereal df_du(doublereal u, doublereal v) {
	       return 1.;
	  }

	  static constexpr doublereal df_dv(doublereal u, doublereal v) {
	       return -int(u / v);
	  }

#ifdef SP_GRAD_DEBUG
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
#ifdef SP_GRAD_DEBUG
	  static void Print(std::ostream& os) {
	       os << "-";
	  }
#endif                
     };


     struct SpGradFabs {
	  static doublereal f(doublereal u) {
	       return fabs(u);
	  }

	  static doublereal df_du(doublereal u) {
	       return std::copysign(1., u);
	  }
                
#ifdef SP_GRAD_DEBUG
	  static void Print(std::ostream& os) {
	       os << "fabs";
	  }
#endif                
     };

     struct SpGradSqrt {
	  static constexpr doublereal f(doublereal u) {
	       return sqrt(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return 1. / (2. * sqrt(u));
	  }
#ifdef SP_GRAD_DEBUG
	  static void Print(std::ostream& os) {
	       os << "sqrt";
	  }
#endif                
     };

     struct SpGradExp {
	  static constexpr doublereal f(doublereal u) {
	       return exp(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return exp(u);
	  }
#ifdef SP_GRAD_DEBUG
	  static void Print(std::ostream& os) {
	       os << "exp";
	  }
#endif                
     };

     struct SpGradLog {
	  static constexpr doublereal f(doublereal u) {
	       return log(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return 1. / u;
	  }
#ifdef SP_GRAD_DEBUG
	  static void Print(std::ostream& os) {
	       os << "log";
	  }
#endif                
     };

     struct SpGradSin {
	  static constexpr doublereal f(doublereal u) {
	       return sin(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return cos(u);
	  }
#ifdef SP_GRAD_DEBUG
	  static void Print(std::ostream& os) {
	       os << "sin";
	  }
#endif                
     };

     struct SpGradCos {
	  static constexpr doublereal f(doublereal u) {
	       return cos(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return -sin(u);
	  }
#ifdef SP_GRAD_DEBUG
	  static void Print(std::ostream& os) {
	       os << "cos";
	  }
#endif                
     };

     struct SpGradTan {
	  static constexpr doublereal f(doublereal u) {
	       return tan(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return 1. + tan(u) * tan(u);
	  }
#ifdef SP_GRAD_DEBUG
	  static void Print(std::ostream& os) {
	       os << "tan";
	  }
#endif                
     };

     struct SpGradSinh {
	  static constexpr doublereal f(doublereal u) {
	       return sinh(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return cosh(u);
	  }
#ifdef SP_GRAD_DEBUG
	  static void Print(std::ostream& os) {
	       os << "sinh";
	  }
#endif                
     };

     struct SpGradCosh {
	  static constexpr doublereal f(doublereal u) {
	       return cosh(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return sinh(u);
	  }
#ifdef SP_GRAD_DEBUG
	  static void Print(std::ostream& os) {
	       os << "cosh";
	  }
#endif                
     };

     struct SpGradTanh {
	  static constexpr doublereal f(doublereal u) {
	       return tanh(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return 1. - tanh(u) * tanh(u);
	  }
#ifdef SP_GRAD_DEBUG
	  static void Print(std::ostream& os) {
	       os << "tanh";
	  }
#endif                
     };

     struct SpGradAsin {
	  static constexpr doublereal f(doublereal u) {
	       return asin(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return 1. / sqrt(1 - u * u);
	  }
#ifdef SP_GRAD_DEBUG
	  static void Print(std::ostream& os) {
	       os << "asin";
	  }
#endif                
     };

     struct SpGradAcos {
	  static constexpr doublereal f(doublereal u) {
	       return acos(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return -1. / sqrt(1 - u * u);
	  }
#ifdef SP_GRAD_DEBUG
	  static void Print(std::ostream& os) {
	       os << "acos";
	  }
#endif                
     };

     struct SpGradAtan {
	  static constexpr doublereal f(doublereal u) {
	       return atan(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return 1. / (1. + u * u);
	  }
#ifdef SP_GRAD_DEBUG
	  static void Print(std::ostream& os) {
	       os << "atan";
	  }
#endif                
     };

     struct SpGradAsinh {
	  static constexpr doublereal f(doublereal u) {
	       return asinh(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return 1. / sqrt(1. + u * u);
	  }
#ifdef SP_GRAD_DEBUG
	  static void Print(std::ostream& os) {
	       os << "asinh";
	  }
#endif                
     };

     struct SpGradAcosh {
	  static constexpr doublereal f(doublereal u) {
	       return acosh(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return 1. / sqrt(u * u - 1.);
	  }
#ifdef SP_GRAD_DEBUG
	  static void Print(std::ostream& os) {
	       os << "acosh";
	  }
#endif                
     };

     struct SpGradAtanh {
	  static constexpr doublereal f(doublereal u) {
	       return atanh(u);
	  }

	  static constexpr doublereal df_du(doublereal u) {
	       return 1. / (1. - u * u);
	  }
#ifdef SP_GRAD_DEBUG
	  static void Print(std::ostream& os) {
	       os << "atanh";
	  }
#endif                
     };

     struct SpGradBoolLessThan {
	  static constexpr bool f(doublereal u, doublereal v) {
	       return u < v;
	  }
     };

     struct SpGradBoolLessEqual {
	  static constexpr bool f(doublereal u, doublereal v) {
	       return u <= v;
	  }
     };

     struct SpGradBoolGreaterThan {
	  static constexpr bool f(doublereal u, doublereal v) {
	       return u > v;
	  }
     };

     struct SpGradBoolGreaterEqual {
	  static constexpr bool f(doublereal u, doublereal v) {
	       return u >= v;
	  }
     };

     struct SpGradBoolEqualTo {
	  static constexpr bool f(doublereal u, doublereal v) {
	       return u == v;
	  }
     };

     struct SpGradBoolNotEqualTo {     
	  static constexpr bool f(doublereal u, doublereal v) {
	       return u != v;
	  }
     };
}

#endif
