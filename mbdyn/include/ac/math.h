/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
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

#if defined(AC_MATH_H) && defined(__CYGWIN__)
#ifndef isfinite
#include <math.h>
#define isfinite(x)     finite(x)
#endif /* !isfinite */
#endif /* AC_MATH_H && __CYGWIN__ */

#ifndef AC_MATH_H
#define AC_MATH_H

/*use cmath only for c++ programs*/
#if defined(HAVE_STD_ISFINITE_IN_CMATH) && defined(__cplusplus)
#include <cmath>
using std::isfinite;
#else
#include <math.h>
#ifndef isfinite
/* Return nonzero value if X is not +-Inf or NaN.  */
#define isfinite(x)     finite(x)
#endif /* isfinite */
#endif /*HAVE_STD_ISFINITE_IN_CMATH && __cplusplus*/

#if !defined(HAVE_MATH_MACROS) && !defined(__USE_XOPEN) && !defined(__USE_BSD)

/* Return nonzero value if X is not +-Inf or NaN.  */
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
extern int finite(double);
extern double copysign(double x, double y);
#ifdef __cplusplus
}
#endif /* __cplusplus */
   
/* Some useful constants.  */
#define M_E            2.7182818284590452354   /* e */
#define M_LOG2E        1.4426950408889634074   /* log_2 e */
#define M_LOG10E       0.43429448190325182765  /* log_10 e */
#define M_LN2          0.69314718055994530942  /* log_e 2 */
#define M_LN10         2.30258509299404568402  /* log_e 10 */
#define M_PI           3.14159265358979323846  /* pi */
#define M_PI_2         1.57079632679489661923  /* pi/2 */
#define M_PI_4         0.78539816339744830962  /* pi/4 */
#define M_1_PI         0.31830988618379067154  /* 1/pi */
#define M_2_PI         0.63661977236758134308  /* 2/pi */
#define M_2_SQRTPI     1.12837916709551257390  /* 2/sqrt(pi) */
#define M_SQRT2        1.41421356237309504880  /* sqrt(2) */
#define M_SQRT1_2      0.70710678118654752440  /* 1/sqrt(2) */

/* The above constants are not adequate for computation using `long double's.
 * Therefore we provide as an extension constants with similar names as a
 * GNU extension.  Provide enough digits for the 128-bit IEEE quad.  */
#define M_El           2.7182818284590452353602874713526625L  /* e */
#define M_LOG2El       1.4426950408889634073599246810018922L  /* log_2 e */
#define M_LOG10El      0.4342944819032518276511289189166051L  /* log_10 e */
#define M_LN2l         0.6931471805599453094172321214581766L  /* log_e 2 */
#define M_LN10l        2.3025850929940456840179914546843642L  /* log_e 10 */
#define M_PIl          3.1415926535897932384626433832795029L  /* pi */
#define M_PI_2l        1.5707963267948966192313216916397514L  /* pi/2 */
#define M_PI_4l        0.7853981633974483096156608458198757L  /* pi/4 */
#define M_1_PIl        0.3183098861837906715377675267450287L  /* 1/pi */
#define M_2_PIl        0.6366197723675813430755350534900574L  /* 2/pi */
#define M_2_SQRTPIl    1.1283791670955125738961589031215452L  /* 2/sqrt(pi) */
#define M_SQRT2l       1.4142135623730950488016887242096981L  /* sqrt(2) */
#define M_SQRT1_2l     0.7071067811865475244008443621048490L  /* 1/sqrt(2) */

#endif /* define fancy constants */
   

#endif /* AC_MATH_H */

