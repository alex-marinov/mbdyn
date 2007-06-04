/* $Header$ */
/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2007
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

/* Return nonzero value if X is not +-Inf or NaN.  */
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
#ifndef HAVE_FINITE
extern int finite(double);
#endif /* HAVE_FINITE */
#ifndef HAVE_COPYSIGN
extern double copysign(double x, double y);
#endif /* HAVE_COPYSIGN */
#ifdef __cplusplus
}
#endif /* __cplusplus */
   
/* Some useful constants.  */
/* Original Copyright statement: */
/* Declarations for math functions.
   Copyright (C) 1991-1993,1995-1999,2001,2002,2004 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */

#ifndef M_E
#define M_E            2.7182818284590452354   /* e */
#endif /* ! M_E */
#ifndef M_LOG2E
#define M_LOG2E        1.4426950408889634074   /* log_2 e */
#endif /* ! M_LOG2E */
#ifndef M_LOG10E
#define M_LOG10E       0.43429448190325182765  /* log_10 e */
#endif /* ! M_LOG10E */
#ifndef M_LN2
#define M_LN2          0.69314718055994530942  /* log_e 2 */
#endif /* ! M_LN2 */
#ifndef M_LN10
#define M_LN10         2.30258509299404568402  /* log_e 10 */
#endif /* ! M_LN10 */
#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#endif /* ! M_PI */
#ifndef M_PI_2
#define M_PI_2         1.57079632679489661923  /* pi/2 */
#endif /* ! M_PI_2 */
#ifndef M_PI_4
#define M_PI_4         0.78539816339744830962  /* pi/4 */
#endif /* ! M_PI_4 */
#ifndef M_1_PI
#define M_1_PI         0.31830988618379067154  /* 1/pi */
#endif /* ! M_1_PI */
#ifndef M_2_PI
#define M_2_PI         0.63661977236758134308  /* 2/pi */
#endif /* ! M_2_PI */
#ifndef M_2_SQRTPI
#define M_2_SQRTPI     1.12837916709551257390  /* 2/sqrt(pi) */
#endif /* ! M_2_SQRTPI */
#ifndef M_SQRT2
#define M_SQRT2        1.41421356237309504880  /* sqrt(2) */
#endif /* ! M_SQRT2 */
#ifndef M_SQRT1_2
#define M_SQRT1_2      0.70710678118654752440  /* 1/sqrt(2) */
#endif /* ! M_SQRT1_2 */

/* The above constants are not adequate for computation using `long double's.
 * Therefore we provide as an extension constants with similar names as a
 * GNU extension.  Provide enough digits for the 128-bit IEEE quad.  */
#ifndef M_El
#define M_El           2.7182818284590452353602874713526625L  /* e */
#endif /* ! M_El */
#ifndef M_LOG2El
#define M_LOG2El       1.4426950408889634073599246810018922L  /* log_2 e */
#endif /* ! M_LOG2El */
#ifndef M_LOG10El
#define M_LOG10El      0.4342944819032518276511289189166051L  /* log_10 e */
#endif /* ! M_LOG10El */
#ifndef M_LN2l
#define M_LN2l         0.6931471805599453094172321214581766L  /* log_e 2 */
#endif /* ! M_LN2l */
#ifndef M_LN10l
#define M_LN10l        2.3025850929940456840179914546843642L  /* log_e 10 */
#endif /* ! M_LN10l */
#ifndef M_PIl
#define M_PIl          3.1415926535897932384626433832795029L  /* pi */
#endif /* ! M_PIl */
#ifndef M_PI_2l
#define M_PI_2l        1.5707963267948966192313216916397514L  /* pi/2 */
#endif /* ! M_PI_2l */
#ifndef M_PI_4l
#define M_PI_4l        0.7853981633974483096156608458198757L  /* pi/4 */
#endif /* ! M_PI_4l */
#ifndef M_1_PIl
#define M_1_PIl        0.3183098861837906715377675267450287L  /* 1/pi */
#endif /* ! M_1_PIl */
#ifndef M_2_PIl
#define M_2_PIl        0.6366197723675813430755350534900574L  /* 2/pi */
#endif /* ! M_2_PIl */
#ifndef M_2_SQRTPIl
#define M_2_SQRTPIl    1.1283791670955125738961589031215452L  /* 2/sqrt(pi) */
#endif /* ! M_2_SQRTPIl */
#ifndef M_SQRT2l
#define M_SQRT2l       1.4142135623730950488016887242096981L  /* sqrt(2) */
#endif /* ! M_SQRT2l */
#ifndef M_SQRT1_2l
#define M_SQRT1_2l     0.7071067811865475244008443621048490L  /* 1/sqrt(2) */
#endif /* ! M_SQRT1_2l */

#endif /* AC_MATH_H */

