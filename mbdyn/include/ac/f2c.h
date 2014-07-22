/* $Header$ */
/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

#ifndef AC_F2C_H
#define AC_F2C_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#if defined(USE_G2C_H)
#include <g2c.h>
#define HAVE_FLAG_T
	
#elif defined(HAVE_F2C_H) /* !HAVE_G2C_H */
#include <f2c.h>
#define HAVE_FLAG_T

#else /* !HAVE_G2C_H && !HAVE_F2C_H */
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif /* HAVE_STDINT_H */
typedef int32_t integer;
typedef float real;
#ifdef MBDYN_SINGLE_PRECISION
#warning "ac/f2c.h: defining doublereal as float"
typedef float doublereal;
#else /* ! MBDYN_SINGLE_PRECISION */
typedef double doublereal;
#endif /* ! MBDYN_SINGLE_PRECISION */
#else /* !HAVE_SYS_TYPES_H */
#if defined(__alpha) || defined(__ia64)
#warning "ac/f2c.h: defining integer as int; edit as required"
typedef int integer;
#else /* !__alpha */
#warning "ac/f2c.h: defining integer as long int; edit as required"
typedef long int integer;
#endif /* !__alpha */
typedef float real;
#ifdef MBDYN_SINGLE_PRECISION
#warning "ac/f2c.h: defining doublereal as float"
typedef float doublereal;
#else /* ! MBDYN_SINGLE_PRECISION */
typedef double doublereal;
#endif /* ! MBDYN_SINGLE_PRECISION */
#endif /* !HAVE_SYS_TYPES_H */

typedef integer logical;
#if 0	/* we define flag somewhere else (we'll get rid of f2c some day!) */
typedef integer flag;
#endif
typedef integer ftnlen;
typedef integer ftnint;
typedef char *address;
typedef struct { doublereal r, i; } doublecomplex;

#endif /* !HAVE_G2C_H && !HAVE_F2C_H */

/*
 * Quoting gcc 3.2.1's cstdlib ...
 *
Get rid of those macros defined in <stdlib.h> in lieu of real functions.
 */
#undef abort
#undef abs
#undef atexit
#undef atof
#undef atoi
#undef atol
#undef bsearch
#undef calloc
#undef div
#undef exit
#undef free
#undef getenv
#undef labs
#undef ldiv
#undef malloc
#undef mblen
#undef mbstowcs
#undef mbtowc
#undef qsort
#undef rand
#undef realloc
#undef srand
#undef strtod
#undef strtol
#undef strtoul
#undef system
#undef wcstombs
#undef wctomb
  
#ifdef USE_UNDERSCORE
#define __FC_DECL__(arg) arg ##_
#else /* USE_UNDERSCORE */
#define __FC_DECL__(arg) arg
#endif /* USE_UNDERSCORE */

#ifdef __cplusplus
#undef max
#undef min
#endif /* __cplusplus */

/*
 * FIXME: non so se e' __FC_DECL__ o solo FALSE_
 * In g2c.h e' cosi':
 */
#define FALSE_ (0)
#define TRUE_ (1)

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* AC_F2C_H */

