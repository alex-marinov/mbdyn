/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
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

#ifndef MYF2C_H
#define MYF2C_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifdef HAVE_F2C_H
#include <f2c.h>
#else /* !HAVE_F2C_H */

#ifndef min
#define min(a,b) ((a) < (b) ? (a) : (b))
#endif /* !min */
#ifndef max
#define max(a,b) ((a) > (b) ? (a) : (b))
#endif /* !max */


#ifdef __alpha
typedef int integer;
typedef float real;
typedef double doublereal;
#else /* !__alpha */
typedef long int integer;
typedef float real;
typedef double doublereal;
#endif /* !__alpha */

typedef integer logical;
typedef integer flag;
typedef integer ftnlen;
typedef integer ftnint;
#warning "Controllare"
typedef char *address;

#endif /* !HAVE_F2C_H */
   
#ifdef USE_UNDERSCORE
#define __FC_DECL__(arg) arg ##_
#else /* USE_UNDERSCORE */
#define __FC_DECL__(arg) arg
#endif /* USE_UNDERSCORE */

extern int finite(double);


/*
 * ??? non so se e' __FC_DECL__ o solo FALSE_
 * In g2c.h e' cosi':
 */
#warning "Controllare"
#define FALSE_ (0)
#define TRUE_ (1)


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* MYF2C_H */
