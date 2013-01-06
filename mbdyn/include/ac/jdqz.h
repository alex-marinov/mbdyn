/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

#ifndef AC_JDQZ_H
#define AC_JDQZ_H

#ifdef USE_JDQZ

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */

#include "ac/f2c.h"

/*
 * Add declarations of jdqz routines used by MBDyn
 */

extern int
__FC_DECL__(jdqz)(
	doublecomplex *alpha,
	doublecomplex *beta,
	doublecomplex *eivec,
	logical *wanted,
	integer *n,
	doublecomplex *target,
	doublereal *eps,
	integer *kmax,
	integer *jmax,
	integer *jmin,
	integer *method,
	integer *m,
	integer *l,
	integer *mxmv,
	integer *maxstep,
	doublereal *lock,
	integer *order,
	integer *testspace,
	doublecomplex *zwork,
	integer *lwork);

/*
The parameters must be of the following data types:
 alpha, beta   double complex array, size jmax
 eivec         two dimensional double complex array, size n × kmax
 wanted        logical, scalar
 n             integer, scalar
 target        double complex, scalar
 eps           double precision, scalar
 kmax          integer, scalar
 jmax          integer, scalar
 jmin          integer, scalar
 method        integer, scalar
 m             integer, scalar
 l             integer, scalar
 mxmv          integer, scalar
 maxstep       integer, scalar
 lock          double precision, scalar
 order         integer, scalar
 testspace     integer, scalar
 zwork         two dimensional double complex array, size n × lwork
 lwork         integer, scalar
*/

#if defined(__cplusplus)
}
#endif /* __cplusplus */

#endif /* USE_JDQZ */

#endif // AC_JDQZ_H
