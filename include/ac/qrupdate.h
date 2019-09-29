/* $Header: /home/reinhard/CVS/mbdyn/include/ac/qrupdate.h,v 2.1 2019/08/05 11:40:08 reinhard Exp $ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2019
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
  AUTHOR: Reinhard Resch <r.resch@a1.net>
  Copyright (C) 2019(-2019) all rights reserved.

  The copyright of this code is transferred
  to Pierangelo Masarati and Paolo Mantegazza
  for use in the software MBDyn as described
  in the GNU Public License version 2.1
*/

#ifndef AC_QRUPDATE_H
#define AC_QRUPDATE_H

#if defined(USE_QRUPDATE)

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */

#include "ac/f2c.h"

/*
  Add declarations of qrupdate routines used by MBDyn

  from function-reference.gz distributed with libqrupdate:

  dqr1up(m,n,k,Q,ldq,R,ldr,u,v,w)
	 purpose:      updates a QR factorization after rank-1 modification
	               i.e., given a m-by-k orthogonal Q and m-by-n upper
	               trapezoidal R, an m-vector u and n-vector v,
	               this subroutine updates Q -> Q1 and R -> R1 so that
	               Q1*R1 = Q*R + u*v', and Q1 is again orthonormal
	               and R1 upper trapezoidal.
	               (real version)
	 arguments:
	 m (in)        number of rows of the matrix Q.
	 n (in)        number of columns of the matrix R.
	 k (in)        number of columns of Q, and rows of R. Must be
	               either k = m (full Q) or k = n < m (economical form).
	 Q (io)        on entry, the orthogonal m-by-k matrix Q.
	               on exit, the updated matrix Q1.
	 ldq (in)      the leading dimension of Q. ldq >= m.
	 R (io)        on entry, the upper trapezoidal m-by-n matrix R..
	               on exit, the updated matrix R1.
	 ldr (in)      the leading dimension of R. ldr >= k.
	 u (io)        the left m-vector. On exit, if k < m, u is destroyed.
	 v (io)        the right n-vector. On exit, v is destroyed.
	 w (out)       a workspace vector of size 2*k
*/
        
extern int
__FC_DECL__(dqr1up)(const integer* m,
                    const integer* n,
                    const integer* k,
                    doublereal* Q,
                    const integer* ldq,
                    doublereal* R,
                    const integer* ldr,
                    doublereal* u,
                    doublereal* v,
                    doublereal* w);
      
#if defined(__cplusplus)
}
#endif /* __cplusplus */

#endif /* USE_QRUPDATE */

#endif // AC_QRUPDATE_H
