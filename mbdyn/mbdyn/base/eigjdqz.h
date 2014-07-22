/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

3.1     User supplied subroutines
The user has to supply three problem dependent routines: one for the mul-
tiplication of a vector with the operator A, one for multiplication with B,
and one for performing the preconditioning operation. The subroutine to
multiply with A must be called AMUL and must have the following header:

      subroutine AMUL( n, q, r )
c...............................................
c...     Subroutine to compute r = Aq
c...............................................
      integer        n
      double complex q(n), r(n)

q is the input vector, r the output vector. n is the dimension of the problem.
The subroutine to multiply with B must be called BMUL and must have the
following header:

        subroutine BMUL( n, q, r )
c...............................................
c...        Subroutine to compute r = Bq
c...............................................
        integer            n
        double complex q(n), r(n)

Finally, the routine to perform the preconditioning operation must be called
PRECON and must have the header

        subroutine PRECON( n, q )
c...............................................
c...        Subroutine to compute q = K^-1 q
c...............................................
        integer            n
        double complex q(n)

The preconditioning matrix should be an approximation of the matrix A -
 B, with  the prechosen target value. Preconditioning within the JDQZ
algorithm is described in section 3.4 of [1]. Preconditioning is not essential
for the correct behavior of the algorithm. It should improve the rate of
convergence, but leaving the vector q untouched should have no influence on
the correctness of the results.

*/

#ifndef EIGJDQZ_H
#define EIGJDQZ_H

#ifdef USE_JDQZ

#include "ac/jdqz.h"
#include "naivemh.h"

extern "C" int
__FC_DECL__(amul)(integer *n, doublecomplex *q, doublecomplex *r);

extern "C" int
__FC_DECL__(bmul)(integer *n, doublecomplex *q, doublecomplex *r);

extern "C" int
__FC_DECL__(precon)(integer *n, doublecomplex *q);

class MBJDQZ {
protected:
	const NaiveMatrixHandler& A;
	const NaiveMatrixHandler& B;

	unsigned cnt;

	void Mul(integer n, doublecomplex *q, doublecomplex *r, const NaiveMatrixHandler& M);

public:
	MBJDQZ(const NaiveMatrixHandler& A, const NaiveMatrixHandler& B);
	virtual ~MBJDQZ(void);

	void AMul(integer n, doublecomplex *q, doublecomplex *r);
	void BMul(integer n, doublecomplex *q, doublecomplex *r);

	unsigned Cnt(void) const;
};

extern MBJDQZ* mbjdqzp;

#endif // USE_JDQZ

#endif // ! EIGJDQZ_H
