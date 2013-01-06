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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

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

#include "eigjdqz.h"

MBJDQZ* mbjdqzp;

extern "C" int
__FC_DECL__(amul)(integer *n, doublecomplex *q, doublecomplex *r)
{
	ASSERT(mbjdqzp != 0);

	mbjdqzp->AMul(*n, q, r);

	return 0;
}

extern "C" int
__FC_DECL__(bmul)(integer *n, doublecomplex *q, doublecomplex *r)
{
	ASSERT(mbjdqzp != 0);

	mbjdqzp->BMul(*n, q, r);

	return 0;
}

extern "C" int
__FC_DECL__(precon)(integer *n, doublecomplex *q)
{
	// leave untouched
	return 0;
}

MBJDQZ::MBJDQZ(const NaiveMatrixHandler& A, const NaiveMatrixHandler& B)
: A(A), B(B), cnt(0)
{
	NO_OP;
}

MBJDQZ::~MBJDQZ(void)
{
	NO_OP;
}

void
MBJDQZ::Mul(integer n, doublecomplex *q, doublecomplex *r, const NaiveMatrixHandler& M)
{
	ASSERT(n == M.iGetNumRows());
	ASSERT(n == M.iGetNumCols());

	for (integer i = 0; i < n; i++) {
		r[i].r = 0.;
		r[i].i = 0.;
	}

	for (NaiveMatrixHandler::const_iterator i = M.begin(); i != M.end(); ++i) {
		r[i->iRow].r += i->dCoef*q[i->iCol].r;
		r[i->iRow].i += i->dCoef*q[i->iCol].i;
	}
}

void
MBJDQZ::AMul(integer n, doublecomplex *q, doublecomplex *r)
{
#define CNT (100)
	if (!(cnt % CNT)) {
		silent_cerr("\r" "cnt=" << cnt);
	}

	cnt++;

	Mul(n, q, r, A);
}

void
MBJDQZ::BMul(integer n, doublecomplex *q, doublecomplex *r)
{
	Mul(n, q, r, B);
}

unsigned
MBJDQZ::Cnt(void) const
{
	return cnt;
}
