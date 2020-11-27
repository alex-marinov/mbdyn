/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

#include "mbconfig.h"

#include <stdlib.h>
#include <stdio.h>
#include "ac/f2c.h"

#include "y12lib.h"

#ifdef USE_Y12

int
main(void)
{
	const integer N = 2;	/* size */
	const integer Z = 3;	/* non-zeroes */

	const integer NN = 3*Z;	/* size of A and SNR */
	doublereal A[9];	/* matrix NN = 3*Z */
	integer SNR[9];		/* column indices NN = 3*Z */

	const integer NN1 = 3*Z;/* size of RNR */
	integer RNR[9];		/* row indices NN1 = 3*Z */

	const integer IHA = N;	/* size of HA */
	integer HA[22];		/* work array IHA*11 */

	doublereal AFLAG[8] = {0., 0., 0., 0., 0., 0., 0., 0.};	/* ? */
	integer IFLAG[10];	/* ? */

	integer IFAIL = 0;	/* diagnostics */

	doublereal PIVOT[2];	/* pivot factors N */
	doublereal B[2];	/* RHS of problem N */

	A[0] = 2.;	/* 1, 1 */
	A[1] = 1.;      /* 2, 1 */
	A[2] = 1.;	/* 2, 2 */

	SNR[0] = 1;
	SNR[1] = 1;
	SNR[2] = 2;

	RNR[0] = 1;
	RNR[1] = 2;
	RNR[2] = 2;

	IFLAG[0] = 0;
	IFLAG[1] = 3;
	IFLAG[2] = 1;	/* total pivoting */
	IFLAG[3] = 0;
	IFLAG[4] = 2;

	__FC_DECL__(y12mbf)((integer *)&N, (integer *)&Z, A, SNR, 
			(integer *)&NN, RNR, 
			(integer *)&NN1, HA, (integer *)&IHA, AFLAG, IFLAG, 
			&IFAIL);


	B[0] = 1.;
	B[1] = 1.;
	printf("AA %d\n", NN1);
	__FC_DECL__(y12mcf)((integer *)&N, (integer *)&Z, A, SNR, 
			(integer *)&NN, RNR, 
			(integer *)&NN1, PIVOT, B, HA, (integer *)&IHA, AFLAG, 
			IFLAG, &IFAIL);

	__FC_DECL__(y12mdf)((integer *)&N, A, (integer *)&NN, B, PIVOT, SNR, 
			HA, 
			(integer *)&IHA, IFLAG, &IFAIL);

	fprintf(stdout, "B = {%e,%e}\n", B[0], B[1]);

	IFLAG[4] = 3;

	B[0] = 2.;
	B[1] = 2.;

	__FC_DECL__(y12mdf)((integer *)&N, A, (integer *)&NN, B, PIVOT, SNR, 
			HA, 
			(integer *)&IHA, IFLAG, &IFAIL);

	fprintf(stdout, "B = {%e,%e}\n", B[0], B[1]);

	B[0] = 2.;
	B[1] = -1.;

	__FC_DECL__(y12mdf)((integer *)&N, A, (integer *)&NN, B, PIVOT, SNR, 
			HA, 
			(integer *)&IHA, IFLAG, &IFAIL);

	fprintf(stdout, "B = {%e,%e}\n", B[0], B[1]);

	return 0;
}

#else /* ! USE_Y12 */

int
main(void)
{
	printf("configure --with-y12 to enable Y12\n");

	return 0;
}

#endif /* ! USE_Y12 */
