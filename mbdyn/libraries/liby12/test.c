#include <stdlib.h>
#include <stdio.h>

#include <myf2c.h>
#include <y12lib.h>

int
main(void)
{
	const integer N = 2;	/* size */
	const integer Z = 3;	/* non-zeroes */

	const integer NN = 3*Z;	/* size of A and SNR */
	doublereal A[NN];	/* matrix */
	integer SNR[NN];	/* column indices */

	const integer NN1 = 3*Z;	/* size of RNR */
	integer RNR[NN1];	/* row indices */

	const integer IHA = N;	/* size of HA */
	integer HA[IHA*11];	/* work array IHA*11 */

	doublereal AFLAG[8];	/* ? */
	integer IFLAG[10];	/* ? */

	integer IFAIL = 0;	/* diagnostics */

	doublereal PIVOT[N];	/* pivot factors */
	doublereal B[N];	/* RHS of problem */

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

	y12mbf_(&N, &Z, A, SNR, &NN, RNR, &NN1, HA, &IHA, AFLAG, IFLAG, &IFAIL);


	B[0] = 1.;
	B[1] = 1.;

	y12mcf_(&N, &Z, A, SNR, &NN, RNR, &NN1, PIVOT, B, HA, &IHA, AFLAG, IFLAG, &IFAIL);

	y12mdf_(&N, A, &NN, B, PIVOT, SNR, HA, &IHA, IFLAG, &IFAIL);

	fprintf(stdout, "B = {%e,%e}\n", B[0], B[1]);

	IFLAG[4] = 3;

	B[0] = 2.;
	B[1] = 2.;

	y12mdf_(&N, A, &NN, B, PIVOT, SNR, HA, &IHA, IFLAG, &IFAIL);

	fprintf(stdout, "B = {%e,%e}\n", B[0], B[1]);

	B[0] = 2.;
	B[1] = -1.;

	y12mdf_(&N, A, &NN, B, PIVOT, SNR, HA, &IHA, IFLAG, &IFAIL);

	fprintf(stdout, "B = {%e,%e}\n", B[0], B[1]);

	return 0;
}
