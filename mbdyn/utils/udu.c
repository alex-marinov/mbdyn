/* $Header$
 * Copyright (C) 2007
 * Pierangelo Masarati
 * GPL'd 
 */

#include <math.h>
#include "udu.h"

/* assume a[] is stored in upper diagonal form:
 *
 *	a[0] == a_11
 *	a[1] == a_12
 *	a[2] == a_22
 *	a[3] == a_13
 *	a[4] == a_23
 *	a[5] == a_33
 *	a[6] == a_14
 *	...
 *	a[(nc-1)*nc/2 + (nr-1)] = a_nr_nc
 */

/* U D U' */

/* u is stored in strictly upper diagonal form,
 * with d on the diagonal, in place */

int
udu_factor(int n, vector *a, vector *udu)
{
	int	c;

	for (c = n; c-- > 0; ) {
		int	r, k;

		/* column index */
		int	ci = c*(c+1)/2;

		for (k = c + 1; k < n; k++) {
			/* column index */
			int	ki = k*(k+1)/2;

			udu->vec[ci + c] -= udu->vec[ki + c]*udu->vec[ki + k]*udu->vec[ki + c];
		}

		/* d_cc must be non-zero (negative allowed) */
		if (fabs(udu->vec[ci + c]) < 1e-7) {
			return -(c + 1);
		}

		/* d_cc = p_cc */
		for (r = c; r-- > 0; ) {
			/* */
			for (k = c + 1; k < n; k++) {
				/* column index */
				int	ki = k*(k+1)/2;

				/* */
				udu->vec[ci + r] -= udu->vec[ki + r]*udu->vec[ki + k]*udu->vec[ki + c];
			}
			udu->vec[ci + r] /= udu->vec[ci + c];
		}
	}

	return 0;
}

/*
 * solves P x = b, with P decomposed as U D U'
 *
 * actually, solves U D U' x = b
 *
 * split in
 *
 * U' x = y
 * D y = z
 * U z = b
 *
 */

int
udu_solve(int n, vector *udu, vector *b)
{
	int	r;

	/* U z = b */
	for (r = n - 1; r-- > 0; ) {
		int	k;

		for (k = r + 1; k < n; k++) {
			b->vec[r] -= udu->vec[k*(k+1)/2 + r]*b->vec[k];
		}
	}

	/* D y = z */
	for (r = 0; r < n; r++) {
		b->vec[r] /= udu->vec[r*(r+1)/2 + r];
	}

	/* U' x = y */
	for (r = 1; r < n; r++) {
		int	k;

		/* column index */
		int	ri = r*(r+1)/2;

		for (k = 0; k < r; k++) {
			b->vec[r] -= udu->vec[ri + k]*b->vec[k];
		}
	}

	return 0;
}

static int matrix2upper_diagonal( matrix *MAT, vector *a){
	
	unsigned i,j;

	for( i=0; i<MAT->Nrow; i++){
		for( j=0; j<=i; j++){
			a->vec[i*(i+1)/2 + j] = MAT->mat[i][j];
		}
	}

	return 0;
}

int LinearSystemV( matrix *A, vector *b, vector *x, vector *aa){


	matrix2upper_diagonal( A, aa);			

	udu_factor(A->Nrow, aa, aa);
	udu_solve(A->Nrow, aa, b);

	vector_copy( x, b, 1.);


	return 0;
}

int LinearSystemM( matrix *A, matrix *b, matrix *x, vector *aa, vector *bb){

	unsigned i,j;

	matrix2upper_diagonal( A, aa);			
	udu_factor(A->Nrow, aa, aa);

	for( i=0; i<b->Ncolumn; i++){
		for( j=0; j<b->Nrow; j++){
			bb->vec[j] = b->mat[j][i];
		}
		udu_solve(A->Nrow, aa, bb);
		for( j=0; j<b->Nrow; j++){
			x->mat[j][i] = bb->vec[j];
		}

	}

	return 0;
}





#ifdef MAIN

#ifndef DIM
#define DIM 2
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int
main(void)
{
	int	i;
	matrix	A, bM ,xM;
	vector 	b, x;
	FILE *fhA, *fhb;

	matrix_init(&A, DIM, DIM);
	matrix_init(&bM, DIM, 5);
	matrix_init(&xM, DIM, 5);
	vector_init(&b, DIM);
	vector_init(&x, DIM);

#ifdef READ
	fhA = fopen("A.dat","r");
	fhb = fopen("b.dat","r");
	matrix_read( &A, fhA, 1 );
	//vector_read( &b, fhb, 1 );
	matrix_read( &bM, fhb, 1 );
#endif
#ifndef READ
	A.mat[0][0] = 1.247;
	A.mat[1][0] = 0.446;
	A.mat[0][1] = 0.446;
	A.mat[1][1] = 3.148;

	b.vec[0] = 1.;
	b.vec[1] = 1.;
#endif
	//vector_null(&x);
	matrix_null(&xM);

	//LinearSystemV( &A, &b, &x);
	LinearSystemM( &A, &bM, &xM);


	matrix_write(&A, stdout, W_M_TEXT);
	matrix_write(&bM, stdout, W_M_TEXT);
	matrix_write(&xM, stdout, W_M_TEXT);
	//vector_write(&b, stdout, W_M_TEXT);
	//vector_write(&x, stdout, W_M_TEXT);
}

#if 0
int
main(void)
{
	int	i;
	matrix	A;
	vector 	b, x;
	REAL	*aa;
	REAL	*bb;
	FILE *fhA, *fhb;

	matrix_init(&A, DIM, DIM);
	vector_init(&b, DIM);
	vector_init(&x, DIM);

#ifdef READ
	fhA = fopen("A.dat","r");
	fhb = fopen("b.dat","r");
	matrix_read( &A, fhA, 1 );
	vector_read( &b, fhb, 1 );
#endif
#ifndef READ
	A.mat[0][0] = 1.247;
	A.mat[1][0] = 0.446;
	A.mat[0][1] = 0.446;
	A.mat[1][1] = 3.148;

	b.vec[0] = 1.;
	b.vec[1] = 1.;
#endif
	vector_null(&x);

	aa = (REAL *)calloc( (A.Nrow*(A.Nrow+1))/2, sizeof(REAL) );
	bb = (REAL *)calloc( A.Nrow, sizeof(REAL) );

	matrix2upper_diagonal( &A, aa);			
	//printf("\n\n %e %e \n %e %e\n\n", aa[0], aa[1], aa[1], aa[2] );
	vector2array( &b, bb);
	//printf("\n\n %e \n %e \n\n", bb[0], bb[1] );
	
	udu_factor(A.Nrow, aa, aa);
	udu_solve(A.Nrow, aa, bb);

	array2vector( &x, bb );

	matrix_write(&A, stdout, W_M_TEXT);
	vector_write(&b, stdout, W_M_TEXT);
	vector_write(&x, stdout, W_M_TEXT);



/*
	memset(aa, 0, sizeof(aa));
	memset(bb, 0, sizeof(bb));

	aa[0] = A.mat[0][0];
	aa[1] = A.mat[0][1];
	aa[2] = A.mat[1][1];

	bb[0] = b.vec[0];
	bb[1] = b.vec[1];

	udu_factor(2, aa, aa);
	udu_solve(2, aa, bb);

	x.vec[0] = bb[0];
	x.vec[1] = bb[1];

	vector_write(&x, stdout, W_M_TEXT);
*/

	return 0;
}
#endif

#endif // MAIN
