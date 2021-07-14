/* $Header$
 * Copyright (C) 2007
 * Pierangelo Masarati
 * GPL'd 
 */

#ifndef UDU_H
#define UDU_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <math.h>
#include "matrix.h"

typedef double REAL;

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

extern int
udu_factor(int n, vector *a, vector *udu);

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

extern int
udu_solve(int n, vector *udu, vector *b);

int LinearSystemV( matrix *, vector *, vector *, vector*);
int LinearSystemM( matrix *, matrix *, matrix *, vector *, vector * );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* ! UDU_H */
