#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define  SPRSPIV

// questi 2 non sono colpa mia, li ha chiesti PierMasa
typedef int     integer;
typedef double  doublereal;

#define BIGINT   (1 << 30)

#define ENULCOL  1000000
#define ENOPIV   2000000

#define HIGH  0xF0000000
#define LOW   0x0FFFFFFF

#define MINPIV   1.0e-8

int naivfct(void *mat_a, integer neq, integer *vet_nzr, void *mat_ri, integer *vet_nzc, void *mat_ci, integer ncd, integer *vet_piv, doublereal minpiv)
{
        struct { doublereal mat[neq][ncd]; } *_mat_a_  = mat_a;
        struct   { integer mat[neq][ncd]; } *_mat_ri_ = mat_ri;
        struct   { integer mat[neq][ncd]; } *_mat_ci_ = mat_ci;

        #define  a(i, k) (_mat_a_->mat)[i][k]
        #define ri(i, k) (_mat_ri_->mat)[i][k]
        #define ci(i, k) (_mat_ci_->mat)[i][k]
        #define   nzr(i) vet_nzr[i]
        #define   nzc(i) vet_nzc[i]
        #define   piv(i) vet_piv[i]

	integer vet_todo[neq];
        #define todo(i) vet_todo[i]

	integer i, j, k, m, pvr, pvc, nr, nc, r;
	doublereal den, mul;

	if (!minpiv) {
		minpiv = MINPIV;
	}
	for (pvr = 0; pvr < neq; pvr++) {
		piv(pvr) = todo(pvr) = -1;
	}
	for (i = 0; i < neq; i++) {
		if (!nzr(i)) { return -(ENULCOL + i); }
		m = BIGINT;	
		nr = nzr(i);
		mul = 0.0;
		for (k = 0; k < nr; k++) {
			r = ri(i, k) & LOW;
			if (todo(r) && fabs(a(r, i)) > mul) {
				mul = fabs(a(r, i));
				m = nzc(pvr = r);	
			}
		}
		if (m == BIGINT) { return -(ENOPIV + i); }
#ifdef SPRSPIV
		mul *= minpiv;
		for (k = 0; k < nr; k++) {
			r = ri(i, k) & LOW;
			if (todo(r) && nzc(r) < m && fabs(a(r, i)) > mul) {
				m = nzc(pvr = r);	
			}
		}
#endif
		piv(i) = pvr;
		todo(pvr) = 0;
		den = a(pvr, i) = 1.0/a(pvr, i);
		nr = nzr(i);
		nc = nzc(pvr);
		for (k = 0; k < nr; k++) {
			if (!todo(r = ri(i, k) & LOW)) { continue; }
			mul = a(r, i) = a(r, i)*den;
			for (j = 0; j < nc; j++) {
				if ((pvc = ci(pvr, j)) <= i) { continue; }
				if (ri(r, pvc) & HIGH) {
					a(r, pvc) -= mul*a(pvr, pvc);
				} else {
					a(r, pvc) = -mul*a(pvr, pvc);
					ri(r, pvc) |= HIGH;
					ri(pvc, nzr(pvc)++) |= r;
					ci(r, nzc(r)++) = pvc;
				}
			}
		}
	}
	return 0;
}

void naivslv(void *mat_a, integer neq, integer *vet_nzc, void *mat_ci, integer ncd, doublereal *vet_rhs, integer *vet_piv)
{
        struct { doublereal mat[neq][ncd]; } *_mat_a_  = mat_a;
        struct     { integer mat[neq][ncd]; } *_mat_ci_ = mat_ci;

        #define  a(i, k) (_mat_a_->mat)[i][k]
        #define ci(i, k) (_mat_ci_->mat)[i][k]
        #define   rhs(i) vet_rhs[i]
        #define   piv(i) vet_piv[i]

	doublereal  vet_fwd[neq];
        #define fwd(i) vet_fwd[i]

	integer i, k, nc, r, c;
	doublereal s;

	fwd(0) = rhs(piv(0));
	for (i = 1; i < neq; i++) {
		nc = nzc(r = piv(i));
		s = rhs(r);
		for (k = 0; k < nc; k++) {
			if ((c = ci(r, k)) < i) {
				s -= a(r, c)*fwd(c);
			}
		}
		fwd(i) = s;
	}

	r = piv(--neq);
	rhs(neq) = fwd(neq)*a(r, neq);
	for (i = neq - 1; i >= 0; i--) {
		nc = nzc(r = piv(i));
		s = fwd(i);
		for (k = 0; k < nc; k++) {
			if ((c = ci(r, k)) > i) {
				s -= a(r, c)*rhs(c);
			}
		}
		rhs(i) = s*a(r, i);
	}

	return;
}

#if 0

/********************************** TEST **************************************/

//#define CPUFRQ  1000000000
#define CPUFRQ  1800000000
static inline unsigned long long rd_CPU_ts(void)
{
	unsigned long long time;
	__asm__ __volatile__( "rdtsc" : "=A" (time));
	return time;
}

#define NEQ       120
#define HALFBAND  12
#define ACTIVCOL  0
#define MINPIV    1.0e-3

static doublereal a[NEQ][NEQ];
//static doublereal a[NEQ][NEQ] = { { 17, 3, 0, 0, -7, 0 }, { -5, 0, 13, 19, 0, 0 }, { 7, 0, -17, 10, 0, 1 }, { 6, -29, 0, -2, 0, 1 }, { 0, -2, 0, 0, 3, 0 }, { 9, -0, 0, 0, 11, 0 } };
//static doublereal a[NEQ][NEQ] = { { 1, 2, 0, 0, 7 }, { 6, 0, 1, 0, 0, }, { 7, 0, 0, 10, 0 }, { 6, 0, 0, -2, 0 }, { 6, -2, 0, 0, 3 } };
//static doublereal a[NEQ][NEQ] = { { 1, 2, 4, 0 }, { 6, 0, 1, 8, }, { 7, 0, 0, 10 }, { 6, -2, 0, 0 } };
//static doublereal a[NEQ][NEQ] = { { 1, 2, 4, 2 }, { 6, 5, 8, 9, }, { 7, 7, 9, 10 }, { 6, -2, 19, 12 } };
//static doublereal a[NEQ][NEQ] = { { 1 }, { 4, 5 }, { 7, 8, 9 } };
//static doublereal a[NEQ][NEQ] = { { 1, 2, 4 }, { 0, 4, 6 }, { 0, 0, 9 } };

int main(void)
{
	long long tf, ts;
	integer nzr[NEQ], ri[NEQ][NEQ], nzc[NEQ], ci[NEQ][NEQ];
	integer i, k, piv[NEQ], retval, nzb, nza;
	doublereal rhs[NEQ][NEQ];

	for (i = 0; i < NEQ; i++) {
		for (k = (i - HALFBAND) < 0 ? 0 : i - HALFBAND; k < ((i + HALFBAND) > NEQ ? NEQ : i + HALFBAND); k++) {
			a[i][k] = 2.0*(((doublereal)rand())/RAND_MAX - 0.5);
		}
	}
	for (i = NEQ - ACTIVCOL; i < NEQ; i++) {
		for (k = 0; k < NEQ; k++) {
			a[k][i] = a[i][k] = 2.0*(((doublereal)rand())/RAND_MAX - 0.5);
		}
	}
	nzb = 0;
	for (i = 0; i < NEQ*NEQ; i++) {
		if (a[0][i]) {
			nzb++;
		}
	}
	for (i = 0; i < NEQ; i++) {
		nzr[i] = nzc[i] = 0;
		for (k = 0; k < NEQ; k++) {
			rhs[i][k] = a[k][i];
		}
		for (k = 0; k < NEQ; k++) {
			if (a[k][i]) {
				ri[i][nzr[i]++] = k;
			}
		}
		for (k = 0; k < NEQ; k++) {
			if (a[i][k]) {
				ci[i][nzc[i]++] = k;
			}
		}
	}

#if 0
printf("ROW\n");
for (i = 0; i < NEQ; i++) {
	printf("%d ", nzr[i]);
	for (k = 0; k < nzr[i]; k++) {
		printf("%d ", ri[i][k]);
	}
	printf("\n");
}
printf("COL\n");
for (i = 0; i < NEQ; i++) {
	printf("%d ", nzc[i]);
	for (k = 0; k < nzc[i]; k++) {
		printf("%d ", ci[i][k]);
	}
	printf("\n");
}
printf("MAT\n");
for (i = 0; i < NEQ; i++) {
	for (k = 0; k < NEQ; k++) {
		printf("%f ", a[i][k]);
	}
	printf("\n");
}
printf("\n");
#endif

	tf = rd_CPU_ts();
	retval = naivfct(a, NEQ, nzr, ri, nzc, ci, NEQ, piv, MINPIV);
	tf = rd_CPU_ts() - tf;
	nza = 0;
	for (i = 0; i < NEQ*NEQ; i++) {
		if (a[0][i]) {
			nza++;
		}
	}

#if 1
printf("RES\n");
	for (i = NEQ - 1; i < NEQ; i++) {
		ts = rd_CPU_ts();
		naivslv(a, NEQ, nzc, ci, NEQ, rhs[i], piv);
		ts = rd_CPU_ts() - ts;
for (k = 0; k < NEQ; k++) {
	printf("%f ", rhs[i][k]);
}
printf("\nRETVAL: %d, FILBEF: %d, FILAFT: %d\n", retval, (int)(100.0*nzb/(float)(NEQ*NEQ)), (int)(100.0*nza/(float)(NEQ*NEQ)));
	}
#endif

printf("\nFACT TIME: %d, SOL TIME: %d - (us)\n", (int)(tf*1000000.0/CPUFRQ), (int)(ts*1000000.0/CPUFRQ));

	return 0;
}

#endif /* 0 */

