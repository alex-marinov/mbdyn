#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define  SPRSPIV

#define BIGINT   (1 << 30)
#define ENULCOL  1000000
#define ENOPIV   2000000

#define HIGH     0x80000000
#define LOW      0x7FFFFFFF

#define MINPIV   1.0e-8

// questi 2 non sono colpa mia, li ha chiesti PierMasa
typedef int     integer;
typedef double  doublereal;

typedef integer** IMAT;
typedef doublereal** RMAT;

int naivfct(RMAT a, integer neq, integer *nzr, IMAT ri, integer *nzc, IMAT ci, integer *piv, doublereal minpiv)
{
	integer todo[neq];
	integer i, j, k, m, pvr, pvc, nr, nc, r;
	integer *pri, *prik, *pci;
	doublereal den, mul;
	doublereal *par, *papvr;

	if (!minpiv) {
		minpiv = MINPIV;
	}
	for (pvr = 0; pvr < neq; pvr++) {
		piv[pvr] = todo[pvr] = -1;
	}
	for (i = 0; i < neq; i++) {
		if (!nzr[i]) { return ENULCOL + i; }
		m = BIGINT;	
		nr = nzr[i];
		mul = 0.0;
		pri = ri[i];
		for (k = 0; k < nr; k++) {
			r = pri[k] & LOW;
			if (todo[r] && fabs(a[r][i]) > mul) {
				mul = fabs(a[r][i]);
				m = nzc[pvr = r];	
			}
		}
		if (m == BIGINT) { return ENOPIV + i; }
#ifdef SPRSPIV
		mul *= minpiv;
		for (k = 0; k < nr; k++) {
			r = pri[k] & LOW;
			if (todo[r] && nzc[r] < m && fabs(a[r][i]) > mul) {
				m = nzc[pvr = r];	
			}
		}
#endif
		piv[i] = pvr;
		todo[pvr] = 0;
		papvr = a[pvr];
		den = papvr[i] = 1.0/papvr[i];
		nr = nzr[i];
		nc = nzc[pvr];
		for (k = 0; k < nr; k++) {
			if (!todo[r = pri[k] & LOW]) { continue; }
			par = a[r];
			mul = par[i] = par[i]*den;
			prik = ri[r];
			pci = ci[pvr];
			for (j = 0; j < nc; j++) {
				if ((pvc = pci[j]) <= i) { continue; }
				if (prik[pvc] & HIGH) {
					par[pvc] -= mul*papvr[pvc];
				} else {
					par[pvc] = -mul*papvr[pvc];
					prik[pvc] |= HIGH;
					ri[pvc][nzr[pvc]++] |= r;
					ci[r][nzc[r]++] = pvc;
				}
			}
		}
	}
	return 0;
}

void naivslv(RMAT a, integer neq, integer *nzc, IMAT ci, doublereal *rhs, integer *piv)
{
	doublereal fwd[neq];

	integer i, k, nc, r, c;
	integer *pci;
	doublereal s;
	doublereal *par;

	fwd[0] = rhs[piv[0]];
	for (i = 1; i < neq; i++) {
		nc = nzc[r = piv[i]];
		s = rhs[r];
		par = a[r];
		pci = ci[r];
		for (k = 0; k < nc; k++) {
			if ((c = pci[k]) < i) {
				s -= par[c]*fwd[c];
			}
		}
		fwd[i] = s;
	}

	r = piv[--neq];
	rhs[neq] = fwd[neq]*a[r][neq];
	for (i = neq - 1; i >= 0; i--) {
		nc = nzc[r = piv[i]];
		s = fwd[i];
		par = a[r];
		pci = ci[r];
		for (k = 0; k < nc; k++) {
			if ((c = pci[k]) > i) {
				s -= par[c]*rhs[c];
			}
		}
		rhs[i] = s*par[i];
	}

	return;
}

