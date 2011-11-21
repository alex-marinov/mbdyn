/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2004-2011
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
neq:		is the matrix size;
a:		is used to store matrix elements; a[row][col] is the element in row
  		row and column col;
nzr and nzc:	are vectors of size neq, used to count the nonzero elements of a particular row or
		column, respectively;
ri and ci:	are neq x neq matrices used to store nonzero element indices; ri[col][i]
		(resp. ci[row][i]), with i < nzr[col] (i < nzc[row]), is the row (column)
		index of one of the nzr[col] (nzc[row]) nonzero elements in
		column col (row row); note that indices in ri[col] and ci[row] are
		not ordered;
nz:		nz[row][col] is true if the element in row row and column col is
		nonzero, false otherwise;
piv:		is a vector of size neq. 


The subroutine naivfct perform the LU factorization, naivslv the back-solve.

*/


#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <math.h>
#include "ac/f2c.h"

#else /* !HAVE_CONFIG_H */
/* to ease compilation outside of MBDyn...
 * replace long and double with the preferred types */
#include <math.h>
typedef long int integer;
typedef double doublereal;
#endif /* !HAVE_CONFIG_H */

#include <stdlib.h>
#include <limits.h>

#include "mthrdslv.h"

#define BIGINT   (1 << 30)

#define MINPIV   1.0e-8

int
naivfct(RMAT a, integer neq, integer *nzr, IMAT ri,
		integer *nzc, IMAT ci, NZMAT nz, 
		integer *piv, doublereal minpiv)
{
	char todo[neq];
	integer i, j, k, pvr, pvc, nr, nc, r;
	integer *pri, *pci;
	char *prik;
	doublereal den, mul, mulpiv, fari;
#if 0
	doublereal fapvr;
#endif
	doublereal *par, *papvr;

	if (neq <= 0 || (unsigned long)neq > NAIVE_MAX) {
		return NAIVE_ERANGE;
	}

	if (!minpiv) {
		minpiv = MINPIV;
	}
	for (pvr = 0; pvr < neq; pvr++) {
		todo[pvr] = 1;
	}
	for (i = 0; i < neq; i++) {
		if (!nzr[i]) { return NAIVE_ENULCOL + i; }
		nc = neq + 1;	
		nr = nzr[i];
		mul = 0.0;
		pri = ri[i];
		pvr = pri[0];
		mulpiv = 0.;
#if 0
		fapvr = 0.;
 		for (k = 0; k < nr; k++) {
 			r = pri[k];
 			if (todo[r]) {
 				fari = fabs(a[r][i]);
 				if (fari > mul) {
 					mulpiv = fari*minpiv;
 					if (nzc[r] <= nc  || mulpiv > fapvr) {
 						nc = nzc[pvr = r];
 						fapvr = fari;
 					}
 					mul = fari;
 				} else if (nzc[r] < nc && fari > mulpiv) {
 					nc = nzc[pvr = r];
 				}
 			}
 		}
#endif
		for (k = 0; k < nr; k++) {
			r = pri[k];
			if (todo[r]) {
				fari = fabs(a[r][i]);
				if (fari > mul) {
					mul = fari;
				}
			}
		}
		mulpiv = mul*minpiv;
		for (k = 0; k < nr; k++) {
			r = pri[k];
			if (todo[r]) {
				fari = fabs(a[r][i]);
				if (fari >= mulpiv && nzc[r] < nc) {
					nc = nzc[pvr = r];
				}
			}
		}
		if (nc == neq + 1) { return NAIVE_ENOPIV + i; }
		if (mulpiv == 0.)  { return NAIVE_ENOPIV + i; }

		piv[i] = pvr;
		todo[pvr] = 0;
		papvr = a[pvr];
		den = papvr[i] = 1.0/papvr[i];

		for (k = 0; k < nr; k++) {
			if (!todo[r = pri[k]]) { continue; }
			par = a[r];
			mul = par[i] = par[i]*den;
			prik = nz[r];
			pci = ci[pvr];
			for (j = 0; j < nc; j++) {
				if ((pvc = pci[j]) <= i) { continue; }
				if (prik[pvc]) {
					par[pvc] -= mul*papvr[pvc];
				} else {
					par[pvc] = -mul*papvr[pvc];
					prik[pvc] = 1;
					ri[pvc][nzr[pvc]++] = r;
					ci[r][nzc[r]++] = pvc;
				}
			}
		}
	}
	return 0;
}

/*
 * to solve A * x = b
 *
 * actually solve P * A * x = P *b
 *
 * compute P * L and P * U such that P * A = P * L * P^-1 * P * U
 *
 * and store P * L, P * U
 *
 * first step: P * L * f = P * b (but actually store P * f)
 *
 * second step: P * U * x = P * f
 */
int
naivslv(RMAT a, integer neq, integer *nzc, IMAT ci, 
		doublereal *rhs, doublereal * sol,
		integer *piv)
{
	doublereal fwd[neq];

	integer i, k, nc, r, c;
	integer *pci;
	doublereal s;
	doublereal *par;

	if (neq <= 0 || (unsigned long)neq > NAIVE_MAX) {
		return NAIVE_ERANGE;
	}

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
	sol[neq] = fwd[neq]*a[r][neq];
	for (i = neq - 1; i >= 0; i--) {
		nc = nzc[r = piv[i]];
		s = fwd[i];
		par = a[r];
		pci = ci[r];
		for (k = 0; k < nc; k++) {
			if ((c = pci[k]) > i) {
				s -= par[c]*sol[c];
			}
		}
		sol[i] = s*par[i];
	}

	return 0;
}

