/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2004-2005
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <stdlib.h>
#include <limits.h>

#include "ac/math.h"
#include "ac/f2c.h"
#include "mthrdslv.h"

#define BIGINT   (1 << 30)

#define MINPIV   1.0e-8

int
naivfct(RMAT a, integer neq, integer *nzr, IMAT ri,
		integer *nzc, IMAT ci, NZMAT nz, 
		integer *piv, doublereal minpiv)
{
	integer todo[neq];
	integer i, j, k, m, pvr, pvc, nr, nc, r;
	integer *pri, *pci;
	char *prik;
	doublereal den, mul, mulpiv, fapvr,fari;
	doublereal *par, *papvr;

	if (!minpiv) {
		minpiv = MINPIV;
	}
	for (pvr = 0; pvr < neq; pvr++) {
		piv[pvr] = todo[pvr] = -1;
	}
	for (i = 0; i < neq; i++) {
		if (!nzr[i]) { return ENULCOL + i; }
		nc = BIGINT;	
		nr = nzr[i];
		mul = 0.0;
		pri = ri[i];
		pvr = pri[0];
		fapvr = 0.;
		mulpiv = 0.;
		for (k = 0; k < nr; k++) {
			r = pri[k];
			if (todo[r]) {
				fari = fabs(a[r][i]);
				if (fari > mul) {
					mul = fari;
					mulpiv = mul*minpiv;
					if (nzc[r] <= nc  || mul > fapvr) {
						nc = nzc[pvr = r];
						fapvr = mul/minpiv;
					}
				} else if (nzc[r] < nc && fari > mulpiv) {
					nc = nzc[pvr = r];
				}
			}
		}
		if (nc == BIGINT) { return ENOPIV + i; }
		//printf("i = %d, pvr = %d\n",i,pvr);
/*// #if PIVMETH == SPRSPIV
// 		mul *= minpiv;
// 		for (k = 0; k < nr; k++) {
// 			r = pri[k];
// 			if (todo[r] && nzc[r] < m && fabs(a[r][i]) > mul) {
// 				m = nzc[pvr = r];	
// 			}
// 		}
// #endif *//* PIVMETH == SPRSPIV */
		piv[i] = pvr;
		todo[pvr] = 0;
		papvr = a[pvr];
		den = papvr[i] = 1.0/papvr[i];
		//nr = nzr[i];
		//nc = nzc[pvr];
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

void
naivslv(RMAT a, integer neq, integer *nzc, IMAT ci, 
		doublereal *rhs, doublereal * sol,
		integer *piv)
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

	return;
}

