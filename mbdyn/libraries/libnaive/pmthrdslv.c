/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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
 * Copyright Morandini (please complete)
 */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#ifdef USE_NAIVE_MULTITHREAD

#define __KERNEL__
#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <sys/poll.h>

struct task_struct;

#include <asm/system.h>
#include <asm/atomic.h>
#include <asm/bitops.h>

#define  SPRSPIV


/* FIXME: from <f2c.h> */
typedef int     integer;
typedef double  doublereal;


#define BIGINT   (1 << 30)
#define ENULCOL  1000000
#define ENOPIV   2000000

#define MINPIV   1.0e-8


int pnaivfct(doublereal** a,
	integer neq,
	integer *nzr, integer** ri,
	integer *nzc, integer** ci,
	char** nz,
	integer *piv,
	integer *todo,
	doublereal minpiv,
	unsigned long *row_locks,
	unsigned long *col_locks,
	int task,
	int ncpu)
{
	integer i, j, k, m, pvr = 0, pvc, nr, nc, r;
	integer *pri, *pci;
	char *pnzk;
	doublereal den, mul;
	doublereal *par, *papvr;
	unsigned long *sync_lock = &row_locks[neq];
	unsigned long *pivot_lock = &row_locks[neq + 1];

	if (minpiv == 0.) {
		minpiv = MINPIV;
	}

/*	atomic_inc((atomic_t *)sync_lock);
 	while (atomic_read((atomic_t *)sync_lock) < ncpu);
*/	
	for (i = 0; i < neq; i++) {
		nr = nzr[i];
		if (nr == 0) {
			return ENULCOL + i; 
		}
		pri = ri[i];
		if (task == 0) {
			m = BIGINT;	
			mul = 0.0;
			for (k = 0; k < nr; k++) {
				r = pri[k];
				doublereal dAri = fabs(a[r][i]);
				if (todo[r] && dAri > mul) {
					mul = dAri;
					pvr = r;
					m = nzc[pvr];	
				}
			}
			if (m == BIGINT) {
				return ENOPIV + i;
			}
#ifdef SPRSPIV
			mul *= minpiv;
			for (k = 0; k < nr; k++) {
				r = pri[k];
				if (todo[r] && nzc[r] < m && fabs(a[r][i]) > mul) {
					pvr = r;
					m = nzc[pvr];
				}
			}
#endif
			todo[pvr] = 0;
			papvr = a[pvr];
			den = 1.0/papvr[i];
			if (pri[i]%ncpu == task) {
				papvr[i] = den;
				wmb();
				atomic_set((atomic_t *)pivot_lock, 1);
			}

			set_wmb(piv[i], pvr);

			while (atomic_read((atomic_t *)pivot_lock) != 1);

		} else {
			while ((pvr = atomic_read((atomic_t *)&piv[i])) < 0);
			papvr = a[pvr];

			if (pri[i]%ncpu == task) {
				den = papvr[i] = 1./papvr[i];
				wmb();
				atomic_set((atomic_t *)pivot_lock, 1);
			} else {
				while (atomic_read((atomic_t *)pivot_lock) != 1);
				den = papvr[i];
			}
		}

		nc  = nzc[pvr];
		pci = ci[pvr];

		for (k = task; k < nr; k += ncpu) {
			r = pri[k];
			if (todo[r] == 0) {
				continue;
			}
			pnzk = nz[r];
			par = a[r];
			mul = par[i] = par[i]*den;
			for (j = 0; j < nc; j++) {
				pvc = pci[j];
				if (pvc <= i) {
					continue;
				}

				if (pnzk[pvc]) {
					par[pvc] -= mul*papvr[pvc];
				} else {
					par[pvc] = -mul*papvr[pvc];
					ci[r][nzc[r]++] = pvc;
					while (atomic_inc_and_test((atomic_t *)&col_locks[pvc]));
					pnzk[pvc] = 1;
					ri[pvc][nzr[pvc]++] = r;
					atomic_set((atomic_t *)&col_locks[pvc], 0); 
				}
			}
		}

		atomic_inc((atomic_t *)&row_locks[i]);
		while (atomic_read((atomic_t *)&row_locks[i]) < ncpu);
	}

	if (task == 0) {
		/* NOTE: same as sync_lock[0] */
		sync_lock[0] = 0;

		/* NOTE: same as pivot_lock[0] */
		pivot_lock[0] = 0;

		wmb();
	}

	return 0;
}

void pnaivslv(doublereal** a,
		integer neq,
		integer *nzc,
		integer** ci, 
		doublereal *rhs,
		integer *piv, 
		doublereal *fwd,
		doublereal *sol,
		unsigned long *locks,
		int task,
		int ncpu)
{

	integer i, k, nc, r, c;
	integer *pci;
	doublereal s;
	doublereal *par;

	if (!task) {
		fwd[0] = rhs[piv[0]];
		wmb();
		atomic_set((atomic_t *)locks, 1);
	}

	for (i = 1; i < neq; i++) {
	//for (i = 1 + task; i < neq; i+=ncpu) {
		if (i%ncpu != task) { continue; }
		//used[task + 2]++;
		nc = nzc[r = piv[i]];
		s = rhs[r];
		par = a[r];
		pci = ci[r];
		for (k = 0; k < nc; k++) {
			c = pci[k];
			if (c < i) {
				while (!atomic_read((atomic_t *)&locks[c]));
				s -= par[c]*fwd[c];
			}
		}
		set_wmb(fwd[i], s);
		//fwd[i] = s;
		atomic_set((atomic_t *)&locks[i], 1);
	}

	atomic_inc((atomic_t *)&locks[neq]);
	while (atomic_read((atomic_t *)&locks[neq]) < ncpu);

	neq--;
	r = piv[neq];
	if (!task) {
		sol[neq] = fwd[neq]*a[r][neq];
		atomic_set((atomic_t *)&locks[neq], 0);
	}
	for (i = neq - 1; i >= 0; i--) {
	//for (i = neq - 1 - task; i >= 0; i-=ncpu) {
		if (i%ncpu != task) { continue; }
		//used[task + 4]++;
		r = piv[i];
		nc = nzc[r];
		s = fwd[i];
		par = a[r];
		pci = ci[r];
		for (k = 0; k < nc; k++) {
			if ((c = pci[k]) > i) {
				while (atomic_read((atomic_t *)&locks[c]));
				s -= par[c]*sol[c];
			}
		}
		set_wmb(sol[i], s*par[i]);
		//sol[i] = s*par[i];
		atomic_set((atomic_t *)&locks[i], 0);
	}

	return;
}

void naivesad(doublereal** a,
	integer** ri,
	integer *nzr,
	integer** ci,
	integer *nzc,
	char ** nz,
	integer neq,
	integer *rowindx,
	integer *colindx,
	doublereal *elmat,
	integer elr,
	integer elc,
	integer eldc)
{
	integer i, k, r, c, er;
	doublereal *par;

	er = 0;
	for (i = 0; i < elr; i++) {
		if ((r = rowindx[i]) >= 0) {
			par = a[r];
			for (k = 0; k < elc; k++) {
				c = colindx[k];
				if (!elmat[er + k] && c >= 0) {
					if (nz[r][c]) {
						par[c] += elmat[er + k];
					} else {
						par[c] = elmat[er + k];
						nz[r][c] = 1;
						ri[c][nzr[c]++] = r;
						ci[r][nzc[r]++] = c;
					}
				}
			}
		}
		er += eldc;
	}
}

void naivepsad(doublereal** ga,
	integer** gri,
	integer *gnzr,
	integer** gci,
	integer *gnzc,
	char ** nz,
	doublereal** a,
	integer** ci,
	integer *nzc,
	integer from,
	integer to)
{
	integer i, r, c, nc;
	doublereal *pgar, *par;

	for (r = from; r <= to; r++) {
		if ((nc = nzc[r])) {
			pgar = ga[r];
			par  = a[r];
			for (i = 0; i < nc; i++) {
				if (nz[r][c = ci[r][i]]) {
					pgar[c] += par[c];
				} else {
					pgar[c] = par[c];
					nz[r][c] = 1;
					gri[c][gnzr[c]++] = r;
					gci[r][gnzc[r]++] = c;
				}
			}
		}
	}
}

#endif /* USE_NAIVE_MULTITHREAD */
