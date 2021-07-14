/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

#include "ldl.h"

/*
C**********************     UPDLDL_ADD       ***********************

      subroutine uldlad(ldl, nrdldl, n, x, z, nrdz, nz, y)

c     Esegue una modifica di rango unitario del problema dei minimi quadrati 
c     in ricorsione all'aggiunta di un vettore, ovvero aggiorna la
c     fattorizzazione LDL della matrice normale B = B + xx' e la matrice
c     dei termini noti Z = Z + xy'. 
c
c     Parametri:    ldl = matrice contenente L e D
c                nrdldl = numero di righe del dimension di LDL
c                     n = ordine di LDL
c                     x = vettore dei coefficienti da 'aggiungere'
c                     z = matrice dei termini noti
c                  nrdz = numero di righe del dimension di Z
c                    nz = numero di termini noti
c                     y = vettore dei termini noti da 'aggiungere'

      implicit none
      integer nrdldl, n, nrdz, nz
      double precision ldl(nrdldl,n), x(n), z(nrdz,nz), y(nz)
      
      double precision t_old, t_new, temp, beta
      integer i, j
*/

int
uldlad_(doublereal *ldl, integer *pnrdldl, integer *pn, doublereal *x, doublereal *z, integer *pnrdz, integer *pnz, doublereal *y)
{
	integer nrdldl = *pnrdldl, n = *pn, nrdz = *pnrdz, nz = *pnz;

	integer i, j;
	doublereal t_old;

	for (i = 0; i < n; i++) {
		for (j = 0; j < nz; j++) {
			z[i + nrdz*j] += x[i]*y[j];
		}
	}

	t_old = 1.;
	for (i = 0; i < n - 1; i++) {
		doublereal temp, t_new, beta;
		doublereal *pldl;

		temp = x[i];
		pldl = &ldl[i + nrdldl*i];
		t_new = t_old + temp*temp*pldl[0];
		beta = pldl[0]*temp/t_new;
		pldl[0] *= t_old/t_new;
		t_old = t_new;
		for (j = i + 1; j < n; j++) {
			pldl++;
			x[j] -= temp*pldl[0];
			pldl[0] += beta*x[j];
		}
	}

	ldl[nrdldl*n - 1] *= t_old/(t_old + x[n-1]*x[n-1]*ldl[nrdldl*n - 1]);

	return 0;
}

/*
C**********************       LDL_SOLVE       ***********************
      subroutine ldlsol(ldl, nrdldl, b, nrdb, n, nvet)

c     Solutore di un sistema fattorizzato LDL', LDL'x = b
c
c     Parametri:    ldl = matqrice contenente L e D
c                nrdldl = numero di righe del dimension di LDL
c                     b = matrice dei termini noti
c                  nrdb = numero di righe del dimension di B
c                     n = ordine di LDL, nonche' dei termini noti
c                  nvet = numero di termini noti (numero di colonne di B)
c
c     Restituisce: la soluzione sovrascritta in B. 
      
      implicit none
      integer nrdldl, nrdb, n, nvet
      double precision ldl(nrdldl,n), b(nrdb, nvet)
      
      double precision s
      integer i, j, k
*/

int
ldlsol_(doublereal *ldl, integer *pnrdldl, doublereal *b, integer *pnrdb, integer *pn, integer *pnvet)
{
	integer nrdldl = *pnrdldl, nrdb = *pnrdb, n = *pn, nvet = *pnvet;

	integer k;

/*
c     risolve il sistema Lv = b con v = DL'x
*/

	for (k = 0; k < nvet; k++) {
		integer i;

		for (i = 1; i < n; i++) {
			integer j;
			doublereal s;

			s = b[i + k*nrdb];
			for (j = 0; j < i - 1; j++) {
				s -= ldl[i + nrdldl*j]*b[j + nrdb*k];
			}
			b[i + nrdb*k] = s;
		}
	}

/*
c      risolve il sistema L'x = inv(D)v
*/

	for (k = 0; k < nvet; k++) {
		integer i;

		b[n - 1 + nrdb*k] *= ldl[nrdldl*n - 1];
		for (i = n - 2; i >= 0; i--) {
			integer j;
			doublereal s;

			s = b[i + nrdb*k]*ldl[i + nrdldl*i];
			for (j = i + 1; j < n; j++) {
				s -= ldl[j + nrdldl*i]*b[j + nrdb*k];
			}
			b[i + nrdb*k] = s;
		}
	}

	return 0;
}

