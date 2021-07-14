/* $Header$ */
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

#ifndef LDL_H
#define LDL_H

#include "ac/f2c.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*
 **********************     UPDLDL_ADD       ***********************
 * 
 * subroutine updldl_add(ldl, nrdldl, n, x, z, nrdz, nz, y)
 * 
 *     Esegue una modifica di rango unitario del problema dei minimi quadrati 
 *     in ricorsione all'aggiunta di un vettore, ovvero aggiorna la
 *     fattorizzazione LDL della matrice normale B = B + xx' e la matrice
 *     dei termini noti Z = Z + xy'. 
 *
 *     Parametri:    ldl = matrice contenente L e D
 *                nrdldl = numero di righe del dimension di LDL
 *                     n = ordine di LDL
 *                     x = vettore dei coefficienti da 'aggiungere'
 *                     z = matrice dei termini noti
 *                  nrdz = numero di righe del dimension di Z
 *                    nz = numero di termini noti
 *                     y = vettore dei termini noti da 'aggiungere'
 */

extern int __FC_DECL__(uldlad) (doublereal* ldl, 
				    integer* nrdldl,
				    integer* n, 
				    doublereal* x, 
				    doublereal* z, 
				    integer* nrdz, 
				    integer* nz, 
				    doublereal* y);

/*
 **********************       LDL_SOLVE       ***********************
 * 
 * subroutine ldl_solve(ldl, nrdldl, b, nrdb, n, nvet)
 * 
 *     Solutore di un sistema fattorizzato LDL', LDL'x = b
 *
 *     Parametri:    ldl = matrice contenente L e D
 *                nrdldl = numero di righe del dimension di LDL
 *                     b = matrice dei termini noti
 *                  nrdb = numero di righe del dimension di B
 *                     n = ordine di LDL, nonche' dei termini noti
 *                  nvet = numero di termini noti (numero di colonne di B)
 *
 *     Restituisce: la soluzione sovrascritta in B. 
 */

extern int __FC_DECL__(ldlsol) (doublereal* ldl,
				   integer* nrdldl,
				   doublereal* b,
				   integer* nrdb,
				   integer* n, 
				   integer* nvet);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* LDL_H */
