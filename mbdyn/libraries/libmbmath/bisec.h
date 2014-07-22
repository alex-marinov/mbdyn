/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

#ifndef BISEC_H
#define BISEC_H

#ifdef __cplusplus

/*
 * algoritmo di bisezione per ricerca efficiente degli angoli
 *
 * v:	array dei valori
 * val:	valore da cercare
 * lb:	indice inferiore
 * ub:	indice superiore
 *
 * restituisce l'indice corrispondente al valore di v
 * che approssima per difetto val;
 * se val > v[ub] restituisce ub;
 * se val < v[lb] restituisce lb - 1.
 * quindi v[ret] <= val < v[ret + 1]
 *
 * Nota: ovviamente si presuppone che i valori di v siano ordinati in modo
 * crescente e strettamente monotono, ovvero v[i] < v[i + 1].
 */
template <class T>
int 
bisec(const T *const v, const T& val, int lb, int ub)
{
	if (val < v[lb]) {
		return lb - 1;
	}
	
	if (val > v[ub]) {
		return ub;
	}
	
	while (ub > lb + 1) {
		int b = (lb + ub)/2;
		if (v[b] > val) {
			ub = b;
		} else if (v[b] < val) {
			lb = b;
		} else {
			return b;
		}
	}
   
   	return lb;
}

#endif /* __cplusplus */

#ifdef __cplusplus
extern "C"
#endif
int bisec_d(doublereal *v, doublereal val, int lb, int ub);

#endif /* BISEC_H */
