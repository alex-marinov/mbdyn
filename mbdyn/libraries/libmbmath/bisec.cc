/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

#include <cmath>

#include "ac/f2c.h"
#include "bisec.h"

/*
 * algoritmo di bisezione per ricerca efficiente degli angoli
 *
 * v:	array dei valori
 * val:	valore da cercare
 * lb:	indice inferiore
 * ub:	indice superiore
 *
 * restituisce l'indice corrispondente al valore di v
 * che approssima per eccesso val;
 * se val > v[ub] restituisce ub + 1;
 * se val < v[lb] restituisce lb.
 * quindi v[ret - 1] <= val <= v[ret]
 *
 * Nota: ovviamente si presuppone che i valori di v siano ordinati in modo
 * crescente e strettamente monotono, ovvero v[i] < v[i + 1].
 */
extern "C" int 
bisec_d(doublereal* v, doublereal val, int lb, int ub)
{
	return bisec<doublereal>(v, val, lb, ub);
}

