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

 /*
  * With the contribution of Runsen Zhang <runsen.zhang@polimi.it>
  * during Google Summer of Code 2020
  */

#ifndef MBDYN_CE_H
#define MBDYN_CE_H

extern "C" {

// opaque pointer to element's data
typedef void * MBDyn_CE_t *;

// creates a new instance of the element, returning an opaque pointer to its data
extern MBDyn_CE_t *
MBDyn_CE_init(void);

// destroy
extern void
MBDyn_CE_destroy(MBDyn_CE_t *);

// add if needed
extern void
MBDyn_CE_AfterPredict(MBDyn_CE_t *);

// add arguments as needed
extern void
MBDyn_CE_Exchange(MBDyn_CE_t *, double *x, double *R, double *f, double *m);

// add if needed
extern void
MBDyn_CE_AfterConvergence(MBDyn_CE_t *);

}

#endif // MBDYN_CE_H
