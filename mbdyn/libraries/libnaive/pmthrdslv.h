/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2004-2013
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

#ifndef pmthrdslv_h
#define pmthrdslv_h

#ifdef USE_NAIVE_MULTITHREAD

#include <atomic_ops.h>
#include "mthrdslv.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

extern int pnaivfct(RMAT a, integer neq, integer *nzr, IMAT ri,
	integer *nzc, IMAT ci, integer* nril, IMAT ril, NZMAT nz,
	integer *piv, integer *todo, doublereal minpiv,
	AO_t *row_locks, volatile AO_TS_t *col_locks, int task, int ncpu);

extern int pnaivslv(RMAT a, integer neq, integer *nzc, IMAT ci, 
	doublereal *rhs, integer *piv, doublereal *fwd, doublereal *sol,
	unsigned long *locks, int task, int ncpu);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* USE_NAIVE_MULTITHREAD */

#endif /* mthrdslv_h */

