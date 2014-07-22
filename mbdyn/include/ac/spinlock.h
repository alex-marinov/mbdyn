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
/*
 * Copyright (C) 2003-2014
 *
 * Marco Morandini
 *
 */
#ifndef ac_spinlock_h
#define ac_spinlock_h

#include <atomic_ops.h>

#if 0
static inline int
mbdyn_compare_and_swap(
	AO_t *val, 
	AO_t newval, 
	AO_t oldval)
{
	return AO_compare_and_swap_full(val, oldval, newval);
};
#endif

static inline AO_TS_VAL_t
mbdyn_test_and_set(
	volatile AO_TS_t *val)
{
	return AO_test_and_set_full(val);
};

	
#endif /* ac_spinlock_h */
