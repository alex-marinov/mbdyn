/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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

#ifndef ac_spinlock_h
#define ac_spinlock_h

/*
 * inspired by <asm/system.h>; requires <sys/types.h>
 */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifdef HAVE_CMPXCHG

#if 1

struct __xchg_dummy { unsigned long a[100]; };
#define __xg(x) ((struct __xchg_dummy *)(x))

/* from <asm/bitops.h> */
#define LOCK_PREFIX "lock ; "

#else

#define __xg(x) (x)
#define LOCK_PREFIX ""

#endif

static inline int8_t
mbdyn_cmpxchgb(int8_t *valptr, int8_t newval, int8_t oldval)
{
	int8_t	prev;

	__asm__ __volatile__( LOCK_PREFIX "cmpxchgb %b1,%2"
			: "=a"(prev)
			: "q"(newval), "m"(*__xg(valptr)), "0"(oldval)
			: "memory");

	return prev;
}

static inline int16_t
mbdyn_cmpxchgw(int16_t *valptr, int16_t newval, int16_t oldval)
{
	int16_t	prev;

	__asm__ __volatile__(LOCK_PREFIX "cmpxchgw %w1,%2"
    			: "=a"(prev)
   			: "q"(newval), "m"(*__xg(valptr)), "0"(oldval)
			: "memory");

	return prev;
}

static inline int32_t
mbdyn_cmpxchgl(int32_t *valptr, int32_t newval, int32_t oldval)
{
	int32_t	prev;

	__asm__ __volatile__(LOCK_PREFIX "cmpxchgl %1,%2"
			: "=a"(prev)
			: "q"(newval), "m"(*__xg(valptr)), "0"(oldval)
			: "memory");

	return prev;
}

#endif /* HAVE_CMPXCHG */

#ifdef __cplusplus
}

static inline bool
mbdyn_compare_and_swap(int8_t *valptr, int8_t newval, int8_t oldval)
{
#if defined(HAVE_CMPXCHG)
	return (mbdyn_cmpxchgb(valptr, newval, oldval) == oldval);
#elif defined(HAVE_COMPARE_AND_SWAP)
	atomic_p	word_addr = (atomic_p *)valptr;
	int		old_val = oldval;
	int		new_val = newval;

	return compare_and_swap(word_addr, &old_val, new_val);
#else
	/* FIXME: provide an alternative ... */
#endif
}

static inline bool
mbdyn_compare_and_swap(int16_t *valptr, int16_t newval, int16_t oldval)
{
#if defined(HAVE_CMPXCHG)
	return (mbdyn_cmpxchgw(valptr, newval, oldval) == oldval);
#elif defined(HAVE_COMPARE_AND_SWAP)
	atomic_p	word_addr = (atomic_p *)valptr;
	int		old_val = oldval;
	int		new_val = newval;

	return compare_and_swap(word_addr, &old_val, new_val);
#else
	/* FIXME: provide an alternative ... */
#endif
}

static inline bool
mbdyn_compare_and_swap(int32_t *valptr, int32_t newval, int32_t oldval)
{
#if defined(HAVE_CMPXCHG)
	return (mbdyn_cmpxchgl(valptr, newval, oldval) == oldval);
#elif defined(HAVE_COMPARE_AND_SWAP)
	atomic_p	word_addr = (atomic_p *)valptr;
	int		old_val = oldval;
	int		new_val = newval;

	return compare_and_swap(word_addr, &old_val, new_val);
#else
	/* FIXME: provide an alternative ... */
#endif
}

#endif /* __cplusplus */
	
#endif /* ac_spinlock_h */
