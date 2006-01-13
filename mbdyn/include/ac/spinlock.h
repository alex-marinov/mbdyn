/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2006
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

#if defined(HAVE_i486_CMPXCHG)

#if 1

struct mbdyn__xchg_dummy { unsigned long a[100]; };
#define mbdyn__xg(x) ((struct mbdyn__xchg_dummy *)(x))

/* from <asm/bitops.h> */
#define LOCK_PREFIX "lock ; "

#else

#define mbdyn__xg(x) (x)
#define LOCK_PREFIX ""

#endif

#define CMPXCHG(f) \
	__asm__ __volatile__( LOCK_PREFIX f \
			: "=a"(prev) \
			: "q"(newval), "m"(*mbdyn__xg(valptr)), "0"(oldval) \
			: "memory")
#define CMPXCHGB CMPXCHG("cmpxchgb %b1,%2")
#define CMPXCHGW CMPXCHG("cmpxchgw %w1,%2")
#define CMPXCHGL CMPXCHG("cmpxchgl %1,%2")

static inline int8_t
mbdyn_cmpxchgb(int8_t *valptr, int8_t newval, int8_t oldval)
{
	int8_t	prev;

	CMPXCHGB;

	return prev;
}

static inline int16_t
mbdyn_cmpxchgw(int16_t *valptr, int16_t newval, int16_t oldval)
{
	int16_t	prev;

	CMPXCHGW;

	return prev;
}

static inline int32_t
mbdyn_cmpxchgl(int32_t *valptr, int32_t newval, int32_t oldval)
{
	int32_t	prev;

	CMPXCHGL;

	return prev;
}

#ifdef NEED_INT_VERSION
static inline int
mbdyn_cmpxchg(int *valptr, int newval, int oldval)
{
	int prev;

#if __INT_MAX__ == 32767
	CMPXCHGW;
#else 
	CMPXCHGL;
#endif

	return prev;
}
#endif /* NEED_INT_VERSION */

#elif defined(HAVE_IA64_CMPXCHG) /* HAVE_IA64_CMPXCHG */

#ifdef HAVE_ASM_SYSTEM_H
#include <asm/system.h>
#endif /* HAVE_ASM_SYSTEM_H */

#define	mbdyn_cmpxchg(ptr, newval, oldval) cmpxchg((ptr), (oldval), (newval))

#endif /* */

#ifdef __cplusplus
}

static inline bool
mbdyn_compare_and_swap(int8_t &val, int8_t newval, int8_t oldval)
{
#if defined(HAVE_i486_CMPXCHG)
	return (mbdyn_cmpxchgb(&val, newval, oldval) == oldval);
#elif defined(HAVE_IA64_CMPXCHG)
	return (mbdyn_cmpxchg(&val, newval, oldval) == oldval);
#elif defined(HAVE_COMPARE_AND_SWAP)
	atomic_p	word_addr = (atomic_p *)&val;
	int		old_val = oldval;
	int		new_val = newval;

	return compare_and_swap(word_addr, &old_val, new_val);
#else
	/* FIXME: provide an alternative ... */
#endif
}

static inline bool
mbdyn_compare_and_swap(int16_t &val, int16_t newval, int16_t oldval)
{
#if defined(HAVE_i486_CMPXCHG)
	return (mbdyn_cmpxchgw(&val, newval, oldval) == oldval);
#elif defined(HAVE_IA64_CMPXCHG)
	return (mbdyn_cmpxchg(&val, newval, oldval) == oldval);
#elif defined(HAVE_COMPARE_AND_SWAP)
	atomic_p	word_addr = (atomic_p *)&val;
	int		old_val = oldval;
	int		new_val = newval;

	return compare_and_swap(word_addr, &old_val, new_val);
#else
	/* FIXME: provide an alternative ... */
#endif
}

static inline bool
mbdyn_compare_and_swap(int32_t &val, int32_t newval, int32_t oldval)
{
#if defined(HAVE_i486_CMPXCHG)
	return (mbdyn_cmpxchgl(&val, newval, oldval) == oldval);
#elif defined(HAVE_IA64_CMPXCHG)
	return (mbdyn_cmpxchg(&val, newval, oldval) == oldval);
#elif defined(HAVE_COMPARE_AND_SWAP)
	atomic_p	word_addr = (atomic_p *)&val;
	int		old_val = oldval;
	int		new_val = newval;

	return compare_and_swap(word_addr, &old_val, new_val);
#else
	/* FIXME: provide an alternative ... */
#endif
}

#ifdef NEED_INT_VERSION
static inline bool
mbdyn_compare_and_swap(int &val, int newval, int oldval)
{
#if defined(HAVE_i486_CMPXCHG)
	return (mbdyn_cmpxchg(&val, newval, oldval) == oldval);
#elif defined(HAVE_IA64_CMPXCHG)
	return (mbdyn_cmpxchg(&val, newval, oldval) == oldval);
#elif defined(HAVE_COMPARE_AND_SWAP)
	atomic_p	word_addr = (atomic_p *)&val;
	int		old_val = oldval;
	int		new_val = newval;

	return compare_and_swap(word_addr, &old_val, new_val);
#else
	/* FIXME: provide an alternative ... */
#endif
}
#endif /* NEED_INT_VERSION */

#endif /* __cplusplus */
	
#endif /* ac_spinlock_h */
