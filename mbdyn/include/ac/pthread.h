/* $Header$ */
/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2007
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

#ifndef AC_PTHREAD_H
#define AC_PTHREAD_H

/* either libc's <pthread.h>, or RTAI's <rtai_usp_posix.h> 
 * #ifdef USE_RTAI, must be included before the pthread_* 
 * funcs are called */

#ifdef USE_RTAI

#ifdef HAVE_RTAI_USP_POSIX_H
#include <rtai_usp_posix.h>
#endif /* HAVE_RTAI_USP_POSIX_H */

/* from <pthread.h> */
#define pthread_create			pthread_create_rt
#define pthread_self			pthread_self_rt
#define pthread_equal			pthread_equal_rt
#define pthread_exit			pthread_exit_rt
#define pthread_join			pthread_join_rt
#define pthread_detach			pthread_detach_rt
#define pthread_attr_init		pthread_attr_init_rt
#define pthread_attr_destroy		pthread_attr_destroy_rt
#define pthread_attr_setdetachstate	pthread_attr_setdetachstate_rt
#define pthread_attr_getdetachstate	pthread_attr_getdetachstate_rt
#define pthread_attr_setschedparam	pthread_attr_setschedparam_rt
#define pthread_attr_getschedparam	pthread_attr_getschedparam_rt
#define pthread_attr_setschedpolicy	pthread_attr_setschedpolicy_rt
#define pthread_attr_getschedpolicy	pthread_attr_getschedpolicy_rt
#define pthread_attr_setinheritsched	pthread_attr_setinheritsched_rt
#define pthread_attr_getinheritsched	pthread_attr_getinheritsched_rt
#define pthread_attr_setscope		pthread_attr_setscope_rt
#define pthread_attr_getscope		pthread_attr_getscope_rt
#define pthread_attr_setguardsize	pthread_attr_setguardsize_rt
#define pthread_attr_getguardsize	pthread_attr_getguardsize_rt
#define pthread_attr_setstackaddr	pthread_attr_setstackaddr_rt
#define pthread_attr_getstackaddr	pthread_attr_getstackaddr_rt
#define pthread_attr_setstack		pthread_attr_setstack_rt
#define pthread_attr_getstack		pthread_attr_getstack_rt
#define pthread_attr_setstacksize	pthread_attr_setstacksize_rt
#define pthread_attr_getstacksize	pthread_attr_getstacksize_rt
#define pthread_setschedparam		pthread_setschedparam_rt
#define pthread_getschedparam		pthread_getschedparam_rt
#define pthread_setconcurrency		pthread_setconcurrency_rt
#define pthread_getconcurrency		pthread_getconcurrency_rt
#define pthread_mutex_init		pthread_mutex_init_rt
#define pthread_mutex_destroy		pthread_mutex_destroy_rt
#define pthread_mutex_trylock		pthread_mutex_trylock_rt
#define pthread_mutex_lock		pthread_mutex_lock_rt
#define pthread_mutex_timedlock		pthread_mutex_timedlock_rt
#define pthread_mutex_unlock		pthread_mutex_unlock_rt
#define pthread_mutexattr_init		pthread_mutexattr_init_rt
#define pthread_mutexattr_destroy	pthread_mutexattr_destroy_rt
#define pthread_mutexattr_setpshared	pthread_mutexattr_setpshared_rt
#define pthread_mutexattr_getpshared	pthread_mutexattr_getpshared_rt
#define pthread_mutexattr_settype	pthread_mutexattr_settype_rt
#define pthread_mutexattr_gettype	pthread_mutexattr_gettype_rt
#define pthread_cond_init		pthread_cond_init_rt
#define pthread_cond_destroy		pthread_cond_destroy_rt
#define pthread_cond_signal		pthread_cond_signal_rt
#define pthread_cond_broadcast		pthread_cond_broadcast_rt
#define pthread_cond_wait		pthread_cond_wait_rt
#define pthread_cond_timedwait		pthread_cond_timedwait_rt
#define pthread_condattr_init		pthread_condattr_init_rt
#define pthread_condattr_destroy	pthread_condattr_destroy_rt
#define pthread_condattr_getpshared	pthread_condattr_getpshared_rt
#define pthread_condattr_setpshared	pthread_condattr_setpshared_rt
#define pthread_rwlock_init		pthread_rwlock_init_rt
#define pthread_rwlock_destroy		pthread_rwlock_destroy_rt
#define pthread_rwlock_rdlock		pthread_rwlock_rdlock_rt
#define pthread_rwlock_tryrdlock	pthread_rwlock_tryrdlock_rt
#define pthread_rwlock_timedrdlock	pthread_rwlock_timedrdlock_rt
#define pthread_rwlock_wrlock		pthread_rwlock_wrlock_rt
#define pthread_rwlock_trywrlock	pthread_rwlock_trywrlock_rt
#define pthread_rwlock_timedwrlock	pthread_rwlock_timedwrlock_rt
#define pthread_rwlock_unlock		pthread_rwlock_unlock_rt
#define pthread_rwlockattr_init		pthread_rwlockattr_init_rt
#define pthread_rwlockattr_destroy	pthread_rwlockattr_destroy_rt
#define pthread_rwlockattr_setpshared	pthread_rwlockattr_setpshared_rt
#define pthread_rwlockattr_getpshared	pthread_rwlockattr_getpshared_rt
#define pthread_spin_init		pthread_spin_init_rt
#define pthread_spin_destroy		pthread_spin_destroy_rt
#define pthread_spin_lock		pthread_spin_lock_rt
#define pthread_spin_trylock		pthread_spin_trylock_rt
#define pthread_spin_unlock		pthread_spin_unlock_rt
/* barriers ?... */
#define pthread_key_create		pthread_key_create_rt
#define pthread_key_delete		pthread_key_delete_rt
#define pthread_setspecific		pthread_setspecific_rt
#define pthread_getspecific		pthread_getspecific_rt
#define pthread_once			pthread_once_rt
#define pthread_setcancelstate		pthread_setcancelstate_rt
#define pthread_setcanceltype		pthread_setcanceltype_rt
#define pthread_cancel			pthread_cancel_rt
#define pthread_testcancel		pthread_testcancel_rt
#define pthread_getcpuclockid		pthread_getcpuclockid_rt
#define pthread_cleanup_push		pthread_cleanup_push_rt
#define pthread_cleanup_pop		pthread_cleanup_pop_rt
#define pthread_atfork			pthread_atfork_rt
/* from <bits/sigthread.h> */
#define pthread_sigmask			pthread_sigmask_rt
#define pthread_kill			pthread_kill_rt
/* from <semaphore.h> */
#define sem_init			sem_init_rt
#define sem_destroy			sem_destroy_rt
#define sem_open			sem_open_rt
#define sem_close			sem_close_rt
#define sem_unlink			sem_unlink_rt
#define sem_wait			sem_wait_rt
#define	sem_timedwait			sem_timedwait_rt
#define	sem_trywait			sem_trywait_rt
#define	sem_post			sem_post_rt
#define	sem_getvalue			sem_getvalue_rt

#else /* !USE_RTAI */

#ifdef HAVE_PTHREAD_H
#include <pthread.h>
#endif /* HAVE_PTHREAD_H */
#ifdef HAVE_SEMAPHORE_H
#include <semaphore.h>
#endif /* HAVE_SEMAPHORE_H */

#endif /* !USE_RTAI */

#endif /* AC_PTHREAD_H */

