/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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
 * #ifdef USE_RTAI, must be included before the mbdyn_pthread_* 
 * funcs are called */

#ifdef USE_RTAI
/* from <pthread.h> */
#define mbdyn_pthread_create			pthread_create_rt
#define mbdyn_pthread_self			pthread_self_rt
#define mbdyn_pthread_equal			pthread_equal_rt
#define mbdyn_pthread_exit			pthread_exit_rt
#define mbdyn_pthread_join			pthread_join_rt
#define mbdyn_pthread_detach			pthread_detach_rt
#define mbdyn_pthread_attr_init			pthread_attr_init_rt
#define mbdyn_pthread_attr_destroy		pthread_attr_destroy_rt
#define mbdyn_pthread_attr_setdetachstate	pthread_attr_setdetachstate_rt
#define mbdyn_pthread_attr_getdetachstate	pthread_attr_getdetachstate_rt
#define mbdyn_pthread_attr_setschedparam	pthread_attr_setschedparam_rt
#define mbdyn_pthread_attr_getschedparam	pthread_attr_getschedparam_rt
#define mbdyn_pthread_attr_setschedpolicy	pthread_attr_setschedpolicy_rt
#define mbdyn_pthread_attr_getschedpolicy	pthread_attr_getschedpolicy_rt
#define mbdyn_pthread_attr_setinheritsched	pthread_attr_setinheritsched_rt
#define mbdyn_pthread_attr_getinheritsched	pthread_attr_getinheritsched_rt
#define mbdyn_pthread_attr_setscope		pthread_attr_setscope_rt
#define mbdyn_pthread_attr_getscope		pthread_attr_getscope_rt
#define mbdyn_pthread_attr_setguardsize		pthread_attr_setguardsize_rt
#define mbdyn_pthread_attr_getguardsize		pthread_attr_getguardsize_rt
#define mbdyn_pthread_attr_setstackaddr		pthread_attr_setstackaddr_rt
#define mbdyn_pthread_attr_getstackaddr		pthread_attr_getstackaddr_rt
#define mbdyn_pthread_attr_setstack		pthread_attr_setstack_rt
#define mbdyn_pthread_attr_getstack		pthread_attr_getstack_rt
#define mbdyn_pthread_attr_setstacksize		pthread_attr_setstacksize_rt
#define mbdyn_pthread_attr_getstacksize		pthread_attr_getstacksize_rt
#define mbdyn_pthread_setschedparam		pthread_setschedparam_rt
#define mbdyn_pthread_getschedparam		pthread_getschedparam_rt
#define mbdyn_pthread_setconcurrency		pthread_setconcurrency_rt
#define mbdyn_pthread_getconcurrency		pthread_getconcurrency_rt
#define mbdyn_pthread_mutex_init		pthread_mutex_init_rt
#define mbdyn_pthread_mutex_destroy		pthread_mutex_destroy_rt
#define mbdyn_pthread_mutex_trylock		pthread_mutex_trylock_rt
#define mbdyn_pthread_mutex_lock		pthread_mutex_lock_rt
#define mbdyn_pthread_mutex_timedlock		pthread_mutex_timedlock_rt
#define mbdyn_pthread_mutex_unlock		pthread_mutex_unlock_rt
#define mbdyn_pthread_mutexattr_init		pthread_mutexattr_init_rt
#define mbdyn_pthread_mutexattr_destroy		pthread_mutexattr_destroy_rt
#define mbdyn_pthread_mutexattr_setpshared	pthread_mutexattr_setpshared_rt
#define mbdyn_pthread_mutexattr_getpshared	pthread_mutexattr_getpshared_rt
#define mbdyn_pthread_mutexattr_settype		pthread_mutexattr_settype_rt
#define mbdyn_pthread_mutexattr_gettype		pthread_mutexattr_gettype_rt
#define mbdyn_pthread_cond_init			pthread_cond_init_rt
#define mbdyn_pthread_cond_destroy		pthread_cond_destroy_rt
#define mbdyn_pthread_cond_signal		pthread_cond_signal_rt
#define mbdyn_pthread_cond_broadcast		pthread_cond_broadcast_rt
#define mbdyn_pthread_cond_wait			pthread_cond_wait_rt
#define mbdyn_pthread_cond_timedwait		pthread_cond_timedwait_rt
#define mbdyn_pthread_condattr_init		pthread_condattr_init_rt
#define mbdyn_pthread_condattr_destroy		pthread_condattr_destroy_rt
#define mbdyn_pthread_condattr_getpshared	pthread_condattr_getpshared_rt
#define mbdyn_pthread_condattr_setpshared	pthread_condattr_setpshared_rt
#define mbdyn_pthread_rwlock_init		pthread_rwlock_init_rt
#define mbdyn_pthread_rwlock_destroy		pthread_rwlock_destroy_rt
#define mbdyn_pthread_rwlock_rdlock		pthread_rwlock_rdlock_rt
#define mbdyn_pthread_rwlock_tryrdlock		pthread_rwlock_tryrdlock_rt
#define mbdyn_pthread_rwlock_timedrdlock	pthread_rwlock_timedrdlock_rt
#define mbdyn_pthread_rwlock_wrlock		pthread_rwlock_wrlock_rt
#define mbdyn_pthread_rwlock_trywrlock		pthread_rwlock_trywrlock_rt
#define mbdyn_pthread_rwlock_timedwrlock	pthread_rwlock_timedwrlock_rt
#define mbdyn_pthread_rwlock_unlock		pthread_rwlock_unlock_rt
#define mbdyn_pthread_rwlockattr_init		pthread_rwlockattr_init_rt
#define mbdyn_pthread_rwlockattr_destroy	pthread_rwlockattr_destroy_rt
#define mbdyn_pthread_rwlockattr_setpshared	pthread_rwlockattr_setpshared_rt
#define mbdyn_pthread_rwlockattr_getpshared	pthread_rwlockattr_getpshared_rt
#define mbdyn_pthread_spin_init			pthread_spin_init_rt
#define mbdyn_pthread_spin_destroy		pthread_spin_destroy_rt
#define mbdyn_pthread_spin_lock			pthread_spin_lock_rt
#define mbdyn_pthread_spin_trylock		pthread_spin_trylock_rt
#define mbdyn_pthread_spin_unlock		pthread_spin_unlock_rt
/* barriers ?... */
#define mbdyn_pthread_key_create		pthread_key_create_rt
#define mbdyn_pthread_key_delete		pthread_key_delete_rt
#define mbdyn_pthread_setspecific		pthread_setspecific_rt
#define mbdyn_pthread_getspecific		pthread_getspecific_rt
#define mbdyn_pthread_once			pthread_once_rt
#define mbdyn_pthread_setcancelstate		pthread_setcancelstate_rt
#define mbdyn_pthread_setcanceltype		pthread_setcanceltype_rt
#define mbdyn_pthread_cancel			pthread_cancel_rt
#define mbdyn_pthread_testcancel		pthread_testcancel_rt
#define mbdyn_pthread_getcpuclockid		pthread_getcpuclockid_rt
#define mbdyn_pthread_cleanup_push		pthread_cleanup_push_rt
#define mbdyn_pthread_cleanup_pop		pthread_cleanup_pop_rt
#define mbdyn_pthread_atfork			pthread_atfork_rt
/* from <bits/sigthread.h> */
#define mbdyn_pthread_sigmask			pthread_sigmask_rt
#define mbdyn_pthread_kill			pthread_kill_rt
#else /* !USE_RTAI */
/* from <pthread.h> */
#define mbdyn_pthread_create			pthread_create
#define mbdyn_pthread_self			pthread_self
#define mbdyn_pthread_equal			pthread_equal
#define mbdyn_pthread_exit			pthread_exit
#define mbdyn_pthread_join			pthread_join
#define mbdyn_pthread_detach			pthread_detach
#define mbdyn_pthread_attr_init			pthread_attr_init
#define mbdyn_pthread_attr_destroy		pthread_attr_destroy
#define mbdyn_pthread_attr_setdetachstate	pthread_attr_setdetachstate
#define mbdyn_pthread_attr_getdetachstate	pthread_attr_getdetachstate
#define mbdyn_pthread_attr_setschedparam	pthread_attr_setschedparam
#define mbdyn_pthread_attr_getschedparam	pthread_attr_getschedparam
#define mbdyn_pthread_attr_setschedpolicy	pthread_attr_setschedpolicy
#define mbdyn_pthread_attr_getschedpolicy	pthread_attr_getschedpolicy
#define mbdyn_pthread_attr_setinheritsched	pthread_attr_setinheritsched
#define mbdyn_pthread_attr_getinheritsched	pthread_attr_getinheritsched
#define mbdyn_pthread_attr_setscope		pthread_attr_setscope
#define mbdyn_pthread_attr_getscope		pthread_attr_getscope
#define mbdyn_pthread_attr_setguardsize		pthread_attr_setguardsize
#define mbdyn_pthread_attr_getguardsize		pthread_attr_getguardsize
#define mbdyn_pthread_attr_setstackaddr		pthread_attr_setstackaddr
#define mbdyn_pthread_attr_getstackaddr		pthread_attr_getstackaddr
#define mbdyn_pthread_attr_setstack		pthread_attr_setstack
#define mbdyn_pthread_attr_getstack		pthread_attr_getstack
#define mbdyn_pthread_attr_setstacksize		pthread_attr_setstacksize
#define mbdyn_pthread_attr_getstacksize		pthread_attr_getstacksize
#define mbdyn_pthread_setschedparam		pthread_setschedparam
#define mbdyn_pthread_getschedparam		pthread_getschedparam
#define mbdyn_pthread_setconcurrency		pthread_setconcurrency
#define mbdyn_pthread_getconcurrency		pthread_getconcurrency
#define mbdyn_pthread_mutex_init		pthread_mutex_init
#define mbdyn_pthread_mutex_destroy		pthread_mutex_destroy
#define mbdyn_pthread_mutex_trylock		pthread_mutex_trylock
#define mbdyn_pthread_mutex_lock		pthread_mutex_lock
#define mbdyn_pthread_mutex_timedlock		pthread_mutex_timedlock
#define mbdyn_pthread_mutex_unlock		pthread_mutex_unlock
#define mbdyn_pthread_mutexattr_init		pthread_mutexattr_init
#define mbdyn_pthread_mutexattr_destroy		pthread_mutexattr_destroy
#define mbdyn_pthread_mutexattr_setpshared	pthread_mutexattr_setpshared
#define mbdyn_pthread_mutexattr_getpshared	pthread_mutexattr_getpshared
#define mbdyn_pthread_mutexattr_settype		pthread_mutexattr_settype
#define mbdyn_pthread_mutexattr_gettype		pthread_mutexattr_gettype
#define mbdyn_pthread_cond_init			pthread_cond_init
#define mbdyn_pthread_cond_destroy		pthread_cond_destroy
#define mbdyn_pthread_cond_signal		pthread_cond_signal
#define mbdyn_pthread_cond_broadcast		pthread_cond_broadcast
#define mbdyn_pthread_cond_wait			pthread_cond_wait
#define mbdyn_pthread_cond_timedwait		pthread_cond_timedwait
#define mbdyn_pthread_condattr_init		pthread_condattr_init
#define mbdyn_pthread_condattr_destroy		pthread_condattr_destroy
#define mbdyn_pthread_condattr_getpshared	pthread_condattr_getpshared
#define mbdyn_pthread_condattr_setpshared	pthread_condattr_setpshared
#define mbdyn_pthread_rwlock_init		pthread_rwlock_init
#define mbdyn_pthread_rwlock_destroy		pthread_rwlock_destroy
#define mbdyn_pthread_rwlock_rdlock		pthread_rwlock_rdlock
#define mbdyn_pthread_rwlock_tryrdlock		pthread_rwlock_tryrdlock
#define mbdyn_pthread_rwlock_timedrdlock	pthread_rwlock_timedrdlock
#define mbdyn_pthread_rwlock_wrlock		pthread_rwlock_wrlock
#define mbdyn_pthread_rwlock_trywrlock		pthread_rwlock_trywrlock
#define mbdyn_pthread_rwlock_timedwrlock	pthread_rwlock_timedwrlock
#define mbdyn_pthread_rwlock_unlock		pthread_rwlock_unlock
#define mbdyn_pthread_rwlockattr_init		pthread_rwlockattr_init
#define mbdyn_pthread_rwlockattr_destroy	pthread_rwlockattr_destroy
#define mbdyn_pthread_rwlockattr_setpshared	pthread_rwlockattr_setpshared
#define mbdyn_pthread_rwlockattr_getpshared	pthread_rwlockattr_getpshared
#define mbdyn_pthread_spin_init			pthread_spin_init
#define mbdyn_pthread_spin_destroy		pthread_spin_destroy
#define mbdyn_pthread_spin_lock			pthread_spin_lock
#define mbdyn_pthread_spin_trylock		pthread_spin_trylock
#define mbdyn_pthread_spin_unlock		pthread_spin_unlock
/* barriers ?... */
#define mbdyn_pthread_key_create		pthread_key_create
#define mbdyn_pthread_key_delete		pthread_key_delete
#define mbdyn_pthread_setspecific		pthread_setspecific
#define mbdyn_pthread_getspecific		pthread_getspecific
#define mbdyn_pthread_once			pthread_once
#define mbdyn_pthread_setcancelstate		pthread_setcancelstate
#define mbdyn_pthread_setcanceltype		pthread_setcanceltype
#define mbdyn_pthread_cancel			pthread_cancel
#define mbdyn_pthread_testcancel		pthread_testcancel
#define mbdyn_pthread_getcpuclockid		pthread_getcpuclockid
#define mbdyn_pthread_cleanup_push		pthread_cleanup_push
#define mbdyn_pthread_cleanup_pop		pthread_cleanup_pop
#define mbdyn_pthread_atfork			pthread_atfork
/* from <bits/sigthread.h> */
#define mbdyn_pthread_sigmask			pthread_sigmask
#define mbdyn_pthread_kill			pthread_kill
#endif /* !USE_RTAI */

#endif /* AC_PTHREAD_H */

