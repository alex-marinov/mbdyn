/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2008
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

/* socket driver */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#ifdef USE_RTAI

/*
 * Mante prefers KEEP_STATIC_INLINE instead of linking liblxrt.a
 */
#define KEEP_STATIC_INLINE

#if defined(HAVE_RTAI_LXRT_H)		/* RTAI 3.0 */
#include <rtai_lxrt.h>
#elif defined(HAVE_RTAI_LXRT_USER_H)	/* up to RTAI 24.1.13 */
#include <rtai_lxrt_user.h>
#endif /* HAVE_RTAI_LXRT_USER_H */
#if defined(HAVE_RTAI_NETRPC_H)		/* RTAI 3.0 */
#include <rtai_netrpc.h>
#elif defined(HAVE_NET_RPC_H)		/* up to RTAI 24.1.13 */
#include <net_rpc.h>
#endif /* HAVE_NET_RPC_H */

#include <assert.h>
#include "mbrtai_utils.h"

int
mbdyn_rt_task_init(const char *name, int priority, int stack_size,
		int max_msg_size, int cpu, void **__task)
{
	assert(name != NULL);
	assert(strlen(name) == 6);
	assert(__task != NULL);
	assert(*__task == NULL);

	*__task = (void *)rt_task_init_schmod(nam2num(name), priority,
			stack_size, max_msg_size, SCHED_FIFO,/*0x2*/ cpu);
	return (*__task == NULL);
}

int
mbdyn_rt_task_delete(void **__task)
{
	RT_TASK		*task;
	int		rc;

	assert(__task != NULL);
	assert(*__task != NULL);

	task = (RT_TASK *)*__task;
	rc = rt_task_delete(task);

	*__task = NULL;

	return /* rc */ 0;
}
void
mbdyn_rt_make_hard_real_time(void)
{
	rt_make_hard_real_time();
}
void
mbdyn_rt_make_soft_real_time(void)
{
	rt_make_soft_real_time();
}
int
mbdyn_rt_task_make_periodic(void *__task, long long __start_time,
		long long __period)
{
	
	
	RT_TASK		*task = (RT_TASK *)__task;
	RTIME		start_time = (RTIME)__start_time;
	RTIME		period = (RTIME)__period;

	assert(__task != NULL);

	return rt_task_make_periodic(task, start_time, period);
}

void
mbdyn_rt_task_wait_period(void)
{
	rt_task_wait_period();
}

void
mbdyn_rt_allow_nonroot_hrt(void)
{
	rt_allow_nonroot_hrt();
}

int
mbdyn_rt_request_port(unsigned long node)
{
	assert(node > 0);

	return rt_request_port(node);
}

void
mbdyn_rt_set_oneshot_mode(void)
{
	rt_set_oneshot_mode();
}

void
mbdyn_rt_set_periodic_mode(void)
{
	rt_set_periodic_mode();
}

int
mbdyn_rt_is_hard_timer_running(void)
{
	return rt_is_hard_timer_running();
}

long long
mbdyn_start_rt_timer(long long __period)
{
	RTIME period = (RTIME)__period;
	return start_rt_timer(period);
}

void
mbdyn_stop_rt_timer(void)
{
	stop_rt_timer();
}

long long
mbdyn_rt_get_time(void)
{	
	return (long long)rt_get_time();
}

long long
mbdyn_count2nano(long long count)
{
	return (long long)count2nano((RTIME)count);
}

long long
mbdyn_nano2count(long long nanos)
{
	return (long long)nano2count((RTIME)nanos);
}

int
mbdyn_rt_mbx_init(const char *name, int size, void **__mbx)
{	
	assert(name);
	assert(strlen(name) <= 6);
	assert(size > 0);
	assert(__mbx != NULL);
	assert(*__mbx == NULL);		/* questa e' bastarda */
	
	*__mbx = (void *)rt_mbx_init(nam2num(name), size);

	return (*__mbx == NULL);
}

int
mbdyn_rt_mbx_delete(void **__mbx)
{
	MBX	*mbx;
	int	rc;

	assert(__mbx != NULL);
	assert(*__mbx != NULL);

	mbx = (MBX *)*__mbx;
	rc = rt_mbx_delete(mbx);

	*__mbx = NULL;

	return rc;
}

int
mbdyn_RT_get_adr(unsigned long node, int port, const char *name, void **__task)
{
	assert(node >= 0);
	assert(node == 0 || port > 0);
	assert(name != NULL);
	assert(strlen(name) == 6);
	assert(__task != NULL);
	assert(*__task == NULL);	/* questa e' bastarda */

	/* non ci va nam2num(name) */
	*__task = (void *)RT_get_adr(node, port, name);

	return (*__task == NULL);
}

int
mbdyn_RT_mbx_send_if(unsigned long node, int port, void *__mbx,
		void *msg, int msg_size)
{
	MBX	*mbx = (MBX *)__mbx;
	
	assert(__mbx != NULL);
	assert(msg != NULL);
	assert(msg_size > 0);
	
	return RT_mbx_send_if(node, port, mbx, msg, msg_size);
	
}

int
mbdyn_RT_mbx_receive_if(unsigned long node, int port, void *__mbx,
		void *msg, int msg_size)
{
	MBX	*mbx = (MBX *)__mbx;

	assert(__mbx != NULL);
	assert(msg != NULL);
	assert(msg_size > 0);

	return RT_mbx_receive_if(node, port, mbx, msg, msg_size);
}

int
mbdyn_rt_task_suspend(void *__task)
{
	RT_TASK		*task = (RT_TASK *)__task;
	
	assert(__task != NULL);

	return rt_task_suspend(task);

}

int
mbdyn_rt_task_resume(void *__task)
{
	RT_TASK		*task = (RT_TASK *)__task;
	
	assert(__task != NULL);

	return rt_task_resume(task);

}
void
mbdyn_rt_sleep(long long count)
{
		rt_sleep((RTIME)count);
}
int
mbdyn_rt_sem_init(char *name, int value, void **__sem)
{	
	assert(strlen(name) == 6);
	assert(__sem != NULL);
	assert(*__sem == NULL);	/* questa e' bastarda */
	*__sem = (void *)rt_sem_init(nam2num(name),value);

	return (*__sem == NULL);
}

int
mbdyn_rt_sem_delete(void **__sem)
{
	SEM	*sem;
	int	rc;

	assert(__sem != NULL);
	assert(*__sem != NULL);

	sem = (SEM *)*__sem;
	rc = rt_sem_delete(sem);

	*__sem = NULL;

	return rc;
}
int
mbdyn_rt_sem_signal(void *__sem)
{
	SEM	*sem = (SEM *)__sem;
	
	return rt_sem_signal(sem);
}
int
mbdyn_rt_sem_wait(void *__sem)
{
	SEM	*sem=(SEM *)__sem;
	
	return rt_sem_wait(sem);
}

void *
mbdyn_rt_receive_if(void *__task, int *msg)
{
	RT_TASK *task = (RT_TASK *)__task;
	
	return (void *)rt_receive_if(task, msg);
}

long long
mbdyn_rt_get_cpu_time_ns(void)
{
	return rt_get_cpu_time_ns();
}

#endif /* USE_RTAI */

