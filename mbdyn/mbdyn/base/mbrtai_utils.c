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

/* socket driver */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

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
rtmbdyn_rt_task_init(const char *name, int priority, int stack_size,
		int max_msg_size, int cpu, void **v_task)
{
	assert(name != NULL);
	assert(strlen(name) == 6);
	assert(v_task != NULL);
	assert(*v_task == NULL);

	*v_task = (void *)rt_task_init_schmod(nam2num(name), priority,
			stack_size, max_msg_size, SCHED_FIFO,/*0x2*/ cpu);
	return (*v_task == NULL);
}

int
rtmbdyn_rt_task_delete(void **v_task)
{
	RT_TASK		*task;
	int		rc;

	assert(v_task != NULL);
	assert(*v_task != NULL);

	task = (RT_TASK *)*v_task;
	rc = rt_task_delete(task);

	*v_task = NULL;

	return /* rc */ 0;
}
void
rtmbdyn_rt_make_hard_real_time(void)
{
	rt_make_hard_real_time();
}
void
rtmbdyn_rt_make_soft_real_time(void)
{
	rt_make_soft_real_time();
}
int
rtmbdyn_rt_task_make_periodic(void *v_task, long long v_start_time,
		long long v_period)
{
	
	
	RT_TASK		*task = (RT_TASK *)v_task;
	RTIME		start_time = (RTIME)v_start_time;
	RTIME		period = (RTIME)v_period;

	assert(v_task != NULL);

	return rt_task_make_periodic(task, start_time, period);
}

void
rtmbdyn_rt_task_wait_period(void)
{
	rt_task_wait_period();
}

void
rtmbdyn_rt_allow_nonroot_hrt(void)
{
	rt_allow_nonroot_hrt();
}

int
rtmbdyn_rt_request_port(unsigned long node)
{
	assert(node > 0);

	return rt_request_port(node);
}

void
rtmbdyn_rt_set_oneshot_mode(void)
{
	rt_set_oneshot_mode();
}

void
rtmbdyn_rt_set_periodic_mode(void)
{
	rt_set_periodic_mode();
}

int
rtmbdyn_rt_is_hard_timer_running(void)
{
	return rt_is_hard_timer_running();
}

long long
rtmbdyn_start_rt_timer(long long v_period)
{
	RTIME period = (RTIME)v_period;
	return start_rt_timer(period);
}

void
rtmbdyn_stop_rt_timer(void)
{
	stop_rt_timer();
}

long long
rtmbdyn_rt_get_time(void)
{	
	return (long long)rt_get_time();
}

long long
rtmbdyn_count2nano(long long count)
{
	return (long long)count2nano((RTIME)count);
}

long long
rtmbdyn_nano2count(long long nanos)
{
	return (long long)nano2count((RTIME)nanos);
}

int
rtmbdyn_rt_mbx_init(const char *name, int size, void **v_mbx)
{	
	assert(name);
	assert(strlen(name) <= 6);
	assert(size > 0);
	assert(v_mbx != NULL);
	assert(*v_mbx == NULL);		/* questa e' bastarda */
	
	*v_mbx = (void *)rt_mbx_init(nam2num(name), size);

	return (*v_mbx == NULL);
}

int
rtmbdyn_rt_mbx_delete(void **v_mbx)
{
	MBX	*mbx;
	int	rc;

	assert(v_mbx != NULL);
	assert(*v_mbx != NULL);

	mbx = (MBX *)*v_mbx;
	rc = rt_mbx_delete(mbx);

	*v_mbx = NULL;

	return rc;
}

int
rtmbdyn_RT_get_adr(unsigned long node, int port, const char *name, void **v_task)
{
	assert(node >= 0);
	assert(node == 0 || port > 0);
	assert(name != NULL);
	assert(strlen(name) == 6);
	assert(v_task != NULL);
	assert(*v_task == NULL);	/* questa e' bastarda */

	/* non ci va nam2num(name) */
	*v_task = (void *)RT_get_adr(node, port, name);

	return (*v_task == NULL);
}

int
rtmbdyn_RT_mbx_send(unsigned long node, int port, void *v_mbx,
		void *msg, int msg_size)
{
	MBX	*mbx = (MBX *)v_mbx;
	
	assert(v_mbx != NULL);
	assert(msg != NULL);
	assert(msg_size > 0);
	
	return RT_mbx_send(node, port, mbx, msg, msg_size);
	
}

int
rtmbdyn_RT_mbx_send_if(unsigned long node, int port, void *v_mbx,
		void *msg, int msg_size)
{
	MBX	*mbx = (MBX *)v_mbx;
	
	assert(v_mbx != NULL);
	assert(msg != NULL);
	assert(msg_size > 0);
	
	return RT_mbx_send_if(node, port, mbx, msg, msg_size);
	
}

int
rtmbdyn_RT_mbx_receive(unsigned long node, int port, void *v_mbx,
		void *msg, int msg_size)
{
	MBX	*mbx = (MBX *)v_mbx;

	assert(v_mbx != NULL);
	assert(msg != NULL);
	assert(msg_size > 0);

	return RT_mbx_receive(node, port, mbx, msg, msg_size);
}

int
rtmbdyn_RT_mbx_receive_if(unsigned long node, int port, void *v_mbx,
		void *msg, int msg_size)
{
	MBX	*mbx = (MBX *)v_mbx;

	assert(v_mbx != NULL);
	assert(msg != NULL);
	assert(msg_size > 0);

	return RT_mbx_receive_if(node, port, mbx, msg, msg_size);
}

int
rtmbdyn_rt_task_suspend(void *v_task)
{
	RT_TASK		*task = (RT_TASK *)v_task;
	
	assert(v_task != NULL);

	return rt_task_suspend(task);

}

int
rtmbdyn_rt_task_resume(void *v_task)
{
	RT_TASK		*task = (RT_TASK *)v_task;
	
	assert(v_task != NULL);

	return rt_task_resume(task);

}
void
rtmbdyn_rt_sleep(long long count)
{
		rt_sleep((RTIME)count);
}
int
rtmbdyn_rt_sem_init(char *name, int value, void **v_sem)
{	
	assert(strlen(name) == 6);
	assert(v_sem != NULL);
	assert(*v_sem == NULL);	/* questa e' bastarda */
	*v_sem = (void *)rt_sem_init(nam2num(name),value);

	return (*v_sem == NULL);
}

int
rtmbdyn_rt_sem_delete(void **v_sem)
{
	SEM	*sem;
	int	rc;

	assert(v_sem != NULL);
	assert(*v_sem != NULL);

	sem = (SEM *)*v_sem;
	rc = rt_sem_delete(sem);

	*v_sem = NULL;

	return rc;
}
int
rtmbdyn_rt_sem_signal(void *v_sem)
{
	SEM	*sem = (SEM *)v_sem;
	
	return rt_sem_signal(sem);
}
int
rtmbdyn_rt_sem_wait(void *v_sem)
{
	SEM	*sem=(SEM *)v_sem;
	
	return rt_sem_wait(sem);
}

void *
rtmbdyn_rt_receive(void *v_task, int *msg)
{
	RT_TASK *task = (RT_TASK *)v_task;
	
	return (void *)rt_receive(task, msg);
}

void *
rtmbdyn_rt_receive_if(void *v_task, int *msg)
{
	RT_TASK *task = (RT_TASK *)v_task;
	
	return (void *)rt_receive_if(task, msg);
}

long long
rtmbdyn_rt_get_cpu_time_ns(void)
{
	return rt_get_cpu_time_ns();
}

#endif /* USE_RTAI */

