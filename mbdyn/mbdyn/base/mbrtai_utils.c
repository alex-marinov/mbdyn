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

/* socket driver */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#ifdef USE_RTAI

/*
 * Mante prefers KEEP_STATIC_INLINE instead of linking liblxrt.a
 */
#define KEEP_STATIC_INLINE
#include <rtai_lxrt_user.h>
#include <net_rpc.h>

//#define NDEBUG
#include <assert.h>

#include <mbrtai_utils.h>

int
mbdyn_rt_task_init(const char *name, int priority, int stack_size,
		int max_msg_size, void **__task)
{
	assert(name != NULL);
	assert(strlen(name) == 6);
	assert(__task != NULL);
	assert(*__task == NULL);

	*__task = (void *)rt_task_init(nam2num(name), priority,
			stack_size, max_msg_size);

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

int
mbdyn_rt_request_port(unsigned long node)
{
	assert(node > 0);

	return rt_request_port(node);
}

long long
mbdyn_rt_get_time(void)
{
	return (long long)rt_get_time();
}

long long
mbdyn_count2nano(long long count)
{
	return (RTIME)count2nano((RTIME)count);
}

long long
mbdyn_nano2count(long long nanos)
{
	return (RTIME)nano2count((RTIME)nanos);
}

int
mbdyn_rt_mbx_init(const char *name, int size, void **__mbx)
{
	assert(name);
	assert(strlen(name) == 6);
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

	return /* rc */ 0;
}

int
mbdyn_RT_get_adr(unsigned long node, int port, const char *name, void **__task)
{
	assert(node > 0);
	assert(port > 0);
	assert(name != NULL);
	assert(strlen(name) == 6);
	assert(__task != NULL);
	assert(*__task == NULL);	/* questa e' bastarda */

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

#endif /* USE_RTAI */

