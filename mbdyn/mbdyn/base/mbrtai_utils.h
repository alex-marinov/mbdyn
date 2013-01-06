/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

#ifndef MBRTAI_UTILS_H
#define MBRTAI_UTILS_H

// NOTE: all calls to RTAI are wrapped by MBDyn calls prefixed by "rtmbdyn_"

#ifdef USE_RTAI

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

extern int rtmbdyn_rt_task_init(const char *name, int priority, int stack_size,
	int max_msg_size, int cpu, void **task);
extern int rtmbdyn_rt_task_delete(void **task);

extern	void rtmbdyn_rt_make_hard_real_time(void);
extern void rtmbdyn_rt_make_soft_real_time(void);

extern void rtmbdyn_rt_allow_nonroot_hrt(void);

extern void rtmbdyn_rt_set_oneshot_mode(void);
extern void rtmbdyn_rt_set_periodic_mode(void);

extern int rtmbdyn_rt_is_hard_timer_running(void);
extern long long rtmbdyn_start_rt_timer(long long period);
extern void rtmbdyn_stop_rt_timer(void);
extern long long rtmbdyn_rt_get_time(void);
extern int rtmbdyn_rt_task_make_periodic(void *task, long long start_time,
	long long period);
extern void rtmbdyn_rt_task_wait_period(void);
extern long long rtmbdyn_count2nano(long long count);
extern long long rtmbdyn_nano2count(long long nanos);

extern int rtmbdyn_rt_request_port(unsigned long node);

extern int rtmbdyn_rt_mbx_init(const char *name, int size, void **mbx);
extern int rtmbdyn_rt_mbx_delete(void **mbx);
extern int rtmbdyn_RT_mbx_send(unsigned long node, int port, void *mbx,
		void *msg, int msg_size);
extern int rtmbdyn_RT_mbx_send_if(unsigned long node, int port, void *mbx,
		void *msg, int msg_size);
extern int rtmbdyn_RT_mbx_receive(unsigned long node, int port, void *mbx,
		void *msg, int msg_size);
extern int rtmbdyn_RT_mbx_receive_if(unsigned long node, int port, void *mbx,
		void *msg, int msg_size);
		
extern int rtmbdyn_RT_get_adr(unsigned long node, int port, const char *name, 
		void **task);
		
extern int rtmbdyn_rt_task_suspend(void *task);
extern int rtmbdyn_rt_task_resume(void *task);
extern void rtmbdyn_rt_sleep(long long count);

extern int rtmbdyn_rt_sem_init(char *name, int value, void **sem);
extern int rtmbdyn_rt_sem_delete(void **sem);
extern int rtmbdyn_rt_sem_signal(void *sem);
extern int rtmbdyn_rt_sem_wait(void *sem);

extern void *rtmbdyn_rt_receive(void *task,int *msg);
extern void *rtmbdyn_rt_receive_if(void *task,int *msg);
extern long long rtmbdyn_rt_get_cpu_time_ns(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* USE_RTAI */

#endif /* MBRTAI_UTILS_H */

