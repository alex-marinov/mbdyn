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

#ifndef RTAI_UTILS_H
#define RTAI_UTILS_H

#ifdef USE_RTAI

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

extern int mbdyn_rt_task_init(const char *name, int priority, int stack_size,
	int max_msg_size, void **__task);
extern int mbdyn_rt_task_delete(void **__task);
extern void mbdyn_rt_set_oneshot_mode(void);
extern void mbdyn_rt_set_periodic_mode(void);
extern int mbdyn_rt_is_hard_timer_running(void);
extern long long mbdyn_start_rt_timer(long long __period);
extern void mbdyn_stop_rt_timer(void);
extern int mbdyn_rt_task_make_periodic(void *__task, long long __start_time,
	long long __period);

extern void mbdyn_rt_task_wait_period(void);
extern void mbdyn_rt_allow_nonroot_hrt(void);

extern long long mbdyn_rt_get_time(void);
extern long long mbdyn_count2nano(long long count);
extern long long mbdyn_nano2count(long long nanos);

extern int mbdyn_rt_request_port(unsigned long node);

extern int mbdyn_rt_mbx_init(const char *name, int size, void **__mbx);
extern int mbdyn_rt_mbx_delete(void **__mbx);
extern int mbdyn_RT_mbx_send_if(unsigned long node, int port, void *__mbx,
		void *msg, int msg_size);
extern int mbdyn_RT_mbx_receive_if(unsigned long node, int port, void *__mbx,
		void *msg, int msg_size);
extern int mbdyn_RT_get_adr(unsigned long node, int port, const char *name, 
		void **__task);
extern int mbdyn_rt_task_suspend(void *__task);
extern void mbdyn_rt_sleep(long long count);
extern int mbdyn_rt_sem_init(char *name, int value, void **__sem);
extern int mbdyn_rt_sem_delete(void **__sem);
extern int mbdyn_rt_sem_signal(void *__sem);
extern int mbdyn_rt_sem_wait(void *__sem);
#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* USE_RTAI */

#endif /* RTAI_UTILS_H */

