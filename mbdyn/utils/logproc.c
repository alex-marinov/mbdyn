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

/*
 * Author: Michele Attolico <attolico@aero.polimi.it>
 */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <stdio.h>
#include <stdlib.h>

#ifdef USE_RTAI

#include <unistd.h>
#include <sys/types.h>
#include <sys/user.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sched.h>

#define KEEP_STATIC_INLINE

#if defined (HAVE_RTAI_LXRT_USER_H)
#include <rtai_lxrt_user.h>
#elif defined (HAVE_RTAI_LXRT_H)
#include <rtai_lxrt.h>
#include <rtai_mbx.h>
#else
#error "No rtai_lxrt_user and no rtai_lxrt within rtai"
#endif

int main(int argc, char* argv[])
{
	RT_TASK *logtask, *mbdtask;
	MBX *mbxlog;
	int i = 0, time = 0;
	char 	*mbdynname = argv[1];
	char	*mbxname = argv[2];
	int 	CpuMap = atoi(argv[3]);
	struct dati{
		int step;
		int time;
	};
	int dim =sizeof(struct dati);
	
	struct dati msg;
	
	
	struct sched_param mysched;

	mysched.sched_priority = 1;

	if( sched_setscheduler( 0, SCHED_FIFO, &mysched ) == -1 ) {
	puts(" ERROR IN SETTING THE SCHEDULER UP");
	perror( "errno" );
	exit( 0 );
 	}       
	
	mlockall(MCL_CURRENT | MCL_FUTURE);
 	
	if (!(logtask = rt_task_init_schmod(nam2num("LTSK"), 20, 0, 0,
			SCHED_FIFO, CpuMap))) {
		printf("CANNOT INIT LOG TASK\n");
		exit(1);
	}
	
	 if (!(mbdtask = rt_get_adr(nam2num(mbdynname)))) {
		printf("CANNOT FIND MBDyn TASK\n");
		exit(1);
	}
	
	rt_task_resume(mbdtask);
	
	if (!( mbxlog= rt_get_adr(nam2num(mbxname)))) {
		printf("CANNOT FIND LOG MBX\n");
		exit(1);
	}

	printf("\nOVERRUNS MONITOR:\n");
	printf("     step    t[micro s]\n");
	while (!rt_mbx_receive(mbxlog, &msg, dim)){
		i++;
		time += msg.time;
		printf("%3d %5d %10d\n", i, msg.step, msg.time);
		
	}
	rt_sleep(nano2count(1000000000));
	rt_task_delete(logtask);
	//sleep(1);
	printf("\n\nOVERRUNS MONITOR:\n");
	printf("Total overruns detected: %d\n",i);
	printf("Mean overruns time: %6.2lf micro s\n",(double)time/(double)i);
	printf("End of overruns monitor.\n");

	return 0;	
}

#else /* !USE_RTAI */
int
main(void)
{
	fprintf(stderr, "need RTAI support\n");
	exit(EXIT_FAILURE);
}
#endif /* !USE_RTAI */

