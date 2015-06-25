/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <stdio.h>
#include <stdlib.h>

#include <unistd.h>
#include <sys/types.h>
#include <sys/user.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sched.h>
#include <errno.h>

#define KEEP_STATIC_INLINE

#if defined (HAVE_RTAI_LXRT_H)
#include <rtai_lxrt.h>
#include <rtai_mbx.h>
#else
#error "No rtai_lxrt_user and no rtai_lxrt within rtai"
#endif

int
main(int argc, char* argv[])
{
	RT_TASK *logtask, *mbdtask;
	MBX *mbxlog;
	int i = 0, time = 0;
	char *mbdynname;
	char *mbxname;
	int CpuMap;
	struct data {
		int step;
		int time;
	};
	int dim = sizeof(struct data);
	struct data msg;
	
	struct sched_param mysched;

	mysched.sched_priority = 1;

	if (argc != 5) {
		fprintf(stderr, "INVALID NUMBER OF ARGUMENTS: ARGC=%d\n"
			"\n"
			"USAGE: logproc <MBDTSK> <MBXLOG> <CPUMAP> <NONROOT>\n",
			argc);
		exit(EXIT_FAILURE);
	}

	mbdynname = argv[1];
	mbxname = argv[2];
	CpuMap = atoi(argv[3]);
	if (strcasecmp(argv[4], "true") == 0) {
		rt_allow_nonroot_hrt();

	} else if (strcasecmp(argv[4], "false" ) != 0) {
		fprintf(stderr, "INVALID VALUE \"%s\" FOR NONROOT\n", argv[4]);
		exit(EXIT_FAILURE);
	}

	if (sched_setscheduler(0, SCHED_FIFO, &mysched) == -1) {
		fputs("ERROR IN SETTING THE SCHEDULER UP", stderr);
		perror("errno");
		exit( 0 );
 	}       
	
	mlockall(MCL_CURRENT | MCL_FUTURE);
 	
	if (!(logtask = rt_task_init_schmod(nam2num("LTSK"), 20, 0, 0,
			SCHED_FIFO, CpuMap)))
	{
		fputs("CANNOT INIT LOG TASK\n", stderr);
		exit(EXIT_FAILURE);
	}
	
	if (!(mbdtask = rt_get_adr(nam2num(mbdynname)))) {
		fprintf(stderr, "CANNOT FIND MBDyn TASK \"%s\"\n", mbdynname);
		exit(EXIT_FAILURE);
	}
	
	rt_task_resume(mbdtask);
	
	if (!(mbxlog = rt_get_adr(nam2num(mbxname)))) {
		fprintf(stderr, "CANNOT FIND LOG MBX\n");
		exit(EXIT_FAILURE);
	}

	printf("\n" "OVERRUNS MONITOR:\n");
	printf("             step     t [ns]\n");
	while (!rt_mbx_receive(mbxlog, &msg, dim)) {
		i++;
		time += msg.time;
		printf("%8d %8d %10d\n", i, msg.step, msg.time);
	}

	rt_sleep(nano2count(1000000000));
	rt_task_delete(logtask);
	printf("\n\n" "OVERRUNS MONITOR:\n");
	printf("Total overruns detected: %d\n", i);
	printf("Mean overruns time: %6.2lf ns\n",
		i ? ((double)time)/((double)i) : 0.);
	printf("End of overruns monitor.\n");

	return 0;	
}

