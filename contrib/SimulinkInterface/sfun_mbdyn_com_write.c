/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2005
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
 * COPYRIGHT (C) 2003-2004
 *
 * Michele Attolico <attolico@aero.polimi.it>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
 */

#undef VERBOSE

#define S_FUNCTION_NAME		sfun_mbdyn_com_write
#define S_FUNCTION_LEVEL	2

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif
#include "simstruc.h"

#define HOST_NAME_PARAM		ssGetSFcnParam(S, 0)
#define MBD_NAME_PARAM		ssGetSFcnParam(S, 1)
#define MBX_NAME_PARAM		ssGetSFcnParam(S, 2)
#define MBX_N_CHN_PARAM		ssGetSFcnParam(S, 3)
#define SAMPLE_TIME_PARAM	ssGetSFcnParam(S, 4)

#define NET_PARAM		ssGetSFcnParam(S, 5)
#define PORT_PARAM		ssGetSFcnParam(S, 6)
#define PATH_PARAM		ssGetSFcnParam(S, 7)

#define NUMBER_OF_PARAMS	(8)

#define WAIT_LOOP  (10000)

#define MBX_N_CHN		((uint_T) mxGetPr(MBX_N_CHN_PARAM)[0])
#define SAMPLE_TIME		((real_T) mxGetPr(SAMPLE_TIME_PARAM)[0])

#define NET			((uint_T) mxGetPr(NET_PARAM)[0])
#define PORT			((uint_T) mxGetPr(PORT_PARAM)[0])

#define NANO_SAMPLE_TIME	((long int)(SAMPLE_TIME*1000000000))

#include <errno.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/param.h>
#include <netdb.h>
#include <netinet/in.h>
#include <arpa/inet.h>

#ifndef MATLAB_MEX_FILE

#define KEEP_STATIC_INLINE
#include <math.h>
#include <rtai_lxrt.h>
#include <rtai_mbx.h>
#include <rtai_netrpc.h>
#include "mbdyn_rtai.h"

bool          MBDynTaskActive[MAX_MBDYN_TASK]={false};
long int 	MBDynNode[MAX_MBDYN_TASK]={1};
unsigned long MBDynName[MAX_MBDYN_TASK]={0xFFFFFFFF};

RT_TASK *rt_HostInterfaceTask;

extern long int 	MBDynNode[];
extern unsigned long	MBDynName[];
extern bool		MBDynTaskActive[];

static char		msg = 't';
extern RT_TASK		*rt_HostInterfaceTask;
struct mbd_t {
	void		*comptr;
	unsigned long	node;
	int		port;
	int		mbdtask_count;
};

#else

#include <string.h>
#include <stdlib.h>
#include <sys/socket.h>
#include <sys/un.h>
#include <sys/poll.h>

#define TIME_OUT	10000	/* milliseconds */

#endif /* MATLAB_MEX_FILE */

#define MDL_CHECK_PARAMETERS
#if defined(MDL_CHECK_PARAMETERS) && defined(MATLAB_MEX_FILE)
static void
mdlCheckParameters(SimStruct *S)
{
	static char_T errMsg[BUFSIZ];

	if (mxGetNumberOfElements(MBX_N_CHN_PARAM) != 1) {
		snprintf(errMsg, sizeof(errMsg),
			"Channel parameter must be a scalar.\n");
		ssSetErrorStatus(S, errMsg);
	}

	if (mxGetNumberOfElements(SAMPLE_TIME_PARAM) != 1) {
		snprintf(errMsg, sizeof(errMsg),
			"Sample time parameter must be a scalar\n");
		ssSetErrorStatus(S, errMsg);
	}

	if (mxGetNumberOfElements(PORT_PARAM) != 1) {
		snprintf(errMsg, sizeof(errMsg),
			"Port parameter must be a scalar\n");
		ssSetErrorStatus(S, errMsg);
	}

	if (NET) {
		if (PORT == 0) {
			snprintf(errMsg, sizeof(errMsg),
				"Port must be defined\n");
			ssSetErrorStatus(S, errMsg);
		}
	}
}
#endif /* MDL_CHECK_PARAMETERS && MATLAB_MEX_FILE */

static void
mdlInitializeSizes(SimStruct *S)
{
	uint_T i;

	ssSetNumSFcnParams(S, NUMBER_OF_PARAMS);
#if defined(MATLAB_MEX_FILE)
	if (ssGetNumSFcnParams(S) == ssGetSFcnParamsCount(S)) {
		mdlCheckParameters(S);
		if (ssGetErrorStatus(S) != NULL) {
			return;
		}

	} else {
		return;
	}
#endif /* MATLAB_MEX_FILE */

	for (i = 0; i < NUMBER_OF_PARAMS; i++) {
		ssSetSFcnParamNotTunable(S, i);
	}
	ssSetNumInputPorts(S, MBX_N_CHN);
	ssSetNumOutputPorts(S, MBX_N_CHN);
	for (i = 0; i < MBX_N_CHN; i++) {
		ssSetOutputPortWidth(S, i, 1);
		ssSetInputPortWidth(S, i, 1);
#ifdef MATLAB_MEX_FILE
		ssSetInputPortRequiredContiguous(S, i, 1);
    		ssSetInputPortDirectFeedThrough(S, i, 1);
#endif /*MATLAB_MEX_FILE*/
	}

	ssSetNumContStates(S, 0);
	ssSetNumDiscStates(S, 0);
	ssSetNumSampleTimes(S, 1);
	ssSetNumPWork(S, 1);

#ifdef MATLAB_MEX_FILE
	ssSetNumIWork(S, 2);
#endif /* MATLAB_MEX_FILE */
}

static void
mdlInitializeSampleTimes(SimStruct *S)
{
	ssSetSampleTime(S, 0, SAMPLE_TIME);
	ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_START
#if defined(MDL_START)
static void
mdlStart(SimStruct *S)
{
#ifndef MATLAB_MEX_FILE
	struct mbd_t		*ptrstr = NULL;
	static char_T		errMsg[BUFSIZ];
	register int		i;
	char			mbx_name[7], mbd_name[7];
	struct in_addr		addr;

	char			host_name[MAXHOSTNAMELEN];
	void			*mbdtask = NULL;
	int			timerflag = 0, count = 0;

	mxGetString(MBD_NAME_PARAM, mbd_name, 7);
	mxGetString(MBX_NAME_PARAM, mbx_name, 7);

	/******************************************************************
	 * alloc struct
	 ******************************************************************/

	ptrstr = (struct mbd_t *)malloc(sizeof(struct mbd_t));
	if (!ptrstr) {
		snprintf(errMsg, sizeof(errMsg),
			"\ncannot alloc memory for struct mbd_t\n");
		ssSetErrorStatus(S, errMsg);
		printf("%s", errMsg);
	}
	ssGetPWork(S)[0] = ptrstr;

	/******************************************************************
	 * converting host_name into node
	 ******************************************************************/
	mxGetString(HOST_NAME_PARAM, host_name, sizeof(host_name));

	inet_aton(host_name, &addr);
	ptrstr->node = addr.s_addr;
	if (ptrstr->node) {
		ptrstr->port = rt_request_hard_port(ptrstr->node);
		if (!ptrstr->port) {
			snprintf(errMsg, sizeof(errMsg),
				"rt_request_hard_port(%s) failed\n",
				ptrstr->node );
			ssSetErrorStatus(S, errMsg);
			printf("%s", errMsg);
			return;
		}

	} else {
		ptrstr->port = 0;
	}

	/* start real-time timer */
	if (!rt_is_hard_timer_running()) {
		rt_set_oneshot_mode();
		start_rt_timer(0);
		timerflag = 1;
	}

	/******************************************************************
	 * find mbdyn task
	 ******************************************************************/
	count = 0;
	printf("\nwrite: MAX_MBDYN_TASK = %d \n", MAX_MBDYN_TASK);
	while (count < MAX_MBDYN_TASK) {
		if (MBDynNode[count] == ptrstr->node
			&& MBDynName[count] == nam2num(mbd_name))
		{
			mbdtask = (void *)RT_get_adr(ptrstr->node,
				ptrstr->port, mbd_name);
			ptrstr->mbdtask_count = count;
			break;
		}

		if (MBDynNode[count] == 1 && MBDynName[count] == 0xFFFFFFFF) {
			i = WAIT_LOOP;
			printf("\nwrite: Connecting to MBDyn task\n");

			while (!(mbdtask = (void *)RT_get_adr(ptrstr->node,
				ptrstr->port,mbd_name)))
			{
				if (!i) {
					snprintf(errMsg, sizeof(errMsg),
						"\ncannot find MBDyn task\n");
					ssSetErrorStatus(S, errMsg);
					printf("%s", errMsg);
					return;
				}
				rt_sleep(nano2count(1000000));
				i--;
			}

			MBDynNode[count] = ptrstr->node;
			MBDynName[count] = nam2num(mbd_name);
			MBDynTaskActive[count] = true;
			ptrstr->mbdtask_count = count;
			break;
		}

		count++;
	}

	if (count == MAX_MBDYN_TASK) {
		snprintf(errMsg, sizeof(errMsg), "\nToo many MBDyn tasks\n");
		ssSetErrorStatus(S, errMsg);
		printf("%s", errMsg);
		return;
	}

	/******************************************************************
	 * find mailbox
	 ******************************************************************/

	i = WAIT_LOOP;

	while (!(ptrstr->comptr = (void *)RT_get_adr(ptrstr->node,
		ptrstr->port, mbx_name)))
	{
		if (!i) {
			snprintf(errMsg, sizeof(errMsg),
				"\ncannot find output mailbox %s\n",
				mbx_name);
			ssSetErrorStatus(S, errMsg);
			printf("%s", errMsg);
			return;
		}
		rt_sleep(nano2count(1000000));
		i--;
	}
	printf("\nsfun_mbdyn_com_write:\n\n");
	printf("Host name: %s\n", host_name);
	printf("node: %lx\n", ptrstr->node);
	printf("port: %d\n", ptrstr->port);
	printf("mbdtask: %p\n", mbdtask);
	printf("mbx name: %s\n", mbx_name);
	printf("mbx: %p\n", ptrstr->comptr);
	printf("mbx channel: %d\n", MBX_N_CHN);
	printf("NANO_SAMPLE_TIME: %ld\n", NANO_SAMPLE_TIME);

	if (timerflag) {
		stop_rt_timer();
	}

#else
	/* inizializza la socket */
	int sock = 0;
	int conn = 0;
	/* salva la socket */
	ssGetIWork(S)[0] = (int_T)sock;
	/* imposta conn a 0 */
	ssGetIWork(S)[1] = (int_T)conn;
#endif /* MATLAB_MEX_FILE */
}
#endif /* MDL_START */

static void
mdlOutputs(SimStruct *S, int_T tid)
{
#ifndef MATLAB_MEX_FILE
	register int	i;
	double		y[MBX_N_CHN];
	struct mbd_t	*ptrstr = (struct mbd_t *)ssGetPWork(S)[0];

	for (i = 0; i < MBX_N_CHN; i++) {
		y[i] = ssGetInputPortRealSignal(S, i)[0];
		ssGetOutputPortRealSignal(S, i)[0] = y[i];
	}

	if (RT_mbx_send_if(ptrstr->node, ptrstr->port, ptrstr->comptr, y,
		sizeof(double)*MBX_N_CHN) < 0)
	{
		MBDynTaskActive[ptrstr->mbdtask_count] = false;
		rt_send(rt_HostInterfaceTask, (int)msg);
	}
#else
	register int	i;
	double		y[MBX_N_CHN];
	int		sock = ssGetIWork(S)[0];
	int		conn = ssGetIWork(S)[1];
	static char_T	errMsg[BUFSIZ];
	int		save_errno;

	if (conn == 0) {
		if (sock == 0) {
			if (NET) {
				/* usa il protocollo tcp/ip */
				struct sockaddr_in	addr;
				char			*host = NULL;
				int			flags;

				/* crea una socket */
				sock = socket(PF_INET, SOCK_STREAM, 0);
				if (sock < 0) {
					snprintf(errMsg, sizeof(errMsg),
						"\nWRITE sfunction: unable to create a INET socket\n");
					ssSetErrorStatus(S, errMsg);
					printf("%s", errMsg);
					return;
				}

				/* imposta la socket come non bloccante */
        			flags = fcntl(sock, F_GETFL, 0);
				if (flags == -1) {
					snprintf(errMsg, sizeof(errMsg),
						"\nWRITE sfunction: unable to get socket flags\n");
					ssSetErrorStatus(S, errMsg);
					printf("%s", errMsg);
					return;
				}

				flags |= O_NONBLOCK;
				if (fcntl(sock, F_SETFL, flags) == -1) {
					snprintf(errMsg, sizeof(errMsg),
						"\nWRITE sfunction: unable to set socket flags\n");
					ssSetErrorStatus(S, errMsg);
					printf("%s", errMsg);
					return;
				}

				/* get host data */
				host = mxArrayToString(HOST_NAME_PARAM);
				if (host == NULL ) {
					snprintf(errMsg, sizeof(errMsg),
						"\nWRITE sfunction: unable to get host parameter\n");
					ssSetErrorStatus(S, errMsg);
					printf("%s", errMsg);
					return;
				}

				addr.sin_family = AF_INET;
				addr.sin_port = htons(PORT);
				if (inet_aton(host, &addr.sin_addr) == 0) {
					snprintf(errMsg, sizeof(errMsg),
						"\nWRITE sfunction: unknown host '%s'\n", host);
					ssSetErrorStatus(S, errMsg);
					printf("%s", errMsg);
					mxFree(host);
					return;
				}
				mxFree(host);

				/* connect */
				if (connect(sock, (struct sockaddr *)&addr, sizeof(addr)) != -1) {
					conn = 1;
					/* reimposta la socket come bloccante */
					flags &= (~O_NONBLOCK);
					if (fcntl(sock, F_SETFL, flags) == -1) {
						snprintf(errMsg, sizeof(errMsg),
							"\nWRITE sfunction: unable to set socket flags\n");
						ssSetErrorStatus(S, errMsg);
						printf("%s", errMsg);
						return;
					}
				}

			} else {
				/* usa le socket local */
				struct sockaddr_un	*addrp;
				char			*path = NULL;
				int			flags;
				size_t			len, size;

				/* crea una socket */
				sock = socket(PF_LOCAL, SOCK_STREAM, 0);
				if (sock < 0) {
					snprintf(errMsg, sizeof(errMsg),
						"\nWRITE sfunction: unable to create a LOCAL socket\n");
					ssSetErrorStatus(S, errMsg);
					printf("%s", errMsg);
					return;
				}

				/* imposta la socket come non bloccante */
	        		flags = fcntl(sock, F_GETFL, 0);
	        		if (flags == -1) {
					snprintf(errMsg, sizeof(errMsg),
						"\nWRITE sfunction: unable to get socket flags\n");
					ssSetErrorStatus(S, errMsg);
					printf("%s", errMsg);
					return;
				}

				flags |= O_NONBLOCK;
				if (fcntl(sock, F_SETFL, flags) == -1) {
					snprintf(errMsg, sizeof(errMsg),
						"\nWRITE sfunction: unable to set socket flags\n");
					ssSetErrorStatus(S, errMsg);
					printf("%s", errMsg);
					return;
				}

				/* get path data */
				path = mxArrayToString(PATH_PARAM);
				if (path == NULL ) {
					snprintf(errMsg, sizeof(errMsg),
						"\nWRITE sfunction: unable to get path parameter\n");
					ssSetErrorStatus(S, errMsg);
					printf("%s", errMsg);
					return;
				}

				len = strlen(path);
				addrp = malloc(sizeof(struct sockaddr_un) + len + 1);
				addrp->sun_family = AF_LOCAL;
				/* Note: sun_path has a fixed length (108 bytes), but it is at the end of sockaddr_un;
				   by malloc'ing it of the correct size (actually, way in excess), we can use paths
				   of arbitrary length */
				strncpy(addrp->sun_path, path, len);
				mxFree(path);

				addrp->sun_path[len] = '\0';
				size = (offsetof (struct sockaddr_un, sun_path) + len + 1);

				/* connect */
				if (connect(sock, (struct sockaddr *)addrp, size) == 0) {
					conn = 1;
					/* reimposta la socket come bloccante */
					flags &= (~O_NONBLOCK);
					if (fcntl(sock, F_SETFL, flags) == -1) {
						snprintf(errMsg, sizeof(errMsg),
							"\nWRITE sfunction: unable to set socket flags\n");
						ssSetErrorStatus(S, errMsg);
						printf("%s", errMsg);
						return;
					}
				}
			}

			/* salva la socket */
			ssGetIWork(S)[0] = (int_T)sock;

		} else {
			/* poll */
			struct pollfd	ufds;
			int		rc;
			int		retries = 0;

retry:;
			ufds.fd = sock;
			ufds.events = (POLLIN | POLLOUT);
			rc = poll(&ufds, 1, TIME_OUT);
			switch (rc) {
			case -1:
				save_errno = errno;
				snprintf(errMsg, sizeof(errMsg),
					"\nWRITE sfunction: POLL error (%d: %s)\n",
					save_errno, strerror(save_errno));
				ssSetErrorStatus(S, errMsg);
				printf("%s", errMsg);
				return;

			case 0:
				snprintf(errMsg, sizeof(errMsg),
					"\nWRITE sfunction: connection timeout reached\n");
				ssSetErrorStatus(S, errMsg);
				printf("%s", errMsg);
				return;

			default :
			{
				int	flags;

				if ((ufds.revents & POLLHUP) && ++retries < 200) {
					usleep(10000);
					goto retry;
				}

				if (ufds.revents & (POLLERR | POLLHUP | POLLNVAL)) {
					snprintf(errMsg, sizeof(errMsg),
						"\nWRITE sfunction: POLL error "
						"(%d ERR=%d HUP=%d NVAL=%d)\n",
						rc,
						ufds.revents & POLLERR,
						ufds.revents & POLLHUP,
						ufds.revents & POLLNVAL);
					ssSetErrorStatus(S, errMsg);
					printf("%s", errMsg);
					return;
				}

				conn = 1;
				/* reimposta la socket come bloccante */
        			flags = fcntl(sock, F_GETFL, 0);
        			if (flags == -1) {
					snprintf(errMsg, sizeof(errMsg),
						"\nWRITE sfunction: unable to get socket flags\n");
					ssSetErrorStatus(S, errMsg);
					printf("%s", errMsg);
					return;
				}
				flags &= (~O_NONBLOCK);
				if (fcntl(sock, F_SETFL, flags) == -1) {
					snprintf(errMsg, sizeof(errMsg),
						"\nWRITE sfunction: unable to set socket flags\n");
					ssSetErrorStatus(S, errMsg);
					printf("%s", errMsg);
					return;
				}
			}
			}
			/* imposta conn a 1 */

		} /* socket == 0 */
		ssGetIWork(S)[1] = (int_T)conn;
	} /* conn == 0 */

	if (conn) {
		/* scrive sulla socket */
		for (i = 0; i < MBX_N_CHN; i++) {
			y[i] = ssGetInputPortRealSignal(S, i)[0];
			ssGetOutputPortRealSignal(S, i)[0] = y[i];
		}
		if (send(sock, y, sizeof(double)*MBX_N_CHN, 0) == -1) {
			snprintf(errMsg, sizeof(errMsg),
				"\nWriteSfunction: Communication closed by host\n");
			ssSetStopRequested(S, 1);
			printf("%s", errMsg);
			return;
		}
	}
#endif /* MATLAB_MEX_FILE */
}

static void
mdlTerminate(SimStruct *S)
{
#ifndef MATLAB_MEX_FILE
	struct 	mbd_t *ptrstr = (struct mbd_t *)ssGetPWork(S)[0];
	void	*mbdtask;
	char 	mbd_name[7];

	printf("\nsfun_mbdyn_com_write, %s:\n", ssGetModelName(S));

	if (ptrstr) {
		/****************************************************************
		 * terminate mbdyn task
		 ****************************************************************/
		if (MBDynTaskActive[ptrstr->mbdtask_count]) {
			mxGetString(MBD_NAME_PARAM, mbd_name, 7);
			mbdtask = (void *)RT_get_adr(ptrstr->node,
				ptrstr->port, mbd_name);
			printf("MBDyn is stopped by write s-function\n");
			RT_send_timed(ptrstr->node, ptrstr->port, mbdtask, 1, NANO_SAMPLE_TIME*2);
			MBDynTaskActive[ptrstr->mbdtask_count] = false;
		}
		/****************************************************************
		 * release port
		 ****************************************************************/

		if (ptrstr->node) {
			rt_release_port(ptrstr->node, ptrstr->port);
		}

		/****************************************************************
		 * free structure memory
		 ****************************************************************/
		free(ptrstr);
	}
#else
	int sock = ssGetIWork(S)[0];
#ifdef VERBOSE
	fprintf(stderr, "WRITE: closing socket...\n");
#endif /* VERBOSE */
	/* chiude la socket */
	shutdown(sock, SHUT_RDWR);
#endif /* MATLAB_MEX_FILE */
}

#ifdef  MATLAB_MEX_FILE
#include "simulink.c"
#else
#include "cg_sfun.h"
#endif

