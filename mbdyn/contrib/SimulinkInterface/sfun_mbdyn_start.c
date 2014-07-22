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

#define S_FUNCTION_NAME  sfun_mbdyn_start
#define S_FUNCTION_LEVEL 2

/*
 * Need to include simstruc.h for the definition of the SimStruct and
 * its associated macro definitions.
 */
#include "simstruc.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <signal.h>


#define N_OF_PARAM 7

#define MBDYN_PATH_NAME_PARAM	ssGetSFcnParam(S, 0)
#define MAX_MBDYN_PATH_NAME_DIM 256
#define MBDYN_PATH_PARAM	ssGetSFcnParam(S, 1)
#define MBDYN_PATH		((uint_T) mxGetPr(MBDYN_PATH_PARAM)[0])

#define FILE_NAME_PARAM		ssGetSFcnParam(S, 2)
#define MAX_FILE_NAME_DIM	256

#define FILE_OUTPUT_PARAM	ssGetSFcnParam(S, 3)
#define MAX_FILE_OUTPUT_DIM	256
#define OUTPUT_PARAM		ssGetSFcnParam(S, 4)
#define OUTPUT			((uint_T) mxGetPr(OUTPUT_PARAM)[0])

#define VERBOSE_PARAM		ssGetSFcnParam(S, 5)
#define VERBOSE			((uint_T) mxGetPr(VERBOSE_PARAM)[0])

#define PEDANTIC_PARAM		ssGetSFcnParam(S, 6)
#define PEDANTIC		((uint_T) mxGetPr(PEDANTIC_PARAM)[0])


static void
mdlInitializeSizes(SimStruct *S)
{
	/* Number of expected parameters */
	ssSetNumSFcnParams(S, N_OF_PARAM);
	if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
		/* Return if # of expected != # of actual parameters */
		printf("S-function params %d do not match expected %d\n",
				ssGetNumSFcnParams(S), ssGetSFcnParamsCount(S));
		return;
	}

	ssSetNumContStates(S, 0);
	ssSetNumDiscStates(S, 0);

	if (!ssSetNumInputPorts(S, 0)) {
		printf("No input ports defined\n");
		return;
	}
	
	/*
	 * Set direct feedthrough flag (1=yes, 0=no).
	 * A port has direct feedthrough if the input is used in either
	 * the mdlOutputs or mdlGetTimeOfNextVarHit functions.
	 * See matlabroot/simulink/src/sfuntmpl_directfeed.txt.
	 */
    
	if (!ssSetNumOutputPorts(S, 0)) {
		printf("No output ports defined\n");
		return;
	}
 
	ssSetNumSampleTimes(S, 1);
	ssSetNumRWork(S, 0);
	ssSetNumIWork(S, 1);
	ssSetNumPWork(S, 0);
	ssSetNumModes(S, 0);
	ssSetNumNonsampledZCs(S, 0);

	ssSetOptions(S, 0);
}

static void
mdlInitializeSampleTimes(SimStruct *S)
{
	ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
	ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_START 
#if defined(MDL_START)
static void
mdlStart(SimStruct *S)
{
#ifdef MATLAB_MEX_FILE
	int pid;

	pid = fork();

	switch (pid) {
	case 0 :
	{
		int		count = 0;
		char		*parameter[8] = { '\0' };
		char		mbdyn_name[MAX_MBDYN_PATH_NAME_DIM] = { '\0' };
		char		file_output[MAX_FILE_OUTPUT_DIM] = { '\0' };
		char		file_name[MAX_FILE_NAME_DIM] = { '\0' };
		char		file_option[] = "-f";
		char		output_option[] = "-o";
		char		verbose_option[] = "-ss";
		char		pedantic_option[] = "-P";
		struct stat	st;
		static char_T	errMsg[BUFSIZ];
		int		executable = 1;

		if (MBDYN_PATH) {
			mxGetString(MBDYN_PATH_NAME_PARAM, mbdyn_name,
					sizeof(mbdyn_name));
		} else {
			strcpy(mbdyn_name, "./mbdyn.sh");
		}
		parameter[count] = mbdyn_name;
		count++;	

		if (stat(mbdyn_name, &st) == -1) {
			int	save_errno = errno;
			char	*msg = strerror(errno);

			if (msg == NULL) {
				snprintf(errMsg, sizeof(errMsg),
					"stat(%s) failed and "
					"strerror(%d) too\n",
					mbdyn_name, save_errno);

			} else {
				snprintf(errMsg, sizeof(errMsg),
					"stat(%s) failed: %d (%s)\n",
					mbdyn_name, save_errno, msg);
			}
			ssSetErrorStatus(S, errMsg);
			fputs(errMsg, stderr);
			exit(EXIT_FAILURE);
		}

		if (st.st_uid == getuid() && !(st.st_mode & S_IXUSR)) {
			/* not executable by owner */
			executable = 0;

		} else if (st.st_gid == getgid() && !(st.st_mode & S_IXGRP)) {
			/* not executable by group */
			executable = 0;

		} else if (!(st.st_mode & S_IXOTH)) {
			/* not executable by other */
			executable = 0;
			
		}

		if (!executable) {
			/* not executable */
			snprintf(errMsg, sizeof(errMsg),
				"program %s is not executable\n",
				mbdyn_name);
			ssSetErrorStatus(S, errMsg);
			fputs(errMsg, stderr);
			exit(EXIT_FAILURE);
		}

		mxGetString(FILE_NAME_PARAM, file_name, sizeof(file_name));
		parameter[count] = file_option;
		count++;		
		parameter[count] = file_name;
		count++;
		
		if (OUTPUT) {
			parameter[count] = output_option;
			count++;
			mxGetString(FILE_OUTPUT_PARAM, file_output,
					sizeof(file_output));
			parameter[count] = file_output;
			count++;
		}
		
		if (!VERBOSE) {
			parameter[count] = verbose_option;
			count++;				
		}
		
		if (PEDANTIC) {
			parameter[count] = pedantic_option;
			count++;				
		}
		parameter[count] = NULL;
		
		if (execv(mbdyn_name, parameter) == -1) {
			char		*msg;
			int		save_errno = errno;

			msg = strerror(save_errno);
			if (msg == NULL) {
				fprintf(stderr, "execv(%s) failed, "
						"and strerror_r too\n",
						parameter[0]);
			} else {
				fprintf(stderr, "execv(%s) failed: %d (%s)\n",
						parameter[0], save_errno, msg);
			}
			exit(EXIT_FAILURE);
		}
	}
	
	case -1:
	{
		char		*msg;
		int		save_errno = errno;
		static char_T	errMsg[BUFSIZ];

		msg = strerror(save_errno);
		if (msg == NULL) {
			snprintf(errMsg, sizeof(errMsg),
				"fork() failed, and strerror too\n");
		} else {
			snprintf(errMsg, sizeof(errMsg),
				"fork() failed: %d (%s)\n",
				save_errno, msg);
		}
		ssSetErrorStatus(S, errMsg);
		break;
	}

	default:
	{
		int sleeptime = 0;
		printf("Simulink start MBDyn task\n", sleeptime);
		if (sleeptime) {
			printf("\tsleeping %d s\n", sleeptime);
			while (sleeptime) {
				sleeptime = sleep(sleeptime);
			}
		}
		break;
	}
	}

	ssGetIWork(S)[0] = (int_T)pid;
#endif /*MATLAB_MEX_FILE*/
}
#endif /*MDL_START*/

static void
mdlOutputs(SimStruct *S, int_T tid)
{
	;
}


static void mdlTerminate(SimStruct *S)
{
#ifdef MATLAB_MEX_FILE
	int pid = (int_T)ssGetIWork(S)[0];
#endif /*MATLAB_MEX_FILE*/
}



/*=============================*
 * Required S-function trailer *
 *=============================*/

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
