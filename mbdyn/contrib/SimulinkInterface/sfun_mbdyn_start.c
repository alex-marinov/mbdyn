/*
COPYRIGHT (C) 2003  Michele Attolico (attolico@aero.polimi.it)

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
*/

#define S_FUNCTION_NAME  sfun_mbdyn_start
#define S_FUNCTION_LEVEL 2

/*
 * Need to include simstruc.h for the definition of the SimStruct and
 * its associated macro definitions.
 */
#include "simstruc.h"
#include <sys/types.h>
#include <unistd.h>
#include <signal.h>


#define N_OF_PARAM 7

#define MBDYN_PATH_NAME_PARAM	ssGetSFcnParam(S,0)
#define MAX_MBDYN_PATH_NAME_DIM 200
#define MBDYN_PATH_PARAM	ssGetSFcnParam(S,1)
#define MBDYN_PATH		((uint_T) mxGetPr(MBDYN_PATH_PARAM)[0])

#define FILE_NAME_PARAM		ssGetSFcnParam(S,2)
#define MAX_FILE_NEME_DIM	100

#define FILE_OUTPUT_PARAM	ssGetSFcnParam(S,3)
#define MAX_FILE_OUTPUT_DIM	20
#define OUTPUT_PARAM		ssGetSFcnParam(S,4)
#define OUTPUT			((uint_T) mxGetPr(OUTPUT_PARAM)[0])

#define VERBOSE_PARAM		ssGetSFcnParam(S,5)
#define VERBOSE			((uint_T) mxGetPr(VERBOSE_PARAM)[0])

#define PEDANTIC_PARAM		ssGetSFcnParam(S,6)
#define PEDANTIC		((uint_T) mxGetPr(PEDANTIC_PARAM)[0])


static void mdlInitializeSizes(SimStruct *S)
{

  	ssSetNumSFcnParams(S, N_OF_PARAM);  /* Number of expected parameters */
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        /* Return if number of expected != number of actual parameters */
        return;
    }

    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);

    if (!ssSetNumInputPorts(S, 0)) return;
    /*
     * Set direct feedthrough flag (1=yes, 0=no).
     * A port has direct feedthrough if the input is used in either
     * the mdlOutputs or mdlGetTimeOfNextVarHit functions.
     * See matlabroot/simulink/src/sfuntmpl_directfeed.txt.
     */
    
    if (!ssSetNumOutputPorts(S, 0)) return;
 
    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 1);
    ssSetNumPWork(S, 0);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);

    ssSetOptions(S, 0);
}

static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);

}


#define MDL_START 
#if defined(MDL_START)
static void mdlStart(SimStruct *S)
{
#ifdef MATLAB_MEX_FILE
	int pid;
	printf("Qui!\n");
	switch(pid = fork()) {
	case 0 :{
		int count = 0;
		char *parameter[8];
		char * mbdyn_name;
		char * file_output;
		char file_name[MAX_FILE_NEME_DIM];
		char file_option[] = "-f";
		char output_option[] = "-o";
		char verbose_option[] = "-ss";
		char pedantic_option[] = "-P";

		fprintf(stderr,"qui\n");

		if(MBDYN_PATH) {
			mbdyn_name = calloc(MAX_MBDYN_PATH_NAME_DIM, sizeof(char));
			mxGetString(MBDYN_PATH_NAME_PARAM, mbdyn_name, MAX_MBDYN_PATH_NAME_DIM);	
			parameter[count] = mbdyn_name;
			count ++;	
		} else {
			mbdyn_name= calloc(8, sizeof(char));
			memcpy(mbdyn_name,"./mbdyn\0",8);
			parameter[count] = mbdyn_name;
			count ++;
		}
		/*char mbdyn_name [] ="../mbdyn-devel1/mbdyn/mbdyn";
		count ++;*/
		mxGetString(FILE_NAME_PARAM, file_name, MAX_FILE_NEME_DIM);
		parameter[count] = file_option;
		count ++;		
		parameter[count] = file_name;
		count ++;
		
		if(OUTPUT) {
			file_output=calloc(MAX_FILE_OUTPUT_DIM, sizeof(char));
			mxGetString(FILE_OUTPUT_PARAM, file_output, MAX_FILE_OUTPUT_DIM);	
			parameter[count] = output_option;
			count ++;
			parameter[count] = file_output;
			count ++;
		}
		
		if(!VERBOSE) {
			parameter[count] = verbose_option;
			count ++;				
		}
		
		if(PEDANTIC) {
			parameter[count] = pedantic_option;
			count ++;				
		}
		parameter[count] = NULL;
		
		if(execv(mbdyn_name, parameter) == -1) {
			fprintf(stderr,"Cannot start child process\n");
			exit(-1);
		}
	}
	case -1: 
		fprintf(stderr,"Cannot fork\n");
		break;

	default:
		{
		int sleeptime = 1;
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

static void mdlOutputs(SimStruct *S, int_T tid)
{
}


static void mdlTerminate(SimStruct *S)
{
#ifdef MATLAB_MEX_FILE
	int pid = (int_T)ssGetIWork(S)[0];
	printf("Qui!\n");
	printf("Qui2!\n");
	/*kill(pid,SIGINT);*/
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
