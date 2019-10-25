/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2010
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
 * Copyright (C) 2010
 *
 * Mattia Mattaboni	<mattaboni@aero.polimi.it>
 */

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/un.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <fcntl.h>
#include <signal.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#include <iostream>
#include <sstream>
#include <vector>

#include "sock.h"
#include "s2s.h"
#include "matrix.h"
#include "GPC.h"

#define MAX_STR_LENGTH 200
#define VERBOSE
#define COMPUTE_TIME

typedef struct GPC_data_struct{
	char IDinputFileName[MAX_STR_LENGTH], NoiseFileName[MAX_STR_LENGTH];
	char MeasuresSocketPath[MAX_STR_LENGTH], ControlsSocketPath[MAX_STR_LENGTH];
	char ComputedControlInputsFileName[MAX_STR_LENGTH], IdentifiedOutputsFileName[MAX_STR_LENGTH];
	char MeasuredOutputsFileName[MAX_STR_LENGTH], ARXParametersFileName[MAX_STR_LENGTH];
	int n, m, p, FlagSimplyProper, s, rhoLength;
	int PreConditioningWindowLength;
	int IDinputLength, NoiseLength;
	int FlagSaveOutputs;
	double mu, delta, rho1, rho2, fc;
	double IdentificationON, IdentificationOFF;
	double ControlON1, ControlOFF1;
	double ControlON2, ControlOFF2;
	double ControlON3, ControlOFF3;
	double ControlON4, ControlOFF4;
	double IDENTIFICATIONINPUT_AMPLITUDE, NOISE_AMPLITUDE;

} GPC_data_struct;

typedef struct timeval time_eval;

int GPCdataRead( GPC_data_struct *, char * );

int
main(int argc, char *argv[])
{
	s2s_t	s2s_measures;
	s2s_t	s2s_controls;
	GPC_data_struct GPCdata;
	ARX_Model ARX_m;
	ARMAX_Model ARMAX_m;
	GPC_Model GPC_m;
	unsigned  iStepCONCounter, NstepXsec;
	matrix IDinput, MeasNoise;
	vector u_k, u_CON_k, Us_ID, y_k;
	double WindowLength;
	vector Mean;
	matrix Y_prev;
	char FileName[MAX_STR_LENGTH];
	FILE *fh, *fh_ComputedControlInputs, *fh_IdentifiedOutputs, *fh_MeasuredOutputs, *fh_ARXParameters;
	#ifdef COMPUTE_TIME
	FILE *fh_time;
	#endif
	
	if( argc == 2){
		strcpy(	FileName, argv[1]);
	} else {
		fprintf( stderr, "DATA FILE MISSING! (usage: ./s2sctrl DataFile)\n");
		exit(EXIT_FAILURE);
	}
		
	
	if (!GPCdataRead(&GPCdata, FileName)){
		fprintf( stderr, "ERROR IN GPC DATA READING");
		exit(EXIT_FAILURE);
	} 
	
#ifdef VERBOSE
	if (GPCdata.n > 0){
		fprintf( stdout, "\nARX MODEL PROPERTIES:\n");
		fprintf( stdout, "- order: %d\n", GPCdata.n);
	} else {
		fprintf( stdout, "\nARMAX MODEL PROPERTIES:\n");
		fprintf( stdout, "- order: %d\n", -GPCdata.n);
	}
	fprintf( stdout, "- inputs number: %d\n", GPCdata.m);
	fprintf( stdout, "- outputs number: %d\n", GPCdata.p);
	if (GPCdata.FlagSimplyProper == 1){
		fprintf( stdout, "- simply proper model\n");
	} else {
		fprintf( stdout, "- strictly proper model\n");
	}
	fprintf( stdout, "- discretization frequency: %le\n", GPCdata.fc);
	fprintf( stdout, "\nMODEL IDENTIFICATION:\n");
	fprintf( stdout, "- forgetting factor: %le\n", GPCdata.mu);
	fprintf( stdout, "- P matrix initialization: P = %leI\n", GPCdata.delta);
	fprintf( stdout, "\nCONTROLLER:\n");
	fprintf( stdout, "- control horizon: %d\n", GPCdata.s);
	fprintf( stdout, "- control penalty function changes linearly from %le to %le \n\tin %d time steps after the controller activation\n", GPCdata.rho1, GPCdata.rho2, GPCdata.rhoLength);
	fprintf( stdout, "\nPRE-CONDITIONING:\n");
	if (GPCdata.PreConditioningWindowLength <= 0){
		fprintf( stdout, "- no measures pre-conditioning\n");
	} else {
		fprintf( stdout, "- measures pre-conditioning: the average value is computed with the last %d time steps\n", GPCdata.PreConditioningWindowLength);
	}
	fprintf( stdout, "\nIDENTIFICATION:\n");
	fprintf( stdout, "- identification is active with the interval from %le s to %le s\n", GPCdata.IdentificationON, GPCdata.IdentificationOFF);
	fprintf( stdout, "\nCONTROL:\n");
	fprintf( stdout, "- control is active with the intervals\n");
	fprintf( stdout, "\t\tfrom %le s to %le s\n", GPCdata.ControlON1, GPCdata.ControlOFF1);
	fprintf( stdout, "\t\tfrom %le s to %le s\n", GPCdata.ControlON2, GPCdata.ControlOFF2);
	fprintf( stdout, "\t\tfrom %le s to %le s\n", GPCdata.ControlON3, GPCdata.ControlOFF3);
	fprintf( stdout, "\t\tfrom %le s to %le s\n", GPCdata.ControlON4, GPCdata.ControlOFF4);
	fprintf( stdout, "\nINPUTS DATA FILES:\n");
	fprintf( stdout, "- identification inputs: %s (length: %d, scale factor:%le)\n", GPCdata.IDinputFileName, GPCdata.IDinputLength, GPCdata.IDENTIFICATIONINPUT_AMPLITUDE);
	fprintf( stdout, "- noise: %s (length: %d, scale factor: %le)\n", GPCdata.NoiseFileName, GPCdata.NoiseLength, GPCdata.NOISE_AMPLITUDE);
	fprintf( stdout, "\nSOCKETS:\n");
	fprintf( stdout, "- measures socket path: %s\n", GPCdata.MeasuresSocketPath);
	fprintf( stdout, "- controls socket path: %s\n", GPCdata.ControlsSocketPath);
	if (GPCdata.FlagSaveOutputs == 0 ){
		fprintf( stdout, "\nOUTPUTS NOT SAVED\n");
	} else {
		fprintf( stdout, "\nOUTPUTS DATA FILES:\n");
		fprintf( stdout, "- computed control inputs: %s\n", GPCdata.ComputedControlInputsFileName);
		fprintf( stdout, "- identified outputs: %s\n", GPCdata.IdentifiedOutputsFileName);
		fprintf( stdout, "- measureded outputs: %s\n", GPCdata.MeasuredOutputsFileName);
		fprintf( stdout, "- ARX parameters: %s\n", GPCdata.ARXParametersFileName);
	}
#endif

	iStepCONCounter = 0;
	
	NstepXsec = GPCdata.fc;
	GPCdata.IdentificationON = GPCdata.IdentificationON*NstepXsec;
	GPCdata.IdentificationOFF = GPCdata.IdentificationOFF*NstepXsec;
	GPCdata.ControlON1 = GPCdata.ControlON1*NstepXsec;
	GPCdata.ControlOFF1 = GPCdata.ControlOFF1*NstepXsec;
	GPCdata.ControlON2 = GPCdata.ControlON2*NstepXsec;
	GPCdata.ControlOFF2 = GPCdata.ControlOFF2*NstepXsec;
	GPCdata.ControlON3 = GPCdata.ControlON3*NstepXsec;
	GPCdata.ControlOFF3 = GPCdata.ControlOFF3*NstepXsec;
	GPCdata.ControlON4 = GPCdata.ControlON4*NstepXsec;
	GPCdata.ControlOFF4 = GPCdata.ControlOFF4*NstepXsec;

	/* identifier and controller initialization*/
	if ( GPCdata.n > 0 ) {
		ARX_Initialize( &ARX_m, GPCdata.n, GPCdata.n, GPCdata.m, GPCdata.p, GPCdata.mu, GPCdata.delta, GPCdata.FlagSimplyProper);
		GPC_Initialize( &GPC_m, GPCdata.n, GPCdata.n, 0 , GPCdata.m, GPCdata.p, GPCdata.s, GPCdata.rho1);
	} else {
		ARMAX_Initialize( &ARMAX_m, -GPCdata.n, -GPCdata.n, -GPCdata.n, GPCdata.m, GPCdata.p, GPCdata.mu, GPCdata.delta, GPCdata.FlagSimplyProper);
		GPC_Initialize( &GPC_m, -GPCdata.n, -GPCdata.n, -GPCdata.n, GPCdata.m, GPCdata.p, GPCdata.s, GPCdata.rho1);
	}

	/* THETA MATRIX initiailization */
	//fh = fopen("thetaGPC_175C.txt","r");
	//matrix_read( &ARX_m.theta, fh, 1);
	//fprintf(stdout,"!!!!!! MATRIX THETA INITIALIZATION !!!!!!\n\n");
	//fclose(fh);
	

	/* Pre-conditioning initialization  */
	WindowLength = GPCdata.PreConditioningWindowLength;
	vector_init( &Mean, GPCdata.p);
	if (WindowLength > 0){
		matrix_init( &Y_prev, unsigned(WindowLength), GPCdata.p);
	}
	
	/* working vectors initialization */
	vector_init( &u_k, GPCdata.m);
	vector_init( &y_k, GPCdata.p);
	vector_init( &u_CON_k, GPCdata.m);
	vector_init( &Us_ID, GPCdata.m*(GPCdata.s));

	/* Inputs files reading */
	fh = fopen(GPCdata.IDinputFileName,"r");
	matrix_init( &IDinput, GPCdata.IDinputLength, GPCdata.m );
	matrix_read( &IDinput, fh, 1);
	fclose(fh);
	fh = fopen(GPCdata.NoiseFileName,"r");
	matrix_init( &MeasNoise, GPCdata.NoiseLength, GPCdata.p );
	matrix_read( &MeasNoise, fh, 1);
	fclose(fh);

	/* sockets initialization */
	s2s_measures.nChannels = GPCdata.p;
	s2s_controls.nChannels = GPCdata.m;
	s2s_measures.path = GPCdata.MeasuresSocketPath;
	s2s_controls.path = GPCdata.ControlsSocketPath;
	s2s_controls.create = 0;
	s2s_measures.create = 0;
	try {
		s2s_measures.prepare();
		s2s_controls.prepare();
	} catch (...) {
		s2s_measures.shutdown();
		s2s_controls.shutdown();
		exit(EXIT_FAILURE);
	}
	s2s_measures.dbuf.resize(s2s_measures.nChannels);
	s2s_controls.dbuf.resize(s2s_controls.nChannels);

	/* OUTPUT FILES */	
	if (GPCdata.FlagSaveOutputs == 1 ){
		fh_ComputedControlInputs = fopen( GPCdata.ComputedControlInputsFileName, "w");
		if (fh_ComputedControlInputs == NULL){
			fprintf( stderr, "FILE NOT FOUND %s\n", GPCdata.ComputedControlInputsFileName);
			exit(EXIT_FAILURE);
		}
		fh_IdentifiedOutputs = fopen( GPCdata.IdentifiedOutputsFileName, "w");
		if (fh_IdentifiedOutputs == NULL){
			fprintf( stderr, "FILE NOT FOUND %s\n", GPCdata.IdentifiedOutputsFileName);
			exit(EXIT_FAILURE);
		}
		fh_MeasuredOutputs = fopen( GPCdata.MeasuredOutputsFileName, "w");
		if (fh_MeasuredOutputs == NULL){
			fprintf( stderr, "FILE NOT FOUND %s\n", GPCdata.MeasuredOutputsFileName);
			exit(EXIT_FAILURE);
		}
		fh_ARXParameters = fopen( GPCdata.ARXParametersFileName, "w");
		if (fh_ARXParameters == NULL){
			fprintf( stderr, "FILE NOT FOUND %s\n", GPCdata.ARXParametersFileName);
			exit(EXIT_FAILURE);
		}
	}	
	#ifdef COMPUTE_TIME
	fh_time = fopen("ComputationalTime.txt","w");
	if (fh_time == NULL){
		fprintf( stderr, "FILE NOT FOUND (ComputationalTime.txt)\n" );
		exit(EXIT_FAILURE);
	}
	double time1, time2, delta_time;
	time_eval t1;
	time_eval t2;
	#endif

	while (true) {
		// read new measures
		int len = s2s_measures.recv(0);

		// check sanity
		switch (len) {
		case -1: {
			int		save_errno = errno;
			const char	*err_msg = strerror(save_errno);
			
			silent_cerr("recv(" << s2s_measures.sock << ",\"" << s2s_measures.buf << "\") "
				"failed (" << save_errno << ": " << err_msg << ")"
				<< std::endl);
		}
		case 0:{
			goto done;}
		default:{
			break;}
		}
		
		/* new measures */
		for( int i=0; i< s2s_measures.nChannels; i++){
			y_k.vec[i] = (s2s_measures.dbuf[i]);
		}
		/* sommo all'ingresso di controllo u_CON_k l'identification input IDinput */	
		for( int i= 0; i < s2s_controls.nChannels ; i++){
			if( iStepCONCounter < GPCdata.IdentificationOFF){
			//if( iStepCONCounter < 8.*NstepXsec ){
				u_k.vec[i] = u_CON_k.vec[i] + GPCdata.IDENTIFICATIONINPUT_AMPLITUDE*IDinput.mat[iStepCONCounter][i];
			} else { 
				u_k.vec[i] = u_CON_k.vec[i];
			}
		}
		//u_k.vec[0] = 0.;
		// send new controls
		for (int i = 0; i < s2s_controls.nChannels ; i++) {
			s2s_controls.dbuf[i] = u_k.vec[i];
		}
		s2s_controls.send(0);

		/* add measurement noise */
		for( int i=0; i< s2s_measures.nChannels; i++){
			y_k.vec[i] += GPCdata.NOISE_AMPLITUDE*MeasNoise.mat[iStepCONCounter][i];
		}
		//printf("AA %e\n",ARX_m.theta.mat[0][0]);
#ifdef VERBOSE
		if( iStepCONCounter%NstepXsec == 0){
			printf("time: %e s\n", (double)iStepCONCounter/(double)NstepXsec);
		}
#endif
		/* MEASUREMENTS PRE-CONDITIONING */
		if (WindowLength > 0){
			/* previous WindowLength measurements */
			for( int i=(unsigned)WindowLength-1; i>0; i--){
				for( int j=0; j<GPCdata.p; j++){
					Y_prev.mat[i][j] = Y_prev.mat[i-1][j];
				}
			}
			for( int i=0; i<GPCdata.p; i++){
				Y_prev.mat[0][i] = y_k.vec[i];
			}
			/* measurement average value on the last WindowLength time steps */
			if( iStepCONCounter < (unsigned)WindowLength ){
				/* intial start-up algorith*/
				for( int i=0; i<GPCdata.p; i++){
					Mean.vec[i] = Mean.vec[i]*((double)iStepCONCounter/((double)iStepCONCounter+1.)) + y_k.vec[i]/((double)iStepCONCounter+1.);
				}
			} else {
				/* regime algorith */
				for( int i=0; i<GPCdata.p; i++){
					Mean.vec[i] = (y_k.vec[i]/WindowLength) + Mean.vec[i] - (Y_prev.mat[unsigned(WindowLength-1)][i]/WindowLength);
				}
			}
		}
			
		/* calcolo l'ingresso di controllo al passo K+1 
		   a partire dalla misura dell'uscita al passo K 
		*/
		/* uscita misurata al passo K tolto il suo valor medio */
		vector_sum( &y_k, &Mean, &y_k, -1. );
		/* aggiorno il modello ARX */
		if ((iStepCONCounter >= GPCdata.IdentificationON) && (iStepCONCounter <= GPCdata.IdentificationOFF))  {
			if ( GPCdata.n > 0 ) {
				ARX_RLS( &ARX_m, &y_k);
			} else {
				ARMAX_ELS( &ARMAX_m, &y_k);
				//ARMAX_RML( &ARMAX_m, &y_k);
			}
		}
		/* aggiorno i vettori con Ym(K) e U(K) che ho alcolato al passo precedente */
		if ( GPCdata.n > 0 ) {
			ARX_UpdatePhi( &ARX_m, &y_k, &u_k);
		} else {
			ARMAX_UpdatePhi( &ARMAX_m, &y_k, &u_k, &ARMAX_m.eps);
		}
		/* aggiorno i vettori con U(K) e Ym(K)*/
		if ( GPCdata.n > 0 ) {
			GPC_VectorUpdate( &GPC_m, &y_k, &u_k, &u_k);
		} else {
			GPC_VectorUpdate( &GPC_m, &y_k, &u_k, &ARMAX_m.eps);
		}
		#ifdef COMPUTE_TIME
		gettimeofday( &t1, 0 );
    		time1 = ( t1.tv_sec + t1.tv_usec*1e-6 );
		#endif
		vector_null(&u_CON_k);
		if (  ((iStepCONCounter >= GPCdata.ControlON1) && (iStepCONCounter <= GPCdata.ControlOFF1)) ||
					((iStepCONCounter >= GPCdata.ControlON2) && (iStepCONCounter <= GPCdata.ControlOFF2)) ||
					((iStepCONCounter >= GPCdata.ControlON3) && (iStepCONCounter <= GPCdata.ControlOFF3)) ||
					((iStepCONCounter >= GPCdata.ControlON4) && (iStepCONCounter <= GPCdata.ControlOFF4)) ){
			/* assemblo le matrici dell'equazione di predizione delle uscite */
			if ( GPCdata.n > 0 ) {
				GPC_PredictionFunction( &GPC_m, &ARX_m.theta, ARX_m.SimplyProper, 0 );
			} else {
				GPC_PredictionFunction( &GPC_m, &ARMAX_m.theta, ARMAX_m.SimplyProper, 0 );
			}
			/* calcolo le variabili di controllo */
			/* A). lambda variabile per evitare una brusca applicazione del controllo */
			if (iStepCONCounter < GPCdata.ControlON1+GPCdata.rhoLength){
				GPC_m.lambda = GPCdata.rho1 + ((GPCdata.rho2-GPCdata.rho1)/GPCdata.rhoLength)*(iStepCONCounter-GPCdata.ControlON1); 
			} else {
				GPC_m.lambda = GPCdata.rho2;
			}
			/* B). considero il contributo dei futuri ingressi di identificazione */
			for( int i=0; i<(GPCdata.s); i++ ){
				for( int j=0; j<GPCdata.m; j++ ){
					if( iStepCONCounter < GPCdata.IdentificationOFF){
					//if( iStepCONCounter < 8.*NstepXsec ){
						Us_ID.vec[(GPCdata.s-i-1)*GPCdata.m+j] = GPCdata.IDENTIFICATIONINPUT_AMPLITUDE*IDinput.mat[iStepCONCounter+i+1][j];
					} else {
						Us_ID.vec[(GPCdata.s-i-1)*GPCdata.m+j] = 0.;
					}
				}
			}
			/* calcolo le future s variabili di controllo */
			//GPC_Control( &GPC_m, &Us_ID );
			GPC_ControlW( &GPC_m, &Us_ID );
			/* REciding horizon */
			for( int i= 0; i < GPCdata.m ; i++){
				u_CON_k.vec[i] = GPC_m.Us.vec[GPC_m.m*(GPC_m.s-1)+i];
			}
		}
	
		#ifdef COMPUTE_TIME
		gettimeofday( &t2, 0 );
    		time2 = ( t2.tv_sec + t2.tv_usec*1e-6 );
		delta_time = time2 - time1;
		fprintf(fh_time, "%le\n", delta_time);
		#endif
		/* incremento il contatore dei passi del controllore */
		iStepCONCounter++;
		if (GPCdata.FlagSaveOutputs == 1 ){
			vector_write( &u_CON_k, fh_ComputedControlInputs, W_M_BIN_ROW);
			vector_write( &y_k, fh_MeasuredOutputs, W_M_BIN_ROW);
			if ( GPCdata.n > 0 ){
				vector_write( &ARX_m.yp, fh_IdentifiedOutputs, W_M_BIN_ROW);
				matrix_write( &ARX_m.theta, fh_ARXParameters, W_M_BIN);
			} else {
				vector_write( &ARMAX_m.yp, fh_IdentifiedOutputs, W_M_BIN_ROW);
				matrix_write( &ARMAX_m.theta, fh_ARXParameters, W_M_BIN);
			}
		}
	}

done:	{
	fh = fopen("thetaGPC.txt","w");
	if ( GPCdata.n > 0 ){
		matrix_write( &ARX_m.theta, stdout, W_M_BIN);
		matrix_write( &ARX_m.theta, fh, W_M_BIN);
	} else {
		matrix_write( &ARMAX_m.theta, stdout, W_M_BIN);
		matrix_write( &ARX_m.theta, fh, W_M_BIN);
	}
	fclose(fh);
	#ifdef COMPUTE_TIME
	fclose(fh_time);
	#endif	
	if (GPCdata.FlagSaveOutputs == 1 ){
		fclose( fh_ComputedControlInputs);
		fclose( fh_IdentifiedOutputs);
		fclose( fh_MeasuredOutputs);
		fclose( fh_ARXParameters);
	}	
	s2s_measures.shutdown();
	s2s_controls.shutdown();
	//fh = fopen("C2.txt","w");
	//matrix_write(&GPC_m.C2, fh, W_M_BIN);
	//fclose(fh);
	if ( GPCdata.n > 0 ){
		ARX_Destroy( &ARX_m);
	} else {
		ARMAX_Destroy( &ARMAX_m);
	}
	GPC_Destroy( &GPC_m);
	vector_destroy( &u_k);
	vector_destroy( &y_k);
	vector_destroy( &u_CON_k);
	vector_destroy( &Us_ID);
	if (WindowLength > 0 ){
		matrix_destroy( &Y_prev);
	}
	vector_destroy( &Mean);
	matrix_destroy( &IDinput);
	matrix_destroy( &MeasNoise);
	printf("THE END\n");
	exit(EXIT_SUCCESS);
	}

}

void check_fgets(char*a){
	if (a == NULL){
        	fprintf( stderr, "ERROR READING FROM FILE\n");
	        exit(EXIT_FAILURE);
	}
	return;
}

int GPCdataRead( GPC_data_struct *data, char *FileName ){

char line[MAX_STR_LENGTH], name[MAX_STR_LENGTH];
FILE *fh, *fh1;
char *a;

/* DATA FILE */
fh = fopen(FileName,"r");
if (fh == NULL){
	fprintf( stderr, "FILE NOT FOUND\n");
	return 0;
}

/* Read MODEL ORDER n */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%d", name, &data->n);
if ( data->n == 0 ){
	fprintf( stderr, "Model order must be different from zero \nerror: n = %d\n", data->n);
	return 0;
} 
/* Read SYSTEM INPUTS NUMBER m */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%d", name, &data->m);
if ( data->m <= 0 ){
	fprintf( stderr, "Inputs number must be positive \nerror: m = %d\n", data->m);
	return 0;
} 
/* Read SYSTEM OUTPUTS NUMBER p */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%d", name, &data->p);
if ( data->p <= 0 ){
	fprintf( stderr, "Outputs number must be positive \nerror: p = %d\n", data->p);
	return 0;
} 
/* Read FLAG SIMPLY PROPER  FlagSimplyProper*/
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%d", name, &data->FlagSimplyProper);
if ( (data->FlagSimplyProper != 0) && (data->FlagSimplyProper != 1) ){
	fprintf( stderr, "FlagSimplyProper must be equal to 1 (simply proper model) or 0 (strictly proper model) \nerror: FlagSimplyProper = %d\n", data->FlagSimplyProper);
	return 0;
} 
/* Read CONTROL FREQUENCY fc */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%le", name, &data->fc);
if ( data->fc <= 0 ){
	fprintf( stderr, "Control frequency must be positive \nerror: fc = %le\n", data->fc);
	return 0;
}
/* Read RLS FORGETTING FACTOR mu */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%le", name, &data->mu);
if ( (data->mu <= 0.) || (data->mu > 1.) ){
	fprintf( stderr, "RLS forgetting factor must be 0<mu<=1 \nerror: mu = %le\n", data->mu);
	return 0;
} 
/* Read DELTA delta */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%le", name, &data->delta);
if ( data->delta <= 0 ){
	fprintf( stderr, "Delta number must be positive \nerror: delta = %le\n", data->delta);
	return 0;
}
/* Read LENGTH OF THE WINDOW USED TO COMPUTE THE MEASURES AVERAGE PreConditioningWindowLength*/
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%d", name, &data->PreConditioningWindowLength);
if ( data->PreConditioningWindowLength <= 0 ){
	fprintf( stderr, "WARNING: No measures pre-conditioning.");
} 
/* Read CONTROL HORIZON LENGTH s */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%d", name, &data->s);
if ( data->s <= 0 ){
	fprintf( stderr, "Control horizon must be positive \nerror: s = %d\n", data->s);
	return 0;
} 
/* Read INITIAL CONTROL PENALTY FUNCTION rho1 */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%le", name, &data->rho1);
if ( data->rho1 <= 0 ){
	fprintf( stderr, "Control penalty function must be positive \nerror: rho1 = %le\n", data->rho1);
	return 0;
}
/* Read FINAL CONTROL PENALTY FUNCTION rho2 */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%le", name, &data->rho2);
if ( data->rho2 <= 0 ){
	fprintf( stderr, "Control penalty function must be positive \nerror: rho2 = %le\n", data->rho2);
	return 0;
}
/* Read LENGTH OF CONTROL PENALTY FUNCTION VARIATION INTERVAL rhoLength */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%d", name, &data->rhoLength);
if ( data->rhoLength <= 0 ){
	fprintf( stderr, "The length of the interval within the penalty function changes must be positive \nerror: rhoLength = %d\n", data->rhoLength);
	return 0;
}
/* Read IDENTIFICATION ON IdentificationON */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%le", name, &data->IdentificationON);
/* Read IDENTIFICATION OFF IdentificationOFF */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%le", name, &data->IdentificationOFF);
/* Read CONTROL ON ControlON1 */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%le", name, &data->ControlON1);
/* Read CONTROL OFF ControlOFF1 */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%le", name, &data->ControlOFF1);
/* Read CONTROL ON ControlON2 */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%le", name, &data->ControlON2);
/* Read CONTROL OFF ControlOFF2 */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%le", name, &data->ControlOFF2);
/* Read CONTROL ON ControlON3 */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%le", name, &data->ControlON3);
/* Read CONTROL OFF ControlOFF3 */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%le", name, &data->ControlOFF3);
/* Read CONTROL ON ControlON4 */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%le", name, &data->ControlON4);
/* Read CONTROL OFF ControlOFF4 */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%le", name, &data->ControlOFF4);
/* Read IDENTIFICATION INPUT FILE NAME IDinputFileName */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%s", name, data->IDinputFileName);
fh1 = fopen(data->IDinputFileName, "r");
if ( fh1 == NULL ){
	fprintf( stderr, "Identification input file not found\n FileName= %s\n", data->IDinputFileName);
	return 0;
}
fclose(fh1);
/* Read IDENTIFICATION INPUT FILE LENGTH IDinputLength */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%d", name, &data->IDinputLength);
if ( data->IDinputLength <= 0 ){
	fprintf( stderr, "The length of the identification input file must be positive \nerror: IDinputLength = %d\n", data->IDinputLength);
	return 0;
}
/* Read IDENTIFICATION INPUT SCALE FACTOR  IDENTIFICATIONINPUT_AMPLITUDE */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%le", name, &data->IDENTIFICATIONINPUT_AMPLITUDE);
/* Read MEASUREMENT NOISE FILE NAME NoiseFileName */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%s", name, data->NoiseFileName);
fh1 = fopen(data->NoiseFileName, "r");
if ( fh1 == NULL ){
	fprintf( stderr, "Noise file not found\n FileName= %s\n", data->NoiseFileName);
	return 0;
}
fclose(fh1);
/* Read NOISE FILE LENGTH NoiseLength */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%d", name, &data->NoiseLength);
if ( data->NoiseLength <= 0 ){
	fprintf( stderr, "The length of the noise file must be positive \nerror: NoiseLength = %d\n", data->NoiseLength);
	return 0;
}
/* Read NOISE SCALE FACTOR  NOISE_AMPLITUDE */
a = fgets(line, MAX_STR_LENGTH, fh);
sscanf(line,"%s%le", name, &data->NOISE_AMPLITUDE);
/* Read MEASURES SOCKET NAME MeasuresSocketPath */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%s", name, data->MeasuresSocketPath);
/* Read CONTROLS SOCKET NAME ControlsSocketPath */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%s", name, data->ControlsSocketPath);
/* Read FLAG SAVe OUTPUTS  FlagSaveOutputs*/
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%d", name, &data->FlagSaveOutputs);
if ( (data->FlagSaveOutputs != 0) && (data->FlagSaveOutputs != 1) ){
	fprintf( stderr, "FlagSaveOutputs must be equal to 1 (save) or 0 (do not save) \nerror: FlagSaveOutputs = %d\n", data->FlagSaveOutputs);
	return 0;
} 
/* Read COMPUTED CONTROL INPUTS FILE NAME ComputedControlInputsFileName */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%s", name, data->ComputedControlInputsFileName);
/* Read IDENTIFIED OUTPUTS FILE NAME IdentifiedOutputsFileName */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%s", name, data->IdentifiedOutputsFileName);
/* Read MEASURED OUTPUTS FILE NAME MeasuredOutputsFileName */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%s", name, data->MeasuredOutputsFileName);
/* Read ARX PARAMETERS FILE NAME ARXParametersFileName */
a = fgets(line, MAX_STR_LENGTH, fh); check_fgets(a);
sscanf(line,"%s%s", name, data->ARXParametersFileName);

fclose(fh);
return 1;
}
 
