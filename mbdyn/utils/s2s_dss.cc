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

typedef struct FCS_data_struct{
	int N, Nvel;
	double fc;
	char matrixA[MAX_STR_LENGTH], matrixB[MAX_STR_LENGTH], matrixC[MAX_STR_LENGTH], matrixD[MAX_STR_LENGTH], velocity[MAX_STR_LENGTH];
	char MeasuresSocketPath[MAX_STR_LENGTH], ControlsSocketPath[MAX_STR_LENGTH];

} FCS_data_struct;

typedef struct timeval time_eval;

int FCSdataRead( FCS_data_struct *, char * );
int SearchInterval( matrix *, double );
int matrix_interpolation( matrix *, matrix *, double , double , double,  matrix * );

int
main(int argc, char *argv[])
{
	s2s_t	s2s_measures;
	s2s_t	s2s_controls;
	FCS_data_struct FCSdata;
	unsigned  iStepCONCounter, NstepXsec;
	matrix matA, matB, matC, matD, vel;
	matrix *MAT_A, *MAT_B, *MAT_C, *MAT_D;
	matrix Ad, Bd, Cd, Dd;
	unsigned i, j, k;
	unsigned NumInputs, NumOutputs;
	double Velocity;
	int IND;
	vector x_k, u_k, tmp_u1, tmp_u2, tmp_x1, tmp_x2, y_k;
	char FileName[MAX_STR_LENGTH];
	FILE *fh;
	#ifdef COMPUTE_TIME
	FILE *fh_time;
	#endif
	
	if( argc == 2){
		strcpy(	FileName, argv[1]);
	} else {
		fprintf( stderr, "DATA FILE MISSING! (usage: ./s2sctrl DataFile)\n");
		exit(EXIT_FAILURE);
	}
		
	
	if (!FCSdataRead(&FCSdata, FileName)){
		fprintf( stderr, "ERROR IN FCS DATA READING");
		exit(EXIT_FAILURE);
	} 
	
#ifdef VERBOSE
	fprintf( stdout, "\nFCS PROPERTIES:\n");
	fprintf( stdout, "- order: %d\n", FCSdata.N);
	fprintf( stdout, "- number of controllers: %d\n", FCSdata.Nvel);
	fprintf( stdout, "- discretization frequency: %le\n", FCSdata.fc);
	fprintf( stdout, "\nMATRICES:\n");
	fprintf( stdout, "- matrixA path: %s\n", FCSdata.matrixA);
	fprintf( stdout, "- matrixB path: %s\n", FCSdata.matrixB);
	fprintf( stdout, "- matrixC path: %s\n", FCSdata.matrixC);
	fprintf( stdout, "- matrixD path: %s\n", FCSdata.matrixD);
	fprintf( stdout, "- velocity path: %s\n", FCSdata.velocity);
	fprintf( stdout, "\nSOCKETS:\n");
	fprintf( stdout, "- measures socket path: %s\n", FCSdata.MeasuresSocketPath);
	fprintf( stdout, "- controls socket path: %s\n", FCSdata.ControlsSocketPath);
#endif

	iStepCONCounter = 0;
	
	NstepXsec = FCSdata.fc;
	
	/* leggo la parametrizzazione nella velocità di volo */
	matrix_init( &vel, FCSdata.Nvel, 1 );
	fh = fopen(FCSdata.velocity, "r");
	matrix_read( &vel, fh, 1);
	fclose(fh);

	/* leggo le matrici del controllore */
	/* MATRIX A */
	matrix_init( &Ad, FCSdata.N, FCSdata.N );
	matrix_init( &matA, FCSdata.Nvel, FCSdata.N*FCSdata.N );
	fh = fopen(FCSdata.matrixA, "r");
	matrix_read( &matA, fh, 1);
	fclose(fh);
	//matrix_write( &matA, stdout, W_M_BIN);
	/* riordino i dati */
        MAT_A = (matrix *)calloc( FCSdata.Nvel, sizeof(matrix) );
	for ( k=0; k<FCSdata.Nvel; k++){
		matrix_init( &MAT_A[k], FCSdata.N, FCSdata.N);
		for ( i=0; i<FCSdata.N; i++){
			for ( j=0; j<FCSdata.N; j++){
				MAT_A[k].mat[j][i] = matA.mat[k][j+i*FCSdata.N];
			}
		}
		//matrix_write( &MAT_A[k], stdout, W_M_BIN);
	}
	/* MATRIX B */
	matrix_init( &Bd, FCSdata.N, 1 );
	matrix_init( &matB, FCSdata.Nvel, FCSdata.N );
	fh = fopen(FCSdata.matrixB, "r");
	matrix_read( &matB, fh, 1);
	fclose(fh);
	//matrix_write( &matA, stdout, W_M_BIN);
	/* riordino i dati */
        MAT_B = (matrix *)calloc( FCSdata.Nvel, sizeof(matrix) );
	for ( k=0; k<FCSdata.Nvel; k++){
		matrix_init( &MAT_B[k], FCSdata.N, 1);
		for ( i=0; i<1; i++){
			for ( j=0; j<FCSdata.N; j++){
				MAT_B[k].mat[j][i] = matB.mat[k][j+i*FCSdata.N];
			}
		}
		//matrix_write( &MAT_B[k], stdout, W_M_BIN);
	}
	/* MATRIX C */
	matrix_init( &Cd, 1, FCSdata.N );
	matrix_init( &matC, FCSdata.Nvel, FCSdata.N);
	fh = fopen(FCSdata.matrixC, "r");
	matrix_read( &matC, fh, 1);
	fclose(fh);
	//matrix_write( &matA, stdout, W_M_BIN);
	/* riordino i dati */
        MAT_C = (matrix *)calloc( FCSdata.Nvel, sizeof(matrix) );
	for ( k=0; k<FCSdata.Nvel; k++){
		matrix_init( &MAT_C[k], 1, FCSdata.N);
		for ( i=0; i<FCSdata.N; i++){
			for ( j=0; j<1; j++){
				MAT_C[k].mat[j][i] = matC.mat[k][i];
			}
		}
		//matrix_write( &MAT_C[k], stdout, W_M_BIN);
	}
	/* MATRIX D */
	matrix_init( &Dd, 1, 1 );
	matrix_init( &matD, FCSdata.Nvel, 1);
	fh = fopen(FCSdata.matrixD, "r");
	matrix_read( &matD, fh, 1);
	fclose(fh);
	//matrix_write( &matA, stdout, W_M_BIN);
	/* riordino i dati */
        MAT_D = (matrix *)calloc( FCSdata.Nvel, sizeof(matrix) );
	for ( k=0; k<FCSdata.Nvel; k++){
		matrix_init( &MAT_D[k], 1, 1);
		for ( i=0; i<1; i++){
			for ( j=0; j<1; j++){
				MAT_D[k].mat[j][i] = matD.mat[k][i];
			}
		}
		//matrix_write( &MAT_D[k], stdout, W_M_BIN);
	}

	/* working vectors initialization */
	vector_init( &x_k, FCSdata.N);
	vector_init( &u_k, 1);
	vector_init( &y_k, 1);
	vector_init( &tmp_u1, 1);
	vector_init( &tmp_u2, 1);
	vector_init( &tmp_x1, FCSdata.N);
	vector_init( &tmp_x2, FCSdata.N);

	NumInputs = 2;
	NumOutputs = 1;

	/* sockets initialization */
	s2s_measures.nChannels = NumInputs;
	s2s_controls.nChannels = NumOutputs;
	s2s_measures.path = FCSdata.MeasuresSocketPath;
	s2s_controls.path = FCSdata.ControlsSocketPath;
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
		Velocity = s2s_measures.dbuf[0]; 
		y_k.vec[0] = s2s_measures.dbuf[1]; 
		
		if( iStepCONCounter%NstepXsec == 0){
			printf("time: %e s\n", (double)iStepCONCounter/(double)NstepXsec);
		}

		/* matrix interpolation */
		IND = SearchInterval( &vel, Velocity);
		if ( IND == -1){
			matrix_copy(&Ad, &MAT_A[0], 1.);
			matrix_copy(&Bd, &MAT_B[0], 1.);
			matrix_copy(&Cd, &MAT_C[0], 1.);
			matrix_copy(&Dd, &MAT_D[0], 1.);
		}
		if ( IND == -2){
			matrix_copy(&Ad, &MAT_A[vel.Nrow-1], 1.);
			matrix_copy(&Bd, &MAT_B[vel.Nrow-1], 1.);
			matrix_copy(&Cd, &MAT_C[vel.Nrow-1], 1.);
			matrix_copy(&Dd, &MAT_D[vel.Nrow-1], 1.);
		}
		if (IND>=0){
			matrix_interpolation( &MAT_A[IND], &MAT_A[IND+1], vel.mat[IND][0] , vel.mat[IND+1][0] , Velocity,  &Ad );
			matrix_interpolation( &MAT_B[IND], &MAT_B[IND+1], vel.mat[IND][0] , vel.mat[IND+1][0] , Velocity,  &Bd );
			matrix_interpolation( &MAT_C[IND], &MAT_C[IND+1], vel.mat[IND][0] , vel.mat[IND+1][0] , Velocity,  &Cd );
			matrix_interpolation( &MAT_D[IND], &MAT_D[IND+1], vel.mat[IND][0] , vel.mat[IND+1][0] , Velocity,  &Dd );
		}

		/* calcolo la nuova variabile di controllo */
		matrix_vector_prod( &Cd, &x_k, &tmp_u1 );
		matrix_vector_prod( &Dd, &y_k, &tmp_u2 );
		matrix_vector_prod( &Ad, &x_k, &tmp_x1 );
		matrix_vector_prod( &Bd, &y_k, &tmp_x2 );
		vector_sum( &tmp_u1, &tmp_u2, &u_k, 1. );
		vector_sum( &tmp_x1, &tmp_x2, &x_k, 1. );

		// send new controls
		for (int i = 0; i < s2s_controls.nChannels ; i++) {
			s2s_controls.dbuf[i] = u_k.vec[i];
		}
		s2s_controls.send(0);
	
		#ifdef COMPUTE_TIME
		gettimeofday( &t2, 0 );
    		time2 = ( t2.tv_sec + t2.tv_usec*1e-6 );
		delta_time = time2 - time1;
		fprintf(fh_time, "%le\n", delta_time);
		#endif
		/* incremento il contatore dei passi del controllore */
		iStepCONCounter++;
	}

done:	{
	#ifdef COMPUTE_TIME
	fclose(fh_time);
	#endif	
	s2s_measures.shutdown();
	s2s_controls.shutdown();

	vector_destroy( &x_k);
	vector_destroy( &u_k);
	vector_destroy( &y_k);
	vector_destroy( &tmp_u1);
	vector_destroy( &tmp_u2);
	vector_destroy( &tmp_x1);
	vector_destroy( &tmp_x2);
	matrix_destroy( &vel );
	
	matrix_destroy( &matA );	
	matrix_destroy( &matB );	
	matrix_destroy( &matC );	
	matrix_destroy( &matD );	
	matrix_destroy( &Ad );	
	matrix_destroy( &Bd );	
	matrix_destroy( &Cd );	
	matrix_destroy( &Dd );	
	for ( i=0; i<FCSdata.Nvel; i++){
		matrix_destroy( &MAT_A[i] );
		matrix_destroy( &MAT_B[i] );
		matrix_destroy( &MAT_C[i] );
		matrix_destroy( &MAT_D[i] );
	}
	free(MAT_A);
	free(MAT_B);
	free(MAT_C);
	free(MAT_D);
	printf("THE END\n");
	exit(EXIT_SUCCESS);
	}
}

int FCSdataRead( FCS_data_struct *data, char *FileName ){

char line[MAX_STR_LENGTH], name[MAX_STR_LENGTH];
FILE *fh, *fh1;
char *a;

/* DATA FILE */
fh = fopen(FileName,"r");
if (fh == NULL){
	fprintf( stderr, "FILE NOT FOUND\n");
	return 0;
}

/* Read FCS CONTROLLER ORDER n */
a = fgets(line, MAX_STR_LENGTH, fh);
sscanf(line,"%s%d", name, &data->N);
if ( data->N == 0 ){
	fprintf( stderr, "Controller order must be different from zero \nerror: n = %d\n", data->N);
	return 0;
} 
/* Read VELOCITY NUMBER m */
a = fgets(line, MAX_STR_LENGTH, fh);
sscanf(line,"%s%d", name, &data->Nvel);
if ( data->Nvel <= 0 ){
	fprintf( stderr, "Velocity number must be positive \nerror: m = %d\n", data->Nvel);
	return 0;
} 
/* Read CONTROL FREQUENCY fc */
a = fgets(line, MAX_STR_LENGTH, fh);
sscanf(line,"%s%le", name, &data->fc);
if ( data->fc <= 0 ){
	fprintf( stderr, "Control frequency must be positive \nerror: fc = %le\n", data->fc);
	return 0;
}
/* Read MATRIX A FileName */
a = fgets(line, MAX_STR_LENGTH, fh);
sscanf(line,"%s%s", name, data->matrixA);
fh1 = fopen(data->matrixA, "r");
if ( fh1 == NULL ){
	fprintf( stderr, "Noise file not found\n FileName= %s\n", data->matrixA);
	return 0;
}
fclose(fh1);
/* Read MATRIX B FileName */
a = fgets(line, MAX_STR_LENGTH, fh);
sscanf(line,"%s%s", name, data->matrixB);
fh1 = fopen(data->matrixB, "r");
if ( fh1 == NULL ){
	fprintf( stderr, "Noise file not found\n FileName= %s\n", data->matrixB);
	return 0;
}
fclose(fh1);
/* Read MATRIX C FileName */
a = fgets(line, MAX_STR_LENGTH, fh);
sscanf(line,"%s%s", name, data->matrixC);
fh1 = fopen(data->matrixC, "r");
if ( fh1 == NULL ){
	fprintf( stderr, "Noise file not found\n FileName= %s\n", data->matrixC);
	return 0;
}
fclose(fh1);
/* Read MATRIX D FileName */
a = fgets(line, MAX_STR_LENGTH, fh);
sscanf(line,"%s%s", name, data->matrixD);
fh1 = fopen(data->matrixD, "r");
if ( fh1 == NULL ){
	fprintf( stderr, "Noise file not found\n FileName= %s\n", data->matrixD);
	return 0;
}
fclose(fh1);
/* Read Velocity FileName */
a = fgets(line, MAX_STR_LENGTH, fh);
sscanf(line,"%s%s", name, data->velocity);
fh1 = fopen(data->velocity, "r");
if ( fh1 == NULL ){
	fprintf( stderr, "Noise file not found\n FileName= %s\n", data->velocity);
	return 0;
}
fclose(fh1);

/* Read MEASURES SOCKET NAME MeasuresSocketPath */
a = fgets(line, MAX_STR_LENGTH, fh);
sscanf(line,"%s%s", name, data->MeasuresSocketPath);
/* Read CONTROLS SOCKET NAME ControlsSocketPath */
a = fgets(line, MAX_STR_LENGTH, fh);
sscanf(line,"%s%s", name, data->ControlsSocketPath);

fclose(fh);
return 1;
}

int SearchInterval( matrix *Vel, double velocity ){

	int IND;

	IND = -1;
	for (int i = 0; i<Vel->Nrow-1; i++){
		if ( (velocity>=Vel->mat[i][0]) && (velocity<=Vel->mat[i+1][0]) ){
			IND = i;
		}
	}
	if (IND == -1){
		if ( velocity<Vel->mat[0][0] ){
			IND = -1;
		}
		if ( velocity>Vel->mat[Vel->Nrow-1][0] ){
			IND = -2;
		}
	}

	return IND;
}

int matrix_interpolation( matrix *mat1, matrix *mat2, double vel1, double vel2, double Velocity,  matrix *matOUT ){

	unsigned i, j;

	for ( i=0; i<matOUT->Nrow; i++){	
		for ( j=0; j<matOUT->Ncolumn; j++){	
			matOUT->mat[i][j] = mat1->mat[i][j] + ((mat2->mat[i][j]-mat1->mat[i][j])/(vel2-vel1))*(Velocity-vel1);
		}
	}

	return 1;
}
