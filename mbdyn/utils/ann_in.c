/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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
 * Copyright (C) 2008
 *
 * Mattia Mattaboni	<mattaboni@aero.polimi.it>
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <time.h>

#include "ann.h"

struct option options[] = {
        { "usage",      0, 0, 'u'  },
        { "mode",       1, 0, 'm'  },
        { "mean",    	1, 0, 'M'  },
        { "variance",  	1, 0, 'V'  },
        { "min",    	1, 0, 'N'  },
        { "max",  	1, 0, 'X'  },
        { "minW",    	1, 0, 'w'  },
        { "maxW",  	1, 0, 'W'  },
        { "file",    	1, 0, 'f'  }, 
        { "fileINPUT",	1, 0, 'i'  },
        { "fileOUTPUT",	1, 0, 'o'  } 
};
void print_usage( void ){

        fprintf( stdout, "\nUSAGE OPTIONS:\n"
                "  -u, --usage\n"
                "       print usage\n"
                "  -f, --file\n"
                "       ann initialized file name (default data/ann.dat)\n"
                "  -i, --fileINPUT\n"
                "       input file name (default data/Input.dat)\n"
                "  -o, --fileOUTPUT\n"
                "       desired output file name (defalut/DOutput.dat)\n"
                "  -m, --mode\n"
                "       scaling input/output mode: 1 mean-variance  	(default)\n"
                "                                  2 min-max\n" 
                "                                  0 no scaling\n" 
                "  -M, --mean\n"
                "       mean value (defaul 0.)\n"
                "  -V, --variance"
                "       variance value (defaul 1.)\n"
                "  -N, --min\n"
                "       minimum value (defaul -1.)\n"
                "  -X, --max\n"
                "       maximum value (defaul 1.)\n"
                "  -w, --minW\n"
                "       random initialization minimum weight value (defaul -0.001)\n"
                "  -W, --maxW\n"
                "       random initialization maximum weight value (defaul 0.001)\n"
         );
        exit(0);

}


int MODE = 1;
float MEAN = 0.;
float VARIANCE = 1.;
float MIN = -1.;
float MAX = 1.;
float minW = -0.001;
float maxW = 0.001;
static char *fileINPUT = "data/Input.dat";
static char *fileOUTPUT = "data/DOutput.dat";
static char *file = "data/ann.dat";

int main( int argc, char **argv ){

	int opt;
        extern char *optarg;
	matrix INPUT, OUTPUT, SF_INPUT, SF_OUTPUT , MAT;
	double mean, var, min, max;
	int N_input, N_output, N_layer, r, N_sample;
	int *N_neuron = NULL;
	double doub;
	int in;
	unsigned i;
	FILE * fh;

        /* 0. Training options */
        do {
		opt = getopt_long( argc, argv, "uf:i:o:m:M:V:N:X:w:W:", options, NULL  );
                switch( opt ) {
                case 'u':       print_usage();
                                break;   
		case 'f':       file = strdup( optarg );
                                break;
		case 'i':       fileINPUT = strdup( optarg );
                                break;
		case 'o':       fileOUTPUT = strdup( optarg );
                                break;
                case 'm':       MODE = atoi( optarg );
                                break;
                case 'M':       MEAN = atof( optarg );
                                break;
                case 'V':       VARIANCE = atof( optarg );
                                break;
                case 'N':       MIN = atof( optarg );
                                break;
                case 'X':       MAX = atof( optarg );
                                break;
                case 'w':       minW = atof( optarg );
                                break;
                case 'W':       maxW = atof( optarg );
                                break;
                default:        break;
                }
        } while (opt >= 0);

	/* initialize random seed: */
  	srand ( time(NULL) );

	/* inizailizzazione del file contenente le informazioni sulla rete neurale */	
	fprintf( stdout, "ARTIFICIAL NEURAL NETWORK INITIALIZATION....\n" );
	fprintf( stdout, "Output save in %s\n\n", file );
	fh = fopen( file, "w" );

	fprintf( stdout, "*** Network architecture ***\n" );
	fprintf( stdout, "\tNetwork inputs' number?\n" );
	fscanf( stdin, "%d", &N_input );
	if( N_input <= 0 ){
                fprintf( stderr, "ERROR: Input number must be greater than zero.\n");
                goto exit_error;
        }
	fprintf( fh, "%d\n", N_input );
	
	fprintf( stdout, "\tNetwork outputs' number?\n" );
	fscanf( stdin, "%d", &N_output );
	if( N_output <= 0 ){
                fprintf( stderr, "ERROR: Output number must be greater than zero.\n");
                goto exit_error;
        }
	fprintf( fh, "%d\n", N_output );

        fprintf( stdout, "\tHidden layers' number?\n" );
        fscanf( stdin, "%d", &N_layer );
        if( N_layer < 0 ){
		fprintf( stderr, "ERROR: Hidden layer number must be not negative.\n");
                goto exit_error; 
        }       
        fprintf( fh, "%d\n", N_layer );

        fprintf( stdout, "\tTime step delay's number?\n" );
        fscanf( stdin, "%d", &r );
        if( r < 0 ){
		fprintf( stderr, "ERROR Timestep delay number must be not negative.\n");
                goto exit_error; 
        }       

        fprintf( fh, "%d\n", r );
	N_neuron = (int *)malloc( (N_layer+2) * sizeof(int) );
	if (N_neuron==NULL) {
		return 1;
        }
        N_neuron[0] = 0;

        for( i=0; i<N_layer; i++ ){
        	fprintf( stdout, "\tNeuron number at %d hidden layer number?\n", i+1 );
        	fscanf( stdin, "%d", &in );
		N_neuron[i+1] = in;
                if( in <= 0 ){
                        fprintf( stderr, "ERROR: Neuron number at %d hidden layer must be greater than zero.\n",i+1);
                        goto exit_error;
                }
        	fprintf( fh, "%d ", in );
        }
        fprintf( stdout, "\tNeuron number at the visible layer number?\n" );
       	fscanf( stdin, "%d", &in );
	N_neuron[N_layer+1] = in;
        if( in < N_output ){
                fprintf( stderr, "ERROR: Neuron number at the visible layer must be not lesser than output number.\n");
                goto exit_error;
        }
        fprintf( fh, "%d ", in );
	N_neuron[0] = N_input + ( r * N_neuron[N_layer+1] );
	
        fprintf( fh, " \n\n\n" );

	fprintf( stdout, "*** Activation function definition***\n" );
	fprintf( stdout, " 1: Hyperbolic tangent [ y = a*tanh(b*x) ]\n" );
	fprintf( stdout, " 2: Linear function [ y = mx +q ]\n" );
	fprintf( stdout, "\tChoose activation function [ 1 2 ]\n" );
	fscanf( stdin, "%d", &in );
	fprintf( fh, "%d\n", in );
        fprintf( stdout, "\tActivation function parameters?\n" );
        fscanf( stdin, "%lf", &doub );
        fprintf( fh, "\n%e ", doub );
        fscanf( stdin, "%lf", &doub );
        fprintf( fh, " %e\n", doub );
        
	fprintf( fh, " \n\n\n" );

	fprintf( stdout, "*** Training parameters***\n" );
        fprintf( stdout, "\tLearning rate?\n" );
        fscanf( stdin, "%lf", &doub );
        if( doub < 0 ){
		fprintf( stderr, "ERROR: Learning rate must be not negative.\n");
                goto exit_error; 
        }       
        fprintf( fh, "%e\n", doub );
	fprintf( stdout, "\nMomentum term?\n" );
        fscanf( stdin, "%lf", &doub );
        fprintf( fh, "%e\n", doub );

        fprintf( fh, " \n\n\n" );

	fprintf( stdout, "*** Random weight matrices initialization ( Wmin = %f; Wmax = %f )***\n", minW, maxW );
	for( i=0; i<N_layer+1; i++){
		matrix_init( &MAT, N_neuron[i], N_neuron[i+1] );
		matrix_random( &MAT, minW, maxW );
		matrix_write( &MAT, fh, W_M_BIN );
		matrix_destroy( &MAT );
	}

	fprintf( stdout, "*** Input/Output Scale factor computing ***\n" );
        N_sample = 0;
        if( ANN_DataRead( &INPUT, &N_sample, fileINPUT ) ){
                fprintf( stderr, "ERRRO: Error in Input data acquisition\n");
                goto exit_error; 
        }	
        if( ANN_DataRead( &OUTPUT, &N_sample, fileOUTPUT ) ){
                fprintf( stderr, "ERROR: Error in Input data acquisition\n");
                goto exit_error; 
        }	
	matrix_init( &SF_INPUT, INPUT.Ncolumn, 2 );
	matrix_init( &SF_OUTPUT, OUTPUT.Ncolumn, 2 );
	for( i=0; i<N_input; i++ ){
		SF_INPUT.mat[i][0] = 1.0;
	}
	for( i=0; i<N_output; i++ ){
		SF_OUTPUT.mat[i][0] = 1.0;
	}
	if( MODE == 1 ){
		for( i=0; i<INPUT.Ncolumn; i++ ){
			mean = mean_value( &INPUT, i );
			var = variance( &INPUT, i );
			fprintf( stdout, "Input number %d (mean value = %e, variance = %e), do you want ot scale it? [1 = yes, 0 = no]\n", i+1, mean, var );
			fscanf( stdin, "%d", &in);
			if( in == 1){
				SF_INPUT.mat[i][0] = sqrt( VARIANCE/var );
				SF_INPUT.mat[i][1] = MEAN - SF_INPUT.mat[i][0]*mean;
			}	
		}
		for( i=0; i<OUTPUT.Ncolumn; i++ ){
			mean = mean_value( &OUTPUT, i );
			var = variance( &OUTPUT, i );
			SF_OUTPUT.mat[i][0] = sqrt( VARIANCE/var );
			SF_OUTPUT.mat[i][1] = MEAN - SF_OUTPUT.mat[i][0]*mean;
		}
	}
	if( MODE == 2 ){
		for( i=0; i<INPUT.Ncolumn; i++ ){
			min = minimum( &INPUT, i );
			max = maximum( &INPUT, i );
			fprintf( stdout, "Input number %d (min value = %e, max = %e), do you want ot scale it? [1 = yes, 0 = no]\n", i+1, min, max );
			fscanf( stdin, "%d", &in);
			if( in == 1){
				SF_INPUT.mat[i][0] = (MAX-MIN)/(max-min);
				SF_INPUT.mat[i][1] = (max*MIN-MAX*min)/(max-min);
			}
		}
		for( i=0; i<OUTPUT.Ncolumn; i++ ){
			min = minimum( &OUTPUT, i );
			max = maximum( &OUTPUT, i );
			SF_OUTPUT.mat[i][0] = (MAX-MIN)/(max-min);
			SF_OUTPUT.mat[i][1] = (max*MIN-MAX*min)/(max-min);
		}
	}

	matrix_write( &SF_INPUT, fh, W_M_BIN );
	matrix_write( &SF_OUTPUT, fh, W_M_BIN );
	
	exit_error:  	
			matrix_destroy(&INPUT);
			matrix_destroy(&OUTPUT);
			matrix_destroy(&SF_INPUT);
			matrix_destroy(&SF_OUTPUT);
			free(N_neuron);


	fclose(fh);
	return 0;

}
