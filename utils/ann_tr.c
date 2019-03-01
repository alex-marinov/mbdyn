/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

#include "ann.h"

/* TRAINING PARAMETERS DEFAULT */
static ann_training_mode_t TRAINING_MODE = ANN_TM_BATCH;
static float TOLL = 0.;
static int MAXITER = 1000;
static int PRINTSTEP = 1;
static int SAVESTEP = 1;
static int verbose = 0;
static char *ANNfile = "data/ann.dat";
static char *INPUTfile = "data/Input.dat";
static char *DOUTPUTfile = "data/DOutput.dat";
static char *SAVEfile = "data/ann_tr.dat";
static char *NN_OUTPUTfile = "data/NNOutput.dat";

static void
print_usage(void)
{
	fprintf(stdout, "\nUSAGE OPTIONS:\n"
		"  -u, --usage\n"
		"	print usage\n"		
		"  -v, --verbose\n"
		"	verbose output\n"		
		"  -m, --mode\n"
		"	training mode: BATCH (default)\n"		
		"	               SEQUENTIAL\n"		
		"  -t, --toll\n"
		"	tollerance ( default 0. )\n"		
		"  -i, --maxiter\n"
		"	max number of training iteration (default) 1000\n"		
		"  -p, --print\n"
		"	printing output step (default 1)\n"		
		"  -s, --save\n"
		"	saving ANN trained step (default 1)\n"		
		"  -A, --ann_init\n"
		"	filename of initialized neural network (default data/ann.dat)\n"		
		"  -I, --input\n"
		"	filename of network training input (default data/Input.dat)\n"		
		"  -O, --des_output\n"
		"	filename of network desired output (default data/DOutput.dat)\n"		
		"  -T, --ann_tr\n"
		"	filename where to save trained neural network (default data/ann_tr.dat)\n"		
		"  -N, --nn_output\n"
		"	filename where to save trained neural network output (default data/NNOutput.dat)\n"		
	);
	exit(0);
}

int
main(int argc, char *argv[])
{
        ANN net = { 0 };
	matrix INPUT, DES_OUTPUT, NN_OUTPUT;
	matrix INPUT2, DES_OUTPUT2;
	vector ERR2;
        int N_sample,i,j,ind;
	int Niter, CNT;
	double err1, err2;
	ANN_vector_matrix W1, W2, dWold, dWnew;
	FILE *fh;
	int opt;
	extern char *optarg;


	/* 0. Training options */
	do {
#ifdef HAVE_GETOPT_LONG
		static struct option options[] = {
			{ "usage",	0, 0, 'u'  },
			{ "verbose",	0, 0, 'v'  },	
			{ "mode",	1, 0, 'm'  },
			{ "toll",	1, 0, 't'  },
			{ "maxiter",	1, 0, 'i'  },
			{ "print",	1, 0, 'p'  },
			{ "save",	1, 0, 's'  },
			{ "ann_init",	1, 0, 'A'  },
			{ "input",	1, 0, 'I'  },
			{ "des_output",	1, 0, 'O'  },
			{ "ann_tr",	1, 0, 'T'  },
			{ "nn_output",	1, 0, 'N'  }
		};
		opt = getopt_long(argc, argv, "uvm:t:i:p:s:A:I:O:T:N:", options, NULL);
#else /* ! HAVE_GETOPT_LONG */
		opt = getopt(argc, argv, "uvm:t:i:p:s:A:I:O:T:N:");
#endif /* ! HAVE_GETOPT_LONG */
		switch (opt) {
		case 'u':	print_usage();
				break;
		case 'v':	verbose = 1;
				break;
		case 'm':
			if (strcasecmp(optarg, "batch") == 0) {
				TRAINING_MODE = ANN_TM_BATCH;
			} else if (strcasecmp(optarg, "sequential") == 0) {
				TRAINING_MODE = ANN_TM_SEQUENTIAL;
			} else {
				fprintf(stderr, "unknown training mode \"%s\" {batch|sequential}\n", optarg);
				return 1;
			}
			break;
				break;
		case 't':	TOLL = atof(optarg);
				break;
		case 'i':	MAXITER = atoi(optarg);
				break;
		case 'p':	PRINTSTEP = atoi(optarg);
				break;
		case 's':	SAVESTEP = atoi(optarg);
				break;
		case 'A':	ANNfile = optarg;
				break;
		case 'I':	INPUTfile = optarg;
				break;
		case 'O':	DOUTPUTfile = optarg;
				break;
		case 'T': 	SAVEfile = optarg;
				break;
		case 'N':	NN_OUTPUTfile = optarg;
				break;
		default:	break;
		}
	} while (opt >= 0);


        /* 1. Artificial Neural Network inizialization*/
        printf("LOADING DATA...\n");
        if (ANN_init(&net, ANNfile)) {
                fprintf(stdout, "Error in ANN initialization\n");
                return 1;
        }
	if (TRAINING_MODE == ANN_TM_BATCH) {    // ADAPTIVE LEARNING RATE
		if (ANN_vector_matrix_init(&W1, net.N_neuron, net.N_layer)) {
			fprintf(stderr, "Initialization error\n");
			return 1;
		}
		if (ANN_vector_matrix_init(&W2, net.N_neuron, net.N_layer)) {
			fprintf(stderr, "Initialization error\n");
			return 1;
		}
		if (ANN_vector_matrix_init(&dWold, net.N_neuron, net.N_layer)) {
			fprintf(stderr, "Initialization error\n");
			return 1;
		}
		if (ANN_vector_matrix_init(&dWnew, net.N_neuron, net.N_layer)) {
			fprintf(stderr, "Initialization error\n");
			return 1;
		}
	}

        /* 2. Trainig data data acquisition*/
        N_sample = 0;
        if (ANN_DataRead(&INPUT, &N_sample, INPUTfile)) {
		fprintf(stderr, "Error in Input data acquisition\n");
                return 1;
        }
        if (ANN_DataRead(&DES_OUTPUT, &N_sample, DOUTPUTfile)) {
		fprintf(stderr, "Error in Output data acquisition\n");
                return 1;
        }
	if (matrix_init(&NN_OUTPUT, N_sample, net.N_output)) {
		fprintf(stderr, "Error in NN_output matrix initialization\n");
		return 1;
	}
	if (matrix_init(&DES_OUTPUT2, N_sample, net.N_output)) {
		fprintf(stderr, "Error in NN_output matrix initialization\n");
		return 1;
	}
	if (matrix_init(&INPUT2, N_sample, net.N_input)) {
		fprintf(stderr, "Error in NN_output matrix initialization\n");
		return 1;
	}
	if (vector_init(&ERR2, MAXITER)) {
		fprintf(stderr, "Error in NN_output matrix initialization\n");
		return 1;
	}

        ANN_write(&net, stdout, ANN_W_A_TEXT);
	fprintf(stdout, "TRAINING....\n");

	Niter = 0;
	err2 = 10000000000.;
	CNT = 0;
	do {
		Niter++;
		err1 = err2;
		if (TRAINING_MODE == ANN_TM_BATCH) {
			if (ANN_vector_matrix_ass(&W2, &W1, net.N_neuron, net.N_layer, 1.)) {
				fprintf(stderr, "Error in ....\n");
			}
			if (ANN_vector_matrix_ass(&W1, &net.W, net.N_neuron, net.N_layer, 1.)) {
				fprintf(stderr, "Error in ....\n");
			}
		}
		ANN_reset(&net);
		
                if (ANN_TrainingEpoch(&net, &INPUT, &DES_OUTPUT, &NN_OUTPUT, N_sample, TRAINING_MODE)) {
                        fprintf(stderr, "Error: ANN_TrainingEpoch@main ppp\n");
                        return 1;
                }
		if (verbose) {
		        ANN_write(&net, stdout, ANN_W_A_TEXT);
		}
		ANN_TotalError(&DES_OUTPUT, &NN_OUTPUT, &err2);

		/* per l'addestramento in modalità BATCH il tasso di apprendimento è
		 * adattativo!!! */
		if (TRAINING_MODE == ANN_TM_BATCH) {
			CNT++;
			while (err2 > err1) {
                                CNT = 0;
                                net.eta = 0.5*net.eta;
				if (verbose) {
					fprintf(stdout, "Network's learning rate decreasing (eta = %lf)\n", net.eta);
				}
				ANN_vector_matrix_ass(&net.W, &W2, net.N_neuron, net.N_layer, 1.);
				ANN_vector_matrix_ass(&dWold, &dWold, net.N_neuron, net.N_layer, 0.5);
				ANN_WeightUpdate(&net, dWold, 1.);
				ANN_vector_matrix_ass(&W1, &net.W, net.N_neuron, net.N_layer, 1.);

				ANN_reset(&net);
                		if (ANN_TrainingEpoch(&net, &INPUT, &DES_OUTPUT, &NN_OUTPUT, N_sample, TRAINING_MODE)) {
		                        fprintf(stderr, "Error: ANN_TrainingEpoch@main\n");
                		        return 1;
                		}
				ANN_TotalError(&DES_OUTPUT, &NN_OUTPUT, &err2);
			}
			ANN_vector_matrix_ass(&dWold, &net.dW, net.N_neuron, net.N_layer, 1.);
                        if (CNT == 20) {
                                net.eta = 1.1*net.eta;
                                if (verbose) {
					fprintf(stdout, "Network's learning rate increasing (eta = %lf)\n", net.eta);
				}
                                CNT = 0;
                        }
		}
		/* randimize training example */
		/*matrix_write(&INPUT, stdout, W_M_BIN);
		matrix_write(&DES_OUTPUT, stdout, W_M_BIN);*/
		matrix_copy( &INPUT2, &INPUT, 1. );
		matrix_copy( &DES_OUTPUT2, &DES_OUTPUT, 1. );
		for ( i=0; i<N_sample; i++ ){
			ind = floor(rand()%(N_sample-i));
			for ( j=0; j< INPUT2.Ncolumn; j++){
				INPUT.mat[i][j] = INPUT2.mat[ind][j];
				INPUT2.mat[ind][j] = INPUT2.mat[N_sample-1-i][j];
			}
			for ( j=0; j< DES_OUTPUT2.Ncolumn; j++){

				DES_OUTPUT.mat[i][j] = DES_OUTPUT2.mat[ind][j];
				DES_OUTPUT2.mat[ind][j] = DES_OUTPUT2.mat[N_sample-1-i][j];
			}
		}
		/*matrix_write(&INPUT, stdout, W_M_BIN);
		matrix_write(&DES_OUTPUT, stdout, W_M_BIN);
		getchar();*/

                if (!(Niter%PRINTSTEP)) {
                        fprintf(stdout, "TRAINING:    iter:%d       ", Niter);
                	fprintf(stdout, "Square error: :%le\n", err2);
                }

                if( !(Niter%SAVESTEP) ){
			fh = fopen(SAVEfile, "w");
                        fprintf( stdout, "SAVING DATA...\n");
                        if (ANN_write( &net, fh, ANN_W_A_BIN)) {
                                fprintf(stderr, "Error in data saving\n");
                                return 1;
                        }
			fclose(fh);
                }
		ERR2.vec[Niter-1] = err2;

        } while ((err2>TOLL) && (Niter<MAXITER));


        fprintf(stdout, "SAVING DATA...\n");
	fh = fopen(SAVEfile, "w");
        if (ANN_write(&net, fh, ANN_W_A_BIN)) {
                fprintf(stderr, "Error: ANN_save@main\n");
                return 1;
        }
	fclose(fh);
	ANN_DataWrite(&NN_OUTPUT, NN_OUTPUTfile);
	fh = fopen("ERR.txt", "w");
	vector_write(&ERR2,fh,W_M_BIN );
	fclose(fh);

        /* dynamic memory free*/
	matrix_destroy(&INPUT);
	matrix_destroy(&DES_OUTPUT);
	matrix_destroy(&INPUT2);
	matrix_destroy(&DES_OUTPUT2);
	matrix_destroy(&NN_OUTPUT);
	vector_destroy(&ERR2);
	if (TRAINING_MODE == ANN_TM_BATCH) {
		for (i = 0; i < net.N_layer + 1; i++) {
			matrix_destroy(&W1[i]);
			matrix_destroy(&W2[i]);
			matrix_destroy(&dWnew[i]);
			matrix_destroy(&dWold[i]);
		}
		free(W1);
		free(W2);
		free(dWnew);
		free(dWold);	
	}

        ANN_destroy(&net);
	fprintf(stdout, "END.......\n");
        return 0;
}
