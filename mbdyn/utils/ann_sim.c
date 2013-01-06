/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

static char *ANNfile = "data/ann.dat";
static char *INPUTfile = "data/Input.dat";
static char *NN_OUTPUTfile = "data/NNOutput.dat";

void
print_usage(void)
{
        fprintf(stdout, "\nUSAGE OPTIONS:\n"
                "  -u, --usage\n"
                "       print usage\n"
                "  -A, --ann\n"
                "       filename of initialized neural network (default data/ann.dat)\n"
                "  -I, --input\n"
                "       filename of network training input (default data/Input.dat)\n"
                "  -N, --nn_output\n"
                "       filename where to save trained neural network output (default data/NNOutput.dat)\n"
        );
        exit(0);

}

int
main(int argc, char *argv[])
{
        ANN net;
	matrix INPUT, NN_OUTPUT;
        int i, j, N_sample;
        int opt;
        extern char *optarg;


        /* 0. Training options */
        do {
#ifdef HAVE_GETOPT_LONG
		static struct option options[] = {
			{ "usage",      0, 0, 'u'  },
			{ "ann",    	1, 0, 'A'  },
			{ "input",      1, 0, 'I'  },
			{ "nn_output",  1, 0, 'N'  }
		};
                opt = getopt_long(argc, argv, "uA:I:T:N:", options, NULL);
#else /* ! HAVE_GETOPT_LONG */
                opt = getopt(argc, argv, "uA:I:T:N:");
#endif /* ! HAVE_GETOPT_LONG */
                switch (opt) {
                case 'u':       print_usage();
                                break;
                case 'A':       ANNfile = optarg;
                                break;
                case 'I':       INPUTfile = optarg;
                                break;
                case 'N':       NN_OUTPUTfile = optarg;
                                break;
                default:        break;
                }
        } while (opt >= 0);

        /* Artificial Neural Network inizialization*/
        printf("LOADING DATA...\n");
        if (ANN_init(&net, ANNfile)) {
		fprintf(stderr, "Initialization error\n");
                return 1;
        }
        /* Input data acquisition*/
        N_sample = 0;
        if (ANN_DataRead(&INPUT, &N_sample, INPUTfile)) {
		fprintf(stderr, "Data input acquisition error\n");
                return 1;
        }

	
	if (matrix_init(&NN_OUTPUT, N_sample, net.N_output)) {
		fprintf(stderr, "MAtrix initailization error\n");
		return 1;
	}
	matrix_write(&INPUT, stdout, W_M_BIN);

        ANN_write(&net, stdout, ANN_W_A_TEXT);
	fprintf(stdout, "SIMULATION....\n");
        for (i = 0; i < N_sample; i++) {
                /* aggiorno il vettore degli ingressi */
		for (j = 0; j < net.N_input; j++) {
                        net.input.vec[j] = INPUT.mat[i][j];
                }
		/* simulo la rete */
                if (ANN_sim(&net, &net.input, &net.output, ANN_FEEDBACK_UPDATE)) {
			fprintf(stderr, "Network simulation error\n");
                        return 1;
                }
		/* aggiorno la matrice delle uscite */
                for (j = 0; j < net.N_output; j++) {
                        NN_OUTPUT.mat[i][j] = net.output.vec[j];
                }
        }

	ANN_DataWrite(&NN_OUTPUT, NN_OUTPUTfile);

        /* dynamic memory free*/
        ANN_destroy(&net);
	matrix_destroy(&INPUT);
	matrix_destroy(&NN_OUTPUT);
        
	return 0;
}

