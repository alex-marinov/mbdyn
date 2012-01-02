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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <string.h>

#include "ann.h"

static int TRAINING_MODE = 1;

void
print_usage(void)
{
	fprintf(stdout, "\nUSAGE OPTIONS:\n"
                "  -u, --usage\n"
                "       print usage\n"
                "  -f, --file\n"
                "       input file name\n"
                "  -m, --mode\n"
                "       scaling mode: 1 mean-variance  	(default)\n"
                "                     2 min-max\n" 
                "  -M, --mean\n"
                "       mean value (defaul 0.)\n"
                "  -V, --variance"
                "       variance value (defaul 1.)\n"
                "  -N, --min\n"
                "       minimum value (defaul -1.)\n"
                "  -X, --max\n"
                "       maximum value (defaul 1.)\n"
	);
        exit(0);
}


static int MODE = 1;
static float MEAN = 0.;
static float VARIANCE = 1.;
static float MIN = -1.;
static float MAX = 1.;
static char *filename = "data.dat";

int
main(int argc, char *argv[])
{
	int opt;
        extern char *optarg;
	matrix MAT, SF;
	int N_sample, i;
	double mean, var, min, max;

        /* 0. Training options */
        do {
#ifdef HAVE_GETOPT_LONG
		static struct option options[] = {
			{ "usage",      0, 0, 'u'  },
			{ "mode",       1, 0, 'm'  },
			{ "mean",    	1, 0, 'M'  },
			{ "variance",  	1, 0, 'V'  },
			{ "min",    	1, 0, 'N'  },
			{ "max",  	1, 0, 'X'  },
			{ "file",    	1, 0, 'f'  } 
		};
		opt = getopt_long(argc, argv, "uf:m:M:V:N:X:", options, NULL);
#else /* ! HAVE_GETOPT_LONG */
		opt = getopt(argc, argv, "uf:m:M:V:N:X:");
#endif /* ! HAVE_GETOPT_LONG */

                switch (opt) {
                case 'u':       print_usage();
                                break; 
		case 'f':       filename = optarg;
                                break;
                case 'm':       MODE = atoi(optarg);
                                break;
                case 'M':       MEAN = atof(optarg);
                                break;
                case 'V':       VARIANCE = atof(optarg);
                                break;
                case 'N':       MIN = atof(optarg);
                                break;
                case 'X':       MAX = atof(optarg);
                                break;
                default:        break;
                }
	} while (opt >= 0);
	
        N_sample = 0;

        if (ANN_DataRead(&MAT, &N_sample, filename)) {
                fprintf(stderr, "Error in Input data acquisition\n");
                return 1; 
        }	

	matrix_init(&SF, MAT.Ncolumn, 2);

	switch (MODE) {
	case 1:
		for (i = 0; i < MAT.Ncolumn; i++) {
			mean = mean_value(&MAT, i);
			var = variance(&MAT, i);
			SF.mat[i][0] = sqrt(VARIANCE/var);
			SF.mat[i][1] = MEAN - SF.mat[i][0]*mean;
		}
		break;

	case 2:
		for (i = 0; i < MAT.Ncolumn; i++) {
			min = minimum(&MAT, i);
			max = maximum(&MAT, i);
			printf("AAAA   %lf   %lf\n", min, max);
			SF.mat[i][0] = (MAX - MIN)/(max - min);
			SF.mat[i][1] = (max*MIN-MAX*min)/(max-min);
		}
		break;
	}

	matrix_write(&SF, stdout, W_M_TEXT);
	
	return 0;
}
