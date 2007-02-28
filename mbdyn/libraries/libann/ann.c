/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2007
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
 * Copyright (C) 2007
 *
 * Mattia Mattaboni	<mattaboni@aero.polimi.it>
 */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <malloc.h>
#include <stdio.h>

#include "ann.h"

int ANN_sim( ANN *net, double *input, double *output ){

	int i, j, k;
	double *IA, *y, *y_p;

	if( !(y_p = (double *)malloc( net->N_neuron[0] * sizeof(double) ) ) ){
		printf("\nError: malloc@ANN_sim(y_p)\n");
		return(0);
	}
	if( !(IA = (double *)malloc( net->N_neuron[0] * sizeof(double) ) ) ){
		printf("\nError: malloc@ANN_sim(IA)\n");
		return(0);
	}
	if( !(y = (double *)malloc( net->N_neuron[0] * sizeof(double) ) ) ){
		printf("\nError: malloc@ANN_sim(y)\n");
		return(0);
	}

	for( i=0;i<net->N_neuron[0];i++ ){
		y_p[i] = input[i];
	}


	for( i=0;i<net->N_layer+1;i++ ){
		if( !( IA = (double *)realloc( IA , net->N_neuron[i+1]*sizeof(double) ) ) ){
			printf("\nError: malloc@ANN_sim(IA)\n");
			return(0);
		}
		if( !( y = (double *)realloc( y , net->N_neuron[i+1]*sizeof(double) ) ) ){
			printf("\nError: malloc@ANN_sim(y)\n");
			return(0);
		}
		for( j=0;j<net->N_neuron[i+1];j++ ){
			IA[j] = 0.;
			for( k=0;k<net->N_neuron[i];k++ ){
				IA[j] += net->W[i][k][j]*y_p[k];
			}

			y[j] = ANN_InternalFunction(IA[j],net);
		}


		if( !( y_p = (double *)realloc( y_p , net->N_neuron[i+1]*sizeof(double) ) ) ){
			printf("\nError: malloc@ANN_sim(y_p)\n");
			return(0);
		}
		for( j=0;j<net->N_neuron[i+1];j++ ){
			y_p[j] = y[j];
		}
	}

	for( i=0;i<net->N_output;i++ ){
		output[i] = y_p[i];
	}

	if( net->r != 0 ){
		for( i=0;i<(net->r-1)*(net->N_neuron[net->N_layer+1]);i++ ){
			net->yD[((net->r)*(net->N_neuron[net->N_layer+1]))-1-i] = net->yD[((net->r-1)*(net->N_neuron[net->N_layer+1]))-1-i];
		}
		for( i=0;i<net->N_neuron[net->N_layer+1];i++ ){
			net->yD[i] = y_p[i];
		}
	}

	free(IA);
	free(y);
	free(y_p);

	return(1);

}
