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

/* Artificial Neural Network structure*/

typedef struct ANN{
	int N_input;     /* network input's number*/
	int N_output;    /* network (visible)output's number*/
	int N_layer;     /* network layer's number*/
	int *N_neuron;   /* neurons number*/
	double alpha;    /* activation function coefficient*/
	double beta;     /* activation function coefficient*/
	double eta;      /* learning rate*/
	double rho;      /* momentum term*/
	int r;           /* number of timestep delays*/
	double ***W;     /* network's synaptic weights*/
	double ***dW;    /* error gradient*/
	double **v;      /* neurons' internal activity*/
	double *****dy;  /* output network gradient*/
	double *yD;
	}ANN;


/*FUNCTIONS' PROTOTYPES:    */

int ANN_sim( ANN * , double * , double * );
