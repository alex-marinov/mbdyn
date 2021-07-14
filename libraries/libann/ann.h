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

#ifndef ANN_H
#define ANN_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "ActivationFunction.h"
#include "matrix.h"

#define ANN_W_A_NONE		(0x0U)
#define ANN_W_A_TEXT		(0x1U)
#define ANN_W_A_BIN		(0x2U)

#define ANN_FEEDBACK_NONE	(0x0U)
#define ANN_FEEDBACK_UPDATE	(0x1U)

/* diagnostics */
typedef enum { 
        ANN_OK = 0,
        ANN_NO_MEMORY,
	ANN_NO_FILE,
	ANN_DATA_ERROR,
        ANN_MATRIX_ERROR,
	ANN_GEN_ERROR
} ann_res_t;

typedef enum {
	ANN_TM_BATCH,
	ANN_TM_SEQUENTIAL
} ann_training_mode_t;

typedef matrix*  ANN_vector_matrix;
typedef vector*  ANN_vector_vector;


/* Artificial Neural Network structure*/
typedef struct ANN {
	/* NEURON'S ACTIVATION FUNCTION */
        w_init_f        	w_init;
        w_destroy_f		w_destroy;
        w_read_f  		w_read;
        w_write_f       	w_write;
        w_eval_f        	w_eval;
        void            	*w_priv;

	/* NETWORK ARCHITECTURE */
        int 			N_input;    	/* network input number*/
        int 			N_output;	/* network (visible) output number*/
        int 			N_layer;    	/* network layer number*/
        int 			*N_neuron;  	/* neuron number*/
        int 			r;          	/* number of time delay*/

	/* TRAINING PARAMETERS */
        double 			eta;     	/* learning rate*/
        double 			rho;     	/* momentum term*/

	/* SYNAPTIC WEIGHTS */
        ANN_vector_matrix 	W;    		/* network synaptic weights*/

	/* JCOBIAN MATRIX */
	matrix			jacobian;

	/* SCALE FACTORS */
	matrix			input_scale;	/* input scale factors */
	matrix			output_scale;	/* output scale factors */

	/* SIMULATION DATA */
        ANN_vector_vector	v;	     	/* neurons' internal activity*/
        ANN_vector_vector	Y_neuron;  	/* neurons' output*/
        vector 			yD;      	/* output */

	/* TRAINING DATA */
        ANN_vector_matrix 	dEdW;		/* error gradient*/
        ANN_vector_matrix	**dy; 		/* output network gradient*/

	/* PRIVATE TRAINING DATA */
	ANN_vector_matrix	*dydW, dW;
	ANN_vector_vector	dXdW, temp, dydV, dEdV, dXdu;
	vector 			input, output, error;
	
	  
} ANN;


/* FUNCTIONS' PROTOTYPES  */
ann_res_t ANN_init( ANN *, const char * );
ann_res_t ANN_destroy( ANN * );
ann_res_t ANN_write( ANN *, FILE *, unsigned );
ann_res_t ANN_sim( ANN *, vector *, vector *, unsigned );
ann_res_t ANN_DataRead( matrix *, int *, char * );
ann_res_t ANN_DataWrite( matrix *, char * );
double ANN_InternalFunction( double , ANN * );
double ANN_InternalFunctionDer( double , ANN * );
ann_res_t ANN_vector_matrix_init( ANN_vector_matrix *, int *, int );
ann_res_t ANN_vector_vector_init( ANN_vector_vector *, int *, int );
ann_res_t ANN_dEdW( ANN *, vector *);
ann_res_t ANN_dXdW( ANN *, int , int , int );
ann_res_t ANN_WeightUpdate( ANN *, ANN_vector_matrix, double );
ann_res_t ANN_TrainingEpoch( ANN *, matrix *, matrix *, matrix *, int, ann_training_mode_t);
ann_res_t ANN_reset( ANN * );
ann_res_t ANN_TotalError( matrix *, matrix *, double *);
ann_res_t ANN_vector_matrix_ass( ANN_vector_matrix *, ANN_vector_matrix *, int* , int, double );
ann_res_t ANN_jacobian_matrix( ANN *, matrix * );

void ANN_error( ann_res_t, char * );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* ANN_H */
