/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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
 * Copyright (C) 2007-2014
 *
 * Mattia Mattaboni	<mattaboni@aero.polimi.it>
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ann.h"

/* Ininzializza ed alloca la memoria per la struttura dati
 * contenente le informazioni relative alla rete neurale */

ann_res_t ANN_init( ANN *net, const char * FileName){

        int i,j;
	int ActFnc;
        FILE *fh;

	memset( net, 0, sizeof( ANN ) );

        if ( !(fh = fopen( FileName, "r" ) ) ){
		fprintf( stderr, "Input file doesn't exist.\n" );
		ANN_error( ANN_NO_FILE, "ANN_init" );
                return ANN_NO_FILE;
        }

        fscanf(fh,"%d",&(net->N_input));
	if( net->N_input <= 0 ){
		fprintf( stderr, "Input number must be greater than zero.\n");
		ANN_error( ANN_DATA_ERROR, "ANN_init" );
		return ANN_DATA_ERROR;
	}
        fscanf(fh,"%d",&(net->N_output));
	if( net->N_output <= 0 ){
		fprintf( stderr, "Output number must be greater than zero.\n");
		ANN_error( ANN_DATA_ERROR, "ANN_init" );
		return ANN_DATA_ERROR;
	}
        fscanf(fh,"%d",&(net->N_layer));
	if( net->N_layer < 0 ){
		fprintf( stderr, "Hidden layer number must be not negative.\n");
		ANN_error( ANN_DATA_ERROR, "ANN_init" );
		return ANN_DATA_ERROR;
	}
        fscanf(fh,"%d",&(net->r));
	if( net->r < 0 ){
		fprintf( stderr, "Timestep delay number must be not negative.\n");
		ANN_error( ANN_DATA_ERROR, "ANN_init" );
		return ANN_DATA_ERROR;
	}

        if( !(net->N_neuron = (int *)malloc( (net->N_layer+2) * sizeof(int) ) ) ){
			ANN_error( ANN_NO_MEMORY, "ANN_init" );
                        return ANN_NO_MEMORY;
        }
        net->N_neuron[0] = 0;
        for( i=0;i<net->N_layer+1;i++ ){
                fscanf(fh,"%d",&(net->N_neuron[i+1]));
		if( net->N_neuron[i+1] <= 0 ){
			fprintf( stderr, "Neuron number at %d layer must be greater than zero.\n",i+1);
			ANN_error( ANN_DATA_ERROR, "ANN_init" );
			return ANN_DATA_ERROR;
		}
        }
        net->N_neuron[0] = net->N_input + ( net->r * net->N_neuron[net->N_layer+1] );

        fscanf( fh, "%d", &ActFnc); 
	switch( ActFnc ){
	case 1:		/* use tanh */
        		net->w_init = w_tanh_init;
		        net->w_destroy = w_tanh_destroy;
		        net->w_read = w_tanh_read;
		        net->w_write = w_tanh_write;
		        net->w_eval = w_tanh_eval;	
			break;
	case 2:		/* use linear activation function */
        		net->w_init = w_linear_init;
		        net->w_destroy = w_linear_destroy;
		        net->w_read = w_linear_read;
		        net->w_write = w_linear_write;
		        net->w_eval = w_linear_eval;	
			break;
	default:	
			fprintf( stderr, "Unknown activation function\n" );
			ANN_error( ANN_DATA_ERROR, "ANN_init" );
			return ANN_DATA_ERROR;
	}
        /* end of switch */

        net->w_priv = NULL;
        if (net->w_init(&net->w_priv) != 0) {
                /* error */
        }

        net->w_read(net->w_priv, fh, W_F_NONE);

        fscanf( fh, "%le", &(net->eta));
	if( net->eta < 0 ){
		fprintf( stderr, "Learning rate must be not negative.\n");
		ANN_error( ANN_DATA_ERROR, "ANN_init" );
		return ANN_DATA_ERROR;
	}
        fscanf( fh, "%le", &(net->rho));

	ANN_vector_matrix_init(&net->W, net->N_neuron, net->N_layer);
	ANN_vector_matrix_init(&net->dEdW, net->N_neuron, net->N_layer);
	ANN_vector_matrix_init(&net->dW, net->N_neuron, net->N_layer);

        for( i=0; i<net->N_layer+1; i++ ){
		matrix_read( &net->W[i], fh, W_M_BIN );
        }
        if( matrix_init( &net->jacobian, net->N_input, net->N_output) ){
		ANN_error( ANN_MATRIX_ERROR, "ANN_init" );
                return ANN_MATRIX_ERROR;
        }
        if( matrix_init( &net->input_scale, net->N_input, 2) ){
		ANN_error( ANN_MATRIX_ERROR, "ANN_init" );
                return ANN_MATRIX_ERROR;
        }
        if( matrix_init( &net->output_scale, net->N_output, 2 ) ){
		ANN_error( ANN_MATRIX_ERROR, "ANN_init" );
                return ANN_MATRIX_ERROR;
        }
	matrix_read( &net->input_scale, fh, W_M_BIN );
	matrix_read( &net->output_scale, fh, W_M_BIN );

	ANN_vector_vector_init(&net->v, net->N_neuron, net->N_layer);
	ANN_vector_vector_init(&net->Y_neuron, net->N_neuron, net->N_layer);
	ANN_vector_vector_init(&net->dXdW, net->N_neuron, net->N_layer);
	ANN_vector_vector_init(&net->temp, net->N_neuron, net->N_layer);
	ANN_vector_vector_init(&net->dydV, net->N_neuron, net->N_layer);
	ANN_vector_vector_init(&net->dEdV, net->N_neuron, net->N_layer);
	ANN_vector_vector_init(&net->dXdu, net->N_neuron, net->N_layer);
	
	if( net->r > 0 ){
        	if( vector_init( &net->yD, net->r*net->N_neuron[net->N_layer+1] ) ){
			ANN_error( ANN_MATRIX_ERROR, "ANN_init" );
               		return ANN_MATRIX_ERROR;
        	}
 	}
        if( vector_init( &net->input, net->N_input ) ){
		ANN_error( ANN_MATRIX_ERROR, "ANN_init" );
                return ANN_MATRIX_ERROR;
        }
        if( vector_init( &net->output, net->N_output ) ){
		ANN_error( ANN_MATRIX_ERROR, "ANN_init" );
                return ANN_MATRIX_ERROR;
        }
        if( vector_init( &net->error, net->N_output ) ){
		ANN_error( ANN_MATRIX_ERROR, "ANN_init" );
                return ANN_MATRIX_ERROR;
        }
	
	
        if( !(net->dydW = (ANN_vector_matrix *)calloc( net->N_neuron[net->N_layer+1], sizeof(ANN_vector_matrix) ) ) ){
                ANN_error( ANN_NO_MEMORY, "ANN_init" );
		return ANN_NO_MEMORY;
        }
	for( i=0; i<net->N_neuron[net->N_layer+1]; i++ ){
		ANN_vector_matrix_init(&net->dydW[i], net->N_neuron, net->N_layer);
	}

        if( !(net->dy = (ANN_vector_matrix **)calloc( net->r, sizeof(ANN_vector_matrix *) ) ) ){
                ANN_error( ANN_NO_MEMORY, "ANN_init" );
		return ANN_NO_MEMORY;
        }
	for( i=0; i<net->r; i++){
        	if( !(net->dy[i] = (ANN_vector_matrix *)calloc( net->N_neuron[net->N_layer+1], sizeof(ANN_vector_matrix) ) ) ){
                	ANN_error( ANN_NO_MEMORY, "ANN_init" );
	                return ANN_NO_MEMORY;
       		}
		for( j=0; j<net->N_neuron[net->N_layer+1]; j++ ){
			ANN_vector_matrix_init( &net->dy[i][j], net->N_neuron, net->N_layer );
		}
	}	

        fclose(fh);

        return ANN_OK;
}

/* distrugge una rete neurale */
ann_res_t ANN_destroy( ANN * net ){

        int i,j,k;

        if( vector_destroy( &net->yD ) ){
		ANN_error( ANN_MATRIX_ERROR, "ANN_destry" );
		return ANN_MATRIX_ERROR;
	}
        if( vector_destroy( &net->input ) ){
		ANN_error( ANN_MATRIX_ERROR, "ANN_destry" );
		return ANN_MATRIX_ERROR;
	}
        if( vector_destroy( &net->output ) ){
		ANN_error( ANN_MATRIX_ERROR, "ANN_destry" );
		return ANN_MATRIX_ERROR;
	}
        if( vector_destroy( &net->error ) ){
		ANN_error( ANN_MATRIX_ERROR, "ANN_destry" );
		return ANN_MATRIX_ERROR;
	}
        if( matrix_destroy( &net->jacobian ) ){
		ANN_error( ANN_MATRIX_ERROR, "ANN_destry" );
		return ANN_MATRIX_ERROR;
	}
        if( matrix_destroy( &net->input_scale ) ){
		ANN_error( ANN_MATRIX_ERROR, "ANN_destry" );
		return ANN_MATRIX_ERROR;
	}
        if( matrix_destroy( &net->output_scale ) ){
		ANN_error( ANN_MATRIX_ERROR, "ANN_destry" );
		return ANN_MATRIX_ERROR;
	}

        for( i=0;i<net->N_layer+1;i++ ){
                if( matrix_destroy( &net->W[i] ) ){
			ANN_error( ANN_MATRIX_ERROR, "ANN_destroy" );
			return ANN_MATRIX_ERROR;
		}
                if( matrix_destroy( &net->dEdW[i] ) ){
			ANN_error( ANN_MATRIX_ERROR, "ANN_destroy" );
			return ANN_MATRIX_ERROR;
		}
                if( matrix_destroy( &net->dW[i] ) ){
			ANN_error( ANN_MATRIX_ERROR, "ANN_destroy" );
			return ANN_MATRIX_ERROR;
		}
        }
        free(net->dEdW);
        free(net->W);
        free(net->dW);

        for( i=0;i<net->N_layer+2;i++ ){
                if( vector_destroy( &net->v[i] ) ){
			ANN_error( ANN_MATRIX_ERROR, "ANN_destroy" );
			return ANN_MATRIX_ERROR;
		}
                if( vector_destroy( &net->Y_neuron[i] ) ){
			ANN_error( ANN_MATRIX_ERROR, "ANN_destroy" );
			return ANN_MATRIX_ERROR;
		}
                if( vector_destroy( &net->dXdW[i] ) ){
			ANN_error( ANN_MATRIX_ERROR, "ANN_destroy" );
			return ANN_MATRIX_ERROR;
		}
                if( vector_destroy( &net->dXdu[i] ) ){
			ANN_error( ANN_MATRIX_ERROR, "ANN_destroy" );
			return ANN_MATRIX_ERROR;
		}
                if( vector_destroy( &net->temp[i] ) ){
			ANN_error( ANN_MATRIX_ERROR, "ANN_destroy" );
			return ANN_MATRIX_ERROR;
		}
                if( vector_destroy( &net->dydV[i] ) ){
			ANN_error( ANN_MATRIX_ERROR, "ANN_destroy" );
			return ANN_MATRIX_ERROR;
		}
                if( vector_destroy( &net->dEdV[i] ) ){
			ANN_error( ANN_MATRIX_ERROR, "ANN_destroy" );
			return ANN_MATRIX_ERROR;
		}
        }
        free(net->v);
        free(net->Y_neuron);
        free(net->dXdW);
        free(net->dXdu);
        free(net->temp);
        free(net->dydV);
        free(net->dEdV);
	
	for( i=0; i<net->N_neuron[net->N_layer+1]; i++ ){
		for( j=0; j<net->N_layer+1; j++ ){
			if( matrix_destroy( &net->dydW[i][j] )) {
				ANN_error( ANN_MATRIX_ERROR, "ANN_destroy" );
				return ANN_MATRIX_ERROR;
			}
		}
		free(net->dydW[i]);
	}	
	free(net->dydW);
	
	for( k=0; k<net->r; k++ ){	
		for( i=0; i<net->N_neuron[net->N_layer+1]; i++ ){
			for( j=0; j<net->N_layer+1; j++ ){
				if( matrix_destroy( &net->dy[k][i][j] )){
					ANN_error( ANN_MATRIX_ERROR, "ANN_destroy" );
					return ANN_MATRIX_ERROR;
				}
			}
			free(net->dy[k][i]);
		}	
		free( net->dy[k] );
	}
	free(net->dy);
        free(net->N_neuron);
        
	if (net->w_destroy(net->w_priv) != 0) {
                /* error */
        }

	return ANN_OK;

}
/* scrive su file o a video una rete neurale */
ann_res_t ANN_write( ANN *net, FILE * fh, unsigned flags ){

        int i;

	if( flags & ANN_W_A_TEXT ){
		fprintf( fh, "ARTIFICIAL NEURAL NETWORK\n");
		fprintf( fh, "Network topology\n");
	        fprintf( fh, "-Input number: %d \n", net->N_input);
       		fprintf( fh, "-Output number: %d \n", net->N_output);
	        fprintf( fh, "-Hidden layers number: %d \n", net->N_layer);

	        for( i=0;i<net->N_layer+1;i++ ){
	               fprintf( fh, "-Neurons number (layer number %d) : %d\n", i+1, net->N_neuron[i+1]);
       	 	}
	        fprintf( fh, "-Time delay number: %d \n", net->r);

		fprintf( fh, "Training parameters\n");
	        fprintf( fh, "-Learning rate: %e \n", net->eta);
       		fprintf( fh, "-Momentum term: %e \n", net->rho);

		fprintf( fh, "Activation function parameters\n");
	        net->w_write(net->w_priv, fh, W_F_TEXT);

		fprintf( fh, "Synaptic weight\n");

	        for( i=0; i<net->N_layer+1; i++ ){
        	        if( i!=net->N_layer )
                	        fprintf( fh, "-Layer number %d :\n",i+1);
                	else
                        	fprintf( fh, "-Visible layer :\n");
	                if( matrix_write( &net->W[i], fh, W_M_TEXT ) ){
				ANN_error( ANN_MATRIX_ERROR, "ANN_write" );
				return ANN_MATRIX_ERROR;
			}
       		}
		fprintf( fh, "Input scale factors\n" );
	        if( matrix_write( &net->input_scale, fh, W_M_TEXT ) ){
			ANN_error( ANN_MATRIX_ERROR, "ANN_write" );
			return ANN_MATRIX_ERROR;
		}
		fprintf( fh, "Output scale factors\n" );
	        if( matrix_write( &net->output_scale, fh, W_M_TEXT ) ){
			ANN_error( ANN_MATRIX_ERROR, "ANN_write" );
			return ANN_MATRIX_ERROR;
		}
	}
	if( flags & ANN_W_A_BIN ){
	        fprintf( fh, " %d \n", net->N_input);
       		fprintf( fh, " %d \n", net->N_output);
	        fprintf( fh, " %d \n", net->N_layer);
	        fprintf( fh, " %d \n", net->r);

	        for( i=0;i<net->N_layer+1;i++ ){
	               fprintf( fh, " %d ", net->N_neuron[i+1]);
       	 	}
		fprintf( fh, "\n\n" );

        	net->w_write(net->w_priv, fh, W_F_BIN);
		fprintf( fh, "\n" );

	        fprintf( fh, " %e \n", net->eta);
       		fprintf( fh, " %e \n", net->rho);
		fprintf( fh, "\n" );

	        for( i=0; i<net->N_layer+1; i++ ){
	                if( matrix_write( &net->W[i], fh, W_M_BIN ) ){
				ANN_error( ANN_MATRIX_ERROR, "ANN_write" );
				return ANN_MATRIX_ERROR;
			}
			fprintf( fh, "\n\n" );
       		}
	        if( matrix_write( &net->input_scale, fh, W_M_BIN ) ){
			ANN_error( ANN_MATRIX_ERROR, "ANN_write" );
			return ANN_MATRIX_ERROR;
		}
	        if( matrix_write( &net->output_scale, fh, W_M_BIN ) ){
			ANN_error( ANN_MATRIX_ERROR, "ANN_write" );
			return ANN_MATRIX_ERROR;
		}
	}
	return ANN_OK;
}

/* calcola l'uscita della rete noti gli ingressi */
ann_res_t ANN_sim( ANN *net , vector *input , vector *output, unsigned flags ){

        int i,j;

	/* costruisco il vettore degli ingressi
	 * incolonnando agli ingressi esterni le retroazioni
	 * delle uscite precedenti */

        for( i=0; i<net->N_input; i++ ){
		net->Y_neuron[0].vec[i] = net->input_scale.mat[i][1] + input->vec[i]*net->input_scale.mat[i][0];
        }
	for( i=0; i<net->N_neuron[net->N_layer+1]*net->r; i++ ){
		net->Y_neuron[0].vec[i+net->N_input] = net->yD.vec[i];
	}

        for( i=0;i<net->N_neuron[0];i++ ){
                net->v[0].vec[i] = net->Y_neuron[0].vec[i];
	}
	
	/* calcolo il vettore delle uscite */
        for( i=0;i<net->N_layer+1;i++ ){
                if( matrixT_vector_prod( &net->W[i], &net->Y_neuron[i] ,&net->v[i+1] ) ){
			ANN_error( ANN_MATRIX_ERROR, "ANN_sim" );
			return ANN_MATRIX_ERROR;
		}
                for( j=0; j<net->N_neuron[i+1];j++ ){
                        net->Y_neuron[i+1].vec[j] = ANN_InternalFunction(net->v[i+1].vec[j], net );
                }
        }

        for( i=0;i<net->N_output;i++ ){
                //output->vec[i] = net->output_scale.mat[i][1] + net->Y_neuron[net->N_layer+1].vec[i]*net->output_scale.mat[i][0];
                output->vec[i] = (net->Y_neuron[net->N_layer+1].vec[i]-net->output_scale.mat[i][1])/net->output_scale.mat[i][0];
        }
	
	/* aggiorno il vettore delle uscite retroazionate */
	if( flags & ANN_FEEDBACK_UPDATE ){
	        for( i=0;i<(net->r-1)*(net->N_neuron[net->N_layer+1]);i++ ){
        	        net->yD.vec[((net->r)*(net->N_neuron[net->N_layer+1]))-1-i] = net->yD.vec[((net->r-1)*(net->N_neuron[net->N_layer+1]))-1-i];
        	}
        	if( net->r != 0 ){
               		for( i=0;i<net->N_neuron[net->N_layer+1];i++ ){
                        	net->yD.vec[i] = net->Y_neuron[net->N_layer+1].vec[i];
	                }
       		}
	}

        return ANN_OK;
}
/* legge da file una matrice */
ann_res_t ANN_DataRead( matrix *MAT, int *N_sample, char *FileName ){

	int Nrow, Ncolumn;
	FILE *fh;

	if( !( fh = fopen( FileName, "r" ) ) ){
		ANN_error( ANN_NO_FILE, "ANN_DataRead" );
		return ANN_NO_FILE;
	}

        fscanf( fh, "%d", &Nrow);
        fscanf( fh, "%d", &Ncolumn);
	if( matrix_init( MAT, Nrow, Ncolumn ) ){
		ANN_error( ANN_MATRIX_ERROR, "ANN_DataRead" );
		return ANN_MATRIX_ERROR;
	}

	if(  matrix_read( MAT, fh, W_M_BIN) ){
		ANN_error( ANN_MATRIX_ERROR, "ANN_DataRead" );
		return ANN_MATRIX_ERROR;
	}

	fclose(fh);

	*N_sample = Nrow;

	return ANN_OK;
}
/* scrive su file una matrice */
ann_res_t ANN_DataWrite( matrix *MAT, char *FileName ){

	FILE *fh;

	if( !( fh = fopen( FileName, "w" ) ) ){
		ANN_error( ANN_NO_FILE, "ANN_DataWrite" );
		return ANN_NO_FILE;
	}

	fprintf( fh, "%d %d", MAT->Nrow, MAT->Ncolumn );
	
	if(  matrix_write(MAT, fh, W_M_BIN) ){
		ANN_error( ANN_MATRIX_ERROR, "ANN_DataWrite" );
		return ANN_MATRIX_ERROR;
	}

	fclose(fh);

	return ANN_OK;
}
/* calcola l'uscita di un neurone nota la sua attività
 * interna */
double ANN_InternalFunction( double v , ANN * net ){

        double y;

        if (net->w_eval(net->w_priv, v, 0, &y) != 0) {
                ANN_error( ANN_GEN_ERROR, "ANN_InternalFunction" );
		return ANN_GEN_ERROR;
        }

        return y;
}

/* calcola la derivata dell'uscita di un neurone nota la 
 * sua attività interna */
double ANN_InternalFunctionDer( double v , ANN * net ){

        double y;

        if (net->w_eval(net->w_priv, v, 1, &y) != 0) {
                ANN_error( ANN_GEN_ERROR, "ANN_InternalFunctionDer" );
		return ANN_GEN_ERROR;
        }
        return y;
}

/* funzione che aggiorna i pesi sinaptici */
ann_res_t ANN_WeightUpdate( ANN *net , ANN_vector_matrix DW , double K){

        int i;

        for( i=0;i<net->N_layer+1;i++ ){
		if( matrix_sum( &net->W[i], &DW[i], &net->W[i], K ) ){
			ANN_error( ANN_MATRIX_ERROR, "ANN_WeightUpdate" );
			return ANN_MATRIX_ERROR;
		}
        }

	return ANN_OK;

}

/* inizializzo una struttura di vettori di matrici associata all'architettura
 * della rete neurale */
ann_res_t ANN_vector_matrix_init( ANN_vector_matrix *vm, int *N_neuron, int N_layer ){

	int i;

        if( !( *vm = (matrix *)calloc( N_layer+1, sizeof(matrix) ) ) ){
		ANN_error( ANN_NO_MEMORY, "ANN_vector_matrix_init" );
                return ANN_NO_MEMORY;
        }
        for( i=0; i<N_layer+1; i++ ){
                if( matrix_init( &(*vm)[i], N_neuron[i], N_neuron[i+1] ) ){
			ANN_error( ANN_MATRIX_ERROR, "ANN_vector_matrix_init" );
                        return ANN_MATRIX_ERROR;
                }
        }

	return ANN_OK;
}
/* inizializzo una struttura di vettori di vettori associata all'architettura
 * della rete neurale */
ann_res_t ANN_vector_vector_init( ANN_vector_vector *vv, int *N_neuron, int N_layer ){

	int i;

        if( !( *vv = (ANN_vector_vector)calloc( N_layer+2, sizeof(vector) ) ) ){
		ANN_error( ANN_NO_MEMORY, "ANN_vector_vector_init" );
                return ANN_NO_MEMORY;
        }
        for( i=0; i<N_layer+2; i++ ){
                if( vector_init( &(*vv)[i], N_neuron[i] ) ){
			ANN_error( ANN_MATRIX_ERROR, "ANN_vector_vector_init" );
                        return ANN_MATRIX_ERROR;
                }
        }

	return ANN_OK;
}


/* calcolo la derivata dell'uscita dell'N-esimo strato di neuroni rispetto
 * al peso I,J del medesimo stato */
ann_res_t ANN_dXdW( ANN * net , int I , int J , int N ){

        int i,j,k;
	
	/* inizializzo con la derivata degli ingressi rispetto al peso I,J
	 * delle strato N */

	/* gli ingressi esterni non dipendono dai pesi */
        i = 0;
        for( j=0; j<net->N_input; j++ ){
                net->dXdW[0].vec[i] = 0.;
                i++;
        }
        for( k=0;k<net->r;k++ ){
                for( j=0;j<net->N_neuron[net->N_layer+1];j++ ){
                        net->dXdW[0].vec[i] = net->dy[k][j][N].mat[I][J];
                        i++;
                }
        }
	/* calcolo la derivata allo strato di interesse a partire da quella degli ingressi*/
        for( i=0;i<N;i++ ){
                if( matrixT_vector_prod( &net->W[i], &net->dXdW[i] ,&net->temp[i+1] ) ){
			ANN_error( ANN_MATRIX_ERROR, "ANN_dXdW" );
			return ANN_MATRIX_ERROR;
		}
                for( j=0; j<net->N_neuron[i+1];j++ ){
                        net->dXdW[i+1].vec[j] = ANN_InternalFunctionDer(net->v[i+1].vec[j], net )*net->temp[i+1].vec[j];
                }
        }

        return ANN_OK;
}


/* funzione che calcola la derivata dell'errore rispetto a tutti pesi della
 * rete */
ann_res_t ANN_dEdW( ANN * net , vector *e ){

        int i,j,k,p,l,q;
        double temp;

        /* Output gradient ( visible layer )*/
        for( k=0;k<net->N_neuron[net->N_layer];k++ ){
                for( l=0;l<net->N_neuron[net->N_layer+1];l++ ){
                        if( ANN_dXdW( net , k , l , net->N_layer ) ){
				ANN_error( ANN_GEN_ERROR, "ANN_dEdW" );
                                return ANN_GEN_ERROR;
                        }
                	if (matrixT_vector_prod( &net->W[net->N_layer], &net->dXdW[net->N_layer] ,&net->temp[net->N_layer+1] ) ){
				ANN_error( ANN_MATRIX_ERROR, "ANN_dEdW" );
				return ANN_MATRIX_ERROR;
			}
                        for( j=0; j<net->N_neuron[net->N_layer+1]; j++ ){
                                net->dydW[j][net->N_layer].mat[k][l] = ANN_InternalFunctionDer(net->v[net->N_layer+1].vec[j], net)*( net->temp[net->N_layer+1].vec[j] + net->Y_neuron[net->N_layer].vec[k]*(l==j) );
                        }

			temp = 0.;
                        for( j=0; j<net->N_output; j++ ){
				temp += -net->dydW[j][net->N_layer].mat[k][l]*e->vec[j];
                        }
			net->dEdW[net->N_layer].mat[k][l] = net->rho*net->dEdW[net->N_layer].mat[k][l] - net->eta*temp;
                }
        }

	/* Output gradient (hidden layer) */
        for( q=0; q<net->N_neuron[net->N_layer+1]; q++ ){
                for( i=0; i<(net->N_neuron[net->N_layer+1] ); i++ ){
                        net->dydV[net->N_layer+1].vec[i] = 0.;
                }
                net->dydV[net->N_layer+1].vec[q] = ANN_InternalFunctionDer(net->v[net->N_layer+1].vec[q],net);

                for( k=0;k<net->N_layer;k++ ){
        		if( matrix_vector_prod( &net->W[net->N_layer-k], &net->dydV[net->N_layer-k+1] ,&net->temp[net->N_layer-k] ) ){	
				ANN_error( ANN_MATRIX_ERROR, "ANN_dEdW" );
				return ANN_MATRIX_ERROR;
			}
                        for( j=0;j<net->N_neuron[net->N_layer-k];j++ ){
                                net->dydV[net->N_layer-k].vec[j] = net->temp[net->N_layer-k].vec[j]*ANN_InternalFunctionDer(net->v[net->N_layer-k].vec[j],net);
                        }
                       for( i=0;i<net->N_neuron[net->N_layer-k-1];i++ ){
                                for( j=0;j<net->N_neuron[net->N_layer-k];j++ ){
                                        if( ANN_dXdW( net ,i ,j ,(net->N_layer-k-1) ) ){
						ANN_error( ANN_GEN_ERROR, "ANN_dEdW" );
		                                return ANN_GEN_ERROR;
                                        }
                                        if( net->N_layer-k-1 != 0 ){
                                                temp = ANN_InternalFunction(net->v[net->N_layer-k-1].vec[i],net);
                                        }
                                        else{
                                                temp = (net->v[0].vec[i]);
                                        }

                                        for( p=0;p<(net->N_neuron[net->N_layer-k-1]);p++ ){
                                                temp += net->W[net->N_layer-k-1].mat[p][j]*net->dXdW[net->N_layer-k-1].vec[p];
                                        }
                                        //net->dydW[q][net->N_layer-1-k].mat[i][j] = net->dydV[net->N_layer+1-k].vec[j]*temp;
                                        net->dydW[q][net->N_layer-1-k].mat[i][j] = net->dydV[net->N_layer-k].vec[j]*temp;
                                }
                        }
                }
        }
	/* calcolo la derivata dell'errore rispetto ai pesi degli
	 * degli strati non visibili */
        for( i=0; i<(net->N_output); i++ ){
                net->dEdV[net->N_layer+1].vec[i] = -e->vec[i]*ANN_InternalFunctionDer(net->v[net->N_layer+1].vec[i],net);
        }
	
	 for( k=0; k<net->N_layer; k++ ){
        	matrix_vector_prod( &net->W[net->N_layer-k], &net->dEdV[net->N_layer-k+1] ,&net->temp[net->N_layer-k] );
                for( j=0;j<net->N_neuron[net->N_layer-k];j++ ){
                	net->dEdV[net->N_layer-k].vec[j] = net->temp[net->N_layer-k].vec[j]*ANN_InternalFunctionDer(net->v[net->N_layer-k].vec[j],net);
                }
                for( i=0;i<net->N_neuron[net->N_layer-k-1];i++ ){
                        for( j=0;j<net->N_neuron[net->N_layer-k];j++ ){
                                if( ANN_dXdW( net , i , j , (net->N_layer-k-1) ) ){
					ANN_error( ANN_GEN_ERROR, "ANN_dEdW" );
		                        return ANN_GEN_ERROR;
                                }

                                if( net->N_layer-k-1 != 0 ){
                                        temp = ANN_InternalFunction(net->v[net->N_layer-k-1].vec[i],net);
                                }
                                else{
                                        temp = (net->v[0].vec[i]);
                                }

                                for( p=0;p<(net->N_neuron[net->N_layer-k-1]);p++ ){
                                        temp += net->W[net->N_layer-k-1].mat[p][j]*net->dXdW[net->N_layer-k-1].vec[p];
                                }
                                //net->dEdW[net->N_layer-1-k].mat[i][j] = net->rho*net->dEdW[net->N_layer-1-k].mat[i][j] - net->eta*net->dEdV[net->N_layer+1-k].vec[j]*temp;
                                net->dEdW[net->N_layer-1-k].mat[i][j] = net->rho*net->dEdW[net->N_layer-1-k].mat[i][j] - net->eta*net->dEdV[net->N_layer-k].vec[j]*temp;
                        }
                }
        }
	/* aggiorno la struttura dati contenente la derivata di tutte le uscite rispetto
	 * a tutti i pesi della rete salvata per gli r passi precedenti */

        for( p=0;p<(net->r-1);p++ ){
		for( i=0; i<net->N_neuron[net->N_layer+1];i++ ){	
			if( ANN_vector_matrix_ass( &net->dy[net->r-1-p][i], &net->dy[net->r-2-p][i], net->N_neuron, net->N_layer, 1. ) ){
				ANN_error( ANN_GEN_ERROR, "ANN_dEDW" );
				return ANN_GEN_ERROR;
			}
		}
        }
	if( net->r != 0 ){
		for( i=0; i<net->N_neuron[net->N_layer+1];i++ ){	
			if( ANN_vector_matrix_ass( &net->dy[0][i], &net->dydW[i], net->N_neuron, net->N_layer, 1. ) ){
				ANN_error( ANN_GEN_ERROR, "ANN_dEdW" );
				return ANN_GEN_ERROR;
			}
		}
	}

	return ANN_OK;
}


/* funzione che esegue un epoca di addestramento in modalità BATCH o 
 * SEQUENTIAL */
ann_res_t ANN_TrainingEpoch( ANN * net, matrix *INPUT, matrix *DES_OUTPUT, matrix *NN_OUTPUT, int N_sample, ann_training_mode_t mode ){

	int i,t;

	for( t=0; t<N_sample; t++ ){

		/* costrisco il vettore degli ingressi al tempo t
		 * partendo dalla matrice degli ingressi INPUT */
		for( i=0; i<net->N_input; i++ ){
			net->input.vec[i] = INPUT->mat[t][i];
		}
		/* simulo la rete per calcolare le uscite al passo t
		 * e il corrispondente errore */
		if( ANN_sim( net, &net->input, &net->output, ANN_FEEDBACK_UPDATE) ){
			ANN_error( ANN_GEN_ERROR, "ANN_TrainigEpoch" );
			return ANN_GEN_ERROR;
		}
		for( i=0; i<net->N_output; i++ ){
			net->error.vec[i] = DES_OUTPUT->mat[t][i]-net->output.vec[i];
			NN_OUTPUT->mat[t][i] = net->output.vec[i];
		}
		/* calcolo la derivata dell'errore rispetto a tutti i pesi
		 * sinaptici */
		if( ANN_dEdW( net, &net->error) ){
			ANN_error( ANN_GEN_ERROR, "ANN_TrainingEpoch" );
			return ANN_GEN_ERROR;
		}
		/* addestramento in modalità BATCH: accumulo la derivata
		 * e la applico solamente alla fine di un epoca di 
		 * addestramento 		 */
		if( mode == ANN_TM_BATCH ){
        		for( i=0;i<net->N_layer+1;i++ ){
				if( matrix_sum( &net->dW[i], &net->dEdW[i], &net->dW[i], 1. ) ){
					ANN_error( ANN_MATRIX_ERROR, "ANN_TrainingEpoch" );
					return ANN_MATRIX_ERROR;
				}
       			 }
		}
		/* addestramento in modalità SEQUENTIAL: applico 
		 * immediatamente la variazione dei pesi */
		if( mode == ANN_TM_SEQUENTIAL ){
			if( ANN_WeightUpdate( net, net->dEdW, 1. ) ){
				ANN_error( ANN_GEN_ERROR, "ANN_TrainingEpoch" );
				return ANN_GEN_ERROR;
			}
		}
	}
	if( mode == ANN_TM_BATCH ){
		if( ANN_WeightUpdate( net, net->dW, 1. ) ){
			ANN_error( ANN_GEN_ERROR, "ANN_trainingEpoch" );
			return ANN_GEN_ERROR;
		}
	}

	return ANN_OK;
}

/* funzione che azzera alcune matrici della rete neurale usate
 * durante l'addestramento */			
ann_res_t ANN_reset( ANN * net){

	int i,j,k;
	
	if( net->r != 0 ) {
		if( vector_null( &net->yD ) ){
			ANN_error( ANN_MATRIX_ERROR, "ANN_reset" );
			return ANN_MATRIX_ERROR;
		}
	}
	for( i=0; i<net->N_layer+1; i++ ){
		if( matrix_null( &net->dW[i] ) ){
			ANN_error( ANN_MATRIX_ERROR, "ANN_reset" );
			return ANN_MATRIX_ERROR;
		}
	}

	for( i=0; i<net->r; i++ ){
		for( j=0; j<net->N_neuron[net->N_layer+1]; j++ ){
			for( k=0; k<net->N_layer+1; k++ ){
				if( matrix_null( &net->dy[i][j][k] ) ){
					ANN_error( ANN_MATRIX_ERROR, "ANN_error" );
					return ANN_MATRIX_ERROR;
				}
			}
		}
	}

	return ANN_OK;
}

/* funzione che calcola l'errore quadratico note la matrici delle uscite
 * della rete e delle uscite desiderate */
ann_res_t ANN_TotalError( matrix *DES_OUTPUT, matrix *NN_OUTPUT, double * err){

	int i,j;

	if( DES_OUTPUT->Nrow != NN_OUTPUT->Nrow || DES_OUTPUT->Ncolumn != NN_OUTPUT->Ncolumn ){
		fprintf( stderr, "Incompatible dimensions\n" );
		ANN_error( ANN_GEN_ERROR, "ANN_TotalError" );
		return ANN_GEN_ERROR;
	}

	*err = 0.;
	for( i=0; i<DES_OUTPUT->Nrow; i++ ){
		for( j=0; j<DES_OUTPUT->Ncolumn; j++ ){
			*err += 0.5*( DES_OUTPUT->mat[i][j]-NN_OUTPUT->mat[i][j] )*( DES_OUTPUT->mat[i][j]-NN_OUTPUT->mat[i][j] );
		}
	}

	return ANN_OK;
}

/* funzione che esegue l'assegnamento tra due variabili di tipo ANN_vector_matrix:
 * vm1 = K*vm2 */ 
ann_res_t ANN_vector_matrix_ass( ANN_vector_matrix *vm1, ANN_vector_matrix *vm2, int* N_neuron, int N_layer, double K ){

	unsigned i,j,k;
	
	for( i=0; i<N_layer+1; i++ ){
		for( j=0; j<N_neuron[i]; j++ ){
			for( k=0; k<N_neuron[i+1]; k++ ){
				//vm1[i]->mat[j][k] = K*vm2[i]->mat[j][k];
				((*vm1)[i]).mat[j][k] = K*((*vm2)[i]).mat[j][k];
				//(&((*vm1)[i]))->mat[j][k] = K*(&((*vm2)[i]))->mat[j][k];
			}
		}
	}
	return ANN_OK;
}


void ANN_error( ann_res_t error, char * string){

        switch(error) {
        case ANN_NO_MEMORY:     fprintf( stderr, "Memory error(@ %s)\n", string );
                                break;
        case ANN_MATRIX_ERROR:  fprintf( stderr, "Error in using matrix library(@ %s)\n", string );
                                break;
        case ANN_NO_FILE:       fprintf( stderr, "Error in file opening(@ %s)\n", string );
                                break;
        case ANN_DATA_ERROR:    fprintf( stderr, "Error in data value(@ %s)\n", string );
                                break;
        case ANN_GEN_ERROR:     fprintf( stderr, "Error(@ %s)\n", string );
                                break;
        default:                break;
        }
}


/* calcola la matrice delle derivate delle uscite rispetto alle
 * derivate degli ingressi */

ann_res_t ANN_jacobian_matrix( ANN *net, matrix *jacobian ){

	unsigned i, j, k;
	
	for( i=0; i<net->N_input; i++ ){
		vector_null( &net->dXdu[0] );
                net->dXdu[0].vec[i] = 1.;
		for( j=0; j<net->N_layer+1; j++ ){
                	if( matrixT_vector_prod( &net->W[j], &net->dXdu[j] ,&net->temp[j+1] ) ){
				ANN_error( ANN_MATRIX_ERROR, "ANN_dXdW" );
				return ANN_MATRIX_ERROR;
			}
                	for( k=0; k<net->N_neuron[j+1];k++ ){
                	 	net->dXdu[j+1].vec[k] = ANN_InternalFunctionDer(net->v[j+1].vec[k], net )*net->temp[j+1].vec[k];
                	}
		}
		for( k=0; k<net->N_output; k++ ){
			//jacobian->mat[i][k] = ( net->output_scale.mat[k][0] * net->input_scale.mat[i][0] )*net->dXdu[net->N_layer+1].vec[k];
			jacobian->mat[i][k] = ( net->input_scale.mat[i][0] / net->output_scale.mat[k][0] )*net->dXdu[net->N_layer+1].vec[k];
		}
	}
	return ANN_OK;
}
			

