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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <matrix.h>

/* inizializza un elemento della classe matrice */
mat_res_t matrix_init( matrix *MAT, unsigned Nrow, unsigned Ncolumn ){

	unsigned i;

	if( Nrow<= 0 || Ncolumn <= 0 ){
		matrix_error( MAT_DIMENSION, "matrix_init" );
		return MAT_DIMENSION;
	}
	
	MAT->Nrow = Nrow;
	MAT->Ncolumn = Ncolumn;
	if( !( MAT->mat = (double **)calloc( Nrow, sizeof(double *) ) ) ){
		matrix_error( MAT_NO_MEMORY, "matrix_init" );
		return MAT_NO_MEMORY;
	}
	for( i=0; i<Nrow; i++ ){
		if( !( MAT->mat[i] = (double *)calloc( Ncolumn, sizeof(double) ) ) ){
			matrix_error( MAT_NO_MEMORY, "matrix_init" );
			return MAT_NO_MEMORY;
		}
	}

	return MAT_OK;
}

/* inizializza un elemento della classe vettore */
mat_res_t vector_init( vector *VEC, unsigned dimension ){
	if( dimension<= 0 ){
		matrix_error( MAT_DIMENSION, "vector_init" );
		return MAT_DIMENSION;
	}
	
	VEC->dimension = dimension;

	if( !( VEC->vec = (double *)calloc( dimension, sizeof(double) ) ) ){
		matrix_error( MAT_NO_MEMORY, "vector_init");
		return MAT_NO_MEMORY;
	}

	return MAT_OK;
}

/* distrugge un elemento della classe matrice*/
mat_res_t matrix_destroy( matrix *MAT ) {

	unsigned i;

	for( i=0; i<MAT->Nrow; i++ ){
		free(MAT->mat[i]);
	}
	free(MAT->mat);

	return MAT_OK;
}

/* distrugge un elemento della classe vettore */
mat_res_t vector_destroy( vector *VEC ) {

	free(VEC->vec);

	return MAT_OK;
}

/* OPERAZIONI TRA MATRICI E VETTORI */

/* azzera una matrice */
mat_res_t matrix_null( matrix *MAT ){

	unsigned i;

	for( i=0; i<MAT->Nrow; i++ ){
		if( !(memset( MAT->mat[i], 0, MAT->Ncolumn*sizeof(double) ) ) ){
			matrix_error( MAT_GEN_ERROR, "matrix_null" );
			return MAT_GEN_ERROR;
		}	       
	}

	return MAT_OK;
}

/* azzera un vettore */
mat_res_t vector_null( vector *VEC ){

	if( !(memset( VEC->vec, 0, VEC->dimension*sizeof(double) ) ) ){
		matrix_error( MAT_GEN_ERROR, "vector_null" );
		return MAT_GEN_ERROR;
	}	       

	return MAT_OK;
}

/* prodotto tra matrici MAT_R = MAT1*MAT2*/
mat_res_t matrix_prod( matrix *MAT1 ,matrix *MAT2, matrix *MAT_R ){

	unsigned i,j,k;
	
	/* controllo dimensionale */
	if( MAT1->Ncolumn != MAT2->Nrow || MAT_R->Nrow != MAT1->Nrow || MAT_R->Ncolumn != MAT2->Ncolumn  ){
		matrix_error( MAT_DIMENSION, "matrix_prod" );
		return MAT_DIMENSION;
	}
	/* azzero la matrice del risultato */
	if( matrix_null(MAT_R) != MAT_OK ){
		matrix_error( MAT_GEN_ERROR, "matrix_prod" );
		return MAT_GEN_ERROR;
	}
	for( i=0; i<MAT1->Nrow; i++ ){
		for( j=0; j<MAT2->Ncolumn; j++ ){
			for( k=0; k<MAT1->Ncolumn; k++ ){
				MAT_R->mat[i][j] += MAT1->mat[i][k]*MAT2->mat[k][j];
			}
		}
	}

	return MAT_OK;
}

/* prodotto matrice vettore VEC_R = MAT*VEC */
mat_res_t matrix_vector_prod( matrix *MAT, vector *VEC, vector *VEC_R){

	unsigned i,j;

	/* controllodimensionale */
	if( MAT->Ncolumn != VEC->dimension || MAT->Nrow != VEC_R->dimension ){
		matrix_error( MAT_DIMENSION, "matrix_vector_prod" );
		return MAT_DIMENSION;
	}
	/* azzero il vettore risultato */
	if( vector_null(VEC_R) != MAT_OK ){
		matrix_error( MAT_GEN_ERROR, "matrix_vector_prod");
		return MAT_GEN_ERROR;
	}

	for( i=0; i<MAT->Nrow; i++ ){
		for( j=0; j<MAT->Ncolumn; j++ ){
			VEC_R->vec[i] += MAT->mat[i][j]*VEC->vec[j];
		}
	}
	
	return MAT_OK;
}	

/* prodotto tra una matrice trasposta ed un vettore VEC_R = MAT'*VEC */
mat_res_t matrixT_vector_prod( matrix *MAT, vector *VEC, vector *VEC_R){

	unsigned i,j;
	
	/* controllo dimensionale */
	if( MAT->Nrow != VEC->dimension || MAT->Ncolumn != VEC_R->dimension ){
		matrix_error( MAT_DIMENSION, "matrixT_vector_prod" );
		return MAT_DIMENSION;
	}
	/* annullo il vettore risulatante */
	if( vector_null(VEC_R) != MAT_OK ){
		matrix_error( MAT_GEN_ERROR, "matrixT_vector_prod" );
		return MAT_GEN_ERROR;
	}

	for( i=0; i<MAT->Nrow; i++ ){
		for( j=0; j<MAT->Ncolumn; j++ ){
			VEC_R->vec[j] += MAT->mat[i][j]*VEC->vec[i];
		}
	}
	
	return MAT_OK;
}

/* somma tra matrici MAT_R = MAT1+K*MAT" */
mat_res_t matrix_sum( matrix *MAT1, matrix *MAT2, matrix *MAT_R, double K ){

	unsigned i,j;

	/* controllo dimensionale */
	if( MAT1->Ncolumn != MAT2->Ncolumn || MAT_R->Ncolumn != MAT1->Ncolumn || MAT_R->Nrow != MAT2->Nrow || MAT1->Nrow != MAT2->Nrow ){
		matrix_error( MAT_DIMENSION, "matrix_sum" );
		return MAT_DIMENSION;
	}
	for( i=0; i<MAT1->Nrow; i++ ){
		for( j=0; j<MAT1->Ncolumn; j++ ){
			MAT_R->mat[i][j] = MAT1->mat[i][j] + K*MAT2->mat[i][j];
		}
	}

	return MAT_OK;
}

/* FUNZIONI ACCESSORIE */

/* scrive a video o du file una matrice */
mat_res_t matrix_write( matrix *MAT, FILE *fh, unsigned flags ){

	unsigned i,j;

	if( flags & W_M_TEXT )		fprintf( fh, "matrix = [\n" );
	if( flags & W_M_BIN )		fprintf( fh, "\n" );
	for( i=0; i<MAT->Nrow; i++ ){
		for( j=0; j<MAT->Ncolumn; j++ ){
			fprintf( fh, "%e ", MAT->mat[i][j] );
		}	
		fprintf( fh, "\n");
	}

	if( flags & W_M_TEXT )		fprintf( fh, "]\n" );
	if( flags & W_M_BIN )		fprintf( fh, "\n" );
	
	return MAT_OK;
}

/* scrive a video o su file un vettore */
mat_res_t vector_write( vector *VEC, FILE *fh, unsigned flags ){

	unsigned i;

	if( flags & W_M_TEXT )		fprintf( fh, "vector = [\n" );
	if( flags & W_M_BIN )		fprintf( fh, "\n" );
	for( i=0; i<VEC->dimension; i++ ){
		fprintf( fh, "%e\n", VEC->vec[i] );
	}

	if( flags & W_M_TEXT )		fprintf( fh, "]\n" );
	if( flags & W_M_BIN )		fprintf( fh, "\n" );
	
	return MAT_OK;
}

/* legge una matrice da file */
mat_res_t matrix_read( matrix *MAT, FILE * fh, unsigned flags ){

	unsigned i,j;

	for( i=0; i<MAT->Nrow; i++ ){
		for( j=0; j<MAT->Ncolumn; j++ ){
			if( fscanf( fh, "%le", &MAT->mat[i][j] ) <= 0 ){
				return MAT_INPUT;
			}
		}
	}		
		
	return MAT_OK;
}

/* legge una vettore da file */
mat_res_t vector_read( vector *VEC, FILE * fh, unsigned flags ){

	unsigned i;

	for( i=0; i<VEC->dimension; i++ ){
		if( fscanf( fh, "%le", &VEC->vec[i] ) <= 0 ){
			return MAT_INPUT;
		}
	}		
		
	return MAT_OK;
}
/* gestione degli errori */
void matrix_error( mat_res_t error, char * string){

	switch(error) {
	case MAT_NO_MEMORY: 	fprintf( stderr, "Memory error( @ %s )\n", string );
			   	break;
	case MAT_DIMENSION:	fprintf( stderr, "Matrix dimension mismatch( @ %s )\n", string );
				break;
	case MAT_INPUT:		fprintf( stderr, "Reading error( @ %s )\n",string );
				break;
	case MAT_GEN_ERROR:	fprintf( stderr, "Error( @ %s )\n",string );
				break;
	default:		break;
	}
}

/* genera una matrice random di numeri compresi tra min e max*/

mat_res_t matrix_random( matrix *MAT, double min, double max ){

	double y;
	unsigned i,j;

	for( i=0; i<MAT->Nrow; i++ ){
		for( j=0; j<MAT->Ncolumn; j++ ){
	        	y = rand();
			y = y/RAND_MAX;
			y = y*(max-min);
        		y += min ;
			MAT->mat[i][j] = y;
		}
	}

	return MAT_OK;
}

/* calcola il valor medio della colonna "column" della matrice "MAT" */
double mean_value( matrix *MAT, int column ){

        int i;
        double mean = 0.;

        for( i=0; i<MAT->Nrow; i++ ){
                mean += MAT->mat[i][column];
        }
        return( mean/MAT->Nrow );
}

/* calcola la varianza della colonna "column" della matrice "MAT" */
double variance( matrix *MAT, int column ){

        int i;
        double mean, var = 0.;

        mean = mean_value( MAT, column);
        for( i=0; i<MAT->Nrow; i++ ){
                var += (MAT->mat[i][column]-mean)*(MAT->mat[i][column]-mean);
        }

        return( var/MAT->Nrow );
}

/* calcola il valor massimo della colonna "column" della matrice "MAT" */
double maximum( matrix *MAT, int column ){

        int i;
        double MAX;
        
        MAX = MAT->mat[0][column];
        for( i=0; i<MAT->Nrow; i++ ){
                if( MAT->mat[i][column] > MAX ){
                        MAX = MAT->mat[i][column];
                }
        }

        return MAX; 
}

/* calcola il valor minimo della colonna "column" della matrice "MAT" */
double minimum( matrix *MAT, int column ){

        int i;
        double MIN;

        MIN = MAT->mat[0][column];
        for( i=0; i<MAT->Nrow; i++ ){
                if( MAT->mat[i][column] < MIN ){
                        MIN = MAT->mat[i][column];
                }
        }

        return MIN; 
}

