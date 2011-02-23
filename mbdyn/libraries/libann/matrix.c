/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2008
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "matrix.h"

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

/* azzera una matrice MAT = zeros*/
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

/* Identity matrix: MAT = K*eye */
mat_res_t matrix_eye( matrix *MAT, double K ){

	unsigned i;

	/* azzero la matrice del risultato */
	if( matrix_null(MAT) != MAT_OK ){
		matrix_error( MAT_GEN_ERROR, "matrix_eye" );
		return MAT_GEN_ERROR;
	}

	for( i=0; i<MAT->Nrow; i++ ){
		MAT->mat[i][i] = K;
	}

	return MAT_OK;
}

/* copia una matrice MAT1 = MAT2*/
mat_res_t matrix_copy( matrix *MAT1, matrix *MAT2, double K ){

	unsigned i, j;

	/* controllo dimensionale */
	if( MAT1->Ncolumn != MAT2->Ncolumn || MAT1->Nrow != MAT2->Nrow){
		matrix_error( MAT_DIMENSION, "matrix_copy" );
		return MAT_DIMENSION;
	}

	for( i=0; i<MAT1->Nrow; i++ ){
		for( j=0; j<MAT1->Ncolumn; j++ ){
			MAT1->mat[i][j] = K*MAT2->mat[i][j];
		}
	}

	return MAT_OK;
}
/* copia un vettore VEC1 = VEC2*/
mat_res_t vector_copy( vector *VEC1, vector *VEC2, double K ){

	unsigned i;

	/* controllo dimensionale */
	if( VEC1->dimension != VEC2->dimension){
		matrix_error( MAT_DIMENSION, "vector_copy" );
		return MAT_DIMENSION;
	}

	for( i=0; i<VEC1->dimension; i++ ){
			VEC1->vec[i] = K*VEC2->vec[i];
	}

	return MAT_OK;
}
/* azzera un vettore VEC = zeros */
mat_res_t vector_null( vector *VEC ){

	if( !(memset( VEC->vec, 0, VEC->dimension*sizeof(double) ) ) ){
		matrix_error( MAT_GEN_ERROR, "vector_null" );
		return MAT_GEN_ERROR;
	}	       

	return MAT_OK;
}

/* prodotto tra matrici MAT_R = K*MAT1*MAT2*/
mat_res_t matrix_prod( matrix *MAT1 ,matrix *MAT2, matrix *MAT_R, double K ){

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
				MAT_R->mat[i][j] += K*MAT1->mat[i][k]*MAT2->mat[k][j];
			}
		}
	}

	return MAT_OK;
}
/* prodotto tra matrici MAT_R = K*MAT1*MAT2*/
mat_res_t matrix_prod_sym( matrix *MAT1 ,matrix *MAT2, matrix *MAT_R, double K ){

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
		for( j=i; j<MAT2->Ncolumn; j++ ){
			for( k=0; k<MAT1->Ncolumn; k++ ){
				MAT_R->mat[i][j] += K*MAT1->mat[i][k]*MAT2->mat[k][j];
			}
			MAT_R->mat[j][i] = MAT_R->mat[i][j];
		}
	}

	return MAT_OK;
}

/* matrice trasposta MAT1 = MAT2^T */
mat_res_t matrix_transpose( matrix *MAT1 ,matrix *MAT2){

	unsigned i,j;

	/* controllo dimensionale */
	if( MAT1->Ncolumn != MAT2->Nrow || MAT1->Nrow != MAT2->Ncolumn ){
		matrix_error( MAT_DIMENSION, "matrix_transpose" );
		return MAT_DIMENSION;
	}

	for( i=0; i<MAT2->Nrow; i++ ){
		for( j=0; j<MAT2->Ncolumn; j++ ){
			MAT1->mat[j][i] = MAT2->mat[i][j];
		}
	}

	return MAT_OK;
}
/* prodotto tra matrice e matrice trasposta MAT_R =K*MAT1^T*MAT2*/
mat_res_t matrix_transpose_prod( matrix *MAT1 ,matrix *MAT2, matrix *MAT_R, double K ){

	unsigned i,j,k;
	
	/* controllo dimensionale */
	if( MAT1->Nrow != MAT2->Nrow || MAT_R->Nrow != MAT1->Ncolumn || MAT_R->Ncolumn != MAT2->Ncolumn  ){
		matrix_error( MAT_DIMENSION, "matrix_transpose_prod" );
		return MAT_DIMENSION;
	}
	/* azzero la matrice del risultato */
	if( matrix_null(MAT_R) != MAT_OK ){
		matrix_error( MAT_GEN_ERROR, "matrix_transpose_prod" );
		return MAT_GEN_ERROR;
	}
	for( i=0; i<MAT1->Ncolumn; i++ ){
		for( j=0; j<MAT2->Ncolumn; j++ ){
			for( k=0; k<MAT1->Nrow; k++ ){
				MAT_R->mat[i][j] += K*MAT1->mat[k][i]*MAT2->mat[k][j];
			}
		}
	}

	return MAT_OK;
}

/* prodotto tra matrice e matrice trasposta MAT_R =K*MAT1*MAT2^T*/
mat_res_t matrix_prod_transpose( matrix *MAT1 ,matrix *MAT2, matrix *MAT_R, double K ){

	unsigned i,j,k;
	
	/* controllo dimensionale */
	if( MAT1->Ncolumn != MAT2->Ncolumn || MAT_R->Nrow != MAT1->Nrow || MAT_R->Ncolumn != MAT2->Nrow  ){
		matrix_error( MAT_DIMENSION, "matrix_prod_transpose" );
		return MAT_DIMENSION;
	}
	/* azzero la matrice del risultato */
	if( matrix_null(MAT_R) != MAT_OK ){
		matrix_error( MAT_GEN_ERROR, "matrix_prod_transpose" );
		return MAT_GEN_ERROR;
	}
	for( i=0; i<MAT1->Nrow; i++ ){
		for( j=0; j<MAT2->Nrow; j++ ){
			for( k=0; k<MAT1->Ncolumn; k++ ){
				MAT_R->mat[i][j] += K*MAT1->mat[i][k]*MAT2->mat[j][k];
			}
		}
	}

	return MAT_OK;
}

/* prodotto scalare RES = VEC1 dot VEC2*/
mat_res_t scalar_prod( vector *VEC1, vector *VEC2, double *RES){

	unsigned i;
	double res;
	
	/* controllodimensionale */
	if( VEC1->dimension != VEC2->dimension ){
		matrix_error( MAT_DIMENSION, "matrix_vector_prod" );
		return MAT_DIMENSION;
	}

	res = 0.;
	for( i=0; i<VEC1->dimension; i++){
		res += VEC1->vec[i]*VEC2->vec[i];
	}

	*RES = res;

	return MAT_OK;

}

/* prodotto vettore*vettore = matrice RES = K*VEC1*VEC2^T */
mat_res_t vector_vector_prod( vector *VEC1, vector *VEC2, matrix *RES, double K){

	unsigned i, j;
	
	/* controllodimensionale */
	if( VEC1->dimension != RES->Nrow || VEC2->dimension != RES->Ncolumn){
		matrix_error( MAT_DIMENSION, "vector_vector_prod" );
		return MAT_DIMENSION;
	}

	for( i=0; i<VEC1->dimension; i++){
		for( j=0; j<VEC2->dimension; j++){
			RES->mat[i][j] = K*VEC1->vec[i]*VEC2->vec[j];
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
	
/* prodotto tra una matrice trasposta ed un vettore VEC_R = MAT^T*VEC */
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

	//for( i=0; i<MAT->Nrow; i++ ){
	//	for( j=0; j<MAT->Ncolumn; j++ ){
	//		VEC_R->vec[j] += MAT->mat[i][j]*VEC->vec[i];
	//	}
	//}
	for( i=0; i<MAT->Ncolumn; i++ ){
		for( j=0; j<MAT->Nrow; j++ ){
			VEC_R->vec[i] += MAT->mat[j][i]*VEC->vec[j];
		}
	}
	
	return MAT_OK;
}

/* somma tra matrici MAT_R = MAT1+K*MAT */
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

/* somma tra matrici MAT_R = MAT1+K*MAT2^T */
mat_res_t matrix_sum_transpose( matrix *MAT1, matrix *MAT2, matrix *MAT_R, double K ){

	unsigned i,j;

	/* controllo dimensionale */
	if( MAT1->Ncolumn != MAT2->Nrow || MAT_R->Ncolumn != MAT1->Ncolumn || MAT_R->Nrow != MAT1->Nrow || MAT1->Nrow != MAT2->Ncolumn ){
		matrix_error( MAT_DIMENSION, "matrix_sum_transpose" );
		return MAT_DIMENSION;
	}
	for( i=0; i<MAT1->Nrow; i++ ){
		for( j=0; j<MAT1->Ncolumn; j++ ){
			MAT_R->mat[i][j] = MAT1->mat[i][j] + K*MAT2->mat[j][i];
		}
	}

	return MAT_OK;
}

/* somma tra vettori VEC_R = VEC1+K*VEC2 */
mat_res_t vector_sum( vector *VEC1, vector *VEC2, vector *VEC_R, double K ){

	unsigned i;

	/* controllo dimensionale */
	if( VEC1->dimension != VEC2->dimension || VEC_R->dimension != VEC1->dimension ){
		matrix_error( MAT_DIMENSION, "vector_sum" );
		return MAT_DIMENSION;
	}
	for( i=0; i<VEC1->dimension; i++ ){
		VEC_R->vec[i] = VEC1->vec[i] + K*VEC2->vec[i];
	}

	return MAT_OK;
}


/* FUNZIONI ACCESSORIE */

/* scrive a video o su file una matrice */
mat_res_t matrix_write( matrix *MAT, FILE *fh, unsigned flags ){

	unsigned i,j;

	if( flags & W_M_TEXT )		fprintf( fh, "matrix = [\n" );
	if( flags & W_M_BIN )		fprintf( fh, "\n" );
	for( i=0; i<MAT->Nrow; i++ ){
		for( j=0; j<MAT->Ncolumn; j++ ){
			fprintf( fh, "%15.16e ", MAT->mat[i][j] );
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
	/*if( flags & W_M_BIN )		fprintf( fh, "\n" );*/
	for( i=0; i<VEC->dimension; i++ ){
		if( flags & W_M_BIN )		fprintf( fh, "\n" );
		fprintf( fh, "%15.16e ", VEC->vec[i] );
	}

	if( flags & W_M_TEXT )		fprintf( fh, "]\n" );
	if( flags & W_M_BIN )		fprintf( fh, "\n" );
	if( flags & W_M_BIN_ROW )	fprintf( fh, "\n" );
	
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
void matrix_error( mat_res_t error, const char * string){

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

mat_res_t vector_random( vector *VEC, double min, double max ){

	double y;
	unsigned i;

	for( i=0; i<VEC->dimension; i++ ){
	       	y = rand();
		y = y/RAND_MAX;
		y = y*(max-min);
        	y += min ;
		VEC->vec[i] = y;
	}
	return MAT_OK;
}

/* calcola il valor medio della colonna "column" della matrice "MAT" */
double mean_value( matrix *MAT, int column ){

        unsigned i;
        double mean = 0.;

        for( i=0; i<MAT->Nrow; i++ ){
                mean += MAT->mat[i][column];
        }
        return( mean/MAT->Nrow );
}

/* calcola la varianza della colonna "column" della matrice "MAT" */
double variance( matrix *MAT, int column ){

        unsigned i;
        double mean, var = 0.;

        mean = mean_value( MAT, column);
        for( i=0; i<MAT->Nrow; i++ ){
                var += (MAT->mat[i][column]-mean)*(MAT->mat[i][column]-mean);
        }

        return( var/MAT->Nrow );
}

/* calcola il valor massimo della colonna "column" della matrice "MAT" */
double maximum( matrix *MAT, int column ){

        unsigned i;
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

        unsigned i;
        double MIN;

        MIN = MAT->mat[0][column];
        for( i=0; i<MAT->Nrow; i++ ){
                if( MAT->mat[i][column] < MIN ){
                        MIN = MAT->mat[i][column];
                }
        }

        return MIN; 
}

/* calcola la traccia di una matrice */
double matrix_trace( matrix *MAT ){

        unsigned i;
        double trace;

        trace = 0.;
        for( i=0; i<MAT->Nrow; i++ ){
		trace += MAT->mat[i][i];
        }

        return trace; 
}
mat_res_t sub_matrix_extract( matrix *BIG, matrix *SUB, unsigned RowIndex, unsigned ColumnIndex){

	unsigned i, j;

	if( (RowIndex+SUB->Nrow > BIG->Nrow) || (ColumnIndex+SUB->Ncolumn > BIG->Ncolumn) ){
		matrix_error( MAT_DIMENSION, "sub_matrix_extract" );
		return MAT_DIMENSION;
	}

	for( i=0; i<SUB->Nrow; i++){
		for( j=0; j<SUB->Ncolumn; j++){
			SUB->mat[i][j] = BIG->mat[RowIndex+i][ColumnIndex+j];
		}
	}

	return MAT_OK;
}

mat_res_t sub_matrix_insert( matrix *BIG, matrix *SUB, unsigned RowIndex, unsigned ColumnIndex){

	unsigned i, j;

	if( (RowIndex+SUB->Nrow > BIG->Nrow) || (ColumnIndex+SUB->Ncolumn > BIG->Ncolumn) ){
		matrix_error( MAT_DIMENSION, "sub_matrix_insert" );
		return MAT_DIMENSION;
	}

	for( i=0; i<SUB->Nrow; i++){
		for( j=0; j<SUB->Ncolumn; j++){
			BIG->mat[RowIndex+i][ColumnIndex+j] = SUB->mat[i][j];
		}
	}

	return MAT_OK;

}
