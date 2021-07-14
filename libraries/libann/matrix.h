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

#ifndef ANN_MATRIX_H
#define ANN_MATRIX_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdio.h>

#define W_M_NONE        (0x00U)
#define W_M_TEXT        (0x01U)
#define W_M_BIN         (0x02U)
#define W_M_BIN_ROW     (0x04U)

/* codici errore */
typedef enum {
        MAT_OK = 0,
	MAT_DIMENSION,
        MAT_NO_MEMORY,
	MAT_INPUT,
	MAT_GEN_ERROR
} mat_res_t;

/* elemento classe matrice */
typedef struct matrix {
	double ** mat;
	unsigned Nrow;
	unsigned Ncolumn;
} matrix;

/* elemento classe vettore */
typedef struct vector {
	double * vec;
	unsigned dimension;
} vector;

/* prototipi funzioni */

mat_res_t matrix_init( matrix *, unsigned , unsigned );
mat_res_t vector_init( vector *, unsigned );
mat_res_t matrix_destroy( matrix * ); 
mat_res_t vector_destroy( vector * );

mat_res_t matrix_null( matrix * );
mat_res_t matrix_eye( matrix *, double );
mat_res_t matrix_copy( matrix *, matrix *, double );
mat_res_t vector_copy( vector *, vector *, double );
mat_res_t vector_null( vector * );
mat_res_t matrix_prod( matrix * , matrix *, matrix *, double );
mat_res_t matrix_prod_sym( matrix * , matrix *, matrix *, double );
mat_res_t matrix_prod_transpose( matrix * , matrix *, matrix *, double );
mat_res_t matrix_transpose_prod( matrix * , matrix *, matrix *, double );
mat_res_t matrix_vector_prod( matrix *, vector *, vector *);
mat_res_t matrixT_vector_prod( matrix *, vector *, vector *);
mat_res_t matrix_sum( matrix *, matrix *, matrix *, double );
mat_res_t matrix_sum_transpose( matrix * , matrix *, matrix *, double );
mat_res_t vector_sum( vector *, vector *, vector *, double );
mat_res_t matrix_transpose( matrix *, matrix *);
mat_res_t scalar_prod( vector *, vector *, double *);
mat_res_t vector_vector_prod( vector *, vector *, matrix *, double);

mat_res_t matrix_write( matrix *, FILE *, unsigned );
mat_res_t vector_write( vector *, FILE *, unsigned );
mat_res_t matrix_read( matrix *, FILE *, unsigned );
mat_res_t vector_read( vector *, FILE *, unsigned );

mat_res_t matrix_random( matrix *, double, double );
mat_res_t vector_random( vector *, double, double );

mat_res_t sub_matrix_extract( matrix *, matrix *, unsigned, unsigned);
mat_res_t sub_matrix_insert( matrix *, matrix *, unsigned, unsigned);

double mean_value( matrix *, int );
double variance( matrix *, int );
double maximum( matrix *, int );
double minimum( matrix *, int  );
double matrix_trace( matrix * );
void matrix_error( mat_res_t, const char * );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* ANN_MATRIX_H */

