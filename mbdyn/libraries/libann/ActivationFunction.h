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

/* ACTIVATION FUNCTION */
#define W_F_NONE        (0x00U)
#define W_F_TEXT        (0x01U)
#define W_F_BIN         (0x02U)
/* activation interface */
typedef int (*w_init_f)(void **);
typedef int (*w_destroy_f)(void *);
typedef int (*w_read_f)(void *, FILE *fd, unsigned flags);
typedef int (*w_write_f)(void *, FILE *fd, unsigned flags);
typedef int (*w_eval_f)(void *, double, int, double *);

/* tanh activation function */
typedef struct w_tanh_t {
        double alpha;
        double beta;
} w_tanh_t;

int w_tanh_init(void ** );
int w_tanh_destroy(void *);
int w_tanh_read(void * , FILE * , unsigned );
int w_tanh_write(void * , FILE * , unsigned );
int w_tanh_eval(void * , double, int, double *);

typedef struct w_linear_t {
        double m;
        double q;
} w_linear_t;

int w_linear_init(void ** );
int w_linear_destroy(void *);
int w_linear_read(void * , FILE * , unsigned );
int w_linear_write(void * , FILE * , unsigned );
int w_linear_eval(void * , double, int, double *);
