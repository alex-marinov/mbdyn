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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ActivationFunction.h"
#include "matrix.h"

/* ACTIVATION FUNCTION 1: hyperbolic tangent */

int
w_tanh_init(void **privp)
{
        *privp = malloc(sizeof(w_tanh_t));

        return (*privp == NULL);
}

int
w_tanh_destroy(void *priv)
{
#if 0
        w_tanh_t        *data = (w_tanh_t *)priv;
#endif

        free(priv);

        return 0;
}

int
w_tanh_read(void *priv, FILE *fh, unsigned flags)
{
        w_tanh_t        *data = (w_tanh_t *)priv;

        fscanf(fh, "%lf", &data->alpha);
        fscanf(fh, "%lf", &data->beta);

        return 0;
}

int
w_tanh_write(void *priv, FILE *fh, unsigned flags)
{
        w_tanh_t        *data = (w_tanh_t *)priv;

        if (flags & W_F_BIN) {
                fprintf(fh, "%d\n%e %e",
                        1,data->alpha, data->beta);

        } else if (flags & W_F_TEXT) {
                fprintf(fh, "tanh alpha=%e beta=%e\n",
                        data->alpha, data->beta);
        }

        return 0;
}

int
w_tanh_eval(void *priv, double in, int order, double *outp)
{
        w_tanh_t        *data = (w_tanh_t *)priv;
        double          y;

        switch (order) {
        case 0:
                y = data->alpha*tanh(data->beta*in);
                break;

        case 1:
                y = tanh(data->beta*in);
                y = data->alpha*data->beta*(1. - y*y);
                break;

        default:
                return 1;
        }

        *outp = y;

        return 0;
}

/* Linear Activation Function */
int
w_linear_init(void **privp)
{
        *privp = malloc(sizeof(w_linear_t));

        return (*privp == NULL);
}

int
w_linear_destroy(void *priv)
{
#if 0
        w_tanh_t        *data = (w_tanh_t *)priv;
#endif

        free(priv);

        return 0;
}

int
w_linear_read(void *priv, FILE *fh, unsigned flags)
{
        w_linear_t        *data = (w_linear_t *)priv;

        fscanf(fh, "%lf", &data->m);
        fscanf(fh, "%lf", &data->q);

        return 0;
}

int
w_linear_write(void *priv, FILE *fh, unsigned flags)
{
        w_linear_t        *data = (w_linear_t *)priv;

        if (flags & W_F_BIN) {
                fprintf(fh, "%d\n%e %e",
                        2,data->m, data->q);

        } else if (flags & W_F_TEXT) {
                fprintf(fh, "linear m=%e q=%e\n",
                        data->m, data->q);
        }

        return 0;
}

int
w_linear_eval(void *priv, double in, int order, double *outp)
{
        w_linear_t        *data = (w_linear_t *)priv;
        double          y;

        switch (order) {
        case 0:
                y = data->m*in + data->q;
                break;

        case 1:
                y = data->m;
                break;

        default:
                return 1;
        }

        *outp = y;

        return 0;
}

