/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2007
 *
 * Marco Morandini	<morandini@aero.polimi.it>
 * Pierangelo Masarati	<masarati@aero.polimi.it>
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
 * AgustaWestland is granted permission to use this file internally
 * in violation of the GNU General Public License (GPL); however,
 * distribution to third-parties must occur according to GPL.
 */

#ifndef NLRHEO_DAMPER_H
#define NLRHEO_DAMPER_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <gsl/gsl_errno.h>
#include "gsl/gsl_interp.h"
#include "gsl/gsl_odeiv.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_linalg.h"

struct sym_params {
	int n_parallelo;
	int *n_serie;			/* [n_parallelo] */
	int max_n_serie;
	double **s_max;			/* [max_n_serie][n_parallelo] */
	double **v_max;			/* [max_n_serie][n_parallelo] */
	int n_elementi;
	int *npti_ks;			/* [n_elementi] */
	int *npti_kv;			/* [n_elementi] */
	int *npti_cs;			/* [n_elementi] */
	int *npti_cv;			/* [n_elementi] */
	int n_variabili_k, n_variabili_c, n_variabili;
	int max_npti_ks, max_npti_kv, max_npti_cs, max_npti_cv;
	double **k_s;			/* [n_elementi][max_npti_ks] */
	double **k_v;			/* [n_elementi][max_npti_kv] */
	double **c_s;			/* [n_elementi][max_npti_cs] */
	double **c_v;			/* [n_elementi][max_npti_cv] */
	gsl_interp ** ik_s, ** ik_v, ** ic_s, ** ic_v;	/* tutti [n_elementi] */
	double *x, *y, *y_dummy;

	/* scale factors */
	double scale_eps;
	double scale_f;
	/* filter angular frequency */
	double a;

	/* per integrazione */
	double sf, si, vf, vi, tf, ti;
	double f, f_s, f_v;

	double F, FDE, FDEPrime;

	/* di servizio */
	gsl_matrix **gsl_C, **gsl_K;
	gsl_vector **gsl_xp, **gsl_x, **gsl_b;
	gsl_permutation **gsl_perm;

	const gsl_odeiv_step_type *T;
	gsl_odeiv_step *stepint;
	gsl_odeiv_evolve * evolve;
	gsl_odeiv_control * control;
	gsl_odeiv_system sys;
	double prev_time, current_time, dt;
	double prev_eps, prev_epsPrime;
};

/* read functions */
extern int nlrheo_get_int(int *i);
extern int nlrheo_get_real(double *d);

extern "C" int nlrheo_init(void *v_nlrheo);
extern "C" int nlrheo_destroy(void *v_nlrheo);
extern "C" int nlrheo_update(void *v_nlrheo,
	double t_curr, double eps, double epsPrime, int do_try);
extern "C" int nlrheo_parse(void *v_nlrheo,
	double scale_eps, double scale_f, double filter);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* ! NLRHEO_DAMPER_H */

