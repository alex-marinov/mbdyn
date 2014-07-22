/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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


#include "nlrheo_damper.h"

#include <math.h>
#include <float.h>

static void
nlrheo_int_compute_kc(double *k, double *c, 
	double * k_partial_s,
	double * k_partial_v,
	double * c_partial_s,
	double * c_partial_v,
	const double s, 
	const double v,
	const int el, 
	const int n_variabili_k,
	const int *const npti_ks,
	const int *const npti_kv,
	const int *const npti_cs,
	const int *const npti_cv,
	const gsl_interp *const ik_s, 
	const gsl_interp *const ik_v, 
	const gsl_interp *const ic_s,
	const gsl_interp *const ic_v,
	const double *const k_s,
	const double *const k_v,
	const double *const c_s,
	const double *const c_v,
	const double *const x)
{
	int i, ii, e;

	double k_v_i[npti_kv[el]];
	double k_v_i_partial_s[npti_kv[el]];
	int somma_k = 0;

	double c_v_i[npti_cv[el]];
	double c_v_i_partial_s[npti_cv[el]];
	int somma_c = n_variabili_k;

	/* calcola k di el dato s,v */
	for (e = 0; e < el; e++) {
		somma_k += npti_ks[e] * npti_kv[e];
	}
	if (ik_s != 0) {
		for (i = 0; i < npti_kv[el]; i++) {
			if (fabs(s) > k_s[npti_ks[el]-1]) {
				k_v_i[i] = x[somma_k + npti_ks[el] - 1];
#if 0
					(x[somma_k + npti_ks[el] - 1] - x[somma_k + npti_ks[el] - 2])
					/ (k_s[npti_ks[el]-1] - k_s[npti_ks[el]-2])
					* (fabs(s) - k_s[npti_ks[el]-1])
					+ x[somma_k + npti_ks[el] - 1];
#endif
				k_v_i_partial_s[i] = 0.;
			} else if (fabs(s) == 0) {
				k_v_i[i] = x[somma_k];
				k_v_i_partial_s[i] = 0.;
			} else {
				double tx[npti_ks[el]];
				for (ii = 0; ii < npti_ks[el]; ii++) {
					tx[ii] = x[somma_k + ii];
				}
				k_v_i[i] = gsl_interp_eval(ik_s, k_s, tx, fabs(s), 0);
				k_v_i_partial_s[i] = gsl_interp_eval_deriv(ik_s, k_s, tx, fabs(s), 0) 
					* copysign(1., s);
			}
			somma_k += npti_ks[el];
		}
	} else {
		for (i = 0; i < npti_kv[el]; i++) {
			k_v_i[i] = x[somma_k + i];
			k_v_i_partial_s[i] = 0.;
		}
	}
	if (ik_v != 0) {
		if (fabs(v) > k_v[npti_kv[el] - 1]) {
			*k = k_v_i[npti_kv[el] - 1];
#if 0
				(k_v_i[npti_kv[el]-1] - k_v_i[npti_kv[el]-2]) 
				/ (k_v[npti_kv[el] - 1] - k_v[npti_kv[el] - 2])
				* (fabs(v) - k_v[npti_kv[el] - 1])
				+ k_v_i[npti_kv[el]-1];
#endif
			*k_partial_s = k_v_i_partial_s[npti_kv[el] - 1];
			*k_partial_v = 0.;
		} else if (fabs(v) == 0) {
			*k = k_v_i[0];
			*k_partial_s = k_v_i_partial_s[0];
			*k_partial_v = 0.;
		} else {
			*k = gsl_interp_eval(ik_v, k_v, k_v_i, fabs(v), 0);
			*k_partial_s = gsl_interp_eval(ik_v, k_v, k_v_i_partial_s, fabs(v), 0);
			*k_partial_v = gsl_interp_eval_deriv(ik_v, k_v, k_v_i, fabs(v), 0)
				* copysign(1., v);
		}
	} else {
		*k = k_v_i[0];
		*k_partial_s = k_v_i_partial_s[0];
		*k_partial_v = 0.;
	}
	
	/* calcola c di el dato s,v */
	for (e = 0; e < el; e++) {
		somma_c += npti_cs[e] * npti_cv[e];
	}
	if (ic_s != 0) {
		for (i = 0; i < npti_cv[el]; i++) {
			if (fabs(s) > c_s[npti_cs[el] - 1]) {
				c_v_i[i] = x[somma_c + npti_cs[el] - 1];
#if 0
					(x[somma_c + npti_cs[el] - 1] - x[somma_c + npti_cs[el] - 2])
					/ (c_s[npti_cs[el] - 1] - c_s[npti_cs[el] - 2])
					* (fabs(s) - c_s[npti_cs[el] - 1])
					+ x[somma_c + npti_cs[el] - 1];
#endif
				c_v_i_partial_s[i] = 0.;
			} else if (fabs(s) == 0) {
				c_v_i[i] = x[somma_c];
				c_v_i_partial_s[i] = 0.;
			} else {
				double tx[npti_cs[el]];
				for (ii = 0; ii < npti_cs[el]; ii++) {
					tx[ii] = x[somma_c + ii];
				}
				c_v_i[i] = gsl_interp_eval(ic_s, c_s, tx, fabs(s), 0);
				c_v_i_partial_s[i] = gsl_interp_eval_deriv(ic_s, c_s, tx, fabs(s), 0)
					* copysign(1., s);
			}
			somma_c += npti_cs[el];
		}
	} else {
		for (i = 0; i < npti_cv[el]; i++) {
			c_v_i[i] = x[somma_c + i];
			c_v_i_partial_s[i] = 0.;
		}
	}
	if (ic_v != 0) {
		if (fabs(v) > c_v[npti_cv[el] - 1]) {
			*c = c_v_i[npti_cv[el] - 1];
#if 0
				(c_v_i[npti_cv[el]-1] - c_v_i[npti_cv[el]-2])
				/ (c_v[npti_cv[el] - 1] - c_v[npti_cv[el] - 2])
				* (fabs(v) - c_v[npti_cv[el] - 1]) +
				+ c_v_i[npti_cv[el]-1];
#endif
			*c_partial_s = c_v_i_partial_s[npti_cv[el] - 1];
			*c_partial_v = 0.;
		} else if (fabs(v) == 0) {
			*c = c_v_i[0];
			*c_partial_s = c_v_i_partial_s[0];
			*c_partial_v = 0.;
		} else {
			*c = gsl_interp_eval(ic_v, c_v, c_v_i, fabs(v), 0);
			*c_partial_s = gsl_interp_eval(ic_v, c_v, c_v_i_partial_s, fabs(v), 0);
			*c_partial_v = gsl_interp_eval_deriv(ic_v, c_v, c_v_i, fabs(v), 0)
				* copysign(1., v);
		}
	} else {
		*c = c_v_i[0];
		*c_partial_s = c_v_i_partial_s[0];
		*c_partial_v = 0.;
	}
#if 0
 	std::cerr << fabs(s) << " " << fabs(v) << " " << k << " " << c << std::endl;
#endif
}

/* from the reference manual of gsl:

25.1 Defining the ODE System
============================

The routines solve the general n-dimensional first-order system,

     dy_i(t)/dt = f_i(t, y_1(t), ..., y_n(t))

for i = 1, \dots, n.  The stepping functions rely on the vector of
derivatives f_i and the Jacobian matrix, J_{ij} = df_i(t,y(t)) / dy_j.
A system of equations is defined using the `gsl_odeiv_system' datatype.

 -- Data Type: gsl_odeiv_system
     This data type defines a general ODE system with arbitrary
     parameters.

    `int (* function) (double t, const double y[], double dydt[], void * params)'
          This function should store the vector elements
          f_i(t,y,params) in the array DYDT, for arguments (T,Y) and
          parameters PARAMS.  The function should return `GSL_SUCCESS'
          if the calculation was completed successfully. Any other
          return value indicates an error.
*/

int
nlrheo_int_func(double t, const double y[], double f[], void *para)
{
	sym_params * pa = (sym_params *)para;
	int displ_hi_freq_low_pass_filter = 1; /* true */
	int vel_hi_freq_low_pass_filter = 1;   /* true */
	double sstatic = y[pa->n_elementi - pa->n_parallelo + 2];
	double vstatic = y[pa->n_elementi - pa->n_parallelo + 2 + 1];
 	double sdynamic = (pa->sf - pa->si) / (pa->tf - pa->ti) * (t - pa->ti) + pa->si; 
	double s;

	double mbdynv;
	double v;
	
	int el = 0;
	int unk = 0;

	int i, ii, row;
	
	if (displ_hi_freq_low_pass_filter) {
		s = y[pa->n_elementi - pa->n_parallelo + 1] - sstatic;
	} else {
		s = sdynamic - sstatic;
	}
	mbdynv = (pa->vf - pa->vi) / (pa->tf - pa->ti) * (t - pa->ti) + pa->vi;
	if (vel_hi_freq_low_pass_filter) {
		v = y[pa->n_elementi - pa->n_parallelo + 2 + 1 + 1 + 1];
	} else {
		v = mbdynv;
	}

	/* Per ogni componente in parallelo */
	pa->f = pa->f_s = pa->f_v = 0.;
	if (pa->nlrheo_t_cur > t) {
		for (i = 0; i < pa->n_elementi - pa->n_parallelo + 2 + 1 + 1 + 2; i++) {
			pa->yp[i] = pa->yp_prev[i];
		}
		pa->nlrheo_t_cur = pa->nlrheo_t_prev;
	} else {
		pa->nlrheo_t_cur = t;
	}
	for (i = 0; i < pa->n_parallelo; i++) {
		double c[pa->n_serie[i]], k[pa->n_serie[i]];
		double c_partial_s[pa->n_serie[i]], c_partial_v[pa->n_serie[i]];
		double k_partial_s[pa->n_serie[i]], k_partial_v[pa->n_serie[i]];
#if 0
		for (int ii = 0; ii < pa->n_serie[i]; ii++) {
			nlrheo_int_compute_kc(k[ii], c[ii], s, v, el + ii, pa->n_variabili_k, 
				pa->npti_ks, pa->npti_kv, pa->npti_cs, pa->npti_cv,
				pa->ik_s[el+ii], pa->ik_v[el+ii], pa->ic_s[el+ii], pa->ic_v[el+ii], 
				pa->k_s[el+ii], pa->k_v[el+ii], pa->c_s[el+ii], pa->c_v[el+ii], 
				pa->x);
				std::cerr << "k[" << ii << "] " << k[ii] << "c[" << ii <<"] " << c[ii] << std::endl;
		}
#endif
		int nincognite = pa->n_serie[i] - 1;
		c[0] = 0;
		if (nincognite > 0) {
			ii = 0;
			nlrheo_int_compute_kc(&k[ii], &c[ii], 
				&k_partial_s[ii],
				&k_partial_v[ii],
				&c_partial_s[ii],
				&c_partial_v[ii],
				y[unk + ii], 
				pa->yp[unk+ii], el + ii, pa->n_variabili_k, 
				pa->npti_ks, pa->npti_kv, pa->npti_cs, pa->npti_cv,
				pa->ik_s[el+ii], pa->ik_v[el+ii], pa->ic_s[el+ii], pa->ic_v[el+ii], 
				pa->k_s[el+ii], pa->k_v[el+ii], pa->c_s[el+ii], pa->c_v[el+ii], 
				pa->x);
		}
		for (ii = 1; ii < nincognite; ii++) {
			nlrheo_int_compute_kc(&k[ii], &c[ii], 
				&k_partial_s[ii],
				&k_partial_v[ii],
				&c_partial_s[ii],
				&c_partial_v[ii],
				y[unk + ii] - y[unk + ii - 1], 
				pa->yp[unk + ii] - pa->yp[unk + ii - 1], 
				el + ii, pa->n_variabili_k, 
				pa->npti_ks, pa->npti_kv, pa->npti_cs, pa->npti_cv,
				pa->ik_s[el+ii], pa->ik_v[el+ii], pa->ic_s[el+ii], pa->ic_v[el+ii], 
				pa->k_s[el+ii], pa->k_v[el+ii], pa->c_s[el+ii], pa->c_v[el+ii], 
				pa->x);
		}
		for (ii = nincognite; ii < pa->n_serie[i]; ii++) {
			nlrheo_int_compute_kc(&k[ii], &c[ii], 
				&k_partial_s[ii],
				&k_partial_v[ii],
				&c_partial_s[ii],
				&c_partial_v[ii],
				s - y[unk + ii], 
				v - pa->yp[unk+ii], 
				el + ii, pa->n_variabili_k, 
				pa->npti_ks, pa->npti_kv, pa->npti_cs, pa->npti_cv,
				pa->ik_s[el+ii], pa->ik_v[el+ii], pa->ic_s[el+ii], pa->ic_v[el+ii], 
				pa->k_s[el+ii], pa->k_v[el+ii], pa->c_s[el+ii], pa->c_v[el+ii], 
				pa->x);
		}

		if (nincognite > 0) {
			double s1, s2, v1, v2;
			int ints;
			gsl_matrix_set_zero(pa->gsl_C[i]);
			gsl_matrix_set_zero(pa->gsl_C_partial_s[i]);
			gsl_matrix_set_zero(pa->gsl_C_partial_v[i]);
			gsl_matrix_set_zero(pa->gsl_K[i]);
			gsl_matrix_set_zero(pa->gsl_K_partial_s[i]);
			gsl_matrix_set_zero(pa->gsl_K_partial_v[i]);
			gsl_vector_set_zero(pa->gsl_xp[i]);
			gsl_vector_set_zero(pa->gsl_xp_partial_s[i]);
			gsl_vector_set_zero(pa->gsl_xp_partial_v[i]);
			gsl_vector_set_zero(pa->gsl_x[i]);
			gsl_vector_set_zero(pa->gsl_x_partial_s[i]);
			gsl_vector_set_zero(pa->gsl_x_partial_v[i]);
			gsl_vector_set_zero(pa->gsl_b[i]);
			gsl_vector_set_zero(pa->gsl_b_partial_s_for_x[i]);
			gsl_vector_set_zero(pa->gsl_b_partial_v_for_x[i]);
			gsl_vector_set_zero(pa->gsl_b_partial_s_for_xp[i]);
			gsl_vector_set_zero(pa->gsl_b_partial_v_for_xp[i]);
			for (row = 0; row < nincognite; row++) {
				gsl_matrix_set(pa->gsl_C[i], row, row, c[row] + c[row + 1]);
				gsl_matrix_set(pa->gsl_K[i], row, row, k[row] + k[row + 1]);
				if (row == 0) {
					s1 = y[unk + row];
					s2 = y[unk + row + 1] - y[unk + row];
					v1 = pa->yp[unk + row];
					v2 = pa->yp[unk + row + 1] - pa->yp[unk + row];
				} else if (row == nincognite - 1) {
					s1 = y[unk + row] - y[unk + row - 1];
					s2 = s - y[unk + row];
					v1 = pa->yp[unk + row] - pa->yp[unk + row - 1];
					v2 = v - pa->yp[unk + row];
				} else {
					s1 = y[unk + row] - y[unk + row - 1];
					s2 = y[unk + row + 1] - y[unk + row];
					v1 = pa->yp[unk + row] - pa->yp[unk + row - 1];
					v2 = pa->yp[unk + row + 1] - pa->yp[unk + row];
				}
				gsl_matrix_set(pa->gsl_C_partial_v[i], row, row, 
					c_partial_v[row] * v1 + c_partial_v[row + 1] * v2);
				gsl_matrix_set(pa->gsl_K_partial_v[i], row, row, 
					k_partial_v[row] * s1 + k_partial_v[row + 1] * s2);
				
				if (row < nincognite - 1) {
					gsl_matrix_set(pa->gsl_C[i], row, row + 1, - c[row + 1]);
					gsl_matrix_set(pa->gsl_K[i], row, row + 1, - k[row + 1]);
					gsl_matrix_set(pa->gsl_C_partial_v[i], row, row + 1, 
						- c_partial_v[row + 1] * v2);
					gsl_matrix_set(pa->gsl_K_partial_v[i], row, row + 1, 
						- k_partial_v[row + 1] * s2);
				} 
#if 0
				if (row == nincognite - 2) {
					gsl_matrix_set(pa->gsl_C_partial_s[i], row, row + 1, - c_partial_s[row + 1]);
					gsl_matrix_set(pa->gsl_K_partial_s[i], row, row + 1, - k_partial_s[row + 1]);
					gsl_matrix_set(pa->gsl_C_partial_v[i], row, row + 1, - c_partial_v[row + 1]);
					gsl_matrix_set(pa->gsl_K_partial_v[i], row, row + 1, - k_partial_v[row + 1]);
				}
#endif
				if (row > 0) {
					gsl_matrix_set(pa->gsl_C[i], row, row - 1, - c[row]);
					gsl_matrix_set(pa->gsl_K[i], row, row - 1, - k[row]);
					gsl_matrix_set(pa->gsl_C_partial_v[i], row, row - 1, 
						- c_partial_v[row] * v1);
					gsl_matrix_set(pa->gsl_K_partial_v[i], row, row - 1, 
						- k_partial_v[row] * s1);
				}
				gsl_vector_set(pa->gsl_b[i], row, 0.);
				gsl_vector_set(pa->gsl_x[i], row, y[unk + row]);
			}
			gsl_vector_set(pa->gsl_b[i], nincognite - 1, c[nincognite] * v + k[nincognite] * s);
#if 0
			gsl_vector_set(pa->gsl_b_partial_s_for_x[i], nincognite - 1, 
				c_partial_s[nincognite] * v + k_partial_s[nincognite] * s + k[nincognite]);
#endif
			gsl_vector_set(pa->gsl_b_partial_s_for_xp[i], nincognite - 1, 
				- c_partial_v[nincognite] * pa->yp[unk + nincognite - 1]
				- k_partial_v[nincognite] * y[unk + nincognite - 1]
				+ c_partial_s[nincognite] * v2 
				+ k_partial_s[nincognite] * s2 
				+ k[nincognite]);
#if 0
			gsl_vector_set(pa->gsl_b_partial_v_for_x[i], nincognite - 1, 
				c_partial_v[nincognite] * v + k_partial_v[nincognite] * s + c[nincognite]);
#endif
			gsl_vector_set(pa->gsl_b_partial_v_for_xp[i], nincognite - 1, 
				c_partial_v[nincognite] * v2 + k_partial_v[nincognite] * s2 
				+ c[nincognite]);
				
			gsl_matrix_add(pa->gsl_C_partial_v[i], pa->gsl_C[i]);
			gsl_matrix_add(pa->gsl_C_partial_v[i], pa->gsl_K_partial_v[i]);

			/*
			 * calcola b -= Kx, b_{/s}_xp -= K_{/s}x , b_{/v}_xp -= K_{/v}x, 
			                    b_{/s}_x  -= K_{/s}x,  b_{/v}_x  -= K_{/v}p, 
			 */
			gsl_blas_dgemv(CblasNoTrans, -1., pa->gsl_K[i], pa->gsl_x[i], 1., pa->gsl_b[i]);
#if 0
			gsl_blas_dgemv(CblasNoTrans, -1., pa->gsl_K_partial_s[i], pa->gsl_x[i], 1., 
				pa->gsl_b_partial_s_for_xp[i]);
			gsl_blas_dgemv(CblasNoTrans, -1., pa->gsl_K_partial_v[i], pa->gsl_x[i], 1., 
				pa->gsl_b_partial_v_for_xp[i]);
			gsl_blas_dgemv(CblasNoTrans, -1., pa->gsl_K_partial_s[i], pa->gsl_x[i], 1., 
				pa->gsl_b_partial_s_for_x[i]);
			gsl_blas_dgemv(CblasNoTrans, -1., pa->gsl_K_partial_v[i], pa->gsl_x[i], 1., 
				pa->gsl_b_partial_v_for_x[i]);
#endif
			/* calcola xp = C^-1 b */
			gsl_linalg_LU_decomp(pa->gsl_C[i], pa->gsl_perm_C[i], &ints);
			gsl_linalg_LU_solve(pa->gsl_C[i], pa->gsl_perm_C[i], pa->gsl_b[i], pa->gsl_xp[i]);
			
			/* calcola xp_s = C_v^-1 b_s e xp_v = C_v^-1 b_v */
			gsl_linalg_LU_decomp(pa->gsl_C_partial_v[i], pa->gsl_perm_C[i], &ints);
			gsl_linalg_LU_solve(pa->gsl_C_partial_v[i], pa->gsl_perm_C[i], 
				pa->gsl_b_partial_s_for_xp[i], pa->gsl_xp_partial_s[i]);
			gsl_linalg_LU_solve(pa->gsl_C_partial_v[i], pa->gsl_perm_C[i], 
				pa->gsl_b_partial_v_for_xp[i], pa->gsl_xp_partial_v[i]);
			
			/* TODO: setta f e par.f */
			for (ii = 0; ii < nincognite; ii++) {
				pa->yp[unk + ii] = f[unk + ii] = gsl_vector_get(pa->gsl_xp[i], ii);
			}
			pa->f += (c[0] * f[unk] + k[0] * y[unk]);
				
#if 0
			gsl_vector_set(pa->gsl_b[i], nincognite - 1, 
				gsl_vector_get(pa->gsl_b[i], nincognite - 1) 
					- k[nincognite] * s
					+ k[nincognite]);
			gsl_linalg_LU_solve(pa->gsl_C[i], pa->gsl_perm[i], pa->gsl_b[i], pa->gsl_xp[i]);
			pa->f_s += (c[0] * gsl_vector_get(pa->gsl_xp[i], 0) + k[0] * y[unk]);	
#endif
			pa->f_s += k_partial_v[0] * y[unk] * gsl_vector_get(pa->gsl_xp_partial_s[i], 0)
				+ c_partial_v[0] * pa->yp[unk] * gsl_vector_get(pa->gsl_xp_partial_s[i], 0)
				+ c[0] * gsl_vector_get(pa->gsl_xp_partial_s[i], 0);
#if 0
			gsl_vector_set(pa->gsl_b[i], nincognite - 1, 
				gsl_vector_get(pa->gsl_b[i], nincognite - 1) 
					- k[nincognite] - c[nincognite] * v + c[nincognite]);
			gsl_linalg_LU_solve(pa->gsl_C[i], pa->gsl_perm[i], pa->gsl_b[i], pa->gsl_xp[i]);
			pa->f_v += (c[0] * gsl_vector_get(pa->gsl_xp[i], 0) + k[0] * y[unk]);	
#endif
			pa->f_v += k_partial_v[0] * y[unk] * gsl_vector_get(pa->gsl_xp_partial_v[i], 0)
				+ c_partial_v[0] * pa->yp[unk] * gsl_vector_get(pa->gsl_xp_partial_v[i], 0)
				+ c[0] * gsl_vector_get(pa->gsl_xp_partial_v[i], 0);

			el += pa->n_serie[i];
			unk += pa->n_serie[i] - 1;
		} else {
			pa->f += c[0] * v + k[0] * s;
			pa->f_s += k[0] + k_partial_s[0] * s + c_partial_s[0] * v;
			pa->f_v += c[0] + k_partial_v[0] * s + c_partial_v[0] * v;
		}
	}
#if 0
	f[pa->n_elementi - pa->n_parallelo] = -y[pa->n_elementi - pa->n_parallelo] * 2. / pa->hi_freq_force_filter_coeff + 
		-y[pa->n_elementi - pa->n_parallelo + 1] / pa->hi_freq_force_filter_coeff / pa->hi_freq_force_filter_coeff + 
		((pa->x[pa->n_variabili] * std::atan(v / pa->x[pa->n_variabili+1]) +
		pa->x[pa->n_variabili+2] * std::atan(v / pa->x[pa->n_variabili+3]))/2.) / pa->hi_freq_force_filter_coeff / pa->hi_freq_force_filter_coeff;
#endif

	/* filtro forze */
#if 0
	f[pa->n_elementi - pa->n_parallelo] = -y[pa->n_elementi - pa->n_parallelo] * 2. / pa->hi_freq_force_filter_coeff + 
		-y[pa->n_elementi - pa->n_parallelo + 1] / pa->hi_freq_force_filter_coeff / pa->hi_freq_force_filter_coeff + 
		(mbdynv) / pa->hi_freq_force_filter_coeff / pa->hi_freq_force_filter_coeff;
	f[pa->n_elementi - pa->n_parallelo + 1] = y[pa->n_elementi - pa->n_parallelo];
#endif

	/* filtro secondo ordine spostamento */
	pa->yp[pa->n_elementi - pa->n_parallelo] = 
	f[pa->n_elementi - pa->n_parallelo] = 
		- sqrt(2.) * pa->hi_freq_force_filter_coeff * y[pa->n_elementi - pa->n_parallelo]
		- pow(pa->hi_freq_force_filter_coeff, 2) * y[pa->n_elementi - pa->n_parallelo + 1]
		+ pow(pa->hi_freq_force_filter_coeff, 2) * sdynamic;
	pa->yp[pa->n_elementi - pa->n_parallelo + 1] = 
	f[pa->n_elementi - pa->n_parallelo + 1] = 
		y[pa->n_elementi - pa->n_parallelo];

	/* filtro primo ordine spostamento */
	pa->yp[pa->n_elementi - pa->n_parallelo + 2] = 
	f[pa->n_elementi - pa->n_parallelo + 2] = 
		-pa->low_freq_displ_filter_coeff * y[pa->n_elementi - pa->n_parallelo + 2] +
		pa->low_freq_displ_filter_coeff * sdynamic; 
	/* filtro primo ordine velocita' */
	pa->yp[pa->n_elementi - pa->n_parallelo + 2 + 1] = 
	f[pa->n_elementi - pa->n_parallelo + 2 + 1] = 
		-pa->low_freq_displ_filter_coeff * y[pa->n_elementi - pa->n_parallelo + 2 + 1] +
		pa->low_freq_displ_filter_coeff * mbdynv; 

	/* filtro secondo ordine velocita' */
	pa->yp[pa->n_elementi - pa->n_parallelo + 2 + 1 + 1] = 
	f[pa->n_elementi - pa->n_parallelo + 2 + 1 + 1] = 
		- sqrt(2.) * pa->hi_freq_force_filter_coeff * y[pa->n_elementi - pa->n_parallelo + 2 + 1 + 1]
		- pow(pa->hi_freq_force_filter_coeff, 2) * y[pa->n_elementi - pa->n_parallelo + 2 + 1 + 1 + 1]
		+ pow(pa->hi_freq_force_filter_coeff, 2) * mbdynv;
	pa->yp[pa->n_elementi - pa->n_parallelo + 2 + 1 + 1 + 1] = 
	f[pa->n_elementi - pa->n_parallelo + 2 + 1 + 1 + 1] = 
		y[pa->n_elementi - pa->n_parallelo + 2 + 1 + 1];

#if 0
	pa->f += y[pa->n_elementi - pa->n_parallelo + 1];
	if (fabs(f[pa->n_elementi - pa->n_parallelo + 2]) > 5.)
#endif

	pa->f += ((pa->x[pa->n_variabili] * atan(vstatic / pa->x[pa->n_variabili+1]) +
		pa->x[pa->n_variabili+2] * atan(vstatic / pa->x[pa->n_variabili+3]))/2.);

	pa->f += sstatic * pa->static_low_freq_stiffness;
	pa->f_s += pa->static_low_freq_stiffness;

	return GSL_SUCCESS;
}

int
nlrheo_init(sym_params *pa)
{
	int i;

	pa->T = gsl_odeiv_step_rkf45;
	pa->prev_time = 0.;
	pa->current_time = 0.;
	pa->dt = 0.;
	pa->prev_s = 0.;
	pa->prev_sPrime = 0.;
	pa->stepint = gsl_odeiv_step_alloc(pa->T, pa->n_elementi - pa->n_parallelo + 2 + 1 + 1 + 2);
	pa->evolve = gsl_odeiv_evolve_alloc(pa->n_elementi - pa->n_parallelo + 2 + 1 + 1 + 2);
	pa->eps_abs = 1.E-2;
	pa->eps_rel = 1.E-2;
	pa->control = gsl_odeiv_control_standard_new(pa->eps_abs, pa->eps_rel,
		1., 1.);
	pa->sys.function = nlrheo_int_func;
	pa->sys.jacobian = 0;
	pa->sys.dimension = 1;
	pa->sys.params = pa;
	pa->y = (double*) malloc(sizeof(double)*(pa->n_elementi - pa->n_parallelo + 2 + 1 + 1 + 2));
	pa->y_dummy = (double*) malloc(sizeof(double)*(pa->n_elementi - pa->n_parallelo + 2 + 1 + 1 + 2));
	pa->yp = (double*) malloc(sizeof(double)*(pa->n_elementi - pa->n_parallelo + 2 + 1 + 1 + 2));
	pa->yp_prev = (double*) malloc(sizeof(double)*(pa->n_elementi - pa->n_parallelo + 2 + 1 + 1 + 2));
	pa->yp_saved = (double*) malloc(sizeof(double)*(pa->n_elementi - pa->n_parallelo + 2 + 1 + 1 + 2));

	pa->f = 0.;
	for (i = 0; i < pa->n_elementi - pa->n_parallelo + 2 + 1 + 1 + 2; i++) {
		pa->y[i] = pa->y_dummy[i] = pa->yp[i] = pa->yp_prev[i] = pa->yp_saved[i] = 0.;
	}

	return 0;
}

int
nlrheo_destroy(sym_params *pa)
{
	int i;
	
	gsl_odeiv_step_free(pa->stepint);
	gsl_odeiv_evolve_free(pa->evolve);
	gsl_odeiv_control_free(pa->control);
	for (i = 0; i < pa->n_parallelo; i++) {
		if (pa->n_serie[i] - 1 > 0) {
			gsl_matrix_free(pa->gsl_C[i]);
			gsl_matrix_free(pa->gsl_C_partial_s[i]);
			gsl_matrix_free(pa->gsl_C_partial_v[i]);
			gsl_matrix_free(pa->gsl_K[i]);
			gsl_matrix_free(pa->gsl_K_partial_s[i]);
			gsl_matrix_free(pa->gsl_K_partial_v[i]);
			gsl_vector_free(pa->gsl_xp[i]);
			gsl_vector_free(pa->gsl_xp_partial_s[i]);
			gsl_vector_free(pa->gsl_xp_partial_v[i]);
			gsl_vector_free(pa->gsl_x[i]);
			gsl_vector_free(pa->gsl_x_partial_s[i]);
			gsl_vector_free(pa->gsl_x_partial_v[i]);
			gsl_vector_free(pa->gsl_b[i]);
			gsl_vector_free(pa->gsl_b_partial_s_for_x[i]);
			gsl_vector_free(pa->gsl_b_partial_v_for_x[i]);
			gsl_vector_free(pa->gsl_b_partial_s_for_xp[i]);
			gsl_vector_free(pa->gsl_b_partial_v_for_xp[i]);
			gsl_permutation_free(pa->gsl_perm_K[i]);
			gsl_permutation_free(pa->gsl_perm_C[i]);
		}
	}
	free(pa->gsl_C);
	free(pa->gsl_C_partial_s);
	free(pa->gsl_C_partial_v);
	free(pa->gsl_K);
	free(pa->gsl_K_partial_s);
	free(pa->gsl_K_partial_v);
	free(pa->gsl_xp);
	free(pa->gsl_xp_partial_s);
	free(pa->gsl_xp_partial_v);
	free(pa->gsl_x);
	free(pa->gsl_x_partial_s);
	free(pa->gsl_x_partial_v);
	free(pa->gsl_b);
	free(pa->gsl_b_partial_s_for_x);
	free(pa->gsl_b_partial_v_for_x);
	free(pa->gsl_b_partial_s_for_xp);
	free(pa->gsl_b_partial_v_for_xp);
	free(pa->gsl_perm_K);
	free(pa->gsl_perm_C);
	free(pa->n_serie);
	
	free(pa->s_max[0]);
	free(pa->s_max);

	free(pa->v_max[0]);
	free(pa->v_max);
		
	for (i = 0; i < pa->n_elementi; i++) {
		if (pa->npti_ks[i] > 1) {
			gsl_interp_free(pa->ik_s[i]);
		}
		if (pa->npti_kv[i] > 1) {
			gsl_interp_free(pa->ik_v[i]);
		}
		if (pa->npti_cs[i] > 1) {
			gsl_interp_free(pa->ic_s[i]);
		}
		if (pa->npti_cv[i] > 1) {
			gsl_interp_free(pa->ic_v[i]);
		}
	}

	free(pa->npti_ks);
	free(pa->npti_kv);
	free(pa->npti_cs);
	free(pa->npti_cv);
	free(pa->ik_s);
	free(pa->ik_v);
	free(pa->ic_s);
	free(pa->ic_v);

	free(pa->k_s[0]);
	free(pa->k_s);
	free(pa->k_v[0]);
	free(pa->k_v);
	free(pa->c_s[0]);
	free(pa->c_s);
	free(pa->c_v[0]);
	free(pa->c_v);
	
	free(pa->x);
	
	free(pa->y);
	free(pa->y_dummy);
	free(pa->yp);
	free(pa->yp_prev);
	free(pa->yp_saved);

	free(pa);

	return 0;
}

int
nlrheo_update2(sym_params *pa,
	double t_curr, double eps, double epsPrime, int do_try)
{
	int i;
	double t;
	
	double *y = pa->y;
	for (i = 0; i < pa->n_elementi - pa->n_parallelo + 2 + 1 + 1 + 2; i++) {
		pa->yp_prev[i] = pa->yp_saved[i];
		pa->yp[i] = pa->yp_saved[i];
		pa->nlrheo_t_cur = pa->nlrheo_t_prev = pa->prev_time;
	}

	pa->tf = t_curr;
	pa->ti = pa->prev_time;

	pa->dt = (pa->tf - pa->ti) / pa->nsubsteps;

	pa->sf = eps * pa->scale_eps;
	pa->si = pa->prev_s;
	pa->vf = epsPrime * pa->scale_eps;
	pa->vi = pa->prev_sPrime;

	if (do_try) {
		for (i = 0; i < pa->n_elementi - pa->n_parallelo + 2 + 1 + 1 + 2; i++) {
			pa->y_dummy[i] = pa->y[i];
		}
		y = pa->y_dummy;
	}
	
	if (pa->dt > 0.) {
		for (t = pa->prev_time; t < pa->tf; ) {
			int rc;
#if 0
			double tt = t;
			/* FIXME: should add a check on the number
			 * of iterations */
#endif
			gsl_odeiv_evolve_reset(pa->evolve);
			rc = gsl_odeiv_evolve_apply(pa->evolve,
				pa->control, pa->stepint,
				&pa->sys, &t, pa->tf, &pa->dt, y);
			if (rc != GSL_SUCCESS) {
				/* error? */
			}
			for (i = 0; i < pa->n_elementi - pa->n_parallelo + 2 + 1 + 1 + 2; i++) {
				pa->yp_prev[i] = pa->yp[i];
			}

			pa->dt = (pa->dt > pa->dtmin ? pa->dt : pa->dtmin); /* max(pa->dt, pa->dtmin) */
			pa->nlrheo_t_prev = pa->nlrheo_t_cur = t;
#if 0
			std::cerr << "### " << t << " " << t - tt
				<< " " << pa->F << " " << pa->FDE
				<< " " << pa->FDEPrime << std::endl;
#endif
		}
#if 0
		std::cerr << "#########" << std::endl;
#endif

		pa->F = pa->f * pa->scale_f;
		pa->FDE = pa->f_s * pa->scale_f * pa->scale_eps;
		pa->FDEPrime = pa->f_v * pa->scale_f * pa->scale_eps;
#if 0
		std::cerr << t_curr << " " << pa->F << " " << pa->FDE << " " << pa->FDEPrime << " "
			<< 15 * pa->scale_f / pa->scale_eps << "\n";
#endif
	}

	if (!do_try) {
		pa->prev_time = pa->current_time = pa->tf;
		pa->prev_s = eps * pa->scale_eps;
		pa->prev_sPrime = epsPrime * pa->scale_eps;
		for (i = 0; i < pa->n_elementi - pa->n_parallelo + 2 + 1 + 1 + 2; i++) {
			pa->yp_saved[i] = pa->yp[i];
		}
	}

	return 0;
}


int
nlrheo_update(sym_params *pa,
	double t_curr, double eps, double epsPrime, int do_try)
{

	double diffeps = 1.E-10;
	double diffepsp = 1.E-6;

	double feps;
	double fepsp;

	nlrheo_update2(pa, t_curr, eps + diffeps, epsPrime, 1);
	feps = pa->f;
	nlrheo_update2(pa, t_curr, eps, epsPrime + diffepsp, 1);
	fepsp = pa->f;
	nlrheo_update2(pa, t_curr, eps, epsPrime, do_try);

#if 0
 	std::cerr << pa->f << " " << feps << " " << fepsp;
 	std::cerr << " " << pa->FDE / pa->scale_f / pa->scale_eps << " " << 
 		(feps - pa->f) / diffeps / pa->scale_eps;
 	std::cerr << " " << pa->FDEPrime / pa->scale_f / pa->scale_eps << " " << 
 		(fepsp - pa->f) / diffepsp / pa->scale_eps << " " << pa->scale_eps << "\n";
#endif

	pa->FDE = (feps - pa->f) / diffeps * pa->scale_f;
	pa->FDEPrime = (fepsp - pa->f) / diffepsp * pa->scale_f;

#if 0
	fprintf(stderr, "=> F=%e FDE=%e FDEPrime=%e\n", pa->f, pa->FDE, pa->FDEPrime);
#endif
	
	return 0;
}


int
nlrheo_parse(sym_params **pap,
	double scale_eps, double scale_f, double hi_filter,
	double lo_filter, double lo_stiffness, int nsubsteps, double dtmin)
{
	sym_params* pa;
	int i, j;

	*pap = 0;

	pa = (sym_params*) calloc(1, sizeof(sym_params));

	pa->scale_eps = scale_eps;
	pa->scale_f = scale_f;
	pa->hi_freq_force_filter_coeff = hi_filter;
	pa->low_freq_displ_filter_coeff = lo_filter;
	pa->static_low_freq_stiffness = lo_stiffness;
	pa->nsubsteps = nsubsteps;
	pa->dtmin = dtmin;

	if (nlrheo_get_int(&pa->n_parallelo)) {
		return -1;
	}

	pa->gsl_C = (gsl_matrix**) malloc(sizeof(gsl_matrix*) * (pa->n_parallelo));
	pa->gsl_C_partial_s = (gsl_matrix**) malloc(sizeof(gsl_matrix*) * (pa->n_parallelo));
	pa->gsl_C_partial_v = (gsl_matrix**) malloc(sizeof(gsl_matrix*) * (pa->n_parallelo));
	pa->gsl_K = (gsl_matrix**) malloc(sizeof(gsl_matrix*) * (pa->n_parallelo));
	pa->gsl_K_partial_s = (gsl_matrix**) malloc(sizeof(gsl_matrix*) * (pa->n_parallelo));
	pa->gsl_K_partial_v = (gsl_matrix**) malloc(sizeof(gsl_matrix*) * (pa->n_parallelo));
	pa->gsl_xp = (gsl_vector**) malloc(sizeof(gsl_vector*) * (pa->n_parallelo));
	pa->gsl_xp_partial_s = (gsl_vector**) malloc(sizeof(gsl_vector*) * (pa->n_parallelo));
	pa->gsl_xp_partial_v = (gsl_vector**) malloc(sizeof(gsl_vector*) * (pa->n_parallelo));
	pa->gsl_x = (gsl_vector**) malloc(sizeof(gsl_vector*) * (pa->n_parallelo));
	pa->gsl_x_partial_s = (gsl_vector**) malloc(sizeof(gsl_vector*) * (pa->n_parallelo));
	pa->gsl_x_partial_v = (gsl_vector**) malloc(sizeof(gsl_vector*) * (pa->n_parallelo));
	pa->gsl_b = (gsl_vector**) malloc(sizeof(gsl_vector*) * (pa->n_parallelo));
	pa->gsl_b_partial_s_for_x = (gsl_vector**) malloc(sizeof(gsl_vector*) * (pa->n_parallelo));
	pa->gsl_b_partial_v_for_x = (gsl_vector**) malloc(sizeof(gsl_vector*) * (pa->n_parallelo));
	pa->gsl_b_partial_s_for_xp = (gsl_vector**) malloc(sizeof(gsl_vector*) * (pa->n_parallelo));
	pa->gsl_b_partial_v_for_xp = (gsl_vector**) malloc(sizeof(gsl_vector*) * (pa->n_parallelo));
	pa->gsl_perm_K = (gsl_permutation**) malloc(sizeof(gsl_permutation*) * (pa->n_parallelo));
	pa->gsl_perm_C = (gsl_permutation**) malloc(sizeof(gsl_permutation*) * (pa->n_parallelo));
	pa->n_serie = (int*) malloc (sizeof(int) * (pa->n_parallelo));
	for (i = 0; i < pa->n_parallelo; i++) {
		if (nlrheo_get_int(&pa->n_serie[i])) {
			return -1;
		}
		if (pa->n_serie[i] - 1 > 0) {
			pa->gsl_C[i] = gsl_matrix_alloc(pa->n_serie[i] - 1, pa->n_serie[i] - 1);
			pa->gsl_C_partial_s[i] = gsl_matrix_alloc(pa->n_serie[i] - 1, pa->n_serie[i] - 1);
			pa->gsl_C_partial_v[i] = gsl_matrix_alloc(pa->n_serie[i] - 1, pa->n_serie[i] - 1);
			pa->gsl_K[i] = gsl_matrix_alloc(pa->n_serie[i] - 1, pa->n_serie[i] - 1);
			pa->gsl_K_partial_s[i] = gsl_matrix_alloc(pa->n_serie[i] - 1, pa->n_serie[i] - 1);
			pa->gsl_K_partial_v[i] = gsl_matrix_alloc(pa->n_serie[i] - 1, pa->n_serie[i] - 1);
			pa->gsl_xp[i] = gsl_vector_alloc(pa->n_serie[i] - 1);
			pa->gsl_xp_partial_s[i] = gsl_vector_alloc(pa->n_serie[i] - 1);
			pa->gsl_xp_partial_v[i] = gsl_vector_alloc(pa->n_serie[i] - 1);
			pa->gsl_x[i] = gsl_vector_alloc(pa->n_serie[i] - 1);
			pa->gsl_x_partial_s[i] = gsl_vector_alloc(pa->n_serie[i] - 1);
			pa->gsl_x_partial_v[i] = gsl_vector_alloc(pa->n_serie[i] - 1);
			pa->gsl_b[i] = gsl_vector_alloc(pa->n_serie[i] - 1);
			pa->gsl_b_partial_s_for_x[i] = gsl_vector_alloc(pa->n_serie[i] - 1);
			pa->gsl_b_partial_v_for_x[i] = gsl_vector_alloc(pa->n_serie[i] - 1);
			pa->gsl_b_partial_s_for_xp[i] = gsl_vector_alloc(pa->n_serie[i] - 1);
			pa->gsl_b_partial_v_for_xp[i] = gsl_vector_alloc(pa->n_serie[i] - 1);
			pa->gsl_perm_K[i] = gsl_permutation_alloc(pa->n_serie[i] - 1);
			pa->gsl_perm_C[i] = gsl_permutation_alloc(pa->n_serie[i] - 1);
		} else {
			pa->gsl_C[i] = 0;
			pa->gsl_C_partial_s[i] = 0;
			pa->gsl_C_partial_v[i] = 0;
			pa->gsl_K[i] = 0;
			pa->gsl_K_partial_s[i] = 0;
			pa->gsl_K_partial_v[i] = 0;
			pa->gsl_xp[i] = 0;
			pa->gsl_xp_partial_s[i] = 0;
			pa->gsl_xp_partial_v[i] = 0;
			pa->gsl_x[i] = 0;
			pa->gsl_x_partial_s[i] = 0;
			pa->gsl_x_partial_v[i] = 0;
			pa->gsl_b[i] = 0;
			pa->gsl_b_partial_s_for_x[i] = 0;
			pa->gsl_b_partial_v_for_x[i] = 0;
			pa->gsl_b_partial_s_for_xp[i] = 0;
			pa->gsl_b_partial_v_for_xp[i] = 0;
			pa->gsl_perm_K[i] = 0;
			pa->gsl_perm_C[i] = 0;
		}
	}
	if (nlrheo_get_int(&pa->max_n_serie)) {
		return -1;
	}
	
	pa->s_max = (double**) malloc (sizeof(double*) * (pa->max_n_serie));
	pa->s_max[0] = (double*) malloc (sizeof(double) * (pa->max_n_serie*pa->n_parallelo));
	for (i = 0; i < pa->max_n_serie; i++) {
		if (i > 0) {
			pa->s_max[i] = pa->s_max[i - 1] + pa->n_parallelo;
		}

		for (j = 0; j < pa->n_parallelo; j++) {
			if (nlrheo_get_real(&pa->s_max[i][j])) {
				return -1;
			}
		}
	}

	pa->v_max = (double**) malloc (sizeof(double*) * (pa->max_n_serie));
	pa->v_max[0] = (double*) malloc (sizeof(double) * (pa->max_n_serie*pa->n_parallelo));
	for (i = 0; i < pa->max_n_serie; i++) {
		if (i > 0) {
			pa->v_max[i] = pa->v_max[i - 1] + pa->n_parallelo;
		}
		for (j = 0; j < pa->n_parallelo; j++) {
			if (nlrheo_get_real(&pa->v_max[i][j])) {
				return -1;
			}
		}
	}

	if (nlrheo_get_int(&pa->n_elementi)) {
		return -1;
	}

	pa->npti_ks = (int*) malloc (sizeof(int) * (pa->n_elementi));
	for (i = 0; i < pa->n_elementi; i++) {
		if (nlrheo_get_int(&pa->npti_ks[i])) {
			return -1;
		}
	}

	pa->npti_kv = (int*) malloc (sizeof(int) * (pa->n_elementi));
	for (i = 0; i < pa->n_elementi; i++) {
		if (nlrheo_get_int(&pa->npti_kv[i])) {
			return -1;
		}
	}

	pa->npti_cs = (int*) malloc (sizeof(int) * (pa->n_elementi));
	for (i = 0; i < pa->n_elementi; i++) {
		if (nlrheo_get_int(&pa->npti_cs[i])) {
			return -1;
		}
	}

	pa->npti_cv = (int*) malloc (sizeof(int) * (pa->n_elementi));
	for (i = 0; i < pa->n_elementi; i++) {
		if (nlrheo_get_int(&pa->npti_cv[i])) {
			return -1;
		}
	}

	if (nlrheo_get_int(&pa->n_variabili_k)) {
		return -1;
	}
	if (nlrheo_get_int(&pa->n_variabili_c)) {
		return -1;
	}
	if (nlrheo_get_int(&pa->n_variabili)) {
		return -1;
	}

	if (nlrheo_get_int(&pa->max_npti_ks)) {
		return -1;
	}
	if (nlrheo_get_int(&pa->max_npti_kv)) {
		return -1;
	}
	if (nlrheo_get_int(&pa->max_npti_cs)) {
		return -1;
	}
	if (nlrheo_get_int(&pa->max_npti_cv)) {
		return -1;
	}

	pa->k_s = (double**) malloc (sizeof(double*) * (pa->n_elementi));
	pa->k_s[0] = (double*) malloc (sizeof(double) * (pa->n_elementi*pa->max_npti_ks));
	for (i = 0; i < pa->n_elementi; i++) {
		if (i > 0) {
			pa->k_s[i] = pa->k_s[i - 1] + pa->max_npti_ks;
		}
		for (j = 0; j < pa->max_npti_ks; j++) {
			if (nlrheo_get_real(&pa->k_s[i][j])) {
				return -1;
			}
		}
	}
	pa->k_v = (double**) malloc (sizeof(double*) * (pa->n_elementi));
	pa->k_v[0] = (double*) malloc (sizeof(double) * (pa->n_elementi*pa->max_npti_kv));
	for (i = 0; i < pa->n_elementi; i++) {
		if (i > 0) {
			pa->k_v[i] = pa->k_v[i - 1] + pa->max_npti_kv;
		}
		for (j = 0; j < pa->max_npti_kv; j++) {
			if (nlrheo_get_real(&pa->k_v[i][j])) {
				return -1;
			}
		}
	}
	pa->c_s = (double**) malloc (sizeof(double*) * (pa->n_elementi));
	pa->c_s[0] = (double*) malloc (sizeof(double) * (pa->n_elementi*pa->max_npti_cs));
	for (i = 0; i < pa->n_elementi; i++) {
		if (i > 0) {
			pa->c_s[i] = pa->c_s[i - 1] + i * pa->max_npti_cs;
		}
		for (j = 0; j < pa->max_npti_cs; j++) {
			if (nlrheo_get_real(&pa->c_s[i][j])) {
				return -1;
			}
		}
	}
	pa->c_v = (double**) malloc (sizeof(double*) * (pa->n_elementi));
	pa->c_v[0] = (double*) malloc (sizeof(double) * (pa->n_elementi*pa->max_npti_cv));
	for (i = 0; i < pa->n_elementi ; i++) {
		if (i > 0) {
			pa->c_v[i] = pa->c_v[i - 1] + pa->max_npti_cv;
		}
		for (j = 0; j < pa->max_npti_cv; j++) {
			if (nlrheo_get_real(&pa->c_v[i][j])) {
				return -1;
			}
		}
	}


	pa->x = (double*) malloc (sizeof(double) * (pa->n_variabili+4));
	for (i = 0; i < pa->n_variabili+4; i++) {
		if (nlrheo_get_real(&pa->x[i])) {
			return -1;
		}
	}

	pa->ik_s = (gsl_interp**) malloc (sizeof(gsl_interp*) * pa->n_elementi);
	pa->ik_v = (gsl_interp**) malloc (sizeof(gsl_interp*) * pa->n_elementi);
	pa->ic_s = (gsl_interp**) malloc (sizeof(gsl_interp*) * pa->n_elementi);
	pa->ic_v = (gsl_interp**) malloc (sizeof(gsl_interp*) * pa->n_elementi);
	for (i = 0; i < pa->n_elementi; i++) {
		if (pa->npti_ks[i] > 1) {
			pa->ik_s[i] = gsl_interp_alloc(gsl_interp_linear, pa->npti_ks[i]);
			gsl_interp_init(pa->ik_s[i], pa->k_s[i], pa->k_s[i], pa->npti_ks[i]);
		} else  {
			pa->ik_s[i] = 0;
		}
		if (pa->npti_kv[i] > 1) {
			pa->ik_v[i] = gsl_interp_alloc(gsl_interp_linear, pa->npti_kv[i]);
			gsl_interp_init(pa->ik_v[i], pa->k_v[i], pa->k_v[i], pa->npti_kv[i]);
		} else  {
			pa->ik_v[i] = 0;
		}
		if (pa->npti_cs[i] > 1) {
			pa->ic_s[i] = gsl_interp_alloc(gsl_interp_linear, pa->npti_cs[i]);
			gsl_interp_init(pa->ic_s[i], pa->c_s[i], pa->c_s[i], pa->npti_cs[i]);
		} else  {
			pa->ic_s[i] = 0;
		}
		if (pa->npti_cv[i] > 1) {
			pa->ic_v[i] = gsl_interp_alloc(gsl_interp_linear, pa->npti_cv[i]);
			gsl_interp_init(pa->ic_v[i], pa->c_v[i], pa->c_v[i], pa->npti_cv[i]);
		} else  {
			pa->ic_v[i] = 0;
		}
	}

	*pap = pa;

	return 0;
}

