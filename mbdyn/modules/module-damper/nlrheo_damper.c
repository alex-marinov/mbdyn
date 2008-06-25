
#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "nlrheo_damper.h"

#include <cmath>
#include <cfloat>
#include <algorithm>
#include <iostream>

static void
nlrheo_int_compute_kc(double &k, double &c, 
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
	// calcola k di el dato s,v
	double k_v_i[npti_kv[el]];
	int somma_k = 0;
	for (int e = 0; e < el; e++) {
		somma_k += npti_ks[e] * npti_kv[e];
	}
	if (ik_s != 0) {
		for (int i = 0; i < npti_kv[el]; i++) {
			if (std::abs(s) > k_s[npti_ks[el]-1]) {
				k_v_i[i] = x[somma_k + npti_ks[el] - 1];
// 					(x[somma_k + npti_ks[el] - 1] - x[somma_k + npti_ks[el] - 2])
// 					/ (k_s[npti_ks[el]-1] - k_s[npti_ks[el]-2])
// 					* (std::abs(s) - k_s[npti_ks[el]-1])
// 					+ x[somma_k + npti_ks[el] - 1]
// 				;
			} else if (std::abs(s) == 0) {
				k_v_i[i] = x[somma_k];
			} else {
				double tx[npti_ks[el]];
				for (int ii = 0; ii < npti_ks[el]; ii++) {
					tx[ii] = x[somma_k + ii];
				}
				k_v_i[i] = gsl_interp_eval(ik_s, k_s, tx, std::abs(s), 0);
			}
			somma_k += npti_ks[el];
		}
	} else {
		for (int i = 0; i < npti_kv[el]; i++) {
			k_v_i[i] = x[somma_k + i];
		}
	}
	if (ik_v != 0) {
		if (std::abs(v) > k_v[npti_kv[el] - 1]) {
			k = k_v_i[npti_kv[el] - 1];
// 				(k_v_i[npti_kv[el]-1] - k_v_i[npti_kv[el]-2]) 
// 				/ (k_v[npti_kv[el] - 1] - k_v[npti_kv[el] - 2])
// 				* (std::abs(v) - k_v[npti_kv[el] - 1])
// 				+ k_v_i[npti_kv[el]-1]
// 			;
		} else if (std::abs(v) == 0) {
			k = k_v_i[0];
		} else {
			k = gsl_interp_eval(ik_v, k_v, k_v_i, std::abs(v), 0);
		}
	} else {
		k = k_v_i[0];
	}
	
	// calcola c di el dato s,v
	double c_v_i[npti_cv[el]];
	int somma_c = n_variabili_k;
	for (int e = 0; e < el; e++) {
		somma_c += npti_cs[e] * npti_cv[e];
	}
	if (ic_s != 0) {
		for (int i = 0; i < npti_cv[el]; i++) {
			if (std::abs(s) > c_s[npti_cs[el] - 1]) {
				c_v_i[i] = x[somma_c + npti_cs[el] - 1];
// 					(x[somma_c + npti_cs[el] - 1] - x[somma_c + npti_cs[el] - 2])
// 					/ (c_s[npti_cs[el] - 1] - c_s[npti_cs[el] - 2])
// 					* (std::abs(s) - c_s[npti_cs[el] - 1])
// 					+ x[somma_c + npti_cs[el] - 1]
// 				;
			} else if (std::abs(s) == 0) {
				c_v_i[i] = x[somma_c];
			}
			double tx[npti_cs[el]];
			for (int ii = 0; ii < npti_cs[el]; ii++) {
				tx[ii] = x[somma_c + ii];
			}
			c_v_i[i] = gsl_interp_eval(ic_s, c_s, tx, std::abs(s), 0);
			somma_c += npti_cs[el];
		}
	} else {
		for (int i = 0; i < npti_cv[el]; i++) {
			c_v_i[i] = x[somma_c + i];
		}
	}
	if (ic_v != 0) {
		if (std::abs(v) > c_v[npti_cv[el] - 1]) {
			c = c_v_i[npti_cv[el] - 1];
// 				(c_v_i[npti_cv[el]-1] - c_v_i[npti_cv[el]-2])
// 				/ (c_v[npti_cv[el] - 1] - c_v[npti_cv[el] - 2])
// 				* (std::abs(v) - c_v[npti_cv[el] - 1]) +
// 				+ c_v_i[npti_cv[el]-1]
// 			;
		} else if (std::abs(v) == 0) {
			c = c_v_i[0];
		} else {
			c = gsl_interp_eval(ic_v, c_v, c_v_i, std::abs(v), 0);
		}
	} else {
		c = c_v_i[0];
	}
// 	std::cerr << std::abs(s) << " " << std::abs(v) << " " << k << " " << c << std::endl;
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

extern "C" int
nlrheo_int_func(double t, const double y[], double f[], void *para)
{
	sym_params & pa = *((sym_params *)para);
	double sstatic = y[pa.n_elementi - pa.n_parallelo + 2];
 	double sdynamic = (pa.sf - pa.si) / (pa.tf - pa.ti) * (t - pa.ti) + pa.si; 
	double s = sdynamic - sstatic;
// 	double s = y[pa.n_elementi - pa.n_parallelo + 1] - sstatic;
	double mbdynv = (pa.vf - pa.vi) / (pa.tf - pa.ti) * (t - pa.ti) + pa.vi;
	double v = mbdynv;//y[pa.n_elementi - pa.n_parallelo];

// 	Per ogni componente in parallelo
	int el = 0;
	int unk = 0;
	pa.f = pa.f_s = pa.f_v = 0.;
	if (pa.nlrheo_t_cur > t) {
		for (int i = 0; i < pa.n_elementi - pa.n_parallelo + 2 + 1; i++) {
			pa.yp[i] = pa.yp_prev[i];
		}
		pa.nlrheo_t_cur = pa.nlrheo_t_prev;
	} else {
		pa.nlrheo_t_cur = t;
	}
	for (int i = 0; i < pa.n_parallelo; i++) {
		double c[pa.n_serie[i]], k[pa.n_serie[i]];
// 		for (int ii = 0; ii < pa.n_serie[i]; ii++) {
// 			nlrheo_int_compute_kc(k[ii], c[ii], s, v, el + ii, pa.n_variabili_k, 
// 				pa.npti_ks, pa.npti_kv, pa.npti_cs, pa.npti_cv,
// 				pa.ik_s[el+ii], pa.ik_v[el+ii], pa.ic_s[el+ii], pa.ic_v[el+ii], 
// 				pa.k_s[el+ii], pa.k_v[el+ii], pa.c_s[el+ii], pa.c_v[el+ii], 
// 				pa.x);
// 				// std::cerr << "k[" << ii << "] " << k[ii] << "c[" << ii <<"] " << c[ii] << std::endl;
// 		}
		int nincognite = pa.n_serie[i] - 1;
		if (nincognite > 0) {
			int ii = 0;
			nlrheo_int_compute_kc(k[ii], c[ii], 
				y[unk + ii], 
				pa.yp[unk+ii], el + ii, pa.n_variabili_k, 
				pa.npti_ks, pa.npti_kv, pa.npti_cs, pa.npti_cv,
				pa.ik_s[el+ii], pa.ik_v[el+ii], pa.ic_s[el+ii], pa.ic_v[el+ii], 
				pa.k_s[el+ii], pa.k_v[el+ii], pa.c_s[el+ii], pa.c_v[el+ii], 
				pa.x);
		}
		for (int ii = 1; ii < nincognite; ii++) {
			nlrheo_int_compute_kc(k[ii], c[ii], 
				y[unk + ii] - y[unk + ii - 1], 
				pa.yp[unk + ii] - pa.yp[unk + ii - 1], 
				el + ii, pa.n_variabili_k, 
				pa.npti_ks, pa.npti_kv, pa.npti_cs, pa.npti_cv,
				pa.ik_s[el+ii], pa.ik_v[el+ii], pa.ic_s[el+ii], pa.ic_v[el+ii], 
				pa.k_s[el+ii], pa.k_v[el+ii], pa.c_s[el+ii], pa.c_v[el+ii], 
				pa.x);
		}
		for (int ii = nincognite; ii < pa.n_serie[i]; ii++) {
			nlrheo_int_compute_kc(k[ii], c[ii], 
				s - y[unk + ii], 
				v - pa.yp[unk+ii], 
				el + ii, pa.n_variabili_k, 
				pa.npti_ks, pa.npti_kv, pa.npti_cs, pa.npti_cv,
				pa.ik_s[el+ii], pa.ik_v[el+ii], pa.ic_s[el+ii], pa.ic_v[el+ii], 
				pa.k_s[el+ii], pa.k_v[el+ii], pa.c_s[el+ii], pa.c_v[el+ii], 
				pa.x);
		}

		if (nincognite > 0) {
			gsl_matrix_set_zero(pa.gsl_C[i]);
			gsl_matrix_set_zero(pa.gsl_K[i]);
			gsl_vector_set_zero(pa.gsl_xp[i]);
			gsl_vector_set_zero(pa.gsl_x[i]);
			gsl_vector_set_zero(pa.gsl_b[i]);
			for (int row = 0; row < nincognite; row++) {
				gsl_matrix_set(pa.gsl_C[i], row, row, c[row] + c[row + 1]);
				gsl_matrix_set(pa.gsl_K[i], row, row, k[row] + k[row + 1]);
				if (row < nincognite - 1) {
					gsl_matrix_set(pa.gsl_C[i], row, row + 1, - c[row + 1]);
					gsl_matrix_set(pa.gsl_K[i], row, row + 1, - k[row + 1]);
				}
				if (row > 0) {
					gsl_matrix_set(pa.gsl_C[i], row, row - 1, - c[row]);
					gsl_matrix_set(pa.gsl_K[i], row, row - 1, - k[row]);
				}
				gsl_vector_set(pa.gsl_b[i], row, 0.);
				gsl_vector_set(pa.gsl_x[i], row, y[unk + row]);
			}
			gsl_vector_set(pa.gsl_b[i], nincognite - 1, c[nincognite] * v + k[nincognite] * s);

			// calcola b -= Kx
			gsl_blas_dgemv(CblasNoTrans, -1., pa.gsl_K[i], pa.gsl_x[i], 1., pa.gsl_b[i]);
			// calcola xp = C^-1 b
			int ints;
			gsl_linalg_LU_decomp(pa.gsl_C[i], pa.gsl_perm[i], &ints);
			gsl_linalg_LU_solve(pa.gsl_C[i], pa.gsl_perm[i], pa.gsl_b[i], pa.gsl_xp[i]);
			// TODO: setta f e par.f
			for (int ii = 0; ii < nincognite; ii++) {
				pa.yp[unk + ii] = f[unk + ii] = gsl_vector_get(pa.gsl_xp[i], ii);
			}
			pa.f += (c[0] * f[unk] + k[0] * y[unk]);
				
			gsl_vector_set(pa.gsl_b[i], nincognite - 1, 
				gsl_vector_get(pa.gsl_b[i], nincognite - 1) 
					- k[nincognite] * s
					+ k[nincognite]);
			gsl_linalg_LU_solve(pa.gsl_C[i], pa.gsl_perm[i], pa.gsl_b[i], pa.gsl_xp[i]);
			pa.f_s += (c[0] * gsl_vector_get(pa.gsl_xp[i], 0) + k[0] * y[unk]);	

			gsl_vector_set(pa.gsl_b[i], nincognite - 1, 
				gsl_vector_get(pa.gsl_b[i], nincognite - 1) 
					- k[nincognite] - c[nincognite] * v + c[nincognite]);
			gsl_linalg_LU_solve(pa.gsl_C[i], pa.gsl_perm[i], pa.gsl_b[i], pa.gsl_xp[i]);
			pa.f_v += (c[0] * gsl_vector_get(pa.gsl_xp[i], 0) + k[0] * y[unk]);	

			el += pa.n_serie[i];
			unk += pa.n_serie[i] - 1;
		} else {
			pa.f += c[0] * v + k[0] * s;
			pa.f_s += k[0];
			pa.f_v += c[0];
		}
	}
// 	f[pa.n_elementi - pa.n_parallelo] = -y[pa.n_elementi - pa.n_parallelo] * 2. / pa.hi_freq_force_filter_coeff + 
// 		-y[pa.n_elementi - pa.n_parallelo + 1] / pa.hi_freq_force_filter_coeff / pa.hi_freq_force_filter_coeff + 
// 		((pa.x[pa.n_variabili] * std::atan(v / pa.x[pa.n_variabili+1]) +
// 		pa.x[pa.n_variabili+2] * std::atan(v / pa.x[pa.n_variabili+3]))/2.) / pa.hi_freq_force_filter_coeff / pa.hi_freq_force_filter_coeff;

//filtro forze
// 	f[pa.n_elementi - pa.n_parallelo] = -y[pa.n_elementi - pa.n_parallelo] * 2. / pa.hi_freq_force_filter_coeff + 
// 		-y[pa.n_elementi - pa.n_parallelo + 1] / pa.hi_freq_force_filter_coeff / pa.hi_freq_force_filter_coeff + 
// 		(mbdynv) / pa.hi_freq_force_filter_coeff / pa.hi_freq_force_filter_coeff;
// 	f[pa.n_elementi - pa.n_parallelo + 1] = y[pa.n_elementi - pa.n_parallelo];
	pa.yp[pa.n_elementi - pa.n_parallelo] = 
	f[pa.n_elementi - pa.n_parallelo] = 
		- std::sqrt(2.) * pa.hi_freq_force_filter_coeff * y[pa.n_elementi - pa.n_parallelo]
		- std::pow(pa.hi_freq_force_filter_coeff, 2) * y[pa.n_elementi - pa.n_parallelo + 1]
		+ std::pow(pa.hi_freq_force_filter_coeff, 2) * sdynamic;
	pa.yp[pa.n_elementi - pa.n_parallelo + 1] = 
	f[pa.n_elementi - pa.n_parallelo + 1] = 
		y[pa.n_elementi - pa.n_parallelo];

	pa.yp[pa.n_elementi - pa.n_parallelo + 2] = 
	f[pa.n_elementi - pa.n_parallelo + 2] = 
		-pa.low_freq_displ_filter_coeff * y[pa.n_elementi - pa.n_parallelo + 2] +
		pa.low_freq_displ_filter_coeff * sdynamic; 
// 	pa.f += y[pa.n_elementi - pa.n_parallelo + 1];
//	if (std::abs(f[pa.n_elementi - pa.n_parallelo + 2]) > 5.)
		pa.f += ((pa.x[pa.n_variabili] * std::atan(f[pa.n_elementi - pa.n_parallelo + 2] / pa.x[pa.n_variabili+1]) +
			pa.x[pa.n_variabili+2] * std::atan(f[pa.n_elementi - pa.n_parallelo + 2] / pa.x[pa.n_variabili+3]))/2.);
	pa.f += sstatic * pa.static_low_freq_stiffness;
	pa.f_s += pa.static_low_freq_stiffness;
	//std::cerr << t << " " << pa.f << std::endl;

	return GSL_SUCCESS;
}

extern "C" int
nlrheo_init(sym_params *nlrheo)
{
	sym_params& pa = *nlrheo;

	pa.T = gsl_odeiv_step_rkf45;
	pa.prev_time = 0.;
	pa.current_time = 0.;
	pa.dt = 0.;
	pa.prev_s = 0.;
	pa.prev_sPrime = 0.;
	pa.stepint = gsl_odeiv_step_alloc(pa.T, pa.n_elementi - pa.n_parallelo + 2 + 1);
	pa.evolve = gsl_odeiv_evolve_alloc(pa.n_elementi - pa.n_parallelo + 2 + 1);
	pa.eps_abs = 1.E-2;
	pa.eps_rel = 1.E-2;
	pa.control = gsl_odeiv_control_standard_new(pa.eps_abs, pa.eps_rel,
		1., 1.);
	pa.sys.function = nlrheo_int_func;
	pa.sys.jacobian = 0;
	pa.sys.dimension = 1;
	pa.sys.params = &pa;
	pa.y = new double[pa.n_elementi - pa.n_parallelo + 2 + 1];
	pa.y_dummy = new double[pa.n_elementi - pa.n_parallelo + 2 + 1];
	pa.yp = new double[pa.n_elementi - pa.n_parallelo + 2 + 1];
	pa.yp_prev = new double[pa.n_elementi - pa.n_parallelo + 2 + 1];
	pa.yp_saved = new double[pa.n_elementi - pa.n_parallelo + 2 + 1];

	pa.f = 0.;
	for (int i = 0; i < pa.n_elementi - pa.n_parallelo + 2 + 1; i++) {
		pa.y[i] = pa.y_dummy[i] = pa.yp[i] = pa.yp_prev[i] = pa.yp_saved[i] = 0.;
	}

	return 0;
}

extern "C" int
nlrheo_destroy(sym_params *nlrheo)
{
	sym_params& pa = *nlrheo;

	gsl_odeiv_step_free(pa.stepint);
	gsl_odeiv_evolve_free(pa.evolve);
	for (int i = 0; i < pa.n_parallelo; i++) {
		if (pa.n_serie[i] - 1 > 0) {
			gsl_matrix_free(pa.gsl_C[i]);
			gsl_matrix_free(pa.gsl_K[i]);
			gsl_vector_free(pa.gsl_xp[i]);
			gsl_vector_free(pa.gsl_x[i]);
			gsl_vector_free(pa.gsl_b[i]);
			gsl_permutation_free(pa.gsl_perm[i]);
		}
	}
	delete[] pa.gsl_C;
	delete[] pa.gsl_K;
	delete[] pa.gsl_xp;
	delete[] pa.gsl_x;
	delete[] pa.gsl_b;
	delete[] pa.gsl_perm;
	
	delete[] pa.s_max[0];
	delete[] pa.s_max;

	delete[] pa.v_max[0];
	delete[] pa.v_max;
		
	for (int i = 0; i < pa.n_elementi; i++) {
		if (pa.npti_ks[i] > 1) {
			gsl_interp_free(pa.ik_s[i]);
		}
		if (pa.npti_kv[i] > 1) {
			gsl_interp_free(pa.ik_v[i]);
		}
		if (pa.npti_cs[i] > 1) {
			gsl_interp_free(pa.ic_s[i]);
		}
		if (pa.npti_cv[i] > 1) {
			gsl_interp_free(pa.ic_v[i]);
		}
	}

	delete[] pa.npti_ks;
	delete[] pa.npti_kv;
	delete[] pa.npti_cs;
	delete[] pa.npti_cv;
	delete[] pa.ik_s;
	delete[] pa.ik_v;
	delete[] pa.ic_s;
	delete[] pa.ic_v;

	delete[] pa.k_s[0];
	delete[] pa.k_s;
	delete[] pa.k_v[0];
	delete[] pa.k_v;
	delete[] pa.c_s[0];
	delete[] pa.c_s;
	delete[] pa.c_v[0];
	delete[] pa.c_v;
	
	delete[] pa.x;
	
	delete[] pa.y;
	delete[] pa.y_dummy;
	delete[] pa.yp;
	delete[] pa.yp_prev;
	delete[] pa.yp_saved;

	return 0;
};
#include <iostream>
extern "C" int
nlrheo_update(sym_params *nlrheo,
	double t_curr, double eps, double epsPrime, int do_try)
{
	sym_params& pa = *nlrheo;

	double *y = pa.y;
	for (int i = 0; i < pa.n_elementi - pa.n_parallelo + 2 + 1; i++) {
		pa.yp[i] = pa.yp_prev[i] = pa.yp_saved[i];
		pa.nlrheo_t_cur = pa.nlrheo_t_prev = pa.prev_time;
	}

	pa.tf = t_curr;
	pa.ti = pa.prev_time;

	pa.dt = (pa.tf - pa.ti) / pa.nsubsteps;

	pa.sf = eps * pa.scale_eps;
	pa.si = pa.prev_s;
	pa.vf = epsPrime * pa.scale_eps;
	pa.vi = pa.prev_sPrime;

	if (do_try) {
		for (int i = 0; i < pa.n_elementi - pa.n_parallelo + 2 + 1; i++) {
			pa.y_dummy[i] = pa.y[i];
		}
		y = pa.y_dummy;
	}
	
	if (pa.dt > 0.) {
		for (double t = pa.prev_time; t < pa.tf; ) {
			// double tt = t;
			// FIXME: should add a check on the number
			// of iterations
			gsl_odeiv_evolve_reset(pa.evolve);
			int rc = gsl_odeiv_evolve_apply(pa.evolve,
				pa.control, pa.stepint,
				&pa.sys, &t, pa.tf, &pa.dt, y);
			if (rc != GSL_SUCCESS) {
				// error?
			}
			for (int i = 0; i < pa.n_elementi - pa.n_parallelo + 2 + 1; i++) {
				pa.yp_prev[i] = pa.yp[i];
			}

			pa.dt = std::max(pa.dt, pa.dtmin);
			pa.nlrheo_t_prev = pa.nlrheo_t_cur = t;
			// std::cerr << "### " << t << " " << t - tt
			// 	<< " " << pa.F << " " << pa.FDE
			// 	<< " " << pa.FDEPrime << std::endl;
		}
		// std::cerr << "#########" << std::endl;

		pa.F = pa.f * pa.scale_f;
		pa.FDE = pa.f_s * pa.scale_f / pa.scale_eps;
		pa.FDEPrime = pa.f_v * pa.scale_f / pa.scale_eps;
	}

	if (!do_try) {
		pa.prev_time = pa.current_time = pa.tf;
		pa.prev_s = eps * pa.scale_eps;
		pa.prev_sPrime = epsPrime * pa.scale_eps;
		for (int i = 0; i < pa.n_elementi - pa.n_parallelo + 2 + 1; i++) {
			pa.yp_saved[i] = pa.yp[i];
		}
	}

	return 0;
}

extern "C" int
nlrheo_parse(sym_params **nlrheop,
	double scale_eps, double scale_f, double hi_filter,
	double lo_filter, double lo_stiffness, int nsubsteps, double dtmin)
{
	*nlrheop = 0;

	sym_params* pap = new sym_params();
	sym_params &pa(*pap);

	pa.scale_eps = scale_eps;
	pa.scale_f = scale_f;
	pa.hi_freq_force_filter_coeff = hi_filter;
	pa.low_freq_displ_filter_coeff = lo_filter;
	pa.static_low_freq_stiffness = lo_stiffness;
	pa.nsubsteps = nsubsteps;
	pa.dtmin = dtmin;

	if (nlrheo_get_int(&pa.n_parallelo)) {
		return -1;
	}

	pa.gsl_C = new gsl_matrix*[pa.n_parallelo];
	pa.gsl_K = new gsl_matrix*[pa.n_parallelo];
	pa.gsl_xp = new gsl_vector*[pa.n_parallelo];
	pa.gsl_x = new gsl_vector*[pa.n_parallelo];
	pa.gsl_b = new gsl_vector*[pa.n_parallelo];
	pa.gsl_perm = new gsl_permutation*[pa.n_parallelo];
	pa.n_serie = new int[pa.n_parallelo];
	for (int i = 0; i < pa.n_parallelo; i++) {
		if (nlrheo_get_int(&pa.n_serie[i])) {
			return -1;
		}
		if (pa.n_serie[i] - 1 > 0) {
			pa.gsl_C[i] = gsl_matrix_alloc(pa.n_serie[i] - 1, pa.n_serie[i] - 1);
			pa.gsl_K[i] = gsl_matrix_alloc(pa.n_serie[i] - 1, pa.n_serie[i] - 1);
			pa.gsl_xp[i] = gsl_vector_alloc(pa.n_serie[i] - 1);
			pa.gsl_x[i] = gsl_vector_alloc(pa.n_serie[i] - 1);
			pa.gsl_b[i] = gsl_vector_alloc(pa.n_serie[i] - 1);
			pa.gsl_perm[i] = gsl_permutation_alloc(pa.n_serie[i] - 1);
		} else {
			pa.gsl_C[i] = 0;
			pa.gsl_K[i] = 0;
			pa.gsl_xp[i] = 0;
			pa.gsl_x[i] = 0;
			pa.gsl_b[i] = 0;
			pa.gsl_perm[i] = 0;
		}
	}
	if (nlrheo_get_int(&pa.max_n_serie)) {
		return -1;
	}
	
	pa.s_max = new double*[pa.max_n_serie];
	pa.s_max[0] = new double[pa.max_n_serie*pa.n_parallelo];
	for (int i = 0; i < pa.max_n_serie; i++) {
		pa.s_max[i] = pa.s_max[0] + i * pa.n_parallelo;
		for (int j = 0; j < pa.n_parallelo; j++) {
			if (nlrheo_get_real(&pa.s_max[i][j])) {
				return -1;
			}
		}
	}

	pa.v_max = new double*[pa.max_n_serie];
	pa.v_max[0] = new double[pa.max_n_serie*pa.n_parallelo];
	for (int i = 0; i < pa.max_n_serie; i++) {
		pa.v_max[i] = pa.v_max[0] + i * pa.n_parallelo;
		for (int j = 0; j < pa.n_parallelo; j++) {
			if (nlrheo_get_real(&pa.v_max[i][j])) {
				return -1;
			}
		}
	}

	if (nlrheo_get_int(&pa.n_elementi)) {
		return -1;
	}

	pa.npti_ks = new int[pa.n_elementi];
	for (int i = 0; i < pa.n_elementi; i++) {
		if (nlrheo_get_int(&pa.npti_ks[i])) {
			return -1;
		}
	}

	pa.npti_kv = new int[pa.n_elementi];
	for (int i = 0; i < pa.n_elementi; i++) {
		if (nlrheo_get_int(&pa.npti_kv[i])) {
			return -1;
		}
	}

	pa.npti_cs = new int[pa.n_elementi];
	for (int i = 0; i < pa.n_elementi; i++) {
		if (nlrheo_get_int(&pa.npti_cs[i])) {
			return -1;
		}
	}

	pa.npti_cv = new int[pa.n_elementi];
	for (int i = 0; i < pa.n_elementi; i++) {
		if (nlrheo_get_int(&pa.npti_cv[i])) {
			return -1;
		}
	}

	if (nlrheo_get_int(&pa.n_variabili_k)) {
		return -1;
	}
	if (nlrheo_get_int(&pa.n_variabili_c)) {
		return -1;
	}
	if (nlrheo_get_int(&pa.n_variabili)) {
		return -1;
	}

	if (nlrheo_get_int(&pa.max_npti_ks)) {
		return -1;
	}
	if (nlrheo_get_int(&pa.max_npti_kv)) {
		return -1;
	}
	if (nlrheo_get_int(&pa.max_npti_cs)) {
		return -1;
	}
	if (nlrheo_get_int(&pa.max_npti_cv)) {
		return -1;
	}

	pa.k_s = new double*[pa.n_elementi];
	pa.k_s[0] = new double[pa.n_elementi*pa.max_npti_ks];
	for (int i = 0; i < pa.n_elementi; i++) {
		pa.k_s[i] = pa.k_s[0] + i * pa.max_npti_ks;
		for (int j = 0; j < pa.max_npti_ks; j++) {
			if (nlrheo_get_real(&pa.k_s[i][j])) {
				return -1;
			}
		}
	}
	pa.k_v = new double*[pa.n_elementi];
	pa.k_v[0] = new double[pa.n_elementi*pa.max_npti_kv];
	for (int i = 0; i < pa.n_elementi; i++) {
		pa.k_v[i] = pa.k_v[0] + i * pa.max_npti_kv;
		for (int j = 0; j < pa.max_npti_kv; j++) {
			if (nlrheo_get_real(&pa.k_v[i][j])) {
				return -1;
			}
		}
	}
	pa.c_s = new double*[pa.n_elementi];
	pa.c_s[0] = new double[pa.n_elementi*pa.max_npti_cs];
	for (int i = 0; i < pa.n_elementi; i++) {
		pa.c_s[i] = pa.c_s[0] + i * pa.max_npti_cs;
		for (int j = 0; j < pa.max_npti_cs; j++) {
			if (nlrheo_get_real(&pa.c_s[i][j])) {
				return -1;
			}
		}
	}
	pa.c_v = new double*[pa.n_elementi];
	pa.c_v[0] = new double[pa.n_elementi*pa.max_npti_cv];
	for (int i = 0; i < pa.n_elementi ; i++) {
		pa.c_v[i] = pa.c_v[0] + i * pa.max_npti_cv;
		for (int j = 0; j < pa.max_npti_cv; j++) {
			if (nlrheo_get_real(&pa.c_v[i][j])) {
				return -1;
			}
		}
	}


	pa.x = new double[pa.n_variabili+4];
	for (int i = 0; i < pa.n_variabili+4; i++) {
		if (nlrheo_get_real(&pa.x[i])) {
			return -1;
		}
	}

	pa.ik_s = new gsl_interp *[pa.n_elementi];
	pa.ik_v = new gsl_interp *[pa.n_elementi];
	pa.ic_s = new gsl_interp *[pa.n_elementi];
	pa.ic_v = new gsl_interp *[pa.n_elementi];
	for (int i = 0; i < pa.n_elementi; i++) {
		if (pa.npti_ks[i] > 1) {
			pa.ik_s[i] = gsl_interp_alloc(gsl_interp_linear, pa.npti_ks[i]);
				gsl_interp_init(pa.ik_s[i], pa.k_s[i], pa.k_s[i], pa.npti_ks[i]);
		} else  {
			pa.ik_s[i] = 0;
		}
		if (pa.npti_kv[i] > 1) {
			pa.ik_v[i] = gsl_interp_alloc(gsl_interp_linear, pa.npti_kv[i]);
				gsl_interp_init(pa.ik_v[i], pa.k_v[i], pa.k_v[i], pa.npti_kv[i]);
		} else  {
			pa.ik_v[i] = 0;
		}
		if (pa.npti_cs[i] > 1) {
			pa.ic_s[i] = gsl_interp_alloc(gsl_interp_linear, pa.npti_cs[i]);
				gsl_interp_init(pa.ic_s[i], pa.c_s[i], pa.c_s[i], pa.npti_cs[i]);
		} else  {
			pa.ic_s[i] = 0;
		}
		if (pa.npti_cv[i] > 1) {
			pa.ic_v[i] = gsl_interp_alloc(gsl_interp_linear, pa.npti_cv[i]);
				gsl_interp_init(pa.ic_v[i], pa.c_v[i], pa.c_v[i], pa.npti_cv[i]);
		} else  {
			pa.ic_v[i] = 0;
		}
	}

	*nlrheop = &pa;

	return 0;
}

