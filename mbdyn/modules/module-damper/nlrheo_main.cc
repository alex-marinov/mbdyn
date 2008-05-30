/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2008
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
 */

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>


#define _GNU_SOURCE 1
#include <fenv.h>
static void __attribute__ ((constructor))
trapfpe ()
{
  /* Enable some exceptions.  At startup all exceptions are masked.  */

  feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}

#include "nlrheo_damper.h"

extern "C" int
nlrheo_int_func(double t, const double y[], double f[], void *para);

std::istream *fdatainp;

int
nlrheo_get_int(int *i)
{
	(*fdatainp) >> *i;
	return 0;
}

int
nlrheo_get_real(double *d)
{
	(*fdatainp) >> *d;
	return 0;
}

int
main(void)
{
	sym_params pa = { 0 };
	std::string file_name_dati_simulazione,
		file_name_variabili_ottimizzazione;
	std::cin >> file_name_dati_simulazione
		>> file_name_variabili_ottimizzazione;

	std::ifstream fdatain(file_name_dati_simulazione.c_str());
	if (!fdatain) {
		std::cerr << "unable to open file "
			<< file_name_dati_simulazione << std::endl;
		throw;
	}

	pa.scale_eps = 1.;
	pa.scale_f = 1.;

	fdatain >> pa.hi_freq_force_filter_coeff; 
	pa.hi_freq_force_filter_coeff = 1./pa.hi_freq_force_filter_coeff;
	fdatain >> pa.low_freq_displ_filter_coeff;
	fdatain >> pa.static_low_freq_stiffness;

	fdatain >> pa.nsubsteps;

	::fdatainp = &fdatain;

	sym_params *pap = 0;
	if (nlrheo_parse(&pap,
		pa.scale_eps, pa.scale_f,
		pa.hi_freq_force_filter_coeff,
		pa.low_freq_displ_filter_coeff,
		pa.static_low_freq_stiffness,
		pa.nsubsteps))
	{
		return -1;
	}
	::fdatainp = 0;

	pa = *pap;

	std::ifstream fxin(file_name_variabili_ottimizzazione.c_str());
	if (!fxin) {
		std::cerr << "unable to open file "
			<< file_name_variabili_ottimizzazione << std::endl;
		throw;
	}

	pa.x = new double[pa.n_variabili+4];
	for (int i = 0; i<pa.n_variabili+4; i++) {
		fxin >> pa.x[i];
	}
	
	pa.ik_s = new gsl_interp *[pa.n_elementi];
	pa.ik_v = new gsl_interp *[pa.n_elementi];
	pa.ic_s = new gsl_interp *[pa.n_elementi];
	pa.ic_v = new gsl_interp *[pa.n_elementi];
	for (int i = 0; i<pa.n_elementi; i++) {
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

	const gsl_odeiv_step_type *T = gsl_odeiv_step_rkf45;
	gsl_odeiv_step *stepint = gsl_odeiv_step_alloc(T, pa.n_elementi - pa.n_parallelo + 2 + 1);
	gsl_odeiv_evolve * evolve = gsl_odeiv_evolve_alloc(pa.n_elementi - pa.n_parallelo + 2 + 1);
	double eps_abs = 1.E-15;
	double eps_rel = 1.E-15;
	gsl_odeiv_control * control = 
		gsl_odeiv_control_standard_new(eps_abs, eps_rel, 1., 1.);
	gsl_odeiv_system sys;
	sys.function = nlrheo_int_func;
	sys.jacobian = NULL;
	sys.dimension = 1;
	sys.params = &pa;
	pa.y = new double[pa.n_elementi - pa.n_parallelo + 2 + 1];
	double t, dt;
	
	pa.f = 0.;
	for (int i = 0; i < pa.n_elementi - pa.n_parallelo + 2 + 1; i++) {
		pa.y[i] = 0.;
	}

	double final_time = 20.;
	double delta_t = 1.e-3;
	dt = delta_t / pa.nsubsteps;
	for (t = 0.; t < final_time;) {
		pa.tf = t + delta_t;
		pa.ti = t;
		if (t < 5.) {
			pa.sf = 1. * pa.tf;
			pa.si = 1. * pa.ti;
			pa.vf = 1.;
			pa.vi = 1.;
		} else if (t < 10.) {
			pa.sf = 5. - 1. * (pa.tf - 5.);
			pa.si = 5. - 1. * (pa.tf - 5.);
			pa.vf = -1.;
			pa.vi = -1.;
		} else {
			pa.sf = 0.;
			pa.si = 0.;
			pa.vf = 0.;
			pa.vi = 0.;
		}
		gsl_odeiv_evolve_apply(evolve, control, stepint,
			&sys, &t, pa.tf, &dt, pa.y);
		std::cout << t << " " << pa.f << std::endl;
	}
   
	return 0;
}
