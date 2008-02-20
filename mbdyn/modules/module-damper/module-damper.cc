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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <cmath>
#include <cfloat>

#include "dataman.h"
#include "constltp.h"
#include "drive_.h"


extern "C" {
#include <gsl/gsl_errno.h>
#include "gsl/gsl_interp.h"
#include "gsl/gsl_odeiv.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_linalg.h"

int func (double t, const double y[], double f[],
	void *params);


struct sym_params {
	int n_parallelo;
	int *n_serie;			//[n_parallelo];
	int max_n_serie;
	double **s_max;			//[max_n_serie][n_parallelo];
	double **v_max;			//[max_n_serie][n_parallelo];
	int n_elementi;
	int *npti_ks;			//[n_elementi];
	int *npti_kv;			//[n_elementi];
	int *npti_cs;			//[n_elementi];
	int *npti_cv;			//[n_elementi];
	int n_variabili_k, n_variabili_c, n_variabili;
	int max_npti_ks, max_npti_kv, max_npti_cs, max_npti_cv;
	double **k_s;			//[n_elementi][max_npti_ks];
	double **k_v;			//[n_elementi][max_npti_kv];
	double **c_s;			//[n_elementi][max_npti_cs];
	double **c_v;			//[n_elementi][max_npti_cv];
	gsl_interp ** ik_s, ** ik_v, ** ic_s, ** ic_v; // tutti [n_elementi]
	double *x, *y;

	double a;

	// per integrazione
	double sf, si, vf, vi, tf, ti;
	double f, f_s, f_v;

	//di servizio
	gsl_matrix **gsl_C, **gsl_K;
	gsl_vector **gsl_xp, **gsl_x, **gsl_b;
	gsl_permutation **gsl_perm;
};

} //extern "C"

void compute_kc(double &k, double &c, 
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
	const double *const x) {
	//calcola k di el dato s,v
	//std::cerr << "e: " << el << std::endl;
	double k_v_i[npti_kv[el]];
	int somma_k = 0;
	for (int e=0; e < el; e++) {
		somma_k += npti_ks[e] * npti_kv[e];
	}
	if (ik_s != 0) {
		for (int i=0; i < npti_kv[el]; i++) {
			//somma_k += npti_ks[el] * i;
			//std::cerr << "el: " << el << " s: " << s << " k_s[npti_ks[el]-1]: " << k_s[npti_ks[el]-1]<< std::endl;
			if (std::abs(s) > k_s[npti_ks[el]-1]) {
				//std::cerr << "1: " << somma_k + npti_ks[el] - 1 << std::endl;
				k_v_i[i] = x[somma_k + npti_ks[el] - 1];
			} else if (std::abs(s) == 0) {
				k_v_i[i] = x[somma_k];
				//std::cerr << "2: " <<  somma_k << std::endl;
			} else {
				double tx[npti_ks[el]];
				for (int ii=0; ii < npti_ks[el]; ii++) {
					tx[ii] = x[somma_k + ii];
					//std::cerr << "3: " <<  somma_k + ii << " " << " " << k_s[ii] << " " << tx[ii] << std::endl;
				}
				k_v_i[i] = gsl_interp_eval(ik_s, k_s, tx, std::abs(s), NULL);
				//std::cerr << "interp[" << i << "]: " << k_v_i[0] << " s: " << std::abs(s) << std::endl;
			}
			somma_k += npti_ks[el];
		}
	} else {
		for (int i=0; i < npti_kv[el]; i++) {
			//std::cerr << somma_k + i << "x: " << k_v_i[i] << std::endl;
			k_v_i[i] = x[somma_k + i];
		}
	}
	//std::cerr << "<<<<<<<<<<<<<<<<<<\n";
	if (ik_v != 0) {
		for (int i=0; i < npti_kv[el]; i++) {
			//std::cerr << "xx " << i << " " << k_v[i] << " " << k_v_i[i] << " " << std::abs(v) << std::endl;
		}
		//std::cerr << "el: " << el << " v: " << s << " k_v[npti_kv[el]-1]: " << k_v[npti_kv[el]-1]<< std::endl;
		if (std::abs(v) > k_v[npti_kv[el] - 1]) {
			k = k_v_i[npti_kv[el]-1];
			//std::cerr << "4: " << k << " " << std::abs(v) << " " << k_v[npti_kv[el]] << std::endl;
		} else if (std::abs(v) == 0) {
			k = k_v_i[0];
			//std::cerr << "5: " << k << std::endl;
		} else {
			k = gsl_interp_eval(ik_v, k_v, k_v_i, v, NULL);
			//std::cerr << "6: " << k << std::endl;
		}
	} else {
		k = k_v_i[0];
	}
	
	//calcola c di el dato s,v
	double c_v_i[npti_cv[el]];
	//std::cerr << "c: " << n_variabili_k << std::endl;
	int somma_c = n_variabili_k;
	for (int e=0; e < el; e++) {
		somma_c += npti_cs[e] * npti_cv[e];
	}
	if (ic_s != 0) {
		for (int i=0; i < npti_cv[el]; i++) {
			//somma_c += npti_cs[el] * i;
			//std::cerr << "el: " << el << " s: " << s << " c_s[npti_cs[el]-1]: " << c_s[npti_cs[el]-1]<< std::endl;
			if (std::abs(s) > c_s[npti_cs[el] - 1]) {
				//std::cerr << somma_c + npti_cs[el] - 1 << std::endl;
				c_v_i[i] = x[somma_c + npti_cs[el] - 1];
			} else if (std::abs(s) == 0) {
				//std::cerr << somma_c << std::endl;
				c_v_i[i] = x[somma_c];
			}
			double tx[npti_cs[el]];
			for (int ii=0; ii < npti_cs[el]; ii++) {
				//std::cerr << somma_c + ii << std::endl;
				tx[ii] = x[somma_c + ii];
			}
			c_v_i[i] = gsl_interp_eval(ic_s, c_s, tx, std::abs(s), NULL);
			somma_c += npti_cs[el];
		}
	} else {
		for (int i=0; i < npti_cv[el]; i++) {
			//std::cerr << somma_c + i << std::endl;
			c_v_i[i] = x[somma_c + i];
		}
	}
	if (ic_v != 0) {
		//std::cerr << "el: " << el << " v: " << v << " c_v[npti_cv[el]-1]: " << c_v[npti_cv[el]-1]<< std::endl;
		if (std::abs(v) > c_v[npti_cv[el] - 1]) {
			c = c_v_i[npti_cv[el]-1];
		} else if (std::abs(v) == 0) {
			c = c_v_i[0];
		} else {
			c = gsl_interp_eval(ic_v, c_v, c_v_i, v, NULL);
		}
	} else {
		c = c_v_i[0];
	}
	//std::cerr << "--------------------------------" << std::endl;
}


int func (double t, const double y[], double f[],
	void *para) {
	sym_params & pa = *(sym_params *)para;
	double s = (pa.sf - pa.si) / (pa.tf - pa.ti) * (t - pa.ti) + pa.si;
	double v = (pa.vf - pa.vi) / (pa.tf - pa.ti) * (t - pa.ti) + pa.vi;

// 	Per ogni componente in parallelo
	int el = 0;
	int unk = 0;
	pa.f = pa.f_s = pa.f_v = 0.;
	for (int i=0; i < pa.n_parallelo; i++) {
		double c[pa.n_serie[i]], k[pa.n_serie[i]];
		for (int ii = 0; ii < pa.n_serie[i]; ii++) {
			compute_kc(k[ii], c[ii], s, v, el + ii, pa.n_variabili_k, 
				pa.npti_ks, pa.npti_kv, pa.npti_cs, pa.npti_cv,
				pa.ik_s[el+ii], pa.ik_v[el+ii], pa.ic_s[el+ii], pa.ic_v[el+ii], 
				pa.k_s[el+ii], pa.k_v[el+ii], pa.c_s[el+ii], pa.c_v[el+ii], 
				pa.x);
				//std::cerr << "k[" << ii << "] " << k[ii] << "c[" << ii <<"] " << c[ii] << std::endl;
		}
		int nincognite = pa.n_serie[i] - 1;

		if (nincognite > 0) {
			gsl_matrix_set_zero(pa.gsl_C[i]);
			gsl_matrix_set_zero(pa.gsl_K[i]);
			gsl_vector_set_zero(pa.gsl_xp[i]);
			gsl_vector_set_zero(pa.gsl_x[i]);
			gsl_vector_set_zero(pa.gsl_b[i]);
			for (int row = 0; row < nincognite - 1; row++) {
				for (int col = 0; col < nincognite - 1; col++) {
					gsl_matrix_set(pa.gsl_C[i], row, col, c[col]);
					gsl_matrix_set(pa.gsl_C[i], row, col+1, -c[col+1]);
					gsl_matrix_set(pa.gsl_K[i], row, col, k[col]);
					gsl_matrix_set(pa.gsl_K[i], row, col+1, -k[col+1]);
				}
				gsl_matrix_set(pa.gsl_C[i], nincognite, row, c[nincognite]);
				gsl_matrix_set(pa.gsl_K[i], nincognite, row, k[nincognite]);
				gsl_vector_set(pa.gsl_b[i], row, 0.);
				gsl_vector_set(pa.gsl_x[i], row, y[unk + row]);
			}
			gsl_matrix_set(pa.gsl_C[i], nincognite - 1, nincognite - 1, c[nincognite - 1] + c[nincognite]);
			gsl_matrix_set(pa.gsl_K[i], nincognite - 1, nincognite - 1, k[nincognite - 1] + k[nincognite]);
			gsl_vector_set(pa.gsl_b[i], nincognite - 1, c[nincognite] * v + k[nincognite] * s);
			gsl_vector_set(pa.gsl_x[i], nincognite - 1, y[unk + nincognite - 1]);

			//calcola b -= Kx
			gsl_blas_dgemv(CblasNoTrans, -1., pa.gsl_K[i], pa.gsl_x[i], 1., pa.gsl_b[i]);
			//calcola xp = C^-1 b
			int s;
			gsl_linalg_LU_decomp(pa.gsl_C[i], pa.gsl_perm[i], &s);
			gsl_linalg_LU_solve(pa.gsl_C[i], pa.gsl_perm[i], pa.gsl_b[i], pa.gsl_xp[i]);
			//TODO: setta f e par.f
			for (int ii=0; ii < nincognite; ii++) {
				f[unk + ii] = gsl_vector_get(pa.gsl_xp[i], ii);
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
	f[pa.n_elementi - pa.n_parallelo] = -y[pa.n_elementi - pa.n_parallelo] * 2. / pa.a + 
		-y[pa.n_elementi - pa.n_parallelo + 1] / pa.a / pa.a + 
		((pa.x[pa.n_variabili] * std::atan(v / pa.x[pa.n_variabili+1]) +
		pa.x[pa.n_variabili+2] * std::atan(v / pa.x[pa.n_variabili+3]))/2. + pa.f) / pa.a / pa.a;
	f[pa.n_elementi - pa.n_parallelo + 1] = y[pa.n_elementi - pa.n_parallelo];
	pa.f = y[pa.n_elementi - pa.n_parallelo + 1];
// 	pa.f += (pa.x[pa.n_variabili] * std::atan(f[0] / pa.x[pa.n_variabili+1]) +
// 		pa.x[pa.n_variabili+2] * std::atan((v-f[0]) / pa.x[pa.n_variabili+3]))/2.;

/*
	pa.f += (pa.x[pa.n_variabili] * std::atan(v / pa.x[pa.n_variabili+1]) +
		pa.x[pa.n_variabili+2] * std::atan(v / pa.x[pa.n_variabili+3]))/2.;
	pa.f_v += (pa.x[pa.n_variabili] * 1./(1. + std::pow(v / pa.x[pa.n_variabili+1], 2)) / pa.x[pa.n_variabili+1] +
		pa.x[pa.n_variabili+2] * 1./(1. + std::pow(v / pa.x[pa.n_variabili+3], 2)) / pa.x[pa.n_variabili+3])/2.;
*/

	//pa.f += pa.x[pa.n_variabili+4] * std::atan(v / pa.x[pa.n_variabili+5]);

// 	//std::cerr << y[0] << " " << f[0] << "\n";
// 	f[0] = (par.k0 * s + c0 * v - ( k1 + par.k0 ) * y[0] ) / (c1 + c0);
// 	//f[0] = (k1 * s + c0 * v - ( k1 + k1 ) * y[0] ) / (c1 + c0);
// 	par.f = k1 * y[0] + c1 * f[0];
// 	//par.f = k1 * y[0] + c1 * f[0];
// 	par.v1 = par.f;
// 	par.v0 = v - par.v1;
	return GSL_SUCCESS;
}


class DamperConstitutiveLaw
: public ConstitutiveLaw<doublereal, doublereal> {
private:
	sym_params & pa;

	const gsl_odeiv_step_type *T;
	gsl_odeiv_step *stepint;
	gsl_odeiv_evolve * evolve;
	gsl_odeiv_control * control;
	gsl_odeiv_system sys;
	DriveCaller *pTime;
	double *y, *y_dummy;
	double prev_time, current_time, dt;
	double scale_eps, scale_f;
public:
	DamperConstitutiveLaw(sym_params & pap, DriveCaller *pT, double l, double f)
	: pa(pap), T(gsl_odeiv_step_rkf45), pTime(pT),
	prev_time(0.), current_time(0.), dt(0.),
	scale_eps(l), scale_f(f) {
		ConstitutiveLaw<doublereal, doublereal>::FDE =  0.; //FIXME dStiffness;

		stepint = gsl_odeiv_step_alloc(T, pa.n_elementi - pa.n_parallelo + 2);
		evolve = gsl_odeiv_evolve_alloc(pa.n_elementi - pa.n_parallelo + 2);
		double eps_abs = 1.E-15;
		double eps_rel = 1.E-15;
		control = gsl_odeiv_control_standard_new(eps_abs, eps_rel, 1., 1.);
		sys.function = func;
		sys.jacobian = NULL;
		sys.dimension = 1;
		sys.params = &pa;
		y = new double[pa.n_elementi - pa.n_parallelo + 2];
		y_dummy = new double[pa.n_elementi - pa.n_parallelo + 2];

		pa.f = 0.;
		for (int i = 0; i < pa.n_elementi - pa.n_parallelo + 2; i++) {
			y[i] = y_dummy[i] = 0.;
		}
	};

	virtual ~DamperConstitutiveLaw(void) {
		gsl_odeiv_step_free(stepint);
		gsl_odeiv_evolve_free(evolve);
		for (int i=0; i<pa.n_parallelo; i++) {
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
		
		for (int i=0; i<pa.n_elementi; i++) {
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
		
		SAFEDELETE(pTime);
		
		delete[] y;
		delete[] y_dummy;

	};


	virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const {
		silent_cerr("DamperConstitutiveLaw1D::pCopy "
			"not implemented.\n" 
			"Please build explicitly "
			"different instances of the constitutive law" << std::endl);
		throw ErrGeneric();
		return 0;
	};

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOELASTIC;
	};

	virtual void Update(const doublereal& Eps, const doublereal& EpsPrime = 0.) {
		double t = prev_time;
// 		std::cerr << "Update: 1" << std::endl;
// 		std::cerr << "\tEps: " << Eps << std::endl;
// 		std::cerr << "\tEpsPrime: " << EpsPrime << std::endl;

		pa.tf = pTime->dGet();
		pa.ti = prev_time;
// 		std::cerr << "\tpa.tf: " << pa.tf << std::endl;
// 		std::cerr << "\tpa.ti: " << pa.ti << std::endl;

		dt = (pa.tf - pa.ti) / 10.;

		pa.sf = Eps * scale_eps;
		pa.si = Epsilon;
		pa.vf = EpsPrime * scale_eps;
		pa.vi = EpsilonPrime;
		
		for (int i=0; i < pa.n_elementi - pa.n_parallelo + 2; i++) {
			y_dummy[i] = y[i];
		}
		if (dt >0.) {
			gsl_odeiv_evolve_apply(evolve, control, stepint,
				&sys, &t, pa.tf, &dt, y_dummy);
// 			std::cerr << "Update: 2" << std::endl;

			F = pa.f * scale_f;
			FDE = pa.f_s * scale_f;
			FDEPrime = pa.f_v * scale_f;
// 			std::cerr << "\tF: " << F << std::endl;
// 			std::cerr << "\tFDE: " << FDE << std::endl;
// 			std::cerr << "\tFDEPrime: " << FDEPrime << std::endl;
		}
	};

	virtual void AfterConvergence(const doublereal& Eps, const doublereal& EpsPrime = 0.) {
		double t = prev_time;
// 		std::cerr << "AfterConvergence: 1" << std::endl;

		pa.tf = pTime->dGet();
		pa.ti = prev_time;

		dt = pa.tf - pa.ti;

		pa.sf = Eps * scale_eps;
		pa.si = Epsilon;
		pa.vf = EpsPrime * scale_eps;
		pa.vi = EpsilonPrime;
		
		if (dt >0.) {
			gsl_odeiv_evolve_apply(evolve, control, stepint,
				&sys, &t, pa.tf, &dt, y);
// 			std::cerr << "AfterConvergence: 2" << std::endl;

			current_time = pa.tf;
			prev_time = current_time;

			F = pa.f * scale_f;
			FDE = pa.f_s * scale_f;
			FDEPrime = pa.f_v * scale_f;
			Epsilon = Eps * scale_eps;
			EpsilonPrime = EpsPrime * scale_eps;
		}
	};
   
	virtual const doublereal& GetF(void) const {
// 		std::cerr << "F: " << F << std::endl;
		return F;
	};

	virtual const doublereal& GetFDE(void) const {
// 		std::cerr << "FDE: " << FDE << std::endl;
		return FDE;
	};

	virtual const doublereal& GetFDEPrime(void) const {
// 		std::cerr << "FDEPrime: " << FDEPrime << std::endl;
		return FDEPrime;
	};

};

/* specific functional object(s) */
struct DamperCLR : public ConstitutiveLawRead<doublereal, doublereal> {
	virtual ConstitutiveLaw<doublereal, doublereal> *
	Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) {
		ConstitutiveLaw<doublereal, doublereal>* pCL = 0;
		
		
		CLType = ConstLawType::VISCOELASTIC;
		
		double scale_eps = 1;
		if (HP.IsKeyWord("scale_eps")) {
			scale_eps = HP.GetReal();
			// check?
		}

		double scale_f = 1;
		if (HP.IsKeyWord("scale_f")) {
			scale_f = HP.GetReal();
			// check?
		}

		sym_params* pap = new sym_params();
		sym_params &pa(*pap);

		// FIXME: what does it mean?
		pa.a = 0.002;
		if (HP.IsKeyWord("filter")) {
			pa.a = HP.GetReal();
			// check?
		}

		pa.n_parallelo = HP.GetInt();
		pa.gsl_C = new gsl_matrix*[pa.n_parallelo];
		pa.gsl_K = new gsl_matrix*[pa.n_parallelo];
		pa.gsl_xp = new gsl_vector*[pa.n_parallelo];
		pa.gsl_x = new gsl_vector*[pa.n_parallelo];
		pa.gsl_b = new gsl_vector*[pa.n_parallelo];
		pa.gsl_perm = new gsl_permutation*[pa.n_parallelo];
		pa.n_serie = new int[pa.n_parallelo];
		for (int i=0; i<pa.n_parallelo; i++) {
			pa.n_serie[i] = HP.GetInt();
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
		pa.max_n_serie = HP.GetInt();
	
		pa.s_max = new double*[pa.max_n_serie];
		pa.s_max[0] = new double[pa.max_n_serie*pa.n_parallelo];
		for (int i=0; i < pa.max_n_serie; i++) {
			pa.s_max[i] = pa.s_max[0] + i * pa.n_parallelo;
			for (int j=0; j<pa.n_parallelo; j++) {
				pa.s_max[i][j] = HP.GetReal();
			}
		}
	
		pa.v_max = new double*[pa.max_n_serie];
		pa.v_max[0] = new double[pa.max_n_serie*pa.n_parallelo];
		for (int i=0; i < pa.max_n_serie; i++) {
			pa.v_max[i] = pa.v_max[0] + i * pa.n_parallelo;
			for (int j=0; j<pa.n_parallelo; j++) {
				pa.v_max[i][j] = HP.GetReal();
			}
		}
	
		pa.n_elementi = HP.GetInt();
	
		pa.npti_ks = new int[pa.n_elementi];
		for (int i=0; i<pa.n_elementi; i++) {
			pa.npti_ks[i] = HP.GetInt();
		}
	
		pa.npti_kv = new int[pa.n_elementi];
		for (int i=0; i<pa.n_elementi; i++) {
			pa.npti_kv[i] = HP.GetInt();
		}
	
		pa.npti_cs = new int[pa.n_elementi];
		for (int i=0; i<pa.n_elementi; i++) {
			pa.npti_cs[i] = HP.GetInt();
		}
	
		pa.npti_cv = new int[pa.n_elementi];
		for (int i=0; i<pa.n_elementi; i++) {
			pa.npti_cv[i] = HP.GetInt();
		}
	
		pa.n_variabili_k = HP.GetInt();
		pa.n_variabili_c = HP.GetInt();
		pa.n_variabili = HP.GetInt();
	
		pa.max_npti_ks = HP.GetInt();
		pa.max_npti_kv = HP.GetInt();
		pa.max_npti_cs = HP.GetInt();
		pa.max_npti_cv = HP.GetInt();
	
		pa.k_s = new double*[pa.n_elementi];
		pa.k_s[0] = new double[pa.n_elementi*pa.max_npti_ks];
		for (int i=0; i < pa.n_elementi; i++) {
			pa.k_s[i] = pa.k_s[0] + i * pa.max_npti_ks;
			for (int j=0; j < pa.max_npti_ks; j++) {
				pa.k_s[i][j] = HP.GetReal();
// 				std::cerr << "letto k_s[" << i << "][" << j << "] = " << pa.k_s[i][j] << "\n";
			}
		}
		pa.k_v = new double*[pa.n_elementi];
		pa.k_v[0] = new double[pa.n_elementi*pa.max_npti_kv];
		for (int i=0; i < pa.n_elementi; i++) {
			pa.k_v[i] = pa.k_v[0] + i * pa.max_npti_kv;
			for (int j=0; j < pa.max_npti_kv; j++) {
				pa.k_v[i][j] = HP.GetReal();
// 				std::cerr << "letto k_v[" << i << "][" << j << "] = " << pa.k_v[i][j] << "\n";
			}
		}
		pa.c_s = new double*[pa.n_elementi];
		pa.c_s[0] = new double[pa.n_elementi*pa.max_npti_cs];
		for (int i=0; i < pa.n_elementi; i++) {
			pa.c_s[i] = pa.c_s[0] + i * pa.max_npti_cs;
			for (int j=0; j < pa.max_npti_cs; j++) {
				pa.c_s[i][j] = HP.GetReal();
// 				std::cerr << "letto c_s[" << i << "][" << j << "] = " << pa.c_s[i][j] << "\n";
			}
		}
		pa.c_v = new double*[pa.n_elementi];
		pa.c_v[0] = new double[pa.n_elementi*pa.max_npti_cv];
		for (int i=0; i < pa.n_elementi ; i++) {
			pa.c_v[i] = pa.c_v[0] + i * pa.max_npti_cv;
			for (int j=0; j < pa.max_npti_cv; j++) {
				pa.c_v[i][j] = HP.GetReal();
// 				std::cerr << "letto c_v[" << i << "][" << j << "] = " << pa.c_v[i][j] << "\n";
			}
		}


		pa.x = new double[pa.n_variabili+4];
		for (int i=0; i<pa.n_variabili+4; i++) {
			pa.x[i] = HP.GetReal();
		}
	
		pa.ik_s = new gsl_interp *[pa.n_elementi];
		pa.ik_v = new gsl_interp *[pa.n_elementi];
		pa.ic_s = new gsl_interp *[pa.n_elementi];
		pa.ic_v = new gsl_interp *[pa.n_elementi];
		for (int i=0; i<pa.n_elementi; i++) {
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

		DriveCaller *pT = 0;
		SAFENEWWITHCONSTRUCTOR(pT,
			TimeDriveCaller,
			TimeDriveCaller(pDM->pGetDrvHdl()));
		SAFENEWWITHCONSTRUCTOR(pCL, 
			DamperConstitutiveLaw, 
			DamperConstitutiveLaw(pa, pT, scale_eps, scale_f));

		return pCL;
	};
};

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
#if 0
	DataManager	*pDM = (DataManager *)pdm;
	MBDynParser	*pHP = (MBDynParser *)php;
#endif

	ConstitutiveLawRead<doublereal, doublereal> *rf1D
		= new DamperCLR;
	if (!SetCL1D("damper", rf1D)) {
		delete rf1D;

		silent_cerr("DamperConstitutiveLaw: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}


