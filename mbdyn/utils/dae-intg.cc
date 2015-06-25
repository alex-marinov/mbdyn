/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#ifdef USE_RUNTIME_LOADING

#include <cstring>
#include <cstdio>
#include <ltdl.h>

#include "ac/getopt.h"
#include "myassert.h"
#include "solman.h"
#include "fullmh.h"
#include "linsol.h"

#include "dae-intg.h"

enum method_t {
	METHOD_UNKNOWN,
	METHOD_MULTISTEP,
	METHOD_HOPE,
	METHOD_CUBIC,
	METHOD_RADAU_II,
	METHOD_CRANK_NICOLSON
};

struct integration_data {
   	doublereal ti;
   	doublereal tf;
   	doublereal dt;
   	doublereal tol;
   	int maxiter;
   	doublereal rho;
};

static int method_multistep(const char*, integration_data*, void*, const char*);
static int method_hope(const char*, integration_data*, void*, const char*);
static int method_cubic(const char*, integration_data*, void*, const char*);
static int method_radau_II(const char*, integration_data*, void*, const char*);
static int method_cn(const char*, integration_data*, void*, const char*);

static void* get_method_data(method_t, const char*);

static int open_module(const char* module);

static struct funcs *ff = NULL;
static bool print_help = false;
static void* p_data = NULL;

int
main(int argc, char *const argv[])
{
   	struct s_method {
      		const char* s;
      		method_t m;
   	} s_m[] = {
		{ "ms",                METHOD_MULTISTEP         },
		{ "hope",              METHOD_HOPE              },
		{ "cubic",             METHOD_CUBIC             },
		{ "radau-II",	       METHOD_RADAU_II		},
		{ "crank-nicolson",    METHOD_CRANK_NICOLSON    },

		{ NULL,                METHOD_UNKNOWN           }
   	};
   	method_t curr_method = METHOD_UNKNOWN;
   	const char* module = "intg-default.so";
   	char* user_defined = NULL;
   	void* method_data = NULL;
   	integration_data d = { 0., 1., 1.e-3, 1.e-6, 10 };

   	/* opzioni */
   	const char optstring[] = "i:lm:t:T:n:r:u:hH";
   	const struct option longopts[] = {
		{ "integrator",     required_argument, NULL, int('i') },
		{ "method-data",    required_argument, NULL, int('m') },
		{ "timing",         required_argument, NULL, int('t') },
		{ "tolerance",      required_argument, NULL, int('T') },
		{ "iterations",     required_argument, NULL, int('n') },
		{ "rho",            required_argument, NULL, int('r') },
		{ "user-defined",   required_argument, NULL, int('u') },
		{ "help",           no_argument,       NULL, int('h') },
		{ "print-help",     no_argument,       NULL, int('H') },

		{ NULL,             no_argument,       NULL, int(0)   }
   	};

   	while (true) {
		int curropt;

		curropt = getopt_long(argc, argv, optstring, longopts, NULL);

		if (curropt == EOF) {
			break;
		}

      		switch(curropt) {
       		case int('?'):
	  		/*
			 * std::cerr << "unknown option -" << char(optopt) << std::endl;
			 */
	  		break;

       		case int(':'):
	  		std::cerr << "missing parameter for option -"
				<< optopt << std::endl;
	  		break;

       		case int('h'):
	  		std::cerr << "usage: int -[imtTnruh] [module]"
				<< std::endl
	    			<< std::endl
	    			<< " -i,--integrator   : integration method"
				<< std::endl
				<< " -l                : list available methods"
				<< std::endl
	    			<< " -m,--method-data  : method-dependent data"
				<< std::endl
	    			<< " -t,--timing       : ti:dt:tf" << std::endl
	    			<< " -T,--tolerance" << std::endl
	    			<< " -n,--maxiter" << std::endl
	    			<< " -r,--rho          : asymptotic radius"
				<< std::endl
	    			<< " -u,--user         : user-defined data"
				<< std::endl
	    			<< std::endl
	    			<< " -h,--help         : print this message "
					"and exit" << std::endl

	    			<< " -H,--print-help   : print module's "
					"help message" << std::endl;
	  		exit(EXIT_SUCCESS);

		case int('H'):
			print_help = true;
			break;

       		case int('i'): {
	  		s_method* p = s_m;
	  		while (p->s != NULL) {
	     			if (strcmp(p->s, optarg) == 0) {
					curr_method = p->m;
					break;
	     			}
	     			p++;
	  		}
	  		if (p->s == NULL) {
	     			std::cerr << "unknown integrator "
					<< optarg << std::endl;
	     			exit(EXIT_FAILURE);
	  		}
	  		break;
       		}

		case 'l':
			for (int i = 0; s_m[i].s; i++) {
				std::cout << "    " << s_m[i].s << std::endl;
			}
	  		break;

       		case int('m'):
	  		method_data = get_method_data(curr_method, optarg);
	  		break;

       		case int('t'): {
			char	*next = optarg;

			d.ti = strtod(next, &next);
			if (next[0] != ':') {
	     			std::cerr << "syntax: ti:dt:tf" << std::endl;
	     			exit(EXIT_FAILURE);
	  		}

			next++;
			d.dt = strtod(next, &next);
			if (next[0] != ':') {
	     			std::cerr << "syntax: ti:dt:tf" << std::endl;
	     			exit(EXIT_FAILURE);
	  		}

			next++;
			d.tf = strtod(next, &next);
			if (next[0] != '\0') {
	     			std::cerr << "syntax: ti:dt:tf" << std::endl;
	     			exit(EXIT_FAILURE);
	  		}

	  		break;
       		}

       		case int ('T'):
	  		d.tol = strtod(optarg, NULL);
	  		break;

       		case int ('n'):
	  		d.maxiter = strtol(optarg, NULL, 10);
	  		break;

       		case int ('r'):
	  		d.rho = strtod(optarg, NULL);
	  		break;

       		case int('u'):
	  		user_defined = optarg;
	  		break;

       		default:
	  		/* std::cerr << "unknown option" << std::endl; */
	  		break;
      		}
   	}

   	if (optind >= argc) {
		std::cerr << "need a module" << std::endl;
		exit(EXIT_FAILURE);
	}

	module = argv[optind];

	open_module(module);

	if (::ff->size == 0) {
		std::cerr << "no \"size\" func in module" << std::endl;
		exit(EXIT_FAILURE);
	}

	if (::ff->func == 0) {
		std::cerr << "no \"func\" func in module" << std::endl;
		exit(EXIT_FAILURE);
	}

	if (::ff->grad == 0) {
		std::cerr << "no \"grad\" func in module" << std::endl;
		exit(EXIT_FAILURE);
	}

   	if (::ff->read && ::ff->read(&p_data, user_defined)) {
		std::cerr << "unable to read(" << user_defined << ") "
			"data for module \"" << module << "\"" << std::endl;
		exit(EXIT_FAILURE);
   	}

	if (print_help) {
		if (::ff->help) {
			::ff->help(p_data, std::cerr);
		} else {
			std::cerr << "help not available for module "
				"\"" << module << "\" (ignored)" << std::endl;
		}
	}

   	int rc = 0;
   	switch (curr_method) {
    	case METHOD_MULTISTEP :
      		rc = method_multistep(module, &d, method_data, user_defined);
      		break;
    	case METHOD_HOPE:
      		rc = method_hope(module, &d, method_data, user_defined);
      		break;
    	case METHOD_CUBIC:
      		rc = method_cubic(module, &d, method_data, user_defined);
      		break;
    	case METHOD_RADAU_II:
      		rc = method_radau_II(module, &d, method_data, user_defined);
      		break;
    	case METHOD_CRANK_NICOLSON:
      		rc = method_cn(module, &d, method_data, user_defined);
      		break;
    	default:
      		std::cerr << "you must select an integration method" << std::endl;
      		exit(EXIT_FAILURE);
   	}

	if (lt_dlexit()) {
		std::cerr << "lt_dlexit() failed" << std::endl;
		exit(EXIT_FAILURE);
	}

   	return 0;
}

static void*
get_method_data(method_t curr_method, const char* optarg)
{
   	switch (curr_method) {
	case METHOD_CUBIC:
	case METHOD_RADAU_II:
		if (strcasecmp(optarg, "linear") == 0) {
			bool *pi = new bool;
			*pi = true;
			return pi;

		} else if (strcasecmp(optarg, "cubic") == 0) {
			bool *pi = new bool;
			*pi = false;
			return pi;

		} else {
      			std::cerr << "unknown data \"" << optarg << "\""
				<< std::endl;
			return NULL;
		}
		break;

    	default:
      		std::cerr << "not implemented yet" << std::endl;
      		exit(EXIT_FAILURE);
   	}

   	return NULL;
}

static int
open_module(const char* module)
{
   	const char* err = NULL;

	lt_dlhandle handle;

	if (lt_dlinit()) {
		std::cerr << "lt_dlinit() failed" << std::endl;
      		exit(EXIT_FAILURE);
   	}

   	if ((handle = lt_dlopen(module)) == NULL) {
      		err = lt_dlerror();
      		std::cerr << "lt_dlopen(\"" << module
			<< "\") returned \"" << err << "\"" << std::endl;
      		exit(EXIT_FAILURE);
   	}

	struct funcs **pf = (struct funcs **)lt_dlsym(handle, "ff");
   	if (pf == NULL) {
      		err = lt_dlerror();
      		std::cerr << "lt_dlsym(\"ff\") returned \"" << err << "\""
			<< std::endl;
      		exit(EXIT_FAILURE);
   	}

	::ff = *pf;
	if (::ff == NULL) {
		std::cerr << "invalid \"ff\" symbol in \"" << module << "\""
			<< std::endl;
      		exit(EXIT_FAILURE);
   	}

   	return 0;
}

static void
flip(VectorHandler** ppX, VectorHandler** ppXP,
	VectorHandler** ppXm1, VectorHandler** ppXPm1,
	VectorHandler** ppXm2, VectorHandler** ppXPm2)
{
   	VectorHandler* p = *ppXm2;
   	*ppXm2 = *ppXm1;
   	*ppXm1 = *ppX;
   	*ppX = p;

   	p = *ppXPm2;
   	*ppXPm2 = *ppXPm1;
   	*ppXPm1 = *ppXP;
   	*ppXP = p;
}

static doublereal
cm0(const doublereal& xi)
{
	return 1 - 3*xi*xi - 2*xi*xi*xi;
}

static doublereal
cm0p(const doublereal& xi)
{
	return - 6.*xi - 6.*xi*xi;
}

static doublereal
cm1(const doublereal& xi)
{
	return 3.*xi*xi + 2.*xi*xi*xi;
}

static doublereal
cm1p(const doublereal& xi)
{
	return 6.*xi + 6.*xi*xi;
}

static doublereal
cn0(const doublereal& xi)
{
	return xi + 2.*xi*xi + xi*xi*xi;
}

static doublereal
cn0p(const doublereal& xi)
{
	return 1. + 4.*xi + 3.*xi*xi;
}

static doublereal
cn1(const doublereal& xi)
{
	return xi*xi + xi*xi*xi;
}

static doublereal
cn1p(const doublereal& xi)
{
	return 2.*xi + 3.*xi*xi;
}

static int
method_multistep(const char* module, integration_data* d,
		 void* method_data, const char* user_defined)
{
   	/* prepara le strutture dati per il calcolo */
   	int size = ::ff->size(p_data);
   	MyVectorHandler v0(size);
   	MyVectorHandler v1(size);
   	MyVectorHandler v2(size);
   	MyVectorHandler v3(size);
   	MyVectorHandler v4(size);
   	MyVectorHandler v5(size);

   	VectorHandler* pX = &v0;
   	VectorHandler* pXP = &v1;
   	VectorHandler* pXm1 = &v2;
   	VectorHandler* pXPm1 = &v3;
   	VectorHandler* pXm2 = &v4;
   	VectorHandler* pXPm2 = &v5;
   	pX->Reset();
   	pXP->Reset();
   	pXm1->Reset();
   	pXPm1->Reset();
   	pXm2->Reset();
   	pXPm2->Reset();

   	doublereal* pd = new doublereal[2 * size*size];
   	doublereal** ppd = new doublereal*[2 * size];
   	FullMatrixHandler J(pd, ppd, size*size, size, size);
	FullMatrixHandler JP(pd+size*size, ppd+size, size*size, size, size);
   	MyVectorHandler R(size);

	// TODO: abstract from a specific solution manager?
	LinSol ls;
	SolutionManager *sm = ls.GetSolutionManager(size);

   	MatrixHandler* Jac = sm->pMatHdl();
   	VectorHandler* Res = sm->pResHdl();
   	VectorHandler* Sol = sm->pSolHdl();

   	/* parametri di integrazione */
   	doublereal ti = d->ti;
   	doublereal dt = d->dt;
   	doublereal tf = d->tf;

   	doublereal tol = d->tol;
   	int maxiter = d->maxiter;

   	/* coefficienti del metodo */
   	doublereal rho = d->rho;
   	doublereal beta =
		(4.*rho*rho-(1.-rho)*(1.-rho))/(4.-(1.-rho)*(1.-rho));
   	doublereal delta = (1.-rho)*(1.-rho)/(2.*(4.-(1.-rho)*(1.-rho)));
   	doublereal a1 = 1.-beta;
   	doublereal a2 = beta;
   	doublereal b0 = delta+.5;
   	doublereal b1 = .5*beta+.5-2.*delta;
   	doublereal b2 = .5*beta+delta;

   	doublereal t = ti;

   	/* inizializza la soluzione */
	if (::ff->init) {
   		::ff->init(p_data, *pX, *pXP);
	}
   	for (int k = 1; k <= size; k++) {
      		doublereal x = (*pX)(k);
      		doublereal xp = (*pXP)(k);
      		(*pXPm1)(k) = xp;
      		(*pXm1)(k) = x-dt*xp;
   	}

   	/* output iniziale */
   	std::cout << ti << " " << 0. << " ";
	if (::ff->out) {
   		::ff->out(p_data, std::cout, *pX, *pXP) << std::endl;
	}

   	flip(&pX, &pXP, &pXm1, &pXPm1, &pXm2, &pXPm2);

   	while (t < tf) {
      		t += dt;

      		/* predict */
      		for (int ir = size; ir > 0; ir--) {
	 		doublereal xm1 = (*pXm1)(ir);
	 		doublereal xPm1 = (*pXPm1)(ir);
	 		doublereal xm2 = (*pXm2)(ir);
	 		doublereal xPm2 = (*pXPm2)(ir);
	 		doublereal xP = cm0p(1.)*xm1
				+ cm1p(1.)*xm2
				+ dt*(cn0p(1.)*xPm1 + cn1p(1.)*xPm2);
			doublereal x = a1*xm1
				+ a2*xm2
				+ dt*(b0*xP + b1*xPm1 + b2*xPm2);
	 		(*pX)(ir) = x;
	 		(*pXP)(ir) = xP;
      		}

      		/* test */
      		int j = 0;
      		doublereal test;
      		doublereal coef = dt*b0;
      		do {
			Jac = sm->pMatHdl();
			Res = sm->pResHdl();
			Sol = sm->pSolHdl();

			::ff->func(p_data, *Res, *pX, *pXP, t);

	 		test = Res->Norm();
			std::cerr << j << " " << test << std::endl;
	 		if (test < tol) {
	    			break;
	 		}
	 		if (++j > maxiter) {
	    			std::cerr << "current iteration " << j
	      				<< " exceeds max iteration number "
					<< maxiter << std::endl;
	    			exit(EXIT_FAILURE);
	 		}

	 		/* correct */
	 		sm->MatrReset();
	 		J.Reset();
	 		JP.Reset();
	 		::ff->grad(p_data, J, JP, *pX, *pXP, t);
	 		for (int ir = 1; ir <= size; ir++) {
	    			for (int ic = 1; ic <= size; ic++) {
	       				(*Jac)(ir, ic) = JP(ir, ic)
						+ coef*J(ir, ic);
	    			}
	 		}
	 		sm->Solve();

	 		/* update */
	 		for (int ir = size; ir > 0; ir--) {
	    			doublereal dxP = (*Sol)(ir);
	    			(*pXP)(ir) += dxP;
	    			(*pX)(ir) += coef*dxP;
	 		}
      		} while (true);

      		/* output */
      		std::cout << t << " " << test << " ";
		if (::ff->out) {
      			::ff->out(p_data, std::cout, *pX, *pXP) << std::endl;
		}

      		flip(&pX, &pXP, &pXm1, &pXPm1, &pXm2, &pXPm2);
	}

	if (::ff->destroy) {
   		::ff->destroy(&p_data);
	}

   	delete[] pd;
   	delete[] ppd;

   	return 0;
}

static int
method_hope(const char* module, integration_data* d,
	    void* method_data, const char* user_defined)
{
   	std::cerr << __FUNCTION__ << "not implemented yet!" << std::endl;
   	exit(EXIT_FAILURE);

   	return 0;
}

static int
method_cubic(const char* module, integration_data* d,
	     void* method_data, const char* user_defined)
{
	bool linear = false;

	if (method_data) {
		bool *pi = (bool *)method_data;

		linear = *pi;
	}

   	/* prepara le strutture dati per il calcolo */
   	int size = ::ff->size(p_data);
   	MyVectorHandler v0(size);
   	MyVectorHandler v1(size);
   	MyVectorHandler v2(size);
   	MyVectorHandler v3(size);
   	MyVectorHandler v4(size);
   	MyVectorHandler v5(size);

   	VectorHandler* pX = &v0;
   	VectorHandler* pXP = &v1;
   	VectorHandler* pXm1 = &v2;
   	VectorHandler* pXPm1 = &v3;
   	VectorHandler* pXm2 = &v4;
   	VectorHandler* pXPm2 = &v5;
   	pX->Reset();
   	pXP->Reset();
   	pXm1->Reset();
   	pXPm1->Reset();
   	pXm2->Reset();
   	pXPm2->Reset();

   	doublereal* pd = new doublereal[4*size*size];
   	doublereal** ppd = new doublereal*[4*size];
   	FullMatrixHandler Jz(pd, ppd, size*size, size, size);
   	FullMatrixHandler J0(pd+size*size, ppd+size, size*size, size, size);
	FullMatrixHandler JPz(pd+2*size*size, ppd+2*size,
			size*size, size, size);
	FullMatrixHandler JP0(pd+3*size*size, ppd+3*size,
			size*size, size, size);
   	MyVectorHandler Xz(size);
   	MyVectorHandler XPz(size);

	// TODO: abstract from a specific solution manager?
	LinSol ls;
	SolutionManager *sm = ls.GetSolutionManager(2*size);

   	MatrixHandler* Jac = sm->pMatHdl();
   	VectorHandler* Res = sm->pResHdl();
   	VectorHandler* Sol = sm->pSolHdl();

   	MyVectorHandler Resz(size, Res->pdGetVec() + size);

   	/* parametri di integrazione */
   	doublereal ti = d->ti;
   	doublereal dt = d->dt;
   	doublereal tf = d->tf;

   	doublereal tol = d->tol;
   	int maxiter = d->maxiter;

   	/* coefficienti del metodo */
   	doublereal rho = d->rho;
   	doublereal z = -rho/(1.+rho);
   	doublereal w1 = (2.+3.*z)/(6.*(1.+z));
   	doublereal wz = -1./(6.*z*(1.+z));
   	doublereal w0 = (1.+3.*z)/(6.*z);

	doublereal jzz = (1.+3.*rho)/(6.*rho*(1.+rho))*dt;
	doublereal jz0 = -1./(6.*rho*(1.+rho)*(1.+rho))*dt;
	doublereal j0z = (1.+rho)*(1.+rho)/(6.*rho)*dt;
	doublereal j00 = (2.*rho-1.)/(6.*rho)*dt;

   	doublereal t = ti;

   	/* inizializza la soluzione */
	if (::ff->init) {
   		::ff->init(p_data, *pX, *pXP);
	}
   	::ff->func(p_data, *Res, *pX, *pXP, t);
   	for (int k = 1; k <= size; k++) {
      		doublereal x = (*pX)(k);
      		doublereal xp = (*pXP)(k);
      		(*pXPm1)(k) = xp;
      		(*pXm1)(k) = x - dt*xp;
   	}

   	/* output iniziale */
   	std::cout << ti << " " << 0. << " ";
	if (::ff->out) {
   		::ff->out(p_data, std::cout, *pX, *pXP) << std::endl;
	}

   	flip(&pX, &pXP, &pXm1, &pXPm1, &pXm2, &pXPm2);

   	while (t < tf) {
      		t += dt;

		if (linear) {

			std::cerr << "linear predictor" << std::endl;

	      		/* predict lineare */
	      		for (int k = 1; k <= size; k++) {
		 		doublereal xm1 = (*pXm1)(k);
		 		doublereal xPm1 = (*pXPm1)(k);
		 		doublereal xP = xPm1;
		 		doublereal xPz = xPm1;
		 		doublereal x = xm1 + dt*xP;
		 		doublereal xz = xm1 + dt*(1.+z)*xP;

		 		(*pX)(k) = x;
		 		(*pXP)(k) = xP;
				Xz(k) = xz;
				XPz(k) = xPz;
      			}

		} else {

      			/* predict cubico */
      			for (int k = 1; k <= size; k++) {
		 		doublereal xm1 = (*pXm1)(k);
		 		doublereal xPm1 = (*pXPm1)(k);
		 		doublereal xm2 = (*pXm2)(k);
		 		doublereal xPm2 = (*pXPm2)(k);
		 		doublereal xP = (cm0p(1.)*xm1 + cm1p(1.)*xm2)/dt
					+ cn0p(1.)*xPm1 + cn1p(1.)*xPm2;
		 		doublereal xPz = (cm0p(1. + z)*xm1 + cm1p(1. + z)*xm2)/dt
					+ cn0p(1. + z)*xPm1 + cn1p(1. + z)*xPm2;
		 		doublereal x = xm1 + dt*(w1*xPm1 + wz*xPz + w0*xP);
		 		doublereal xz = cm0(z)*x + cm1(z)*xm1 + dt*(cn0(z)*xP + cn1(z)*xPm1);

		 		(*pX)(k) = x;
		 		(*pXP)(k) = xP;
				Xz(k) = xz;
				XPz(k) = xPz;
			}
      		}


      		/* test */
      		int j = 0;
      		doublereal test;
      		do {

			Jac = sm->pMatHdl();
			Res = sm->pResHdl();
			Sol = sm->pSolHdl();

			Resz.Attach(size, Res->pdGetVec() + size);

	 		Res->Reset();
	 		::ff->func(p_data, *Res, *pX, *pXP, t);
			/* Res->Reset() zeros Resz as well */
	 		::ff->func(p_data, Resz, Xz, XPz, t + z*dt);

			/* Res->Norm() computes Resz norm as well */
	 		test = Res->Norm();
	 		if (test < tol) {
	    			break;
	 		}
	 		if (++j > maxiter) {
	    			std::cerr << "current iteration " << j
	      				<< " exceeds max iteration number "
					<< maxiter << std::endl;
	    			exit(EXIT_FAILURE);
	 		}

	 		/* correct */
	 		sm->MatrReset();
	 		Jz.Reset();
	 		JPz.Reset();
	 		J0.Reset();
	 		JP0.Reset();
	 		::ff->grad(p_data, Jz, JPz, Xz, XPz, t + z*dt);
	 		::ff->grad(p_data, J0, JP0, *pX, *pXP, t);

			for (int ir = 1; ir <= size; ir++) {
				for (int ic = 1; ic <= size; ic++) {
					(*Jac)(ir, ic)
						= j00*J0(ir, ic)
						+ JP0(ir, ic);
					(*Jac)(ir, size + ic)
						= j0z*J0(ir, ic);
					(*Jac)(size + ir, ic)
						= jz0*Jz(ir, ic);
					(*Jac)(size + ir, size + ic)
						= jzz*Jz(ir, ic)
						+ JPz(ir, ic);
				}
			}

	 		sm->Solve();

	 		/* update */
	 		for (int ir = size; ir > 0; ir--) {
	    			doublereal dxP0 = (*Sol)(ir);
	    			doublereal dxPz = (*Sol)(size + ir);
	    			(*pXP)(ir) += dxP0;
				XPz(ir) += dxPz;
	    			(*pX)(ir) += dt*(wz*dxPz + w0*dxP0);
				Xz(ir) += dt*(cm0(z)*(wz*dxPz + w0*dxP0)+cn0(z)*dxP0);
	 		}
     	 	} while (true);

      		/* output */
      		std::cout << t << " " << test << " ";
		if (::ff->out) {
      			::ff->out(p_data, std::cout, *pX, *pXP) << std::endl;
		}

      		flip(&pX, &pXP, &pXm1, &pXPm1, &pXm2, &pXPm2);
   	}

	if (::ff->destroy) {
   		::ff->destroy(&p_data);
	}
   	delete[] pd;
   	delete[] ppd;

   	return 0;
}

static doublereal
radau_II_cm0(const doublereal& xi)
{
	return 1 - xi*xi;
}

static doublereal
radau_II_cm0p(const doublereal& xi)
{
	return - 2.*xi;
}

static doublereal
radau_II_cm1(const doublereal& xi)
{
	return xi*xi;
}

static doublereal
radau_II_cm1p(const doublereal& xi)
{
	return 2.*xi;
}

static doublereal
radau_II_cn0(const doublereal& xi)
{
	return xi + xi*xi;
}

static doublereal
radau_II_cn0p(const doublereal& xi)
{
	return 1. + 2.*xi;
}

static doublereal
radau_II_cn1(const doublereal& xi)
{
	return 0.;
}

static doublereal
radau_II_cn1p(const doublereal& xi)
{
	return 0.;
}

static int
method_radau_II(const char* module, integration_data* d,
	     void* method_data, const char* user_defined)
{
	bool linear = false;

	if (method_data) {
		bool *pi = (bool *)method_data;

		linear = *pi;
	}

   	/* prepara le strutture dati per il calcolo */
   	int size = ::ff->size(p_data);
   	MyVectorHandler v0(size);
   	MyVectorHandler v1(size);
   	MyVectorHandler v2(size);
   	MyVectorHandler v3(size);
   	MyVectorHandler v4(size);
   	MyVectorHandler v5(size);

   	VectorHandler* pX = &v0;
   	VectorHandler* pXP = &v1;
   	VectorHandler* pXm1 = &v2;
   	VectorHandler* pXPm1 = &v3;
   	VectorHandler* pXm2 = &v4;
   	VectorHandler* pXPm2 = &v5;
   	pX->Reset();
   	pXP->Reset();
   	pXm1->Reset();
   	pXPm1->Reset();
   	pXm2->Reset();
   	pXPm2->Reset();

   	doublereal* pd = new doublereal[4*size*size];
   	doublereal** ppd = new doublereal*[4*size];
   	FullMatrixHandler Jz(pd, ppd, size*size, size, size);
   	FullMatrixHandler J0(pd+size*size, ppd+size, size*size, size, size);
	FullMatrixHandler JPz(pd+2*size*size, ppd+2*size,
			size*size, size, size);
	FullMatrixHandler JP0(pd+3*size*size, ppd+3*size,
			size*size, size, size);
   	MyVectorHandler Xz(size);
   	MyVectorHandler XPz(size);

	// TODO: abstract from a specific solution manager?
	LinSol ls;
	SolutionManager *sm = ls.GetSolutionManager(2*size);

   	MatrixHandler* Jac = sm->pMatHdl();
   	VectorHandler* Res = sm->pResHdl();
   	VectorHandler* Sol = sm->pSolHdl();

   	MyVectorHandler Resz(size, Res->pdGetVec() + size);

   	/* parametri di integrazione */
   	doublereal ti = d->ti;
   	doublereal dt = d->dt;
   	doublereal tf = d->tf;

   	doublereal tol = d->tol;
   	int maxiter = d->maxiter;

   	/* coefficienti del metodo */
   	doublereal z = -2./3.;
   	doublereal w1 = 0.;
   	doublereal wz = -1./(2.*z);
   	doublereal w0 = (1.+2.*z)/(2.*z);
   	doublereal m0 = radau_II_cm0(z);
   	doublereal m1 = radau_II_cm1(z);
   	doublereal n0 = radau_II_cn0(z);
   	doublereal n1 = radau_II_cn1(z);

	doublereal jzz = 5./12.*dt;
	doublereal jz0 = -1./12.*dt;
	doublereal j0z = 3./4.*dt;
	doublereal j00 = 1./4.*dt;

   	doublereal t = ti;

   	/* inizializza la soluzione */
	if (::ff->init) {
   		::ff->init(p_data, *pX, *pXP);
	}
   	::ff->func(p_data, *Res, *pX, *pXP, t);
   	for (int k = 1; k <= size; k++) {
      		doublereal x = (*pX)(k);
      		doublereal xp = (*pXP)(k);
      		(*pXPm1)(k) += xp;
      		(*pXm1)(k) += x - dt*xp;
   	}

   	/* output iniziale */
   	std::cout << ti << " " << 0. << " ";
	if (::ff->out) {
   		::ff->out(p_data, std::cout, *pX, *pXP) << std::endl;
	}

   	flip(&pX, &pXP, &pXm1, &pXPm1, &pXm2, &pXPm2);

   	while (t < tf) {
      		t += dt;

		if (linear) {

			std::cerr << "linear predictor" << std::endl;

	      		/* predict lineare */
	      		for (int k = 1; k <= size; k++) {
		 		doublereal xm1 = (*pXm1)(k);
		 		doublereal xPm1 = (*pXPm1)(k);
		 		doublereal xP = xPm1;
		 		doublereal xPz = xPm1;
		 		doublereal x = xm1 + dt*xP;
		 		doublereal xz = xm1 + dt*(1.+z)*xP;

		 		(*pX)(k) = x;
		 		(*pXP)(k) = xP;
				Xz(k) = xz;
				XPz(k) = xPz;
      			}

		} else {

      			/* predict cubico */
      			for (int k = 1; k <= size; k++) {
		 		doublereal xm1 = (*pXm1)(k);
		 		doublereal xPm1 = (*pXPm1)(k);
		 		doublereal xm2 = (*pXm2)(k);
		 		doublereal xPm2 = (*pXPm2)(k);
		 		doublereal xP = (radau_II_cm0p(1.)*xm1 + radau_II_cm1p(1.)*xm2)/dt
					+ radau_II_cn0p(1.)*xPm1 + radau_II_cn1p(1.)*xPm2;
		 		doublereal xPz = (radau_II_cm0p(1. + z)*xm1 + radau_II_cm1p(1. + z)*xm2)/dt
					+ radau_II_cn0p(1. + z)*xPm1 + radau_II_cn1p(1. + z)*xPm2;
		 		doublereal x = xm1 + dt*(w1*xPm1 + wz*xPz + w0*xP);
		 		doublereal xz = m0*x + m1*xm1 + dt*(n0*xP + n1*xPm1);

		 		(*pX)(k) = x;
		 		(*pXP)(k) = xP;
				Xz(k) = xz;
				XPz(k) = xPz;
			}
      		}


      		/* test */
      		int j = 0;
      		doublereal test;
      		do {

			Jac = sm->pMatHdl();
			Res = sm->pResHdl();
			Sol = sm->pSolHdl();

			Resz.Attach(size, Res->pdGetVec() + size);

	 		Res->Reset();
	 		::ff->func(p_data, *Res, *pX, *pXP, t);
			/* Res->Reset() zeros Resz as well */
	 		::ff->func(p_data, Resz, Xz, XPz, t + z*dt);

			/* Res->Norm() computes Resz norm as well */
	 		test = Res->Norm();
	 		if (test < tol) {
	    			break;
	 		}
	 		if (++j > maxiter) {
	    			std::cerr << "current iteration " << j
	      				<< " exceeds max iteration number "
					<< maxiter << std::endl;
	    			exit(EXIT_FAILURE);
	 		}

	 		/* correct */
	 		sm->MatrReset();
	 		Jz.Reset();
	 		JPz.Reset();
	 		J0.Reset();
	 		JP0.Reset();
	 		::ff->grad(p_data, Jz, JPz, Xz, XPz, t + z*dt);
	 		::ff->grad(p_data, J0, JP0, *pX, *pXP, t);

			for (int ir = 1; ir <= size; ir++) {
				for (int ic = 1; ic <= size; ic++) {
					(*Jac)(ir, ic)
						= j00*J0(ir, ic)
						+ JP0(ir, ic);
					(*Jac)(ir, size + ic)
						= j0z*J0(ir, ic);
					(*Jac)(size + ir, ic)
						= jz0*Jz(ir, ic);
					(*Jac)(size + ir, size + ic)
						= jzz*Jz(ir, ic)
						+ JPz(ir, ic);
				}
			}

	 		sm->Solve();

	 		/* update */
	 		for (int ir = size; ir > 0; ir--) {
	    			doublereal dxP0 = (*Sol)(ir);
	    			doublereal dxPz = (*Sol)(size + ir);
	    			(*pXP)(ir) += dxP0;
				XPz(ir) += dxPz;
	    			(*pX)(ir) += dt*(wz*dxPz + w0*dxP0);
				Xz(ir) += dt*(m0*(wz*dxPz + w0*dxP0)+n0*dxP0);
	 		}
     	 	} while (true);

      		/* output */
      		std::cout << t << " " << test << " ";
		if (::ff->out) {
      			::ff->out(p_data, std::cout, *pX, *pXP) << std::endl;
		}

      		flip(&pX, &pXP, &pXm1, &pXPm1, &pXm2, &pXPm2);
	}

	if (::ff->destroy) {
   		::ff->destroy(&p_data);
	}
   	delete[] pd;
   	delete[] ppd;

   	return 0;
}

static int
method_cn(const char* module, integration_data* d,
	  void* method_data, const char* user_defined)
{
   	std::cerr << __FUNCTION__ << "not implemented yet!" << std::endl;
   	exit(EXIT_FAILURE);

   	return 0;
}

#else // ! USE_RUNTIME_LOADING

#include <iostream>
#include <stdlib.h>

int
main(void)
{
	std::cerr << "Need dynamic load capabilities" << std::endl;
   	exit(EXIT_FAILURE);
}

#endif // ! USE_RUNTIME_LOADING

