/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ltdl.h>

#include "ac/getopt.h"
#include <iostream>

#include "myassert.h"
#include "solman.h"
#include "linsol.h"
#include "fullmh.h"
#include "harwrap.h"
#include "mschwrap.h"

#include "intg.h"

struct integration_data {
   	doublereal ti;
   	doublereal tf;
   	doublereal dt;
   	doublereal tol;
   	int maxiter;
   	doublereal rho;
};

static int open_module(const char* module);

static int method_multistep(const char*, integration_data*, void*, const char*);
static int method_hope(const char*, integration_data*, void*, const char*);
static int method_cubic(const char*, integration_data*, void*, const char*);
static int method_cn(const char*, integration_data*, void*, const char*);

void* get_method_data(int, const char*);

static struct funcs *ff = NULL;
static bool print_help = false;
static void* p_data = NULL;

int
main(int argn, char *const argv[])
{
   	enum {
      		METHOD_UNKNOWN,
		METHOD_MULTISTEP,
		METHOD_HOPE,
		METHOD_CUBIC,
		METHOD_CRANK_NICOLSON
   	};

   	struct s_method {
      		const char* s;
      		int m;
   	} s_m[] = {
		{ "ms",                METHOD_MULTISTEP         },
		{ "hope",              METHOD_HOPE              },
		{ "cubic",             METHOD_CUBIC             },
		{ "crank-nicolson",    METHOD_CRANK_NICOLSON    },

		{ 0,                METHOD_UNKNOWN           }
   	};
   	int curr_method = METHOD_UNKNOWN;
   	char* module = 0;
   	char* user_defined = 0;
   	void* method_data = 0;
   	integration_data d = { 0., 1., 1.e-3, 1.e-6, 10 };

   	/* opzioni */
   	const char optstring[] = "i:m:t:T:n:r:u:h";
   	const struct option longopts[] = {
		{ "integrator",		required_argument,	0, int('i') },
		{ "method-data",	required_argument,	0, int('m') },
		{ "timing",		required_argument,	0, int('t') },
		{ "tolerance",		required_argument,	0, int('T') },
		{ "iterations",		required_argument,	0, int('n') },
		{ "rho",		required_argument,	0, int('r') },
		{ "user-defined",	required_argument,	0, int('u') },
		{ "help",		no_argument,		0, int('h') },
		{ "print-help",		no_argument,		0, int('H') },

		{ 0,			no_argument,		0, int(0)   }
   	};

   	while (true) {
		int curropt;

		curropt = getopt_long(argn, argv, optstring, longopts, 0);

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
	  		std::cerr << "usage: int -[imtTnruh] [module]" << std::endl
	    			<< std::endl
	    			<< " -i,--integrator   : integration method" << std::endl
	    			<< " -m,--method-data  : method-dependent data" << std::endl
	    			<< " -t,--timing       : ti:dt:tf" << std::endl
	    			<< " -T,--tolerance" << std::endl
	    			<< " -n,--maxiter" << std::endl
	    			<< " -r,--rho          : asymptotic radius" << std::endl
	    			<< " -u,--user         : user-defined data" << std::endl
	    			<< std::endl
	    			<< " -h,--help         : print this message and exit" << std::endl;
	  		exit(EXIT_SUCCESS);

		case int('H'):
			print_help = true;
			break;

       		case int('i'): {
	  		s_method* p = s_m;
	  		while (p->s != 0) {
	     			if (strcmp(p->s, optarg) == 0) {
					curr_method = p->m;
					break;
	     			}
	     			p++;
	  		}
	  		if (p->s == 0) {
	     			std::cerr << "unknown integrator "
					<< optarg << std::endl;
	     			exit(EXIT_FAILURE);
	  		}
	  		break;
       		}

       		case int('m'):
	  		method_data = get_method_data(curr_method, optarg);
	  		break;

       		case int('t'): {
	  		char* s = new char[strlen(optarg)+1];
	  		strcpy(s, optarg);
	  		char* p = std::strrchr(s, ':');
	  		if (p == 0) {
	     			std::cerr << "syntax: ti:dt:tf" << std::endl;
	     			exit(EXIT_FAILURE);
	  		}
	  		d.tf = atof(p+1);
	  		p[0] = '\0';
	  		p = std::strrchr(s, ':');
	  		if (p == 0) {
	     			std::cerr << "syntax: ti:dt:tf" << std::endl;
	     			exit(EXIT_FAILURE);
	  		}
	  		d.dt = atof(p+1);
	  		p[0] = '\0';
	  		d.ti = atof(s);
	  		delete[] s;
	  		break;
       		}

       		case int ('T'):
	  		d.tol = atof(optarg);
	  		break;

       		case int ('n'):
	  		d.maxiter = atoi(optarg);
	  		break;

       		case int ('r'):
	  		d.rho = atof(optarg);
	  		break;

       		case int('u'):
	  		user_defined = optarg;
	  		break;

       		default:
	  		/* std::cerr << "unknown option" << std::endl; */
	  		break;
      		}
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
    	case METHOD_CRANK_NICOLSON:
      		rc = method_cn(module, &d, method_data, user_defined);
      		break;
    	default:
      		std::cerr << "you must select an integration method" << std::endl;
      		exit(EXIT_FAILURE);
   	}

   	return 0;
}

void *
get_method_data(int curr_method, const char* optarg)
{
	switch (curr_method) {
	default:
		std::cerr << "not implemented yet" << std::endl;
		exit(EXIT_FAILURE);
	}

	return 0;
}

static int
open_module(const char* module)
{
   	const char* err = 0;

	lt_dlhandle handle;

	if (lt_dlinit()) {
		std::cerr << "lt_dlinit() failed" << std::endl;
      		exit(EXIT_FAILURE);
   	}

   	handle = lt_dlopen(module);
	if (handle == 0) {
      		err = lt_dlerror();
      		std::cerr << "lt_dlopen(\"" << module
			<< "\") returned \"" << err << "\"" << std::endl;
      		exit(EXIT_FAILURE);
   	}

	struct funcs **pf = (struct funcs **)lt_dlsym(handle, "ff");
   	if (pf == 0) {
      		err = lt_dlerror();
      		std::cerr << "lt_dlsym(\"ff\") returned \"" << err << "\""
			<< std::endl;
      		exit(EXIT_FAILURE);
   	}

	::ff = *pf;
	if (::ff == 0) {
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

static int
method_multistep(const char* module, integration_data* d,
	void* method_data, const char* user_defined)
{
	// prepara i dati
	void* p_data = 0;
	::ff->read(&p_data, user_defined);

	// prepara le strutture dati per il calcolo
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

	FullMatrixHandler J(size);
	MyVectorHandler R(size);

	// TODO: abstract from a specific solution manager?
	LinSol ls;
	SolutionManager *sm = ls.GetSolutionManager(size);

	MatrixHandler& Jac = *sm->pMatHdl();
	VectorHandler& Res = *sm->pResHdl();
	VectorHandler& Sol = *sm->pSolHdl();

	// parametri di integrazione
	doublereal ti = d->ti;
	doublereal dt = d->dt;
	doublereal tf = d->tf;

	doublereal tol = d->tol;
	int maxiter = d->maxiter;

	// coefficienti del metodo
	doublereal rho = d->rho;
	doublereal beta = (4.*rho*rho-(1.-rho)*(1.-rho))/(4.-(1.-rho)*(1.-rho));
	doublereal delta = (1.-rho)*(1.-rho)/(2.*(4.-(1.-rho)*(1.-rho)));
	doublereal a1 = 1.-beta;
	doublereal a2 = beta;
	doublereal b0 = delta+.5;
	doublereal b1 = .5*beta+.5-2.*delta;
	doublereal b2 = .5*beta+delta;

	doublereal t = ti;

	// inizializza la soluzione
	::ff->init(p_data, *pX);
	::ff->func(p_data, *pXP, *pX, t);
	for (int k = 1; k <= size; k++) {
		doublereal x = pX->operator()(k);
		doublereal xp = pXP->operator()(k);
		pXPm1->PutCoef(k, xp);
		pXm1->PutCoef(k, x-dt*xp);
	}

	// output iniziale
	std::cout << ti << " " << 0. << " ";
	::ff->out(p_data, std::cout, *pX, *pXP) << std::endl;

	flip(&pX, &pXP, &pXm1, &pXPm1, &pXm2, &pXPm2);

	while (t < tf) {
		t += dt;
		// predict
		for (int k = size; k > 0; k--) {
			doublereal xm1 = pXm1->operator()(k);
			doublereal xPm1 = pXPm1->operator()(k);
			doublereal xm2 = pXm2->operator()(k);
			doublereal xPm2 = pXPm2->operator()(k);
			doublereal x = -4.*xm1+5.*xm2+dt*(4.*xPm1+2.*xPm2);
			pX->PutCoef(k, x);
			R.PutCoef(k, a1*xm1+a2*xm2+dt*(b1*xPm1+b2*xPm2));
		}

		// test
		int j = 0;
		doublereal test;
		doublereal coef = dt*b0;
		do {
			::ff->func(p_data, *pXP, *pX, t);
			for (int k = 1; k <= size; k++) {
				doublereal x = pX->operator()(k);
				doublereal xP = pXP->operator()(k);
				Res.PutCoef(k, R(k)-x+coef*xP);
			}

			test = Res.Norm();
			if (test < tol) {
				break;
			}
			if (++j > maxiter) {
				std::cerr << "current iteration " << j
					<< " exceeds max iteration number "
					<< maxiter << std::endl;
					exit(EXIT_FAILURE);
			}

			// correct
			sm->MatrReset();
			J.Reset();
			::ff->grad(p_data, J, *pX, t);
			for (int k = 1; k <= size; k++) {
				for (int l = 1; l <= size; l++) {
					Jac.PutCoef(k, l, -coef*J(k, l));
				}
				Jac.IncCoef(k, k, 1.);
			}
			sm->Solve();

			// update
			for (int k = size; k > 0; k--) {
				doublereal dx = Sol(k);
				pX->IncCoef(k, dx);
			}
		} while (true);

		// output
		std::cout << t << " " << test << " ";
			::ff->out(p_data, std::cout, *pX, *pXP) << std::endl;

		flip(&pX, &pXP, &pXm1, &pXPm1, &pXm2, &pXPm2);
	}

	::ff->destroy(&p_data);

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
	// prepara i dati
	void* p_data = 0;
	::ff->read(&p_data, user_defined);

	// prepara le strutture dati per il calcolo
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

	FullMatrixHandler Jz(size);
	FullMatrixHandler J0(size);
	MyVectorHandler Xz(size);
	MyVectorHandler XPz(size);

	// TODO: abstract from a specific solution manager?
	LinSol ls;
	SolutionManager *sm = ls.GetSolutionManager(size);

	MatrixHandler& Jac = *sm->pMatHdl();
	VectorHandler& Res = *sm->pResHdl();
	VectorHandler& Sol = *sm->pSolHdl();

	// paramteri di integrazione
	doublereal ti = d->ti;
	doublereal dt = d->dt;
	doublereal tf = d->tf;

	doublereal tol = d->tol;
	int maxiter = d->maxiter;

	// coefficienti del metodo
	doublereal rho = d->rho;
	doublereal z = -rho/(1.+rho);
	doublereal w1 = (2.+3.*z)/(6.*(1.+z));
	doublereal wz = -1./(6.*z*(1.+z));
	doublereal w0 = (1.+3.*z)/(6.*z);
	doublereal m0 = 1.-z*z*(3.+2.*z);
	doublereal m1 = z*z*(3.+2.*z);
	doublereal n0 = z*(1.+z)*(1.+z);
	doublereal n1 = z*z*(1.+z);

	doublereal t = ti;

	// inizializza la soluzione
	::ff->init(p_data, *pX);
	::ff->func(p_data, *pXP, *pX, t);
	for (int k = 1; k <= size; k++) {
		doublereal x = pX->operator()(k);
		doublereal xp = pXP->operator()(k);
		pXPm1->PutCoef(k, xp);
		pXm1->PutCoef(k, x-dt*xp);
	}

	// output iniziale
	std::cout << ti << " " << 0. << " ";
	::ff->out(p_data, std::cout, *pX, *pXP) << std::endl;

	flip(&pX, &pXP, &pXm1, &pXPm1, &pXm2, &pXPm2);

	while (t < tf) {
		t += dt;
		// predict
		for (int k = 1; k <= size; k++) {
			doublereal xm1 = pXm1->operator()(k);
			doublereal xPm1 = pXPm1->operator()(k);
			doublereal xm2 = pXm2->operator()(k);
			doublereal xPm2 = pXPm2->operator()(k);
			doublereal x = -4.*xm1+5.*xm2+dt*(4.*xPm1+2.*xPm2);
			pX->PutCoef(k, x);
		}

		// test
		int j = 0;
		doublereal test;
		do {
			pXP->Reset();
			::ff->func(p_data, *pXP, *pX, t);
			for (int k = 1; k <= size; k++) {
				doublereal x = pX->operator()(k);
				doublereal xP = pXP->operator()(k);
				doublereal xm1 = pXm1->operator()(k);
				doublereal xPm1 = pXPm1->operator()(k);
				doublereal xz = m0*x+m1*xm1+dt*(n0*xP+n1*xPm1);
				Xz.PutCoef(k, xz);
			}
			XPz.Reset();
			::ff->func(p_data, XPz, Xz, t+z*dt);
			for (int k = 1; k <= size; k++) {
				doublereal d = dt*(
					w1*pXPm1->operator()(k)
					+ wz*XPz(k)
					+ w0*pXP->operator()(k)
				) - (
					pX->operator()(k)
					- pXm1->operator()(k)
				);
				Res.PutCoef(k, d);
			}

			test = Res.Norm();
			if (test < tol) {
				break;
			}
			if (++j > maxiter) {
				std::cerr << "current iteration " << j
					<< " exceeds max iteration number "
					<< maxiter << std::endl;
				exit(EXIT_FAILURE);
			}

			// correct
			sm->MatrReset();
			Jz.Reset();
			J0.Reset();
			::ff->grad(p_data, Jz, Xz, t+z*dt);
			::ff->grad(p_data, J0, *pX, t);
			for (int k = 1; k <= size; k++) {
				for (int l = 1; l <= size; l++) {
					doublereal d = 0.;
					for (int m = 1; m <= size; m++) {
						d += Jz(k, m)*J0(m, l);
					}
					d = -dt*(wz*(Jz(k, l)+dt*n0*d)+w0*J0(k, l));
					Jac.PutCoef(k, l, d);
				}
				Jac.IncCoef(k, k, 1.);
			}
			sm->Solve();

			// update
			for (int k = size; k > 0; k--) {
				doublereal dx = Sol(k);
				pX->IncCoef(k, dx);
			}
		} while (true);

		// output
		std::cout << t << " " << test << " ";
		::ff->out(p_data, std::cout, *pX, *pXP) << std::endl;

		flip(&pX, &pXP, &pXm1, &pXPm1, &pXm2, &pXPm2);
	}

	::ff->destroy(&p_data);

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
