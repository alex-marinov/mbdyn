/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#if defined(HAVE_LTDL_H) || defined(HAVE_DLFCN_H)

/* FIXME: temporary hack ... */
//#undef HAVE_LTDL_H

#include <string.h>
#ifdef HAVE_LTDL_H
#include <ltdl.h>
#elif defined(HAVE_DLFCN_H)
#include <dlfcn.h>
#else /* !HAVE_LTDL_H && !HAVE_DLFCN_H */
#error "no dynamic linking headers, sorry"
#endif /* !HAVE_LTDL_H && HAVE_DLFCN_H */

#include <ac/getopt.h>

#include <myassert.h>
#include <solman.h>
#include <fullmh.h>
#include <y12wrap.h>
#include <harwrap.h>
#include <mschwrap.h>
#include <umfpackwrap.h>

struct integration_data {
   	doublereal ti;
   	doublereal tf;
   	doublereal dt;
   	doublereal tol;
   	int maxiter;
   	doublereal rho;
};

int method_multistep(const char*, integration_data*, void*, const char*);
int method_hope(const char*, integration_data*, void*, const char*);
int method_cubic(const char*, integration_data*, void*, const char*);
int method_cn(const char*, integration_data*, void*, const char*);

void* get_method_data(int, const char*);

int 
main(int argn, char *const argv[])
{
   	enum {
      		METHOD_UNKNOWN,
		METHOD_MULTISTEP,
		METHOD_HOPE,
		METHOD_CUBIC,
		METHOD_CRANK_NICHOLSON
   	};
	
   	struct s_method {
      		const char* s;
      		int m;
   	} s_m[] = {
		{ "ms",                METHOD_MULTISTEP         },
		{ "hope",              METHOD_HOPE              },
		{ "cubic",             METHOD_CUBIC             },
		{ "crank-nicholson",   METHOD_CRANK_NICHOLSON   },
	
		{ NULL,                METHOD_UNKNOWN           }
   	};
   	int curr_method = METHOD_UNKNOWN;   
   	char* module = "intg-default.so";
   	char* user_defined = NULL;
   	void* method_data = NULL;
   	integration_data d = { 0., 1., 1.e-3, 1.e-6, 10 };

   	/* opzioni */
   	const char optstring[] = "i:m:t:T:n:r:u:h";
   	const struct option longopts[] = { 
		{ "integrator",     required_argument, NULL, int('i') },
		{ "method-data",    required_argument, NULL, int('m') },
		{ "timing",         required_argument, NULL, int('t') },
		{ "tolerance",      required_argument, NULL, int('T') },
		{ "iterations",     required_argument, NULL, int('n') },
		{ "rho",            required_argument, NULL, int('r') },
		{ "user-defined",   required_argument, NULL, int('u') },
		{ "help",           no_argument,       NULL, int('h') },
	
		{ NULL,             no_argument,       NULL, int(0)   }  
   	};

   	while (true) {
		int curropt;
		
		curropt = getopt_long(argn, argv, optstring, longopts, NULL);
		
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
	    			<< " -i,--integrator   : integration method"
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
	    			<< " -h,--help         : print this message"
				" and exit" << std::endl;
	  		exit(EXIT_SUCCESS);

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
   
   	if (argv[optind] != NULL) {
      		module = argv[optind];
      		/* std::cerr << "using module " << module << std::endl; */
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
    	case METHOD_CRANK_NICHOLSON:
      		rc = method_cn(module, &d, method_data, user_defined);
      		break;
    	default:
      		std::cerr << "you must select an integration method" << std::endl;
      		exit(EXIT_FAILURE);
   	}
   
#ifdef HAVE_LTDL_H
	if (lt_dlexit()) {
		std::cerr << "lt_dlexit() failed" << std::endl;
		exit(EXIT_FAILURE);
	}
#endif /* HAVE_LTDL_H */
	
   	return 0;
}

void* 
get_method_data(int curr_method, const char* optarg)
{
   	switch (curr_method) {
    	default:
      		std::cerr << "not implemented yet" << std::endl;
      		exit(EXIT_FAILURE);
   	}
   
   	return NULL;
}

#include <dae-intg.h>
static struct funcs *ff = NULL;

int
open_module(const char* module) 
{
   	const char* err = NULL;

#ifdef HAVE_LTDL_H
	lt_dlhandle handle;

	if (lt_dlinit()) {
		std::cerr << "lt_dlinit() failed" << std::endl;
      		exit(EXIT_FAILURE);
   	}

   	if ((handle = lt_dlopen(module)) == NULL) {
      		err = lt_dlerror();
      		std::cerr << "lt_dlopen(\"" << module << "\") returned \"" << err
			<< "\"" << std::endl;
      		exit(EXIT_FAILURE);
   	}

	struct funcs **pf = (struct funcs **)lt_dlsym(handle, "ff");
   	if (pf == NULL) {
      		err = lt_dlerror();
      		std::cerr << "lt_dlsym(\"ff\") returned \"" << err << "\""
			<< std::endl;
      		exit(EXIT_FAILURE);
   	}

#elif defined(HAVE_DLFCN_H)
   	void* handle = NULL;

   	if ((handle = dlopen(module, RTLD_NOW)) == NULL) {
      		err = dlerror();
      		std::cerr << "dlopen(\"" << module << "\") returned \"" << err
			<< "\"" << std::endl;
      		exit(EXIT_FAILURE);
   	}
   
	struct funcs **pf = (struct funcs **)dlsym(handle, "ff");
   	if (pf == NULL) {
      		err = dlerror();
      		std::cerr << "dlsym(\"ff\") returned \"" << err << "\""
			<< std::endl;
      		exit(EXIT_FAILURE);
   	}
#endif /* !HAVE_LTDL_H && HAVE_DLFCN_H */

	::ff = *pf;
	if (::ff == NULL) {
		std::cerr << "invalid \"ff\" symbol in \"" << module << "\""
			<< std::endl;
      		exit(EXIT_FAILURE);
   	}

   	return 0;
}

void
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

int
method_multistep(const char* module, integration_data* d, 
		 void* method_data, const char* user_defined)
{
   	open_module(module);
 
   	// prepara i dati
   	void* p_data = NULL;
   	(*::ff->read)(&p_data, user_defined);
   
   	// prepara le strutture dati per il calcolo
   	int size = (*::ff->size)(p_data);
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

#if defined(USE_UMFPACK)
	UmfpackSparseLUSolutionManager sm(size);
#elif defined(USE_Y12)
	Y12SparseLUSolutionManager sm(size, size*(size+1)+1, 1.);
#elif defined(USE_HARWELL)
   	HarwellSparseLUSolutionManager sm(size, size*(size+1)+1, 1.);
#elif defined(USE_MESCHACH)
	MeschachSparseLUSolutionManager sm(size, size*(size+1)+1, 1.);
#endif
   
   	MatrixHandler& Jac = *sm.pMatHdl();
   	VectorHandler& Res = *sm.pResHdl();   
   	VectorHandler& Sol = *sm.pSolHdl();   

   	// paramteri di integrazione
   	doublereal ti = d->ti;
   	doublereal dt = d->dt;
   	doublereal tf = d->tf;
   
   	doublereal tol = d->tol;
   	int maxiter = d->maxiter;
   
   	// coefficienti del metodo
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
   
   	// inizializza la soluzione
   	(*::ff->init)(p_data, *pX, *pXP);
   	for (int k = 1; k <= size; k++) {
      		doublereal x = pX->dGetCoef(k);
      		doublereal xp = pXP->dGetCoef(k);
      		pXPm1->fPutCoef(k, xp);
      		pXm1->fPutCoef(k, x-dt*xp);
   	}
   
   	// output iniziale
   	std::cout << ti << " " << 0. << " ";
   	(*::ff->out)(p_data, std::cout, *pX, *pXP) << std::endl;
   
   	flip(&pX, &pXP, &pXm1, &pXPm1, &pXm2, &pXPm2);
   
   	while (t < tf) {
      		t += dt;
		
      		// predict
      		for (int ir = size; ir > 0; ir--) {
	 		doublereal xm1 = pXm1->dGetCoef(ir);
	 		doublereal xPm1 = pXPm1->dGetCoef(ir);
	 		doublereal xm2 = pXm2->dGetCoef(ir);
	 		doublereal xPm2 = pXPm2->dGetCoef(ir);
	 		doublereal xP = cm0p(1.)*xm1+cm1p(1.)*xm2+dt*(cn0p(1.)*xPm1+cn1p(1.)*xPm2);
			doublereal x = a1*xm1 + a2*xm2 + dt*(b0*xP + b1*xPm1 + b2*xPm2);
	 		pX->fPutCoef(ir, x);
	 		pXP->fPutCoef(ir, xP);
      		}
      
      		// test
      		int j = 0;
      		doublereal test;
      		doublereal coef = dt*b0;
      		do {
	 		(*::ff->func)(p_data, Res, *pX, *pXP, t);

	 		test = Res.Norm();
	 		if (test < tol) {
	    			break;
	 		}
	 		if (++j > maxiter) {
	    			std::cerr << "current iteration " << j 
	      				<< " exceedes max iteration number "
					<< maxiter << std::endl;
	    			exit(EXIT_FAILURE);
	 		}
	 
	 		// correct
	 		sm.MatrInit();
	 		J.Init();
	 		(*::ff->grad)(p_data, J, JP, *pX, *pXP, t);
	 		for (int ir = 1; ir <= size; ir++) {
	    			for (int ic = 1; ic <= size; ic++) {
	       				Jac.fPutCoef(ir, ic, JP.dGetCoef(ir, ic) + coef*J.dGetCoef(ir, ic));
	    			}
	 		}
	 		sm.Solve();
			
	 		// update
	 		for (int ir = size; ir > 0; ir--) {
	    			doublereal dxP = Sol.dGetCoef(ir);	
	    			pXP->fIncCoef(ir, dxP);
	    			pX->fIncCoef(ir, coef*dxP);
	 		}
      		} while (true);
      
      		// output
      		std::cout << t << " " << test << " ";
      		(*::ff->out)(p_data, std::cout, *pX, *pXP) << std::endl;
      
      		flip(&pX, &pXP, &pXm1, &pXPm1, &pXm2, &pXPm2);            
   	}
   
   	(*::ff->destroy)(&p_data);
   	delete[] pd;
   	delete[] ppd;

   	return 0;
}

int
method_hope(const char* module, integration_data* d, 
	    void* method_data, const char* user_defined)
{
   	open_module(module);
   
   	std::cerr << __FUNCTION__ << "not implemented yet!" << std::endl;
   	exit(EXIT_FAILURE);
   
   	return 0;
}

int
method_cubic(const char* module, integration_data* d, 
	     void* method_data, const char* user_defined)
{
   	open_module(module);

   	// prepara i dati
   	void* p_data = NULL;
   	(*::ff->read)(&p_data, user_defined);
   
   	// prepara le strutture dati per il calcolo
   	int size = (*::ff->size)(p_data);
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
   
#if defined(USE_UMFPACK)
	UmfpackSparseLUSolutionManager sm(2*size);
#elif defined(USE_Y12)
	Y12SparseLUSolutionManager sm(2*size, 2*size*(2*size+1)+1, 1.);
#elif defined(USE_HARWELL)
   	HarwellSparseLUSolutionManager sm(2*size, 2*size*(2*size+1)+1, 1.);
#elif defined(USE_MESCHACH)
	MeschachSparseLUSolutionManager sm(2*size, 2*size*(2*size+1)+1, 1.);
#endif
   
   	MatrixHandler& Jac = *sm.pMatHdl();
   	VectorHandler& Res = *sm.pResHdl();   
   	VectorHandler& Sol = *sm.pSolHdl();   

   	MyVectorHandler Resz(size, Res.pdGetVec() + size);

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
   	doublereal m0 = cm0(z);
   	doublereal m1 = cm1(z);
   	doublereal n0 = cn0(z);
   	doublereal n1 = cn1(z);
   	doublereal m0p = cm0p(z);
   	doublereal m1p = cm1p(z);
   	doublereal n0p = cn0p(z);
   	doublereal n1p = cn1p(z);

	doublereal jzz = (1.+3.*rho)/(6.*rho*(1.+rho))*dt;
	doublereal jz0 = -1./(6.*rho*(1.+rho)*(1.+rho))*dt;
	doublereal j0z = (1.+rho)*(1.+rho)/(6.*rho)*dt;
	doublereal j00 = (2.*rho-1.)/(6.*rho)*dt;

   	doublereal t = ti;
   
   	// inizializza la soluzione
   	(*::ff->init)(p_data, *pX, *pXP);
   	(*::ff->func)(p_data, Res, *pX, *pXP, t);
   	for (int k = 1; k <= size; k++) {
      		doublereal x = pX->dGetCoef(k);
      		doublereal xp = pXP->dGetCoef(k);
      		pXPm1->fPutCoef(k, xp);
      		pXm1->fPutCoef(k, x-dt*xp);
   	}
   
   	// output iniziale
   	std::cout << ti << " " << 0. << " ";
   	(*::ff->out)(p_data, std::cout, *pX, *pXP) << std::endl;
   
   	flip(&pX, &pXP, &pXm1, &pXPm1, &pXm2, &pXPm2);
   
   	while (t < tf) {
      		t += dt;
		
      		// predict cubico
      		for (int k = 1; k <= size; k++) {
	 		doublereal xm1 = pXm1->dGetCoef(k);
	 		doublereal xPm1 = pXPm1->dGetCoef(k);
	 		doublereal xm2 = pXm2->dGetCoef(k);
	 		doublereal xPm2 = pXPm2->dGetCoef(k);
	 		doublereal xP = cm0p(1.)*xm1+cm1p(1.)*xm2
				+dt*(cn0p(1.)*xPm1+cn1p(1.)*xPm2);
	 		doublereal xPz = cm0p(1.+z)*xm1+cm1p(1.+z)*xm2
				+dt*(cn0p(1.+z)*xPm1+cn1p(1.+z)*xPm2);
	 		doublereal x = xm1 + dt*(w1*xPm1+wz*xPz+w0*xP);
	 		doublereal xz = m0*x+m1*xm1+dt*(n0*xP+n1*xPm1);

	 		pX->fPutCoef(k, x);
	 		pXP->fPutCoef(k, xP);
			Xz.fPutCoef(k, xz);
			XPz.fPutCoef(k, xPz);
      		}
      
      		// test
      		int j = 0;
      		doublereal test;      
      		do {
	 		Res.Reset();
	 		(*::ff->func)(p_data, Res, *pX, *pXP, t);
	 		Resz.Reset();
	 		(*::ff->func)(p_data, Resz, Xz, XPz, t+z*dt);

	 		test = Res.Norm() + Resz.Norm();
	 		if (test < tol) {
	    			break;
	 		}
	 		if (++j > maxiter) {
	    			std::cerr << "current iteration " << j 
	      				<< " exceedes max iteration number "
					<< maxiter << std::endl;
	    			exit(EXIT_FAILURE);
	 		}
	 
	 		// correct
	 		sm.MatrInit();
	 		Jz.Init();
	 		JPz.Init();
	 		J0.Init();
	 		JP0.Init();
	 		(*::ff->grad)(p_data, Jz, JPz, Xz, XPz, t+z*dt);
	 		(*::ff->grad)(p_data, J0, JP0, *pX, *pXP, t);

			for (int ir = 1; ir <= size; ir++) {
				for (int ic = 1; ic <= size; ic++) {
					Jac.fIncCoef(ir, ic, j00*J0.dGetCoef(ir, ic) + JP0.dGetCoef(ir, ic));
					Jac.fIncCoef(ir, size+ic, j0z*J0.dGetCoef(ir, ic));
					Jac.fIncCoef(size+ir, ic, jz0*Jz.dGetCoef(ir, ic));
					Jac.fIncCoef(size+ir, size+ic, jzz*Jz.dGetCoef(ir, ic) + JPz.dGetCoef(ir, ic));
				}
			}

	 		sm.Solve();
	 
	 		// update
	 		for (int ir = size; ir > 0; ir--) {
	    			doublereal dxP0 = Sol.dGetCoef(ir);
	    			pXP->fIncCoef(ir, dxP0);
	    			pX->fIncCoef(ir, dt*w0*dxP0);
	    			doublereal dxPz = Sol.dGetCoef(size+ir);
				XPz.fIncCoef(ir, dxPz);
				Xz.fIncCoef(ir, dt*(m0*(wz*dxPz+w0*dxP0)+n0*dxP0));
	 		}
     	 	} while (true);
      
      		// output
      		std::cout << t << " " << test << " ";
      		(*::ff->out)(p_data, std::cout, *pX, *pXP) << std::endl;
      
      		flip(&pX, &pXP, &pXm1, &pXPm1, &pXm2, &pXPm2);            
   	}
   
   	(*::ff->destroy)(&p_data);
   	delete[] pd;
   	delete[] ppd;

   	return 0;
}

int
method_cn(const char* module, integration_data* d,
	  void* method_data, const char* user_defined)
{
   	open_module(module);
   
   	std::cerr << __FUNCTION__ << "not implemented yet!" << std::endl;
   	exit(EXIT_FAILURE);

   	return 0;
}

#else /* defined(HAVE_LTDL_H) || defined(HAVE_DLFCN_H) */

int
main(void)
{
	std::cerr << "Need dynamic load capabilities" << atd::endl;
   	exit(EXIT_FAILURE);
}

#endif /* defined(HAVE_LTDL_H) || defined(HAVE_DLFCN_H) */

