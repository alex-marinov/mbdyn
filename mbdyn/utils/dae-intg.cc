#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <myassert.h>
#include <solman.h>
#include <fullmh.h>
#include <y12wrap.h>
#include <harwrap.h>
#include <mschwrap.h>

#include <getopt.h>
#include <string.h>
#include <dlfcn.h>

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

   	while (1) {
		int curropt;
		
		curropt = getopt_long(argn, argv, optstring, longopts, NULL);
		
		if (curropt == EOF) {
			break;
		}
		
      		switch(curropt) {
       		case int('?'):
	  		/* 
			 * cerr << "unknown option -" << char(optopt) << endl;
			 */
	  		break;
			
       		case int(':'):
	  		cerr << "missing parameter for option -" 
				<< optopt << endl;
	  		break;
			
       		case int('h'):
	  		cerr << "usage: int -[imtTnruh] [module]" << endl
	    			<< endl
	    			<< " -i,--integrator   : integration method"
				<< endl
	    			<< " -m,--method-data  : method-dependent data"
				<< endl
	    			<< " -t,--timing       : ti:dt:tf" << endl
	    			<< " -T,--tolerance" << endl
	    			<< " -n,--maxiter" << endl
	    			<< " -r,--rho          : asymptotic radius"
				<< endl
	    			<< " -u,--user         : user-defined data"
				<< endl
	    			<< endl
	    			<< " -h,--help         : print this message"
				" and exit" << endl;
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
	     			cerr << "unknown integrator " 
					<< optarg << endl;
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
	  		char* p = strrchr(s, ':');
	  		if (p == NULL) {
	     			cerr << "syntax: ti:dt:tf" << endl;
	     			exit(EXIT_FAILURE);
	  		}
	  		d.tf = atof(p+1);
	  		p[0] = '\0';
	  		p = strrchr(s, ':');
	  		if (p == NULL) {
	     			cerr << "syntax: ti:dt:tf" << endl;
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
	  		user_defined = new char[strlen(optarg)+1];
	  		strcpy(user_defined, optarg);
	  		break;
			
       		default:
	  		/* cerr << "unknown option" << endl; */
	  		break;
      		}
   	}
   
   	if (argv[optind] != NULL) {
      		module = argv[optind];
      		/* cerr << "using module " << module << endl; */
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
      		cerr << "you must select an integration method" << endl;
      		exit(EXIT_FAILURE);
   	}
   
   	if (user_defined != NULL) {
      		delete[] user_defined;
      		user_defined = NULL;
   	}
   	return 0;
}

void* 
get_method_data(int curr_method, const char* optarg)
{
   	switch (curr_method) {
    	default:
      		cerr << "not implemented yet" << endl;
      		exit(EXIT_FAILURE);
   	}
   
   	return NULL;
}

#include <dae-intg.h>
static funcs *ff;

int
open_module(const char* module) 
{
   	void* handle = NULL;
   	const char* err = NULL;
   	if ((handle = dlopen(module, RTLD_NOW)) == NULL) {
      		err = dlerror();
      		cerr << "dlopen(\"" << module << "\") returned \"" << err
			<< "\"" << endl;
      		exit(EXIT_FAILURE);
   	}
   
   	if ((::ff = (funcs *)dlsym(handle, "ff")) == NULL) {
      		err = dlerror();
      		cerr << "dlsym(\"ff\") returned \"" << err << "\""
			<< endl;
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
	FullMatrixHandler Jp(pd+size*size, ppd+size, size*size, size, size);
   	MyVectorHandler R(size);

#if defined(USE_Y12)
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
   	(*::ff->init)(p_data, *pX);
   	(*::ff->func)(p_data, *pXP, *pX, t);
   	for (int k = 1; k <= size; k++) {
      		doublereal x = pX->dGetCoef(k);
      		doublereal xp = pXP->dGetCoef(k);
      		pXPm1->fPutCoef(k, xp);
      		pXm1->fPutCoef(k, x-dt*xp);
   	}
   
   	// output iniziale
   	cout << ti << " " << 0. << " ";
   	(*::ff->out)(p_data, cout, *pX, *pXP) << endl;
   
   	flip(&pX, &pXP, &pXm1, &pXPm1, &pXm2, &pXPm2);
   
   	while (t < tf) {
      		t += dt;
		
      		// predict
      		for (int k = size; k > 0; k--) {
	 		doublereal xm1 = pXm1->dGetCoef(k);
	 		doublereal xPm1 = pXPm1->dGetCoef(k);
	 		doublereal xm2 = pXm2->dGetCoef(k);
	 		doublereal xPm2 = pXPm2->dGetCoef(k);
	 		doublereal x = -4.*xm1+5.*xm2+dt*(4.*xPm1+2.*xPm2);
	 		pX->fPutCoef(k, x);
	 		R.fPutCoef(k, a1*xm1+a2*xm2+dt*(b1*xPm1+b2*xPm2));
      		}
      
      		// test
      		int j = 0;
      		doublereal test;
      		doublereal coef = dt*b0;
      		do {
	 		(*::ff->func)(p_data, *pXP, *pX, t);
	 		for (int k = 1; k <= size; k++) {
	    			doublereal x = pX->dGetCoef(k);
	    			doublereal xP = pXP->dGetCoef(k);
	    			Res.fPutCoef(k, R.dGetCoef(k)-x+coef*xP);
	 		}

	 		test = Res.Norm();
	 		if (test < tol) {
	    			break;
	 		}
	 		if (++j > maxiter) {
	    			cerr << "current iteration " << j 
	      				<< " exceedes max iteration number "
					<< maxiter << endl;
	    			exit(EXIT_FAILURE);
	 		}
	 
	 		// correct
	 		sm.MatrInit();
	 		J.Init();
	 		(*::ff->grad)(p_data, J, Jp, *pX, t);
	 		for (int k = 1; k <= size; k++) {
	    			for (int l = 1; l <= size; l++) {
	       				Jac.fPutCoef(k, l,
						     -coef*J.dGetCoef(k, l));
	    			}
	    			Jac.fIncCoef(k, k, 1.);
	 		}
	 		sm.Solve();
			
	 		// update
	 		for (int k = size; k > 0; k--) {
	    			doublereal dx = Sol.dGetCoef(k);	
	    			pX->fIncCoef(k, dx);
	 		}
      		} while (1);
      
      		// output
      		cout << t << " " << test << " ";
      		(*::ff->out)(p_data, cout, *pX, *pXP) << endl;
      
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
   
   	cerr << __FUNCTION__ << "not implemented yet!" << endl;
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
	FullMatrixHandler Jpz(pd+2*size*size,
			      ppd+2*size, size*size, size, size);
	FullMatrixHandler Jp0(pd+3*size*size,
			      ppd+3*size, size*size, size, size);
   	MyVectorHandler Xz(size);
   	MyVectorHandler XPz(size);
   
#if defined(USE_Y12)
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
   	(*::ff->init)(p_data, *pX);
   	(*::ff->func)(p_data, *pXP, *pX, t);
   	for (int k = 1; k <= size; k++) {
      		doublereal x = pX->dGetCoef(k);
      		doublereal xp = pXP->dGetCoef(k);
      		pXPm1->fPutCoef(k, xp);
      		pXm1->fPutCoef(k, x-dt*xp);
   	}
   
   	// output iniziale
   	cout << ti << " " << 0. << " ";
   	(*::ff->out)(p_data, cout, *pX, *pXP) << endl;
   
   	flip(&pX, &pXP, &pXm1, &pXPm1, &pXm2, &pXPm2);
   
   	while (t < tf) {
      		t += dt;
		
      		// predict
      		for (int k = 1; k <= size; k++) {
	 		doublereal xm1 = pXm1->dGetCoef(k);
	 		doublereal xPm1 = pXPm1->dGetCoef(k);
	 		doublereal xm2 = pXm2->dGetCoef(k);
	 		doublereal xPm2 = pXPm2->dGetCoef(k);
	 		doublereal x = -4.*xm1+5.*xm2+dt*(4.*xPm1+2.*xPm2);
	 		pX->fPutCoef(k, x);
      		}
      
      		// test
      		int j = 0;
      		doublereal test;      
      		do {
	 		pXP->Reset();
	 		(*::ff->func)(p_data, *pXP, *pX, t);
	 		for (int k = 1; k <= size; k++) {
	    			doublereal x = pX->dGetCoef(k);
	    			doublereal xP = pXP->dGetCoef(k);
	    			doublereal xm1 = pXm1->dGetCoef(k);
	    			doublereal xPm1 = pXPm1->dGetCoef(k);
	    			doublereal xz = m0*x+m1*xm1+dt*(n0*xP+n1*xPm1);
	    			Xz.fPutCoef(k, xz);
	 		}
	 		XPz.Reset();
	 		(*::ff->func)(p_data, XPz, Xz, t+z*dt);
	 		for (int k = 1; k <= size; k++) {
	    			doublereal d = dt*(
			       		w1*pXPm1->dGetCoef(k)
			       		+wz*XPz.dGetCoef(k)
			       		+w0*pXP->dGetCoef(k)
			       		)-(
				  	pX->dGetCoef(k)
				  	-pXm1->dGetCoef(k)
				  	);
	    			Res.fPutCoef(k, d);
	 		}

	 		test = Res.Norm();
	 		if (test < tol) {
	    			break;
	 		}
	 		if (++j > maxiter) {
	    			cerr << "current iteration " << j 
	      				<< " exceedes max iteration number "
					<< maxiter << endl;
	    			exit(EXIT_FAILURE);
	 		}
	 
	 		// correct
	 		sm.MatrInit();
	 		Jz.Init();
	 		J0.Init();
	 		(*::ff->grad)(p_data, Jz, Jpz, Xz, t+z*dt);
	 		(*::ff->grad)(p_data, J0, Jp0, *pX, t);	 
	 		for (int k = 1; k <= size; k++) {
	    			for (int l = 1; l <= size; l++) {
	       				doublereal d = 0.;
	       				for (int m = 1; m <= size; m++) {
		  				d += Jz.dGetCoef(k, m)
							*J0.dGetCoef(m, l);
	       				}
	       				d = -dt*(wz*(Jz.dGetCoef(k, l)+dt*n0*d)
						+w0*J0.dGetCoef(k, l));
	       				Jac.fPutCoef(k, l, d);
	    			}
	    			Jac.fIncCoef(k, k, 1.);
	 		}
	 		sm.Solve();
	 
	 		// update
	 		for (int k = size; k > 0; k--) {
	    			doublereal dx = Sol.dGetCoef(k);
	    			pX->fIncCoef(k, dx);
	 		}
     	 	} while (1);
      
      		// output
      		cout << t << " " << test << " ";
      		(*::ff->out)(p_data, cout, *pX, *pXP) << endl;
      
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
   
   	cerr << __FUNCTION__ << "not implemented yet!" << endl;
   	exit(EXIT_FAILURE);

   	return 0;
}

