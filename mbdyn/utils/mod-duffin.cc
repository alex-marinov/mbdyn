#include <myassert.h>
#include <solman.h>
#include <harwrap.h>
#include <fstream.h>
#include <math.h>

struct private_data {
   doublereal m;
   doublereal c;
   doublereal k;
   doublereal alpha;
   doublereal f;
   doublereal omega;
   doublereal x[2];
};

int read(void** pp, const char* user_defined)
{
   *pp = (void*)new private_data;
   private_data* pd = (private_data*)*pp;
   
   if (user_defined != NULL) {
      // cerr << "opening file \"" << user_defined << "\"" << endl;
      ifstream in(user_defined);
      if (!in) {
	 cerr << "unable to open file \"" << user_defined << "\"" << endl;
	 exit(EXIT_FAILURE);
      }
      in >> pd->m >> pd->c >> pd->k >> pd->alpha >> pd->f >> pd->omega
	>> pd->x[0] >> pd->x[1];
   } else {
      pd->m = 1.;
      pd->c = 1.e-2;
      pd->k = 1.;
      pd->alpha = 0.;
      pd->f = 1.;
      pd->omega = 1.;
      pd->x[0] = 0.;
      pd->x[1] = 0.;
   }
   
   // cerr << "m=" << pd->m << ", c=" << pd->c << ", k=" << pd->k << endl
   //   << "alpha=" << pd->alpha << ", f=" << pd->f << ", omega=" << pd->omega << endl
   //   << "x={" << pd->x[0] << "," << pd->x[1] << "}" << endl;
   
   return 0;
}

int size(void* p)
{
   // private_data* pd = (private_data*)p;
   return 2;
}

int init(void* p, VectorHandler& X)
{
   private_data* pd = (private_data*)p;
   X.Reset(0.);
   for (int i = 1; i <= size(p); i++) {      
      X.fPutCoef(i, pd->x[i-1]); /* posiz. iniziale */
   }
   return 0;
}

int jac(void* p, MatrixHandler& J, const VectorHandler& X, const doublereal& t)
{
   private_data* pd = (private_data*)p;
   doublereal x = X.dGetCoef(1);
   J.fPutCoef(1, 2, 1.);
   J.fPutCoef(2, 1, -(pd->k+3.*x*x*pd->alpha)/pd->m);
   J.fPutCoef(2, 2, -pd->c/pd->m);
   return 0;
}

int res(void* p, VectorHandler& R, const VectorHandler& X, const doublereal& t)
{
   private_data* pd = (private_data*)p;
   doublereal x = X.dGetCoef(1);
   doublereal v = X.dGetCoef(2);
   R.fPutCoef(1, v);
   R.fPutCoef(2, (pd->f*sin(pd->omega*t)-((pd->k+pd->alpha*x*x)*x+pd->c*v))/pd->m);
   return 0;
}

ostream& out(void* p, ostream& o, 
	     const VectorHandler& X, const VectorHandler& XP)
{
   private_data* pd = (private_data*)p;
   doublereal x = X.dGetCoef(1);
   doublereal v = X.dGetCoef(2);
   
   doublereal E = .5*(pd->m*v*v+(pd->k+.5*pd->alpha*x*x)*x*x);
		      
   return o << X.dGetCoef(1) << " " << X.dGetCoef(2)
     << " " << XP.dGetCoef(1) << " " << XP.dGetCoef(2) << " " << E;
}

int destroy(void** p)
{
   // private_data* pd = (private_data*)p;
   delete *p;
   *p = NULL;
   return 0;
}
