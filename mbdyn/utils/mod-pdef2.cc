#include <myassert.h>
#include <solman.h>
#include <fstream.h>
#include <math.h>

struct private_data {
   doublereal m;
   doublereal c;
   doublereal k;
   doublereal l;
   doublereal g;
   doublereal x[4];
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
      in >> pd->m >> pd->c >> pd->k >> pd->l >> pd->g
	>> pd->x[0] >> pd->x[1] >> pd->x[2] >> pd->x[3];
   } else {
      pd->m = 1.;
      pd->c = 0.;
      pd->k = 1.;     
      pd->l = 1.;
      pd->g = 9.81;     
      pd->x[0] = 0.;
      pd->x[1] = 0.;
      pd->x[2] = 1.;
      pd->x[3] = 0.;
   }
   
   doublereal theta = pd->x[0];
   doublereal uu = pd->x[1];
   doublereal thetap = pd->x[2];
   doublereal uup = pd->x[3];
   
   pd->x[0] = (pd->l+uu)*sin(theta);
   pd->x[1] = -(pd->l+uu)*cos(theta);
   pd->x[2] = thetap*(pd->l+uu)*cos(theta)+uup*sin(theta);
   pd->x[3] = thetap*(pd->l+uu)*sin(theta)-uup*cos(theta);   
   
   cerr << "m=" << pd->m << ", c=" << pd->c << ", k=" << pd->k << endl
     << "l=" << pd->l << ", g=" << pd->g << endl
     << "x={" << pd->x[0] << "," << pd->x[1] << "," 
     << pd->x[2] << "," << pd->x[3] << "}" << endl;
   
   return 0;
}

int size(void* p)
{
   // private_data* pd = (private_data*)p;
   return 4;
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

int grad(void* p, MatrixHandler& J, const VectorHandler& X, const doublereal& t)
{
   private_data* pd = (private_data*)p;
   
   doublereal x = X.dGetCoef(1);
   doublereal y = X.dGetCoef(2);
   doublereal u = X.dGetCoef(3);
   doublereal v = X.dGetCoef(4);
   
   doublereal l = sqrt(x*x+y*y);
   doublereal m = pd->m;
   // doublereal c = pd->c;
   doublereal k = pd->k;
   doublereal l0 = pd->l; 
   doublereal g = pd->g;
   
   J.fPutCoef(1, 3, 1.);
   J.fPutCoef(2, 4, 1.);
   J.fPutCoef(3, 1, -k/m*(1.-l0/l*(1.-(x*x)/(l*l))));
   J.fPutCoef(3, 2, -k/m*x*y*l0/(l*l*l));   
   J.fPutCoef(4, 1, -k/m*x*y*l0/(l*l*l));
   J.fPutCoef(4, 2, -k/m*(1.-l0/l*(1.-(y*y)/(l*l))));

   return 0;
}

int func(void* p, VectorHandler& R, const VectorHandler& X, const doublereal& t)
{
   private_data* pd = (private_data*)p;
   
   doublereal x = X.dGetCoef(1);
   doublereal y = X.dGetCoef(2);
   doublereal u = X.dGetCoef(3);
   doublereal v = X.dGetCoef(4);
   
   doublereal m = pd->m;
   // doublereal c = pd->c;
   doublereal k = pd->k;
   doublereal l0 = pd->l;
   doublereal l = sqrt(x*x+y*y);
   doublereal g = pd->g;

   R.fPutCoef(1, u);
   R.fPutCoef(2, v);
   R.fPutCoef(3, -k/m*x*(1.-l0/l));
   R.fPutCoef(4, -k/m*y*(1.-l0/l)-g);

   return 0;
}

ostream& out(void* p, ostream& o, 
	     const VectorHandler& X, const VectorHandler& XP)
{
   private_data* pd = (private_data*)p;
  
   doublereal x = X.dGetCoef(1);
   doublereal y = X.dGetCoef(2);
   doublereal u = X.dGetCoef(3);
   doublereal v = X.dGetCoef(4);   
   doublereal m = pd->m;
   doublereal k = pd->k;
   doublereal l0 = pd->l;
   doublereal l = sqrt(x*x+y*y);
   doublereal theta = atan2(x, -y);
   doublereal g = pd->g;
   
   doublereal E = .5*m*(u*u+v*v)+m*g*y+.5*k*(l-l0)*(l-l0);
  
   
   return o << theta << " " << (l-l0)
     << " " << X.dGetCoef(3) << " " << X.dGetCoef(4)
     << " " << XP.dGetCoef(1) << " " << XP.dGetCoef(2)
     << " " << XP.dGetCoef(3) << " " << XP.dGetCoef(4)
     << " " << x << " " << y << " " << E;
}

int destroy(void** p)
{
   // private_data* pd = (private_data*)p;
   delete *p;
   *p = NULL;
   return 0;
}
