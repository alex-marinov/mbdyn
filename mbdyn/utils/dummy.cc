#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <myassert.h>
#include <solman.h>
#include <harwrap.h>

struct private_data {
   int i;
};

int read(void** pp, const char* user_defined)
{
   *pp = (void*)new private_data;
   return 0;
}

int init(void* p, VectorHandler& X, VectorHandler& XP)
{
   // private_data* pd = (private_data*)p;
   X.Reset(0.);
   XP.Reset(0.);
   return 0;
}

int size(void* p)
{
   // private_data* pd = (private_data*)p;
   return 1;
}

int jac(void* p, MatrixHandler& JP, MatrixHandler& J,
	const VectorHandler& X, const VectorHandler& XP)
{
   // private_data* pd = (private_data*)p;
   J.fPutCoef(1, 1, 1.);
   return 0;
}

int res(void* p, VectorHandler& R, const doublereal&, 
	const VectorHandler& X, const VectorHandler& XP)
{
   // private_data* pd = (private_data*)p;
   R.fPutCoef(1, -1.*X.dGetCoef(1)-1.*XP.dGetCoef(1));
   return 0;
}

ostream& out(void* p, ostream& o, 
	     const VectorHandler& X, const VectorHandler& XP)
{
   // private_data* pd = (private_data*)p;
   return o;
}

int destroy(void** p)
{
   // private_data* pd = (private_data*)p;
   delete *p;
   *p = NULL;
   return 0;
}
