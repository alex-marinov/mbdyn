#ifndef DAE_INTG_H
#define DAE_INTG_H

typedef int (*pread)(void**, const char*);
typedef int (*pinit)(void*, VectorHandler&);
typedef int (*psize)(void*);
typedef int (*pgrad)(void*, MatrixHandler&, MatrixHandler&, 
		     const VectorHandler&, const doublereal&);
typedef int (*pfunc)(void*, VectorHandler&,
		     const VectorHandler&, const doublereal&);
typedef ostream& (*pout)(void*, ostream&,
			 const VectorHandler&, const VectorHandler&);
typedef int (*pdestroy)(void**);

typedef struct _funcs {
	pread read;
	pinit init;
	psize size;
	pgrad grad;
	pfunc func;
	pout out;
	pdestroy destroy;
} funcs;

#endif /* DAE_INTG_H */

