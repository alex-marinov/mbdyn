/* metodi di interpolazione meshless */

#ifndef MLS_H
#define MLS_H

#ifdef USE_ANN

#if defined(HAVE_ANN_H)
#include <ANN.h>
#elif defined(HAVE_ANN_ANN_H)
#include <ANN/ANN.h>
#else
// should never happen
#error "need ANN.h"
#endif

#include "fullmh.h"
#include "spmapmh.h"
#include "geomdata.h"

class Weight
{
public:
	virtual ~Weight(void) {};
	virtual void SetWeight(MyVectorHandler& r_in, SparseMatrixHandler& w_out) = 0;
};

class RBF0: public Weight
{
public:
	inline void SetWeight(MyVectorHandler& r, SparseMatrixHandler& w_out);
};

class RBF2: public Weight
{
public:
	inline void SetWeight(MyVectorHandler& r, SparseMatrixHandler& w_out);
};

class RBF4: public Weight
{
public:
	inline void SetWeight(MyVectorHandler& r, SparseMatrixHandler& w_out);
};

class POrder
{
public:
	virtual ~POrder(void) {};
	virtual void
	SetP(double* PosFem, double** PosMb, int* id,
		unsigned int k, MyVectorHandler& p, FullMatrixHandler& P) = 0;
};
        
class Linear: public POrder
{
public:
	/* 1 x y z */
 	void
	SetP(double* PosFem, double** PosMb, int* id,
		unsigned int k, MyVectorHandler& p, FullMatrixHandler& P);
};

class Quadratic: public POrder
{
public:
	/* 1 x y z x^2 xy xz y^2 yz z^2 */
	void
	SetP(double* PosFem, double** PosMb, int* id,
		unsigned int k, MyVectorHandler& p, FullMatrixHandler& P);
};

class InterpMethod
{
public:
	virtual ~InterpMethod(void) {};

          // virtual unsigned int MaxNds(void) = 0;
	virtual void
	Interpolate(const GeometryData&, const GeometryData&, SpMapMatrixHandler*) = 0;
	  
	virtual void
	Interpolate_Adj(const GeometryData&, const GeometryData&,
		SpMapMatrixHandler*, const std::vector<Vec3>&) = 0;
};

/* Mean Least Square Polinomial */
class MLSP: public InterpMethod
{
	unsigned int N;
	FullMatrixHandler pP;
	SpMapMatrixHandler pW;
	MyVectorHandler pp;    
    
	Weight* pWg;
	POrder* pOr;
    
	ANNkd_tree* pTree;
    
public:
    	MLSP(unsigned int NodesN, int WgType, bool LinQuad);
	void Interpolate(const GeometryData&, const GeometryData&, SpMapMatrixHandler*);
	void Interpolate_Adj(const GeometryData&, const GeometryData&, SpMapMatrixHandler*,const std::vector<Vec3>&);
};

#endif // USE_ANN

#endif // MLS_H
