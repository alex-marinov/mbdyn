/* metodi di interpolazione meshless */

#ifndef MLS_H
#define MLS_H

#include <stroutput.h>
#include <spmapmh.h>
#include <ANN.h>

struct GeometryData;

class Weight
{
    public:
	Weight(void) {};
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
	virtual void SetP(double* PosFem, double** PosMb, int* id, unsigned int k, MyVectorHandler& p, FullMatrixHandler& P) = 0;
    
};
        
class Linear: public POrder
{
    public:
    
    Linear(void) {};
	    
        	/* 1 x y z */
 	void SetP(double* PosFem, double** PosMb, int* id, unsigned int k, MyVectorHandler& p, FullMatrixHandler& P);
    	
};

class Quadratic: public POrder
{
    public:
        Quadratic(void) {};
    
		    /* 1 x y z x^2 xy xz y^2 yz z^2 */
	void SetP(double* PosFem, double** PosMb, int* id, unsigned int k, MyVectorHandler& p, FullMatrixHandler& P);
    	
};

class InterpMethod
{

  public:

          //virtual unsigned int MaxNds(void) = 0;
	  
	  virtual void Interpolate(const GeometryData&, const GeometryData&, SpMapMatrixHandler*) = 0;
	  
	  virtual void Interpolate_Adj(const GeometryData&, const GeometryData&, SpMapMatrixHandler*,const std::vector<Vec3>&) = 0;

	  virtual ~InterpMethod(void){};
	 	
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

#endif // MLS_H
