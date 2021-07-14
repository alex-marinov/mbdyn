/* metodi di interpolazione meshless */

#ifndef INTERP_HH
#define INTERP_HH

#include <map>
#include <vector>

#include "fullmh.h"
#include "spmapmh.h"
#include "vh.h"
#if defined(HAVE_ANN_H)
#include <ANN.h>
#elif defined(HAVE_ANN_ANN_H)
#include <ANN/ANN.h>
#else
// should never happen
#error "need ANN.h"
#endif

#include "node3.h"

class Weight;
class POrder;

class InterpMethod 
{

  public:
   
	//virtual unsigned int MaxNds(void) = 0;
	
	virtual void Interpolate(std::vector<Node>&, std::vector<Node>&, SpMapMatrixHandler*) = 0;
	
	virtual ~InterpMethod(void){};
};

/* Mean Least Square Polinomial */
class MLSP: public InterpMethod
{
     	
    unsigned int N;
    unsigned int MaxN;
    unsigned int nadj;
    FullMatrixHandler pP;
    SpMapMatrixHandler pW;
    MyVectorHandler pp;
    
    Weight* pWg;
    POrder* pOr;
    
    ANNkd_tree* pTree; 
	
    public:
    
    	MLSP();
        
	MLSP(unsigned int NodesN, int WgType, bool LinQuad, int Nadj);
	
	void ShowMe();
	
	void Interpolate(std::vector<Node>&, std::vector<Node>&, SpMapMatrixHandler*);
				
};

#endif //INTERP_HH
