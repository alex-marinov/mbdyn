#ifndef JacSubMatrix_hh
#define JacSubMatrix_hh

#include <vector>

#include "ac/f2c.h"
#include "submat.h"

class ExpandableRowVector {
private:
	std::vector<doublereal> x;
	std::vector<ExpandableRowVector*> xm;
	std::vector<integer> idx;
	ExpandableRowVector & operator = (const ExpandableRowVector &); // not to be implemented
	ExpandableRowVector (const ExpandableRowVector &); // not to be implemented
public:
	ExpandableRowVector();
	ExpandableRowVector(const integer n);
	virtual ~ExpandableRowVector();
	void ReDim(const integer n);
	void Zero();
	void Reset();
	void Link(const integer i, ExpandableRowVector* xp);
	void Set(doublereal xx, integer i);
	void SetIdx(integer i, integer iidx);
	void Set(doublereal xx, integer i, integer iidx);
	void Add(doublereal xx, integer i);
	void Sub(doublereal xx, integer i);
	void Add(SubVectorHandler& WorkVec, const doublereal c = 1.);
	void Sub(SubVectorHandler& WorkVec, const doublereal c = 1.);
	void Add(FullSubMatrixHandler& WM, const integer eq, const doublereal c = 1.);
	void Sub(FullSubMatrixHandler& WM, const integer eq, const doublereal c = 1.);
};


#endif //JacSubMatrix_hh
