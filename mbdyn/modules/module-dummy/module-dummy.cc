#include <loadable.h>


/* funzioni di default */
void* read(LoadableElem*,
	   DataManager*,
	   MBDynParser&,
	   const DriveHandler*)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);
   return NULL;
}

unsigned int i_get_num_dof(const LoadableElem* pEl)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);
   return 0;
}

DofOrder::Order set_dof(const LoadableElem*, unsigned int i)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);
   return DofOrder::UNKNOWN;
}

void output(const LoadableElem* pEl, OutputHandler& OH)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl); 
}

ostream& restart(const LoadableElem* pEl, ostream& out)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);
   return out << "not implemented yet;" << endl;
}

void work_space_dim(const LoadableElem* pEl, 
		    integer* piNumRows, 
		    integer* piNumCols)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);
   *piNumRows = 0;
   *piNumCols = 0;
}

VariableSubMatrixHandler& 
ass_jac(LoadableElem* pEl, 
	VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);
   WorkMat.SetNullMatrix();
   return WorkMat;
}

void
ass_eig(LoadableElem* pEl, 
	VariableSubMatrixHandler& WorkMatA,
	VariableSubMatrixHandler& WorkMatB,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);
   WorkMatA.SetNullMatrix();  
   WorkMatB.SetNullMatrix();  
}

SubVectorHandler& 
ass_res(LoadableElem* pEl, 
	SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);
   WorkVec.Resize(0);
   return WorkVec;
}

void before_predict(const LoadableElem* pEl, 
		    VectorHandler& X,
		    VectorHandler& XP,
		    VectorHandler& XPrev,
		    VectorHandler& XPPrev)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);  
}

void after_predict(const LoadableElem* pEl, 
		   VectorHandler& X,
		   VectorHandler& XP)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);  
}

void update(LoadableElem* pEl, 
	    const VectorHandler& X,
	    const VectorHandler& XP)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);  
}

unsigned int i_get_initial_num_dof(const LoadableElem* pEl)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);
   return 0;
}

void initial_work_space_dim(const LoadableElem* pEl, 
			    integer* piNumRows, 
			    integer* piNumCols)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);
   *piNumRows = 0;
   *piNumCols = 0;   
}

VariableSubMatrixHandler& 
initial_ass_jac(LoadableElem* pEl, 
		VariableSubMatrixHandler& WorkMat, 
		const VectorHandler& XCurr)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);
   WorkMat.SetNullMatrix();
   return WorkMat;
}

SubVectorHandler& 
initial_ass_res(LoadableElem* pEl, 
		SubVectorHandler& WorkVec, 
		const VectorHandler& XCurr)
{  
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);
   WorkVec.Resize(0);
   return WorkVec;
}

void set_value(const LoadableElem* pEl, VectorHandler& X, VectorHandler& XP)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);
}
   
void set_initial_value(const LoadableElem* pEl, VectorHandler& X)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);
}

unsigned int i_get_num_priv_data(const LoadableElem* pEl)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);
   return 0;
}

doublereal d_get_priv_data(const LoadableElem* pEl, unsigned int i)
{
   cerr << "You shouldn't be here!" << endl;
   THROW(ErrGeneric());
}

void destroy(LoadableElem* pEl)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);
}
