#include <mbconfig.h>
#include <loadable.h>


/* funzioni di default */
static void *
read(LoadableElem*, DataManager*, MBDynParser&, const DriveHandler*)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);
   return NULL;
}

static unsigned int
i_get_num_dof(const LoadableElem* pEl)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);
   return 0;
}

static DofOrder::Order
set_dof(const LoadableElem*, unsigned int i)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);
   return DofOrder::UNKNOWN;
}

static void 
output(const LoadableElem* pEl, OutputHandler& OH)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl); 
}

static ostream&
restart(const LoadableElem* pEl, ostream& out)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);
   return out << "not implemented yet;" << endl;
}

static void
work_space_dim(const LoadableElem* pEl, integer* piNumRows, integer* piNumCols)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);
   *piNumRows = 0;
   *piNumCols = 0;
}

static VariableSubMatrixHandler& 
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

static void
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

static SubVectorHandler& 
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

static void
before_predict(const LoadableElem* pEl, 
		VectorHandler& X,
		VectorHandler& XP,
		VectorHandler& XPrev,
		VectorHandler& XPPrev)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);  
}

static void
after_predict(const LoadableElem* pEl, VectorHandler& X, VectorHandler& XP)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);  
}

static void
update(LoadableElem* pEl, const VectorHandler& X, const VectorHandler& XP)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);  
}

static unsigned int
i_get_initial_num_dof(const LoadableElem* pEl)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);
   return 0;
}

static void
initial_work_space_dim(const LoadableElem* pEl, 
		integer* piNumRows, 
		integer* piNumCols)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);
   *piNumRows = 0;
   *piNumCols = 0;   
}

static VariableSubMatrixHandler& 
initial_ass_jac(LoadableElem* pEl, 
		VariableSubMatrixHandler& WorkMat, 
		const VectorHandler& XCurr)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);
   WorkMat.SetNullMatrix();
   return WorkMat;
}

static SubVectorHandler& 
initial_ass_res(LoadableElem* pEl, 
		SubVectorHandler& WorkVec, 
		const VectorHandler& XCurr)
{  
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);
   WorkVec.Resize(0);
   return WorkVec;
}

static void
set_value(const LoadableElem* pEl, VectorHandler& X, VectorHandler& XP)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);
}
   
static void
set_initial_value(const LoadableElem* pEl, VectorHandler& X)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);
}

static unsigned int
i_get_num_priv_data(const LoadableElem* pEl)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);
   return 0;
}

static doublereal
d_get_priv_data(const LoadableElem* pEl, unsigned int i)
{
   cerr << "You shouldn't be here!" << endl;
   THROW(ErrGeneric());
}

static void
destroy(LoadableElem* pEl)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << endl);
}

static struct LoadableCalls lc = {
	read,
	i_get_num_dof,
	set_dof,
	output,
	restart,
	work_space_dim,
	ass_jac,
	ass_eig,
	ass_res,
	before_predict,
	after_predict,
	update,
	i_get_initial_num_dof,
	initial_work_space_dim,
	initial_ass_jac,
	initial_ass_res,
	set_value,
	set_initial_value,
	i_get_num_priv_data,
	d_get_priv_data,
	destroy
};

extern "C" void *calls = &lc;

