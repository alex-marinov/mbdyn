#include <loadable.h>

struct module_template {
   // user-defined struct
   int i; 
};

/* funzioni di default */
void* read(LoadableElem* pEl,
	   DataManager* pDM,
	   MBDynParser& HP,
	   const DriveHandler* pDH)
{
   DEBUGCOUTFNAME("read");
   
   // allocation of user-defined struct
   module_template* p = NULL;
   SAFENEW(p, module_template, EMmm);

   // read data
   
   return (void *)p;
}

unsigned int i_get_num_dof(const LoadableElem* pEl)
{
   DEBUGCOUTFNAME("i_get_num_dof");
   return 0;
}

DofOrder::Order set_dof(const LoadableElem*, unsigned int i)
{
   DEBUGCOUTFNAME("set_dof");
   return DofOrder::UNKNOWN;
}

void output(const LoadableElem* pEl, OutputHandler& OH)
{
   DEBUGCOUTFNAME("output");
}

ostream& restart(const LoadableElem* pEl, ostream& out)
{
   DEBUGCOUTFNAME("restart");
   return out << "not implemented yet;" << endl;
}

void work_space_dim(const LoadableElem* pEl, 
		    integer* piNumRows, 
		    integer* piNumCols)
{
   DEBUGCOUTFNAME("work_space_dim");
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
   DEBUGCOUTFNAME("ass_jac");
   integer iNumRows = 0;
   integer iNumCols = 0;
   pEl->WorkSpaceDim(&iNumRows, &iNumCols);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.ResizeInit(iNumRows, iNumCols, 0.);
   
   module_template* p = (module_template *)pEl->pGetData();
   
   // set sub-matrix indices and coefs
   
   return WorkMat;
}

void
ass_eig(LoadableElem* pEl, 
	VariableSubMatrixHandler& WorkMatA,
	VariableSubMatrixHandler& WorkMatB,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{  
   DEBUGCOUTFNAME("ass_eig");
   integer iNumRows = 0;
   integer iNumCols = 0;
   pEl->WorkSpaceDim(&iNumRows, &iNumCols);
   
   FullSubMatrixHandler& WMA = WorkMatA.SetFull();
   WMA.ResizeInit(iNumRows, iNumCols, 0.);
   
   FullSubMatrixHandler& WMB = WorkMatB.SetFull();
   WMB.ResizeInit(iNumRows, iNumCols, 0.);
   
   module_template* p = (module_template *)pEl->pGetData();
   
   // set sub-matrix indices and coefs
}

SubVectorHandler& 
ass_res(LoadableElem* pEl, 
	SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
   DEBUGCOUTFNAME("ass_res");
   integer iNumRows = 0;
   integer iNumCols = 0;
   pEl->WorkSpaceDim(&iNumRows, &iNumCols);
   
   WorkVec.Resize(iNumRows);
   
   module_template* p = (module_template *)pEl->pGetData(); 

   // set sub-vector indices and coefs
   
   return WorkVec;
}

void before_predict(const LoadableElem* pEl, 
		    VectorHandler& X,
		    VectorHandler& XP,
		    VectorHandler& XPrev,
		    VectorHandler& XPPrev)
{
   DEBUGCOUTFNAME("before_predict");
}

void after_predict(const LoadableElem* pEl, 
		   VectorHandler& X,
		   VectorHandler& XP)
{
   DEBUGCOUTFNAME("after_predict");
}

void update(LoadableElem* pEl, 
	    const VectorHandler& X,
	    const VectorHandler& XP)
{
   DEBUGCOUTFNAME("update");
}

unsigned int i_get_initial_num_dof(const LoadableElem* pEl)
{
   DEBUGCOUTFNAME("i_get_initial_num_dof");
   return 0;
}

void initial_work_space_dim(const LoadableElem* pEl, 
			    integer* piNumRows, 
			    integer* piNumCols)
{
   DEBUGCOUTFNAME("initial_work_space_dim");
   *piNumRows = 0;
   *piNumCols = 0;   
}

VariableSubMatrixHandler& 
initial_ass_jac(LoadableElem* pEl, 
		VariableSubMatrixHandler& WorkMat, 
		const VectorHandler& XCurr)
{
   DEBUGCOUTFNAME("initial_ass_jac");
   integer iNumRows = 0;
   integer iNumCols = 0;
   pEl->InitialWorkSpaceDim(&iNumRows, &iNumCols);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.ResizeInit(iNumRows, iNumCols, 0.);
   
   module_template* p = (module_template *)pEl->pGetData();
   
   // set sub-matrix indices and coefs
   
   return WorkMat;
}

SubVectorHandler& 
initial_ass_res(LoadableElem* pEl, 
		SubVectorHandler& WorkVec, 
		const VectorHandler& XCurr)
{
   DEBUGCOUTFNAME("initial_ass_res");
   integer iNumRows = 0;
   integer iNumCols = 0;
   pEl->WorkSpaceDim(&iNumRows, &iNumCols);
   
   WorkVec.Resize(iNumRows);
   
   module_template* p = (module_template *)pEl->pGetData(); 

   // set sub-vector indices and coefs
   
   return WorkVec;
}

void set_value(const LoadableElem* pEl, VectorHandler& X, VectorHandler& XP)
{
   DEBUGCOUTFNAME("set_value");
}
   
void set_initial_value(const LoadableElem* pEl, VectorHandler& X)
{
   DEBUGCOUTFNAME("set_initial_value");
}

unsigned int i_get_num_priv_data(const LoadableElem* pEl)
{
   DEBUGCOUTFNAME("i_get_num_priv_data");
   return 0;
}

doublereal d_get_priv_data(const LoadableElem* pEl, unsigned int i)
{
   DEBUGCOUTFNAME("d_get_priv_data");
   ASSERT(pEl->iGetNumPrivData() > 0);
   if (i > pEl->iGetNumPrivData()) {
      cerr << "Module-template Elem: illegal private data index " << i << endl;      
      THROW(ErrGeneric());
   }
   
   // return i-th priv data
   return 0.;
}

void destroy(LoadableElem* pEl)
{
   DEBUGCOUTFNAME("destroy");
   module_template* p = (module_template *)pEl->pGetData();
   SAFEDELETE(p, EMmm);
}
