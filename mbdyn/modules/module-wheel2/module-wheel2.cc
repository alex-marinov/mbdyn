#include <mbconfig.h>

#include <loadable.h>

/*
 *	3 corpi: 
 *	-	Axle
 *	-	Wheel
 *	-	Ground
 *
 *	Note: 
 *	-	Axle e Wheel sono collegati da un giunto che consente
 *		solo rotazione relativa attorno ad un asse (axle)
 *	-	Ground supporta un piano definito dalla normale e da
 *		un punto (potra' essere reso deformabile e "non piano"
 *		in futuro)
 */

struct module_wheel {
	StructNode *pAxle;
	StructNode *pWheel;
	StructNode *pGround;

	Vec3 GroundPosition;
	Vec3 GroundDirection;
};

/* funzioni di default */
static void*
read(LoadableElem* pEl,
	   DataManager* pDM,
	   MBDynParser& HP,
	   const DriveHandler* pDH)
{
	DEBUGCOUTFNAME("read");
	
	/* allocation of user-defined struct */
	module_wheel* p = NULL;
	SAFENEW(p, module_wheel, EMmm);

	p->pAxle = (StructNode *)pDM->ReadNode(HP, Node::STRUCTURAL);
	
	return (void *)p;
}

unsigned int
i_get_num_dof(const LoadableElem* pEl)
{
	DEBUGCOUTFNAME("i_get_num_dof");
	return 0;
}

static void
output(const LoadableElem* pEl, OutputHandler& OH)
{
	DEBUGCOUTFNAME("output");
	
	if (pEl->fToBeOutput()) {
		
	}
}

ostream&
restart(const LoadableElem* pEl, ostream& out)
{
	DEBUGCOUTFNAME("restart");
	return out << "not implemented yet;" << endl;
}

static void
work_space_dim(const LoadableElem* pEl, 
		    integer* piNumRows, 
		    integer* piNumCols)
{
	DEBUGCOUTFNAME("work_space_dim");
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
	DEBUGCOUTFNAME("ass_jac");   
	WorkMat.SetNullMatrix();
	
	return WorkMat;
}

static SubVectorHandler& 
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
	WorkVec.Reset(0.);
	
	module_wheel* p = (module_wheel *)pEl->pGetData();

	return WorkVec;
}

void
before_predict(const LoadableElem* pEl, 
		VectorHandler& X,
		VectorHandler& XP,
		VectorHandler& XPrev,
		VectorHandler& XPPrev)
{
   DEBUGCOUTFNAME("before_predict");
}

void
after_predict(const LoadableElem* pEl, 
		VectorHandler& X,
		VectorHandler& XP)
{
	DEBUGCOUTFNAME("after_predict");
}

void
update(LoadableElem* pEl, 
		const VectorHandler& X,
		const VectorHandler& XP)
{
	DEBUGCOUTFNAME("update");
}

unsigned int
i_get_initial_num_dof(const LoadableElem* pEl)
{
	DEBUGCOUTFNAME("i_get_initial_num_dof");
	return 0;
}

void
initial_work_space_dim(const LoadableElem* pEl, 
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
	
	module_wheel* p = (module_wheel *)pEl->pGetData();
   
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
	
	module_wheel* p = (module_wheel *)pEl->pGetData(); 
	
	// set sub-vector indices and coefs
   
	return WorkVec;
}

void
set_value(const LoadableElem* pEl, VectorHandler& X, VectorHandler& XP)
{
	DEBUGCOUTFNAME("set_value");
}
   
void
set_initial_value(const LoadableElem* pEl, VectorHandler& X)
{
	DEBUGCOUTFNAME("set_initial_value");
}

unsigned int
i_get_num_priv_data(const LoadableElem* pEl)
{
	DEBUGCOUTFNAME("i_get_num_priv_data");
	return 0;
}

doublereal d_get_priv_data(const LoadableElem* pEl, unsigned int i)
{
	DEBUGCOUTFNAME("d_get_priv_data");
	ASSERT(pEl->iGetNumPrivData() > 0);
	if (i > pEl->iGetNumPrivData()) {
		cerr << "Module-template Elem: illegal private data index "
			<< i << endl;      
		THROW(ErrGeneric());
	}
	
	// return i-th priv data

	return 0.;
}

static void
destroy(LoadableElem* pEl)
{
	DEBUGCOUTFNAME("destroy");
	module_wheel* p = (module_wheel *)pEl->pGetData();
	
	SAFEDELETE(p, EMmm);
}

static struct
LoadableCalls lc = {
	read, /* */
	NULL /* i_get_num_dof */ ,
	NULL /* set_dof */ ,
	output, /* */
	NULL /* restart */ ,
	work_space_dim, /* */
	ass_jac, /* */
	NULL /* ass_eig */ ,
	ass_res, /* */
	NULL /* before_predict */ ,
	NULL /* after_predict */ ,
	NULL /* update */ ,
	NULL /* i_get_initial_num_dof */ ,
	NULL /* initial_work_space_dim */ ,
	NULL /* initial_ass_jac */ ,
	NULL /* initial_ass_res */ ,
	NULL /* set_value */ ,
	NULL /* set_initial_value */ ,
	NULL /* i_get_num_priv_data */ ,
	NULL /* d_get_priv_data */ ,
	destroy
};

extern "C" {
void *calls = &lc;
}


