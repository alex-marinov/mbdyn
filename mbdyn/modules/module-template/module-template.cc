/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 * 
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <mbconfig.h>

#include <loadable.h>

/*
 * user-defined struct
 */
struct module_template {
	int i; 
};

/* default funcs */
static void *
read(
		LoadableElem* pEl,
		DataManager* pDM,
		MBDynParser& HP,
		const DriveHandler* pDH
)
{
	DEBUGCOUTFNAME("read");
	
	/*
	 * allocation of user-defined struct
	 */
	module_template* p = NULL;
	SAFENEW(p, module_template);
	
	/*
	 * read data
	 */
	if (HP.IsKeyWord("help")) {
		std::cout << "Module; template" << std::endl;
	}
	
	return (void *)p;
}

static unsigned int
i_get_num_dof(const LoadableElem* pEl)
{
	DEBUGCOUTFNAME("i_get_num_dof");
	return 0;
}

static DofOrder::Order
set_dof(const LoadableElem*, unsigned int i)
{
	DEBUGCOUTFNAME("set_dof");
	return DofOrder::UNKNOWN;
}

static void
output(const LoadableElem* pEl, OutputHandler& OH)
{
	DEBUGCOUTFNAME("output");
}

static std::ostream&
restart(const LoadableElem* pEl, std::ostream& out)
{
	DEBUGCOUTFNAME("restart");
	return out << "not implemented yet;" << std::endl;
}

static void
work_space_dim(const LoadableElem* pEl, integer* piNumRows, integer* piNumCols)
{
	DEBUGCOUTFNAME("work_space_dim");
	*piNumRows = 0;
	*piNumCols = 0;
}

static VariableSubMatrixHandler& 
ass_jac(
		LoadableElem* pEl, 
		VariableSubMatrixHandler& WorkMat,
		doublereal dCoef, 
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr
)
{  
	DEBUGCOUTFNAME("ass_jac");
	integer iNumRows = 0;
	integer iNumCols = 0;
	pEl->WorkSpaceDim(&iNumRows, &iNumCols);
	
	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeInit(iNumRows, iNumCols, 0.);

#if 0
	module_template* p = (module_template *)pEl->pGetData();
#endif /* 0 */
	
	/*
	 * set sub-matrix indices and coefs
	 */
	
	return WorkMat;
}

static void
ass_mats(
		LoadableElem* pEl, 
		VariableSubMatrixHandler& WorkMatA,
		VariableSubMatrixHandler& WorkMatB,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr
)
{  
	DEBUGCOUTFNAME("ass_mats");
	integer iNumRows = 0;
	integer iNumCols = 0;
	pEl->WorkSpaceDim(&iNumRows, &iNumCols);
	
	FullSubMatrixHandler& WMA = WorkMatA.SetFull();
	WMA.ResizeInit(iNumRows, iNumCols, 0.);
	
	FullSubMatrixHandler& WMB = WorkMatB.SetFull();
	WMB.ResizeInit(iNumRows, iNumCols, 0.);

#if 0
	module_template* p = (module_template *)pEl->pGetData();
#endif /* 0 */
	
	/*
	 * set sub-matrix indices and coefs
	 */
}

static SubVectorHandler& 
ass_res(
		LoadableElem* pEl, 
		SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr
)
{
	DEBUGCOUTFNAME("ass_res");
	integer iNumRows = 0;
	integer iNumCols = 0;
	pEl->WorkSpaceDim(&iNumRows, &iNumCols);
	
	WorkVec.Resize(iNumRows);

#if 0
	module_template* p = (module_template *)pEl->pGetData(); 
#endif /* 0 */
	
	/*
	 * set sub-vector indices and coefs
	 */
	
	return WorkVec;
}

static void
before_predict(
		const LoadableElem* pEl, 
		VectorHandler& X,
		VectorHandler& XP,
		VectorHandler& XPrev,
		VectorHandler& XPPrev
)
{
	DEBUGCOUTFNAME("before_predict");
}

static void
after_predict(
		const LoadableElem* pEl, 
		VectorHandler& X,
		VectorHandler& XP
)
{
	DEBUGCOUTFNAME("after_predict");
}

static void
update(
		LoadableElem* pEl, 
		const VectorHandler& X,
		const VectorHandler& XP
)
{
	DEBUGCOUTFNAME("update");
}

static unsigned int
i_get_initial_num_dof(const LoadableElem* pEl)
{
	DEBUGCOUTFNAME("i_get_initial_num_dof");
	return 0;
}

static void
initial_work_space_dim(
		const LoadableElem* pEl, 
		integer* piNumRows, 
		integer* piNumCols
)
{
	DEBUGCOUTFNAME("initial_work_space_dim");
	*piNumRows = 0;
	*piNumCols = 0;   
}

static VariableSubMatrixHandler& 
initial_ass_jac(
		LoadableElem* pEl, 
		VariableSubMatrixHandler& WorkMat, 
		const VectorHandler& XCurr
)
{
	DEBUGCOUTFNAME("initial_ass_jac");
	integer iNumRows = 0;
	integer iNumCols = 0;
	pEl->InitialWorkSpaceDim(&iNumRows, &iNumCols);
	
	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeInit(iNumRows, iNumCols, 0.);
	
#if 0
	module_template* p = (module_template *)pEl->pGetData();
#endif /* 0 */
	
	/*
	 * set sub-matrix indices and coefs
	 */
	
	return WorkMat;
}

static SubVectorHandler& 
initial_ass_res(
		LoadableElem* pEl, 
		SubVectorHandler& WorkVec, 
		const VectorHandler& XCurr
)
{
	DEBUGCOUTFNAME("initial_ass_res");
	integer iNumRows = 0;
	integer iNumCols = 0;
	pEl->WorkSpaceDim(&iNumRows, &iNumCols);
	
	WorkVec.Resize(iNumRows);
	
#if 0
	module_template* p = (module_template *)pEl->pGetData(); 
#endif /* 0 */
	
	/*
	 * set sub-vector indices and coefs
	 */
	
	return WorkVec;
}

static void
set_value(const LoadableElem* pEl, VectorHandler& X, VectorHandler& XP)
{
	DEBUGCOUTFNAME("set_value");
}
   
static void
set_initial_value(const LoadableElem* pEl, VectorHandler& X)
{
	DEBUGCOUTFNAME("set_initial_value");
}

static unsigned int
i_get_num_priv_data(const LoadableElem* pEl)
{
	DEBUGCOUTFNAME("i_get_num_priv_data");
	return 0;
}

static doublereal
d_get_priv_data(const LoadableElem* pEl, unsigned int i)
{
	DEBUGCOUTFNAME("d_get_priv_data");
	ASSERT(pEl->iGetNumPrivData() > 0);
	if (i > pEl->iGetNumPrivData()) {
		std::cerr << "Module-template Elem: illegal private data index "
			<< i << std::endl;      
		THROW(ErrGeneric());
	}
   
	/*
	 * return i-th priv data
	 */
	return 0.;
}

static void
destroy(LoadableElem* pEl)
{
	DEBUGCOUTFNAME("destroy");

	module_template* p = (module_template *)pEl->pGetData();

	/*
	 * delete private data
	 */
	
	SAFEDELETE(p);
}

static int
i_get_num_connected_nodes(const LoadableElem* pEl)
{
	DEBUGCOUTFNAME("i_get_num_connected_nodes");
	return 0;
}

static void
get_connected_nodes(const LoadableElem* pEl, 
		int& NumNodes, Node::Type* NdTyp, unsigned int* NdLabels)
{
	DEBUGCOUTFNAME("get_connected_nodes");

#if 0
	module_template* p = (module_template *)pEl->pGetData();
#endif /* 0 */

	/*
	 * set args according to element connections
	 */
	NumNodes = i_get_num_connected_nodes(pEl);
}

static struct
LoadableCalls lc = {
	"template",
	read,
	i_get_num_dof,
	set_dof,
	output,
	restart,
	work_space_dim,
	ass_jac,
	ass_mats,
	ass_res,
	before_predict,
	after_predict,
	update,
	NULL, /* after_convergence */
	i_get_initial_num_dof,
	initial_work_space_dim,
	initial_ass_jac,
	initial_ass_res,
	set_value,
	set_initial_value,
	i_get_num_priv_data,
	d_get_priv_data,
	i_get_num_connected_nodes,
	get_connected_nodes,
	destroy
};

extern "C" {
void *calls = &lc;
}

