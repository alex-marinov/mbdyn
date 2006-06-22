/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2006
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h> 		/* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "loadable.h"
#include "module-aerodyn.h"

module_aerodyn_t *module_aerodyn = NULL;

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

	if (::module_aerodyn != NULL) {
		silent_cerr("Another module-aerodyn might be running; error"
			<< std::endl);
		throw ErrGeneric();
	}
	
	/*
	 * allocation of user-defined struct
	 */
	module_aerodyn_t* p = NULL;
	SAFENEW(p, module_aerodyn_t);
	
	/*
	 * read data
	 */
	if (HP.IsKeyWord("help")) {
		/* NOTE: add help message */
		silent_cout("Module: AeroDyn" << std::endl);
	}

	/* hub node */
	unsigned int uNode = (unsigned int)HP.GetInt();
       
	DEBUGCOUT("Hub: node " << uNode << std::endl);
       
	/* verifica di esistenza del nodo */
	p->pHub = pDM->pFindStructNode(uNode);
	if (p->pHub  == NULL) {
		silent_cerr("Aerodyn(" << pEl->GetLabel() << "): "
			"StructuralNode(" << uNode << ") not defined "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric();
	}

	/*
	 * Initialize AeroDyn package
	 *
	 * FIXME: does it return an error code?
	 */
	__FC_DECL__(adinputgate)();

	::module_aerodyn = p;

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
	module_aerodyn_t* p = (module_aerodyn_t *)pEl->pGetData(); 
#endif /* 0 */
	
	/*
	 * set sub-vector indices and coefs
	 */
	F_LOGICAL	FirstLoop;
	F_INTEGER	JElem;
	F_REAL		DFN, DFT, PMA;

	__FC_DECL__(aerofrcintrface)(&FirstLoop, &JElem, &DFN, &DFT, &PMA);
	
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

static void 
after_convergence(const LoadableElem* /* pEl */ ,
		const VectorHandler& /* X */ ,
		const VectorHandler& /* XP */ )
{
	DEBUGCOUTFNAME("after_convergence");
}

static unsigned int
i_get_initial_num_dof(const LoadableElem* pEl)
{
	DEBUGCOUTFNAME("i_get_initial_num_dof");
	return 0;
}

static void
set_value(const LoadableElem* pEl, DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
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

static void
destroy(LoadableElem* pEl)
{
	DEBUGCOUTFNAME("destroy");

	module_aerodyn_t* p = (module_aerodyn_t *)pEl->pGetData();

	/*
	 * delete private data
	 */
	::module_aerodyn = NULL;
	
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
	module_aerodyn_t* p = (module_aerodyn_t *)pEl->pGetData();
#endif /* 0 */

	/*
	 * set args according to element connections
	 */
	NumNodes = i_get_num_connected_nodes(pEl);
}

static struct
LoadableCalls lc = {
	LOADABLE_VERSION_SET(1, 2, 0),

	"AeroDyn",
	"1.1",
	"Dipartimento di Ingegneria Aerospaziale, Politecnico di Milano",
	"wind tower aerodynamics; uses AeroDyn package 12.3\n"
		"\tas distributed by the National Renewable Energy Laboratory,\n"
		"\thttp://www.nrel.gov",

	read,
	i_get_num_dof,
	set_dof,
	output,
	restart,
	work_space_dim,
	NULL /* ass_jac */ ,
	NULL /* ass_mats */ ,
	ass_res,
	before_predict,
	after_predict,
	update,
	after_convergence,
	i_get_initial_num_dof,
	NULL /* initial_work_space_dim */ ,
	NULL /* initial_ass_jac */ ,
	NULL /* initial_ass_res */ ,
	set_value,
	set_initial_value,
	i_get_num_priv_data,
	NULL /* i_get_priv_data_idx */ ,
	NULL /* d_get_priv_data */ ,
	i_get_num_connected_nodes,
	get_connected_nodes,
	destroy,
	NULL,
	NULL,
	NULL,
	NULL
};

extern "C" {
void *calls = &lc;
}

