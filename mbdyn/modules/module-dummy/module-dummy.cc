/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2011
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
#include "loadable.h"


/* funzioni di default */
static void *
read(LoadableElem*, DataManager*, MBDynParser&)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << std::endl);
   return NULL;
}

static unsigned int
i_get_num_dof(const LoadableElem* pEl)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << std::endl);
   return 0;
}

static DofOrder::Order
set_dof(const LoadableElem*, unsigned int i)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << std::endl);
   return DofOrder::UNKNOWN;
}

static void 
output(const LoadableElem* pEl, OutputHandler& OH)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << std::endl); 
}

static std::ostream&
restart(const LoadableElem* pEl, std::ostream& out)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << std::endl);
   return out << "not implemented yet;" << std::endl;
}

static void
work_space_dim(const LoadableElem* pEl, integer* piNumRows, integer* piNumCols)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << std::endl);
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
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << std::endl);
   WorkMat.SetNullMatrix();
   return WorkMat;
}

static void
ass_mats(LoadableElem* pEl, 
	VariableSubMatrixHandler& WorkMatA,
	VariableSubMatrixHandler& WorkMatB,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << std::endl);
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
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << std::endl);
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
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << std::endl);  
}

static void
after_predict(const LoadableElem* pEl, VectorHandler& X, VectorHandler& XP)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << std::endl);  
}

static void
update(LoadableElem* pEl, const VectorHandler& X, const VectorHandler& XP)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << std::endl);  
}

static unsigned int
i_get_initial_num_dof(const LoadableElem* pEl)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << std::endl);
   return 0;
}

static void
initial_work_space_dim(const LoadableElem* pEl, 
		integer* piNumRows, 
		integer* piNumCols)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << std::endl);
   *piNumRows = 0;
   *piNumCols = 0;   
}

static VariableSubMatrixHandler& 
initial_ass_jac(LoadableElem* pEl, 
		VariableSubMatrixHandler& WorkMat, 
		const VectorHandler& XCurr)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << std::endl);
   WorkMat.SetNullMatrix();
   return WorkMat;
}

static SubVectorHandler& 
initial_ass_res(LoadableElem* pEl, 
		SubVectorHandler& WorkVec, 
		const VectorHandler& XCurr)
{  
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << std::endl);
   WorkVec.Resize(0);
   return WorkVec;
}

static void
set_value(const LoadableElem* pEl, DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << std::endl);
}
   
static void
set_initial_value(const LoadableElem* pEl, VectorHandler& X)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << std::endl);
}

static unsigned int
i_get_num_priv_data(const LoadableElem* pEl)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << std::endl);
   return 0;
}

static unsigned int
i_get_priv_data_idx(const LoadableElem* pEl, const char *s)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << std::endl);
   return 0;
}

static doublereal
d_get_priv_data(const LoadableElem* pEl, unsigned int i)
{
   std::cerr << "You shouldn't be here!" << std::endl;
   throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

static void
destroy(LoadableElem* pEl)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << std::endl);
}

static int
i_get_num_connected_nodes(const LoadableElem* pEl)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << std::endl);
   return 0;
}

static void
get_connected_nodes(const LoadableElem* pEl, 
		int& NumNodes, Node::Type* NdTyp, unsigned int* NdLabels)
{
   DEBUGCOUT("Dummy Elem: " << __PRETTY_FUNCTION__ << std::endl);
   NumNodes = i_get_num_connected_nodes(pEl);
}

static struct LoadableCalls lc = {
	LOADABLE_VERSION_SET(1, 5, 0),

	"dummy",
	"1.2",
	"Dipartimento di Ingegneria Aerospaziale, Politecnico di Milano",
	"Dummy module --- does nothing, used to test the infrastructure",
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
	NULL /* after_convergence */ ,
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
	destroy,
	NULL,
	NULL,
	NULL,
	NULL
};

extern "C" void *calls = &lc;

