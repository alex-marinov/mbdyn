/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 * Paolo Mantegazza     <mantegazza@aero.polimi.it>
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <iostream>
#include <limits>
#include <cfloat>
#include <limits>

#include "dataman.h"
#include "loadable.h"

/*
 * Usage:
 *
 *	loadable: <label>, <module name>, help [ , <module data> ] ;
 */

struct module_marble {
	StructNode *pNode;

	doublereal dRadius;
	Vec3 e3R;

	Vec3 F;
};

static const Vec3 e1(1., 0., 0.);
static const Vec3 e2(0., 1., 0.);
static const Vec3 e3(0., 0., 1.);

/* funzioni di default */
static void*
marble_read(LoadableElem* pEl,
		DataManager* pDM,
		MBDynParser& HP)
{
	DEBUGCOUTFNAME("marble_read");
	
	/* allocation of user-defined struct */
	module_marble* p = NULL;
	SAFENEW(p, module_marble);

	/*
	 * help
	 */
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	marble							\n"
"Author: 	Pierangelo Masarati <masarati@aero.polimi.it>		\n"
"Organization:	Dipartimento di Ingegneria Aerospaziale			\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it				\n"
"									\n"
"	All rights reserved						\n"
"									\n"
"     -	Input:								\n"
"		<marble structural node label> ,			\n"
"		<radius>						\n"
"									\n"
"     -	Output:								\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	p->pNode = (StructNode *)pDM->ReadNode(HP, Node::STRUCTURAL);
	p->dRadius = HP.GetReal();
	p->e3R = e3*p->dRadius;

	return (void *)p;
}

static unsigned int
marble_i_get_num_dof(const LoadableElem* pEl)
{
	DEBUGCOUTFNAME("marble_i_get_num_dof");
	return 3;
}

static DofOrder::Order
marble_set_dof(const LoadableElem*, unsigned int i)
{
	DEBUGCOUTFNAME("marble_set_dof");
	ASSERT(i >= 0 && i < 2);
	return DofOrder::ALGEBRAIC; 
}

static void
marble_output(const LoadableElem* pEl, OutputHandler& OH)
{
	DEBUGCOUTFNAME("marble_output");
	
	if (pEl->bToBeOutput()) {
		module_marble* p = (module_marble *)pEl->pGetData();      

		std::ostream& out = OH.Loadable();

		out << std::setw(8) << pEl->GetLabel()	/* 1:	label */
			<< " " << p->F
			<< std::endl;
	}
}

#if 0
static std::ostream&
restart(const LoadableElem* pEl, std::ostream& out)
{
	DEBUGCOUTFNAME("restart");
	return out << "not implemented yet;" << std::endl;
}
#endif

static void
marble_work_space_dim(const LoadableElem* pEl, 
	integer* piNumRows, 
	integer* piNumCols)
{
	DEBUGCOUTFNAME("marble_work_space_dim");
	*piNumRows = 9;
	*piNumCols = 9;
}

static VariableSubMatrixHandler& 
marble_ass_jac(LoadableElem* pEl, 
	VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{  
	DEBUGCOUTFNAME("marble_ass_jac");   

	module_marble* p = (module_marble *)pEl->pGetData();

	integer iNumRows = 0;
	integer iNumCols = 0;
	pEl->WorkSpaceDim(&iNumRows, &iNumCols);
   
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	WM.Resize(iNumRows, iNumCols);

	integer iFirstPosIndex = p->pNode->iGetFirstPositionIndex();
	integer iFirstMomIndex = p->pNode->iGetFirstMomentumIndex();
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iFirstMomIndex + iCnt);
		WM.PutColIndex(iCnt, iFirstPosIndex + iCnt);
	}

	integer iFirstIndex = pEl->iGetFirstIndex();
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(6 + iCnt, iFirstIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iFirstIndex + iCnt);
	}

	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutCoef(iCnt, 6 + iCnt, 1.);

		WM.PutCoef(6 + iCnt, iCnt, 1.);
	}

	WM.PutCoef(4, 6 + 2, p->dRadius);
	WM.PutCoef(5, 6 + 1, -p->dRadius);

	WM.PutCoef(6 + 1, 5, -p->dRadius);
	WM.PutCoef(6 + 2, 4, p->dRadius);

	const Vec3& W = p->pNode->GetWCurr();

	WM.PutCoef(6 + 1, 4, p->dRadius*W(3)*dCoef);
	WM.PutCoef(6 + 2, 5, p->dRadius*W(3)*dCoef);

	WM.PutCoef(6 + 1, 6, -p->dRadius*W(1)*dCoef);
	WM.PutCoef(6 + 2, 6, -p->dRadius*W(2)*dCoef);

	return WorkMat;
}

static SubVectorHandler& 
marble_ass_res(LoadableElem* pEl, 
	SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("marble_ass_res");
	
	module_marble* p = (module_marble *)pEl->pGetData();

	integer iNumRows = 0;
	integer iNumCols = 0;
	pEl->WorkSpaceDim(&iNumRows, &iNumCols);
   
	WorkVec.ResizeReset(iNumRows);

	integer iFirstMomIndex = p->pNode->iGetFirstMomentumIndex();
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstMomIndex + iCnt);
	}

	integer iFirstIndex = pEl->iGetFirstIndex();
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(6 + iCnt, iFirstIndex + iCnt);
	}

	p->F = Vec3(XCurr, iFirstIndex + 1);

	const Vec3& X = p->pNode->GetXCurr();
	Vec3 V = p->pNode->GetVCurr() - p->pNode->GetWCurr().Cross(p->e3R);

	WorkVec.Sub(1, p->F);
	WorkVec.Add(4, p->e3R.Cross(p->F));

	WorkVec.PutCoef(7, -V(1));
	WorkVec.PutCoef(8, -V(2));
	WorkVec.PutCoef(9, (p->dRadius - X(3))/dCoef);

	return WorkVec;
}

static unsigned int
marble_i_get_initial_num_dof(const LoadableElem* pEl)
{
	DEBUGCOUTFNAME("marble_i_get_initial_num_dof");
	return 0;
}

static void
marble_destroy(LoadableElem* pEl)
{
	DEBUGCOUTFNAME("marble_destroy");
	module_marble* p = (module_marble *)pEl->pGetData();
	
	SAFEDELETE(p);
}

static int
marble_i_get_num_connected_nodes(const LoadableElem* /* pEl */ )
{
	DEBUGCOUTFNAME("marble_i_get_num_connected_nodes");
	
	return 1;
}

static void
marble_get_connected_nodes(const LoadableElem* pEl, 
		std::vector<const Node *>& connectedNodes)
{
	DEBUGCOUTFNAME("marble_get_connected_nodes");
	module_marble *p = (module_marble *)pEl->pGetData();
	
	connectedNodes.resize(1);

	connectedNodes[0] = p->pNode;
}

#ifdef MBDYN_MODULE
static
#endif /* MBDYN_MODULE */
struct
LoadableCalls module_marble_lc = {
	LOADABLE_VERSION_SET(1, 5, 0),

	"marble",
	"0.1.0",
	"Dipartimento di Ingegneria Aerospaziale, Politecnico di Milano",
	"marble model for icecube robot\n"
		"\tcontact Pierangelo Masarati <masarati@aero.polimi.it>",

	marble_read /* */ ,
	marble_i_get_num_dof /* */ ,
	marble_set_dof /* */ ,
	marble_output /* */ ,
	NULL /* marble_restart */ ,
	marble_work_space_dim /* */ ,
	marble_ass_jac /* */ ,
	NULL /* marble_ass_mats */ ,
	marble_ass_res /* */ ,
	NULL /* marble_before_predict */ ,
	NULL /* marble_after_predict */ ,
	NULL /* marble_update */ ,
	NULL /* marble_after_convergence */ ,
	marble_i_get_initial_num_dof,
	NULL /* marble_initial_work_space_dim */ ,
	NULL /* marble_initial_ass_jac */ ,
	NULL /* marble_initial_ass_res */ ,
	NULL /* marble_set_value */ ,
	NULL /* marble_set_initial_value */ ,
	NULL /* marble_i_get_num_priv_data */ ,
	NULL /* marble_i_get_priv_data_idx */ ,
	NULL /* marble_d_get_priv_data */ ,
	marble_i_get_num_connected_nodes,
	marble_get_connected_nodes,
	marble_destroy,
	NULL,
	NULL,
	NULL,
	NULL
};

extern "C" {

void *calls = &module_marble_lc;

#ifndef STATIC_MODULES
extern "C" int
module_init(const char *s, void *dm, void *)
{
	DataManager *pDM = (DataManager *)dm;

	if (pDM == 0) {
		silent_cerr("module-marble: DataManager unavailable (module_init() called too early?)" << std::endl);
		return 1;
	}
	
	pDM->SetLoadableElemModule(module_marble_lc.name, &module_marble_lc);

	return 0;
}
#endif /* !STATIC_MODULES */

}

