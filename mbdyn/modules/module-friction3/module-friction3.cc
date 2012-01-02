/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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

#include "loadable.h"
#include "friction.h"

/*
 * Usage:
 *
 *	loadable: <label>, <module name>, help [ , <module data> ] ;
 */


struct module_friction {
	BasicShapeCoefficient * Sh_c;
	BasicFriction * fc;
	doublereal mass;
	ScalarDof pos;
	ScalarDof mom;
	bool firstUpdate;
};

/* funzioni di default */
static void*
read(LoadableElem* pEl,
		DataManager* pDM,
		MBDynParser& HP)
{
	DEBUGCOUTFNAME("read");
	
	/* allocation of user-defined struct */
	module_friction* p = NULL;
	SAFENEW(p, module_friction);

	/*
	 * help
	 */
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	friction3						\n"
"Author: 	Marco Morandini <morandini@aero.polimi.it>		\n"
"Organization:	Dipartimento di Ingegneria Aerospaziale			\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it				\n"
"									\n"
"	All rights reserved						\n"
"									\n"
"	EXPERIMENTAL: used to test friction model			\n"
"									\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	p->firstUpdate = true;

	p->mom = ReadScalarDof(pDM, HP, 0);
	p->pos = ReadScalarDof(pDM, HP, 0);
	if (p->mom.pNode->GetDofType(0) != DofOrder::DIFFERENTIAL) {
		silent_cerr("Friction3(" << pEl->GetLabel() << "): "
			"need a differential node" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	if (p->pos.pNode->GetDofType(0) != DofOrder::DIFFERENTIAL) {
		silent_cerr("Friction3(" << pEl->GetLabel() << "): "
			"need a differential node" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	p->mass = HP.GetReal();

	if (HP.IsKeyWord("friction")) {
		p->fc = ParseFriction(HP, pDM);
		p->Sh_c = ParseShapeCoefficient(HP);
	} else {
		silent_cerr("Friction3(" << pEl->GetLabel() << "): "
			"perhaps you need a friction model; "
			"anyway, I'll go ahead" << std::endl);
		p->fc = 0;
	}

	return (void *)p;
}

static unsigned int
i_get_num_dof(const LoadableElem* pEl)
{
	DEBUGCOUTFNAME("i_get_num_dof");
	module_friction* p = (module_friction*)pEl->pGetData();
	if (p->fc) {
		return p->fc->iGetNumDof();
	}
	return 0;
}

DofOrder::Order
set_dof(const LoadableElem*pEl, unsigned int i)
{
	module_friction* p = (module_friction*)pEl->pGetData();
	if (p->fc) {
		return p->fc->GetDofType(i);
	}
	return DofOrder::ALGEBRAIC;
}

static void
output(const LoadableElem* pEl, OutputHandler& OH)
{
	DEBUGCOUTFNAME("output");
	
	if (pEl->fToBeOutput()) {
		module_friction* p = (module_friction *)pEl->pGetData();      
		std::ostream& out = OH.Loadable();

		out << std::setw(8) << pEl->GetLabel()
			<< " " << p->pos.pNode->dGetX()
			<< " " << p->pos.pNode->dGetXPrime()
			<< " " << p->mom.pNode->dGetX()
			<< " " << p->mom.pNode->dGetXPrime();
		if (p->fc) {
			out << " " << p->fc->fc();
		}
		out << std::endl;
	}
}

static void
work_space_dim(const LoadableElem* pEl, 
		    integer* piNumRows, 
		    integer* piNumCols)
{
	DEBUGCOUTFNAME("work_space_dim");
	module_friction* p = (module_friction*)pEl->pGetData();
	*piNumRows = 2;
	*piNumCols = 2;
	if (p->fc) {
		*piNumRows += p->fc->iGetNumDof();
		*piNumCols += p->fc->iGetNumDof();
	}
}

static VariableSubMatrixHandler& 
ass_jac(LoadableElem* pEl, 
	VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{  
	DEBUGCOUTFNAME("ass_jac");   
	FullSubMatrixHandler& WM = WorkMat.SetFull();
	
	module_friction* p = (module_friction*)pEl->pGetData();

	/*
	 * Dimensiono la matrice
	 */
	integer iNumRows = 0;
	integer iNumCols = 0;
	pEl->WorkSpaceDim(&iNumRows, &iNumCols);

	WM.ResizeReset(iNumRows, iNumCols);

	integer iFirstMomentumIndex = p->mom.pNode->iGetFirstRowIndex();
	integer iFirstPositionIndex = p->pos.pNode->iGetFirstRowIndex();
	integer iFirstReactionIndex = pEl->iGetFirstIndex();

	WM.PutRowIndex(1, iFirstMomentumIndex+1);
	WM.PutRowIndex(2, iFirstPositionIndex+1);
	WM.PutColIndex(1, iFirstMomentumIndex+1);
	WM.PutColIndex(2, iFirstPositionIndex+1);
	if (p->fc) {
	for (unsigned int i=1; i<=p->fc->iGetNumDof(); i++) {
		WM.PutRowIndex(2+i,iFirstReactionIndex+i);
		WM.PutColIndex(2+i,iFirstReactionIndex+i);
	}
	}


/*
 * 	dynamic system
 */
	// dot x matrix
	WM.PutCoef(1, 1, 1);
	WM.PutCoef(2, 2, -1);
	// x matrix
	WM.PutCoef(2, 1, dCoef/p->mass);
	
	if (p->fc) {
		doublereal f = p->fc->fc();
		doublereal shc = p->Sh_c->Sh_c();
		doublereal v = XPrimeCurr(iFirstPositionIndex+1);
		doublereal modF = p->mass;
		//doublereal F = shc*modF*f;
		
		ExpandableRowVector dfc;
		ExpandableRowVector dF;
		ExpandableRowVector dv;
		dF.ReDim(0);
		dv.ReDim(1);
		dv.Set(1,1,2);
		p->fc->AssJac(WorkMat,dfc,2,iFirstReactionIndex,dCoef,modF,v,
			XCurr,XPrimeCurr,dF,dv);
		ExpandableRowVector dFr;
		ExpandableRowVector dShc;
		p->Sh_c->dSh_c(dShc,f,modF,v,dfc,dF,dv);
		dFr.ReDim(2);
		dFr.Set(shc,1); dFr.Link(1,&dF);
		dFr.Set(F[0],2); dFr.Link(2,&dShc);
		dFr.Add(WorkMat,1,1.);
	}

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
	
	module_friction* p = (module_friction*)pEl->pGetData();

	/*
	 * Dimensiono il vettore
	 */
	integer iNumRows = 0;
	integer iNumCols = 0;
	pEl->WorkSpaceDim(&iNumRows, &iNumCols);
   
	WorkVec.ResizeReset(iNumRows);

	integer iFirstMomentumIndex = p->mom.pNode->iGetFirstRowIndex();
	integer iFirstPositionIndex = p->pos.pNode->iGetFirstRowIndex();
	integer iFirstReactionIndex = pEl->iGetFirstIndex();
	
	WorkVec.PutRowIndex(1,iFirstMomentumIndex+1);
	WorkVec.PutRowIndex(2,iFirstPositionIndex+1);
	if (p->fc) {
	for (unsigned int i=1; i<=p->fc->iGetNumDof(); i++) {
		WorkVec.PutRowIndex(2+i,iFirstReactionIndex+i);
	}
	}
	
	WorkVec.PutCoef(1,-XPrimeCurr(iFirstMomentumIndex+1));
	WorkVec.PutCoef(2,
		XPrimeCurr(iFirstPositionIndex+1));
	WorkVec.IncCoef(2,
		-XCurr(iFirstMomentumIndex+1)/p->mass);
	if (p->fc) {
		bool ChangeJac(false);
		doublereal v = XPrimeCurr(iFirstPositionIndex+1);
		try {
			p->fc->AssRes(WorkVec,2,iFirstReactionIndex,p->mass,
				v,XCurr,XPrimeCurr);
		}
		catch (Elem::ChangedEquationStructure) {
			ChangeJac = true;
		}
		doublereal f = p->fc->fc();
		doublereal shc = p->Sh_c->Sh_c(f,p->mass,v);
		WorkVec.Add(1,-shc*p->mass);
		if (ChangeJac) {
			throw Elem::ChangedEquationStructure();
		}
	}

	return WorkVec;
}

static void
before_predict(const LoadableElem* pEl, 
		VectorHandler& X,
		VectorHandler& XP,
		VectorHandler& XPrev,
		VectorHandler& XPPrev)
{
   DEBUGCOUTFNAME("before_predict");
}

static void
after_predict(const LoadableElem* pEl, 
		VectorHandler& X,
		VectorHandler& XP)
{
	DEBUGCOUTFNAME("after_predict");
}

static void
update(LoadableElem* pEl, 
		const VectorHandler& X,
		const VectorHandler& XP)
{
	DEBUGCOUTFNAME("update");
}

void
after_convergence(const LoadableElem* pEl, 
		const VectorHandler& X,
		const VectorHandler& XP)
{
	DEBUGCOUTFNAME("after_convergence");
	module_friction* p = (module_friction*)pEl->pGetData();
	if (p->fc) {
		integer iFirstPositionIndex = p->pos.pNode->iGetFirstRowIndex();
		doublereal v = XP(iFirstPositionIndex+1);
		p->fc->AfterConvergence(p->mass,v,X,XP,pEl->iGetFirstIndex());
	}
}

static void
destroy(LoadableElem* pEl)
{
	DEBUGCOUTFNAME("destroy");
	module_friction* p = (module_friction *)pEl->pGetData();
	
	SAFEDELETE(p);
}

static void
set_value(const LoadableElem* pEl, DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
	DEBUGCOUTFNAME("set_value");
	module_friction* p = (module_friction *)pEl->pGetData();
	if (p->fc) {
		p->fc->SetValue(X,XP,pEl->iGetFirstIndex());
	}
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
	NumNodes = i_get_num_connected_nodes(pEl);
}

static struct
LoadableCalls lc = {
	LOADABLE_VERSION_SET(1, 3, 0),

	"friction2",
	"1.0",
	"Dipartimento di Ingegneria Aerospaziale, Politecnico di Milano",
	"eXperimental friction model (2nd try)",

	read, /* */
	i_get_num_dof,
	set_dof,
	output, /* */
	NULL /* restart */ ,
	work_space_dim, /* */
	ass_jac, /* */
	NULL /* ass_mats */ ,
	ass_res, /* */
	before_predict,
	after_predict /* */ ,
	update /* */ ,
	after_convergence,
	NULL /* i_get_initial_num_dof */ ,
	NULL /* initial_work_space_dim */ ,
	NULL /* initial_ass_jac */ ,
	NULL /* initial_ass_res */ ,
	set_value /* */ ,
	NULL /* set_initial_value */ ,
	NULL /* i_get_num_priv_data */ ,
	NULL /* i_get_priv_data_idx */ ,
	NULL /* d_get_priv_data */ ,
	i_get_num_connected_nodes,
	get_connected_nodes,
	destroy
};

extern "C" {
void *calls = &lc;
}

