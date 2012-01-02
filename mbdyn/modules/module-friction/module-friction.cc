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

#if 0
static StepFriction<doublereal> sf(1.);
static TanhFriction<doublereal> tf(1., 1.);
#endif

struct module_friction {
	bool			firstUpdate;
	Friction<doublereal> 	*f;
	ScalarDof		dof;
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
"Module: 	friction						\n"
"Author: 	Stefania Gualdi <gualdi@aero.polimi.it>			\n"
"		Marco Morandini <morandini@aero.polimi.it>		\n"
"		Pierangelo Masarati <masarati@aero.polimi.it>		\n"
"Organization:	Dipartimento di Ingegneria Aerospaziale			\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it				\n"
"									\n"
"	All rights reserved						\n"
"									\n"
"	EXPERIMENTAL: used to test different friction models		\n"
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

	p->dof = ReadScalarDof(pDM, HP, 0);
	if (p->dof.pNode->GetDofType(0) != DofOrder::DIFFERENTIAL) {
		silent_cerr("Friction(" << pEl->GetLabel() << "): "
			"need a differential node" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	p->f = NULL;
	if (HP.IsKeyWord("step")) {
		doublereal mu = HP.GetReal();

		typedef StepFriction<doublereal> S;
		SAFENEWWITHCONSTRUCTOR(p->f, S, S(mu));

	} else if (HP.IsKeyWord("tanh")) {
		doublereal mu = HP.GetReal();
		doublereal vrif = HP.GetReal();

		typedef TanhFriction<doublereal> S;
		SAFENEWWITHCONSTRUCTOR(p->f, S, S(mu, vrif));

	} else if (HP.IsKeyWord("discretestate")) {
		doublereal maxForce = HP.GetReal();
		doublereal stiffness = HP.GetReal();
		doublereal damping = HP.GetReal();
		doublereal velTreshold = HP.GetReal();

		DiscStateFriction<doublereal>::State state = DiscStateFriction<doublereal>::Stick;
		if (HP.IsKeyWord("initialstate")) {
			if (HP.IsKeyWord("stick")) {
				state = DiscStateFriction<doublereal>::Stick;
			} else if (HP.IsKeyWord("slip")) {
				state = DiscStateFriction<doublereal>::Slip;
			} else {
				silent_cerr("Friction(" << pEl->GetLabel() << "): "
					"unknown state for discrete state friction model"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
		
		typedef DiscStateFriction<doublereal> S;
		SAFENEWWITHCONSTRUCTOR(p->f, S, S(maxForce, stiffness, 
					damping, velTreshold, state));

	} else {
		silent_cerr("Friction(" << pEl->GetLabel() << "): "
			"unknown friction model" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return (void *)p;
}

static unsigned int
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
		module_friction* p = (module_friction *)pEl->pGetData();      
		std::ostream& out = OH.Loadable();

		out << std::setw(8) << pEl->GetLabel() << " " << p->f->F() << std::endl;
	}
}

static std::ostream&
restart(const LoadableElem* pEl, std::ostream& out)
{
	DEBUGCOUTFNAME("restart");
	return out << "not implemented yet;" << std::endl;
}

static void
work_space_dim(const LoadableElem* pEl, 
		    integer* piNumRows, 
		    integer* piNumCols)
{
	DEBUGCOUTFNAME("work_space_dim");
	*piNumRows = 1;
	*piNumCols = 1;
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
	
	module_friction* p = (module_friction*)pEl->pGetData();

	/*
	 * Dimensiono il vettore
	 */
	integer iNumRows = 0;
	integer iNumCols = 0;
	pEl->WorkSpaceDim(&iNumRows, &iNumCols);
   
	WorkVec.Resize(iNumRows);
	WorkVec.Reset(0.);

	integer iIndex = p->dof.pNode->iGetFirstRowIndex() + 1;

	//std::cerr << "firstUpdate: " << p->firstUpdate;

	if (p->firstUpdate) {
		p->firstUpdate = false;
	} else {
		p->f->Update(1., p->dof.pNode->dGetX(), p->dof.pNode->dGetXPrime());
	}

	//std::cerr << "; F=" << p->f->F() << std::endl;

	WorkVec.PutItem(1, iIndex, p->f->F());

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

	module_friction* p = (module_friction*)pEl->pGetData();

	p->firstUpdate = true;
	p->f->Update(1., p->dof.pNode->dGetX(), p->dof.pNode->dGetXPrime(), Friction<doublereal>::FIRST);

	// std::cerr << "update; F=" << p->f->F() << std::endl;
}

static void
update(LoadableElem* pEl, 
		const VectorHandler& X,
		const VectorHandler& XP)
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
initial_work_space_dim(const LoadableElem* pEl, 
		integer* piNumRows, 
		integer* piNumCols)
{
	DEBUGCOUTFNAME("initial_work_space_dim");
	*piNumRows = 0;
	*piNumCols = 0;   
}

static VariableSubMatrixHandler& 
initial_ass_jac(LoadableElem* pEl, 
		VariableSubMatrixHandler& WorkMat, 
		const VectorHandler& XCurr)
{
	DEBUGCOUTFNAME("initial_ass_jac");
	integer iNumRows = 0;
	integer iNumCols = 0;
	pEl->InitialWorkSpaceDim(&iNumRows, &iNumCols);
   
	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeReset(iNumRows, iNumCols, 0.);

	return WorkMat;
}

static SubVectorHandler& 
initial_ass_res(LoadableElem* pEl, 
		SubVectorHandler& WorkVec, 
		const VectorHandler& XCurr)
{
	DEBUGCOUTFNAME("initial_ass_res");
	integer iNumRows = 0;
	integer iNumCols = 0;
	pEl->WorkSpaceDim(&iNumRows, &iNumCols);
	
	WorkVec.Resize(iNumRows);

	/* set sub-vector indices and coefs */
   
	return WorkVec;
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
	module_friction* p = (module_friction *)pEl->pGetData();
	
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
	NumNodes = i_get_num_connected_nodes(pEl);
}

static struct
LoadableCalls lc = {
	LOADABLE_VERSION_SET(1, 3, 0),

	"friction",
	"1.1",
	"Dipartimento di Ingegneria Aerospaziale, Politecnico di Milano",
	"eXperimental friction model",

	read, /* */
	NULL /* i_get_num_dof */ ,
	NULL /* set_dof */ ,
	output, /* */
	NULL /* restart */ ,
	work_space_dim, /* */
	ass_jac, /* */
	NULL /* ass_mats */ ,
	ass_res, /* */
	NULL /* before_predict */ ,
	after_predict /* */ ,
	NULL /* update */ ,
	NULL /* after_convergence */ ,
	NULL /* i_get_initial_num_dof */ ,
	NULL /* initial_work_space_dim */ ,
	NULL /* initial_ass_jac */ ,
	NULL /* initial_ass_res */ ,
	NULL /* set_value */ ,
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

