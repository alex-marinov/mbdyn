/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
 *
 * Pierangelo Masarati	<masarati@aero.polimi.it>
 * Paolo Mantegazza	<mantegazza@aero.polimi.it>
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
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#ifdef HAVE_LOADABLE
#include <unistd.h>
#include <loadable.h>
#include <dataman.h>

/* funzioni di default */
static unsigned int 
__i_get_num_dof(const LoadableElem* /* pEl */ )
{
   	return 0;
}

static DofOrder::Order 
__set_dof(const LoadableElem*, unsigned int /* i */ )
{
   	std::cerr << "You shouldn't be here!" << std::endl;
   	THROW(ErrGeneric());
#ifndef USE_EXCEPTIONS
   	return DofOrder::UNKNOWN;
#endif /* USE_EXCEPTIONS */
}

static void
__output(const LoadableElem* /* pEl */ , OutputHandler& /* OH */ )
{
   	NO_OP;
}

static std::ostream& 
__restart(const LoadableElem* pEl , std::ostream& out)
{
   	return out << "loadable: " << pEl->GetLabel() 
		<< ", not implemented yet;" << std::endl;
}

static void 
__work_space_dim(const LoadableElem* /* pEl */ ,
		 integer* piNumRows, 
		 integer* piNumCols)
{
   	*piNumRows = 0;
   	*piNumCols = 0;
}

static VariableSubMatrixHandler& 
__ass_jac(LoadableElem* /* pEl */ ,
	  VariableSubMatrixHandler& WorkMat,
	  doublereal /* dCoef */ ,
	  const VectorHandler& /* XCurr */ ,
	  const VectorHandler& /* XPrimeCurr */ )
{
   	WorkMat.SetNullMatrix();
   	return WorkMat;
}

static void 
__ass_mats(LoadableElem* /* pEl */ ,
	  VariableSubMatrixHandler& WorkMatA,
	  VariableSubMatrixHandler& WorkMatB,
	  const VectorHandler& /* XCurr */ ,
	  const VectorHandler& /* XPrimeCurr */ )
{
   	WorkMatA.SetNullMatrix();
   	WorkMatB.SetNullMatrix(); 
}

static SubVectorHandler& 
__ass_res(LoadableElem* /* pEl */ ,
	  SubVectorHandler& WorkVec,
	  doublereal /* dCoef */ ,
	  const VectorHandler& /* XCurr */ ,
	  const VectorHandler& /* XPrimeCurr */ )
{
   	WorkVec.Resize(0);
   	return WorkVec;
}

static void 
__before_predict(const LoadableElem* /* pEl */ ,
		 VectorHandler& /* X */ ,
		 VectorHandler& /* XP */ ,
		 VectorHandler& /* XPrev */ ,
		 VectorHandler& /* XPPrev */ )
{
   	NO_OP;
}

static void 
__after_predict(const LoadableElem* /* pEl */ ,
		VectorHandler& /* X */ ,
		VectorHandler& /* XP */ )
{
   	NO_OP;
}

static void 
__update(LoadableElem* /* pEl */ ,
	 const VectorHandler& /* X */ ,
	 const VectorHandler& /* XP */ )
{
   	NO_OP;
}

static void 
__after_convergence(const LoadableElem* /* pEl */ ,
		const VectorHandler& /* X */ ,
		const VectorHandler& /* XP */ )
{
   	NO_OP;
}

static unsigned int 
__i_get_initial_num_dof(const LoadableElem* /* pEl */ )
{
   	return 0;
}

static void 
__initial_work_space_dim(const LoadableElem* /* pEl */ ,
			 integer* piNumRows, 
			 integer* piNumCols)
{
   	*piNumRows = 0;
   	*piNumCols = 0;   
}

static VariableSubMatrixHandler& 
__initial_ass_jac(LoadableElem* /* pEl */ ,
		  VariableSubMatrixHandler& WorkMat, 
		  const VectorHandler& /* XCurr */ )
{
   	WorkMat.SetNullMatrix();
   	return WorkMat;
}

static SubVectorHandler& 
__initial_ass_res(LoadableElem* /* pEl */ ,
		  SubVectorHandler& WorkVec, 
		  const VectorHandler& /* XCurr */ )
{  
   	WorkVec.Resize(0);
   	return WorkVec;
}

static void 
__set_value(const LoadableElem* /* pEl */ , 
	    VectorHandler& /* X */ ,
	    VectorHandler& /* XP */ )
{
   	NO_OP;
}

static void 
__set_initial_value(const LoadableElem* /* pEl */ , VectorHandler& /* X */ )
{
   	NO_OP;
}

static unsigned int 
__i_get_num_priv_data(const LoadableElem* /* pEl */ )
{
   	return 0;
}

static unsigned int 
__i_get_priv_data_idx(const LoadableElem* /* pEl */ , const char *s)
{
   	std::cerr << "You shouldn't be here!" << std::endl;
   	THROW(ErrGeneric());
#ifndef USE_EXCEPTIONS
        return 0;
#endif /* USE_EXCEPTIONS */
}

static doublereal 
__d_get_priv_data(const LoadableElem* /* pEl */ , unsigned int /* i */ )
{
   	std::cerr << "You shouldn't be here!" << std::endl;
   	THROW(ErrGeneric());
#ifndef USE_EXCEPTIONS
        return 0.;
#endif /* USE_EXCEPTIONS */
}

static int
__i_get_num_connected_nodes(const LoadableElem* pEl)
{
	std::cerr << psElemNames[Elem::LOADABLE] << "(" << pEl->GetLabel()
		<< ") cannot be used in parallel environment" << std::endl;
	THROW(ErrGeneric());
}

static void
__get_connected_nodes(const LoadableElem* pEl, 
		int& NumNodes, 
		Node::Type* /* NdTyps */ , 
		unsigned int* /* NdLabels */ )
{
	NumNodes = __i_get_num_connected_nodes(pEl);
}

static void 
__destroy(LoadableElem* /* pEl */ )
{
   	NO_OP;
}

/* metodi della classe */

LoadableElem::LoadableElem(unsigned int uLabel, 
			   const DofOwner* pDO, 
			   DataManager* pDM, 
			   MBDynParser& HP)
: Elem(uLabel, Elem::LOADABLE, flag(0)),
#ifdef USE_STRUCT_NODES
InitialAssemblyElem(uLabel, Elem::LOADABLE, flag(0)),
#ifdef USE_AERODYNAMIC_ELEMS
AerodynamicElem(uLabel, AerodynamicElem::AERODYNAMICLOADABLE, flag(0)),
#endif /* USE_AERODYNAMIC_ELEMS */
ElemGravityOwner(uLabel, Elem::LOADABLE, flag(0)),
#endif /* USE_STRUCT_NODES */
ElemWithDofs(uLabel, Elem::LOADABLE, pDO, flag(0)),
priv_data(NULL),
module_name(NULL),
handle(NULL),
calls(NULL),
needsAirProperties(false)
{
   	ASSERT(pDM != NULL);

	GetCalls(HP);
	BindCalls(pDM, HP);
}

LoadableElem::LoadableElem(unsigned int uLabel, 
			   const DofOwner* pDO, 
			   LoadableCalls *c,
			   DataManager* pDM, 
			   MBDynParser& HP)
: Elem(uLabel, Elem::LOADABLE, flag(0)),
#ifdef USE_STRUCT_NODES
InitialAssemblyElem(uLabel, Elem::LOADABLE, flag(0)),
#ifdef USE_AERODYNAMIC_ELEMS
AerodynamicElem(uLabel, AerodynamicElem::AERODYNAMICLOADABLE, flag(0)),
#endif /* USE_AERODYNAMIC_ELEMS */
ElemGravityOwner(uLabel, Elem::LOADABLE, flag(0)),
#endif /* USE_STRUCT_NODES */
ElemWithDofs(uLabel, Elem::LOADABLE, pDO, flag(0)),
priv_data(NULL),
module_name(NULL),
handle(NULL),
calls(c),
needsAirProperties(false)
{
   	ASSERT(pDM != NULL);

	BindCalls(pDM, HP);
}

void
LoadableElem::GetCalls(MBDynParser& HP)
{
   	/* nome del modulo */
   	const char* s = HP.GetFileName();
	if (s == NULL) {
		std::cerr << "Loadable(" << GetLabel()
			<< "): unable to get module name" << std::endl;
		THROW(ErrGeneric());
	}

   	SAFESTRDUP(module_name, s);
#ifdef HAVE_LTDL_H
	handle = lt_dlopenext(module_name);
#elif defined(HAVE_DLFCN_H)
   	handle = dlopen(module_name, RTLD_NOW /* RTLD_LAZY */ );
#endif /* !HAVE_LTDL_H && HAVE_DLFCN_H */

	if (handle == NULL) {
      		std::cerr << "Loadable(" << uLabel 
			<< "): unable to open module <" << module_name 
#ifdef HAVE_LTDL_H
			<< "> (lt_dlopen returns \"" << lt_dlerror() << "\")"
#elif defined(HAVE_DLFCN_H)
			<< "> (dlopen returns \"" << dlerror() << "\")" 
#endif /* !HAVE_LTDL_H && HAVE_DLFCN_H */
			<< std::endl;
      		THROW(ErrGeneric());
   	}

	/* default LoadableCalls struct */
	const char *data_name = "calls";

	LoadableCalls **tmpcalls = NULL;
	if (HP.IsKeyWord("name")) {
		data_name = HP.GetStringWithDelims();
	}
   	DEBUGCOUT("binding to data \"" << data_name
     		<< "\" (must be def'd!)" << std::endl);
#ifdef HAVE_LTDL_H
	tmpcalls = (LoadableCalls **)lt_dlsym(handle, data_name);
#elif defined(HAVE_DLFCN_H)
	tmpcalls = (LoadableCalls **)dlsym(handle, data_name);
#endif /* !HAVE_LTDL_H && HAVE_DLFCN_H */
	
   	if (tmpcalls == NULL) {
      		const char* err =
#ifdef HAVE_LTDL_H
			lt_dlerror();
#elif defined(HAVE_DLFCN_H)
			dlerror();
#endif /* !HAVE_LTDL_H && HAVE_DLFCN_H */
      		if (err == NULL) {
	 		std::cerr << "Loadable(" << uLabel 
	   			<< "): data \"" << data_name
	   			<< "\" must be defined in module <" 
				<< module_name << ">" << std::endl;
      		} else {
	 		std::cerr << "Loadable(" << uLabel
	   			<< "): error while binding to data \"" 
				<< data_name
	   			<< "\" in module <" << module_name
	   			<< "> (\"" << err 
				<< "\")" << std::endl;
      		}
      		THROW(ErrGeneric());
   	}

	calls = *tmpcalls;
}

void
LoadableElem::BindCalls(DataManager* pDM, MBDynParser& HP)
{
	if (calls->loadable_version != LOADABLE_VERSION) {
		std::cerr << "Loadable(" << uLabel
			<< "): incompatible version; need "
			<< LOADABLE_VERSION_OUT(LOADABLE_VERSION)
			<< ", got " 
			<< LOADABLE_VERSION_OUT(calls->loadable_version) 
			<< std::endl;
		THROW(ErrGeneric());
	}
   
	if (calls->read == NULL) {
		std::cerr << "Loadable(" << uLabel
			<< "): function \"read\" must be defined in module <"
			<< module_name << "> data" << std::endl;
		THROW(ErrGeneric());
	}

	if (calls->name == NULL) {
		calls->name = module_name;
	}

	if (calls->version == NULL) {
		calls->version = "(undefined)";
	}

	if (calls->vendor == NULL) {
		calls->vendor = "(undefined)";
	}

	if (calls->description == NULL) {
		calls->description = "";
	}

	silent_cout("Loadable(" << uLabel << "): " << calls->name
			<< " version " << calls->version << std::endl
			<< "\tvendor: " << calls->vendor << std::endl
			<< "\tdescription: " << calls->description
			<< std::endl);

	/*
	 * Mette i default ove servono
	 */
	if (calls->i_get_num_dof == NULL) {
		calls->i_get_num_dof = __i_get_num_dof;
	}

	if (calls->set_dof == NULL) {
		calls->set_dof = __set_dof;
	}

	if (calls->output == NULL) {
		calls->output = __output;
	}

	if (calls->restart == NULL) {
		calls->restart = __restart;
	}

	if (calls->work_space_dim == NULL) {
		calls->work_space_dim = __work_space_dim;
	}

	if (calls->ass_jac == NULL) {
		calls->ass_jac = __ass_jac;
	}

	if (calls->ass_mats == NULL) {
		calls->ass_mats = __ass_mats;
	}

	if (calls->ass_res == NULL) {
		calls->ass_res = __ass_res;
	}

	if (calls->before_predict == NULL) {
		calls->before_predict = __before_predict;
	}

	if (calls->after_predict == NULL) {
		calls->after_predict = __after_predict;
	}

	if (calls->update == NULL) {
		calls->update = __update;
	}

	if (calls->after_convergence == NULL) {
		calls->after_convergence = __after_convergence;
	}

	if (calls->i_get_initial_num_dof == NULL) {
		calls->i_get_initial_num_dof = __i_get_initial_num_dof;
	}

	if (calls->initial_work_space_dim == NULL) {
		calls->initial_work_space_dim = __initial_work_space_dim;
	}

	if (calls->initial_ass_jac == NULL) {
		calls->initial_ass_jac = __initial_ass_jac;
	}

	if (calls->initial_ass_res == NULL) {
		calls->initial_ass_res = __initial_ass_res;
	}

	if (calls->set_value == NULL) {
		calls->set_value = __set_value;
	}

	if (calls->set_initial_value == NULL) {
		calls->set_initial_value = __set_initial_value;
	}

	if (calls->i_get_num_priv_data == NULL) {
		calls->i_get_num_priv_data = __i_get_num_priv_data;
	}

	if (calls->i_get_priv_data_idx == NULL) {
		calls->i_get_priv_data_idx = __i_get_priv_data_idx;
	}

	if (calls->d_get_priv_data == NULL) {
		calls->d_get_priv_data = __d_get_priv_data;
	}

	if (calls->i_get_num_connected_nodes == NULL) {
		calls->i_get_num_connected_nodes = __i_get_num_connected_nodes;
	}

	if (calls->get_connected_nodes == NULL) {
		calls->get_connected_nodes = __get_connected_nodes;
	}

	if (calls->destroy == NULL) {
		calls->destroy = __destroy;
	}
	
   	priv_data = (*calls->read)(this, pDM, HP, pDM->pGetDrvHdl());
   	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE)); 
}

LoadableElem::~LoadableElem(void)
{
	ASSERT(calls->destroy != NULL);
   	(*calls->destroy)(this);
   
#if !defined(HAVE_LTDL_H) && defined(HAVE_DLFCN_H)
   	if (handle != NULL) {
   		if (dlclose(handle) != 0) {
			std::cerr << "unable to dlclose module \"" 
				<< module_name << "\"" << std::endl;
			THROW(ErrGeneric());
		}
	}
#endif /* !HAVE_LTDL_H && HAVE_DLFCN_H */
	
   	SAFEDELETEARR(module_name);
}

Elem::Type 
LoadableElem::GetElemType(void) const
{
   	return Elem::LOADABLE;
}

unsigned int 
LoadableElem::iGetNumDof(void) const
{
	ASSERT(calls->i_get_num_dof != NULL);
   	return (*calls->i_get_num_dof)(this);
}

DofOrder::Order 
LoadableElem::GetDofType(unsigned int i) const
{
   	ASSERT(i < iGetNumDof());
	ASSERT(calls->set_dof != NULL);
   	return (*calls->set_dof)(this, i);
}

void 
LoadableElem::Output(OutputHandler& OH) const
{
	ASSERT(calls->output != NULL);
   	(*calls->output)(this, OH);
}

std::ostream& 
LoadableElem::Restart(std::ostream& out) const
{
	ASSERT(calls->restart != NULL);
   	out << "    loadable: " << GetLabel() << ", \"" 
		<< module_name << "\", ";
   	return (*calls->restart)(this, out) << ';' << std::endl;
}

void 
LoadableElem::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	ASSERT(calls->work_space_dim != NULL);
   	(*calls->work_space_dim)(this, piNumRows, piNumCols);
}

VariableSubMatrixHandler& 
LoadableElem::AssJac(VariableSubMatrixHandler& WorkMat,
		     doublereal dCoef, 
		     const VectorHandler& XCurr,
		     const VectorHandler& XPCurr)
{
	ASSERT(calls->ass_jac != NULL);
   	return (*calls->ass_jac)(this, WorkMat, dCoef, XCurr, XPCurr);
}

void
LoadableElem::AssMats(VariableSubMatrixHandler& WorkMatA,
		     VariableSubMatrixHandler& WorkMatB,
		     const VectorHandler& XCurr,
		     const VectorHandler& XPCurr)
{
   	ASSERT(calls->ass_mats != NULL);
   	(*calls->ass_mats)(this, WorkMatA, WorkMatB, XCurr, XPCurr);
}

SubVectorHandler& 
LoadableElem::AssRes(SubVectorHandler& WorkVec,
		     doublereal dCoef,
		     const VectorHandler& XCurr, 
		     const VectorHandler& XPCurr)
{
   	ASSERT(calls->ass_res != NULL);
   	return (*calls->ass_res)(this, WorkVec, 
					    dCoef, XCurr, XPCurr);
}

void 
LoadableElem::BeforePredict(VectorHandler& X,
			    VectorHandler& XP,
			    VectorHandler& XPrev,
			    VectorHandler& XPPrev) const
{
   	ASSERT(calls->before_predict != NULL);
   	(*calls->before_predict)(this, X, XP, XPrev, XPPrev);
}

void 
LoadableElem::AfterPredict(VectorHandler& X,
			   VectorHandler& XP)
{
   	ASSERT(calls->after_predict != NULL);
   	(*calls->after_predict)(this, X, XP);
}

void 
LoadableElem::Update(const VectorHandler& XCurr, 
		     const VectorHandler& XPrimeCurr)
{
   	ASSERT(calls->update != NULL);
   	(*calls->update)(this, XCurr, XPrimeCurr);
}

void 
LoadableElem::AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP)
{
   	ASSERT(calls->after_convergence != NULL);
   	(*calls->after_convergence)(this, X, XP);
}

#ifdef USE_STRUCT_NODES
unsigned int 
LoadableElem::iGetInitialNumDof(void) const
{
   	ASSERT(calls->i_get_initial_num_dof != NULL);
   	return (*calls->i_get_initial_num_dof)(this);
}

void 
LoadableElem::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
   	ASSERT(calls->initial_work_space_dim != NULL);
   	(*calls->initial_work_space_dim)(this, piNumRows, piNumCols);
}

VariableSubMatrixHandler& 
LoadableElem::InitialAssJac(VariableSubMatrixHandler& WorkMat,
			    const VectorHandler& XCurr)
{
   	ASSERT(calls->initial_ass_jac != NULL);
   	return (*calls->initial_ass_jac)(this, WorkMat, XCurr);
}

SubVectorHandler& 
LoadableElem::InitialAssRes(SubVectorHandler& WorkVec, 
			    const VectorHandler& XCurr)
{
   	ASSERT(calls->initial_ass_res != NULL);
   	return (*calls->initial_ass_res)(this, WorkVec, XCurr);
}

void 
LoadableElem::SetInitialValue(VectorHandler& X) const
{   
   	ASSERT(calls->set_initial_value != NULL);
   	(*calls->set_initial_value)(this, X);
}
#endif /* USE_STRUCT_NODES */

void 
LoadableElem::SetValue(VectorHandler& X, VectorHandler& XP) const
{
   	ASSERT(calls->set_value != NULL);
   	(*calls->set_value)(this, X, XP);
}

unsigned int 
LoadableElem::iGetNumPrivData(void) const
{
   	ASSERT(calls->i_get_num_priv_data != NULL);
   	return (*calls->i_get_num_priv_data)(this);
}

unsigned int 
LoadableElem::iGetPrivDataIdx(const char *s) const
{
   	ASSERT(calls->i_get_priv_data_idx != NULL);
   	return (*calls->i_get_priv_data_idx)(this, s);
}

doublereal 
LoadableElem::dGetPrivData(unsigned int i) const
{
   	ASSERT(calls->d_get_priv_data != NULL);
   	return (*calls->d_get_priv_data)(this, i);
}

int
LoadableElem::GetNumConnectedNodes(void) const
{
	ASSERT(calls->i_get_num_connected_nodes != NULL);
	return (*calls->i_get_num_connected_nodes)(this);
}

void
LoadableElem::GetConnectedNodes(int& NumNodes, 
		Node::Type* NdTyps, 
		unsigned int* NdLabels)
{
	ASSERT(calls->get_connected_nodes != NULL);
	return (*calls->get_connected_nodes)(this, NumNodes, NdTyps, NdLabels);
}

Elem* 
ReadLoadable(DataManager* pDM, 
	     MBDynParser& HP, 
	     const DofOwner* pDO, 
	     unsigned int uLabel)
{
   	Elem* pEl = NULL;
   
   	SAFENEWWITHCONSTRUCTOR(pEl,
			       LoadableElem,
			       LoadableElem(uLabel, pDO, pDM, HP));
   	return pEl;
}

bool
LoadableElem::NeedsAirProperties(void) const
{
	return needsAirProperties;
}

void
LoadableElem::NeedsAirProperties(bool yesno)
{
	needsAirProperties = yesno;
}

#endif /* HAVE_LOADABLE */

