/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <unistd.h>

#include "dataman.h"
#include "loadable.h"

/* funzioni di default */
static unsigned int
int_i_get_num_dof(const LoadableElem* /* pEl */ )
{
   	return 0;
}

static DofOrder::Order
int_set_dof(const LoadableElem*, unsigned int /* i */ )
{
   	silent_cerr("You shouldn't be here!" << std::endl);
   	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

static void
int_output(const LoadableElem* /* pEl */ , OutputHandler& /* OH */ )
{
   	NO_OP;
}

static std::ostream&
int_restart(const LoadableElem* pEl , std::ostream& out)
{
   	return out << "loadable: " << pEl->GetLabel()
		<< ", not implemented yet;" << std::endl;
}

static void
int_work_space_dim(const LoadableElem* /* pEl */ ,
		 integer* piNumRows,
		 integer* piNumCols)
{
   	*piNumRows = 0;
   	*piNumCols = 0;
}

static VariableSubMatrixHandler&
int_ass_jac(LoadableElem* /* pEl */ ,
	  VariableSubMatrixHandler& WorkMat,
	  doublereal /* dCoef */ ,
	  const VectorHandler& /* XCurr */ ,
	  const VectorHandler& /* XPrimeCurr */ )
{
   	WorkMat.SetNullMatrix();
   	return WorkMat;
}

static void
int_ass_mats(LoadableElem* /* pEl */ ,
	  VariableSubMatrixHandler& WorkMatA,
	  VariableSubMatrixHandler& WorkMatB,
	  const VectorHandler& /* XCurr */ ,
	  const VectorHandler& /* XPrimeCurr */ )
{
   	WorkMatA.SetNullMatrix();
   	WorkMatB.SetNullMatrix();
}

static SubVectorHandler&
int_ass_res(LoadableElem* /* pEl */ ,
	  SubVectorHandler& WorkVec,
	  doublereal /* dCoef */ ,
	  const VectorHandler& /* XCurr */ ,
	  const VectorHandler& /* XPrimeCurr */ )
{
   	WorkVec.Resize(0);
   	return WorkVec;
}

static void
int_before_predict(const LoadableElem* /* pEl */ ,
		 VectorHandler& /* X */ ,
		 VectorHandler& /* XP */ ,
		 VectorHandler& /* XPrev */ ,
		 VectorHandler& /* XPPrev */ )
{
   	NO_OP;
}

static void
int_after_predict(const LoadableElem* /* pEl */ ,
		VectorHandler& /* X */ ,
		VectorHandler& /* XP */ )
{
   	NO_OP;
}

static void
int_update(LoadableElem* /* pEl */ ,
	 const VectorHandler& /* X */ ,
	 const VectorHandler& /* XP */ )
{
   	NO_OP;
}

static void
int_after_convergence(const LoadableElem* /* pEl */ ,
		const VectorHandler& /* X */ ,
		const VectorHandler& /* XP */ )
{
   	NO_OP;
}

static unsigned int
int_i_get_initial_num_dof(const LoadableElem* /* pEl */ )
{
   	return 0;
}

static void
int_initial_work_space_dim(const LoadableElem* /* pEl */ ,
			 integer* piNumRows,
			 integer* piNumCols)
{
   	*piNumRows = 0;
   	*piNumCols = 0;
}

static VariableSubMatrixHandler&
int_initial_ass_jac(LoadableElem* /* pEl */ ,
		  VariableSubMatrixHandler& WorkMat,
		  const VectorHandler& /* XCurr */ )
{
   	WorkMat.SetNullMatrix();
   	return WorkMat;
}

static SubVectorHandler&
int_initial_ass_res(LoadableElem* /* pEl */ ,
		  SubVectorHandler& WorkVec,
		  const VectorHandler& /* XCurr */ )
{
   	WorkVec.Resize(0);
   	return WorkVec;
}

static void
int_set_value(const LoadableElem* /* pEl */ ,
		DataManager *pDM,
		VectorHandler& /* X */ ,
		VectorHandler& /* XP */ ,
		SimulationEntity::Hints *ph)
{
   	NO_OP;
}

static void
int_set_initial_value(const LoadableElem* /* pEl */ , VectorHandler& /* X */ )
{
   	NO_OP;
}

static unsigned int
int_i_get_num_priv_data(const LoadableElem* /* pEl */ )
{
   	return 0;
}

static unsigned int
int_i_get_priv_data_idx(const LoadableElem* /* pEl */ , const char *s)
{
   	silent_cerr("You shouldn't be here!" << std::endl);
   	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

static doublereal
int_d_get_priv_data(const LoadableElem* /* pEl */ , unsigned int /* i */ )
{
   	silent_cerr("You shouldn't be here!" << std::endl);
   	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

static int
int_i_get_num_connected_nodes(const LoadableElem* pEl)
{
	return 0;
}

static void
int_get_connected_nodes(const LoadableElem* pEl, std::vector<const Node *>& connectedNodes)
{
	connectedNodes.resize(0);
}

static void
int_destroy(LoadableElem* /* pEl */ )
{
   	NO_OP;
}

/* metodi della classe */

LoadableElem::LoadableElem(unsigned int uLabel,
			   const DofOwner* pDO,
			   DataManager* pDM,
			   MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO),
priv_data(0),
module_name(0),
#ifdef USE_RUNTIME_LOADING
handle(0),
#endif // USE_RUNTIME_LOADING
calls(0)
{
   	ASSERT(pDM != 0);

	GetCalls(HP);
	BindCalls(pDM, HP);
}

LoadableElem::LoadableElem(unsigned int uLabel,
			   const DofOwner* pDO,
			   const LoadableCalls *c,
			   DataManager* pDM,
			   MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO),
priv_data(0),
module_name(0),
#ifdef USE_RUNTIME_LOADING
handle(0),
#endif // USE_RUNTIME_LOADING
calls(const_cast<LoadableCalls *>(c))
{
   	ASSERT(pDM != NULL);

	BindCalls(pDM, HP);
}

void
LoadableElem::GetCalls(MBDynParser& HP)
{
#ifdef USE_RUNTIME_LOADING
   	/* nome del modulo */
   	const char* s = HP.GetFileName();
	if (s == NULL) {
		silent_cerr("Loadable(" << GetLabel()
			<< "): unable to get module name" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

   	SAFESTRDUP(module_name, s);
	handle = lt_dlopenext(module_name);

	if (handle == NULL) {
		const char *err = lt_dlerror();
		if (err == 0) {
			err = "";
		}

      		silent_cerr("Loadable(" << uLabel
			<< "): unable to open module <" << module_name
			<< "> (" << err << ") at line " << HP.GetLineData()
			<< std::endl);
      		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   	}

	/* default LoadableCalls struct */
	const char *data_name = "calls";

	LoadableCalls **tmpcalls = NULL;
	if (HP.IsKeyWord("name")) {
		data_name = HP.GetStringWithDelims();
	}
   	DEBUGCOUT("binding to data \"" << data_name
     		<< "\" (must be def'd!)" << std::endl);
	tmpcalls = (LoadableCalls **)lt_dlsym(handle, data_name);

   	if (tmpcalls == NULL) {
      		const char* err = lt_dlerror();
      		if (err == NULL) {
	 		silent_cerr("Loadable(" << uLabel
	   			<< "): data \"" << data_name
	   			<< "\" must be defined in module <"
				<< module_name << ">" << std::endl);
      		} else {
	 		silent_cerr("Loadable(" << uLabel
	   			<< "): error while binding to data \""
				<< data_name
	   			<< "\" in module <" << module_name
	   			<< "> (\"" << err
				<< "\")" << std::endl);
      		}
      		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   	}

	calls = *tmpcalls;
#else // !USE_RUNTIME_LOADING
	silent_cerr("LoadableElem(" << GetLabel() << ") GetCalls: "
		"should not be called when --disable-runtime-loading" << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif // ! USE_RUNTIME_LOADING
}

void
LoadableElem::BindCalls(DataManager* pDM, MBDynParser& HP)
{
	if (calls->loadable_version != LOADABLE_VERSION) {
		silent_cerr("Loadable(" << uLabel
			<< "): incompatible version; need "
			<< LOADABLE_VERSION_OUT(LOADABLE_VERSION)
			<< ", got "
			<< LOADABLE_VERSION_OUT(calls->loadable_version)
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (calls->read == NULL) {
		silent_cerr("Loadable(" << uLabel
			<< "): function \"read\" must be defined in module <"
			<< module_name << "> data" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
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
		calls->i_get_num_dof = int_i_get_num_dof;
	}

	if (calls->set_dof == NULL) {
		calls->set_dof = int_set_dof;
	}

	if (calls->output == NULL) {
		calls->output = int_output;
	}

	if (calls->restart == NULL) {
		calls->restart = int_restart;
	}

	if (calls->work_space_dim == NULL) {
		calls->work_space_dim = int_work_space_dim;
	}

	if (calls->ass_jac == NULL) {
		calls->ass_jac = int_ass_jac;
	}

	if (calls->ass_mats == NULL) {
		calls->ass_mats = int_ass_mats;
	}

	if (calls->ass_res == NULL) {
		calls->ass_res = int_ass_res;
	}

	if (calls->before_predict == NULL) {
		calls->before_predict = int_before_predict;
	}

	if (calls->after_predict == NULL) {
		calls->after_predict = int_after_predict;
	}

	if (calls->update == NULL) {
		calls->update = int_update;
	}

	if (calls->after_convergence == NULL) {
		calls->after_convergence = int_after_convergence;
	}

	if (calls->i_get_initial_num_dof == NULL) {
		calls->i_get_initial_num_dof = int_i_get_initial_num_dof;
	}

	if (calls->initial_work_space_dim == NULL) {
		calls->initial_work_space_dim = int_initial_work_space_dim;
	}

	if (calls->initial_ass_jac == NULL) {
		calls->initial_ass_jac = int_initial_ass_jac;
	}

	if (calls->initial_ass_res == NULL) {
		calls->initial_ass_res = int_initial_ass_res;
	}

	if (calls->set_value == NULL) {
		calls->set_value = int_set_value;
	}

	if (calls->set_initial_value == NULL) {
		calls->set_initial_value = int_set_initial_value;
	}

	if (calls->i_get_num_priv_data == NULL) {
		calls->i_get_num_priv_data = int_i_get_num_priv_data;
	}

	if (calls->i_get_priv_data_idx == NULL) {
		calls->i_get_priv_data_idx = int_i_get_priv_data_idx;
	}

	if (calls->d_get_priv_data == NULL) {
		calls->d_get_priv_data = int_d_get_priv_data;
	}

	if (calls->i_get_num_connected_nodes == NULL) {
		calls->i_get_num_connected_nodes = int_i_get_num_connected_nodes;
	}

	if (calls->get_connected_nodes == NULL) {
		calls->get_connected_nodes = int_get_connected_nodes;
	}

	if (calls->destroy == NULL) {
		calls->destroy = int_destroy;
	}

   	priv_data = calls->read(this, pDM, HP);
   	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
}

LoadableElem::~LoadableElem(void)
{
	ASSERT(calls->destroy != NULL);
   	calls->destroy(this);

#ifdef USE_RUNTIME_LOADING
   	if (handle != NULL) {
   		if (lt_dlclose(handle) != 0) {
			silent_cerr("unable to close module "
				"\"" << module_name << "\"" << std::endl);
		}
	}
#endif // USE_RUNTIME_LOADING

   	SAFEDELETEARR(module_name);
}

unsigned int
LoadableElem::iGetNumDof(void) const
{
	ASSERT(calls->i_get_num_dof != NULL);
   	return calls->i_get_num_dof(this);
}

DofOrder::Order
LoadableElem::GetDofType(unsigned int i) const
{
   	ASSERT(i < iGetNumDof());
	ASSERT(calls->set_dof != NULL);
   	return calls->set_dof(this, i);
}

void
LoadableElem::Output(OutputHandler& OH) const
{
	ASSERT(calls->output != NULL);
   	calls->output(this, OH);
}

std::ostream&
LoadableElem::Restart(std::ostream& out) const
{
	ASSERT(calls->restart != NULL);
   	out << "    loadable: " << GetLabel() << ", \""
		<< module_name << "\", ";
   	return calls->restart(this, out) << ';' << std::endl;
}

void
LoadableElem::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	ASSERT(calls->work_space_dim != NULL);
   	calls->work_space_dim(this, piNumRows, piNumCols);
}

VariableSubMatrixHandler&
LoadableElem::AssJac(VariableSubMatrixHandler& WorkMat,
		     doublereal dCoef,
		     const VectorHandler& XCurr,
		     const VectorHandler& XPCurr)
{
	ASSERT(calls->ass_jac != NULL);
   	return calls->ass_jac(this, WorkMat, dCoef, XCurr, XPCurr);
}

void
LoadableElem::AssMats(VariableSubMatrixHandler& WorkMatA,
		     VariableSubMatrixHandler& WorkMatB,
		     const VectorHandler& XCurr,
		     const VectorHandler& XPCurr)
{
   	ASSERT(calls->ass_mats != NULL);
   	calls->ass_mats(this, WorkMatA, WorkMatB, XCurr, XPCurr);
}

SubVectorHandler&
LoadableElem::AssRes(SubVectorHandler& WorkVec,
		     doublereal dCoef,
		     const VectorHandler& XCurr,
		     const VectorHandler& XPCurr)
{
   	ASSERT(calls->ass_res != NULL);
   	return calls->ass_res(this, WorkVec,
					    dCoef, XCurr, XPCurr);
}

void
LoadableElem::BeforePredict(VectorHandler& X,
			    VectorHandler& XP,
			    VectorHandler& XPrev,
			    VectorHandler& XPPrev) const
{
   	ASSERT(calls->before_predict != NULL);
   	calls->before_predict(this, X, XP, XPrev, XPPrev);
}

void
LoadableElem::AfterPredict(VectorHandler& X,
			   VectorHandler& XP)
{
   	ASSERT(calls->after_predict != NULL);
   	calls->after_predict(this, X, XP);
}

void
LoadableElem::Update(const VectorHandler& XCurr,
		     const VectorHandler& XPrimeCurr)
{
   	ASSERT(calls->update != NULL);
   	calls->update(this, XCurr, XPrimeCurr);
}

void
LoadableElem::AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP)
{
   	ASSERT(calls->after_convergence != NULL);
   	calls->after_convergence(this, X, XP);
}

unsigned int
LoadableElem::iGetInitialNumDof(void) const
{
   	ASSERT(calls->i_get_initial_num_dof != NULL);
   	return calls->i_get_initial_num_dof(this);
}

void
LoadableElem::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
   	ASSERT(calls->initial_work_space_dim != NULL);
   	calls->initial_work_space_dim(this, piNumRows, piNumCols);
}

VariableSubMatrixHandler&
LoadableElem::InitialAssJac(VariableSubMatrixHandler& WorkMat,
			    const VectorHandler& XCurr)
{
   	ASSERT(calls->initial_ass_jac != NULL);
   	return calls->initial_ass_jac(this, WorkMat, XCurr);
}

SubVectorHandler&
LoadableElem::InitialAssRes(SubVectorHandler& WorkVec,
			    const VectorHandler& XCurr)
{
   	ASSERT(calls->initial_ass_res != NULL);
   	return calls->initial_ass_res(this, WorkVec, XCurr);
}

void
LoadableElem::SetInitialValue(VectorHandler& X)
{
   	ASSERT(calls->set_initial_value != NULL);
   	calls->set_initial_value(this, X);
}

void
LoadableElem::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
   	ASSERT(calls->set_value != NULL);
   	calls->set_value(this, pDM, X, XP, ph);
}

unsigned int
LoadableElem::iGetNumPrivData(void) const
{
   	ASSERT(calls->i_get_num_priv_data != NULL);
   	return calls->i_get_num_priv_data(this);
}

unsigned int
LoadableElem::iGetPrivDataIdx(const char *s) const
{
   	ASSERT(calls->i_get_priv_data_idx != NULL);
   	return calls->i_get_priv_data_idx(this, s);
}

doublereal
LoadableElem::dGetPrivData(unsigned int i) const
{
   	ASSERT(calls->d_get_priv_data != NULL);
   	return calls->d_get_priv_data(this, i);
}

int
LoadableElem::GetNumConnectedNodes(void) const
{
	ASSERT(calls->i_get_num_connected_nodes != NULL);
	return calls->i_get_num_connected_nodes(this);
}

void
LoadableElem::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
	ASSERT(calls->get_connected_nodes != NULL);
	return calls->get_connected_nodes(this, connectedNodes);
}

UserDefinedElem*
LoadableElemRead::Read(unsigned int uLabel, const DofOwner* pDO,
	DataManager* pDM, MBDynParser& HP) const
{
	UserDefinedElem* pEl = NULL;

	if (HP.IsKeyWord("reference")) {
		const char *s = HP.GetStringWithDelims();

		const LoadableCalls *c = pDM->GetLoadableElemModule(s);

		if (c == 0) {
			silent_cerr("Loadable(" << uLabel << "): "
				"unable to find loadable element module \""
				<< s << "\" at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		SAFENEWWITHCONSTRUCTOR(pEl,
				LoadableElem,
				LoadableElem(uLabel, pDO, c, pDM, HP));

	} else {
		SAFENEWWITHCONSTRUCTOR(pEl,
				LoadableElem,
				LoadableElem(uLabel, pDO, pDM, HP));
	}

   	return pEl;
}
