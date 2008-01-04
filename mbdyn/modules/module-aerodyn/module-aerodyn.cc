/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2008
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

// Called by AeroDyn

/*
 * Rotor parameters - called once per time step.
 */
int
__FC_DECL__(getrotorparams)(
	F_REAL *Omega,
	F_REAL *gamma,
	F_REAL *VHUB,
	F_REAL *tau)
{
	return 0;
}

/*
 * Blade parameters - called once for each blade at each time step.
 */
int
__FC_DECL__(getbladeparams)(F_REAL *psi)
{
	return 0;
}

/*
 * Element parameters - called once for each element at each time step.
 */
int
__FC_DECL__(getelemparams)(
	F_INTEGER *MulTabLoc,
	F_REAL *phi,
	F_REAL *radius,
	F_REAL *XGRND,
	F_REAL *YGRND,
	F_REAL *ZGRND)
{
	return 0;
}

/*
 * Compute VT, VN{W,E} based on VX, VY, VZ of the wind.
 */
int
__FC_DECL__(getvnvt)(
	F_REAL *VX,
	F_REAL *VY,
	F_REAL *VZ,
	F_REAL *VT,
	F_REAL *VNW,
	F_REAL *VNE)
{
	return 0;
}

/*
 * Write an error message to the appropriate stream
 *
 * FIXME: the "msg" and "level" arrays should be reset by the caller
 * before writing the message, otherwise they're not '\0' terminated
 */
int
__FC_DECL__(usrmes)(
	F_LOGICAL *Logical,
	F_CHAR msg[],
	F_INTEGER *code,
	F_CHAR level[])
{
	return 0;
}

module_aerodyn_t *module_aerodyn = NULL;

/* default funcs */
static void *
read(
		LoadableElem *pEl,
		DataManager* pDM,
		MBDynParser& HP
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
		silent_cout(
"Module: AeroDyn" << std::endl
<< std::endl
<< "This is the MBDyn interface to AeroDyn, the aerodynamic routines" << std::endl
<< "developed by NREL <http://ww.nrel.gov/> to model the aerodynamic" << std::endl
<< "forces acting on wind turbines" << std::endl
);
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
	char Version[26];
	snprintf(Version, sizeof(Version), "(%s)", VERSION);
	for (unsigned i = strlen(Version); i < sizeof(Version); i++) {
		Version[i] = ' ';
	}

	// number of blades
	F_INTEGER NBlades = HP.GetInt();
	if (NBlades <= 0) {
		silent_cerr("Aerodyn(" << pEl->GetLabel() << "): "
			"invalid number of blades " << NBlades
			<< " at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric();
	}

	__FC_DECL__(mbdyn_init)(Version, &NBlades);

#if 0
	// does not work as expected...
	const char *tmp = HP.GetStringWithDelims();
	if (tmp == 0) {
		silent_cerr("unable to get input file "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric();
	}

	char input_file[80];
	if (strlen(tmp) >= sizeof(input_file)) {
		silent_cerr("Aerodyn(" << pEl->GetLabel() << "): "
			"file name \"" << tmp << "\" too long "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric();
	}

	strncpy(input_file, tmp, sizeof(input_file));
	for (unsigned i = strlen(input_file); i < sizeof(input_file); i++) {
		input_file[i] = ' ';
	}
	int rc = __FC_DECL__(ad_inputgate)(input_file);
#else
	// the input file name must be "aerodyn.ipt"
	int rc = __FC_DECL__(adinputgate)();
#endif
	if (rc != 0) {
		silent_cerr("Aerodyn(" << pEl->GetLabel() << "): "
			"initialization failed "
			"(err=" << rc << ")" << std::endl);
		throw ErrGeneric();
	}

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
		std::vector<const Node *>& connectedNodes)
{
	DEBUGCOUTFNAME("get_connected_nodes");

#if 0
	module_aerodyn_t* p = (module_aerodyn_t *)pEl->pGetData();
#endif /* 0 */

	/*
	 * set args according to element connections
	 */
	connectedNodes.resize(i_get_num_connected_nodes(pEl));
}

static struct
LoadableCalls lc = {
	LOADABLE_VERSION_SET(1, 5, 0),

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

