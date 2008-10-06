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
/* module-aerodyn
 * Authors: Fanzhong Meng, Pierangelo Masarati
 *
 * Copyright (C) 2008
 *
 * Fanzhong Meng <f.meng@tudelft.nl>
 * 
 * Faculty of Aerospace Engineering - Delft University of Technology
 * Kluyverweg 1, 2629HS Delft, the Netherlands
 * http://www.tudelft.nl
 *
 */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h> 		/* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "Rot.hh"

#include "loadable.h"
#include "drive_.h"
#include "module-aerodyn.h"

static module_aerodyn_t *module_aerodyn;

// Called by AeroDyn

/*
 * Rotor parameters - called once per time step.
 */
/* This comment is added by Fanzhong MENG 9th.Feb.2008
 *
 * I think in these three functions we should add the codes
 * that can deal with the communication with MBDyn to get
 * the structure information which is need for the
 * Aerodynamic calculation from libAeroDyn.a. 
 */


int
__FC_DECL__(getrotorparams)(
	F_REAL *Omega,
	F_REAL *gamma,
	F_REAL *VHUB,
	F_REAL *tau)
{
    	/*
	 * This comment is added by Fanzhong MENG 10th.Feb.2008
	 * NOTE: Add code here to get rotorspeed from MBDyn.  
	 * Rotor can be running clockwise or ccw.  But the
	 * value of *Omega must be positive [rad/s]
	 */
	/*
	 * PM 2008-09-06:
	 * nacelle & rotor axis is node axis 3, positive opposite
	 * to wind direction
	 */

	const Mat3x3& nacelle_R = ::module_aerodyn->pNacelle->GetRCurr();
	const Vec3& hub_Omega = ::module_aerodyn->pHub->GetWCurr();
	Vec3 rotation_axis = nacelle_R.GetVec(3);

	*Omega = std::fabs(hub_Omega*rotation_axis);

    	/* 
	 * This comment is added by Fanzhong MENG 10th.Feb.2008
	 * NOTE: Add code here to get Yawangle from MBDyn.  
	 * Yaw angle of the shaft relative to the ground reference.
	 * Positive value is clockwise when viewed from above of nacelle.
	 * Variable name is *gamma. [rad].
	 */
	*gamma = std::atan2(-rotation_axis(2), -rotation_axis(1));

    	/*
	 * This comment is added by Fanzhong MENG 10th.Feb.2008
	 * NOTE: Add code here to get tilt angle from MBDyn.
	 * Tilt angle of the shaft relative to the ground reference.
	 * Positive value is Tilt UP
	 * Variable name is *tau. [rad].
	 */
	*tau = std::atan2(rotation_axis(3),
		std::sqrt(rotation_axis(1)*rotation_axis(1) + rotation_axis(2)*rotation_axis(2)));

    	/*
	 * This comment is added by Fanzhong MENG 10th.Feb.2008
	 * NOTE: Add code here to get Hub Tangential velocity from MBDyn.
	 * Positive value is in the same direction with Yaw Angle.
	 * Variable name is *VHUB.[m/s or feet/s]
	 */ 
	*VHUB = hub_Omega(3)*::module_aerodyn->Hub_Tower_xy_distance;

#ifdef MODULE_AERODYN_DEBUG
	silent_cerr("getrotorparams: "
		"Omega=" << *Omega << " "
		"gamma=" << (*gamma)*180/M_PI << " "
		"tau=" << (*tau)*180/M_PI << " "
		"VHUB=" << *VHUB << std::endl);
#endif // MODULE_AERODYN_DEBUG

	return 0;
}

/*
 * Blade parameters - called once for each blade at each time step.
 */
int
__FC_DECL__(getbladeparams)(F_REAL *psi)
{
    	/*
	 * This comment is added by Fanzhong MENG 10th.Feb.2008
	 * NOTE: Add code here too get blade Azimuth angle.
	 * The 6 o'clock position is positive value
	 * Variable name is *psi. [rad].
	 */
	const Mat3x3& R_nacelle = ::module_aerodyn->pNacelle->GetRCurr();
	const Mat3x3& R_hub = ::module_aerodyn->pHub->GetRCurr();
	const Mat3x3& R_blade_root = ::module_aerodyn->bladeR[::module_aerodyn->c_blade - 1];
	Mat3x3 R = (R_hub*R_blade_root).MulTM(R_nacelle);
	Vec3 Phi(RotManip::VecRot(R));
	*psi = -Phi(3);

	// unwrap psi (0 <= psi <= 360)
	while (*psi < 0) {
	    *psi += 2*M_PI;
	}

#ifdef MODULE_AERODYN_DEBUG
	silent_cerr("getbladeparams: blade=" << ::module_aerodyn->c_blade
		<< " psi=" << (*psi)*180/M_PI << std::endl);
#endif // MODULE_AERODYN_DEBUG

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
    	/*
	 * This comment is added by Fanzhong MENG 10th.Feb.2008
	 * Get coordinates of blade element relative to the undeflected
	 * tower top reference frame. At the Hub height.
	 * The variable name is *XGRND *YGRND *ZGRND
	 */
    	/*
	 * Get blade elements coordinates X,Y,Z.
	 */ 
	Vec3 x = ::module_aerodyn->nodes[::module_aerodyn->elem].pNode->GetXCurr()
		- ::module_aerodyn->pNacelle->GetXCurr(); // Maybe we have to modified here. Here we need
						   	  // hub height coordinates.
	*XGRND = x(1);
	*YGRND = x(2);
	*ZGRND = x(3);
	/*
	 * Get blade elements radius. It is the perpendicular distance 
	 * from the rotor axis to the element aerodynamic reference point.
	 */ 
	const Vec3& X_node = ::module_aerodyn->nodes[::module_aerodyn->elem].pNode->GetXCurr();
	const Vec3& X_hub = ::module_aerodyn->pHub->GetXCurr();
	const Mat3x3& R_nacelle = ::module_aerodyn->pNacelle->GetRCurr();
	Vec3 d = R_nacelle.MulTV(X_node - X_hub);
	d(3) = 0.;
	*radius = d.Norm();

	/*
	 * Get the blade elements local pitch angle. Here X-axis is in the blade
	 * spanwise direction.
	 */ 
	unsigned blade = ::module_aerodyn->elem/::module_aerodyn->nelems;
	const Mat3x3& R_node = ::module_aerodyn->nodes[::module_aerodyn->elem].pNode->GetRCurr(); //Pitch bug is not fixed!
	const Mat3x3& R_hub = ::module_aerodyn->pHub->GetRCurr();
	const Mat3x3& R_blade_root = ::module_aerodyn->bladeR[blade];
	Mat3x3 R = (R_hub*R_blade_root).MulTM(R_node);
	Vec3 Phi(RotManip::VecRot(R.Transpose()));

	// Pitch Now = Pitch + Twist 
	*phi = Phi(1) + ::module_aerodyn->nodes[::module_aerodyn->elem].dBuiltInTwist;
	::module_aerodyn->nodes[::module_aerodyn->elem].PITNOW = *phi;

    	/* 
	 * This comment is added by Fanzhong MENG 10th.Feb.2008
	 * When the Multi-airfoil table option is open, we should also 
	 * switch on this option. By default, this value should be set to ZERO.
	 * I don't know how this Multi-airfoil table work,yet.
	 */ 
	*MulTabLoc = 0;

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
    	/*
	 * Fanzhong MENG 18th.Feb.2008.
	 * This function returns velocities of the wind
	 * and the element in the ground reference frame.
	 * VX,VY,VZ are provided by AeroDyn for MBDyn.
	 */
        /*
	 * This "getvnvt" function is correct now!
         */
	/*
	 * velocity of the element related to Ground reference frame.
	 */
	const Vec3& ve_G = ::module_aerodyn->nodes[::module_aerodyn->elem].pNode->GetVCurr();

	/*
	 * velocity of the wind related to Ground reference frame.
	 * The "-" sign change the wind velocity from Aerodyn reference frame to MBDyn 
	 *  Global reference frame. AeroDyn global reference frame can be found on page 35 of Aerodyn user's guide--Figure C1
	 */
	Vec3 vw_G = Vec3(-(*VX), -(*VY), *VZ);

	/*
	 * Transform wind vector from ground coordinate system (vw_G) to
	 * blade Node reference frame (vw_El)
	 */
	const Mat3x3& R_node = ::module_aerodyn->nodes[::module_aerodyn->elem].pNode->GetRCurr();
	/*
	 * 10th of June. Modified by fanzhong meng.
	 * get orientation matrix of node twist with respect to node without pitch.
	 */
	const Mat3x3& R_node_twist = ::module_aerodyn->nodes[::module_aerodyn->elem].Ra; 

	
	Vec3 vw_El = R_node_twist.MulTV( R_node.MulTV(vw_G) );
	Vec3 ve_El = R_node_twist.MulTV( R_node.MulTV(ve_G) );

	/*
	 * Calculate SIN and COS value of current pitch angle.
	 */

	doublereal PitchNow = ::module_aerodyn->nodes[::module_aerodyn->elem].PITNOW; 
	doublereal CPitchNow = std::cos(PitchNow);
	doublereal SPitchNow = std::sin(PitchNow);
	
	/*
	 * Play around with these SIN and COS value to get
	 * *VT,*VNW and *VNE.
	 */

	doublereal vt_wind = -vw_El(2)*CPitchNow - vw_El(3)*SPitchNow;
	doublereal vt_element = ve_El(2)*CPitchNow + ve_El(3)*SPitchNow;

	*VT = vt_wind + vt_element;
	*VNW = -vw_El(2)*SPitchNow + vw_El(3)*CPitchNow;
	/* 
	 * Element Normal velocity opposites to the Wind Normal velocity. 
	 * You can find is on page 35 of Aerodyn user's guide--Figure C2.
	 */
	*VNE = -(-ve_El(2)*SPitchNow + ve_El(3)*CPitchNow);

	return 0;
}

/*
 * Write an error message to the appropriate stream
 * By Fanzhong Meng: usrmes is an Adams built in function.
 * We have to rewrite this function. Because libAeroDyn.a library-GenSub.f90 will call this function,
 * so we can not simply remove this function.
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
#if 0
	// can't work, unless we know the size of the arrays!
	silent_cerr("module-aerodyn: msg=\"" << msg << "\" "
		"code=" << *code
		<< " level=" << level
		<< std::endl);
#endif
	silent_cerr("module-aerodyn: diagnostics from AeroDyn, "
		"code=" << *code << std::endl);

	return 0;
}

/* default funcs */
static void *
module_aerodyn_read(
		LoadableElem *pEl,
		DataManager* pDM,
		MBDynParser& HP
)
{
	DEBUGCOUTFNAME("module_aerodyn_read");

	if (::module_aerodyn != NULL) {
		silent_cerr("Another module-aerodyn might be running; error"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
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
		<< "Authors: Fanzhong Meng, Pierangelo Masarati" << std::endl
		<< std::endl
		<< "This is the MBDyn interface to AeroDyn, the aerodynamic routines" << std::endl
		<< "developed by NREL <http://www.nrel.gov/> to model the aerodynamic" << std::endl
		<< "forces acting on wind turbines" << std::endl
		<< std::endl
		<< "usage:" << std::endl
		<< "Nacelle node; requirements:" << std::endl
		<< "\t\t- axis 3 is the shaft axis" << std::endl
		<< "\t\t- axis 3 in wind direction" << std::endl
		<< "\t<nacelle node label> ," << std::endl
		<< "\t<hub node label> ," << std::endl
		<< "\t<pylon top-hub xy distance> ," << std::endl
		<< "\t<hub radius> ," << std::endl
		<< "\t<number of blades> ," << std::endl
		<< "\t<number of elements per blade> ," << std::endl
		<< "\t# for each blade... axis 1 must be the blade spanwise direction" << std::endl
		<< "\t\t<i-th blade root orientation matrix> ," << std::endl
		<< "\t\t# for each blade element..." << std::endl
		<< "\t\t\t<i-th blade j-th node label> ," << std::endl
		<< "\t\t\t[ orientation , <i-th blade j-th node orientation> , ]" << std::endl
		<< "\t[ output file name , \" <file name> \" ]" << std::endl
		);
	}

	/* read nacelle node */
	p->pNacelle = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
	if (p->pNacelle == 0) {
		silent_cerr("Aerodyn(" << pEl->GetLabel() << "): "
			"nacelle node not defined "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* read hub node */
	p->pHub = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
	if (p->pHub == 0) {
		silent_cerr("Aerodyn(" << pEl->GetLabel() << "): "
			"hub node not defined "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	p->Hub_Tower_xy_distance = HP.GetReal();

	p->r_hub = HP.GetReal();

	/*
	 * Initialize AeroDyn package
	 *
	 * FIXME: does it return an error code?
	 */
	char Version[26 + 1];
	snprintf(Version, sizeof(Version), "(%s)", VERSION);
	for (unsigned i = strlen(Version); i < sizeof(Version); i++) {
		Version[i] = ' ';
	}
	Version[sizeof(Version) - 1] = '\0';

	// number of blades
	F_INTEGER NBlades = HP.GetInt();
	if (NBlades <= 0) {
		silent_cerr("Aerodyn(" << pEl->GetLabel() << "): "
			"invalid number of blades " << NBlades
			<< " at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// number of elements per blade
	F_INTEGER NElems = HP.GetInt();
	if (NElems <= 0) {
		silent_cerr("Aerodyn(" << pEl->GetLabel() << "): "
			"invalid number of elements per blade " << NElems
			<< " at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	p->nblades = NBlades;
	p->nelems = NElems;

	// get blade root orientation(without pitch angle).

    	p->nodes.resize(NBlades*NElems);
    	p->bladeR.resize(NBlades);
	ReferenceFrame rf(p->pHub);
	for (unsigned e = 0; e < unsigned(NBlades*NElems); e++) {
		if ((e % NElems) == 0) {
		    	/*
			 * Get the orientation matrix of blade root with respect to hub reference frame
			 */ 
			p->bladeR[e/NElems] = HP.GetRotRel(rf);
#if 0
		std::cerr << "Root[" << e/NElems << "] Rotation Matrix:" << p->bladeR[e/NElems] << std::endl;
		std::cerr << "Hub Rotation Matrix:" << p->pHub->GetRCurr() << std::endl;
		std::cerr << "Nacelle Rotation Matrix:" << p->pNacelle->GetRCurr() << std::endl<< std::endl;
#endif
		}

		p->nodes[e].pNode = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));

		// (optional) aerodynamics offset with respect to the node,
		// constant in the reference frame of the node
		if (HP.IsKeyWord("position")) {
			p->nodes[e].f = HP.GetPosRel(ReferenceFrame(p->nodes[e].pNode));

		} else {
			p->nodes[e].f = Zero3;
		}

		// (optional) relative orientation between the aerodynamics
		// and the node
		if (HP.IsKeyWord("orientation")) {
			p->nodes[e].Ra = HP.GetRotRel(ReferenceFrame(p->nodes[e].pNode));
			// built-in twist: the component about blade axis 1
			// (x) of the relative orientation between
			// the aerodynamics and the node
			Vec3 Twist(RotManip::VecRot(p->nodes[e].Ra));
			p->nodes[e].dBuiltInTwist = -Twist(1);

		} else {
			p->nodes[e].Ra = Eye3;
			p->nodes[e].dBuiltInTwist = 0.;
		}
	}

	if (HP.IsKeyWord("output" "file" "name")) {
		const char *ofname = HP.GetFileName();
		if (ofname == 0) {
			silent_cerr("Aerodyn(" << pEl->GetLabel() << "): "
				"unable to get file name "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		p->ofname = ofname;
		p->out.open(ofname);
		if (!p->out) {
			silent_cerr("Aerodyn(" << pEl->GetLabel() << "): "
				"unable to open file \"" << ofname << "\" "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	__FC_DECL__(mbdyn_init)(Version, &NBlades);

	int rc;
	if (HP.IsKeyWord("input" "file" "name")) {
		const char *input_file_name = HP.GetStringWithDelims();
		if (input_file_name == 0) {
			silent_cerr("unable to get input file name "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		F_INTEGER input_file_name_len = strlen(input_file_name);
		rc = __FC_DECL__(mbdyn_ad_inputgate)((F_CHAR *)input_file_name, &input_file_name_len);

	} else {
		rc = __FC_DECL__(adinputgate)();
	}

	if (rc != 0) {
		silent_cerr("Aerodyn(" << pEl->GetLabel() << "): "
			"initialization failed "
			"(err=" << rc << ")" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	(void)__FC_DECL__(mbdyn_true)(&p->FirstLoop);

	// Calculate the Tip and Hub loss constants for AeroDyn.
	for (int e = 0; e < p->nelems; e++) {
	    p->elem = e;
	    const Vec3& X_node = p->nodes[p->elem].pNode->GetXCurr();
	    const Vec3& X_hub = p->pHub->GetXCurr();
	    const Mat3x3& R_nacelle = p->pNacelle->GetRCurr();
	    Vec3 d = R_nacelle.MulTV(X_node - X_hub);
//	    d(3) = 0.;
	    p->rlocal = (F_REAL)d.Norm();
	    p->c_elem = e + 1;
	    __FC_DECL__(mbdyn_get_tl_const)(&p->rlocal, &p->c_elem);
    	    __FC_DECL__(mbdyn_get_hl_const)(&p->rlocal, &p->c_elem, &p->r_hub);
	}


	// FIXME: only needed when Beddoes and also Dynamic Inflow model are enabled 
	p->Time.Set(new TimeDriveCaller(pDM->pGetDrvHdl()));
	p->dCurTime = p->Time.dGet();
	(void)__FC_DECL__(mbdyn_sim_time)(&p->dCurTime);

	::module_aerodyn = p;

	return (void *)p;
}

static unsigned int
module_aerodyn_i_get_num_dof(const LoadableElem* pEl)
{
	DEBUGCOUTFNAME("module_aerodyn_i_get_num_dof");
	return 0;
}

static DofOrder::Order
module_aerodyn_set_dof(const LoadableElem*, unsigned int i)
{
	DEBUGCOUTFNAME("module_aerodyn_set_dof");
	return DofOrder::UNKNOWN;
}

static void
module_aerodyn_output(const LoadableElem* pEl, OutputHandler& OH)
{
	DEBUGCOUTFNAME("module_aerodyn_output");

	module_aerodyn_t* p = (module_aerodyn_t *)pEl->pGetData(); 
	
	if (p->out) {
		unsigned c = 0;
		for (unsigned b = 0; b < unsigned(p->nblades); b++) {
			for (unsigned e = 0; e < unsigned(p->nelems); e++, c++) {
				p->out
					<< b + 1 << '.' << e + 1
					<< " " << p->nodes[c].F
					<< " " << p->nodes[c].M
					<< std::endl;
			}
		}
	}
}

static std::ostream&
module_aerodyn_restart(const LoadableElem* pEl, std::ostream& out)
{
	DEBUGCOUTFNAME("module_aerodyn_restart");
	return out << "not implemented yet;" << std::endl;
}

static void
module_aerodyn_work_space_dim(const LoadableElem* pEl, integer* piNumRows, integer* piNumCols)
{
	DEBUGCOUTFNAME("module_aerodyn_work_space_dim");

	module_aerodyn_t* p = (module_aerodyn_t *)pEl->pGetData(); 
	
	*piNumRows = 6*p->nblades*p->nelems;
	*piNumCols = 1;
}

static SubVectorHandler& 
module_aerodyn_ass_res(
		LoadableElem* pEl, 
		SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr
)
{
	DEBUGCOUTFNAME("module_aerodyn_ass_res");
	integer iNumRows = 0;
	integer iNumCols = 0;
	pEl->WorkSpaceDim(&iNumRows, &iNumCols);
	
	WorkVec.ResizeReset(iNumRows);

	module_aerodyn_t* p = (module_aerodyn_t *)pEl->pGetData(); 
	
	/*
	 * set sub-vector indices and coefs
	 */

	if (p->bFirst) {
		int c_blade = -1;
		p->c_blade = 0;
		Mat3x3 BladeR;
		Vec3 BladeAxis;

		for (int e = 0; e < p->nelems*p->nblades; e++) {
			/*
			 * get the current blade number.
			 */ 
			if ((e % p->nelems) == 0) {
				c_blade++;
				p->c_blade++;
				BladeR = p->pHub->GetRCurr()*p->bladeR[c_blade];
				BladeAxis = BladeR.GetVec(1);
			}

			/*
			 * get force/moment
			 */
			F_REAL		DFN, DFT, PMA;
	
			/*
			 * current element; use as index to access array
			 * when inside callbacks
			 */
			p->elem = e;

			p->c_elem = e%p->nelems + 1;
			(void)__FC_DECL__(mbdyn_com_data)(&p->c_blade, &p->c_elem);

                        /*
	                 * Get blade elements radius. It is the perpendicular distance 
	                 * from the rotor axis to the element aerodynamic reference point.
	                 */ 

			/*
			 * forces and moments are in the blade reference frame,
			 * not in the element's.
			 */
			__FC_DECL__(aerofrcintrface)(&p->FirstLoop, &p->c_elem, &DFN, &DFT, &PMA);

#ifdef MODULE_AERODYN_DEBUG
		        silent_cerr("aerodyn[" << e << "] in rotation plane: DFN=" << DFN << " DFT=" << DFT << " PMA=" << PMA << std::endl);
#endif // MODULE_AERODYN_DEBUG

			/*
			 * turn force/moment into the global frame,
			 * the moment with respect to the node
			 */
			p->nodes[e].F = BladeR*Vec3(0., DFT, DFN);
			Vec3 Arm(p->nodes[e].pNode->GetRCurr()*p->nodes[e].f);
			p->nodes[e].M = BladeAxis*PMA + Arm.Cross(p->nodes[e].F);

#if 0 // def MODULE_AERODYN_DEBUG
	        	silent_cerr("aerodyn[" << e << "] in BEM frame: "
				<< F_BEM << " " << PMA << " " << 0.0 << " " << 0.0 << std::endl
				<< "aerodyn[" << e << "] in global frame: "
				<< p->nodes[e].F << " " << p->nodes[e].M
				<< std::endl << std::endl);
#endif // MODULE_AERODYN_DEBUG
		}

#ifdef MODULE_AERODYN_DEBUG
 	       silent_cerr(std::endl);
#endif // MODULE_AERODYN_DEBUG

		p->bFirst = false;
	}

	for (unsigned e = 0; e < unsigned(p->nblades*p->nelems); e++) {
		/*
		 * set indices where force/moment need to be put
		 */
		integer iFirstIndex = p->nodes[e].pNode->iGetFirstMomentumIndex();
		for (int i = 1; i <= 6; i++) {
			WorkVec.PutRowIndex(6*e + i, iFirstIndex + i);
		}

		/*
		 * add force/moment to residual, after rotating them
		 * into the global frame
		 */
		/*
		 * The force/moment is already in the global reference frame
	         */	 
		WorkVec.Add(6*e + 1, p->nodes[e].F);
		WorkVec.Add(6*e + 4, p->nodes[e].M);

	}

	/* 
	 * make sure next time FirstLoop will be false
	 */
	(void)__FC_DECL__(mbdyn_false)(&p->FirstLoop);
	
	return WorkVec;
}

#if 0
static void
module_aerodyn_before_predict(const LoadableElem* pEl, 
	VectorHandler& X, VectorHandler& XP,
	VectorHandler& XPrev, VectorHandler& XPPrev)
{
	DEBUGCOUTFNAME("module_aerodyn_before_predict");
}
#endif

static void
module_aerodyn_after_predict(const LoadableElem* pEl, 
	VectorHandler& X, VectorHandler& XP)
{
    	/* NOTE: think about what should be done here!*/
	DEBUGCOUTFNAME("module_aerodyn_after_predict");

	module_aerodyn_t* p = (module_aerodyn_t *)pEl->pGetData(); 

	/*
	 * Get the current simulation time and time step!
	 * MENG: 19 June 2008.
	 */ 
	p->dDT = p->Time.dGet() - p->dOldTime;
	p->dCurTime = p->Time.dGet();
	(void)__FC_DECL__(mbdyn_sim_time)(&p->dCurTime);
	(void)__FC_DECL__(mbdyn_time_step)(&p->dDT);
}

#if 0
static void
module_aerodyn_update(LoadableElem* pEl, 
	const VectorHandler& X, const VectorHandler& XP)
{
	DEBUGCOUTFNAME("module_aerodyn_update");
}
#endif

static void 
module_aerodyn_after_convergence(const LoadableElem* pEl,
	const VectorHandler& /* X */ , const VectorHandler& /* XP */ )
{
    	/* NOTE: think about what should be done here!*/
	DEBUGCOUTFNAME("module_aerodyn_after_convergence");

	module_aerodyn_t* p = (module_aerodyn_t *)pEl->pGetData(); 

	p->bFirst = true;
	p->dOldTime = p->Time.dGet();
}

static unsigned int
module_aerodyn_i_get_initial_num_dof(const LoadableElem* pEl)
{
	DEBUGCOUTFNAME("module_aerodyn_i_get_initial_num_dof");
	return 0;
}

static void
module_aerodyn_set_value(const LoadableElem* pEl, DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
	DEBUGCOUTFNAME("module_aerodyn_set_value");

	module_aerodyn_t* p = (module_aerodyn_t *)pEl->pGetData(); 

	p->bFirst = true;
}

#if 0   
static void
module_aerodyn_set_initial_value(const LoadableElem* pEl, VectorHandler& X)
{
	DEBUGCOUTFNAME("module_aerodyn_set_initial_value");
}
#endif

static unsigned int
module_aerodyn_i_get_num_priv_data(const LoadableElem* pEl)
{
	DEBUGCOUTFNAME("module_aerodyn_i_get_num_priv_data");
	return 0;
}

static void
module_aerodyn_destroy(LoadableElem* pEl)
{
	DEBUGCOUTFNAME("module_aerodyn_destroy");

	module_aerodyn_t* p = (module_aerodyn_t *)pEl->pGetData();

	/*
	 * delete private data
	 */
	::module_aerodyn = NULL;
	
	SAFEDELETE(p);
}

static int
module_aerodyn_i_get_num_connected_nodes(const LoadableElem* pEl)
{
	DEBUGCOUTFNAME("module_aerodyn_i_get_num_connected_nodes");
	return 0;
}

static void
module_aerodyn_get_connected_nodes(const LoadableElem* pEl, 
		std::vector<const Node *>& connectedNodes)
{
	DEBUGCOUTFNAME("module_aerodyn_get_connected_nodes");

#if 0
	module_aerodyn_t* p = (module_aerodyn_t *)pEl->pGetData();
#endif // 0

	/*
	 * set args according to element connections
	 */
	connectedNodes.resize(module_aerodyn_i_get_num_connected_nodes(pEl));
}

static struct
LoadableCalls lc = {
	LOADABLE_VERSION_SET(1, 5, 0),

	"AeroDyn",
	"1.1",
	"Dipartimento di Ingegneria Aerospaziale, Politecnico di Milano",
	"wind turbine aerodynamics; uses AeroDyn package 12.58\n"
		"\tas distributed by the National Renewable Energy Laboratory,\n"
		"\thttp://www.nrel.gov",

	module_aerodyn_read,
	module_aerodyn_i_get_num_dof,
	module_aerodyn_set_dof,
	module_aerodyn_output,
	module_aerodyn_restart,
	module_aerodyn_work_space_dim,
	NULL /* module_aerodyn_ass_jac */ ,
	NULL /* module_aerodyn_ass_mats */ ,
	module_aerodyn_ass_res,
	NULL /* module_aerodyn_before_predict */ ,
	module_aerodyn_after_predict,
	NULL /* module_aerodyn_update */ ,
	module_aerodyn_after_convergence,
	module_aerodyn_i_get_initial_num_dof,
	NULL /* module_aerodyn_initial_work_space_dim */ ,
	NULL /* module_aerodyn_initial_ass_jac */ ,
	NULL /* module_aerodyn_initial_ass_res */ ,
	module_aerodyn_set_value,
	NULL /* module_aerodyn_set_initial_value */ ,
	module_aerodyn_i_get_num_priv_data,
	NULL /* module_aerodyn_i_get_priv_data_idx */ ,
	NULL /* module_aerodyn_d_get_priv_data */ ,
	module_aerodyn_i_get_num_connected_nodes,
	module_aerodyn_get_connected_nodes,
	module_aerodyn_destroy,
	NULL,
	NULL,
	NULL,
	NULL
};

extern "C" {
void *calls = &lc;
}

