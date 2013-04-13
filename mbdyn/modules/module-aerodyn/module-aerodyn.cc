/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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
 * Copyright (C) 2008-2013
 *
 * Fanzhong Meng <f.meng@tudelft.nl>
 * 
 * Faculty of Aerospace Engineering - Delft University of Technology
 * Kluyverweg 1, 2629HS Delft, the Netherlands
 * http://www.tudelft.nl
 *
 */
/*
 * AeroDynModule: interface between MBDyn and NREL's AeroDyn.
 *
 * Author: Pierangelo Masarati
 *
 * Copyright (C) 2008-2013 All rights reserved
 *
 * Refactoring of module-aerodyn using the UserDefinedElem interface.
 * Based on module-aerodyn.cc by Fanzhong Meng and Pierangelo Masarati.
 */

#include "mbconfig.h" 		/* This goes first in every *.c,*.cc file */

#include "Rot.hh"

#include "dataman.h"
#include "userelem.h"
//#include "loadable.h"
#include "drive_.h"
//#include "module-aerodyn.h"

// uncomment to enable debug output
// #define MODULE_AERODYN_DEBUG

// keep this consistent with AeroDyn build
// #define USE_DOUBLE_PRECISION
#define USE_SINGLE_PRECISION

#include "NREL_AeroDyn.h"

class AeroDynModule
: virtual public Elem,
public UserDefinedElem
{
private:

	struct AeroNode {
		StructNode	*pNode;
		Vec3		f;	// offset of the aero point wrt./ the node,
					// constant in the node's reference frame
		Mat3x3		Ra;	// aerodynamic orientation of the aero point
					// wrt./ the node

		doublereal	dBuiltInTwist;

		Vec3		F;	// Force acting on Node
		Vec3		M;	// Moment acting on Node, with respect
					// to node's position

 		doublereal      FN;    // Normal Force on each blade element (Not being used now!).
 		doublereal      FT;    // Tangental force on each blade element.(Not being used now!)
 		doublereal      AM;    // Aerodynamic moment on each blade element.(Not being used now!)
 
		doublereal	PITNOW; // Node pitch angle
	};

	/**
	 * Nacelle node; requirements:
	 * - axis 3 is the shaft axis
	 * - axis 3 in wind direction
	 */
	StructNode	*pNacelle;
	StructNode	*pHub;

	integer		nblades; // the number of blades.
	integer		nelems;  // the number of elements per blade.

	doublereal	Hub_Tower_xy_distance;

	/*
	 * node data
	 */
	std::vector<AeroNode>	nodes;		// nodes
	std::vector<Mat3x3>	bladeR;		// orientation matrix of each blade root in the hub reference frame

	std::string ofname;
	mutable std::ofstream out;

	/*
	 * Total aerodynamic data
	 */
	Vec3		TF;     // Total Force on the rotor in the absolute frame
	Vec3		TM;     // Total Moment on the rotor in the absolute frame.
	Vec3		TF_h;     // Total Force on the rotor in the hub frame
	Vec3		TM_h;     // Total Moment on the rotor in the hub frame.
	doublereal      Thrust;   // Rotor Thrust.
	doublereal      Torque;   // Rotor Torque.
	doublereal      Rotor_speed;   // Rotor angular velocity.
	/* 
	 * internal states to access the Variables which is defined in the 
	 * common MODULEs of AeroDyn
	 */
	F_LOGICAL	FirstLoop;
	F_INTEGER	elem;    // use to identify the current element in the interface module!
	F_INTEGER	c_elem;  // use to identify the current element in AeroDyn!
	F_INTEGER	c_blade; // use to identify the current blade in AeroDyn!
	F_REAL		rlocal;  // use to identify the current element position 
	F_REAL		r_hub;   // use to identify the hub radius.
	F_REAL		r_rotor; // use to identify the rotor radius.

	bool bFirst;
	DriveOwner	Time;		// time drive
	doublereal	dOldTime;	// old time
	doublereal      dCurTime;       // current time
	F_REAL          dDT;		// time step

	DriveOwner FSF;

public:
	// interface for AeroDyn callbacks

	// module_aerodyn->pNacelle
	const StructNode *pGetNacelleNode(void) const;
	// module_aerodyn->pHub
	const StructNode *pGetHubNode(void) const;
	// module_aerodyn->Rotor_speed
	void SetRotorSpeed(doublereal Omega);
	// module_aerodyn->Hub_Tower_xy_distance;
	doublereal dGetHubTowerXYDistance(void) const;
	// module_aerodyn->nblades
	F_INTEGER iGetNumBlades(void) const;
	// module_aerodyn->c_blade
	F_INTEGER iGetCurrBlade(void) const;
	// module_aerodyn->bladeR
	const Mat3x3& GetCurrBladeR(void) const;
	// module_aerodyn->nelems
	F_INTEGER iGetNumBladeElems(void) const;
	// module_aerodyn->elem
	F_INTEGER iGetCurrBladeElem(void) const;
	// module_aerodyn->nodes
	const StructNode *pGetCurrBladeNode(void) const;
	// module_aerodyn->nodes[::module_aerodyn->elem].dBuiltInTwist
	doublereal dGetCurrBladeNodeBuiltinTwist(void) const;
	// module_aerodyn->nodes[::module_aerodyn->elem].PITNOW
	void SetCurrBladeNodePITNOW(doublereal PITNOW);
	doublereal dGetCurrBladeNodePITNOW(void) const;
	// module_aerodyn->nodes[::module_aerodyn->elem].Ra
        const Mat3x3& GetCurrBladeNodeRa(void) const;

public:
	AeroDynModule(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	~AeroDynModule(void);

	unsigned int iGetNumDof(void) const;
	DofOrder::Order GetDofType(unsigned int i) const;
	void Output(OutputHandler& OH) const;
	std::ostream& Restart(std::ostream& out) const;
	void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
	VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef, 
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr);
	void AfterPredict(VectorHandler& X, VectorHandler& XP);
	void AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP);
	unsigned int iGetInitialNumDof(void) const;
	void InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat, 
		const VectorHandler& XCurr);
   	SubVectorHandler& 
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
	void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints* h = 0);
	unsigned int iGetNumPrivData(void) const;
#if 0
	unsigned int iGetPrivDataIdx(const char *s) const;
	doublereal dGetPrivData(unsigned int i) const;
#endif
	int GetNumConnectedNodes(void) const;
	void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
};

// static handler for AeroDyn callbacks
// only one element per simulation can be active (AeroDyn limitation)
static AeroDynModule *module_aerodyn;

AeroDynModule::AeroDynModule(
	unsigned uLabel,
	const DofOwner *pDO,
	DataManager* pDM,
	MBDynParser& HP)
: Elem(uLabel, 0),
UserDefinedElem(uLabel, pDO),
bFirst(true)
{
	if (HP.IsKeyWord("help")) {
		/* NOTE: add help message */
		silent_cout(
		"Module: AeroDyn" << std::endl
		<< std::endl
		<< "Author: Fanzhong Meng, Pierangelo Masarati" << std::endl
		<< std::endl
		<< "This is the MBDyn interface to AeroDyn, the aerodynamic routines" << std::endl
		<< "developed by NREL <http://www.nrel.gov/> to model the aerodynamic" << std::endl
		<< "forces acting on wind turbines" << std::endl
		<< std::endl
		<< "usage:" << std::endl
		<< "#Nacelle node; requirements:" << std::endl
		<< "#\t\t- axis 3 is the shaft axis, pointing downstream the wind" << std::endl
		<< "\t<nacelle node label> ," << std::endl
		<< "\t<hub node label> ," << std::endl
		<< "\t<pylon top-hub xy distance> ," << std::endl
		<< "\t<hub radius> ," << std::endl
		<< "\t<rotor radius> ," << std::endl
		<< "\t<number of blades> ," << std::endl
		<< "\t<number of elements per blade> ," << std::endl
		<< "\t# for each blade..." << std::endl
		<< "\t#\t- axis 1 in the blade spanwise direction, towards blade tip" << std::endl
		<< "\t#\t- axis 2 in coord direction, towards leading edge" << std::endl
		<< "\t#\t- axis 3 in thickness direction, towards low pressure side" << std::endl
		<< "\t\t<i-th blade root orientation matrix> ," << std::endl
		<< "\t\t# for each blade element..." << std::endl
		<< "\t\t\t<i-th blade j-th node label>" << std::endl
		<< "\t\t\t[ , orientation , <i-th blade j-th node orientation> ]" << std::endl
		<< "\t[ , force scale factor , (DriveCaller)<factor> ]" << std::endl
		<< "\t[ , output file name , \"<file name>\" ]" << std::endl
		<< "\t[ , AeroDyn input file name , \"<file name>\" ]" << std::endl
		<< "\t[ , element file name , \"<file name>\" ]" << std::endl
		);
	}

	if (::module_aerodyn != 0) {
		silent_cerr("AeroDynModule::AeroDynModule(" << uLabel << "): "
			"AeroDyn already in use" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* read nacelle node */
	pNacelle = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
	if (pNacelle == 0) {
		silent_cerr("AeroDynModule(" << GetLabel() << "): "
			"nacelle node not defined "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* read hub node */
	pHub = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
	if (pHub == 0) {
		silent_cerr("AeroDynModule(" << GetLabel() << "): "
			"hub node not defined "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	Hub_Tower_xy_distance = HP.GetReal();

	r_hub = HP.GetReal();
	
	r_rotor = HP.GetReal();

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
		silent_cerr("AeroDynModule(" << GetLabel() << "): "
			"invalid number of blades " << NBlades
			<< " at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// number of elements per blade
	F_INTEGER NElems = HP.GetInt();
	if (NElems <= 0) {
		silent_cerr("AeroDynModule(" << GetLabel() << "): "
			"invalid number of elements per blade " << NElems
			<< " at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	nblades = NBlades;
	nelems = NElems;

	// get blade root orientation(without pitch angle).

    	nodes.resize(NBlades*NElems);
    	bladeR.resize(NBlades);
	ReferenceFrame rf(pHub);
	for (elem = 0; elem < NBlades*NElems; elem++) {
		if ((elem % NElems) == 0) {
		    	/*
			 * Get the orientation matrix of blade root with respect to hub reference frame
			 */ 
			bladeR[elem/NElems] = HP.GetRotRel(rf);
#if 0
		std::cerr << "Root[" << elem/NElems << "] Rotation Matrix:" << bladeR[elem/NElems] << std::endl;
		std::cerr << "Hub Rotation Matrix:" << pHub->GetRCurr() << std::endl;
		std::cerr << "Nacelle Rotation Matrix:" << pNacelle->GetRCurr() << std::endl<< std::endl;
#endif
		}

		nodes[elem].pNode = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));

		// (optional) aerodynamics offset with respect to the node,
		// constant in the reference frame of the node
		if (HP.IsKeyWord("position")) {
			nodes[elem].f = HP.GetPosRel(ReferenceFrame(nodes[elem].pNode));

		} else {
			nodes[elem].f = Zero3;
		}

		// (optional) relative orientation between the aerodynamics
		// and the node
		if (HP.IsKeyWord("orientation")) {
			nodes[elem].Ra = HP.GetRotRel(ReferenceFrame(nodes[elem].pNode));
			// built-in twist: the component about blade axis 1
			// (x) of the relative orientation between
			// the aerodynamics and the node
			Vec3 Twist(RotManip::VecRot(nodes[elem].Ra));
			nodes[elem].dBuiltInTwist = -Twist(1);

		} else {
			nodes[elem].Ra = Eye3;
			nodes[elem].dBuiltInTwist = 0.;
		}
	}

	if (HP.IsKeyWord("force" "scale" "factor")) {
		FSF.Set(HP.GetDriveCaller());

	} else {
		FSF.Set(new OneDriveCaller);
	}

	if (HP.IsKeyWord("output" "file" "name")) {
		const char *ofname = HP.GetFileName();
		if (ofname == 0) {
			silent_cerr("AeroDynModule(" << GetLabel() << "): "
				"unable to get file name "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		ofname = ofname;
		out.open(ofname);
		if (!out) {
			silent_cerr("AeroDynModule(" << GetLabel() << "): "
				"unable to open file \"" << ofname << "\" "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		} else {
			out << std::setw(16) << "Time s"
			    << std::setw(16) << "Thrust kN"
			    << std::setw(16) << "Torque kNm"
			    << std::setw(16) << "R.Speed Rad/s"
			    << std::setw(16) << "Power kW"
			    << std::endl;
		}
	}

	__FC_DECL__(mbdyn_init)(Version, &NBlades, &r_rotor);

	std::string input_file_name;
	std::string elem_file_name;

	if (HP.IsKeyWord("input" "file" "name")) {
		const char *fname = HP.GetStringWithDelims();
		if (fname == 0) {
			silent_cerr("unable to get input file name "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		input_file_name = fname;
	}

	if (HP.IsKeyWord("element" "file" "name")) {
		const char *fname = HP.GetStringWithDelims();
		if (fname == 0) {
			silent_cerr("unable to get element file name "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		elem_file_name = fname;
	}

	F_INTEGER input_file_name_len = input_file_name.size();
	F_INTEGER elem_file_name_len = elem_file_name.size();
	int rc = __FC_DECL__(mbdyn_ad_inputgate)((F_CHAR *)input_file_name.c_str(), &input_file_name_len,
			(F_CHAR *)elem_file_name.c_str(), &elem_file_name_len);
	if (rc != 0) {
		silent_cerr("AeroDynModule(" << GetLabel() << "): "
			"initialization failed "
			"(err=" << rc << ")" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	(void)__FC_DECL__(mbdyn_true)(&FirstLoop);

	// Calculate the Tip and Hub loss constants for AeroDyn.
	for (elem = 0; elem < nelems; elem++) {
		const Vec3& X_node = nodes[elem].pNode->GetXCurr();
		const Vec3& X_hub = pHub->GetXCurr();
		const Mat3x3& R_nacelle = pNacelle->GetRCurr();
		Vec3 d = R_nacelle.MulTV(X_node - X_hub);
		rlocal = (F_REAL)d.Norm();
		c_elem = elem + 1;
		__FC_DECL__(mbdyn_get_tl_const)(&rlocal, &c_elem);
		__FC_DECL__(mbdyn_get_hl_const)(&rlocal, &c_elem, &r_hub);
	}

	// FIXME: only needed when Beddoes and also Dynamic Inflow model are enabled 
	Time.Set(new TimeDriveCaller(pDM->pGetDrvHdl()));
	dCurTime = Time.dGet();
	(void)__FC_DECL__(mbdyn_sim_time)(&dCurTime);

	::module_aerodyn = this;
}

AeroDynModule::~AeroDynModule(void)
{
	NO_OP;
}

unsigned int
AeroDynModule::iGetNumDof(void) const
{
	return 0;
}

DofOrder::Order
AeroDynModule::GetDofType(unsigned int i) const
{
	return DofOrder::UNKNOWN;
}

void
AeroDynModule::Output(OutputHandler& OH) const
{
	if (out) {
		out << std::scientific
		    << std::setw(16) << dCurTime 
		    << std::setw(16) << Thrust/1000.0
		    << std::setw(16) << Torque/1000.0
		    << std::setw(16) << Rotor_speed
		    << std::setw(16) << Torque*Rotor_speed/1000.0
		    << std::endl; 
	}
}

std::ostream&
AeroDynModule::Restart(std::ostream& out) const
{
	return out << "not implemented yet;" << std::endl;
}

void
AeroDynModule::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 6*nblades*nelems;
	*piNumCols = 1;
}

VariableSubMatrixHandler& 
AeroDynModule::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	// should do something useful
	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
AeroDynModule::AssRes(
	SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	
	WorkVec.ResizeReset(iNumRows);

	/*
	 * set sub-vector indices and coefs
	 */

	if (bFirst) {
		c_blade = 0;
		Mat3x3 BladeR;
		Vec3 BladeAxis;

		// force scale factor to "ramp up" loads
		doublereal dFSF = FSF.dGet();

		// sanity check
		ASSERT(dFSF >= 0.);
		ASSERT(dFSF <= 1.);

		for (elem = 0; elem < nelems*nblades; elem++) {
			/*
			 * get the current blade number.
			 */ 
			if ((elem % nelems) == 0) {
				c_blade++;
				BladeR = pHub->GetRCurr()*bladeR[c_blade];
				BladeAxis = BladeR.GetVec(1);
			}

			/*
			 * get force/moment
			 */
			F_REAL		DFN, DFT, PMA;
	
			c_elem = elem%nelems + 1;
			(void)__FC_DECL__(mbdyn_com_data)(&c_blade, &c_elem);

                        /*
	                 * Get blade elements radius. It is the perpendicular distance 
	                 * from the rotor axis to the element aerodynamic reference point.
	                 */ 

			/*
			 * forces and moments are in the blade reference frame,
			 * not in the element's.
			 */
			__FC_DECL__(aerofrcintrface)(&FirstLoop, &c_elem, &DFN, &DFT, &PMA);
			
			// ramp forces in the very first second up to ease convergence
			DFN *= dFSF;
			DFT *= dFSF;
			PMA *= dFSF;

 			nodes[elem].FN = DFN;
 			nodes[elem].FT = DFT;
 			nodes[elem].AM = PMA;

#ifdef MODULE_AERODYN_DEBUG
		        silent_cerr("aerodyn[" << elem << "] in rotation plane: DFN=" << DFN << " DFT=" << DFT << " PMA=" << PMA << std::endl);
#endif // MODULE_AERODYN_DEBUG
			/*
			 * turn force/moment into the node frame (blade element frame used by Aerodyn to calculation forces)
			 * (passing thru local element frame),
			 * including offset
			 */

			doublereal c = std::cos(nodes[elem].PITNOW);
			doublereal s = std::sin(nodes[elem].PITNOW);
			nodes[elem].F = Vec3(0., DFT*c - DFN*s, DFT*s + DFN*c);
			nodes[elem].M = nodes[elem].Ra.GetVec(1)*PMA + nodes[elem].f.Cross(nodes[elem].F);


#ifdef MODULE_AERODYN_DEBUG

 	        silent_cerr("Moment on Node[" << elem << "]: " << nodes[elem].f.Cross(nodes[elem].F) << std::endl);

#endif // MODULE_AERODYN_DEBUG

#ifdef MODULE_AERODYN_DEBUG
			silent_cerr(std::setprecision(3) << std::fixed
			          << " aerodyn[" << elem << "]"
			          << " F = " << nodes[elem].F 
			          << " M = " << nodes[elem].M
			          << " FN = " << nodes[elem].FN
			          << " FN = " << nodes[elem].FT
			          << " AM = " << nodes[elem].AM
			          << std::endl << std::endl);
#endif // MODULE_AERODYN_DEBUG

		}

#ifdef MODULE_AERODYN_DEBUG
 	       silent_cerr(std::endl);
#endif // MODULE_AERODYN_DEBUG

		bFirst = false;
	}


	TF = Vec3(Zero3);
	TM = Vec3(Zero3);
	for (elem = 0; elem < nblades*nelems; elem++) {
		/*
		 * set indices where force/moment need to be put
		 */
		integer iFirstIndex = nodes[elem].pNode->iGetFirstMomentumIndex();
		for (int i = 1; i <= 6; i++) {
			WorkVec.PutRowIndex(6*elem + i, iFirstIndex + i);
		}

		/*
		 * add force/moment to residual, after rotating them
		 * into the global frame
		 */
		/*
		 * Transfer the force/moment to global frame. 
	         */	 
  		WorkVec.Add(6*elem + 1, nodes[elem].pNode->GetRCurr()*(nodes[elem].Ra*nodes[elem].F));
  		WorkVec.Add(6*elem + 4, nodes[elem].pNode->GetRCurr()*(nodes[elem].Ra*nodes[elem].M));

		/*
		 * calculate the force and moment contributions on the rotor in absolute frame.
	         */
		const Vec3& Xp = nodes[elem].pNode->GetXCurr();
		const Vec3& Xh = pHub->GetXCurr();

		TF = TF + nodes[elem].pNode->GetRCurr()*nodes[elem].Ra*nodes[elem].F;
		TM = TM + nodes[elem].pNode->GetRCurr()*nodes[elem].Ra*nodes[elem].M
			+ (Xp-Xh).Cross(nodes[elem].pNode->GetRCurr()*nodes[elem].Ra*nodes[elem].F);
  	}

	/*
	 * Transfer the force/moment to hub reference frame.
	 */
	const Mat3x3& Rh = pHub->GetRCurr();
	TF_h = Rh.MulTV(TF);
	TM_h = Rh.MulTV(TM);

	/*
	 * The 3rd component of TF_h and TM_h are the rotor Thrust and rotor Torque.
	 */
	Thrust = TF_h(3);
	Torque = TM_h(3);
  
	/* 
	 * make sure next time FirstLoop will be false
	 */
	if (FirstLoop) {
		(void)__FC_DECL__(mbdyn_false)(&FirstLoop);
	}
	
	return WorkVec;
}

#if 0
void
module_aerodyn_before_predict(const LoadableElem* pEl, 
	VectorHandler& X, VectorHandler& XP,
	VectorHandler& XPrev, VectorHandler& XPPrev)
{
	DEBUGCOUTFNAME("module_aerodyn_before_predict");
}
#endif

void
AeroDynModule::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	bFirst = true;

	/*
	 * Get the current simulation time and time step!
	 * MENG: 19 June 2008.
	 */ 
	dDT = Time.dGet() - dOldTime;
	dCurTime = Time.dGet();
	(void)__FC_DECL__(mbdyn_sim_time)(&dCurTime);
	(void)__FC_DECL__(mbdyn_time_step)(&dDT);
}

#if 0
void
module_aerodyn_update(LoadableElem* pEl, 
	const VectorHandler& X, const VectorHandler& XP)
{
	DEBUGCOUTFNAME("module_aerodyn_update");
}
#endif

void 
AeroDynModule::AfterConvergence(
	const VectorHandler& X,
	const VectorHandler& XP)
{
	dOldTime = Time.dGet();
}

unsigned int
AeroDynModule::iGetInitialNumDof(void) const
{
	return 0;
}

void
AeroDynModule::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 6*nblades*nelems;
	*piNumCols = 1;
}

VariableSubMatrixHandler&
AeroDynModule::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
AeroDynModule::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}

void
AeroDynModule::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	bFirst = true;
}

unsigned int
AeroDynModule::iGetNumPrivData(void) const
{
	return 0;
}

int
AeroDynModule::GetNumConnectedNodes() const
{
	return nodes.size() + 2;
}

void
AeroDynModule::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	/*
	 * set args according to element connections
	 */
	connectedNodes.resize(nodes.size() + 2);

	unsigned n;
	for (n = 0; n < nodes.size(); n++) {
		connectedNodes[n] = nodes[n].pNode;
	}

	connectedNodes[n++] = pNacelle;
	connectedNodes[n++] = pHub;
}

// module_aerodyn->pNacelle
const StructNode *
AeroDynModule::pGetNacelleNode(void) const
{
	return pNacelle;
}

// module_aerodyn->pHub
const StructNode *
AeroDynModule::pGetHubNode(void) const
{
	return pHub;
}

// module_aerodyn->Rotor_speed
void
AeroDynModule::SetRotorSpeed(doublereal Omega)
{
	Rotor_speed = Omega;
}

// module_aerodyn->Hub_Tower_xy_distance;
doublereal
AeroDynModule::dGetHubTowerXYDistance(void) const
{
	return Hub_Tower_xy_distance;
}

// module_aerodyn->nblades
F_INTEGER
AeroDynModule::iGetNumBlades(void) const
{
	return nblades;
}

// module_aerodyn->c_blade
F_INTEGER
AeroDynModule::iGetCurrBlade(void) const
{
	ASSERT(c_blade >= 0);
	ASSERT(c_blade < nblades);

	return c_blade;
}

// module_aerodyn->bladeR
const Mat3x3&
AeroDynModule::GetCurrBladeR(void) const
{
	ASSERT(c_blade >= 0);
	ASSERT(c_blade < nblades);

	return bladeR[c_blade - 1];
}

// module_aerodyn->nelems
F_INTEGER
AeroDynModule::iGetNumBladeElems(void) const
{
	return nelems;
}

// module_aerodyn->elem
F_INTEGER
AeroDynModule::iGetCurrBladeElem(void) const
{
	ASSERT(elem >= 0);
	ASSERT(elem < nblades*nelems);

	return elem;
}

// module_aerodyn->nodes
const StructNode *
AeroDynModule::pGetCurrBladeNode(void) const
{
	ASSERT(elem >= 0);
	ASSERT(elem < nblades*nelems);

	return nodes[elem].pNode;
}

// module_aerodyn->nodes[::module_aerodyn->elem].dBuiltInTwist
doublereal
AeroDynModule::dGetCurrBladeNodeBuiltinTwist(void) const
{
	ASSERT(elem >= 0);
	ASSERT(elem < nblades*nelems);

	return nodes[elem].dBuiltInTwist;
}

// module_aerodyn->nodes[::module_aerodyn->elem].PITNOW
void
AeroDynModule::SetCurrBladeNodePITNOW(doublereal PITNOW)
{
	ASSERT(elem >= 0);
	ASSERT(elem < nblades*nelems);

	nodes[elem].PITNOW = PITNOW;
}

doublereal
AeroDynModule::dGetCurrBladeNodePITNOW(void) const
{
	ASSERT(elem >= 0);
	ASSERT(elem < nblades*nelems);

	return nodes[elem].PITNOW;
}

// module_aerodyn->nodes[::module_aerodyn->elem].Ra
const Mat3x3&
AeroDynModule::GetCurrBladeNodeRa(void) const
{
	ASSERT(elem >= 0);
	ASSERT(elem < nblades*nelems);

	return nodes[elem].Ra;
}

// Functions called by AeroDyn

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

	const Mat3x3& nacelle_R = ::module_aerodyn->pGetNacelleNode()->GetRCurr();
	const Vec3& hub_Omega = ::module_aerodyn->pGetHubNode()->GetWCurr();
	Vec3 rotation_axis = nacelle_R.GetVec(3);

	*Omega = std::abs(hub_Omega*rotation_axis);
	module_aerodyn->SetRotorSpeed(*Omega);

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
	// *VHUB = hub_Omega(3)*::module_aerodyn->dGetHubTowerXYDistance();
	/* See Emails by Fanzhong Meng August 27, 2010 */
	const Vec3& nacelle_Omega = ::module_aerodyn->pGetNacelleNode()->GetWCurr();
	*VHUB = nacelle_Omega(3)*::module_aerodyn->dGetHubTowerXYDistance();
	
	__FC_DECL__(elemout)();

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
	const Mat3x3& R_nacelle = ::module_aerodyn->pGetNacelleNode()->GetRCurr();
	const Mat3x3& R_hub = ::module_aerodyn->pGetHubNode()->GetRCurr();
	const Mat3x3& R_blade_root = ::module_aerodyn->GetCurrBladeR();
	Mat3x3 R = (R_hub*R_blade_root).MulTM(R_nacelle);
	Vec3 Phi(RotManip::VecRot(R));
	*psi = -Phi(3);

	// unwrap psi (0 <= psi <= 360)
	while (*psi < 0) {
	    *psi += 2*M_PI;
	}

#ifdef MODULE_AERODYN_DEBUG
	silent_cerr("getbladeparams: blade=" << ::module_aerodyn->iGetCurrBlade()
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
	Vec3 x = ::module_aerodyn->pGetCurrBladeNode()->GetXCurr()
		- ::module_aerodyn->pGetNacelleNode()->GetXCurr();
	// Maybe we have to modified here. Here we need
	// hub height coordinates.
	*XGRND = x(1);
	*YGRND = x(2);
	*ZGRND = x(3);
	/*
	 * Get blade elements radius. It is the perpendicular distance 
	 * from the rotor axis to the element aerodynamic reference point.
	 */ 
	const Vec3& X_node = ::module_aerodyn->pGetCurrBladeNode()->GetXCurr();
	const Vec3& X_hub = ::module_aerodyn->pGetHubNode()->GetXCurr();
	const Mat3x3& R_nacelle = ::module_aerodyn->pGetNacelleNode()->GetRCurr();
	Vec3 d = R_nacelle.MulTV(X_node - X_hub);
	d(3) = 0.;
	*radius = d.Norm();

	/*
	 * Get the blade elements local pitch angle. Here X-axis is in the blade
	 * spanwise direction.
	 */ 
	// unsigned blade = ::module_aerodyn->iGetCurrBlade();
	const Mat3x3& R_node = ::module_aerodyn->pGetCurrBladeNode()->GetRCurr(); 
	const Mat3x3& R_hub = ::module_aerodyn->pGetHubNode()->GetRCurr();
	const Mat3x3& R_blade_root = ::module_aerodyn->GetCurrBladeR();
	Mat3x3 R = (R_hub*R_blade_root).MulTM(R_node);
	Vec3 Phi(RotManip::VecRot(R.Transpose()));

	// Pitch Now = Pitch + Twist 
	*phi = Phi(1) + ::module_aerodyn->dGetCurrBladeNodeBuiltinTwist();
	::module_aerodyn->SetCurrBladeNodePITNOW(*phi);

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
	const Vec3& ve_G = ::module_aerodyn->pGetCurrBladeNode()->GetVCurr();

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
	const Mat3x3& R_node = ::module_aerodyn->pGetCurrBladeNode()->GetRCurr();
	/*
	 * 10th of June. Modified by fanzhong meng.
	 * get orientation matrix of node twist with respect to node without pitch.
	 */
	const Mat3x3& R_node_twist = ::module_aerodyn->GetCurrBladeNodeRa(); 

	
	Vec3 vw_El = R_node_twist.MulTV( R_node.MulTV(vw_G) );
	Vec3 ve_El = R_node_twist.MulTV( R_node.MulTV(ve_G) );

	/*
	 * Calculate SIN and COS value of current pitch angle.
	 */

	doublereal PitchNow = ::module_aerodyn->dGetCurrBladeNodePITNOW(); 
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

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	UserDefinedElemRead *rf = new UDERead<AeroDynModule>;

	if (!SetUDE("aerodyn", rf)) {
		delete rf;

		silent_cerr("module-aerodyn: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

