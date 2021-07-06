/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <fstream>
#include <iomanip>
#include <cfloat>
#include <vector>
#include "elem.h"
#include "strnode.h"
#include "dataman.h"
#include "userelem.h"
#include "module-uvlm.h"
#include "mbdyn_uvlm.h"


UvlmInterfaceBaseElem::UvlmInterfaceBaseElem(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO),
m_PDM(pDM),
MBDyn_UVLM_Model_Converged(pDM),
bMBDyn_UVLM_FirstSend(true),
bMBDyn_UVLM_Model_DoStepDynamics(true)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
			"									\n"
			"Module: 	module-chrono-interface					\n"
			"Organization:	Dipartimento di Scienze e Tecnologie Aerospaziali	\n"
			"		Politecnico di Milano					\n"
			"		http://www.aero.polimi.it/				\n"
			"Author: Shubhaditya Burela, Andrea Zanoni, and Pierangelo Masarati	\n"
			"Supported: Google Summer of Code 2021						\n"
			"Description: This module dynamically links the MBDyn solver to the UVLM C++ library(https://github.com/ImperialCollegeLondon/UVLM)\n"
			"Usage: \n"
			"<elem_type> ::= user defined \n"
			"<user_defined_type> ::=UvlmInterface \n"
			"<ChronoInterface_arglist> ::= <cosimulation_platform>, \n"
			"                              <coupling_variables>, \n"
			"                              <coupling_bodies>, \n"
			"                              <other_settings>; \n"
			"1. <cosimulation_platform> ::= coupling, {none | <loose_coupling> | <tight_coupling> }, \n"
			"                               <loose_coupling> ::= loose, [{embedded | jacobian | gauss}], \n"
			"                               <tight_coupling> ::= tight, (int) <num_iterations>, {tolerance, (real) <tolerance>}, \n"
			"2. <coupling_variables>    ::= [force type, {reaction | contact}], \n"
			"                               [length scale, (real) <length_scale>], [mass scale, (real) <math_scale>], \n"
			"3. <coupling_bodies>       ::= nodes numeber, (int) <num_coupling_nodes>, \n"
			"                               <mbdyn_nodes_info>, <chrono_bodies_info>, [<coupling_constraint>], [<chbody_output>] \n"
			"                               ...\n"
			"                               <chrono_ground>\n"
			"                               <mbdyn_nodes_info>     ::= mbdyn node, (int) <node_label>,\n"
			"                                                          [offset, (Vec3) <offset_coupling_point>, rotation, (Mat3x3) <relative_orientation>], \n"
			"                               <chorno_bodies_info>   ::= chrono body, (int) <chbody_label>, \n"
			"                                                          [offset, (Vec3) <offset_coupling_point>], \n"
			"                               <coupling_constraints> ::= [position constraint, (int) <bool_x>, (int) <bool_y>, (int) <bool_z>], \n"
			"                                                          [rotation constraint, (int) <bool_x>, (int) <bool_y>, (int) <bool_z>], \n"
			"                               <chbody_output>        ::= [output chbody, {yes | no}], \n"
			"                               <chrono_ground>        ::= ground, (int) <chbody_label>, \n"
			"                                                          [position, (Vec3) <absolute position>], \n"
			"                                                          [orientation, (Mat3x3) <absolute_rotation>], [<chbody_output>],\n"
			"4. <other_settings>        ::= [output all chbodies], [coupling forces filename, (filename) <filename>], [verbose]. \n"
			"-------------------------------------------------start simulation-------------------------------------------------\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	// Read information from the script
	

	//--------------- 1. <cosimulation_platform> ------------------------
	MBDyn_UVLM_CouplingIter_Count = 0;
	MBDyn_UVLM_Coupling_Tol = 1.0e-3; // by default, tolerance == 1.0e-6;
	MBDyn_UVLM_CouplingType = MBDyn_UVLM_COUPLING::COUPLING_NONE;
	if (HP.IsKeyWord("coupling")) {
		if (HP.IsKeyWord("none"))
		{
			MBDyn_UVLM_CouplingType = MBDyn_UVLM_COUPLING::COUPLING_NONE;
			MBDyn_UVLM_CouplingIter_Max = 0;
		}
		//--------------- <loose_coupling> ---------------
		else if (HP.IsKeyWord("loose"))
		{
			MBDyn_UVLM_CouplingType = MBDyn_UVLM_COUPLING::COUPLING_LOOSE;
			// By default, using embedded scheme.
			MBDyn_UVLM_CouplingType_loose = MBDyn_UVLM_COUPLING_LOOSE::LOOSE_EMBEDDED;
			MBDyn_UVLM_CouplingIter_Max = 1;
			std::cout << "Embedded-loose method is chosen" << std::endl;
			if (HP.IsKeyWord("embedded"))
			{
				MBDyn_UVLM_CouplingType_loose = MBDyn_UVLM_COUPLING_LOOSE::LOOSE_EMBEDDED;
			}
			else if (HP.IsKeyWord("jacobian"))
			{
				MBDyn_UVLM_CouplingType_loose = MBDyn_UVLM_COUPLING_LOOSE::LOOSE_JACOBIAN;
			}
			else if (HP.IsKeyWord("gauss"))
			{
				MBDyn_UVLM_CouplingType_loose = MBDyn_UVLM_COUPLING_LOOSE::LOOSE_GAUSS;
			}
			else
			{
				silent_cerr("The type of loose coupling scheme at line " << HP.GetLineData() << " is unknown \n");
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
		//--------------- <tight_coupling> ---------------
		else if (HP.IsKeyWord("tight"))
		{
			MBDyn_UVLM_CouplingType = MBDyn_UVLM_COUPLING::COUPLING_TIGHT;
			MBDyn_UVLM_CouplingType_loose = MBDyn_UVLM_COUPLING_LOOSE::TIGHT; // meaningless. 
			MBDyn_UVLM_CouplingIter_Max = HP.GetInt();
			if (HP.IsKeyWord("tolerance"))
			{
				MBDyn_UVLM_Coupling_Tol = HP.GetReal();
			}
			else
			{
				MBDyn_UVLM_Coupling_Tol = 1.0e-3; //pDM -> GetSolver()->pGetStepIntegrator()->GetIntegratorDTol();
			}
		}
		else
		{
			silent_cerr("Unknown coupling type at line " << HP.GetLineData() << "\n");
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
	else
	{
		silent_cerr("Unknown keyword at line " << HP.GetLineData() << "\n"); //- the first keywork must be "coupling"
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	//--------------- 2. <coupling_variables> ------------------------
	// read the force type: transformed contact forces or reaction forces acting on the UVLM coupling bodies
	MBDyn_UVLM_ForceType = MBDyn_UVLM_FORCETYPE::REACTION_FORCE;           //- by default, using the reaction force 
	if (HP.IsKeyWord("force" "type"))
	{
		if (HP.IsKeyWord("reaction"))
		{
			MBDyn_UVLM_ForceType = MBDyn_UVLM_FORCETYPE::REACTION_FORCE; //- reaction force is imposed to objects in MBDyn
		}
		else if (HP.IsKeyWord("contact"))
		{
			MBDyn_UVLM_ForceType = MBDyn_UVLM_FORCETYPE::CONTACT_FORCE; //- contact force is imposed to objects in MBDyn
		}
		else
		{
			silent_cerr("Force type at line " << HP.GetLineData() << " is not implemeted yet \n");
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
	//--------------- read the scale type-----------------------------	
	MBDyn_UVLM_Scale[0] = 1.0; //- default value for length scale
	MBDyn_UVLM_Scale[1] = 1.0; //- default value for mass scale
	MBDyn_UVLM_Scale[2] = 1.0; //- default value for force scaleforce scale
	MBDyn_UVLM_Scale[3] = 1.0; //- default value for torque scalemoment scale
	if (HP.IsKeyWord("length" "scale"))
	{
		MBDyn_UVLM_Scale[0] = HP.GetReal();		   //- length scale
		MBDyn_UVLM_Scale[1] = MBDyn_UVLM_Scale[0]; //- if not mentioned, mass scale = length scale;
	}
	if (HP.IsKeyWord("mass" "scale"))
	{
		MBDyn_UVLM_Scale[1] = HP.GetReal();          //- mass scale
	}
	MBDyn_UVLM_Scale[2] = MBDyn_UVLM_Scale[1] * MBDyn_UVLM_Scale[0]; //- force scale
	MBDyn_UVLM_Scale[3] = MBDyn_UVLM_Scale[2] * MBDyn_UVLM_Scale[0]; //- moment scale


}

void
UvlmInterfaceBaseElem::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler& 
UvlmInterfaceBaseElem::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr){

	pedantic_cout("\tMBDyn::AssJac()" << std::endl);
	WorkMat.SetNullMatrix();

	return WorkMat;
}


unsigned int
UvlmInterfaceBaseElem::iGetNumPrivData(void) const
{
	return 0;
}

int
UvlmInterfaceBaseElem::iGetNumConnectedNodes(void) const
{
	return 0;
}

void
UvlmInterfaceBaseElem::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(0);
}


std::ostream&
UvlmInterfaceBaseElem::Restart(std::ostream& out) const
{
	return out << "# UvlmInterface: is doing now" << std::endl;
}

unsigned int
UvlmInterfaceBaseElem::iGetInitialNumDof(void) const
{
	return 0;
}

void 
UvlmInterfaceBaseElem::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
UvlmInterfaceBaseElem::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr){

	pedantic_cout("\tInitialAssJac()" << std::endl);
	WorkMat.SetNullMatrix();
	switch (MBDyn_UVLM_CouplingType)
	{
	case MBDyn_UVLM_COUPLING::COUPLING_NONE:
		break; //- do nothing
	case MBDyn_UVLM_COUPLING::COUPLING_TIGHT: //- do nothing
		break;
	case MBDyn_UVLM_COUPLING::COUPLING_LOOSE: //- do nothing
	case MBDyn_UVLM_COUPLING::COUPLING_STSTAGGERED: //- do nothing
	default:
		break;
	}
	return WorkMat;
}

SubVectorHandler& 
UvlmInterfaceBaseElem::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr){

	pedantic_cout("\tInitialAssRes()" << std::endl);
	WorkVec.ResizeReset(0);
	switch (MBDyn_UVLM_CouplingType)
	{
	case MBDyn_UVLM_COUPLING::COUPLING_NONE:
		break; //- do nothing
	case MBDyn_UVLM_COUPLING::COUPLING_TIGHT: //- do nothing
		break;
	case MBDyn_UVLM_COUPLING::COUPLING_LOOSE:	   //- do nothing
	case MBDyn_UVLM_COUPLING::COUPLING_STSTAGGERED: //- do nothing
	default:
		break;
	}
	return WorkVec;
}

void 
UvlmInterfaceBaseElem::MBDyn_UVLM_SendDataToBuf_Prev() {

	/* formulation for calculation
			mbdynuvlm_x = x + mbdynuvlm_f
			mbdynuvlm_R = R
			mbdynuvlm_v = xp + mbdynuvlm_w cross mbdynuvlm_f
			mbdynuvlm_w = w
			mbdynuvlm_a = xpp + mbdynuvlm_wp cross mbdynuvlm_f + mbdynuvlm_w cross mbdynuvlm_w cross mbdynuvlm_f
			mbdynuvlm_wp = wp
	*/

	for (unsigned i = 0; i < MBDyn_UVLM_NodesNum; i++) {
		const MBDYN_UVLM_POINTDATA& mbdynuvlm_point = MBDyn_UVLM_Nodes[i];
		// Rotation and position
		const Mat3x3& mbdynuvlm_R = mbdynuvlm_point.pMBDyn_UVLM_Node->GetRPrev();
		const Mat3x3& mbdynuvlm_Rh = MBDyn_UVLM_Nodes[i].MBDyn_UVLM_RhM;   // Relative orientation of the marker in local ref.
		Mat3x3 mbdynuvlm_R_marker = mbdynuvlm_R * mbdynuvlm_Rh;
		Vec3 mbdynuvlm_f = mbdynuvlm_R * mbdynuvlm_point.MBDyn_UVLM_Offset;
		Vec3 mbdynuvlm_x = mbdynuvlm_point.pMBDyn_UVLM_Node->GetXPrev() + mbdynuvlm_f;
		// Angular velocity and velocity
		const Vec3& mbdynuvlm_w = mbdynuvlm_point.pMBDyn_UVLM_Node->GetWPrev();
		Vec3 mbdynuvlm_wCrossf = mbdynuvlm_w.Cross(mbdynuvlm_f);
		Vec3 mbdynuvlm_v = mbdynuvlm_point.pMBDyn_UVLM_Node->GetVPrev() + mbdynuvlm_wCrossf;
		// Angular accelaration and accelaration
		const Vec3 &mbdynuvlm_wp = mbdynuvlm_point.pMBDyn_UVLM_Node->GetWPPrev();
		Vec3 mbdynuvlm_a = mbdynuvlm_point.pMBDyn_UVLM_Node->GetXPPPrev() + mbdynuvlm_wp.Cross(mbdynuvlm_f) + 
			mbdynuvlm_w.Cross(mbdynuvlm_wCrossf);

		double mbdynuvlm_tempvec3_x[3];
		double mbdynuvlm_tempvec3_v[3];
		double mbdynuvlm_tempvec3_a[3];
		double mbdynuvlm_tempvec3_w[3];
		double mbdynuvlm_tempvec3_wp[3];
		MBDyn_UVLM_Vec3D(mbdynuvlm_x, mbdynuvlm_tempvec3_x, MBDyn_CE_CEScale[0]);
		MBDyn_UVLM_Vec3D(mbdynuvlm_v, mbdynuvlm_tempvec3_v, MBDyn_CE_CEScale[0]);
		MBDyn_UVLM_Vec3D(mbdynuvlm_a, mbdynuvlm_tempvec3_a, MBDyn_CE_CEScale[0]);
		MBDyn_UVLM_Vec3D(mbdynuvlm_w, mbdynuvlm_tempvec3_w, 1.0);
		MBDyn_UVLM_Vec3D(mbdynuvlm_wp, mbdynuvlm_tempvec3_wp, 1.0);
		double mbdynuvlm_tempmat3x3_R[9];
		MBDyn_UVLM_Mat3x3D(mbdynuvlm_R_marker, mbdynuvlm_tempmat3x3_R);


		memcpy(&pMBDyn_UVLM_CouplingKinematic_x[3 * i], mbdynuvlm_tempvec3_x, 3 * sizeof(double));
		memcpy(&pMBDyn_UVLM_CouplingKinematic_R[9 * i], mbdynuvlm_tempmat3x3_R, 9 * sizeof(double));
		memcpy(&pMBDyn_UVLM_CouplingKinematic_xp[3 * i], mbdynuvlm_tempvec3_v, 3 * sizeof(double));
		memcpy(&pMBDyn_UVLM_CouplingKinematic_omega[3 * i], mbdynuvlm_tempvec3_w, 3 * sizeof(double));
		memcpy(&pMBDyn_UVLM_CouplingKinematic_xpp[3 * i], mbdynuvlm_tempvec3_a, 3 * sizeof(double));
		memcpy(&pMBDyn_UVLM_CouplingKinematic_omegap[3 * i], mbdynuvlm_tempvec3_wp, 3 * sizeof(double));
	}
}

void 
UvlmInterfaceBaseElem::MBDyn_UVLM_SendDataToBuf_Curr() {

	/* formulation for calculation
			mbdynuvlm_x = x + mbdynuvlm_f
			mbdynuvlm_R = R
			mbdynuvlm_v = xp + mbdynuvlm_w cross mbdynuvlm_f
			mbdynuvlm_w = w
			mbdynuvlm_a = xpp + mbdynuvlm_wp cross mbdynuvlm_f + mbdynuvlm_w cross mbdynuvlm_w cross mbdynuvlm_f
			mbdynuvlm_wp = wp
	*/
	for (unsigned i = 0; i < MBDyn_CE_NodesNum; i++)
	{
		const MBDYN_UVLM_POINTDATA& mbdynuvlm_point = MBDyn_UVLM_Nodes[i];
		//- Rotation and position
		const Mat3x3 & mbdynuvlm_R = mbdynuvlm_point.pMBDyn_UVLM_Node->GetRCurr();
		const Mat3x3 & mbdynce_Rh = MBDyn_UVLM_Nodes[i].MBDyn_UVLM_RhM;             //- relative orientation of the marker in local ref.
		Mat3x3 mbdynuvlm_R_marker = mbdynuvlm_R * mbdynuvlm_Rh; //- absolute orientation of the marker in local ref.
		Vec3 mbdynuvlm_f = mbdynuvlm_R * mbdynuvlm_point.MBDyn_UVLM_Offset;
		Vec3 mbdynuvlm_x = mbdynuvlm_point.pMBDyn_UVLM_Node->GetXCurr() + mbdynuvlm_f;
		//- Angular velocity and velocity
		const Vec3 &mbdynuvlm_w = mbdynuvlm_point.pMBDyn_UVLM_Node->GetWCurr();
		Vec3 mbdynuvlm_wCrossf = mbdynuvlm_w.Cross(mbdynuvlm_f);
		Vec3 mbdynuvlm_v = mbdynuvlm_point.pMBDyn_UVLM_Node->GetVCurr() + mbdynuvlm_wCrossf;
		//- Angular acceleration and acceleration
		const Vec3 &mbdynuvlm_wp = mbdynuvlm_point.pMBDyn_UVLM_Node->GetWPCurr();
		Vec3 mbdynuvlm_a = mbdynuvlm_point.pMBDyn_UVLM_Node->GetXPPCurr() + mbdynuvlm_wp.Cross(mbdynuvlm_f) + mbdynuvlm_w.Cross(mbdynuvlm_wCrossf);

		double mbdynuvlm_tempvec3_x[3];
		double mbdynuvlm_tempvec3_v[3];
		double mbdynuvlm_tempvec3_a[3];
		double mbdynuvlm_tempvec3_w[3];
		double mbdynuvlm_tempvec3_wp[3];
		MBDyn_UVLM_Vec3D(mbdynuvlm_x, mbdynuvlm_tempvec3_x, MBDyn_UVLM_Scale[0]);
		MBDyn_UVLM_Vec3D(mbdynuvlm_v, mbdynuvlm_tempvec3_v, MBDyn_UVLM_Scale[0]);
		MBDyn_UVLM_Vec3D(mbdynuvlm_a, mbdynuvlm_tempvec3_a, MBDyn_UVLM_Scale[0]);
		MBDyn_UVLM_Vec3D(mbdynuvlm_w, mbdynuvlm_tempvec3_w, 1.0);
		MBDyn_UVLM_Vec3D(mbdynuvlm_wp, mbdynuvlm_tempvec3_wp, 1.0);
		double mbdynuvlm_tempmat3x3_R[9];
		MBDyn_UVLM_Mat3x3D(mbdynuvlm_R_marker, mbdynuvlm_tempmat3x3_R);


		memcpy(&pMBDyn_UVLM_CouplingKinematic_x[3 * i], mbdynuvlm_tempvec3_x, 3 * sizeof(double));
		memcpy(&pMBDyn_UVLM_CouplingKinematic_R[9 * i], mbdynuvlm_tempmat3x3_R, 9 * sizeof(double));
		memcpy(&pMBDyn_UVLM_CouplingKinematic_xp[3 * i], mbdynuvlm_tempvec3_v, 3 * sizeof(double));
		memcpy(&pMBDyn_UVLM_CouplingKinematic_omega[3 * i], mbdynuvlm_tempvec3_w, 3 * sizeof(double));
		memcpy(&pMBDyn_UVLM_CouplingKinematic_xpp[3 * i], mbdynuvlm_tempvec3_a, 3 * sizeof(double));
		memcpy(&pMBDyn_UVLM_CouplingKinematic_omegap[3 * i], mbdynuvlm_tempvec3_wp, 3 * sizeof(double));
	}
}


/**************************Private function definations**************************/
void 
UvlmInterfaceBaseElem::MBDyn_UVLM_Vec3D(const Vec3& mbdynuvlm_Vec3, double* mbdynuvlm_temp, 
	double MBDyn_UVLM_LengthScale) const {

	mbdynuvlm_temp[0] = MBDyn_UVLM_LengthScale * static_cast<double>(*(mbdynuvlm_Vec3.pGetVec()));
	mbdynuvlm_temp[1] = MBDyn_UVLM_LengthScale * static_cast<double>(*(mbdynuvlm_Vec3.pGetVec() + 1));
	mbdynuvlm_temp[2] = MBDyn_UVLM_LengthScale * static_cast<double>(*(mbdynuvlm_Vec3.pGetVec() + 2));
}

void 
UvlmInterfaceBaseElem::MBDyn_UVLM_Mat3x3D(const Mat3x3& mbdynuvlm_Mat3x3, double* mbdynuvlm_temp) const {

	for (unsigned i = 0; i < 9; i++) {
		mbdynuvlm_temp[i] = static_cast<double>(mbdynuvlm_Mat3x3.pGetMat()[i]);
	}
}
/********************************************************************************/


extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	UserDefinedElemRead *rf = new UDERead<UvlmInterfaceBaseElem>;

	if (!SetUDE("UVLM", rf)) {
		delete rf;

		silent_cerr("module-UVLM: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

