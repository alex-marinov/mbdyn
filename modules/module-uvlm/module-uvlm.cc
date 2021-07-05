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
"Module: 	template2						\n"
"Author: 	Pierangelo Masarati <masarati@aero.polimi.it>		\n"
"Organization:	Dipartimento di Ingegneria Aerospaziale			\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it/				\n"
"									\n"
"	All rights reserved						\n"
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
			mbdynce_x = x + mbdynce_f
			mbdynce_R = R
			mbdynce_v = xp + mbdynce_w cross mbdynce_f
			mbdynce_w = w
			mbdynce_a = xpp + mbdynce_wp cross mbdynce_f + mbdynce_w cross mbdynce_w cross mbdynce_f
			mbdynce_wp = wp
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

