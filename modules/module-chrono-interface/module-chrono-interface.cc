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

 /*
  * With the contribution of Runsen Zhang <runsen.zhang@polimi.it>
  * during Google Summer of Code 2020
  */


#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <iostream>
#include <cfloat>
#include <vector>

#include "elem.h"
#include "strnode.h"
#include "dataman.h"
#include "userelem.h"

#include "module-chrono-interface.h"
#include "mbdyn_ce.h"

ChronoInterfaceBaseElem::ChronoInterfaceBaseElem(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO),
MBDyn_CE_CEModel_Converged(pDM),
bMBDyn_CE_FirstSend(true),
bMBDyn_CE_CEModel_DoStepDynamics(true)
{
	// Read element: obtain information from MBDyn script
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	module-chrono-interface						\n"
"Organization:	Dipartimento di Ingegneria Aerospaziale			\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it/				\n"
"									\n"
"	All rights reserved						\n"
"user defined, ID, ChronoInterface, chrono, tight, \n"
"coordinate, Vec3, Mat3x3... \n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}
	
	/* read information from script - start*/
	// read the coupling type
    MBDyn_CE_CouplingIter_Count=0; 
	MBDyn_CE_CouplingType = ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_NONE;
	if (HP.IsKeyWord("coupling")){
		if (HP.IsKeyWord("none")){
			MBDyn_CE_CouplingType = ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_NONE;
			MBDyn_CE_CouplingIter_Max = 1;
		}
		else if (HP.IsKeyWord("loose")){
			MBDyn_CE_CouplingType = ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_LOOSE;
			MBDyn_CE_CouplingIter_Max = 1;
			std::cout << "loose coupling C::E interface: not implemented" << std::endl;
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		else if (HP.IsKeyWord("tight")){
			MBDyn_CE_CouplingType = ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_TIGHT;
			MBDyn_CE_CouplingIter_Max = HP.GetInt();
		}
		else{
			MBDyn_CE_CouplingType = HP.GetInt();
			if (MBDyn_CE_CouplingType > 1)
			{
				MBDyn_CE_CouplingIter_Max = HP.GetInt();
				std::cout << "multirate coupling C::E interface: not implemented" << std::endl;
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
	}

	// read the Units relationship between MBDyn and C::E
	MBDyn_CE_CEScale[0] = 1.0; // default value for length scale
	MBDyn_CE_CEScale[1] = 1.0; // default value for mass scale
	
	if (HP.IsKeyWord("length" "scale"))
	{
		MBDyn_CE_CEScale[0] = HP.GetReal(); // legnth scale: obtained by input
		MBDyn_CE_CEScale[1] = MBDyn_CE_CEScale[0]; // by default, mass scale = length scale;
	}
	if (HP.IsKeyWord("mass" "scale"))
	{
		MBDyn_CE_CEScale[1] = HP.GetReal(); // mass scale: obtained by input
	}
	MBDyn_CE_CEScale[2] = MBDyn_CE_CEScale[1] * MBDyn_CE_CEScale[0]; // force scale
	MBDyn_CE_CEScale[3] = MBDyn_CE_CEScale[2] * MBDyn_CE_CEScale[0]; // moment scale
	std::cout << "read C::E element\n";

	// read the coupling node information
	int MBDyn_CE_NodesNum;
	if (HP.IsKeyWord("nodes" "number"))
	{
		MBDyn_CE_NodesNum = HP.GetInt();
	}	
	if (MBDyn_CE_NodesNum <= 0) // for coupling none ?
	{
		silent_cerr("ChronoInterface(" << uLabel << "): illegal node number " << MBDyn_CE_NodesNum <<
			" at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	MBDyn_CE_Nodes.resize(MBDyn_CE_NodesNum);
	MBDyn_CE_CEModel_Label.resize(MBDyn_CE_NodesNum + 1);
	for (int i = 0; i < MBDyn_CE_NodesNum; i++)
	{
		MBDyn_CE_Nodes[i].pMBDyn_CE_Node = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP); // nodes in MBDyn model
		ReferenceFrame pNode_RF(MBDyn_CE_Nodes[i].pMBDyn_CE_Node); //get its coordinate
		if (HP.IsKeyWord("offset")){
			MBDyn_CE_Nodes[i].MBDyn_CE_Offset = HP.GetPosRel(pNode_RF); // return offset in the node coordinate
		}
		else {
			MBDyn_CE_Nodes[i].MBDyn_CE_Offset= Zero3;
		}
		MBDyn_CE_Nodes[i].MBDyn_CE_F = Zero3;
		MBDyn_CE_Nodes[i].MBDyn_CE_M = Zero3;
		MBDyn_CE_Nodes[i].MBDyn_CE_uLabel = MBDyn_CE_Nodes[i].pMBDyn_CE_Node->GetLabel();
		// get coupling bodies in C::E model, bodies' ID;
		// nodes' motion is imposed between C::E bodies and C::E Ground.
		// the order should be consistent with the nodes; one node corresponds to one body
		MBDyn_CE_Nodes[i].MBDyn_CE_CEBody_Label = HP.GetInt(); // coupling bodies in C::E model, bodies' ID
		MBDyn_CE_CEModel_Label[i].MBDyn_CE_CEBody_Label = MBDyn_CE_Nodes[i].MBDyn_CE_CEBody_Label; // coupling bodies in C::E model, bodies' ID
	}

	// read C::E ground information
	// read the chrono reference::to do!!!!
	// default case
	MBDyn_CE_CEModel_Label[MBDyn_CE_NodesNum].MBDyn_CE_CEBody_Label = 0; // by default: ID of C::E ground in C::E model =0
	double mbdyn_ce_ref_x[3] = {0.0, 0.0, 0.0};
	double mbdyn_ce_ref_R[9];
	for (unsigned i = 0; i < 3;i++){
		mbdyn_ce_ref_x[i] = MBDyn_CE_CEScale[0] * mbdyn_ce_ref_x[i];
	}
	for (unsigned i = 0; i < 9; i++) // default: two system has the same coordinate
	{
		mbdyn_ce_ref_R[i] = MBDyn_CE_CEScale[0] * 0.0;
	}
	mbdyn_ce_ref_R[0] = MBDyn_CE_CEScale[0]* 1.0;
	mbdyn_ce_ref_R[4] = MBDyn_CE_CEScale[0]* 1.0;
	mbdyn_ce_ref_R[8] = MBDyn_CE_CEScale[0]* 1.0;
	// or read from inputs
	if (HP.IsKeyWord("ground"))
	{
		MBDyn_CE_CEModel_Label[MBDyn_CE_NodesNum].MBDyn_CE_CEBody_Label = HP.GetInt(); // the ID of the ground body in C::E
		Vec3 mbdyn_ce_ref_x_Vec3 = HP.GetVec3(); // get C::E ground position information, with respect to the ref GLOBAL
		Mat3x3 mbdyn_ce_ref_R_Mat3x3=HP.GetMatR2vec(); // get C::E ground orietation information
		for (unsigned i = 0; i < 3; i++)
		{
			mbdyn_ce_ref_x[i] = (mbdyn_ce_ref_x_Vec3.pGetVec())[i];
		}
		for (unsigned i = 0; i < 9; i++)
		{
			mbdyn_ce_ref_R[i] = (mbdyn_ce_ref_R_Mat3x3.pGetMat())[i];
		}
	}
	/* read information from script - end*/

	/* initial public vectors (containers for the coupling variables) - start*/
	// allocate space for coupling variables
	MBDyn_CE_CouplingSize[0] = MBDyn_CE_NodesNum * (3 + 9 + 3 + 3 + 3 + 3) + 12; //kinematic motion + 12 (for the coordinate in chrono)
	MBDyn_CE_CouplingSize[1] = MBDyn_CE_NodesNum*(3 + 3); //dynamic variables

	MBDyn_CE_CouplingKinematic.resize(MBDyn_CE_CouplingSize[0], 0.0);
	MBDyn_CE_CouplingDynamic.resize(MBDyn_CE_CouplingSize[1], 0.0);

	pMBDyn_CE_CouplingKinematic_x = &MBDyn_CE_CouplingKinematic[0];
	pMBDyn_CE_CouplingKinematic_R = &MBDyn_CE_CouplingKinematic[3*MBDyn_CE_NodesNum];
	pMBDyn_CE_CouplingKinematic_xp = &MBDyn_CE_CouplingKinematic[12*MBDyn_CE_NodesNum];
	pMBDyn_CE_CouplingKinematic_omega = &MBDyn_CE_CouplingKinematic[15*MBDyn_CE_NodesNum];
	pMBDyn_CE_CouplingKinematic_xpp = &MBDyn_CE_CouplingKinematic[18*MBDyn_CE_NodesNum];
	pMBDyn_CE_CouplingKinematic_omegap = &MBDyn_CE_CouplingKinematic[21*MBDyn_CE_NodesNum];

	pMBDyn_CE_CEFrame = &MBDyn_CE_CouplingKinematic[24 * MBDyn_CE_NodesNum];
	memcpy(pMBDyn_CE_CEFrame, mbdyn_ce_ref_x, 3 * sizeof(double));
	memcpy((pMBDyn_CE_CEFrame+3), mbdyn_ce_ref_R, 9 * sizeof(double));

	pMBDyn_CE_CouplingDynamic_f = &MBDyn_CE_CouplingDynamic[0];
	pMBDyn_CE_CouplingDynamic_m = &MBDyn_CE_CouplingDynamic[3*MBDyn_CE_NodesNum];
	/* initial public vectors - end*/


	/* initial chrono::engine system - start*/
	// intial CE model, and allocate space for reloading C::E model data
	// the default graivity is: (0.0,-9.81,0.0)
	pMBDyn_CE_CEModel = MBDyn_CE_CEModel_Init(MBDyn_CE_CEModel_Data,pMBDyn_CE_CEFrame, MBDyn_CE_CEScale, MBDyn_CE_CEModel_Label, MBDyn_CE_CouplingType);
	/* initial chrono::engine system - end*/
}

ChronoInterfaceBaseElem::~ChronoInterfaceBaseElem(void)
{
	// destroy private data
	MBDyn_CE_CEModel_Destroy(pMBDyn_CE_CEModel);
}

void 
ChronoInterfaceBaseElem::SetValue(DataManager *pDM,
									   VectorHandler &X, VectorHandler &XP,
									   SimulationEntity::Hints *h)
{
	std::cout << "\t\tSetValue()\n";
	/*switch (MBDyn_CE_CouplingType)
	{
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_NONE:
		//break;
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_TIGHT:
		//break;
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_LOOSE:
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_STSTAGGERED:
	default:
		break;
	}*/
}

void 
ChronoInterfaceBaseElem::Update(const VectorHandler &XCurr,
                        const VectorHandler &XprimeCurr)
{
	// a regular step
	// 1. mbdyn writes kinematic coupling variables to buffer;
	// 2. C::E models reload data;
	// 3. C::E models read the coupling data from buffer;
	if (bMBDyn_CE_CEModel_DoStepDynamics || bMBDyn_CE_FirstSend)
	{
		MBDyn_CE_SendDataToBuf();
		if(MBDyn_CE_CouplingType >=1) // tight coupling or multirate coupling
		{
			if(MBDyn_CE_CEModel_DataReload(pMBDyn_CE_CEModel, MBDyn_CE_CEModel_Data))
			{
				silent_cerr("ChronoInterface(" << uLabel << ") data reloading process is wrong " << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
		MBDyn_CE_CEModel_RecvFromBuf(pMBDyn_CE_CEModel, MBDyn_CE_CouplingKinematic, MBDyn_CE_NodesNum);
		MBDyn_CE_CEModel_DoStepDynamics(pMBDyn_CE_CEModel, 1.0e-3);
		// write coupling force to MBDyn_CE_CouplingDynamic
		MBDyn_CE_CEModel_SendToBuf(pMBDyn_CE_CEModel, MBDyn_CE_CouplingDynamic,pMBDyn_CE_CEFrame, MBDyn_CE_NodesNum, MBDyn_CE_CEScale);
	}
	/*switch (MBDyn_CE_CouplingType)
	{
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_NONE:
		std::cout << "\t\tUpdate()\n";
		break; // do nothing
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_TIGHT:
		if (bMBDyn_CE_CEModel_DoStepDynamics){
			MBDyn_CE_SendDataToBuf()
			if (MBDyn_CE_CEModel_DataReload(pMBDyn_CE_CEModel, MBDyn_CE_CEModel_Data))
			{
				silent_cerr("ChronoInterface(" << uLabel << ") data reloading process is wrong " << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			MBDyn_CE_CEModel_RecvFromBuf(pMBDyn_CE_CEModel, MBDyn_CE_CouplingKinematic, MBDyn_CE_NodesNum);
			MBDyn_CE_CEModel_DoStepDynamics(pMBDyn_CE_CEModel, 1.0e-3);
			// write coupling force to MBDyn_CE_CouplingDynamic
			MBDyn_CE_CEModel_SendToBuf(pMBDyn_CE_CEModel, MBDyn_CE_CouplingDynamic, pMBDyn_CE_CEFrame, MBDyn_CE_NodesNum, MBDyn_CE_CEScale);
		}
		break;
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_LOOSE:	   // to do
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_STSTAGGERED: // to do
	default:
		break; 
	}*/
}

void 
ChronoInterfaceBaseElem::AfterConvergence(const VectorHandler &X,
                                  const VectorHandler &XP)
{
	switch (MBDyn_CE_CouplingType)
	{
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_NONE:
		std::cout << "\t\tAfterConvergence()\n";
		break; // do nothing
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_TIGHT: // to do 
		break;
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_LOOSE: // to do
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_STSTAGGERED: // to do
	default:
		break; 
	}
	return;
}

void 
ChronoInterfaceBaseElem::AfterPredict(VectorHandler &X,
                              VectorHandler &XP)
{
	MBDyn_CE_CouplingIter_Count = 0;
	// After prediction, mark next residual as first
	bMBDyn_CE_FirstSend = true;
	// After prediction, tell CEModel that can do time integration and send back data
	bMBDyn_CE_CEModel_DoStepDynamics = true;

	// 1. mbdyn writes kinematic coupling variables to buffer;
	// 2. C::E models reload data;
	// 3. C::E models read the coupling data from buffer; 
	switch (MBDyn_CE_CouplingType)
	{
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_NONE:
		std::cout << "\t\tAfterPredict()\n";
		MBDyn_CE_CEModel_DoStepDynamics(pMBDyn_CE_CEModel, 1.0e-3); // how to get current time step size;	
		bMBDyn_CE_FirstSend = false;
		bMBDyn_CE_CEModel_DoStepDynamics = false; // only do time integration once
		MBDyn_CE_CEModel_Converged.Set(Converged::State::CONVERGED);
		break; // do nothing
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_TIGHT:
		if(MBDyn_CE_CEModel_DataSave(pMBDyn_CE_CEModel, MBDyn_CE_CEModel_Data))
		{
			silent_cerr("ChronoInterface(" << uLabel << ") data saving process is wrong " << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		Update(X, XP); // do step dynamics at the first time.
		bMBDyn_CE_FirstSend = false;
		bMBDyn_CE_CEModel_DoStepDynamics = true;
		break;
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_LOOSE: // to do
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_STSTAGGERED: // to do
	default: // multirate  to do 
		break; 
	}

}


void
ChronoInterfaceBaseElem::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler& 
ChronoInterfaceBaseElem::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering C::E-interface::AssJac()" << std::endl);
	WorkMat.SetNullMatrix();
	/*switch (MBDyn_CE_CouplingType)
	{
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_NONE:
		std::cout << "\t\tAssJac()\n";
		break; // do nothing
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_TIGHT: // do nothing
		break;
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_LOOSE: // do nothing
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_STSTAGGERED: // do nothing
	default:
		break; 
	}*/
	return WorkMat;
}

SubVectorHandler& 
ChronoInterfaceBaseElem::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	// C::E model do time integration(Update the force)
	/*switch (MBDyn_CE_CouplingType)
	{
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_NONE:
		std::cout << "\t\tAssRes()\n";
		break; // do nothing
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_TIGHT:	
		break;
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_LOOSE:		// to do
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_STSTAGGERED: // to do
	default:
		break; 
	}*/
	
	if (MBDyn_CE_CouplingType == MBDyn_CE_COUPLING::COUPLING_NONE)
	{
		std::cout << "\t\tAssRes()\n";
		WorkVec.Resize(0);
		return WorkVec;
	}

	else if (MBDyn_CE_CouplingType == MBDyn_CE_COUPLING::COUPLING_TIGHT)
	{
		// read coupling force from the buffer;
		MBDyn_CE_RecvDataFromBuf();
		const int iOffset = 6;
		WorkVec.ResizeReset(iOffset * MBDyn_CE_NodesNum);
		for (unsigned i = 0; i < MBDyn_CE_NodesNum; i++) {
			const MBDYN_CE_POINTDATA& point = MBDyn_CE_Nodes[i];

			integer iFirstIndex = point.pMBDyn_CE_Node->iGetFirstMomentumIndex();
			for (int r = 1; r <= iOffset; r++) {
				WorkVec.PutRowIndex(i*iOffset + r, iFirstIndex + r);
			}

			WorkVec.Add(i*iOffset + 1, point.MBDyn_CE_F);
			WorkVec.Add(i*iOffset + 4, point.MBDyn_CE_M + (point.pMBDyn_CE_Node->GetRCurr()*point.MBDyn_CE_Offset).Cross(point.MBDyn_CE_F)); // AssJac
		} 
	}
	return WorkVec;
}

void 
ChronoInterfaceBaseElem::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
ChronoInterfaceBaseElem::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering C::E-interface::InitialAssJac()" << std::endl);
	WorkMat.SetNullMatrix();
	switch (MBDyn_CE_CouplingType)
	{
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_NONE:
		std::cout << "\t\tInitialAssJac()\n";
		break; // do nothing
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_TIGHT: // do nothing
		break;
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_LOOSE: // do nothing
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_STSTAGGERED: // do nothing
	default:
		break; 
	}
	return WorkMat;
}

SubVectorHandler& 
ChronoInterfaceBaseElem::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	WorkVec.ResizeReset(0);
	switch (MBDyn_CE_CouplingType)
	{
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_NONE:
		std::cout << "\t\tInitialAssRes()\n";
		break; // do nothing
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_TIGHT: // do nothing
		break;
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_LOOSE:	   // do nothing
	case ChronoInterfaceBaseElem::MBDyn_CE_COUPLING::COUPLING_STSTAGGERED: // do nothing
	default:
		break; 
	}
	return WorkVec;
}

void
ChronoInterfaceBaseElem::MBDyn_CE_SendDataToBuf()
{
	/*
				mbdynce_x = x + mbdynce_f
				mbdynce_R = R
				mbdynce_v = xp + mbdynce_w cross mbdynce_f
				mbdynce_w = w
				mbdynce_a = xpp + mbdynce_wp cross mbdynce_f + mbdynce_w cross mbdynce_w cross mbdynce_f
				mbdynce_wp = wp
	*/
	for (unsigned mbdynce_i = 0; mbdynce_i < MBDyn_CE_NodesNum; mbdynce_i++)
	{
		const MBDYN_CE_POINTDATA& mbdynce_point = MBDyn_CE_Nodes[mbdynce_i];
		// rotation and position
		const Mat3x3 & mbdynce_R = mbdynce_point.pMBDyn_CE_Node->GetRCurr();
		Vec3 mbdynce_f = mbdynce_R * mbdynce_point.MBDyn_CE_Offset;
		Vec3 mbdynce_x = mbdynce_point.pMBDyn_CE_Node->GetXCurr() + mbdynce_f;
		// angular velocity and velocity
		const Vec3 &mbdynce_w = mbdynce_point.pMBDyn_CE_Node->GetWCurr();
		Vec3 mbdynce_wCrossf = mbdynce_w.Cross(mbdynce_f);
		Vec3 mbdynce_v = mbdynce_point.pMBDyn_CE_Node->GetVCurr() + mbdynce_wCrossf;
		// angular acceleration and acceleration
		const Vec3 &mbdynce_wp = mbdynce_point.pMBDyn_CE_Node->GetWPCurr();
		Vec3 mbdynce_a = mbdynce_point.pMBDyn_CE_Node->GetXPPCurr() + mbdynce_wp.Cross(mbdynce_f) + mbdynce_w.Cross(mbdynce_wCrossf);

		double mbdynce_tempvec3_x[3];
		double mbdynce_tempvec3_v[3];
		double mbdynce_tempvec3_a[3];
		double mbdynce_tempvec3_w[3];
		double mbdynce_tempvec3_wp[3];
		MBDyn_CE_Vec3D(mbdynce_x, mbdynce_tempvec3_x, MBDyn_CE_CEScale[0]);
		MBDyn_CE_Vec3D(mbdynce_v, mbdynce_tempvec3_v, MBDyn_CE_CEScale[0]);
		MBDyn_CE_Vec3D(mbdynce_a, mbdynce_tempvec3_a, MBDyn_CE_CEScale[0]);
		MBDyn_CE_Vec3D(mbdynce_w, mbdynce_tempvec3_w, 1.0);
		MBDyn_CE_Vec3D(mbdynce_wp, mbdynce_tempvec3_wp, 1.0);

		double mbdynce_tempmat3x3_R[9];
		MBDyn_CE_Mat3x3D(mbdynce_R, mbdynce_tempmat3x3_R);

		memcpy(&pMBDyn_CE_CouplingKinematic_x[3 * mbdynce_i], mbdynce_tempvec3_x, 3 * sizeof(double));
		memcpy(&pMBDyn_CE_CouplingKinematic_R[9* mbdynce_i], mbdynce_tempmat3x3_R, 9 * sizeof(double));
		memcpy(&pMBDyn_CE_CouplingKinematic_xp[3 * mbdynce_i], mbdynce_tempvec3_v, 3 * sizeof(double));
		memcpy(&pMBDyn_CE_CouplingKinematic_omega[3 * mbdynce_i], mbdynce_tempvec3_w, 3 * sizeof(double));
		memcpy(&pMBDyn_CE_CouplingKinematic_xpp[3 * mbdynce_i], mbdynce_tempvec3_a, 3 * sizeof(double));
		memcpy(&pMBDyn_CE_CouplingKinematic_omegap[3 * mbdynce_i], mbdynce_tempvec3_wp, 3 * sizeof(double));
	}
}

// Read the data from buffer, and write them to the Vec3 MBDyn_CE_F and Vec3 MBDyn_CE_M;
void 
ChronoInterfaceBaseElem::MBDyn_CE_RecvDataFromBuf()
{
	for (unsigned mbdynce_i = 0; mbdynce_i < MBDyn_CE_NodesNum; mbdynce_i++)
	{
		// read from the buffer
		MBDYN_CE_POINTDATA& mbdynce_point = MBDyn_CE_Nodes[mbdynce_i];
		mbdynce_point.MBDyn_CE_F = Vec3(pMBDyn_CE_CouplingDynamic_f[3*mbdynce_i],pMBDyn_CE_CouplingDynamic_f[3*mbdynce_i+1],pMBDyn_CE_CouplingDynamic_f[3*mbdynce_i+2]);
		mbdynce_point.MBDyn_CE_M = Vec3(pMBDyn_CE_CouplingDynamic_m[3*mbdynce_i],pMBDyn_CE_CouplingDynamic_m[3*mbdynce_i+1],pMBDyn_CE_CouplingDynamic_m[3*mbdynce_i+2]);
	}
	MBDyn_CE_CouplingIter_Count++;
	std::cout << "\t\t" << MBDyn_CE_CouplingIter_Count << "\n";
	MBDyn_CE_CEModel_Converged.Set(Converged::State::NOT_CONVERGED);
	bMBDyn_CE_CEModel_DoStepDynamics = true;
	if (MBDyn_CE_CouplingIter_Count >= MBDyn_CE_CouplingIter_Max)  // reach the max coupling iterations or convergence
	{
		MBDyn_CE_CEModel_Converged.Set(Converged::State::CONVERGED);
		bMBDyn_CE_CEModel_DoStepDynamics = false; // C::E don't do time integration until next step;
	}
}


void
ChronoInterfaceBaseElem::Output(OutputHandler& OH) const
{
	// should do something useful
	NO_OP;
}

std::ostream&
ChronoInterfaceBaseElem::Restart(std::ostream& out) const
{
	return out << "# ModuleChronoInterface: is doing now" << std::endl;
}

/* private functions: start*/
void 
ChronoInterfaceBaseElem::MBDyn_CE_Vec3D(const Vec3& mbdynce_Vec3, double* mbdynce_temp, double MBDyn_CE_CELengthScale)
{
	mbdynce_temp[0] = MBDyn_CE_CELengthScale*static_cast<double>(*(mbdynce_Vec3.pGetVec()));
	mbdynce_temp[1] = MBDyn_CE_CELengthScale*static_cast<double>(*(mbdynce_Vec3.pGetVec() + 1));
	mbdynce_temp[2] = MBDyn_CE_CELengthScale*static_cast<double>(*(mbdynce_Vec3.pGetVec() + 2));
}

void 
ChronoInterfaceBaseElem::MBDyn_CE_Mat3x3D(const Mat3x3& mbdynce_Mat3x3, double *mbdynce_temp)
{
	for (unsigned i = 0; i < 8;i++)
	{
		*(mbdynce_temp + i) = static_cast<double>(*(mbdynce_Mat3x3.pGetMat() + i));
	}
}
/* private functions: end*/

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	UserDefinedElemRead *rf = new UDERead<ChronoInterfaceBaseElem>; // or new ChronoInterfaceElemRead;
	std::cout << "create your C::E models:"
			  << "\n";
	if (!SetUDE("ChronoInterface", rf)) {
		delete rf;

		silent_cerr("module-Chrono-Interface: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}