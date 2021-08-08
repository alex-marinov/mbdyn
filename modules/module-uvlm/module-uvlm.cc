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
#include "matvec3.h"
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


	// Read information from the script    [START]
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

	//--------------- 3. <coupling_bodies> ------------------------
	//---------- default coupling data
	//----- nodes number
	MBDyn_UVLM_NodesNum = 0;
	//----- UVLM ground reference
	std::vector<double> mbdyn_uvlm_ref_x(3, 0.0);
	std::vector<double> mbdyn_uvlm_ref_R(9, 0.0);
	for (unsigned i = 0; i < 3; i++) //- default: two system has the same coordinate
	{
		mbdyn_uvlm_ref_x[i] = MBDyn_UVLM_Scale[0] * mbdyn_uvlm_ref_x[i];
	}
	for (unsigned i = 0; i < 9; i++)
	{
		mbdyn_uvlm_ref_R[i] = 0.0;
	}
	mbdyn_uvlm_ref_R[0] = 1.0;
	mbdyn_uvlm_ref_R[4] = 1.0;
	mbdyn_uvlm_ref_R[8] = 1.0;

	//---------- values obtained from scripts
	if (MBDyn_UVLM_CouplingType >= -1) {		//- if coupled

		//---- Total number of coupling nodes 
		MBDyn_UVLM_NodesNum = 0;
		//---- read the ID of the aerodynamic beam 
		unsigned int uBeam = (unsigned int)HP.GetInt();

		DEBUGLCOUT(MYDEBUG_INPUT, "Linked to beam: " << uBeam << std::endl);

		Elem* p = pDM->pFindElem(Elem::BEAM, uBeam);
		if (p == 0) {
			silent_cerr("Beam3(" << uBeam << ") not defined "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		Beam *pBeam = dynamic_cast<Beam *>(p);
		if (pBeam == 0) {
			silent_cerr("Beam(" << uBeam << ") is not a Beam3 "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		ASSERT(pBeam != 0);

		/* Node 1: */

		/* Offset of the aerodynamic body w.r.t. the node */
		const StructNode* pNode1 = pBeam->pGetNode(1);
		
		ReferenceFrame RF(pNode1);
		Vec3 f1(HP.GetPosRel(RF));
		DEBUGLCOUT(MYDEBUG_INPUT, "Node 1 offset: " << f1 << std::endl);

		Mat3x3 Ra1(HP.GetRotRel(RF));
		DEBUGLCOUT(MYDEBUG_INPUT,
			"Node 1 rotation matrix: " << std::endl << Ra1 << std::endl);
		MBDyn_UVLM_NodesNum += 1;
		MBDyn_UVLM_AeroNodes.push_back(pNode1);

		/* Node 2: */

		/* Offset of the aerodynamic body w.r.t. the node */
		const StructNode* pNode2 = pBeam->pGetNode(2);

		RF = ReferenceFrame(pNode2);
		Vec3 f2(HP.GetPosRel(RF));
		DEBUGLCOUT(MYDEBUG_INPUT, "Node 2 offset: " << f2 << std::endl);

		Mat3x3 Ra2(HP.GetRotRel(RF));
		DEBUGLCOUT(MYDEBUG_INPUT,
			"Node 2 rotation matrix: " << std::endl << Ra2 << std::endl);
		MBDyn_UVLM_NodesNum += 1;
		MBDyn_UVLM_AeroNodes.push_back(pNode2);

		/* Node 3: */

		/* Offset of the aerodynamic body w.r.t. the node */
		const StructNode* pNode3 = pBeam->pGetNode(3);

		RF = ReferenceFrame(pNode3);
		Vec3 f3(HP.GetPosRel(RF));
		DEBUGLCOUT(MYDEBUG_INPUT, "Node 3 offset: " << f3 << std::endl);

		Mat3x3 Ra3(HP.GetRotRel(RF));
		DEBUGLCOUT(MYDEBUG_INPUT,
			"Node 3 rotation matrix: " << std::endl << Ra3 << std::endl);
		MBDyn_UVLM_NodesNum += 1;
		MBDyn_UVLM_AeroNodes.push_back(pNode3);

		//----- preparing for reading information of coupling nodes
		MBDyn_UVLM_Nodes.resize(MBDyn_UVLM_NodesNum);
		MBDyn_UVLM_Model_Label.resize(MBDyn_UVLM_NodesNum + 1);                //- node number + 1 (UVLM Ground);
		MBDyn_UVLM_Model_Label[MBDyn_UVLM_NodesNum].MBDyn_UVLM_Body_Label = 0; //- ID of UVLM ground in UVLM model = 0;
		MBDyn_UVLM_Model_Label[MBDyn_UVLM_NodesNum].bMBDyn_UVLM_Body_Output = false; //- Output: no
		//--------reading information of coupling nodes
		for (int i = 0; i < MBDyn_UVLM_NodesNum; ++i) {

			//---------------  <mbdyn_nodes_info> ---------------
			MBDyn_UVLM_Nodes[i].pMBDyn_UVLM_Node = MBDyn_UVLM_AeroNodes[i];        //- pointer to nodes in MBDyn model;
			const DynamicStructNode *temp_pMBDyn_UVLM_Node = dynamic_cast<const DynamicStructNode *>(MBDyn_UVLM_Nodes[i].pMBDyn_UVLM_Node);
			if (temp_pMBDyn_UVLM_Node == NULL) {     //- not a dynamic structural node
				silent_cerr("The node" << MBDyn_UVLM_Nodes[i].pMBDyn_UVLM_Node->GetLabel() << "is not a dynamic structural node" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			const_cast<DynamicStructNode *>(temp_pMBDyn_UVLM_Node)->ComputeAccelerations(true);  //- Also ouput acceleration information;

			//----- coupling node offset
			ReferenceFrame pNode_RF(MBDyn_UVLM_Nodes[i].pMBDyn_UVLM_Node);  //-get node reference
			if (HP.IsKeyWord("offset")) {
				MBDyn_UVLM_Nodes[i].MBDyn_UVLM_Offset = HP.GetPosRel(pNode_RF);  //- return offset in node reference
			}
			else {
				MBDyn_UVLM_Nodes[i].MBDyn_UVLM_Offset = Zero3;
			}
			//----- orientation of the marker in MBDyn. By default, MBDyn_UVLM_RhM == eye;
			//----- orientation of the marker (MBDyn_UVLM_Body_Rh[9]) in UVLM is calculated by this matrix MBDyn_UVLM_RhM;
			MBDyn_UVLM_Nodes[i].MBDyn_UVLM_RhM = Eye3;
			if (HP.IsKeyWord("rotation"))
			{
				MBDyn_UVLM_Nodes[i].MBDyn_UVLM_RhM = HP.GetRotRel(pNode_RF);
			}
			Mat3x3 mbdyn_UVLM_RhM_abs;
			mbdyn_UVLM_RhM_abs = pNode_RF.GetR() * (MBDyn_UVLM_Nodes[i].MBDyn_UVLM_RhM);
			for (unsigned j = 0; j < 9; j++)
			{
				MBDyn_UVLM_Model_Label[i].MBDyn_UVLM_Body_Rh[j] = mbdyn_UVLM_RhM_abs.pGetMat()[j];
				std::cout << "\t\trotation matrix:\t" << MBDyn_UVLM_Model_Label[i].MBDyn_UVLM_Body_Rh[j] << "\n";
			}
			//--------- Forces associated without and with the time derivatives (i.e. forces and dynamic forces)
			MBDyn_UVLM_Nodes[i].MBDyn_UVLM_F = Zero3;
			MBDyn_UVLM_Nodes[i].MBDyn_UVLM_DF = Zero3;
			MBDyn_UVLM_Nodes[i].MBDyn_UVLM_M = Zero3;
			MBDyn_UVLM_Nodes[i].MBDyn_UVLM_uLabel = MBDyn_UVLM_Nodes[i].pMBDyn_UVLM_Node->GetLabel();  //- node label (MBDyn)

			// ----------------------  <UVLM_bodies_info>--------------------
			if (HP.IsKeyWord("UVLM" "body")) {
				MBDyn_UVLM_Nodes[i].MBDyn_UVLM_Body_Label = HP.GetInt();                                       //- bodies' ID, recorded in MBDyn;
				MBDyn_UVLM_Model_Label[i].MBDyn_UVLM_Body_Label = MBDyn_UVLM_Nodes[i].MBDyn_UVLM_Body_Label;   //- bodies' ID, used in UVLM;
			}

			// ------------- offset relative to UVLM body, expressed in body local ref.----
			if (HP.IsKeyWord("offset")) {
				for (unsigned j = 0; j < 3; ++j) {
					MBDyn_UVLM_Model_Label[i].MBDyn_UVLM_Body_Offset[j] = HP.GetReal() * MBDyn_UVLM_Scale[0];    //- return offset in UVLM body reference;
				}
			}
			else {
				for (unsigned j = 0; j < 3; ++j) {
					MBDyn_UVLM_Model_Label[i].MBDyn_UVLM_Body_Offset[j] = 0.0;                                  //- return offset in UVLM body reference;
				}
			}

			//----- position and rotation constraints;
			//----- through it, although mbdyn still pass all postion and rotation to UVLM,
			//----- UVLM body will be constrained to the node frame according to the settings.
			for (unsigned j = 0; j < 6; ++j) {  //---- default case: all positions and rotation are imposed; 
				MBDyn_UVLM_Model_Label[i].bMBDyn_UVLM_Constraint[j] = true;
			}
			if (HP.IsKeyWord("position" "constraint")) {
				for (unsigned j = 0; j < 3; ++j) {
					bool bActive = HP.GetInt();
					MBDyn_UVLM_Model_Label[i].bMBDyn_UVLM_Constraint[j] = bActive;
				}
			}
			if (HP.IsKeyWord("orientation" "constraint")) {
				for (unsigned j = 3; j < 6; ++j) {
					bool bActive = HP.GetInt();
					MBDyn_UVLM_Model_Label[i].bMBDyn_UVLM_Constraint[j] = bActive;
				}
			}


			//------- Write UVLM coupling bodies data to file, or not;
			MBDyn_UVLM_Model_Label[i].bMBDyn_UVLM_Body_Output = false;  //-  by default: no output;
			if (HP.IsKeyWord("output" "uvlmbody")) {
				if (HP.IsKeyWord("yes")) {
					MBDyn_UVLM_Model_Label[i].bMBDyn_UVLM_Body_Output = true;
				}
				else if (HP.IsKeyWord("no")) {
					// NO_OP;
				}
				else {
					silent_cerr("UVLMInterface(" << uLabel << "): lacks keyword at line " << HP.GetLineData() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
			}
		}
		//-------------- <uvlm_ground_info> ----------------------------------
		if (HP.IsKeyWord("ground")) {
			MBDyn_UVLM_Model_Label[MBDyn_UVLM_NodesNum].MBDyn_UVLM_Body_Label = HP.GetInt(); //- the ID of the ground body in UVLM
			Vec3 mbdyn_uvlm_ref_x_Vec3;
			Mat3x3 mbdyn_uvlm_ref_R_Mat3x3;
			ReferenceFrame mbdyn_uvlm_ref(MBDyn_UVLM_Nodes[0].pMBDyn_UVLM_Node);  //- ground ref (keyword 'Node' and 'Local' denotes the ref of first node)
			if (HP.IsKeyWord("position")) {
				mbdyn_uvlm_ref_x_Vec3 = HP.GetPosAbs(mbdyn_uvlm_ref);   //- get UVLM ground position information, with respect to the ref GLOBAL
			}
			if (HP.IsKeyWord("orientation")) {
				mbdyn_uvlm_ref_R_Mat3x3 = HP.GetRotAbs(mbdyn_uvlm_ref);  //- get UVLM ground orientation information
			}

			//- Write UVLM ground's motion, or not.
			MBDyn_UVLM_Model_Label[MBDyn_UVLM_NodesNum].bMBDyn_UVLM_Body_Output = false;   // by default: no output
			if (HP.IsKeyWord("output" "uvlmbody")) {
				if (HP.IsKeyWord("yes")) {
					MBDyn_UVLM_Model_Label[MBDyn_UVLM_NodesNum].bMBDyn_UVLM_Body_Output = true;
				}
				else if (HP.IsKeyWord("no")) {
					// NO_OP;
				}
				else {
					silent_cerr("UVLMInterface(" << uLabel << "): lacks keyword at line " << HP.GetLineData() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
			}
			//- save the ref Vec3 in double [] (can also use functions: void MBDyn_UVLM_Vec3D() and  void MBDyn_UVLM_Mat3x3D());
			for (unsigned i = 0; i < 3; ++i) {
				mbdyn_uvlm_ref_x[i] = (mbdyn_uvlm_ref_x_Vec3.pGetVec())[i] * MBDyn_UVLM_Scale[0];
			}
			for (unsigned i = 0; i < 9; ++i) {
				//- matrix in mbdyn:
				//- [ m_11,  m_12,  m_13 ]
				//- [ m_21,  m_22,  m_23 ]
				//- [ m_31,  m_32,  m_33 ]
				//- R[0]=m_11, R[1]=m_21, R[2]=m_31, R[3]=m_12....
				mbdyn_uvlm_ref_R[i] = (mbdyn_uvlm_ref_R_Mat3x3.pGetMat())[i];
			}
		}
	}
	//--------------- 4. <other_settings> ------------------------
	//----- whether to output kinematics of all bodies
	MBDyn_UVLM_OutputType = MBDYN_UVLM_OUTPUTTYPE::MBDYN_UVLM_OUTPUT_SELECTEDCOUPLINGBODIES;  //- by default: output selected coupling bodies in UVLM
	if (HP.IsKeyWord("output" "all" "uvlmbodies")) {
		MBDyn_UVLM_OutputType = MBDYN_UVLM_OUTPUTTYPE::MBDYN_UVLM_OUTPUT_ALLBODIES;           //- output all bodies in uvlm
	}
	//---------- output information about coupling forces;
	if (HP.IsKeyWord("coupling" "forces" "filename"))
	{
		const char *MBDyn_UVLM_Output_Filename = HP.GetFileName();
		if (MBDyn_UVLM_Output_Filename == NULL)
		{
			silent_cerr("ChronoInterface(" << uLabel << "): unable to get file name for coupling forces " <<
				" at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		MBDyn_UVLM_out_forces.open(MBDyn_UVLM_Output_Filename);
		if (!MBDyn_UVLM_out_forces.is_open())
		{
			silent_cerr("ChronoInterface(" << uLabel << "): unable to open" << MBDyn_UVLM_Output_Filename <<
				" at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
	//----- printf information in screen ? 
	bMBDyn_UVLM_Verbose = false;
	if (HP.IsKeyWord("verbose"))
	{
		bMBDyn_UVLM_Verbose = true;
	}
	// Read information from the script    [END]


	// Initialization (std::vector<>) - [START]

	// Initialize the constant values for the UVLM part
	MBDyn_UVLM_StepUVLM_settings->gamma_dot_filtering = ;   // Filtering parameter for the Welch filter for the Gamma_dot estimation. Used when ``unsteady_force_contribution`` is ``on``.'
	MBDyn_UVLM_StepUVLM_settings->gust_intensity = ;
	MBDyn_UVLM_StepUVLM_settings->gust_length = ;
	MBDyn_UVLM_StepUVLM_settings->gust_offset = ;
	MBDyn_UVLM_StepUVLM_settings->gust_shape = ;
	MBDyn_UVLM_StepUVLM_settings->n_time_steps = ;
	MBDyn_UVLM_StepUVLM_settings->relative_motion = "off";   // "on" when there is NOT free flight
	MBDyn_UVLM_StepUVLM_settings->span = ;
	MBDyn_UVLM_StepUVLM_settings->velocity_field_generator = ;
	
	MBDyn_UVLM_Aerogrid_settings->aligned_grid = true;
	MBDyn_UVLM_Aerogrid_settings->freestream_dir = ;
	MBDyn_UVLM_Aerogrid_settings->mstar = ;  // Number of chordwise panels
	MBDyn_UVLM_Aerogrid_settings->unsteady = false;
	MBDyn_UVLM_Aerogrid_settings->wake_shape_generator = "StraightWake";

	MBDyn_UVLM_StraightWake_settings->dt = ;  // Times step;
	MBDyn_UVLM_StraightWake_settings->dx1 = -1.0;  // Size of the first wake panel
	MBDyn_UVLM_StraightWake_settings->dxmax = -1.0; // Maximum panel size
	MBDyn_UVLM_StraightWake_settings->ndx1 = 1; // Number of panels with size "dx1"
	MBDyn_UVLM_StraightWake_settings->r = 1.0; // Growth rate after "ndx1" panels
	MBDyn_UVLM_StraightWake_settings->u_inf = ;
	MBDyn_UVLM_StraightWake_settings->u_inf_direction = ;

	MBDyn_UVLM_UVMopts->cfl1 = true;  // 'If it is ``True``, it assumes that the discretisation complies with CFL=1'
	MBDyn_UVLM_UVMopts->convection_scheme = 3;  // '``0``: fixed wake, ' \'``2``: convected with background flow;' \'``3``: full force-free wake'
	MBDyn_UVLM_UVMopts->convect_wake = true;
	MBDyn_UVLM_UVMopts->dt = ;
	MBDyn_UVLM_UVMopts->filter_method = 0;  // Method to filter the points: no filter (0) moving average(2)'
	MBDyn_UVLM_UVMopts->ImageMethod = false;
	MBDyn_UVLM_UVMopts->interp_coords = 0; // Coordinates to use for wake description: cartesian(0) or cylindrical_z(1)'
	MBDyn_UVLM_UVMopts->interp_method = 0; // Method of interpolation: linear(0), parabolic(1), splines(2), slerp around z(3), slerp around yaw_slerp(4)'
	MBDyn_UVLM_UVMopts->iterative_precond = false;
	MBDyn_UVLM_UVMopts->iterative_solver = false;
	MBDyn_UVLM_UVMopts->iterative_tol = 1.0e-4;
	MBDyn_UVLM_UVMopts->NumCores = 0;  // Number of cores to use in UVLM solver
	MBDyn_UVLM_UVMopts->NumSurfaces = ;  // Number of surfaces
	MBDyn_UVLM_UVMopts->quasi_steady = false;  // Use quasi-steady approximation in UVLM
	MBDyn_UVLM_UVMopts->vortex_radius = 1.0e-6;  //Distance between points below which induction is not computed
	MBDyn_UVLM_UVMopts->vortex_radius_wake_ind = 1.0e-6; // Distance between points below which induction is not computed in the wake convection
	MBDyn_UVLM_UVMopts->yaw_slerp = 0;  // Yaw angle in radians to be used when interp_metod == 4

	MBDyn_UVLM_FlightConditions->c_ref = ;
	MBDyn_UVLM_FlightConditions->rho = ;
	MBDyn_UVLM_FlightConditions->uinf = ;
	MBDyn_UVLM_FlightConditions->uinf_direction = ;

	//---------- allocate space for coupling variables (std::vectors for the coupling variables, motion and force)
	//- kinematic motion + 12 (for the global coordinate in uvlm)
	MBDyn_UVLM_CouplingSize.Size_Kinematic = MBDyn_UVLM_NodesNum * (3 + 9 + 3 + 3 + 3 + 3) + 12; //- motion
	MBDyn_UVLM_CouplingSize.Size_Dynamic = MBDyn_UVLM_NodesNum * (3 + 3 + 3); //- dynamic variables

	MBDyn_UVLM_CouplingKinematic.resize(MBDyn_UVLM_CouplingSize.Size_Kinematic, 0.0);
	MBDyn_UVLM_CouplingDynamic.resize(MBDyn_UVLM_CouplingSize.Size_Dynamic, 0.0);
	MBDyn_UVLM_CouplingDynamic_pre.resize(MBDyn_UVLM_CouplingSize.Size_Dynamic, 0.0); //- vector for last iteration
	//- pointer to motion
	/*
	std::copy(MBDyn_UVLM_CouplingKinematic.begin(), MBDyn_UVLM_CouplingKinematic.begin() + 3 * MBDyn_UVLM_NodesNum, pMBDyn_UVLM_CouplingKinematic_x);
	std::copy(MBDyn_UVLM_CouplingKinematic.begin() + 3 * MBDyn_UVLM_NodesNum, MBDyn_UVLM_CouplingKinematic.begin() + 12 * MBDyn_UVLM_NodesNum, pMBDyn_UVLM_CouplingKinematic_R);
	std::copy(MBDyn_UVLM_CouplingKinematic.begin() + 12 * MBDyn_UVLM_NodesNum, MBDyn_UVLM_CouplingKinematic.begin() + 15 * MBDyn_UVLM_NodesNum, pMBDyn_UVLM_CouplingKinematic_xp);
	std::copy(MBDyn_UVLM_CouplingKinematic.begin() + 15 * MBDyn_UVLM_NodesNum, MBDyn_UVLM_CouplingKinematic.begin() + 18 * MBDyn_UVLM_NodesNum, pMBDyn_UVLM_CouplingKinematic_omega);
	std::copy(MBDyn_UVLM_CouplingKinematic.begin() + 18 * MBDyn_UVLM_NodesNum, MBDyn_UVLM_CouplingKinematic.begin() + 21 * MBDyn_UVLM_NodesNum, pMBDyn_UVLM_CouplingKinematic_xpp);
	std::copy(MBDyn_UVLM_CouplingKinematic.begin() + 21 * MBDyn_UVLM_NodesNum, MBDyn_UVLM_CouplingKinematic.begin() + 24 * MBDyn_UVLM_NodesNum, pMBDyn_UVLM_CouplingKinematic_omegap);
	*/

	//- UVLM ground coordinate  
	MBDyn_UVLM_CouplingKinematic[24 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_ref_x[0];
	MBDyn_UVLM_CouplingKinematic[24 * MBDyn_UVLM_NodesNum + 1] = mbdyn_uvlm_ref_x[1];
	MBDyn_UVLM_CouplingKinematic[24 * MBDyn_UVLM_NodesNum + 2] = mbdyn_uvlm_ref_x[2];

	MBDyn_UVLM_CouplingKinematic[24 * MBDyn_UVLM_NodesNum + 3] = mbdyn_uvlm_ref_R[0];
	MBDyn_UVLM_CouplingKinematic[24 * MBDyn_UVLM_NodesNum + 4] = mbdyn_uvlm_ref_R[1];
	MBDyn_UVLM_CouplingKinematic[24 * MBDyn_UVLM_NodesNum + 5] = mbdyn_uvlm_ref_R[2];
	MBDyn_UVLM_CouplingKinematic[24 * MBDyn_UVLM_NodesNum + 6] = mbdyn_uvlm_ref_R[3];
	MBDyn_UVLM_CouplingKinematic[24 * MBDyn_UVLM_NodesNum + 7] = mbdyn_uvlm_ref_R[4];
	MBDyn_UVLM_CouplingKinematic[24 * MBDyn_UVLM_NodesNum + 8] = mbdyn_uvlm_ref_R[5];
	MBDyn_UVLM_CouplingKinematic[24 * MBDyn_UVLM_NodesNum + 9] = mbdyn_uvlm_ref_R[6];
	MBDyn_UVLM_CouplingKinematic[24 * MBDyn_UVLM_NodesNum + 10] = mbdyn_uvlm_ref_R[7];
	MBDyn_UVLM_CouplingKinematic[24 * MBDyn_UVLM_NodesNum + 11] = mbdyn_uvlm_ref_R[8];

	/*
	//- pointer to force
	pMBDyn_UVLM_CouplingDynamic_f = &MBDyn_UVLM_CouplingDynamic[0];
	pMBDyn_UVLM_CouplingDynamic_m = &MBDyn_UVLM_CouplingDynamic[3 * MBDyn_UVLM_NodesNum];
	//- pointer to force of last iteration (not time step)
	pMBDyn_UVLM_CouplingDynamic_f_pre = &MBDyn_UVLM_CouplingDynamic_pre[0];
	pMBDyn_UVLM_CouplingDynamic_m_pre = &MBDyn_UVLM_CouplingDynamic_pre[3 * MBDyn_UVLM_NodesNum];
	*/

	// initialization (std::vector<>) - [END]


	// Initial UVLM system - [START]
	//---------- initial UVLM model, and allocate space for reloading UVLM model data
	MBDyn_UVLM_Model_Init(MBDyn_UVLM_StepUVLM_settings, MBDyn_UVLM_Aerogrid_settings, MBDyn_UVLM_StraightWake_settings, 
		MBDyn_UVLM_UVMopts, MBDyn_UVLM_FlightConditions, MBDyn_UVLM_Beam_inputs, MBDyn_UVLM_Aero_inputs,
		MBDyn_UVLM_StraightWake, MBDyn_UVLM_SteadyVelocityField, MBDyn_UVLM_UvlmLibVar, VMoptions, UVMoptions, 
		FlightConditions, MBDyn_UVLM_NodesNum);
	// Initial UVLM system - [END]

}

UvlmInterfaceBaseElem::~UvlmInterfaceBaseElem() {

	NO_OP;
}

void
UvlmInterfaceBaseElem::SetValue(DataManager *pDM, VectorHandler &X, VectorHandler &XP,
	SimulationEntity::Hints *h) {

	pedantic_cerr("\t MBDyn::SetValue() \n");
	//- SetValue function doesn't call MBDyn_UVLM_Model_DoStepDynamics()
	//- just save data if using tight coupling scheme for the first step
	MBDyn_UVLM_CouplingIter_Count = MBDyn_UVLM_CouplingIter_Max;
	bMBDyn_UVLM_FirstSend = false;
	bMBDyn_UVLM_Model_DoStepDynamics = false;

	//---------- set gravity 
	//- the default graivity is: (0.0,-9.81,0.0)
	Vec3 mbdyn_uvlm_mbdyn_gravity_vec3;
	std::vector<double> mbdyn_uvlm_mbdyn_gravity(3);
	Vec3 mbdyn_uvlm_mbdyn_vec3 = Zero3; // arm of the gravity. is not used.
	bool bmbdyn_uvlm_gravity;
	GravityOwner::bGetGravity(mbdyn_uvlm_mbdyn_vec3, mbdyn_uvlm_mbdyn_gravity_vec3);
	MBDyn_UVLM_Vec3D(mbdyn_uvlm_mbdyn_gravity_vec3, mbdyn_uvlm_mbdyn_gravity, MBDyn_UVLM_Scale[0]); //- transfer Vec3 to double vector
	std::cout << "\t\tgravity is: " << mbdyn_uvlm_mbdyn_gravity[0] << "\t" << mbdyn_uvlm_mbdyn_gravity[1] << "\t" << mbdyn_uvlm_mbdyn_gravity[2] << "\n"; //- print in screen

	//-------------- Initial send data to Buffer
	MBDyn_UVLM_SendDataToBuf_Curr(); //- send initial data to BUF, where UVLM model can read infromation.

	switch (MBDyn_UVLM_CouplingType) {    //- only tight scheme needs to save the data
	
	case MBDyn_UVLM_COUPLING::COUPLING_NONE:  // - do nothing
		break;
	case MBDyn_UVLM_COUPLING::COUPLING_TIGHT:  //- Save the data
		/*
		!
		!
		!
		!
		!
		*/
		break;
	case MBDyn_UVLM_COUPLING::COUPLING_LOOSE:  // - do nothing
		break;
	case MBDyn_UVLM_COUPLING::COUPLING_STSTAGGERED:  //-  to do
	default:  //-  multirate to do
		break;
	}
}

void
UvlmInterfaceBaseElem::Update(const VectorHandler &XCurr, const VectorHandler &XprimeCurr) {

	//---------- A regular step to update UVLM model
	//---------- 1. mbdyn writes kinematic coupling variables to buffer;
	//---------- 2. UVLM models reload data (only tight coupling scheme);
	//---------- 3. UVLM models read the coupling data from buffer;
	//---------- 4. UVLM models do integration (one step);
	//---------- 5. MBDyn receives data from buffer;
	pedantic_cout("\t MBDyn::Update() \n");
	if (bMBDyn_UVLM_Model_DoStepDynamics || bMBDyn_UVLM_FirstSend) {
		if (MBDyn_UVLM_CouplingType == MBDyn_UVLM_COUPLING::COUPLING_TIGHT) {
			//---------- 1. mbdyn writes kinematic coupling variables to buffer;
			MBDyn_UVLM_SendDataToBuf_Curr();
			//---------- 2. UVLM models reload data (tight coupling scheme);
			/*
			!
			!
			!
			!
			!
			*/
			//---------- steps 3 and steps 4;
			//---------- 3. UVLM models read the coupling data from buffer;
			//---------- 4. UVLM models do integration (one step);
			MBDyn_UVLM_UpdateUVLMModel();
			//---------- 5. MBDyn receives data from Buf;
			MBDyn_UVLM_RecvDataFromBuf();
			//---------- other settings
			bMBDyn_UVLM_FirstSend = false;
			bMBDyn_UVLM_Model_DoStepDynamics = true;
			MBDyn_UVLM_Model_Converged.Set(Converged::State::NOT_CONVERGED);
		}
		else if (MBDyn_UVLM_CouplingType == MBDyn_UVLM_COUPLING::COUPLING_LOOSE) {
			switch (MBDyn_UVLM_CouplingType_loose) {

				case MBDyn_UVLM_COUPLING_LOOSE::LOOSE_EMBEDDED: {
					//---------- 1. mbdyn writes kinematic coupling variables to buffer;
					MBDyn_UVLM_SendDataToBuf_Curr();
					//---------- steps 3 and steps 4;
					//---------- 3. UVLM models read the coupling data from buffer;
					//---------- 4. UVLM models do integration (one step);
					MBDyn_UVLM_UpdateUVLMModel();
					//---------- 5. MBDyn receives data from Buf;
					MBDyn_UVLM_RecvDataFromBuf();
					break;
				}
				case MBDyn_UVLM_COUPLING_LOOSE::LOOSE_JACOBIAN:{ //- x^(P)_k+1=x_k..., f^{P}_{k+1}=f^{P}
					//---------- 1. mbdyn writes kinematic coupling variables to buffer;
					MBDyn_UVLM_SendDataToBuf_Prev();
					//---------- 5. MBDyn receives data from Buf before integration in UVLM. data of last step.
					MBDyn_UVLM_RecvDataFromBuf();
					//---------- steps 3 and steps 4;
					//---------- 3. UVLM models read the coupling data from buffer;
					//---------- 4. UVLM models do integration (one step);
					MBDyn_UVLM_UpdateUVLMModel();
					break;
				}
				case MBDyn_UVLM_COUPLING_LOOSE::LOOSE_GAUSS: //- to do
					break;
				default:
					silent_cerr("unkown loose coupling scheme" << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					break;
			}
			bMBDyn_UVLM_Model_DoStepDynamics = false;
			bMBDyn_UVLM_FirstSend = false;
			MBDyn_UVLM_Model_Converged.Set(Converged::State::CONVERGED);
		}
	}
}

void
UvlmInterfaceBaseElem::BeforePredict(VectorHandler &X, VectorHandler &XP, VectorHandler &XPrev, VectorHandler &XPPrev) const {

	pedantic_cout("\t MBDyn::BeforePredict() \n");
}

void 
UvlmInterfaceBaseElem::AfterConvergence(const VectorHandler &X, const VectorHandler &XP) {
	
	pedantic_cout("\t MBDyn::AfterConvergence() \n");
	switch (MBDyn_UVLM_CouplingType) {

	case MBDyn_UVLM_COUPLING::COUPLING_NONE: //- do nothing
		break;
	case MBDyn_UVLM_COUPLING::COUPLING_TIGHT: //- do nothing
		break;
	case MBDyn_UVLM_COUPLING::COUPLING_LOOSE: //- to do
	case MBDyn_UVLM_COUPLING::COUPLING_STSTAGGERED: //- to do
	default:
		break;
	}
	if (bMBDyn_UVLM_Verbose) {
		MBDyn_UVLM_MBDynPrint();
	}
	return;
}

void 
UvlmInterfaceBaseElem::AfterPredict(VectorHandler &X, VectorHandler &XP) {

	//- first send data after predict in each co-simulation (tight);
	//- after predict saves UVLM data and call MBDyn_UVLM_Model_DoStepDynamics();
	MBDyn_UVLM_CouplingIter_Count = 0;
	bMBDyn_UVLM_FirstSend = true;
	bMBDyn_UVLM_Model_DoStepDynamics = true;

	switch (MBDyn_UVLM_CouplingType) {

	case MBDyn_UVLM_COUPLING::COUPLING_NONE:
		pedantic_cout("\t MBDyn::AfterPredict() - no_coupling \n");
		time_step = m_PDM->pGetDrvHdl()->dGetTimeStep();  //-  get time step

		// Perform one step dynamics
		MBDyn_UVLM_Model_DoStepDynamics(MBDyn_UVLM_StepUVLM_settings, MBDyn_UVLM_Aerogrid_settings, MBDyn_UVLM_StraightWake_settings,
			MBDyn_UVLM_UVMopts, MBDyn_UVLM_FlightConditions, MBDyn_UVLM_Beam_inputs, MBDyn_UVLM_Aero_inputs,
			MBDyn_UVLM_StraightWake, MBDyn_UVLM_SteadyVelocityField, MBDyn_UVLM_Aerogrid,
			MBDyn_UVLM_UvlmLibVar, UVMoptions, FlightConditions, time_step);

		bMBDyn_UVLM_FirstSend = false;
		bMBDyn_UVLM_Model_DoStepDynamics = false;   //- only do time integration once
		MBDyn_UVLM_Model_Converged.Set(Converged::State::CONVERGED);
		break;
	case MBDyn_UVLM_COUPLING::COUPLING_TIGHT:
		pedantic_cout("\t MBDyn::AfterPredict() -tight_embedded \n");
		/*
		!
		!
		!
		!
		!
		!
		*/
		Update(X, XP);  // A regular step in UVLM
		break;
	case MBDyn_UVLM_COUPLING::COUPLING_LOOSE: {

		switch (MBDyn_UVLM_CouplingType_loose) {

		case MBDyn_UVLM_COUPLING_LOOSE::TIGHT:
			break;  //- do nothing (the work is finished in case MBDyn_UVLM_COUPLING::COUPLING_TIGHT)
		case MBDyn_UVLM_COUPLING_LOOSE::LOOSE_EMBEDDED:
			pedantic_cout("\t MBDyn::AfterPredict() -Loose_embedded \n");
			Update(X, XP);  //- A regular step in UVLM
			break;
		case MBDyn_UVLM_COUPLING_LOOSE::LOOSE_JACOBIAN:
			pedantic_cout("\t MBDyn::AfterPredict() -Loose_Jacobian \n");
			Update(X, XP);  //- A regular step in UVLM
			break;
		case MBDyn_UVLM_COUPLING_LOOSE::LOOSE_GAUSS:  //- to do
		default: //- to do
			break;
		}

		break;
	}
	case MBDyn_UVLM_COUPLING::COUPLING_STSTAGGERED: //- to do
	default: //- multirate to do
		break;
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
	const VectorHandler& XPrimeCurr) {

	pedantic_cout("\tMBDyn::AssJac()" << std::endl);
	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler&
UvlmInterfaceBaseElem::AssRes(SubVectorHandler& WorkVec, 
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr) {

	if (MBDyn_UVLM_CouplingType == MBDyn_UVLM_COUPLING::COUPLING_NONE)
	{
		pedantic_cout("\tMBDyn::AssRes()\n");
		WorkVec.Resize(0);
		return WorkVec;
	}
	else if (MBDyn_UVLM_CouplingType >= -1) {   //- coupling case

		//---------- A regular step in MBDyn
		//---------- 1. detect the convergence of coupling data;
		//---------- 2. add forces to Res;
		pedantic_cout("\t MBDyn::AssRes() \n");

		//---------- 1. detect the convergence of coupling data;
		if ((bMBDyn_UVLM_Model_DoStepDynamics) && (MBDyn_UVLM_CouplingType == MBDyn_UVLM_COUPLING::COUPLING_TIGHT)) {  //- tight coupling case
			pedantic_cout("\t \t Ass_Coupling_Error() \n");
			pedantic_cout("\t \t iterations: " << MBDyn_UVLM_CouplingIter_Count << "\n");
			MBDyn_UVLM_Model_Converged.Set(Converged::State::NOT_CONVERGED);
			bMBDyn_UVLM_Model_DoStepDynamics = true;
			if (MBDyn_UVLM_CouplingIter_Count > 0) {  //- iters 2nd
				double mbdyn_uvlm_error;
				mbdyn_uvlm_error = MBDyn_UVLM_CalculateError();
				if (mbdyn_uvlm_error < MBDyn_UVLM_Coupling_Tol || MBDyn_UVLM_CouplingIter_Count >= MBDyn_UVLM_CouplingIter_Max) {
					MBDyn_UVLM_Model_Converged.Set(Converged::State::CONVERGED);
					bMBDyn_UVLM_Model_DoStepDynamics = false;    //- UVLM doesn't need to do time integration until next step;
					pedantic_cout("\t \t Coupling error: " << mbdyn_uvlm_error << "\n");
				}
				else {
					pedantic_cout("\t \t Coupling error: " << mbdyn_uvlm_error << "\n");
				}
			}
			MBDyn_UVLM_CouplingIter_Count++;
		}
		//---------- 2. add forces to Nodes;
		const int iOffset = 6;
		WorkVec.ResizeReset(iOffset * MBDyn_UVLM_NodesNum);
		for (unsigned i = 0; i < MBDyn_UVLM_NodesNum; ++i) {
			const MBDYN_UVLM_POINTDATA& point = MBDyn_UVLM_Nodes[i];

			integer iFirstIndex = point.pMBDyn_UVLM_Node->iGetFirstMomentumIndex();
			for (int r = 1; r <= iOffset; ++r) {
				WorkVec.PutRowIndex(i*iOffset + r, iFirstIndex + r);
			}

			WorkVec.Add(i*iOffset + 1, point.MBDyn_UVLM_F);
			WorkVec.Add(i*iOffset + 4, point.MBDyn_UVLM_DF);
			WorkVec.Add(i*iOffset + 7, point.MBDyn_UVLM_M + (point.pMBDyn_UVLM_Node->GetRCurr()*point.MBDyn_UVLM_Offset).Cross(point.MBDyn_UVLM_F));
			//- save the force of last iteration (only used in tight co-simulation)
			// Force components 
			MBDyn_UVLM_CouplingDynamic_pre[3 * i] = MBDyn_UVLM_CouplingDynamic[3 * i];
			MBDyn_UVLM_CouplingDynamic_pre[3 * i + 1] = MBDyn_UVLM_CouplingDynamic[3 * i + 1];
			MBDyn_UVLM_CouplingDynamic_pre[3 * i + 2] = MBDyn_UVLM_CouplingDynamic[3 * i + 2];
			// Dynamic force components
			MBDyn_UVLM_CouplingDynamic_pre[3 * i + 3 * MBDyn_UVLM_NodesNum] = MBDyn_UVLM_CouplingDynamic[3 * i + 3 * MBDyn_UVLM_NodesNum];
			MBDyn_UVLM_CouplingDynamic_pre[3 * i + 1 + 3 * MBDyn_UVLM_NodesNum] = MBDyn_UVLM_CouplingDynamic[3 * i + 1 + 3 * MBDyn_UVLM_NodesNum];
			MBDyn_UVLM_CouplingDynamic_pre[3 * i + 2 + 3 * MBDyn_UVLM_NodesNum] = MBDyn_UVLM_CouplingDynamic[3 * i + 2 + 3 * MBDyn_UVLM_NodesNum];
			//- Moment components
			MBDyn_UVLM_CouplingDynamic_pre[3 * i + 6 * MBDyn_UVLM_NodesNum] = MBDyn_UVLM_CouplingDynamic[3 * i + 6 * MBDyn_UVLM_NodesNum];
			MBDyn_UVLM_CouplingDynamic_pre[3 * i + 1 + 6 * MBDyn_UVLM_NodesNum] = MBDyn_UVLM_CouplingDynamic[3 * i + 1 + 6 * MBDyn_UVLM_NodesNum];
			MBDyn_UVLM_CouplingDynamic_pre[3 * i + 2 + 6 * MBDyn_UVLM_NodesNum] = MBDyn_UVLM_CouplingDynamic[3 * i + 2 + 6 * MBDyn_UVLM_NodesNum];
		}
	}
	return WorkVec;
}

void
UvlmInterfaceBaseElem::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
UvlmInterfaceBaseElem::InitialAssJac(
	VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr) {

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
	const VectorHandler& XCurr) {

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
			mbdyn_uvlm_x = x + mbdyn_uvlm_f
			mbdyn_uvlm_R = R
			mbdyn_uvlm_v = xp + mbdyn_uvlm_w cross mbdyn_uvlm_f
			mbdyn_uvlm_w = w
			mbdyn_uvlm_a = xpp + mbdyn_uvlm_wp cross mbdyn_uvlm_f + mbdyn_uvlm_w cross mbdyn_uvlm_w cross mbdyn_uvlm_f
			mbdyn_uvlm_wp = wp
	*/

	for (unsigned i = 0; i < MBDyn_UVLM_NodesNum; ++i) {
		const MBDYN_UVLM_POINTDATA& mbdyn_uvlm_point = MBDyn_UVLM_Nodes[i];
		// Rotation and position
		const Mat3x3& mbdyn_uvlm_R = mbdyn_uvlm_point.pMBDyn_UVLM_Node->GetRPrev();
		const Mat3x3& mbdyn_uvlm_Rh = MBDyn_UVLM_Nodes[i].MBDyn_UVLM_RhM;   // Relative orientation of the marker in local ref.
		Mat3x3 mbdyn_uvlm_R_marker = mbdyn_uvlm_R * mbdyn_uvlm_Rh;
		Vec3 mbdyn_uvlm_f = mbdyn_uvlm_R * mbdyn_uvlm_point.MBDyn_UVLM_Offset;
		Vec3 mbdyn_uvlm_x = mbdyn_uvlm_point.pMBDyn_UVLM_Node->GetXPrev() + mbdyn_uvlm_f;
		// Angular velocity and velocity
		const Vec3& mbdyn_uvlm_w = mbdyn_uvlm_point.pMBDyn_UVLM_Node->GetWPrev();
		Vec3 mbdyn_uvlm_wCrossf = mbdyn_uvlm_w.Cross(mbdyn_uvlm_f);
		Vec3 mbdyn_uvlm_v = mbdyn_uvlm_point.pMBDyn_UVLM_Node->GetVPrev() + mbdyn_uvlm_wCrossf;
		// Angular accelaration and accelaration
		const Vec3 &mbdyn_uvlm_wp = mbdyn_uvlm_point.pMBDyn_UVLM_Node->GetWPPrev();
		Vec3 mbdyn_uvlm_a = mbdyn_uvlm_point.pMBDyn_UVLM_Node->GetXPPPrev() + mbdyn_uvlm_wp.Cross(mbdyn_uvlm_f) +
			mbdyn_uvlm_w.Cross(mbdyn_uvlm_wCrossf);

		std::vector<double> mbdyn_uvlm_tempvec3_x(3);
		std::vector<double> mbdyn_uvlm_tempvec3_v(3);
		std::vector<double> mbdyn_uvlm_tempvec3_a(3);
		std::vector<double> mbdyn_uvlm_tempvec3_w(3);
		std::vector<double> mbdyn_uvlm_tempvec3_wp(3);
		MBDyn_UVLM_Vec3D(mbdyn_uvlm_x, mbdyn_uvlm_tempvec3_x, MBDyn_UVLM_Scale[0]);
		MBDyn_UVLM_Vec3D(mbdyn_uvlm_v, mbdyn_uvlm_tempvec3_v, MBDyn_UVLM_Scale[0]);
		MBDyn_UVLM_Vec3D(mbdyn_uvlm_a, mbdyn_uvlm_tempvec3_a, MBDyn_UVLM_Scale[0]);
		MBDyn_UVLM_Vec3D(mbdyn_uvlm_w, mbdyn_uvlm_tempvec3_w, 1.0);
		MBDyn_UVLM_Vec3D(mbdyn_uvlm_wp, mbdyn_uvlm_tempvec3_wp, 1.0);
		std::vector<double> mbdyn_uvlm_tempmat3x3_R(9);
		MBDyn_UVLM_Mat3x3D(mbdyn_uvlm_R_marker, mbdyn_uvlm_tempmat3x3_R);


		MBDyn_UVLM_CouplingKinematic[3 * i] = mbdyn_uvlm_tempvec3_x[0];
		MBDyn_UVLM_CouplingKinematic[3 * i + 1] = mbdyn_uvlm_tempvec3_x[1];
		MBDyn_UVLM_CouplingKinematic[3 * i + 2] = mbdyn_uvlm_tempvec3_x[2];

		MBDyn_UVLM_CouplingKinematic[9 * i + 3 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempmat3x3_R[0];
		MBDyn_UVLM_CouplingKinematic[9 * i + 1 + 3 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempmat3x3_R[1];
		MBDyn_UVLM_CouplingKinematic[9 * i + 2 + 3 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempmat3x3_R[2];
		MBDyn_UVLM_CouplingKinematic[9 * i + 3 + 3 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempmat3x3_R[3];
		MBDyn_UVLM_CouplingKinematic[9 * i + 4 + 3 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempmat3x3_R[4];
		MBDyn_UVLM_CouplingKinematic[9 * i + 5 + 3 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempmat3x3_R[5];
		MBDyn_UVLM_CouplingKinematic[9 * i + 6 + 3 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempmat3x3_R[6];
		MBDyn_UVLM_CouplingKinematic[9 * i + 7 + 3 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempmat3x3_R[7];
		MBDyn_UVLM_CouplingKinematic[9 * i + 8 + 3 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempmat3x3_R[8];

		MBDyn_UVLM_CouplingKinematic[3 * i + 12 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempvec3_v[0];
		MBDyn_UVLM_CouplingKinematic[3 * i + 1 + 12 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempvec3_v[1];
		MBDyn_UVLM_CouplingKinematic[3 * i + 2 + 12 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempvec3_v[2];

		MBDyn_UVLM_CouplingKinematic[3 * i + 15 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempvec3_w[0];
		MBDyn_UVLM_CouplingKinematic[3 * i + 1 + 15 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempvec3_w[1];
		MBDyn_UVLM_CouplingKinematic[3 * i + 2 + 15 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempvec3_w[2];

		MBDyn_UVLM_CouplingKinematic[3 * i + 18 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempvec3_a[0];
		MBDyn_UVLM_CouplingKinematic[3 * i + 1 + 18 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempvec3_a[1];
		MBDyn_UVLM_CouplingKinematic[3 * i + 2 + 18 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempvec3_a[2];

		MBDyn_UVLM_CouplingKinematic[3 * i + 21 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempvec3_wp[0];
		MBDyn_UVLM_CouplingKinematic[3 * i + 1 + 21 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempvec3_wp[1];
		MBDyn_UVLM_CouplingKinematic[3 * i + 2 + 21 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempvec3_wp[2];

	}
}

void
UvlmInterfaceBaseElem::MBDyn_UVLM_SendDataToBuf_Curr() {

	/* formulation for calculation
			mbdyn_uvlm_x = x + mbdyn_uvlm_f
			mbdyn_uvlm_R = R
			mbdyn_uvlm_v = xp + mbdyn_uvlm_w cross mbdyn_uvlm_f
			mbdyn_uvlm_w = w
			mbdyn_uvlm_a = xpp + mbdyn_uvlm_wp cross mbdyn_uvlm_f + mbdyn_uvlm_w cross mbdyn_uvlm_w cross mbdyn_uvlm_f
			mbdyn_uvlm_wp = wp
	*/
	for (unsigned i = 0; i < MBDyn_UVLM_NodesNum; ++i)
	{
		const MBDYN_UVLM_POINTDATA& mbdyn_uvlm_point = MBDyn_UVLM_Nodes[i];
		//- Rotation and position
		const Mat3x3 & mbdyn_uvlm_R = mbdyn_uvlm_point.pMBDyn_UVLM_Node->GetRCurr();
		const Mat3x3 & mbdyn_uvlm_Rh = MBDyn_UVLM_Nodes[i].MBDyn_UVLM_RhM;             //- relative orientation of the marker in local ref.
		Mat3x3 mbdyn_uvlm_R_marker = mbdyn_uvlm_R * mbdyn_uvlm_Rh; //- absolute orientation of the marker in local ref.
		Vec3 mbdyn_uvlm_f = mbdyn_uvlm_R * mbdyn_uvlm_point.MBDyn_UVLM_Offset;
		Vec3 mbdyn_uvlm_x = mbdyn_uvlm_point.pMBDyn_UVLM_Node->GetXCurr() + mbdyn_uvlm_f;
		//- Angular velocity and velocity
		const Vec3 &mbdyn_uvlm_w = mbdyn_uvlm_point.pMBDyn_UVLM_Node->GetWCurr();
		Vec3 mbdyn_uvlm_wCrossf = mbdyn_uvlm_w.Cross(mbdyn_uvlm_f);
		Vec3 mbdyn_uvlm_v = mbdyn_uvlm_point.pMBDyn_UVLM_Node->GetVCurr() + mbdyn_uvlm_wCrossf;
		//- Angular acceleration and acceleration
		const Vec3 &mbdyn_uvlm_wp = mbdyn_uvlm_point.pMBDyn_UVLM_Node->GetWPCurr();
		Vec3 mbdyn_uvlm_a = mbdyn_uvlm_point.pMBDyn_UVLM_Node->GetXPPCurr() + mbdyn_uvlm_wp.Cross(mbdyn_uvlm_f) + mbdyn_uvlm_w.Cross(mbdyn_uvlm_wCrossf);

		std::vector<double> mbdyn_uvlm_tempvec3_x(3);
		std::vector<double> mbdyn_uvlm_tempvec3_v(3);
		std::vector<double> mbdyn_uvlm_tempvec3_a(3);
		std::vector<double> mbdyn_uvlm_tempvec3_w(3);
		std::vector<double> mbdyn_uvlm_tempvec3_wp(3);
		MBDyn_UVLM_Vec3D(mbdyn_uvlm_x, mbdyn_uvlm_tempvec3_x, MBDyn_UVLM_Scale[0]);
		MBDyn_UVLM_Vec3D(mbdyn_uvlm_v, mbdyn_uvlm_tempvec3_v, MBDyn_UVLM_Scale[0]);
		MBDyn_UVLM_Vec3D(mbdyn_uvlm_a, mbdyn_uvlm_tempvec3_a, MBDyn_UVLM_Scale[0]);
		MBDyn_UVLM_Vec3D(mbdyn_uvlm_w, mbdyn_uvlm_tempvec3_w, 1.0);
		MBDyn_UVLM_Vec3D(mbdyn_uvlm_wp, mbdyn_uvlm_tempvec3_wp, 1.0);
		std::vector<double> mbdyn_uvlm_tempmat3x3_R(9);
		MBDyn_UVLM_Mat3x3D(mbdyn_uvlm_R_marker, mbdyn_uvlm_tempmat3x3_R);

		
		MBDyn_UVLM_CouplingKinematic[3 * i] = mbdyn_uvlm_tempvec3_x[0];
		MBDyn_UVLM_CouplingKinematic[3 * i + 1] = mbdyn_uvlm_tempvec3_x[1];
		MBDyn_UVLM_CouplingKinematic[3 * i + 2] = mbdyn_uvlm_tempvec3_x[2];

		MBDyn_UVLM_CouplingKinematic[9 * i + 3 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempmat3x3_R[0];
		MBDyn_UVLM_CouplingKinematic[9 * i + 1 + 3 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempmat3x3_R[1];
		MBDyn_UVLM_CouplingKinematic[9 * i + 2 + 3 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempmat3x3_R[2];
		MBDyn_UVLM_CouplingKinematic[9 * i + 3 + 3 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempmat3x3_R[3];
		MBDyn_UVLM_CouplingKinematic[9 * i + 4 + 3 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempmat3x3_R[4];
		MBDyn_UVLM_CouplingKinematic[9 * i + 5 + 3 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempmat3x3_R[5];
		MBDyn_UVLM_CouplingKinematic[9 * i + 6 + 3 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempmat3x3_R[6];
		MBDyn_UVLM_CouplingKinematic[9 * i + 7 + 3 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempmat3x3_R[7];
		MBDyn_UVLM_CouplingKinematic[9 * i + 8 + 3 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempmat3x3_R[8];

		MBDyn_UVLM_CouplingKinematic[3 * i + 12 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempvec3_v[0];
		MBDyn_UVLM_CouplingKinematic[3 * i + 1 + 12 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempvec3_v[1];
		MBDyn_UVLM_CouplingKinematic[3 * i + 2 + 12 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempvec3_v[2];

		MBDyn_UVLM_CouplingKinematic[3 * i + 15 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempvec3_w[0];
		MBDyn_UVLM_CouplingKinematic[3 * i + 1 + 15 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempvec3_w[1];
		MBDyn_UVLM_CouplingKinematic[3 * i + 2 + 15 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempvec3_w[2];

		MBDyn_UVLM_CouplingKinematic[3 * i + 18 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempvec3_a[0];
		MBDyn_UVLM_CouplingKinematic[3 * i + 1 + 18 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempvec3_a[1];
		MBDyn_UVLM_CouplingKinematic[3 * i + 2 + 18 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempvec3_a[2];

		MBDyn_UVLM_CouplingKinematic[3 * i + 21 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempvec3_wp[0];
		MBDyn_UVLM_CouplingKinematic[3 * i + 1 + 21 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempvec3_wp[1];
		MBDyn_UVLM_CouplingKinematic[3 * i + 2 + 21 * MBDyn_UVLM_NodesNum] = mbdyn_uvlm_tempvec3_wp[2];
	}
}

//- intepolation of pos, vel, and acc. (used in Jacobian/Gauss scheme) to do
void
UvlmInterfaceBaseElem::MBDyn_UVLM_KinematicData_Interpolate()
{
	NO_OP;
}

//- a regular UVLM integration step.
void
UvlmInterfaceBaseElem::MBDyn_UVLM_UpdateUVLMModel() {

	//---------- A regular step to update UVLM model (this function executate steps 3 and 4)
	//---------- 3. UVLM models read the coupling data from buffer;
	//---------- 4. UVLM models do integration (one step);

	//---------- 3. UVLM models read the coupling data from buffer;
	time_step = m_PDM->pGetDrvHdl()->dGetTimeStep();  //- get time step
	MBDyn_UVLM_Model_RecvFromBuf(MBDyn_UVLM_Aerogrid, MBDyn_UVLM_CouplingKinematic, MBDyn_UVLM_NodesNum);

	//---------- 4. UVLM models do integration (one step);
	MBDyn_UVLM_Model_DoStepDynamics(MBDyn_UVLM_StepUVLM_settings, MBDyn_UVLM_Aerogrid_settings, MBDyn_UVLM_StraightWake_settings,
		MBDyn_UVLM_UVMopts, MBDyn_UVLM_FlightConditions, MBDyn_UVLM_Beam_inputs, MBDyn_UVLM_Aero_inputs,
		MBDyn_UVLM_StraightWake, MBDyn_UVLM_SteadyVelocityField, MBDyn_UVLM_Aerogrid,
		MBDyn_UVLM_UvlmLibVar, UVMoptions, FlightConditions, time_step);

	pedantic_cout("\t MBDyn_UVLM::Update a regular step\n");
}

//- Read the data from buffer, and write them to the Vec3 MBDyn_UVLM_F, Vec3 MBDyn_UVLM_DF and Vec3 MBDyn_UVLM_M;
void
UvlmInterfaceBaseElem::MBDyn_UVLM_RecvDataFromBuf() {

	//---------- C::E models sends data to the buffer;
	/*
	!
	!
	!
	!
	!
	!
	*/
	for (unsigned i = 0; i < MBDyn_UVLM_NodesNum; ++i) {

		//- read from the buffer
		MBDYN_UVLM_POINTDATA& mbdyn_uvlm_point = MBDyn_UVLM_Nodes[i];
		mbdyn_uvlm_point.MBDyn_UVLM_F = Vec3(MBDyn_UVLM_CouplingDynamic[3 * i], MBDyn_UVLM_CouplingDynamic[3 * i + 1], MBDyn_UVLM_CouplingDynamic[3 * i + 2]);
		mbdyn_uvlm_point.MBDyn_UVLM_DF = Vec3(MBDyn_UVLM_CouplingDynamic[3 * i + 3 * MBDyn_UVLM_NodesNum], MBDyn_UVLM_CouplingDynamic[3 * i + 1 + 3 * MBDyn_UVLM_NodesNum], MBDyn_UVLM_CouplingDynamic[3 * i + 2 + 3 * MBDyn_UVLM_NodesNum]);
		mbdyn_uvlm_point.MBDyn_UVLM_M = Vec3(MBDyn_UVLM_CouplingDynamic[3 * i + 6 * MBDyn_UVLM_NodesNum], MBDyn_UVLM_CouplingDynamic[3 * i + 1 + 6 * MBDyn_UVLM_NodesNum], MBDyn_UVLM_CouplingDynamic[3 * i + 2 + 6 * MBDyn_UVLM_NodesNum]);
	}
}

void
UvlmInterfaceBaseElem::Output(OutputHandler& OH) const {

	pedantic_cout("\t MBDyn::Output() \n");
	//- Output the motion of UVLM bodies in LOADABLE files;
	if (!OH.IsOpen(OutputHandler::OutFiles::LOADABLE)) {
		OH.Open(OutputHandler::OutFiles::LOADABLE);
	}
	std::ostream &mbdyn_uvlm_body_output = OH.Loadable();
	/*
	!
	!
	!
	!
	!
	!
	*/
	if (MBDyn_UVLM_out_forces) {
		for (unsigned i = 0; i < MBDyn_UVLM_NodesNum; ++i) {
			const MBDYN_UVLM_POINTDATA &point = MBDyn_UVLM_Nodes[i];
			const MBDYN_UVLM_MODELDATA &uvlm_body = MBDyn_UVLM_Model_Label[i];
			Vec3 temp_moment = point.MBDyn_UVLM_M + (point.pMBDyn_UVLM_Node->GetRCurr() * point.MBDyn_UVLM_Offset).Cross(point.MBDyn_UVLM_F);
			if (uvlm_body.bMBDyn_UVLM_Body_Output || MBDyn_UVLM_OutputType == MBDYN_UVLM_OUTPUTTYPE::MBDYN_UVLM_OUTPUT_ALLBODIES) {

				MBDyn_UVLM_out_forces << std::setw(8) << point.MBDyn_UVLM_uLabel
					<< std::setw(16) << MBDyn_UVLM_CouplingDynamic[3 * i]
					<< std::setw(16) << MBDyn_UVLM_CouplingDynamic[3 * i + 1]
					<< std::setw(16) << MBDyn_UVLM_CouplingDynamic[3 * i + 2]
					<< std::setw(16) << MBDyn_UVLM_CouplingDynamic[3 * i + 3 * MBDyn_UVLM_NodesNum]         
					<< std::setw(16) << MBDyn_UVLM_CouplingDynamic[3 * i + 1 + 3 * MBDyn_UVLM_NodesNum]     
					<< std::setw(16) << MBDyn_UVLM_CouplingDynamic[3 * i + 2 + 3 * MBDyn_UVLM_NodesNum]     
					<< std::setw(16) << temp_moment.pGetVec()[0]                                            // moment in the res.
					<< std::setw(16) << temp_moment.pGetVec()[1]                                            // moment in the res.
					<< std::setw(16) << temp_moment.pGetVec()[2]                                            // moment in the res.
					<< std::setw(16) << MBDyn_UVLM_CouplingDynamic[3 * i + 6 * MBDyn_UVLM_NodesNum]         // received moment.
					<< std::setw(16) << MBDyn_UVLM_CouplingDynamic[3 * i + 1 + 6 * MBDyn_UVLM_NodesNum]     // received moment.
					<< std::setw(16) << MBDyn_UVLM_CouplingDynamic[3 * i + 2 + 6 * MBDyn_UVLM_NodesNum]     // received moment.
					<< std::endl;
			}
		}
	}
}

std::ostream&
UvlmInterfaceBaseElem::Restart(std::ostream& out) const
{
	return out << "# UvlmInterface: is doing now" << std::endl;
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

unsigned int
UvlmInterfaceBaseElem::iGetInitialNumDof(void) const
{
	return 0;
}


// private functions: [START]

void
UvlmInterfaceBaseElem::MBDyn_UVLM_Vec3D(const Vec3& mbdyn_uvlm_Vec3, std::vector<double>& mbdyn_uvlm_temp,
	double MBDyn_UVLM_LengthScale) const {

	mbdyn_uvlm_temp[0] = MBDyn_UVLM_LengthScale * static_cast<double>(*(mbdyn_uvlm_Vec3.pGetVec()));
	mbdyn_uvlm_temp[1] = MBDyn_UVLM_LengthScale * static_cast<double>(*(mbdyn_uvlm_Vec3.pGetVec() + 1));
	mbdyn_uvlm_temp[2] = MBDyn_UVLM_LengthScale * static_cast<double>(*(mbdyn_uvlm_Vec3.pGetVec() + 2));
}

void
UvlmInterfaceBaseElem::MBDyn_UVLM_Mat3x3D(const Mat3x3& mbdyn_uvlm_Mat3x3, std::vector<double>& mbdyn_uvlm_temp) const {

	for (unsigned i = 0; i < 9; i++) {
		mbdyn_uvlm_temp[i] = static_cast<double>(mbdyn_uvlm_Mat3x3.pGetMat()[i]);
	}
}

double
UvlmInterfaceBaseElem::MBDyn_UVLM_CalculateError() {    //- calculate the error of coupling force.

	double mbdyn_uvlm_temp_error = 0.0;  //- using Euclidean norm
	for (unsigned i = 0; i < MBDyn_UVLM_NodesNum; ++i) {
		for (unsigned j = 0; j < 3; ++j) {
			mbdyn_uvlm_temp_error += pow(MBDyn_UVLM_CouplingDynamic_pre[3 * i + j] - MBDyn_UVLM_CouplingDynamic[3 * i + j], 2);
			mbdyn_uvlm_temp_error += pow(MBDyn_UVLM_CouplingDynamic_pre[3 * i + j + 3 * MBDyn_UVLM_NodesNum] - MBDyn_UVLM_CouplingDynamic[3 * i + j + 3 * MBDyn_UVLM_NodesNum], 2);
			mbdyn_uvlm_temp_error += pow(MBDyn_UVLM_CouplingDynamic_pre[3 * i + j + 6 * MBDyn_UVLM_NodesNum] - MBDyn_UVLM_CouplingDynamic[3 * i + j + 6 * MBDyn_UVLM_NodesNum], 2);
		}
	}
	return sqrt(mbdyn_uvlm_temp_error);
}

void
UvlmInterfaceBaseElem::MBDyn_UVLM_MBDynPrint() const {

	//- print the data calculated in MBDyn
	for (unsigned i = 0; i < MBDyn_UVLM_NodesNum; ++i) {

		double time;
		time = m_PDM->dGetTime();
		const MBDYN_UVLM_POINTDATA &mbdyn_uvlm_point = MBDyn_UVLM_Nodes[i];
		std::cout << "time in MBDyn: " << time << "\n\t"
			<< "coupling node" << mbdyn_uvlm_point.MBDyn_UVLM_uLabel << "\n"
			<< "\t\t forces obtained: " << MBDyn_UVLM_CouplingDynamic[3 * i] << "\t" << MBDyn_UVLM_CouplingDynamic[3 * i + 1] << "\t" << MBDyn_UVLM_CouplingDynamic[3 * i + 2]
			<< "\n"
			<< "\t \t torques obtained: " << MBDyn_UVLM_CouplingDynamic[3 * i + 6 * MBDyn_UVLM_NodesNum] << "\t" << MBDyn_UVLM_CouplingDynamic[3 * i + 1 + 6 * MBDyn_UVLM_NodesNum] << "\t" << MBDyn_UVLM_CouplingDynamic[3 * i + 2 + 6 * MBDyn_UVLM_NodesNum]
			<< "\n";

		//- current states
		std::vector<double> mbdyn_uvlm_tempvec3_x(3);
		std::vector<double> mbdyn_uvlm_tempvec3_v(3);
		std::vector<double> mbdyn_uvlm_tempvec3_a(3);
		std::vector<double> mbdyn_uvlm_tempvec3_w(3);
		std::vector<double> mbdyn_uvlm_tempvec3_wp(3);
		std::vector<double> mbdyn_uvlm_tempmat3x3_R(9);

		//- states of last step
		std::vector<double> mbdyn_uvlm_tempvec3_xPrev(3);
		std::vector<double> mbdyn_uvlm_tempvec3_vPrev(3);
		std::vector<double> mbdyn_uvlm_tempvec3_aPrev(3);
		std::vector<double> mbdyn_uvlm_tempvec3_wPrev(3);
		std::vector<double> mbdyn_uvlm_tempvec3_wpPrev(3);
		std::vector<double> mbdyn_uvlm_tempmat3x3_RPrev(9);
		{
			//- rotation and position
			const Mat3x3 &mbdyn_uvlm_R = mbdyn_uvlm_point.pMBDyn_UVLM_Node->GetRCurr();
			Vec3 mbdyn_uvlm_f = mbdyn_uvlm_R * mbdyn_uvlm_point.MBDyn_UVLM_Offset;
			Vec3 mbdyn_uvlm_x = mbdyn_uvlm_point.pMBDyn_UVLM_Node->GetXCurr() + mbdyn_uvlm_f;

			//- angular velocity and velocity
			const Vec3 &mbdyn_uvlm_w = mbdyn_uvlm_point.pMBDyn_UVLM_Node->GetWCurr();
			Vec3 mbdyn_uvlm_wCrossf = mbdyn_uvlm_w.Cross(mbdyn_uvlm_f);
			Vec3 mbdyn_uvlm_v = mbdyn_uvlm_point.pMBDyn_UVLM_Node->GetVCurr() + mbdyn_uvlm_wCrossf;

			//- angular acceleration and acceleration
			const Vec3 &mbdyn_uvlm_wp = mbdyn_uvlm_point.pMBDyn_UVLM_Node->GetWPCurr();
			Vec3 mbdyn_uvlm_a = mbdyn_uvlm_point.pMBDyn_UVLM_Node->GetXPPCurr() + mbdyn_uvlm_wp.Cross(mbdyn_uvlm_f) + mbdyn_uvlm_w.Cross(mbdyn_uvlm_wCrossf);
		
			MBDyn_UVLM_Vec3D(mbdyn_uvlm_x, mbdyn_uvlm_tempvec3_x, MBDyn_UVLM_Scale[0]);
			MBDyn_UVLM_Vec3D(mbdyn_uvlm_v, mbdyn_uvlm_tempvec3_v, MBDyn_UVLM_Scale[0]);
			MBDyn_UVLM_Vec3D(mbdyn_uvlm_a, mbdyn_uvlm_tempvec3_a, MBDyn_UVLM_Scale[0]);
			MBDyn_UVLM_Vec3D(mbdyn_uvlm_w, mbdyn_uvlm_tempvec3_w, 1.0);
			MBDyn_UVLM_Vec3D(mbdyn_uvlm_wp, mbdyn_uvlm_tempvec3_wp, 1.0);
			MBDyn_UVLM_Mat3x3D(mbdyn_uvlm_R, mbdyn_uvlm_tempmat3x3_R);
			std::cout << "\t\tpos: " << mbdyn_uvlm_tempvec3_x[0] << "\t" << mbdyn_uvlm_tempvec3_x[1] << "\t" << mbdyn_uvlm_tempvec3_x[2] << "\n"
				<< "\t\tvel: " << mbdyn_uvlm_tempvec3_v[0] << "\t" << mbdyn_uvlm_tempvec3_v[1] << "\t" << mbdyn_uvlm_tempvec3_v[2] << "\n"
				<< "\t\tacc: " << mbdyn_uvlm_tempvec3_a[0] << "\t" << mbdyn_uvlm_tempvec3_a[1] << "\t" << mbdyn_uvlm_tempvec3_a[2] << "\n";
		}

		{
			//- data of last step
			//- rotation and position
			const Mat3x3 &mbdyn_uvlm_R = mbdyn_uvlm_point.pMBDyn_UVLM_Node->GetRPrev();
			Vec3 mbdyn_uvlm_f = mbdyn_uvlm_R * mbdyn_uvlm_point.MBDyn_UVLM_Offset;
			Vec3 mbdyn_uvlm_x = mbdyn_uvlm_point.pMBDyn_UVLM_Node->GetXPrev() + mbdyn_uvlm_f;
			
			//- angular velocity and velocity
			const Vec3 &mbdyn_uvlm_w = mbdyn_uvlm_point.pMBDyn_UVLM_Node->GetWPrev();
			Vec3 mbdyn_uvlm_wCrossf = mbdyn_uvlm_w.Cross(mbdyn_uvlm_f);
			Vec3 mbdyn_uvlm_v = mbdyn_uvlm_point.pMBDyn_UVLM_Node->GetVPrev() + mbdyn_uvlm_wCrossf;

			//- angular acceleration and acceleration
			const Vec3 &mbdyn_uvlm_wp = mbdyn_uvlm_point.pMBDyn_UVLM_Node->GetWPPrev();
			Vec3 mbdyn_uvlm_a = mbdyn_uvlm_point.pMBDyn_UVLM_Node->GetXPPPrev() + mbdyn_uvlm_wp.Cross(mbdyn_uvlm_f) + mbdyn_uvlm_w.Cross(mbdyn_uvlm_wCrossf);

			MBDyn_UVLM_Vec3D(mbdyn_uvlm_x, mbdyn_uvlm_tempvec3_xPrev, MBDyn_UVLM_Scale[0]);
			MBDyn_UVLM_Vec3D(mbdyn_uvlm_v, mbdyn_uvlm_tempvec3_vPrev, MBDyn_UVLM_Scale[0]);
			MBDyn_UVLM_Vec3D(mbdyn_uvlm_a, mbdyn_uvlm_tempvec3_aPrev, MBDyn_UVLM_Scale[0]);
			MBDyn_UVLM_Vec3D(mbdyn_uvlm_w, mbdyn_uvlm_tempvec3_wPrev, 1.0);
			MBDyn_UVLM_Vec3D(mbdyn_uvlm_wp, mbdyn_uvlm_tempvec3_wpPrev, 1.0);
			MBDyn_UVLM_Mat3x3D(mbdyn_uvlm_R, mbdyn_uvlm_tempmat3x3_RPrev);
			std::cout << "\t\t pos_prev: " << mbdyn_uvlm_tempvec3_xPrev[0] << "\t" << mbdyn_uvlm_tempvec3_xPrev[1] << "\t" << mbdyn_uvlm_tempvec3_xPrev[2] << "\n"
				<< "\t\t vel_prev: " << mbdyn_uvlm_tempvec3_vPrev[0] << "\t" << mbdyn_uvlm_tempvec3_vPrev[1] << "\t" << mbdyn_uvlm_tempvec3_vPrev[2] << "\n"
				<< "\t\t acc_prev: " << mbdyn_uvlm_tempvec3_aPrev[0] << "\t" << mbdyn_uvlm_tempvec3_aPrev[1] << "\t" << mbdyn_uvlm_tempvec3_aPrev[2] << "\n";
		}

		//- predicted states
		double mbdyn_uvlm_temp_vPred;
		double e1 = 12000.0, e2 = 12000.0, f1 = 8.0, f2 = 5.0;
		mbdyn_uvlm_temp_vPred = e1 * mbdyn_uvlm_tempvec3_x[0] + e2 * mbdyn_uvlm_tempvec3_xPrev[0] + f1 * mbdyn_uvlm_tempvec3_v[0] + f2 * mbdyn_uvlm_tempvec3_vPrev[0];
		std::cout << "\t\t vel_pred: " << mbdyn_uvlm_temp_vPred << "\t"
			<< "x-xPrev: " << mbdyn_uvlm_tempvec3_x[0] - mbdyn_uvlm_tempvec3_xPrev[0] << "\n";
	}
}

// private functions: [END]


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

