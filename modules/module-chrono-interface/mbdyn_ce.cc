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
// standard library
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include <cmath>
#include <unistd.h>
#include <sys/time.h> //for both Unix/Windowns <sys/time.h> only for Unix

#include "mbdyn_ce.h"
#include "chrono/ChConfig.h"
#include "chrono_parallel/ChDataManager.h" // for simulation of parallel system, data_manager
#include "chrono_parallel/solver/ChIterativeSolverParallel.h"
#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono/physics/ChLinkMotionImposed.h" //for 3-D dimension
#include "chrono/motion_functions/ChFunctionPosition_line.h"
#include "chrono/motion_functions/ChFunctionPosition_setpoint.h"
#include "chrono/motion_functions/ChFunctionRotation_setpoint.h"
#include "chrono/motion_functions/ChFunctionRotation_spline.h"
#include "chrono/geometry/ChLineSegment.h"
#include "chrono_thirdparty/filesystem/path.h"

#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

using namespace chrono;
using namespace chrono::collision;

clock_t startTime,endTime;
struct timeval start_time, end_time;

extern "C" void
MBDyn_CE_CEModel_Create(ChSystemParallelNSC *pMBDyn_CE_CEModel);

// extern "C" void
// MBDyn_CE_CEModel_Create(ChSystem *pMBDyn_CE_CEModel);

extern "C" pMBDyn_CE_CEModel_t MBDyn_CE_CEModel_Init
(std::vector<double> & MBDyn_CE_CEModel_Data, //- all body infor in C::E Model for saving data and reloading data to restart C::E step (used in tight co-simulation)
const int MBDyn_CE_CEMotorType, //- Motor type: by position, velocity, or acceleration...
const double* pMBDyn_CE_CEFrame, const double* MBDyn_CE_CEScale,
std::vector<MBDYN_CE_CEMODELDATA> & MBDyn_CE_CEModel_Label, //- IDs of coupling bodies and motors in C::E model, and output cmd
const int& MBDyn_CE_CouplingType) //- Coupling type
{
	gettimeofday(&start_time, NULL);//get the start time
	startTime = clock();
	std::cout << "Initial MBDyn_CE_CEModel pointer:\n";
	ChSystemParallelNSC *pMBDyn_CE_CEModel = new ChSystemParallelNSC; //- By now, only support ChSystemParallel.
	//--------------- create C::E model (defined by user)
	MBDyn_CE_CEModel_Create(pMBDyn_CE_CEModel);
	if(pMBDyn_CE_CEModel==NULL)
	{
		std::cout << "\t\tInitial MBDyn_CE_CEModel pointer fails\n";
		//- MBDyn_CE_CEModel_Destroy(pMBDyn_CE_CEModel);
		pMBDyn_CE_CEModel = NULL;
		return NULL;
	}
	else
	{
		//--------------- initialize information
		//---------- initial the frame MBDyn in C::E;
		//- r = (-CEF1,-CEF2,-CEF3);
		//- R= [CEF3,  CEF6,  CEF9 ]^T           [CEF3,   CEF4,   CEF5]
		//-    [CEF4,  CEF7,  CEF10]     ====    [CEF6,   CEF7,   CEF8]
		//-    [CEF5,  CEF8,  CEF11];            [CEF9,   CEF10,  CEF11]
		ChVector<> mbdynce_temp_frameMBDyn_pos(-pMBDyn_CE_CEFrame[0], -pMBDyn_CE_CEFrame[1], -pMBDyn_CE_CEFrame[2]);
		ChVector<> mbdynce_temp_frameMBDyn_rot_axis_X(pMBDyn_CE_CEFrame[3], pMBDyn_CE_CEFrame[6], pMBDyn_CE_CEFrame[9]);
		ChVector<> mbdynce_temp_frameMBDyn_rot_axis_Y(pMBDyn_CE_CEFrame[4], pMBDyn_CE_CEFrame[7], pMBDyn_CE_CEFrame[10]);
		ChVector<> mbdynce_temp_frameMBDyn_rot_axis_Z(pMBDyn_CE_CEFrame[5], pMBDyn_CE_CEFrame[8], pMBDyn_CE_CEFrame[11]);
		ChMatrix33<> mbdynce_temp_frameMBDyn_rot(mbdynce_temp_frameMBDyn_rot_axis_X, mbdynce_temp_frameMBDyn_rot_axis_Y, mbdynce_temp_frameMBDyn_rot_axis_Z);
		ChFrame<> mbdynce_temp_frameMBDyn(mbdynce_temp_frameMBDyn_pos, mbdynce_temp_frameMBDyn_rot);
		//---------- initial the gravity of C::E model
		pMBDyn_CE_CEModel->Set_G_acc(mbdynce_temp_frameMBDyn.TransformDirectionLocalToParent(MBDyn_CE_CEScale[0]*pMBDyn_CE_CEModel->Get_G_acc()));
		std::cout << pMBDyn_CE_CEModel->Get_G_acc()<<"\n";
		//---------- initial motor for coupling bodies.
		if (MBDyn_CE_CouplingType >= -1) //- coupling
		{
			unsigned mbdynce_temp_bodies_num = MBDyn_CE_CEModel_Label.size() - 1; //- num of coupling bodies
			unsigned mbdynce_temp_ground_id = MBDyn_CE_CEModel_Label[mbdynce_temp_bodies_num].MBDyn_CE_CEBody_Label; //- ground ID is set in the last element
			MBDyn_CE_CEModel_Label[mbdynce_temp_bodies_num].MBDyn_CE_CEMotor_Label = 0;								 //- the last label is invalid
			auto mbdynce_temp_ground = pMBDyn_CE_CEModel->SearchBodyID(mbdynce_temp_ground_id);
			if (mbdynce_temp_ground == NULL)
			{
				std::cout << "\t\tCannot find C::E ground\n";
				delete pMBDyn_CE_CEModel;
				pMBDyn_CE_CEModel = NULL;
				return NULL;
			}
			//----- create motor for each coupling body
			for (unsigned i = 0; i < mbdynce_temp_bodies_num; i++)
			{
				unsigned body_i_id = MBDyn_CE_CEModel_Label[i].MBDyn_CE_CEBody_Label;
				auto body_i = pMBDyn_CE_CEModel->SearchBodyID(body_i_id);
				if (body_i == NULL)
				{
					std::cout << "\t\tCannot find C::E body " << body_i_id << "\n";
					delete pMBDyn_CE_CEModel;
					pMBDyn_CE_CEModel = NULL;
					return NULL;
				}
				auto motor3d_body_i = std::make_shared<ChLinkMotionImposed>();
				//--------------- initialize offset
				//---------- Marker in MBDyn >>>> Marker in C::E;
				//- r = (off1, off2,  off3); (offset in C::E bodies' local ref)
				//- Rh_MBDynG: (joint orientation in MBDyn frame)
				//-    [CEF0,  CEF3,  CEF6]
				//-    [CEF1,  CEF4,  CEF7]
				//-    [CEF2,  CEF5,  CEF8];
				ChVector<> temp_mbdyn_CE_offset(MBDyn_CE_CEModel_Label[i].MBDyn_CE_CEBody_Offset[0],
												MBDyn_CE_CEModel_Label[i].MBDyn_CE_CEBody_Offset[1],
												MBDyn_CE_CEModel_Label[i].MBDyn_CE_CEBody_Offset[2]);
				ChVector<> temp_mbdyn_CE_Rh_x(MBDyn_CE_CEModel_Label[i].MBDyn_CE_CEBody_Rh[0], MBDyn_CE_CEModel_Label[i].MBDyn_CE_CEBody_Rh[1], MBDyn_CE_CEModel_Label[i].MBDyn_CE_CEBody_Rh[2]);
				ChVector<> temp_mbdyn_CE_Rh_y(MBDyn_CE_CEModel_Label[i].MBDyn_CE_CEBody_Rh[3], MBDyn_CE_CEModel_Label[i].MBDyn_CE_CEBody_Rh[4], MBDyn_CE_CEModel_Label[i].MBDyn_CE_CEBody_Rh[5]);
				ChVector<> temp_mbdyn_CE_Rh_z(MBDyn_CE_CEModel_Label[i].MBDyn_CE_CEBody_Rh[6], MBDyn_CE_CEModel_Label[i].MBDyn_CE_CEBody_Rh[7], MBDyn_CE_CEModel_Label[i].MBDyn_CE_CEBody_Rh[8]);
				ChMatrix33<> temp_mbdyn_CE_Rh_MBDyn(temp_mbdyn_CE_Rh_x, temp_mbdyn_CE_Rh_y, temp_mbdyn_CE_Rh_z);
				ChQuaternion<> temp_mbdyn_CE_Rh_G(temp_mbdyn_CE_Rh_MBDyn.Get_A_quaternion()>> mbdynce_temp_frameMBDyn);
				ChQuaternion<> temp_mbdyn_CE_Rh_local(body_i->GetRot().GetConjugate()*temp_mbdyn_CE_Rh_G);
				std::cout << "temp_mbdyn_CE_Rh_G " << temp_mbdyn_CE_Rh_G << "\n";
				motor3d_body_i->Initialize(body_i,
										   mbdynce_temp_ground,
										   true,													//- connecting frames are described in local ref.
										   ChFrame<>(temp_mbdyn_CE_offset, temp_mbdyn_CE_Rh_local), //- By default: using the mass center and the body orientation
										   ChFrame<>(ChVector<>(0.0, 0.0, 0.0), ChQuaternion<>(1.0, 0.0, 0.0, 0.0)));
				pMBDyn_CE_CEModel->Add(motor3d_body_i);

				if (MBDyn_CE_CEMotorType == MBDyn_CE_CEMOTORTYPE::VELOCITY) //- coupling by velocity
				{
					// the setpoint function can be used for co-simulation by position, or velocity, or both(
					// if both position and velocity are imposed to C::E objects, there are small errors in the resulting 
					// position and velocity of coupled bodies.)					).
					// here, we use it to achieve co-simulation by velocity, which means that the resulting velocity of the
					// coupled bodies are totally equal to the inputting values.
					auto motor3d_function_pos = std::make_shared<ChFunctionPosition_setpoint>();;
					motor3d_function_pos->SetMode(ChFunctionPosition_setpoint::eChSetpointMode::OVERRIDE);
					auto motor3d_function_rot = std::make_shared<ChFunctionRotation_setpoint>();
					motor3d_function_rot->SetMode(ChFunctionRotation_setpoint::eChSetpointMode::OVERRIDE);
					motor3d_body_i->SetPositionFunction(motor3d_function_pos);
					motor3d_body_i->SetRotationFunction(motor3d_function_rot);
					std::cout<<"motor type is velocity (using function setpoint)\n";
				}
				else if (MBDyn_CE_CEMotorType == MBDyn_CE_CEMOTORTYPE::POSITION) //- coupling by position
				{
					//- cosim using spline/line interpolation for position
					auto motor3d_function_pos = std::make_shared<ChFunctionPosition_line>();
					auto motor3d_function_rot = std::make_shared<ChFunctionRotation_spline>();
					motor3d_body_i->SetPositionFunction(motor3d_function_pos);
					motor3d_body_i->SetRotationFunction(motor3d_function_rot);
					std::cout<<"motor type is position (using function spline)\n";
				}	
				else if (MBDyn_CE_CEMotorType == MBDyn_CE_CEMOTORTYPE::ACCELERATION) //- coupling by acceleration
				{
					//- the motor is the same as that for motortype::velocity
					auto motor3d_function_pos = std::make_shared<ChFunctionPosition_setpoint>();;
					motor3d_function_pos->SetMode(ChFunctionPosition_setpoint::eChSetpointMode::OVERRIDE);
					auto motor3d_function_rot = std::make_shared<ChFunctionRotation_setpoint>();
					motor3d_function_rot->SetMode(ChFunctionRotation_setpoint::eChSetpointMode::OVERRIDE);
					motor3d_body_i->SetPositionFunction(motor3d_function_pos);
					motor3d_body_i->SetRotationFunction(motor3d_function_rot);
					std::cout<<"motor type is acceleration (using function setpoint)\n";
				}
				else if (MBDyn_CE_CEMotorType == MBDyn_CE_CEMOTORTYPE::ORIGIN) //- coupling by origin formulation
				{
					//- the motor is the same as that for motortype::velocity
					auto motor3d_function_pos = std::make_shared<ChFunctionPosition_setpoint>();;
					motor3d_function_pos->SetMode(ChFunctionPosition_setpoint::eChSetpointMode::OVERRIDE);
					auto motor3d_function_rot = std::make_shared<ChFunctionRotation_setpoint>();
					motor3d_function_rot->SetMode(ChFunctionRotation_setpoint::eChSetpointMode::OVERRIDE);
					motor3d_body_i->SetPositionFunction(motor3d_function_pos);
					motor3d_body_i->SetRotationFunction(motor3d_function_rot);
					std::cout<<"motor type is origin (using function setpoint)\n";
				}
				bool *pmbdynce_temp_bConstraint = &MBDyn_CE_CEModel_Label[i].bMBDyn_CE_Constraint[0];
				motor3d_body_i->SetConstrainedCoords(pmbdynce_temp_bConstraint[0],pmbdynce_temp_bConstraint[1],pmbdynce_temp_bConstraint[2],pmbdynce_temp_bConstraint[3],pmbdynce_temp_bConstraint[4],pmbdynce_temp_bConstraint[5]);
				std::cout << "\tpostion constraints are:\t"
						  << pmbdynce_temp_bConstraint[0] << "\t" << pmbdynce_temp_bConstraint[1] << "\t" << pmbdynce_temp_bConstraint[2] << "\n"
						  << "\trotation constraints are:\t"
						  << pmbdynce_temp_bConstraint[3] << "\t" << pmbdynce_temp_bConstraint[4] << "\t" << pmbdynce_temp_bConstraint[5] << "\n";
				MBDyn_CE_CEModel_Label[i].MBDyn_CE_CEMotor_Label = motor3d_body_i->GetIdentifier();
				//----- output motor information
				if (true)
				{
					std::cout << "\tC::E motor\t" << i + 1 << "\tID:\t" << MBDyn_CE_CEModel_Label[i].MBDyn_CE_CEMotor_Label
							  << "\tbody 1:\t" << body_i->GetIdentifier() << "\tbody 2:\t" << mbdynce_temp_ground->GetIdentifier() << "\n"
							  << "\tFrame1_pos:\t" << motor3d_body_i->GetFrame1().GetPos() << "\n"
							  << "\tFrame1_rot:\t" << motor3d_body_i->GetFrame1().GetRot() << "\n"
							  << "\tRelative_matrix_mbdyn (Ground ref):\t" << temp_mbdyn_CE_Rh_MBDyn << "\n"
							  << "\tAbsolute_rot:\t" << temp_mbdyn_CE_Rh_G<< "\n"
							  << "\tRelative_rot:\t" << temp_mbdyn_CE_Rh_local << "\n";							  
				}
			}
		}
		else
		{
			std::cout << "Coupling none in C::E model.\n";
		}

		//---------- allocate space for C::E_Model_Data;
		unsigned int bodies_size, system_size;
		bodies_size = pMBDyn_CE_CEModel->Get_bodylist().size();
		switch (MBDyn_CE_CouplingType)
		{
			case MBDyn_CE_COUPLING::COUPLING_TIGHT:
			{
				system_size = (3 * 3 + 3 * 4) * bodies_size + 1; // +1: save Chtime; system_size=body_size+1(time)
				MBDyn_CE_CEModel_Data.resize(system_size, 0.0);
				break;
			}
			case MBDyn_CE_COUPLING::COUPLING_NONE: //- for these three cases, it
			case MBDyn_CE_COUPLING::COUPLING_STSTAGGERED:
			case MBDyn_CE_COUPLING::COUPLING_LOOSE:
			default:
			{
				system_size = 0; //- save nothing
				MBDyn_CE_CEModel_Data.resize(system_size, 0.0);
			}
		}

		//---------- print model information after adding motor
		std::cout << "there is the coupling C::E model \n";
		std::cout << "num of links:\t" << pMBDyn_CE_CEModel->Get_linklist().size() << "\n";
		std::cout << "num of other physicslist:\t" << pMBDyn_CE_CEModel->Get_otherphysicslist().size() << "\n";
		std::cout << "num of rigid bodies:\t" << pMBDyn_CE_CEModel->Get_bodylist().size() << "\n";
		std::cout << "num of speed motor:\t" << pMBDyn_CE_CEModel->data_manager->num_linmotors << "\n";
		//---------- print detailed information about the C::E model
		if(true)
		{
			std::cout << "Detailed information of the C::E model:\n";
			for (unsigned i = 0; i < bodies_size; i++)
			{
				std::cout << "Body " << pMBDyn_CE_CEModel->Get_bodylist()[i]->GetIdentifier() << "\n";
				std::cout << "\tpos: " << pMBDyn_CE_CEModel->Get_bodylist()[i]->GetPos() << "\n";
				std::cout << "\tpos_dt: " << pMBDyn_CE_CEModel->Get_bodylist()[i]->GetPos_dt() << "\n";
				std::cout << "\tpos_dtdt: " << pMBDyn_CE_CEModel->Get_bodylist()[i]->GetPos_dtdt() << "\n";
				std::cout << "\trot: " << pMBDyn_CE_CEModel->Get_bodylist()[i]->GetRot() << "\n";
				std::cout << "\tWvel_par: " << pMBDyn_CE_CEModel->Get_bodylist()[i]->GetWvel_par() << "\n";
				std::cout << "\tWacc_par: " << pMBDyn_CE_CEModel->Get_bodylist()[i]->GetWacc_par() << "\n";
			}
			for (unsigned i = 0; i < pMBDyn_CE_CEModel->Get_linklist().size(); i++)
			{
				std::cout << "Link " << pMBDyn_CE_CEModel->Get_linklist()[i]->GetIdentifier() << "\n";
				std::cout << "\tFrame2G_pos: " << pMBDyn_CE_CEModel->Get_linklist()[i]->GetLinkAbsoluteCoords().pos << "\n"; //- Frame2G
				std::cout << "\tFrame2G_rot: " << pMBDyn_CE_CEModel->Get_linklist()[i]->GetLinkAbsoluteCoords().rot << "\n";
			}
		}
#ifdef CHRONO_OPENGL
//#if 0
		opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
		gl_window.Initialize(1280, 720, "Pendulum", pMBDyn_CE_CEModel);
		gl_window.SetCamera(ChVector<>(-20.89465, 104.58118, 223.05379), 
							ChVector<>(-20.94352, 104.29030, 222.09828), ChVector<>(0, 1.0, 0), 5.0f);
		gl_window.SetRenderMode(opengl::WIREFRAME);
#endif
		return pMBDyn_CE_CEModel;
	}
}

extern "C" void
MBDyn_CE_CEModel_Destroy(pMBDyn_CE_CEModel_t &pMBDyn_CE_CEModel)
{
	gettimeofday(&end_time, NULL);//the end time
	double time_diff_sec = (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_usec - start_time.tv_usec) / 1.0e6;
	endTime = clock();

	std::cout << "time consuming:\t" << time_diff_sec << "s\t" << "CPU time:\t" << (double)(endTime - startTime) / CLOCKS_PER_SEC <<"s"
			 << "\n";
	if(pMBDyn_CE_CEModel!=NULL)
	{
		// must convert to the correct type
		// void * pointer is transfored to type ChsystemParallelNSC *.(downcast)
		// ChSystemParallelNSC *temp_pointer = dynamic_cast<ChSystemParallelNSC *> (pMBDyn_CE_CEModel);
		ChSystemParallelNSC *temp_pointer = (ChSystemParallelNSC *)pMBDyn_CE_CEModel;
		if (temp_pointer!=NULL)
			delete temp_pointer;
		temp_pointer = NULL;
		pMBDyn_CE_CEModel = NULL;
		std::cout << "destroy the CE_model...\n";
	}
}

extern "C" 
bool MBDyn_CE_CEModel_InitCheck(pMBDyn_CE_CEModel_t pMBDyn_CE_CEModel,
                                  const std::vector<double> &MBDyn_CE_CouplingKinematic,
                                  const unsigned &MBDyn_CE_NodesNum,
                                  const std::vector<MBDYN_CE_CEMODELDATA> &MBDyn_CE_CEModel_Label,
								  const double * pMBDyn_CE_Gravity)
{	
	if (pMBDyn_CE_CEModel==NULL)
	{
		std::cout << "Error: the C::E model pointer is NULL.\n";
		return false;
	}
	std::cout << "CE_models MBDyn_CE_CEModel_InitCheck():\n";
	ChSystemParallelNSC *tempsys = (ChSystemParallelNSC *)pMBDyn_CE_CEModel;
	//---------- 1. obtain the data from MBDyn;
	//---------- 2. transfer it to the coordinate in C::E;
	//---------- 3. check the consistency;

	//---------- 1. obtain the data from MBDyn;
	//----- MBDyn Global ref, expressed in C::E Ground Ref.
	const double *pmbdynce_tempvec3_x = &MBDyn_CE_CouplingKinematic[0];
	const double *pmbdynce_tempemat3x3_R = &MBDyn_CE_CouplingKinematic[3 * MBDyn_CE_NodesNum];
	const double *pmbdynce_tempvec3_xp = &MBDyn_CE_CouplingKinematic[12 * MBDyn_CE_NodesNum];
	const double *pmbdynce_tempvec3_omega = &MBDyn_CE_CouplingKinematic[15 * MBDyn_CE_NodesNum];
	const double *pmbdynce_tempvec3_xpp = &MBDyn_CE_CouplingKinematic[18 * MBDyn_CE_NodesNum];
	const double *pmbdynce_tempvec3_omegap = &MBDyn_CE_CouplingKinematic[21 * MBDyn_CE_NodesNum];
	const double *pmbdynce_temp_frame = &MBDyn_CE_CouplingKinematic[24 * MBDyn_CE_NodesNum];
	ChVector<> mbdynce_temp_frameMBDyn_pos(-pmbdynce_temp_frame[0], -pmbdynce_temp_frame[1], -pmbdynce_temp_frame[2]);
	ChVector<> mbdynce_temp_frameMBDyn_rot_axis_X(pmbdynce_temp_frame[3], pmbdynce_temp_frame[6], pmbdynce_temp_frame[9]);
	ChVector<> mbdynce_temp_frameMBDyn_rot_axis_Y(pmbdynce_temp_frame[4], pmbdynce_temp_frame[7], pmbdynce_temp_frame[10]);
	ChVector<> mbdynce_temp_frameMBDyn_rot_axis_Z(pmbdynce_temp_frame[5], pmbdynce_temp_frame[8], pmbdynce_temp_frame[11]);
	ChMatrix33<> mbdynce_temp_frameMBDyn_rot(mbdynce_temp_frameMBDyn_rot_axis_X,mbdynce_temp_frameMBDyn_rot_axis_Y,mbdynce_temp_frameMBDyn_rot_axis_Z);
	ChFrame<> mbdynce_temp_frameMBDyn(mbdynce_temp_frameMBDyn_pos, mbdynce_temp_frameMBDyn_rot);
	//----- Set the gravity of C::E model
	tempsys->Set_G_acc(mbdynce_temp_frameMBDyn.TransformDirectionLocalToParent(ChVector<>(pMBDyn_CE_Gravity[0],pMBDyn_CE_Gravity[1],pMBDyn_CE_Gravity[2])));
	std::cout <<"\t\tC::E gravity is set as: "<< tempsys->Get_G_acc() << "\n";
	
	//---------- 2. transfer it to the coordinate in C::E; (only check the consistence of position)
	bool bmbdynce_temp_check=true;
	for (unsigned i = 0; i < MBDyn_CE_NodesNum; i++)
	{
		//----- 2.1 MBDyn data
		ChVector<> mbdynce_tempmbdyn_pos = ChVector<>(pmbdynce_tempvec3_x[3 * i], pmbdynce_tempvec3_x[3 * i + 1], pmbdynce_tempvec3_x[3 * i + 2]) >> mbdynce_temp_frameMBDyn;
		//- ChVector<> mbdynce_tempmbdyn_pos_dt = mbdynce_temp_frameMBDyn.TransformDirectionLocalToParent(ChVector<>(pmbdynce_tempvec3_xp[3 * i], pmbdynce_tempvec3_xp[3 * i + 1], pmbdynce_tempvec3_xp[3 * i + 2]));
		ChMatrix33<> mbdynce_tempmbdyn_R1(ChVector<>(pmbdynce_tempemat3x3_R[9 * i], pmbdynce_tempemat3x3_R[9 * i + 1], pmbdynce_tempemat3x3_R[9 * i + 2]),
										  ChVector<>(pmbdynce_tempemat3x3_R[9 * i + 3], pmbdynce_tempemat3x3_R[9 * i + 4], pmbdynce_tempemat3x3_R[9 * i + 5]),
										  ChVector<>(pmbdynce_tempemat3x3_R[9 * i + 6], pmbdynce_tempemat3x3_R[9 * i + 7], pmbdynce_tempemat3x3_R[9 * i + 8])); // three column vectors
		ChQuaternion<> mbdynce_tempmbdyn_rot(mbdynce_tempmbdyn_R1.Get_A_quaternion() >> mbdynce_temp_frameMBDyn);
		//- ChVector<> mbdynce_tempmbdyn_Wvel_par=mbdynce_temp_frameMBDyn.TransformDirectionLocalToParent(ChVector<>(pmbdynce_tempvec3_omega[3 * i], pmbdynce_tempvec3_omega[3 * i + 1], pmbdynce_tempvec3_omega[3 * i + 2]));

		//----- 2.2 C::E data
		unsigned mbdynce_tempce_body_id = MBDyn_CE_CEModel_Label[i].MBDyn_CE_CEBody_Label; //- body ID
		auto mbdynce_tempce_body_i = tempsys->SearchBodyID(mbdynce_tempce_body_id);//- the corresponding body
		unsigned motor3d_motor_i_id = MBDyn_CE_CEModel_Label[i].MBDyn_CE_CEMotor_Label;
		auto motor3d_motor_i = std::dynamic_pointer_cast<ChLinkMotionImposed>(tempsys->SearchLinkID(motor3d_motor_i_id)); //- the corresponding motor
		ChFrame<> mbdynce_tempframe1b1_start = motor3d_motor_i->GetFrame1();
		ChFrame<> mbdynce_tempframe1G_start = mbdynce_tempframe1b1_start >> *(motor3d_motor_i->GetBody1());//- at initial, Frame_Node==Frame1G.
		//- ChVector<> mbdynce_tempce_pos = mbdynce_tempce_body_i->GetPos();
		ChVector<> mbdynce_tempce_pos = mbdynce_tempframe1G_start.GetPos();
		//- ChVector<> mbdynce_tempce_pos_dt = mbdynce_tempce_body_i->GetPos_dt();
		//- ChQuaternion<> mbdynce_tempce_rot = mbdynce_tempce_body_i->GetRot();
		ChQuaternion<> mbdynce_tempce_rot = mbdynce_tempframe1G_start.GetRot();
		//- ChVector<> mbdynce_tempce_Wvel_par = mbdynce_tempce_body_i->GetWvel_par(); //- angular velocity in global ref.

		//---------- 3. check the consistency;
		if (((mbdynce_tempmbdyn_pos - mbdynce_tempce_pos).Length()) > 1.0e-6)
		{
			bmbdynce_temp_check = false;
			std::cout << "\t\tpos of coupling body " << mbdynce_tempce_body_id << " doesn't agree with that in MBDyn\n";
			std::cout << "\t\tpos in MBDyn\t" << mbdynce_tempmbdyn_pos << "\n";
			std::cout << "\t\tpos in C::E\t" << mbdynce_tempce_pos << "\n";
			std::cout << "\t\tdifference between MBDyn and C::E\t" << (mbdynce_tempmbdyn_pos - mbdynce_tempce_pos).Length() << "\n";
			break;
		}
		if (mbdynce_tempmbdyn_rot.Dot(mbdynce_tempce_rot)<0)
		{
			mbdynce_tempmbdyn_rot = -mbdynce_tempmbdyn_rot;
		}
		if (((mbdynce_tempmbdyn_rot - mbdynce_tempce_rot).Length()) > 1.0e-6)
		{
			bmbdynce_temp_check = false;
			std::cout << "\t\trot of coupling body " << mbdynce_tempce_body_id << " doesn't agree with that in MBDyn\n";
			std::cout << "\t\trot in MBDyn\t" << mbdynce_tempmbdyn_R1.Get_A_quaternion() << "\n";
			std::cout << "\t\trot in C::E\t" << mbdynce_tempce_rot << "\n";
			std::cout << "\t\tFrame1 in C::E\t" << motor3d_motor_i->GetFrame1() << "\n";
			std::cout << "\t\tFrameBody1 in C::E\t" << *(motor3d_motor_i->GetBody1()) << "\n";
			std::cout << "\t\tdifference between MBDyn and C::E\t" << (mbdynce_tempmbdyn_rot - mbdynce_tempce_rot).Length() << "\n";
			break;
		}
		/*if (mbdynce_tempmbdyn_pos_dt!=mbdynce_tempce_pos_dt)
		{
			bmbdynce_temp_check = false;
			std::cout << "\t\tpos_dt of coupling body " << mbdynce_tempce_body_id << " doesn't agree with that in MBDyn\n";
			std::cout << "\t\tpos_dt in MBDyn\t" << mbdynce_tempmbdyn_pos_dt << "\n";
			std::cout << "\t\tpos_dt in C::E\t" << mbdynce_tempce_pos_dt << "\n";
			break;
		}
		if (mbdynce_tempmbdyn_Wvel_par!=mbdynce_tempce_Wvel_par)
		{
			bmbdynce_temp_check = false;
			std::cout << "\t\tWvel_par of coupling body " << mbdynce_tempce_body_id << " doesn't agree with that in MBDyn\n";
			std::cout << "\t\tWvel_par in MBDyn\t" << mbdynce_tempmbdyn_Wvel_par << "\n";
			std::cout << "\t\tWvel_par in C::E\t" << mbdynce_tempce_Wvel_par << "\n";
			break;
		}*/
	}
	if (bmbdynce_temp_check)
	{
		std::cout << "the C::E model is well established.\n\n";
	}
	return bmbdynce_temp_check;
}

//- save CEModel at current step for reloading them in the tight coupling scheme
//- (before advance())
extern "C" int
MBDyn_CE_CEModel_DataSave(pMBDyn_CE_CEModel_t pMBDyn_CE_CEModel, 
                        std::vector<double> & MBDyn_CE_CEModel_Data)
{
	if (pMBDyn_CE_CEModel==NULL)
	{
		std::cout << "Error: the C::E model pointer is NULL.\n";
		return 1;
	}
	ChSystemParallelNSC *tempsys = (ChSystemParallelNSC *)pMBDyn_CE_CEModel; 
	//- std::cout << "\t\tCE_models DataSave():\n";
	unsigned int tempsys_bodies_size = tempsys->Get_bodylist().size();
	unsigned int tempsys_size = (3 * 3 + 3 * 4) * tempsys_bodies_size + 1; //- +1 Chtime: Sys_size=body_size + 1(for time);
	unsigned int vector_size = MBDyn_CE_CEModel_Data.size();
	if (tempsys_size != vector_size)
	{
		std::cout << "Error: the vector to save data is not consistent with the C::E model:\n"; //- how to safely exit MBDyn?
		return 1;
	}

	//- save data
//#pragma omp parallel for
	for (unsigned int i = 0; i < tempsys_bodies_size; i++)
	{
		const ChVector<> &body_pos = tempsys->Get_bodylist()[i]->GetPos();		 //- 3 components
		const ChQuaternion<> &body_rot = tempsys->Get_bodylist()[i]->GetRot();	 //- 4 components
		const ChVector<> &body_pos_dt = tempsys->Get_bodylist()[i]->GetPos_dt(); //- 3 components
		const ChQuaternion<> &body_rot_dt = tempsys->Get_bodylist()[i]->GetRot_dt(); //- 4 components
		const ChVector<> &body_pos_dtdt = tempsys->Get_bodylist()[i]->GetPos_dtdt(); //- 3 components
		const ChQuaternion<> &body_rot_dtdt = tempsys->Get_bodylist()[i]->GetRot_dtdt(); //- 4 components

		unsigned int i_pos = (3 * 3 + 3 * 4) * i; 
		MBDyn_CE_CEModel_Data[i_pos] = body_pos.x();
		MBDyn_CE_CEModel_Data[i_pos+1] = body_pos.y();
		MBDyn_CE_CEModel_Data[i_pos+2] = body_pos.z();

		unsigned int i_rot = i_pos + 3;
		MBDyn_CE_CEModel_Data[i_rot] = body_rot.e0();
		MBDyn_CE_CEModel_Data[i_rot+1] = body_rot.e1();
		MBDyn_CE_CEModel_Data[i_rot+2] = body_rot.e2();
		MBDyn_CE_CEModel_Data[i_rot+3] = body_rot.e3();

		unsigned int i_pos_dt = i_rot + 4;
		MBDyn_CE_CEModel_Data[i_pos_dt] = body_pos_dt.x();
		MBDyn_CE_CEModel_Data[i_pos_dt+1] = body_pos_dt.y();
		MBDyn_CE_CEModel_Data[i_pos_dt+2] = body_pos_dt.z();

		unsigned int i_rot_dt = i_pos_dt + 3;
		MBDyn_CE_CEModel_Data[i_rot_dt] = body_rot_dt.e0();
		MBDyn_CE_CEModel_Data[i_rot_dt+1] = body_rot_dt.e1();
		MBDyn_CE_CEModel_Data[i_rot_dt+2] = body_rot_dt.e2();
		MBDyn_CE_CEModel_Data[i_rot_dt+3] = body_rot_dt.e3();

		unsigned int i_pos_dtdt = i_rot_dt + 4;
		MBDyn_CE_CEModel_Data[i_pos_dtdt] = body_pos_dtdt.x();
		MBDyn_CE_CEModel_Data[i_pos_dtdt+1] = body_pos_dtdt.y();
		MBDyn_CE_CEModel_Data[i_pos_dtdt+2] = body_pos_dtdt.z();

		unsigned int i_rot_dtdt = i_pos_dtdt + 3;
		MBDyn_CE_CEModel_Data[i_rot_dtdt] = body_rot_dtdt.e0();
		MBDyn_CE_CEModel_Data[i_rot_dtdt+1] = body_rot_dtdt.e1();
		MBDyn_CE_CEModel_Data[i_rot_dtdt+2] = body_rot_dtdt.e2();
		MBDyn_CE_CEModel_Data[i_rot_dtdt+3] = body_rot_dtdt.e3();
	}
	MBDyn_CE_CEModel_Data[tempsys_size-1] = tempsys->GetChTime();
	return 0;
}


//- reload data in the tight coupling scheme at each iteration
extern "C" int
MBDyn_CE_CEModel_DataReload(pMBDyn_CE_CEModel_t pMBDyn_CE_CEModel, 
                        std::vector<double> & MBDyn_CE_CEModel_Data)
{
	if (pMBDyn_CE_CEModel==NULL)
	{
		std::cout << "Error: the C::E model pointer is NULL.\n";
		return 1;
	}
	//- std::cout << "\t\tCE_models DataReload():\n";
	ChSystemParallelNSC *tempsys = (ChSystemParallelNSC *)pMBDyn_CE_CEModel;
	unsigned int tempsys_bodies_size = tempsys->Get_bodylist().size();
	unsigned int tempsys_size = (3 * 3 + 3 * 4) * tempsys_bodies_size + 1; //- +1 Chtime: Sys_size=body_size + 1(for time);
	unsigned int vector_size = MBDyn_CE_CEModel_Data.size();
	if (tempsys_size != vector_size)
	{
		std::cout << "Error: the vector used to save data is not consistent with the C::E model:\n"; //- how to safely exit MBDyn?
		return 1;
	}
	double mbdynce_ce_time = MBDyn_CE_CEModel_Data[tempsys_size - 1];
	tempsys->SetChTime(mbdynce_ce_time);
//- #pragma omp parallel for
	for (int i = 0; i < tempsys_bodies_size; i++)
	{
		unsigned int i_pos = (3 * 3 + 3 * 4) * i;
		tempsys->Get_bodylist()[i]->SetPos(ChVector<>(MBDyn_CE_CEModel_Data[i_pos],
										 MBDyn_CE_CEModel_Data[i_pos + 1],
										 MBDyn_CE_CEModel_Data[i_pos + 2]));
		unsigned int i_rot = i_pos + 3;
		tempsys->Get_bodylist()[i]->SetRot(ChQuaternion<>(MBDyn_CE_CEModel_Data[i_rot],
										 MBDyn_CE_CEModel_Data[i_rot + 1],
										 MBDyn_CE_CEModel_Data[i_rot + 2],
										 MBDyn_CE_CEModel_Data[i_rot + 3]));

		unsigned int i_pos_dt = i_rot + 4;
		tempsys->Get_bodylist()[i]->SetPos_dt(ChVector<>(MBDyn_CE_CEModel_Data[i_pos_dt],
													  MBDyn_CE_CEModel_Data[i_pos_dt + 1],
													  MBDyn_CE_CEModel_Data[i_pos_dt + 2]));

		unsigned int i_rot_dt = i_pos_dt + 3;
		tempsys->Get_bodylist()[i]->SetRot_dt(ChQuaternion<>(MBDyn_CE_CEModel_Data[i_rot_dt],
														  MBDyn_CE_CEModel_Data[i_rot_dt + 1],
														  MBDyn_CE_CEModel_Data[i_rot_dt + 2],
														  MBDyn_CE_CEModel_Data[i_rot_dt + 3]));
		
		unsigned int i_pos_dtdt = i_rot_dt + 4;
		tempsys->Get_bodylist()[i]->SetPos_dtdt(ChVector<>(MBDyn_CE_CEModel_Data[i_pos_dtdt],
													  MBDyn_CE_CEModel_Data[i_pos_dtdt + 1],
													  MBDyn_CE_CEModel_Data[i_pos_dtdt + 2]));

		unsigned int i_rot_dtdt = i_pos_dtdt + 3;
		tempsys->Get_bodylist()[i]->SetRot_dtdt(ChQuaternion<>(MBDyn_CE_CEModel_Data[i_rot_dtdt],
														  MBDyn_CE_CEModel_Data[i_rot_dtdt + 1],
														  MBDyn_CE_CEModel_Data[i_rot_dtdt + 2],
														  MBDyn_CE_CEModel_Data[i_rot_dtdt + 3]));
		tempsys->Get_bodylist()[i]->Update(mbdynce_ce_time);
	}
	return 0;
}

extern "C" int
MBDyn_CE_CEModel_DoStepDynamics(pMBDyn_CE_CEModel_t pMBDyn_CE_CEModel, double time_step, bool bMBDyn_CE_Verbose)
{
	if (pMBDyn_CE_CEModel==NULL)
	{
		std::cout << "Error: the C::E model pointer is NULL.\n";
		return 1;
	}
	//- std::cout << "\t\tCE_models DoStepDynamics():\n";
	ChSystemParallelNSC *tempsys = (ChSystemParallelNSC *)pMBDyn_CE_CEModel;
	if (bMBDyn_CE_Verbose)
	{
		std::cout << "\tC::E model before integration:\n";
		for (unsigned i = 0; i < tempsys->Get_bodylist().size(); i++)
		{
			std::cout << "\t\tBody " << tempsys->Get_bodylist()[i]->GetIdentifier() << "\n";
			std::cout << "\t\t\tpos: " << tempsys->Get_bodylist()[i]->GetPos() << "\n";
			std::cout << "\t\t\tpos_dt: " << tempsys->Get_bodylist()[i]->GetPos_dt() << "\n";
			std::cout << "\t\t\tpos_dtdt: " << tempsys->Get_bodylist()[i]->GetPos_dtdt() << "\n";
			std::cout << "\t\t\trot: " << tempsys->Get_bodylist()[i]->GetRot().Q_to_Euler123()/CH_C_PI*180.0 << "\n";
			std::cout << "\t\t\tWvel_par: " << tempsys->Get_bodylist()[i]->GetWvel_par() << "\n";
			std::cout << "\t\t\tWacc_par: " << tempsys->Get_bodylist()[i]->GetWacc_par() << "\n";
		}
		for (unsigned i = 0; i < tempsys->Get_linklist().size(); i++)
		{
			std::cout << "\t\tLink " << tempsys->Get_linklist()[i]->GetIdentifier() << "\n";
			auto motor_3d=std::dynamic_pointer_cast<ChLinkMotionImposed>(tempsys->Get_linklist()[i]);
			if (motor_3d)
			{
				//- some information about the motion of marker on body 1
				//- start velocity in C::E Ground
				ChVector<> frame1_pos_loc, frame1_pos_par, frame1_pos_dt_par, frame1_pos_dtdt_par;				
				ChVector<> frame1_Wvel_loc, frame1_Wvel_par, frame1_Wacc_loc, frame1_Wacc_par;
				ChQuaternion<> frame1_q_loc, frame1_q_par;
				frame1_pos_loc = motor_3d->GetFrame1().GetPos();
				frame1_pos_par = motor_3d->GetBody1()->TransformDirectionLocalToParent(frame1_pos_loc);
				frame1_q_loc = motor_3d->GetFrame1().GetRot();
				frame1_q_par = motor_3d->GetBody1()->GetRot()*frame1_q_loc;
				frame1_Wvel_loc = motor_3d->GetBody1()->GetWvel_loc();
				frame1_Wvel_par = motor_3d->GetBody1()->GetWvel_par();
				frame1_Wacc_loc = motor_3d->GetBody1()->GetWacc_loc();
				frame1_Wacc_par = motor_3d->GetBody1()->GetWacc_par();
				frame1_pos_dt_par.Cross(frame1_Wvel_par, frame1_pos_par);
				frame1_pos_dtdt_par.Cross(frame1_Wvel_par, frame1_pos_dt_par);
				frame1_pos_dt_par += motor_3d->GetBody1()->GetPos_dt();
				frame1_pos_dtdt_par += motor_3d->GetBody1()->GetPos_dtdt() + frame1_Wacc_par.Cross(frame1_pos_par);
				frame1_pos_par += motor_3d->GetBody1()->GetPos();

				std::cout << "\t\t\tframe1_pos_par : " << frame1_pos_par  << "\n";
				std::cout << "\t\t\tframe1_pos_dt_par : " << frame1_pos_dt_par  << "\n";
				std::cout << "\t\t\tframe1_pos_dtdt_par : " << frame1_pos_dtdt_par << "\n";
				std::cout << "\t\t\tframe1_q_par : " << frame1_q_par  << "\n";
				std::cout << "\t\t\tframe1_Wvel_par : " << frame1_Wvel_par << "\n";
				std::cout << "\t\t\tframe1_Wacc_par : " << frame1_Wacc_par << "\n";
				
				std::cout << "\t\t\tframe1pos: " << motor_3d->GetFrame1().GetPos() << "\n"; //- marker on body 1 (local ref body 1)
				std::cout << "\t\t\tframe1rot: " << motor_3d->GetFrame1().GetRot()<< "\n";
				std::cout << "\t\t\tframe2pos: " << motor_3d->GetFrame2().GetPos() << "\n"; //- marker on body 2 (local ref body 2)
				std::cout << "\t\t\tframe2rot: " << motor_3d->GetFrame2().GetRot() << "\n";
				std::cout << "\t\t\tframeM2pos: " << motor_3d->GetFrameM2().GetPos() << "\n"; //- marker from mbdyn (local ref body 2)
				std::cout << "\t\t\tframeM2rot: " << motor_3d->GetFrameM2().GetRot()<< "\n";
				std::cout << "\t\t\tlinkpos_abs: " << motor_3d->GetLinkAbsoluteCoords().pos << "\n";
				std::cout << "\t\t\tlinkrot_abs: " << motor_3d->GetLinkAbsoluteCoords().rot.Q_to_Euler123()/CH_C_PI*180.0 << "\n";
				std::cout << "\t\t\tbody1pos: " << motor_3d->GetBody1()->GetPos() << "\n";
				std::cout << "\t\t\tbody1rot: " << motor_3d->GetBody1()->GetRot().Q_to_Euler123()/CH_C_PI*180.0 << "\n";
				std::cout << "\t\t\tbody1pos_dt: " << motor_3d->GetBody1()->GetPos_dt() << "\n";
				std::cout << "\t\t\tbody1Wvel_par: " << motor_3d->GetBody1()->GetWvel_par() << "\n";
				std::cout << "\t\t\tbody1pos_dtdt: " << motor_3d->GetBody1()->GetPos_dtdt() << "\n";
				std::cout << "\t\t\tbody1Wacc_par: " << motor_3d->GetBody1()->GetWacc_par() << "\n";
				std::cout << "\t\t\tbody2pos: " << motor_3d->GetBody2()->GetPos() << "\n";
				std::cout << "\t\t\tbody2rot: " << motor_3d->GetBody2()->GetRot().Q_to_Euler123()/CH_C_PI*180.0 << "\n";
				std::cout << "\t\t\tbody2pos_dt: " << motor_3d->GetBody2()->GetPos_dt() << "\n";
				std::cout << "\t\t\tbody2Wvel_par: " << motor_3d->GetBody2()->GetWvel_par() << "\n";
				std::cout << "\t\t\tbody2pos_dtdt: " << motor_3d->GetBody2()->GetPos_dtdt() << "\n";
				std::cout << "\t\t\tbody2Wacc_par: " << motor_3d->GetBody2()->GetWacc_par() << "\n";
			}
		}
	}
#ifdef CHRONO_OPENGL 
//#if 0
	opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
	if (!gl_window.Active())
	{
		std::cout << "GL_WINDOW is close.\n";
		return 1;
	}
	gl_window.DoStepDynamics(time_step);						
	gl_window.Render();
#else
	tempsys->DoStepDynamics(time_step);
#endif
	tempsys->CalculateContactForces();
	if (bMBDyn_CE_Verbose)
	{
		std::cout << "\t\ttime: " << tempsys->GetChTime() << " s\n";
		std::cout << "\t\ttime step: " << time_step << " s\n";
	}
	//- print detailed information about the C::E model
	if (bMBDyn_CE_Verbose)
	{
		std::cout << "\tC::E model after integration:\n";
		for (unsigned i = 0; i < tempsys->Get_bodylist().size(); i++)
		{
			std::cout << "\t\tBody " << tempsys->Get_bodylist()[i]->GetIdentifier() << "\n";
			std::cout << "\t\t\tpos: " << tempsys->Get_bodylist()[i]->GetPos() << "\n";
			std::cout << "\t\t\tpos_dt: " << tempsys->Get_bodylist()[i]->GetPos_dt() << "\n";
			std::cout << "\t\t\tpos_dtdt: " << tempsys->Get_bodylist()[i]->GetPos_dtdt() << "\n";
			std::cout << "\t\t\trot: " << tempsys->Get_bodylist()[i]->GetRot().Q_to_Euler123()/CH_C_PI*180.0 << "\n";
			std::cout << "\t\t\tWvel_par: " << tempsys->Get_bodylist()[i]->GetWvel_par() << "\n";
			std::cout << "\t\t\tWacc_par: " << tempsys->Get_bodylist()[i]->GetWacc_par() << "\n";
		}
		for (unsigned i = 0; i < tempsys->Get_linklist().size(); i++)
		{
			std::cout << "\t\tLink " << tempsys->Get_linklist()[i]->GetIdentifier() << "\n";
			auto motor_3d=std::dynamic_pointer_cast<ChLinkMotionImposed>(tempsys->Get_linklist()[i]);
			if (motor_3d)
			{
				//- some information about the motion of marker on body 1
				//- start velocity in C::E Ground
				ChVector<> frame1_pos_loc, frame1_pos_par, frame1_pos_par2, frame1_pos_dt_par, frame1_pos_dt_par2, frame1_pos_dtdt_par, frame1_pos_dtdt_par2;				
				ChVector<> frame1_Wvel_loc, frame1_Wvel_par, frame1_Wacc_loc, frame1_Wacc_par;
				ChQuaternion<> frame1_q_loc, frame1_q_par;
				frame1_pos_loc = motor_3d->GetFrame1().GetPos();
				frame1_pos_par = motor_3d->GetBody1()->TransformDirectionLocalToParent(frame1_pos_loc);
				
				frame1_q_loc = motor_3d->GetFrame1().GetRot();
				frame1_q_par = motor_3d->GetBody1()->GetRot()*frame1_q_loc;
				frame1_Wvel_loc = motor_3d->GetBody1()->GetWvel_loc();
				frame1_Wvel_par = motor_3d->GetBody1()->GetWvel_par();
				frame1_Wacc_loc = motor_3d->GetBody1()->GetWacc_loc();
				frame1_Wacc_par = motor_3d->GetBody1()->GetWacc_par();
				frame1_pos_dt_par.Cross(frame1_Wvel_par, frame1_pos_par);
				frame1_pos_dtdt_par.Cross(frame1_Wvel_par, frame1_pos_dt_par);
				frame1_pos_dt_par += motor_3d->GetBody1()->GetPos_dt();
				frame1_pos_dtdt_par += motor_3d->GetBody1()->GetPos_dtdt() + frame1_Wacc_par.Cross(frame1_pos_par);
				frame1_pos_par += motor_3d->GetBody1()->GetPos();

				frame1_pos_par2 = motor_3d->GetBody1()->TransformDirectionLocalToParent(frame1_pos_loc);
				frame1_pos_dt_par2.Cross(frame1_Wvel_par, frame1_pos_par2);
				frame1_pos_dtdt_par2.Cross(frame1_Wvel_par, frame1_pos_dt_par2);
				frame1_pos_dt_par2 += motor_3d->GetBody1()->GetPos_dt();
				frame1_pos_dtdt_par2 += motor_3d->GetBody1()->GetPos_dtdt() + frame1_Wacc_par.Cross(frame1_pos_par2);
				frame1_pos_par2 += motor_3d->GetBody1()->GetPos();

				std::cout << "\t\t\tframe1_pos_par : " << frame1_pos_par  << "\n";
				std::cout << "\t\t\tframe1_pos_dt_par : " << frame1_pos_dt_par  << "\n";
				std::cout << "\t\t\tframe1_pos_dtdt_par : " << frame1_pos_dtdt_par << "\n";
				std::cout << "\t\t\tframe1_pos_par2 : " << frame1_pos_par2  << "\n";
				std::cout << "\t\t\tframe1_pos_dt_par2 : " << frame1_pos_dt_par2  << "\n";
				std::cout << "\t\t\tframe1_pos_dtdt_par2 : " << frame1_pos_dtdt_par2 << "\n";
				std::cout << "\t\t\tframe1_q_par : " << frame1_q_par  << "\n";
				std::cout << "\t\t\tframe1_Wvel_par : " << frame1_Wvel_par << "\n";
				std::cout << "\t\t\tframe1_Wacc_par : " << frame1_Wacc_par << "\n";

				std::cout << "\t\t\tframe1pos: " << motor_3d->GetFrame1().GetPos() << "\n"; //- marker on body 1 (local ref body 1)
				std::cout << "\t\t\tframe1rot: " << motor_3d->GetFrame1().GetRot() << "\n";
				std::cout << "\t\t\tframe2pos: " << motor_3d->GetFrame2().GetPos() << "\n"; //- marker on body 2 (local ref body 2)
				std::cout << "\t\t\tframe2rot: " << motor_3d->GetFrame2().GetRot() << "\n";
				std::cout << "\t\t\tframeM2pos: " << motor_3d->GetFrameM2().GetPos() << "\n"; //- marker from mbdyn (local ref body 2)
				std::cout << "\t\t\tframeM2rot: " << motor_3d->GetFrameM2().GetRot() << "\n";
				std::cout << "\t\t\tlinkpos_abs: " << motor_3d->GetLinkAbsoluteCoords().pos << "\n";
				std::cout << "\t\t\tlinkrot_abs: " << motor_3d->GetLinkAbsoluteCoords().rot << "\n";
				std::cout << "\t\t\tbody1pos: " << motor_3d->GetBody1()->GetPos() << "\n";
				std::cout << "\t\t\tbody1rot: " << motor_3d->GetBody1()->GetRot() << "\n";
				std::cout << "\t\t\tbody1pos_dt: " << motor_3d->GetBody1()->GetPos_dt() << "\n";
				std::cout << "\t\t\tbody1Wvel_par: " << motor_3d->GetBody1()->GetWvel_par() << "\n";
				std::cout << "\t\t\tbody1Wvel_loc: " << motor_3d->GetBody1()->GetWvel_loc() << "\n";
				std::cout << "\t\t\tbody1pos_dtdt: " << motor_3d->GetBody1()->GetPos_dtdt() << "\n";
				std::cout << "\t\t\tbody1Wacc_par: " << motor_3d->GetBody1()->GetWacc_par() << "\n";
				std::cout << "\t\t\tbody1Wacc_loc: " << motor_3d->GetBody1()->GetWacc_loc() << "\n";
				std::cout << "\t\t\tbody2pos: " << motor_3d->GetBody2()->GetPos() << "\n";
				std::cout << "\t\t\tbody2rot: " << motor_3d->GetBody2()->GetRot()<< "\n";
				std::cout << "\t\t\tbody2pos_dt: " << motor_3d->GetBody2()->GetPos_dt() << "\n";
				std::cout << "\t\t\tbody2Wvel_par: " << motor_3d->GetBody2()->GetWvel_par() << "\n";
				std::cout << "\t\t\tbody2pos_dtdt: " << motor_3d->GetBody2()->GetPos_dtdt() << "\n";
				std::cout << "\t\t\tbody2Wacc_par: " << motor_3d->GetBody2()->GetWacc_par() << "\n";
			}
		}
	}
	return 0;
}

//- C::E models receive coupling motion from the buffer
extern "C" int 
MBDyn_CE_CEModel_RecvFromBuf(pMBDyn_CE_CEModel_t pMBDyn_CE_CEModel,
const std::vector<double>& MBDyn_CE_CouplingKinematic, 
const unsigned& MBDyn_CE_NodesNum,
const std::vector<MBDYN_CE_CEMODELDATA> & MBDyn_CE_CEModel_Label,
const int MBDyn_CE_CEMotorType,
double* MBDyn_CE_CEMotorType_coeff, //- only used for motor type: acceleration
double time_step,
bool bMBDyn_CE_Verbose)
{
	if (pMBDyn_CE_CEModel==NULL)
	{
		std::cout << "Error: the C::E model pointer is NULL.\n";
		return 1;
	}
	//- std::cout << "\t\tCE_models RecvFromMBDyn():\n";
	ChSystemParallelNSC *tempsys = (ChSystemParallelNSC *)pMBDyn_CE_CEModel;
	//---------- 1. obtain the data;
	//---------- 2. transfer it to the coordinate in C::E;
	//---------- 3. update motor functions;

	//---------- 1. obtain the data;
	const double *pmbdynce_tempvec3_x = &MBDyn_CE_CouplingKinematic[0];
	const double *pmbdynce_tempemat3x3_R = &MBDyn_CE_CouplingKinematic[3 * MBDyn_CE_NodesNum];
	const double *pmbdynce_tempvec3_xp = &MBDyn_CE_CouplingKinematic[12 * MBDyn_CE_NodesNum];
	const double *pmbdynce_tempvec3_omega = &MBDyn_CE_CouplingKinematic[15 * MBDyn_CE_NodesNum];
	const double *pmbdynce_tempvec3_xpp = &MBDyn_CE_CouplingKinematic[18 * MBDyn_CE_NodesNum];
	const double *pmbdynce_tempvec3_omegap = &MBDyn_CE_CouplingKinematic[21 * MBDyn_CE_NodesNum];
	const double *pmbdynce_temp_frame = &MBDyn_CE_CouplingKinematic[24 * MBDyn_CE_NodesNum];
	ChVector<> mbdynce_temp_frameMBDyn_pos(-pmbdynce_temp_frame[0], -pmbdynce_temp_frame[1], -pmbdynce_temp_frame[2]);
	ChVector<> mbdynce_temp_frameMBDyn_rot_axis_X(pmbdynce_temp_frame[3], pmbdynce_temp_frame[6], pmbdynce_temp_frame[9]);
	ChVector<> mbdynce_temp_frameMBDyn_rot_axis_Y(pmbdynce_temp_frame[4], pmbdynce_temp_frame[7], pmbdynce_temp_frame[10]);
	ChVector<> mbdynce_temp_frameMBDyn_rot_axis_Z(pmbdynce_temp_frame[5], pmbdynce_temp_frame[8], pmbdynce_temp_frame[11]);
	ChMatrix33<> mbdynce_temp_frameMBDyn_rot(mbdynce_temp_frameMBDyn_rot_axis_X,mbdynce_temp_frameMBDyn_rot_axis_Y,mbdynce_temp_frameMBDyn_rot_axis_Z);
	ChFrame<> mbdynce_temp_frameMBDyn(mbdynce_temp_frameMBDyn_pos, mbdynce_temp_frameMBDyn_rot);
	double time = tempsys->GetChTime();
	//---------- 2. transfer it to the coordinate in C::E, and 3. update motor functions
	for (unsigned i = 0; i < MBDyn_CE_NodesNum;i++)
	{
		//----- 2.1 coordinate transformation
		//- the C::E ground keeps still relative to MBDyn ground.
		ChVector<> mbdynce_tempmbdyn_pos = ChVector<>(pmbdynce_tempvec3_x[3 * i], pmbdynce_tempvec3_x[3 * i + 1], pmbdynce_tempvec3_x[3 * i + 2]) >> mbdynce_temp_frameMBDyn;
		ChVector<> mbdynce_tempmbdyn_pos_dt = mbdynce_temp_frameMBDyn.TransformDirectionLocalToParent(ChVector<>(pmbdynce_tempvec3_xp[3 * i], pmbdynce_tempvec3_xp[3 * i + 1], pmbdynce_tempvec3_xp[3 * i + 2]));
		ChVector<> mbdynce_tempmbdyn_pos_dtdt = mbdynce_temp_frameMBDyn.TransformDirectionLocalToParent(ChVector<>(pmbdynce_tempvec3_xpp[3 * i], pmbdynce_tempvec3_xpp[3 * i + 1], pmbdynce_tempvec3_xpp[3 * i + 2]));
		ChMatrix33<> mbdynce_tempmbdyn_R1(ChVector<>(pmbdynce_tempemat3x3_R[9 * i], pmbdynce_tempemat3x3_R[9 * i + 1], pmbdynce_tempemat3x3_R[9 * i + 2]),
										  ChVector<>(pmbdynce_tempemat3x3_R[9 * i + 3], pmbdynce_tempemat3x3_R[9 * i + 4], pmbdynce_tempemat3x3_R[9 * i + 5]),
										  ChVector<>(pmbdynce_tempemat3x3_R[9 * i + 6], pmbdynce_tempemat3x3_R[9 * i + 7], pmbdynce_tempemat3x3_R[9 * i + 8])); // three column vectors
		ChMatrix33<> mbdynce_tempmbdyn_R(mbdynce_tempmbdyn_R1.Get_A_quaternion() >> mbdynce_temp_frameMBDyn);
		ChVector<> mbdynce_tempmbdyn_Wvel_par = mbdynce_temp_frameMBDyn.TransformDirectionLocalToParent(ChVector<>(pmbdynce_tempvec3_omega[3 * i], pmbdynce_tempvec3_omega[3 * i + 1], pmbdynce_tempvec3_omega[3 * i + 2]));
		ChVector<> mbdynce_tempmbdyn_Wacc_par = mbdynce_temp_frameMBDyn.TransformDirectionLocalToParent(ChVector<>(pmbdynce_tempvec3_omegap[3 * i], pmbdynce_tempvec3_omegap[3 * i + 1], pmbdynce_tempvec3_omegap[3 * i + 2]));
		ChFrame<> mbdynce_tempMBDynG_end(mbdynce_tempmbdyn_pos, mbdynce_tempmbdyn_R); //- pos and rot in ground ref.
		//----- 2.2 create motor functions
		//- find motor i
		unsigned motor3d_motor_i_id = MBDyn_CE_CEModel_Label[i].MBDyn_CE_CEMotor_Label;
		auto motor3d_motor_i = std::dynamic_pointer_cast<ChLinkMotionImposed>(tempsys->SearchLinkID(motor3d_motor_i_id));
		//- frame_start
		if (motor3d_motor_i!=NULL)
		{
			auto motor3d_function_pos_base = motor3d_motor_i->GetPositionFunction();
			auto motor3d_function_rot_base = motor3d_motor_i->GetRotationFunction();	
			
			if (motor3d_function_pos_base != NULL & motor3d_function_rot_base != NULL)
			{
				ChFrame<> mbdynce_tempframe1b1_start, mbdynce_tempframe1G_start, mbdynce_tempframeM2_start, mbdynce_tempframeM2_end;
				ChFrame<> mbdynce_tempframe2G(((motor3d_motor_i->GetFrame2()) >> *motor3d_motor_i->GetBody2()));
				ChVector<> mbdynce_tempframeM2_end_pos_dt, mbdynce_tempframeM2_end_pos_dtdt;
				ChVector<> mbdynce_tempframeM2_end_Wvel_loc, mbdynce_tempframeM2_end_Wacc_loc; // angular velocity and acceleration in frame M loc

				//- start/end position in frame 2 x_{B,n} = x_{2,n} + G_{n}*r^{0}_{B}
				mbdynce_tempframe1b1_start = motor3d_motor_i->GetFrame1();
				mbdynce_tempframe1G_start = mbdynce_tempframe1b1_start >> *(motor3d_motor_i->GetBody1());			
				// mbdynce_tempframe1G_start.SetPos(motor3d_motor_i->GetBody1()->GetPos()+mbdynce_temp_offset);
				mbdynce_tempframeM2_start = mbdynce_tempframe1G_start >> (mbdynce_tempframe2G.GetInverse()); // expressed in Frame 2
				mbdynce_tempframeM2_end = mbdynce_tempMBDynG_end >> (mbdynce_tempframe2G.GetInverse()); // pos and rot of the node expressed in Frame 2
				
				//- start velocity in C::E Ground
				ChVector<> mbdynce_tempframeG_start_pos_dt, mbdynce_temp_start_vel_rel;
				//ChVector<> mbdynce_temp_offset(motor3d_motor_i->GetBody1()->TransformDirectionLocalToParent(mbdynce_tempframe1b1_start.GetPos()));
				//------ Chrono uses the orientation of last step, when calculating the velocity of the marker. v_{B,n}= G_{n}*W^{l}_{n}\times G_{n-1}*r^{0}_{B} + v_{n}
				ChQuaternion<> FrameBody1_prev_rot; //- if tight scheme is used, how to save the motor function ?
				if (time<time_step) //- in the fisrt step, using the G_{-1} = G_{0} for v_{B,0}
				{
					FrameBody1_prev_rot = motor3d_motor_i->GetBody1()->GetRot();
				}
				else //- from step n->n+1, using G_{n-1} for v_{B,n}.
				{
					FrameBody1_prev_rot.Cross(motor3d_function_rot_base->Get_q(time), mbdynce_tempframe1b1_start.GetRot().GetConjugate());
				}				
				ChVector<> mbdynce_temp_offset(FrameBody1_prev_rot.Rotate(mbdynce_tempframe1b1_start.GetPos()));
				mbdynce_temp_start_vel_rel.Cross(motor3d_motor_i->GetBody1()->GetWvel_par(), mbdynce_temp_offset);
				mbdynce_tempframeG_start_pos_dt = motor3d_motor_i->GetBody1()->GetPos_dt() + mbdynce_temp_start_vel_rel;
				//- start acc in C::E Ground a_{B,n} = a_{2,n} + G_{n}*W^{l}_{n}\times G_{n}*W^{l}_{n}\times G_{n-1}*r^{0}_{B} +G_{n} Wacc^{l}_{n}\times G_{n-1}*r^{0}_{B}
				ChVector<> mbdynce_temp_start_acc, mbdynce_temp_start_acc_rel;
				mbdynce_temp_start_acc_rel.Cross(motor3d_motor_i->GetBody1()->GetWacc_par(), mbdynce_temp_offset);
				mbdynce_temp_start_acc = motor3d_motor_i->GetBody1()->GetPos_dtdt() + mbdynce_temp_start_acc_rel+motor3d_motor_i->GetBody1()->GetWvel_par().Cross(mbdynce_temp_start_vel_rel);

				if (bMBDyn_CE_Verbose)
				{
					// data receive from MBDyn
					std::cout << "\tData receives from MBDyn\n\t\tmotor ID: " << motor3d_motor_i_id << " \n";//motions start at: \n";
					std::cout << "\t\t\tpos in frame C::E Ground: " << mbdynce_tempmbdyn_pos << "\n";
					std::cout << "\t\t\tpos_dt in frame C::E Ground: " << mbdynce_tempmbdyn_pos_dt << "\n";
					//std::cout << "\t\t\tpos_dt in frame 2: " << mbdynce_tempframeM2_end_pos_dt << "\n";
					std::cout << "\t\t\tpos_dtdt in frame C::E Ground: " << mbdynce_tempmbdyn_pos_dtdt << "\n";
					std::cout << "\t\t\trot in frame C::E Ground: " << mbdynce_tempMBDynG_end.GetRot()<< "\n";
					std::cout << "\t\t\trot_dt in frame C::E Ground: " << mbdynce_tempmbdyn_Wvel_par << "\n";
					std::cout << "\t\t\trot_dtdt in frame C::E Ground: " << mbdynce_tempmbdyn_Wacc_par << "\n";
					std::cout << "\t\t\tmbdynce_temp_offset: " << mbdynce_temp_offset << "\n";
					std::cout << "\t\t\tframeMG_prev: " << motor3d_function_rot_base->Get_q(time) << "\n";
					std::cout << "\t\t\tmbdynce_tempframe1G_start_pos: " << mbdynce_tempframe1G_start.GetPos() << "\n";
					std::cout << "\t\t\tmbdynce_tempframe1G_start_rot: " << mbdynce_tempframeM2_start.GetRot()<< "\n";
					std::cout << "\t\t\tmbdynce_tempframe1G_start_posdt: " << mbdynce_tempframeG_start_pos_dt << "\n";
					std::cout << "\t\t\tmbdynce_tempframe1G_start_posdt_loc: " << motor3d_function_rot_base->Get_q(time).RotateBack(mbdynce_tempframeG_start_pos_dt) << "\n";
					std::cout << "\t\t\tmbdynce_tempframe1G_start_rotdt: " << motor3d_motor_i->GetBody1()->GetWvel_par() << "\n";
					std::cout << "\t\t\tmbdynce_tempframe1G_start_rotdt_loc: " << motor3d_motor_i->GetBody1()->GetWvel_loc() << "\n";
					std::cout << "\t\t\tmbdynce_tempframe1G_start_posdtdt: " << mbdynce_temp_start_acc << "\n";
					std::cout << "\t\t\tmbdynce_tempframe1G_start_rotdtdt: " << motor3d_motor_i->GetBody1()->GetWacc_par() << "\n";
				}

				if (MBDyn_CE_CEMotorType==MBDyn_CE_CEMOTORTYPE::ACCELERATION) //- accelerations
				{
					// using setpoint functions for position and rotation
					auto motor3d_function_pos = std::dynamic_pointer_cast<ChFunctionPosition_setpoint>(motor3d_function_pos_base);
					auto motor3d_function_rot = std::dynamic_pointer_cast<ChFunctionRotation_setpoint>(motor3d_function_rot_base);
					
					//- the velocity includes three part, 1. {(x_k+1-x_k)/h}, 2. {xp_k+1}, 3. {xp_k+h*xpp_k+1}
					//- 1. {(x_k+1-x_k)/h} velocity by position difference
					ChVector<> mbdynce_accmotor_posdt1, mbdynce_accmotor_posdt2, mbdynce_accmotor_posdt3;
					ChVector<> mbdynce_accmotor_w_loc1, mbdynce_accmotor_w_loc2, mbdynce_accmotor_w_loc3;
					//- position function using line/spline interpolation, in frame 2
					auto mbdynce_temp_pos_func = std::make_shared<ChFunctionPosition_line>(); //- a position function
					auto mbdynce_temp_pos_line = chrono_types::make_shared<geometry::ChLineSegment>(mbdynce_tempframeM2_start.GetPos(), mbdynce_tempframeM2_end.GetPos());
					mbdynce_temp_pos_func->SetLine(mbdynce_temp_pos_line);
					//- chrono_types::make_shared<>: a more safety case in C::E
					mbdynce_temp_pos_func->SetSpaceFunction(chrono_types::make_shared<ChFunction_Ramp>(-time / time_step, 1 / time_step));
					//- rotation function using Slerp interpolation, in frame 2
					auto mbdynce_temp_rot_func = std::make_shared<ChFunctionRotation_spline>();
					if (mbdynce_tempframeM2_start.GetRot().Dot(mbdynce_tempframeM2_end.GetRot()) < 0)
					{
						mbdynce_tempframeM2_end.GetRot() = -mbdynce_tempframeM2_end.GetRot();
					}
					std::vector<ChQuaternion<>> mbdynce_temp_rot_spline = {{mbdynce_tempframeM2_start.GetRot()}, {mbdynce_tempframeM2_end.GetRot()}};
					mbdynce_temp_rot_func->SetupData(1, mbdynce_temp_rot_spline);
					mbdynce_temp_rot_func->SetSpaceFunction(chrono_types::make_shared<ChFunction_Ramp>(-time / time_step, 1 / time_step));
					mbdynce_accmotor_posdt1 = mbdynce_temp_pos_func->Get_p_ds(time); //- constant velocity, expressed in frame 2
					mbdynce_accmotor_w_loc1 = mbdynce_temp_rot_func->Get_w_loc(time); //- constant angular velocity, expressed in local ref: {mbdynce_tempframeM2_start.GetRot()}
					
					//- 2. {xp_k+1} velocity obtained from coupled velocity.
					//- linear velocity and acceleration expressed in frame 2 (setpoint function);
					mbdynce_tempframeM2_end_pos_dt = mbdynce_tempframe2G.TransformDirectionParentToLocal(mbdynce_tempmbdyn_pos_dt);
					mbdynce_tempframeM2_end_pos_dtdt = mbdynce_tempframe2G.TransformDirectionParentToLocal(mbdynce_tempmbdyn_pos_dtdt);
					//- angular veloctiy/acceleration expressed in frame M (in real: mbdynce_tempframeM2_start=mbdynce_tempframeM2_end) (setpoint function)
					mbdynce_tempframeM2_end_Wvel_loc = mbdynce_tempframe2G.TransformDirectionParentToLocal(mbdynce_tempmbdyn_Wvel_par); //- angular velocity: expressed in frame 2
					mbdynce_tempframeM2_end_Wacc_loc = mbdynce_tempframe2G.TransformDirectionParentToLocal(mbdynce_tempmbdyn_Wacc_par);
					mbdynce_tempframeM2_end_Wvel_loc = mbdynce_tempframeM2_start.TransformDirectionParentToLocal(mbdynce_tempframeM2_end_Wvel_loc); //- angular velocity: expressed in frame M
					mbdynce_tempframeM2_end_Wacc_loc = mbdynce_tempframeM2_start.TransformDirectionParentToLocal(mbdynce_tempframeM2_end_Wacc_loc); //- should be expressed in the get_q(s) frame, namely, frame M.
					mbdynce_accmotor_posdt2 = mbdynce_tempframeM2_end_pos_dt;
					mbdynce_accmotor_w_loc2 = mbdynce_tempframeM2_end_Wvel_loc;
					//- 3. {xp_k+h*xpp_k+1} velocity obtained from accleration*h.
					//ChVector<> mbdynce_tempframeG_start_pos_dt, mbdynce_tempframeG_start_Wvel_par;
					//mbdynce_tempframeG_start_pos_dt = motor3d_motor_i->GetBody1()->GetPos_dt(); //- pos_dt expressed in C::E global ref.
					/*ChVector<> mbdynce_tempframeG_start_Wvel_par;
					mbdynce_tempframeG_start_Wvel_par = motor3d_motor_i->GetBody1()->GetWvel_par(); //- Wvel expressed in C::E global ref. 
					ChVector<> mbdynce_tempframeM2_start_pos_dt, mbdynce_tempframeM2_start_Wvel_loc;
					mbdynce_tempframeM2_start_pos_dt = mbdynce_tempframe2G.TransformDirectionParentToLocal(mbdynce_tempframeG_start_pos_dt); //- expressed in frame 2.
					mbdynce_tempframeM2_start_Wvel_loc = mbdynce_tempframe2G.TransformDirectionParentToLocal(mbdynce_tempframeG_start_Wvel_par);		//- expressed in frame 2.
					mbdynce_tempframeM2_start_Wvel_loc = mbdynce_tempframeM2_start.TransformDirectionParentToLocal(mbdynce_tempframeM2_start_Wvel_loc); //- expressed in frame loc
					mbdynce_accmotor_posdt3 = mbdynce_tempframeM2_start_pos_dt + time_step * mbdynce_tempframeM2_end_pos_dtdt;		
					mbdynce_accmotor_w_loc3 = mbdynce_tempframeM2_start_Wvel_loc + time_step * mbdynce_tempframeM2_end_Wacc_loc;*/
					ChVector<> mbdynce_tempframeG_start_Wvel_par;
					mbdynce_tempframeG_start_Wvel_par = motor3d_motor_i->GetBody1()->GetWvel_par(); //- Wvel expressed in C::E global ref.
					ChVector<> mbdynce_accmotor_posdt3_G, mbdynce_accmotor_w_par3;
					mbdynce_accmotor_posdt3_G = mbdynce_tempframeG_start_pos_dt + time_step * mbdynce_tempmbdyn_pos_dtdt;
					mbdynce_accmotor_w_par3 = mbdynce_tempframeG_start_Wvel_par + time_step * mbdynce_tempmbdyn_Wacc_par;
					mbdynce_accmotor_posdt3 = mbdynce_tempframe2G.TransformDirectionParentToLocal(mbdynce_accmotor_posdt3_G); //- expressed in frame 2.
					mbdynce_accmotor_w_loc3 = mbdynce_tempframe2G.TransformDirectionParentToLocal(mbdynce_accmotor_w_par3);		//- expressed in frame 2.
					mbdynce_accmotor_w_loc3 = mbdynce_tempframeM2_start.TransformDirectionParentToLocal(mbdynce_accmotor_w_loc3); //- expressed in frame loc.
					if (bMBDyn_CE_Verbose)
					{
						std::cout << "\t\tmbdynce_tempframeM2_start.GetPos():\t" << mbdynce_tempframeM2_start.GetPos() << "\n";
						std::cout << "\t\tmbdynce_tempframeM2_end.GetPos():\t" << mbdynce_tempframeM2_end.GetPos() << "\n";
						std::cout << "\t\tmbdynce_tempframeM2_start.GetRot()\t" << mbdynce_tempframeM2_start.GetRot() << "\n";
						std::cout << "\t\tmbdynce_tempframeM2_end.GetRot()\t" << mbdynce_tempframeM2_end.GetRot() << "\n";
						std::cout << "\t\tvel1: " << mbdynce_accmotor_posdt1 << "\t"<< mbdynce_accmotor_w_loc1 << "\n";
						std::cout << "\t\tmbdynce_tempmbdyn_pos_dt:\t" << mbdynce_tempmbdyn_pos_dt << "\n";
						std::cout << "\t\tmbdynce_tempframeM2_end_pos_dt :\t" << mbdynce_tempframeM2_end_pos_dt << "\n";
						std::cout << "\t\tmbdynce_tempmbdyn_Wvel_par:\t" << mbdynce_tempmbdyn_Wvel_par << "\n";
						std::cout << "\t\tmbdynce_tempframeM2_end_Wvel_loc:\t" << mbdynce_tempframeM2_end_Wvel_loc << "\n";
						std::cout << "\t\tvel2: " << mbdynce_accmotor_posdt2 << "\t" << mbdynce_accmotor_w_loc2 << "\n";
						std::cout << "\t\tmbdynce_tempframeG_start_pos_dt:\t" << mbdynce_tempframeG_start_pos_dt << "\n";
						//std::cout << "\t\tmbdynce_tempframeM2_start_pos_dt:\t" << mbdynce_tempframeM2_start_pos_dt << "\n";
						std::cout << "\t\tmbdynce_tempframeM2_end_pos_dtdt:\t" << mbdynce_tempframeM2_end_pos_dtdt << "\n";
						std::cout << "\t\tmbdynce_tempframeG_start_Wvel_par:\t" << mbdynce_tempframeG_start_Wvel_par << "\n";
						//std::cout << "\t\tmbdynce_tempframeM2_start_Wvel_loc:\t" << mbdynce_tempframeM2_start_Wvel_loc << "\n";
						std::cout << "\t\tmbdynce_tempframeM2_end_pos_dtdt:\t" << mbdynce_tempframeM2_end_Wacc_loc << "\n";
						std::cout << "\t\tvel3: " << mbdynce_accmotor_posdt3 << "\t" << mbdynce_accmotor_w_loc3 << "\n";
					}
					//- average velocit of velocity 1, 2 and 3
					mbdynce_tempframeM2_end_pos_dt = time_step * MBDyn_CE_CEMotorType_coeff[0] * mbdynce_accmotor_posdt1 + MBDyn_CE_CEMotorType_coeff[1] * mbdynce_accmotor_posdt2 + MBDyn_CE_CEMotorType_coeff[2] / time_step * mbdynce_accmotor_posdt3;
					mbdynce_tempframeM2_end_Wvel_loc = time_step * MBDyn_CE_CEMotorType_coeff[0] * mbdynce_accmotor_w_loc1 + MBDyn_CE_CEMotorType_coeff[1] * mbdynce_accmotor_w_loc2 + MBDyn_CE_CEMotorType_coeff[2] / time_step * mbdynce_accmotor_w_loc3;
					double coeff = time_step / (MBDyn_CE_CEMotorType_coeff[0] * time_step * time_step + MBDyn_CE_CEMotorType_coeff[1] * time_step + MBDyn_CE_CEMotorType_coeff[2]);
					mbdynce_tempframeM2_end_pos_dt *= coeff;
					mbdynce_tempframeM2_end_Wvel_loc *= coeff;

					//- position function using setpoint functions
					motor3d_function_pos->SetSetpointAndDerivatives(mbdynce_tempframeM2_start.GetPos(), mbdynce_tempframeM2_end_pos_dt, mbdynce_tempframeM2_end_pos_dtdt);
					motor3d_function_rot->SetSetpointAndDerivatives(mbdynce_tempframeM2_start.GetRot(), mbdynce_tempframeM2_end_Wvel_loc, mbdynce_tempframeM2_end_Wacc_loc);
				}
				else if (MBDyn_CE_CEMotorType == MBDyn_CE_CEMOTORTYPE::VELOCITY) //- velocity
				{
					//- when using setpoint function to impose velocity, 
					//- end position/rotation == start position/rotation.
					mbdynce_tempframeM2_end = mbdynce_tempframeM2_start;

					//- end linear veloctiy/acceleration;
					//- expressed in frame 2 (setpoint function);
					mbdynce_tempframeM2_end_pos_dt = mbdynce_tempframe2G.TransformDirectionParentToLocal(mbdynce_tempmbdyn_pos_dt);
					mbdynce_tempframeM2_end_pos_dtdt = mbdynce_tempframe2G.TransformDirectionParentToLocal(mbdynce_tempmbdyn_pos_dtdt);
					//- end angular veloctiy/acceleration;
					//- expressed in frame M (setpoint function)
					mbdynce_tempframeM2_end_Wvel_loc = mbdynce_tempframe2G.TransformDirectionParentToLocal(mbdynce_tempmbdyn_Wvel_par); //- angular velocity: expressed in frame 2
					mbdynce_tempframeM2_end_Wacc_loc = mbdynce_tempframe2G.TransformDirectionParentToLocal(mbdynce_tempmbdyn_Wacc_par);
					mbdynce_tempframeM2_end_Wvel_loc = mbdynce_tempframeM2_end.TransformDirectionParentToLocal(mbdynce_tempframeM2_end_Wvel_loc); //- angular velocity: expressed in frame M
					mbdynce_tempframeM2_end_Wacc_loc = mbdynce_tempframeM2_end.TransformDirectionParentToLocal(mbdynce_tempframeM2_end_Wacc_loc); //- should be expressed in the get_q(s) frame, namely, frame M.	

					// get the pointer to functions
					auto motor3d_function_pos = std::dynamic_pointer_cast<ChFunctionPosition_setpoint>(motor3d_function_pos_base);
					auto motor3d_function_rot = std::dynamic_pointer_cast<ChFunctionRotation_setpoint>(motor3d_function_rot_base);

					//- position function using setpoint functions
					motor3d_function_pos->SetSetpointAndDerivatives(mbdynce_tempframeM2_end.GetPos(), mbdynce_tempframeM2_end_pos_dt, mbdynce_tempframeM2_end_pos_dtdt);
					motor3d_function_rot->SetSetpointAndDerivatives(mbdynce_tempframeM2_end.GetRot(), mbdynce_tempframeM2_end_Wvel_loc, mbdynce_tempframeM2_end_Wacc_loc);
					//- std::cout << "\t\t\tpos_dt in frame 2 (pos_function): " << motor3d_function_pos->Get_p_ds(time)<< "\n";
					//- std::cout << "\t\t\trot_dt in frame M (rot_function): " << motor3d_function_rot->Get_w_loc(time)<< "\n";
				}
				else if (MBDyn_CE_CEMotorType == MBDyn_CE_CEMOTORTYPE::POSITION) //- position
				{
					auto motor3d_function_pos = std::dynamic_pointer_cast<ChFunctionPosition_line>(motor3d_function_pos_base);
					auto motor3d_function_rot = std::dynamic_pointer_cast<ChFunctionRotation_spline>(motor3d_function_rot_base);
					//- position function using line/spline interpolation
					auto mbdynce_temp_pos_line = chrono_types::make_shared<geometry::ChLineSegment>(mbdynce_tempframeM2_start.GetPos(), mbdynce_tempframeM2_end.GetPos());
					motor3d_function_pos->SetLine(mbdynce_temp_pos_line);
					//- chrono_types::make_shared<>: a more safety case in C::E
					motor3d_function_pos->SetSpaceFunction(chrono_types::make_shared<ChFunction_Ramp>(-time / time_step, 1 / time_step));
					//- rotation function
					if (mbdynce_tempframeM2_start.GetRot().Dot(mbdynce_tempframeM2_end.GetRot())<0)
					{
						mbdynce_tempframeM2_end.GetRot() = -mbdynce_tempframeM2_end.GetRot();
					}
					std::vector<ChQuaternion<>> mbdynce_temp_rot_spline = {{mbdynce_tempframeM2_start.GetRot()}, {mbdynce_tempframeM2_end.GetRot()}};
					motor3d_function_rot->SetupData(1, mbdynce_temp_rot_spline);
					motor3d_function_rot->SetSpaceFunction(chrono_types::make_shared<ChFunction_Ramp>(-time / time_step, 1 / time_step));
				}	
				else if (MBDyn_CE_CEMotorType == MBDyn_CE_CEMOTORTYPE::ORIGIN) //- origin
				{
					//- when using setpoint function to impose origin, 
					//- end position/rotation == position/rotation (last step).

					//- rotation: quanternion at the same direction
					if (mbdynce_tempframeM2_start.GetRot().Dot(mbdynce_tempframeM2_end.GetRot())<0)
					{
						mbdynce_tempframeM2_end.GetRot() = -mbdynce_tempframeM2_end.GetRot();
					}

					//- end linear veloctiy/acceleration;
					//- expressed in frame 2 (setpoint function);
					mbdynce_tempframeM2_end_pos_dt = mbdynce_tempframe2G.TransformDirectionParentToLocal(mbdynce_tempmbdyn_pos_dt);
					mbdynce_tempframeM2_end_pos_dtdt = mbdynce_tempframe2G.TransformDirectionParentToLocal(mbdynce_tempmbdyn_pos_dtdt);
					//- end angular veloctiy/acceleration;
					//- expressed in frame M (setpoint function)
					mbdynce_tempframeM2_end_Wvel_loc = mbdynce_tempframe2G.TransformDirectionParentToLocal(mbdynce_tempmbdyn_Wvel_par); //- angular velocity: expressed in frame 2
					mbdynce_tempframeM2_end_Wacc_loc = mbdynce_tempframe2G.TransformDirectionParentToLocal(mbdynce_tempmbdyn_Wacc_par);
					mbdynce_tempframeM2_end_Wvel_loc = mbdynce_tempframeM2_end.TransformDirectionParentToLocal(mbdynce_tempframeM2_end_Wvel_loc); //- angular velocity: expressed in frame M
					mbdynce_tempframeM2_end_Wacc_loc = mbdynce_tempframeM2_end.TransformDirectionParentToLocal(mbdynce_tempframeM2_end_Wacc_loc); //- should be expressed in the get_q(s) frame, namely, frame M.	

					// get the pointer to functions
					auto motor3d_function_pos = std::dynamic_pointer_cast<ChFunctionPosition_setpoint>(motor3d_function_pos_base);
					auto motor3d_function_rot = std::dynamic_pointer_cast<ChFunctionRotation_setpoint>(motor3d_function_rot_base);

					//- position function using setpoint functions
					motor3d_function_pos->SetSetpointAndDerivatives(mbdynce_tempframeM2_end.GetPos(), mbdynce_tempframeM2_end_pos_dt, mbdynce_tempframeM2_end_pos_dtdt);
					motor3d_function_rot->SetSetpointAndDerivatives(mbdynce_tempframeM2_end.GetRot(), mbdynce_tempframeM2_end_Wvel_loc, mbdynce_tempframeM2_end_Wacc_loc);
					//- std::cout << "\t\t\tpos_dt in frame 2 (pos_function): " << motor3d_function_pos->Get_p_ds(time)<< "\n";
					//- std::cout << "\t\t\trot_dt in frame M (rot_function): " << motor3d_function_rot->Get_w_loc(time)<< "\n";
				}
				if(bMBDyn_CE_Verbose)
				{
					std::cout << "\t\tC::E model motor " << motor3d_motor_i_id << " function start at (expressed in Frame 2): \n";
					//std::cout << "\t\t\tpos in C::E: " << mbdynce_tempframeM2_start.GetPos() << "\n";
					//std::cout << "\t\t\treceived from MBDyn: " << mbdynce_tempframeM2_end.GetPos() << "\n";
					std::cout << "\t\t\tpos: " << motor3d_motor_i->GetPositionFunction()->Get_p(time) << "\n";
					std::cout << "\t\t\trot: " << motor3d_motor_i->GetRotationFunction()->Get_q(time) << "\n";
					std::cout << "\t\t\tpos_dt: " << motor3d_motor_i->GetPositionFunction()->Get_p_ds(time) << "\n";
					std::cout << "\t\t\tWvel_loc: " << motor3d_motor_i->GetRotationFunction()->Get_w_loc(time) << "\n";
					std::cout << "\t\tC::E model motor " << motor3d_motor_i_id << " function end at (expressed in Frame 2): \n";
					std::cout << "\t\t\tpos: " << motor3d_motor_i->GetPositionFunction()->Get_p(time + time_step) << "\n";
					std::cout << "\t\t\trot: " << motor3d_motor_i->GetRotationFunction()->Get_q(time + time_step) << "\n";
					std::cout << "\t\t\tpos_dt: " << motor3d_motor_i->GetPositionFunction()->Get_p_ds(time + time_step) << "\n";
					std::cout << "\t\t\tWvel_loc: " << motor3d_motor_i->GetRotationFunction()->Get_w_loc(time + time_step) << "\n";
				}
			}
			else
			{
				std::cout << "Error: functions of motor3d_motor_i" << motor3d_motor_i_id << " are null. \n";
				return 1;
			}
		}
		else
		{
			std::cout << "Error: motor3d_motor_i" << motor3d_motor_i_id << " is null. \n";
			return 1;
		}
	}
	return 0;
}

//- C::E models send coupling forces to the buffer
extern "C" int 
MBDyn_CE_CEModel_SendToBuf(pMBDyn_CE_CEModel_t pMBDyn_CE_CEModel, std::vector<double> &MBDyn_CE_CouplingDynamic, 
							double* pMBDyn_CE_CEFrame, const unsigned& MBDyn_CE_NodesNum, const double* MBDyn_CE_CEScale,
							const std::vector<MBDYN_CE_CEMODELDATA> & MBDyn_CE_CEModel_Label,
							const int MBDyn_CE_CEForceType, 
							bool bMBDyn_CE_Verbose)
{
	if (pMBDyn_CE_CEModel==NULL)
	{
		std::cout << "\t\tCE_models SendToMBDyn() fails:\n";
		return 1;
	}
	//- std::cout << "\t\tCE_models SendToMBDyn():\n";
	ChSystemParallelNSC *tempsys = (ChSystemParallelNSC *)pMBDyn_CE_CEModel;
	//- obtain the transform matrix
	ChVector<> mbdynce_temp_frameMBDyn_pos(-pMBDyn_CE_CEFrame[0], -pMBDyn_CE_CEFrame[1], -pMBDyn_CE_CEFrame[2]);
	ChVector<> mbdynce_temp_frameMBDyn_rot_axis_X(pMBDyn_CE_CEFrame[3], pMBDyn_CE_CEFrame[6], pMBDyn_CE_CEFrame[9]);
	ChVector<> mbdynce_temp_frameMBDyn_rot_axis_Y(pMBDyn_CE_CEFrame[4], pMBDyn_CE_CEFrame[7], pMBDyn_CE_CEFrame[10]);
	ChVector<> mbdynce_temp_frameMBDyn_rot_axis_Z(pMBDyn_CE_CEFrame[5], pMBDyn_CE_CEFrame[8], pMBDyn_CE_CEFrame[11]);
	ChMatrix33<> mbdynce_temp_frameMBDyn_rot(mbdynce_temp_frameMBDyn_rot_axis_X,mbdynce_temp_frameMBDyn_rot_axis_Y,mbdynce_temp_frameMBDyn_rot_axis_Z);
	ChFrame<> mbdynce_temp_frameMBDyn(mbdynce_temp_frameMBDyn_pos, mbdynce_temp_frameMBDyn_rot);
	//- write exchange force/torque to the buffer
	ChVector<> mbdynce_ce_force_G(0.0, 0.0, 0.0);
	ChVector<> mbdynce_ce_torque_G(0.0, 0.0, 0.0);
	ChVector<> mbdynce_mbdyn_force = mbdynce_temp_frameMBDyn.TransformDirectionParentToLocal(mbdynce_ce_force_G);
	ChVector<> mbdynce_mbdyn_torque = mbdynce_temp_frameMBDyn.TransformDirectionParentToLocal(mbdynce_ce_force_G);
	double* pmbdynce_tempvec3_f = &MBDyn_CE_CouplingDynamic[0];
	double* pmbdynce_tempvec3_m = &MBDyn_CE_CouplingDynamic[3*MBDyn_CE_NodesNum];
	for (unsigned i = 0; i < MBDyn_CE_NodesNum; i++)
	{
		//- obtain the exchange force/torque from CE model; 
		unsigned motor_i_id = MBDyn_CE_CEModel_Label[i].MBDyn_CE_CEMotor_Label;
		auto motor_i = std::dynamic_pointer_cast<ChLinkMotionImposed>(tempsys->SearchLinkID(motor_i_id));
		if (motor_i == NULL)
		{
			std::cout << "\t\tcannot read coupling motor\n";
			return 1;
		}
		else
		{
			if (MBDyn_CE_CEForceType==MBDyn_CE_CEFORCETYPE::REACTION_FORCE)
			{
				ChCoordsys<> mbdynce_temp_frameMG_calcu(motor_i->GetLinkAbsoluteCoords());
				//- get the reaction force at the global frame and MBDyn frame, apply to the node(body 2)
				mbdynce_ce_force_G = mbdynce_temp_frameMG_calcu.TransformDirectionLocalToParent(motor_i->Get_react_force());
				mbdynce_ce_torque_G = mbdynce_temp_frameMG_calcu.TransformDirectionLocalToParent(motor_i->Get_react_torque());
				mbdynce_mbdyn_force = mbdynce_temp_frameMBDyn.TransformDirectionParentToLocal(mbdynce_ce_force_G);
				mbdynce_mbdyn_torque = mbdynce_temp_frameMBDyn.TransformDirectionParentToLocal(mbdynce_ce_torque_G);
			}
			else if (MBDyn_CE_CEForceType==MBDyn_CE_CEFORCETYPE::CONTACT_FORCE)
			{
				//- get the contact force
				ChCoordsys<> mbdynce_temp_frameMG_calcu(motor_i->GetLinkAbsoluteCoords());
				auto temp_body = tempsys->SearchBodyID(MBDyn_CE_CEModel_Label[i].MBDyn_CE_CEBody_Label);
				ChVector<> contact_forc(tempsys->GetBodyContactForce(temp_body).x, tempsys->GetBodyContactForce(temp_body).y, tempsys->GetBodyContactForce(temp_body).z);
				ChVector<> contact_torq(tempsys->GetBodyContactTorque(temp_body).x, tempsys->GetBodyContactTorque(temp_body).y, tempsys->GetBodyContactTorque(temp_body).z);
				mbdynce_ce_force_G=mbdynce_temp_frameMG_calcu.TransformDirectionLocalToParent(contact_forc);
				mbdynce_ce_torque_G=mbdynce_temp_frameMG_calcu.TransformDirectionLocalToParent(contact_torq);
				mbdynce_mbdyn_force = mbdynce_temp_frameMBDyn.TransformDirectionParentToLocal(mbdynce_ce_force_G);
				mbdynce_mbdyn_force = mbdynce_temp_frameMBDyn.TransformDirectionParentToLocal(mbdynce_ce_torque_G);
			}
		}
		//- write in the buffer
		double mbdynce_tempvec3_f[3]={mbdynce_mbdyn_force.x()/MBDyn_CE_CEScale[2],
										mbdynce_mbdyn_force.y()/MBDyn_CE_CEScale[2],
										mbdynce_mbdyn_force.z()/MBDyn_CE_CEScale[2]};
		double mbdynce_tempvec3_m[3]={mbdynce_mbdyn_torque.x()/(MBDyn_CE_CEScale[3]),
										mbdynce_mbdyn_torque.y()/(MBDyn_CE_CEScale[3]),
										mbdynce_mbdyn_torque.z()/(MBDyn_CE_CEScale[3])};
		memcpy(&pmbdynce_tempvec3_f[3*i], &mbdynce_tempvec3_f[0], 3 * sizeof(double));
		memcpy(&pmbdynce_tempvec3_m[3*i], &mbdynce_tempvec3_m[0], 3 * sizeof(double));
		if (bMBDyn_CE_Verbose)
		{
			std::cout << "\t\tcoupling forces i " << motor_i_id << "\n";
			std::cout << "\t\t\tcoupling forces_f: " << mbdynce_tempvec3_f[0] << "\t" << mbdynce_tempvec3_f[1] << "\t" << mbdynce_tempvec3_f[2] << "\n";
			std::cout << "\t\t\tcoupling forces_q: " << mbdynce_tempvec3_m[0] << "\t" << mbdynce_tempvec3_m[1] << "\t" << mbdynce_tempvec3_m[2] << "\n";
		}
	}
	return 0;
}

//- write data to files (.usr)
extern "C" int 
MBDyn_CE_CEModel_WriteToFiles(pMBDyn_CE_CEModel_t pMBDyn_CE_CEModel,
					const std::vector<MBDYN_CE_CEMODELDATA> & MBDyn_CE_CEModel_Label,
					double *pMBDyn_CE_CEFrame,  
					const double* MBDyn_CE_CEScale,
					std::ostream & out)
{
	if (pMBDyn_CE_CEModel==NULL)
	{
		std::cout << "Error: the C::E model pointer is NULL.\n";
		return 1;
	}
	ChSystemParallelNSC *tempsys = (ChSystemParallelNSC *)pMBDyn_CE_CEModel; // static_cast?
	//- std::cout << "\t\tCE_models DataSave():\n";
	unsigned int tempsys_couplingbodies_size = MBDyn_CE_CEModel_Label.size();
	double time = tempsys->GetChTime();
	//- obtain the transform matrix
	ChVector<> mbdynce_temp_frameMBDyn_pos(-pMBDyn_CE_CEFrame[0], -pMBDyn_CE_CEFrame[1], -pMBDyn_CE_CEFrame[2]);
	ChVector<> mbdynce_temp_frameMBDyn_rot_axis_X(pMBDyn_CE_CEFrame[3], pMBDyn_CE_CEFrame[6], pMBDyn_CE_CEFrame[9]);
	ChVector<> mbdynce_temp_frameMBDyn_rot_axis_Y(pMBDyn_CE_CEFrame[4], pMBDyn_CE_CEFrame[7], pMBDyn_CE_CEFrame[10]);
	ChVector<> mbdynce_temp_frameMBDyn_rot_axis_Z(pMBDyn_CE_CEFrame[5], pMBDyn_CE_CEFrame[8], pMBDyn_CE_CEFrame[11]);
	ChMatrix33<> mbdynce_temp_frameMBDyn_rot(mbdynce_temp_frameMBDyn_rot_axis_X,mbdynce_temp_frameMBDyn_rot_axis_Y,mbdynce_temp_frameMBDyn_rot_axis_Z);
	ChFrame<> mbdynce_temp_frameMBDyn(mbdynce_temp_frameMBDyn_pos, mbdynce_temp_frameMBDyn_rot);
	//- write data to files
//- #pragma omp parallel for
	/*for (unsigned int i = 0; i < tempsys_couplingbodies_size; i++)
	{
		if (MBDyn_CE_CEModel_Label[i].bMBDyn_CE_CEBody_Output)
		{
			auto temp_body = tempsys->SearchBodyID(MBDyn_CE_CEModel_Label[i].MBDyn_CE_CEBody_Label);
			const ChVector<>& body_pos = temp_body->GetPos(); // 3
			const ChQuaternion<> &body_rot = temp_body->GetRot();   // 4
			const ChVector<> &body_pos_dt = temp_body->GetPos_dt(); // 3
			const ChVector<> &body_Wvel_par = temp_body->GetWvel_par(); // 3
			const ChVector<> &body_pos_dtdt = temp_body->GetPos_dtdt(); // 3
			const ChVector<> &body_Wacc_par = temp_body->GetWacc_par(); // 3

			ChVector<> out_body_pos = mbdynce_temp_frameMBDyn.TransformParentToLocal(body_pos);
			ChQuaternion<> out_body_rot_temp = body_rot >> (mbdynce_temp_frameMBDyn.GetInverse());
			ChVector<> out_body_rot = out_body_rot_temp.Q_to_Euler123();
			ChVector<> out_body_pos_dt = mbdynce_temp_frameMBDyn.TransformDirectionParentToLocal(body_pos_dt);
			ChVector<> out_body_Wvel_par = mbdynce_temp_frameMBDyn.TransformDirectionParentToLocal(body_Wvel_par);
			ChVector<> out_body_pos_dtdt = mbdynce_temp_frameMBDyn.TransformDirectionParentToLocal(body_pos_dtdt);
			ChVector<> out_body_Wacc_par = mbdynce_temp_frameMBDyn.TransformDirectionParentToLocal(body_Wacc_par);
			double degtorad = 180.0 / CH_C_PI;

			out << MBDyn_CE_CEModel_Label[i].MBDyn_CE_CEBody_Label << " "
				<< out_body_pos.x() / MBDyn_CE_CEScale[0] << " " << out_body_pos.y() / MBDyn_CE_CEScale[0] << " " << out_body_pos.z() / MBDyn_CE_CEScale[0] << " "
				<< out_body_rot.z() * degtorad << " " << out_body_rot.y() * degtorad << " " << out_body_rot.x() * degtorad << " "
				<< out_body_pos_dt.x() / MBDyn_CE_CEScale[0] << " " << out_body_pos_dt.y() / MBDyn_CE_CEScale[0] << " " << out_body_pos_dt.z() / MBDyn_CE_CEScale[0] << " "
				<< out_body_Wvel_par.x() << " " << out_body_Wvel_par.y() << " " << out_body_Wvel_par.z() << " "
				<< out_body_pos_dtdt.x() / MBDyn_CE_CEScale[0] << " " << out_body_pos_dtdt.y() / MBDyn_CE_CEScale[0] << " " << out_body_pos_dtdt.z() / MBDyn_CE_CEScale[0] << " "
				<< out_body_Wacc_par.x() << " " << out_body_Wacc_par.y() << " " << out_body_Wacc_par.z() << "\n";
		}
	}*/
	// output all bodies in C::E
	for (unsigned i = 0; i < tempsys->Get_bodylist().size(); i++)
	{
		auto temp_body = tempsys->Get_bodylist()[i];
		const ChVector<> &body_pos = temp_body->GetPos();			// 3
		const ChQuaternion<> &body_rot = temp_body->GetRot();		// 4
		const ChVector<> &body_pos_dt = temp_body->GetPos_dt();		// 3
		const ChVector<> &body_Wvel_par = temp_body->GetWvel_par(); // 3
		const ChVector<> &body_pos_dtdt = temp_body->GetPos_dtdt(); // 3
		const ChVector<> &body_Wacc_par = temp_body->GetWacc_par(); // 3

		ChVector<> out_body_pos = mbdynce_temp_frameMBDyn.TransformParentToLocal(body_pos);
		ChQuaternion<> out_body_rot_temp = body_rot >> (mbdynce_temp_frameMBDyn.GetInverse());
		ChVector<> out_body_rot = out_body_rot_temp.Q_to_Euler123();
		ChVector<> out_body_pos_dt = mbdynce_temp_frameMBDyn.TransformDirectionParentToLocal(body_pos_dt);
		ChVector<> out_body_Wvel_par = mbdynce_temp_frameMBDyn.TransformDirectionParentToLocal(body_Wvel_par);
		ChVector<> out_body_pos_dtdt = mbdynce_temp_frameMBDyn.TransformDirectionParentToLocal(body_pos_dtdt);
		ChVector<> out_body_Wacc_par = mbdynce_temp_frameMBDyn.TransformDirectionParentToLocal(body_Wacc_par);
		double degtorad = 180.0 / CH_C_PI;

		out << temp_body->GetIdentifier() << " "
			<< out_body_pos.x() / MBDyn_CE_CEScale[0] << " " << out_body_pos.y() / MBDyn_CE_CEScale[0] << " " << out_body_pos.z() / MBDyn_CE_CEScale[0] << " "
			<< out_body_rot.z() * degtorad << " " << out_body_rot.y() * degtorad << " " << out_body_rot.x() * degtorad << " "
			<< out_body_pos_dt.x() / MBDyn_CE_CEScale[0] << " " << out_body_pos_dt.y() / MBDyn_CE_CEScale[0] << " " << out_body_pos_dt.z() / MBDyn_CE_CEScale[0] << " "
			<< out_body_Wvel_par.x() << " " << out_body_Wvel_par.y() << " " << out_body_Wvel_par.z() << " "
			<< out_body_pos_dtdt.x() / MBDyn_CE_CEScale[0] << " " << out_body_pos_dtdt.y() / MBDyn_CE_CEScale[0] << " " << out_body_pos_dtdt.z() / MBDyn_CE_CEScale[0] << " "
			<< out_body_Wacc_par.x() << " " << out_body_Wacc_par.y() << " " << out_body_Wacc_par.z() << "\n";
	}
	return 0;
}
