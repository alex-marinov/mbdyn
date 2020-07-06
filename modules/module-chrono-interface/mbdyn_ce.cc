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

#include "mbdyn_ce.h"
#include "chrono/ChConfig.h"
#include "chrono_parallel/ChDataManager.h" // for simulation of parallel system, data_manager
#include "chrono_parallel/solver/ChIterativeSolverParallel.h"
#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono/physics/ChLinkMotionImposed.h" //for 3-D dimension
#include "chrono/motion_functions/ChFunctionPosition_line.h"
#include "chrono/motion_functions/ChFunctionRotation_spline.h"

using namespace chrono;
using namespace chrono::collision;

extern "C" void
MBDyn_CE_CEModel_Create(ChSystemParallelNSC *pMBDyn_CE_CEModel);

extern "C" pMBDyn_CE_CEModel_t MBDyn_CE_CEModel_Init
(std::vector<double> & MBDyn_CE_CEModel_Data,
const double* pMBDyn_CE_CEFrame, const double* MBDyn_CE_CEScale,
std::vector<MBDYN_CE_CEMODELDATA> & MBDyn_CE_CEModel_Label,
const int& MBDyn_CE_CouplingType)
{
	std::cout << "Initial MBDyn_CE_CEModel pointer:\n";
	ChSystemParallelNSC *pMBDyn_CE_CEModel = new ChSystemParallelNSC;
	if(pMBDyn_CE_CEModel==NULL)
	{
		std::cout << "\t\tInitial MBDyn_CE_CEModel pointer fails\n";
	}
	MBDyn_CE_CEModel_Create(pMBDyn_CE_CEModel);

	// initial the ground coordinate in C::E;
	// r = (-CEF1,-CEF2,-CEF3);
	// R= [CEF3,  CEF4,  CEF5 ]^T           [CEF3,  CEF6,  CEF9]
	//    [CEF6,  CEF7,  CEF8 ]     ====    [CEF4,  CEF7,  CEF10]
	//    [CEF9,  CEF10, CEF11];            [CEF5,  CEF8,  CEF11]
	ChVector<> mbdynce_temp_frameMBDyn_pos(-pMBDyn_CE_CEFrame[0], -pMBDyn_CE_CEFrame[1], -pMBDyn_CE_CEFrame[2]);
	ChVector<> mbdynce_temp_frameMBDyn_rot_axis_X(pMBDyn_CE_CEFrame[3], pMBDyn_CE_CEFrame[4], pMBDyn_CE_CEFrame[5]);
	ChVector<> mbdynce_temp_frameMBDyn_rot_axis_Y(pMBDyn_CE_CEFrame[6], pMBDyn_CE_CEFrame[7], pMBDyn_CE_CEFrame[8]);
	ChVector<> mbdynce_temp_frameMBDyn_rot_axis_Z(pMBDyn_CE_CEFrame[9], pMBDyn_CE_CEFrame[10], pMBDyn_CE_CEFrame[11]);
	ChMatrix33<> mbdynce_temp_frameMBDyn_rot(mbdynce_temp_frameMBDyn_rot_axis_X,mbdynce_temp_frameMBDyn_rot_axis_Y,mbdynce_temp_frameMBDyn_rot_axis_Z);
	ChFrame<> mbdynce_temp_frameMBDyn(mbdynce_temp_frameMBDyn_pos, mbdynce_temp_frameMBDyn_rot);
	// initial the gravity of C::E model
	pMBDyn_CE_CEModel->Set_G_acc(mbdynce_temp_frameMBDyn.TransformDirectionLocalToParent(ChVector<>(0.0,-9.81*MBDyn_CE_CEScale[0],0.0)));

	// initial motor for coupling bodies.
	if (MBDyn_CE_CouplingType>=-2)//coupling
	{
		unsigned mbdynce_temp_bodies_num = MBDyn_CE_CEModel_Label.size() - 1;
		unsigned mbdynce_temp_ground_id = MBDyn_CE_CEModel_Label[mbdynce_temp_bodies_num].MBDyn_CE_CEBody_Label; // ground ID is set in the last element
		auto mbdynce_temp_ground = pMBDyn_CE_CEModel->SearchBodyID(mbdynce_temp_ground_id);
		for (unsigned i = 0; i < mbdynce_temp_bodies_num;i++)
		{
			unsigned body_i_id = MBDyn_CE_CEModel_Label[i].MBDyn_CE_CEBody_Label;
			auto body_i = pMBDyn_CE_CEModel->SearchBodyID(body_i_id);
			auto motor3d_body_i = std::make_shared<ChLinkMotionImposed>();
			motor3d_body_i->Initialize(mbdynce_temp_ground,
									   body_i,
									   false, //connecting frames are described in ground ref.
									   ChFrame<>(ChVector<>(0.0, 0.0, 0.0)), // how to determine these two frames::TO DO !!!!!!
									   ChFrame<>(ChVector<>(0.0, 0.0, 0.0)));
			pMBDyn_CE_CEModel->Add(motor3d_body_i);
			auto motor3d_function_pos = std::make_shared<ChFunctionPosition_line>(); // impose veloctiy:: TO DO !!!!!!!
			auto motor3d_function_rot = std::make_shared<ChFunctionRotation_spline>();
			motor3d_body_i->SetPositionFunction(motor3d_function_pos);
			motor3d_body_i->SetRotationFunction(motor3d_function_rot);
			MBDyn_CE_CEModel_Label[i].MBDyn_CE_CEMotor_Label = motor3d_body_i->GetIdentifier();
			std::cout << "C::E motor" << i + 1 << "ID:\t" << motor3d_body_i->GetIdentifier() << "\n";
		}
		std::cout << "C::E ground ID is \t"<<mbdynce_temp_ground->GetIdentifier()<<"\n";
	}
	else
	{
		std::cout << "Coupling none in C::E model.\n";
	}
	// allocate space for C::E_Model_Data;
	unsigned int bodies_size = pMBDyn_CE_CEModel->Get_bodylist().size();
	unsigned int system_size = (3 * 3 + 3 * 4) * bodies_size + 1; // +1: save Chtime; system_size=body_size+1(time)
	MBDyn_CE_CEModel_Data.resize(system_size,0.0);
	return pMBDyn_CE_CEModel;
}

extern "C" void
MBDyn_CE_CEModel_Destroy(pMBDyn_CE_CEModel_t pMBDyn_CE_CEModel)
{
	std::cout << "destroy the CE_model...\n";
	if(pMBDyn_CE_CEModel!=NULL)
	{
		// must convert to the correct type
		// delete  (int *) pMBDyn_CE_CEModel;
		delete  (ChSystemParallelNSC *) pMBDyn_CE_CEModel;
		pMBDyn_CE_CEModel = NULL;
	}
}

// save CEModel at current step for reloading them in the tight coupling scheme
// (before advance())
extern "C" int
MBDyn_CE_CEModel_DataSave(pMBDyn_CE_CEModel_t pMBDyn_CE_CEModel, 
                        std::vector<double> & MBDyn_CE_CEModel_Data)
{
	ChSystemParallelNSC *tempsys = (ChSystemParallelNSC *)pMBDyn_CE_CEModel;
	unsigned int tempsys_bodies_size = tempsys->Get_bodylist().size();
	unsigned int tempsys_size = (3 * 3 + 3 * 4) * tempsys_bodies_size + 1; // +1 Chtime: Sys_size=body_size + 1(for time);
	unsigned int vector_size = MBDyn_CE_CEModel_Data.size();
	if (tempsys_size != vector_size)
	{
		std::cout << "Error: the vector to save data is not consistent with the C::E model:\n"; // how to safely exit MBDyn?
		return 1;
	}

	// save data
#pragma omp parallel for
	for (unsigned int i = 0; i < tempsys_bodies_size; i++)
	{
		ChVector<>& body_pos = tempsys->Get_bodylist()[i]->GetPos(); // 3
		ChQuaternion<> &body_rot = tempsys->Get_bodylist()[i]->GetRot();   // 4
		ChVector<> &body_pos_dt = tempsys->Get_bodylist()[i]->GetPos_dt(); // 3
		ChQuaternion<> &body_rot_dt = tempsys->Get_bodylist()[i]->GetRot_dt(); // 4
		ChVector<> &body_pos_dtdt = tempsys->Get_bodylist()[i]->GetPos_dtdt(); // 3
		ChQuaternion<> &body_rot_dtdt = tempsys->Get_bodylist()[i]->GetRot_dtdt(); //4

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


// reload data in the tight coupling scheme at each iteration
extern "C" int
MBDyn_CE_CEModel_DataReload(pMBDyn_CE_CEModel_t pMBDyn_CE_CEModel, 
                        std::vector<double> & MBDyn_CE_CEModel_Data)
{
	ChSystemParallelNSC *tempsys = (ChSystemParallelNSC *)pMBDyn_CE_CEModel;
	unsigned int tempsys_bodies_size = tempsys->Get_bodylist().size();
	unsigned int tempsys_size = (3 * 3 + 3 * 4) * tempsys_bodies_size + 1; // +1 Chtime: Sys_size=body_size + 1(for time);
	unsigned int vector_size = MBDyn_CE_CEModel_Data.size();
	if (tempsys_size != vector_size)
	{
		std::cout << "Error: the vector to save data is not consistent with the C::E model:\n"; // how to safely exit MBDyn?
		return 1;
	}
	double MBDyn_CE_CEModel_time = MBDyn_CE_CEModel_Data[tempsys_size - 1];
	tempsys->SetChTime(MBDyn_CE_CEModel_time);
#pragma omp parallel for
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
		tempsys->Get_bodylist()[i]->SetRot(ChQuaternion<>(MBDyn_CE_CEModel_Data[i_rot_dt],
														  MBDyn_CE_CEModel_Data[i_rot_dt + 1],
														  MBDyn_CE_CEModel_Data[i_rot_dt + 2],
														  MBDyn_CE_CEModel_Data[i_rot_dt + 3]));
		
		unsigned int i_pos_dtdt = i_rot_dt + 4;
		tempsys->Get_bodylist()[i]->SetPos_dt(ChVector<>(MBDyn_CE_CEModel_Data[i_pos_dtdt],
													  MBDyn_CE_CEModel_Data[i_pos_dtdt + 1],
													  MBDyn_CE_CEModel_Data[i_pos_dtdt + 2]));

		unsigned int i_rot_dtdt = i_pos_dtdt + 3;
		tempsys->Get_bodylist()[i]->SetRot(ChQuaternion<>(MBDyn_CE_CEModel_Data[i_rot_dtdt],
														  MBDyn_CE_CEModel_Data[i_rot_dtdt + 1],
														  MBDyn_CE_CEModel_Data[i_rot_dtdt + 2],
														  MBDyn_CE_CEModel_Data[i_rot_dtdt + 3]));
		tempsys->Get_bodylist()[i]->Update(MBDyn_CE_CEModel_time);
	}
	return 0;
}

extern "C" void
MBDyn_CE_CEModel_DoStepDynamics(pMBDyn_CE_CEModel_t pMBDyn_CE_CEModel, double time_step)
{
	std::cout << "\t\tCE_models DoStepDynamics():\n";
	if (pMBDyn_CE_CEModel!=NULL)
	{   // it's not a good idea to do this convert
		// it's dangerous, how to avoid?
		ChSystemParallelNSC* tempsys=(ChSystemParallelNSC *) pMBDyn_CE_CEModel;
		tempsys->DoStepDynamics(time_step);
		tempsys->CalculateContactForces();
		std::cout << "\t\t" << tempsys->GetChTime() << "\n";
		//pedantic_cout();
	}
}

// C::E models receive coupling motion from the buffer
void 
MBDyn_CE_CEModel_RecvFromBuf(pMBDyn_CE_CEModel_t pMBDyn_CE_CEModel, std::vector<double>& MBDyn_CE_CouplingKinematic, const unsigned& MBDyn_CE_NodesNum)
{
	ChSystemParallelNSC *mbdynce_tempsys = (ChSystemParallelNSC *)pMBDyn_CE_CEModel;
	// 1. obtain the data;
	// 2. transfer it to the coordinate in C::E;
	// 3. apply to the model;
	
	// 1. obtain the data;
	double *pmbdynce_tempvec3_x = &MBDyn_CE_CouplingKinematic[0];
	double *pmbdynce_tempemat3x3_R = &MBDyn_CE_CouplingKinematic[3 * MBDyn_CE_NodesNum];
	double *pmbdynce_tempvec3_xp = &MBDyn_CE_CouplingKinematic[12 * MBDyn_CE_NodesNum];
	double *pmbdynce_tempvec3_omega = &MBDyn_CE_CouplingKinematic[15 * MBDyn_CE_NodesNum];
	double *pmbdynce_tempvec3_xpp = &MBDyn_CE_CouplingKinematic[18 * MBDyn_CE_NodesNum];
	double *pmbdynce_tempvec3_omegap = &MBDyn_CE_CouplingKinematic[21 * MBDyn_CE_NodesNum];
	double *pmbdynce_temp_frame = &MBDyn_CE_CouplingKinematic[24 * MBDyn_CE_NodesNum];
	ChVector<> mbdynce_temp_frameMBDyn_pos(-(*pmbdynce_temp_frame), -(*(pmbdynce_temp_frame + 1)), -(*(pmbdynce_temp_frame + 2)));
	ChVector<> mbdynce_temp_frameMBDyn_rot_axis_X(*(pmbdynce_temp_frame + 3), *(pmbdynce_temp_frame + 4), *(pmbdynce_temp_frame + 5));
	ChVector<> mbdynce_temp_frameMBDyn_rot_axis_Y(*(pmbdynce_temp_frame + 6), *(pmbdynce_temp_frame + 7), *(pmbdynce_temp_frame + 8));
	ChVector<> mbdynce_temp_frameMBDyn_rot_axis_Z(*(pmbdynce_temp_frame + 9), *(pmbdynce_temp_frame + 10), *(pmbdynce_temp_frame + 11));
	ChMatrix33<> mbdynce_temp_frameMBDyn_rot(mbdynce_temp_frameMBDyn_rot_axis_X,mbdynce_temp_frameMBDyn_rot_axis_Y,mbdynce_temp_frameMBDyn_rot_axis_Z);
	ChFrame<> mbdynce_temp_frameMBDyn(mbdynce_temp_frameMBDyn_pos, mbdynce_temp_frameMBDyn_rot);

	// 2. transfer it to the coordinate in C::E, and create motor functions
	for (unsigned i = 0; i < MBDyn_CE_NodesNum;i++)
	{
		// 2.1 coordinate transformation
		// Currently, the codes only use the pos and rotation.
		ChVector<> mbdynce_tempmbdyn_pos = ChVector<>(pmbdynce_tempvec3_x[3 * i], pmbdynce_tempvec3_x[3 * i + 1], pmbdynce_tempvec3_x[3 * i + 2]) >> mbdynce_temp_frameMBDyn;
		ChMatrix33<> mbdynce_tempmbdyn_R1(ChVector<>(pmbdynce_tempemat3x3_R[9 * i], pmbdynce_tempemat3x3_R[9 * i + 3], pmbdynce_tempemat3x3_R[9 * i + 6]),
										  ChVector<>(pmbdynce_tempemat3x3_R[9 * i + 1], pmbdynce_tempemat3x3_R[9 * i + 4], pmbdynce_tempemat3x3_R[9 * i + 7]),
										  ChVector<>(pmbdynce_tempemat3x3_R[9 * i + 2], pmbdynce_tempemat3x3_R[9 * i + 5], pmbdynce_tempemat3x3_R[9 * i + 8])); // three column vectors
		ChMatrix33<> mbdynce_tempmbdyn_R(mbdynce_tempmbdyn_R1.Get_A_quaternion() >> mbdynce_temp_frameMBDyn);
		ChFrame<> mbdynce_tempframe_end(mbdynce_tempmbdyn_pos,mbdynce_tempmbdyn_R);

		// 2.2 create motor functions
		// TO DO
	}
}

// C::E models send coupling forces to the buffer
void 
MBDyn_CE_CEModel_SendToBuf(pMBDyn_CE_CEModel_t pMBDyn_CE_CEModel, std::vector<double> &MBDyn_CE_CouplingDynamic, 
							double* pMBDyn_CE_CEFrame, const unsigned& MBDyn_CE_NodesNum, const double* MBDyn_CE_CEScale)
{
	// obtain the transform matrix
	ChVector<> mbdynce_temp_frameMBDyn_pos(-pMBDyn_CE_CEFrame[0], -pMBDyn_CE_CEFrame[1], -pMBDyn_CE_CEFrame[2]);
	ChVector<> mbdynce_temp_frameMBDyn_rot_axis_X(pMBDyn_CE_CEFrame[3], pMBDyn_CE_CEFrame[4], pMBDyn_CE_CEFrame[5]);
	ChVector<> mbdynce_temp_frameMBDyn_rot_axis_Y(pMBDyn_CE_CEFrame[6], pMBDyn_CE_CEFrame[7], pMBDyn_CE_CEFrame[8]);
	ChVector<> mbdynce_temp_frameMBDyn_rot_axis_Z(pMBDyn_CE_CEFrame[9], pMBDyn_CE_CEFrame[10], pMBDyn_CE_CEFrame[11]);
	ChMatrix33<> mbdynce_temp_frameMBDyn_rot(mbdynce_temp_frameMBDyn_rot_axis_X,mbdynce_temp_frameMBDyn_rot_axis_Y,mbdynce_temp_frameMBDyn_rot_axis_Z);
	ChFrame<> mbdynce_temp_frameMBDyn(mbdynce_temp_frameMBDyn_pos, mbdynce_temp_frameMBDyn_rot);

	
	// write exchange force/torque to the buffer
	ChVector<> mbdynce_ce_force_G(0.0, 0.0, 0.0);
	ChVector<> mbdynce_ce_torque_G(0.0, 0.0, 0.0);
	ChVector<> mbdynce_mbdyn_force = mbdynce_temp_frameMBDyn.TransformDirectionParentToLocal(mbdynce_ce_force_G);
	ChVector<> mbdynce_mbdyn_torque = mbdynce_temp_frameMBDyn.TransformDirectionParentToLocal(mbdynce_ce_force_G);
	double* pmbdynce_tempvec3_f = &MBDyn_CE_CouplingDynamic[0];
	double* pmbdynce_tempvec3_m = &MBDyn_CE_CouplingDynamic[3*MBDyn_CE_NodesNum];
	for (unsigned i = 0; i < MBDyn_CE_NodesNum; i++)
	{
		// obtain the exchange force/torque from CE model; 
		// TO DO ...


		// write in the buffer
		double mbdynce_tempvec3_f[3]={mbdynce_mbdyn_force.x()/MBDyn_CE_CEScale[2],
										mbdynce_mbdyn_force.y()/MBDyn_CE_CEScale[2],
										mbdynce_mbdyn_force.z()/MBDyn_CE_CEScale[2]};
		double mbdynce_tempvec3_m[3]={mbdynce_mbdyn_torque.x()/(MBDyn_CE_CEScale[3]),
										mbdynce_mbdyn_torque.y()/(MBDyn_CE_CEScale[3]),
										mbdynce_mbdyn_torque.z()/(MBDyn_CE_CEScale[3])};
		memcpy(&pmbdynce_tempvec3_f[3*i], &mbdynce_tempvec3_f[0], 3 * sizeof(double));
		memcpy(&pmbdynce_tempvec3_m[3*i], &mbdynce_tempvec3_m[0], 3 * sizeof(double));
	}
}