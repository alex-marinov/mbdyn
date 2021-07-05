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
 
#ifndef MODULE_UVLM_H
#define MODULE_UVLM_H

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <iostream>
#include <cfloat>
#include <vector>
#include <map>
#include <string>

#include "dataman.h"
#include "userelem.h"
#include "elem.h"
#include "strnode.h"
#include "dataman.h"
#include "converged.h"

#include "mbdyn_uvlm.h"


/* We introduce a start class UVLM*/

class UvlmInterfaceBaseElem
: virtual public Elem, public UserDefinedElem {
private:
	DataManager *m_PDM;
	// A function that transfer doublereal Vec3 to double Vec3 (data between MBDyn structural part and the UVLM aero part should be of the same type)
	void MBDyn_UVLM_Vec3D(const Vec3& mbdynuvlm_Vec3, double *mbdynuvlm_temp, double MBDyn_UVLM_LengthScale) const;
	// A function that transfer doublereal Mat3x3 to double Mat3x3 (data between MBDyn structural part and the UVLM aero part should be of the same type)
	void MBDyn_UVLM_Mat3x3D(const Mat3x3& mbdynuvlm_Mat3x3, double *mbdynuvlm_temp) const;
	double MBDyn_UVLM_calculateError();  // Calculate the error of coupling forces
	// A function that prints the coupling data on the screen
	void MBDyn_UVLM_MBDynPrint() const;
	
protected:
	std::vector<double> MBDyn_UVLM_CouplingKinematic;                       //- for coupling motion
	std::vector<double> MBDyn_UVLM_CouplingDynamic;                         //- for coupling forces
	std::vector<double> MBDyn_UVLM_CouplingDynamic_pre;                     //- for coupling forces in last iterations.
	double *pMBDyn_UVLM_CouplingKinematic_x = NULL;                         //- consistent with the external struc force element
	double *pMBDyn_UVLM_CouplingKinematic_R = NULL;
	double *pMBDyn_UVLM_CouplingKinematic_xp = NULL;
	double *pMBDyn_UVLM_CouplingKinematic_omega = NULL;
	double *pMBDyn_UVLM_CouplingKinematic_xpp = NULL;
	double *pMBDyn_UVLM_CouplingKinematic_omegap = NULL;
	double *pMBDyn_UVLM_CouplingDynamic_f = NULL;
	double *pMBDyn_UVLM_CouplingDynamic_m = NULL;
	double *pMBDyn_UVLM_CouplingDynamic_f_pre = NULL;
	double *pMBDyn_UVLM_CouplingDynamic_m_pre = NULL;
	struct {
		unsigned Size_Kinematic;
		unsigned Size_Dynamic;
	} MBDyn_UVLM_CouplingSize;
	//- some parameters about the convergence
	unsigned MBDyn_UVLM_CouplingIter_Max;
	unsigned MBDyn_UVLM_CouplingIter_Count;
	double MBDyn_UVLM_Coupling_Tol;

protected:
	Converged MBDyn_UVLM_Model_Converged;                            //- denote whether the coupling variables are converged
	bool bMBDyn_UVLM_Model_DoStepDynamics;                           //- detect whether UVLM model is needed to be simulated and sends back data
	bool bMBDyn_UVLM_FirstSend;                                      //- whether the current residual is the first or not..
	bool bMBDyn_UVLM_Verbose;                                        //- whether UVLM codes print the solution process at each iteration.
	int MBDyn_UVLM_OutputType;                                       //- type of outputs
	mutable std::ofstream MBDyn_UVLM_out_forces;                     //- ofstream for outputing coupling forces.

public:
	double MBDyn_UVLM_Scale[4];

	//- Coupling nodes information
	struct MBDYN_UVLM_POINTDATA {
		unsigned MBDyn_UVLM_uLabel;
//		unsigned MBDyn_CE_CEBody_Label;                              //- coupling bodies in C::E model
		const StructNode *pMBDyn_UVLM_Node;
		Vec3 MBDyn_UVLM_Offset;                                      //- offset of the marker in MBDyn. By default, MBDyn_CE_Offset == null;
		Mat3x3 MBDyn_UVLM_RhM;                                       //- orientation of the marker in MBDyn. By default, Rh_M == eye; bool constraints are also defined in this orientation,
		Vec3 MBDyn_UVLM_F;
		Vec3 MBDyn_UVLM_M;
	};
	
protected:
	double time_step;
	std::vector<MBDYN_UVLM_POINTDATA> MBDyn_UVLM_Nodes;              //- Nodes info in MBDyn
	unsigned MBDyn_UVLM_NodesNum;

public:
	int MBDyn_UVLM_CouplingType;
	int MBDyn_UVLM_CouplingType_loose;
	int MBDyn_UVLM_ForceType;
	
public:
	std::map<std::string, std::variant<std::string, int, double>> UVLM_STEPUVLM_settings;
	std::map<std::string, std::variant<std::string, int, double>> UVLM_STATICUVLM_settings;
	struct MBDyn_UVLM_VMopts{
		bool ImageMethod;
		bool Steady;
		bool horseshoe;
		bool KJMeth;
		bool NewAIC;
		double DelTime;
		bool Rollup;
		unsigned int NumCores;
		unsigned int NumSurfaces;
		double dt;
		unsigned int n_rollup;
		double rollup_tolerance;
		unsigned int rollup_aic_refresh;
		bool iterative_solver;
		double iterative_tol;
		bool iterative_precond;
		bool cfl1;
		double vortex_radius;
		double vortex_radius_wake_ind;
		std::vector<double> rbm_vel_g;
	};
	struct MBDyn_UVLM_UVMopts{
		double dt;
		unsigned int NumCores;
		unsigned int NumSurfaces;
		unsigned int convection_scheme;
		bool ImageMethod;
		bool iterative_solver;
		double iterative_tol;
		bool iterative_precond;
		bool convect_wake;
		bool cfl1;
		double vortex_radius;
		double vortex_radius_wake_ind;
		unsigned int interp_coords;
		unsigned int filter_method;
		unsigned int interp_method;
		double yaw_slerp;
		bool quasi_steady;
	};
	struct MBDyn_UVLM_FlightConditions{
		double uinf;
		std::vector<double> uinf_direction;
		double rho;
		double c_ref;
	};
	
	// decalra the variables for raw inputs
	
	
	
	
	

public:

    int MBDyn_UVLM_CouplingType;
    int MBDyn_UVLM_CouplingType_loose;
    int MBDyn_UVLM_ForceType;


//  Constructor for the UvlmInterfaceBaseElem 
	UvlmInterfaceBaseElem(unsigned uLabel,  // Label
		const DofOwner *pDO,
		DataManager* pDM,           // Information for solvers (nodes, elements, solver, ....)
		MBDynParser& HP);           // Parse data from MBDyn input file.
		
//  Destructor
	virtual ~UvlmInterfaceBaseElem(void);
	
//  Functions that introduce member functions to handle the simulation
	virtual void SetValue(DataManager *pDM, 
		VectorHandler& X, 
		VectorHandler& XP,
		SimulationEntity::Hints *ph = 0);
	virtual void Update(const VectorHandler &XCurr,
		const VectorHandler &XprimeCurr);
	virtual void AfterConvergence(const VectorHandler &X,
		const VectorHandler &XP);
	virtual void AfterPredict(VectorHandler &X,
		VectorHandler &XP);
	virtual void BeforePredict(VectorHandler &X, VectorHandler &XP, VectorHandler &XPrev, VectorHandler & XPPrev) const;
//	unsigned int iGetNumPrivData(void) const;
	
//  Functions for the element, which set Jac and Res
//  Intial Assembly
	virtual void InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   	VariableSubMatrixHandler& InitialAssJac(VariableSubMatrixHandler& WorkMat, 
		const VectorHandler& XCurr);
   	SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec, 
		const VectorHandler& XCurr);
//  Assembly
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
	VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef, 
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
	SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr);
	
//	Functions for coupling variables
	void MBDyn_UVLM_UpdateUVLMModel();                    //- update a regular step in C::E
	void MBDyn_UVLM_SendDataToBuf_Curr();                 //- write current kinematic variables to the vector.
	void MBDyn_UVLM_SendDataToBuf_Prev();                 //- write previous kinematic variables to the vector.
	void MBDyn_UVLM_KinematicData_Interpolate();          //- intepolates the kinematic for peer (to do)
	void MBDyn_UVLM_RecvDataFromBuf();                    //- read the dynamic variables from the Buf.

//  Miscellaneous member functions
	virtual void Output(OutputHandler& OH) const;
	std::ostream& Restart(std::ostream& out) const;
	virtual unsigned int iGetInitialNumDof(void) const;
	int iGetNumConnectedNodes(void) const;
	void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	
};

#endif