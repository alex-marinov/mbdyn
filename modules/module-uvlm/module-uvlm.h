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


#include "mbdyn_uvlm.h"


/* We introduce a start class UVLM*/

class UvlmInterfaceBaseElem
: virtual public Elem, public UserDefinedElem {
private:
	DataManager *m_PDM;
	
	// Some other functions which are needed
	
	// A function that will print the coupling data on screen
	void MBDyn_UVLMPrint() const;
	
protected:
	// Coupling variables to be defined (see chrono interface code)
	
protected:
	double time_step;
	
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
	UvlmInterfaceBaseElem(unsigned uLabel, 
		const DofOwner *pDO,
		DataManager* pDM,           // Information for solvers (nodes, elements, solver, ....)
		MBDynParser& HP);           // Parse data from MBDyn input file.
		
//  Destructor
	virtual ~UvlmInterfaceBaseElem(void);
	
//  Functions that introduce member functions to handle the simulation
	void SetValue(DataManager *pDM, 
		VectorHandler& X, 
		VectorHandler& XP,
		SimulationEntity::Hints *ph);
	unsigned int iGetNumPrivData(void) const;
	
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
	
//  Miscellaneous member functions
	virtual void Output(OutputHandler& OH) const;
	std::ostream& Restart(std::ostream& out) const;
	virtual unsigned int iGetInitialNumDof(void) const;
	int iGetNumConnectedNodes(void) const;
	void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	
};

#endif