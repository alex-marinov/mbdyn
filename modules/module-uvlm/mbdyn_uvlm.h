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



#ifndef MBDYN_UVLM_H
#define MBDYN_UVLM_H

#include <vector>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include <cmath>


#include "dataman.h"
#include "aeroelem.h"
#include "UVLM-master/src/uvlm.h"

extern "C" {

	//- IDs of coupling bodies UVLM
	//- a set of IDs
	//- Body_Label: body label in UVLM;
	//- bool output: output the body motion in UVLM;
	//- constraint (MotionImposed) in UVLM:
	//- Node in MBDyn  >>>>> Body Ground (Marker Ground);
	//- Body in UVLM >>>>> Body i (Marker i);
	//- constraint connects Marker i and Marker Ground;
	//- constraints are described in Marker i local ref. 

	struct MBDYN_UVLM_MODELDATA {
		unsigned MBDyn_UVLM_Body_Label;
		bool bMBDyn_UVLM_Body_Output;
		bool bMBDyn_UVLM_Constraint[6];              //- position and rotation constraints, expresed in node local ref.
		double MBDyn_UVLM_Body_Offset[3];            //- relative position between constraints (Marker i) and UVLM body, expressed in UVLM body ref.
		double MBDyn_UVLM_Body_Rh[9];                //- relative rotation between constraints and MBDyn Node, expressed in node ref.
	};

	struct StepUVLM_settings {
		int gamma_dot_filtering;
		std::string velocity_field_generator;
		std::string gust_shape;
		double gust_length;
		double gust_intensity;
		double gust_offset;
		double span;
		bool relative_motion;
		unsigned int n_time_steps;
	};

	struct Aerogrid_settings {
		bool unsteady;
		bool aligned_grid;
		std::vector<double> freestream_dir;
		int mstar;
		std::string wake_shape_generator;  // this should contain a string which will decide which struct to use out of starightwake and helicoidal wake
	};

	struct StraightWake_settings {
		double u_inf;
		std::vector<double> u_inf_direction;
		double dt;
		double dx1;
		int ndx1;
		double r;
		double dxmax;
	};

	struct UVMopts {
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

	struct FlightConditions {
		double uinf;
		std::vector<double> uinf_direction;
		double rho;
		double c_ref;
	};


	class AeroTimeStepInfo
	{
		typedef std::vector<double> vector_1d;
		typedef std::vector<vector_1d> vector_2d;
		typedef std::vector<vector_2d> vector_3d;
		typedef std::vector<vector_3d> vector_4d;

	public:
		vector_4d zeta;
		vector_4d zeta_dot;
		vector_4d normals;
		vector_4d forces;
		vector_4d dynamic_forces;
		vector_4d zeta_star;
		vector_4d u_ext;
		vector_4d u_ext_star;
		vector_3d gamma;
		vector_3d gamma_star;
		vector_3d gamma_dot;
		vector_3d dist_to_orig;

		AeroTimeStepInfo();
		void initialize(std::vector<std::pair<unsigned int, unsigned int>>&,
			std::vector<std::pair<unsigned int, unsigned int>>&);
		~AeroTimeStepInfo();
	};


	class Beam_inputs {

	public:
		int num_node_elem;
		int num_elem;
		int num_node;
		std::vector<std::vector<int>> connectivities;
		std::vector<std::vector<int>> reordered_connectivities; 
		std::vector<std::vector<std::vector<double>>> frame_of_reference_delta;

		Beam_inputs();
		void Beam_inputs_generate(int);
		~Beam_inputs();
	};


	class Aero_inputs {

	public:
		std::map<int, std::vector<std::pair<double, double>>> airfoils; // This is an airfoil group which contains airfoil data for each and every part corresponding to map id
		std::vector<std::vector<double>> chords;    // array with the chords of every airfoil given in an element/node basis
		std::vector<std::vector<double>> twist;     // Has the twist angle in radians. It is implemented as rotation around the loxal x axis
		std::vector<std::vector<double>> sweep;     // Has the twist angle in radians. It is implemented as rotation around the local z axis
		std::vector<std::vector<int>> airfoil_distribution_input; // Airfoil distribution. Contains the indices of the airfoils that are present in the airfoils group.
		std::vector<int> surface_distribution_input; // It contains the index of the surface that element belongs to. Surfaces need to be continuous.
		std::vector<int> surface_m;   // number of chordwise pannels for each surface
		std::string m_distribution;   // is a string with the chordwie panel distribution generally "uniform"
		std::vector<bool> aero_node_input;    // Aerodyanmic node definition. if indicates if that node has a lifting surface attached to it
		std::vector<std::vector<double>> elastic_axis;  // Indicates the elastic axis location with respect to the leading edge as a fraction of the chord of that rib. Note that the elastic axis is already determined, as the beam is fixed now, so this setting controls the location of the lifting surface with respect to the beam

		Aero_inputs();
		void Aero_inputs_generate(int);
		~Aero_inputs();

	private:
		double naca(double, double, double);
		std::vector<std::pair<double, double>>& generate_naca_camber(double, double);
	};


	class Aerogrid {

	public:

		//- Aerodynamic data at each time step
		AeroTimeStepInfo aero_ini_info;
		std::vector<AeroTimeStepInfo> aero_timestep_info;
		std::vector<int> surface_distribution;
		std::vector<int> surface_m;
		std::vector<std::pair<unsigned int, unsigned int>> dimensions;
		std::vector<std::pair<unsigned int, unsigned int>> dimensions_star;
		std::map<int, std::vector<double>> airfoil_db;
		std::vector<std::vector<std::pair<int, int>>> struct2aero_mapping;
		std::vector<std::vector<int>> aero2struct_mapping;

		//- Kinematic data required at each time step
		std::vector<std::vector<double>> node_displacements;
		std::vector<std::vector<double>> node_displacements_der;
		std::vector<std::vector<double>> node_CRV;
		std::vector<std::vector<double>> node_CRV_der;
		std::vector<std::vector<std::vector<double>>> node_cga;


		Aerogrid();
		void generate(Aero_inputs&, Beam_inputs&, const Aerogrid_settings*, double);
		void output_info(Aero_inputs&);
		void calculate_dimensions(Aero_inputs&, Beam_inputs&, const Aerogrid_settings*);
		void add_timestep();
		void generate_zeta_timestep_info(Aero_inputs&, Beam_inputs&, AeroTimeStepInfo&, const Aerogrid_settings*);
		void generate_zeta(Aero_inputs&, Beam_inputs&, const Aerogrid_settings*, double);
		void generate_mapping(Aero_inputs&, Beam_inputs&);
		void compute_gamma_dot(Aero_inputs&, double, double);
		void generate_strip(Aero_inputs&, const Aerogrid_settings*, AeroTimeStepInfo&);
		~Aerogrid();

	private:
		int n_node = 0;
		int n_elem = 0;
		int n_surf = 0;
		int n_aero_node = 0;
		int n_control_surfaces = 0;

		struct node_info_dict {
			int i_node;
			int i_local_node;
			int i_surf;
			int i_n;
			double chord;
			double eaxis;
			double twist;
			double sweep;
			int M;
			std::string M_distribution;
			int airfoil;
			std::vector<double> beam_coord;
			std::vector<double> pos_dot;
			std::vector<double> beam_psi;
			std::vector<double> psi_dot;
			std::vector<double> for_delta;
			std::vector<std::vector<double>> cga;
		} node_info;
	};


	class UvlmLibVar {

	public:

		unsigned int** p_dimensions = NULL;
		unsigned int** p_dimensions_star = NULL;
		unsigned int i_iter;
		double** p_uext = NULL;
		double** p_uext_star = NULL;
		double** p_zeta = NULL;
		double** p_zeta_star = NULL;
		double** p_zeta_dot = NULL;
		double** p_zeta_star_dot = NULL;
		double* p_rbm_vel = NULL;
		double* p_centre_rot = NULL;
		double** p_gamma = NULL;
		double** p_gamma_star = NULL;
		double** p_gamma_dot = NULL;
		double** p_dist_to_orig = NULL;
		double** p_normals = NULL;
		double** p_forces = NULL;
		double** p_dynamic_forces = NULL;

		UvlmLibVar();
		void UvlmLibVar_generate(Aerogrid&, double);
		void UvlmLibVar_save(Aerogrid&, double);
		~UvlmLibVar();
	};


	class StraightWake {

	public:
		StraightWake();
		void StraightWake_initialize(const StraightWake_settings*);
		void StraightWake_generate(AeroTimeStepInfo&);
		~StraightWake();

	private:
		double _u_inf;                                  // Freestream velocity magnitude
		std::vector<double> _u_inf_direction;           // x, y, z relative components of the freestream velocity
		double _dt;                                     // Time step
		double _dx1;                                    // Size of the first wake panel
		int _ndx1;                                      // Number of panels with size "dx1"
		double _r;                                      // Growth rate after "ndx1" panels
		double _dxmax;                                  // Maximum panel size

		double StraightWake_get_deltax(int, double, double, double, double);
	};


	class SteadyVelocityField {

	public:
		SteadyVelocityField();
		void SteadyVelocityField_initialize(const FlightConditions*);
		void SteadyVelocityField_generate(AeroTimeStepInfo&, const UVMopts*);
		~SteadyVelocityField();

	private:

		double u_inf;
		std::vector<double> u_inf_direction;

	};

	// Coupling cases: 
	enum MBDyn_UVLM_COUPLING
	{
		COUPLING_NONE = -2,
		COUPLING_STSTAGGERED = -1,   // Not implemented right now
		COUPLING_LOOSE = 0,          // Loose coupling
		COUPLING_TIGHT = 1           // Tight coupling
	};

	// Loose coupling:
	enum MBDyn_UVLM_COUPLING_LOOSE
	{
		TIGHT = 0,                   // meaningless
		LOOSE_EMBEDDED = 1,          // Embedded function method
		LOOSE_JACOBIAN = 2,          // Jacobian/parallel type
		LOOSE_GAUSS = 3              // Gauss/serial type
	};

	// Type of force:
	enum MBDyn_UVLM_FORCETYPE
	{
		REACTION_FORCE = 0,  //- using reaction forces as coupling forces (coupled by constraints)
		CONTACT_FORCE = 1,   //- using contact forces as coupling forces (coupled by anytype, but using coupling forces acting on the coupling bodies)
	};

	enum MBDYN_UVLM_OUTPUTTYPE
	{
		MBDYN_UVLM_OUTPUT_ALLBODIES = 0,              //- outputing kinematic data of all bodies in UVLM
		MBDYN_UVLM_OUTPUT_SELECTEDCOUPLINGBODIES = 1, //- outputing kinematic data of selected coupling bodies in UVLM
	};

	//- create the Init function
	void
		MBDyn_UVLM_Model_Init(const StepUVLM_settings* MBDyn_UVLM_StepUVLM_settings,
			const Aerogrid_settings* MBDyn_UVLM_Aerogrid_settings,
			const StraightWake_settings* MBDyn_UVLM_StraightWake_settings,
			const UVMopts* MBDyn_UVLM_UVMopts,
			const FlightConditions* MBDyn_UVLM_FlightConditions,
			Beam_inputs& MBDyn_UVLM_Beam_inputs,
			Aero_inputs& MBDyn_UVLM_Aero_inputs,
			StraightWake& MBDyn_UVLM_StraightWake,
			SteadyVelocityField& MBDyn_UVLM_SteadyVelocityField,
			UvlmLibVar& MBDyn_UVLM_UvlmLibVar,
			UVLM::Types::VMopts& VMoptions,
			UVLM::Types::UVMopts& UVMoptions,
			UVLM::Types::FlightConditions& FlightConditions,
			unsigned MBDyn_UVLM_NodesNum);

	//- Performs a step of Integration
	void
		MBDyn_UVLM_Model_DoStepDynamics(const StepUVLM_settings* MBDyn_UVLM_StepUVLM_settings,
			const Aerogrid_settings* MBDyn_UVLM_Aerogrid_settings,
			const StraightWake_settings* MBDyn_UVLM_StraightWake_settings,
			const UVMopts* MBDyn_UVLM_UVMopts,
			const FlightConditions* MBDyn_UVLM_FlightConditions,
			Beam_inputs& MBDyn_UVLM_Beam_inputs,
			Aero_inputs& MBDyn_UVLM_Aero_inputs,
			StraightWake& MBDyn_UVLM_StraightWake,
			SteadyVelocityField& MBDyn_UVLM_SteadyVelocityField,
			Aerogrid& MBDyn_UVLM_Aerogrid,
			UvlmLibVar& MBDyn_UVLM_UvlmLibVar,
			UVLM::Types::UVMopts& UVMoptions,
			UVLM::Types::FlightConditions& FlightConditions,
			double time_step);

	//- UVLM models receive coupling motion from the buffer
	void
		MBDyn_UVLM_Model_RecvFromBuf(Aerogrid& MBDyn_UVLM_Aerogrid, 
			std::vector<double> &MBDyn_CE_CouplingKinematic,
			const unsigned &MBDyn_CE_NodesNum);

	//- UVLM models send the force information to the buffer 
	void
		MBDyn_UVLM_Model_SendToBuf(Aerogrid& MBDyn_UVLM_Aerogrid,
			std::vector<double>& MBDyn_UVLM_CouplingDynamic,
			const unsigned &MBDyn_UVLM_NodesNum,
			double time_step);

}
#endif //- MBDYN_UVLM_H
