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



extern "C" {

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
	std::string control_surface_generator;
	// input array or matrix for the control surface generator needs to be defined
	std::string wake_shape_generator;  // this should contain a string which will decide which struct to use out of starightwake and helicoidal wake
};

struct StraighWake_settings {
	double u_inf;
	std::vector<double> u_inf_direction;
	double dt;
	double dx1;
	int ndx1;
	double r;
	double dxmax;
};

struct MBDyn_UVLM_UVMopts {
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

struct MBDyn_UVLM_FlightConditions {
	double uinf;
	std::vector<double> uinf_direction;
	double rho;
	double c_ref;
};
	

class MBDyn_UVLM_AeroTimeStepInfo
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

	MBDyn_UVLM_AeroTimeStepInfo(std::vector<std::pair<int, int>>&,
		std::vector<std::pair<int, int>>&);
	~MBDyn_UVLM_AeroTimeStepInfo();

};



class Aero_inputs {

public:
	// have to add the airfoil also
	std::vector<std::vector<double>> chords;    // array with the chords of every airfoil given in an element/node basis
	std::vector<std::vector<double>> twist;     // Has the twist angle in radians. It is implemented as rotation around the loxal x axis
	std::vector<std::vector<double>> sweep;     // Has the twist angle in radians. It is implemented as rotation around the local z axis
	std::vector<std::vector<int>> airfoil_distribution_input; // Airfoil distribution. Contains the indices of the airfoils that are present in the airfoils group.
	std::vector<int> surface_distribution_input; // It contains the index of the surface that element belongs to. Surfaces need to be continuous.
	std::vector<int> surface_m;   // number of chordwise pannels for each surface
	std::string m_distribution;   // is a string with the chordwie panel distribution generally "uniform"
	std::vector<bool> aero_node_input;    // Aerodyanmic node definition. if indicates if that node has a lifting surface attached to it
	std::vector<std::vector<double>> elastic_axis;  // Indicates the elastic axis location wit respect to the leading edge as a fraction of the chord of that rib. Note that the elastic axis is already determined, as the beam is fixed now, so this setting controls the location of the lifting surface with respect to the beam
	std::vector<std::vector<int>> control_surface;  // Integer array containing -1 if that section has no control surface associated to it, and 0,1,2... if the section belongs to the conrol surface 0, 1, 2....
	std::vector<int> control_surface_type;         // contains 0 if the control surface deflection is static and 1 if it is dynamic
	std::vector<int> control_surface_chord;  // It is an integer array with the number of panels belonging to the control surface. For example , if M=4 and we want our control surface to ne 0.25c then we need to put 1
	std::vector<int> control_surface_hinge_coord;  // only necessary for lifting surfaces that are deflected as a whole, like some horizontal tails in some aircrafts. we let it be 0 if we are not modelling this

	Aero_inputs(int, int, int, int, int);
	void Aero_inputs_generate(int, int, int, int, int);

	// throw in some ostream functions too for some fancy stuff
	~Aero_inputs() {};

private:
	int _num_elem;      // Total number of elements
	int _num_node_elem;  // number of nodes per element
	int _num_surfaces;   // numer of surfaces
	int _num_node;    // total number of nodes
	int _num_control_surfaces;   // number of control surfaces
};



class Aerogrid {

public:

	MBDyn_UVLM_AeroTimeStepInfo aero_ini_info;
	std::vector<MBDyn_UVLM_AeroTimeStepInfo> aero_timestep_info;
	std::vector<int> surface_distribution;
	std::vector<int> surface_m;
	std::vector<std::pair<int, int>> dimensions;
	std::vector<std::pair<int, int>> dimensions_star;
	std::vector<std::vector<std::pair<int, int>>> struct2aero_mapping;
	std::vector<std::vector<int>> aero2struct_mapping;

	Aerogrid();
	void generate(Aero_inputs*, Aerogrid_settings*, double);
	void output_info(Aero_inputs*);
	void calculate_dimensions(Aero_inputs*, Aerogrid_settings*);
	void add_timestep();
	void generate_zeta_timestep_info(Aero_inputs*, MBDyn_UVLM_AeroTimeStepInfo&, Aerogrid_settings*);
	void generate_zeta(Aero_inputs*, Aerogrid_settings*, double);
	void generate_mapping(Aero_inputs*);
	void compute_gamma_dot(Aero_inputs*, double, MBDyn_UVLM_AeroTimeStepInfo*, std::vector<MBDyn_UVLM_AeroTimeStepInfo*>&);
	void generate_strip();


	~Aerogrid();

private:
	int n_node = 0;
	int n_elem = 0;
	int n_surf = 0;
	int n_aero_node = 0;
	int n_control_surfaces = 0;
};


class StraightWake {

public:
	StraightWake();
	void StraightWake_initialize(StraighWake_settings*);
	void StraightWake_generate(MBDyn_UVLM_AeroTimeStepInfo*);
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




































}
#endif //- MBDYN_UVLM_H
