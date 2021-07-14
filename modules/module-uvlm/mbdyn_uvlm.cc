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

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include <cmath>
#include <unistd.h>
#include <map>
#include <vector>
#include <variant>
#include <algorithm>


#include "mbdyn_uvlm.h"


/***************************************************************/
// Functions private to this .cc file
void
write_2files();
/***************************************************************/


/***************************************************************/
// Functions used in coupling with MBDyn and UVLM subsystem
extern "C" void MBDyn_UVLM_Init()
{
	
	
	
	
	
	
}
/***************************************************************/





MBDyn_UVLM_AeroTimeStepInfo::MBDyn_UVLM_AeroTimeStepInfo
	(std::vector<std::pair<int, int>>& dimensions,
	std::vector<std::pair<int, int>>& dimensions_star){
		
	unsigned int number_of_surfaces = dimensions.size();

	zeta.resize(number_of_surfaces);
	for (int i = 0; i < number_of_surfaces; ++i)
	{
		zeta[i].resize(3);
		for (int j = 0; j < 3; ++j) {
			zeta[i][j].resize(dimensions[i].first + 1);
			for (int k = 0; k < dimensions[i].first + 1; ++k) {
				zeta[i][j][k].resize(dimensions[i].second + 1);
			}
		}
	}

	zeta_dot.resize(number_of_surfaces);
	for (int i = 0; i < number_of_surfaces; ++i)
	{
		zeta_dot[i].resize(3);
		for (int j = 0; j < 3; ++j) {
			zeta_dot[i][j].resize(dimensions[i].first + 1);
			for (int k = 0; k < dimensions[i].first + 1; ++k) {
				zeta_dot[i][j][k].resize(dimensions[i].second + 1);
			}
		}
	}
		
	normals.resize(number_of_surfaces);
	for (int i = 0; i < number_of_surfaces; ++i)
	{
		normals[i].resize(3);
		for (int j = 0; j < 3; ++j) {
			normals[i][j].resize(dimensions[i].first);
			for (int k = 0; k < dimensions[i].first; ++k) {
				normals[i][j][k].resize(dimensions[i].second);
			}
		}
	}

	forces.resize(number_of_surfaces);
	for (int i = 0; i < number_of_surfaces; ++i)
	{
		forces[i].resize(6);
		for (int j = 0; j < 6; ++j) {
			forces[i][j].resize(dimensions[i].first + 1);
			for (int k = 0; k < dimensions[i].first + 1; ++k) {
				forces[i][j][k].resize(dimensions[i].second + 1);
			}
		}
	}

	dynamic_forces.resize(number_of_surfaces);
	for (int i = 0; i < number_of_surfaces; ++i)
	{
		dynamic_forces[i].resize(6);
		for (int j = 0; j < 6; ++j) {
			dynamic_forces[i][j].resize(dimensions[i].first + 1);
			for (int k = 0; k < dimensions[i].first + 1; ++k) {
				dynamic_forces[i][j][k].resize(dimensions[i].second + 1);
			}
		}
	}

	zeta_star.resize(number_of_surfaces);
	for (int i = 0; i < number_of_surfaces; ++i)
	{
		zeta_star[i].resize(3);
		for (int j = 0; j < 3; ++j) {
			zeta_star[i][j].resize(dimensions_star[i].first + 1);
			for (int k = 0; k < dimensions_star[i].first + 1; ++k) {
				zeta_star[i][j][k].resize(dimensions_star[i].second + 1);
			}
		}
	}

	u_ext.resize(number_of_surfaces);
	for (int i = 0; i < number_of_surfaces; ++i)
	{
		u_ext[i].resize(3);
		for (int j = 0; j < 3; ++j) {
			u_ext[i][j].resize(dimensions[i].first + 1);
			for (int k = 0; k < dimensions[i].first + 1; ++k) {
				u_ext[i][j][k].resize(dimensions[i].second + 1);
			}
		}
	}

	u_ext_star.resize(number_of_surfaces);
	for (int i = 0; i < number_of_surfaces; ++i)
	{
		u_ext_star[i].resize(3);
		for (int j = 0; j < 3; ++j) {
			u_ext_star[i][j].resize(dimensions_star[i].first + 1);
			for (int k = 0; k < dimensions_star[i].first + 1; ++k) {
				u_ext_star[i][j][k].resize(dimensions_star[i].second + 1);
			}
		}
	}

	gamma.resize(number_of_surfaces);
	for (int i = 0; i < number_of_surfaces; ++i)
	{
		gamma[i].resize(dimensions[i].first);
		for (int j = 0; j < dimensions[i].first; ++j) {
			gamma[i][j].resize(dimensions[i].second);
		}
	}

	gamma_star.resize(number_of_surfaces);
	for (int i = 0; i < number_of_surfaces; ++i)
	{
		gamma_star[i].resize(dimensions_star[i].first);
		for (int j = 0; j < dimensions_star[i].first; ++j) {
			gamma_star[i][j].resize(dimensions_star[i].second);
		}
	}

	gamma_dot.resize(number_of_surfaces);
	for (int i = 0; i < number_of_surfaces; ++i)
	{
		gamma_dot[i].resize(dimensions[i].first);
		for (int j = 0; j < dimensions[i].first; ++j) {
			gamma_dot[i][j].resize(dimensions[i].second);
		}
	}

	dist_to_orig.resize(number_of_surfaces);
	for (int i = 0; i < number_of_surfaces; ++i)
	{
		dist_to_orig[i].resize(dimensions[i].first + 1);
		for (int j = 0; j < dimensions[i].first + 1; ++j) {
			dist_to_orig[i][j].resize(dimensions[i].second + 1);
		}
	}

}

MBDyn_UVLM_AeroTimeStepInfo::~MBDyn_UVLM_AeroTimeStepInfo(){
	
	NO_OP;
}


Aero_inputs::Aero_inputs(int num_elem, int num_node_elem, int num_surfaces,
	int num_node, int num_control_surfaces)
	:_num_elem(num_elem), _num_node_elem(num_node_elem), _num_surfaces(num_surfaces),
	_num_node(num_node), _num_control_surfaces(num_control_surfaces)
{
	chords.resize(_num_elem);
	twist.resize(_num_elem);
	sweep.resize(_num_elem);
	airfoil_distribution_input.resize(_num_elem);
	elastic_axis.resize(_num_elem);
	control_surface.resize(_num_elem);
	for (int i = 0; i < _num_elem; ++i) {
		chords[i].resize(_num_node_elem);
		twist[i].resize(_num_node_elem);
		sweep[i].resize(_num_node_elem);
		airfoil_distribution_input[i].resize(_num_node_elem);
		elastic_axis[i].resize(_num_node_elem);
		control_surface[i].resize(_num_node_elem);
	}

	surface_distribution_input.resize(_num_elem);
	surface_m.resize(_num_surfaces);
	aero_node_input.resize(_num_node);
	control_surface_type.resize(_num_control_surfaces);
	control_surface_chord.resize(_num_control_surfaces);
	control_surface_hinge_coord.resize(_num_control_surfaces);
}


void Aero_inputs::Aero_inputs_generate(int num_elem, int num_node_elem, int num_surfaces,
	int num_node, int num_control_surfaces) {

	// this function will generate all the necessary arrays and matrices for the input of the aero side of the solver;

}





void Aerogrid::generate(Aero_inputs* aero_inputs, Aerogrid_settings* aerogrid_settings,
	double ts) {


	// number of total nodes (structural + aero&struct)
	n_node = aero_inputs->aero_node_input.size();

	// number of elements
	n_elem = aero_inputs->surface_distribution_input.size();

	// Surface distribution
	surface_distribution = aero_inputs->surface_distribution_input;

	// number of surfaces
	n_surf = aero_inputs->surface_m.size();

	// number of chordwise panels
	surface_m = aero_inputs->surface_m;

	// number of aero nodes
	n_aero_node = aero_inputs->aero_node_input.size();

	// get N per surface
	calculate_dimensions(aero_inputs, aerogrid_settings);

	// Write grid info on the screen
	output_info(aero_inputs);

	// allocating the initial storage
	aero_ini_info = MBDyn_UVLM_AeroTimeStepInfo(dimensions, dimensions_star);

	
	// load airfoils db
	/*****************************/


	// pending

	/****************************/
	add_timestep();
	generate_mapping(aero_inputs);
	generate_zeta(aero_inputs, aerogrid_settings, ts);
}

void Aerogrid::calculate_dimensions(Aero_inputs* aero_inputs, Aerogrid_settings* aerogrid_settings) {

	dimensions.resize(n_surf);
	dimensions_star.resize(n_surf);

	for (int i = 0; i < n_surf; ++i) {
		// adding M values
		dimensions[i].first = surface_m[i];
	}
	// Count N values (actually the count result will be N+1)
	std::vector<std::vector<int>> nodes_in_surface;
	for (int i_elem = 0; i_elem < /*beam.num_elem*/; ++i_elem) {
		int nodes = ; // beam.elements[i_elem].global_connectivities
		int i_surf = aero_inputs->surface_distribution_input[i_elem];
		if (i_surf < 0) {
			continue;
		}
		for (int i_global_node = 0; i_global_node < nodes; ++i_global_node) {
			if (std::find(nodes_in_surface[i_surf].begin(), 
				nodes_in_surface[i_surf].end(), i_global_node) 
				!= nodes_in_surface[i_surf].end()) {
				continue;
			}
			else {
				nodes_in_surface[i_surf].push_back(i_global_node);
			}
			if (aero_inputs->aero_node_input[i_global_node]) {
				dimensions[i_surf].second += 1;
			}
		}
	}
	// Accounting for N+1 nodes --> N panels
	for (auto& itr : dimensions) {
		itr.second = itr.second - 1;
	}

	dimensions_star = dimensions;
	for (int i_surf = 0; i_surf < n_surf; ++i_surf) {
		dimensions_star[i_surf].first = aerogrid_settings->mstar;
	}

}

void Aerogrid::output_info(Aero_inputs* aero_inputs) {

	int _n_surf = aero_inputs->surface_m.size(); // should i use this defination of _n_surf or should I use the updated value of the private variable of this class(DOUBT###############)
	std::cout << "The aerodynamic grid contains " << _n_surf << " surfaces" << std::endl;

	for (int i_surf = 0; i_surf < _n_surf; ++i_surf) {
		std::cout << "Surface : " << i_surf << " ; " << "M : " << 
			dimensions[i_surf].first << " ; " << "N : " << 
			dimensions[i_surf].second << std::endl;
		std::cout << "Wake : " << i_surf << " ; " << "M : " <<
			dimensions_star[i_surf].first << " ; " << "N : " <<
			dimensions_star[i_surf].second << std::endl;
	}

	int total_bound_panels = 0;
	int total_wake_panels = 0;
	for (auto itr : dimensions) {
		total_bound_panels = total_bound_panels + itr.first * itr.second;
	}
	for (auto itr : dimensions_star) {
		total_wake_panels = total_wake_panels + itr.first * itr.second;
	}
	std::cout << "In total : " << total_bound_panels << " bound panels" << std::endl;
	std::cout << "In total : " << total_wake_panels << " wake panels" << std::endl;
	std::cout << "Total number of panels : " << total_bound_panels + total_wake_panels << std::endl;
}

void Aerogrid::add_timestep() {

	if (aero_timestep_info.empty()) {
		aero_timestep_info.push_back(aero_ini_info);
	}
	else {
		aero_timestep_info.push_back(aero_timestep_info.back());    // doubtful**************************************
	}
}

void Aerogrid::generate_zeta_timestep_info(Aero_inputs* aero_inputs, MBDyn_UVLM_AeroTimeStepInfo& aero_tstep, 
	Aerogrid_settings* aerogrid_settings) {

	int it = ; // beam.timestep_info - 1
	std::vector<std::vector<int>> global_node_in_surface;
	global_node_in_surface.resize(n_surf);

	// One surface per element
	for (int i_elem = 0; i_elem < n_elem; ++i_elem) {
		int i_surf = aero_inputs->surface_distribution_input[i_elem];
		// Check if we have to generate a surface here
		if (i_surf == -1) {
			continue;
		}

		for (int i_local_node = 0; i_local_node < /* beam.elements[i_elem].global_connectivities*/; ++i_local_node) {
			int i_global_node = ; // beam.elements[i_elem].global_connectivities[i_local_node]
			if (aero_inputs->aero_node_input[i_global_node] <= 0) {
				continue;
			}
			if (std::find(global_node_in_surface[i_surf].begin(),
				global_node_in_surface[i_surf].end(), i_global_node)
				!= global_node_in_surface[i_surf].end()) {
				continue;
			}
			else {
				global_node_in_surface[i_surf].push_back(i_global_node);
			}
			// find the i_surf and i_n data from the mapping
			int i_n = -1;
			int ii_surf = -1;
			for (int i = 0; i < struct2aero_mapping[i_global_node].size(); ++i) {
				i_n = struct2aero_mapping[i_global_node][i].second;
				ii_surf = struct2aero_mapping[i_global_node][i].first;
				if (ii_surf == i_surf) {
					break;
				}
			}
			if (i_n == -1 || ii_surf == -1) {
				std::cout << "Something failed with the mapping in aerogrid" << std::endl;
			}


			/// ////// Here define the node_info and fill it. Alternatively it is better to define the node info as a private data member
		}
	}
}

void Aerogrid::generate_zeta(Aero_inputs* aero_inputs, Aerogrid_settings* aerogrid_settings,
	double ts) {
	generate_zeta_timestep_info(aero_inputs, aero_timestep_info[ts], aerogrid_settings);
}

void Aerogrid::generate_mapping(Aero_inputs* aero_inputs) {

	struct2aero_mapping.resize(n_node);
	std::vector<int> surf_n_counter(n_surf, 0);
	std::vector<std::vector<int>> nodes_in_surface;
	nodes_in_surface.resize(n_surf);

	for (int i_elem = 0; i_elem < n_elem; ++i_elem) {
		int i_surf = aero_inputs->surface_distribution_input[i_elem];
		if (i_surf == -1) {
			continue;
		}
		for (int i_global_node = 0; i_global_node </*beam.elements[i_elem].reordered_global_connectivities*/; ++i_global_node) {
			if (aero_inputs->aero_node_input[i_global_node] <= 0) {
				continue;
			}
			if (std::find(nodes_in_surface[i_surf].begin(),
				nodes_in_surface[i_surf].end(), i_global_node)
				!= nodes_in_surface[i_surf].end()) {
				continue;
			}
			else {
				nodes_in_surface[i_surf].push_back(i_global_node);
				surf_n_counter[i_surf] += 1;
			}
			int i_n = surf_n_counter[i_surf] - 1;
			struct2aero_mapping[i_global_node].push_back(std::make_pair(i_surf, i_n));
		}
	}

	nodes_in_surface.clear();

	nodes_in_surface.resize(n_surf);
	aero2struct_mapping.resize(n_surf);
	for (int i_surf = 0; i_surf < n_surf; ++i_surf) {
		aero2struct_mapping[i_surf].resize(surf_n_counter[i_surf]);
	}
	for (int i_elem = 0; i_elem < n_elem; ++i_elem) {
		for (int i_global_node = 0; i_global_node </*beam.elements[i_elem].global_connectivities*/; ++i_global_node) {
			for (int i = 0; i < struct2aero_mapping[i_global_node].size(); ++i) {
				int i_surf = struct2aero_mapping[i_global_node][i].first;
				int i_n = struct2aero_mapping[i_global_node][i].second;
				if (std::find(nodes_in_surface[i_surf].begin(),
					nodes_in_surface[i_surf].end(), i_global_node)
					!= nodes_in_surface[i_surf].end()) {
					continue;
				}
				else {
					nodes_in_surface[i_surf].push_back(i_global_node);
				}
				aero2struct_mapping[i_surf][i_n] = i_global_node;
			}
		}
	}
}

void Aerogrid::compute_gamma_dot(Aero_inputs* aero_inputs, double dt, MBDyn_UVLM_AeroTimeStepInfo* tstep, 
	std::vector<MBDyn_UVLM_AeroTimeStepInfo*>& previous_tsteps) {

	int n_surf = aero_inputs->surface_m.size();

	if (previous_tsteps.size() == 0) {
		for (int i_surf = 0; i_surf < n_surf; ++i_surf) {
			for (int j = 0; j < tstep->gamma_dot[i_surf].size(); ++j) {
				for (int k = 0; k < tstep->gamma_dot[i_surf][j].size(); ++k) {
					tstep->gamma_dot[i_surf][j][k] = 0.0;
				}
			}
		}
	}
	else {
		for (int i_surf = 0; i_surf < n_surf; ++i_surf) {
			for (int j = 0; j < tstep->gamma_dot[i_surf].size(); ++j) {
				for (int k = 0; k < tstep->gamma_dot[i_surf][j].size(); ++k) {
					tstep->gamma_dot[i_surf][j][k] = (tstep->gamma[i_surf][j][k] - previous_tsteps[/*-2 in python sense*/]->gamma[i_surf][j][k]) / dt;
				}
			}
		}
	}
}

void Aerogrid::generate_strip() {

}



StraightWake::StraightWake() {
	NO_OP;
}

void StraightWake::StraightWake_initialize(StraighWake_settings* straightwakesettings) {

	
	_u_inf = straightwakesettings->u_inf;
	_u_inf_direction = straightwakesettings->u_inf_direction;
	_dt = straightwakesettings->dt;

	if (straightwakesettings->dx1 == -1) {
		_dx1 = _u_inf * _dt;
	}
	else {
		_dx1 = straightwakesettings->dx1;
	}

	_ndx1 = straightwakesettings->ndx1;
	_r = straightwakesettings->r;

	if (straightwakesettings->dxmax == -1) {
		_dxmax = _dx1;
	}
	else {
		_dxmax = straightwakesettings->dxmax;
	}
}

void StraightWake::StraightWake_generate(MBDyn_UVLM_AeroTimeStepInfo* aero_tstep) {

	int nsurf = aero_tstep->zeta.size();
	for (int isurf = 0; isurf < nsurf; ++isurf) {
		int M = aero_tstep->zeta_star[isurf][0].size();
		int N = aero_tstep->zeta_star[isurf][0][0].size();
		for (int j = 0; j < N; ++j) {
			for (int k = 0; k < 3; ++k) {
				aero_tstep->zeta_star[isurf][k][0][j] = aero_tstep->zeta[isurf][k][M][j];
			}
			for (int i = 1; i < M; ++i) {
				double deltax = StraightWake_get_deltax(i, _dx1, _ndx1, _r, _dxmax);
				for (int k = 0; k < 3; ++k) {
					aero_tstep->zeta_star[isurf][k][i][j] = aero_tstep->zeta_star[isurf][k][i - 1][j] + deltax * _u_inf_direction[k];
				}
			}
		}
		for (int j = 0; j < aero_tstep->gamma[isurf].size(); ++j) {
			for (int i = 0; i < aero_tstep->gamma[isurf][0].size(); ++i) {
				aero_tstep->gamma[isurf][j][i] *= 0;
				aero_tstep->gamma_star[isurf][j][i] *= 0;
			}
		}
	}

	for (int isurf = 0; isurf < nsurf; ++isurf) {
		int M = aero_tstep->zeta_star[isurf][0].size();
		int N = aero_tstep->zeta_star[isurf][0][0].size();
		for (int j = 0; j < aero_tstep->dist_to_orig[isurf].size(); ++j) {
			for (int i = 0; i < aero_tstep->dist_to_orig[isurf][0].size(); ++i) {
				aero_tstep->dist_to_orig[isurf][j][i] = 0;
			}
		}
		for (int j = 0; j < N; ++j) {
			for (int i = 1; i < M; ++i) {
				double norm = 0;
				for (int k = 0; k < 3; ++k) {
					norm = norm + (aero_tstep->zeta_star[isurf][k][i][j] - aero_tstep->zeta_star[isurf][k][i - 1][j]) *
						(aero_tstep->zeta_star[isurf][k][i][j] - aero_tstep->zeta_star[isurf][k][i - 1][j]);
				}
				norm = sqrt(norm);
				aero_tstep->dist_to_orig[isurf][i][j] = aero_tstep->dist_to_orig[isurf][i - 1][j] + norm;
			}
			for (int i = 0; i < M; ++i) {
				aero_tstep->dist_to_orig[isurf][i][j] /= aero_tstep->dist_to_orig[isurf][M][j];
			}
		}
	}
}

StraightWake::~StraightWake() {
	NO_OP;
}

double StraightWake::StraightWake_get_deltax(int i, double dx1, double ndx1, double r, double dxmax) {

	double deltax;
	if (i < ndx1 + 1) {
		deltax = dx1;
	}
	else {
		deltax = dx1 * pow(r, i - ndx1);
	}
	deltax = std::min(deltax, dxmax);

	return deltax;
}


