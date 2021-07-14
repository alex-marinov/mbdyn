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

std::vector<double>
QuadInterpol1D(std::vector<double>& x, std::vector<double>& y, std::vector<double>& xp) {

	std::vector<double> yp(xp.size(), 0);
	// This function implements Lagrange interpolation method
	for (int k = 0; k < xp.size(); ++k) {
		for (int i = 0; i < y.size(); ++i) {
			int p = 1;
			for (int j = 0; j < x.size(); ++j) {
				if (i != j) {
					p = p * (xp[k] - x[j]) / (x[i] - x[j]);
				}
			}
			yp[k] = yp[k] + p * y[i];
		}
	}

	return yp;
}

void
rotation3d_x(std::vector<std::vector<double>>& Mat, double angle) {

	double c = cos(angle);
	double s = sin(angle);

	Mat[0][0] = 1.0;
	Mat[0][1] = 0.0;
	Mat[0][2] = 0.0;
	Mat[1][0] = 0.0;
	Mat[1][1] = c;
	Mat[1][2] = -s;
	Mat[2][0] = 0.0;
	Mat[2][1] = s;
	Mat[2][2] = c;
}

void
rotation3d_z(std::vector<std::vector<double>>& Mat, double angle) {

	double c = cos(angle);
	double s = sin(angle);

	Mat[0][0] = c;
	Mat[0][1] = -s;
	Mat[0][2] = 0.0;
	Mat[1][0] = s;
	Mat[1][1] = c;
	Mat[1][2] = 0.0;
	Mat[2][0] = 0.0;
	Mat[2][1] = 0.0;
	Mat[2][2] = 1.0;
}

void
skew(std::vector<std::vector<double>>& Mat, std::vector<double>& vec) {

	if (vec.size() != 3) {
		std::cout << "The input vector is not 3D" << std::endl;
	}
	else {
		Mat[1][2] = -vec[0];
		Mat[2][0] = -vec[1];
		Mat[0][1] = -vec[2];
		Mat[2][1] = vec[0];
		Mat[0][2] = vec[1];
		Mat[1][0] = vec[2];
	}
}

void
MatMatMul(std::vector<std::vector<double>>& Mat_1, std::vector<std::vector<double>>& Mat_2,
	std::vector<std::vector<double>>& Mat_3) {

	for (int i = 0; i < Mat_1.size(); i++) {
		for (int j = 0; j < Mat_1.size(); j++) {
			Mat_3[i][j] = 0;
			for (int k = 0; k < Mat_1.size(); k++) {
				Mat_3[i][j] += Mat_1[i][k] * Mat_2[k][j];
			}
		}
	}
}

std::vector<double>
MatVecMul(std::vector<std::vector<double>>& Mat, std::vector<double>& Vec) {

	std::vector<double> Res_Vec(3, 0);

	Res_Vec[0] = Mat[0][0] * Vec[0] + Mat[0][1] * Vec[1] + Mat[0][2] * Vec[2];
	Res_Vec[1] = Mat[1][0] * Vec[0] + Mat[1][1] * Vec[1] + Mat[1][2] * Vec[2];
	Res_Vec[2] = Mat[2][0] * Vec[0] + Mat[2][1] * Vec[1] + Mat[2][2] * Vec[2];

	return Res_Vec;
}

void
crv2rotation(std::vector<std::vector<double>>& Mat, std::vector<double>& psi) {

	double norm_psi = 0;
	for (const auto& itr : psi) {
		norm_psi += itr * itr;
	}
	norm_psi = sqrt(norm_psi);

	std::vector<std::vector<double>> eye(3, std::vector<double>(3, 0));
	eye[0][0] = 1.0;
	eye[1][1] = 1.0;
	eye[2][2] = 1.0;

	std::vector<std::vector<double>> Matmul(3, std::vector<double>(3, 0));

	if (norm_psi < 1.0e-15) {
		std::vector<std::vector<double>> skew_psi(3, std::vector<double>(3, 0));
		skew(skew_psi, psi);
		MatMatMul(skew_psi, skew_psi, Matmul);
		for (int i = 0; i < Mat.size(); ++i) {
			for (int j = 0; j < Mat[0].size(); ++j) {
				Mat[i][j] = eye[i][j] + skew_psi[i][j] + 0.5 * Matmul[i][j];
			}
		}
	}
	else {
		std::vector<double> normal;
		normal = psi;
		std::transform(normal.begin(), normal.end(), normal.begin(), [norm_psi](auto& c) {return c / norm_psi; });
		std::vector<std::vector<double>> skew_normal(3, std::vector<double>(3, 0));
		skew(skew_normal, normal);

		MatMatMul(skew_normal, skew_normal, Matmul);
		for (int i = 0; i < Mat.size(); ++i) {
			for (int j = 0; j < Mat[0].size(); ++j) {
				Mat[i][j] = eye[i][j] + sin(norm_psi) * skew_normal[i][j] + (1.0 - cos(norm_psi)) * Matmul[i][j];
			}
		}
	}
}

void
crv2tan(std::vector<std::vector<double>>& Mat, std::vector<double>& psi) {

	double norm_psi = 0;
	for (const auto& itr : psi) {
		norm_psi += itr * itr;
	}
	norm_psi = sqrt(norm_psi);
	std::vector<std::vector<double>> psi_skew(3, std::vector<double>(3, 0));
	skew(psi_skew, psi);

	std::vector<std::vector<double>> eye(3, std::vector<double>(3, 0));
	eye[0][0] = 1.0;
	eye[1][1] = 1.0;
	eye[2][2] = 1.0;
	std::vector<std::vector<double>> Matmul(3, std::vector<double>(3, 0));
	MatMatMul(psi_skew, psi_skew, Matmul);
	double eps = 1.0e-8;
	if (norm_psi < eps) {
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				Mat[i][j] = eye[i][j] - 0.5*psi_skew[i][j] + (1.0 / 6.0)*Matmul[i][j];
			}
		}
	}
	else {
		double k1 = (cos(norm_psi) - 1.0) / (norm_psi*norm_psi);
		double k2 = (1.0 - sin(norm_psi) / norm_psi) / (norm_psi*norm_psi);
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				Mat[i][j] = eye[i][j] + k1 * psi_skew[i][j] + k2 * Matmul[i][j];
			}
		}
	}
}

void
crv_dot2omega(std::vector<double>& crv, std::vector<double>& crv_dot, std::vector<double>& omega_a) {

	std::vector<std::vector<double>> Mat(3, std::vector<double>(3, 0));
	crv2tan(Mat, crv);
	std::vector<std::vector<double>> Mat_T(3, std::vector<double>(3, 0));
	Mat_T[0][0] = Mat[0][0];
	Mat_T[0][1] = Mat[1][0];
	Mat_T[0][2] = Mat[2][0];
	Mat_T[1][0] = Mat[0][1];
	Mat_T[1][1] = Mat[1][1];
	Mat_T[1][2] = Mat[2][1];
	Mat_T[2][0] = Mat[0][2];
	Mat_T[2][1] = Mat[1][2];
	Mat_T[2][2] = Mat[2][2];
	omega_a = MatVecMul(Mat_T, crv_dot);
}

double
angle_between_vectors_sign(std::vector<double>& vec_a, std::vector<double>& vec_b, std::vector<double>& orientation_in) {

	std::vector<double> cross_prod;
	cross_prod[0] = vec_a[1] * vec_b[2] - vec_a[2] * vec_b[1];
	cross_prod[1] = vec_a[2] * vec_b[0] - vec_a[0] * vec_b[2];
	cross_prod[2] = vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0];

	double dot_prod = 0;
	for (int i = 0; i < vec_a.size(); i++) {
		dot_prod += vec_a[i] * vec_b[i];
	}

	double norm = 0;
	for (const auto& itr : cross_prod) {
		norm += itr * itr;
	}
	norm = sqrt(norm);

	double angle = atan2(norm, dot_prod);

	dot_prod = 0;
	for (int i = 0; i < orientation_in.size(); ++i) {
		dot_prod += orientation_in[i] * cross_prod[i];
	}

	if (dot_prod < 0) {
		angle *= -1;
	}

	return angle;
}


/***************************************************************/


/***************************************************************/
// Functions used in coupling with MBDyn and UVLM subsystem
extern "C" void MBDyn_UVLM_Init()
{






}
/***************************************************************/





MBDyn_UVLM_AeroTimeStepInfo::MBDyn_UVLM_AeroTimeStepInfo(std::vector<std::pair<int, int>>& dimensions,
	std::vector<std::pair<int, int>>& dimensions_star) {

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

MBDyn_UVLM_AeroTimeStepInfo::~MBDyn_UVLM_AeroTimeStepInfo() {
	NO_OP;
}



Beam_inputs::Beam_inputs() {
	NO_OP;
}

void Beam_inputs::Beam_inputs_generate() {

	// do remember to input the raw inputs from MBDyn

	connectivities.resize(num_elem);
	reordered_connectivities.resize(num_elem);
	frame_of_reference_delta.resize(num_elem);
	for (int i = 0; i < num_elem; ++i) {
		connectivities[i].resize(num_node_elem);
		reordered_connectivities[i].resize(num_node_elem);
		frame_of_reference_delta[i].resize(num_node_elem);
		for (int j = 0; j < 3; ++j) {
			frame_of_reference_delta[i][j].resize(3);
		}
	}

	// right wing (beam 0)
	int working_elem = 0;
	for (int ielem = 0; ielem < num_elem / 2; ++ielem) {
		for (int inode = 0; inode < num_node_elem; ++inode) {
			frame_of_reference_delta[ielem][inode][0] = -1;
			frame_of_reference_delta[ielem][inode][1] = 0;
			frame_of_reference_delta[ielem][inode][2] = 0;
		}
	}
	// connectivity
	for (int ielem = 0; ielem < num_elem / 2; ++ielem) {
		connectivities[ielem][0] = 0 + ielem * (num_node_elem - 1);
		connectivities[ielem][1] = 2 + ielem * (num_node_elem - 1);
		connectivities[ielem][2] = 1 + ielem * (num_node_elem - 1);
		working_elem++;
	}

	// left wing (beam 1)
	for (int ielem = working_elem; ielem < num_elem; ++ielem) {
		for (int inode = 0; inode < num_node_elem; ++inode) {
			frame_of_reference_delta[ielem][inode][0] = -1;
			frame_of_reference_delta[ielem][inode][1] = 0;
			frame_of_reference_delta[ielem][inode][2] = 0;
		}
	}
	// connectivity
	for (int ielem = working_elem; ielem < num_elem; ++ielem) {
		connectivities[ielem][0] = 0 + ielem * (num_node_elem - 1) + 1;
		connectivities[ielem][1] = 2 + ielem * (num_node_elem - 1) + 1;
		connectivities[ielem][2] = 1 + ielem * (num_node_elem - 1) + 1;
		working_elem++;
	}

	reordered_connectivities = connectivities;
}

Beam_inputs::~Beam_inputs() {
	NO_OP;
}


Aero_inputs::Aero_inputs() {
	NO_OP;
}

void
Aero_inputs::Aero_inputs_generate() {

	// do remeber to input the values of the private variables of this class
	// along with the necessary raw inputs (see generate_planarwing in SHARPy)

	surface_distribution_input.resize(_num_elem);
	surface_m.resize(_num_surfaces);
	m_distribution = "uniform";
	aero_node_input.resize(_num_node);

	airfoil_distribution_input.resize(_num_elem);
	twist.resize(_num_elem);
	chords.resize(_num_elem);
	elastic_axis.resize(_num_elem);
	sweep.resize(_num_elem);
	for (int i = 0; i < _num_elem; ++i) {
		chords[i].resize(_num_node_elem);
		twist[i].resize(_num_node_elem);
		sweep[i].resize(_num_node_elem);
		airfoil_distribution_input[i].resize(_num_node_elem);
		elastic_axis[i].resize(_num_node_elem);
	}

	int working_elem = 0;
	int working_node = 0;
	// right wing (surface 0 , beam 0)
	int i_surf = 0;
	for (int i = 0; i < _num_elem / 2; ++i) {
		for (int j = 0; j < _num_node_elem; ++j) {
			airfoil_distribution_input[i][j] = 0;
			chords[i][j] = ;  // input from the MBDyn input
			elastic_axis[i][j] = ; // input from the MBDyn input
			twist[i][j] = ; // input from MBDyn
			sweep[i][j] = ; // input from MBDyn
		}
		surface_distribution_input[i] = i_surf;
		working_elem++;
	}
	surface_m[i_surf] = ; // number of chordwise panelling for each surface
	for (int i = 0; i < _num_node / 2; ++i) {
		aero_node_input[i] = true;
		working_node++;
	}

	// left wing (surface 1, beam 1)
	i_surf = 1;
	for (int i = working_elem; i < _num_elem; ++i) {
		for (int j = 0; j < _num_node_elem; ++j) {
			airfoil_distribution_input[i][j] = 0;
			chords[i][j] = ;  // input from the MBDyn input
			elastic_axis[i][j] = ; // input from the MBDyn input
			twist[i][j] = ; // input from MBDyn
			sweep[i][j] = ; // input from MBDyn
		}
		surface_distribution_input[i] = i_surf;
		working_elem++;
	}
	surface_m[i_surf] = ; // number of chordwise panelling for each surface
	for (int i = working_node; i < _num_node; ++i) {
		aero_node_input[i] = true;
		working_node++;
	}

	// Generate data for a single airfoil. 
	// In order to create data for more airfoils we just have to change the vaue of M and P.
	double P = 0;
	double M = 0;
	airfoils[0] = generate_naca_camber(P, M);

}

double
Aero_inputs::naca(double x, double m, double p) {

	if (x < 1.0e-6) {
		return 0.0;
	}
	else if (x < p) {
		return m / (p*p)*(2 * p*x - x * x);
	}
	else if (x > p && x < 1.0e-6) {
		return m / ((1 - p)*(1 - p))*(1 - 2 * p + 2 * p*x - x * x);
	}
}

std::vector<std::pair<double, double>>&
Aero_inputs::generate_naca_camber(double P, double M) {

	double m = M * 1.0e-2;
	double p = P * 1.0e-1;

	std::vector<std::pair<double, double>> airfoil_data(1000, std::make_pair(0, 0));
	for (int i = 0; i < airfoil_data.size(); ++i) {
		airfoil_data[i].first = i / (1000 - 1);
		airfoil_data[i].second = naca(airfoil_data[i].first, m, p);
	}

	return airfoil_data;
}

Aero_inputs::~Aero_inputs() {
	NO_OP;
}





void Aerogrid::generate(Aero_inputs* aero_inputs, Beam_inputs* beam_inputs, Aerogrid_settings* aerogrid_settings,
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
	calculate_dimensions(aero_inputs, beam_inputs, aerogrid_settings);

	// Write grid info on the screen
	output_info(aero_inputs);

	// allocating the initial storage
	aero_ini_info = MBDyn_UVLM_AeroTimeStepInfo(dimensions, dimensions_star);

	add_timestep();
	generate_mapping(aero_inputs, beam_inputs);
	generate_zeta(aero_inputs, beam_inputs, aerogrid_settings, ts);
}

void Aerogrid::calculate_dimensions(Aero_inputs* aero_inputs, Beam_inputs* beam_inputs,
	Aerogrid_settings* aerogrid_settings) {

	dimensions.resize(n_surf);
	dimensions_star.resize(n_surf);

	for (int i = 0; i < n_surf; ++i) {
		// adding M values
		dimensions[i].first = surface_m[i];
	}
	// Count N values (actually the count result will be N+1)
	std::vector<std::vector<int>> nodes_in_surface;
	for (int i_elem = 0; i_elem < beam_inputs->num_elem; ++i_elem) {
		std::vector<int> nodes;
		nodes = beam_inputs->connectivities[i_elem];
		int i_surf = aero_inputs->surface_distribution_input[i_elem];
		if (i_surf < 0) {
			continue;
		}
		for (int i_global_node = 0; i_global_node < nodes.size(); ++i_global_node) {
			if (std::find(nodes_in_surface[i_surf].begin(),
				nodes_in_surface[i_surf].end(), nodes[i_global_node])
				!= nodes_in_surface[i_surf].end()) {
				continue;
			}
			else {
				nodes_in_surface[i_surf].push_back(nodes[i_global_node]);
			}
			if (aero_inputs->aero_node_input[nodes[i_global_node]]) {
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

void Aerogrid::generate_zeta_timestep_info(Aero_inputs* aero_inputs, Beam_inputs* beam_inputs, MBDyn_UVLM_AeroTimeStepInfo& aero_tstep,
	Aerogrid_settings* aerogrid_settings) {

	std::vector<std::vector<int>> global_node_in_surface;
	global_node_in_surface.resize(n_surf);

	// One surface per element
	for (int i_elem = 0; i_elem < n_elem; ++i_elem) {
		int i_surf = aero_inputs->surface_distribution_input[i_elem];
		// Check if we have to generate a surface here
		if (i_surf == -1) {
			continue;
		}

		for (int i_local_node = 0; i_local_node < beam_inputs->connectivities[i_elem].size(); ++i_local_node) {
			int i_global_node = beam_inputs->connectivities[i_elem][i_local_node];
			if (aero_inputs->aero_node_input[i_global_node] == false) {
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


			node_info.i_node = i_global_node;
			node_info.i_local_node = i_local_node;
			node_info.i_surf = i_surf;
			node_info.i_n = i_n;
			node_info.chord = aero_inputs->chords[i_elem][i_local_node];
			node_info.eaxis = aero_inputs->elastic_axis[i_elem][i_local_node];
			node_info.twist = aero_inputs->twist[i_elem][i_local_node];
			node_info.sweep = aero_inputs->sweep[i_elem][i_local_node];
			node_info.M = dimensions[i_surf].first;
			node_info.M_distribution = aero_inputs->m_distribution;
			node_info.airfoil = aero_inputs->airfoil_distribution_input[i_elem][i_local_node];
			//node_info.beam_coord = ;
			//node_info.pos_dot = ;
			//node_info.beam_psi = ;
			//node_info.psi_dot = ;
			node_info.for_delta = beam_inputs->frame_of_reference_delta[i_elem][i_local_node];
			//node_info.cga = ;

			generate_strip(aero_inputs, aerogrid_settings, aero_tstep);
		}
	}
}

void Aerogrid::generate_zeta(Aero_inputs* aero_inputs, Beam_inputs* beam_inputs, Aerogrid_settings* aerogrid_settings,
	double ts) {
	generate_zeta_timestep_info(aero_inputs, beam_inputs, aero_timestep_info[ts], aerogrid_settings);
}

void Aerogrid::generate_mapping(Aero_inputs* aero_inputs, Beam_inputs* beam_inputs) {

	struct2aero_mapping.resize(n_node);
	std::vector<int> surf_n_counter(n_surf, 0);
	std::vector<std::vector<int>> nodes_in_surface;
	nodes_in_surface.resize(n_surf);

	std::vector<int> global_nodes;
	for (int i_elem = 0; i_elem < n_elem; ++i_elem) {
		int i_surf = aero_inputs->surface_distribution_input[i_elem];
		if (i_surf == -1) {
			continue;
		}
		global_nodes = beam_inputs->reordered_connectivities[i_elem];
		for (int i_global_node = 0; i_global_node < global_nodes.size(); ++i_global_node) {
			if (aero_inputs->aero_node_input[global_nodes[i_global_node]] == false) {
				continue;
			}
			if (std::find(nodes_in_surface[i_surf].begin(),
				nodes_in_surface[i_surf].end(), global_nodes[i_global_node])
				!= nodes_in_surface[i_surf].end()) {
				continue;
			}
			else {
				nodes_in_surface[i_surf].push_back(global_nodes[i_global_node]);
				surf_n_counter[i_surf] += 1;
			}
			int i_n = surf_n_counter[i_surf] - 1;
			struct2aero_mapping[global_nodes[i_global_node]].push_back(std::make_pair(i_surf, i_n));
		}
	}

	nodes_in_surface.clear();

	nodes_in_surface.resize(n_surf);
	aero2struct_mapping.resize(n_surf);
	for (int i_surf = 0; i_surf < n_surf; ++i_surf) {
		aero2struct_mapping[i_surf].resize(surf_n_counter[i_surf]);
	}
	for (int i_elem = 0; i_elem < n_elem; ++i_elem) {
		global_nodes = beam_inputs->connectivities[i_elem];
		for (int i_global_node = 0; i_global_node < global_nodes.size(); ++i_global_node) {
			for (int i = 0; i < struct2aero_mapping[global_nodes[i_global_node]].size(); ++i) {
				int i_surf = struct2aero_mapping[global_nodes[i_global_node]][i].first;
				int i_n = struct2aero_mapping[global_nodes[i_global_node]][i].second;

				if (std::find(nodes_in_surface[i_surf].begin(),
					nodes_in_surface[i_surf].end(), global_nodes[i_global_node])
					!= nodes_in_surface[i_surf].end()) {
					continue;
				}
				else {
					nodes_in_surface[i_surf].push_back(global_nodes[i_global_node]);
				}
				aero2struct_mapping[i_surf][i_n] = global_nodes[i_global_node];
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

void Aerogrid::generate_strip(Aero_inputs* aero_inputs, Aerogrid_settings* aerogrid_settings, MBDyn_UVLM_AeroTimeStepInfo& aero_tstep) {

	std::vector<std::vector<double>> strip_coordinates_a_frame(3, std::vector<double>(node_info.M + 1));
	std::vector<std::vector<double>> strip_coordinates_b_frame(3, std::vector<double>(node_info.M + 1));
	std::vector<std::vector<double>> zeta_dot_a_frame(3, std::vector<double>(node_info.M + 1));

	// airfoil coordinates
	// we are going to store everything in the x-z plane of the b
	// FoR, so that the transformation Cab rotates everything in place.
	if (node_info.M_distribution == "uniform") {
		for (int i = 0; i < strip_coordinates_b_frame[0].size(); ++i) {
			strip_coordinates_b_frame[1][i] = i / (strip_coordinates_b_frame[0].size() - 1);
		}
	}
	else if (node_info.M_distribution == "1-cos") {
		for (int i = 0; i < strip_coordinates_b_frame[0].size(); ++i) {
			double arg = i / (strip_coordinates_b_frame[0].size() - 1);
			strip_coordinates_b_frame[1][i] = 0.5*(1 - cos(arg * M_PI));
		}
	}
	else {
		std::cout << "M_distribution is not implemented" << std::endl;
	}

	// load airfoils db
	std::vector<std::pair<double, double>> airfoils_coords;
	airfoils_coords = aero_inputs->airfoils[node_info.airfoil];
	std::vector<double> airfoils_coords_x_c(airfoils_coords.size());
	std::vector<double> airfoils_coords_y_c(airfoils_coords.size());
	for (int i = 0; i < airfoils_coords.size(); ++i) {
		airfoils_coords_x_c[i] = airfoils_coords[i].first;
		airfoils_coords_y_c[i] = airfoils_coords[i].second;
	}
	strip_coordinates_b_frame[2] = QuadInterpol1D(airfoils_coords_x_c, airfoils_coords_y_c, strip_coordinates_b_frame[1]);


	// Elastic axis correction
	for (int i_M = 0; i_M < node_info.M + 1; ++i_M) {
		strip_coordinates_b_frame[1][i_M] -= node_info.eaxis;
	}

	// Chord scaling
	double cons = node_info.chord;
	for (int row = 0; row < strip_coordinates_b_frame.size(); ++row) {
		std::transform(strip_coordinates_b_frame[row].begin(), strip_coordinates_b_frame[row].end(),
			strip_coordinates_b_frame[row].begin(), [cons](auto& c) {return cons * c; });
	}

	// Twist transformation (rotation around x_b axis)
	std::vector<std::vector<double>> Ctwist(3, std::vector<double>(3, 0));
	if (abs(node_info.twist) > 1.0e-6) {
		rotation3d_x(Ctwist, node_info.twist);
	}
	else {
		Ctwist[0][0] = 1.0;
		Ctwist[1][1] = 1.0;
		Ctwist[2][2] = 1.0;
	}

	// Cab transformation
	std::vector<std::vector<double>> Cab(3, std::vector<double>(3, 0));
	crv2rotation(Cab, node_info.beam_psi);

	std::vector<double> vec_a(Cab.size(), 0);
	std::vector<double> vec_b(Cab.size(), 0);
	for (int i = 0; i < Cab.size(); ++i) {
		vec_a[i] = Cab[i][1];
		vec_b[i] = Cab[i][2];
	}
	double rot_angle = angle_between_vectors_sign(vec_a, vec_b, aerogrid_settings->freestream_dir);

	double dot_prod = 0;
	for (int i = 0; i < Cab.size(); ++i) {
		dot_prod += aerogrid_settings->freestream_dir[i] * Cab[i][1];
	}
	if (dot_prod >= 0) {
		rot_angle += 0;
	}
	else {
		rot_angle += -2.0*M_PI;
	}
	std::vector<std::vector<double>> Crot(3, std::vector<double>(3, 0));
	rotation3d_z(Crot, -rot_angle);

	std::vector<std::vector<double>> c_sweep(3, std::vector<double>(3, 0));
	if (abs(node_info.sweep) > 1.0e-6) {
		rotation3d_z(c_sweep, node_info.sweep);
	}
	else {
		c_sweep[0][0] = 1.0;
		c_sweep[1][1] = 1.0;
		c_sweep[2][2] = 1.0;
	}

	// Transformation from beam to beam prime (with sweep and twist)
	std::vector<double> MatVecMul_Res(3, 0);
	std::vector<double> strip_coordinates_b_frame_vec(strip_coordinates_b_frame.size(), 0);
	for (int i_M = 0; i_M < node_info.M + 1; ++i_M) {
		for (int i = 0; i < strip_coordinates_b_frame.size(); ++i) {
			strip_coordinates_b_frame_vec[i] = strip_coordinates_b_frame[i][i_M];
		}
		MatVecMul_Res = MatVecMul(Ctwist, strip_coordinates_b_frame_vec);
		MatVecMul_Res = MatVecMul(Crot, MatVecMul_Res);
		MatVecMul_Res = MatVecMul(c_sweep, MatVecMul_Res);
		strip_coordinates_b_frame[0][i_M] = MatVecMul_Res[0];
		strip_coordinates_b_frame[1][i_M] = MatVecMul_Res[1];
		strip_coordinates_b_frame[2][i_M] = MatVecMul_Res[2];


		MatVecMul_Res = MatVecMul(Cab, strip_coordinates_b_frame_vec);
		strip_coordinates_a_frame[0][i_M] = MatVecMul_Res[0];
		strip_coordinates_a_frame[1][i_M] = MatVecMul_Res[1];
		strip_coordinates_a_frame[2][i_M] = MatVecMul_Res[2];
	}

	// Zeta_dot
	// Velocity due to pos_dot
	for (int i_M = 0; i_M < node_info.M + 1; ++i_M) {
		zeta_dot_a_frame[0][i_M] += node_info.pos_dot[0];
		zeta_dot_a_frame[1][i_M] += node_info.pos_dot[1];
		zeta_dot_a_frame[2][i_M] += node_info.pos_dot[2];
	}
	// Velocity due to psi_dot
	std::vector<double> omega_a(3, 0);
	crv_dot2omega(node_info.beam_psi, node_info.psi_dot, omega_a);
	std::vector<std::vector<double>> skew_omega_a(3, std::vector<double>(3, 0));
	skew(skew_omega_a, omega_a);
	std::vector<double> strip_coordinates_a_frame_vec(strip_coordinates_a_frame.size(), 0);
	for (int i_M = 0; i_M < node_info.M + 1; ++i_M) {
		for (int i = 0; i < strip_coordinates_a_frame.size(); ++i) {
			strip_coordinates_a_frame_vec[i] = strip_coordinates_a_frame[i][i_M];
		}
		MatVecMul_Res = MatVecMul(skew_omega_a, strip_coordinates_a_frame_vec);

		zeta_dot_a_frame[0][i_M] += MatVecMul_Res[0];
		zeta_dot_a_frame[1][i_M] += MatVecMul_Res[1];
		zeta_dot_a_frame[2][i_M] += MatVecMul_Res[2];
	}

	// Add node coords
	for (int i_M = 0; i_M < node_info.M + 1; ++i_M) {
		strip_coordinates_a_frame[0][i_M] += node_info.beam_coord[0];
		strip_coordinates_a_frame[1][i_M] += node_info.beam_coord[1];
		strip_coordinates_a_frame[2][i_M] += node_info.beam_coord[2];
	}

	// Add quarter-chord displacement
	std::vector<double> delta_c(strip_coordinates_a_frame.size(), 0);
	delta_c[0] = (strip_coordinates_a_frame[0][strip_coordinates_a_frame[0].back()] - strip_coordinates_a_frame[0][0]) / node_info.M;
	delta_c[1] = (strip_coordinates_a_frame[1][strip_coordinates_a_frame[1].back()] - strip_coordinates_a_frame[1][0]) / node_info.M;
	delta_c[2] = (strip_coordinates_a_frame[2][strip_coordinates_a_frame[2].back()] - strip_coordinates_a_frame[2][0]) / node_info.M;
	if (node_info.M_distribution == "uniform") {
		for (int i_M = 0; i_M < node_info.M + 1; ++i_M) {
			strip_coordinates_a_frame[0][i_M] += 0.25*delta_c[0];
			strip_coordinates_a_frame[1][i_M] += 0.25*delta_c[1];
			strip_coordinates_a_frame[2][i_M] += 0.25*delta_c[2];
		}
	}
	else {
		std::cout << "No quarter chord disp of grid for non uniform grids" << std::endl;
	}

	// Rotation from a to g
	std::vector<double> zeta_dot_a_frame_vec(zeta_dot_a_frame.size(), 0);
	for (int i_M = 0; i_M < node_info.M + 1; ++i_M) {
		for (int i = 0; i < strip_coordinates_a_frame.size(); ++i) {
			strip_coordinates_a_frame_vec[i] = strip_coordinates_a_frame[i][i_M];
			zeta_dot_a_frame_vec[i] = zeta_dot_a_frame[i][i_M];
		}
		MatVecMul_Res = MatVecMul(node_info.cga, strip_coordinates_a_frame_vec);
		strip_coordinates_a_frame[0][i_M] = MatVecMul_Res[0];
		strip_coordinates_a_frame[1][i_M] = MatVecMul_Res[1];
		strip_coordinates_a_frame[2][i_M] = MatVecMul_Res[2];

		MatVecMul_Res = MatVecMul(node_info.cga, zeta_dot_a_frame_vec);
		zeta_dot_a_frame[0][i_M] = MatVecMul_Res[0];
		zeta_dot_a_frame[1][i_M] = MatVecMul_Res[1];
		zeta_dot_a_frame[2][i_M] = MatVecMul_Res[2];
	}


	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < node_info.M + 1; ++j) {
			aero_tstep.zeta[node_info.i_surf][i][j][node_info.i_n] = strip_coordinates_a_frame[i][j];
			aero_tstep.zeta_dot[node_info.i_surf][i][j][node_info.i_n] = zeta_dot_a_frame[i][j];
		}
	}

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


