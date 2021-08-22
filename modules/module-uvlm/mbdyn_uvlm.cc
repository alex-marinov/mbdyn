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
#include <sys/time.h>
#include <unistd.h>
#include <dlfcn.h>


#include "mbdyn_uvlm.h"


clock_t startTime, endTime;
struct timeval start_time, end_time;


 // Functions private to this .cc file [START]

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

std::vector<double>
matrix2skewvec(std::vector<std::vector<double>>& Mat) {

	std::vector<double> vec(3, 0);

	vec[0] = Mat[2][1] - Mat[1][2];
	vec[1] = Mat[0][2] - Mat[2][0];
	vec[2] = Mat[1][0] - Mat[0][1];

	return vec;
}

std::vector<double>
rotation2quat(std::vector<std::vector<double>>& Mat) {

	std::vector<std::vector<double>> s(4, std::vector<double>(4, 0));
	
	s[0][0] = 1.0 + Mat[0][0] + Mat[1][1] + Mat[2][2];
	std::vector<double> vec(3, 0);
	vec = matrix2skewvec(Mat);
	s[0][1] = vec[0];
	s[0][2] = vec[1];
	s[0][3] = vec[2];

	s[1][0] = Mat[2][1] - Mat[1][2];
	s[1][1] = 1.0 + Mat[0][0] - Mat[1][1] - Mat[2][2];
	s[1][2] = Mat[0][1] + Mat[1][0];
	s[1][3] = Mat[0][2] + Mat[2][0];

	s[2][0] = Mat[0][2] - Mat[2][0];
	s[2][1] = Mat[1][0] + Mat[0][1];
	s[2][2] = 1.0 - Mat[0][0] + Mat[1][1] - Mat[2][2];
	s[2][3] = Mat[1][2] + Mat[2][1];

	s[3][0] = Mat[1][0] - Mat[0][1];
	s[3][1] = Mat[0][2] + Mat[2][0];
	s[3][2] = Mat[1][2] + Mat[2][1];
	s[3][3] = 1.0 - Mat[0][0] - Mat[1][1] + Mat[2][2];

	std::vector<double> diag_s(4, 0);
	diag_s[0] = s[0][0];
	diag_s[1] = s[1][1];
	diag_s[2] = s[2][2];
	diag_s[3] = s[3][3];

	double smax = *std::max_element(diag_s.begin(), diag_s.end());
	int ismax = std::distance(diag_s.begin(), std::max_element(diag_s.begin(), diag_s.end()));
	// compute quaternion angles
	std::vector<double> quat(4, 0);
	quat[ismax] = 0.5 * sqrt(smax);
	for (int i = 0; i < 4; i++) {
		if (i == ismax) {
			continue;
		}
		quat[i] = 0.25 * s[ismax][i] / quat[ismax];
	}

	// quat bound
	if (quat[0] < 0) {
		quat[0] *= -1;
		quat[1] *= -1;
		quat[2] *= -1;
		quat[3] *= -1;
	}

	return quat;
}

std::vector<double>
quat2crv(std::vector<double>& quat) {

	double crv_norm = 2.0*acos(std::max(-1.0, std::min(quat[0], 1.0)));

	// normal vector
	std::vector<double> psi(3, 0);
	if (abs(crv_norm) < 1.0e-15) {
		psi = { 0.0, 0.0, 0.0 };
	}
	else {
		psi[0] = crv_norm * quat[1] / sin(crv_norm*0.5);
		psi[1] = crv_norm * quat[2] / sin(crv_norm*0.5);
		psi[2] = crv_norm * quat[3] / sin(crv_norm*0.5);
	}

	return psi;
}


std::vector<double>
rotation2crv(std::vector<std::vector<double>>& Mat) {

	std::vector<double> quat;
	std::vector<double> psi;

	quat = rotation2crv(Mat);
	psi = quat2crv(quat);

	// crv bounds
	double norm_ini = 0;
	for (const auto& itr : psi) {
		norm_ini += itr * itr;
	}
	norm_ini = sqrt(norm_ini);

	// forces the norm to be in [-pi, pi]
	double norm = norm_ini - 2.0*M_PI * int(norm_ini / (2 * M_PI));

	if (norm == 0.0) {
		psi[0] *= 0.0;
		psi[1] *= 0.0;
		psi[2] *= 0.0;
	}
	else {
		if (norm > M_PI) {
			norm -= 2.0*M_PI;
		}
		else if (norm < -M_PI) {
			norm += 2.0*M_PI;
		}
		psi[0] *= (norm / norm_ini);
		psi[1] *= (norm / norm_ini);
		psi[2] *= (norm / norm_ini);
	}

	return psi;
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
angle_between_vectors_sign(std::vector<double>& vec_a, std::vector<double>& vec_b, const std::vector<double>& orientation_in) {

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

std::vector<double>
flatten(std::vector<std::vector<double>>& Mat) {

	std::vector<double> ret;
	for (auto &v : Mat)
		ret.insert(ret.end(), v.begin(), v.end());

	return ret;
}


/*  [START] */
// Functions used in coupling with MBDyn and UVLM subsystem
extern "C" void 
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
	UVLM::Types::FlightConditions& Flight_Conditions,
	unsigned MBDyn_UVLM_NodesNum) {

	gettimeofday(&start_time, NULL);  // get the start time
	startTime = clock();

	// Initialize the StraightWake
	MBDyn_UVLM_StraightWake.StraightWake_initialize(MBDyn_UVLM_StraightWake_settings);

	// Initialize the Velocity field
	MBDyn_UVLM_SteadyVelocityField.SteadyVelocityField_initialize(MBDyn_UVLM_FlightConditions);

	// Create variables for the initialization function of the Uvlm library
	Aerogrid aerogrid;
	double time_step = 0;  // Initial iteration
	aerogrid.generate(MBDyn_UVLM_Aero_inputs, MBDyn_UVLM_Beam_inputs, MBDyn_UVLM_Aerogrid_settings, time_step);
	MBDyn_UVLM_StraightWake.StraightWake_generate(aerogrid.aero_timestep_info[0]);
	MBDyn_UVLM_SteadyVelocityField.SteadyVelocityField_generate(aerogrid.aero_timestep_info[0], MBDyn_UVLM_UVMopts);

	
	// constructing the VMopts as input
	VMoptions.cfl1 = MBDyn_UVLM_UVMopts->cfl1;
	VMoptions.DelTime = 1.0;
	VMoptions.dt = MBDyn_UVLM_UVMopts->dt;
	VMoptions.horseshoe = true;
	VMoptions.ImageMethod = MBDyn_UVLM_UVMopts->ImageMethod;
	VMoptions.iterative_precond = MBDyn_UVLM_UVMopts->iterative_precond;
	VMoptions.iterative_solver = MBDyn_UVLM_UVMopts->iterative_solver;
	VMoptions.iterative_tol = MBDyn_UVLM_UVMopts->iterative_tol;
	VMoptions.KJMeth = false;
	VMoptions.NewAIC = false;
	VMoptions.NumCores = MBDyn_UVLM_UVMopts->NumCores;
	VMoptions.NumSurfaces = MBDyn_UVLM_UVMopts->NumSurfaces;
	VMoptions.n_rollup = 0;
	VMoptions.Rollup = false;
	VMoptions.rollup_aic_refresh = 1;
	VMoptions.rollup_tolerance = 1.0e-5;
	VMoptions.Steady = true;
	VMoptions.vortex_radius = MBDyn_UVLM_UVMopts->vortex_radius;
	VMoptions.vortex_radius_wake_ind = MBDyn_UVLM_UVMopts->vortex_radius_wake_ind;

	// constructing the UVMopts as the input
	UVMoptions.cfl1 = MBDyn_UVLM_UVMopts->cfl1;
	UVMoptions.convection_scheme = MBDyn_UVLM_UVMopts->convection_scheme;
	UVMoptions.convect_wake = MBDyn_UVLM_UVMopts->convect_wake;
	UVMoptions.dt = MBDyn_UVLM_UVMopts->dt;
	UVMoptions.filter_method = MBDyn_UVLM_UVMopts->filter_method;
	UVMoptions.ImageMethod = MBDyn_UVLM_UVMopts->ImageMethod;
	UVMoptions.interp_coords = MBDyn_UVLM_UVMopts->interp_coords;
	UVMoptions.interp_method = MBDyn_UVLM_UVMopts->interp_method;
	UVMoptions.iterative_precond = MBDyn_UVLM_UVMopts->iterative_precond;
	UVMoptions.iterative_solver = MBDyn_UVLM_UVMopts->iterative_solver;
	UVMoptions.iterative_tol = MBDyn_UVLM_UVMopts->iterative_tol;
	UVMoptions.NumCores = MBDyn_UVLM_UVMopts->NumCores;
	UVMoptions.NumSurfaces = MBDyn_UVLM_UVMopts->NumSurfaces;
	UVMoptions.quasi_steady = MBDyn_UVLM_UVMopts->quasi_steady;
	UVMoptions.vortex_radius = MBDyn_UVLM_UVMopts->vortex_radius;
	UVMoptions.vortex_radius_wake_ind = MBDyn_UVLM_UVMopts->vortex_radius_wake_ind;
	UVMoptions.yaw_slerp = MBDyn_UVLM_UVMopts->yaw_slerp;

	// constructing the FlightConditions as input
	Flight_Conditions.c_ref = MBDyn_UVLM_FlightConditions->c_ref;
	Flight_Conditions.rho = MBDyn_UVLM_FlightConditions->rho;
	Flight_Conditions.uinf = MBDyn_UVLM_FlightConditions->uinf;
	Flight_Conditions.uinf_direction[0] = MBDyn_UVLM_FlightConditions->uinf_direction[0];
	Flight_Conditions.uinf_direction[1] = MBDyn_UVLM_FlightConditions->uinf_direction[1];
	Flight_Conditions.uinf_direction[2] = MBDyn_UVLM_FlightConditions->uinf_direction[2];
	

	// generate the pointers for the UVLM variables 
	MBDyn_UVLM_UvlmLibVar.UvlmLibVar_generate(aerogrid, time_step);

	// call the initializer function from the UVLM library
	init_UVLM(VMoptions, Flight_Conditions, MBDyn_UVLM_UvlmLibVar.p_dimensions, MBDyn_UVLM_UvlmLibVar.p_dimensions_star, 
		MBDyn_UVLM_UvlmLibVar.p_uext, MBDyn_UVLM_UvlmLibVar.p_zeta, MBDyn_UVLM_UvlmLibVar.p_zeta_star, MBDyn_UVLM_UvlmLibVar.p_zeta_dot, 
		MBDyn_UVLM_UvlmLibVar.p_zeta_star_dot, MBDyn_UVLM_UvlmLibVar.p_rbm_vel, MBDyn_UVLM_UvlmLibVar.p_gamma, MBDyn_UVLM_UvlmLibVar.p_gamma_star, 
		MBDyn_UVLM_UvlmLibVar.p_normals, MBDyn_UVLM_UvlmLibVar.p_forces);

}


extern "C" void
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
	UVLM::Types::FlightConditions& Flight_Conditions,
	double time_step) {

	if (MBDyn_UVLM_Aerogrid.aero_timestep_info.size() == 0) {    
		// This is for the very first time step (here we will use the values obtained from the init function 
		// and will pass these in the UVLM solver)
		run_UVLM(UVMoptions, Flight_Conditions, MBDyn_UVLM_UvlmLibVar.p_dimensions, MBDyn_UVLM_UvlmLibVar.p_dimensions_star, 
			MBDyn_UVLM_UvlmLibVar.i_iter, MBDyn_UVLM_UvlmLibVar.p_uext, MBDyn_UVLM_UvlmLibVar.p_uext_star, 
			MBDyn_UVLM_UvlmLibVar.p_zeta, MBDyn_UVLM_UvlmLibVar.p_zeta_star, MBDyn_UVLM_UvlmLibVar.p_zeta_dot, 
			MBDyn_UVLM_UvlmLibVar.p_rbm_vel, MBDyn_UVLM_UvlmLibVar.p_centre_rot, MBDyn_UVLM_UvlmLibVar.p_gamma, 
			MBDyn_UVLM_UvlmLibVar.p_gamma_star, MBDyn_UVLM_UvlmLibVar.p_dist_to_orig, MBDyn_UVLM_UvlmLibVar.p_normals, 
			MBDyn_UVLM_UvlmLibVar.p_forces, MBDyn_UVLM_UvlmLibVar.p_dynamic_forces);
	}
	else {
		// Generate the Grid (This function creates the grid for each and every time step)
		MBDyn_UVLM_Aerogrid.generate(MBDyn_UVLM_Aero_inputs, MBDyn_UVLM_Beam_inputs, MBDyn_UVLM_Aerogrid_settings, time_step);

		// Wake shape generator
		MBDyn_UVLM_StraightWake.StraightWake_generate(MBDyn_UVLM_Aerogrid.aero_timestep_info[time_step]);

		// Generate u_ext and u_ext_star
		MBDyn_UVLM_SteadyVelocityField.SteadyVelocityField_generate(MBDyn_UVLM_Aerogrid.aero_timestep_info[time_step], MBDyn_UVLM_UVMopts);

		// Generate the pointers for the UVLM variables
		MBDyn_UVLM_UvlmLibVar.UvlmLibVar_generate(MBDyn_UVLM_Aerogrid, time_step);

		// call the solver 
		run_UVLM(UVMoptions, Flight_Conditions, MBDyn_UVLM_UvlmLibVar.p_dimensions, MBDyn_UVLM_UvlmLibVar.p_dimensions_star, 
			MBDyn_UVLM_UvlmLibVar.i_iter, MBDyn_UVLM_UvlmLibVar.p_uext, MBDyn_UVLM_UvlmLibVar.p_uext_star, 
			MBDyn_UVLM_UvlmLibVar.p_zeta, MBDyn_UVLM_UvlmLibVar.p_zeta_star, MBDyn_UVLM_UvlmLibVar.p_zeta_dot,
			MBDyn_UVLM_UvlmLibVar.p_rbm_vel, MBDyn_UVLM_UvlmLibVar.p_centre_rot, MBDyn_UVLM_UvlmLibVar.p_gamma, 
			MBDyn_UVLM_UvlmLibVar.p_gamma_star, MBDyn_UVLM_UvlmLibVar.p_dist_to_orig, MBDyn_UVLM_UvlmLibVar.p_normals, 
			MBDyn_UVLM_UvlmLibVar.p_forces, MBDyn_UVLM_UvlmLibVar.p_dynamic_forces);
	}


	// Unsteady Contribution
	int n_surf = MBDyn_UVLM_Aerogrid.aero_timestep_info[time_step].gamma.size();
	if (!MBDyn_UVLM_UVMopts->quasi_steady) {

		//1. Copy the value of pointer gamma obtained from the "run_UVLM" function into the aerogrid timestep info gamma.
		for (int i_surf = 0; i_surf < n_surf; i_surf++) {
			for (int i = 0; i < MBDyn_UVLM_Aerogrid.aero_timestep_info[time_step].gamma[i_surf].size(); i++) {
				int N = MBDyn_UVLM_Aerogrid.aero_timestep_info[time_step].gamma[i_surf][i].size();  // Number of columns in the core matrix
				for (int j = 0; j < MBDyn_UVLM_Aerogrid.aero_timestep_info[time_step].gamma[i_surf][i].size(); j++) {
					MBDyn_UVLM_Aerogrid.aero_timestep_info[time_step].gamma[i_surf][i][j] = MBDyn_UVLM_UvlmLibVar.p_gamma[i_surf][i*N + j];
				}
			}
		}
		double dt = MBDyn_UVLM_UVMopts->dt;

		//2. pass the aerogrid into the compute_gamma_dot function and use these gamma values to compute the gamma_dot
		MBDyn_UVLM_Aerogrid.compute_gamma_dot(MBDyn_UVLM_Aero_inputs, dt, time_step);
		
		//3. gamma_dot_filtering is yet to be implemented (but not now)
		// This part is for the time being ommitted and will later be added in the future version of the code

		//4. copy the computed gamma_dot value into the pointer gamma_dot.
		std::vector<double> flatten_vec_gamma_dot;
		for (unsigned i = 0; i < n_surf; i++) {
			unsigned int NumOfPanels = (MBDyn_UVLM_Aerogrid.dimensions[i].first)*(MBDyn_UVLM_Aerogrid.dimensions[i].second);
			flatten_vec_gamma_dot = flatten(MBDyn_UVLM_Aerogrid.aero_timestep_info[time_step].gamma_dot[i]);
			for (unsigned j = 0; j < NumOfPanels; j++) {
				MBDyn_UVLM_UvlmLibVar.p_gamma_dot[i][j] = flatten_vec_gamma_dot[j];
			}
		}

		//5. Call the unsteady forces function from the UVLM lib
		calculate_unsteady_forces(UVMoptions, Flight_Conditions, MBDyn_UVLM_UvlmLibVar.p_dimensions, MBDyn_UVLM_UvlmLibVar.p_dimensions_star,
			MBDyn_UVLM_UvlmLibVar.p_zeta, MBDyn_UVLM_UvlmLibVar.p_zeta_star, MBDyn_UVLM_UvlmLibVar.p_rbm_vel, MBDyn_UVLM_UvlmLibVar.p_gamma,
			MBDyn_UVLM_UvlmLibVar.p_gamma_star, MBDyn_UVLM_UvlmLibVar.p_gamma_dot, MBDyn_UVLM_UvlmLibVar.p_normals,
			MBDyn_UVLM_UvlmLibVar.p_dynamic_forces);
	}
	else {
		for (unsigned i = 0; i < n_surf; i++) {
			unsigned int NumOfPanels = (MBDyn_UVLM_Aerogrid.dimensions[i].first)*(MBDyn_UVLM_Aerogrid.dimensions[i].second);
			for (unsigned j = 0; j < NumOfPanels; j++) {
				MBDyn_UVLM_UvlmLibVar.p_gamma_dot[i][j] = 0.0;
			}
		}
	}

	// Save the updated values of the UVLM pointer variables into the current time step "aero_timestep_info"
	MBDyn_UVLM_UvlmLibVar.UvlmLibVar_save(MBDyn_UVLM_Aerogrid, time_step);
}


extern "C" void
MBDyn_UVLM_Model_RecvFromBuf(Aerogrid& MBDyn_UVLM_Aerogrid,
	std::vector<double> &MBDyn_UVLM_CouplingKinematic,
	const unsigned &MBDyn_UVLM_NodesNum) {

	// Resize the kinematic data members of the Aerogrid class
	MBDyn_UVLM_Aerogrid.node_displacements.resize(MBDyn_UVLM_NodesNum);
	MBDyn_UVLM_Aerogrid.node_displacements_der.resize(MBDyn_UVLM_NodesNum);
	MBDyn_UVLM_Aerogrid.node_CRV.resize(MBDyn_UVLM_NodesNum);
	MBDyn_UVLM_Aerogrid.node_CRV_der.resize(MBDyn_UVLM_NodesNum);
	for (int i = 0; i < MBDyn_UVLM_NodesNum; ++i){
		MBDyn_UVLM_Aerogrid.node_displacements[i].resize(3);
		MBDyn_UVLM_Aerogrid.node_displacements_der[i].resize(3);
		MBDyn_UVLM_Aerogrid.node_CRV[i].resize(3);
		MBDyn_UVLM_Aerogrid.node_CRV_der[i].resize(3);
	}

	std::vector<std::vector<double>> Cag(3, std::vector<double>(3, 0));
	// Obtain the kinematic data
	for (unsigned i = 0; i < MBDyn_UVLM_NodesNum; ++i) {
		// The coordinates and the derivative of the coordinates are all in the "G" frame of reference. We need to 
		// convert them into "A" frame of reference. For this we multiply these by the rotation matrix "Cag".
		Cag[0][0] = MBDyn_UVLM_CouplingKinematic[9 * i + 3 * MBDyn_UVLM_NodesNum];
		Cag[1][0] = MBDyn_UVLM_CouplingKinematic[9 * i + 1 + 3 * MBDyn_UVLM_NodesNum];
		Cag[2][0] = MBDyn_UVLM_CouplingKinematic[9 * i + 2 + 3 * MBDyn_UVLM_NodesNum];
		Cag[0][1] = MBDyn_UVLM_CouplingKinematic[9 * i + 3 + 3 * MBDyn_UVLM_NodesNum];
		Cag[1][1] = MBDyn_UVLM_CouplingKinematic[9 * i + 4 + 3 * MBDyn_UVLM_NodesNum];
		Cag[2][1] = MBDyn_UVLM_CouplingKinematic[9 * i + 5 + 3 * MBDyn_UVLM_NodesNum];
		Cag[0][2] = MBDyn_UVLM_CouplingKinematic[9 * i + 6 + 3 * MBDyn_UVLM_NodesNum];
		Cag[1][2] = MBDyn_UVLM_CouplingKinematic[9 * i + 7 + 3 * MBDyn_UVLM_NodesNum];
		Cag[2][2] = MBDyn_UVLM_CouplingKinematic[9 * i + 8 + 3 * MBDyn_UVLM_NodesNum];

		//- Filling up the coordinates of the structural beam nodes 
		MBDyn_UVLM_Aerogrid.node_displacements[i] = MatVecMul(Cag, { MBDyn_UVLM_CouplingKinematic[3 * i],
			MBDyn_UVLM_CouplingKinematic[3 * i + 1], MBDyn_UVLM_CouplingKinematic[3 * i + 2] });

		//- Filling up the derivatives of the coordinates of the structural beam nodes
		MBDyn_UVLM_Aerogrid.node_displacements_der[i] = MatVecMul(Cag, { MBDyn_UVLM_CouplingKinematic[3 * i + 12 * MBDyn_UVLM_NodesNum],
			MBDyn_UVLM_CouplingKinematic[3 * i + 1 + 12 * MBDyn_UVLM_NodesNum], MBDyn_UVLM_CouplingKinematic[3 * i + 2 + 12 * MBDyn_UVLM_NodesNum] });

		//- Filling up the Cga matrix for each node 
		MBDyn_UVLM_Aerogrid.node_cga[i][0][0] = Cag[0][0];
		MBDyn_UVLM_Aerogrid.node_cga[i][0][1] = Cag[1][0];
		MBDyn_UVLM_Aerogrid.node_cga[i][0][2] = Cag[2][0];
		MBDyn_UVLM_Aerogrid.node_cga[i][1][0] = Cag[0][1];
		MBDyn_UVLM_Aerogrid.node_cga[i][1][1] = Cag[1][1];
		MBDyn_UVLM_Aerogrid.node_cga[i][1][2] = Cag[2][1];
		MBDyn_UVLM_Aerogrid.node_cga[i][2][0] = Cag[0][2];
		MBDyn_UVLM_Aerogrid.node_cga[i][2][1] = Cag[1][2];
		MBDyn_UVLM_Aerogrid.node_cga[i][2][2] = Cag[2][2];

		//- Filling up the Cartesian rotation vector for each node of each element
		MBDyn_UVLM_Aerogrid.node_CRV[i] = rotation2crv(Cag);

		//- Filling up the derivatives of the cartesian rotation vector of each node of each element
		MBDyn_UVLM_Aerogrid.node_CRV_der[i] = { 0.0, 0.0, 0.0 };
	}
}

extern "C" void
MBDyn_UVLM_Model_SendToBuf(Aerogrid& MBDyn_UVLM_Aerogrid,
	std::vector<double>& MBDyn_UVLM_CouplingDynamic,
	const unsigned &MBDyn_UVLM_NodesNum,
	double time_step) {

	// The values returned are the forces, dynamic forces and moments at each node .
	for (unsigned i = 0; i < MBDyn_UVLM_NodesNum; ++i) {
		if (MBDyn_UVLM_Aerogrid.aero_timestep_info.size() == 0) {
			MBDyn_UVLM_CouplingDynamic[3 * i] = 0.0;
			MBDyn_UVLM_CouplingDynamic[3 * i + 1] = 0.0;
			MBDyn_UVLM_CouplingDynamic[3 * i + 2] = 0.0;
			MBDyn_UVLM_CouplingDynamic[3 * i + 3 * MBDyn_UVLM_NodesNum] = 0.0;
			MBDyn_UVLM_CouplingDynamic[3 * i + 1 + 3 * MBDyn_UVLM_NodesNum] = 0.0;
			MBDyn_UVLM_CouplingDynamic[3 * i + 2 + 3 * MBDyn_UVLM_NodesNum] = 0.0;
			MBDyn_UVLM_CouplingDynamic[3 * i + 6 * MBDyn_UVLM_NodesNum] = 0.0;
			MBDyn_UVLM_CouplingDynamic[3 * i + 1 + 6 * MBDyn_UVLM_NodesNum] = 0.0;
			MBDyn_UVLM_CouplingDynamic[3 * i + 2 + 6 * MBDyn_UVLM_NodesNum] = 0.0;
		}
		else {
			MBDyn_UVLM_CouplingDynamic[3 * i] = 1.0;
			MBDyn_UVLM_CouplingDynamic[3 * i + 1] = 1.0;
			MBDyn_UVLM_CouplingDynamic[3 * i + 2] = 1.0;
			MBDyn_UVLM_CouplingDynamic[3 * i + 3 * MBDyn_UVLM_NodesNum] = 1.0;
			MBDyn_UVLM_CouplingDynamic[3 * i + 1 + 3 * MBDyn_UVLM_NodesNum] = 1.0;
			MBDyn_UVLM_CouplingDynamic[3 * i + 2 + 3 * MBDyn_UVLM_NodesNum] = 1.0;
			MBDyn_UVLM_CouplingDynamic[3 * i + 6 * MBDyn_UVLM_NodesNum] = 1.0;
			MBDyn_UVLM_CouplingDynamic[3 * i + 1 + 6 * MBDyn_UVLM_NodesNum] = 1.0;
			MBDyn_UVLM_CouplingDynamic[3 * i + 2 + 6 * MBDyn_UVLM_NodesNum] = 1.0;
		}
	}
}



/* Class member function definitions [START] */ 
AeroTimeStepInfo::AeroTimeStepInfo() {
	NO_OP;
}

void
AeroTimeStepInfo::initialize(std::vector<std::pair<unsigned int, unsigned int>>& dimensions,
	std::vector<std::pair<unsigned int, unsigned int>>& dimensions_star) {

	unsigned int number_of_surfaces = dimensions.size();

	zeta.resize(number_of_surfaces);
	for (int i = 0; i < number_of_surfaces; ++i){
		zeta[i].resize(3);
		for (int j = 0; j < 3; ++j) {
			zeta[i][j].resize(dimensions[i].first + 1);
			for (int k = 0; k < dimensions[i].first + 1; ++k) {
				zeta[i][j][k].resize(dimensions[i].second + 1);
			}
		}
	}

	zeta_dot.resize(number_of_surfaces);
	for (int i = 0; i < number_of_surfaces; ++i){
		zeta_dot[i].resize(3);
		for (int j = 0; j < 3; ++j) {
			zeta_dot[i][j].resize(dimensions[i].first + 1);
			for (int k = 0; k < dimensions[i].first + 1; ++k) {
				zeta_dot[i][j][k].resize(dimensions[i].second + 1);
			}
		}
	}

	normals.resize(number_of_surfaces);
	for (int i = 0; i < number_of_surfaces; ++i){
		normals[i].resize(3);
		for (int j = 0; j < 3; ++j) {
			normals[i][j].resize(dimensions[i].first);
			for (int k = 0; k < dimensions[i].first; ++k) {
				normals[i][j][k].resize(dimensions[i].second);
			}
		}
	}

	forces.resize(number_of_surfaces);
	for (int i = 0; i < number_of_surfaces; ++i){
		forces[i].resize(6);
		for (int j = 0; j < 6; ++j) {
			forces[i][j].resize(dimensions[i].first + 1);
			for (int k = 0; k < dimensions[i].first + 1; ++k) {
				forces[i][j][k].resize(dimensions[i].second + 1);
			}
		}
	}

	dynamic_forces.resize(number_of_surfaces);
	for (int i = 0; i < number_of_surfaces; ++i){
		dynamic_forces[i].resize(6);
		for (int j = 0; j < 6; ++j) {
			dynamic_forces[i][j].resize(dimensions[i].first + 1);
			for (int k = 0; k < dimensions[i].first + 1; ++k) {
				dynamic_forces[i][j][k].resize(dimensions[i].second + 1);
			}
		}
	}

	zeta_star.resize(number_of_surfaces);
	for (int i = 0; i < number_of_surfaces; ++i){
		zeta_star[i].resize(3);
		for (int j = 0; j < 3; ++j) {
			zeta_star[i][j].resize(dimensions_star[i].first + 1);
			for (int k = 0; k < dimensions_star[i].first + 1; ++k) {
				zeta_star[i][j][k].resize(dimensions_star[i].second + 1);
			}
		}
	}

	u_ext.resize(number_of_surfaces);
	for (int i = 0; i < number_of_surfaces; ++i){
		u_ext[i].resize(3);
		for (int j = 0; j < 3; ++j) {
			u_ext[i][j].resize(dimensions[i].first + 1);
			for (int k = 0; k < dimensions[i].first + 1; ++k) {
				u_ext[i][j][k].resize(dimensions[i].second + 1);
			}
		}
	}

	u_ext_star.resize(number_of_surfaces);
	for (int i = 0; i < number_of_surfaces; ++i){
		u_ext_star[i].resize(3);
		for (int j = 0; j < 3; ++j) {
			u_ext_star[i][j].resize(dimensions_star[i].first + 1);
			for (int k = 0; k < dimensions_star[i].first + 1; ++k) {
				u_ext_star[i][j][k].resize(dimensions_star[i].second + 1);
			}
		}
	}

	gamma.resize(number_of_surfaces);
	for (int i = 0; i < number_of_surfaces; ++i){
		gamma[i].resize(dimensions[i].first);
		for (int j = 0; j < dimensions[i].first; ++j) {
			gamma[i][j].resize(dimensions[i].second);
		}
	}

	gamma_star.resize(number_of_surfaces);
	for (int i = 0; i < number_of_surfaces; ++i){
		gamma_star[i].resize(dimensions_star[i].first);
		for (int j = 0; j < dimensions_star[i].first; ++j) {
			gamma_star[i][j].resize(dimensions_star[i].second);
		}
	}

	gamma_dot.resize(number_of_surfaces);
	for (int i = 0; i < number_of_surfaces; ++i){
		gamma_dot[i].resize(dimensions[i].first);
		for (int j = 0; j < dimensions[i].first; ++j) {
			gamma_dot[i][j].resize(dimensions[i].second);
		}
	}

	dist_to_orig.resize(number_of_surfaces);
	for (int i = 0; i < number_of_surfaces; ++i){
		dist_to_orig[i].resize(dimensions[i].first + 1);
		for (int j = 0; j < dimensions[i].first + 1; ++j) {
			dist_to_orig[i][j].resize(dimensions[i].second + 1);
		}
	}
}

AeroTimeStepInfo::~AeroTimeStepInfo() {
	NO_OP;
}



Beam_inputs::Beam_inputs() {
	NO_OP;
}

Beam_inputs::~Beam_inputs() {
	NO_OP;
}



Aero_inputs::Aero_inputs() {
	NO_OP;
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




Aerogrid::Aerogrid() {

	NO_OP;
}

void 
Aerogrid::generate(Aero_inputs& aero_inputs, Beam_inputs& beam_inputs, const Aerogrid_settings* aerogrid_settings,
	double ts) {

	// number of total nodes (structural + aero&struct)
	n_node = aero_inputs.aero_node_input.size();

	// number of elements
	n_elem = aero_inputs.surface_distribution_input.size();

	// Surface distribution
	surface_distribution = aero_inputs.surface_distribution_input;

	// number of surfaces
	n_surf = aero_inputs.surface_m.size();

	// number of chordwise panels
	surface_m = aero_inputs.surface_m;

	// number of aero nodes
	n_aero_node = aero_inputs.aero_node_input.size();

	// get N per surface
	calculate_dimensions(aero_inputs, beam_inputs, aerogrid_settings);

	// Write grid info on the screen
	output_info(aero_inputs);

	// allocating the initial storage
	aero_ini_info.initialize(dimensions, dimensions_star);

	add_timestep();
	generate_mapping(aero_inputs, beam_inputs);
	generate_zeta(aero_inputs, beam_inputs, aerogrid_settings, ts);
}

void 
Aerogrid::calculate_dimensions(Aero_inputs& aero_inputs, Beam_inputs& beam_inputs,
	const Aerogrid_settings* aerogrid_settings) {

	dimensions.resize(n_surf);
	dimensions_star.resize(n_surf);

	for (int i = 0; i < n_surf; ++i) {
		// adding M values
		dimensions[i].first = surface_m[i];
	}
	// Count N values (actually the count result will be N+1)
	std::vector<std::vector<int>> nodes_in_surface;
	for (int i_elem = 0; i_elem < beam_inputs.num_elem; ++i_elem) {
		std::vector<int> nodes;
		nodes = beam_inputs.connectivities[i_elem];
		int i_surf = aero_inputs.surface_distribution_input[i_elem];
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
			if (aero_inputs.aero_node_input[nodes[i_global_node]]) {
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

void 
Aerogrid::output_info(Aero_inputs& aero_inputs) {

	int _n_surf = aero_inputs.surface_m.size(); 
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

void 
Aerogrid::add_timestep() {

	if (aero_timestep_info.empty()) {
		aero_timestep_info.push_back(aero_ini_info);
	}
	else {
		aero_timestep_info.push_back(aero_timestep_info.back());
	}
}

void 
Aerogrid::generate_zeta_timestep_info(Aero_inputs& aero_inputs, Beam_inputs& beam_inputs, AeroTimeStepInfo& aero_tstep,
	const Aerogrid_settings* aerogrid_settings) {

	std::vector<std::vector<int>> global_node_in_surface;
	global_node_in_surface.resize(n_surf);

	// One surface per element
	for (int i_elem = 0; i_elem < n_elem; ++i_elem) {
		int i_surf = aero_inputs.surface_distribution_input[i_elem];
		// Check if we have to generate a surface here
		if (i_surf == -1) {
			continue;
		}

		for (int i_local_node = 0; i_local_node < beam_inputs.connectivities[i_elem].size(); ++i_local_node) {
			int i_global_node = beam_inputs.connectivities[i_elem][i_local_node];
			if (aero_inputs.aero_node_input[i_global_node] == false) {
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
			node_info.chord = aero_inputs.chords[i_elem][i_local_node];
			node_info.eaxis = aero_inputs.elastic_axis[i_elem][i_local_node];
			node_info.twist = aero_inputs.twist[i_elem][i_local_node];
			node_info.sweep = aero_inputs.sweep[i_elem][i_local_node];
			node_info.M = dimensions[i_surf].first;
			node_info.M_distribution = aero_inputs.m_distribution;
			node_info.airfoil = aero_inputs.airfoil_distribution_input[i_elem][i_local_node];
			node_info.beam_coord = node_displacements[i_global_node];
			node_info.pos_dot = node_displacements_der[i_global_node];
			node_info.beam_psi = node_CRV[beam_inputs.connectivities[i_elem][i_local_node]];
			node_info.psi_dot = node_CRV_der[beam_inputs.connectivities[i_elem][i_local_node]];
			node_info.for_delta = beam_inputs.frame_of_reference_delta[i_elem][i_local_node];
			node_info.cga = node_cga[i_global_node];

			generate_strip(aero_inputs, aerogrid_settings, aero_tstep);
		}
	}
}

void 
Aerogrid::generate_zeta(Aero_inputs& aero_inputs, Beam_inputs& beam_inputs, const Aerogrid_settings* aerogrid_settings,
	double ts) {
	generate_zeta_timestep_info(aero_inputs, beam_inputs, aero_timestep_info[ts], aerogrid_settings);
}

void 
Aerogrid::generate_mapping(Aero_inputs& aero_inputs, Beam_inputs& beam_inputs) {

	struct2aero_mapping.resize(n_node);
	std::vector<int> surf_n_counter(n_surf, 0);
	std::vector<std::vector<int>> nodes_in_surface;
	nodes_in_surface.resize(n_surf);

	std::vector<int> global_nodes;
	for (int i_elem = 0; i_elem < n_elem; ++i_elem) {
		int i_surf = aero_inputs.surface_distribution_input[i_elem];
		if (i_surf == -1) {
			continue;
		}
		global_nodes = beam_inputs.reordered_connectivities[i_elem];
		for (int i_global_node = 0; i_global_node < global_nodes.size(); ++i_global_node) {
			if (aero_inputs.aero_node_input[global_nodes[i_global_node]] == false) {
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
		global_nodes = beam_inputs.connectivities[i_elem];
		for (int i_global_node = 0; i_global_node < global_nodes.size(); ++i_global_node) {
			for (int i = 0; i < struct2aero_mapping[global_nodes[i_global_node]].size(); ++i) {  // this loop is if a single node belongs to multiple surfaces
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

void 
Aerogrid::compute_gamma_dot(Aero_inputs& aero_inputs, double dt, double time_step) {

	// Computes the temporal derivative of circulation (gamma) using finite differences.
	// It will use a first order approximation for the first evaluation and second order for the subsequent evaluations.

	int n_surf = aero_inputs.surface_m.size();

	if (aero_timestep_info.size() <= 2) {
		for (int i_surf = 0; i_surf < n_surf; ++i_surf) {
			for (int j = 0; j < aero_timestep_info[time_step].gamma_dot[i_surf].size(); ++j) {
				for (int k = 0; k < aero_timestep_info[time_step].gamma_dot[i_surf][j].size(); ++k) {
					aero_timestep_info[time_step].gamma_dot[i_surf][j][k] = 0.0;
				}
			}
		}
	}
	else {
		for (int i_surf = 0; i_surf < n_surf; ++i_surf) {
			for (int j = 0; j < aero_timestep_info[time_step].gamma_dot[i_surf].size(); ++j) {
				for (int k = 0; k < aero_timestep_info[time_step].gamma_dot[i_surf][j].size(); ++k) {
					aero_timestep_info[time_step].gamma_dot[i_surf][j][k] = (aero_timestep_info[time_step].gamma[i_surf][j][k] - 
						aero_timestep_info[time_step - 2].gamma[i_surf][j][k]) / dt;
				}
			}
		}
	}
}

void 
Aerogrid::generate_strip(Aero_inputs& aero_inputs, const Aerogrid_settings* aerogrid_settings, AeroTimeStepInfo& aero_tstep) {

	// Returns a strip of panels in "G" frame of reference, it has to be then rotated to simulate angles of attack, etc
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
	airfoils_coords = aero_inputs.airfoils[node_info.airfoil];
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

	// Rotation from A to G frame of reference
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

Aerogrid::~Aerogrid() {

	NO_OP;
}



UvlmLibVar::UvlmLibVar() {

	NO_OP;
}

void
UvlmLibVar::UvlmLibVar_generate(Aerogrid& aerogrid, double time_step) {

	// constructing the "p_dimensions" pointer
	p_dimensions = new unsigned int*[aerogrid.dimensions.size()];
	for (unsigned i = 0; i < aerogrid.dimensions.size(); i++) {
		p_dimensions[i] = new unsigned int[2];
		p_dimensions[i][0] = aerogrid.dimensions[i].first;
		p_dimensions[i][1] = aerogrid.dimensions[i].second;
	}

	// constructing the "p_dimensions_star" pointer
	p_dimensions_star = new unsigned int*[aerogrid.dimensions_star.size()];
	for (unsigned i = 0; i < aerogrid.dimensions_star.size(); i++) {
		p_dimensions_star[i] = new unsigned int[2];
		p_dimensions_star[i][0] = aerogrid.dimensions_star[i].first;
		p_dimensions_star[i][1] = aerogrid.dimensions_star[i].second;
	}

	// constructing the "i_iter"
	i_iter = static_cast<unsigned int>(time_step);

	// constructing the "p_uext" pointer
	int n_surf = aerogrid.aero_timestep_info[time_step].u_ext.size();
	p_uext = new double*[n_surf * 3];
	std::vector<double> flatten_vec_x;
	std::vector<double> flatten_vec_y;
	std::vector<double> flatten_vec_z;
	for (unsigned i = 0; i < n_surf; i++) {
		unsigned int NumOfGridVertices = (aerogrid.dimensions[i].first + 1)*(aerogrid.dimensions[i].second + 1);
		p_uext[3 * i] = new double[NumOfGridVertices];
		p_uext[3 * i + 1] = new double[NumOfGridVertices];
		p_uext[3 * i + 2] = new double[NumOfGridVertices];
		flatten_vec_x = flatten(aerogrid.aero_timestep_info[time_step].u_ext[i][0]);
		flatten_vec_y = flatten(aerogrid.aero_timestep_info[time_step].u_ext[i][1]);
		flatten_vec_z = flatten(aerogrid.aero_timestep_info[time_step].u_ext[i][2]);
		for (unsigned j = 0; j < NumOfGridVertices; j++) {
			p_uext[3 * i][j] = flatten_vec_x[j];
			p_uext[3 * i + 1][j] = flatten_vec_y[j];
			p_uext[3 * i + 2][j] = flatten_vec_z[j];
		}
	}

	// constructing the "p_uext_star" pointer
	n_surf = aerogrid.aero_timestep_info[time_step].u_ext_star.size();
	p_uext_star = new double*[n_surf * 3];
	for (unsigned i = 0; i < n_surf; i++) {
		unsigned int NumOfWakeGridVertices = (aerogrid.dimensions_star[i].first + 1)*(aerogrid.dimensions_star[i].second + 1);
		p_uext_star[3 * i] = new double[NumOfWakeGridVertices];
		p_uext_star[3 * i + 1] = new double[NumOfWakeGridVertices];
		p_uext_star[3 * i + 2] = new double[NumOfWakeGridVertices];
		flatten_vec_x = flatten(aerogrid.aero_timestep_info[time_step].u_ext_star[i][0]);
		flatten_vec_y = flatten(aerogrid.aero_timestep_info[time_step].u_ext_star[i][1]);
		flatten_vec_z = flatten(aerogrid.aero_timestep_info[time_step].u_ext_star[i][2]);
		for (unsigned j = 0; j < NumOfWakeGridVertices; j++) {
			p_uext_star[3 * i][j] = flatten_vec_x[j];
			p_uext_star[3 * i + 1][j] = flatten_vec_y[j];
			p_uext_star[3 * i + 2][j] = flatten_vec_z[j];
		}
	}

	// constructing the "p_zeta" pointer
	n_surf = aerogrid.aero_timestep_info[time_step].zeta.size();
	p_zeta = new double*[n_surf * 3];
	for (unsigned i = 0; i < n_surf; i++) {
		unsigned int NumOfGridVertices = (aerogrid.dimensions[i].first + 1)*(aerogrid.dimensions[i].second + 1);
		p_zeta[3 * i] = new double[NumOfGridVertices];
		p_zeta[3 * i + 1] = new double[NumOfGridVertices];
		p_zeta[3 * i + 2] = new double[NumOfGridVertices];
		flatten_vec_x = flatten(aerogrid.aero_timestep_info[time_step].zeta[i][0]);
		flatten_vec_y = flatten(aerogrid.aero_timestep_info[time_step].zeta[i][1]);
		flatten_vec_z = flatten(aerogrid.aero_timestep_info[time_step].zeta[i][2]);
		for (unsigned j = 0; j < NumOfGridVertices; j++) {
			p_zeta[3 * i][j] = flatten_vec_x[j];
			p_zeta[3 * i + 1][j] = flatten_vec_y[j];
			p_zeta[3 * i + 2][j] = flatten_vec_z[j];
		}
	}

	// constructing the "p_zeta_star" pointer
	n_surf = aerogrid.aero_timestep_info[time_step].zeta_star.size();
	p_zeta_star = new double*[n_surf * 3];
	for (unsigned i = 0; i < n_surf; i++) {
		unsigned int NumOfWakeGridVertices = (aerogrid.dimensions_star[i].first + 1)*(aerogrid.dimensions_star[i].second + 1);
		p_zeta_star[3 * i] = new double[NumOfWakeGridVertices];
		p_zeta_star[3 * i + 1] = new double[NumOfWakeGridVertices];
		p_zeta_star[3 * i + 2] = new double[NumOfWakeGridVertices];
		flatten_vec_x = flatten(aerogrid.aero_timestep_info[time_step].zeta_star[i][0]);
		flatten_vec_y = flatten(aerogrid.aero_timestep_info[time_step].zeta_star[i][1]);
		flatten_vec_z = flatten(aerogrid.aero_timestep_info[time_step].zeta_star[i][2]);
		for (unsigned j = 0; j < NumOfWakeGridVertices; j++) {
			p_zeta_star[3 * i][j] = flatten_vec_x[j];
			p_zeta_star[3 * i + 1][j] = flatten_vec_y[j];
			p_zeta_star[3 * i + 2][j] = flatten_vec_z[j];
		}
	}

	// constructing the "p_zeta_dot" pointer
	n_surf = aerogrid.aero_timestep_info[time_step].zeta_dot.size();
	p_zeta_dot = new double*[n_surf * 3];
	for (unsigned i = 0; i < n_surf; i++) {
		unsigned int NumOfGridVertices = (aerogrid.dimensions[i].first + 1)*(aerogrid.dimensions[i].second + 1);
		p_zeta_dot[3 * i] = new double[NumOfGridVertices];
		p_zeta_dot[3 * i + 1] = new double[NumOfGridVertices];
		p_zeta_dot[3 * i + 2] = new double[NumOfGridVertices];
		flatten_vec_x = flatten(aerogrid.aero_timestep_info[time_step].zeta_dot[i][0]);
		flatten_vec_y = flatten(aerogrid.aero_timestep_info[time_step].zeta_dot[i][1]);
		flatten_vec_z = flatten(aerogrid.aero_timestep_info[time_step].zeta_dot[i][2]);
		for (unsigned j = 0; j < NumOfGridVertices; j++) {
			p_zeta_dot[3 * i][j] = flatten_vec_x[j];
			p_zeta_dot[3 * i + 1][j] = flatten_vec_y[j];
			p_zeta_dot[3 * i + 2][j] = flatten_vec_z[j];
		}
	}

	// constructing the "p_zeta_star_dot" pointer
	n_surf = aerogrid.aero_timestep_info[time_step].zeta_star.size();
	p_zeta_star_dot = new double*[n_surf * 3];
	for (unsigned i = 0; i < n_surf; i++) {
		unsigned int NumOfGridVertices = (aerogrid.dimensions_star[i].first + 1)*(aerogrid.dimensions_star[i].second + 1);
		p_zeta_star_dot[3 * i] = new double[NumOfGridVertices];
		p_zeta_star_dot[3 * i + 1] = new double[NumOfGridVertices];
		p_zeta_star_dot[3 * i + 2] = new double[NumOfGridVertices];
		for (unsigned j = 0; j < NumOfGridVertices; j++) {
			p_zeta_star_dot[3 * i][j] = 0.0;
			p_zeta_star_dot[3 * i + 1][j] = 0.0;
			p_zeta_star_dot[3 * i + 2][j] = 0.0;
		}
	}

	// constructing the "p_rbm_vel" pointer
	p_rbm_vel = new double[6];
	for (unsigned i = 0; i < 6; i++) {
		p_rbm_vel[i] = 0.0;     // for now it is considered to be 0 (has to be included if the dynamic simulation is required)
	}

	// constructing the "p_centre_rot" pointer
	p_centre_rot = new double[3];
	p_centre_rot[0] = 0.0;
	p_centre_rot[1] = 0.0;
	p_centre_rot[2] = 0.0;

	// constructing the "p_gamma" pointer
	n_surf = aerogrid.aero_timestep_info[time_step].gamma.size();
	p_gamma = new double*[n_surf];
	std::vector<double> flatten_vec_gamma;
	for (unsigned i = 0; i < n_surf; i++) {
		unsigned int NumOfPanels = (aerogrid.dimensions[i].first)*(aerogrid.dimensions[i].second);
		p_gamma[i] = new double[NumOfPanels];
		flatten_vec_gamma = flatten(aerogrid.aero_timestep_info[time_step].gamma[i]);
		for (unsigned j = 0; j < NumOfPanels; j++) {
			p_gamma[i][j] = flatten_vec_gamma[j];
		}
	}

	// constructing the "p_gamma_star" pointer
	n_surf = aerogrid.aero_timestep_info[time_step].gamma_star.size();
	p_gamma_star = new double*[n_surf];
	for (unsigned i = 0; i < n_surf; i++) {
		unsigned int NumOfWakePanels = (aerogrid.dimensions_star[i].first)*(aerogrid.dimensions_star[i].second);
		p_gamma_star[i] = new double[NumOfWakePanels];
		flatten_vec_gamma = flatten(aerogrid.aero_timestep_info[time_step].gamma_star[i]);
		for (unsigned j = 0; j < NumOfWakePanels; j++) {
			p_gamma_star[i][j] = flatten_vec_gamma[j];
		}
	}

	// constructing the "p_gamma_dot" pointer
	n_surf = aerogrid.aero_timestep_info[time_step].gamma_dot.size();
	p_gamma_dot = new double*[n_surf];
	for (unsigned i = 0; i < n_surf; i++) {
		unsigned int NumOfPanels = (aerogrid.dimensions[i].first)*(aerogrid.dimensions[i].second);
		p_gamma_dot[i] = new double[NumOfPanels];
		flatten_vec_gamma = flatten(aerogrid.aero_timestep_info[time_step].gamma_dot[i]);
		for (unsigned j = 0; j < NumOfPanels; j++) {
			p_gamma_dot[i][j] = flatten_vec_gamma[j];
		}
	}

	// constructing the "p_dist_to_orig" pointer
	// Distance from the trailing edge of the wake vertices
	n_surf = aerogrid.aero_timestep_info[time_step].dist_to_orig.size();
	p_dist_to_orig = new double*[n_surf];
	std::vector<double> flatten_vec_dist_to_orig;
	for (unsigned i = 0; i < n_surf; i++) {
		unsigned int NumOfWakeGridVertices = (aerogrid.dimensions_star[i].first + 1)*(aerogrid.dimensions_star[i].second + 1);
		p_dist_to_orig[i] = new double[NumOfWakeGridVertices];
		flatten_vec_dist_to_orig = flatten(aerogrid.aero_timestep_info[time_step].dist_to_orig[i]);
		for (unsigned j = 0; j < NumOfWakeGridVertices; j++) {
			p_dist_to_orig[i][j] = flatten_vec_dist_to_orig[j];
		}
	}

	// constructing the "p_normals" pointer
	n_surf = aerogrid.aero_timestep_info[time_step].normals.size();
	p_normals = new double*[n_surf * 3];
	for (unsigned i = 0; i < n_surf; i++) {
		unsigned int NumOfPanels = (aerogrid.dimensions[i].first)*(aerogrid.dimensions[i].second);
		p_normals[3 * i] = new double[NumOfPanels];
		p_normals[3 * i + 1] = new double[NumOfPanels];
		p_normals[3 * i + 2] = new double[NumOfPanels];
		flatten_vec_x = flatten(aerogrid.aero_timestep_info[time_step].normals[i][0]);
		flatten_vec_y = flatten(aerogrid.aero_timestep_info[time_step].normals[i][1]);
		flatten_vec_z = flatten(aerogrid.aero_timestep_info[time_step].normals[i][2]);
		for (unsigned j = 0; j < NumOfPanels; j++) {
			p_normals[3 * i][j] = flatten_vec_x[j];
			p_normals[3 * i + 1][j] = flatten_vec_y[j];
			p_normals[3 * i + 2][j] = flatten_vec_z[j];
		}
	}

	// constructing the "p_forces" pointer
	n_surf = aerogrid.aero_timestep_info[time_step].forces.size();
	p_forces = new double*[n_surf * 6];
	std::vector<double> flatten_vec_x_prime;
	std::vector<double> flatten_vec_y_prime;
	std::vector<double> flatten_vec_z_prime;
	for (unsigned i = 0; i < n_surf; i++) {
		unsigned int NumOfGridVertices = (aerogrid.dimensions[i].first + 1)*(aerogrid.dimensions[i].second + 1);
		p_forces[3 * i] = new double[NumOfGridVertices];
		p_forces[3 * i + 1] = new double[NumOfGridVertices];
		p_forces[3 * i + 2] = new double[NumOfGridVertices];
		p_forces[3 * i + 3] = new double[NumOfGridVertices];
		p_forces[3 * i + 4] = new double[NumOfGridVertices];
		p_forces[3 * i + 5] = new double[NumOfGridVertices];
		flatten_vec_x = flatten(aerogrid.aero_timestep_info[time_step].forces[i][0]);
		flatten_vec_y = flatten(aerogrid.aero_timestep_info[time_step].forces[i][1]);
		flatten_vec_z = flatten(aerogrid.aero_timestep_info[time_step].forces[i][2]);
		flatten_vec_x_prime = flatten(aerogrid.aero_timestep_info[time_step].forces[i][3]);
		flatten_vec_y_prime = flatten(aerogrid.aero_timestep_info[time_step].forces[i][4]);
		flatten_vec_z_prime = flatten(aerogrid.aero_timestep_info[time_step].forces[i][5]);
		for (unsigned j = 0; j < NumOfGridVertices; j++) {
			p_forces[3 * i][j] = flatten_vec_x[j];
			p_forces[3 * i + 1][j] = flatten_vec_y[j];
			p_forces[3 * i + 2][j] = flatten_vec_z[j];
			p_forces[3 * i + 3][j] = flatten_vec_x_prime[j];
			p_forces[3 * i + 4][j] = flatten_vec_y_prime[j];
			p_forces[3 * i + 5][j] = flatten_vec_z_prime[j];
		}
	}

	// constructing the "p_dynamic_forces" pointer
	n_surf = aerogrid.aero_timestep_info[time_step].dynamic_forces.size();
	p_dynamic_forces = new double*[n_surf * 6];
	for (unsigned i = 0; i < n_surf; i++) {
		unsigned int NumOfGridVertices = (aerogrid.dimensions[i].first + 1)*(aerogrid.dimensions[i].second + 1);
		p_dynamic_forces[3 * i] = new double[NumOfGridVertices];
		p_dynamic_forces[3 * i + 1] = new double[NumOfGridVertices];
		p_dynamic_forces[3 * i + 2] = new double[NumOfGridVertices];
		p_dynamic_forces[3 * i + 3] = new double[NumOfGridVertices];
		p_dynamic_forces[3 * i + 4] = new double[NumOfGridVertices];
		p_dynamic_forces[3 * i + 5] = new double[NumOfGridVertices];
		flatten_vec_x = flatten(aerogrid.aero_timestep_info[time_step].dynamic_forces[i][0]);
		flatten_vec_y = flatten(aerogrid.aero_timestep_info[time_step].dynamic_forces[i][1]);
		flatten_vec_z = flatten(aerogrid.aero_timestep_info[time_step].dynamic_forces[i][2]);
		flatten_vec_x_prime = flatten(aerogrid.aero_timestep_info[time_step].dynamic_forces[i][3]);
		flatten_vec_y_prime = flatten(aerogrid.aero_timestep_info[time_step].dynamic_forces[i][4]);
		flatten_vec_z_prime = flatten(aerogrid.aero_timestep_info[time_step].dynamic_forces[i][5]);
		for (unsigned j = 0; j < NumOfGridVertices; j++) {
			p_dynamic_forces[3 * i][j] = flatten_vec_x[j];
			p_dynamic_forces[3 * i + 1][j] = flatten_vec_y[j];
			p_dynamic_forces[3 * i + 2][j] = flatten_vec_z[j];
			p_dynamic_forces[3 * i + 3][j] = flatten_vec_x_prime[j];
			p_dynamic_forces[3 * i + 4][j] = flatten_vec_y_prime[j];
			p_dynamic_forces[3 * i + 5][j] = flatten_vec_z_prime[j];
		}
	}
}

void
UvlmLibVar::UvlmLibVar_save(Aerogrid& aerogrid, double time_step) {

	// saving the values of "p_uext"
	int n_surf = aerogrid.aero_timestep_info[time_step].u_ext.size();
	for (int i_surf = 0; i_surf < n_surf; i_surf++) {
		for (int i = 0; i < aerogrid.aero_timestep_info[time_step].u_ext[i_surf].size(); i++) {
			for (int j = 0; j < aerogrid.aero_timestep_info[time_step].u_ext[i_surf][i].size(); j++) {
				int N = aerogrid.aero_timestep_info[time_step].u_ext[i_surf][i][j].size();  // Number of columns in the core matrix
				for (int k = 0; k < aerogrid.aero_timestep_info[time_step].u_ext[i_surf][i][j].size(); k++) {
					aerogrid.aero_timestep_info[time_step].u_ext[i_surf][i][j][k] = p_uext[i + i_surf * n_surf][j*N + k];
				}
			}
		}
	}

	// saving the value of "p_uext_star" and "p_zeta_star"
	n_surf = aerogrid.aero_timestep_info[time_step].u_ext_star.size();
	for (int i_surf = 0; i_surf < n_surf; i_surf++) {
		for (int i = 0; i < aerogrid.aero_timestep_info[time_step].u_ext_star[i_surf].size(); i++) {
			for (int j = 0; j < aerogrid.aero_timestep_info[time_step].u_ext_star[i_surf][i].size(); j++) {
				int N = aerogrid.aero_timestep_info[time_step].u_ext_star[i_surf][i][j].size();  // Number of columns in the core matrix
				for (int k = 0; k < aerogrid.aero_timestep_info[time_step].u_ext_star[i_surf][i][j].size(); k++) {
					aerogrid.aero_timestep_info[time_step].u_ext_star[i_surf][i][j][k] = p_uext_star[i + i_surf * n_surf][j*N + k];
					aerogrid.aero_timestep_info[time_step].zeta_star[i_surf][i][j][k] = p_zeta_star[i + i_surf * n_surf][j*N + k];
				}
			}
		}
	}

	// saving the value of "p_zeta" and "p_zeta_dot"
	n_surf = aerogrid.aero_timestep_info[time_step].zeta.size();
	for (int i_surf = 0; i_surf < n_surf; i_surf++) {
		for (int i = 0; i < aerogrid.aero_timestep_info[time_step].zeta[i_surf].size(); i++) {
			for (int j = 0; j < aerogrid.aero_timestep_info[time_step].zeta[i_surf][i].size(); j++) {
				int N = aerogrid.aero_timestep_info[time_step].zeta[i_surf][i][j].size();  // Number of columns in the core matrix
				for (int k = 0; k < aerogrid.aero_timestep_info[time_step].zeta[i_surf][i][j].size(); k++) {
					aerogrid.aero_timestep_info[time_step].zeta[i_surf][i][j][k] = p_zeta[i + i_surf * n_surf][j*N + k];
					aerogrid.aero_timestep_info[time_step].zeta_dot[i_surf][i][j][k] = p_zeta_dot[i + i_surf * n_surf][j*N + k];
				}
			}
		}
	}

	// saving the values of "p_gamma" and "p_gamma_dot"
	n_surf = aerogrid.aero_timestep_info[time_step].gamma.size();
	for (int i_surf = 0; i_surf < n_surf; i_surf++) {
		for (int i = 0; i < aerogrid.aero_timestep_info[time_step].gamma[i_surf].size(); i++) {
			int N = aerogrid.aero_timestep_info[time_step].gamma[i_surf][i].size();  // Number of columns in the core matrix
			for (int j = 0; j < aerogrid.aero_timestep_info[time_step].gamma[i_surf][i].size(); j++) {
				aerogrid.aero_timestep_info[time_step].gamma[i_surf][i][j] = p_gamma[i_surf][i*N + j];
				aerogrid.aero_timestep_info[time_step].gamma_dot[i_surf][i][j] = p_gamma_dot[i_surf][i*N + j];
			}
		}
	}

	// saving the values of "p_gamma_star"
	n_surf = aerogrid.aero_timestep_info[time_step].gamma_star.size();
	for (int i_surf = 0; i_surf < n_surf; i_surf++) {
		for (int i = 0; i < aerogrid.aero_timestep_info[time_step].gamma_star[i_surf].size(); i++) {
			int N = aerogrid.aero_timestep_info[time_step].gamma_star[i_surf][i].size();  // Number of columns in the core matrix
			for (int j = 0; j < aerogrid.aero_timestep_info[time_step].gamma_star[i_surf][i].size(); j++) {
				aerogrid.aero_timestep_info[time_step].gamma_star[i_surf][i][j] = p_gamma_star[i_surf][i*N + j];
			}
		}
	}

	// saving the values of "p_dist_to_orig"
	n_surf = aerogrid.aero_timestep_info[time_step].dist_to_orig.size();
	for (int i_surf = 0; i_surf < n_surf; i_surf++) {
		for (int i = 0; i < aerogrid.aero_timestep_info[time_step].dist_to_orig[i_surf].size(); i++) {
			int N = aerogrid.aero_timestep_info[time_step].dist_to_orig[i_surf][i].size();  // Number of columns in the core matrix
			for (int j = 0; j < aerogrid.aero_timestep_info[time_step].dist_to_orig[i_surf][i].size(); j++) {
				aerogrid.aero_timestep_info[time_step].dist_to_orig[i_surf][i][j] = p_dist_to_orig[i_surf][i*N + j];
			}
		}
	}

	// saving the values of "p_normals"
	n_surf = aerogrid.aero_timestep_info[time_step].normals.size();
	for (int i_surf = 0; i_surf < n_surf; i_surf++) {
		for (int i = 0; i < aerogrid.aero_timestep_info[time_step].normals[i_surf].size(); i++) {
			for (int j = 0; j < aerogrid.aero_timestep_info[time_step].normals[i_surf][i].size(); j++) {
				int N = aerogrid.aero_timestep_info[time_step].normals[i_surf][i][j].size();  // Number of columns in the core matrix
				for (int k = 0; k < aerogrid.aero_timestep_info[time_step].normals[i_surf][i][j].size(); k++) {
					aerogrid.aero_timestep_info[time_step].normals[i_surf][i][j][k] = p_normals[i + i_surf * n_surf][j*N + k];
				}
			}
		}
	}

	// saving the "p_forces" and the "p_dynamic_forces"
	n_surf = aerogrid.aero_timestep_info[time_step].forces.size();
	for (int i_surf = 0; i_surf < n_surf; i_surf++) {
		for (int i = 0; i < aerogrid.aero_timestep_info[time_step].forces[i_surf].size(); i++) {
			for (int j = 0; j < aerogrid.aero_timestep_info[time_step].forces[i_surf][i].size(); j++) {

				int N = aerogrid.aero_timestep_info[time_step].forces[i_surf][i][j].size();  // Number of columns in the core matrix
				for (int k = 0; k < aerogrid.aero_timestep_info[time_step].forces[i_surf][i][j].size(); k++) {
					aerogrid.aero_timestep_info[time_step].forces[i_surf][i][j][k] = p_forces[i + i_surf * n_surf][j*N + k];
					aerogrid.aero_timestep_info[time_step].dynamic_forces[i_surf][i][j][k] = p_dynamic_forces[i + i_surf * n_surf][j*N + k];
				}
			}
		}
	}
}

UvlmLibVar::~UvlmLibVar() {

	NO_OP;
}



StraightWake::StraightWake() {
	NO_OP;
}

void 
StraightWake::StraightWake_initialize(const StraightWake_settings* straightwakesettings) {


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

void 
StraightWake::StraightWake_generate(AeroTimeStepInfo& aero_tstep) {

	// The object creates a straight wake shedding from the trailing edge based on
	// the time step "dt", the incoming velocity magnitude "u_inf" and
	// direction "u_inf_direction"

	int nsurf = aero_tstep.zeta.size();
	for (int isurf = 0; isurf < nsurf; ++isurf) {
		int M = aero_tstep.zeta_star[isurf][0].size();
		int N = aero_tstep.zeta_star[isurf][0][0].size();
		for (int j = 0; j < N; ++j) {
			for (int k = 0; k < 3; ++k) {
				aero_tstep.zeta_star[isurf][k][0][j] = aero_tstep.zeta[isurf][k][M][j];
			}
			for (int i = 1; i < M; ++i) {
				double deltax = StraightWake_get_deltax(i, _dx1, _ndx1, _r, _dxmax);
				for (int k = 0; k < 3; ++k) {
					aero_tstep.zeta_star[isurf][k][i][j] = aero_tstep.zeta_star[isurf][k][i - 1][j] + deltax * _u_inf_direction[k];
				}
			}
		}
		for (int j = 0; j < aero_tstep.gamma[isurf].size(); ++j) {
			for (int i = 0; i < aero_tstep.gamma[isurf][0].size(); ++i) {
				aero_tstep.gamma[isurf][j][i] *= 0;
				aero_tstep.gamma_star[isurf][j][i] *= 0;
			}
		}
	}

	for (int isurf = 0; isurf < nsurf; ++isurf) {
		int M = aero_tstep.zeta_star[isurf][0].size();
		int N = aero_tstep.zeta_star[isurf][0][0].size();
		for (int j = 0; j < aero_tstep.dist_to_orig[isurf].size(); ++j) {
			for (int i = 0; i < aero_tstep.dist_to_orig[isurf][0].size(); ++i) {
				aero_tstep.dist_to_orig[isurf][j][i] = 0;
			}
		}
		for (int j = 0; j < N; ++j) {
			for (int i = 1; i < M; ++i) {
				double norm = 0;
				for (int k = 0; k < 3; ++k) {
					norm = norm + (aero_tstep.zeta_star[isurf][k][i][j] - aero_tstep.zeta_star[isurf][k][i - 1][j]) *
						(aero_tstep.zeta_star[isurf][k][i][j] - aero_tstep.zeta_star[isurf][k][i - 1][j]);
				}
				norm = sqrt(norm);
				aero_tstep.dist_to_orig[isurf][i][j] = aero_tstep.dist_to_orig[isurf][i - 1][j] + norm;
			}
			for (int i = 0; i < M; ++i) {
				aero_tstep.dist_to_orig[isurf][i][j] /= aero_tstep.dist_to_orig[isurf][M][j];
			}
		}
	}
}

StraightWake::~StraightWake() {
	NO_OP;
}

double 
StraightWake::StraightWake_get_deltax(int i, double dx1, double ndx1, double r, double dxmax) {

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



SteadyVelocityField::SteadyVelocityField() {

	NO_OP;
}

void
SteadyVelocityField::SteadyVelocityField_initialize(const FlightConditions* flight_conditions) {

	u_inf = flight_conditions->uinf;
	u_inf_direction = flight_conditions->uinf_direction;
}

void 
SteadyVelocityField::SteadyVelocityField_generate(AeroTimeStepInfo& aero_tstep, const UVMopts* uvmopts) {

	// The object creates a steady velocity field with the velocity and flow direction specified by the user.

	for (int i_surf = 0; i_surf < aero_tstep.zeta.size(); ++i_surf) {
		for (int k = 0; k < aero_tstep.zeta[i_surf].size(); ++k) {
			for (int i = 0; i < aero_tstep.zeta[i_surf][k].size(); ++i) {
				for (int j = 0; j < aero_tstep.zeta[i_surf][k][i].size(); ++j) {
					aero_tstep.u_ext[i_surf][k][i][j] += u_inf * u_inf_direction[k];
				}
			}
		}
	}

	if (((uvmopts->convection_scheme > 1) && uvmopts->convect_wake) || !uvmopts->cfl1) {
		for (int i_surf = 0; i_surf < aero_tstep.zeta_star.size(); ++i_surf) {
			for (int k = 0; k < aero_tstep.zeta_star[i_surf].size(); ++k) {
				for (int i = 0; i < aero_tstep.zeta_star[i_surf][k].size(); ++i) {
					for (int j = 0; j < aero_tstep.zeta_star[i_surf][k][i].size(); ++j) {
						aero_tstep.u_ext_star[i_surf][k][i][j] += u_inf * u_inf_direction[k];
					}
				}
			}
		}
	}
}

SteadyVelocityField::~SteadyVelocityField() {

	NO_OP;
}




/*
std::vector<double> pMBDyn_UVLM_CouplingKinematic_x;                    //- consistent with the external struct force element
std::vector<double> pMBDyn_UVLM_CouplingKinematic_R;
std::vector<double> pMBDyn_UVLM_CouplingKinematic_xp;
std::vector<double> pMBDyn_UVLM_CouplingKinematic_omega;
std::vector<double> pMBDyn_UVLM_CouplingKinematic_xpp;
std::vector<double> pMBDyn_UVLM_CouplingKinematic_omegap;
std::vector<double> pMBDyn_UVLM_Frame;                                       //- the position [3] and the orietation [9] of UVLM ground coordinate
std::vector<double> pMBDyn_UVLM_CouplingDynamic_f;
std::vector<double> pMBDyn_UVLM_CouplingDynamic_m;
std::vector<double> pMBDyn_UVLM_CouplingDynamic_f_pre;
std::vector<double> pMBDyn_UVLM_CouplingDynamic_m_pre;
*/