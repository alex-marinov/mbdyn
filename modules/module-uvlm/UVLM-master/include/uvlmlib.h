#include "cpp_interface.h"
#include <fenv.h>


DLLEXPORT void run_VLM
(
	const UVLM::Types::VMopts& options,
	const UVLM::Types::FlightConditions& flightconditions,
	unsigned int** p_dimensions,
	unsigned int** p_dimensions_star,
	double** p_zeta,
	double** p_zeta_star,
	double** p_zeta_dot,
	double** p_u_ext,
	double** p_gamma,
	double** p_gamma_star,
	double** p_forces,
	double* p_rbm_vel_g,
	double* p_centre_rot_g
);

DLLEXPORT void init_UVLM
(
	const UVLM::Types::VMopts& options,
	const UVLM::Types::FlightConditions& flightconditions,
	unsigned int** p_dimensions,
	unsigned int** p_dimensions_star,
	double** p_uext,
	double** p_zeta,
	double** p_zeta_star,
	double** p_zeta_dot,
	double** p_zeta_star_dot,
	double*  p_rbm_vel,
	double** p_gamma,
	double** p_gamma_star,
	double** p_normals,
	double** p_forces
);



DLLEXPORT void run_UVLM
(
	const UVLM::Types::UVMopts& options,
	const UVLM::Types::FlightConditions& flightconditions,
	unsigned int** p_dimensions,
	unsigned int** p_dimensions_star,
	unsigned int i_iter,
	double** p_uext,
	double** p_uext_star,
	double** p_zeta,
	double** p_zeta_star,
	double** p_zeta_dot,
	double*  p_rbm_vel,
	double*  p_centre_rot,
	double** p_gamma,
	double** p_gamma_star,
	double** p_dist_to_orig,
	// double** p_previous_gamma,
	double** p_normals,
	double** p_forces,
	double** p_dynamic_forces
);


DLLEXPORT void calculate_unsteady_forces
(
	const UVLM::Types::UVMopts& options,
	const UVLM::Types::FlightConditions& flightconditions,
	unsigned int** p_dimensions,
	unsigned int** p_dimensions_star,
	double** p_zeta,
	double** p_zeta_star,
	double*  p_rbm_vel,
	double** p_gamma,
	double** p_gamma_star,
	double** p_gamma_dot,
	double** p_normals,
	double** p_dynamic_forces
);


DLLEXPORT void UVLM_check_incidence_angle
(
	uint& n_surf,
	unsigned int** p_dimensions,
	double** p_uext,
	double** p_zeta,
	double** p_zeta_dot,
	double** p_normals,
	double*  p_rbm_vel,
	double** p_incidence_angle
);


DLLEXPORT void total_induced_velocity_at_points
(
	const UVLM::Types::UVMopts& options,
	unsigned int** p_dimensions,
	unsigned int** p_dimensions_star,
	double** p_zeta,
	double** p_zeta_star,
	double** p_gamma,
	double** p_gamma_star,
	double* p_target_triads,
	double* p_uout,
	unsigned int npoints
);


DLLEXPORT void multisurface
(
	const UVLM::Types::UVMopts& options,
	unsigned int** p_dimensions,
	unsigned int** p_dimensions_target,
	unsigned int** p_dimensions_uout,
	double** p_zeta,
	double** p_gamma,
	double** p_target_surface,
	double** p_uout
	// double  vortex_radius
);



// linear UVLM interface

DLLEXPORT void call_der_biot_panel(double p_DerP[9],
	double p_DerVertices[36],// 4x9
	double p_zetaP[3],
	double p_ZetaPanel[12],
	const double& gamma,
	double& vortex_radius);



DLLEXPORT void call_biot_panel(double p_vel[3],
	double p_zetaP[3],
	double p_ZetaPanel[12],
	const double& gamma,
	double& vortex_radius);



DLLEXPORT void call_dvinddzeta(double p_DerC[9],
	double p_DerV[],
	double p_zetaC[3],
	double p_ZetaIn[],
	double p_GammaIn[],
	int& M_in,
	int& N_in,
	bool& IsBound,
	int& M_in_bound, // M of bound surf associated
	double& vortex_radius);



DLLEXPORT void call_aic3(double p_AIC3[],
	double p_zetaC[3],
	double p_ZetaIn[],
	int& M_in,
	int& N_in,
	double& vortex_radius);



DLLEXPORT void call_ind_vel(
	double p_vel[3],
	double p_zetaC[3],
	double p_ZetaIn[],
	double p_GammaIn[],
	int& M_in,
	int& N_in,
	double& vortex_radius);

