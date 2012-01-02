# $Header$
# MBDyn (C) is a multibody analysis code. 
# http://www.mbdyn.org
# 
# Copyright (C) 1996-2012
# 
# Pierangelo Masarati	<masarati@aero.polimi.it>
# Paolo Mantegazza	<mantegazza@aero.polimi.it>
# 
# Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
# via La Masa, 34 - 20156 Milano, Italy
# http://www.aero.polimi.it
# 
# Changing this copyright notice is forbidden.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation (version 2 of the License).
# 
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# Author:	Pierangelo Masarati	<masarati@aero.polimi.it>
#
# Use:		awk -f abs2rel.awk -v<option>=<value> [...] <file>.mov
#
# options:
#		option name	expected value
#		RefNode		label of reference node (required)
#		InputMode	{euler123|vector|matrix} (default: euler123)
#		OutputMode	{euler123|vector|matrix} (default: euler123)
#		RefOnly		{0|1} (default: 0): re-position & re-orient only
#
# undocumented options:
#		NumBlades	(int) <nb>
#		BladeLabelOff	(int) <off>
#		RotorAxis	{1,2,3}|{x,y,z}

# R matrix from Euler 123 angles
function euler2R(Alpha, Beta, Gamma, R,   dCosAlpha, dSinAlpha, dCosBeta, dSinBeta, dCosGamma, dSinGamma) {
	dCosAlpha = cos(Alpha*AngleScale);
	dSinAlpha = sin(Alpha*AngleScale);
	dCosBeta = cos(Beta*AngleScale);
	dSinBeta = sin(Beta*AngleScale);
	dCosGamma = cos(Gamma*AngleScale);
	dSinGamma = sin(Gamma*AngleScale);

	R[1, 1] = dCosBeta*dCosGamma;
	R[2, 1] = dCosAlpha*dSinGamma+dSinAlpha*dSinBeta*dCosGamma;
	R[3, 1] = dSinAlpha*dSinGamma-dCosAlpha*dSinBeta*dCosGamma;

	R[1, 2] = -dCosBeta*dSinGamma;
	R[2, 2] = dCosAlpha*dCosGamma-dSinAlpha*dSinBeta*dSinGamma;
	R[3, 2] = dSinAlpha*dCosGamma+dCosAlpha*dSinBeta*dSinGamma;

	R[1, 3] = dSinBeta;
	R[2, 3] = -dSinAlpha*dCosBeta;
	R[3, 3] = dCosAlpha*dCosBeta;
}

# R matrix from orientation vector
function vec2R(Phi1, Phi2, Phi3, R,   dPhi, dCosPhi, dSinPhi)
{
	dPhi = sqrt(Phi1*Phi1 + Phi2*Phi2 + Phi3*Phi3);

	dCosPhi = cos(dPhi);
	dSinPhi = sin(dPhi);

	if (dPhi > 1.e-15) {
		Phi1 /= dPhi;
		Phi2 /= dPhi;
		Phi3 /= dPhi;

		R[1, 1] = 1. - (1. - dCosPhi)*(Phi2*Phi2 + Phi3*Phi3);
		R[1, 2] = -dSinPhi*Phi3 + (1. - dCosPhi)*Phi1*Phi2;
		R[1, 3] = dSinPhi*Phi2 + (1. - dCosPhi)*Phi1*Phi3;

		R[2, 1] = dSinPhi*Phi3 + (1. - dCosPhi)*Phi2*Phi1;
		R[2, 2] = 1. - (1. - dCosPhi)*(Phi3*Phi3 + Phi1*Phi1);
		R[2, 3] = -dSinPhi*Phi1 + (1. - dCosPhi)*Phi2*Phi3;

		R[3, 1] = -dSinPhi*Phi2 + (1. - dCosPhi)*Phi3*Phi1;
		R[3, 2] = dSinPhi*Phi1 + (1. - dCosPhi)*Phi3*Phi2;
		R[3, 3] = 1. - (1. - dCosPhi)*(Phi1*Phi1 + Phi2*Phi2);

	} else {
		R[1, 1] = 1.;
		R[1, 2] = 0.;
		R[1, 3] = 0.;

		R[2, 1] = 0.;
		R[2, 2] = 1.;
		R[2, 3] = 0.;

		R[3, 1] = 0.;
		R[3, 2] = 0.;
		R[3, 3] = 1.;
	}
}

# Matrix R to Euler 123 params
function R2euler(R, e,   dAlpha, dCosAlpha, dSinAlpha)
{
	dAlpha = atan2(-R[2, 3], R[3, 3]);
	dCosAlpha = cos(dAlpha);
	dSinAlpha = sin(dAlpha);

	e[1] = dAlpha/AngleScale;
	e[2] = atan2(R[1, 3], dCosAlpha*R[3, 3] - dSinAlpha*R[2, 3])/AngleScale;
	e[3] = atan2(dCosAlpha*R[2, 1]+dSinAlpha*R[3, 1],
		dCosAlpha*R[2, 2]+dSinAlpha*R[3, 2])/AngleScale;
}

# Matrix R to orientation vector
function R2vec(R, Phi,   dPhi, dCosPhi, dSinPhi)
{
	dCosPhi = (R[1, 1] + R[2, 2] + R[3, 3] - 1.)/2.;

	Phi[1] = (R[3, 2] - R[2, 3])/2.;
	Phi[2] = (R[1, 3] - R[3, 1])/2.;
	Phi[3] = (R[2, 1] - R[1, 2])/2.;

	dSinPhi = sqrt(Phi[1]*Phi[1] + Phi[2]*Phi[2] + Phi[3]*Phi[3]);

	dPhi = atan2(dSinPhi, dCosPhi);

	# above 1e-7 they start differ in double precision
	if (dSinPhi > 1.e-7) {
		Phi[1] *= dPhi/dSinPhi;
		Phi[2] *= dPhi/dSinPhi;
		Phi[3] *= dPhi/dSinPhi;
	}
}

# ...
function mat3_mul_vec3(m, v, r)
{
	r[1] = m[1, 1]*v[1] + m[1, 2]*v[2] + m[1, 3]*v[3];
	r[2] = m[2, 1]*v[1] + m[2, 2]*v[2] + m[2, 3]*v[3];
	r[3] = m[3, 1]*v[1] + m[3, 2]*v[2] + m[3, 3]*v[3];
}

# ...
function mat3T_mul_vec3(m, v, r)
{
	r[1] = m[1, 1]*v[1] + m[2, 1]*v[2] + m[3, 1]*v[3];
	r[2] = m[1, 2]*v[1] + m[2, 2]*v[2] + m[3, 2]*v[3];
	r[3] = m[1, 3]*v[1] + m[2, 3]*v[2] + m[3, 3]*v[3];
}

# ...
function mat3_mul_mat3(m1, m2, r)
{
	r[1, 1] = m1[1, 1]*m2[1, 1] + m1[1, 2]*m2[2, 1] + m1[1, 3]*m2[3, 1];
	r[2, 1] = m1[2, 1]*m2[1, 1] + m1[2, 2]*m2[2, 1] + m1[2, 3]*m2[3, 1];
	r[3, 1] = m1[3, 1]*m2[1, 1] + m1[3, 2]*m2[2, 1] + m1[3, 3]*m2[3, 1];

	r[1, 2] = m1[1, 1]*m2[1, 2] + m1[1, 2]*m2[2, 2] + m1[1, 3]*m2[3, 2];
	r[2, 2] = m1[2, 1]*m2[1, 2] + m1[2, 2]*m2[2, 2] + m1[2, 3]*m2[3, 2];
	r[3, 2] = m1[3, 1]*m2[1, 2] + m1[3, 2]*m2[2, 2] + m1[3, 3]*m2[3, 2];

	r[1, 3] = m1[1, 1]*m2[1, 3] + m1[1, 2]*m2[2, 3] + m1[1, 3]*m2[3, 3];
	r[2, 3] = m1[2, 1]*m2[1, 3] + m1[2, 2]*m2[2, 3] + m1[2, 3]*m2[3, 3];
	r[3, 3] = m1[3, 1]*m2[1, 3] + m1[3, 2]*m2[2, 3] + m1[3, 3]*m2[3, 3];
}

# ...
function mat3T_mul_mat3(m1, m2, r)
{
	r[1, 1] = m1[1, 1]*m2[1, 1] + m1[2, 1]*m2[2, 1] + m1[3, 1]*m2[3, 1];
	r[2, 1] = m1[1, 2]*m2[1, 1] + m1[2, 2]*m2[2, 1] + m1[3, 2]*m2[3, 1];
	r[3, 1] = m1[1, 3]*m2[1, 1] + m1[2, 3]*m2[2, 1] + m1[3, 3]*m2[3, 1];

	r[1, 2] = m1[1, 1]*m2[1, 2] + m1[2, 1]*m2[2, 2] + m1[3, 1]*m2[3, 2];
	r[2, 2] = m1[1, 2]*m2[1, 2] + m1[2, 2]*m2[2, 2] + m1[3, 2]*m2[3, 2];
	r[3, 2] = m1[1, 3]*m2[1, 2] + m1[2, 3]*m2[2, 2] + m1[3, 3]*m2[3, 2];

	r[1, 3] = m1[1, 1]*m2[1, 3] + m1[2, 1]*m2[2, 3] + m1[3, 1]*m2[3, 3];
	r[2, 3] = m1[1, 2]*m2[1, 3] + m1[2, 2]*m2[2, 3] + m1[3, 2]*m2[3, 3];
	r[3, 3] = m1[1, 3]*m2[1, 3] + m1[2, 3]*m2[2, 3] + m1[3, 3]*m2[3, 3];
}

# ...
function vec3_cross_vec3(v1, v2, r)
{
	r[1] = v1[2]*v2[3] - v1[3]*v2[2];
	r[2] = v1[3]*v2[1] - v1[1]*v2[3];
	r[3] = v1[1]*v2[2] - v1[2]*v2[1];
}

# Modify reference configuration
function prepare(X, R, V, W, acc, A, WP)
{
	# add user-defined modifications to reference configuration
	# operate on A, WP only if acc != 0
	return;

	# 45 deg rotation about z axis ...
	euler2R(0., 0., 45., Rtmp);

	RR[1, 1] = R[1, 1];
	RR[2, 1] = R[2, 1];
	RR[3, 1] = R[3, 1];

	RR[1, 2] = R[1, 2];
	RR[2, 2] = R[2, 2];
	RR[3, 2] = R[3, 2];

	RR[1, 3] = R[1, 3];
	RR[2, 3] = R[2, 3];
	RR[3, 3] = R[3, 3];

	mat3_mul_mat3(RR, Rtmp, R);
}

# Skip undesired labels
function skip_label(l)
{
	return 0;
}

# Output
function output() {
	for (node = 0; node < num_labels; node++) {
		if (label[node] == RefNode) {
			printf("%8d %13.6e %13.6e %13.6e",
				label[node],
				0., 0., 0.);
			if (omode == 2) {
				printf(" %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e ",
					0., 0., 0.,
					0., 0., 0.,
					0., 0., 0.);
			} else {
				printf(" %13.6e %13.6e %13.6e",
					0., 0., 0.);
			}
			if (RefOnly) {
				# velocity
				Vr[1] = Vel[node, 1];
				Vr[2] = Vel[node, 2];
				Vr[3] = Vel[node, 3];
	
				mat3T_mul_vec3(R0, Vr, V);
	
				# angular velocity
				Wr[1] = Ome[node, 1];
				Wr[2] = Ome[node, 2];
				Wr[3] = Ome[node, 3];
	
				mat3T_mul_vec3(R0, Wr, W);
	
				if (accel[node]) {
					# acceleration
					Ar[1] = Acc[node, 1];
					Ar[2] = Acc[node, 2];
					Ar[3] = Acc[node, 3];
	
					mat3T_mul_vec3(R0, Ar, A);
	
					# angular acceleration
					WPr[1] = OmeP[node, 1];
					WPr[2] = OmeP[node, 2];
					WPr[3] = OmeP[node, 3];
	
					mat3T_mul_vec3(R0, WPr, WP);
				}

			} else {
				printf(" %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e ",
					V[1], V[2], V[3],
					W[1], W[2], W[3]);
				if (acc) {
					printf(" %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e",
						A[1], A[2], A[3],
						WP[1], WP[2], WP[3]);
				}
			}
			printf("\n");
			continue;
		}

		if (NumBlades > 0) {
			# Note: blades must be defined with local axis 1
			# along the blade axis, axis 3 along the flapwise
			# direction, axis 2 according to a right triad
			# (chordwise trail to lead for counter-clockwise
			# rotors)
			# Blade #1 is assumed to be aligned with the reference
			# node; higher blade numbers lag behind
			#
			#         ^            blade #4  |
			#         | y                    |
			#         |                      |
			#         |     x                |  blade #1
			#	  +------>      ---------o---------
			#                      blade #3  |
			#                                |
			#                                |
			#                                | blade #2
			#
			ib = int(label[node]/BladeLabelOff);
			if (ib > 0 && ib <= NumBlades) {
				psi = -2*3.14159265358979323846*(ib - 1)/NumBlades;

				theta[1] = 0.;
				theta[2] = 0.;
				theta[3] = 0.;
				theta[RotorAxis] = psi;
				vec2R(theta[1], theta[2], theta[3], Rpsi);
				mat3_mul_mat3(R0, Rpsi, Rb)
			} else {
				Rb[1, 1] = R0[1, 1];
				Rb[1, 2] = R0[1, 2];
				Rb[1, 3] = R0[1, 3];

				Rb[2, 1] = R0[2, 1];
				Rb[2, 2] = R0[2, 2];
				Rb[2, 3] = R0[2, 3];

				Rb[3, 1] = R0[3, 1];
				Rb[3, 2] = R0[3, 2];
				Rb[3, 3] = R0[3, 3];
			}
		}

		# position
		Xr[1] = Pos[node, 1] - X0[1];
		Xr[2] = Pos[node, 2] - X0[2];
		Xr[3] = Pos[node, 3] - X0[3];

		if (NumBlades > 0) {
			mat3T_mul_vec3(Rb, Xr, X);
		} else {
			mat3T_mul_vec3(R0, Xr, X);
		}

		# orientation
		if (imode == 0) {
			euler2R(The[node, 1], The[node, 2], The[node, 3], Rn);

		} else if (imode == 1) {
			vec2R(The[node, 1], The[node, 2], The[node, 3], Rn);

		} else if (imode == 2) {
			Rn[1, 1] = The[node, 1];
			Rn[1, 2] = The[node, 2];
			Rn[1, 3] = The[node, 3];

			Rn[2, 1] = The[node, 4];
			Rn[2, 2] = The[node, 5];
			Rn[2, 3] = The[node, 6];

			Rn[3, 1] = The[node, 7];
			Rn[3, 2] = The[node, 8];
			Rn[3, 3] = The[node, 9];
		}

		if (NumBlades > 0) {
			mat3T_mul_mat3(Rb, Rn, R);
		} else {
			mat3T_mul_mat3(R0, Rn, R);
		}

		if (omode == 0) {
			R2euler(R, e);

		} else if (omode == 1) {
			R2vec(R, e);
		}

		if (RefOnly) {
			# velocity
			Vr[1] = Vel[node, 1];
			Vr[2] = Vel[node, 2];
			Vr[3] = Vel[node, 3];

			if (NumBlades > 0) {
				mat3T_mul_vec3(Rb, Vr, V);
			} else {
				mat3T_mul_vec3(R0, Vr, V);
			}

			# angular velocity
			Wr[1] = Ome[node, 1];
			Wr[2] = Ome[node, 2];
			Wr[3] = Ome[node, 3];

			if (NumBlades > 0) {
				mat3T_mul_vec3(Rb, Wr, W);
			} else {
				mat3T_mul_vec3(R0, Wr, W);
			}

			if (accel[node]) {
				# acceleration
				Ar[1] = Acc[node, 1];
				Ar[2] = Acc[node, 2];
				Ar[3] = Acc[node, 3];

				if (NumBlades > 0) {
					mat3T_mul_vec3(Rb, Ar, A);
				} else {
					mat3T_mul_vec3(R0, Ar, A);
				}

				# angular acceleration
				WPr[1] = OmeP[node, 1];
				WPr[2] = OmeP[node, 2];
				WPr[3] = OmeP[node, 3];

				if (NumBlades > 0) {
					mat3T_mul_vec3(Rb, WPr, WP);
				} else {
					mat3T_mul_vec3(R0, WPr, WP);
				}
			}

		} else {
			# velocity
			vec3_cross_vec3(W0, Xr, Vtmp);

			Vr[1] = Vel[node, 1] - V0[1] - Vtmp[1];
			Vr[2] = Vel[node, 2] - V0[2] - Vtmp[2];
			Vr[3] = Vel[node, 3] - V0[3] - Vtmp[3];

			if (NumBlades > 0) {
				mat3T_mul_vec3(Rb, Vr, V);
			} else {
				mat3T_mul_vec3(R0, Vr, V);
			}

			# angular velocity
			Wr[1] = Ome[node, 1] - W0[1];
			Wr[2] = Ome[node, 2] - W0[2];
			Wr[3] = Ome[node, 3] - W0[3];

			if (NumBlades > 0) {
				mat3T_mul_vec3(Rb, Wr, W);
			} else {
				mat3T_mul_vec3(R0, Wr, W);
			}

			if (accel[node]) {
				# acceleration
				vec3_cross_vec3(WP0, Xr, ATmp);

				Ar[1] = Acc[node, 1] - A0[1] - ATmp[1];
				Ar[2] = Acc[node, 2] - A0[2] - ATmp[2];
				Ar[3] = Acc[node, 3] - A0[3] - ATmp[3];

				vec3_cross_vec3(W0, Xr, ATmp1);
				vec3_cross_vec3(W0, ATmp1, ATmp);

				Ar[1] -= ATmp[1]
				Ar[2] -= ATmp[2]
				Ar[3] -= ATmp[3]

				vec3_cross_vec3(W0, Vr, ATmp);

				Ar[1] -= 2*ATmp[1]
				Ar[2] -= 2*ATmp[2]
				Ar[3] -= 2*ATmp[3]
			
				if (NumBlades > 0) {
					mat3T_mul_vec3(Rb, Ar, A);
				} else {
					mat3T_mul_vec3(R0, Ar, A);
				}

				# angular acceleration
				vec3_cross_vec3(W0, Wr, WPTmp);

				WPr[1] = OmeP[node, 1] - WP0[1] - WPTmp[1];
				WPr[2] = OmeP[node, 2] - WP0[2] - WPTmp[2];
				WPr[3] = OmeP[node, 3] - WP0[3] - WPTmp[3];

				if (NumBlades > 0) {
					mat3T_mul_vec3(Rb, WPr, WP);
				} else {
					mat3T_mul_vec3(R0, WPr, WP);
				}
			}
		}

		printf("%8d %13.6e %13.6e %13.6e", label[node], X[1], X[2], X[3]);

		if (omode == 2) {
			printf(" %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e",
				R[1, 1], R[1, 2], R[1, 3],
				R[2, 1], R[2, 2], R[2, 3],
				R[3, 1], R[3, 2], R[3, 3]);
		} else {
			printf(" %13.6e %13.6e %13.6e", e[1], e[2], e[3]);
		}
	
		printf(" %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e",
			V[1], V[2], V[3], W[1], W[2], W[3]);
	
		if (accel[node]) {
			printf(" %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e",
				A[1], A[2], A[3], WP[1], WP[2], WP[3]);
		}
	
		printf("\n");
	}
}

# Initialize vars
BEGIN {
	FirstLabel = -1;
	FirstStep = 1;
	# SkipSteps = 0;
	step = 0;
	node = 0;

	if (InputMode == "euler123") {
		imode = 0;

	} else if (InputMode == "vector") {
		imode = 1;

	} else if (InputMode == "matrix") {
		imode = 2;

	} else if (InputMode != 0) {
		printf("Unknown InputMode=%s\n", InputMode);
		exit;
	}

	if (OutputMode == "euler123") {
		omode = 0;

	} else if (OutputMode == "vector") {
		omode = 1;

	} else if (OutputMode == "matrix") {
		omode = 2;

	} else if (OutputMode != 0) {
		printf("Unknown OutputMode=%s\n", OutputMode);
		exit;
	}

	if (NumBlades < 0) {
		printf("NumBlades must be 0 (not a rotor) or positive\n");
		exit;
	}

	if (NumBlades > 0) {
		if (BladeLabelOff == 0) {
			printf("Need BladeLabelOffset\n");
			exit;

		} else if (int(BladeLabelOff) != BladeLabelOff) {
			printf("BladeLabelOffset must be an integer\n");
			exit;

		} else if (BladeLabelOff <= 0) {
			printf("BladeLabelOffset must be positive\n");
			exit;
		}

		if (RotorAxis == "x") {
			RotorAxis = 0 + 1;
		} else if (RotorAxis == "y") {
			RotorAxis = 0 + 2;
		} else if (RotorAxis == "z") {
			RotorAxis = 0 + 3;
		} else if (RotorAxis != 1 && RotorAxis != 2 && RotorAxis != 3) {
			printf("Invalid RotorAxis\n");
			exit;
		}
	}

	deg2rad = 0.017453293;
	rad2deg = 57.29578;
	AngleScale = deg2rad;
}

{
	# every time a new step starts, output the previous one
	if ($1 == FirstLabel) {
		if (num_labels == 0) {
			num_labels = node;
		}

		if (step >= SkipSteps) {
			output();
		}
		step++;

		node = 0;
	}

	# get the first label, which marks the beginning of a new time step
	if (FirstStep) {
		FirstLabel = $1;
		FirstStep = 0;
	}

	if (step < SkipSteps) {
		next;
	}

	# specially store the reference node configuration
	if ($1 == RefNode) {
		X0[1] = $2;
		X0[2] = $3;
		X0[3] = $4;
		if (imode == 0) {
			euler2R($5, $6, $7, R0);
			offset = 7;

		} else if (imode == 1) {
			vec2R($5, $6, $7, R0);
			offset = 7;

		} else if (imode == 2) {
			R0[1, 1] = $5;
			R0[1, 2] = $6;
			R0[1, 3] = $7;

			R0[2, 1] = $8;
			R0[2, 2] = $9;
			R0[2, 3] = $10;

			R0[3, 1] = $11;
			R0[3, 2] = $12;
			R0[3, 3] = $13;

			offset = 13;
		}
		V0[1] = $(offset + 1);
		V0[2] = $(offset + 2);
		V0[3] = $(offset + 3);

		W0[1] = $(offset + 4);
		W0[2] = $(offset + 5);
		W0[3] = $(offset + 6);

		acc = 0;

		if (NF > offset + 7) {
			A0[1] = $(offset + 7);
			A0[2] = $(offset + 8);
			A0[3] = $(offset + 9);

			WP0[1] = $(offset + 10);
			WP0[2] = $(offset + 11);
			WP0[3] = $(offset + 12);

			acc = 1;
		}

		prepare(X0, R0, V0, W0, acc, A0, WP0);
	}

	# skip undesired labels (modify the body of the function at will)
	if (skip_label($1)) {
		next;
	}

	# store the configuration of the current node
	label[node] = $1;

	Pos[node, 1] = $2;
	Pos[node, 2] = $3;
	Pos[node, 3] = $4;

	if (imode == 0 || imode == 1) {
		The[node, 1] = $5;
		The[node, 2] = $6;
		The[node, 3] = $7;

		offset = 7;

	} else if (imode == 2) {
		The[node, 1] = $5;
		The[node, 2] = $6;
		The[node, 3] = $7;

		The[node, 4] = $8;
		The[node, 5] = $9;
		The[node, 6] = $10;

		The[node, 7] = $11;
		The[node, 8] = $12;
		The[node, 9] = $13;

		offset = 13;
	}

	Vel[node, 1] = $(offset + 1);
	Vel[node, 2] = $(offset + 2);
	Vel[node, 3] = $(offset + 3);

	Ome[node, 1] = $(offset + 4);
	Ome[node, 2] = $(offset + 5);
	Ome[node, 3] = $(offset + 6);

	if (NF > offset + 7) {
		accel[node] = 1;

		Acc[node, 1] = $(offset + 7);
		Acc[node, 2] = $(offset + 8);
		Acc[node, 3] = $(offset + 9);

		OmeP[node, 1] = $(offset + 10);
		OmeP[node, 2] = $(offset + 11);
		OmeP[node, 3] = $(offset + 12);
	}

	node++;
}

END {
	if (num_labels == 0) {
		num_labels = node;
	}

	output();
}
