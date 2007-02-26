# $Header$
# MBDyn (C) is a multibody analysis code. 
# http://www.mbdyn.org
# 
# Copyright (C) 1996-2007
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
# R matrix
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

# Euler params
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
function prepare(X, R, V, W)
{
	# add user-defined modifications to reference configuration
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
			printf("%8d %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n",
				label[node],
				0., 0., 0., 0., 0., 0.,
				0., 0., 0., 0., 0., 0.);
			continue;
		}

		# position
		Xr[1] = Pos[node, 1] - X0[1];
		Xr[2] = Pos[node, 2] - X0[2];
		Xr[3] = Pos[node, 3] - X0[3];

		mat3T_mul_vec3(R0, Xr, X);

		# orientation
		euler2R(The[node, 1], The[node, 2], The[node, 3], Rn);

		mat3T_mul_mat3(R0, Rn, R);

		R2euler(R, e);

		# Velocity
		Wn[1] = Ome[node, 1];
		Wn[2] = Ome[node, 2];
		Wn[3] = Ome[node, 3];

		vec3_cross_vec3(Wn, Xr, Vtmp);

		Vtmp[1] = Vel[node, 1] - V0[1] - Vtmp[1];
		Vtmp[2] = Vel[node, 2] - V0[2] - Vtmp[2];
		Vtmp[3] = Vel[node, 3] - V0[3] - Vtmp[3];

		mat3T_mul_vec3(R0, Vtmp, V);

		# angular velocity
		Wtmp[1] = W0[1] - Ome[node, 1];
		Wtmp[2] = W0[2] - Ome[node, 2];
		Wtmp[3] = W0[3] - Ome[node, 3];

		mat3T_mul_vec3(R0, Wtmp, W);

		printf("%8d %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n",
			label[node],
			X[1], X[2], X[3], e[1], e[2], e[3],
			V[1], V[2], V[3], W[1], W[2], W[3]);
	}
}

# Initialize vars
BEGIN {
	FirstLabel = -1;
	FirstStep = 1;
	# SkipSteps = 0;
	step = 0;
	node = 0;

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
		euler2R($5, $6, $7, R0);
		V0[1] = $8;
		V0[2] = $9;
		V0[3] = $10;
		W0[1] = $11;
		W0[2] = $12;
		W0[3] = $13;

		prepare(X0, R0, V0, W0);
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

	The[node, 1] = $5;
	The[node, 2] = $6;
	The[node, 3] = $7;

	Vel[node, 1] = $8;
	Vel[node, 2] = $9;
	Vel[node, 3] = $10;

	Ome[node, 1] = $11;
	Ome[node, 2] = $12;
	Ome[node, 3] = $13;

	node++;
}

END {
	output();
}
