# MBDyn (C) is a multibody analysis code. 
# $Header$
# http://www.mbdyn.org
# 
# Copyright (C) 1996-2015
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

# Skip undesired labels
function skip_label(l)
{
	return 0;
}

# Condition
function condition(d) {
	# return 1 if the condition of periodicity is met;
	# plus, put in d[1] and d[2] the weights of the previous
	# and the current step for linear interpolation of the
	# periodicity time, and in d[3] the label of a node
	# that is considered as reference

	# this condition applies to the ADYN model of rotor + pylon:
	# when node 5000 (the shaft) completes a round about its
	# axis (approximate, in the global frame) a period is
	# completed.
	# (Note: the node is oriented with axis 3 in direction -1,
	# axis 1 in direction 3 and axis 2 in direction 2 when the
	# model is assembled; the condition is represented by axis
	# 1 of the node crossing the global plane 1-2 coming from 
	# below
	euler2R(The[5000, 1], The[5000, 2], The[5000, 3], Rtmp);
	euler2R(OldThe[5000, 1], OldThe[5000, 2], OldThe[5000, 3], OldRtmp);
	if (Rtmp[2, 1] >= 0 && OldRtmp[2, 1] < 0) {
		dd = Rtmp[2, 1] - OldRtmp[2, 1];
		d[1] = Rtmp[2, 1]/dd;
		d[2] = - OldRtmp[2, 1]/dd;
		d[3] = 5000;
		return 1;
	}

	# default: no periodicity
	return 0;
}

# Output
function output() {
	if (condition(d)) {
		for (node in label) {

			if (node == d[3]) {
	
				printf("%8d %e %e %e %e %e %e %e %e %e %e %e %e        %e %e %e %d\n",
					node,
					d[2]*Pos[node, 1] + d[1]*OldPos[node, 1],
					d[2]*Pos[node, 2] + d[1]*OldPos[node, 2],
					d[2]*Pos[node, 3] + d[1]*OldPos[node, 3],
					d[2]*The[node, 1] + d[1]*OldThe[node, 1],
					d[2]*The[node, 2] + d[1]*OldThe[node, 2],
					d[2]*The[node, 3] + d[1]*OldThe[node, 3],
					d[2]*Vel[node, 1] + d[1]*OldVel[node, 1],
					d[2]*Vel[node, 2] + d[1]*OldVel[node, 2],
					d[2]*Vel[node, 3] + d[1]*OldVel[node, 3],
					d[2]*Ome[node, 1] + d[1]*OldOme[node, 1],
					d[2]*Ome[node, 2] + d[1]*OldOme[node, 2],
					d[2]*Ome[node, 3] + d[1]*OldOme[node, 3],
					d[1], d[2], d[1] + d[2], step);
			} else {

				printf("%8d %e %e %e %e %e %e %e %e %e %e %e %e\n",
					node,
					d[2]*Pos[node, 1] + d[1]*OldPos[node, 1],
					d[2]*Pos[node, 2] + d[1]*OldPos[node, 2],
					d[2]*Pos[node, 3] + d[1]*OldPos[node, 3],
					d[2]*The[node, 1] + d[1]*OldThe[node, 1],
					d[2]*The[node, 2] + d[1]*OldThe[node, 2],
					d[2]*The[node, 3] + d[1]*OldThe[node, 3],
					d[2]*Vel[node, 1] + d[1]*OldVel[node, 1],
					d[2]*Vel[node, 2] + d[1]*OldVel[node, 2],
					d[2]*Vel[node, 3] + d[1]*OldVel[node, 3],
					d[2]*Ome[node, 1] + d[1]*OldOme[node, 1],
					d[2]*Ome[node, 2] + d[1]*OldOme[node, 2],
					d[2]*Ome[node, 3] + d[1]*OldOme[node, 3]);
			}
		}
	}

	for (node in label) {
		OldPos[node, 1] = Pos[node, 1];
		OldPos[node, 2] = Pos[node, 2];
		OldPos[node, 3] = Pos[node, 3];

		OldThe[node, 1] = The[node, 1];
		OldThe[node, 2] = The[node, 2];
		OldThe[node, 3] = The[node, 3];

		OldVel[node, 1] = Vel[node, 1];
		OldVel[node, 2] = Vel[node, 2];
		OldVel[node, 3] = Vel[node, 3];

		OldOme[node, 1] = Ome[node, 1];
		OldOme[node, 2] = Ome[node, 2];
		OldOme[node, 3] = Ome[node, 3];
	}
}

# Initialize vars
BEGIN {
	FirstLabel = -1;
	FirstStep = 1;
	# SkipSteps = 0;
	step = 0;

	deg2rad = 0.017453293;
	rad2deg = 57.29578;
	AngleScale = deg2rad;
}

{
	# every time a new step starts, output the previous one
	if ($1 == FirstLabel) {
		if (step >= SkipSteps) {
			output();
		}
		step++;
	}

	# get the first label, which marks the beginning of a new time step
	if (FirstStep) {
		FirstLabel = $1;
		FirstStep = 0;
	}

	if (step < SkipSteps) {
		next;
	}

	# skip undesired labels (modify the body of the function at will)
	if (skip_label($1)) {
		next;
	}

	# store the configuration of the current node
	label[$1] = $1;

	Pos[$1, 1] = $2;
	Pos[$1, 2] = $3;
	Pos[$1, 3] = $4;

	The[$1, 1] = $5;
	The[$1, 2] = $6;
	The[$1, 3] = $7;

	Vel[$1, 1] = $8;
	Vel[$1, 2] = $9;
	Vel[$1, 3] = $10;

	Ome[$1, 1] = $11;
	Ome[$1, 2] = $12;
	Ome[$1, 3] = $13;
}

END {
	output();
}
