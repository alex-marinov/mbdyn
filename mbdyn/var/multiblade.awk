# $Header$
# MBDyn (C) is a multibody analysis code. 
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
# Use:		awk -f multiblabe.awk -v<option>=<value> [...] <file>.mov
#
# options:
#		option name	expected value
#		NumBlades	number of blades (required)
#		BladeSeq	sequence of blades ([counter]clockwise)
#		HubNode		label of reference rotating node (required)
#		ShaftNode	label of reference non-rotating node (required)
#		RotorAxis	{x|y|z} (required)
#		BladeLabelOff	offset of blade labels (required)
#		BladeLabelDelta	delta of blade labels (default: BladeLabelOff)
#		PsiOff		azimuth offset (default: 0)
#		InputMode	{euler123|vector|matrix} (default: euler123)
#
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
	for (node2 = 0; node2 < num_labels/NumBlades; node2++) {
		X[node2, 0, 1] = 0;
		X[node2, 0, 2] = 0;
		X[node2, 0, 3] = 0;

		T[node2, 0, 1] = 0;
		T[node2, 0, 2] = 0;
		T[node2, 0, 3] = 0;

		V[node2, 0, 1] = 0;
		V[node2, 0, 2] = 0;
		V[node2, 0, 3] = 0;

		W[node2, 0, 1] = 0;
		W[node2, 0, 2] = 0;
		W[node2, 0, 3] = 0;

		for (n = 1; n <= NumCyc; n++) {
			X[node2, n "c", 1] = 0;
			X[node2, n "c", 2] = 0;
			X[node2, n "c", 3] = 0;

			T[node2, n "c", 1] = 0;
			T[node2, n "c", 2] = 0;
			T[node2, n "c", 3] = 0;

			V[node2, n "c", 1] = 0;
			V[node2, n "c", 2] = 0;
			V[node2, n "c", 3] = 0;

			W[node2, n "c", 1] = 0;
			W[node2, n "c", 2] = 0;
			W[node2, n "c", 3] = 0;

			X[node2, n "s", 1] = 0;
			X[node2, n "s", 2] = 0;
			X[node2, n "s", 3] = 0;

			T[node2, n "s", 1] = 0;
			T[node2, n "s", 2] = 0;
			T[node2, n "s", 3] = 0;

			V[node2, n "s", 1] = 0;
			V[node2, n "s", 2] = 0;
			V[node2, n "s", 3] = 0;

			W[node2, n "s", 1] = 0;
			W[node2, n "s", 2] = 0;
			W[node2, n "s", 3] = 0;
		}

		if (Even != 0) {
			X[node2, Even, 1] = 0;
			X[node2, Even, 2] = 0;
			X[node2, Even, 3] = 0;

			T[node2, Even, 1] = 0;
			T[node2, Even, 2] = 0;
			T[node2, Even, 3] = 0;

			V[node2, Even, 1] = 0;
			V[node2, Even, 2] = 0;
			V[node2, Even, 3] = 0;

			W[node2, Even, 1] = 0;
			W[node2, Even, 2] = 0;
			W[node2, Even, 3] = 0;
		}
	}

	mat3T_mul_mat3(Rs, R0, DR);
	R2vec(DR, Psi);
	psi = Psi[RotorAxis];

	CurrBladeLabelMax = BladeLabelOff;
	psi_m = psi
	m1m = 1;

	blade = 0;
	theta[1] = 0.;
	theta[2] = 0.;
	theta[3] = 0.;
	theta[RotorAxis] = PsiOff;

	for (node = 0; node < num_labels; node++) {
		if (label[node] >= CurrBladeLabelMax) {
			blade++;
			CurrBladeLabelMax += BladeLabelDelta;
			psi_m += PsiDelta;
			cospsi_m = cos(psi_m);
			sinpsi_m = sin(psi_m);
			m1m *= -1;

			theta[RotorAxis] += PsiDelta;

			# printf("### blade=%d theta[RotorAxis]=%e\n", blade, theta[RotorAxis]);

			vec2R(theta[1], theta[2], theta[3], Rpsi);
			mat3_mul_mat3(R0, Rpsi, Rb)
		}

		node2 = node%(num_labels/NumBlades);

		# position
		Xr[1] = Pos[node, 1] - X0[1];
		Xr[2] = Pos[node, 2] - X0[2];
		Xr[3] = Pos[node, 3] - X0[3];

		mat3T_mul_vec3(Rb, Xr, X);

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

		mat3T_mul_mat3(Rb, Rn, R);
		R2vec(R, T);

		# printf("### label=%d X={%e %e %e}\n", label[node], X[1], X[2], X[3]);
		# printf("### label=%d theta={%e %e %e}\n", label[node], T[1], T[2], T[3]);

		# velocity
		vec3_cross_vec3(W0, Xr, Vtmp);

		Vr[1] = Vel[node, 1] - V0[1] - Vtmp[1];
		Vr[2] = Vel[node, 2] - V0[2] - Vtmp[2];
		Vr[3] = Vel[node, 3] - V0[3] - Vtmp[3];

		mat3T_mul_vec3(Rb, Vr, V);

		# angular velocity
		Wr[1] = Ome[node, 1] - W0[1];
		Wr[2] = Ome[node, 2] - W0[2];
		Wr[3] = Ome[node, 3] - W0[3];

		mat3T_mul_vec3(Rb, Wr, W);

		X[node2, 0, 1] += X[1];
		X[node2, 0, 2] += X[2];
		X[node2, 0, 3] += X[3];

		T[node2, 0, 1] += T[1];
		T[node2, 0, 2] += T[2];
		T[node2, 0, 3] += T[3];

		V[node2, 0, 1] += V[1];
		V[node2, 0, 2] += V[2];
		V[node2, 0, 3] += V[3];

		W[node2, 0, 1] += W[1];
		W[node2, 0, 2] += W[2];
		W[node2, 0, 3] += W[3];

		for (n = 1; n <= NumCyc; n++) {
			X[node2, n "c", 1] += cospsi_m*X[1];
			X[node2, n "c", 2] += cospsi_m*X[2];
			X[node2, n "c", 3] += cospsi_m*X[3];

			T[node2, n "c", 1] += cospsi_m*T[1];
			T[node2, n "c", 2] += cospsi_m*T[2];
			T[node2, n "c", 3] += cospsi_m*T[3];

			V[node2, n "c", 1] += cospsi_m*V[1];
			V[node2, n "c", 2] += cospsi_m*V[2];
			V[node2, n "c", 3] += cospsi_m*V[3];

			W[node2, n "c", 1] += cospsi_m*W[1];
			W[node2, n "c", 2] += cospsi_m*W[2];
			W[node2, n "c", 3] += cospsi_m*W[3];

			X[node2, n "s", 1] += sinpsi_m*X[1];
			X[node2, n "s", 2] += sinpsi_m*X[2];
			X[node2, n "s", 3] += sinpsi_m*X[3];

			T[node2, n "s", 1] += sinpsi_m*T[1];
			T[node2, n "s", 2] += sinpsi_m*T[2];
			T[node2, n "s", 3] += sinpsi_m*T[3];

			V[node2, n "s", 1] += sinpsi_m*V[1];
			V[node2, n "s", 2] += sinpsi_m*V[2];
			V[node2, n "s", 3] += sinpsi_m*V[3];

			W[node2, n "s", 1] += sinpsi_m*W[1];
			W[node2, n "s", 2] += sinpsi_m*W[2];
			W[node2, n "s", 3] += sinpsi_m*W[3];
		}

		if (Even != 0) {
			X[node2, Even, 1] += m1m*X[1];
			X[node2, Even, 2] += m1m*X[2];
			X[node2, Even, 3] += m1m*X[3];

			T[node2, Even, 1] += m1m*T[1];
			T[node2, Even, 2] += m1m*T[2];
			T[node2, Even, 3] += m1m*T[3];

			V[node2, Even, 1] += m1m*V[1];
			V[node2, Even, 2] += m1m*V[2];
			V[node2, Even, 3] += m1m*V[3];

			W[node2, Even, 1] += m1m*W[1];
			W[node2, Even, 2] += m1m*W[2];
			W[node2, Even, 3] += m1m*W[3];
		}
	}

	for (node2 = 0; node2 < num_labels/NumBlades; node2++) {
		X[node2, 0, 1] /= NumBlades;
		X[node2, 0, 2] /= NumBlades;
		X[node2, 0, 3] /= NumBlades;

		T[node2, 0, 1] /= NumBlades;
		T[node2, 0, 2] /= NumBlades;
		T[node2, 0, 3] /= NumBlades;

		V[node2, 0, 1] /= NumBlades;
		V[node2, 0, 2] /= NumBlades;
		V[node2, 0, 3] /= NumBlades;

		W[node2, 0, 1] /= NumBlades;
		W[node2, 0, 2] /= NumBlades;
		W[node2, 0, 3] /= NumBlades;

		printf("%8d:0", label[node2]);
		printf(" %13.6e %13.6e %13.6e", X[node2, 0, 1], X[node2, 0, 2], X[node2, 0, 3]);
		printf(" %13.6e %13.6e %13.6e", T[node2, 0, 1], T[node2, 0, 2], T[node2, 0, 3]);
		printf(" %13.6e %13.6e %13.6e", V[node2, 0, 1], V[node2, 0, 2], V[node2, 0, 3]);
		printf(" %13.6e %13.6e %13.6e", W[node2, 0, 1], W[node2, 0, 2], W[node2, 0, 3]);
		printf("\n");

		for (n = 1; n <= NumCyc; n++) {
			X[node2, n "c", 1] *= 2/NumBlades;
			X[node2, n "c", 2] *= 2/NumBlades;
			X[node2, n "c", 3] *= 2/NumBlades;

			T[node2, n "c", 1] *= 2/NumBlades;
			T[node2, n "c", 2] *= 2/NumBlades;
			T[node2, n "c", 3] *= 2/NumBlades;

			V[node2, n "c", 1] *= 2/NumBlades;
			V[node2, n "c", 2] *= 2/NumBlades;
			V[node2, n "c", 3] *= 2/NumBlades;

			W[node2, n "c", 1] *= 2/NumBlades;
			W[node2, n "c", 2] *= 2/NumBlades;
			W[node2, n "c", 3] *= 2/NumBlades;

			printf("%8d:%dc", label[node2], n);
			printf(" %13.6e %13.6e %13.6e", X[node2, n "c", 1], X[node2, n "c", 2], X[node2, n "c", 3]);
			printf(" %13.6e %13.6e %13.6e", T[node2, n "c", 1], T[node2, n "c", 2], T[node2, n "c", 3]);
			printf(" %13.6e %13.6e %13.6e", V[node2, n "c", 1], V[node2, n "c", 2], V[node2, n "c", 3]);
			printf(" %13.6e %13.6e %13.6e", W[node2, n "c", 1], W[node2, n "c", 2], W[node2, n "c", 3]);
			printf("\n");

			X[node2, n "s", 1] *= 2/NumBlades;
			X[node2, n "s", 2] *= 2/NumBlades;
			X[node2, n "s", 3] *= 2/NumBlades;

			T[node2, n "s", 1] *= 2/NumBlades;
			T[node2, n "s", 2] *= 2/NumBlades;
			T[node2, n "s", 3] *= 2/NumBlades;

			V[node2, n "s", 1] *= 2/NumBlades;
			V[node2, n "s", 2] *= 2/NumBlades;
			V[node2, n "s", 3] *= 2/NumBlades;

			W[node2, n "s", 1] *= 2/NumBlades;
			W[node2, n "s", 2] *= 2/NumBlades;
			W[node2, n "s", 3] *= 2/NumBlades;

			printf("%8d:%ds", label[node2], n);
			printf(" %13.6e %13.6e %13.6e", X[node2, n "s", 1], X[node2, n "s", 2], X[node2, n "s", 3]);
			printf(" %13.6e %13.6e %13.6e", T[node2, n "s", 1], T[node2, n "s", 2], T[node2, n "s", 3]);
			printf(" %13.6e %13.6e %13.6e", V[node2, n "s", 1], V[node2, n "s", 2], V[node2, n "s", 3]);
			printf(" %13.6e %13.6e %13.6e", W[node2, n "s", 1], W[node2, n "s", 2], W[node2, n "s", 3]);
			printf("\n");
		}

		if (Even != 0) {
			X[node2, Even, 1] /= NumBlades;
			X[node2, Even, 2] /= NumBlades;
			X[node2, Even, 3] /= NumBlades;

			T[node2, Even, 1] /= NumBlades;
			T[node2, Even, 2] /= NumBlades;
			T[node2, Even, 3] /= NumBlades;

			V[node2, Even, 1] /= NumBlades;
			V[node2, Even, 2] /= NumBlades;
			V[node2, Even, 3] /= NumBlades;

			W[node2, Even, 1] /= NumBlades;
			W[node2, Even, 2] /= NumBlades;
			W[node2, Even, 3] /= NumBlades;

			printf("%8d:%d", label[node2], Even);
			printf(" %13.6e %13.6e %13.6e", X[node2, Even, 1], X[node2, Even, 2], X[node2, Even, 3]);
			printf(" %13.6e %13.6e %13.6e", T[node2, Even, 1], T[node2, Even, 2], T[node2, Even, 3]);
			printf(" %13.6e %13.6e %13.6e", V[node2, Even, 1], V[node2, Even, 2], V[node2, Even, 3]);
			printf(" %13.6e %13.6e %13.6e", W[node2, Even, 1], W[node2, Even, 2], W[node2, Even, 3]);
			printf("\n");
		}
	}
}

# Initialize vars
BEGIN {
	FirstLabel = -1;
	FirstStep = 1;
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

	if (HubNode == "" || int(HubNode) != HubNode || HubNode < 0) {
		printf("Need HubNode\n");
		exit;
	}

	if (ShaftNode == "" || int(ShaftNode) != ShaftNode || ShaftNode < 0) {
		printf("Need ShaftNode\n");
		exit;
	}

	if (HubNode == ShaftNode) {
		printf("HubNode and ShaftNode must differ\n");
		exit;
	}

	if (NumBlades == "" || int(NumBlades) != NumBlades || NumBlades < 3) {
		printf("NumBlades must be >= 3\n");
		exit;
	}

	if (BladeSeq == "" || BladeSeq == "counterclockwise" || BladeSeq == 1) {
		BladeSeq = 0 + 1;

	} else if (BladeSeq == "clockwise" || BladeSeq == -1) {
		BladeSeq = 0 - 1;

	} else {
		printf("Invalid BladeSeq\n");
		exit;
	}

	PsiDelta = BladeSeq*2*3.14159265358979323846/NumBlades;

	if (int(NumBlades/2)*2 == NumBlades) {
		Even = NumBlades/2;
		NumCyc = NumBlades/2 - 1;

	} else {
		Even = 0;
		NumCyc = (NumBlades - 1)/2;
	}

	if (BladeLabelOff == "" || int(BladeLabelOff) != BladeLabelOff || BladeLabelOff < 0) {
		printf("BladeLabelOff must be integer >= 0\n");
		exit;
	}

	if (BladeLabelDelta == "") {
		BladeLabelDelta = BladeLabelOff;

	} else if (int(BladeLabelDelta) != BladeLabelDelta || BladeLabelDelta <= 0) {
		printf("BladeLabelDelta must be integer > 0\n");
		exit;
	}

	BladeLabelMax = BladeLabelOff + NumBlades*BladeLabelDelta;

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

	deg2rad = 3.14159265358979323846/180;
	rad2deg = 180/3.14159265358979323846;
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
	if ($1 == HubNode) {
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

		prepare(X0, R0, V0, W0);

	} else if ($1 == ShaftNode) {
		if (imode == 0) {
			euler2R($5, $6, $7, Rs);

		} else if (imode == 1) {
			vec2R($5, $6, $7, Rs);

		} else if (imode == 2) {
			Rs[1, 1] = $5;
			Rs[1, 2] = $6;
			Rs[1, 3] = $7;

			Rs[2, 1] = $8;
			Rs[2, 2] = $9;
			Rs[2, 3] = $10;

			Rs[3, 1] = $11;
			Rs[3, 2] = $12;
			Rs[3, 3] = $13;
		}
	}

	if ($1 < BladeLabelOff || $1 >= BladeLabelMax) {
		next;
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

	node++;
}

END {
	if (num_labels == 0) {
		num_labels = node;
	}

	output();
}
