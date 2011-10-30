/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 2010
 *
 * Marco Morandini	<morandini@aero.polimi.it>
 * Riccardo Vescovini	<vescovini@aero.polimi.it>
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

/*
 * Inspired by
 * Wojciech Witkowski
 * "4-Node combined shell element with semi-EAS-ANS strain interpolations
 * in 6-parameter shell theories with drilling degrees of freedom"
 * Comput Mech (2009) 43:307­319 DOI 10.1007/s00466-008-0307-x
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "constltp_impl.h"
#include "tpldrive_impl.h"
#include "shelleas.h"
#include "mynewmem.h"

#include "shell.hc"


/***************************************
* 
*
****************************************/
 
doublereal Shell4EAS::xi_i[Shell4EAS::NUMIP][2] = {
	{-1. / std::sqrt(3.), -1. / std::sqrt(3.)},
	{ 1. / std::sqrt(3.), -1. / std::sqrt(3.)},
	{ 1. / std::sqrt(3.),  1. / std::sqrt(3.)},
	{-1. / std::sqrt(3.),  1. / std::sqrt(3.)}


// 	{-0.774596669241483,  -0.774596669241483},
// 	{-0.774596669241483,                  0.},
// 	{-0.774596669241483,   0.774596669241483},
// 	{                0.,  -0.774596669241483},
// 	{                0.,                  0.},
// 	{                0.,   0.774596669241483},
// 	{ 0.774596669241483,  -0.774596669241483},
// 	{ 0.774596669241483,                  0.},
// 	{ 0.774596669241483,   0.774596669241483}

};

doublereal Shell4EAS::w_i[Shell4EAS::NUMIP] = {
	1.,
	1.,
	1.,
	1.

// 	0.555555555555556 * 0.555555555555556,
// 	0.555555555555556 * 0.888888888888889,
// 	0.555555555555556 * 0.555555555555556,
// 	0.888888888888889 * 0.555555555555556,
// 	0.888888888888889 * 0.888888888888889,
// 	0.888888888888889 * 0.555555555555556,
// 	0.555555555555556 * 0.555555555555556,
// 	0.555555555555556 * 0.888888888888889,
// 	0.555555555555556 * 0.555555555555556
	
};

doublereal Shell4EAS::xi_n[Shell4EAS::NUMNODES][2] = {
	{ 1.,  1.},
	{-1.,  1.},
	{-1., -1.},
	{ 1., -1.}
};

doublereal Shell4EAS::xi_0[2] = {0., 0.};

void
Shell4EAS::UpdateNodalAndAveragePosAndOrientation(void)
{
	Mat3x3 T_avg(Zero3x3);
	Mat3x3 Tn[NUMNODES];
	Mat3x3 R_tilde_n[NUMNODES];
	for (integer i = 0; i < NUMNODES; i++) {
// 		xa[i] = pNode[i]->GetXCurr()+pNode[i]->GetRCurr()*f[i];
		xa[i] = pNode[i]->GetXCurr();
		Tn[i] = pNode[i]->GetRCurr() * iTa[i];
// 		T_a_avg += Ta[i];
		T_avg += pNode[i]->GetRRef() * iTa[i];
	}
	T_avg /= 4.;
	//FIXME; alternative solution: polar decomposition of T_a_avg = T0*U
	T_overline = RotManip::Rot(RotManip::VecRot(T_avg));
	for (integer i = 0; i < NUMNODES; i++) {
		R_tilde_n[i] = T_overline.MulTM(Tn[i]);
		phi_tilde_n[i] = RotManip::VecRot(R_tilde_n[i]);
	}
}

void
Shell4EAS::ComputeInitialNodeAndIptOrientation(void)
{
	for (integer i = 0; i < NUMNODES; i++) {
		xa[i] = pNode[i]->GetXCurr();
	}
	for (integer i = 0; i < NUMNODES; i++) {
		Vec3 t1 = InterpDeriv_xi1(xa, xi_n[i]);
		t1 = t1 / t1.Norm();
		Vec3 t2 = InterpDeriv_xi2(xa, xi_n[i]);
		t2 = t2 / t2.Norm();
		Vec3 t3 = t1.Cross(t2);
		t3 = t3 / t3.Norm();
		t2 = t3.Cross(t1);
		iTa[i] = (pNode[i]->GetRCurr()).MulTM(Mat3x3(t1, t2, t3));
	}
	for (integer i = 0; i < NUMIP; i++) {
		iTa_i[i] = Eye3;
	}
	UpdateNodalAndAveragePosAndOrientation();
	InterpolateOrientation();
	for (integer i = 0; i < NUMIP; i++) {
		Vec3 t1 = InterpDeriv_xi1(xa, xi_i[i]);
		t1 = t1 / t1.Norm();
		Vec3 t2 = InterpDeriv_xi2(xa, xi_i[i]);
		t2 = t2 / t2.Norm();
		Vec3 t3 = t1.Cross(t2);
		t3 = t3 / t3.Norm();
		t2 = t3.Cross(t1);
		iTa_i[i] = (T_i[i]).MulTM(Mat3x3(t1, t2, t3));
	}
	InterpolateOrientation();
}

void
Shell4EAS::InterpolateOrientation(void)
{
	Mat3x3 DRot_I_phi_tilde_n_MT_T_overline[NUMNODES];
	for (integer n = 0; n < NUMNODES; n++) {
		DRot_I_phi_tilde_n_MT_T_overline[n] = 
			RotManip::DRot_I(phi_tilde_n[n]).MulMT(T_overline);
	}
	Mat3x3 Ri, Gammai;
	for (integer i = 0; i < NUMIP; i++) {
		phi_tilde_i[i] = Interp(phi_tilde_n, xi_i[i]);
		RotManip::RotAndDRot(phi_tilde_i[i], Ri, Gammai);
		T_i[i] = T_overline * Ri * iTa_i[i];
		Mat3x3 T_overline_Gamma_tilde_i(T_overline * Gammai);
		for (int n = 0; n < NUMNODES; n++) {
			Phi_Delta_i[i][n] = T_overline_Gamma_tilde_i * 
				DRot_I_phi_tilde_n_MT_T_overline[n];
		}
	}
	Vec3 phi_tilde_0 = Interp(phi_tilde_n, xi_0);
	T_0 = T_overline * RotManip::Rot(phi_tilde_0);
}

// void
// Shell4EAS::ComputeIPSEPRotations(void)
// {
// 	for (integer i = 0; i < NUMIP; i++) {
// 		Q_i[i] = T_i[i].MulMT(T_0_i[i]);
// 	}
// 	for (integer i = 0; i < NUMSSEP; i++) {
// 		Q_A[i] = T_A[i].MulMT(T_0_A[i]);
// 	}
// }

void
Shell4EAS::ComputeIPCurvature(void)
{
	Mat3x3 Gamma_I_n_MT_T_overline[NUMNODES];
	for (integer n = 0; n < NUMNODES; n++) {
		Gamma_I_n_MT_T_overline[n] = RotManip::DRot_I(phi_tilde_n[n]).MulMT(T_overline);
	}
	for (integer i = 0; i < NUMIP; i++) {
		Vec3 phi_tilde_1_i(Zero3);
		Vec3 phi_tilde_2_i(Zero3);
		InterpDeriv(phi_tilde_n, L_alpha_beta_i[i], phi_tilde_1_i, phi_tilde_2_i);
		Mat3x3 T_overlineGamma_tilde_i(T_overline * RotManip::DRot(phi_tilde_i[i]));
		k_1_i[i] = T_overlineGamma_tilde_i * phi_tilde_1_i;
		k_2_i[i] = T_overlineGamma_tilde_i * phi_tilde_2_i;
		Mat3x3 tmp1 = T_overline * RotManip::Elle(phi_tilde_i[i], phi_tilde_1_i);
		Mat3x3 tmp2 = T_overline * RotManip::Elle(phi_tilde_i[i], phi_tilde_2_i);
		for (int n = 0; n < NUMNODES; n++) {
			Kappa_delta_i_1[i][n] = tmp1 * Gamma_I_n_MT_T_overline[n] * LI[n](xi_i[i]) + 
				Phi_Delta_i[i][n] * L_alpha_beta_i[i](n + 1, 1);
			Kappa_delta_i_2[i][n] = tmp2 * Gamma_I_n_MT_T_overline[n] * LI[n](xi_i[i]) + 
				Phi_Delta_i[i][n] * L_alpha_beta_i[i](n + 1, 2);
		}
	}
}

Shell4EAS::Shell4EAS(unsigned int uL,
	const DofOwner* pDO,
	const StructNode* pN1, const StructNode* pN2,
	const StructNode* pN3, const StructNode* pN4,
#if 0 // TODO: offset
 	const Vec3& f1, const Vec3& f2,
 	const Vec3& f3, const Vec3& f4,
#endif
	const Mat3x3& R1, const Mat3x3& R2,
	const Mat3x3& R3, const Mat3x3& R4,
#ifdef USE_CL_IN_SHELL
	const ConstitutiveLaw<vh, fmh>** pDTmp, 
#else // ! USE_CL_IN_SHELL
	const fmh& pDTmp,
	const vh& PreStress,
#endif // ! USE_CL_IN_SHELL
	flag fOut)
: 
Elem(uL, fOut), 
Shell(uL, pDO, fOut),

S_alpha_beta_0(2, 2),

L_alpha_beta_i(NUMIP, fmh(4, 2) ),

B_overline_i(NUMIP, fmh(12, 24) ),

P_i(NUMIP, fmh(12, iGetNumDof()) ),

beta(iGetNumDof()),
epsilon_hat(12),
epsilon(12),

#ifndef USE_CL_IN_SHELL
bPreStress(PreStress.Norm() > 0.),
PreStress(PreStress),
#endif // ! USE_CL_IN_SHELL

DRef(NUMIP, fmh(12, 12)),
stress_i(NUMIP, vh(12))
{
#ifdef USE_CL_IN_SHELL
	for (integer i = 0; i < NUMIP; i++) {
		pD[i] = 0;
		SAFENEWWITHCONSTRUCTOR(pD[i],
			ConstitutiveLawOwnerType,
			ConstitutiveLawOwnerType(pDTmp[i]));	
	}
#else // ! USE_CL_IN_SHELL
	for (unsigned i = 0; i < NUMIP; i++) {
		DRef[i] = pDTmp;
	}
#endif // ! USE_CL_IN_SHELL

	pNode[NODE1] = pN1;
	pNode[NODE2] = pN2;
	pNode[NODE3] = pN3;
	pNode[NODE4] = pN4;
// 	f[NODE1] = f1;
// 	f[NODE2] = f2;
// 	f[NODE3] = f3;
// 	f[NODE4] = f4;
	ComputeInitialNodeAndIptOrientation();
// 	UpdateNodalAndAveragePosAndOrientation();
// 	InterpolateOrientation();
	// copy ref values
	T0_overline = T_overline;
	T_0_0 = T_0;
	for (integer i = 0; i < NUMNODES; i++) {
		xa_0[i] = xa[i];
	}
	for (integer i = 0; i < NUMIP; i++) {
		T_0_i[i] = T_i[i];
	}


	fmh M_0(6, 6);
	fmh M_0_Inv(6, 6);
	{
		Vec3 x_1 = InterpDeriv_xi1(xa, xi_0);
		Vec3 x_2 = InterpDeriv_xi2(xa, xi_0);
		S_alpha_beta_0(1, 1) = T_0_0.GetCol(1) * x_1;
		S_alpha_beta_0(2, 1) = T_0_0.GetCol(2) * x_1;
		S_alpha_beta_0(1, 2) = T_0_0.GetCol(1) * x_2;
		S_alpha_beta_0(2, 2) = T_0_0.GetCol(2) * x_2;
		alpha_0 = S_alpha_beta_0(1, 1) * S_alpha_beta_0(2, 2) -
			S_alpha_beta_0(1, 2) * S_alpha_beta_0(2, 1);

		M_0(1, 1) = S_alpha_beta_0(1, 1) * S_alpha_beta_0(1, 1);
		M_0(1, 2) = S_alpha_beta_0(1, 2) * S_alpha_beta_0(1, 2);
		M_0(1, 3) = S_alpha_beta_0(1, 2) * S_alpha_beta_0(1, 1);
		M_0(1, 4) = S_alpha_beta_0(1, 1) * S_alpha_beta_0(1, 2);	

		M_0(2, 1) = S_alpha_beta_0(2, 1) * S_alpha_beta_0(2, 1);
		M_0(2, 2) = S_alpha_beta_0(2, 2) * S_alpha_beta_0(2, 2);
		M_0(2, 3) = S_alpha_beta_0(2, 2) * S_alpha_beta_0(2, 1);
		M_0(2, 4) = S_alpha_beta_0(2, 1) * S_alpha_beta_0(2, 2);
	
		M_0(3, 1) = S_alpha_beta_0(1, 1) * S_alpha_beta_0(2, 1);
		M_0(3, 2) = S_alpha_beta_0(1, 2) * S_alpha_beta_0(2, 2);
		M_0(3, 3) = S_alpha_beta_0(1, 1) * S_alpha_beta_0(2, 2);
		M_0(3, 4) = S_alpha_beta_0(1, 2) * S_alpha_beta_0(2, 1);

		M_0(4, 1) = S_alpha_beta_0(1, 1) * S_alpha_beta_0(2, 1);
		M_0(4, 2) = S_alpha_beta_0(1, 2) * S_alpha_beta_0(2, 2);
		M_0(4, 3) = S_alpha_beta_0(1, 2) * S_alpha_beta_0(2, 1);
		M_0(4, 4) = S_alpha_beta_0(1, 1) * S_alpha_beta_0(2, 2);
		
		M_0(5, 5) = S_alpha_beta_0(1, 1);
		M_0(5, 6) = S_alpha_beta_0(1, 2);
		M_0(6, 5) = S_alpha_beta_0(2, 1);
		M_0(6, 6) = S_alpha_beta_0(2, 2);
		
		InvBlockDiagonal4_2x4_2(M_0, M_0_Inv);
	}

	for (integer i = 0; i < NUMIP; i++) {
		fmh L_alpha_B_i(4, 2);
		fmh S_alpha_beta_i(2, 2);
		Vec3 x_1 = InterpDeriv_xi1(xa, xi_i[i]);
		Vec3 x_2 = InterpDeriv_xi2(xa, xi_i[i]);
		S_alpha_beta_i(1, 1) = T_0_i[i].GetCol(1) * x_1;
		S_alpha_beta_i(2, 1) = T_0_i[i].GetCol(2) * x_1;
		S_alpha_beta_i(1, 2) = T_0_i[i].GetCol(1) * x_2;
		S_alpha_beta_i(2, 2) = T_0_i[i].GetCol(2) * x_2;
		// alpha_i = det(S_alpha_beta_i)
		alpha_i[i] = S_alpha_beta_i(1, 1) * S_alpha_beta_i(2, 2) -
			S_alpha_beta_i(1, 2) * S_alpha_beta_i(2, 1);

		// xi_i_i = S_alpha_beta_i^{-1}
		FullMatrixHandler xi_i_i(2, 2);
		Inv2x2(S_alpha_beta_i, xi_i_i);
		
		for (integer n = 0; n < NUMNODES; n++) {
			for (integer ii = 0; ii < 2; ii++) {
				L_alpha_B_i(n + 1, ii + 1) = LI_J[n][ii](xi_i[i]);
			}
		}
		
		L_alpha_B_i.MatMatMul(L_alpha_beta_i[i], xi_i_i); 

		fmh H(6, iGetNumDof());
		doublereal t = xi_i[i][0] * xi_i[i][1];
		H(1, 1) = xi_i[i][0];
		H(1, 2) = t;

		H(2, 3) = xi_i[i][1];
		H(2, 4) = t;
		
		H(3, 5) = xi_i[i][0];
		H(3, 6) = t;
		
		H(4, 7) = xi_i[i][1];
		H(4, 6) = t;
		
		doublereal a = 2. * xi_i[i][0] * (xi_i[i][1] * xi_i[i][1] - 1.);
		doublereal b = 2. * xi_i[i][1] * (xi_i[i][0] * xi_i[i][0] - 1.);
		H(5, 8) = a;
		H(5, 9) = b;
		H(5, 10) = a * b;
		H(6, 11) = b;
		H(6, 12) = a;
		H(6, 13) = a * b;
// 		doublereal a = 2. * xi_i[i][0] * (xi_i[i][1] * xi_i[i][1] - 1.);
// 		doublereal b = 2. * xi_i[i][1] * (xi_i[i][0] * xi_i[i][0] - 1.);
// 		H(5, 9) = a;
// 		H(5, 10) = b;
// 		H(5, 11) = a * b;
// 		H(6, 12) = b;
// 		H(6, 13) = a;
// 		H(6, 14) = a * b;
		fmh Perm(12, 6), tmpP(6, iGetNumDof());
			// 1, 5, 4, 2, 3, 6
			Perm(1, 1) = 1.;
			Perm(2, 3) = 1.;
			Perm(3, 5) = 1.;
			Perm(4, 4) = 1.;
			Perm(5, 2) = 1.;
			Perm(6, 6) = 1.;
		
		M_0_Inv.MatTMatMul(tmpP, H);
		Perm.MatMatMul(P_i[i], tmpP);
		P_i[i].ScalarMul(alpha_0 / alpha_i[i]);
	}
	// save initial axial values
	ComputeIPCurvature();
	for (integer i = 0; i < NUMIP; i++) {
		InterpDeriv(xa, L_alpha_beta_i[i], y_i_1[i], y_i_2[i]);
		eps_tilde_1_0_i[i] = T_i[i].MulTV(y_i_1[i]);
		eps_tilde_2_0_i[i] = T_i[i].MulTV(y_i_2[i]);
		k_tilde_1_0_i[i] = T_i[i].MulTV(k_1_i[i]);
		k_tilde_2_0_i[i] = T_i[i].MulTV(k_2_i[i]);
	}
}

Shell4EAS::~Shell4EAS(void)
{
#ifdef USE_CL_IN_SHELL
	for (integer i = 0; i < NUMIP; i++) {
		ASSERT(pD[i] != 0);
		if (pD[i] != 0) {
			SAFEDELETE(pD[i]);
		}
	}
#endif // USE_CL_IN_SHELL
}

SubVectorHandler& Shell4EAS::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera e setta gli indici */
	for (integer i = 0; i < 4; i++) {
		integer iNodeFirstMomIndex = pNode[i]->iGetFirstMomentumIndex();
		for (int iCnt = 1; iCnt <= 6; iCnt++) {
			WorkVec.PutRowIndex(iCnt + 6 * i, iNodeFirstMomIndex + iCnt);
		}
	}

	integer iFirstReactionIndex = iGetFirstIndex();
	for (unsigned int iCnt = 1; iCnt <= iGetNumDof(); iCnt++) {
		WorkVec.PutRowIndex(24 + iCnt, iFirstReactionIndex + iCnt);
	}

	UpdateNodalAndAveragePosAndOrientation();
	InterpolateOrientation();
// 	ComputeIPSEPRotations();
	for (unsigned int i = 1; i <= iGetNumDof(); i++) {
		beta(i) = XCurr(iFirstReactionIndex + i);
	}


	ComputeIPCurvature();
	for (integer i = 0; i < NUMIP; i++) {
		InterpDeriv(xa, L_alpha_beta_i[i], y_i_1[i], y_i_2[i]);
		eps_tilde_1_i[i] = T_i[i].MulTV(y_i_1[i]) - eps_tilde_1_0_i[i];
		eps_tilde_2_i[i] = T_i[i].MulTV(y_i_2[i]) - eps_tilde_2_0_i[i];
		k_tilde_1_i[i] = T_i[i].MulTV(k_1_i[i]) - k_tilde_1_0_i[i];
		k_tilde_2_i[i] = T_i[i].MulTV(k_2_i[i]) - k_tilde_2_0_i[i];
		// parte variabile di B_overline_i
		for (integer n = 0; n < NUMNODES; n++) {
			Mat3x3 Phi_Delta_i_n_LI_i = Phi_Delta_i[i][n] * LI[n](xi_i[i]);

			// delta epsilon_tilde_1_i
			B_overline_i[i].PutT(1, 1 + 6 * n, T_i[i] * L_alpha_beta_i[i](n + 1, 1));
			B_overline_i[i].Put(1, 4 + 6 * n, 
				T_i[i].MulTM(Mat3x3(MatCross, y_i_1[i])) * Phi_Delta_i_n_LI_i
			);

			// delta epsilon_tilde_2_i
			B_overline_i[i].PutT(4, 1 + 6 * n, T_i[i] * L_alpha_beta_i[i](n + 1, 2));
			B_overline_i[i].Put(4, 4 + 6 * n, 
				T_i[i].MulTM(Mat3x3(MatCross, y_i_2[i])) * Phi_Delta_i_n_LI_i
			);

			// delta k_tilde_1_i
			Vec3 phi_tilde_1_i(Zero3);
			Vec3 phi_tilde_2_i(Zero3);
			InterpDeriv(phi_tilde_n, L_alpha_beta_i[i], phi_tilde_1_i, phi_tilde_2_i);
			B_overline_i[i].Put(7, 4 + 6 * n, 
				T_i[i].MulTM(Mat3x3(MatCross, k_1_i[i])) * Phi_Delta_i_n_LI_i
				+
				T_i[i].MulTM(Kappa_delta_i_1[i][n])
			);

			// delta k_tilde_2_i
			B_overline_i[i].Put(10, 4 + 6 * n, 
				T_i[i].MulTM(Mat3x3(MatCross, k_2_i[i])) * Phi_Delta_i_n_LI_i
				+
				T_i[i].MulTM(Kappa_delta_i_2[i][n])
			);
		}
	}
	

	/* Calcola le azioni interne */
	for (integer i = 0; i < NUMIP; i++) {
		epsilon.Put(1, eps_tilde_1_i[i]);
		epsilon.Put(4, eps_tilde_2_i[i]);
		epsilon.Put(7, k_tilde_1_i[i]);
		epsilon.Put(10, k_tilde_2_i[i]);
		// TODO: recupera epsilon_hat con l'ordine giusto per qua
		P_i[i].MatVecMul(epsilon_hat, beta);
// 		{
// 			for (int off = 0; off < 4; off++) {
// 				Vec3 e1(epsilon_hat(1 + 3 * off), epsilon_hat(2 + 3 * off), epsilon_hat(3 + 3 *	off));
// 				e1 = T_i[i].MulTV(e1);
// 				epsilon_hat(1 + 3 * off) = e1(1);
// 				epsilon_hat(2 + 3 * off) = e1(2);
// 				epsilon_hat(3 + 3 * off) = e1(3);
// 			}
// 		}

		epsilon += epsilon_hat;
#ifdef USE_CL_IN_SHELL
		pD[i]->Update(epsilon);
		stress_i[i] = pD[i]->GetF();
#else // ! USE_CL_IN_SHELL
		DRef[i].MatVecMul(stress_i[i], epsilon);
		if (bPreStress) {
			stress_i[i] += PreStress;
		}
#endif // ! USE_CL_IN_SHELL
	}
		
	
	//Residuo
	//forze di volume
	MyVectorHandler rd(24);
	MyVectorHandler rbeta(iGetNumDof());
	for (integer i = 0; i < NUMIP; i++) {
		B_overline_i[i].MatTVecMul(rd, stress_i[i]);
		P_i[i].MatTVecMul(rbeta, stress_i[i]);
		
		AssembleVector(WorkVec,  1, rd, -alpha_i[i] * w_i[i]);
		AssembleVector(WorkVec, 25, rbeta, -alpha_i[i] * w_i[i]);
	}

	return WorkVec;
}
#include <iostream>
// Jacobian matrix assembly
VariableSubMatrixHandler&
Shell4EAS::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera e setta gli indici */
	integer iFirstReactionIndex = iGetFirstIndex();
	for (integer i = 0; i < 4; i++) {
		integer iNodeFirstMomIndex = pNode[i]->iGetFirstMomentumIndex();
		integer iNodeFirstPosIndex = pNode[i]->iGetFirstPositionIndex();
		for (int iCnt = 1; iCnt <= 6; iCnt++) {
			WM.PutRowIndex(iCnt + 6 * i, iNodeFirstMomIndex + iCnt);
			WM.PutColIndex(iCnt + 6 * i, iNodeFirstPosIndex + iCnt);
		}
	}
	for (unsigned int iCnt = 1; iCnt <= iGetNumDof(); iCnt++) {
		WM.PutRowIndex(24 + iCnt, iFirstReactionIndex + iCnt);
		WM.PutColIndex(24 + iCnt, iFirstReactionIndex + iCnt);
	}


///////////////////////////////////////
///////////////////////////////////////
///////////////////////////////////////
///////////////////////////////////////
// 	for (integer ip = 0; ip < NUMIP; ip++) {
// /*	Finite difference check of B_overline matrix */
// 		Vec3 eps1, eps2, k1, k2, y1, y2;
// 
// 		fmh B_overline_base(12, 24);
// 		B_overline_base = B_overline_i[ip];
// 		fmh pippo(12, 24);
// 		fmh pippodiff(12, 24);
// 	
// 		MyVectorHandler basedef(12);
// 		MyVectorHandler incdef(12);
// 	
// 		MyVectorHandler inc(50);
// 		MyVectorHandler xprime(50);
// 		inc.Reset();
// 		xprime.Reset();
// 
// 		std::cerr << "-----------   B_overline[" << ip << "]-----------" << std::endl;
// 		doublereal ddd = 0.000001;
// 		doublereal maxdiff = 0.;
// 		std::cerr << "\nxxxxxxxxxxxxxxx\n" << std::endl;
// 		std::cerr << std::fixed << std::setprecision(6) << B_overline_i[ip] << std::endl;
// 		for (integer n = 0; n < NUMNODES; n++) {
// 			Vec3 xbase = xa[n];
// 			Mat3x3 RDelta = pNode[n]->GetRCurr().MulMT(pNode[n]->GetRRef());
// 			Vec3 gbase(CGR_Rot::Param, RDelta); 
// 			for (integer gdl = 1; gdl <= 3; gdl ++) {
// 				inc.PutCoef(pNode[n]->iGetFirstIndex() + 0 + gdl, xbase(gdl));
// 				inc.PutCoef(pNode[n]->iGetFirstIndex() + 3 + gdl, gbase(gdl));
// 			}
// 		}
// 
// 		for (integer n = 0; n < NUMNODES; n++) {
// 			for (integer gdl = 1; gdl <= 6; gdl ++) {
// 				if (n == 1 && gdl == 1) {
// 					std::cerr << "Eccomi!\n";
// 					for (integer nn=0; nn < 4; nn++) {
// 						std::cerr << std::setprecision(12) << "node  " << nn
// 							<< " " << xa[nn] << "\n";
// 					}
// 					std::cerr << std::setprecision(12) << "etilde 1 " << eps_tilde_1_i[ip] << "\n";
// 					std::cerr << std::setprecision(12) << "etilde 2 " << eps_tilde_2_i[ip] << "\n";
// 					std::cerr << std::setprecision(12) << "y_i_1 " << y_i_1[ip] << "\n";
// 					std::cerr << std::setprecision(12) << "y_i_2 " << y_i_2[ip] << "\n";
// 					std::cerr << std::setprecision(12) << "T y_i_1 " << T_i[ip].MulTV(y_i_1[ip]) << "\n";
// 					std::cerr << std::setprecision(12) << "T y_i_2 " << T_i[ip].MulTV(y_i_2[ip]) << "\n";
// 					std::cerr << std::setprecision(12) << "T " << T_i[ip] << "\n";
// 				}
// 				doublereal saved_value = inc(pNode[n]->iGetFirstIndex() + gdl);
// 				inc.IncCoef(pNode[n]->iGetFirstIndex() + gdl, ddd);
// 				const_cast<StructNode*>(pNode[n])->Update(inc, xprime);
// 				UpdateNodalAndAveragePosAndOrientation();
// 				InterpolateOrientation();
// 				ComputeIPCurvature();
// 				InterpDeriv(xa, L_alpha_beta_i[ip], y1, y2);
// 				eps1 = T_i[ip].MulTV(y1) - Vec3(1., 0., 0.);
// 				eps2 = T_i[ip].MulTV(y2) - Vec3(0., 1., 0.);
// 				k1 = T_i[ip].MulTV(k_1_i[ip]) - T_0_i[ip].MulTV(k_1_0_i[ip]);
// 				k2 = T_i[ip].MulTV(k_2_i[ip]) - T_0_i[ip].MulTV(k_2_0_i[ip]);
// 				if (n == 1 && gdl == 1) {
// 					std::cerr << std::setprecision(12) << "e1       " << eps1 << "\n";
// 					std::cerr << std::setprecision(12) << "e2       " << eps2 << "\n";
// 					std::cerr << std::setprecision(12) << "y1    " << y1 << "\n";
// 					std::cerr << std::setprecision(12) << "y2    " << y2 << "\n";
// 					std::cerr << std::setprecision(12) << "T y_1   " << T_i[ip].MulTV(y1) << "\n";
// 					std::cerr << std::setprecision(12) << "T y_2   " << T_i[ip].MulTV(y2) << "\n";
// 					std::cerr << std::setprecision(12) << T_i[ip] << "\n";
// 					std::cerr << ",,,,,,,,,,,,,,,,,\n";
// 					
// 				}
// 				// Backup changes
// 				inc.PutCoef(pNode[n]->iGetFirstPositionIndex() + gdl, saved_value);
// 				const_cast<StructNode*>(pNode[n])->Update(inc, xprime);
// 				UpdateNodalAndAveragePosAndOrientation();
// 				InterpolateOrientation();
// 				ComputeIPCurvature();
// 				InterpDeriv(xa, L_alpha_beta_i[ip], y1, y2);
// 				for (integer comp = 1; comp <= 3; comp++) {
// 					doublereal d1 = (eps1(comp) - eps_tilde_1_i[ip](comp)) / ddd;
// 					doublereal d2 = (eps2(comp) - eps_tilde_2_i[ip](comp)) / ddd;
// 					doublereal d3 = (k1(comp) - k_tilde_1_i[ip](comp)) / ddd;
// 					doublereal d4 = (k2(comp) - k_tilde_2_i[ip](comp)) / ddd;
// 					d1 = std::abs(d1) > 1.E-100 ? d1 : 0.;
// 					d2 = std::abs(d2) > 1.E-100 ? d2 : 0.;
// 					d3 = std::abs(d3) > 1.E-100 ? d3 : 0.;
// 					d4 = std::abs(d4) > 1.E-100 ? d4 : 0.;
// 					pippo(    comp, gdl + 6 * n) =  d1;
// 					pippo(3 + comp, gdl + 6 * n) =  d2;
// 					pippo(6 + comp, gdl + 6 * n) =  d3;
// 					pippo(9 + comp, gdl + 6 * n) =  d4;
// 					
// 					pippodiff(comp, gdl + 6 * n) = 
// 						std::abs(pippo(comp, gdl + 6 * n) - 
// 						B_overline_i[ip](comp, gdl + 6 * n)) > 1.E-30?
// 						pippo(comp, gdl + 6 * n) - 
// 						B_overline_i[ip](comp, gdl + 6 * n) : 0.;
// 					pippodiff(3 + comp, gdl + 6 * n) = 
// 						std::abs(pippo(3 + comp, gdl + 6 * n) - 
// 						B_overline_i[ip](3 + comp, gdl + 6 * n) ) > 1.E-30?
// 						pippo(3 + comp, gdl + 6 * n) - 
// 						B_overline_i[ip](3 + comp, gdl + 6 * n) : 0.;
// 					pippodiff(6 + comp, gdl + 6 * n) = 
// 						std::abs(pippo(6 + comp, gdl + 6 * n) - 
// 						B_overline_i[ip](6 + comp, gdl + 6 * n)) > 1.E-30?
// 						pippo(6 + comp, gdl + 6 * n) - 
// 						B_overline_i[ip](6 + comp, gdl + 6 * n) : 0.;
// 					pippodiff(9 + comp, gdl + 6 * n) = 
// 						std::abs(pippo(9 + comp, gdl + 6 * n) - 
// 						B_overline_i[ip](9 + comp, gdl + 6 * n)) > 1.E-30?
// 						pippo(9 + comp, gdl + 6 * n) - 
// 						B_overline_i[ip](9 + comp, gdl + 6 * n) : 0.;
// 					maxdiff = std::max(maxdiff, std::abs(pippodiff(0 + comp, gdl + 6 * n)));
// 					maxdiff = std::max(maxdiff, std::abs(pippodiff(3 + comp, gdl + 6 * n)));
// 					maxdiff = std::max(maxdiff, std::abs(pippodiff(6 + comp, gdl + 6 * n)));
// 					maxdiff = std::max(maxdiff, std::abs(pippodiff(9 + comp, gdl + 6 * n)));
// 				}	
// 			}
// // 			check rotation linearization
// 			for (integer gdl = 4; gdl <= 6; gdl ++) {
// 				//Mat3x3 RDelta = pNode[n]->GetRCurr().MulMT(pNode[n]->GetRRef());
// 				//Mat3x3 Gamma = RotManip::DRot(RotManip::VecRot(RDelta));
// 				Mat3x3 Gamma = Mat3x3(CGR_Rot::MatG, pNode[n]->GetgCurr());
// 				doublereal saved_value = inc(pNode[n]->iGetFirstIndex() + gdl);
// 				Mat3x3 SavedT = T_i[ip];
// 				inc.IncCoef(pNode[n]->iGetFirstIndex() + gdl, ddd);
// 				const_cast<StructNode*>(pNode[n])->Update(inc, xprime);
// 				UpdateNodalAndAveragePosAndOrientation();
// 				InterpolateOrientation();
// 				ComputeIPCurvature();
// 				Vec3 phi_delta = ((T_i[ip]-SavedT).MulMT(SavedT)).Ax() / ddd;
// 				// Backup changes
// 				inc.PutCoef(pNode[n]->iGetFirstPositionIndex() + gdl, saved_value);
// 				const_cast<StructNode*>(pNode[n])->Update(inc, xprime);
// 				UpdateNodalAndAveragePosAndOrientation();
// 				InterpolateOrientation();
// 				ComputeIPCurvature();
// 				InterpDeriv(xa, L_alpha_beta_i[ip], y1, y2);
// 				std::cerr << std::fixed << std::setprecision(12) << "gdl "
// 					<< gdl << " " << phi_delta << " ; "
// 					<< (Phi_Delta_i[ip][n] * Gamma ).GetCol(gdl - 3) * LI[n](xi_i[ip]) << "\n";
// 			}
// 		}
// // 		for (integer j = 1; j <= pJac->iGetNumCols(); j++) {
// // 			pippo.PutCoef(j, i, 
// // 				std::abs(incsol(j)) >
// // 				1.E-100?incsol(j):0.);
// // 			pippodiff.PutCoef(j, i, 
// // 				std::abs(pippo(j,i)-pJac->operator()(j,i)) >
// // 				1.E-30?pippo(j,i)-pJac->operator()(j,i):0.);
// // 			maxdiff = std::max(maxdiff, std::abs(pippodiff(j, i)));
// // 			// use finite difference Jacobian matrix
// // 			pJac->PutCoef(j, i, incsol(j));
// // 		}
// 	std::cerr << "\n---------------\n" << std::endl;
// 	std::cerr << std::fixed << std::setprecision(6) << pippo << std::endl;
// 	std::cerr << "\n===============\n" << std::endl;
// 	std::cerr << std::fixed << std::setprecision(6) << pippodiff << std::endl;
// 	std::cerr << "\nZZZZZZZZZZZZZZZ\n" << std::endl;
// 	std::cerr << "MAXDIFF: " << std::fixed << std::setprecision(12) << maxdiff << std::endl;
// 	std::cerr << "\nXXXXXXXXXXXXXXX\n" << std::endl;
// 		
// 	}
//////////////////////////////////////
//////////////////////////////////////
//////////////////////////////////////
//////////////////////////////////////






	// tangente

	FullMatrixHandler Km(24, 24);
	FullMatrixHandler K_beta_q(iGetNumDof(), 24);
	FullMatrixHandler K_beta_beta(iGetNumDof(), iGetNumDof());

	FullMatrixHandler Ktm(24, 12);
// 	FullMatrixHandler Kttm(24, 24);
	
	FullMatrixHandler Ktbetaq(iGetNumDof(), 12);
	
	FullMatrixHandler C(12, 12);
	for (integer i = 0; i < NUMIP; i++) {

#ifdef USE_CL_IN_SHELL
		C = pD[i]->GetFDE();
#else // ! USE_CL_IN_SHELL
		C = DRef[i];
#endif // ! USE_CL_IN_SHELL

		B_overline_i[i].MatTMatMul(Ktm, C);
// 		B_overline_i[i].MatTMatMul(Ktm, C);
// 		FullMatrixHandler Gamma(24, 24);
// 		for (integer n = 0; n < NUMNODES; n++) {
// 			Mat3x3 nGamma = Mat3x3(CGR_Rot::MatG, pNode[n]->GetgCurr());
// 			InsertMatrix(Gamma, 1 + 6 * n, 1 + 6 * n, Mat3x3(1.));
// 			InsertMatrix(Gamma, 4 + 6 * n, 4 + 6 * n, nGamma);
// 			InsertMatrix(Gamma, 4 + 6 * n, 4 + 6 * n, Mat3x3(1.));
// 		}
		

		Ktm.MatMatMul(Km, B_overline_i[i]);
// 		Ktm.MatMatMul(Kttm, B_overline_i[i]);
// 		Kttm.MatMatMul(Km, Gamma);
		
		P_i[i].MatTMatMul(Ktbetaq, C);
		Ktbetaq.MatMatMul(K_beta_q, B_overline_i[i]);
		
		Ktbetaq.MatMatMul(K_beta_beta, P_i[i]);
		
		
// 		{
// 			Vec3 n1, n2, m1, m2;
// 			ExtractVec3(n1, stress_i[i],  1);
// 			ExtractVec3(n2, stress_i[i],  4);
// 			ExtractVec3(m1, stress_i[i],  7);
// 			ExtractVec3(m2, stress_i[i], 10);
// 			n1(3) = 0.;
// 			n2(3) = 0.;
// 		
// 			Vec3 Tn1 = T_i[i] * n1;
// 			Vec3 Tn2 = T_i[i] * n2;
// 			Vec3 Tm1 = T_i[i] * m1;
// 			Vec3 Tm2 = T_i[i] * m2;
// 			for (int n = 0; n < NUMNODES; n++) {
// 				for (int m = 0; m < NUMNODES; m++) {
// 					// forze
// 					WM.Add(4 + 6 * n, 4 + 6 * m, 
// 						// 1
// 						Phi_Delta_i[i][n].MulTM(Mat3x3(y_i_1[i])) *
// 							Mat3x3(Tn1) * Phi_Delta_i[i][m] * 
// 							LI[m](xi_i[i]) * LI[n](xi_i[i])
// 							* alpha_i[i] * w_i[i] * dCoef
// 						+
// 						// 2
// 						Phi_Delta_i[i][n].MulTM(Mat3x3(y_i_2[i])) *
// 							Mat3x3(Tn2) * Phi_Delta_i[i][m] * 
// 							LI[m](xi_i[i]) * LI[n](xi_i[i])
// 							* alpha_i[i] * w_i[i] * dCoef
// 					);
// 					WM.Add(4 + 6 * n, 1 + 6 * m,
// 						// 1
// 						Phi_Delta_i[i][n].MulTM(Mat3x3(Tn1)) * 
// 							L_alpha_beta_i[i](m + 1, 1) * LI[n](xi_i[i])
// 							* alpha_i[i] * w_i[i] * dCoef
// 						+
// 						// 2
// 						Phi_Delta_i[i][n].MulTM(Mat3x3(Tn2)) * 
// 							L_alpha_beta_i[i](m + 1, 2) * LI[n](xi_i[i])
// 							* alpha_i[i] * w_i[i] * dCoef
// 					);
// 					WM.Sub(1 + 6 * n, 4 + 6 * m,
// 						// 1
// 						Mat3x3(Tn1) * Phi_Delta_i[i][m] * 
// 							LI[m](xi_i[i]) * L_alpha_beta_i[i](n + 1, 1)
// 							* alpha_i[i] * w_i[i] * dCoef
// 						+
// 						// 2
// 						Mat3x3(Tn2) * Phi_Delta_i[i][m] * 
// 							LI[m](xi_i[i]) * L_alpha_beta_i[i](n + 1, 2)
// 							* alpha_i[i] * w_i[i] * dCoef
// 					);
// 					// pezzo in Elle
// 					WM.Add(4 + 6 * n, 4 + 6 * m,
// 						// 1
// 						-Phi_Delta_i[i][n] * 
// 							RotManip::Elle(
// 								-phi_tilde_i[i], 
// 								T_overline * (Tn1 * Mat3x3(y_i_1[i]))
// 							) * Phi_Delta_i[i][m]
// 							* LI[n](xi_i[i]) * LI[m](xi_i[i])
// 							* alpha_i[i] * w_i[i] * dCoef
// 						// 2
// 						-
// 						Phi_Delta_i[i][n] * 
// 							RotManip::Elle(
// 								-phi_tilde_i[i], 
// 								T_overline * (Tn2 * Mat3x3(y_i_2[i]))
// 							) * Phi_Delta_i[i][m]
// 							* LI[n](xi_i[i]) * LI[m](xi_i[i])
// 							* alpha_i[i] * w_i[i] * dCoef
// 						// 1
// 						+
// 						T_overline.MulMT(RotManip::DRot_IT(phi_tilde_n[n])) *
// 							RotManip::Elle(
// 								-phi_tilde_n[n], 
// 								RotManip::DRot_IT(phi_tilde_n[n]).MulMT(
// 									RotManip::DRot(phi_tilde_i[i])
// 								).MulMT(T_overline) * Tn1 * Mat3x3(y_i_1[i])
// 							) * Phi_Delta_i[i][m]
// 							* LI[n](xi_i[i]) * LI[m](xi_i[i])
// 							* alpha_i[i] * w_i[i] * dCoef
// 						// 2
// 						+
// 						T_overline.MulMT(RotManip::DRot_IT(phi_tilde_n[n])) *
// 							RotManip::Elle(
// 								-phi_tilde_n[n], 
// 								RotManip::DRot_IT(phi_tilde_n[n]).MulMT(
// 									RotManip::DRot(phi_tilde_i[i])
// 								).MulMT(T_overline) * Tn2 * Mat3x3(y_i_2[i])
// 							) * Phi_Delta_i[i][m]
// 							* LI[n](xi_i[i]) * LI[m](xi_i[i])
// 							* alpha_i[i] * w_i[i] * dCoef
// 					);
// 					
// 					// momenti
// 					WM.Add(4 + 6 * n, 4 + 6 * m, 
// 						// 1
// 						Phi_Delta_i[i][n].MulTM(Mat3x3(k_1_i[i])) * Mat3x3(Tm1) *
// 							Phi_Delta_i[i][m] * LI[m](xi_i[i]) * LI[n](xi_i[i])
// 							* alpha_i[i] * w_i[i] * dCoef
// 						+
// 						Phi_Delta_i[i][n].MulTM(Mat3x3(Tm1)) * 
// 							Phi_Delta_i[i][m] * L_alpha_beta_i[i](m + 1, 1) 
// 							* LI[n](xi_i[i])
// 							* alpha_i[i] * w_i[i] * dCoef
// 						
// 						+
// 						
// 						// 2
// 						Phi_Delta_i[i][n].MulTM(Mat3x3(k_2_i[i])) * Mat3x3(Tm2) *
// 							Phi_Delta_i[i][m] * LI[m](xi_i[i]) 
// 							* LI[n](xi_i[i])
// 							* alpha_i[i] * w_i[i] * dCoef
// 						+
// 						Phi_Delta_i[i][n].MulTM(Mat3x3(Tm2)) * 
// 							Phi_Delta_i[i][m] * L_alpha_beta_i[i](m + 1, 2)	
// 							* LI[n](xi_i[i])
// 							* alpha_i[i] * w_i[i] * dCoef
// 					);
// 					// pezzo in Elle
// 					WM.Add(4 + 6 * n, 4 + 6 * m, 
// 						// 1
// 						-
// 						Phi_Delta_i[i][m].MulTMT(
// 							RotManip::Elle(
// 								phi_tilde_i[i],
// 								T_overline.MulTV(Tm1)
// 							)
// 						) *
// 						RotManip::DRot_I(phi_tilde_n[n]).MulMT(
// 							T_overline
// 						)
// 						* LI[m](xi_i[i]) * L_alpha_beta_i[i](n + 1, 1)
// 						* alpha_i[i] * w_i[i] * dCoef
// 						-
// 						T_overline * RotManip::DRot_IT(phi_tilde_n[m]) *
// 						RotManip::Elle(
// 							phi_tilde_i[i],
// 							T_overline.MulTV(Tm1)
// 						) * Phi_Delta_i[i][n]
// 						* LI[n](xi_i[i]) * L_alpha_beta_i[i](m + 1, 1)
// 						* alpha_i[i] * w_i[i] * dCoef
// 						// 2
// 						-
// 						Phi_Delta_i[i][m].MulTMT(
// 							RotManip::Elle(
// 								phi_tilde_i[i],
// 								T_overline.MulTV(Tm2)
// 							)
// 						) *
// 						RotManip::DRot_I(phi_tilde_n[n]).MulMT(
// 							T_overline
// 						)
// 						* LI[m](xi_i[i]) * L_alpha_beta_i[i](n + 1, 2)
// 						* alpha_i[i] * w_i[i] * dCoef
// 						-
// 						T_overline * RotManip::DRot_IT(phi_tilde_n[m]) *
// 						RotManip::Elle(
// 							phi_tilde_i[i],
// 							T_overline.MulTV(Tm2)
// 						) * Phi_Delta_i[i][n]
// 						* LI[n](xi_i[i]) * L_alpha_beta_i[i](m + 1, 2)
// 						* alpha_i[i] * w_i[i] * dCoef
// 					);
// 				}
// 				WM.Add(4 + 6 * n, 4 + 6 * n, 
// 					// 1
// 					T_overline * RotManip::DRot_IT(phi_tilde_n[n]) *
// 						RotManip::Elle(
// 							phi_tilde_n[n],
// 							RotManip::DRot_IT(phi_tilde_n[n]).MulMT(
// 								RotManip::DRot(phi_tilde_i[i])
// 							).MulMT(T_overline) * Tm1
// 						) * RotManip::DRot_I(phi_tilde_n[n])
// 					* L_alpha_beta_i[i](n + 1, 1)
// 					* alpha_i[i] * w_i[i] * dCoef
// 					// 2
// 					+
// 					T_overline * RotManip::DRot_IT(phi_tilde_n[n]) *
// 						RotManip::Elle(
// 							phi_tilde_n[n],
// 							RotManip::DRot_IT(phi_tilde_n[n]).MulMT(
// 								RotManip::DRot(phi_tilde_i[i])
// 							).MulMT(T_overline) * Tm2
// 						) * RotManip::DRot_I(phi_tilde_n[n])
// 					* L_alpha_beta_i[i](n + 1, 2)
// 					* alpha_i[i] * w_i[i] * dCoef
// 				);
// 			}
// 		}


//////////////////////////////////////////////
// 
// VERSIONE DI PIERO
// 
//////////////////////////////////////////////
		{
			Vec3 n1, n2, m1, m2;
			ExtractVec3(n1, stress_i[i],  1);
			ExtractVec3(n2, stress_i[i],  4);
			ExtractVec3(m1, stress_i[i],  7);
			ExtractVec3(m2, stress_i[i], 10);
// 			n1(3) = 0.;
// 			n2(3) = 0.;
		
			Vec3 Tn1 = T_i[i] * n1;
			Vec3 Tn2 = T_i[i] * n2;
			Vec3 Tm1 = T_i[i] * m1;
			Vec3 Tm2 = T_i[i] * m2;
			for (int n = 0; n < NUMNODES; n++) {
				Mat3x3 ToverGammaITn(T_overline * RotManip::DRot_IT(phi_tilde_n[n]));
				for (int m = 0; m < NUMNODES; m++) {
					Mat3x3 ToverGammaITm(T_overline * RotManip::DRot_IT(phi_tilde_n[m]));
					// forze
					WM.Add(4 + 6 * n, 4 + 6 * m, 
						// 1
						Phi_Delta_i[i][n].MulTM(Mat3x3(MatCross, y_i_1[i])) *
							Tn1.Cross(Phi_Delta_i[i][m]) * 
							LI[m](xi_i[i]) * LI[n](xi_i[i])
							* alpha_i[i] * w_i[i] * dCoef
						+
						// 2
						Phi_Delta_i[i][n].MulTM(Mat3x3(MatCross, y_i_2[i])) *
							Tn2.Cross(Phi_Delta_i[i][m]) * 
							LI[m](xi_i[i]) * LI[n](xi_i[i])
							* alpha_i[i] * w_i[i] * dCoef
					);
					WM.Add(4 + 6 * n, 1 + 6 * m,
						// 1
						Phi_Delta_i[i][n].MulTM(Mat3x3(MatCross, Tn1)) * 
							L_alpha_beta_i[i](m + 1, 1) * LI[n](xi_i[i])
							* alpha_i[i] * dCoef
						+
						// 2
						Phi_Delta_i[i][n].MulTM(Mat3x3(MatCross, Tn2)) * 
							L_alpha_beta_i[i](m + 1, 2) * LI[n](xi_i[i])
							* alpha_i[i] * w_i[i] * w_i[i] * dCoef
					);
					WM.Sub(1 + 6 * n, 4 + 6 * m,
						// 1
						Tn1.Cross(Phi_Delta_i[i][m]) * 
							LI[m](xi_i[i]) * L_alpha_beta_i[i](n + 1, 1)
							* alpha_i[i] * w_i[i] * dCoef
						+
						// 2
						Tn2.Cross(Phi_Delta_i[i][m]) * 
							LI[m](xi_i[i]) * L_alpha_beta_i[i](n + 1, 2)
							* alpha_i[i] * w_i[i] * dCoef
					);
					// pezzo in Elle
					WM.Add(4 + 6 * n, 4 + 6 * m,
						// 1
						-ToverGammaITn * 
							RotManip::Elle(
								-phi_tilde_i[i], 
								T_overline.MulTM(Mat3x3(MatCross, Tn1)) * y_i_1[i]
							).MulMT(ToverGammaITm)
							* LI[n](xi_i[i]) * LI[m](xi_i[i])
							* alpha_i[i] * w_i[i] * dCoef
						// 2
						-ToverGammaITn * 
							RotManip::Elle(
								-phi_tilde_i[i], 
								T_overline.MulTM(Mat3x3(MatCross, Tn2)) * y_i_2[i]
							).MulMT(ToverGammaITm)
							* LI[n](xi_i[i]) * LI[m](xi_i[i])
							* alpha_i[i] * w_i[i] * dCoef
					);
					
					// momenti
					WM.Add(4 + 6 * n, 4 + 6 * m, 
						// 1
						Phi_Delta_i[i][n].MulTM(Mat3x3(MatCross, k_1_i[i])) *
							Tm1.Cross(Phi_Delta_i[i][m]) * LI[m](xi_i[i]) * LI[n](xi_i[i])
							* alpha_i[i] * w_i[i] * dCoef
						+
						// 2
						Phi_Delta_i[i][n].MulTM(Mat3x3(MatCross, k_2_i[i])) *
					       		Tm2.Cross(Phi_Delta_i[i][m]) * LI[m](xi_i[i]) 
							* LI[n](xi_i[i])
							* alpha_i[i] * w_i[i] * dCoef
						+
						// 1
						Phi_Delta_i[i][n].MulTM(Mat3x3(MatCross, Tm1)) * 
							Kappa_delta_i_1[i][m]
							* LI[n](xi_i[i])
							* alpha_i[i] * dCoef
						// 2
						+
						Phi_Delta_i[i][n].MulTM(Mat3x3(MatCross, Tm2)) * 
							Kappa_delta_i_2[i][m]
							* LI[n](xi_i[i])
							* alpha_i[i] * w_i[i] * dCoef
					);
					// pezzo in Elle
					WM.Add(4 + 6 * n, 4 + 6 * m, 
						// 1
						-
						Phi_Delta_i[i][m].MulTMT(
							RotManip::Elle(
								phi_tilde_i[i],
								T_overline.MulTV(Tm1)
							)
						).MulMT(ToverGammaITn)
						* LI[m](xi_i[i]) * L_alpha_beta_i[i](n + 1, 1)
						* alpha_i[i] * w_i[i] * dCoef
						// 2
						-
						Phi_Delta_i[i][m].MulTMT(
							RotManip::Elle(
								phi_tilde_i[i],
								T_overline.MulTV(Tm2)
							)
						).MulMT(ToverGammaITn)
						* LI[m](xi_i[i]) * L_alpha_beta_i[i](n + 1, 2)
						* alpha_i[i] * w_i[i] * dCoef
						// 1
						-
						ToverGammaITm *
							RotManip::Elle(
								phi_tilde_i[i],
								T_overline.MulTV(Tm1)
							)
						* Phi_Delta_i[i][n]
						* LI[n](xi_i[i]) * L_alpha_beta_i[i](m + 1, 1)
						* alpha_i[i] * w_i[i] * dCoef
						// 2
						-
						ToverGammaITm *
							RotManip::Elle(
								phi_tilde_i[i],
								T_overline.MulTV(Tm2)
							)
						* Phi_Delta_i[i][n]
						* LI[n](xi_i[i]) * L_alpha_beta_i[i](m + 1, 2)
						* alpha_i[i] * w_i[i] * dCoef
					);
				}
				//pezzo in Elle sono su n per le forze
				WM.Add(4 + 6 * n, 4 + 6 * n,
					// 1
					ToverGammaITn *
						RotManip::Elle(
							-phi_tilde_n[n], 
							RotManip::DRot_IT(phi_tilde_n[n]).MulMT(
								RotManip::DRot(phi_tilde_i[i])
							).MulMT(T_overline) * Tn1.Cross(y_i_1[i])
						).MulMT(ToverGammaITn)
						* LI[n](xi_i[i])
						* alpha_i[i] * w_i[i] * dCoef
					// 2
					+
					ToverGammaITn *
						RotManip::Elle(
							-phi_tilde_n[n], 
							RotManip::DRot_IT(phi_tilde_n[n]).MulMT(
								RotManip::DRot(phi_tilde_i[i])
							).MulMT(T_overline) * Tn2.Cross(y_i_2[i])
						).MulMT(ToverGammaITn)
						* LI[n](xi_i[i])
						* alpha_i[i] * w_i[i] * dCoef
				);
				//pezzo in Elle sono su n per i momenti
				WM.Add(4 + 6 * n, 4 + 6 * n, 
					// 1
					ToverGammaITn * 
						RotManip::Elle(
							phi_tilde_n[n],
							RotManip::DRot_IT(phi_tilde_n[n]).MulMT(
								RotManip::DRot(phi_tilde_i[i])
							).MulMT(T_overline) * Tm1
						).MulMT(ToverGammaITn)
					* L_alpha_beta_i[i](n + 1, 1)
					* alpha_i[i] * w_i[i] * dCoef
					// 2
					+
					ToverGammaITn * 
						RotManip::Elle(
							phi_tilde_n[n],
							RotManip::DRot_IT(phi_tilde_n[n]).MulMT(
								RotManip::DRot(phi_tilde_i[i])
							).MulMT(T_overline) * Tm2
						).MulMT(ToverGammaITn)
					* L_alpha_beta_i[i](n + 1, 2)
					* alpha_i[i] * w_i[i] * dCoef
				);
			}
		}


		
// 		std::cerr << "Kg:\n" << std::fixed << std::setprecision(12) << Kg << std::endl;
		
		
// 		std::cerr << "Km:\n" << std::fixed << std::setprecision(12) << Km << std::endl;
		// AssembleMatrix(WM, 1, 1, Km, alpha_i[i] * w_i[i] * dCoef);

		// AssembleTransposeMatrix(WM, 1, 25, K_beta_q, alpha_i[i] * w_i[i]);
		// AssembleMatrix(WM, 25, 1, K_beta_q, alpha_i[i] * w_i[i] * dCoef);
		// AssembleMatrix(WM, 25, 25, K_beta_beta, alpha_i[i] * w_i[i]);

		WM.Add(1, 1, Km, alpha_i[i] * w_i[i] * dCoef);

		WM.AddT(1, 25, K_beta_q, alpha_i[i] * w_i[i]);
		WM.Add(25, 1, K_beta_q, alpha_i[i] * w_i[i] * dCoef);
		WM.Add(25, 25, K_beta_beta, alpha_i[i] * w_i[i]);
	}
	return WorkMat;

}

// Contribution to restart file
std::ostream&
Shell4EAS::Restart(std::ostream& out) const
{
	return out << "# not implemented yet" << std::endl;
}

// Initial settings
void
Shell4EAS::SetValue(DataManager *pDM,
	VectorHandler& /* X */ , VectorHandler& /* XP */ ,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

VariableSubMatrixHandler&
Shell4EAS::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	WorkMat.SetNullMatrix();
	return WorkMat;
}

SubVectorHandler&
Shell4EAS::InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr)
{
	WorkVec.Resize(0);
	return WorkVec;
}

// Access to nodes
const StructNode*
Shell4EAS::pGetNode(unsigned int i) const
{
	ASSERT(i >= 1 && i <= 4);
	switch (i) {
	case 1:
	case 2:
	case 3:
	case 4:
		return pNode[i - 1];
	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

void
Shell4EAS::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		if (OH.UseText(OutputHandler::PLATES)) {
			std::ostream& out = OH.Plates();
			out << std::setw(8) << GetLabel();
				// TODO: complete
			for (integer i = 0; i < NUMIP; i++) {
				for (integer r = 1; r <= 12; r++) {
					out << " " << stress_i[i](r);
				}
			}
			out << std::endl;
		}
	}
}


typedef LinearElasticGenericConstitutiveLaw<Shell::vh, Shell::fmh> LEGCLShell;

Elem*
ReadShell4EAS(DataManager* pDM,
	MBDynParser& HP,
	const DofOwner* pDO,
	unsigned int uLabel)
{
	silent_cout("Shell4EAS(" << uLabel << "): warning, \"shell4eas\" is deprecated; use \"shell4easans\" instead" << std::endl);

	const StructNode* pN[4];
	Mat3x3 R[4];
	for (unsigned i = 0; i < 4; i++) {
		pN[i] = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);
		ReferenceFrame RF(pN[i]);
		if (HP.IsKeyWord("orientation")) {
			R[i] = HP.GetRotRel(RF);
		} else {
			R[i] = Eye3;
		}
	}

	/* Prestress and prestrain */
	Shell::vh PreStress(12);
	PreStress.Reset();

#ifdef USE_CL_IN_SHELL
	const ConstitutiveLaw<Shell::vh, Shell::fmh> *pD[4];

	TplDriveCaller<Shell::vh>* pTplDC = new ZeroTplDriveCaller<Shell::vh>;

	Shell::fmh S(12, 12);
	for (unsigned ir = 1; ir <= 12; ir++) {
		for (unsigned ic = 1; ic <= 12; ic++) {
			S(ir, ic) = HP.GetReal();
		}
	}

	pD[0] = 0;
	SAFENEWWITHCONSTRUCTOR(pD[0], LEGCLShell, LEGCLShell(pTplDC, PreStress, S));

	for (unsigned i = 1; i < 4; i++) {
		pD[i] = pD[0]->Copy();
	}
#else // ! USE_CL_IN_SHELL
	Shell::fmh pD(12, 12);

	if (ReadShellConstLaw(HP, pD, PreStress)) {
		silent_cerr("Shell(" << uLabel << "): unable to read constitutive law" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
#endif // ! USE_CL_IN_SHELL

	flag fOut = pDM->fReadOutput(HP, Elem::PLATE);

	Elem *pEl = 0;
	SAFENEWWITHCONSTRUCTOR(pEl, Shell4EAS,
		Shell4EAS(uLabel, pDO,
			pN[0], pN[1], pN[2], pN[3],
			R[0], R[1], R[2], R[3],
#ifdef USE_CL_IN_SHELL
#else // ! USE_CL_IN_SHELL
			pD, PreStress,
#endif // ! USE_CL_IN_SHELL
			fOut));

	std::ostream& out = pDM->GetLogFile();
	out << "shell4: " << uLabel
		<< " " << pN[0]->GetLabel()
		<< " " << pN[1]->GetLabel()
		<< " " << pN[2]->GetLabel()
		<< " " << pN[3]->GetLabel()
		<< std::endl;

	return pEl;
}

