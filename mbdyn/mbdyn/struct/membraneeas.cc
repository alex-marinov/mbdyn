/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 2011-2015
 *
 * Marco Morandini	<morandini@aero.polimi.it>
 * Riccardo Vescovini	<vescovini@aero.polimi.it>
 * Tommaso Solcia	<solcia@aero.polimi.it>
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
#include "membraneeas.h"
#include "mynewmem.h"

#include "shell.hc"


/***************************************
* 
*
****************************************/
 
doublereal Membrane4EAS::xi_i[Membrane4EAS::NUMIP][2] = {
	{-1. / std::sqrt(3.), -1. / std::sqrt(3.)},
	{ 1. / std::sqrt(3.), -1. / std::sqrt(3.)},
	{ 1. / std::sqrt(3.),  1. / std::sqrt(3.)},
	{-1. / std::sqrt(3.),  1. / std::sqrt(3.)}
};

doublereal Membrane4EAS::w_i[Membrane4EAS::NUMIP] = {
	1.,
	1.,
	1.,
	1.
};

doublereal Membrane4EAS::xi_n[Membrane4EAS::NUMNODES][2] = {
	{ 1.,  1.},
	{-1.,  1.},
	{-1., -1.},
	{ 1., -1.}
};

doublereal Membrane4EAS::xi_0[2] = {0., 0.};

void
Membrane4EAS::UpdateNodalAndAveragePos(void)
{
	for (integer i = 0; i < NUMNODES; i++) {
		xa[i] = pNode[i]->GetXCurr();
	}
}

void
Membrane4EAS::ComputeInitialIptOrientation(void)
{
	UpdateNodalAndAveragePos();
	for (integer i = 0; i < NUMIP; i++) {
		Vec3 t1 = InterpDeriv_xi1(xa, xi_i[i]);
		t1 = t1 / t1.Norm();
		Vec3 t2 = InterpDeriv_xi2(xa, xi_i[i]);
		t2 = t2 / t2.Norm();
		Vec3 t3 = t1.Cross(t2);
		t3 = t3 / t3.Norm();
		t2 = t3.Cross(t1);
		T_i[i] = Mat3x3(t1, t2, t3);
	}
	Vec3 t1 = InterpDeriv_xi1(xa, xi_0);
	t1 = t1 / t1.Norm();
	Vec3 t2 = InterpDeriv_xi2(xa, xi_0);
	t2 = t2 / t2.Norm();
	Vec3 t3 = t1.Cross(t2);
	t3 = t3 / t3.Norm();
	t2 = t3.Cross(t1);
	T_0 = Mat3x3(t1, t2, t3);
}

Membrane4EAS::Membrane4EAS(unsigned int uL,
	const DofOwner* pDO,
	const StructDispNode* pN1, const StructDispNode* pN2,
	const StructDispNode* pN3, const StructDispNode* pN4,
#if 0 // TODO: offset
 	const Vec3& f1, const Vec3& f2,
 	const Vec3& f3, const Vec3& f4,
#endif
#ifdef USE_CL_IN_MEMBRANE
	const ConstitutiveLaw<vh, fmh>** pDTmp, 
#else // ! USE_CL_IN_MEMBRANE
	const fmh& pDTmp,
	const vh& PreStress,
#endif // ! USE_CL_IN_MEMBRANE
	flag fOut)
: 
Elem(uL, fOut), 
Membrane(uL, pDO, fOut),

S_alpha_beta_0(2, 2),               // 2x2

L_alpha_beta_i(NUMIP, fmh(4, 2) ),  // 4x2

B_overline_i(NUMIP, fmh(3, 12) ),   // 12x24

P_i(NUMIP, fmh(3, iGetNumDof()) ),  // 12x(iGetNumDof)

beta(iGetNumDof()),                 // 1x(iGetNumDof)
epsilon_hat(3),                     // 1x12
epsilon(3),                         // 1x12

#ifndef USE_CL_IN_MEMBRANE
bPreStress(PreStress.Norm() > 0.),
PreStress(PreStress),
#endif // ! USE_CL_IN_MEMBRANE

DRef(NUMIP, fmh(3, 3)),            // 12x12
stress_i(NUMIP, vh(3))              // 1x12
{
#ifdef USE_CL_IN_MEMBRANE
	for (integer i = 0; i < NUMIP; i++) {
		pD[i] = 0;
		SAFENEWWITHCONSTRUCTOR(pD[i],
			ConstitutiveLawOwnerType,
			ConstitutiveLawOwnerType(pDTmp[i]));	
	}
#else // ! USE_CL_IN_MEMBRANE
	for (unsigned i = 0; i < NUMIP; i++) {
		DRef[i] = pDTmp;
	}
#endif // ! USE_CL_IN_MEMBRANE

	pNode[NODE1] = pN1;
	pNode[NODE2] = pN2;
	pNode[NODE3] = pN3;
	pNode[NODE4] = pN4;
// 	f[NODE1] = f1;
// 	f[NODE2] = f2;
// 	f[NODE3] = f3;
// 	f[NODE4] = f4;
	ComputeInitialIptOrientation();
// 	UpdateNodalAndAveragePos();
	// copy ref values
//	T0_overline = T_overline;
	T_0_0 = T_0;
	for (integer i = 0; i < NUMNODES; i++) {
		xa_0[i] = xa[i];
	}
	for (integer i = 0; i < NUMIP; i++) {
		T_0_i[i] = T_i[i];
	}


	fmh M_0(3, 3);
	fmh M_0_Inv(3, 3);
	// doublereal pd_M_0[9], *ppd_M_0[3]; fmh M_0(pd_M_0, ppd_M_0, 9, 3, 3);
	// doublereal pd_M_0_Inv[9], *ppd_M_0_Inv[3]; fmh M_0_Inv(pd_M_0_Inv, ppd_M_0_Inv, 9, 3, 3);
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

		M_0(2, 1) = S_alpha_beta_0(2, 1) * S_alpha_beta_0(2, 1);
		M_0(2, 2) = S_alpha_beta_0(2, 2) * S_alpha_beta_0(2, 2);
		M_0(2, 3) = S_alpha_beta_0(2, 2) * S_alpha_beta_0(2, 1);
	
		M_0(3, 1) = S_alpha_beta_0(1, 1) * S_alpha_beta_0(2, 1);
		M_0(3, 2) = S_alpha_beta_0(1, 2) * S_alpha_beta_0(2, 2);
		M_0(3, 3) = S_alpha_beta_0(1, 1) * S_alpha_beta_0(2, 2) +
		            S_alpha_beta_0(1, 2) * S_alpha_beta_0(2, 1);

		// M_0(4, 1) = S_alpha_beta_0(1, 1) * S_alpha_beta_0(2, 1);
		// M_0(4, 2) = S_alpha_beta_0(1, 2) * S_alpha_beta_0(2, 2);
		// M_0(4, 3) = S_alpha_beta_0(1, 2) * S_alpha_beta_0(2, 1);
		// M_0(4, 4) = S_alpha_beta_0(1, 1) * S_alpha_beta_0(2, 2);
		
		// M_0(5, 5) = S_alpha_beta_0(1, 1);
		// M_0(5, 6) = S_alpha_beta_0(1, 2);
		// M_0(6, 5) = S_alpha_beta_0(2, 1);
		// M_0(6, 6) = S_alpha_beta_0(2, 2);
		
		Inv3x3(M_0, M_0_Inv);
	}

	for (integer i = 0; i < NUMIP; i++) {
		fmh L_alpha_B_i(4, 2);
		fmh S_alpha_beta_i(2, 2);
		// doublereal pd_L_alpha_B_i[8], *ppd_L_alpha_B_i[2]; fmh L_alpha_B_i(pd_L_alpha_B_i, ppd_L_alpha_B_i, 8, 4, 2);
		// doublereal pd_S_alpha_beta_i[4], *ppd_S_alpha_beta_i[2]; fmh S_alpha_beta_i(pd_S_alpha_beta_i, ppd_S_alpha_beta_i, 4, 2, 2);
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

		fmh H(3, iGetNumDof());
		// integer i_H = iGetNumDof(); doublereal pd_H[3*i_H], *ppd_H[i_H]; fmh H(pd_H, ppd_H, 3*i_H, 3, i_H);
		doublereal t = xi_i[i][0] * xi_i[i][1];
		H(1, 1) = xi_i[i][0];
		H(1, 2) = t;

		H(2, 3) = xi_i[i][1];
		H(2, 4) = t;
		
		H(3, 5) = xi_i[i][0];
		H(3, 6) = t;
		
		H(3, 7) = xi_i[i][1];
//		H(4, 6) = t;
		
		fmh Perm(3, 3), tmpP(3, iGetNumDof());    // Perm (12x6), tmpP(6xiGetNumDof)
		// doublereal pd_Perm[9], *ppd_Perm[3]; fmh Perm(pd_Perm, ppd_Perm, 9, 3, 3);
		// integer i_tmpP = iGetNumDof(); doublereal pd_tmpP[3*i_tmpP], *ppd_tmpP[i_tmpP]; fmh tmpP(pd_tmpP, ppd_tmpP, 3*i_tmpP, 3, i_tmpP);
// FIXME: CHECK ORDER !!!! apparently no need for permutation
		Perm(1, 1) = 1.;
		Perm(2, 2) = 1.;
		Perm(3, 3) = 1.;
		
		M_0_Inv.MatTMatMul(tmpP, H);
		Perm.MatMatMul(P_i[i], tmpP);
		P_i[i].ScalarMul(alpha_0 / alpha_i[i]);
	}
	// save initial axial values
	for (integer i = 0; i < NUMIP; i++) {
		InterpDeriv(xa, L_alpha_beta_i[i], y_i_1[i], y_i_2[i]);
		fmh F(3, 2);
		// doublereal pd_F[6], *ppd_F[2]; fmh(pd_F, ppd_F, 6, 3, 2);
		F.Put(1, 1, y_i_1[i]);
		F.Put(1, 2, y_i_2[i]);
		fmh FTF(2, 2);
		// doublereal pd_FTF[4], *ppd_FTF[2]; fmh FTF(pd_FTF, ppd_FTF, 4, 2, 2);
		F.MatTMatMul(FTF, F);
		eps_tilde_0_i[i](1) = 0.5*( FTF(1, 1) - 1 );
		eps_tilde_0_i[i](2) = 0.5*( FTF(2, 2) - 1 );
		eps_tilde_0_i[i](3) = FTF(1, 2);
/*
		eps_tilde_1_0_i[i] = T_i[i].MulTV(y_i_1[i]);
		eps_tilde_2_0_i[i] = T_i[i].MulTV(y_i_2[i]);
*/
	}
}

Membrane4EAS::~Membrane4EAS(void)
{
#ifdef USE_CL_IN_MEMBRANE
	for (integer i = 0; i < NUMIP; i++) {
		ASSERT(pD[i] != 0);
		if (pD[i] != 0) {
			SAFEDELETE(pD[i]);
		}
	}
#endif // USE_CL_IN_MEMBRANE
}

SubVectorHandler& Membrane4EAS::AssRes(SubVectorHandler& WorkVec,
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
		for (int iCnt = 1; iCnt <= 3; iCnt++) {
			WorkVec.PutRowIndex(iCnt + 3 * i, iNodeFirstMomIndex + iCnt);
		}
	}

	integer iFirstReactionIndex = iGetFirstIndex();
	for (unsigned int iCnt = 1; iCnt <= iGetNumDof(); iCnt++) {
		WorkVec.PutRowIndex(12 + iCnt, iFirstReactionIndex + iCnt);
	}

	UpdateNodalAndAveragePos();
	for (unsigned int i = 1; i <= iGetNumDof(); i++) {
		beta(i) = XCurr(iFirstReactionIndex + i);
	}


	for (integer i = 0; i < NUMIP; i++) {
		InterpDeriv(xa, L_alpha_beta_i[i], y_i_1[i], y_i_2[i]);
		fmh F(3, 2);
		F.Put(1, 1, y_i_1[i]);
		F.Put(1, 2, y_i_2[i]);
		fmh FTF(2, 2);
		F.MatTMatMul(FTF, F);
		eps_tilde_i[i](1) = 0.5*( FTF(1, 1) - 1 );
		eps_tilde_i[i](2) = 0.5*( FTF(2, 2) - 1 );
		eps_tilde_i[i](3) = FTF(1, 2);
		eps_tilde_i[i] -= eps_tilde_0_i[i];
/*
		eps_tilde_1_i[i] = T_i[i].MulTV(y_i_1[i]) - eps_tilde_1_0_i[i];
		eps_tilde_2_i[i] = T_i[i].MulTV(y_i_2[i]) - eps_tilde_2_0_i[i];
*/

		// parte variabile di B_overline_i
		for (integer n = 0; n < NUMNODES; n++) {

			fmh TMP_1(3, 3), TMP_2(3, 3);
			// TMP_1.PutT(1, 1, T_i[i] * L_alpha_beta_i[i](n + 1, 1));
			// TMP_2.PutT(1, 1, T_i[i] * L_alpha_beta_i[i](n + 1, 2));

			B_overline_i[i].PutCoef(1, 1 + 3 * n, F(1, 1) * L_alpha_beta_i[i](n + 1, 1) );
			B_overline_i[i].PutCoef(1, 2 + 3 * n, F(2, 1) * L_alpha_beta_i[i](n + 1, 1) );
			B_overline_i[i].PutCoef(1, 3 + 3 * n, F(3, 1) * L_alpha_beta_i[i](n + 1, 1) );
			B_overline_i[i].PutCoef(2, 1 + 3 * n, F(1, 2) * L_alpha_beta_i[i](n + 1, 2) );
			B_overline_i[i].PutCoef(2, 2 + 3 * n, F(2, 2) * L_alpha_beta_i[i](n + 1, 2) );
			B_overline_i[i].PutCoef(2, 3 + 3 * n, F(3, 2) * L_alpha_beta_i[i](n + 1, 2) );
			B_overline_i[i].PutCoef(3, 1 + 3 * n, F(1, 1) * L_alpha_beta_i[i](n + 1, 2) + F(1, 2) * L_alpha_beta_i[i](n + 1, 1) );
			B_overline_i[i].PutCoef(3, 2 + 3 * n, F(2, 1) * L_alpha_beta_i[i](n + 1, 2) + F(2, 2) * L_alpha_beta_i[i](n + 1, 1) );
			B_overline_i[i].PutCoef(3, 3 + 3 * n, F(3, 1) * L_alpha_beta_i[i](n + 1, 2) + F(3, 2) * L_alpha_beta_i[i](n + 1, 1) );

/* 
			// delta epsilon_tilde_1_i
			B_overline_i[i].PutT(1, 1 + 3 * n, T_i[i] * L_alpha_beta_i[i](n + 1, 1));
			// delta epsilon_tilde_2_i
			// B_overline_i[i].PutT(4, 1 + 3 * n, T_i[i] * L_alpha_beta_i[i](n + 1, 2));
*/

		}

	}
	

	/* Calcola le azioni interne */
	for (integer i = 0; i < NUMIP; i++) {
		epsilon.Put(1, eps_tilde_i[i]);
		// TODO: recupera epsilon_hat con l'ordine giusto per qua
		P_i[i].MatVecMul(epsilon_hat, beta);





// FIXME: add epsilon_hat to enhance strains
		epsilon += epsilon_hat;



#ifdef USE_CL_IN_MEMBRANE
		pD[i]->Update(epsilon);
		stress_i[i] = pD[i]->GetF();
#else // ! USE_CL_IN_MEMBRANE
		DRef[i].MatVecMul(stress_i[i], epsilon);
		if (bPreStress) {
			stress_i[i] += PreStress;
		}
#endif // ! USE_CL_IN_MEMBRANE
	}
		
	
	//Residuo
	//forze di volume
	MyVectorHandler rd(12);
	MyVectorHandler rbeta(iGetNumDof());
	for (integer i = 0; i < NUMIP; i++) {
		B_overline_i[i].MatTVecMul(rd, stress_i[i]);
		P_i[i].MatTVecMul(rbeta, stress_i[i]);
		
		AssembleVector(WorkVec,  1, rd, -alpha_i[i] * w_i[i]);
		AssembleVector(WorkVec, 13, rbeta, -alpha_i[i] * w_i[i]);
	}

	return WorkVec;
}

// Jacobian matrix assembly
VariableSubMatrixHandler&
Membrane4EAS::AssJac(VariableSubMatrixHandler& WorkMat,
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
		integer iNodeFirstPosIndex = pNode[i]->iGetFirstPositionIndex();
		integer iNodeFirstMomIndex = pNode[i]->iGetFirstMomentumIndex();
		for (int iCnt = 1; iCnt <= 3; iCnt++) {
			WM.PutRowIndex(iCnt + 3 * i, iNodeFirstMomIndex + iCnt);
			WM.PutColIndex(iCnt + 3 * i, iNodeFirstPosIndex + iCnt);
		}
	}
	for (unsigned int iCnt = 1; iCnt <= iGetNumDof(); iCnt++) {
		WM.PutRowIndex(12 + iCnt, iFirstReactionIndex + iCnt);
		WM.PutColIndex(12 + iCnt, iFirstReactionIndex + iCnt);
	}

	// tangente
	FullMatrixHandler Km(12, 12);
	FullMatrixHandler K_beta_q(iGetNumDof(), 12);
	FullMatrixHandler K_beta_beta(iGetNumDof(), iGetNumDof());

	FullMatrixHandler Ktm(12, 3);
	
	FullMatrixHandler Ktbetaq(iGetNumDof(), 3);
	
	FullMatrixHandler C(3, 3);
	for (integer i = 0; i < NUMIP; i++) {

#ifdef USE_CL_IN_SHELL
		C = pD[i]->GetFDE();
#else // ! USE_CL_IN_SHELL
		C = DRef[i];
#endif // ! USE_CL_IN_SHELL

		B_overline_i[i].MatTMatMul(Ktm, C);

		Ktm.MatMatMul(Km, B_overline_i[i]);
		
		P_i[i].MatTMatMul(Ktbetaq, C);
		Ktbetaq.MatMatMul(K_beta_q, B_overline_i[i]);
		
		Ktbetaq.MatMatMul(K_beta_beta, P_i[i]);
		
		{
			Vec3 n1;
			ExtractVec3(n1, stress_i[i], 1);

			Vec3 Tn1 = T_i[i] * n1;
			Vec3 Tn2 = T_i[i] * n1;
			Vec3 Tm1 = T_i[i] * n1;
			Vec3 Tm2 = T_i[i] * n1;
			for (int n = 0; n < NUMNODES; n++) {
				for (int m = 0; m < NUMNODES; m++) {
				}
			}
		}

		WM.Add(1, 1, Km, alpha_i[i] * w_i[i] * dCoef);

		WM.AddT(1, 13, K_beta_q, alpha_i[i] * w_i[i]);
		WM.Add(13, 1, K_beta_q, alpha_i[i] * w_i[i] * dCoef);
		WM.Add(13, 13, K_beta_beta, alpha_i[i] * w_i[i]);
		
		for (integer l = 0; l < 3; l++) {
			doublereal Fl_1[4], Fl_2[4];
			for (integer n = 0; n < 4; n++) {
				Fl_1[n] = L_alpha_beta_i[i](n + 1, 1);
				Fl_2[n] = L_alpha_beta_i[i](n + 1, 2);
			}
			for (integer n = 0; n < 4; n++) {
				for (integer m = 0; m < 4; m++) {
					WM.IncCoef(3 * n + l + 1, 3 * m + l + 1, Fl_1[n] * stress_i[i](1) * Fl_1[m] * alpha_i[i] * w_i[i] * dCoef);
					WM.IncCoef(3 * n + l + 1, 3 * m + l + 1, Fl_2[n] * stress_i[i](2) * Fl_2[m] * alpha_i[i] * w_i[i] * dCoef);
					WM.IncCoef(3 * n + l + 1, 3 * m + l + 1, (Fl_1[n] * stress_i[i](3) * Fl_2[m] + Fl_2[n] * stress_i[i](3) * Fl_1[m])* alpha_i[i] * w_i[i] * dCoef);
				}
			}
		}

	}

	return WorkMat;

}

// Contribution to restart file
std::ostream&
Membrane4EAS::Restart(std::ostream& out) const
{
	return out << "# not implemented yet" << std::endl;
}

// Initial settings
void
Membrane4EAS::SetValue(DataManager *pDM,
	VectorHandler& /* X */ , VectorHandler& /* XP */ ,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

VariableSubMatrixHandler&
Membrane4EAS::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	WorkMat.SetNullMatrix();
	return WorkMat;
}

SubVectorHandler&
Membrane4EAS::InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr)
{
	WorkVec.Resize(0);
	return WorkVec;
}

// Access to nodes
const StructDispNode*
Membrane4EAS::pGetNode(unsigned int i) const
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
Membrane4EAS::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		if (OH.UseText(OutputHandler::PLATES)) {
			std::ostream& out = OH.Plates();
			out << std::setw(8) << GetLabel();
				// TODO: complete
			for (integer i = 0; i < NUMIP; i++) {
				for (integer r = 1; r <= 3; r++) {
					out << " " << stress_i[i](r);
				}
			}
			out << std::endl;
		}
	}
}


typedef LinearElasticGenericConstitutiveLaw<Membrane::vh, Membrane::fmh> LEGCLMembrane;

Elem*
ReadMembrane4EAS(DataManager* pDM,
	MBDynParser& HP,
	const DofOwner* pDO,
	unsigned int uLabel)
{

	const StructDispNode* pN[4];
	for (unsigned i = 0; i < 4; i++) {
		pN[i] = pDM->ReadNode<const StructDispNode, Node::STRUCTURAL>(HP);
	}

	/* Prestress and prestrain */
	Membrane::vh PreStress(3);
	PreStress.Reset();

#ifdef USE_CL_IN_MEMBRANE
	const ConstitutiveLaw<Membrane::vh, Membrane::fmh> *pD[4];

	TplDriveCaller<Membrane::vh>* pTplDC = new ZeroTplDriveCaller<Membrane::vh>;

	Membrane::fmh S(3, 3);
	for (unsigned ir = 1; ir <= 3; ir++) {
		for (unsigned ic = 1; ic <= 3; ic++) {
			S(ir, ic) = HP.GetReal();
		}
	}

	pD[0] = 0;
	SAFENEWWITHCONSTRUCTOR(pD[0], LEGCLMembrane, LEGCLMembrane(pTplDC, PreStress, S));

	for (unsigned i = 1; i < 4; i++) {
		pD[i] = pD[0]->Copy();
	}
#else // ! USE_CL_IN_MEMBRANE
	Membrane::fmh pD(3, 3);

	if (ReadMembraneConstLaw(HP, pD, PreStress)) {
		silent_cerr("Membrane(" << uLabel << "): unable to read constitutive law" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
#endif // ! USE_CL_IN_MEMBRANE

	flag fOut = pDM->fReadOutput(HP, Elem::PLATE);

	Elem *pEl = 0;
	SAFENEWWITHCONSTRUCTOR(pEl, Membrane4EAS,
		Membrane4EAS(uLabel, pDO,
			pN[0], pN[1], pN[2], pN[3],
#ifdef USE_CL_IN_MEMBRANE
#else // ! USE_CL_IN_MEMBRANE
			pD, PreStress,
#endif // ! USE_CL_IN_MEMBRANE
			fOut));

	std::ostream& out = pDM->GetLogFile();
	out << "membrane4: " << uLabel
		<< " " << pN[0]->GetLabel()
		<< " " << pN[1]->GetLabel()
		<< " " << pN[2]->GetLabel()
		<< " " << pN[3]->GetLabel()
		<< std::endl;

	return pEl;
}

