/* $Header$ */
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "drive.h"
#include "elec.h"
#include "strnode.h"
#include "elecnode.h"
#include "motor.h"

/*
 * Electric motor: an internal couple between two structural nodes
 * whose value is represented by an internal state that is a current,
 * according to equations:

	M = M0(Phi_e) - (Gain + M1(Phi_e)) * i
	
	C1 = -M
	
	C2 = M

	i1 = - i

	i2 = i
	
	    d i
	L * --- + R * i = - Gain * (Omega2 - Omega1) + V2 - V1
	    d t
 
 * In order to take into account the so called cogging torque or ripple torque 
 * of DC motors, the following data can be provided:
 *   M0(Phi_e) ... variation of the motor torque independent of current
 *   M1(Phi_e) ... variation of the motor torque proportional to current  
 *   p ... number of terminal pairs
 *
 *   Phi_m ... mechanical angle of rotation between rotor and stator
 *   Phi_e = p * Phi_m ... electric angle between rotor field and stator field
 */

const Motor::PrivData Motor::rgPrivData[iNumPrivData] = {
	{1, "omega"},
	{2, "M"},
	{3, "Pmech"},
	{4, "U"},
	{5, "i"},
	{6, "iP"},
	{7, "Pel"},
	{8, "Phimech"},
	{9, "Phiel"},
	{10, "M0"},
	{11, "M1"},
	{12, "R"}
};

Motor::Motor(const unsigned int uL, const DofOwner* pD, 
		const StructNode* pN1, const StructNode* pN2,
		const ElectricNode* pV1, const ElectricNode* pV2,
		const Mat3x3& Rn, doublereal dG,
		const doublereal dL, DriveCaller* dR, const doublereal i0,
		integer p,
		const DriveCaller* pM0,
		const DriveCaller* pM1,
		const flag fOut)
: Elem(uL, fOut), 
Electric(uL, pD, fOut),
pStrNode1(pN1), pStrNode2(pN2), pVoltage1(pV1), pVoltage2(pV2),
Rn(Rn), dGain(dG), dL(dL), dR(dR), p(p), M0(pM0), M1(pM1), M(i0 * dG), i(i0)
{
	const doublereal dU  = dGetVoltage();
	const doublereal omega = dGetOmega();

	iP = (dU - dGain * omega - dR->dGet(fabs(i)) * i) / dL;

	Phi_m = dGetPhiMechanical();
	Phi_e = dGetPhiElectric(Phi_m);
}

Motor::~Motor(void)
{
	NO_OP;
}

Electric::Type Motor::GetElectricType(void) const
{
	return Electric::MOTOR;
}

/* Contributo al file di restart */
std::ostream&
Motor::Restart(std::ostream& out) const
{
	out << "electric: " << GetLabel()
		<< ", motor, "
		<< pStrNode1->GetLabel() << ", "
		<< "orientation, " << Rn << ", "
		<< pStrNode2->GetLabel() << ", "
		<< pVoltage1->GetLabel() << ", "
		<< pVoltage2->GetLabel() << ", "
		<< dGain << ", "
		<< dL << ", ";

	dR.pGetDriveCaller()->Restart(out);

	out << ", initial current, " << i << ", M0, ";

	M0.pGetDriveCaller()->Restart(out);

	out << ", M1, ";

	M1.pGetDriveCaller()->Restart(out);

	out << ", terminal pairs, " << p << ";" << std::endl;

	return out;                
}
   
unsigned int
Motor::iGetNumDof(void) const
{
	return 1;
}

DofOrder::Order
Motor::GetDofType(unsigned int i) const
{
	return DofOrder::DIFFERENTIAL;
}

void
Motor::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 9;
	*piNumCols = 9;
}
      
VariableSubMatrixHandler&
Motor::AssJac(VariableSubMatrixHandler& WorkMat,
		const doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr)
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	integer iStrNode1FirstPosIdx = pStrNode1->iGetFirstPositionIndex() + 3;
	integer iStrNode2FirstPosIdx = pStrNode2->iGetFirstPositionIndex() + 3;
	integer iStrNode1FirstMomIdx = pStrNode1->iGetFirstMomentumIndex() + 3;
	integer iStrNode2FirstMomIdx = pStrNode2->iGetFirstMomentumIndex() + 3;
	integer iElecNode1FirstIndex = pVoltage1->iGetFirstRowIndex() + 1;
	integer iElecNode2FirstIndex = pVoltage2->iGetFirstRowIndex() + 1;
	integer iFirstIndex = iGetFirstIndex() + 1;

	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iStrNode1FirstMomIdx + iCnt);
		WM.PutRowIndex(3+iCnt, iStrNode2FirstMomIdx + iCnt);

		WM.PutColIndex(iCnt, iStrNode1FirstPosIdx + iCnt);
		WM.PutColIndex(3+iCnt, iStrNode2FirstPosIdx + iCnt);
	}

	WM.PutRowIndex(7, iElecNode1FirstIndex);
	WM.PutRowIndex(8, iElecNode2FirstIndex);
	WM.PutRowIndex(9, iFirstIndex);

	WM.PutColIndex(7, iElecNode1FirstIndex);
	WM.PutColIndex(8, iElecNode2FirstIndex);
	WM.PutColIndex(9, iFirstIndex);

	const doublereal i = XCurr(iFirstIndex);
	const Mat3x3& R1 = pStrNode1->GetRCurr();
	const Mat3x3& R2 = pStrNode2->GetRCurr();
	const Mat3x3& R1_0 = pStrNode1->GetRRef();
	const Mat3x3& R2_0 = pStrNode2->GetRRef();
	const Vec3& omega1 = pStrNode1->GetWCurr();
	const Vec3& omega2 = pStrNode2->GetWCurr();
	const Vec3& omega1_0 = pStrNode1->GetWRef();
	const Vec3& omega2_0 = pStrNode2->GetWRef();
	const Vec3& gP1 = pStrNode1->GetgPCurr();
	const Vec3& gP2 = pStrNode2->GetgPCurr();
	const Mat3x3 DeltaR = Rn.MulTM(R1.MulTM(R2));
	const Vec3 dDeltaR21_dg1_T = -R2.GetCol(1).Cross(R1_0 * Rn.GetCol(2));
	const Vec3 dDeltaR21_dg2_T = R2_0.GetCol(1).Cross(R1 * Rn.GetCol(2));
	const Vec3 dDeltaR11_dg1_T = -R2.GetCol(1).Cross(R1_0 * Rn.GetCol(1));
	const Vec3 dDeltaR11_dg2_T = R2_0.GetCol(1).Cross(R1 * Rn.GetCol(1));
	const doublereal a0 = DeltaR(1, 1) * DeltaR(1, 1) + DeltaR(2, 1) * DeltaR(2, 1);
	const Vec3 dPhi_dg1_T = (dDeltaR21_dg1_T * DeltaR(1, 1) - dDeltaR11_dg1_T * DeltaR(2, 1)) / a0;
	const Vec3 dPhi_dg2_T = (dDeltaR21_dg2_T * DeltaR(1, 1) - dDeltaR11_dg2_T * DeltaR(2, 1)) / a0;
	const doublereal dM0_dPhi = M0.dGetP(Phi_e) * p;
	const doublereal dM1_dPhi = M1.dGetP(Phi_e) * p;
	const Vec3 dM0_dg1_T = dPhi_dg1_T * dM0_dPhi;
	const Vec3 dM0_dg2_T = dPhi_dg2_T * dM0_dPhi;
	const Vec3 dM1_dg1_T = dPhi_dg1_T * dM1_dPhi;
	const Vec3 dM1_dg2_T = dPhi_dg2_T * dM1_dPhi;
	const Vec3 dM_dg1_T = dM0_dg1_T + dM1_dg1_T * i;
	const Vec3 dM_dg2_T = dM0_dg2_T + dM1_dg2_T * i;
	const Mat3x3 dC1_dg1 = Mat3x3(MatCross, R1_0 * Rn.GetCol(3) * M) - (R1 * Rn.GetCol(3)).Tens(dM_dg1_T);
	const Mat3x3 dC1_dg2 = (-(R1 * Rn.GetCol(3))).Tens(dM_dg2_T);
	const Vec3 dC1_di = -(R1 * Rn.GetCol(3) * (dGain + M1.dGet(Phi_e)));
	const Vec3 dOmega_dg1_T = (omega1 - omega2).Cross(R1_0 * Rn.GetCol(3)) - (gP1 * 0.5 + omega1_0).Cross(R1 * Rn.GetCol(3));
	const doublereal di1_di = -1;
	const doublereal di2_di = 1;
	const Vec3 dOmega_dgP1_T = -(R1 * Rn.GetCol(3));
	const Vec3 dOmega_dg2_T = (gP2 * 0.5 + omega2_0).Cross(R1 * Rn.GetCol(3));
	const Vec3 dOmega_dgP2_T = R1 * Rn.GetCol(3);
	const Vec3 dfi_dg1_T = dOmega_dg1_T * (-dGain);
	const Vec3 dfi_dgP1_T = dOmega_dgP1_T * (-dGain);
	const Vec3 dfi_dg2_T = dOmega_dg2_T * (-dGain);
	const Vec3 dfi_dgP2_T = dOmega_dgP2_T * (-dGain);
	const doublereal dfi_du1 = -1;
	const doublereal dfi_du2 = 1;
	const doublereal abs_i = fabs(i);
	const doublereal dfi_di = -dR.dGet(abs_i) - dR.dGetP(abs_i) * i;
	const doublereal dfi_diP = -dL;

	WM.Put(1, 1, dC1_dg1 * (-dCoef));
	WM.Put(1, 4, dC1_dg2 * (-dCoef));
	WM.Put(4, 1, dC1_dg1 * dCoef);
	WM.Put(4, 4, dC1_dg2 * dCoef);
	WM.Put(1, 9, dC1_di * (-dCoef));
	WM.Put(4, 9, dC1_di * dCoef);
	WM.PutCoef(7, 9, -dCoef * di1_di);
	WM.PutCoef(8, 9, -dCoef * di2_di);
	WM.PutT(9, 1, -dfi_dgP1_T - dfi_dg1_T * dCoef);
	WM.PutT(9, 4, -dfi_dgP2_T - dfi_dg2_T * dCoef);
	WM.PutCoef(9, 7, -dCoef * dfi_du1);
	WM.PutCoef(9, 8, -dCoef * dfi_du2);
	WM.PutCoef(9, 9, -dfi_diP - dCoef * dfi_di);

	return WorkMat;
}

SubVectorHandler&
Motor::AssRes(SubVectorHandler& WorkVec,
		const doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr)
{
	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	integer iStrNode1FirstIndex = pStrNode1->iGetFirstMomentumIndex() + 3;
	integer iStrNode2FirstIndex = pStrNode2->iGetFirstMomentumIndex() + 3;
	integer iElecNode1FirstIndex = pVoltage1->iGetFirstRowIndex() + 1;
	integer iElecNode2FirstIndex = pVoltage2->iGetFirstRowIndex() + 1;
	integer iFirstIndex = iGetFirstIndex() + 1;

	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iStrNode1FirstIndex + iCnt);
		WorkVec.PutRowIndex(3+iCnt, iStrNode2FirstIndex + iCnt);
	}

	WorkVec.PutRowIndex(7, iElecNode1FirstIndex);
	WorkVec.PutRowIndex(8, iElecNode2FirstIndex);
	WorkVec.PutRowIndex(9, iFirstIndex);

	Phi_m = dGetPhiMechanical();
	Phi_e = dGetPhiElectric(Phi_m);        
	i = XCurr(iFirstIndex);
	iP = XPrimeCurr(iFirstIndex);
	M = M0.dGet(Phi_e) + (dGain + M1.dGet(Phi_e)) * i;

	const Mat3x3& R1 = pStrNode1->GetRCurr();
	const Vec3 n(R1 * Rn.GetCol(3));
	const doublereal dU  = dGetVoltage();
	const Vec3 C(n * M);
	const doublereal Omega = dGetOmega(n);

	WorkVec.Put(1, -C);
	WorkVec.Put(4, C);
	WorkVec.PutCoef(7, -i);
	WorkVec.PutCoef(8, i);
	WorkVec.PutCoef(9, dU - dGain * Omega - dL * iP - dR.dGet(fabs(i)) * i);
      
	return WorkVec;
}

doublereal Motor::dGetVoltage() const {
	return pVoltage2->dGetX() - pVoltage1->dGetX();
}

doublereal Motor::dGetOmega(const Vec3& TmpDir) const {
	return  TmpDir * (pStrNode2->GetWCurr() - pStrNode1->GetWCurr());
}

Vec3 Motor::GetAxisOfRotation() const {
    return pStrNode1->GetRCurr() * Rn.GetCol(3);
}

doublereal Motor::dGetOmega() const {
	return dGetOmega(GetAxisOfRotation());
}

doublereal Motor::dGetPhiMechanical() const {
	const Mat3x3& R1 = pStrNode1->GetRCurr();
	const Mat3x3& R2 = pStrNode2->GetRCurr();
	const Vec3 DeltaR_e1 = Rn.MulTV(R1.MulTV(R2.GetCol(1)));

	return atan2(DeltaR_e1(2), DeltaR_e1(1));
}

doublereal Motor::dGetPhiElectric(doublereal Phi_m) const {
	doublereal Phi_e = fmod(p * Phi_m, 2 * M_PI);

	if (Phi_e < 0.) {
		Phi_e += 2 * M_PI;
	}

	ASSERT(Phi_e >= 0 && Phi_e <= 2 * M_PI);

	return Phi_e;
}

void
Motor::SetInitialValue(VectorHandler& X)
{
	const integer iFirstIndex = iGetFirstIndex() + 1;
 
	X(iFirstIndex) = i;
}

void
Motor::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints* h)
{
	const integer iFirstIndex = iGetFirstIndex() + 1;

	X(iFirstIndex) = i;
	XP(iFirstIndex) = iP;
}

/* *******PER IL SOLUTORE PARALLELO******** */        
/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
 * utile per l'assemblaggio della matrice di connessione fra i dofs */
void Motor::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
	connectedNodes.resize(4);
	connectedNodes[0] = pStrNode1;
	connectedNodes[1] = pStrNode2;
	connectedNodes[2] = pVoltage1;
	connectedNodes[3] = pVoltage2;
};
/* ************************************************ */

unsigned int Motor::iGetNumPrivData(void) const
{
	return iNumPrivData;
}

unsigned int Motor::iGetPrivDataIdx(const char *s) const
{
	for (int i = 0; i < iNumPrivData; ++i ) {
		if (0 == strcmp(rgPrivData[i].name, s)) {
			return rgPrivData[i].index;
		}
	}

	return 0;
}

doublereal Motor::dGetPrivData(unsigned int iIndex) const
{
	switch (iIndex) {
	case 1:
		return dGetOmega();
	case 2:
		return M;
	case 3:
		return dGetOmega() * M;
	case 4:
		return dGetVoltage();
	case 5:
		return i;
	case 6:
		return iP;
	case 7:
		return dGetVoltage() * i;
	case 8:
		return dGetPhiMechanical();
	case 9:
		return dGetPhiElectric(dGetPhiMechanical());
	case 10:
		return M0.dGet(dGetPhiElectric(dGetPhiMechanical()));
	case 11:
		return M1.dGet(dGetPhiElectric(dGetPhiMechanical())) * i;
	case 12:
		return dR.dGet(fabs(i));
	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}
