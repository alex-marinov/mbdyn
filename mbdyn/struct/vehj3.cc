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

/* Cerniera deformabile */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "dataman.h"
#include "vehj3.h"

#include "matvecexp.h"
#include "Rot.hh"

// helper for efficiency
static void
MultRMRtGammam1(Mat6x6& M, const Mat3x3& R, const Mat3x3& Gammam1)
{
	Mat3x3 GM1RT(Gammam1.MulMT(R));

	M.PutMat11(R*(M.GetMat11().MulMT(R)));
	M.PutMat12(R*(M.GetMat12()*GM1RT));
	M.PutMat21(R*(M.GetMat21().MulMT(R)));
	M.PutMat22(R*(M.GetMat22()*GM1RT));
}

/* DeformableJoint - begin */

/* Costruttore non banale */
DeformableJoint::DeformableJoint(unsigned int uL,
	const DofOwner* pDO,
	const ConstitutiveLaw6D* pCL,
	const StructNode* pN1,
	const StructNode* pN2,
	const Vec3& tilde_f1,
	const Vec3& tilde_f2,
	const Mat3x3& tilde_R1h,
	const Mat3x3& tilde_R2h,
	const OrientationDescription& od,
	flag fOut)
: Elem(uL, fOut),
Joint(uL, pDO, fOut),
ConstitutiveLaw6DOwner(pCL),
pNode1(pN1), pNode2(pN2),
tilde_f1(tilde_f1), tilde_f2(tilde_f2),
tilde_R1h(tilde_R1h), tilde_R2h(tilde_R2h),
od(od),
tilde_k(Zero6), tilde_kPrime(Zero6),
bFirstRes(false)
{
	ASSERT(pNode1 != NULL);
	ASSERT(pNode2 != NULL);
	ASSERT(pNode1->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pNode2->GetNodeType() == Node::STRUCTURAL);

	R1h = pNode1->GetRRef()*tilde_R1h;
}


/* Distruttore */
DeformableJoint::~DeformableJoint(void)
{
	NO_OP;
}


/* Contributo al file di restart */
std::ostream&
DeformableJoint::Restart(std::ostream& out) const
{
	Joint::Restart(out) << ", deformable joint, "
		<< pNode1->GetLabel() << ", reference, node, ",
	tilde_f1.Write(out, ", ") << ", hinge, reference, node, 1, ",
	(tilde_R1h.GetVec(1)).Write(out, ", ")
		<< ", 2, ", (tilde_R1h.GetVec(2)).Write(out, ", ") << ", "
		<< pNode2->GetLabel() << ", reference, node, ",
	tilde_f2.Write(out, ", ") << ", hinge, reference, node, 1, ",
	(tilde_R2h.GetVec(1)).Write(out, ", ")
		<< ", 2, ", (tilde_R2h.GetVec(2)).Write(out, ", ") << ", ";
	return pGetConstLaw()->Restart(out) << ';' << std::endl;
}

void
DeformableJoint::OutputPrepare(OutputHandler &OH)
{
	if (bToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			std::string name;
			OutputPrepare_int("Deformable joint", OH, name);
			Var_tilde_d = OH.CreateVar<Vec3>(name + "d",
				MBUnits::Dimensions::Length,
				"relative position in local frame (x, y, z)");
			Var_tilde_dPrime = OH.CreateVar<Vec3>(name + "dPrime",
				MBUnits::Dimensions::Velocity,
				"relative linear velocity in local frame (x, y, z)");
			Var_d = OH.CreateVar<Vec3>(name + "D",
				MBUnits::Dimensions::Length,
				"relative position in global frame (x, y, z)");
			Var_dPrime = OH.CreateVar<Vec3>(name + "DPrime",
				MBUnits::Dimensions::Velocity,
				"relative linear velocity in global frame (x, y, z)");
			Var_Phi = OH.CreateRotationVar(name, "", od, 
				"relative orientation, in joint reference frame");
			Var_Omega = OH.CreateVar<Vec3>(name + "Omega",
				MBUnits::Dimensions::AngularVelocity,
				"local relative angular velocity (x, y, z)");
		}
#endif // USE_NETCDF
	}
}

void
DeformableJoint::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {
		Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
		Mat3x3 R2h(pNode2->GetRCurr()*tilde_R2h);
		Mat3x3 R(R1h.MulTM(R2h));
		Vec3 F(GetF().GetVec1());
		Vec3 M(GetF().GetVec2());
		Vec3 E;

		// angular strain
		switch (od) {
		case EULER_123:
			E = MatR2EulerAngles123(R)*dRaDegr;
			break;

		case EULER_313:
			E = MatR2EulerAngles313(R)*dRaDegr;
			break;

		case EULER_321:
			E = MatR2EulerAngles321(R)*dRaDegr;
			break;

		case ORIENTATION_VECTOR:
			E = RotManip::VecRot(R);
			break;

		case ORIENTATION_MATRIX:
			break;

		default:
			/* impossible */
			break;
		}


#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			Joint::NetCDFOutput(OH, R1h*F, R1h*M, F, M);
			switch (od) {
			case EULER_123:
			case EULER_313:
			case EULER_321:
			case ORIENTATION_VECTOR:
				OH.WriteNcVar(Var_Phi, E);
				break;
			case ORIENTATION_MATRIX:
				OH.WriteNcVar(Var_Phi, R);
				break;
			default:
				/* impossible */
				break;
			}
		}
#endif // USE_NETCDF
		if (OH.UseText(OutputHandler::JOINTS)) {
			Joint::Output(OH.Joints(), "DeformableJoint", GetLabel(),
				F, M, R1h*F, R1h*M);

			// linear strain
			OH.Joints() << " " << tilde_k.GetVec1() << " ";
			
			// Angular strain
			switch (od) {
			case EULER_123:
			case EULER_313:
			case EULER_321:
			case ORIENTATION_VECTOR:
				OH.Joints() << E;
				break;
			case ORIENTATION_MATRIX:
				OH.Joints() << R;
				break;
			default:
				/* impossible */
				break;
			}

			if (GetConstLawType() & ConstLawType::VISCOUS) {
				OH.Joints() << " " << tilde_kPrime;
			}

			OH.Joints() << std::endl;
		}
	}
}

void
DeformableJoint::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	if (ph) {
		for (unsigned i = 0; i < ph->size(); i++) {
			Joint::JointHint *pjh = dynamic_cast<Joint::JointHint *>((*ph)[i]);

			if (pjh) {
				if (dynamic_cast<Joint::OffsetHint<1> *>(pjh)) {
					Mat3x3 R1t(pNode1->GetRCurr().Transpose());
					Vec3 f2(pNode2->GetRCurr()*tilde_f2);
  	 
					tilde_f1 = R1t*(pNode2->GetXCurr() + f2 - pNode1->GetXCurr());
	
				} else if (dynamic_cast<Joint::OffsetHint<2> *>(pjh)) {
					Mat3x3 R2t(pNode2->GetRCurr().Transpose());
					Vec3 f1(pNode1->GetRCurr()*tilde_f1);
  	 
					tilde_f2 = R2t*(pNode1->GetXCurr() + f1 - pNode2->GetXCurr());
	
				} else if (dynamic_cast<Joint::HingeHint<1> *>(pjh)) {
					tilde_R1h = pNode1->GetRCurr().Transpose()*pNode2->GetRCurr()*tilde_R2h;
	
				} else if (dynamic_cast<Joint::HingeHint<2> *>(pjh)) {
					tilde_R2h = pNode2->GetRCurr().Transpose()*pNode1->GetRCurr()*tilde_R1h;
	
				} else if (dynamic_cast<Joint::ReactionsHint *>(pjh)) {
					/* TODO */
				}

				continue;
			}

			/* else, pass to constitutive law */
			ConstitutiveLaw6DOwner::SetValue(pDM, X, XP, ph);
		}
	}
}

/* inverse dynamics capable element */
bool
DeformableJoint::bInverseDynamics(void) const
{
	return true;
}

Hint *
DeformableJoint::ParseHint(DataManager *pDM, const char *s) const
{
	if (strncasecmp(s, "offset{" /*}*/ , STRLENOF("offset{" /*}*/ )) == 0)
	{
		s += STRLENOF("offset{" /*}*/ );

		if (strcmp(&s[1], /*{*/ "}") != 0) {
			return 0;
		}

		switch (s[0]) {
		case '1':
			return new Joint::OffsetHint<1>;

		case '2':
			return new Joint::OffsetHint<2>;
		}

	} else if (strncasecmp(s, "hinge{" /*}*/, STRLENOF("hinge{" /*}*/)) == 0) {
		s += STRLENOF("hinge{" /*}*/);

		if (strcmp(&s[1], /*{*/ "}") != 0) {
			return 0;
		}

		switch (s[0]) {
		case '1':
			return new Joint::HingeHint<1>;

		case '2':
			return new Joint::HingeHint<2>;
		}
	}

	return ConstitutiveLaw6DOwner::ParseHint(pDM, s);
}

unsigned int
DeformableJoint::iGetNumPrivData(void) const
{
	return 18 + ConstitutiveLaw6DOwner::iGetNumPrivData();
}

unsigned int
DeformableJoint::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);

	unsigned idx = 0;

	switch (s[0]) {
	case 'd':
		break;

	case 'r':
		idx += 3;
		break;

	case 'v':
		idx += 6;
		break;

	case 'w':
		idx += 9;
		break;

	case 'F':
		idx += 12;
		break;

	case 'M':
		idx += 15;
		break;

	default:
	{
		size_t l = STRLENOF("constitutiveLaw.");
		if (strncmp(s, "constitutiveLaw.", l) == 0) {
			idx = ConstitutiveLaw6DOwner::iGetPrivDataIdx(&s[l]);
			if (idx > 0) {
				return 18 + idx;
			}
		}
		return 0;
	}
	}

	switch (s[1]) {
	case 'x':
		idx += 1;
		break;
	case 'y':
		idx += 2;
		break;
	case 'z':
		idx += 3;
		break;
	default:
		return 0;
	}

	if (s[2] != '\0') {
		return 0;
	}

	return idx;
}

doublereal
DeformableJoint::dGetPrivData(unsigned int i) const
{
	ASSERT(i > 0);

	ASSERT(i <= iGetNumPrivData());

	switch (i) {
	case 1:
	case 2:
	case 3:
	{
		Vec3 f1(pNode1->GetRCurr()*tilde_f1);
		Vec3 f2(pNode2->GetRCurr()*tilde_f2);
		Mat3x3 R1hT((pNode1->GetRCurr()*tilde_R1h).Transpose());
		Vec3 tilde_d(R1hT*(pNode2->GetXCurr() + f2 - pNode1->GetXCurr() - f1));

		return tilde_d(i);
	}

	case 4:
	case 5:
	case 6:
	{
		Mat3x3 R1hT((pNode1->GetRCurr()*tilde_R1h).Transpose());
		Mat3x3 R2h(pNode2->GetRCurr()*tilde_R2h);

		Vec3 tilde_Theta(RotManip::VecRot(R1hT*R2h));

		return tilde_Theta(i - 3);
	}

	case 7:
	case 8:
	case 9:
	{
		Vec3 f2(pNode2->GetRCurr()*tilde_f2);
		Mat3x3 R1hT(pNode1->GetRCurr().Transpose());
		Vec3 tilde_dPrime(R1hT*(pNode2->GetVCurr() - pNode1->GetVCurr()
					+ (pNode2->GetXCurr() - pNode1->GetXCurr()).Cross(pNode1->GetWCurr())
					- f2.Cross(pNode2->GetWCurr() - pNode1->GetWCurr())));

		return tilde_dPrime(i - 6);
	}

	case 10:
	case 11:
	case 12:
	{
		Mat3x3 R1hT((pNode1->GetRCurr()*tilde_R1h).Transpose());
		Vec3 tilde_Omega = R1hT*(pNode2->GetWCurr() - pNode1->GetWCurr());

		return tilde_Omega(i - 9);
	}

	case 13:
	case 14:
	case 15:
	case 16:
	case 17:
	case 18:
		return GetF()(i - 12);

	default:
		return ConstitutiveLaw6DOwner::dGetPrivData(i - 18);
	}
}

void
DeformableJoint::AssMatCommon(FullSubMatrixHandler& WM, doublereal dCoef)
{
	d2 = pNode2->GetRCurr()*tilde_f2;
	d1 = pNode2->GetXCurr() + d2 - pNode1->GetXCurr();

	Vec3 FTmp(F.GetVec1()*dCoef);
	Mat3x3 FCross(MatCross, FTmp);
	Mat3x3 MCross(MatCross, F.GetVec2()*dCoef);

	WM.Add(1, 4, FCross);
	WM.Sub(4, 1, FCross);
	WM.Add(4, 6 + 1, FCross);
	WM.Sub(6 + 1, 4, FCross);

	Mat3x3 MTmp(MatCrossCross, FTmp, d2);

	WM.Add(6 + 4, 6 + 4, MTmp);
	WM.Sub(4, 6 + 4, MTmp);

	MTmp = Mat3x3(MatCrossCross, d2, FTmp) + MCross;

	WM.Sub(6 + 4, 4, MTmp);

	MTmp = Mat3x3(MatCrossCross, d1, FTmp) + MCross;

	WM.Add(4, 4, MTmp);
}

void
DeformableJoint::AssMatElastic(FullSubMatrixHandler& WM, doublereal dCoef,
	const Mat6x6& FDE)
{
	Mat3x3 F_d = FDE.GetMat11()*dCoef;
	Mat3x3 F_theta = FDE.GetMat12()*dCoef;
	Mat3x3 M_d = FDE.GetMat21()*dCoef;
	Mat3x3 M_theta = FDE.GetMat22()*dCoef;

	/* D11 */
	WM.Add(1, 1, F_d);
	WM.Sub(1, 6 + 1, F_d);
	WM.Sub(6 + 1, 1, F_d);
	WM.Add(6 + 1, 6 + 1, F_d);

	/* D11 * [d1 x] */
	Mat3x3 FTmp = F_d*Mat3x3(MatCross, d1) - F_theta;

	WM.Sub(1, 4, FTmp);
	WM.Add(6 + 1, 4, FTmp);

	/* */
	Mat3x3 MTmp(M_d*Mat3x3(MatCross, d1) - M_theta);

	WM.Sub(4, 4, d1.Cross(FTmp) + MTmp);
	WM.Add(6 + 4, 4, d2.Cross(FTmp) + MTmp);

	/* D11 * [d2 x] */
	FTmp = F_d*Mat3x3(MatCross, d2) - F_theta;

	WM.Add(1, 6 + 4, FTmp);
	WM.Sub(6 + 1, 6 + 4, FTmp);

	/* */
	MTmp = M_d*Mat3x3(MatCross, d2) - M_theta;

	WM.Add(4, 6 + 4, d1.Cross(FTmp) + MTmp);
	WM.Sub(6 + 4, 6 + 4, d2.Cross(FTmp) + MTmp);

	/* [d1 x] * D11 */
	FTmp = d1.Cross(F_d) + M_d;

	WM.Add(4, 1, FTmp);
	WM.Sub(4, 6 + 1, FTmp);

	/* [d2 x] * D11 */
	FTmp = d2.Cross(F_d) + M_d;

	WM.Sub(6 + 4, 1, FTmp);
	WM.Add(6 + 4, 6 + 1, FTmp);
}


void
DeformableJoint::AssMatViscous(FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB,
	doublereal dCoef, const Mat6x6& FDEPrime)
{
	const Mat3x3& F_dPrime = FDEPrime.GetMat11();
	const Mat3x3& F_thetaPrime = FDEPrime.GetMat12();
	const Mat3x3& M_dPrime = FDEPrime.GetMat21();
	const Mat3x3& M_thetaPrime = FDEPrime.GetMat22();

	const Vec3& Omega1 = pNode1->GetWCurr();
	const Vec3& Omega2 = pNode2->GetWCurr();

	/* D11 */
	WMB.Add(1, 1, F_dPrime);
	WMB.Sub(1, 6 + 1, F_dPrime);
	WMB.Sub(6 + 1, 1, F_dPrime);
	WMB.Add(6 + 1, 6 + 1, F_dPrime);

	/* D11 * [d1 x] */
	Mat3x3 FTmp = F_dPrime*Mat3x3(MatCross, d1) - F_thetaPrime;

	WMB.Sub(1, 4, FTmp);
	WMB.Add(6 + 1, 4, FTmp);

	/* */
	Mat3x3 MTmp(M_dPrime*Mat3x3(MatCross, d1) - M_thetaPrime);

	WMB.Sub(4, 4, d1.Cross(FTmp) + MTmp);
	WMB.Add(6 + 4, 4, d2.Cross(FTmp) + MTmp);

	/* D11 * [d2 x] */
	FTmp = F_dPrime*Mat3x3(MatCross, d2) - F_thetaPrime;

	WMB.Add(1, 6 + 4, FTmp);
	WMB.Sub(6 + 1, 6 + 4, FTmp);

	/* */
	MTmp = M_dPrime*Mat3x3(MatCross, d2) - M_thetaPrime;

	WMB.Add(4, 6 + 4, d1.Cross(FTmp) + MTmp);
	WMB.Sub(6 + 4, 6 + 4, d2.Cross(FTmp) + MTmp);

	/* [d1 x] * D11 */
	FTmp = d1.Cross(F_dPrime) + M_dPrime;

	WMB.Add(4, 1, FTmp);
	WMB.Sub(4, 6 + 1, FTmp);

	/* [d2 x] * D11 */
	FTmp = d2.Cross(F_dPrime) + M_dPrime;

	WMB.Sub(6 + 4, 1, FTmp);
	WMB.Add(6 + 4, 6 + 1, FTmp);

	// ~~~ o ~~~ o ~~~ o ~~~

	MTmp = F_dPrime*Mat3x3(MatCross, Omega1*dCoef);

	WMA.Sub(1, 1, MTmp);
	WMA.Add(6 + 1, 1, MTmp);
	WMA.Add(1, 6 + 1, MTmp);
	WMA.Sub(6 + 1, 6 + 1, MTmp);

	Mat3x3 A1dP(d1.Cross(F_dPrime) + M_dPrime);
	Mat3x3 A2dP(d2.Cross(F_dPrime) + M_dPrime);
	
	Mat3x3 A1tPw((d1.Cross(F_thetaPrime) + M_thetaPrime)*Mat3x3(MatCross, Omega2*dCoef));
	Mat3x3 A2tPw((d2.Cross(F_thetaPrime) + M_thetaPrime)*Mat3x3(MatCross, Omega2*dCoef));

	Mat3x3 D1(Mat3x3(MatCross, d1Prime*dCoef) - Mat3x3(MatCrossCross, Omega1, d1*dCoef));
	Mat3x3 D2(Mat3x3(MatCross, d2Prime*dCoef) - Mat3x3(MatCrossCross, Omega1, d2*dCoef));

	MTmp = A1dP*Mat3x3(MatCross, Omega1*dCoef);

	WMA.Sub(4, 1, MTmp);
	WMA.Add(4, 6 + 1, MTmp);

	MTmp = A2dP*Mat3x3(MatCross, Omega1*dCoef);

	WMA.Add(6 + 4, 1, MTmp);
	WMA.Sub(6 + 4, 6 + 1, MTmp);

	FTmp = F_thetaPrime*Mat3x3(MatCross, Omega2*dCoef);

	MTmp = F_dPrime*D1 + FTmp;

	WMA.Sub(1, 4, MTmp);
	WMA.Add(6 + 1, 4, MTmp);

	MTmp = F_dPrime*D2 + FTmp;

	WMA.Add(1, 6 + 4, MTmp);
	WMA.Sub(6 + 1, 6 + 4, MTmp);

	WMA.Sub(4, 4, A1dP*D1 + A1tPw);

	WMA.Add(6 + 4, 4, A2dP*D1 + A2tPw);

	WMA.Add(4, 6 + 4, A1dP*D2 + A1tPw);

	WMA.Sub(6 + 4, 6 + 4, A2dP*D2 + A2tPw);
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
DeformableJoint::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	AssMats(WM, WM, dCoef);

	return WorkMat;
}

/* Jacobian matrix assembly - all but Elastic */
void
DeformableJoint::AssMats(VariableSubMatrixHandler& WorkMatA,
	VariableSubMatrixHandler& WorkMatB,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	FullSubMatrixHandler& WMA = WorkMatA.SetFull();
	FullSubMatrixHandler& WMB = WorkMatB.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WMA.ResizeReset(iNumRows, iNumCols);
	WMB.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WMA.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WMA.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WMA.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
		WMA.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);

		WMB.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WMB.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WMB.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
		WMB.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	AssMats(WMA, WMB, 1.);
}

/* assemblaggio residuo */
SubVectorHandler&
DeformableJoint::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera gli indici */
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
	}

	AssVec(WorkVec);

	return WorkVec;
}

/* Inverse Dynamics Residual Assembly */
SubVectorHandler&
DeformableJoint::AssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */,
	const VectorHandler& /* XPrimeCurr */, 
	const VectorHandler& /* XPrimePrimeCurr */, 
	InverseDynamics::Order iOrder)
{	
	ASSERT(iOrder == InverseDynamics::INVERSE_DYNAMICS
		|| (iOrder == InverseDynamics::POSITION && bIsErgonomy()));

	bFirstRes = false;

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera gli indici */
	integer iNode1FirstMomIndex = pNode1->iGetFirstPositionIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstPositionIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
	}

	AssVec(WorkVec);

	return WorkVec;
}

/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
DeformableJoint::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ )
{
	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	AssVec(WorkVec);

	return WorkVec;
}

const MBUnits::Dimensions
DeformableJoint::GetEquationDimension(integer index) const {
	// DOF == 0
	return MBUnits::Dimensions::UnknownDimension;
}
/* DeformableJoint - end */


/* ElasticJoint - begin */

ElasticJoint::ElasticJoint(unsigned int uL,
	const DofOwner* pDO,
	const ConstitutiveLaw6D* pCL,
	const StructNode* pN1,
	const StructNode* pN2,
	const Vec3& tilde_f1,
	const Vec3& tilde_f2,
	const Mat3x3& tilde_R1h,
	const Mat3x3& tilde_R2h,
	const OrientationDescription& od,
	flag fOut)
: Elem(uL, fOut),
DeformableJoint(uL, pDO, pCL, pN1, pN2, tilde_f1, tilde_f2, tilde_R1h, tilde_R2h, od, fOut),
ThetaRef(Zero3)
{
	/*
	 * Chiede la matrice tangente di riferimento
	 * e la porta nel sistema globale
	 */
	FDE = MultRMRt(ConstitutiveLaw6DOwner::GetFDE(), R1h);
}

ElasticJoint::~ElasticJoint(void)
{
	NO_OP;
}

void
ElasticJoint::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
	ConstitutiveLaw6DOwner::AfterConvergence(tilde_k);
}

/* assemblaggio jacobiano */
void
ElasticJoint::AssMats(VariableSubMatrixHandler& WorkMatA,
	VariableSubMatrixHandler& WorkMatB,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	FullSubMatrixHandler& WMA = WorkMatA.SetFull();
	WorkMatB.SetNullMatrix();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WMA.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WMA.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WMA.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WMA.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
		WMA.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	AssMats(WMA, WMA, 1.);
}

void
ElasticJoint::AssMats(FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& /* WMB */ ,
	doublereal dCoef)
{
	AssMatCommon(WMA, dCoef);
	AssMatElastic(WMA, dCoef, FDE);
}

/* Inverse Dynamics Jacobian matrix assembly */
VariableSubMatrixHandler&
ElasticJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr)
{
	ASSERT(bIsErgonomy());

	// HACK?  Need to call AfterPredict() here to update MDE and so
	// const_cast because AfterPredict() does not modify X, XP
	AfterPredict(const_cast<VectorHandler&>(XCurr), const_cast<VectorHandler&>(XCurr));

	return DeformableJoint::AssJac(WorkMat, 1., XCurr, XCurr);
}

/* Inverse Dynamics update */
void
ElasticJoint::Update(const VectorHandler& XCurr, InverseDynamics::Order iOrder)
{
	NO_OP;
}

/* Inverse Dynamics after convergence */
void
ElasticJoint::AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP, const VectorHandler& XPP)
{
	ConstitutiveLaw6DOwner::AfterConvergence(tilde_k);
}
 
void
ElasticJoint::AssVec(SubVectorHandler& WorkVec)
{
	if (bFirstRes) {
		bFirstRes = false;

	} else {
		R1h = pNode1->GetRCurr()*tilde_R1h;

		d2 = pNode2->GetRCurr()*tilde_f2;
		d1 = pNode2->GetXCurr() + d2 - pNode1->GetXCurr();

		Mat3x3 R2h(pNode2->GetRCurr()*tilde_R2h);
		Vec3 f1(pNode1->GetRCurr()*tilde_f1);

		tilde_k = Vec6(R1h.MulTV(d1 - f1), RotManip::VecRot(R1h.MulTM(R2h)));

		ConstitutiveLaw6DOwner::Update(tilde_k);
	}

	F = MultRV(GetF(), R1h);

	WorkVec.Add(1, F.GetVec1());
	WorkVec.Add(4, d1.Cross(F.GetVec1()) + F.GetVec2());
	WorkVec.Sub(6 + 1, F.GetVec1());
	WorkVec.Sub(6 + 4, d2.Cross(F.GetVec1()) + F.GetVec2());
}

void
ElasticJoint::AfterPredict(VectorHandler& /* X */ ,
	VectorHandler& /* XP */ )
{
	/* Calcola le deformazioni, aggiorna il legame costitutivo
	 * e crea la FDE */

	/* Recupera i dati */
	R1h = pNode1->GetRRef()*tilde_R1h;

	/* Calcola la deformazione corrente nel sistema locale (nodo a) */
	ThetaRef = RotManip::VecRot(R1h.MulTM(pNode2->GetRRef()*tilde_R2h));

	/* Calcola l'inversa di Gamma di ThetaRef */
	Mat3x3 GammaRefm1 = RotManip::DRot_I(ThetaRef);

	d2 = pNode2->GetRRef()*tilde_f2;
	d1 = pNode2->GetXCurr() + d2 - pNode1->GetXCurr();
	Vec3 f1(pNode1->GetRRef()*tilde_f1);

	tilde_k = Vec6(R1h.MulTV(d1 - f1), ThetaRef);

	/* Aggiorna il legame costitutivo */
	ConstitutiveLaw6DOwner::Update(tilde_k);

	/* Chiede la matrice tangente di riferimento e la porta
	 * nel sistema globale */
	FDE = ConstitutiveLaw6DOwner::GetFDE();
        MultRMRtGammam1(FDE, R1h, GammaRefm1);

	bFirstRes = true;
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
ElasticJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& /* XCurr */ )
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	AssMats(WM, WM, 1.);

	return WorkMat;
}

/* ElasticJoint - end */


/* ElasticJointInv - begin */

ElasticJointInv::ElasticJointInv(unsigned int uL,
	const DofOwner* pDO,
	const ConstitutiveLaw6D* pCL,
	const StructNode* pN1,
	const StructNode* pN2,
	const Vec3& tilde_f1,
	const Vec3& tilde_f2,
	const Mat3x3& tilde_R1h,
	const Mat3x3& tilde_R2h,
	const OrientationDescription& od,
	flag fOut)
: Elem(uL, fOut),
DeformableJoint(uL, pDO, pCL, pN1, pN2, tilde_f1, tilde_f2, tilde_R1h, tilde_R2h, od, fOut),
ThetaRef(Zero3)
{
	/*
	 * Chiede la matrice tangente di riferimento
	 * e la porta nel sistema globale
	 */
	FDE = MultRMRt(ConstitutiveLaw6DOwner::GetFDE(), R1h);
}

ElasticJointInv::~ElasticJointInv(void)
{
	NO_OP;
}

void
ElasticJointInv::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
	ConstitutiveLaw6DOwner::AfterConvergence(tilde_k);
}

/* assemblaggio jacobiano */
void
ElasticJointInv::AssMats(VariableSubMatrixHandler& WorkMatA,
	VariableSubMatrixHandler& WorkMatB,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	FullSubMatrixHandler& WMA = WorkMatA.SetFull();
	WorkMatB.SetNullMatrix();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WMA.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WMA.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WMA.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WMA.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
		WMA.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	AssMats(WMA, WMA, 1.);
}

void
ElasticJointInv::AssMats(FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& /* WMB */ ,
	doublereal dCoef)
{
	AssMatCommon(WMA, dCoef);
	AssMatElastic(WMA, dCoef, FDE);
}

void
ElasticJointInv::AssVec(SubVectorHandler& WorkVec)
{
	R1h = pNode1->GetRRef()*tilde_R1h;

	Mat3x3 R2h(pNode2->GetRRef()*tilde_R2h);
	Vec3 ThetaCurr = RotManip::VecRot(R1h.Transpose()*R2h);
	Mat3x3 tilde_R = RotManip::Rot(ThetaCurr/2.);
	Mat3x3 hat_R(R1h*tilde_R);

	Mat3x3 hat_I = hat_R*((Eye3 + tilde_R).Inv()).MulMT(hat_R);

	d1 = pNode1->GetRCurr()*tilde_f1;
	d2 = pNode2->GetRCurr()*tilde_f2;
	Vec3 d(pNode2->GetXCurr() + d2 - pNode1->GetXCurr() - d1);

	if (bFirstRes) {
		bFirstRes = false;

	} else {
		tilde_k = Vec6(hat_R.MulTV(d), ThetaCurr);

		ConstitutiveLaw6DOwner::Update(tilde_k);
	}

	F = MultRV(GetF(), R1h);

	Vec3 dCrossF(d.Cross(F.GetVec1()));

	WorkVec.Add(1, F.GetVec1());
	WorkVec.Add(4, d1.Cross(F.GetVec1()) + hat_I*dCrossF + F.GetVec2());
	WorkVec.Sub(6 + 1, F.GetVec1());
	WorkVec.Sub(6 + 4, d2.Cross(F.GetVec1()) - hat_I.MulTV(dCrossF) + F.GetVec2());
}

void
ElasticJointInv::AfterPredict(VectorHandler& /* X */ ,
	VectorHandler& /* XP */ )
{
	/* Calcola le deformazioni, aggiorna il legame costitutivo
	 * e crea la FDE */

	/* Recupera i dati */
	R1h = pNode1->GetRRef()*tilde_R1h;

	/* Calcola la deformazione corrente nel sistema locale (nodo a) */
	ThetaRef = RotManip::VecRot(R1h.MulTM(pNode2->GetRRef()*tilde_R2h));

	/* Calcola l'inversa di Gamma di ThetaRef */
	Mat3x3 GammaRefm1 = RotManip::DRot_I(ThetaRef);

	Mat3x3 hat_R(R1h*RotManip::Rot(ThetaRef/2.));

	d1 = pNode1->GetRCurr()*tilde_f1;
	d2 = pNode2->GetRRef()*tilde_f2;
	Vec3 d = pNode2->GetXCurr() + d2 - pNode1->GetXCurr() - d1;

	tilde_k = Vec6(hat_R.MulTV(d), ThetaRef);

	/* Aggiorna il legame costitutivo */
	ConstitutiveLaw6DOwner::Update(tilde_k);

	/* Chiede la matrice tangente di riferimento e la porta
	 * nel sistema globale */
	FDE = ConstitutiveLaw6DOwner::GetFDE();
        MultRMRtGammam1(FDE, R1h, GammaRefm1);

	bFirstRes = true;
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
ElasticJointInv::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& /* XCurr */ )
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	AssMats(WM, WM, 1.);

	return WorkMat;
}

/* ElasticJointInv - end */


/* ViscousJoint - begin */

ViscousJoint::ViscousJoint(unsigned int uL,
	const DofOwner* pDO,
	const ConstitutiveLaw6D* pCL,
	const StructNode* pN1,
	const StructNode* pN2,
	const Vec3& tilde_f1,
	const Vec3& tilde_f2,
	const Mat3x3& tilde_R1h,
	const Mat3x3& tilde_R2h,
	const OrientationDescription& od,
	flag fOut)
: Elem(uL, fOut),
DeformableJoint(uL, pDO, pCL, pN1, pN2, tilde_f1, tilde_f2, tilde_R1h, tilde_R2h, od, fOut)
{
	/*
	 * Chiede la matrice tangente di riferimento
	 * e la porta nel sistema globale
	 */
	FDEPrime = MultRMRt(ConstitutiveLaw6DOwner::GetFDEPrime(), R1h);
}

ViscousJoint::~ViscousJoint(void)
{
	NO_OP;
}

void
ViscousJoint::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
	ConstitutiveLaw6DOwner::AfterConvergence(tilde_k, tilde_kPrime);
}

void
ViscousJoint::AssMats(FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB, doublereal dCoef)
{
	AssMatCommon(WMA, dCoef);
	AssMatViscous(WMA, WMB, dCoef, FDEPrime);
}

void
ViscousJoint::AssVec(SubVectorHandler& WorkVec)
{
	if (bFirstRes) {
		bFirstRes = false;

	} else {
		R1h = pNode1->GetRCurr()*tilde_R1h;

		d2 = pNode2->GetRCurr()*tilde_f2;
		d1 = pNode2->GetXCurr() + d2 - pNode1->GetXCurr();

		d2Prime = pNode2->GetWCurr().Cross(d2);
		d1Prime = pNode2->GetVCurr() + d2Prime - pNode1->GetVCurr();

		tilde_kPrime = Vec6(R1h.MulTV(d1Prime - pNode1->GetWCurr().Cross(d1)),
			R1h.MulTV(pNode2->GetWCurr() - pNode1->GetWCurr()));

		ConstitutiveLaw6DOwner::Update(Zero6, tilde_kPrime);
	}

	F = MultRV(ConstitutiveLaw6DOwner::GetF(), R1h);

	WorkVec.Add(1, F.GetVec1());
	WorkVec.Add(4, d1.Cross(F.GetVec1()) + F.GetVec2());
	WorkVec.Sub(6 + 1, F.GetVec1());
	WorkVec.Sub(6 + 4, d2.Cross(F.GetVec1()) + F.GetVec2());
}

void
ViscousJoint::AfterPredict(VectorHandler& /* X */ ,
	VectorHandler& /* XP */ )
{
	/* Calcola le deformazioni, aggiorna il legame costitutivo
	 * e crea la FDE */

	/* Recupera i dati */
	R1h = pNode1->GetRRef()*tilde_R1h;

	d2 = pNode2->GetRCurr()*tilde_f2;
	d1 = pNode2->GetXCurr() + d2 - pNode1->GetXCurr();

	d2Prime = pNode2->GetWCurr().Cross(d2);
	d1Prime = pNode2->GetVCurr() + d2Prime - pNode1->GetVCurr();

	tilde_kPrime = Vec6(R1h.MulTV(d1Prime - pNode1->GetWCurr().Cross(d1)),
		R1h.MulTV(pNode2->GetWCurr() - pNode1->GetWCurr()));

	ConstitutiveLaw6DOwner::Update(Zero6, tilde_kPrime);

	/* FIXME: we need to be able to regenerate FDE
	 * if the constitutive law throws ChangedEquationStructure */
	FDEPrime = MultRMRt(ConstitutiveLaw6DOwner::GetFDEPrime(), R1h);

	bFirstRes = true;
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
ViscousJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& /* XCurr */ )
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	AssMats(WM, WM, 1.);

	return WorkMat;
}


/* ViscousJoint - end */


/* ViscoElasticJoint - begin */

ViscoElasticJoint::ViscoElasticJoint(unsigned int uL,
	const DofOwner* pDO,
	const ConstitutiveLaw6D* pCL,
	const StructNode* pN1,
	const StructNode* pN2,
	const Vec3& tilde_f1,
	const Vec3& tilde_f2,
	const Mat3x3& tilde_R1h,
	const Mat3x3& tilde_R2h,
	const OrientationDescription& od,
	flag fOut)
: Elem(uL, fOut),
DeformableJoint(uL, pDO, pCL, pN1, pN2, tilde_f1, tilde_f2, tilde_R1h, tilde_R2h, od, fOut),
ThetaRef(Zero3)
{
	/*
	 * Chiede la matrice tangente di riferimento
	 * e la porta nel sistema globale
	 */
	FDE = MultRMRt(ConstitutiveLaw6DOwner::GetFDE(), R1h);
	FDEPrime = MultRMRt(ConstitutiveLaw6DOwner::GetFDEPrime(), R1h);
}

ViscoElasticJoint::~ViscoElasticJoint(void)
{
	NO_OP;
}

void
ViscoElasticJoint::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
	ConstitutiveLaw6DOwner::AfterConvergence(tilde_k, tilde_kPrime);
}

void
ViscoElasticJoint::AssMats(FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB, doublereal dCoef)
{
	AssMatCommon(WMA, dCoef);
	AssMatElastic(WMA, dCoef, FDE);
	AssMatViscous(WMA, WMB, dCoef, FDEPrime);
}

void
ViscoElasticJoint::AssVec(SubVectorHandler& WorkVec)
{
	if (bFirstRes) {
		bFirstRes = false;

	} else {
		R1h = pNode1->GetRCurr()*tilde_R1h;

		d2 = pNode2->GetRCurr()*tilde_f2;
		d1 = pNode2->GetXCurr() + d2 - pNode1->GetXCurr();

		d2Prime = pNode2->GetWCurr().Cross(d2);
		d1Prime = pNode2->GetVCurr() + d2Prime - pNode1->GetVCurr();

		Mat3x3 R2h(pNode2->GetRCurr()*tilde_R2h);
		Vec3 f1(pNode1->GetRCurr()*tilde_f1);

		tilde_k = Vec6(R1h.MulTV(d1 - f1), RotManip::VecRot(R1h.MulTM(R2h)));

		tilde_kPrime = Vec6(R1h.MulTV(d1Prime - pNode1->GetWCurr().Cross(d1)),
			R1h.MulTV(pNode2->GetWCurr() - pNode1->GetWCurr()));

		ConstitutiveLaw6DOwner::Update(tilde_k, tilde_kPrime);
	}

	F = MultRV(ConstitutiveLaw6DOwner::GetF(), R1h);

	WorkVec.Add(1, F.GetVec1());
	WorkVec.Add(4, d1.Cross(F.GetVec1()) + F.GetVec2());
	WorkVec.Sub(6 + 1, F.GetVec1());
	WorkVec.Sub(6 + 4, d2.Cross(F.GetVec1()) + F.GetVec2());
}

void
ViscoElasticJoint::AfterPredict(VectorHandler& /* X */ ,
	VectorHandler& /* XP */ )
{
	/* Calcola le deformazioni, aggiorna il legame costitutivo
	 * e crea la FDE */

	/* Recupera i dati */
	R1h = pNode1->GetRRef()*tilde_R1h;

	d2 = pNode2->GetRCurr()*tilde_f2;
	d1 = pNode2->GetXCurr() + d2 - pNode1->GetXCurr();

	d2Prime = pNode2->GetWCurr().Cross(d2);
	d1Prime = pNode2->GetVCurr() + d2Prime - pNode1->GetVCurr();

	/* Calcola la deformazione corrente nel sistema locale (nodo a) */
	ThetaRef = RotManip::VecRot(R1h.MulTM(pNode2->GetRRef()*tilde_R2h));

	/* Calcola l'inversa di Gamma di ThetaRef */
	Mat3x3 GammaRefm1 = RotManip::DRot_I(ThetaRef);

	Vec3 f1(pNode1->GetRRef()*tilde_f1);

	tilde_k = Vec6(R1h.MulTV(d1 - f1), ThetaRef);

	tilde_kPrime = Vec6(R1h.MulTV(d1Prime - pNode1->GetWCurr().Cross(d1)),
		R1h.MulTV(pNode2->GetWCurr() - pNode1->GetWCurr()));

	ConstitutiveLaw6DOwner::Update(tilde_k, tilde_kPrime);

	/* Chiede la matrice tangente di riferimento e la porta
	 * nel sistema globale */
	FDE = ConstitutiveLaw6DOwner::GetFDE();
        MultRMRtGammam1(FDE, R1h, GammaRefm1);

	/* FIXME: we need to be able to regenerate FDE
	 * if the constitutive law throws ChangedEquationStructure */
	FDEPrime = MultRMRt(ConstitutiveLaw6DOwner::GetFDEPrime(), R1h);

	bFirstRes = true;
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
ViscoElasticJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& /* XCurr */ )
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	AssMats(WM, WM, 1.);

	return WorkMat;
}


/* ViscoElasticJoint - end */

