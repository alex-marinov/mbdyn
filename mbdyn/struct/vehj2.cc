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
#include "vehj2.h"
#include "Rot.hh"

/* DeformableDispJoint - begin */

/* Costruttore non banale */
DeformableDispJoint::DeformableDispJoint(unsigned int uL,
	const DofOwner* pDO,
	const ConstitutiveLaw3D* pCL,
	const StructNode* pN1,
	const StructNode* pN2,
	const Vec3& tilde_f1,
	const Vec3& tilde_f2,
	const Mat3x3& tilde_R1h,
	const Mat3x3& tilde_R2h,
	flag fOut)
: Elem(uL, fOut),
Joint(uL, pDO, fOut),
ConstitutiveLaw3DOwner(pCL),
pNode1(pN1), pNode2(pN2),
tilde_f1(tilde_f1), tilde_f2(tilde_f2),
tilde_R1h(tilde_R1h), tilde_R2h(tilde_R2h),
tilde_R1hT_tilde_f1(tilde_R1h.Transpose()*tilde_f1),
tilde_d(Zero3), tilde_dPrime(Zero3),
bFirstRes(false), F(Zero3)
# ifdef USE_NETCDFC // netcdfcxx4 has non-pointer vars...
,
Var_tilde_d(0),
Var_tilde_dPrime(0),
Var_d(0),
Var_dPrime(0)
#endif // USE_NETCDFC
{
	ASSERT(pNode1 != NULL);
	ASSERT(pNode2 != NULL);
	ASSERT(pNode1->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pNode2->GetNodeType() == Node::STRUCTURAL);
}

/* Distruttore */
DeformableDispJoint::~DeformableDispJoint(void)
{
	NO_OP;
}

/* assemblaggio jacobiano - F */
void
DeformableDispJoint::AssMatF(FullSubMatrixHandler& WMA,
	const Vec3& d1, const Vec3& d2, doublereal dCoef)
{
	/* force */
	Vec3 FTmp(F*dCoef);

	/* - [ F x ] * dCoef */
	Mat3x3 MTmp(MatCross, FTmp);
	WMA.Add(1, 4, MTmp);
	WMA.Sub(6 + 1, 4, MTmp);
	WMA.Sub(4, 1, MTmp);
	WMA.Add(4, 6 + 1, MTmp);

	/* [ di x ] [ F x ] * dCoef */
	WMA.Add(4, 4, Mat3x3(MatCrossCross, d1, FTmp));
	WMA.Sub(6 + 4, 4, Mat3x3(MatCrossCross, d2, FTmp));

	/* [ F x ] [ d2 x ] * dCoef */
	MTmp = Mat3x3(MatCrossCross, FTmp, d2);
	WMA.Sub(4, 6 + 4, MTmp);
	WMA.Add(6 + 4, 6 + 4, MTmp);
}

/* assemblaggio jacobiano - FDE */
void
DeformableDispJoint::AssMatFDE(FullSubMatrixHandler& WMA,
	const Vec3& d1, const Vec3& d2, doublereal dCoef)
{
	/* F/d */
	Mat3x3 DTmp(FDE*dCoef);

	/* Force equations */

	/* delta x1 */
	WMA.Add(1, 1, DTmp);
	WMA.Sub(6 + 1, 1, DTmp);

	/* delta x2 */
	WMA.Sub(1, 6 + 1, DTmp);
	WMA.Add(6 + 1, 6 + 1, DTmp);

	/* delta g1 */
	Mat3x3 MTmp(DTmp*Mat3x3(MatCross, d1));
	WMA.Sub(1, 4, MTmp);
	WMA.Add(6 + 1, 4, MTmp);

	/* delta g2 */
	MTmp = DTmp*Mat3x3(MatCross, d2);
	WMA.Add(1, 6 + 4, MTmp);
	WMA.Sub(6 + 1, 6 + 4, MTmp);

	/* Moment equation on node 1 */
	/* d1 x F/d */
	MTmp = d1.Cross(DTmp);

	/* delta x1 */
	WMA.Add(4, 1, MTmp);

	/* delta x2 */
	WMA.Sub(4, 6 + 1, MTmp);

	/* delta g1 */
	WMA.Sub(4, 4, MTmp*Mat3x3(MatCross, d1));

	/* delta g2 */
	WMA.Add(4, 6 + 4, MTmp*Mat3x3(MatCross, d2));

	/* Moment equation on node 2 */
	MTmp = d2.Cross(DTmp);

	/* delta x1 */
	WMA.Sub(6 + 4, 1, MTmp);

	/* delta x2 */
	WMA.Add(6 + 4, 6 + 1, MTmp);

	/* delta g1 */
	WMA.Add(6 + 4, 4, MTmp*Mat3x3(MatCross, d1));

	/* delta g2 */
	WMA.Sub(6 + 4, 6 + 4, MTmp*Mat3x3(MatCross, d2));
}

/* assemblaggio jacobiano - FDEPrime */
void
DeformableDispJoint::AssMatFDEPrime(FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB,
	const Vec3& d1, const Vec3& d2, doublereal dCoef)
{
	/* F/dot{d} */
	WMB.Add(1, 1, FDEPrime);
	WMB.Sub(6 + 1, 1, FDEPrime);
	WMB.Sub(1, 6 + 1, FDEPrime);
	WMB.Add(6 + 1, 6 + 1, FDEPrime);

	Mat3x3 MTmp(d1.Cross(FDEPrime));

	WMB.Add(4, 1, MTmp);
	WMB.Sub(4, 6 + 1, MTmp);

	MTmp = d2.Cross(FDEPrime);

	WMB.Sub(6 + 4, 1, MTmp);
	WMB.Add(6 + 4, 6 + 1, MTmp);

	/* F/dot{d} * [ d2 x ] */
	MTmp = FDEPrime*Mat3x3(MatCross, d2);

	WMB.Add(1, 6 + 4, MTmp);
	WMB.Sub(6 + 1, 6 + 4, MTmp);

	WMB.Add(4, 6 + 4, d1.Cross(MTmp));
	WMB.Sub(6 + 4, 6 + 4, d2.Cross(MTmp));

	/* F/dot{d} * [ d1 x ] */
	MTmp = FDEPrime*Mat3x3(MatCross, d1);

	WMB.Sub(1, 4, MTmp);
	WMB.Add(6 + 1, 4, MTmp);

	WMB.Sub(4, 4, d1.Cross(MTmp));
	WMB.Add(6 + 4, 4, d2.Cross(MTmp));

	/* F/dot{d} * [ ( w1 * dCoef ) x ] */
	Mat3x3 CTmp = FDEPrime*Mat3x3(MatCross, pNode1->GetWCurr()*dCoef);

	WMA.Sub(1, 1, CTmp);
	WMA.Add(6 + 1, 1, CTmp);
	WMA.Add(1, 6 + 1, CTmp);
	WMA.Sub(6 + 1, 6 + 1, CTmp);

	MTmp = d1.Cross(CTmp);

	WMA.Sub(4, 1, MTmp);
	WMA.Add(4, 6 + 1, MTmp);

	MTmp = d2.Cross(CTmp);

	WMA.Add(6 + 4, 1, MTmp);
	WMA.Sub(6 + 4, 6 + 1, MTmp);

	/* F/dot{d} * ( [ d1Prime x ] - [ w1 x ] [ d1 x ] ) * dCoef */
	Vec3 d1Prime(pNode2->GetVCurr()
		+ pNode2->GetWCurr().Cross(d2)
		- pNode1->GetVCurr());
	MTmp = FDEPrime*(Mat3x3(MatCross, d1Prime*dCoef) - Mat3x3(MatCrossCross, pNode1->GetWCurr(), d1*dCoef));
	WMA.Sub(1, 4, MTmp);
	WMA.Add(6 + 1, 4, MTmp);

	WMA.Sub(4, 4, d1.Cross(MTmp));
	WMA.Add(6 + 4, 4, d2.Cross(MTmp));

	/* F/dot{d} * ( [ ( w2 x d2 ) x ] - [ w1 x ] [ d2 x ] ) * dCoef */
	Vec3 d2Prime(pNode2->GetWCurr().Cross(d2));
	MTmp = FDEPrime*(Mat3x3(MatCross, d2Prime*dCoef)
		- Mat3x3(MatCrossCross, pNode1->GetWCurr(), d2*dCoef));
	WMA.Add(1, 6 + 4, MTmp);
	WMA.Sub(6 + 1, 6 + 4, MTmp);

	WMA.Add(4, 6 + 4, d1.Cross(MTmp));
	WMA.Sub(6 + 4, 6 + 4, d2.Cross(MTmp));
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
DeformableDispJoint::AssJac(VariableSubMatrixHandler& WorkMat,
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

/* Jacobian matrix assembly - all but the purely elastic */
void
DeformableDispJoint::AssMats(VariableSubMatrixHandler& WorkMatA,
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

/* Contributo al file di restart */
std::ostream&
DeformableDispJoint::Restart(std::ostream& out) const
{
	Joint::Restart(out) << ", deformable displacement joint, "
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
DeformableDispJoint::OutputPrepare(OutputHandler& OH)
{
	if (bToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			std::string name;
			OutputPrepare_int("deformable displacement", OH, name);

			Var_tilde_d = OH.CreateVar<Vec3>(name + "d",
					OutputHandler::Dimensions::Length,
					"relative position in local frame (x, y, z)");
			Var_tilde_dPrime = OH.CreateVar<Vec3>(name + "dPrime",
					OutputHandler::Dimensions::Velocity,
					"relative linear velocity in local frame (x, y, z)");
			Var_d = OH.CreateVar<Vec3>(name + "D",
					OutputHandler::Dimensions::Length,
					"relative position in global frame (x, y, z)");
			Var_dPrime = OH.CreateVar<Vec3>(name + "DPrime",
					OutputHandler::Dimensions::Velocity,
					"relative linear velocity in global frame (x, y, z)");
		}
#endif // USE_NETCDF
	}
}

void
DeformableDispJoint::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			Joint::NetCDFOutput(OH, GetF(), Zero3, pNode1->GetRCurr()*(tilde_R1h*GetF()), Zero3);
			OH.WriteNcVar(Var_tilde_d, tilde_d);
			OH.WriteNcVar(Var_d, (pNode1->GetRCurr()*(tilde_R1h*tilde_d)));
			if (GetConstLawType() & ConstLawType::VISCOUS) {
				OH.WriteNcVar(Var_tilde_dPrime, tilde_dPrime);
				OH.WriteNcVar(Var_dPrime, (pNode1->GetRCurr()*(tilde_R1h*tilde_dPrime)));
			} else {
				OH.WriteNcVar(Var_tilde_dPrime, Zero3);
				OH.WriteNcVar(Var_dPrime, Zero3);
			}
		}
#endif // USE_NETCDF

		if (OH.UseText(OutputHandler::JOINTS)) {
			Joint::Output(OH.Joints(), "DeformableDispJoint", GetLabel(),
					GetF(), Zero3,
				pNode1->GetRCurr()*(tilde_R1h*GetF()), Zero3)
				<< " " << tilde_d;
			if (GetConstLawType() & ConstLawType::VISCOUS) {
				OH.Joints() << " " << tilde_dPrime;
			}
			OH.Joints() << std::endl;
		}
	}
}

void
DeformableDispJoint::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	if (ph) {
		for (unsigned i = 0; i < ph->size(); i++) {
			Joint::JointHint *pjh = dynamic_cast<Joint::JointHint *>((*ph)[i]);

			if (pjh) {
				if (dynamic_cast<Joint::OffsetHint<1> *>(pjh)) {
					Mat3x3 R1T(pNode1->GetRCurr().Transpose());
					Vec3 f2(pNode2->GetRCurr()*tilde_f2);
  	 
					tilde_f1 = R1T*(pNode2->GetXCurr() + f2 - pNode1->GetXCurr());
	
				} else if (dynamic_cast<Joint::OffsetHint<2> *>(pjh)) {
					Mat3x3 R2T(pNode2->GetRCurr().Transpose());
					Vec3 f1(pNode1->GetRCurr()*tilde_f1);
  	 
					tilde_f2 = R2T*(pNode1->GetXCurr() + f1 - pNode2->GetXCurr());
	
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
			ConstitutiveLaw3DOwner::SetValue(pDM, X, XP, ph);
		}
	}
}

Hint *
DeformableDispJoint::ParseHint(DataManager *pDM, const char *s) const
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

		if (strcmp(&s[1], /* { */ "}") != 0) {
			return 0;
		}

		switch (s[0]) {
		case '1':
			return new Joint::HingeHint<1>;

		case '2':
			return new Joint::HingeHint<2>;
		}
	}

	return ConstitutiveLaw3DOwner::ParseHint(pDM, s);
}

/* inverse dynamics capable element */
bool
DeformableDispJoint::bInverseDynamics(void) const
{
	return true;
}

unsigned int
DeformableDispJoint::iGetNumPrivData(void) const
{
	return 9 + ConstitutiveLaw3DOwner::iGetNumPrivData();
}

unsigned int
DeformableDispJoint::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);

	unsigned idx = 0;

	switch (s[0]) {
	case 'd':
		break;

	case 'v':
		idx += 3;
		break;

	case 'F':
		idx += 6;
		break;

	default:
	{
		size_t l = STRLENOF("constitutiveLaw.");
		if (strncmp(s, "constitutiveLaw.", l) == 0) {
			idx = ConstitutiveLaw3DOwner::iGetPrivDataIdx(&s[l]);
			if (idx > 0) {
				return 9 + idx;
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
DeformableDispJoint::dGetPrivData(unsigned int i) const
{
	ASSERT(i > 0);

	ASSERT(i <= iGetNumPrivData());

	switch (i) {
	case 1:
	case 2:
	case 3:
	{
		/* FIXME: allows simplifications by using only the column
		 * and the components of tilde_R1h that is actually required */
		Vec3 d2(pNode2->GetRCurr()*tilde_f2);
		Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
		Vec3 tilde_d(R1h.MulTV(pNode2->GetXCurr() + d2 - pNode1->GetXCurr()) - tilde_R1hT_tilde_f1);
		return tilde_d(i);
	}

	case 4:
	case 5:
	case 6:
	{
		Vec3 d2(pNode2->GetRCurr()*tilde_f2);
		Mat3x3 R1h(pNode1->GetRCurr());
		Vec3 d1(pNode2->GetXCurr() + d2 - pNode1->GetXCurr());
		Vec3 d1Prime(pNode2->GetVCurr() + pNode2->GetWCurr().Cross(d2) - pNode1->GetVCurr());
		Vec3 tilde_dPrime(R1h.MulTV(d1Prime - pNode1->GetWCurr().Cross(d1)));

		return tilde_dPrime(i - 3);
	}

	case 7:
	case 8:
	case 9:
		return GetF()(i - 6);

	default:
		return ConstitutiveLaw3DOwner::dGetPrivData(i - 9);
	}
}

/* DeformableDispJoint - end */


/* ElasticDispJoint - begin */

ElasticDispJoint::ElasticDispJoint(unsigned int uL,
	const DofOwner* pDO,
	const ConstitutiveLaw3D* pCL,
	const StructNode* pN1,
	const StructNode* pN2,
	const Vec3& tilde_f1,
	const Vec3& tilde_f2,
	const Mat3x3& tilde_R1h,
	const Mat3x3& tilde_R2h,
	flag fOut)
: Elem(uL, fOut),
DeformableDispJoint(uL, pDO, pCL, pN1, pN2, tilde_f1, tilde_f2, tilde_R1h, tilde_R2h, fOut)
{
	/*
	 * Chiede la matrice tangente di riferimento
	 * e la porta nel sistema globale
	 */
	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);
	FDE = R1h*GetFDE().MulMT(R1h);
}

ElasticDispJoint::~ElasticDispJoint(void)
{
	NO_OP;
}

void
ElasticDispJoint::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
	ConstitutiveLaw3DOwner::AfterConvergence(tilde_d);
}

/* assemblaggio jacobiano */
void
ElasticDispJoint::AssMats(VariableSubMatrixHandler& WorkMatA,
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

	// NOTE: the second instance of WMA is ignored (see below)
	AssMats(WMA, WMA, 1.);
}

void
ElasticDispJoint::AssMats(FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& /* WMB */ ,
	doublereal dCoef)
{
	Vec3 d2(pNode2->GetRCurr()*tilde_f2);
	Vec3 d1(pNode2->GetXCurr() + d2 - pNode1->GetXCurr());

	AssMatF(WMA, d1, d2, dCoef);
	AssMatFDE(WMA, d1, d2, dCoef);
}

/* assemblaggio residuo */
SubVectorHandler&
ElasticDispJoint::AssRes(SubVectorHandler& WorkVec,
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

/* Inverse Dynamics Jacobian matrix assembly */
VariableSubMatrixHandler&
ElasticDispJoint::AssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	ASSERT(bIsErgonomy());

	// HACK?  Need to call AfterPredict() here to update MDE and so
	// const_cast because AfterPredict() does not modify X, XP
	AfterPredict(const_cast<VectorHandler&>(XCurr), const_cast<VectorHandler&>(XCurr));

	return DeformableDispJoint::AssJac(WorkMat, 1., XCurr, XCurr);
}

/* Inverse Dynamics Residual Assembly */
SubVectorHandler&
ElasticDispJoint::AssRes(SubVectorHandler& WorkVec,
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

/* Inverse Dynamics update */
void
ElasticDispJoint::Update(const VectorHandler& XCurr, InverseDynamics::Order iOrder)
{
	NO_OP;
}

/* Inverse Dynamics after convergence */
void
ElasticDispJoint::AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP, const VectorHandler& XPP)
{
	ConstitutiveLaw3DOwner::AfterConvergence(tilde_d);
}

void
ElasticDispJoint::AssVec(SubVectorHandler& WorkVec)
{
	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
	Vec3 d2(pNode2->GetRCurr()*tilde_f2);
	Vec3 d1(pNode2->GetXCurr() + d2 - pNode1->GetXCurr());

	if (bFirstRes) {
		bFirstRes = false;

	} else {
		tilde_d = R1h.MulTV(d1) - tilde_R1hT_tilde_f1;

		ConstitutiveLaw3DOwner::Update(tilde_d);
	}

	F = R1h*GetF();

	WorkVec.Add(1, F);
	WorkVec.Add(4, d1.Cross(F));
	WorkVec.Sub(6 + 1, F);
	WorkVec.Sub(6 + 4, d2.Cross(F));
}

void
ElasticDispJoint::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
	Vec3 d2(pNode2->GetRCurr()*tilde_f2);

	tilde_d = R1h.MulTV(pNode2->GetXCurr() + d2 - pNode1->GetXCurr())
		- tilde_R1hT_tilde_f1;

	ConstitutiveLaw3DOwner::Update(tilde_d);

	/* FIXME: we need to be able to regenerate FDE
	 * if the constitutive law throws ChangedEquationStructure */
	FDE = R1h*GetFDE().MulMT(R1h);

	bFirstRes = true;
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
ElasticDispJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
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

/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
ElasticDispJoint::InitialAssRes(SubVectorHandler& WorkVec,
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

/* ElasticDispJoint - end */


/* ElasticDispJointInv - begin */

ElasticDispJointInv::ElasticDispJointInv(unsigned int uL,
	const DofOwner* pDO,
	const ConstitutiveLaw3D* pCL,
	const StructNode* pN1,
	const StructNode* pN2,
	const Vec3& tilde_f1,
	const Vec3& tilde_f2,
	const Mat3x3& tilde_R1h,
	const Mat3x3& tilde_R2h,
	flag fOut)
: Elem(uL, fOut),
DeformableDispJoint(uL, pDO, pCL, pN1, pN2, tilde_f1, tilde_f2, tilde_R1h, tilde_R2h, fOut)
{
	/*
	 * Chiede la matrice tangente di riferimento
	 * e la porta nel sistema globale
	 */
	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);
	FDE = R1h*GetFDE().MulMT(R1h);
}

ElasticDispJointInv::~ElasticDispJointInv(void)
{
	NO_OP;
}

void
ElasticDispJointInv::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
	ConstitutiveLaw3DOwner::AfterConvergence(tilde_d);
}

/* assemblaggio jacobiano */
void
ElasticDispJointInv::AssMats(VariableSubMatrixHandler& WorkMatA,
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

	// NOTE: the second instance of WMA is ignored (see below)
	AssMats(WMA, WMA, 1.);
}

void
ElasticDispJointInv::AssMats(FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& /* WMB */ ,
	doublereal dCoef)
{
	Vec3 d2(pNode2->GetRCurr()*tilde_f2);
	Vec3 d1(pNode2->GetXCurr() + d2 - pNode1->GetXCurr());

	AssMatF(WMA, d1, d2, dCoef);
	AssMatFDE(WMA, d1, d2, dCoef);
}

/* assemblaggio residuo */
SubVectorHandler&
ElasticDispJointInv::AssRes(SubVectorHandler& WorkVec,
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
ElasticDispJointInv::AssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */,
	const VectorHandler& /* XPrimeCurr */, 
	const VectorHandler& /* XPrimePrimeCurr */, 
	InverseDynamics::Order iOrder)
{	
	ASSERT(iOrder == InverseDynamics::INVERSE_DYNAMICS);

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

void
ElasticDispJointInv::AssVec(SubVectorHandler& WorkVec)
{
	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);
	Mat3x3 R2h(pNode2->GetRRef()*tilde_R2h);
	Vec3 ThetaCurr = RotManip::VecRot(R1h.MulTM(R2h));
	Mat3x3 tilde_R = RotManip::Rot(ThetaCurr/2.);
	Mat3x3 hat_R(R1h*tilde_R);

	Mat3x3 hat_I = hat_R*((Eye3 + tilde_R).Inv()).MulMT(hat_R);

	Vec3 f1(pNode1->GetRCurr()*tilde_f1);
	Vec3 f2(pNode2->GetRCurr()*tilde_f2);
	Vec3 d(pNode2->GetXCurr() + f2 - pNode1->GetXCurr() - f1);

	if (bFirstRes) {
		bFirstRes = false;

	} else {
		tilde_d = hat_R.MulTV(d);

		ConstitutiveLaw3DOwner::Update(tilde_d);
	}

	F = hat_R*GetF();

	Vec3 dCrossF(d.Cross(F));
	WorkVec.Add(1, F);
	WorkVec.Add(4, f1.Cross(F) + hat_I*dCrossF);
	WorkVec.Sub(6 + 1, F);
	WorkVec.Sub(6 + 4, f2.Cross(F) - hat_I.MulTV(dCrossF));
}

void
ElasticDispJointInv::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);
	Mat3x3 R2h(pNode2->GetRRef()*tilde_R2h);
	Vec3 ThetaCurr = RotManip::VecRot(R1h.MulTM(R2h));
	Mat3x3 hat_R(R1h*RotManip::Rot(ThetaCurr/2.));

	Vec3 f1(pNode1->GetRCurr()*tilde_f1);
	Vec3 f2(pNode2->GetRCurr()*tilde_f2);
	tilde_d = hat_R.MulTV(pNode2->GetXCurr() + f2 - pNode1->GetXCurr() - f1);

	ConstitutiveLaw3DOwner::Update(tilde_d);

	/* FIXME: we need to be able to regenerate FDE
	 * if the constitutive law throws ChangedEquationStructure */
	FDE = hat_R*GetFDE().MulMT(hat_R);

	bFirstRes = true;
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
ElasticDispJointInv::InitialAssJac(VariableSubMatrixHandler& WorkMat,
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

/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
ElasticDispJointInv::InitialAssRes(SubVectorHandler& WorkVec,
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

/* ElasticDispJointInv - end */


/* ViscousDispJoint - begin */

ViscousDispJoint::ViscousDispJoint(unsigned int uL,
	const DofOwner* pDO,
	const ConstitutiveLaw3D* pCL,
	const StructNode* pN1,
	const StructNode* pN2,
	const Vec3& tilde_f1,
	const Vec3& tilde_f2,
	const Mat3x3& tilde_R1h,
	const Mat3x3& tilde_R2h,
	flag fOut)
: Elem(uL, fOut),
DeformableDispJoint(uL, pDO, pCL, pN1, pN2, tilde_f1, tilde_f2, tilde_R1h, tilde_R2h, fOut)
{
	/*
	 * Chiede la matrice tangente di riferimento
	 * e la porta nel sistema globale
	 */
	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);
	FDEPrime = R1h*GetFDEPrime()*R1h.Transpose();
}

ViscousDispJoint::~ViscousDispJoint(void)
{
	NO_OP;
}

void
ViscousDispJoint::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
	ConstitutiveLaw3DOwner::AfterConvergence(Zero3, tilde_dPrime);
}

void
ViscousDispJoint::AssMats(FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB,
	doublereal dCoef)
{
	Vec3 d2(pNode2->GetRCurr()*tilde_f2);
	Vec3 d1(pNode2->GetXCurr() + d2 - pNode1->GetXCurr());

	AssMatF(WMA, d1, d2, dCoef);
	AssMatFDEPrime(WMA, WMB, d1, d2, dCoef);
}

/* assemblaggio residuo */
SubVectorHandler&
ViscousDispJoint::AssRes(SubVectorHandler& WorkVec,
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
ViscousDispJoint::AssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */,
	const VectorHandler& /* XPrimeCurr */, 
	const VectorHandler& /* XPrimePrimeCurr */, 
	InverseDynamics::Order iOrder)
{	
	ASSERT(iOrder == InverseDynamics::INVERSE_DYNAMICS);

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

/* Inverse Dynamics update */
void
ViscousDispJoint::Update(const VectorHandler& XCurr, InverseDynamics::Order iOrder)
{
	if (iOrder != InverseDynamics::INVERSE_DYNAMICS) {
		// issue "default" warning message
		Joint::Update(XCurr, iOrder);
	}
}

/* Inverse Dynamics after convergence */
void
ViscousDispJoint::AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP, const VectorHandler& XPP)
{
	ConstitutiveLaw3DOwner::AfterConvergence(tilde_d);
}

void
ViscousDispJoint::AssVec(SubVectorHandler& WorkVec)
{
	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
	Vec3 d2(pNode2->GetRCurr()*tilde_f2);
	Vec3 d1(pNode2->GetXCurr() + d2 - pNode1->GetXCurr());

	if (bFirstRes) {
		bFirstRes = false;

	} else {
		Vec3 d1Prime(pNode2->GetVCurr()
			+ pNode2->GetWCurr().Cross(d2)
			- pNode1->GetVCurr());

		tilde_dPrime = R1h.Transpose()*(d1Prime
			- pNode1->GetWCurr().Cross(d1));

		ConstitutiveLaw3DOwner::Update(Zero3, tilde_dPrime);
	}

	Vec3 F(R1h*GetF());

	WorkVec.Add(1, F);
	WorkVec.Add(4, d1.Cross(F));
	WorkVec.Sub(6 + 1, F);
	WorkVec.Sub(6 + 4, d2.Cross(F));
}

void
ViscousDispJoint::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
	Mat3x3 R1hT(R1h.Transpose());
	Vec3 d2(pNode2->GetRCurr()*tilde_f2);
	Vec3 d1(pNode2->GetXCurr() + d2 - pNode1->GetXCurr());
	Vec3 d1Prime(pNode2->GetVCurr()
		+ pNode2->GetWCurr().Cross(d2)
		- pNode1->GetVCurr());

	tilde_dPrime = R1h.Transpose()*(d1Prime
		- pNode1->GetWCurr().Cross(d1));

	ConstitutiveLaw3DOwner::Update(Zero3, tilde_dPrime);

	/* FIXME: we need to be able to regenerate FDE
	 * if the constitutive law throws ChangedEquationStructure */
	FDEPrime = R1h*GetFDEPrime()*R1hT;

	bFirstRes = true;
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
ViscousDispJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
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
	integer iNode1FirstVelIndex = iNode1FirstPosIndex + 6;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstVelIndex = iNode2FirstPosIndex + 6;

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode1FirstVelIndex + iCnt);

		WM.PutRowIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutColIndex(12 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutColIndex(18 + iCnt, iNode2FirstVelIndex + iCnt);
	}

	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);
	Vec3 Omega1(pNode1->GetWRef());
	Vec3 Omega2(pNode2->GetWRef());

	Vec3 F(R1h*GetF());
	Mat3x3 FDEPrime(R1h*GetFDEPrime()*R1h.Transpose());
	Mat3x3 Tmp(Mat3x3(MatCross, F) - FDEPrime*Mat3x3(MatCross, Omega2 - Omega1));

	WM.Add(1, 1, Tmp);
	WM.Sub(4, 1, Tmp);

	WM.Add(1, 4, FDEPrime);
	WM.Add(4, 6 + 4, FDEPrime);

	WM.Sub(1, 6 + 4, FDEPrime);
	WM.Sub(4, 4, FDEPrime);

	return WorkMat;
}

/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
ViscousDispJoint::InitialAssRes(SubVectorHandler& WorkVec,
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

	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
	Vec3 Omega1(pNode1->GetWCurr());
	Vec3 Omega2(pNode2->GetWCurr());

	Vec3 g1(pNode1->GetgCurr());
	Vec3 g2(pNode2->GetgCurr());

	/* Aggiornamento: tilde_d += R1h^T(G(g2)*g2-G(g1)*g1) */
	tilde_d = R1h.Transpose()*(g2*(4./(4. + g2.Dot()))
			- g1*(4./(4. + g1.Dot())));
	tilde_dPrime = R1h.Transpose()*(Omega2 - Omega1);
	ConstitutiveLaw3DOwner::Update(tilde_d, tilde_dPrime);

	Vec3 F(R1h*GetF());

	WorkVec.Add(1, F);
	WorkVec.Sub(4, F);

	return WorkVec;
}

/* ViscousDispJoint - end */


/* ViscoElasticDispJoint - begin */

ViscoElasticDispJoint::ViscoElasticDispJoint(unsigned int uL,
	const DofOwner* pDO,
	const ConstitutiveLaw3D* pCL,
	const StructNode* pN1,
	const StructNode* pN2,
	const Vec3& tilde_f1,
	const Vec3& tilde_f2,
	const Mat3x3& tilde_R1h,
	const Mat3x3& tilde_R2h,
	flag fOut)
: Elem(uL, fOut),
DeformableDispJoint(uL, pDO, pCL, pN1, pN2, tilde_f1, tilde_f2, tilde_R1h, tilde_R2h, fOut)
{
	/*
	 * Chiede la matrice tangente di riferimento
	 * e la porta nel sistema globale
	 */
	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);
	Mat3x3 R1hT(R1h.Transpose());
	FDE = R1h*GetFDE()*R1hT;
	FDEPrime = R1h*GetFDEPrime()*R1hT;
}

ViscoElasticDispJoint::~ViscoElasticDispJoint(void)
{
	NO_OP;
}

void
ViscoElasticDispJoint::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
	ConstitutiveLaw3DOwner::AfterConvergence(tilde_d, tilde_dPrime);
}

/* assemblaggio jacobiano */
void
ViscoElasticDispJoint::AssMats(FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB,
	doublereal dCoef)
{
	Vec3 d2(pNode2->GetRCurr()*tilde_f2);
	Vec3 d1(pNode2->GetXCurr() + d2 - pNode1->GetXCurr());

	AssMatF(WMA, d1, d2, dCoef);
	AssMatFDE(WMA, d1, d2, dCoef);
	AssMatFDEPrime(WMA, WMB, d1, d2, dCoef);
}

/* assemblaggio residuo */
SubVectorHandler&
ViscoElasticDispJoint::AssRes(SubVectorHandler& WorkVec,
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
ViscoElasticDispJoint::AssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */,
	const VectorHandler& /* XPrimeCurr */, 
	const VectorHandler& /* XPrimePrimeCurr */, 
	InverseDynamics::Order iOrder)
{	
	ASSERT(iOrder == InverseDynamics::INVERSE_DYNAMICS);

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

/* Inverse Dynamics update */
void
ViscoElasticDispJoint::Update(const VectorHandler& XCurr, InverseDynamics::Order iOrder)
{
	if (iOrder != InverseDynamics::INVERSE_DYNAMICS) {
		// issue "default" warning message
		Joint::Update(XCurr, iOrder);
	}
}

/* Inverse Dynamics after convergence */
void
ViscoElasticDispJoint::AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP, const VectorHandler& XPP)
{
	ConstitutiveLaw3DOwner::AfterConvergence(tilde_d);
}

void
ViscoElasticDispJoint::AssVec(SubVectorHandler& WorkVec)
{
	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
	Vec3 d2(pNode2->GetRCurr()*tilde_f2);
	Vec3 d1(pNode2->GetXCurr() + d2 - pNode1->GetXCurr());

	if (bFirstRes) {
		bFirstRes = false;

	} else {
		Mat3x3 R1hT(R1h.Transpose());
		Vec3 d1Prime(pNode2->GetVCurr()
			+ pNode2->GetWCurr().Cross(d2)
			- pNode1->GetVCurr());

		tilde_d = R1hT*d1 - tilde_R1hT_tilde_f1;
		tilde_dPrime = R1hT*(d1Prime - pNode1->GetWCurr().Cross(d1));

		ConstitutiveLaw3DOwner::Update(tilde_d, tilde_dPrime);
	}

	Vec3 F(R1h*GetF());

	WorkVec.Add(1, F);
	WorkVec.Add(4, d1.Cross(F));
	WorkVec.Sub(6 + 1, F);
	WorkVec.Sub(6 + 4, d2.Cross(F));
}

void
ViscoElasticDispJoint::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
	Mat3x3 R1hT(R1h.Transpose());
	Vec3 d2(pNode2->GetRCurr()*tilde_f2);
	Vec3 d1(pNode2->GetXCurr() + d2 - pNode1->GetXCurr());
	Vec3 d1Prime(pNode2->GetVCurr()
		+ pNode2->GetWCurr().Cross(d2)
		- pNode1->GetVCurr());
	
	tilde_d = R1hT*d1 - tilde_R1hT_tilde_f1;
	tilde_dPrime = R1hT*(d1Prime - pNode1->GetWCurr().Cross(d1));

	ConstitutiveLaw3DOwner::Update(tilde_d, tilde_dPrime);

	/* FIXME: we need to be able to regenerate FDE and FDEPrime
	 * if the constitutive law throws ChangedEquationStructure */
	FDE = R1h*GetFDE()*R1hT;
	FDEPrime = R1h*GetFDEPrime()*R1hT;

	bFirstRes = true;
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
ViscoElasticDispJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
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
	integer iNode1FirstVelIndex = iNode1FirstPosIndex + 6;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstVelIndex = iNode2FirstPosIndex + 6;

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode1FirstVelIndex + iCnt);

		WM.PutRowIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutColIndex(12 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutColIndex(18 + iCnt, iNode2FirstVelIndex + iCnt);
	}

	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);
	Vec3 Omega1(pNode1->GetWRef());
	Vec3 Omega2(pNode2->GetWRef());

	Vec3 F(R1h*GetF());
	Mat3x3 FDE(R1h*GetFDE().MulMT(R1h));
	Mat3x3 FDEPrime(R1h*GetFDEPrime().MulMT(R1h));

	Mat3x3 Tmp(Mat3x3(MatCross, F) - FDEPrime*Mat3x3(MatCross, Omega2 - Omega1) + FDE);

	// FIXME: check and rewrite

	WM.Add(1, 1, Tmp);
	WM.Sub(4, 1, Tmp);

	WM.Add(4, 6 + 1, FDE);
	WM.Sub(1, 6 + 1, FDE);

	WM.Add(1, 4, FDEPrime);
	WM.Add(4, 6 + 4, FDEPrime);

	WM.Sub(1, 6 + 4, FDEPrime);
	WM.Sub(4, 4, FDEPrime);

	return WorkMat;
}

/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
ViscoElasticDispJoint::InitialAssRes(SubVectorHandler& WorkVec,
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

	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
	Vec3 d2(pNode2->GetRCurr()*tilde_f2);
	Vec3 d1(pNode2->GetXCurr() + d2 - pNode1->GetXCurr());

	Mat3x3 R1hT(R1h.Transpose());
	Vec3 d1Prime(pNode2->GetVCurr()
		+ pNode2->GetWCurr().Cross(d2)
		- pNode1->GetVCurr());

	tilde_d = R1hT*d1 - tilde_R1hT_tilde_f1;
	tilde_dPrime = R1hT*(d1Prime - pNode1->GetWCurr().Cross(d1));

	ConstitutiveLaw3DOwner::Update(tilde_d, tilde_dPrime);

	Vec3 F(R1h*GetF());

	WorkVec.Add(1, F);
	WorkVec.Add(4, d1.Cross(F));
	WorkVec.Sub(6 + 1, F);
	WorkVec.Sub(6 + 4, d2.Cross(F));

	return WorkVec;
}

/* ViscoElasticDispJoint - end */
