/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2011
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
#include "vehj.h"

#include "matvecexp.h"
#include "Rot.hh"

/* DeformableHingeJoint - begin */

/* Costruttore non banale */
DeformableHingeJoint::DeformableHingeJoint(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw3D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Mat3x3& tilde_R1h,
		const Mat3x3& tilde_R2h,
		const OrientationDescription& od,
		flag fOut)
: Elem(uL, fOut),
Joint(uL, pDO, fOut),
ConstitutiveLaw3DOwner(pCL),
pNode1(pN1),
pNode2(pN2),
tilde_R1h(tilde_R1h),
tilde_R2h(tilde_R2h),
od(od),
bFirstRes(false)
{
	ASSERT(pNode1 != NULL);
	ASSERT(pNode2 != NULL);
	ASSERT(pNode1->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pNode2->GetNodeType() == Node::STRUCTURAL);
}

/* Distruttore */
DeformableHingeJoint::~DeformableHingeJoint(void)
{
	NO_OP;
}

/* assemblaggio jacobiano */
/* Used by all "attached" hinges */
void
DeformableHingeJoint::AssMatM(FullSubMatrixHandler& WMA,
	doublereal dCoef)
{
	/* M was updated by AssRes */
	Mat3x3 MTmp(MatCross, M*dCoef);

	WMA.Add(1, 1, MTmp);
	WMA.Sub(4, 1, MTmp);
}

/* Used by all "invariant" hinges */
void
DeformableHingeJoint::AssMatMInv(FullSubMatrixHandler& WMA, doublereal dCoef)
{
	/* M was updated by AssRes;
	 * hat_I and hat_IT were updated by AfterPredict */
	Vec3 MCoef(M*dCoef);
	Mat3x3 MTmp(MCoef.Cross(hat_IT));

	WMA.Add(1, 1, MTmp);
	WMA.Sub(4, 1, MTmp);

	MTmp = MCoef.Cross(hat_I);

	WMA.Add(1, 4, MTmp);
	WMA.Sub(4, 4, MTmp);
}

/* Used by all hinges */
void
DeformableHingeJoint::AssMatMDE(FullSubMatrixHandler& WMA,
	doublereal dCoef)
{
	Mat3x3 MTmp(MDE*dCoef);

	WMA.Add(1, 1, MTmp);
	WMA.Sub(1, 4, MTmp);
	WMA.Sub(4, 1, MTmp);
	WMA.Add(4, 4, MTmp);
}

/* Used by all "attached" hinges */
void
DeformableHingeJoint::AssMatMDEPrime(FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB, doublereal dCoef)
{
	WMB.Add(1, 1, MDEPrime);
	WMB.Sub(1, 4, MDEPrime);
	WMB.Sub(4, 1, MDEPrime);
	WMB.Add(4, 4, MDEPrime);

	Mat3x3 MTmp(MDEPrime*Mat3x3(MatCross, pNode2->GetWCurr()*dCoef));
	WMA.Sub(1, 1, MTmp);
	WMA.Add(1, 4, MTmp);
	WMA.Add(4, 1, MTmp);
	WMA.Sub(4, 4, MTmp);
}

/* Used by all "invariant" hinges */
void
DeformableHingeJoint::AssMatMDEPrimeInv(FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB, doublereal dCoef)
{
	/* hat_I and hat_IT were updated by AfterPredict */
	WMB.Add(1, 1, MDEPrime);
	WMB.Sub(1, 4, MDEPrime);
	WMB.Sub(4, 1, MDEPrime);
	WMB.Add(4, 4, MDEPrime);

	Vec3 W1(pNode1->GetWCurr()*dCoef);
	Vec3 W2(pNode2->GetWCurr()*dCoef);
	Mat3x3 MTmp(MDEPrime*(W2.Cross(hat_IT) + W1.Cross(hat_I)));
	WMA.Sub(1, 1, MTmp);
	WMA.Add(1, 4, MTmp);
	WMA.Add(4, 1, MTmp);
	WMA.Sub(4, 4, MTmp);
}

/* Contributo al file di restart */
std::ostream&
DeformableHingeJoint::Restart(std::ostream& out) const
{
	/* FIXME: does not work for invariant hinge */
	Joint::Restart(out) << ", deformable hinge, "
		<< pNode1->GetLabel() << ", hinge, reference, node, 1, ",
		(tilde_R1h.GetVec(1)).Write(out, ", ")
		<< ", 2, ", (tilde_R1h.GetVec(2)).Write(out, ", ") << ", "
		<< pNode2->GetLabel() << ", hinge, reference, node, 1, ",
		(tilde_R2h.GetVec(1)).Write(out, ", ")
		<< ", 2, ", (tilde_R2h.GetVec(2)).Write(out, ", ") << ", ";
	return pGetConstLaw()->Restart(out) << ';' << std::endl;
}

void
DeformableHingeJoint::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
		Mat3x3 R2h(pNode2->GetRCurr()*tilde_R2h);
		Mat3x3 R(R1h.MulTM(R2h));

		Vec3 v(GetF());
		Joint::Output(OH.Joints(), "DeformableHinge", GetLabel(),
			Zero3, v, Zero3, R1h*v) << " ";

		switch (od) {
		case EULER_123:
			OH.Joints() << MatR2EulerAngles123(R)*dRaDegr;
			break;

		case EULER_313:
			OH.Joints() << MatR2EulerAngles313(R)*dRaDegr;
			break;

		case EULER_321:
			OH.Joints() << MatR2EulerAngles321(R)*dRaDegr;
			break;

		case ORIENTATION_VECTOR:
			OH.Joints() << RotManip::VecRot(R);
			break;

		case ORIENTATION_MATRIX:
			OH.Joints() << R;
			break;

		default:
			/* impossible */
			break;
		}

		if (GetConstLawType() & ConstLawType::VISCOUS) {
			OH.Joints() << " " << R1h.MulTV(pNode2->GetWCurr() - pNode1->GetWCurr());
		}

		OH.Joints() << std::endl;
	}
}

void
DeformableHingeJoint::OutputInv(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
		Mat3x3 R2h(pNode2->GetRCurr()*tilde_R2h);
		Mat3x3 R(R1h.MulTM(R2h));
		Mat3x3 hat_R(R1h*RotManip::Rot(RotManip::VecRot(R)/2.));

		Vec3 v(GetF());
		Joint::Output(OH.Joints(), "DeformableHinge", GetLabel(),
			Zero3, v, Zero3, hat_R*v) << " ";

		switch (od) {
		case EULER_123:
			OH.Joints() << MatR2EulerAngles123(R)*dRaDegr;
			break;

		case EULER_313:
			OH.Joints() << MatR2EulerAngles313(R)*dRaDegr;
			break;

		case EULER_321:
			OH.Joints() << MatR2EulerAngles321(R)*dRaDegr;
			break;

		case ORIENTATION_VECTOR:
			OH.Joints() << RotManip::VecRot(R);
			break;

		case ORIENTATION_MATRIX:
			OH.Joints() << R;
			break;

		default:
			/* impossible */
			break;
		}

		if (GetConstLawType() & ConstLawType::VISCOUS) {
			OH.Joints() << " " << hat_R.MulTV(pNode2->GetWCurr() - pNode1->GetWCurr());
			}
		OH.Joints() << std::endl;
	}
}

void
DeformableHingeJoint::AfterPredict(VectorHandler& /* X */ ,
		VectorHandler& /* XP */ )
{
	bFirstRes = true;

	AfterPredict();
}

void
DeformableHingeJoint::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
	if (ph) {
		for (unsigned i = 0; i < ph->size(); i++) {
			Joint::JointHint *pjh = dynamic_cast<Joint::JointHint *>((*ph)[i]);

			if (pjh) {
				if (dynamic_cast<Joint::HingeHint<1> *>(pjh)) {
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

/* inverse dynamics capable element */
bool
DeformableHingeJoint::bInverseDynamics(void) const
{
	return true;
}

void
DeformableHingeJoint::SetInitialValue(VectorHandler& /* X */ )
{
	AfterPredict();
}

Hint *
DeformableHingeJoint::ParseHint(DataManager *pDM, const char *s) const
{
	if (strncasecmp(s, "hinge{" /*}*/, STRLENOF("hinge{" /*}*/)) == 0) {
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

	return ConstitutiveLaw3DOwner::ParseHint(pDM, s);
}

unsigned int
DeformableHingeJoint::iGetNumPrivData(void) const
{
	return 9 + ConstitutiveLaw3DOwner::iGetNumPrivData();
}

unsigned int
DeformableHingeJoint::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);

	unsigned idx = 0;

	switch (s[0]) {
	case 'r':
		break;

	case 'w':
		idx += 3;
		break;

	case 'M':
		idx += 6;
		break;

	default:
	{
		const size_t l = STRLENOF("constitutiveLaw.");
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
DeformableHingeJoint::dGetPrivData(unsigned int i) const
{
	ASSERT(i > 0);

	ASSERT(i <= iGetNumPrivData());

	switch (i) {
	case 1:
	case 2:
	case 3:
	{
		/* NOTE: this is correct also in the invariant case */
		Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
		Mat3x3 R2h(pNode2->GetRCurr()*tilde_R2h);

		Vec3 v(RotManip::VecRot(R1h.MulTM(R2h)));

		return v(i);
	}

	case 4:
	case 5:
	case 6:
	{
		Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
		Vec3 w(R1h.MulTV(pNode2->GetWCurr() - pNode1->GetWCurr()));

		return w(i - 3);
	}

	case 7:
	case 8:
	case 9:
		return GetF()(i - 6);

	default:
		return ConstitutiveLaw3DOwner::dGetPrivData(i - 9);
	}
}

doublereal
DeformableHingeJoint::dGetPrivDataInv(unsigned int i) const
{
	ASSERT(i > 0);

	ASSERT(i <= iGetNumPrivData());

	switch (i) {
	case 4:
	case 5:
	case 6:
	{
		Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
		Mat3x3 R2h(pNode2->GetRCurr()*tilde_R2h);
		Mat3x3 hat_R(R1h*RotManip::Rot(RotManip::VecRot(R1h.MulTM(R2h)/2.)));
		Vec3 w(hat_R.MulTV(pNode2->GetWCurr() - pNode1->GetWCurr()));

		return w(i - 3);
	}

	default:
		/* fall back to regular one */
		return DeformableHingeJoint::dGetPrivData(i);
	}
}

/* DeformableHingeJoint - end */


/* ElasticHingeJoint - begin */

ElasticHingeJoint::ElasticHingeJoint(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw3D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Mat3x3& tilde_R1h,
		const Mat3x3& tilde_R2h,
		const OrientationDescription& od,
		flag fOut)
: Elem(uL, fOut),
DeformableHingeJoint(uL, pDO, pCL, pN1, pN2, tilde_R1h, tilde_R2h, od, fOut),
ThetaRef(Zero3)
{
	NO_OP;
}

ElasticHingeJoint::~ElasticHingeJoint(void)
{
	NO_OP;
}

void
ElasticHingeJoint::AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP)
{
	ConstitutiveLaw3DOwner::AfterConvergence(ThetaCurr);
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
ElasticHingeJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering ElasticHingeJoint::AssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex() + 3;
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex() + 3;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex() + 3;
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex() + 3;

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(3 + iCnt, iNode2FirstMomIndex + iCnt);
		WM.PutColIndex(3 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	AssMat(WM, dCoef);

	return WorkMat;
}

/* assemblaggio jacobiano */
void
ElasticHingeJoint::AssMats(VariableSubMatrixHandler& WorkMatA,
		VariableSubMatrixHandler& WorkMatB,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering ElasticHingeJoint::AssJac()" << std::endl);

	FullSubMatrixHandler& WMA = WorkMatA.SetFull();
	WorkMatB.SetNullMatrix();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WMA.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex() + 3;
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex() + 3;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex() + 3;
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex() + 3;

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WMA.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WMA.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WMA.PutRowIndex(3 + iCnt, iNode2FirstMomIndex + iCnt);
		WMA.PutColIndex(3 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	AssMat(WMA, 1.);
}

void
ElasticHingeJoint::AfterPredict(void)
{
	/* Calcola le deformazioni, aggiorna il legame costitutivo
	 * e crea la MDE */

	/* Recupera i dati */
	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);
	Mat3x3 R2h(pNode2->GetRRef()*tilde_R2h);

	/* Calcola la deformazione corrente nel sistema locale (nodo a) */
	ThetaCurr = ThetaRef = RotManip::VecRot(R1h.MulTM(R2h));

	/* Aggiorna il legame costitutivo */
	ConstitutiveLaw3DOwner::Update(ThetaRef);

	/* Calcola l'inversa di Gamma di ThetaRef */
	Mat3x3 GammaRefm1 = RotManip::DRot_I(ThetaRef);

	/* Chiede la matrice tangente di riferimento e la porta
	 * nel sistema globale */
	MDE = R1h*ConstitutiveLaw3DOwner::GetFDE()*GammaRefm1.MulMT(R1h);
}

void
ElasticHingeJoint::AssMat(FullSubMatrixHandler& WMA, doublereal dCoef)
{
#if 0
	// Calcola l'inversa di Gamma di ThetaRef
	Mat3x3 GammaCurrm1 = RotManip::DRot_I(ThetaCurr);

	// Chiede la matrice tangente di riferimento e la porta
	// nel sistema globale
	MDE = R1h*ConstitutiveLaw3DOwner::GetFDE()*GammaCurrm1.MulMT(R1h);
#endif

	AssMatM(WMA, dCoef);
	AssMatMDE(WMA, dCoef);
}

/* assemblaggio residuo */
SubVectorHandler&
ElasticHingeJoint::AssRes(SubVectorHandler& WorkVec,
		doublereal /* dCoef */ ,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering ElasticHingeJoint::AssRes()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera gli indici */
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex() + 3;
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex() + 3;

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WorkVec.PutRowIndex(3 + iCnt, iNode2FirstMomIndex + iCnt);
	}

	AssVec(WorkVec);

	return WorkVec;
}

/* Inverse Dynamics Residual Assembly */
SubVectorHandler&
ElasticHingeJoint::AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& /* XCurr */,
		const VectorHandler& /* XPrimeCurr */, 
		const VectorHandler& /* XPrimePrimeCurr */, 
		InverseDynamics::Order iOrder)
{
	DEBUGCOUT("Entering ElasticHingeJoint::AssRes()" << std::endl);

	ASSERT(iOrder == InverseDynamics::INVERSE_DYNAMICS);
	
	/* There is no need to call AfterPredict, everything is done in AssVec*/
	bFirstRes = false;

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera gli indici */
	integer iNode1FirstMomIndex = pNode1->iGetFirstPositionIndex() + 3;
	integer iNode2FirstMomIndex = pNode2->iGetFirstPositionIndex() + 3;

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WorkVec.PutRowIndex(3 + iCnt, iNode2FirstMomIndex + iCnt);
	}

	AssVec(WorkVec);

	return WorkVec;
}


void
ElasticHingeJoint::AssVec(SubVectorHandler& WorkVec)
{
	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);

	if (bFirstRes) {
		bFirstRes = false;

	} else {
		Mat3x3 R2h(pNode2->GetRCurr()*tilde_R2h);
		ThetaCurr = RotManip::VecRot(R1h.MulTM(R2h));
		ConstitutiveLaw3DOwner::Update(ThetaCurr);
	}

	/* Couple attached to node 1 */
	M = R1h*GetF();

	WorkVec.Add(1, M);
	WorkVec.Sub(4, M);
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
ElasticHingeJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& /* XCurr */ )
{
	DEBUGCOUT("Entering ElasticHingeJoint::InitialAssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex() + 3;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex() + 3;

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(3 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutColIndex(3 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	AssMat(WM, 1.);

	return WorkMat;
}

/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
ElasticHingeJoint::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& /* XCurr */ )
{
	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex() + 3;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex() + 3;

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WorkVec.PutRowIndex(3 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	AssVec(WorkVec);

	return WorkVec;
}

/* ElasticHingeJoint - end */


/* ElasticHingeJointInv - begin */

void
ElasticHingeJointInv::AfterPredict(void)
{
	/* Recupera i dati */
	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);
	Mat3x3 R2h(pNode2->GetRRef()*tilde_R2h);
	ThetaCurr = ThetaRef = RotManip::VecRot(R1h.MulTM(R2h));

	/* Aggiorna il legame costitutivo */
	ConstitutiveLaw3DOwner::Update(ThetaRef);

	Mat3x3 tilde_R(RotManip::Rot(ThetaRef/2.));
	Mat3x3 hat_R(R1h*tilde_R);

	/* Calcola l'inversa di Gamma di ThetaRef */
	Mat3x3 GammaRefm1 = RotManip::DRot_I(ThetaRef);

	/* Chiede la matrice tangente di riferimento e la porta
	 * nel sistema globale */
	MDE = hat_R*ConstitutiveLaw3DOwner::GetFDE()*GammaRefm1.MulMT(R1h);

	hat_I = hat_R*(Eye3 + tilde_R).Inv().MulMT(hat_R);
	hat_IT = hat_I.Transpose();
}

ElasticHingeJointInv::ElasticHingeJointInv(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw3D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Mat3x3& tilde_R1h,
		const Mat3x3& tilde_R2h,
		const OrientationDescription& od,
		flag fOut)
: Elem(uL, fOut),
ElasticHingeJoint(uL, pDO, pCL, pN1, pN2, tilde_R1h, tilde_R2h, od, fOut)
{
	NO_OP;
}

ElasticHingeJointInv::~ElasticHingeJointInv(void)
{
	NO_OP;
}

void
ElasticHingeJointInv::AssMatM(FullSubMatrixHandler& WMA, doublereal dCoef)
{
	DeformableHingeJoint::AssMatMInv(WMA, dCoef);
}

void
ElasticHingeJointInv::AssVec(SubVectorHandler& WorkVec)
{
	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);

	if (bFirstRes) {
		/* ThetaCurr and the constitutive laws were updated
		 * by AfterPredict */
		bFirstRes = false;

	} else {
		Mat3x3 R2h(pNode2->GetRCurr()*tilde_R2h);
		ThetaCurr = RotManip::VecRot(R1h.MulTM(R2h));

		/* aggiorna il legame costitutivo */
		ConstitutiveLaw3DOwner::Update(ThetaCurr);
	}

	Mat3x3 hat_R(R1h*RotManip::Rot(ThetaCurr/2.));
	M = hat_R*GetF();

	WorkVec.Add(1, M);
	WorkVec.Sub(4, M);
}

void
ElasticHingeJointInv::Output(OutputHandler& OH) const
{
	DeformableHingeJoint::OutputInv(OH);
}

doublereal
ElasticHingeJointInv::dGetPrivData(unsigned int i) const
{
	return DeformableHingeJoint::dGetPrivDataInv(i);
}

/* ElasticHingeJointInv - end */


/* ViscousHingeJoint - begin */

ViscousHingeJoint::ViscousHingeJoint(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw3D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Mat3x3& tilde_R1h,
		const Mat3x3& tilde_R2h,
		const OrientationDescription& od,
		flag fOut)
: Elem(uL, fOut),
DeformableHingeJoint(uL, pDO, pCL, pN1, pN2, tilde_R1h, tilde_R2h, od, fOut)
{
	NO_OP;
}

ViscousHingeJoint::~ViscousHingeJoint(void)
{
	NO_OP;
}

void
ViscousHingeJoint::AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP)
{
	ConstitutiveLaw3DOwner::AfterConvergence(Zero3, Omega);
}

void
ViscousHingeJoint::AfterPredict(void)
{
	/* Calcola le deformazioni, aggiorna il legame costitutivo
	 * e crea la MDE */

	/* Recupera i dati */
	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);

	/* Aggiorna il legame costitutivo */
	Omega = R1h.MulTV(pNode2->GetWRef() - pNode1->GetWRef());
	ConstitutiveLaw3DOwner::Update(Zero3, Omega);

	/* Chiede la matrice tangente di riferimento e la porta
	 * nel sistema globale */
	MDEPrime = R1h*GetFDEPrime().MulMT(R1h);
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
ViscousHingeJoint::AssJac(VariableSubMatrixHandler& WorkMat,
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
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex() + 3;
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex() + 3;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex() + 3;
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex() + 3;

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(3 + iCnt, iNode2FirstMomIndex + iCnt);
		WM.PutColIndex(3 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	AssMats(WM, WM, dCoef);

	return WorkMat;
}

/* assemblaggio jacobiano */
void
ViscousHingeJoint::AssMats(VariableSubMatrixHandler& WorkMatA,
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
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex() + 3;
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex() + 3;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex() + 3;
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex() + 3;

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WMA.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WMA.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WMA.PutRowIndex(3 + iCnt, iNode2FirstMomIndex + iCnt);
		WMA.PutColIndex(3 + iCnt, iNode2FirstPosIndex + iCnt);

		WMB.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WMB.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WMB.PutRowIndex(3 + iCnt, iNode2FirstMomIndex + iCnt);
		WMB.PutColIndex(3 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	AssMats(WMA, WMA, 1.);
}

/* assemblaggio jacobiano */
void
ViscousHingeJoint::AssMats(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB,
		doublereal dCoef)
{
	AssMatM(WMA, dCoef);
	AssMatMDEPrime(WMA, WMB, dCoef);
}

/* assemblaggio residuo */
SubVectorHandler&
ViscousHingeJoint::AssRes(SubVectorHandler& WorkVec,
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
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex() + 3;
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex() + 3;

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WorkVec.PutRowIndex(3 + iCnt, iNode2FirstMomIndex + iCnt);
	}

	AssVec(WorkVec);

	return WorkVec;
}

/* Inverse Dynamics Residual Assembly */
SubVectorHandler&
ViscousHingeJoint::AssRes(SubVectorHandler& WorkVec,
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
	integer iNode1FirstMomIndex = pNode1->iGetFirstPositionIndex() + 3;
	integer iNode2FirstMomIndex = pNode2->iGetFirstPositionIndex() + 3;

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WorkVec.PutRowIndex(3 + iCnt, iNode2FirstMomIndex + iCnt);
	}

	AssVec(WorkVec);

	return WorkVec;
}

/* assemblaggio residuo */
void
ViscousHingeJoint::AssVec(SubVectorHandler& WorkVec)
{
	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);

	if (bFirstRes) {
		bFirstRes = false;

	} else {
		Omega = R1h.MulTV(pNode2->GetWCurr() - pNode1->GetWCurr());

		ConstitutiveLaw3DOwner::Update(Zero3, Omega);
	}

	M = R1h*ConstitutiveLaw3DOwner::GetF();

	WorkVec.Add(1, M);
	WorkVec.Sub(4, M);
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
ViscousHingeJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& /* XCurr */ )
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex() + 3;
	integer iNode1FirstVelIndex = iNode1FirstPosIndex + 6;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex() + 3;
	integer iNode2FirstVelIndex = iNode2FirstPosIndex + 6;

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(3 + iCnt, iNode1FirstVelIndex + iCnt);

		WM.PutRowIndex(3 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutColIndex(9 + iCnt, iNode2FirstVelIndex + iCnt);
	}

	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);

	const Vec3& W2(pNode2->GetWCurr());

	MDEPrime = R1h*ConstitutiveLaw3DOwner::GetFDEPrime().MulMT(R1h);

	Mat3x3 Tmp(MDEPrime*Mat3x3(MatCross, W2));
	WM.Add(4, 7, Tmp);
	WM.Sub(1, 7, Tmp);

	Tmp += Mat3x3(MatCross, R1h*ConstitutiveLaw3DOwner::GetF());
	WM.Add(1, 1, Tmp);
	WM.Sub(4, 1, Tmp);

	WM.Add(1, 4, MDEPrime);
	WM.Add(4, 10, MDEPrime);

	WM.Sub(1, 10, MDEPrime);
	WM.Sub(4, 4, MDEPrime);

	return WorkMat;
}

/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
ViscousHingeJoint::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& /* XCurr */ )
{
	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex() + 3;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex() + 3;

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WorkVec.PutRowIndex(3 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);

	if (bFirstRes) {
		bFirstRes = false;

	} else {
		Omega = R1h.MulTV(pNode2->GetWCurr() - pNode1->GetWCurr());

		ConstitutiveLaw3DOwner::Update(Zero3, Omega);
	}

	Vec3 M(R1h*ConstitutiveLaw3DOwner::GetF());

	WorkVec.Add(1, M);
	WorkVec.Sub(4, M);

	return WorkVec;
}

/* ViscousHingeJoint - end */


/* ViscousHingeJointInv - begin */

ViscousHingeJointInv::ViscousHingeJointInv(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw3D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Mat3x3& tilde_R1h,
		const Mat3x3& tilde_R2h,
		const OrientationDescription& od,
		flag fOut)
: Elem(uL, fOut),
ViscousHingeJoint(uL, pDO, pCL, pN1, pN2, tilde_R1h, tilde_R2h, od, fOut)
{
	NO_OP;
}

ViscousHingeJointInv::~ViscousHingeJointInv(void)
{
	NO_OP;
}

void
ViscousHingeJointInv::AfterPredict(void)
{
	/* Calcola le deformazioni, aggiorna il legame costitutivo
	 * e crea la MDE */

	/* Recupera i dati */
	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);
	Mat3x3 R2h(pNode2->GetRRef()*tilde_R2h);
	Vec3 tilde_Theta(RotManip::VecRot(R1h.MulTM(R2h))/2.);
	Mat3x3 tilde_R(RotManip::Rot(tilde_Theta));
	Mat3x3 hat_R(R1h*tilde_R);

	/* Aggiorna il legame costitutivo */
	Omega = hat_R.MulTV(pNode2->GetWRef() - pNode1->GetWRef());
	ConstitutiveLaw3DOwner::Update(Zero3, Omega);

	/* Chiede la matrice tangente di riferimento e la porta
	 * nel sistema globale */
	/* FIXME: Jacobian matrix not implemented yet */
	MDEPrime = hat_R*GetFDEPrime().MulMT(hat_R);

	hat_I = hat_R*(Eye3 + tilde_R).Inv().MulMT(hat_R);
	hat_IT = hat_I.Transpose();
}

void
ViscousHingeJointInv::AssMatM(FullSubMatrixHandler& WMA, doublereal dCoef)
{
	DeformableHingeJoint::AssMatMInv(WMA, dCoef);
}

void
ViscousHingeJointInv::AssMatMDEPrime(FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB, doublereal dCoef)
{
	DeformableHingeJoint::AssMatMDEPrimeInv(WMA, WMB, dCoef);
}

/* assemblaggio residuo */
void
ViscousHingeJointInv::AssVec(SubVectorHandler& WorkVec)
{
	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
	Mat3x3 R2h(pNode2->GetRCurr()*tilde_R2h);
	Vec3 tilde_Theta(RotManip::VecRot(R1h.MulTM(R2h))/2.);
	Mat3x3 hat_R(R1h*RotManip::Rot(tilde_Theta));

	if (bFirstRes) {
		bFirstRes = false;

	} else {
		/* velocita' relativa nel riferimento intermedio */
		Omega = hat_R.MulTV(pNode2->GetWCurr() - pNode1->GetWCurr());

		/* aggiorna la legge costitutiva */
		ConstitutiveLaw3DOwner::Update(Zero3, Omega);
	}

	M = hat_R*GetF();

	WorkVec.Add(1, M);
	WorkVec.Sub(4, M);
}

void
ViscousHingeJointInv::Output(OutputHandler& OH) const
{
	DeformableHingeJoint::OutputInv(OH);
}

doublereal
ViscousHingeJointInv::dGetPrivData(unsigned int i) const
{
	return DeformableHingeJoint::dGetPrivDataInv(i);
}

/* ViscousHingeJointInv - end */


/* ViscoElasticHingeJoint - begin */

ViscoElasticHingeJoint::ViscoElasticHingeJoint(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw3D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Mat3x3& tilde_R1h,
		const Mat3x3& tilde_R2h,
		const OrientationDescription& od,
		flag fOut)
: Elem(uL, fOut),
DeformableHingeJoint(uL, pDO, pCL, pN1, pN2, tilde_R1h, tilde_R2h, od, fOut),
ThetaRef(Zero3)
{
	NO_OP;
}

ViscoElasticHingeJoint::~ViscoElasticHingeJoint(void)
{
	NO_OP;
}

void
ViscoElasticHingeJoint::AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP)
{
	ConstitutiveLaw3DOwner::AfterConvergence(ThetaCurr, Omega);
}

void
ViscoElasticHingeJoint::AfterPredict(void)
{
	/* Computes strains, updates constitutive law and generates
	 * MDE and MDEPrime */

	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);
	Mat3x3 R2h(pNode2->GetRRef()*tilde_R2h);

	/* Current strain in material reference frame (node 1) */
	ThetaCurr = ThetaRef = RotManip::VecRot(R1h.MulTM(R2h));

	/* Relative angular velocity */
	Omega = R1h.MulTV(pNode2->GetWRef() - pNode1->GetWRef());

	/* Updates constitutive law */
	ConstitutiveLaw3DOwner::Update(ThetaRef, Omega);

	/* don't repeat the above operations during AssRes */
	bFirstRes = true;

	/* Inverse of Gamma(ThetaRef) */
	Mat3x3 GammaRefm1 = RotManip::DRot_I(ThetaRef);

	/* Tangent matrices are updated and projected in the global
	 * reference frame; they won't change during the solution
	 * of the current time step, according to the updated-updated
	 * approach */
	MDE = R1h*ConstitutiveLaw3DOwner::GetFDE()*GammaRefm1.MulMT(R1h);
	MDEPrime = R1h*GetFDEPrime().MulMT(R1h);
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
ViscoElasticHingeJoint::AssJac(VariableSubMatrixHandler& WorkMat,
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
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex() + 3;
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex() + 3;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex() + 3;
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex() + 3;

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(3 + iCnt, iNode2FirstMomIndex + iCnt);
		WM.PutColIndex(3 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	AssMats(WM, WM, dCoef);

	return WorkMat;
}

/* assemblaggio jacobiano */
void
ViscoElasticHingeJoint::AssMats(VariableSubMatrixHandler& WorkMatA,
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
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex() + 3;
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex() + 3;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex() + 3;
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex() + 3;

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WMA.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WMA.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WMA.PutRowIndex(3 + iCnt, iNode2FirstMomIndex + iCnt);
		WMA.PutColIndex(3 + iCnt, iNode2FirstPosIndex + iCnt);

		WMB.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WMB.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WMB.PutRowIndex(3 + iCnt, iNode2FirstMomIndex + iCnt);
		WMB.PutColIndex(3 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	AssMats(WMA, WMB, 1.);
}

/* assemblaggio jacobiano */
void
ViscoElasticHingeJoint::AssMats(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB,
		doublereal dCoef)
{
	AssMatM(WMA, dCoef);
	AssMatMDE(WMA, dCoef);
	AssMatMDEPrime(WMA, WMB, dCoef);
}

/* assemblaggio residuo */
SubVectorHandler&
ViscoElasticHingeJoint::AssRes(SubVectorHandler& WorkVec,
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
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex() + 3;
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex() + 3;

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WorkVec.PutRowIndex(3 + iCnt, iNode2FirstMomIndex + iCnt);
	}

	AssVec(WorkVec);

	return WorkVec;
}

/* Inverse Dynamics Residual Assembly */
SubVectorHandler&
ViscoElasticHingeJoint::AssRes(SubVectorHandler& WorkVec,
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
	integer iNode1FirstMomIndex = pNode1->iGetFirstPositionIndex() + 3;
	integer iNode2FirstMomIndex = pNode2->iGetFirstPositionIndex() + 3;

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WorkVec.PutRowIndex(3 + iCnt, iNode2FirstMomIndex + iCnt);
	}

	AssVec(WorkVec);

	return WorkVec;
}

/* assemblaggio residuo */
void
ViscoElasticHingeJoint::AssVec(SubVectorHandler& WorkVec)
{
	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);

	if (bFirstRes) {
		bFirstRes = false;

	} else {
		Mat3x3 R2h(pNode2->GetRCurr()*tilde_R2h);

		/* orientazione intermedia */
		ThetaCurr = RotManip::VecRot(R1h.MulTM(R2h));

		/* velocita' relativa nel riferimento intermedio */
		Omega = R1h.MulTV(pNode2->GetWCurr() - pNode1->GetWCurr());

		/* aggiorna il legame costitutivo */
		ConstitutiveLaw3DOwner::Update(ThetaCurr, Omega);
	}

	M = R1h*ConstitutiveLaw3DOwner::GetF();

	WorkVec.Add(1, M);
	WorkVec.Sub(4, M);
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
ViscoElasticHingeJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& /* XCurr */ )
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex() + 3;
	integer iNode1FirstVelIndex = iNode1FirstPosIndex + 6;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex() + 3;
	integer iNode2FirstVelIndex = iNode2FirstPosIndex + 6;

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(3 + iCnt, iNode1FirstVelIndex + iCnt);

		WM.PutRowIndex(3 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutColIndex(9 + iCnt, iNode2FirstVelIndex + iCnt);
	}

	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
	const Vec3& W2(pNode2->GetWCurr());

	MDE = R1h*ConstitutiveLaw3DOwner::GetFDE().MulMT(R1h);
	MDEPrime = R1h*ConstitutiveLaw3DOwner::GetFDEPrime().MulMT(R1h);

	Mat3x3 Tmp(MDE + MDEPrime*Mat3x3(MatCross, W2));
	WM.Add(4, 7, Tmp);
	WM.Sub(1, 7, Tmp);

	Tmp += Mat3x3(MatCross, R1h*ConstitutiveLaw3DOwner::GetF());
	WM.Add(1, 1, Tmp);
	WM.Sub(4, 1, Tmp);

	WM.Add(1, 4, MDEPrime);
	WM.Add(4, 10, MDEPrime);

	WM.Sub(1, 10, MDEPrime);
	WM.Sub(4, 4, MDEPrime);

	return WorkMat;
}

/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
ViscoElasticHingeJoint::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& /* XCurr */ )
{
	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex() + 3;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex() + 3;

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WorkVec.PutRowIndex(3 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);

	if (bFirstRes) {
		bFirstRes = false;

	} else {
		Mat3x3 R2h(pNode2->GetRCurr()*tilde_R2h);

		ThetaCurr = RotManip::VecRot(R1h.MulTM(R2h));
		Omega = R1h.MulTV(pNode2->GetWCurr() - pNode1->GetWCurr());

		ConstitutiveLaw3DOwner::Update(ThetaCurr, Omega);
	}

	M = R1h*ConstitutiveLaw3DOwner::GetF();

	WorkVec.Add(1, M);
	WorkVec.Sub(4, M);

	return WorkVec;
}

/* ViscoElasticHingeJoint - end */

/* ViscoElasticHingeJointInv - begin */

ViscoElasticHingeJointInv::ViscoElasticHingeJointInv(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw3D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Mat3x3& tilde_R1h,
		const Mat3x3& tilde_R2h,
		const OrientationDescription& od,
		flag fOut)
: Elem(uL, fOut),
ViscoElasticHingeJoint(uL, pDO, pCL, pN1, pN2, tilde_R1h, tilde_R2h, od, fOut)
{
	NO_OP;
}

ViscoElasticHingeJointInv::~ViscoElasticHingeJointInv(void)
{
	NO_OP;
}

void
ViscoElasticHingeJointInv::AfterPredict(void)
{
	/* Calcola le deformazioni, aggiorna il legame costitutivo
	 * e crea le MDE e MDEPrime */

	/* Recupera i dati */
	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);
	Mat3x3 R2h(pNode2->GetRRef()*tilde_R2h);
	ThetaCurr = ThetaRef = RotManip::VecRot(R1h.MulTM(R2h));
	Mat3x3 tilde_R(RotManip::Rot(ThetaRef/2.));
	Mat3x3 hat_R(R1h*tilde_R);

	/* velocita' relativa nel sistema intermedio */
	Omega = hat_R.MulTV(pNode2->GetWRef() - pNode1->GetWRef());

	/* Aggiorna il legame costitutivo */
	ConstitutiveLaw3DOwner::Update(ThetaRef, Omega);

	/* Chiede la matrice tangente di riferimento e la porta
	 * nel sistema globale */

	/* Calcola l'inversa di Gamma di ThetaRef */
	Mat3x3 GammaRefm1 = RotManip::DRot_I(ThetaRef);

	/* Chiede la matrice tangente di riferimento e la porta
	 * nel sistema globale */
	MDE = hat_R*ConstitutiveLaw3DOwner::GetFDE()*GammaRefm1.MulMT(R1h);
	MDEPrime = hat_R*ConstitutiveLaw3DOwner::GetFDEPrime().MulMT(hat_R);

	hat_I = hat_R*(Eye3 + tilde_R).Inv().MulMT(hat_R);
	hat_IT = hat_I.Transpose();
}

/* NOTE: duplicate of ElasticHingeJointInv and ViscousHingeJointInv */
void
ViscoElasticHingeJointInv::AssMatM(FullSubMatrixHandler& WMA, doublereal dCoef)
{
	DeformableHingeJoint::AssMatMInv(WMA, dCoef);
}

/* NOTE: duplicate of ViscousHingeJointInv */
void
ViscoElasticHingeJointInv::AssMatMDEPrime(FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB, doublereal dCoef)
{
	DeformableHingeJoint::AssMatMDEPrimeInv(WMA, WMB, dCoef);
}

/* assemblaggio residuo */
void
ViscoElasticHingeJointInv::AssVec(SubVectorHandler& WorkVec)
{
	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
	Mat3x3 hat_R;

	if (bFirstRes) {
		bFirstRes = false;
		hat_R = R1h*RotManip::Rot(ThetaCurr/2.);

	} else {
		Mat3x3 R2h(pNode2->GetRCurr()*tilde_R2h);
		ThetaCurr = RotManip::VecRot(R1h.MulTM(R2h));

		/* orientazione intermedia */
		hat_R = R1h*RotManip::Rot(ThetaCurr/2.);

		/* velocita' relativa nell'orientazione intermedia */
		Omega = hat_R.MulTV(pNode2->GetWCurr() - pNode1->GetWCurr());

		/* aggiorna il legame costitutivo */
		ConstitutiveLaw3DOwner::Update(ThetaCurr, Omega);
	}

	M = hat_R*GetF();

	WorkVec.Add(1, M);
	WorkVec.Sub(4, M);
}

void
ViscoElasticHingeJointInv::Output(OutputHandler& OH) const
{
	DeformableHingeJoint::OutputInv(OH);
}

doublereal
ViscoElasticHingeJointInv::dGetPrivData(unsigned int i) const
{
	return DeformableHingeJoint::dGetPrivDataInv(i);
}

/* ViscoElasticHingeJointInv - end */


/* InvAngularConstitutiveLaw - begin */

class InvAngularConstitutiveLaw
: public ConstitutiveLaw<Vec3, Mat3x3> {
private:
	doublereal dXi;
	ConstitutiveLaw<Vec3, Mat3x3> *pCL;

public:
	InvAngularConstitutiveLaw(const doublereal& dxi,
		ConstitutiveLaw<Vec3, Mat3x3> *pcl);
	virtual ~InvAngularConstitutiveLaw(void);

	ConstLawType::Type GetConstLawType(void) const;

	virtual ConstitutiveLaw<Vec3, Mat3x3>* pCopy(void) const;
	virtual std::ostream& Restart(std::ostream& out) const;

	virtual void Update(const Vec3& Eps, const Vec3& /* EpsPrime */  = Zero3);
};


InvAngularConstitutiveLaw::InvAngularConstitutiveLaw(const doublereal& dxi,
	ConstitutiveLaw<Vec3, Mat3x3> *pcl)
: dXi(dxi), pCL(pcl)
{
	NO_OP;
}

InvAngularConstitutiveLaw::~InvAngularConstitutiveLaw(void)
{
	if (pCL) {
		delete pCL;
	}
}

ConstLawType::Type
InvAngularConstitutiveLaw::GetConstLawType(void) const
{
	return ConstLawType::Type(ConstLawType::ELASTIC | pCL->GetConstLawType());
}

ConstitutiveLaw<Vec3, Mat3x3>*
InvAngularConstitutiveLaw::pCopy(void) const
{
	ConstitutiveLaw<Vec3, Mat3x3>* pcl = NULL;

	typedef InvAngularConstitutiveLaw cl;
	SAFENEWWITHCONSTRUCTOR(pcl, cl, cl(dXi, pCL->pCopy()));
	return pcl;
}

std::ostream&
InvAngularConstitutiveLaw::Restart(std::ostream& out) const
{
	out << "invariant, " << dXi << ", ";
	return pCL->Restart(out);
}

void
InvAngularConstitutiveLaw::Update(const Vec3& Eps, const Vec3& EpsPrime)
{
	// Save Eps, just in case
	ConstitutiveLaw<Vec3, Mat3x3>::Epsilon = Eps;

	// Theta_xi = Theta * xi
	Vec3 Tx(Eps*dXi);

	// R_xi = exp(xi*Theta x)
	Mat3x3 Rx(RotManip::Rot(Tx));

	// If the underlying CL needs EpsPrime, re-orient it accordingly
	Vec3 EP;
	if (pCL->GetConstLawType() & ConstLawType::VISCOUS) {
		EP = Rx.MulTV(EpsPrime);
	}

	// EP is passed uninitialized if the underlying CL has no VISCOUS
	pCL->Update(Eps, EP);

	// F = R_xi * K * Theta
	ConstitutiveLaw<Vec3, Mat3x3>::F =  Rx*pCL->GetF();

	// Gamma_xi
	Mat3x3 Gx(RotManip::DRot(Tx));

	// re-orientation of moment
	ConstitutiveLaw<Vec3, Mat3x3>::FDE = Mat3x3(MatCross, ConstitutiveLaw<Vec3, Mat3x3>::F*(-dXi))*Gx;

	// stiffness, if any
	if (pCL->GetConstLawType() & ConstLawType::ELASTIC) {
		ConstitutiveLaw<Vec3, Mat3x3>::FDE += Rx*pCL->GetFDE();
	}

	// viscosity, if any
	if (pCL->GetConstLawType() & ConstLawType::VISCOUS) {
		// NOTE: RDE will be incomplete; it'll miss the
		// delta hatR^T * (w_b - w_a) term in case of VISCOUS
		ConstitutiveLaw<Vec3, Mat3x3>::FDEPrime = Rx*pCL->GetFDEPrime().MulMT(Rx);
	}
}

/* InvAngularConstitutiveLaw - end */


/* InvAngularCLR - begin */

ConstitutiveLaw<Vec3, Mat3x3> *
InvAngularCLR::Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType)
{
	ConstitutiveLaw<Vec3, Mat3x3>* pCL = 0;

	doublereal dXi = HP.GetReal();

	ConstitutiveLaw<Vec3, Mat3x3>* pCL2 = HP.GetConstLaw3D(CLType);
	CLType = ConstLawType::Type(ConstLawType::ELASTIC | CLType);

	typedef InvAngularConstitutiveLaw L;
	SAFENEWWITHCONSTRUCTOR(pCL, L, L(dXi, pCL2));

	return pCL;
}

/* InvAngularCLR - end */
