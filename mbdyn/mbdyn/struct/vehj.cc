/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2006
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

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "dataman.h"
#include "vehj.h"

#include "matvecexp.h"
#include "Rot.hh"

/* define to use invariant residual */
#undef MBDYN_DEFORMABLE_INVARIANT

/* DeformableHingeJoint - begin */

/* Costruttore non banale */
DeformableHingeJoint::DeformableHingeJoint(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw3D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Mat3x3& tilde_R1h,
		const Mat3x3& tilde_R2h,
		flag fOut)
: Elem(uL, fOut),
Joint(uL, pDO, fOut),
ConstitutiveLaw3DOwner(pCL),
pNode1(pN1), pNode2(pN2), tilde_R1h(tilde_R1h), tilde_R2h(tilde_R2h), bFirstRes(true)
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

/* Contributo al file di restart */
std::ostream&
DeformableHingeJoint::Restart(std::ostream& out) const
{
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
		Mat3x3 R1 = pNode1->GetRCurr()*tilde_R1h;
		Mat3x3 R2 = pNode2->GetRCurr()*tilde_R2h;

		Vec3 d(MatR2EulerAngles(R1.Transpose()*R2)*dRaDegr);
		Vec3 v(GetF());
		Joint::Output(OH.Joints(), "DeformableHinge", GetLabel(),
			Zero3, v, Zero3, R1*v)
			<< " " << d << std::endl;
	}
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

Hint *
DeformableHingeJoint::ParseHint(DataManager *pDM, const char *s) const
{
	if (strncasecmp(s, "hinge{" /* } */, sizeof("hinge{" /* } */) - 1) == 0) {
		s += sizeof("hinge{" /* } */) - 1;

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
		size_t l = sizeof("constitutiveLaw.") - 1;
		if (strncmp(s, "constitutiveLaw.", l) == 0) {
			return 9 + ConstitutiveLaw3DOwner::iGetPrivDataIdx(s + l);
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
		Mat3x3 R1hT((pNode1->GetRCurr()*tilde_R1h).Transpose());
		Mat3x3 R2h(pNode2->GetRCurr()*tilde_R2h);

		Vec3 v(RotManip::VecRot(R1hT*R2h));

		return v(i);
	}

	case 4:
	case 5:
	case 6:
	{
		Mat3x3 R1hT((pNode1->GetRCurr()*tilde_R1h).Transpose());
		Vec3 w = R1hT*(pNode2->GetWCurr() - pNode1->GetWCurr());

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

/* DeformableHingeJoint - end */


/* ElasticHingeJoint - begin */

ElasticHingeJoint::ElasticHingeJoint(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw3D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Mat3x3& tilde_R1h,
		const Mat3x3& tilde_R2h,
		flag fOut)
: Elem(uL, fOut),
DeformableHingeJoint(uL, pDO, pCL, pN1, pN2, tilde_R1h, tilde_R2h, fOut),
tilde_ThetaRef(0.), FDE(0.)
{
	/*
	 * Chiede la matrice tangente di riferimento
	 * e la porta nel sistema globale
	 */
	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);
	FDE = R1h*GetFDE()*R1h.Transpose();
}

ElasticHingeJoint::~ElasticHingeJoint(void)
{
	NO_OP;
}

void
ElasticHingeJoint::AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP)
{
	ConstitutiveLaw3DOwner::AfterConvergence(tilde_ThetaCurr);
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
ElasticHingeJoint::AfterPredict(VectorHandler& /* X */ ,
		VectorHandler& /* XP */ )
{
	/* Calcola le deformazioni, aggiorna il legame costitutivo
	 * e crea la FDE */

	/* Recupera i dati */
	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);
	Mat3x3 R2h(pNode2->GetRRef()*tilde_R2h);
	Mat3x3 R1hT(R1h.Transpose());

	/* Calcola la deformazione corrente nel sistema locale (nodo a) */
	tilde_ThetaCurr = tilde_ThetaRef = RotManip::VecRot(R1hT*R2h);

	/* Calcola l'inversa di Gamma di tilde_ThetaRef */
	Mat3x3 GammaRefm1 = RotManip::DRot_I(tilde_ThetaRef);

	/* Aggiorna il legame costitutivo */
	ConstitutiveLaw3DOwner::Update(tilde_ThetaRef);

	/* Chiede la matrice tangente di riferimento e la porta
	 * nel sistema globale */
	FDE = R1h*ConstitutiveLaw3DOwner::GetFDE()*GammaRefm1*R1hT;

	bFirstRes = true;
}

void
ElasticHingeJoint::AssMat(FullSubMatrixHandler& WM, doublereal dCoef)
{
	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);

	Mat3x3 FDETmp = FDE*dCoef;

	WM.Add(4, 4, FDETmp);
	WM.Sub(1, 4, FDETmp);

	FDETmp += Mat3x3(R1h*(GetF()*dCoef));
	WM.Add(1, 1, FDETmp);
	WM.Sub(4, 1, FDETmp);
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

void
ElasticHingeJoint::AssVec(SubVectorHandler& WorkVec)
{
	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);

	if (bFirstRes) {
		bFirstRes = false;

	} else {
		Mat3x3 R2h(pNode2->GetRCurr()*tilde_R2h);
		tilde_ThetaCurr = RotManip::VecRot(R1h.Transpose()*R2h);
		ConstitutiveLaw3DOwner::Update(tilde_ThetaCurr);
	}

#ifndef MBDYN_DEFORMABLE_INVARIANT
	/* Couple attached to node 1 */
	Vec3 M(R1h*GetF());
#else /* MBDYN_DEFORMABLE_INVARIANT */
	/* Couple in intermediate reference frame (invariant) */
	Vec3 M(R1h*(RotManip::Rot(tilde_ThetaCurr/2.)*GetF()));
#endif /* MBDYN_DEFORMABLE_INVARIANT */

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


/* ViscousHingeJoint - begin */

ViscousHingeJoint::ViscousHingeJoint(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw3D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Mat3x3& tilde_R1h,
		const Mat3x3& tilde_R2h,
		flag fOut)
: Elem(uL, fOut),
DeformableHingeJoint(uL, pDO, pCL, pN1, pN2, tilde_R1h, tilde_R2h, fOut)
{
	Mat3x3 R1h(pNode1->GetRRef());
	FDEPrime = R1h*ConstitutiveLaw3DOwner::GetFDEPrime()*R1h.Transpose();
}

ViscousHingeJoint::~ViscousHingeJoint(void)
{
	NO_OP;
}

void
ViscousHingeJoint::AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP)
{
	ConstitutiveLaw3DOwner::AfterConvergence(Zero3, tilde_Omega);
}

void
ViscousHingeJoint::AfterPredict(VectorHandler& /* X */ ,
		VectorHandler& /* XP */ )
{
	/* Calcola le deformazioni, aggiorna il legame costitutivo
	 * e crea la FDE */

	/* Recupera i dati */
	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);
	Mat3x3 R1hT(R1h.Transpose());

	/* Aggiorna il legame costitutivo */
	tilde_Omega = R1hT*(pNode2->GetWRef() - pNode1->GetWRef());
	ConstitutiveLaw3DOwner::Update(Zero3, tilde_Omega);

	/* Chiede la matrice tangente di riferimento e la porta
	 * nel sistema globale */
	FDEPrime = R1h*GetFDEPrime()*R1hT;

	bFirstRes = true;
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
	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);
	Vec3 W2(pNode1->GetWRef());

	WMB.Add(4, 4, FDEPrime);
	WMB.Sub(1, 4, FDEPrime);
	WMB.Add(1, 1, FDEPrime);
	WMB.Sub(4, 1, FDEPrime);

	Mat3x3 Tmp(FDEPrime*Mat3x3(W2*dCoef));
	WMA.Sub(4, 4, Tmp);
	WMA.Add(1, 4, Tmp);

	Tmp -= Mat3x3(R1h*(ConstitutiveLaw3DOwner::GetF()*dCoef));
	WMA.Sub(1, 1, Tmp);
	WMA.Add(4, 1, Tmp);
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

	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);

	if (bFirstRes) {
		bFirstRes = false;

	} else {
		Mat3x3 R1hT(R1h.Transpose());

		tilde_Omega = R1hT*(pNode2->GetWCurr() - pNode1->GetWCurr());

		ConstitutiveLaw3DOwner::Update(Zero3, tilde_Omega);
	}

#ifndef MBDYN_DEFORMABLE_INVARIANT
	Vec3 M(R1h*ConstitutiveLaw3DOwner::GetF());
#else /* MBDYN_DEFORMABLE_INVARIANT */
	Vec3 M(R1h*(RotManip::Rot(RotManip::VecRot(R1h.Transpose()*pNode2->GetRCurr()*tilde_R2h)/2.)*GetF()));
#endif /* MBDYN_DEFORMABLE_INVARIANT */

	WorkVec.Add(1, M);
	WorkVec.Sub(4, M);

	return WorkVec;
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

	Vec3 W2(pNode2->GetWCurr());

	FDEPrime = R1h*ConstitutiveLaw3DOwner::GetFDEPrime()*R1h.Transpose();

	Mat3x3 Tmp(FDEPrime*Mat3x3(W2));
	WM.Add(4, 7, Tmp);
	WM.Sub(1, 7, Tmp);

	Tmp += Mat3x3(R1h*ConstitutiveLaw3DOwner::GetF());
	WM.Add(1, 1, Tmp);
	WM.Sub(4, 1, Tmp);

	WM.Add(1, 4, FDEPrime);
	WM.Add(4, 10, FDEPrime);

	WM.Sub(1, 10, FDEPrime);
	WM.Sub(4, 4, FDEPrime);

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
		Mat3x3 R1hT(R1h.Transpose());

		tilde_Omega = R1hT*(pNode2->GetWCurr() - pNode1->GetWCurr());

		ConstitutiveLaw3DOwner::Update(0., tilde_Omega);
	}

	Vec3 M(R1h*ConstitutiveLaw3DOwner::GetF());

	WorkVec.Add(1, M);
	WorkVec.Sub(4, M);

	return WorkVec;
}

/* ViscousHingeJoint - end */


/* ViscoElasticHingeJoint - begin */

ViscoElasticHingeJoint::ViscoElasticHingeJoint(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw3D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Mat3x3& tilde_R1h,
		const Mat3x3& tilde_R2h,
		flag fOut)
: Elem(uL, fOut),
DeformableHingeJoint(uL, pDO, pCL, pN1, pN2, tilde_R1h, tilde_R2h, fOut),
tilde_ThetaRef(0.)
{
	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);
	Mat3x3 R1hT(R1h.Transpose());
	FDE = R1h*ConstitutiveLaw3DOwner::GetFDE()*R1hT;
	FDEPrime = R1h*ConstitutiveLaw3DOwner::GetFDEPrime()*R1hT;
}

ViscoElasticHingeJoint::~ViscoElasticHingeJoint(void)
{
	NO_OP;
}

void
ViscoElasticHingeJoint::AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP)
{
	ConstitutiveLaw3DOwner::AfterConvergence(tilde_ThetaCurr, tilde_Omega);
}

void
ViscoElasticHingeJoint::AfterPredict(VectorHandler& /* X */ ,
		VectorHandler& /* XP */ )
{
	/* Calcola le deformazioni, aggiorna il legame costitutivo
	 * e crea la FDE */

	/* Recupera i dati */
	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);
	Mat3x3 R2h(pNode2->GetRRef()*tilde_R2h);
	Mat3x3 R1hT(R1h.Transpose());

	/* Calcola la deformazione corrente nel sistema locale (nodo a) */
	tilde_ThetaCurr = tilde_ThetaRef = RotManip::VecRot(R1hT*R2h);

	/* Calcola l'inversa di Gamma di tilde_ThetaRef */
	Mat3x3 GammaRefm1 = RotManip::DRot_I(tilde_ThetaRef);

	/* Aggiorna il legame costitutivo */
	tilde_Omega = R1hT*(pNode2->GetWRef() - pNode1->GetWRef());

	/* Aggiorna il legame costitutivo */
	ConstitutiveLaw3DOwner::Update(tilde_ThetaRef, tilde_Omega);

	/* Chiede la matrice tangente di riferimento e la porta
	 * nel sistema globale */
	FDE = R1h*ConstitutiveLaw3DOwner::GetFDE()*GammaRefm1*R1hT;
	FDEPrime = R1h*GetFDEPrime()*R1hT;

	bFirstRes = true;
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
	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);
	Vec3 W2(pNode1->GetWRef());

	WMB.Add(1, 1, FDEPrime);
	WMB.Sub(1, 4, FDEPrime);
	WMB.Sub(4, 1, FDEPrime);
	WMB.Add(4, 4, FDEPrime);

	Mat3x3 Tmp(FDE*dCoef);
	Tmp -= FDEPrime*Mat3x3(W2*dCoef);
	WMA.Add(4, 4, Tmp);
	WMA.Sub(1, 4, Tmp);

	Tmp += Mat3x3(R1h*(ConstitutiveLaw3DOwner::GetF()*dCoef));
	WMA.Add(1, 1, Tmp);
	WMA.Sub(4, 1, Tmp);
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

	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);

	if (bFirstRes) {
		bFirstRes = false;

	} else {
		Mat3x3 R1hT(R1h.Transpose());
		Mat3x3 R2h(pNode2->GetRCurr()*tilde_R2h);

		tilde_ThetaCurr = RotManip::VecRot(R1hT*R2h);
		tilde_Omega = R1hT*(pNode2->GetWCurr() - pNode1->GetWCurr());

		ConstitutiveLaw3DOwner::Update(tilde_ThetaCurr, tilde_Omega);
	}

#ifndef MBDYN_DEFORMABLE_INVARIANT
	Vec3 M(R1h*ConstitutiveLaw3DOwner::GetF());
#else /* MBDYN_DEFORMABLE_INVARIANT */
	Vec3 M(R1h*(RotManip::Rot(tilde_ThetaCurr/2.)*GetF()));
#endif /* MBDYN_DEFORMABLE_INVARIANT */

	WorkVec.Add(1, M);
	WorkVec.Sub(4, M);

	return WorkVec;
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
	Mat3x3 R1hT(R1h.Transpose());
	Vec3 W2(pNode2->GetWCurr());

	FDE = R1h*ConstitutiveLaw3DOwner::GetFDE()*R1hT;
	FDEPrime = R1h*ConstitutiveLaw3DOwner::GetFDEPrime()*R1hT;

	Mat3x3 Tmp(FDE + FDEPrime*Mat3x3(W2));
	WM.Add(4, 7, Tmp);
	WM.Sub(1, 7, Tmp);

	Tmp += Mat3x3(R1h*ConstitutiveLaw3DOwner::GetF());
	WM.Add(1, 1, Tmp);
	WM.Sub(4, 1, Tmp);

	WM.Add(1, 4, FDEPrime);
	WM.Add(4, 10, FDEPrime);

	WM.Sub(1, 10, FDEPrime);
	WM.Sub(4, 4, FDEPrime);

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
		Mat3x3 R1hT(R1h.Transpose());
		Mat3x3 R2h(pNode2->GetRCurr()*tilde_R2h);

		tilde_ThetaCurr = RotManip::VecRot(R1hT*R2h);
		tilde_Omega = R1hT*(pNode2->GetWCurr() - pNode1->GetWCurr());

		ConstitutiveLaw3DOwner::Update(tilde_ThetaCurr, tilde_Omega);
	}

	Vec3 F(R1h*ConstitutiveLaw3DOwner::GetF());

	WorkVec.Add(1, F);
	WorkVec.Sub(4, F);

	return WorkVec;
}

/* ViscoElasticHingeJoint - end */

