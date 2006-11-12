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
#include "vehj3.h"

#include "matvecexp.h"
#include "Rot.hh"

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
		flag fOut)
: Elem(uL, fOut),
Joint(uL, pDO, fOut),
ConstitutiveLaw6DOwner(pCL),
pNode1(pN1), pNode2(pN2),
tilde_f1(tilde_f1), tilde_f2(tilde_f2),
tilde_R1h(tilde_R1h), tilde_R2h(tilde_R2h),
tilde_k(0.), tilde_kPrime(0.), bFirstRes(true)
{
	ASSERT(pNode1 != NULL);
	ASSERT(pNode2 != NULL);
	ASSERT(pNode1->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pNode2->GetNodeType() == Node::STRUCTURAL);
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
DeformableJoint::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
      		Vec3 F(GetF().GetVec1());
		Vec3 M(GetF().GetVec2());
		Joint::Output(OH.Joints(), "DeformableJoint", GetLabel(),
				F, M, R1h*F, R1h*M) << std::endl;
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

Hint *
DeformableJoint::ParseHint(DataManager *pDM, const char *s) const
{
	if (strncasecmp(s, "offset{" /* } */ , sizeof("offset{" /* } */ ) - 1) == 0)
	{
		s += sizeof("offset{" /* } */ ) - 1;

		if (strcmp(&s[1], /* { */ "}") != 0) {
			return 0;
		}

		switch (s[0]) {
		case '1':
			return new Joint::OffsetHint<1>;

		case '2':
			return new Joint::OffsetHint<2>;
		}

	} else if (strncasecmp(s, "hinge{" /* } */, sizeof("hinge{" /* } */) - 1) == 0) {
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
		size_t l = sizeof("constitutiveLaw.") - 1;
		if (strncmp(s, "constitutiveLaw.", l) == 0) {
			return 18 + ConstitutiveLaw6DOwner::iGetPrivDataIdx(s + l);
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
		flag fOut)
: Elem(uL, fOut),
DeformableJoint(uL, pDO, pCL, pN1, pN2, tilde_f1, tilde_f2, tilde_R1h, tilde_R2h, fOut),
ThetaRef(0.)
{
	/*
	 * Chiede la matrice tangente di riferimento
	 * e la porta nel sistema globale
	 */
	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);
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
VariableSubMatrixHandler&
ElasticJoint::AssJac(VariableSubMatrixHandler& WorkMat,
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
		WM.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
		WM.PutRowIndex(6 + iCnt, iNode2FirstMomIndex+iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex+iCnt);
	}

	AssMat(WM, dCoef);

	return WorkMat;
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
		WMA.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
		WMA.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
		WMA.PutRowIndex(6 + iCnt, iNode2FirstMomIndex+iCnt);
		WMA.PutColIndex(6 + iCnt, iNode2FirstPosIndex+iCnt);
	}

	AssMat(WMA, 1.);
}

void
ElasticJoint::AssMat(FullSubMatrixHandler& WM, doublereal dCoef)
{
	Vec3 f1(pNode1->GetRRef()*tilde_f1);
	Vec3 f2(pNode2->GetRRef()*tilde_f2);
	Vec3 d1(pNode2->GetXCurr() + f2 - pNode1->GetXCurr());

	/* D11 */
	Mat3x3 DTmp = FDE.GetMat11()*dCoef;

	WM.Add(1, 1, DTmp);
	WM.Sub(1, 6 + 1, DTmp);
	WM.Sub(6 + 1, 1, DTmp);
	WM.Add(6 + 1, 6 + 1, DTmp);

	/* D11 * [d1 x] */
	Mat3x3 FTmp = DTmp*Mat3x3(d1);

	WM.Sub(1, 4, FTmp);
	WM.Add(6 + 1, 4, FTmp);

	/* D11 * [f2 x] */
	FTmp = DTmp*Mat3x3(f2);

	WM.Add(1, 6 + 4, FTmp);
	WM.Sub(6 + 1, 6 + 4, FTmp);

	/* [f1 x] * D11 */
	FTmp = Mat3x3(f1)*DTmp;

	WM.Add(4, 1, FTmp);
	WM.Sub(4, 6 + 1, FTmp);

	WM.Sub(4, 4, FTmp*Mat3x3(d1));
	WM.Add(4, 6 + 4, FTmp*Mat3x3(f2));

	/* [f2 x] * D11 */
	FTmp = Mat3x3(f2)*DTmp;

	WM.Sub(6 + 4, 1, FTmp);
	WM.Add(6 + 4, 6 + 1, FTmp);

	WM.Add(6 + 4, 4, FTmp*Mat3x3(d1));
	WM.Sub(6 + 4, 6 + 4, FTmp*Mat3x3(f2));


	/* D12 */
	DTmp = FDE.GetMat12()*dCoef;

	WM.Sub(1, 6 + 4, DTmp);
	WM.Add(6 + 1, 6 + 4, DTmp);

	/* [f1 x] * D12 */
	WM.Sub(4, 6 + 4, Mat3x3(f1)*DTmp);

	/* [f2 x] * D12 */
	WM.Add(6 + 4, 6 + 4, Mat3x3(f2)*DTmp);


	/* D12 + [F x] */
	DTmp += Mat3x3(F.GetVec1()*dCoef);

	WM.Add(1, 4, DTmp);
	WM.Sub(6 + 1, 4, DTmp);

	/* [f1 x] * D12 */
	WM.Add(4, 4, Mat3x3(f1)*DTmp);

	/* [f2 x] * D12 */
	WM.Sub(6 + 4, 4, Mat3x3(f2)*DTmp);


	/* D21 */
	DTmp = FDE.GetMat21()*dCoef;

	WM.Add(4, 1, DTmp);
	WM.Sub(4, 6 + 1, DTmp);
	WM.Sub(6 + 4, 1, DTmp);
	WM.Add(6 + 4, 6 + 1, DTmp);

	/* D21 * [d1 x] */
	FTmp = DTmp*Mat3x3(f1);

	WM.Sub(4, 4, FTmp);
	WM.Add(6 + 4, 4, FTmp);

	/* D21 * [f2 x] */
	FTmp = DTmp*Mat3x3(f2);

	WM.Add(4, 6 + 4, FTmp);
	WM.Sub(6 + 4, 6 + 4, FTmp);


	/* D22 */
	DTmp = FDE.GetMat22()*dCoef;

	WM.Sub(4, 6 + 4, DTmp);
	WM.Add(6 + 4, 6 + 4, DTmp);

	/* D22 + [M x] */
	DTmp += Mat3x3(F.GetVec2()*dCoef);

	WM.Add(4, 4, DTmp);
	WM.Sub(6 + 4, 4, DTmp);


	/* [F x] [f1 x] */
	WM.Sub(4, 4, Mat3x3((F.GetVec1())*dCoef, f1));

	/* [F x] [f2 x] */
	WM.Add(6 + 4, 6 + 4, Mat3x3((F.GetVec1())*dCoef, f2));
}

/* assemblaggio residuo */
SubVectorHandler&
ElasticJoint::AssRes(SubVectorHandler& WorkVec,
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
		WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstMomIndex+iCnt);
	}

	AssVec(WorkVec);

	return WorkVec;
}


void
ElasticJoint::AssVec(SubVectorHandler& WorkVec)
{
	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
	Vec3 f1(pNode1->GetRCurr()*tilde_f1);
	Vec3 f2(pNode2->GetRCurr()*tilde_f2);

	if (bFirstRes) {
		bFirstRes = false;

	} else {
		Mat3x3 R1hT(R1h.Transpose());
		Mat3x3 R2h(pNode2->GetRCurr()*tilde_R2h);
		Vec3 tilde_d(R1hT*(pNode2->GetXCurr() + f2 - pNode1->GetXCurr() - f1));

		tilde_k = Vec6(tilde_d, RotManip::VecRot(R1hT*R2h));

		ConstitutiveLaw6DOwner::Update(tilde_k);
	}

	F = MultRV(GetF(), R1h);

	WorkVec.Add(1, F.GetVec1());
	WorkVec.Add(4, f1.Cross(F.GetVec1()) + F.GetVec2());
	WorkVec.Sub(6 + 1, F.GetVec1());
	WorkVec.Sub(6 + 4, f2.Cross(F.GetVec1()) + F.GetVec2());
}

void
ElasticJoint::AfterPredict(VectorHandler& /* X */ ,
		VectorHandler& /* XP */ )
{
	/* Calcola le deformazioni, aggiorna il legame costitutivo
	 * e crea la FDE */

	/* Recupera i dati */
	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);
	Mat3x3 R2h(pNode2->GetRRef()*tilde_R2h);
	Mat3x3 R1hT(R1h.Transpose());

	/* Calcola la deformazione corrente nel sistema locale (nodo a) */
	ThetaRef = RotManip::VecRot(R1hT*R2h);

	/* Calcola l'inversa di Gamma di ThetaRef */
	Mat3x3 GammaRefm1 = RotManip::DRot_I(ThetaRef);

	Vec3 f1(pNode1->GetRRef()*tilde_f1);
	Vec3 f2(pNode2->GetRRef()*tilde_f2);
	Vec3 tilde_d(R1hT*(pNode2->GetXCurr() + f2 - pNode1->GetXCurr() - f1));

	tilde_k = Vec6(tilde_d, ThetaRef);

	/* Aggiorna il legame costitutivo */
	ConstitutiveLaw6DOwner::Update(tilde_k);

	/* Chiede la matrice tangente di riferimento e la porta
	 * nel sistema globale */
	/* FIXME: horrible */
        FDE = MultRMRt(ConstitutiveLaw6DOwner::GetFDE()*Mat6x6(Eye3, Zero3x3, Zero3x3, GammaRefm1), R1h);

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
		WM.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
		WM.PutRowIndex(6 + iCnt, iNode2FirstPosIndex+iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex+iCnt);
	}

	AssMat(WM, 1.);

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
ElasticJoint::InitialAssRes(SubVectorHandler& WorkVec,
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
		WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstPosIndex+iCnt);
	}

	AssVec(WorkVec);

	return WorkVec;
}

