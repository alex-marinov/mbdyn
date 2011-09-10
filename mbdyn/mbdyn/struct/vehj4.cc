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
#include "vehj4.h"

#include "matvecexp.h"
#include "Rot.hh"

/* DeformableAxialJoint - begin */

/* assemblaggio jacobiano */
/* Used by all "attached" hinges */
void
DeformableAxialJoint::AssMatM(FullSubMatrixHandler& WMA,
	doublereal dCoef)
{
	/* M was updated by AssRes */
	Vec3 M(pNode1->GetRCurr()*(tilde_R1h.GetVec(3)*(dM*dCoef)));
	Mat3x3 MTmp(MatCross, M);

	WMA.Add(1, 1, MTmp);
	WMA.Sub(4, 1, MTmp);
}

/* Used by all hinges */
void
DeformableAxialJoint::AssMatMDE(FullSubMatrixHandler& WMA,
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
DeformableAxialJoint::AssMatMDEPrime(FullSubMatrixHandler& WMA,
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

/* Costruttore non banale */
DeformableAxialJoint::DeformableAxialJoint(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw1D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Mat3x3& tilde_R1h,
		const Mat3x3& tilde_R2h,
		flag fOut)
: Elem(uL, fOut),
Joint(uL, pDO, fOut),
ConstitutiveLaw1DOwner(pCL),
pNode1(pN1),
pNode2(pN2),
tilde_R1h(tilde_R1h),
tilde_R2h(tilde_R2h),
bFirstRes(false)
{
	ASSERT(pNode1 != NULL);
	ASSERT(pNode2 != NULL);
	ASSERT(pNode1->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pNode2->GetNodeType() == Node::STRUCTURAL);
}

/* Distruttore */
DeformableAxialJoint::~DeformableAxialJoint(void)
{
	NO_OP;
}

/* Contributo al file di restart */
std::ostream&
DeformableAxialJoint::Restart(std::ostream& out) const
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
DeformableAxialJoint::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
		Mat3x3 R2h(pNode2->GetRCurr()*tilde_R2h);
		Mat3x3 R(R1h.MulTM(R2h));

		Vec3 v(0., 0., GetF());
		Joint::Output(OH.Joints(), "DeformableHinge", GetLabel(),
			Zero3, v, Zero3, R1h*v) << " " << RotManip::VecRot(R)(3);

		if (GetConstLawType() & ConstLawType::VISCOUS) {
			OH.Joints() << " " << R1h.GetVec(3).Dot(pNode2->GetWCurr() - pNode1->GetWCurr());
		}

		OH.Joints() << std::endl;
	}
}

void
DeformableAxialJoint::AfterPredict(VectorHandler& /* X */ ,
		VectorHandler& /* XP */ )
{
	bFirstRes = true;

	AfterPredict();
}

void
DeformableAxialJoint::SetValue(DataManager *pDM,
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
			ConstitutiveLaw1DOwner::SetValue(pDM, X, XP, ph);
		}
	}
}

void
DeformableAxialJoint::SetInitialValue(VectorHandler& /* X */ )
{
	AfterPredict();
}

Hint *
DeformableAxialJoint::ParseHint(DataManager *pDM, const char *s) const
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

	return ConstitutiveLaw1DOwner::ParseHint(pDM, s);
}

unsigned int
DeformableAxialJoint::iGetNumPrivData(void) const
{
	return 3 + ConstitutiveLaw1DOwner::iGetNumPrivData();
}

unsigned int
DeformableAxialJoint::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);

	if (s[1] == '\0') {
		switch (s[0]) {
		case 'r':
			return 1;

		case 'w':
			return 2;

		case 'M':
			return 3;
		}
	}

	const size_t l = STRLENOF("constitutiveLaw.");
	if (strncmp(s, "constitutiveLaw.", l) == 0) {
		unsigned idx = ConstitutiveLaw1DOwner::iGetPrivDataIdx(&s[l]);
		if (idx > 0) {
			return 3 + idx;
		}
	}

	return 0;
}

doublereal
DeformableAxialJoint::dGetPrivData(unsigned int i) const
{
	ASSERT(i > 0);

	ASSERT(i <= iGetNumPrivData());

	switch (i) {
	case 1:
	{
		/* NOTE: this is correct also in the invariant case */
		Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
		Mat3x3 R2h(pNode2->GetRCurr()*tilde_R2h);

		Vec3 v(RotManip::VecRot(R1h.MulTM(R2h)));

		return v(3);
	}

	case 2:
	{
		Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
		Vec3 w(R1h.MulTV(pNode2->GetWCurr() - pNode1->GetWCurr()));

		return w(3);
	}

	case 3:
		return GetF();

	default:
		return ConstitutiveLaw1DOwner::dGetPrivData(i - 3);
	}
}

/* DeformableAxialJoint - end */


/* ElasticAxialJoint - begin */

ElasticAxialJoint::ElasticAxialJoint(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw1D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Mat3x3& tilde_R1h,
		const Mat3x3& tilde_R2h,
		flag fOut)
: Elem(uL, fOut),
DeformableAxialJoint(uL, pDO, pCL, pN1, pN2, tilde_R1h, tilde_R2h, fOut),
dThetaRef(0.)
{
	NO_OP;
}

ElasticAxialJoint::~ElasticAxialJoint(void)
{
	NO_OP;
}

void
ElasticAxialJoint::AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP)
{
	ConstitutiveLaw1DOwner::AfterConvergence(dThetaCurr);
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
ElasticAxialJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering ElasticAxialJoint::AssJac()" << std::endl);

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
ElasticAxialJoint::AssMats(VariableSubMatrixHandler& WorkMatA,
		VariableSubMatrixHandler& WorkMatB,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering ElasticAxialJoint::AssJac()" << std::endl);

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
ElasticAxialJoint::AfterPredict(void)
{
	/* Calcola le deformazioni, aggiorna il legame costitutivo
	 * e crea la MDE */

	/* Recupera i dati */
	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);
	Mat3x3 R2h(pNode2->GetRRef()*tilde_R2h);

	/* Calcola la deformazione corrente nel sistema locale (nodo a) */
	dThetaCurr = dThetaRef = RotManip::VecRot(R1h.MulTM(R2h))(3);

	/* Aggiorna il legame costitutivo */
	ConstitutiveLaw1DOwner::Update(dThetaRef);

	/* Chiede la matrice tangente di riferimento e la porta
	 * nel sistema globale */
	Vec3 e1z(R1h.GetVec(3));
	MDE = e1z.Tens(e1z*ConstitutiveLaw1DOwner::GetFDE());
}

void
ElasticAxialJoint::AssMat(FullSubMatrixHandler& WMA, doublereal dCoef)
{
#if 0
	// Calcola l'inversa di Gamma di ThetaRef
	Mat3x3 GammaCurrm1 = RotManip::DRot_I(ThetaCurr);

	// Chiede la matrice tangente di riferimento e la porta
	// nel sistema globale
	MDE = R1h*ConstitutiveLaw1DOwner::GetFDE()*GammaCurrm1.MulMT(R1h);
#endif

	AssMatM(WMA, dCoef);
	AssMatMDE(WMA, dCoef);
}

/* assemblaggio residuo */
SubVectorHandler&
ElasticAxialJoint::AssRes(SubVectorHandler& WorkVec,
		doublereal /* dCoef */ ,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering ElasticAxialJoint::AssRes()" << std::endl);

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
ElasticAxialJoint::AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& /* XCurr */,
		const VectorHandler& /* XPrimeCurr */, 
		const VectorHandler& /* XPrimePrimeCurr */, 
		int iOrder)
{
	DEBUGCOUT("Entering ElasticAxialJoint::AssRes()" << std::endl);

	ASSERT(iOrder = -1);
	
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
ElasticAxialJoint::AssVec(SubVectorHandler& WorkVec)
{
	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);

	if (bFirstRes) {
		bFirstRes = false;

	} else {
		Mat3x3 R2h(pNode2->GetRCurr()*tilde_R2h);
		dThetaCurr = RotManip::VecRot(R1h.MulTM(R2h))(3);
		ConstitutiveLaw1DOwner::Update(dThetaCurr);
	}

	/* Couple attached to node 1 */
	M = R1h.GetVec(3)*GetF();

	WorkVec.Add(1, M);
	WorkVec.Sub(4, M);
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
ElasticAxialJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& /* XCurr */ )
{
	DEBUGCOUT("Entering ElasticAxialJoint::InitialAssJac()" << std::endl);

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
ElasticAxialJoint::InitialAssRes(SubVectorHandler& WorkVec,
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

/* ElasticAxialJoint - end */


/* ViscousAxialJoint - begin */

ViscousAxialJoint::ViscousAxialJoint(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw1D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Mat3x3& tilde_R1h,
		const Mat3x3& tilde_R2h,
		flag fOut)
: Elem(uL, fOut),
DeformableAxialJoint(uL, pDO, pCL, pN1, pN2, tilde_R1h, tilde_R2h, fOut)
{
	NO_OP;
}

ViscousAxialJoint::~ViscousAxialJoint(void)
{
	NO_OP;
}

void
ViscousAxialJoint::AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP)
{
	ConstitutiveLaw1DOwner::AfterConvergence(0, dOmega);
}

void
ViscousAxialJoint::AfterPredict(void)
{
	/* Calcola le deformazioni, aggiorna il legame costitutivo
	 * e crea la MDE */

	/* Aggiorna il legame costitutivo */
	Vec3 e1z(pNode1->GetRRef()*tilde_R1h.GetVec(3));
	dOmega = e1z.Dot(pNode2->GetWRef() - pNode1->GetWRef());
	ConstitutiveLaw1DOwner::Update(0., dOmega);

	/* Chiede la matrice tangente di riferimento e la porta
	 * nel sistema globale */
	MDEPrime = e1z.Tens(e1z*GetFDEPrime());
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
ViscousAxialJoint::AssJac(VariableSubMatrixHandler& WorkMat,
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
ViscousAxialJoint::AssMats(VariableSubMatrixHandler& WorkMatA,
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
ViscousAxialJoint::AssMats(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB,
		doublereal dCoef)
{
	AssMatM(WMA, dCoef);
	AssMatMDEPrime(WMA, WMB, dCoef);
}

/* assemblaggio residuo */
SubVectorHandler&
ViscousAxialJoint::AssRes(SubVectorHandler& WorkVec,
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
ViscousAxialJoint::AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& /* XCurr */,
		const VectorHandler& /* XPrimeCurr */, 
		const VectorHandler& /* XPrimePrimeCurr */, 
		int iOrder)
{
	ASSERT(iOrder = -1);

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
ViscousAxialJoint::AssVec(SubVectorHandler& WorkVec)
{
	Vec3 e1z(pNode1->GetRCurr()*tilde_R1h.GetVec(3));

	if (bFirstRes) {
		bFirstRes = false;

	} else {
		dOmega = e1z.Dot(pNode2->GetWCurr() - pNode1->GetWCurr());

		ConstitutiveLaw1DOwner::Update(0., dOmega);
	}

	M = e1z*ConstitutiveLaw1DOwner::GetF();

	WorkVec.Add(1, M);
	WorkVec.Sub(4, M);
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
ViscousAxialJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
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

	Vec3 e1z(pNode1->GetRCurr()*tilde_R1h.GetVec(3));

	const Vec3& W2(pNode2->GetWCurr());

	MDEPrime = e1z.Tens(e1z*ConstitutiveLaw1DOwner::GetFDEPrime());

	Mat3x3 Tmp(MDEPrime*Mat3x3(MatCross, W2));
	WM.Add(4, 7, Tmp);
	WM.Sub(1, 7, Tmp);

	Tmp += Mat3x3(MatCross, e1z*ConstitutiveLaw1DOwner::GetF());
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
ViscousAxialJoint::InitialAssRes(SubVectorHandler& WorkVec,
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

	Vec3 e1z(pNode1->GetRCurr()*tilde_R1h.GetVec(3));

	if (bFirstRes) {
		bFirstRes = false;

	} else {
		dOmega = e1z.Dot(pNode2->GetWCurr() - pNode1->GetWCurr());

		ConstitutiveLaw1DOwner::Update(0., dOmega);
	}

	Vec3 M(e1z*ConstitutiveLaw1DOwner::GetF());

	WorkVec.Add(1, M);
	WorkVec.Sub(4, M);

	return WorkVec;
}

/* ViscousAxialJoint - end */


/* ViscoElasticAxialJoint - begin */

ViscoElasticAxialJoint::ViscoElasticAxialJoint(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw1D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Mat3x3& tilde_R1h,
		const Mat3x3& tilde_R2h,
		flag fOut)
: Elem(uL, fOut),
DeformableAxialJoint(uL, pDO, pCL, pN1, pN2, tilde_R1h, tilde_R2h, fOut),
dThetaRef(0.)
{
	NO_OP;
}

ViscoElasticAxialJoint::~ViscoElasticAxialJoint(void)
{
	NO_OP;
}

void
ViscoElasticAxialJoint::AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP)
{
	ConstitutiveLaw1DOwner::AfterConvergence(dThetaCurr, dOmega);
}

void
ViscoElasticAxialJoint::AfterPredict(void)
{
	/* Computes strains, updates constitutive law and generates
	 * MDE and MDEPrime */

	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);
	Mat3x3 R2h(pNode2->GetRRef()*tilde_R2h);

	/* Current strain in material reference frame (node 1) */
	dThetaCurr = dThetaRef = RotManip::VecRot(R1h.MulTM(R2h))(3);

	/* Relative angular velocity */
	Vec3 e1z(R1h.GetVec(3));
	dOmega = e1z.Dot(pNode2->GetWRef() - pNode1->GetWRef());

	/* Updates constitutive law */
	ConstitutiveLaw1DOwner::Update(dThetaRef, dOmega);

	/* don't repeat the above operations during AssRes */
	bFirstRes = true;

	/* Tangent matrices are updated and projected in the global
	 * reference frame; they won't change during the solution
	 * of the current time step, according to the updated-updated
	 * approach */
	MDE = e1z.Tens(e1z*ConstitutiveLaw1DOwner::GetFDE());
	MDEPrime = e1z.Tens(e1z*GetFDEPrime());
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
ViscoElasticAxialJoint::AssJac(VariableSubMatrixHandler& WorkMat,
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
ViscoElasticAxialJoint::AssMats(VariableSubMatrixHandler& WorkMatA,
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
ViscoElasticAxialJoint::AssMats(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB,
		doublereal dCoef)
{
	AssMatM(WMA, dCoef);
	AssMatMDE(WMA, dCoef);
	AssMatMDEPrime(WMA, WMB, dCoef);
}

/* assemblaggio residuo */
SubVectorHandler&
ViscoElasticAxialJoint::AssRes(SubVectorHandler& WorkVec,
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
ViscoElasticAxialJoint::AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& /* XCurr */,
		const VectorHandler& /* XPrimeCurr */, 
		const VectorHandler& /* XPrimePrimeCurr */, 
		int iOrder)
{
	ASSERT(iOrder = -1);

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
ViscoElasticAxialJoint::AssVec(SubVectorHandler& WorkVec)
{
	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
	Vec3 e1z(R1h.GetVec(3));

	if (bFirstRes) {
		bFirstRes = false;

	} else {
		Mat3x3 R2h(pNode2->GetRCurr()*tilde_R2h);

		/* orientazione intermedia */
		dThetaCurr = RotManip::VecRot(R1h.MulTM(R2h))(3);

		/* velocita' relativa nel riferimento intermedio */
		dOmega = e1z.Dot(pNode2->GetWCurr() - pNode1->GetWCurr());

		/* aggiorna il legame costitutivo */
		ConstitutiveLaw1DOwner::Update(dThetaCurr, dOmega);
	}

	M = e1z*ConstitutiveLaw1DOwner::GetF();

	WorkVec.Add(1, M);
	WorkVec.Sub(4, M);
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
ViscoElasticAxialJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
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

	Vec3 e1z(pNode1->GetRCurr()*tilde_R1h.GetVec(3));
	const Vec3& W2(pNode2->GetWCurr());

	MDE = e1z.Tens(e1z*ConstitutiveLaw1DOwner::GetFDE());
	MDEPrime = e1z.Tens(e1z*ConstitutiveLaw1DOwner::GetFDEPrime());

	Mat3x3 Tmp(MDE + MDEPrime*Mat3x3(MatCross, W2));
	WM.Add(4, 7, Tmp);
	WM.Sub(1, 7, Tmp);

	Tmp += Mat3x3(MatCross, e1z*ConstitutiveLaw1DOwner::GetF());
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
ViscoElasticAxialJoint::InitialAssRes(SubVectorHandler& WorkVec,
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
	Vec3 e1z(R1h.GetVec(3));

	if (bFirstRes) {
		bFirstRes = false;

	} else {
		Mat3x3 R2h(pNode2->GetRCurr()*tilde_R2h);

		dThetaCurr = RotManip::VecRot(R1h.MulTM(R2h))(3);
		dOmega = e1z.Dot(pNode2->GetWCurr() - pNode1->GetWCurr());

		ConstitutiveLaw1DOwner::Update(dThetaCurr, dOmega);
	}

	M = e1z*ConstitutiveLaw1DOwner::GetF();

	WorkVec.Add(1, M);
	WorkVec.Sub(4, M);

	return WorkVec;
}

/* ViscoElasticAxialJoint - end */

/* InvAngularCLR - end */
