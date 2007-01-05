/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2007
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

/* Cerniera pilotata */

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "dataman.h"
#include "gimbal.h"
#include "Rot.hh"

/* GimbalRotationJoint - begin */

/* Costruttore non banale */
GimbalRotationJoint::GimbalRotationJoint(unsigned int uL,
				 const DofOwner* pDO,
				 const StructNode* pN1,
				 const StructNode* pN2,
				 const Mat3x3& R1,
				 const Mat3x3& R2,
				 flag fOut)
: Elem(uL, fOut),
Joint(uL, pDO, fOut),
pNode1(pN1), pNode2(pN2), R1h(R1), R2h(R2),
M(0.), dTheta(0.), dPhi(0.)
{
	ASSERT(pNode1 != NULL);
	ASSERT(pNode2 != NULL);
	ASSERT(pNode1->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pNode2->GetNodeType() == Node::STRUCTURAL);
}


/* Distruttore */
GimbalRotationJoint::~GimbalRotationJoint(void)
{
	NO_OP;
}


/* Contributo al file di restart */
std::ostream&
GimbalRotationJoint::Restart(std::ostream& out) const
{
	Joint::Restart(out) << ", gimbal rotation, "
		<< pNode1->GetLabel() << ", reference, node, 1, ",
		(R1h.GetVec(1)).Write(out, ", ")
		<< ", 2, ", (R1h.GetVec(2)).Write(out, ", ") << ", "
		<< pNode2->GetLabel() << ", reference, node, 1, ",
		(R2h.GetVec(1)).Write(out, ", ")
		<< ", 2, ", (R2h.GetVec(2)).Write(out, ", ") << ';' << std::endl;
	return out;
}


void
GimbalRotationJoint::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		Mat3x3 Ra(pNode1->GetRCurr());

		Vec3 d(MatR2EulerAngles(pNode1->GetRCurr().Transpose()*pNode2->GetRCurr()));
		Joint::Output(OH.Joints(), "Gimbal", GetLabel(),
				Zero3, M, Zero3, Ra*M)
			<< " " << dTheta << " " << dPhi << std::endl;
	}
}


/* assemblaggio jacobiano */
VariableSubMatrixHandler&
GimbalRotationJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering GimbalRotationJoint::AssJac()" << std::endl);

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
	integer iFirstReactionIndex = iGetFirstIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(3 + iCnt, iNode2FirstMomIndex + iCnt);
		WM.PutColIndex(3 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	for (int iCnt = 1; iCnt <= 5; iCnt++) {
		WM.PutRowIndex(6 + iCnt, iFirstReactionIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iFirstReactionIndex + iCnt);
	}

	AssMat(WM, dCoef);

	return WorkMat;
}


void
GimbalRotationJoint::AfterPredict(VectorHandler& /* X */ ,
		VectorHandler& /* XP */ )
{
	NO_OP;
}


void
GimbalRotationJoint::AssMat(FullSubMatrixHandler& WM, doublereal dCoef)
{
	Mat3x3 Ra(pNode1->GetRRef()*R1h);
	Mat3x3 RaT(Ra.Transpose());

	doublereal dCosTheta = cos(dTheta);
	doublereal dSinTheta = sin(dTheta);
	doublereal dCosPhi = cos(dPhi);
	doublereal dSinPhi = sin(dPhi);

	/* coppie */
	/* termini in Delta lambda */
	WM.Add(1, 6 + 1, Ra);
	WM.Sub(3 + 1, 6 + 1, Ra);

	/* termini in Delta g_a */
	Mat3x3 MTmp(Ra*(M*dCoef));
	WM.Sub(1, 1, MTmp);
	WM.Add(3 + 1, 1, MTmp);

	/* equazioni di vincolo */
	/* termini in Delta g_a, Delta g_b */
	MTmp = RaT*dCoef;
	WM.Add(6 + 1, 1, MTmp);
	WM.Sub(6 + 1, 3 + 1, MTmp);

	/* termini in Delta theta */
	WM.IncCoef(6 + 1, 9 + 1, dSinTheta*dSinPhi);
	WM.IncCoef(6 + 2, 9 + 1, 1. + dCosPhi);
	WM.IncCoef(6 + 3, 9 + 1, dCosTheta*dSinPhi);

	/* termini in Delta phi */
	WM.IncCoef(6 + 1, 9 + 2, dCosTheta);
	WM.DecCoef(6 + 3, 9 + 2, dSinTheta);

	/* equazione in theta */
	/* termini in Delta lambda */
	WM.IncCoef(9 + 1, 6 + 1, dSinTheta*dSinPhi);
	WM.IncCoef(9 + 1, 6 + 2, 1. + dCosPhi);
	WM.IncCoef(9 + 1, 6 + 3, dCosTheta*dSinPhi);

	/* termine in Delta theta */
	WM.IncCoef(9 + 1, 9 + 1, dSinPhi*(dCosTheta*M(1) - dSinTheta*M(3)));

	/* termine in Delta phi */
	WM.IncCoef(9 + 1, 9 + 2, dCosPhi*(dSinTheta*M(1) + dCosTheta*M(3)) - dSinPhi*M(2));

	/* equazione in phi */
	/* termini in Delta lambda */
	WM.IncCoef(9 + 2, 6 + 1, dCosTheta);
	WM.DecCoef(9 + 2, 6 + 3, dSinTheta);

	/* termine in Delta theta */
	WM.DecCoef(9 + 2, 9 + 1, dSinTheta*M(1) + dCosTheta*M(3));
}


/* assemblaggio residuo */
SubVectorHandler&
GimbalRotationJoint::AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering GimbalRotationJoint::AssRes()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera gli indici */
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex() + 3;
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex() + 3;
	integer iFirstReactionIndex = iGetFirstIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WorkVec.PutRowIndex(3 + iCnt, iNode2FirstMomIndex + iCnt);
	}

	for (int iCnt = 1; iCnt <= 5; iCnt++) {
		WorkVec.PutRowIndex(6 + iCnt, iFirstReactionIndex + iCnt);
	}

	M = Vec3(XCurr, iFirstReactionIndex + 1);
	dTheta = XCurr(iFirstReactionIndex + 3 + 1);
	dPhi = XCurr(iFirstReactionIndex + 3 + 2);

	AssVec(WorkVec, dCoef);

	return WorkVec;
}


void
GimbalRotationJoint::AssVec(SubVectorHandler& WorkVec, doublereal dCoef)
{
	Mat3x3 Ra(pNode1->GetRCurr()*R1h);
	Mat3x3 Rb(pNode2->GetRCurr()*R2h);

	Mat3x3 ExpTheta(RotManip::Rot(Vec3(0., dTheta, 0.)));
	Mat3x3 ExpPhi(RotManip::Rot(Vec3(dPhi, 0., 0.)));

	doublereal dCosTheta = cos(dTheta);
	doublereal dSinTheta = sin(dTheta);
	doublereal dCosPhi = cos(dPhi);
	doublereal dSinPhi = sin(dPhi);

	Vec3 MTmp(Ra*M);

	WorkVec.Sub(1, MTmp);
	WorkVec.Add(3 + 1, MTmp);

	WorkVec.Add(6 + 1, RotManip::VecRot(Ra.Transpose()*Rb) - RotManip::VecRot(ExpTheta*ExpPhi*ExpTheta));

	WorkVec.DecCoef(9 + 1, dSinTheta*dSinPhi*M(1) + (1. + dCosPhi)*M(2) + dCosTheta*dSinPhi*M(3));
	WorkVec.DecCoef(9 + 2, dCosTheta*M(1) - dSinTheta*M(3));
}


/*
 * FIXME: the initial assembly is flawed:
 * there is no constraint derivative...
 */

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
GimbalRotationJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering GimbalRotationJoint::InitialAssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex() + 3;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex() + 3;
	integer iFirstReactionIndex = iGetFirstIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(3 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutColIndex(3 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	for (int iCnt = 1; iCnt <= 5; iCnt++) {
		WM.PutRowIndex(6 + iCnt, iFirstReactionIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iFirstReactionIndex + iCnt);
	}

	AssMat(WM, 1.);

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
GimbalRotationJoint::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering GimbalRotationJoint::InitialAssRes()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex() + 3;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex() + 3;
	integer iFirstReactionIndex = iGetFirstIndex();

	/* Setta gli indici del vettore */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WorkVec.PutRowIndex(3 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	for (int iCnt = 1; iCnt <= 5; iCnt++) {
		WorkVec.PutRowIndex(6 + iCnt, iFirstReactionIndex + iCnt);
	}

	AssVec(WorkVec, 1);

	return WorkVec;
}

/* Dati privati (aggiungere magari le reazioni vincolari) */
unsigned int
GimbalRotationJoint::iGetNumPrivData(void) const
{
	return 5;
}

unsigned int
GimbalRotationJoint::iGetPrivDataIdx(const char *s) const
{
	if (strncmp(s, "lambda[", sizeof("lambda[") - 1) == 0) {
		s += sizeof("lambda[") - 1;
		if (s[1] != ']') {
			return 0;
		}

		switch (s[0]) {
		case '1':
		case '2':
		case '3':
			return s[0] - '0';

		default:
			return 0;
		}

	} else if (strcmp(s, "theta") == 0) {
		return 4;

	} else if (strcmp(s, "phi") == 0) {
		return 5;
	}

	return 0;
}

doublereal
GimbalRotationJoint::dGetPrivData(unsigned int i) const
{
	switch (i) {
	case 1:
	case 2:
	case 3:
		return M(i);

	case 4:
		return dTheta;

	case 5:
		return dPhi;
	}

	throw ErrGeneric();
}

/* GimbalRotationJoint - end */
