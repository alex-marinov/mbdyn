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

/* Forze */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cfloat>

#include "dataman.h"
#include "strforce.h"
#include "strforce_impl.h"
#include "tpldrive_impl.h"

/* StructuralForce - begin */

/* Costruttore */
StructuralForce::StructuralForce(unsigned int uL,
	const StructNode* pN,
	const TplDriveCaller<Vec3>* pDC,
	flag fOut)
: Elem(uL, fOut),
Force(uL, fOut),
f(pDC),
pNode(pN)
{
	ASSERT(pNode != NULL);
	ASSERT(pNode->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pDC != NULL);
}


StructuralForce::~StructuralForce(void)
{
	NO_OP;
}

void
StructuralForce::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
	connectedNodes.resize(1);
	connectedNodes[0] = pNode;
}

/* StructuralForce - end */


/* AbsoluteForce - begin */

/* Costruttore non banale */

AbsoluteForce::AbsoluteForce(unsigned int uL,
	const StructNode* pN,
	const TplDriveCaller<Vec3>* pDC,
	const Vec3& TmpArm,
	flag fOut)
: Elem(uL, fOut),
StructuralForce(uL, pN, pDC, fOut),
Arm(TmpArm)
{
	NO_OP;
}


AbsoluteForce::~AbsoluteForce(void)
{
	NO_OP;
}

void
AbsoluteForce::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 6;
	*piNumCols = 3;
}

void
AbsoluteForce::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 12;
	*piNumCols = 6;
}

/* Contributo al file di restart */
std::ostream&
AbsoluteForce::Restart(std::ostream& out) const
{
	Force::Restart(out) << ", absolute, "
		<< pNode->GetLabel() << ", position, reference, node, ",
		Arm.Write(out, ", ") << ", ";
	return f.pGetDriveCaller()->Restart(out) << ';' << std::endl;
}

VariableSubMatrixHandler&
AbsoluteForce::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering AbsoluteForce::AssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	WM.ResizeReset(3, 3);

	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex() + 3;
	integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex() + 3;
	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iFirstMomentumIndex + iCnt);
		WM.PutColIndex(iCnt, iFirstPositionIndex + iCnt);
	}

	/* Dati */
	Vec3 TmpArm(pNode->GetRRef()*Arm);

	/* |    F/\   |           |   F  |
	 * |          | Delta_g = |      |
	 * | (d/\F)/\ |           | d/\F |
	 */

	WM.Sub(1, 1, Mat3x3(MatCrossCross, f.Get(), TmpArm*dCoef));

	return WorkMat;
}


/* Assembla il residuo */
SubVectorHandler&
AbsoluteForce::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering AbsoluteForce::AssRes()" << std::endl);

	integer iNumRows;
	integer iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici delle incognite del nodo */
	integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex + iCnt);
	}

	const Mat3x3& R(pNode->GetRCurr());
	Vec3 F(f.Get());
	Vec3 M((R*Arm).Cross(F));

	WorkVec.Add(1, F);
	WorkVec.Add(4, M);

	return WorkVec;
}

/* Inverse Dynamics*/
SubVectorHandler&
AbsoluteForce::AssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ ,
	const VectorHandler& /* XPrimePrimeCurr */ ,
	int iOrder)
{
	DEBUGCOUT("Entering AbsoluteForce::AssRes()" << std::endl);

	integer iNumRows;
	integer iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici delle incognite del nodo */
	integer iFirstMomentumIndex = pNode->iGetFirstPositionIndex();
	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex + iCnt);
	}

	const Mat3x3& R(pNode->GetRCurr());
	Vec3 F(f.Get());
	Vec3 M((R*Arm).Cross(F));

	WorkVec.Add(1, F);
	WorkVec.Add(4, M);

	return WorkVec;
}

void
AbsoluteForce::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		OH.Forces()
			<< GetLabel()
			<< " " << pNode->GetLabel()
			<< " " << f.Get()
			<< " " << pNode->GetXCurr() + pNode->GetRCurr()*Arm
			<< std::endl;
	}
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
AbsoluteForce::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUT("Entering AbsoluteForce::InitialAssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	WM.ResizeReset(6, 6);

	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex() + 3;
	integer iFirstVelocityIndex = iFirstPositionIndex + 6;
	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
		WM.PutRowIndex(3+iCnt, iFirstVelocityIndex + iCnt);
		WM.PutColIndex(iCnt, iFirstPositionIndex + iCnt);
		WM.PutColIndex(3+iCnt, iFirstVelocityIndex + iCnt);
	}

	/* Dati */
	Vec3 TmpArm(pNode->GetRRef()*Arm);
	Vec3 TmpDir = f.Get();
	const Vec3& Omega(pNode->GetWRef());

	/* |    F/\   |           |   F  |
	 * |          | Delta_g = |      |
	 * | (d/\F)/\ |           | d/\F |
	 */

	Mat3x3 MTmp(MatCrossCross, TmpDir, TmpArm);

	WM.Sub(1, 1, MTmp);
	WM.Sub(4, 1, Mat3x3(MatCrossCross, TmpDir, Omega)*Mat3x3(MatCross, TmpArm));
	WM.Sub(4, 4, MTmp);

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
AbsoluteForce::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUT("Entering AbsoluteForce::InitialAssRes()" << std::endl);

	integer iNumRows;
	integer iNumCols;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici delle incognite del nodo */
	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	integer iFirstVelocityIndex = iFirstPositionIndex + 6;
	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
		WorkVec.PutRowIndex(6+iCnt, iFirstVelocityIndex + iCnt);
	}

	/* Dati */
	const Mat3x3& R(pNode->GetRCurr());
	Vec3 TmpDir(f.Get());
	Vec3 TmpArm(R*Arm);
	const Vec3& Omega(pNode->GetWCurr());

	WorkVec.Add(1, TmpDir);
	WorkVec.Add(4, TmpArm.Cross(TmpDir));
	/* In 7 non c'e' nulla */
	WorkVec.Add(10, (Omega.Cross(TmpArm)).Cross(TmpDir));

	return WorkVec;
}

/* AbsoluteForce - end */


/* FollowerForce - begin */

/* Costruttore non banale */

FollowerForce::FollowerForce(unsigned int uL, const StructNode* pN,
	const TplDriveCaller<Vec3>* pDC,
	const Vec3& TmpArm,
	flag fOut)
: Elem(uL, fOut),
StructuralForce(uL, pN, pDC, fOut),
Arm(TmpArm)
{
	NO_OP;
}


FollowerForce::~FollowerForce(void)
{
	NO_OP;
}


void
FollowerForce::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 6;
	*piNumCols = 3;
}


void
FollowerForce::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 12;
	*piNumCols = 6;
}


/* Contributo al file di restart */
std::ostream&
FollowerForce::Restart(std::ostream& out) const
{
	Force::Restart(out) << ", follower, "
		<< pNode->GetLabel()
		<< ", position, reference, node, ",
		Arm.Write(out, ", ") << ", ";
	return f.pGetDriveCaller()->Restart(out) << ';' << std::endl;
}


VariableSubMatrixHandler&
FollowerForce::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering FollowerForce::AssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	integer iFirstRotationIndex = pNode->iGetFirstPositionIndex() + 3;
	integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iFirstMomentumIndex + iCnt);     /* forza */
		WM.PutRowIndex(3 + iCnt, iFirstMomentumIndex + 3 + iCnt); /* coppia */
		WM.PutColIndex(iCnt, iFirstRotationIndex + iCnt);     /* rotazione */
	}

	/* Dati */
	const Mat3x3& R(pNode->GetRRef());
	Vec3 TmpDir(R*(f.Get()*dCoef));
	Vec3 TmpArm(R*Arm);

	/* |    F/\   |           |   F  |
	 * |          | Delta_g = |      |
	 * | (d/\F)/\ |           | d/\F |
	 */

	WM.Add(1, 1, Mat3x3(MatCross, TmpDir));
	WM.Add(4, 1, Mat3x3(MatCross, TmpArm.Cross(TmpDir)));

	return WorkMat;
}


/* Assembla il residuo */
SubVectorHandler&
FollowerForce::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering FollowerForce::AssRes()" << std::endl);

	integer iNumRows;
	integer iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici delle incognite del nodo */
	integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex + iCnt);
	}

	/* Dati */
	const Mat3x3& R(pNode->GetRCurr());
	Vec3 TmpDir = f.Get();
	Vec3 F(R*TmpDir);
	Vec3 M(R*Arm.Cross(TmpDir));

	WorkVec.Add(1, F);
	WorkVec.Add(4, M);

	return WorkVec;
}


/* Inverse Dynamics*/
SubVectorHandler&
FollowerForce::AssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ ,
	const VectorHandler& /* XPrimePrimeCurr */ ,
	int iOrder)
{
	DEBUGCOUT("Entering FollowerForce::AssRes()" << std::endl);

	integer iNumRows;
	integer iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici delle incognite del nodo */
	integer iFirstMomentumIndex = pNode->iGetFirstPositionIndex();
	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex+iCnt);
	}

	/* Dati */
	const Mat3x3& R(pNode->GetRCurr());
	Vec3 TmpDir = f.Get();
	Vec3 F(R*TmpDir);
	Vec3 M(R*Arm.Cross(TmpDir));

	WorkVec.Add(1, F);
	WorkVec.Add(4, M);

	return WorkVec;
}


void
FollowerForce::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		OH.Forces()
			<< GetLabel()
			<< " " << pNode->GetLabel()
			<< " " << pNode->GetRCurr()*f.Get()
			<< " " << pNode->GetXCurr() + pNode->GetRCurr()*Arm
			<< std::endl;
	}
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
FollowerForce::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUT("Entering FollowerForce::InitialAssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	WM.ResizeReset(12, 6);

	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	integer iFirstVelocityIndex = iFirstPositionIndex + 6;
	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
		WM.PutRowIndex(3 + iCnt, iFirstPositionIndex + 3 + iCnt);
		WM.PutRowIndex(6 + iCnt, iFirstVelocityIndex + iCnt);
		WM.PutRowIndex(9 + iCnt, iFirstVelocityIndex + 3 + iCnt);
		WM.PutColIndex(iCnt, iFirstPositionIndex + 3 + iCnt);
		WM.PutColIndex(3 + iCnt, iFirstVelocityIndex + 3 + iCnt);
	}

	/* Dati */
	const Mat3x3& R(pNode->GetRRef());
	Vec3 TmpArm(R*Arm);
	Vec3 TmpDir = R*f.Get();
	const Vec3& Omega(pNode->GetWRef());

	/* |    F/\   |           |   F  |
	 * |          | Delta_g = |      |
	 * | (d/\F)/\ |           | d/\F |
	 */

	WM.Add(1, 1, Mat3x3(MatCross, TmpDir));
	WM.Add(4, 1, Mat3x3(MatCross, TmpArm.Cross(TmpDir)));
	WM.Add(7, 1, Mat3x3(MatCrossCross, Omega, TmpDir));
	WM.Add(7, 4, Mat3x3(MatCross, TmpDir));
	WM.Add(10, 1, Mat3x3(MatCrossCross, Omega, TmpArm.Cross(TmpDir)));
	WM.Add(10, 4, Mat3x3(MatCross, TmpArm.Cross(TmpDir)));

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
FollowerForce::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUT("Entering FollowerForce::InitialAssRes()" << std::endl);

	integer iNumRows;
	integer iNumCols;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici delle incognite del nodo */
	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	integer iFirstVelocityIndex = iFirstPositionIndex + 6;
	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iFirstVelocityIndex + iCnt);
	}

	/* Dati */
	const Mat3x3& R(pNode->GetRCurr());
	Vec3 TmpDir(R*f.Get());
	Vec3 TmpArm(R*Arm);
	const Vec3& Omega(pNode->GetWCurr());

	WorkVec.Add(1, TmpDir);
	WorkVec.Add(4, TmpArm.Cross(TmpDir));
	WorkVec.Add(7, Omega.Cross(TmpDir));
	WorkVec.Add(10, (Omega.Cross(TmpArm)).Cross(TmpDir)
		+ TmpArm.Cross(Omega.Cross(TmpDir)));

	return WorkVec;
}

/* FollowerForce - end */


/* AbsoluteCouple - begin */

/* Costruttore non banale */

AbsoluteCouple::AbsoluteCouple(unsigned int uL, const StructNode* pN,
	const TplDriveCaller<Vec3>* pDC,
	flag fOut)
: Elem(uL, fOut),
StructuralForce(uL, pN, pDC, fOut)
{
	NO_OP;
}


AbsoluteCouple::~AbsoluteCouple(void)
{
	NO_OP;
}


void
AbsoluteCouple::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 3;
	*piNumCols = 1;
}


void
AbsoluteCouple::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 3;
	*piNumCols = 1;
}


/* Contributo al file di restart */
std::ostream&
AbsoluteCouple::Restart(std::ostream& out) const
{
	out << "  couple: " << GetLabel() << ", absolute, "
		<< pNode->GetLabel() << ", ";
	return f.pGetDriveCaller()->Restart(out) << ';' << std::endl;
}


/* Assembla il residuo */
SubVectorHandler&
AbsoluteCouple::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering AbsoluteCouple::AssRes()" << std::endl);

	integer iNumRows;
	integer iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici delle incognite del nodo */
	integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex()+3;
	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex+iCnt);
	}

	WorkVec.Add(1, f.Get());

	return WorkVec;
}

/* Inverse Dynamics*/
SubVectorHandler&
AbsoluteCouple::AssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ ,
	const VectorHandler& /* XPrimePrimeCurr */ ,
	int iOrder)
{
	DEBUGCOUT("Entering AbsoluteCouple::AssRes()" << std::endl);

	integer iNumRows;
	integer iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici delle incognite del nodo */
	integer iFirstMomentumIndex = pNode->iGetFirstPositionIndex()+3;
	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex+iCnt);
	}

	WorkVec.Add(1, f.Get());

	return WorkVec;
}


void
AbsoluteCouple::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		OH.Forces()
			<< GetLabel()
			<< " " << pNode->GetLabel()
			<< " " << f.Get()
			<< std::endl;
	}
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
AbsoluteCouple::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUT("Entering AbsoluteCouple::InitialAssRes()" << std::endl);

	WorkVec.Resize(3);

	/* Indici delle incognite del nodo */
	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex() + 3;
	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
	}

	WorkVec.Add(1, f.Get());

	return WorkVec;
}

/* AbsoluteCouple - end */


/* FollowerCouple - begin */

FollowerCouple::FollowerCouple(unsigned int uL, const StructNode* pN,
	const TplDriveCaller<Vec3>* pDC,
	flag fOut)
: Elem(uL, fOut),
StructuralForce(uL, pN, pDC, fOut)
{
	NO_OP;
}


FollowerCouple::~FollowerCouple(void)
{
	NO_OP;
}


void
FollowerCouple::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 3;
	*piNumCols = 3;
}


void
FollowerCouple::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 6;
	*piNumCols = 6;
}


/* Contributo al file di restart */
std::ostream&
FollowerCouple::Restart(std::ostream& out) const
{
	out << "  couple: " << GetLabel() << ", follower, "
		<< pNode->GetLabel() << ", ";
	return f.pGetDriveCaller()->Restart(out) << ';' << std::endl;
}


VariableSubMatrixHandler&
FollowerCouple::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering FollowerCouple::AssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	integer iFirstRotationIndex = pNode->iGetFirstPositionIndex() + 3;
	integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex() + 3;
	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iFirstMomentumIndex + iCnt);    /* coppia */
		WM.PutColIndex(iCnt, iFirstRotationIndex + iCnt);    /* rotazione */
	}

	/* Dati */
	const Mat3x3& R(pNode->GetRRef());
	Mat3x3 MWedge(MatCross, R*(f.Get()*dCoef));

	/* | M /\| Delta_g = | M | */

	WM.Add(1, 1, MWedge);

	return WorkMat;
}


/* Assembla il residuo */
SubVectorHandler&
FollowerCouple::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering FollowerCouple::AssRes()" << std::endl);

	integer iNumRows;
	integer iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici delle incognite del nodo */
	integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex() + 3;
	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex + iCnt);
	}

	/* Dati */
	const Mat3x3& R(pNode->GetRCurr());
	WorkVec.Add(1, R*f.Get());

	return WorkVec;
}


/* Inverse Dynamics*/
SubVectorHandler&
FollowerCouple::AssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ ,
	const VectorHandler& /* XPrimePrimeCurr */ ,
	int iOrder)
{
	DEBUGCOUT("Entering FollowerCouple::AssRes()" << std::endl);

	integer iNumRows;
	integer iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici delle incognite del nodo */
	integer iFirstMomentumIndex = pNode->iGetFirstPositionIndex()+3;
	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex+iCnt);
	}

	/* Dati */
	Mat3x3 R(pNode->GetRCurr());
	WorkVec.Add(1, R*f.Get());

	return WorkVec;
}

void
FollowerCouple::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		OH.Forces()
			<< GetLabel()
			<< " " << pNode->GetLabel()
			<< " " << pNode->GetRCurr()*f.Get()
			<< std::endl;
	}
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
FollowerCouple::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUT("Entering FollowerCouple::InitialAssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	WM.ResizeReset(6, 6);

	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex() + 3;
	integer iFirstVelocityIndex = iFirstPositionIndex + 6;
	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
		WM.PutRowIndex(3 + iCnt, iFirstVelocityIndex + iCnt);
		WM.PutColIndex(iCnt, iFirstPositionIndex + iCnt);
		WM.PutColIndex(3 + iCnt, iFirstVelocityIndex + iCnt);
	}

	/* Dati */
	Vec3 TmpDir(pNode->GetRRef()*f.Get());
	const Vec3& Omega(pNode->GetWRef());

	/* |    F/\   |           |   F  |
	 * |          | Delta_g = |      |
	 * | (d/\F)/\ |           | d/\F |
	 */

	WM.Add(1, 1, Mat3x3(MatCross, TmpDir));
	WM.Add(4, 1, Mat3x3(MatCrossCross, Omega, TmpDir));
	WM.Add(4, 4, Mat3x3(MatCross, TmpDir));

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
FollowerCouple::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUT("Entering FollowerCouple::InitialAssRes()" << std::endl);

	WorkVec.ResizeReset(6);

	/* Indici delle incognite del nodo */
	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex() + 3;
	integer iFirstVelocityIndex = iFirstPositionIndex + 6;
	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iFirstVelocityIndex + iCnt);
	}

	/* Dati */
	const Mat3x3& R(pNode->GetRCurr());
	Vec3 TmpDir(R*f.Get());
	const Vec3& Omega(pNode->GetWCurr());

	WorkVec.Add(1, TmpDir);
	WorkVec.Add(4, Omega.Cross(TmpDir));

	return WorkVec;
}

/* FollowerCouple - end */


/* StructuralInternalForce - begin */

/* Costruttore */
StructuralInternalForce::StructuralInternalForce(unsigned int uL,
	const StructNode* pN1, const StructNode* pN2,
	const TplDriveCaller<Vec3>* pDC,
	flag fOut)
: Elem(uL, fOut),
Force(uL, fOut),
f(pDC),
pNode1(pN1), pNode2(pN2)
{
	ASSERT(pNode1 != NULL);
	ASSERT(pNode1->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pNode2 != NULL);
	ASSERT(pNode2->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pDC != NULL);
}


StructuralInternalForce::~StructuralInternalForce(void)
{
	NO_OP;
}

void
StructuralInternalForce::GetConnectedNodes(
	std::vector<const Node *>& connectedNodes) const {
	connectedNodes.resize(2);
	connectedNodes[0] = pNode1;
	connectedNodes[1] = pNode2;
}

/* StructuralInternalForce - end */


/* AbsoluteInternalForce - begin */

/* Costruttore non banale */

AbsoluteInternalForce::AbsoluteInternalForce(unsigned int uL,
	const StructNode* pN1, const StructNode* pN2,
	const TplDriveCaller<Vec3>* pDC,
	const Vec3& TmpArm1, const Vec3& TmpArm2,
	flag fOut)
: Elem(uL, fOut),
StructuralInternalForce(uL, pN1, pN2, pDC, fOut),
Arm1(TmpArm1), Arm2(TmpArm2)
{
	NO_OP;
}


AbsoluteInternalForce::~AbsoluteInternalForce(void)
{
	NO_OP;
}


void
AbsoluteInternalForce::WorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 12;
	*piNumCols = 6;
}


void
AbsoluteInternalForce::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 24;
	*piNumCols = 12;
}


/* Contributo al file di restart */
std::ostream&
AbsoluteInternalForce::Restart(std::ostream& out) const
{
	Force::Restart(out) << ", absolute internal, "
		<< pNode1->GetLabel() << ", position, reference, node, ",
		Arm1.Write(out, ", ") << ", "
		<< pNode2->GetLabel() << ", position, reference, node, ",
		Arm2.Write(out, ", ") << ", ";
	return f.pGetDriveCaller()->Restart(out) << ';' << std::endl;
}


VariableSubMatrixHandler&
AbsoluteInternalForce::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering AbsoluteInternalForce::AssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	WM.ResizeReset(6, 6);

	integer iFirstPositionIndex1 = pNode1->iGetFirstPositionIndex() + 3;
	integer iFirstMomentumIndex1 = pNode1->iGetFirstMomentumIndex() + 3;

	integer iFirstPositionIndex2 = pNode2->iGetFirstPositionIndex() + 3;
	integer iFirstMomentumIndex2 = pNode2->iGetFirstMomentumIndex() + 3;

	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iFirstMomentumIndex1 + iCnt);
		WM.PutColIndex(iCnt, iFirstPositionIndex1 + iCnt);

		WM.PutRowIndex(3 + iCnt, iFirstMomentumIndex2 + iCnt);
		WM.PutColIndex(3 + iCnt, iFirstPositionIndex2 + iCnt);
	}

	/* Dati */
	Vec3 TmpArm1(pNode1->GetRRef()*Arm1);
	Vec3 TmpArm2(pNode2->GetRRef()*Arm2);
	Vec3 TmpDir = f.Get()*dCoef;

	/* |    F/\   |           |   F  |
	 * |          | Delta_g = |      |
	 * | (d/\F)/\ |           | d/\F |
	 */

	WM.Sub(1, 1, Mat3x3(MatCrossCross, TmpDir, TmpArm1));
	WM.Add(4, 4, Mat3x3(MatCrossCross, TmpDir, TmpArm2));

	return WorkMat;
}


/* Assembla il residuo */
SubVectorHandler&
AbsoluteInternalForce::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering AbsoluteInternalForce::AssRes()" << std::endl);

	integer iNumRows;
	integer iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici delle incognite del nodo */
	integer iFirstMomentumIndex1 = pNode1->iGetFirstMomentumIndex();
	integer iFirstMomentumIndex2 = pNode2->iGetFirstMomentumIndex();
	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex1 + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iFirstMomentumIndex2 + iCnt);
	}

	/* Dati */
	Vec3 F(f.Get());
	Vec3 M1((pNode1->GetRCurr()*Arm1).Cross(F));
	Vec3 M2(F.Cross(pNode2->GetRCurr()*Arm1));	/* - x2 /\ F */

	WorkVec.Add(1, F);
	WorkVec.Add(4, M1);
	WorkVec.Sub(7, F);
	WorkVec.Sub(10, M2);

	return WorkVec;
}

/* Inverse Dynamics*/
SubVectorHandler&
AbsoluteInternalForce::AssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ ,
	const VectorHandler& /* XPrimePrimeCurr */ ,
	int iOrder)
{
	DEBUGCOUT("Entering AbsoluteInternalForce::AssRes()" << std::endl);

	integer iNumRows;
	integer iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici delle incognite del nodo */
	integer iFirstMomentumIndex1 = pNode1->iGetFirstPositionIndex();
	integer iFirstMomentumIndex2 = pNode2->iGetFirstPositionIndex();
	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex1 + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iFirstMomentumIndex2 + iCnt);
	}

	/* Dati */
	Vec3 F(f.Get());
	Vec3 M1((pNode1->GetRCurr()*Arm1).Cross(F));
	Vec3 M2(F.Cross(pNode2->GetRCurr()*Arm1));	/* - x2 /\ F */

	WorkVec.Add(1, F);
	WorkVec.Add(4, M1);
	WorkVec.Sub(7, F);
	WorkVec.Sub(10, M2);

	return WorkVec;
}


void
AbsoluteInternalForce::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		Vec3 F(f.Get());
		OH.Forces()
			<< GetLabel()
			<< " " << pNode1->GetLabel()
			<< " " << F
			<< " " << pNode1->GetXCurr() + pNode1->GetRCurr()*Arm1
			<< " " << pNode2->GetLabel()
			<< " " << -F
			<< " " << pNode2->GetXCurr() + pNode2->GetRCurr()*Arm2
			<< std::endl;
	}
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
AbsoluteInternalForce::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUT("Entering AbsoluteInternalForce::InitialAssJac()"
		<< std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	WM.ResizeReset(12, 12);

	integer iFirstPositionIndex1 = pNode1->iGetFirstPositionIndex() + 3;
	integer iFirstVelocityIndex1 = iFirstPositionIndex1 + 6;

	integer iFirstPositionIndex2 = pNode2->iGetFirstPositionIndex() + 3;
	integer iFirstVelocityIndex2 = iFirstPositionIndex2 + 6;

	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iFirstPositionIndex1 + iCnt);
		WM.PutRowIndex(3 + iCnt, iFirstVelocityIndex1 + iCnt);

		WM.PutColIndex(iCnt, iFirstPositionIndex1 + iCnt);
		WM.PutColIndex(3 + iCnt, iFirstVelocityIndex1 + iCnt);

		WM.PutRowIndex(6 + iCnt, iFirstPositionIndex2 + iCnt);
		WM.PutRowIndex(9 + iCnt, iFirstVelocityIndex2 + iCnt);

		WM.PutColIndex(6 + iCnt, iFirstPositionIndex2 + iCnt);
		WM.PutColIndex(9 + iCnt, iFirstVelocityIndex2 + iCnt);
	}

	/* Dati */
	Vec3 TmpArm1(pNode1->GetRRef()*Arm1);
	Vec3 TmpArm2(pNode2->GetRRef()*Arm2);
	Vec3 TmpDir = f.Get();
	Vec3 Omega1(pNode1->GetWRef());
	Vec3 Omega2(pNode2->GetWRef());

	/* |    F/\   |           |   F  |
	 * |          | Delta_g = |      |
	 * | (d/\F)/\ |           | d/\F |
	 */

	Mat3x3 MTmp(MatCrossCross, TmpDir, TmpArm1);
	WM.Sub(1, 1, MTmp);
	WM.Sub(4, 1, Mat3x3(MatCrossCross, TmpDir, Omega1)*Mat3x3(MatCross, TmpArm1));
	WM.Sub(4, 4, MTmp);

	MTmp = Mat3x3(MatCrossCross, TmpDir, TmpArm2);
	WM.Add(7, 7, MTmp);
	WM.Add(10, 7, Mat3x3(MatCrossCross, TmpDir, Omega2)*Mat3x3(MatCross, TmpArm2));
	WM.Add(10, 10, MTmp);

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
AbsoluteInternalForce::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUT("Entering AbsoluteInternalForce::InitialAssRes()"
		<< std::endl);

	integer iNumRows;
	integer iNumCols;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici delle incognite del nodo */
	integer iFirstPositionIndex1 = pNode1->iGetFirstPositionIndex();
	integer iFirstVelocityIndex1 = iFirstPositionIndex1 + 6;

	integer iFirstPositionIndex2 = pNode2->iGetFirstPositionIndex();
	integer iFirstVelocityIndex2 = iFirstPositionIndex2 + 6;

	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstPositionIndex1 + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iFirstVelocityIndex1 + iCnt);

		WorkVec.PutRowIndex(12 + iCnt, iFirstPositionIndex2 + iCnt);
		WorkVec.PutRowIndex(18 + iCnt, iFirstVelocityIndex2 + iCnt);
	}

	/* Dati */
	Vec3 TmpDir(f.Get());
	Vec3 TmpArm1(pNode1->GetRCurr()*Arm1);
	Vec3 TmpArm2(pNode2->GetRCurr()*Arm2);
	const Vec3& Omega1(pNode1->GetWCurr());
	const Vec3& Omega2(pNode2->GetWCurr());

	WorkVec.Add(1, TmpDir);
	WorkVec.Add(4, TmpArm1.Cross(TmpDir));

	/* In 7 non c'e' nulla */
	WorkVec.Add(10, (Omega1.Cross(TmpArm1)).Cross(TmpDir));

	WorkVec.Sub(7, TmpDir);
	WorkVec.Sub(10, TmpArm2.Cross(TmpDir));

	/* In 7 non c'e' nulla */
	WorkVec.Sub(16, (Omega2.Cross(TmpArm2)).Cross(TmpDir));

	return WorkVec;
}

/* AbsoluteInternalForce - end */


/* FollowerInternalForce - begin */

/* Costruttore non banale */

FollowerInternalForce::FollowerInternalForce(unsigned int uL,
	const StructNode* pN1, const StructNode* pN2,
	const TplDriveCaller<Vec3>* pDC,
	const Vec3& TmpArm1, const Vec3& TmpArm2,
	flag fOut)
: Elem(uL, fOut),
StructuralInternalForce(uL, pN1, pN2, pDC, fOut),
Arm1(TmpArm1), Arm2(TmpArm2)
{
	NO_OP;
}


FollowerInternalForce::~FollowerInternalForce(void)
{
	NO_OP;
}


void
FollowerInternalForce::WorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 12;
	*piNumCols = 6;
}


void
FollowerInternalForce::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 24;
	*piNumCols = 12;
}


/* Contributo al file di restart */
std::ostream&
FollowerInternalForce::Restart(std::ostream& out) const
{
	Force::Restart(out) << ", follower internal, "
		<< pNode1->GetLabel()
		<< ", position, reference, node, ",
		Arm1.Write(out, ", ") << ", "
		<< pNode2->GetLabel()
		<< ", position, reference, node, ",
		Arm2.Write(out, ", ") << ", ";
	return f.pGetDriveCaller()->Restart(out) << ';' << std::endl;
}


VariableSubMatrixHandler&
FollowerInternalForce::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering FollowerInternalForce::AssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	integer iFirstRotationIndex1 = pNode1->iGetFirstPositionIndex() + 3;
	integer iFirstMomentumIndex1 = pNode1->iGetFirstMomentumIndex();
	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iFirstMomentumIndex1 + iCnt);     /* forza */
		WM.PutRowIndex(3 + iCnt, iFirstMomentumIndex1 + 3 + iCnt); /* coppia */
		WM.PutColIndex(iCnt, iFirstRotationIndex1 + iCnt);     /* rotazione */

		WM.PutRowIndex(6 + iCnt, iFirstMomentumIndex1 + iCnt);   /* forza */
		WM.PutRowIndex(9 + iCnt, iFirstMomentumIndex1 + 3 + iCnt); /* coppia */
		WM.PutColIndex(3 + iCnt, iFirstRotationIndex1 + iCnt);   /* rotazione */
	}

	/* Dati */
	Vec3 TmpDir(pNode1->GetRRef()*(f.Get()*dCoef));
	Vec3 TmpArm1(pNode1->GetRRef()*Arm1);
	Vec3 TmpArm2(pNode2->GetRRef()*Arm2);

	/* |    F/\       0    |             |   F   |
	 * |                   | Delta_g_1 = |       |
	 * | (d1/\F)/\    0    |             | d1/\F |
	 * |                   | Delta_g_2 = |       |
	 * | -F/\d2/\  d2/\F/\ |             | d2/\F |
	 */

	WM.Add(1, 1, Mat3x3(MatCross, TmpDir));
	WM.Add(4, 1, Mat3x3(MatCross, TmpArm1.Cross(TmpDir)));
	WM.Sub(7, 1, Mat3x3(MatCross, TmpDir));
	WM.Sub(7, 1, Mat3x3(MatCrossCross, TmpArm2, TmpDir));
	WM.Add(7, 4, Mat3x3(MatCrossCross, TmpDir, TmpArm2));

	return WorkMat;
}


/* Assembla il residuo */
SubVectorHandler&
FollowerInternalForce::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering FollowerInternalForce::AssRes()" << std::endl);

	integer iNumRows;
	integer iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici delle incognite del nodo */
	integer iFirstMomentumIndex1 = pNode1->iGetFirstMomentumIndex();
	integer iFirstMomentumIndex2 = pNode2->iGetFirstMomentumIndex();
	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex1 + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iFirstMomentumIndex2 + iCnt);
	}

	/* Dati */
	Vec3 TmpDir(f.Get());
	Vec3 F(pNode1->GetRCurr()*TmpDir);
	Vec3 M1(pNode1->GetRCurr()*Arm1.Cross(TmpDir));
	Vec3 M2(F.Cross(pNode2->GetRCurr()*Arm2));

	WorkVec.Add(1, F);
	WorkVec.Add(4, M1);
	WorkVec.Sub(7, F);
	WorkVec.Add(10, M2);

	return WorkVec;
}

/* Inverse Dynamics*/
SubVectorHandler&
FollowerInternalForce::AssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ ,
	const VectorHandler& /* XPrimePrimeCurr */ ,
	int iOrder)
{
	DEBUGCOUT("Entering FollowerInternalForce::AssRes()" << std::endl);

	integer iNumRows;
	integer iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici delle incognite del nodo */
	integer iFirstMomentumIndex1 = pNode1->iGetFirstPositionIndex();
	integer iFirstMomentumIndex2 = pNode2->iGetFirstPositionIndex();
	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex1 + iCnt);
		WorkVec.PutRowIndex(6+iCnt, iFirstMomentumIndex2 + iCnt);
	}

	/* Dati */
	Vec3 TmpDir(f.Get());
	Vec3 F(pNode1->GetRCurr()*TmpDir);
	Vec3 M1(pNode1->GetRCurr()*Arm1.Cross(TmpDir));
	Vec3 M2(F.Cross(pNode2->GetRCurr()*Arm2));

	WorkVec.Add(1, F);
	WorkVec.Add(4, M1);
	WorkVec.Sub(7, F);
	WorkVec.Add(10, M2);

	return WorkVec;
}


void
FollowerInternalForce::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		Vec3 F(pNode1->GetRCurr()*f.Get());
		OH.Forces()
			<< GetLabel()
			<< " " << pNode1->GetLabel()
			<< " " << F
			<< " " << pNode1->GetXCurr() + pNode1->GetRCurr()*Arm1
			<< " " << pNode2->GetLabel()
			<< " " << -F
			<< " " << pNode2->GetXCurr() + pNode2->GetRCurr()*Arm2
			<< std::endl;
	}
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
FollowerInternalForce::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUT("Entering FollowerInternalForce::InitialAssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	integer iNumRows;
	integer iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);

	/* Dimensiona e resetta la matrice di lavoro */
	WM.ResizeReset(iNumRows, iNumCols);

	integer iFirstPositionIndex1 = pNode1->iGetFirstPositionIndex();
	integer iFirstVelocityIndex1 = iFirstPositionIndex1 + 6;

	integer iFirstPositionIndex2 = pNode2->iGetFirstPositionIndex();
	integer iFirstVelocityIndex2 = iFirstPositionIndex2 + 6;

	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iFirstPositionIndex1 + iCnt);
		WM.PutRowIndex(3 + iCnt, iFirstPositionIndex1 + 3 + iCnt);
		WM.PutRowIndex(6 + iCnt, iFirstVelocityIndex1 + iCnt);
		WM.PutRowIndex(9 + iCnt, iFirstVelocityIndex1 + 3 + iCnt);

		WM.PutColIndex(iCnt, iFirstPositionIndex1 +3+iCnt);
		WM.PutColIndex(3 + iCnt, iFirstVelocityIndex1 + 3 + iCnt);

		WM.PutRowIndex(12 + iCnt, iFirstPositionIndex2 + iCnt);
		WM.PutRowIndex(15 + iCnt, iFirstPositionIndex2 + 3 + iCnt);
		WM.PutRowIndex(18 + iCnt, iFirstVelocityIndex2 + iCnt);
		WM.PutRowIndex(21 + iCnt, iFirstVelocityIndex2 + 3 + iCnt);

		WM.PutColIndex(12 + iCnt, iFirstPositionIndex2 + 3 + iCnt);
		WM.PutColIndex(15 + iCnt, iFirstVelocityIndex2 + 3 + iCnt);
	}

	/* Dati */
	Vec3 TmpArm1(pNode1->GetRRef()*Arm1);
	Vec3 TmpArm2(pNode2->GetRRef()*Arm2);
	Vec3 TmpDir = pNode1->GetRRef()*f.Get();
	const Vec3& Omega1(pNode1->GetWRef());
	const Vec3& Omega2(pNode2->GetWRef());

	WM.Add(1, 1, Mat3x3(MatCross, TmpDir));
	WM.Add(4, 1, Mat3x3(MatCross, TmpArm1.Cross(TmpDir)));

	WM.Add(7, 1, Mat3x3(MatCrossCross, Omega1, TmpDir));
	WM.Add(7, 4, Mat3x3(MatCross, TmpDir));
	WM.Add(10, 1, Mat3x3(MatCrossCross, Omega1, TmpArm1.Cross(TmpDir)));
	WM.Add(10, 4, Mat3x3(MatCross, TmpArm1.Cross(TmpDir)));

	WM.Sub(13, 1, Mat3x3(MatCross, TmpDir));
	WM.Add(16, 1, Mat3x3(MatCrossCross, TmpArm2, TmpDir));
	WM.Sub(16, 7, Mat3x3(MatCrossCross, TmpDir, TmpArm2));

	WM.Sub(19, 1, Mat3x3(MatCrossCross, Omega1, TmpDir));
	WM.Sub(19, 4, Mat3x3(MatCross, TmpDir));

	WM.Add(22, 1, Mat3x3(MatCrossCross, TmpArm2, Omega1)*Mat3x3(MatCross, TmpDir)
		- Mat3x3(MatCrossCross, Omega2.Cross(TmpArm2), TmpDir));
	WM.Add(22, 4, Mat3x3(MatCrossCross, TmpArm2, TmpDir));
	WM.Add(22, 7, Mat3x3(MatCrossCross, TmpDir, Omega2)*Mat3x3(MatCross, TmpArm2)
		- Mat3x3(MatCrossCross, Omega1.Cross(TmpDir), TmpArm2));
	WM.Add(22, 10, Mat3x3(MatCrossCross, TmpDir, TmpArm2));

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
FollowerInternalForce::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUT("Entering FollowerInternalForce::InitialAssRes()" << std::endl);

	integer iNumRows;
	integer iNumCols;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici delle incognite del nodo */
	integer iFirstPositionIndex1 = pNode1->iGetFirstPositionIndex();
	integer iFirstVelocityIndex1 = iFirstPositionIndex1 + 6;

	integer iFirstPositionIndex2 = pNode2->iGetFirstPositionIndex();
	integer iFirstVelocityIndex2 = iFirstPositionIndex2 + 6;

	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstPositionIndex1 + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iFirstVelocityIndex1 + iCnt);

		WorkVec.PutRowIndex(12 + iCnt, iFirstPositionIndex2 + iCnt);
		WorkVec.PutRowIndex(18 + iCnt, iFirstVelocityIndex2 + iCnt);
	}

	/* Dati */
	Vec3 TmpDir(pNode1->GetRCurr()*f.Get());
	Vec3 TmpArm1(pNode1->GetRCurr()*Arm1);
	Vec3 TmpArm2(pNode2->GetRCurr()*Arm2);
	const Vec3& Omega1(pNode1->GetWCurr());
	const Vec3& Omega2(pNode2->GetWCurr());

	WorkVec.Add(1, TmpDir);
	WorkVec.Add(4, TmpArm1.Cross(TmpDir));
	WorkVec.Add(7, Omega1.Cross(TmpDir));
	WorkVec.Add(10, (Omega1.Cross(TmpArm1)).Cross(TmpDir)
		+ TmpArm1.Cross(Omega1.Cross(TmpDir)));

	WorkVec.Sub(13, TmpDir);
	WorkVec.Sub(16, TmpArm2.Cross(TmpDir));
	WorkVec.Sub(19, Omega1.Cross(TmpDir));
	WorkVec.Sub(22, (Omega2.Cross(TmpArm2)).Cross(TmpDir)
		+ TmpArm2.Cross(Omega1.Cross(TmpDir)));

	return WorkVec;
}

/* FollowerInternalForce - end */


/* AbsoluteInternalCouple - begin */

/* Costruttore non banale */

AbsoluteInternalCouple::AbsoluteInternalCouple(unsigned int uL,
	const StructNode* pN1, const StructNode* pN2,
	const TplDriveCaller<Vec3>* pDC,
	flag fOut)
: Elem(uL, fOut),
StructuralInternalForce(uL, pN1, pN2, pDC, fOut)
{
	NO_OP;
}


AbsoluteInternalCouple::~AbsoluteInternalCouple(void)
{
	NO_OP;
}

void
AbsoluteInternalCouple::WorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 6;
	*piNumCols = 1;
}


void
AbsoluteInternalCouple::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 6;
	*piNumCols = 1;
}


/* Contributo al file di restart */
std::ostream&
AbsoluteInternalCouple::Restart(std::ostream& out) const
{
	out << "  couple: " << GetLabel() << ", absolute internal, "
		<< pNode1->GetLabel() << ", "
		<< pNode2->GetLabel() << ", ";
	return f.pGetDriveCaller()->Restart(out) << ';' << std::endl;
}


/* Assembla il residuo */
SubVectorHandler&
AbsoluteInternalCouple::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering AbsoluteInternalCouple::AssRes()" << std::endl);

	integer iNumRows;
	integer iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici delle incognite del nodo */
	integer iFirstMomentumIndex1 = pNode1->iGetFirstMomentumIndex() + 3;
	integer iFirstMomentumIndex2 = pNode2->iGetFirstMomentumIndex() + 3;
	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex1 + iCnt);
		WorkVec.PutRowIndex(3 + iCnt, iFirstMomentumIndex2 + iCnt);
	}

	/* Dati */
	Vec3 F(f.Get());

	WorkVec.Add(1, F);
	WorkVec.Sub(4, F);

	return WorkVec;
}

/* Inverse Dynamics*/
SubVectorHandler&
AbsoluteInternalCouple::AssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ ,
	const VectorHandler& /* XPrimePrimeCurr */ ,
	int iOrder)
{
	DEBUGCOUT("Entering AbsoluteInternalCouple::AssRes()" << std::endl);

	integer iNumRows;
	integer iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici delle incognite del nodo */
	integer iFirstMomentumIndex1 = pNode1->iGetFirstPositionIndex() + 3;
	integer iFirstMomentumIndex2 = pNode2->iGetFirstPositionIndex() + 3;
	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex1 + iCnt);
		WorkVec.PutRowIndex(3 + iCnt, iFirstMomentumIndex2 + iCnt);
	}

	/* Dati */
	Vec3 F(f.Get());

	WorkVec.Add(1, F);
	WorkVec.Sub(4, F);

	return WorkVec;
}


void
AbsoluteInternalCouple::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		Vec3 F(f.Get());
		OH.Forces()
			<< GetLabel()
			<< " " << pNode1->GetLabel()
			<< " " << F
			<< " " << pNode2->GetLabel()
			<< " " << -F
			<< std::endl;
	}
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
AbsoluteInternalCouple::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUT("Entering AbsoluteInternalCouple::InitialAssRes()" << std::endl);

	WorkVec.ResizeReset(6);

	/* Indici delle incognite del nodo */
	integer iFirstPositionIndex1 = pNode1->iGetFirstPositionIndex() + 3;
	integer iFirstPositionIndex2 = pNode2->iGetFirstPositionIndex() + 3;
	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstPositionIndex1 + iCnt);
		WorkVec.PutRowIndex(3 + iCnt, iFirstPositionIndex2 + iCnt);
	}

	/* Dati */
	Vec3 F(f.Get());

	WorkVec.Add(1, F);
	WorkVec.Sub(4, F);

	return WorkVec;
}

/* AbsoluteInternalCouple - end */


/* FollowerInternalCouple - begin */

FollowerInternalCouple::FollowerInternalCouple(unsigned int uL,
	const StructNode* pN1, const StructNode* pN2,
	const TplDriveCaller<Vec3>* pDC,
	flag fOut)
: Elem(uL, fOut),
StructuralInternalForce(uL, pN1, pN2, pDC, fOut)
{
	NO_OP;
}


FollowerInternalCouple::~FollowerInternalCouple(void)
{
	NO_OP;
}

void
FollowerInternalCouple::WorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 6;
	*piNumCols = 3;
}


void
FollowerInternalCouple::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 12;
	*piNumCols = 12;
}


/* Contributo al file di restart */
std::ostream&
FollowerInternalCouple::Restart(std::ostream& out) const
{
	out << "  couple: " << GetLabel() << ", follower internal, "
		<< pNode1->GetLabel() << ", "
		<< pNode2->GetLabel() << ", ";
	return f.pGetDriveCaller()->Restart(out) << ';' << std::endl;
}


VariableSubMatrixHandler&
FollowerInternalCouple::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering FollowerInternalCouple::AssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	integer iFirstRotationIndex1 = pNode1->iGetFirstPositionIndex() + 3;
	integer iFirstMomentumIndex1 = pNode1->iGetFirstMomentumIndex() + 3;
	integer iFirstMomentumIndex2 = pNode2->iGetFirstMomentumIndex() + 3;
	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iFirstMomentumIndex1 + iCnt);    /* coppia */
		WM.PutColIndex(iCnt, iFirstRotationIndex1 + iCnt);    /* rotazione */

		WM.PutRowIndex(3 + iCnt, iFirstMomentumIndex2 + iCnt);    /* coppia */
	}

	/* Dati */
	Mat3x3 MWedge(MatCross, pNode1->GetRRef()*(f.Get()*dCoef));

	/* | M /\| Delta_g = | M | */

	WM.Add(1, 1, MWedge);
	WM.Sub(4, 1, MWedge);

	return WorkMat;
}


/* Assembla il residuo */
SubVectorHandler&
FollowerInternalCouple::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering FollowerInternalCouple::AssRes()" << std::endl);

	integer iNumRows;
	integer iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici delle incognite del nodo */
	integer iFirstMomentumIndex1 = pNode1->iGetFirstMomentumIndex() + 3;
	integer iFirstMomentumIndex2 = pNode2->iGetFirstMomentumIndex() + 3;
	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex1 + iCnt);
		WorkVec.PutRowIndex(3 + iCnt, iFirstMomentumIndex2 + iCnt);
	}

	/* Dati */
	Vec3 M(pNode1->GetRCurr()*f.Get());

	WorkVec.Add(1, M);
	WorkVec.Sub(4, M);

	return WorkVec;
}

/* Inverse Dynamics*/
SubVectorHandler&
FollowerInternalCouple::AssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ ,
	const VectorHandler& /* XPrimePrimeCurr */ ,
	int iOrder)
{
	DEBUGCOUT("Entering FollowerInternalCouple::AssRes()" << std::endl);

	integer iNumRows;
	integer iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici delle incognite del nodo */
	integer iFirstMomentumIndex1 = pNode1->iGetFirstPositionIndex() + 3;
	integer iFirstMomentumIndex2 = pNode2->iGetFirstPositionIndex() + 3;
	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex1 + iCnt);
		WorkVec.PutRowIndex(3 + iCnt, iFirstMomentumIndex2 + iCnt);
	}

	/* Dati */
	Vec3 M(pNode1->GetRCurr()*f.Get());

	WorkVec.Add(1, M);
	WorkVec.Sub(4, M);

	return WorkVec;
}

void
FollowerInternalCouple::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		Vec3 F(pNode1->GetRCurr()*f.Get());
		OH.Forces()
			<< GetLabel()
			<< " " << pNode1->GetLabel()
			<< " " << F
			<< " " << pNode2->GetLabel()
			<< " " << -F
			<< std::endl;
	}
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
FollowerInternalCouple::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUT("Entering FollowerInternalCouple::InitialAssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	WM.ResizeReset(12, 6);

	integer iFirstPositionIndex1 = pNode1->iGetFirstPositionIndex() + 3;
	integer iFirstVelocityIndex1 = iFirstPositionIndex1 + 6;
	integer iFirstPositionIndex2 = pNode2->iGetFirstPositionIndex() + 3;
	integer iFirstVelocityIndex2 = iFirstPositionIndex2 + 6;
	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iFirstPositionIndex1 + iCnt);
		WM.PutRowIndex(3 + iCnt, iFirstVelocityIndex1 + iCnt);

		WM.PutColIndex(iCnt, iFirstPositionIndex1 + iCnt);
		WM.PutColIndex(3 + iCnt, iFirstVelocityIndex1 + iCnt);

		WM.PutRowIndex(6 + iCnt, iFirstPositionIndex2 + iCnt);
		WM.PutRowIndex(9 + iCnt, iFirstVelocityIndex2 + iCnt);
	}

	/* Dati */
	Vec3 TmpDir(pNode1->GetRRef()*f.Get());
	const Vec3& Omega1(pNode1->GetWRef());

	/* |    F/\   |           |   F  |
	 * |          | Delta_g = |      |
	 * | (d/\F)/\ |           | d/\F |
	 */

	WM.Add(1, 1, Mat3x3(MatCross, TmpDir));
	WM.Add(4, 1, Mat3x3(MatCrossCross, Omega1, TmpDir));
	WM.Add(4, 4, Mat3x3(MatCross, TmpDir));

	WM.Sub(7, 1, Mat3x3(MatCross, TmpDir));
	WM.Sub(10, 1, Mat3x3(MatCrossCross, Omega1, TmpDir));
	WM.Sub(10, 4, Mat3x3(MatCross, TmpDir));

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
FollowerInternalCouple::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUT("Entering FollowerInternalCouple::InitialAssRes()" << std::endl);

	WorkVec.ResizeReset(12);

	/* Indici delle incognite del nodo */
	integer iFirstPositionIndex1 = pNode1->iGetFirstPositionIndex() + 3;
	integer iFirstVelocityIndex1 = iFirstPositionIndex1 + 6;

	integer iFirstPositionIndex2 = pNode2->iGetFirstPositionIndex() + 3;
	integer iFirstVelocityIndex2 = iFirstPositionIndex2 + 6;

	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstPositionIndex1+iCnt);
		WorkVec.PutRowIndex(3 + iCnt, iFirstVelocityIndex1 + iCnt);

		WorkVec.PutRowIndex(6 + iCnt, iFirstPositionIndex2 + iCnt);
		WorkVec.PutRowIndex(9 + iCnt, iFirstVelocityIndex2 + iCnt);
	}

	/* Dati */
	Vec3 TmpDir(pNode1->GetRCurr()*f.Get());
	const Vec3& Omega1(pNode1->GetWCurr());

	WorkVec.Add(1, TmpDir);
	WorkVec.Add(4, Omega1.Cross(TmpDir));

	WorkVec.Sub(7, TmpDir);
	WorkVec.Sub(10, Omega1.Cross(TmpDir));

	return WorkVec;
}

/* FollowerInternalCouple - end */

Elem *
ReadStructuralForce(DataManager* pDM,
	MBDynParser& HP,
	unsigned int uLabel,
	bool bCouple,
	bool bFollower,
	bool bInternal)
{
	Elem *pEl = 0;
	const char *sType = bCouple ? "Couple" : "Force";

	if (bCouple && bInternal) {
		silent_cerr(sType << "(" << uLabel << ") "
			"line " << HP.GetLineData() << ": "
			"warning, the syntax changed; "
			"you may safely ignore this warning if you used "
			"the syntax documented for MBDyn >= 1.3.7"
			<< std::endl);
	}

	/* nodo collegato */
	StructNode* pNode = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
	ReferenceFrame rf(pNode);
	Vec3 Arm(Zero3);

	// FIXME: legacy...
	Vec3 Dir(Zero3);
	bool bLegacy(false);

	/* distanza dal nodo (vettore di 3 elementi) (solo se e' una forza) */
	if (HP.IsKeyWord("position")) {
		Arm = HP.GetPosRel(rf);
		DEBUGCOUT("Arm is supplied" << std::endl);

	} else {
		if (!bCouple) {
			silent_cerr(sType << "(" << uLabel << "): "
				"\"position\" keyword expected "
				"at line " << HP.GetLineData() << "; "
				"still using deprecated syntax?"
				<< std::endl);

			if (bFollower) {
				try {
					Dir = HP.GetUnitVecRel(rf);
				} catch (ErrNullNorm) {
					silent_cerr(sType << "(" << uLabel << ") has null direction" << std::endl);
					throw ErrNullNorm(MBDYN_EXCEPT_ARGS);
				}

			} else {
				try {
					Dir = HP.GetUnitVecAbs(rf);
				} catch (ErrNullNorm) {
					silent_cerr(sType << "(" << uLabel << ") has null direction" << std::endl);
					throw ErrNullNorm(MBDYN_EXCEPT_ARGS);
				}
			}

			Arm = HP.GetPosRel(rf);

			bLegacy = true;
		}
	}

	StructNode *pNode2 = 0;
	Vec3 Arm2(Zero3);
	if (bInternal) {
		/* nodo collegato */
		pNode2 = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
		ReferenceFrame rf2(pNode2);

		/* distanza dal nodo (vettore di 3 elementi) ( solo se e' una forza) */
		if (HP.IsKeyWord("position")) {
			Arm2 = HP.GetPosRel(rf2);
			DEBUGCOUT("Node 2 arm is supplied" << std::endl);

		} else if (bLegacy) {
			Arm2 = HP.GetPosRel(rf2);
		}
	}

	TplDriveCaller<Vec3>* pDC = 0;
	if (bLegacy) {
		pDC = DC2TDC(HP.GetDriveCaller(), Dir);

	} else {
		if (bFollower) {
			pDC = ReadDCVecRel(pDM, HP, rf);

		} else {
			pDC = ReadDCVecAbs(pDM, HP, rf);
		}
	}

	flag fOut = pDM->fReadOutput(HP, Elem::FORCE);

	/* Alloca la forza */
	if (!bCouple) {
		if (!bFollower) {
			if (!bInternal) {
				SAFENEWWITHCONSTRUCTOR(pEl,
					AbsoluteForce,
					AbsoluteForce(uLabel, pNode, pDC, Arm, fOut));

			} else {
				SAFENEWWITHCONSTRUCTOR(pEl,
					AbsoluteInternalForce,
					AbsoluteInternalForce(uLabel, pNode, pNode2, pDC, Arm, Arm2, fOut));
			}

		} else {
			if (!bInternal) {
				SAFENEWWITHCONSTRUCTOR(pEl,
					FollowerForce,
					FollowerForce(uLabel, pNode, pDC, Arm, fOut));

			} else {
				SAFENEWWITHCONSTRUCTOR(pEl,
					FollowerInternalForce,
					FollowerInternalForce(uLabel, pNode, pNode2, pDC, Arm, Arm2, fOut));
			}
		}

	} else {
		if (!bFollower) {
			if (!bInternal) {
				SAFENEWWITHCONSTRUCTOR(pEl,
					AbsoluteCouple,
					AbsoluteCouple(uLabel, pNode, pDC, fOut));

			} else {
				SAFENEWWITHCONSTRUCTOR(pEl,
					AbsoluteInternalCouple,
					AbsoluteInternalCouple(uLabel, pNode, pNode2, pDC, fOut));
			}

		} else {
			if (!bInternal) {
				SAFENEWWITHCONSTRUCTOR(pEl,
					FollowerCouple,
					FollowerCouple(uLabel, pNode, pDC, fOut));

			} else {
				SAFENEWWITHCONSTRUCTOR(pEl,
					FollowerInternalCouple,
					FollowerInternalCouple(uLabel, pNode, pNode2, pDC, fOut));
			}
		}
	}

	std::ostream& os = pDM->GetLogFile();

	os << "structural";
	if (bInternal) {
		os << " internal";
	}
	if (bFollower) {
		os << " follower";
	} else {
		os << " absolute";
	}
	if (bCouple) {
		os << " couple";
	} else {
		os << " force";
	}
	os << ": " << uLabel << ' ' << pNode->GetLabel()
		<< ' ' << Arm;

	if (pNode2 != 0) {
		os << ' ' << pNode2->GetLabel() << ' ' << Arm2;
	}

	os << std::endl;

	return pEl;
}

