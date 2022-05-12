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

/* Vincoli generali */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <limits>

#ifdef MBDYN_X_DISTANCE_JOINT

#include "distance.h"
#include "hint.h"
#include "hint_impl.h"

/* DistanceJoint - begin */

/* Costruttore non banale */
DistanceJoint::DistanceJoint(unsigned int uL, const DofOwner* pDO,
		const StructDispNode* pN1, const StructDispNode* pN2,
		const DriveCaller* pDC, flag fOut)
: Elem(uL, fOut),
Joint(uL, pDO, fOut),
DriveOwner(pDC),
pNode1(pN1), pNode2(pN2), Vec(Zero3), 
dAlpha(0.)
{
	NO_OP;
}

/* Distruttore banale - ci pensa il DriveOwner a distruggere il DriveCaller */
DistanceJoint::~DistanceJoint(void)
{
	NO_OP;
}

void
DistanceJoint::Abort(void)
{
	silent_cerr("DistanceJoint(" << GetLabel() << "): distance is null"
		<< std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

void
DistanceJoint::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& /* XP */ ,
	SimulationEntity::Hints *ph)
{
	if (ph) {
		for (unsigned i = 0; i < ph->size(); i++) {
			DriveHint *pdh = dynamic_cast<DriveHint *>((*ph)[i]);

			if (pdh) {
				pedantic_cout("DistanceJoint(" << uLabel << "): "
					"creating drive from hint[" << i << "]..." << std::endl);

				DriveCaller *pDC = pdh->pCreateDrive(pDM);
				if (pDC == 0) {
					silent_cerr("DistanceJoint(" << uLabel << "): "
						"unable to create drive after hint "
						"#" << i << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				DriveOwner::Set(pDC);
				continue;
			}
		}
	}

	doublereal dDistance = pGetDriveCaller()->dGet();

	/* Setta a 1 dAlpha, che e' indipendente dalle altre variabili
	 * in caso di distanza nulla */
	if (fabs(dDistance) <= std::numeric_limits<doublereal>::epsilon()) {
		Abort();
	}

	/* Scrive la direzione della distanza. Se e' stata ottenuta con
	 * l'assemblaggio iniziale bene, se no' la calcola */
	Vec = pNode2->GetXCurr() - pNode1->GetXCurr();
	doublereal d = Vec.Dot();
	if (d <= std::numeric_limits<doublereal>::epsilon()) {
    		silent_cerr("DistanceJoint(" << uLabel << ") "
			"linked to nodes " << pNode1->GetLabel()
			<< " and " << pNode2->GetLabel() << ": "
			"nodes are coincident;" << std::endl
	  		<< "initial joint assembly is recommended"
			<< std::endl);
		throw ErrNullNorm(MBDYN_EXCEPT_ARGS);
	}
}

/* Dati privati */
unsigned int
DistanceJoint::iGetNumPrivData(void) const
{
	return 1;
}

unsigned int
DistanceJoint::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);

	if (strcmp(s, "d") == 0) {
		return 1;
	}

	return 0;
}

doublereal
DistanceJoint::dGetPrivData(unsigned int i) const
{
	ASSERT(i == 1);

	if (i == 1) {
#if 0
		/* FIXME: could simply use */
		return dDistance;
#endif
		return dGet();
	}

	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

/* Contributo al file di restart */
std::ostream& DistanceJoint::Restart(std::ostream& out) const
{
	Joint::Restart(out) << ", distance, "
		<< pNode1->GetLabel() << ", "
		<< pNode2->GetLabel() << ", ",
		pGetDriveCaller()->Restart(out) << ';' << std::endl;
	return out;
}


/* Assemblaggio matrici */
void
DistanceJoint::AssMat(FullSubMatrixHandler& WorkMatA,
		 FullSubMatrixHandler& WorkMatB,
		 doublereal dCoef,
		 const VectorHandler& XCurr,
		 const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering DistanceJoint::AssMat()" << std::endl);

	doublereal dd = dAlpha*dCoef;
	doublereal dv = Vec.Norm();
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		doublereal d = Vec.dGet(iCnt);

		/* Constraint equation */

		/* C: - U^T/|U| Delta x_1 */
		WorkMatB.DecCoef(6 + 1, 0 + iCnt, d/dv);

		/* C: U^T/|U| Delta x_2 */
		WorkMatB.IncCoef(6 + 1, 3 + iCnt, d/dv);

		/* Equilibrium */

		/* F1: - U Delta lambda */
		WorkMatB.DecCoef(0 + iCnt, 6 + 1, d);

		/* F2: U Delta lambda */
		WorkMatB.IncCoef(3 + iCnt, 6 + 1, d);

		/* F1: lambda Delta x_1 */
		WorkMatA.IncCoef(0 + iCnt, 0 + iCnt, d*dd);

		/* F1: - lambda Delta x_2 */
		WorkMatA.DecCoef(0 + iCnt, 3 + iCnt, d*dd);

		/* F2: - lambda Delta x_1 */
		WorkMatA.DecCoef(3 + iCnt, 0 + iCnt, d*dd);

		/* F2: lambda Delta x_2 */
		WorkMatA.IncCoef(3 + iCnt, 3 + iCnt, d*dd);
	}
}


/* Assemblaggio jacobiano */
VariableSubMatrixHandler&
DistanceJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering DistanceJoint::AssJac()" << std::endl);

 	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeReset(7, 7);

	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
	integer iFirstReactionIndex = iGetFirstIndex();

	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(0 + iCnt, iNode1FirstMomIndex + iCnt);
		WM.PutColIndex(0 + iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(3 + iCnt, iNode2FirstMomIndex + iCnt);
		WM.PutColIndex(3 + iCnt, iNode2FirstPosIndex + iCnt);
	}
	WM.PutRowIndex(6 + 1, iFirstReactionIndex + 1);
	WM.PutColIndex(6 + 1, iFirstReactionIndex + 1);

	AssMat(WM, WM, dCoef, XCurr, XPrimeCurr);

	return WorkMat;
}


/* Assemblaggio jacobiano */
void
DistanceJoint::AssMats(VariableSubMatrixHandler& WorkMatA,
		VariableSubMatrixHandler& WorkMatB,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering DistanceJoint::AssMats()" << std::endl);

 	FullSubMatrixHandler& WMA = WorkMatA.SetFull();
	WMA.ResizeReset(7, 7);

 	FullSubMatrixHandler& WMB = WorkMatB.SetFull();
	WMB.ResizeReset(7, 7);

	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
	integer iFirstReactionIndex = iGetFirstIndex();

	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WMA.PutRowIndex(0 + iCnt, iNode1FirstMomIndex + iCnt);
		WMA.PutColIndex(0 + iCnt, iNode1FirstPosIndex + iCnt);
		WMA.PutRowIndex(3 + iCnt, iNode2FirstMomIndex + iCnt);
		WMA.PutColIndex(3 + iCnt, iNode2FirstPosIndex + iCnt);

		WMB.PutRowIndex(0 + iCnt, iNode1FirstMomIndex + iCnt);
		WMB.PutColIndex(0 + iCnt, iNode1FirstPosIndex + iCnt);
		WMB.PutRowIndex(3 + iCnt, iNode2FirstMomIndex + iCnt);
		WMB.PutColIndex(3 + iCnt, iNode2FirstPosIndex + iCnt);
	}
	WMA.PutRowIndex(6 + 1, iFirstReactionIndex + 1);
	WMA.PutColIndex(6 + 1, iFirstReactionIndex + 1);

	WMB.PutRowIndex(6 + 1, iFirstReactionIndex + 1);
	WMB.PutColIndex(6 + 1, iFirstReactionIndex + 1);

	AssMat(WMA, WMB, 1., XCurr, XPrimeCurr);
}


/* Assemblaggio residuo */
SubVectorHandler&
DistanceJoint::AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering DistanceJoint::AssRes()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
	integer iFirstReactionIndex = iGetFirstIndex();

	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		/* Indici del nodo 1 */
		WorkVec.PutRowIndex(0 + iCnt, iNode1FirstMomIndex + iCnt);

		/* Indici del nodo 2 */
		WorkVec.PutRowIndex(3 + iCnt, iNode2FirstMomIndex + iCnt);
	}

	/* Indici del vincolo */
	WorkVec.PutRowIndex(6 + 1, iFirstReactionIndex + 1);

	Vec = pNode2->GetXCurr() - pNode1->GetXCurr();

	/* Aggiorna i dati propri */
	dAlpha = XCurr(iFirstReactionIndex + 1);

	dDistance = pGetDriveCaller()->dGet();

	/* Distanza nulla */
	if (fabs(dDistance) <= std::numeric_limits<doublereal>::epsilon()) {
		Abort();
	}

	WorkVec.IncCoef(6 + 1, (dDistance - Vec.Norm())/dCoef);

	Vec3 Tmp(Vec*dAlpha);
	WorkVec.Add(0 + 1, Tmp);
	WorkVec.Sub(3 + 1, Tmp);

	return WorkVec;
}

void
DistanceJoint::OutputPrepare(OutputHandler& OH)
{
	if (bToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseBinary(OutputHandler::JOINTS)) {
			std::string name;
			OutputPrepare_int("distance", OH, name);
			
			Var_V = OH.CreateVar<Vec3>(name + "V", "-",
				"constrained distance direction unit vector (x, y, z)");
			
			Var_d = OH.CreateVar<doublereal>(name + "d", "m",
				"constrained distance magnitude");
		}
#endif // USE_NETCDF
	}
}

void
DistanceJoint::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {
		
		if (OH.UseText(OutputHandler::JOINTS)) {
			Joint::Output(OH.Joints(), "Distance", GetLabel(),
					Vec3(dAlpha, 0., 0.), Zero3, Vec*dAlpha, Zero3)
				<< " " << Vec/dDistance << " " << dDistance
				<< std::endl;
		}

#ifdef USE_NETCDF
		if (OH.UseBinary(OutputHandler::JOINTS)) {
			Joint::NetCDFOutput(OH, Vec3(dAlpha, 0., 0.), Zero3, Vec*dAlpha, Zero3);
			OH.WriteNcVar(Var_V, Vec/dDistance);
			OH.WriteNcVar(Var_d, dDistance);
		}
#endif // USE_NETCDF

	}
}

/* Nota: vanno modificati in analogia al DistanceWithOffsetJoint */

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
DistanceJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering DistanceJoint::InitialAssJac()" << std::endl);

	WorkMat.SetNullMatrix();

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
DistanceJoint::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering DistanceJoint::InitialAssRes()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	WorkVec.Resize(0);

	return WorkVec;
}


/* DistanceJoint - end */


/* DistanceJointWithOffset - begin */

/* Costruttore non banale */
DistanceJointWithOffset::DistanceJointWithOffset(unsigned int uL,
		const DofOwner* pDO,
		const StructNode* pN1,
		const StructNode* pN2,
		const Vec3& f1Tmp,
		const Vec3& f2Tmp,
		const DriveCaller* pDC,
		flag fOut)
: Elem(uL, fOut),
DistanceJoint(uL, pDO, pN1, pN2, pDC, fOut),
f1(f1Tmp), f2(f2Tmp)
{
	NO_OP;
}

/* Distruttore banale - ci pensa il DriveOwner a distruggere il DriveCaller */
DistanceJointWithOffset::~DistanceJointWithOffset(void)
{
	NO_OP;
}

/* Contributo al file di restart */
std::ostream&
DistanceJointWithOffset::Restart(std::ostream& out) const
{
	Joint::Restart(out) << ", distance with offset, "
		<< pNode1->GetLabel()
		<< ", reference, node, ",
		f1.Write(out, ", ") << ", "
		<< pNode2->GetLabel()
		<< ", reference, node, ",
		f2.Write(out, ", ") << ", ";
		return pGetDriveCaller()->Restart(out) << ';' << std::endl;
}

void
DistanceJointWithOffset::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& /* XP */ ,
	SimulationEntity::Hints *ph)
{
	if (ph) {
		for (unsigned i = 0; i < ph->size(); i++) {
			pedantic_cout("DistanceJointWithOffset(" << uLabel << "): "
				"creating drive from hint..." << std::endl);

			DriveHint *pdh = dynamic_cast<DriveHint *>((*ph)[i]);

			if (pdh) {
				DriveCaller *pDC = pdh->pCreateDrive(pDM);
				if (pDC == 0) {
					silent_cerr("DistanceJointWithOffset(" << uLabel << "): "
						"unable to create drive after hint "
						"#" << i << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				DriveOwner::Set(pDC);
				continue;
			}
		}
	}

	doublereal dDistance = pGetDriveCaller()->dGet();

	/* Setta a 1 dAlpha, che e' indipendente dalle altre variabili
	 * in caso di distanza nulla */
	if (fabs(dDistance) <= std::numeric_limits<doublereal>::epsilon()) {
		Abort();
	}

	 Vec = pNode2->GetXCurr() + dynamic_cast<const StructNode *>(pNode2)->GetRCurr()*f2
		- pNode1->GetXCurr() - dynamic_cast<const StructNode *>(pNode1)->GetRCurr()*f1;
	doublereal d = Vec.Dot();
	if (d <= std::numeric_limits<doublereal>::epsilon()) {
		silent_cerr("DistanceJoint(" << GetLabel() << ") "
			"linked to nodes " << pNode1->GetLabel()
			<< " and " << pNode2->GetLabel() << ": "
			"nodes are coincident;" << std::endl
			<< "this is no longer supported" << std::endl);
		throw ErrNullNorm(MBDYN_EXCEPT_ARGS);
	}
}

/* Assemblaggio matrici */
void
DistanceJointWithOffset::AssMat(FullSubMatrixHandler& WorkMatA,
		FullSubMatrixHandler& WorkMatB,
		doublereal dCoef,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering DistanceJointWithOffset::AssMat()" << std::endl);

	Vec3 f1Tmp(dynamic_cast<const StructNode *>(pNode1)->GetRRef()*f1);
	Vec3 f2Tmp(dynamic_cast<const StructNode *>(pNode2)->GetRRef()*f2);

	Vec3 f1u(f1Tmp.Cross(Vec));
	Vec3 f2u(f2Tmp.Cross(Vec));

	doublereal dd = dAlpha*dCoef;
	doublereal dv = Vec.Norm();
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		doublereal d = Vec.dGet(iCnt);

		/* Constraint equation */

		/* C: - U^T/|U| Delta x_1 */
		WorkMatB.DecCoef(12 + 1, 0 + iCnt, d/dv);

		/* C: U^T/|U| Delta x_2 */
		WorkMatB.IncCoef(12 + 1, 6 + iCnt, d/dv);

		/* Equilibrium */

		/* F1: - U Delta lambda */
		WorkMatB.DecCoef(0 + iCnt, 12 + 1, d);

		/* F2: U Delta lambda */
		WorkMatB.IncCoef(6 + iCnt, 12 + 1, d);

		/* F1: lambda Delta x_1 */
		WorkMatA.IncCoef(0 + iCnt, 0 + iCnt, dd);

		/* F1: - lambda Delta x_2 */
		WorkMatA.DecCoef(0 + iCnt, 6 + iCnt, dd);

		/* F2: - lambda Delta x_1 */
		WorkMatA.DecCoef(6 + iCnt, 0 + iCnt, dd);

		/* F2: lambda Delta x_2 */
		WorkMatA.IncCoef(6 + iCnt, 6 + iCnt, dd);

		d = f1u.dGet(iCnt);

		/* C: - (f1 Cross U)^T/|U| Delta g_1 */
		WorkMatB.DecCoef(12 + 1, 3 + iCnt, d/dv);

		/* M1: - (f1 Cross U) Delta lambda */
		WorkMatB.DecCoef(3 + iCnt, 12 + 1, d);

		d = f2u.dGet(iCnt);

		/* C: (f2 Cross U)^T/|U| Delta g_2 */
		WorkMatB.IncCoef(12 + 1, 9 + iCnt, d/dv);

		/* M2: (f2 Cross U) Delta lambda */
		WorkMatB.IncCoef(9 + iCnt, 12 + 1, d);
	}

	Mat3x3 Tmp;

	Tmp = Mat3x3(MatCross, f1Tmp*dd);

	/* F1: - lambda f1 Cross Delta g_1 */
	WorkMatA.Sub(0 + 1, 3 + 1, Tmp);

	/* F2: lambda f1 Cross Delta g_1 */
	WorkMatA.Add(6 + 1, 3 + 1, Tmp);

	/* M1: lambda f1 Cross Delta x_1 */
	WorkMatA.Add(3 + 1, 0 + 1, Tmp);

	/* M1: - lambda f1 Cross Delta x_2 */
	WorkMatA.Sub(3 + 1, 6 + 1, Tmp);

	Tmp = Mat3x3(MatCross, f2Tmp*dd);

	/* F1: lambda f2 Cross Delta g_2 */
	WorkMatA.Add(0 + 1, 9 + 1, Tmp);

	/* F2: - lambda f2 Cross Delta g_2 */
	WorkMatA.Sub(6 + 1, 9 + 1, Tmp);

	/* M2: - lambda f2 Cross Delta x_1 */
	WorkMatA.Sub(9 + 1, 0 + 1, Tmp);

	/* M2: lambda f2 Cross Delta x_2 */
	WorkMatA.Add(9 + 1, 6 + 1, Tmp);

	Tmp = Mat3x3(MatCrossCross, f1Tmp, f2Tmp*dd);

	/* M1: lambda f1 Cross f2 Cross Delta g_2 */
	WorkMatA.Add(3 + 1, 9 + 1, Tmp);

	Tmp = Mat3x3(MatCrossCross, f2Tmp, f1Tmp*dd);

	/* M2: lambda f2 Cross f1 Cross Delta g_1 */
	WorkMatA.Add(9 + 1, 3 + 1, Tmp);

	Tmp = Mat3x3(MatCrossCross, Vec + f1Tmp, f1Tmp*dd);

	/* M1: - lambda (x2 + f2 - x1) Cross f1 Cross Delta g_1 */
	WorkMatA.Sub(3 + 1, 3 + 1, Tmp);

	Tmp = Mat3x3(MatCrossCross, Vec - f2Tmp, f2Tmp*dd);

	/* M2: - lambda (x2 - x1 - f1) Cross f2 Cross Delta g_2 */
	WorkMatA.Add(9 + 1, 9 + 1, Tmp);
}

/* Assemblaggio jacobiano */
VariableSubMatrixHandler&
DistanceJointWithOffset::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering DistanceJointWithOffset::AssJac()" << std::endl);

 	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeReset(13, 13);

	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
	integer iFirstReactionIndex = iGetFirstIndex();

	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(0 + iCnt, iNode1FirstMomIndex + iCnt);
		WM.PutColIndex(0 + iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}
	WM.PutRowIndex(12 + 1, iFirstReactionIndex + 1);
	WM.PutColIndex(12 + 1, iFirstReactionIndex + 1);

	AssMat(WM, WM, dCoef, XCurr, XPrimeCurr);

	return WorkMat;
}

/* Assemblaggio matrici */
void
DistanceJointWithOffset::AssMats(VariableSubMatrixHandler& WorkMatA,
		VariableSubMatrixHandler& WorkMatB,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering DistanceJointWithOffset::AssMats()" << std::endl);

 	FullSubMatrixHandler& WMA = WorkMatA.SetFull();
	WMA.ResizeReset(13, 13);

 	FullSubMatrixHandler& WMB = WorkMatB.SetFull();
	WMB.ResizeReset(13, 13);

	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
	integer iFirstReactionIndex = iGetFirstIndex();

	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WMA.PutRowIndex(0 + iCnt, iNode1FirstMomIndex + iCnt);
		WMA.PutColIndex(0 + iCnt, iNode1FirstPosIndex + iCnt);
		WMA.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
		WMA.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);

		WMB.PutRowIndex(0 + iCnt, iNode1FirstMomIndex + iCnt);
		WMB.PutColIndex(0 + iCnt, iNode1FirstPosIndex + iCnt);
		WMB.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
		WMB.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}
	WMA.PutRowIndex(12 + 1, iFirstReactionIndex + 1);
	WMA.PutColIndex(12 + 1, iFirstReactionIndex + 1);

	WMB.PutRowIndex(12 + 1, iFirstReactionIndex + 1);
	WMB.PutColIndex(12 + 1, iFirstReactionIndex + 1);

	AssMat(WMA, WMB, 1., XCurr, XPrimeCurr);
}

/* Assemblaggio residuo */
SubVectorHandler&
DistanceJointWithOffset::AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering DistanceJointWithOffset::AssRes()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
	integer iFirstReactionIndex = iGetFirstIndex();

	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		/* Indici del nodo 1 */
		WorkVec.PutRowIndex(0 + iCnt, iNode1FirstMomIndex + iCnt);

		/* Indici del nodo 2 */
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
	}

	/* Indici del vincolo */
	WorkVec.PutRowIndex(12 + 1, iFirstReactionIndex + 1);

	Vec3 f1Tmp(dynamic_cast<const StructNode *>(pNode1)->GetRCurr()*f1);
	Vec3 f2Tmp(dynamic_cast<const StructNode *>(pNode2)->GetRCurr()*f2);

	/* x2 + f2 - x1 - f1 */
	Vec = pNode2->GetXCurr() + f2Tmp
		- pNode1->GetXCurr() - f1Tmp;

	/* Aggiorna i dati propri */
	dAlpha = XCurr(iFirstReactionIndex + 1);

	dDistance = pGetDriveCaller()->dGet();

	/* Distanza nulla */
	if (fabs(dDistance) <= std::numeric_limits<doublereal>::epsilon()) {
		Abort();
	}

	WorkVec.IncCoef(12 + 1, (dDistance - Vec.Norm())/dCoef);

	Vec3 Tmp(Vec*dAlpha);
	WorkVec.Add(0 + 1, Tmp);
	WorkVec.Add(3 + 1, f1Tmp.Cross(Tmp));
	WorkVec.Sub(6 + 1, Tmp);
	WorkVec.Sub(9 + 1, f2Tmp.Cross(Tmp));

	return WorkVec;
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
DistanceJointWithOffset::InitialAssJac(VariableSubMatrixHandler& WorkMat,
				       const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering DistanceJointWithOffset::InitialAssJac()"
			<< std::endl);

	WorkMat.SetNullMatrix();

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
DistanceJointWithOffset::InitialAssRes(SubVectorHandler& WorkVec,
				       const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering DistanceJointWithOffset::InitialAssRes()"
			<< std::endl);

	WorkVec.Resize(0);

	return WorkVec;
}

#endif

/* DistanceJointWithOffset - end */
