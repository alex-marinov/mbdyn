/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <distance.h>

/* DistanceJoint - begin */

/* Costruttore non banale */
DistanceJoint::DistanceJoint(unsigned int uL, const DofOwner* pDO,
		const StructNode* pN1, const StructNode* pN2,
		const DriveCaller* pDC, flag fOut)
: Elem(uL, Elem::JOINT, fOut), 
Joint(uL, Joint::DISTANCE, pDO, fOut),
DriveOwner(pDC),
pNode1(pN1), pNode2(pN2), Vec(0.), dAlpha(0.)
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
	std::cerr << "Joint(" << GetLabel() << "): distance is null"
		<< std::endl;
	THROW(ErrGeneric());
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
		return dGet();
	}

	THROW(ErrGeneric());
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


/* Assemblaggio jacobiano */
VariableSubMatrixHandler& 
DistanceJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering DistanceJoint::AssJac()" << std::endl);
 
	/* Dimensiona e resetta la matrice di lavoro */
	doublereal dDistance = pGetDriveCaller()->dGet();
   
	/* Distanza nulla */
	if (fabs(dDistance) <= DBL_EPSILON) {
		/* abort? */
	}

 	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeInit(7, 7, 0.);
 
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
	integer iFirstReactionIndex = iGetFirstIndex();

	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
		WM.fPutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
		WM.fPutRowIndex(3+iCnt, iNode2FirstMomIndex+iCnt);
		WM.fPutColIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
	}
	WM.fPutRowIndex(6+1, iFirstReactionIndex+1);
	WM.fPutColIndex(6+1, iFirstReactionIndex+1);

	doublereal dd = (dAlpha/dDistance)*dCoef;
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		doublereal d = Vec.dGet(iCnt);

		/* Constraint equation */

		/* C: U^T Delta x_1 */
		WM.fIncCoef(6+1, iCnt, d);

		/* C: - U^T Delta x_2 */
		WM.fDecCoef(6+1, 3+iCnt, d);

		/* Equilibrium */

		/* F1: U Delta lambda */
		WM.fIncCoef(iCnt, 6+1, d);

		/* F2: - U Delta lambda */
		WM.fDecCoef(3+iCnt, 6+1, d);

		/* F1: - lambda/d Delta x_1 */
		WM.fDecCoef(iCnt, iCnt, dd);

		/* F1: lambda/d Delta x_2 */
		WM.fIncCoef(iCnt, 3+iCnt, dd);

		/* F2: lambda/d Delta x_1 */
		WM.fIncCoef(3+iCnt, iCnt, dd);

		/* F2: - lambda/d Delta x_2 */
		WM.fDecCoef(3+iCnt, 3+iCnt, dd);
	}
      
	return WorkMat;
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
	WorkVec.Resize(iNumRows);
	WorkVec.Reset(0.);
        
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();  
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
	integer iFirstReactionIndex = iGetFirstIndex();
   
	for (int iCnt = 1; iCnt <= 3; iCnt++) {      
		/* Indici del nodo 1 */
		WorkVec.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);

		/* Indici del nodo 2 */
		WorkVec.fPutRowIndex(3+iCnt, iNode2FirstMomIndex+iCnt);
	}
	
	/* Indici del vincolo */
	WorkVec.fPutRowIndex(6+1, iFirstReactionIndex+1);

	Vec3 x1(pNode1->GetXCurr());
	Vec3 x2(pNode2->GetXCurr());

	/* Aggiorna i dati propri */
	dAlpha = XCurr.dGetCoef(iFirstReactionIndex+1);

	doublereal dDistance = pGetDriveCaller()->dGet();
   
	/* Distanza nulla */
	if (fabs(dDistance) <= DBL_EPSILON) {	
		Abort();
	}

	Vec = (x2 - x1)/dDistance;
	WorkVec.fIncCoef(6+1, (Vec.Norm()-1.)*dDistance/dCoef);

	Vec3 Tmp(Vec*dAlpha);
	WorkVec.Sub(1, Tmp);
	WorkVec.Add(3+1, Tmp);	

	return WorkVec;
}

void
DistanceJoint::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {      
		doublereal d = dGet();

		Joint::Output(OH.Joints(), "Distance", GetLabel(),
	    			Vec3(dAlpha, 0., 0.), Zero3, Vec*dAlpha, Zero3)
			<< " " << Vec << " " << d << std::endl;
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
	WorkVec.Reset(0.);
	
	return WorkVec;
}

#ifdef USE_ADAMS
void 
DistanceJoint::GetAdamsDummyPart(unsigned int part, 
		Vec3& x, 
		Mat3x3& R) const 
{
	ASSERT(part == 1);
	x = pNode1->GetXCurr();
	R = pNode1->GetRCurr();
}

std::ostream& 
DistanceJoint::WriteAdamsDummyPartCmd(std::ostream& out,
		unsigned int part, 
		unsigned int firstId) const
{
	Vec3 x1 = pNode1->GetXCurr();
	Vec3 x2 = pNode2->GetXCurr();
     
	Vec3 v1 = x2-x1;
	doublereal l = v1.Norm();
	v1 /= l;

	Mat3x3 Rx(Eye3-v1.Tens(v1));
	int index = 1;
	if (fabs(v1.dGet(2)) < fabs(v1.dGet(index))) {
		index = 2;
	}
	if (fabs(v1.dGet(3)) < fabs(v1.dGet(index))) {
		index = 3;
	}

	Vec3 v2(Rx.GetVec(index));
	v2 /= v2.Norm();

	Vec3 e(MatR2EulerAngles(MatR2vec(1, v1, 2, v2))*dRaDegr);

	return out
		<< psAdamsElemCode[GetElemType()] << "_" << GetLabel()
		<< "_" << part << std::endl
		<< firstId << " " << x1 << " "
		<< MatR2EulerAngles(pNode1->GetRCurr())*dRaDegr << " "
		<< x1 << " " << e << " "
		<< l << " " << 0. << " " << 0. << " "
		<< Zero3 << std::endl;
}
#endif /* USE_ADAMS */

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
: Elem(uL, Elem::JOINT, fOut), 
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

/* Assemblaggio jacobiano */
VariableSubMatrixHandler& 
DistanceJointWithOffset::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering DistanceJointWithOffset::AssJac()" << std::endl);
 
	/* Dimensiona e resetta la matrice di lavoro */
	doublereal dDistance = pGetDriveCaller()->dGet();
   
	/* Distanza nulla */
	if (fabs(dDistance) <= DBL_EPSILON) {
		/* abort? */
	}

 	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeInit(13, 13, 0.);
 
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
	integer iFirstReactionIndex = iGetFirstIndex();

	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
		WM.fPutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
		WM.fPutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
		WM.fPutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
	}
	WM.fPutRowIndex(12+1, iFirstReactionIndex+1);
	WM.fPutColIndex(12+1, iFirstReactionIndex+1);

	Vec3 f1Tmp(pNode1->GetRRef()*f1);
	Vec3 f2Tmp(pNode2->GetRRef()*f2);

	Vec3 f1u(f1Tmp.Cross(Vec));
	Vec3 f2u(f2Tmp.Cross(Vec));

	doublereal dd = (dAlpha/dDistance)*dCoef;
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		doublereal d = Vec.dGet(iCnt);

		/* Constraint equation */

		/* C: U^T Delta x_1 */
		WM.fIncCoef(12+1, iCnt, d);

		/* C: - U^T Delta x_2 */
		WM.fDecCoef(12+1, 6+iCnt, d);

		/* Equilibrium */

		/* F1: U Delta lambda */
		WM.fIncCoef(iCnt, 12+1, d);

		/* F2: - U Delta lambda */
		WM.fDecCoef(6+iCnt, 12+1, d);

		/* F1: - lambda/d Delta x_1 */
		WM.fDecCoef(iCnt, iCnt, dd);

		/* F1: lambda/d Delta x_2 */
		WM.fIncCoef(iCnt, 3+iCnt, dd);

		/* F2: lambda/d Delta x_1 */
		WM.fIncCoef(6+iCnt, iCnt, dd);

		/* F2: - lambda/d Delta x_2 */
		WM.fDecCoef(6+iCnt, 6+iCnt, dd);

		d = f1u.dGet(iCnt);

		/* C: - (f1 Cross U)^T Delta g_1 */
		WM.fDecCoef(12+1, 3+iCnt, d);

		/* M1: (f1 Cross u) Delta lambda */
		WM.fIncCoef(3+iCnt, 12+1, d);

		d = f2u.dGet(iCnt);

		/* C: (f2 Cross U)^T Delta g_2 */
		WM.fIncCoef(12+1, 9+iCnt, d);

		/* M2: - (f2 Cross u) Delta lambda */
		WM.fDecCoef(9+iCnt, 12+1, d);
	}

	Mat3x3 Tmp;

	Tmp = Mat3x3(f1Tmp*dd);

	/* F1: lambda/d f1 Cross Delta g_1 */
	WM.Add(1, 3+1, Tmp);

	/* F2: lambda/d f1 Cross Delta g_1 */
	WM.Sub(6+1, 3+1, Tmp);

	/* M1: - lambda/d f1 Cross Delta x_1 */
	WM.Sub(3+1, 1, Tmp);
	
	/* M1: lambda/d f1 Cross Delta x_2 */
	WM.Add(3+1, 6+1, Tmp);
	
	Tmp = Mat3x3(f2Tmp*dd);

	/* F1: lambda/d f2 Cross Delta g_2 */
	WM.Sub(1, 9+1, Tmp);

	/* F2: lambda/d f2 Cross Delta g_2 */
	WM.Add(6+1, 9+1, Tmp);

	/* M1: - lambda/d f2 Cross Delta x_1 */
	WM.Add(9+1, 1, Tmp);
	
	/* M1: lambda/d f2 Cross Delta x_2 */
	WM.Sub(9+1, 6+1, Tmp);

	Tmp = Mat3x3(f1Tmp, f2Tmp*dd);
	
	/* M1: lambda/d f1 Cross f2 Cross Delta g_2 */
	WM.Sub(3+1, 9+1, Tmp);

	Tmp = Mat3x3(f2Tmp, f1Tmp*dd);
	
	/* M2: lambda/d f2 Cross f1 Cross Delta g_1 */
	WM.Add(9+1, 3+1, Tmp);

	Tmp = Mat3x3(Vec+f1Tmp/dDistance, f1Tmp*dd);

	/* M1: - lambda/d (x2 + f2 - x1) Cross f1 Cross Delta g_1 */
	WM.Add(3+1, 3+1, Tmp);
      
	Tmp = Mat3x3(Vec-f2Tmp/dDistance, f2Tmp*dd);

	/* M2: - lambda/d (x2 - x1 - f1) Cross f2 Cross Delta g_2 */
	WM.Sub(9+1, 9+1, Tmp);
      
	return WorkMat;
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
	WorkVec.Resize(iNumRows);
	WorkVec.Reset(0.);
        
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();  
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
	integer iFirstReactionIndex = iGetFirstIndex();
   
	for (int iCnt = 1; iCnt <= 6; iCnt++) {      
		/* Indici del nodo 1 */
		WorkVec.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);

		/* Indici del nodo 2 */
		WorkVec.fPutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
	}
	
	/* Indici del vincolo */
	WorkVec.fPutRowIndex(12+1, iFirstReactionIndex+1);

	Vec3 x1(pNode1->GetXCurr());
	Vec3 x2(pNode2->GetXCurr());

	Vec3 f1Tmp(pNode1->GetRRef()*f1);
	Vec3 f2Tmp(pNode2->GetRRef()*f2);

	/* Aggiorna i dati propri */
	dAlpha = XCurr.dGetCoef(iFirstReactionIndex+1);

	doublereal dDistance = pGetDriveCaller()->dGet();
   
	/* Distanza nulla */
	if (fabs(dDistance) <= DBL_EPSILON) {	
		Abort();
	}

	Vec = (x2 + f2Tmp - x1 - f1Tmp)/dDistance;
	WorkVec.fIncCoef(12+1, (Vec.Norm()-1.)*dDistance/dCoef);

	Vec3 Tmp(Vec*dAlpha);
	WorkVec.Sub(1, Tmp);
	WorkVec.Sub(3+1, f1Tmp.Cross(Tmp));
	WorkVec.Add(6+1, Tmp);	
	WorkVec.Add(9+1, f2Tmp.Cross(Tmp));

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


#ifdef USE_ADAMS
void 
DistanceJointWithOffset::GetAdamsDummyPart(unsigned int part,
		Vec3& x, 
		Mat3x3& R) const 
{
	ASSERT(part == 1);
	x = pNode1->GetXCurr()+pNode1->GetRCurr()*f1;
	
	Vec3 x2 = pNode2->GetXCurr()+pNode2->GetRCurr()*f2;

	Vec3 v1 = x2-x;
	doublereal l = v1.Norm();
	v1 /= l;

	Mat3x3 Rx(Eye3-v1.Tens(v1));
	int index = 1;
	if (fabs(v1.dGet(2)) < fabs(v1.dGet(index))) {
		index = 2;
	}
	if (fabs(v1.dGet(3)) < fabs(v1.dGet(index))) {
		index = 3;
	}

	Vec3 v2(Rx.GetVec(index));
	v2 /= v2.Norm();

	R = MatR2vec(1, v1, 2, v2);

	// R = pNode1->GetRCurr();
}

std::ostream& 
DistanceJointWithOffset::WriteAdamsDummyPartCmd(std::ostream& out,
		unsigned int part, 
		unsigned int firstId) const
{
	Vec3 x1 = pNode1->GetXCurr()+pNode1->GetRCurr()*f1;
	Vec3 x2 = pNode2->GetXCurr()+pNode2->GetRCurr()*f2;

	Vec3 v1 = x2-x1; 
	doublereal l = v1.Norm();
	v1 /= l;

	Mat3x3 Rx(Eye3-v1.Tens(v1));
	int index = 1;
	if (fabs(v1.dGet(2)) < fabs(v1.dGet(index))) {
		index = 2;
	}
	if (fabs(v1.dGet(3)) < fabs(v1.dGet(index))) {
		index = 3;
	}

	Vec3 v2(Rx.GetVec(index));
	v2 /= v2.Norm();

	Vec3 e(MatR2EulerAngles(MatR2vec(1, v1, 2, v2))*dRaDegr);

	return out 
		<< psAdamsElemCode[GetElemType()] << "_" << GetLabel()
		<< "_" << part << std::endl
		<< firstId << " " << x1 << " "
		<< e /* MatR2EulerAngles(pNode1->GetRCurr())*dRaDegr */ << " "
		<< x1 << " "
		<< e << " "
		<< l << " " << 0. << " " << 0. << " "
		<< Zero3 << std::endl;
}
#endif /* USE_ADAMS */

/* DistanceJointWithOffset - end */

