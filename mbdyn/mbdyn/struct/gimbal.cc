/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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

/* GimbalJoint - begin */

/* Costruttore non banale */
GimbalJoint::GimbalJoint(unsigned int uL,			      
				 const DofOwner* pDO, 
				 const StructNode* pN1, 
				 const StructNode* pN2,
				 const Mat3x3& R1,
				 const Mat3x3& R2,
				 flag fOut)
: Elem(uL, Elem::JOINT, fOut), 
Joint(uL, Joint::GIMBAL, pDO, fOut), 
pNode1(pN1), pNode2(pN2), R1h(R1), R2h(R2), 
M(0.), dTheta(0.), dPhi(0.)
{
	ASSERT(pNode1 != NULL);
	ASSERT(pNode2 != NULL);
	ASSERT(pNode1->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pNode2->GetNodeType() == Node::STRUCTURAL);
}

   
/* Distruttore */
GimbalJoint::~GimbalJoint(void)
{
	NO_OP;
}


/* Contributo al file di restart */
std::ostream&
GimbalJoint::Restart(std::ostream& out) const
{
	Joint::Restart(out) << ", gimbal hinge, "
		<< pNode1->GetLabel() << ", reference, node, 1, ",
		(R1h.GetVec(1)).Write(out, ", ")
		<< ", 2, ", (R1h.GetVec(2)).Write(out, ", ") << ", "
		<< pNode2->GetLabel() << ", reference, node, 1, ",
		(R2h.GetVec(1)).Write(out, ", ")
		<< ", 2, ", (R2h.GetVec(2)).Write(out, ", ") << ';' << std::endl;
	return out;
}


void
GimbalJoint::Output(OutputHandler& OH) const
{   
	if (fToBeOutput()) {
		Mat3x3 Ra(pNode1->GetRCurr());

		Vec3 d(MatR2EulerAngles(pNode1->GetRCurr().Transpose()*pNode2->GetRCurr()));
		Joint::Output(OH.Joints(), "Gimbal", GetLabel(),
				Zero3, M, Zero3, Ra*M) << std::endl;
		/* FIXME */
	}
}


/* assemblaggio jacobiano */
VariableSubMatrixHandler& 
GimbalJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef, 
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering GimbalJoint::AssJac()" << std::endl);

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
GimbalJoint::AfterPredict(VectorHandler& /* X */ ,
		VectorHandler& /* XP */ )
{
	NO_OP;
}


void
GimbalJoint::AssMat(FullSubMatrixHandler& WM, doublereal dCoef)   
{
	Mat3x3 Ra(pNode1->GetRRef()*R1h);
	Mat3x3 Rb(pNode2->GetRRef()*R2h);
	Mat3x3 RaT(Ra.Transpose());

	Mat3x3 ExpTheta(RotManip::Rot(Vec3(0., dTheta, 0.)));
	Mat3x3 ExpPhi(RotManip::Rot(Vec3(dPhi, 0., 0.)));

	doublereal dCosTheta = cos(dTheta);
	doublereal dSinTheta = sin(dTheta);
	doublereal dCosPhi = cos(dPhi);
	doublereal dSinPhi = sin(dPhi);

	Mat3x3 MTmp(RaT*dCoef);

	/* equazioni di vincolo */
	/* termini in Delta g_a, Delta g_b */
	WM.Add(6 + 1, 1, MTmp);
	WM.Sub(6 + 1, 3 + 1, MTmp);

	/* termini in Delta theta */
	WM.IncCoef(6 + 1, 9 + 1, dSinTheta*dSinPhi);
	WM.IncCoef(6 + 2, 9 + 1, 2.);
	WM.IncCoef(6 + 3, 9 + 1, (dCosTheta - 1.)*dSinPhi);

	/* termini in Delta phi */
	WM.IncCoef(6 + 1, 9 + 2, dCosTheta);
	WM.DecCoef(6 + 3, 9 + 2, dSinTheta);

	/* coppie */
	/* termini in Delta lambda */
	WM.Sub(1, 6 + 1, Ra);
	WM.Add(3 + 1, 6 + 1, Ra);

	/* termini in Delta g_a, Delta g_b */
	MTmp = Mat3x3(Ra*(M*dCoef));
	WM.Add(1, 1, MTmp);
	WM.Sub(3 + 1, 1, MTmp);

	/* equazione in theta */
	/* termine in Delta theta */
	WM.IncCoef(9 + 1, 9 + 1, dSinPhi*(dCosTheta*M(1) - dSinTheta*M(3)));

	/* termine in Delta phi */
	WM.IncCoef(9 + 1, 9 + 2, dCosPhi*(dSinTheta*M(1) + (dCosTheta - 1.)*M(3)));

	/* termini in Delta lambda */
	WM.IncCoef(9 + 1, 6 + 1, dSinTheta*dSinPhi);
	WM.IncCoef(9 + 1, 6 + 2, 2.);
	WM.IncCoef(9 + 1, 6 + 3, (dCosTheta - 1.)*dSinPhi);

	/* equazione in phi */
	/* termine in Delta theta */
	WM.DecCoef(9 + 2, 9 + 1, dSinTheta*M(1) + dCosTheta*M(3));

	/* termini in Delta lambda */
	WM.IncCoef(9 + 2, 6 + 1, dCosTheta);
	WM.DecCoef(9 + 2, 6 + 3, dSinTheta);
}


/* assemblaggio residuo */
SubVectorHandler& 
GimbalJoint::AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering GimbalJoint::AssRes()" << std::endl);   

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
GimbalJoint::AssVec(SubVectorHandler& WorkVec, doublereal dCoef)
{
	Mat3x3 Ra(pNode1->GetRCurr()*R1h);
	Mat3x3 Rb(pNode2->GetRCurr()*R2h);

	Mat3x3 ExpTheta(RotManip::Rot(Vec3(0., dTheta, 0.)));
	Mat3x3 ExpPhi(RotManip::Rot(Vec3(dPhi, 0., 0.)));

	doublereal dCosTheta = cos(dTheta);
	doublereal dSinTheta = sin(dTheta);
	doublereal dSinPhi = sin(dPhi);

	Vec3 MTmp(Ra*M);

	WorkVec.Add(1, MTmp);
	WorkVec.Sub(3 + 1, MTmp);

	WorkVec.Add(6 + 1, RotManip::VecRot(Ra.Transpose()*Rb) - RotManip::VecRot(ExpTheta*ExpPhi*ExpTheta));

	WorkVec.DecCoef(9 + 1, Vec3(dSinTheta*dSinPhi, 2., (dCosTheta - 1.)*dSinPhi)*M);
	WorkVec.DecCoef(9 + 2, Vec3(dCosTheta, 0., -dSinTheta)*M);
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
GimbalJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering GimbalJoint::InitialAssJac()" << std::endl);

	WorkMat.SetNullMatrix();

#if 0
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex()+3;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex()+3;
	integer iFirstReactionIndex = iGetFirstIndex();
	integer iNode1FirstVelIndex = iNode1FirstPosIndex+6;   
	integer iNode2FirstVelIndex = iNode2FirstPosIndex+6;   
	integer iReactionPrimeIndex = iFirstReactionIndex+3;   

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
		WM.PutRowIndex(3+iCnt, iNode1FirstVelIndex+iCnt);
		WM.PutColIndex(3+iCnt, iNode1FirstVelIndex+iCnt);
		WM.PutRowIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
		WM.PutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
		WM.PutRowIndex(9+iCnt, iNode2FirstVelIndex+iCnt);
		WM.PutColIndex(9+iCnt, iNode2FirstVelIndex+iCnt);
		WM.PutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
		WM.PutColIndex(12+iCnt, iFirstReactionIndex+iCnt);
		WM.PutRowIndex(15+iCnt, iReactionPrimeIndex+iCnt);
		WM.PutColIndex(15+iCnt, iReactionPrimeIndex+iCnt);
	}

	Mat3x3 Ra(pNode1->GetRRef()*R1h);
	Mat3x3 RaT(Ra.Transpose());
	Vec3 Wa(pNode1->GetWRef());
	Vec3 Wb(pNode2->GetWRef());

	Mat3x3 MTmp(M);
	Mat3x3 MPrimeTmp(Ra*Vec3(XCurr, iReactionPrimeIndex+1));

	WM.Add(1, 1, MTmp);
	WM.Add(4, 4, MTmp);   
	WM.Sub(7, 1, MTmp);
	WM.Sub(10, 4, MTmp);

	MTmp = Mat3x3(Wa)*MTmp+MPrimeTmp;
	WM.Add(4, 1, MTmp);
	WM.Sub(10, 1, MTmp);

	WM.Add(7, 13, Ra);
	WM.Add(10, 16, Ra);
	WM.Sub(1, 13, Ra);
	WM.Sub(4, 16, Ra);

	MTmp = Mat3x3(Wa)*Ra;
	WM.Add(10, 13, MTmp);
	WM.Sub(4, 13, MTmp);

	WM.Add(13, 7, RaT);
	WM.Add(16, 10, RaT);
	WM.Sub(13, 1, RaT);
	WM.Sub(16, 4, RaT);
	WM.Add(16, 1, RaT*Mat3x3(Wb));
	WM.Sub(16, 7, RaT*Mat3x3(Wa));
#endif

	return WorkMat;
}

					   
/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
GimbalJoint::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering GimbalJoint::InitialAssRes()" << std::endl);   

	WorkVec.ResizeReset(0);
	
#if 0
	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex()+3;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex()+3;
	integer iFirstReactionIndex = iGetFirstIndex();
	integer iNode1FirstVelIndex = iNode1FirstPosIndex+6;   
	integer iNode2FirstVelIndex = iNode2FirstPosIndex+6;   
	integer iReactionPrimeIndex = iFirstReactionIndex+3;   

	/* Setta gli indici del vettore */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
		WorkVec.PutRowIndex(3+iCnt, iNode1FirstVelIndex+iCnt);
		WorkVec.PutRowIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
		WorkVec.PutRowIndex(9+iCnt, iNode2FirstVelIndex+iCnt);
		WorkVec.PutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
		WorkVec.PutRowIndex(15+iCnt, iReactionPrimeIndex+iCnt);
	}

	Mat3x3 R1(pNode1->GetRCurr()*R1h);
	Mat3x3 R1T(R1.Transpose());
	Vec3 Wa(pNode1->GetWCurr());
	Vec3 Wb(pNode2->GetWCurr());

	M = Vec3(XCurr, iFirstReactionIndex+1);
	Vec3 MPrime = Vec3(XCurr, iReactionPrimeIndex+1);

	Vec3 MPrimeTmp(Wa.Cross(R1*M)+R1*MPrime);

	Mat3x3 R2(pNode2->GetRCurr()*R2h);
	ThetaCurr = RotManip::VecRot(R1.Transpose()*R2);

	Vec3 ThetaPrime = R1T*(Wa.Cross(ThetaCurr)-Wb+Wa);

	WorkVec.Add(1, MPrimeTmp);
	WorkVec.Add(4, MPrimeTmp);
	WorkVec.Sub(7, MPrimeTmp);
	WorkVec.Sub(10, MPrimeTmp);
	WorkVec.Add(13, Get()-ThetaCurr);
	WorkVec.Add(16, ThetaPrime);   
#endif

	return WorkVec;
}
					   
/* GimbalJoint - end */
