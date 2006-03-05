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

/* Cerniera pilotata */

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "dataman.h"
#include "drvdisp.h"
#include "hint_impl.h"

/* DriveDisplacementJoint - begin */

/* Costruttore non banale */
DriveDisplacementJoint::DriveDisplacementJoint(unsigned int uL,			      
				 const DofOwner* pDO, 
				 const TplDriveCaller<Vec3>* pDC,
				 const StructNode* pN1, 
				 const StructNode* pN2,
				 const Vec3& f1,
				 const Vec3& f2,
				 flag fOut)
: Elem(uL, Elem::JOINT, fOut), 
Joint(uL, Joint::DRIVEHINGE, pDO, fOut), 
TplDriveOwner<Vec3>(pDC),
pNode1(pN1), pNode2(pN2), f1(f1), f2(f2), 
R1Ref(Eye3),
RRef(Eye3),
f2Ref(0.),
dRef(0.),
F(0.)
{
	ASSERT(pNode1 != NULL);
	ASSERT(pNode2 != NULL);
	ASSERT(pNode1->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pNode2->GetNodeType() == Node::STRUCTURAL);
}

   
/* Distruttore */
DriveDisplacementJoint::~DriveDisplacementJoint(void)
{
	NO_OP;
}


/* Contributo al file di restart */
std::ostream&
DriveDisplacementJoint::Restart(std::ostream& out) const
{
	Joint::Restart(out) << ", drive displacement, "
		<< pNode1->GetLabel() << ", reference, node, ",
		f1.Write(out, ", ") << ", "
		<< pNode2->GetLabel() << ", reference, node, ",
		f2.Write(out, ", ") << ", ",
		pGetDriveCaller()->Restart(out) << ';' << std::endl;
	return out;
}


void
DriveDisplacementJoint::Output(OutputHandler& OH) const
{   
	if (fToBeOutput()) {
		Vec3 d(pNode2->GetXCurr() + pNode2->GetRCurr()*f2
			- pNode1->GetXCurr() - pNode1->GetRCurr()*f1);
		Joint::Output(OH.Joints(), "DriveHinge", GetLabel(),
				pNode1->GetRCurr().Transpose()*F, Zero3, F, Zero3)
			<< " " << d << std::endl;
	}
}

void
DriveDisplacementJoint::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
	if (ph) {
		for (unsigned i = 0; i < ph->size(); i++) {
			Joint::JointHint *pjh = dynamic_cast<Joint::JointHint *>((*ph)[i]);

			if (pjh) {
				if (dynamic_cast<Joint::OffsetHint<1> *>(pjh)) {
					/* TODO */

				} else if (dynamic_cast<Joint::OffsetHint<2> *>(pjh)) {
					/* TODO */

				} else if (dynamic_cast<Joint::ReactionsHint *>(pjh)) {
					/* TODO */
				}
				continue;
			}

			TplDriveHint<Vec3> *pdh = dynamic_cast<TplDriveHint<Vec3> *>((*ph)[i]);

			if (pdh) {
				pedantic_cout("DriveDisplacementJoint(" << uLabel << "): "
					"creating drive from hint[" << i << "]..." << std::endl);

				TplDriveCaller<Vec3> *pDC = pdh->pCreateDrive(pDM);
				if (pDC == 0) {
					silent_cerr("DriveDisplacementJoint(" << uLabel << "): "
						"unable to create drive "
						"after hint #" << i << std::endl);
					throw ErrGeneric();
				}
				
				TplDriveOwner<Vec3>::Set(pDC);
				continue;
			}
		}
	}
}

Hint *
DriveDisplacementJoint::ParseHint(DataManager *pDM, const char *s) const
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
	}

	/* take care of "drive" hint... */
	return SimulationEntity::ParseHint(pDM, s);
}

std::ostream&
DriveDisplacementJoint::DescribeDof(std::ostream& out,
		char *prefix, bool bInitial, int i) const
{
	integer iIndex = iGetFirstIndex();

	if (i >= 0) {
		silent_cerr("DriveDisplacementJoint(" << GetLabel() << "): "
			"DescribeDof(" << i << ") "
			"not implemented yet" << std::endl);
		throw ErrGeneric();
	}

	out
		<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
			"reaction forces [fx,fy,fz]" << std::endl;

	if (bInitial) {
		iIndex += 3;
		out
			<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
				"reaction force derivatives [fPx,fPy,fPz]" << std::endl;
	}

	return out;
}

std::ostream&
DriveDisplacementJoint::DescribeEq(std::ostream& out,
		char *prefix, bool bInitial, int i) const
{
	integer iIndex = iGetFirstIndex();

	if (i >= 0) {
		silent_cerr("DriveDisplacementJoint(" << GetLabel() << "): "
			"DescribeEq(" << i << ") "
			"not implemented yet" << std::endl);
		throw ErrGeneric();
	}

	out
		<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
			"position constraints [px1=px2,py1=py2,pz1=pz2]" << std::endl;

	if (bInitial) {
		iIndex += 3;
		out
			<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
				"velocity constraints [vx1=vx2,vy1=vy2,vz1=vz2]" << std::endl;
	}

	return out;

}
   
/* Dati privati (aggiungere magari le reazioni vincolari) */
unsigned int
DriveDisplacementJoint::iGetNumPrivData(void) const
{
	return 3;
};

unsigned int
DriveDisplacementJoint::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);
	ASSERT(s[0] != '\0');

	if (s[0] != 'p' || s[2] != '\0') {
		return 0;
	}

	switch (s[1]) {
	case 'x':
		return 1;
	case 'y':
		return 2;
	case 'z':
		return 3;
	}

	return 0;
}

doublereal
DriveDisplacementJoint::dGetPrivData(unsigned int i) const
{
	ASSERT(i >= 1 && i <= 3);
	return Get().dGet(i);
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler& 
DriveDisplacementJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef, 
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering DriveDisplacementJoint::AssJac()" << std::endl);

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
	integer iFirstReactionIndex = iGetFirstIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);	
		WM.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(12 + iCnt, iFirstReactionIndex + iCnt);
		WM.PutColIndex(12 + iCnt, iFirstReactionIndex + iCnt);
	}
   
	AssMat(WM, dCoef);

	return WorkMat;
}


void
DriveDisplacementJoint::AfterPredict(VectorHandler& /* X */ ,
		VectorHandler& /* XP */ )
{
	/* Recupera i dati */
	f2Ref = pNode2->GetRRef()*f2;
	dRef = pNode1->GetRRef()*(f1 + Get());
}


void
DriveDisplacementJoint::AssMat(FullSubMatrixHandler& WM, doublereal dCoef)   
{
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		/* node 1 force */
		WM.DecCoef(iCnt, 12 + iCnt, 1.);
		/* node 2 force */
		WM.IncCoef(6 + iCnt, 12 + iCnt, 1.);
		/* node 1 constraint */
		WM.DecCoef(12 + iCnt, iCnt, 1.);
		/* node 2 constraint */
		WM.IncCoef(12 + iCnt, 6 + iCnt, 1.);
	}

	Mat3x3 MTmp(dRef);

	/* node 1 moment */
	WM.Sub(3 + 1, 12 + 1, MTmp);
	/* node 1 constraint */
	WM.Add(12 + 1, 3 + 1, MTmp);

	MTmp = Mat3x3(f2Ref);
	/* node 2 moment */
	WM.Add(9 + 1, 12 + 1, MTmp);
	/* node 2 constraint */
	WM.Sub(12 + 1, 9 + 1, MTmp);

	MTmp = Mat3x3(F, dRef*dCoef);
	/* node 1 moment */
	WM.Sub(3 + 1, 3 + 1, MTmp);
	
	MTmp = Mat3x3(F, f2Ref*dCoef);
	/* node 2 moment */
	WM.Add(9 + 1, 9 + 1, MTmp);
}


/* assemblaggio residuo */
SubVectorHandler& 
DriveDisplacementJoint::AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering DriveDisplacementJoint::AssRes()" << std::endl);   

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera gli indici */
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
	integer iFirstReactionIndex = iGetFirstIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
	}

	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(12 + iCnt, iFirstReactionIndex + iCnt);
	}

	F = Vec3(XCurr, iFirstReactionIndex + 1);

	AssVec(WorkVec, dCoef);

	return WorkVec;
}


void
DriveDisplacementJoint::AssVec(SubVectorHandler& WorkVec, doublereal dCoef)
{   
	Mat3x3 R1(pNode1->GetRCurr());

	Vec3 f2Tmp = pNode2->GetRCurr()*f2;
	Vec3 d = pNode1->GetRCurr()*(f1 + Get());

	WorkVec.Add(1, F);
	WorkVec.Add(3 + 1, d.Cross(F));
	WorkVec.Sub(6 + 1, F);
	WorkVec.Sub(9 + 1, f2Tmp.Cross(F));
	
	if (dCoef != 0.) {
		WorkVec.Sub(12 + 1, (pNode2->GetXCurr() + f2Tmp - pNode1->GetXCurr() - d)/dCoef);
	}
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
DriveDisplacementJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering DriveDisplacementJoint::InitialAssJac()" << std::endl);

	WorkMat.SetNullMatrix();

#if 0
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iFirstReactionIndex = iGetFirstIndex();
	integer iNode1FirstVelIndex = iNode1FirstPosIndex + 6;
	integer iNode2FirstVelIndex = iNode2FirstPosIndex + 6;   
	integer iReactionPrimeIndex = iFirstReactionIndex + 3;   

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNode1FirstVelIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode1FirstVelIndex + iCnt);
		WM.PutRowIndex(12 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutColIndex(12 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutRowIndex(18 + iCnt, iNode2FirstVelIndex + iCnt);
		WM.PutColIndex(18 + iCnt, iNode2FirstVelIndex + iCnt);
		WM.PutRowIndex(24 + iCnt, iFirstReactionIndex + iCnt);
		WM.PutColIndex(24 + iCnt, iFirstReactionIndex + iCnt);
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
DriveDisplacementJoint::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering DriveDisplacementJoint::InitialAssRes()" << std::endl);   

	WorkVec.ResizeReset(0);

#if 0
	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iFirstReactionIndex = iGetFirstIndex();
	integer iNode1FirstVelIndex = iNode1FirstPosIndex + 6;
	integer iNode2FirstVelIndex = iNode2FirstPosIndex + 6;
	integer iReactionPrimeIndex = iFirstReactionIndex + 3;

	/* Setta gli indici del vettore */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode1FirstVelIndex + iCnt);
		WorkVec.PutRowIndex(12 + iCnt, iNode2FirstPosIndex + iCnt);
		WorkVec.PutRowIndex(18 + iCnt, iNode2FirstVelIndex + iCnt);
		WorkVec.PutRowIndex(24 + iCnt, iFirstReactionIndex + iCnt);
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
					   
/* DriveDisplacementJoint - end */
