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

/* Cerniera pilotata */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

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
: Elem(uL, fOut), 
Joint(uL, pDO, fOut), 
TplDriveOwner<Vec3>(pDC),
pNode1(pN1), pNode2(pN2), f1(f1), f2(f2), 
R1Ref(Eye3),
RRef(Eye3),
f2Ref(Zero3),
dRef(Zero3),
F(Zero3)
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
DriveDisplacementJoint::OutputPrepare(OutputHandler &OH)
{
	if (bToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseBinary(OutputHandler::JOINTS)) {
			std::string name;
			OutputPrepare_int("Drive displacement", OH, name);

			Var_d = OH.GetBinaryFile()->CreateVar<Vec3>(name + "d",
				MBUnits::Dimensions::Length,
				"imposed relative displacement, global frame (x, y, z)");
		}
#endif // USE_NETCDF
	}
}

void
DriveDisplacementJoint::Output(OutputHandler& OH) const
{   
	if (bToBeOutput()) {
		Vec3 d(pNode2->GetXCurr() + pNode2->GetRCurr()*f2
			- pNode1->GetXCurr() - pNode1->GetRCurr()*f1);


		if (OH.UseText(OutputHandler::JOINTS)) {
			Joint::Output(OH.Joints(), "DriveDisplacementJoint", GetLabel(),
					pNode1->GetRCurr().Transpose()*F, Zero3, F, Zero3)
				<< " " << d << std::endl;
		}

#ifdef USE_NETCDF
		if (OH.UseBinary(OutputHandler::JOINTS)) {
			Joint::NetCDFOutput(OH, F, Zero3, F, Zero3);
			OH.WriteVar(Var_d, d);
		}
#endif // USE_NETCDF

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
					Mat3x3 R1t(pNode1->GetRCurr().Transpose());
					Vec3 fTmp2(pNode2->GetRCurr()*f2);
   
					f1 = R1t*(pNode2->GetXCurr() + fTmp2 - pNode1->GetXCurr() - Get());

				} else if (dynamic_cast<Joint::OffsetHint<2> *>(pjh)) {
					Mat3x3 R2t(pNode2->GetRCurr().Transpose());
					Vec3 fTmp1(pNode1->GetRCurr()*f1);
   
					f2 = R2t*(pNode1->GetXCurr() + fTmp1 - pNode2->GetXCurr() + Get());

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
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
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
	if (strncasecmp(s, "offset{" /*}*/ , STRLENOF("offset{" /*}*/ )) == 0)
	{
		s += STRLENOF("offset{" /*}*/ );

		if (strcmp(&s[1], /*{*/ "}") != 0) {
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
		const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

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

static const char xyz[] = "xyz";
static const char *dof[] = { "reaction force f", "reaction force derivative fP" };
static const char *eq[] = { "displacement constraint P", "displacement constraint derivative v" };

void
DriveDisplacementJoint::DescribeDof(std::vector<std::string>& desc,
	bool bInitial, int i) const
{
	int iend = 1;
	if (i == -1) {
		if (bInitial) {
			iend = 6;

		} else {
			iend = 3;
		}
	}
	desc.resize(iend);

	std::ostringstream os;
	os << "DriveDisplacementJoint(" << GetLabel() << ")";

	if (i == -1) {
		std::string name = os.str();
		for (i = 0; i < iend; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": " << dof[i/3] << xyz[i%3];

			desc[i] = os.str();
		}

	} else {
		os << ": " << dof[i/3] << xyz[i%3];
		desc[0] = os.str();
	}
}

std::ostream&
DriveDisplacementJoint::DescribeEq(std::ostream& out,
		const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

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
   
void
DriveDisplacementJoint::DescribeEq(std::vector<std::string>& desc,
	bool bInitial, int i) const
{
	int iend = 1;
	if (i == -1) {
		if (bInitial) {
			iend = 6;

		} else {
			iend = 3;
		}
	}
	desc.resize(iend);

	std::ostringstream os;
	os << "DriveDisplacementJoint(" << GetLabel() << ")";

	if (i == -1) {
		std::string name = os.str();
		for (i = 0; i < iend; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": " << eq[i/3] << xyz[i%3];

			desc[i] = os.str();
		}

	} else {
		os << ": " << eq[i/3] << xyz[i%3];
		desc[0] = os.str();
	}
}

/* Dati privati (aggiungere magari le reazioni vincolari) */
unsigned int
DriveDisplacementJoint::iGetNumPrivData(void) const
{
	return 6;
};

unsigned int
DriveDisplacementJoint::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);
	ASSERT(s[0] != '\0');

	if (s[2] != '\0') {
		return 0;
	}

	unsigned int idx = 0;

	switch (s[0]) {
	case 'f':
		idx += 3;
		/* fallthru */
	case 'd':
		break;

	default:
		return 0;
	}

	switch (s[1]) {
	case 'x':
		return idx + 1;
	case 'y':
		return idx + 2;
	case 'z':
		return idx + 3;
	}

	return 0;
}

doublereal
DriveDisplacementJoint::dGetPrivData(unsigned int i) const
{
	ASSERT(i >= 1 && i <= 6);

	switch (i) {
	case 1:
	case 2:
	case 3:
		return Get().dGet(i);

	case 4:
	case 5:
	case 6:
		return F(i - 3);
	}

	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
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

	Mat3x3 MTmp(MatCross, dRef);

	/* node 1 moment */
	WM.Sub(3 + 1, 12 + 1, MTmp);
	/* node 1 constraint */
	WM.Add(12 + 1, 3 + 1, MTmp);

	MTmp = Mat3x3(MatCross, f2Ref);
	/* node 2 moment */
	WM.Add(9 + 1, 12 + 1, MTmp);
	/* node 2 constraint */
	WM.Sub(12 + 1, 9 + 1, MTmp);

	MTmp = Mat3x3(MatCrossCross, F, dRef*dCoef);
	/* node 1 moment */
	WM.Sub(3 + 1, 3 + 1, MTmp);
	
	MTmp = Mat3x3(MatCrossCross, F, f2Ref*dCoef);
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
	
	ASSERT(dCoef != 0.);
	WorkVec.Sub(12 + 1, (pNode2->GetXCurr() + f2Tmp - pNode1->GetXCurr() - d)/dCoef);
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
DriveDisplacementJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering DriveDisplacementJoint::InitialAssJac()" << std::endl);

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
	integer iFirstReactionPrimeIndex = iFirstReactionIndex + 3;   

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

	Vec3 FPrime = Vec3(XCurr, iFirstReactionPrimeIndex + 1);

	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		/* node 1 force */
		WM.DecCoef(iCnt, 24 + iCnt, 1.);
		/* node 2 force */
		WM.IncCoef(12 + iCnt, 24 + iCnt, 1.);

		/* node 1 force derivative */
		WM.DecCoef(6 + iCnt, 27 + iCnt, 1.);
		/* node 2 force derivative */
		WM.IncCoef(18 + iCnt, 27 + iCnt, 1.);

		/* node 1 constraint */
		WM.DecCoef(24 + iCnt, iCnt, 1.);
		/* node 2 constraint */
		WM.IncCoef(24 + iCnt, 12 + iCnt, 1.);

		/* node 1 constraint derivative */
		WM.DecCoef(27 + iCnt, 6 + iCnt, 1.);
		/* node 2 constraint derivative */
		WM.IncCoef(27 + iCnt, 18 + iCnt, 1.);
	}

	Mat3x3 MTmp(MatCross, dRef);

	/* node 1 moment */
	WM.Sub(3 + 1, 24 + 1, MTmp);
	/* node 1 moment derivative */
	WM.Sub(9 + 1, 27 + 1, MTmp);
	/* node 1 constraint */
	WM.Add(24 + 1, 3 + 1, MTmp);
	/* node 1 constraint derivative */
	WM.Add(27 + 1, 9 + 1, MTmp);

	/* in case the drive is differentiable... */
	if (bIsDifferentiable()) {
		WM.Add(27 + 1, 3 + 1, Mat3x3(MatCross, pNode1->GetRCurr()*GetP()));
	}

	MTmp = Mat3x3(MatCross, f2Ref);
	/* node 2 moment */
	WM.Add(15 + 1, 24 + 1, MTmp);
	/* node 2 moment derivatives */
	WM.Add(21 + 1, 27 + 1, MTmp);
	/* node 2 constraint */
	WM.Sub(24 + 1, 15 + 1, MTmp);
	/* node 2 constraint derivative */
	WM.Sub(27 + 1, 21 + 1, MTmp);

	MTmp = Mat3x3(MatCrossCross, F, dRef);
	/* node 1 moment */
	WM.Sub(3 + 1, 3 + 1, MTmp);
	
	MTmp = Mat3x3(MatCrossCross, FPrime, dRef);
	/* node 1 moment derivative */
	WM.Sub(9 + 1, 9 + 1, MTmp);
	
	MTmp = Mat3x3(MatCrossCross, F, f2Ref);
	/* node 2 moment */
	WM.Add(15 + 1, 15 + 1, MTmp);

	MTmp = Mat3x3(MatCrossCross, FPrime, f2Ref);
	/* node 2 moment derivative */
	WM.Add(21 + 1, 21 + 1, MTmp);

	MTmp = Mat3x3(MatCrossCross, pNode1->GetWRef(), dRef);
	/* node 2 constraint derivative */
	WM.Sub(27 + 1, 3 + 1, MTmp);

	MTmp = Mat3x3(MatCrossCross, pNode2->GetWRef(), f2Ref);
	/* node 2 constraint derivative */
	WM.Add(27 + 1, 15 + 1, MTmp);

	return WorkMat;
}

					   
/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
DriveDisplacementJoint::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering DriveDisplacementJoint::InitialAssRes()" << std::endl);   

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

	F = Vec3(XCurr, iFirstReactionIndex + 1);
	Vec3 FPrime = Vec3(XCurr, iReactionPrimeIndex + 1);

	dRef = pNode1->GetRCurr()*(f1 + Get());
	f2Ref = pNode2->GetRCurr()*f2;

	WorkVec.Add(1, F);
	WorkVec.Add(3 + 1, dRef.Cross(F));
	WorkVec.Add(6 + 1, FPrime);
	WorkVec.Add(9 + 1, dRef.Cross(FPrime));

	WorkVec.Sub(12 + 1, F);
	WorkVec.Sub(15 + 1, f2Ref.Cross(F));
	WorkVec.Sub(18 + 1, FPrime);
	WorkVec.Sub(21 + 1, f2Ref.Cross(FPrime));

	WorkVec.Add(24 + 1, pNode1->GetXCurr() + dRef - pNode2->GetXCurr() - f2Ref);

	/* in case the drive is differentiable... */
	Vec3 PhiPrime = pNode1->GetVCurr() + pNode1->GetWCurr().Cross(dRef)
			- pNode2->GetVCurr() + pNode2->GetWCurr().Cross(f2Ref);
	if (bIsDifferentiable()) {
		PhiPrime += pNode1->GetRCurr()*GetP();
	}
	WorkVec.Add(27 + 1, PhiPrime);

	return WorkVec;
}

const MBUnits::Dimensions
DriveDisplacementJoint::GetEquationDimension(integer index) const {
	// DOF == 3
	MBUnits::Dimensions dimension = MBUnits::Dimensions::UnknownDimension;

	switch (index)
	{
	case 1:
		dimension = MBUnits::Dimensions::Length;
		break;
	case 2:
		dimension = MBUnits::Dimensions::Length;
		break;
	case 3:
		dimension = MBUnits::Dimensions::Length;
		break;
	}

	return dimension;
}
					   
/* DriveDisplacementJoint - end */

/* DriveDisplacementPinJoint - begin */

/* Costruttore non banale */
DriveDisplacementPinJoint::DriveDisplacementPinJoint(unsigned int uL,			      
				 const DofOwner* pDO, 
				 const TplDriveCaller<Vec3>* pDC,
				 const StructNode* pN,
				 const Vec3& f,
				 const Vec3& x,
				 flag fOut)
: Elem(uL, fOut), 
Joint(uL, pDO, fOut), 
TplDriveOwner<Vec3>(pDC),
pNode(pN), f(f), x(x),
fRef(Zero3),
dRef(Zero3),
F(Zero3)
{
	ASSERT(pNode != NULL);
	ASSERT(pNode->GetNodeType() == Node::STRUCTURAL);
}

   
/* Distruttore */
DriveDisplacementPinJoint::~DriveDisplacementPinJoint(void)
{
	NO_OP;
}


/* Contributo al file di restart */
std::ostream&
DriveDisplacementPinJoint::Restart(std::ostream& out) const
{
	Joint::Restart(out) << ", drive displacement pin, "
		<< pNode->GetLabel() << ", reference, node, ",
		f.Write(out, ", ") << ", ",
		x.Write(out, ", ") << ", ",
		pGetDriveCaller()->Restart(out) << ';' << std::endl;
	return out;
}

void
DriveDisplacementPinJoint::OutputPrepare(OutputHandler &OH)
{
	if (bToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseBinary(OutputHandler::JOINTS)) {
			std::string name;
			OutputPrepare_int("Drive displacement pin", OH, name);

			Var_d = OH.GetBinaryFile()->CreateVar<Vec3>(name + "d",
				MBUnits::Dimensions::Length,
				"imposed relative displacement (x, y, z)");
		}
#endif // USE_NETCDF
	}
}

void
DriveDisplacementPinJoint::Output(OutputHandler& OH) const
{   
	if (bToBeOutput()) {
		Vec3 d(pNode->GetXCurr() + pNode->GetRCurr()*f - x);

		if (OH.UseText(OutputHandler::JOINTS)) {
			Joint::Output(OH.Joints(), "DriveDisplacementPinJoint", GetLabel(),
					F, Zero3, F, Zero3)
				<< " " << d << std::endl;
		}

#ifdef USE_NETCDF
		if (OH.UseBinary(OutputHandler::JOINTS)) {
			Joint::NetCDFOutput(OH, F, Zero3, F, Zero3);
			OH.WriteVar(Var_d, d);
		}
#endif // USE_NETCDF

	}
}

void
DriveDisplacementPinJoint::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
	if (ph) {
		for (unsigned i = 0; i < ph->size(); i++) {
			Joint::JointHint *pjh = dynamic_cast<Joint::JointHint *>((*ph)[i]);

			if (pjh) {
				if (dynamic_cast<Joint::OffsetHint<1> *>(pjh)) {
					Mat3x3 Rt(pNode->GetRCurr().Transpose());
   
					f = Rt*(x + Get() - pNode->GetXCurr());

				} else if (dynamic_cast<Joint::OffsetHint<0> *>(pjh)) {
					Vec3 fTmp(pNode->GetRCurr()*f);
   
					x = pNode->GetXCurr() + fTmp - Get();

				} else if (dynamic_cast<Joint::ReactionsHint *>(pjh)) {
					/* TODO */
				}
				continue;
			}

			TplDriveHint<Vec3> *pdh = dynamic_cast<TplDriveHint<Vec3> *>((*ph)[i]);

			if (pdh) {
				pedantic_cout("DriveDisplacementPinJoint(" << uLabel << "): "
					"creating drive from hint[" << i << "]..." << std::endl);

				TplDriveCaller<Vec3> *pDC = pdh->pCreateDrive(pDM);
				if (pDC == 0) {
					silent_cerr("DriveDisplacementPinJoint(" << uLabel << "): "
						"unable to create drive "
						"after hint #" << i << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				
				TplDriveOwner<Vec3>::Set(pDC);
				continue;
			}
		}
	}
}

Hint *
DriveDisplacementPinJoint::ParseHint(DataManager *pDM, const char *s) const
{
	if (strncasecmp(s, "offset{" /*}*/ , STRLENOF("offset{" /*}*/ )) == 0)
	{
		s += STRLENOF("offset{" /*}*/ );

		if (strcmp(&s[1], /*{*/ "}") != 0) {
			return 0;
		}

		switch (s[0]) {
		case '1':
			return new Joint::OffsetHint<1>;

		case '0':
			return new Joint::OffsetHint<0>;
		}
	}

	/* take care of "drive" hint... */
	return SimulationEntity::ParseHint(pDM, s);
}

std::ostream&
DriveDisplacementPinJoint::DescribeDof(std::ostream& out,
		const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

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

void
DriveDisplacementPinJoint::DescribeDof(std::vector<std::string>& desc,
	bool bInitial, int i) const
{
	int iend = 1;
	if (i == -1) {
		if (bInitial) {
			iend = 6;

		} else {
			iend = 3;
		}
	}
	desc.resize(iend);

	std::ostringstream os;
	os << "DriveDisplacementPinJoint(" << GetLabel() << ")";

	if (i == -1) {
		std::string name = os.str();
		for (i = 0; i < iend; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": " << dof[i/3] << xyz[i%3];

			desc[i] = os.str();
		}

	} else {
		os << ": " << dof[i/3] << xyz[i%3];
		desc[0] = os.str();
	}
}

std::ostream&
DriveDisplacementPinJoint::DescribeEq(std::ostream& out,
		const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

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
   
void
DriveDisplacementPinJoint::DescribeEq(std::vector<std::string>& desc,
	bool bInitial, int i) const
{
	int iend = 1;
	if (i == -1) {
		if (bInitial) {
			iend = 6;

		} else {
			iend = 3;
		}
	}
	desc.resize(iend);

	std::ostringstream os;
	os << "DriveDisplacementPinJoint(" << GetLabel() << ")";

	if (i == -1) {
		std::string name = os.str();
		for (i = 0; i < iend; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": " << eq[i/3] << xyz[i%3];

			desc[i] = os.str();
		}

	} else {
		os << ": " << eq[i/3] << xyz[i%3];
		desc[0] = os.str();
	}
}

/* Dati privati (aggiungere magari le reazioni vincolari) */
unsigned int
DriveDisplacementPinJoint::iGetNumPrivData(void) const
{
	return 6;
};

unsigned int
DriveDisplacementPinJoint::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);
	ASSERT(s[0] != '\0');

	if (s[2] != '\0') {
		return 0;
	}

	unsigned int idx = 0;

	switch (s[0]) {
	case 'f':
		idx += 3;
		/* fallthru */
	case 'd':
		break;

	default:
		return 0;
	}

	switch (s[1]) {
	case 'x':
		return idx + 1;
	case 'y':
		return idx + 2;
	case 'z':
		return idx + 3;
	}

	return 0;
}

doublereal
DriveDisplacementPinJoint::dGetPrivData(unsigned int i) const
{
	ASSERT(i >= 1 && i <= 6);

	switch (i) {
	case 1:
	case 2:
	case 3:
		return Get().dGet(i);

	case 4:
	case 5:
	case 6:
		return F(i - 3);
	}

	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler& 
DriveDisplacementPinJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef, 
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering DriveDisplacementPinJoint::AssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNodeFirstPosIndex = pNode->iGetFirstPositionIndex();
	integer iNodeFirstMomIndex = pNode->iGetFirstMomentumIndex();
	integer iFirstReactionIndex = iGetFirstIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNodeFirstMomIndex + iCnt);
		WM.PutColIndex(iCnt, iNodeFirstPosIndex + iCnt);	
	}
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(6 + iCnt, iFirstReactionIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iFirstReactionIndex + iCnt);
	}
   
	AssMat(WM, dCoef);

	return WorkMat;
}


void
DriveDisplacementPinJoint::AfterPredict(VectorHandler& /* X */ ,
		VectorHandler& /* XP */ )
{
	/* Recupera i dati */
	fRef = pNode->GetRRef()*f;
	dRef = Get();
}


void
DriveDisplacementPinJoint::AssMat(FullSubMatrixHandler& WM, doublereal dCoef)   
{
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		/* node force */
		WM.IncCoef(iCnt, 6 + iCnt, 1.);

		/* node constraint */
		WM.IncCoef(6 + iCnt, iCnt, 1.);
	}

	Mat3x3 MTmp(MatCross, fRef);

	/* node moment */
	WM.Add(3 + 1, 6 + 1, MTmp);

	/* node constraint */
	WM.Sub(6 + 1, 3 + 1, MTmp);

	MTmp = Mat3x3(MatCrossCross, F, fRef*dCoef);

	/* node moment */
	WM.Add(3 + 1, 3 + 1, MTmp);
}


/* assemblaggio residuo */
SubVectorHandler& 
DriveDisplacementPinJoint::AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering DriveDisplacementPinJoint::AssRes()" << std::endl);   

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera gli indici */
	integer iNodeFirstMomIndex = pNode->iGetFirstMomentumIndex();
	integer iFirstReactionIndex = iGetFirstIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNodeFirstMomIndex + iCnt);
	}

	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(6 + iCnt, iFirstReactionIndex + iCnt);
	}

	F = Vec3(XCurr, iFirstReactionIndex + 1);

	AssVec(WorkVec, dCoef);

	return WorkVec;
}


void
DriveDisplacementPinJoint::AssVec(SubVectorHandler& WorkVec, doublereal dCoef)
{   
	Vec3 fTmp = pNode->GetRCurr()*f;
	Vec3 d = Get();

	WorkVec.Sub(1, F);
	WorkVec.Sub(3 + 1, fTmp.Cross(F));
	
	ASSERT(dCoef != 0.);
	WorkVec.Sub(6 + 1, (pNode->GetXCurr() + fTmp - x - d)/dCoef);
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
DriveDisplacementPinJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering DriveDisplacementPinJoint::InitialAssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNodeFirstPosIndex = pNode->iGetFirstPositionIndex();
	integer iFirstReactionIndex = iGetFirstIndex();
	integer iNodeFirstVelIndex = iNodeFirstPosIndex + 6;
	integer iReactionPrimeIndex = iFirstReactionIndex + 3;   

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNodeFirstPosIndex + iCnt);
		WM.PutColIndex(iCnt, iNodeFirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNodeFirstVelIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNodeFirstVelIndex + iCnt);
		WM.PutRowIndex(12 + iCnt, iFirstReactionIndex + iCnt);
		WM.PutColIndex(12 + iCnt, iFirstReactionIndex + iCnt);
	}

	Vec3 FPrime = Vec3(XCurr, iReactionPrimeIndex + 1);

	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		/* node 2 force */
		WM.IncCoef(iCnt, 12 + iCnt, 1.);

		/* node 2 force derivative */
		WM.IncCoef(6 + iCnt, 15 + iCnt, 1.);

		/* node 2 constraint */
		WM.IncCoef(12 + iCnt, iCnt, 1.);

		/* node 2 constraint derivative */
		WM.IncCoef(15 + iCnt, 6 + iCnt, 1.);
	}

	Mat3x3 MTmp(MatCross, fRef);
	/* node 2 moment */
	WM.Add(3 + 1, 12 + 1, MTmp);
	/* node 2 moment derivatives */
	WM.Add(9 + 1, 15 + 1, MTmp);
	/* node 2 constraint */
	WM.Sub(12 + 1, 3 + 1, MTmp);
	/* node 2 constraint derivative */
	WM.Sub(15 + 1, 9 + 1, MTmp);

	MTmp = Mat3x3(MatCrossCross, F, fRef);
	/* node 2 moment */
	WM.Add(3 + 1, 3 + 1, MTmp);

	MTmp = Mat3x3(MatCrossCross, FPrime, fRef);
	/* node 2 moment derivative */
	WM.Add(9 + 1, 9 + 1, MTmp);

	MTmp = Mat3x3(MatCrossCross, pNode->GetWRef(), fRef);
	/* node 2 constraint derivative */
	WM.Add(15 + 1, 3 + 1, MTmp);

	return WorkMat;
}

					   
/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
DriveDisplacementPinJoint::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering DriveDisplacementPinJoint::InitialAssRes()" << std::endl);   

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera gli indici */
	integer iNodeFirstPosIndex = pNode->iGetFirstPositionIndex();
	integer iFirstReactionIndex = iGetFirstIndex();
	integer iNodeFirstVelIndex = iNodeFirstPosIndex + 6;
	integer iFirstReactionPrimeIndex = iFirstReactionIndex + 3;

	/* Setta gli indici del vettore */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNodeFirstPosIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNodeFirstVelIndex + iCnt);
		WorkVec.PutRowIndex(12 + iCnt, iFirstReactionIndex + iCnt);
	}

	F = Vec3(XCurr, iFirstReactionIndex + 1);
	Vec3 FPrime = Vec3(XCurr, iFirstReactionPrimeIndex + 1);

	dRef = Get();
	fRef = pNode->GetRCurr()*f;

	WorkVec.Sub(1, F);
	WorkVec.Sub(3 + 1, fRef.Cross(F));
	WorkVec.Sub(6 + 1, FPrime);
	WorkVec.Sub(9 + 1, fRef.Cross(FPrime));

	WorkVec.Add(12 + 1, dRef - pNode->GetXCurr() - fRef);

	/* in case the drive is differentiable... */
	Vec3 PhiPrime = pNode->GetVCurr() + pNode->GetWCurr().Cross(fRef);
	if (bIsDifferentiable()) {
		PhiPrime -= GetP();
	}
	WorkVec.Sub(15 + 1, PhiPrime);

	return WorkVec;
}

const MBUnits::Dimensions
DriveDisplacementPinJoint::GetEquationDimension(integer index) const {
	// DOF == 3
	MBUnits::Dimensions dimension = MBUnits::Dimensions::UnknownDimension;

	switch (index)
	{
	case 1:
		dimension = MBUnits::Dimensions::Length;
		break;
	case 2:
		dimension = MBUnits::Dimensions::Length;
		break;
	case 3:
		dimension = MBUnits::Dimensions::Length;
		break;
	}

	return dimension;
}
					   
/* DriveDisplacementPinJoint - end */
