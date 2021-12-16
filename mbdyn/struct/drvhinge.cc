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
#include "drvhinge.h"
#include "Rot.hh"
#include "hint_impl.h"

/* DriveHingeJoint - begin */

/* Costruttore non banale */
DriveHingeJoint::DriveHingeJoint(unsigned int uL,
				 const DofOwner* pDO,
				 const TplDriveCaller<Vec3>* pDC,
				 const StructNode* pN1,
				 const StructNode* pN2,
				 const Mat3x3& R1,
				 const Mat3x3& R2,
				 flag fOut)
: Elem(uL, fOut),
Joint(uL, pDO, fOut),
TplDriveOwner<Vec3>(pDC),
pNode1(pN1), pNode2(pN2), R1h(R1), R2h(R2),
R1Ref(Eye3),
RRef(Eye3),
ThetaRef(Zero3),
ThetaCurr(Zero3),
M(Zero3),
bFirstRes(false)
{
	ASSERT(pNode1 != NULL);
	ASSERT(pNode2 != NULL);
	ASSERT(pNode1->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pNode2->GetNodeType() == Node::STRUCTURAL);
}


/* Distruttore */
DriveHingeJoint::~DriveHingeJoint(void)
{
	NO_OP;
}


/* Contributo al file di restart */
std::ostream&
DriveHingeJoint::Restart(std::ostream& out) const
{
	Joint::Restart(out) << ", drive hinge, "
		<< pNode1->GetLabel() << ", reference, node, 1, ",
		(R1h.GetVec(1)).Write(out, ", ")
		<< ", 2, ", (R1h.GetVec(2)).Write(out, ", ") << ", "
		<< pNode2->GetLabel() << ", reference, node, 1, ",
		(R2h.GetVec(1)).Write(out, ", ")
		<< ", 2, ", (R2h.GetVec(2)).Write(out, ", ") << ", ",
		pGetDriveCaller()->Restart(out) << ';' << std::endl;
	return out;
}

void
DriveHingeJoint::OutputPrepare(OutputHandler &OH)
{
	if (bToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			std::string name;
			OutputPrepare_int("Drive Hinge", OH, name);

			Var_Phi = OH.CreateVar<Vec3>(name + "Theta",
				OutputHandler::Dimensions::rad,
				"Relative orientation");
		}
#endif // USE_NETCDF
	}
}

void
DriveHingeJoint::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {
		Mat3x3 R1(pNode1->GetRCurr()*R1h);
		Vec3 d(MatR2EulerAngles(R1.Transpose()*(pNode2->GetRCurr()*R2h)));
		
		if (OH.UseText(OutputHandler::JOINTS)) {
			Joint::Output(OH.Joints(), "DriveHinge", GetLabel(),
					Zero3, M, Zero3, R1*M)
				<< " " << d*dRaDegr
				<< " " << ThetaCurr << std::endl;
		}

#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			Joint::NetCDFOutput(OH, Zero3, M, Zero3, R1*M);
			OH.WriteNcVar(Var_Phi, ThetaCurr);
		}
#endif // USE_NETCDF

	}
}

void
DriveHingeJoint::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
	if (ph) {
		for (unsigned i = 0; i < ph->size(); i++) {
			Joint::JointHint *pjh = dynamic_cast<Joint::JointHint *>((*ph)[i]);

			if (pjh) {
				if (dynamic_cast<Joint::HingeHint<1> *>(pjh)) {
					(Mat3x3&)R1h = pNode1->GetRCurr().Transpose()*pNode2->GetRCurr()*R2h;

				} else if (dynamic_cast<Joint::HingeHint<2> *>(pjh)) {
					(Mat3x3&)R2h = pNode2->GetRCurr().Transpose()*pNode1->GetRCurr()*R1h;

				} else if (dynamic_cast<Joint::ReactionsHint *>(pjh)) {
					/* TODO */
				}
				continue;
			}

			TplDriveHint<Vec3> *pdh = dynamic_cast<TplDriveHint<Vec3> *>((*ph)[i]);

			if (pdh) {
				pedantic_cout("DriveHingeJoint(" << uLabel << "): "
					"creating drive from hint[" << i << "]..." << std::endl);

				TplDriveCaller<Vec3> *pDC = pdh->pCreateDrive(pDM);
				if (pDC == 0) {
					silent_cerr("DriveHingeJoint(" << uLabel << "): "
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
DriveHingeJoint::ParseHint(DataManager *pDM, const char *s) const
{
	if (strncasecmp(s, "hinge{" /*}*/ , STRLENOF("hinge{" /*}*/ )) == 0)
	{
		s += STRLENOF("hinge{" /*}*/ );

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

	/* take care of "drive" hint... */
	return SimulationEntity::ParseHint(pDM, s);
}

std::ostream&
DriveHingeJoint::DescribeDof(std::ostream& out,
	const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
			"reaction couples [mx,my,mz]" << std::endl;

	if (bInitial) {
		iIndex += 3;
		out
			<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
				"reaction couple derivatives [mPx,mPy,mPz]" << std::endl;
	}

	return out;
}

static const char xyz[] = "xyz";
static const char *dof[] = { "reaction couple m", "reaction couple derivative mP" };
static const char *eq[] = { "orientation constraint g", "orientation constraint derivative w" };

void
DriveHingeJoint::DescribeDof(std::vector<std::string>& desc,
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
	os << "DriveHingeJoint(" << GetLabel() << ")";

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
DriveHingeJoint::DescribeEq(std::ostream& out,
	const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
			"orientation constraints [gx1=gx2,gy1=gy2,gz1=gz2]" << std::endl;

	if (bInitial) {
		iIndex += 3;
		out
			<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
				"angular velocity constraints [wx1=wx2,wy1=wy2,wz1=wz2]" << std::endl;
	}

	return out;

}

void
DriveHingeJoint::DescribeEq(std::vector<std::string>& desc,
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
	os << "DriveHingeJoint(" << GetLabel() << ")";

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
DriveHingeJoint::iGetNumPrivData(void) const
{
	return 6;
};

unsigned int
DriveHingeJoint::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);
	ASSERT(s[0] != '\0');

	unsigned int idx = 0;

	switch (s[0]) {
	case 'M':
		idx += 3;
		/* fallthru */
	case 'r':
		break;

	default:
		return 0;
	}

	if (s[1] == '\0' || s[2] != '\0') {
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
DriveHingeJoint::dGetPrivData(unsigned int i) const
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
		return M(i - 3);

	default:
		ASSERT(0);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
DriveHingeJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering DriveHingeJoint::AssJac()" << std::endl);

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
		WM.PutRowIndex(6 + iCnt, iFirstReactionIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iFirstReactionIndex + iCnt);
	}

	AssMat(WM, dCoef);

	return WorkMat;
}


void
DriveHingeJoint::AfterPredict(VectorHandler& /* X */ ,
		VectorHandler& /* XP */ )
{
	/* Recupera i dati */
	R1Ref = pNode1->GetRRef()*R1h;
	Mat3x3 R1T = R1Ref.Transpose();
	Mat3x3 RD(R1T*pNode2->GetRRef()*R2h);

	/* Calcola la deformazione corrente nel sistema locale (nodo a) */
	ThetaCurr = ThetaRef = RotManip::VecRot(RD);

	/* Calcola l'inversa di Gamma di ThetaRef */
	Mat3x3 GammaRefm1 = RotManip::DRot_I(ThetaRef);

	/* Contributo alla linearizzazione ... */
	RRef = GammaRefm1*R1T;

	/* Flag di aggiornamento dopo la predizione */
	bFirstRes = true;
}


void
DriveHingeJoint::AssMat(FullSubMatrixHandler& WM, doublereal dCoef)
{
	Mat3x3 MCross(MatCross, R1Ref*(M*dCoef));

	WM.Add(1, 1, MCross);
	WM.Sub(4, 1, MCross);

	WM.Sub(1, 7, R1Ref);
	WM.Add(4, 7, R1Ref);

	WM.Sub(7, 1, RRef);
	WM.Add(7, 4, RRef);
}


/* assemblaggio residuo */
SubVectorHandler&
DriveHingeJoint::AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering DriveHingeJoint::AssRes()" << std::endl);

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
		WorkVec.PutRowIndex(6 + iCnt, iFirstReactionIndex + iCnt);
	}

	M = Vec3(XCurr, iFirstReactionIndex+1);

	AssVec(WorkVec, dCoef);

	return WorkVec;
}


void
DriveHingeJoint::AssVec(SubVectorHandler& WorkVec, doublereal dCoef)
{
	Mat3x3 R1(pNode1->GetRCurr()*R1h);

	if (bFirstRes) {
		/* La rotazione e' gia' stata aggiornata da AfterPredict */
		bFirstRes = false;

	} else {
		Mat3x3 R2(pNode2->GetRCurr()*R2h);
		ThetaCurr = RotManip::VecRot(R1.Transpose()*R2);
	}

	Vec3 MTmp(R1*M);

	WorkVec.Add(1, MTmp);
	WorkVec.Sub(3 + 1, MTmp);

	ASSERT(dCoef != 0.);
	WorkVec.Add(6 + 1, (Get() - ThetaCurr)/dCoef);
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
DriveHingeJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering DriveHingeJoint::InitialAssJac()" << std::endl);

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
	integer iNode1FirstVelIndex = iNode1FirstPosIndex + 6;
	integer iNode2FirstVelIndex = iNode2FirstPosIndex + 6;
	integer iReactionPrimeIndex = iFirstReactionIndex + 3;

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(3 + iCnt, iNode1FirstVelIndex + iCnt);
		WM.PutColIndex(3 + iCnt, iNode1FirstVelIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutRowIndex(9 + iCnt, iNode2FirstVelIndex + iCnt);
		WM.PutColIndex(9 + iCnt, iNode2FirstVelIndex + iCnt);
		WM.PutRowIndex(12 + iCnt, iFirstReactionIndex + iCnt);
		WM.PutColIndex(12 + iCnt, iFirstReactionIndex + iCnt);
		WM.PutRowIndex(15 + iCnt, iReactionPrimeIndex + iCnt);
		WM.PutColIndex(15 + iCnt, iReactionPrimeIndex + iCnt);
	}

	Mat3x3 Ra(pNode1->GetRRef()*R1h);
	Mat3x3 RaT(Ra.Transpose());
	Vec3 Wa(pNode1->GetWRef());
	Vec3 Wb(pNode2->GetWRef());

	Mat3x3 MTmp(MatCross, M);
	Mat3x3 MPrimeTmp(MatCross, Ra*Vec3(XCurr, iReactionPrimeIndex + 1));

	WM.Add(1, 1, MTmp);
	WM.Add(3 + 1, 3 + 1, MTmp);
	WM.Sub(6 + 1, 1, MTmp);
	WM.Sub(9 + 1, 3 + 1, MTmp);

	MTmp = Wa.Cross(MTmp) + MPrimeTmp;
	WM.Add(3 + 1, 1, MTmp);
	WM.Sub(9 + 1, 1, MTmp);

	WM.Add(6 + 1, 12 + 1, Ra);
	WM.Add(9 + 1, 15 + 1, Ra);
	WM.Sub(1, 12 + 1, Ra);
	WM.Sub(3 + 1, 15 + 1, Ra);

	MTmp = Wa.Cross(Ra);
	WM.Add(9 + 1, 12 + 1, MTmp);
	WM.Sub(3 + 1, 12 + 1, MTmp);

	WM.Add(12 + 1, 6 + 1, RaT);
	WM.Sub(12 + 1, 1, RaT);
	WM.Sub(15 + 1, 3 + 1, RaT);
	WM.Add(15 + 1, 9 + 1, RaT);
	WM.Add(15 + 1, 1, RaT*Mat3x3(MatCross, Wb - Wa));

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
DriveHingeJoint::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering DriveHingeJoint::InitialAssRes()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex() + 3;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex() + 3;
	integer iFirstReactionIndex = iGetFirstIndex();
	integer iNode1FirstVelIndex = iNode1FirstPosIndex + 6;
	integer iNode2FirstVelIndex = iNode2FirstPosIndex + 6;
	integer iReactionPrimeIndex = iFirstReactionIndex + 3;

	/* Setta gli indici del vettore */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WorkVec.PutRowIndex(3 + iCnt, iNode1FirstVelIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
		WorkVec.PutRowIndex(9 + iCnt, iNode2FirstVelIndex + iCnt);
		WorkVec.PutRowIndex(12 + iCnt, iFirstReactionIndex + iCnt);
		WorkVec.PutRowIndex(15 + iCnt, iReactionPrimeIndex + iCnt);
	}

	Mat3x3 R1(pNode1->GetRCurr()*R1h);
	Mat3x3 R1T(R1.Transpose());
	Vec3 Wa(pNode1->GetWCurr());
	Vec3 Wb(pNode2->GetWCurr());

	M = Vec3(XCurr, iFirstReactionIndex+1);
	Vec3 MPrime = Vec3(XCurr, iReactionPrimeIndex+1);

	Vec3 MTmp(R1*M);
	Vec3 MPrimeTmp(Wa.Cross(MTmp) + R1*MPrime);

	Mat3x3 R2(pNode2->GetRCurr()*R2h);
	ThetaCurr = RotManip::VecRot(R1T*R2);

	Vec3 ThetaPrime = R1T*(Wb-Wa);

	WorkVec.Add(1, MTmp);
	WorkVec.Add(3 + 1, MPrimeTmp);
	WorkVec.Sub(6 + 1, MTmp);
	WorkVec.Sub(9 + 1, MPrimeTmp);
	WorkVec.Add(12 + 1, Get() - ThetaCurr);
	if (bIsDifferentiable()) {
		ThetaPrime -= GetP();
	}
	WorkVec.Sub(15 + 1, ThetaPrime);

	return WorkVec;
}

const OutputHandler::Dimensions
DriveHingeJoint::GetEquationDimension(integer index) const {
	// DOF == 6
	OutputHandler::Dimensions dimension = OutputHandler::Dimensions::UnknownDimension;

	switch (index)
	{
	case 1:
		dimension = OutputHandler::Dimensions::rad;
		break;
	case 2:
		dimension = OutputHandler::Dimensions::rad;
		break;
	case 3:
		dimension = OutputHandler::Dimensions::rad;
		break;
	case 4:
		dimension = OutputHandler::Dimensions::AngularVelocity;
		break;
	case 5:
		dimension = OutputHandler::Dimensions::AngularVelocity;
		break;
	case 6:
		dimension = OutputHandler::Dimensions::AngularVelocity;
		break;
	}

	return dimension;
}
/* DriveHingeJoint - end */
