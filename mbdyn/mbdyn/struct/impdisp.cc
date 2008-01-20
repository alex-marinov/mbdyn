/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2008
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
#include "impdisp.h"
#include "hint_impl.h"

/* ImposedDisplacementJoint - begin */

/* Costruttore non banale */
ImposedDisplacementJoint::ImposedDisplacementJoint(unsigned int uL,
	const DofOwner* pDO,
	const DriveCaller* pDC,
	const StructNode* pN1,
	const StructNode* pN2,
	const Vec3& f1,
	const Vec3& f2,
	const Vec3& e1,
	flag fOut)
: Elem(uL, fOut),
Joint(uL, pDO, fOut),
DriveOwner(pDC),
pNode1(pN1), pNode2(pN2), f1(f1), f2(f2), e1(e1),
e1xf1(e1*f1),
f2Ref(0.),
dRef(0.),
e1Ref(0.),
F(0.)
{
	ASSERT(pNode1 != NULL);
	ASSERT(pNode2 != NULL);
	ASSERT(pNode1->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pNode2->GetNodeType() == Node::STRUCTURAL);
}


/* Distruttore */
ImposedDisplacementJoint::~ImposedDisplacementJoint(void)
{
	NO_OP;
}


/* Contributo al file di restart */
std::ostream&
ImposedDisplacementJoint::Restart(std::ostream& out) const
{
	Joint::Restart(out) << ", imposed displacement, "
		<< pNode1->GetLabel() << ", "
		"reference, node, ", f1.Write(out, ", ") << ", "
		<< pNode2->GetLabel() << ", "
		"reference, node, ", f2.Write(out, ", ") << ", "
		"reference, node, ", e1.Write(out, ", ") << ", ",
		pGetDriveCaller()->Restart(out) << ';' << std::endl;
	return out;
}


void
ImposedDisplacementJoint::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		Vec3 d(pNode2->GetXCurr() + pNode2->GetRCurr()*f2
			- pNode1->GetXCurr() - pNode1->GetRCurr()*f1);
		Vec3 FTmp(e1*F);
		Joint::Output(OH.Joints(), "ImposedDisplacementJoint", GetLabel(),
				FTmp, Zero3, pNode1->GetRCurr().Transpose()*FTmp, Zero3)
			<< " " << d << std::endl;
	}
}

void
ImposedDisplacementJoint::SetValue(DataManager *pDM,
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

					f1 = R1t*(pNode2->GetXCurr() + fTmp2 - pNode1->GetXCurr()) - e1*dGet();

				} else if (dynamic_cast<Joint::OffsetHint<2> *>(pjh)) {
					Mat3x3 R2t(pNode2->GetRCurr().Transpose());
					Vec3 fTmp1(pNode1->GetRCurr()*f1);

					f2 = R2t*(pNode1->GetXCurr() + fTmp1 - pNode2->GetXCurr()) + e1*dGet();

				} else if (dynamic_cast<Joint::ReactionsHint *>(pjh)) {
					/* TODO */
				}
				continue;
			}

			DriveHint *pdh = dynamic_cast<DriveHint *>((*ph)[i]);

			if (pdh) {
				pedantic_cout("ImposedDisplacementJoint(" << uLabel << "): "
					"creating drive from hint[" << i << "]" << std::endl);

				DriveCaller *pDC = pdh->pCreateDrive(pDM);
				if (pDC == 0) {
					silent_cerr("ImposedDisplacementJoint(" << uLabel << "): "
						"unable to create drive "
						"after hint[" << i << "]" << std::endl);
					throw ErrGeneric();
				}

				DriveOwner::Set(pDC);
				continue;
			}
		}
	}
}

Hint *
ImposedDisplacementJoint::ParseHint(DataManager *pDM, const char *s) const
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
ImposedDisplacementJoint::DescribeDof(std::ostream& out,
	const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << ": "
			"reaction force" << std::endl;

	if (bInitial) {
		iIndex += 1;
		out
			<< prefix << iIndex + 1 << ": "
				"reaction force derivative" << std::endl;
	}

	return out;
}

static const char *dof[] = {
	"reaction force f",
	"reaction force derivative fP"
};
static const char *eq[] = {
	"position constraint P",
	"position constraint derivative v"
};

void
ImposedDisplacementJoint::DescribeDof(std::vector<std::string>& desc,
	bool bInitial, int i) const
{
	int nself = 1;
	
	if (i == -1) {
		if (bInitial) {
			nself *= 2;
		}
	}
	desc.resize(nself);

	std::ostringstream os;
	os << "ImposedDisplacementJoint(" << GetLabel() << ")";

	if (i == -1) {
		std::string name(os.str());

		for (int i = 0; i < nself; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": " << dof[i];
			desc[i] = os.str();
		}

	} else {
		if (i < 0 || i >= nself) {
			// error
			throw ErrGeneric();
		}

		os << ": " << dof[i];
		desc[0] = os.str();
	}
}

std::ostream&
ImposedDisplacementJoint::DescribeEq(std::ostream& out,
	const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << ": "
			"position constraint" << std::endl;

	if (bInitial) {
		iIndex += 1;
		out
			<< prefix << iIndex + 1 << ": "
				"velocity constraint" << std::endl;
	}

	return out;

}

void
ImposedDisplacementJoint::DescribeEq(std::vector<std::string>& desc,
	bool bInitial, int i) const
{
	int nself = 1;
	
	if (i == -1) {
		if (bInitial) {
			nself *= 2;
		}
	}
	desc.resize(nself);

	std::ostringstream os;
	os << "ImposedDisplacementJoint(" << GetLabel() << ")";

	if (i == -1) {
		std::string name(os.str());

		for (int i = 0; i < nself; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": " << eq[i];
			desc[i] = os.str();
		}

	} else {
		if (i < 0 || i >= nself) {
			// error
			throw ErrGeneric();
		}

		os << ": " << eq[i];
		desc[0] = os.str();
	}
}

/* Dati privati (aggiungere magari le reazioni vincolari) */
unsigned int
ImposedDisplacementJoint::iGetNumPrivData(void) const
{
	return 2;
};

unsigned int
ImposedDisplacementJoint::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);
	ASSERT(s[0] != '\0');

	if (s[1] != '\0') {
		return 0;
	}

	switch (s[0]) {
	case 'd':
		return 1;

	case 'f':
		return 2;
	}

	return 0;
}

doublereal
ImposedDisplacementJoint::dGetPrivData(unsigned int i) const
{
	ASSERT(i >= 1 && i <= 2);

	switch (i) {
	case 1:
		return dGet();

	case 2:
		return F;
	}

	throw ErrGeneric();
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
ImposedDisplacementJoint::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering ImposedDisplacementJoint::AssJac()" << std::endl);

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
	WM.PutRowIndex(12 + 1, iFirstReactionIndex + 1);
	WM.PutColIndex(12 + 1, iFirstReactionIndex + 1);

	AssMat(WM, dCoef);

	return WorkMat;
}


void
ImposedDisplacementJoint::AfterPredict(VectorHandler& /* X */ ,
		VectorHandler& /* XP */ )
{
	/* Recupera i dati */
	f2Ref = pNode2->GetRRef()*f2;
	dRef = pNode2->GetXCurr() + f2Ref - pNode1->GetXCurr();
	e1Ref = pNode1->GetRRef()*e1;
}


void
ImposedDisplacementJoint::AssMat(FullSubMatrixHandler& WM, doublereal dCoef)
{
	Vec3 dxe1(dRef.Cross(e1Ref));
	Vec3 f2xe1(f2Ref.Cross(e1Ref));

	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		doublereal d = e1Ref(iCnt);

		/* node 1 force */
		WM.DecCoef(iCnt, 12 + 1, d);
		/* node 2 force */
		WM.IncCoef(6 + iCnt, 12 + 1, d);

		/* node 1 constraint */
		WM.DecCoef(12 + 1, iCnt, d);
		/* node 2 constraint */
		WM.IncCoef(12 + 1, 6 + iCnt, d);

		d = dxe1(iCnt);

		/* node 1 moment */
		WM.DecCoef(3 + iCnt, 12 + 1, d);

		/* node 1 constraint */
		WM.DecCoef(12 + 1, 3 + iCnt, d);

		d = f2xe1(iCnt);

		/* node 2 moment */
		WM.IncCoef(9 + iCnt, 12 + 1, d);

		/* node 2 constraint */
		WM.IncCoef(12 + 1, 9 + iCnt, d);
	}

	Vec3 FTmp(e1Ref*(F*dCoef));
	Mat3x3 MTmp(FTmp);

	/* node 1 force */
	WM.Add(1, 3 + 1, MTmp);

	/* node 2 force */
	WM.Sub(6 + 1, 3 + 1, MTmp);

	/* node 1 moment */
	WM.Sub(3 + 1, 1, MTmp);
	WM.Add(3 + 1, 6 + 1, MTmp);

	MTmp = Mat3x3(dRef, FTmp);

	/* node 1 moment */
	WM.Add(3 + 1, 3 + 1, MTmp);

	MTmp = Mat3x3(f2Ref, FTmp);

	/* node 2 moment */
	WM.Sub(9 + 1, 3 + 1, MTmp);

	MTmp = Mat3x3(FTmp, f2Ref);

	/* node 2 moment */
	WM.Sub(3 + 1, 9 + 1, MTmp);
	WM.Add(9 + 1, 9 + 1, MTmp);
}


/* assemblaggio residuo */
SubVectorHandler&
ImposedDisplacementJoint::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering ImposedDisplacementJoint::AssRes()" << std::endl);

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

	WorkVec.PutRowIndex(12 + 1, iFirstReactionIndex + 1);

	F = XCurr(iFirstReactionIndex + 1);

	AssVec(WorkVec, dCoef);

	return WorkVec;
}


void
ImposedDisplacementJoint::AssVec(SubVectorHandler& WorkVec, doublereal dCoef)
{
	Mat3x3 R1(pNode1->GetRCurr());

	Vec3 f2Tmp(pNode2->GetRCurr()*f2);
	Vec3 dTmp(pNode2->GetXCurr() + f2Tmp - pNode1->GetXCurr());
	Vec3 e1Tmp(pNode1->GetRCurr()*e1);

	Vec3 FTmp(e1Tmp*F);

	WorkVec.Add(1, FTmp);
	WorkVec.Add(3 + 1, dTmp.Cross(FTmp));
	WorkVec.Sub(6 + 1, FTmp);
	WorkVec.Sub(9 + 1, f2Tmp.Cross(FTmp));

	if (dCoef != 0.) {
		WorkVec.DecCoef(12 + 1, (e1Tmp*dTmp - (e1xf1 + dGet()))/dCoef);
	}
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
ImposedDisplacementJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering ImposedDisplacementJoint::InitialAssJac()" << std::endl);

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
	integer iFirstReactionPrimeIndex = iFirstReactionIndex + 1;

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
	}

	WM.PutRowIndex(24 + 1, iFirstReactionIndex + 1);
	WM.PutColIndex(24 + 1, iFirstReactionIndex + 1);
	WM.PutRowIndex(24 + 2, iFirstReactionPrimeIndex + 1);
	WM.PutColIndex(24 + 2, iFirstReactionPrimeIndex + 1);

	doublereal FPrime = XCurr(iFirstReactionPrimeIndex + 1);

	Vec3 dxe1(dRef.Cross(e1Ref));
	Vec3 f2xe1(f2Ref.Cross(e1Ref));
	Vec3 w1xe1(pNode1->GetWCurr().Cross(e1Ref));
	Vec3 dPrime(pNode2->GetVCurr()
		+ pNode2->GetWCurr().Cross(f2Ref)
		- pNode1->GetVCurr());
	Vec3 t1(e1Ref.Cross(dPrime - pNode1->GetWCurr().Cross(dRef)));
	Vec3 t2(f2Ref.Cross(e1Ref.Cross(pNode2->GetWCurr() - pNode1->GetWCurr())));

	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		doublereal d = e1Ref(iCnt);

		/* node 1 force */
		WM.DecCoef(iCnt, 24 + 1, d);
		/* node 2 force */
		WM.IncCoef(12 + iCnt, 24 + 1, d);

		/* node 1 force derivative */
		WM.DecCoef(6 + iCnt, 24 + 2, d);
		/* node 2 force derivative */
		WM.IncCoef(18 + iCnt, 24 + 2, d);

		/* node 1 constraint */
		WM.DecCoef(24 + 1, iCnt, d);
		/* node 2 constraint */
		WM.IncCoef(24 + 1, 12 + iCnt, d);

		/* node 1 constraint derivative */
		WM.DecCoef(24 + 2, 6 + iCnt, d);
		/* node 2 constraint derivative */
		WM.IncCoef(24 + 2, 18 + iCnt, d);

		d = dxe1(iCnt);

		/* node 1 moment */
		WM.DecCoef(3 + iCnt, 24 + 1, d);

		/* node 1 moment derivative */
		WM.DecCoef(9 + iCnt, 24 + 2, d);

		/* node 1 constraint */
		WM.DecCoef(24 + 1, 3 + iCnt, d);

		/* node 1 constraint derivative */
		WM.DecCoef(24 + 2, 9 + iCnt, d);

		d = f2xe1(iCnt);

		/* node 2 moment */
		WM.DecCoef(15 + iCnt, 24 + 1, d);

		/* node 2 moment derivative */
		WM.DecCoef(21 + iCnt, 24 + 2, d);

		/* node 2 constraint */
		WM.DecCoef(24 + 1, 15 + iCnt, d);

		/* node 2 constraint derivative */
		WM.DecCoef(24 + 2, 21 + iCnt, d);

		d = w1xe1(iCnt);

		/* node 1 constraint derivative */
		WM.DecCoef(24 + 2, iCnt, d);

		/* node 2 constraint derivative */
		WM.IncCoef(24 + 2, 12 + iCnt, d);

		d = t1(iCnt);

		/* node 1 constraint derivative */
		WM.IncCoef(24 + 2, 3 + iCnt, d);

		d = t2(iCnt);

		/* node 2 constraint derivative */
		WM.IncCoef(24 + 2, 9 + iCnt, d);
	}

	Vec3 FTmp(e1Ref*F);
	Vec3 FPrimeTmp(e1Ref*FPrime);

	Mat3x3 MTmp(FTmp);

	/* node 1 force */
	WM.Add(1, 3 + 1, MTmp);

	/* node 2 force */
	WM.Sub(12 + 1, 3 + 1, MTmp);

	/* node 1 moment */
	WM.Sub(3 + 1, 1, MTmp);
	WM.Add(3 + 1, 12 + 1, MTmp);

	MTmp = Mat3x3(FTmp, f2Ref);

	/* node 1 moment */
	WM.Sub(3 + 1, 15 + 1, MTmp);

	/* node 2 moment */
	WM.Add(15 + 1, 15 + 1, MTmp);

	MTmp = Mat3x3(dRef, FTmp);

	/* node 1 moment */
	WM.Add(3 + 1, 3 + 1, MTmp);

	MTmp = Mat3x3(f2Ref, FTmp);

	/* node 1 moment */
	WM.Sub(15 + 1, 3 + 1, MTmp);



	MTmp = Mat3x3(FPrimeTmp);

	/* node 1 force derivative */
	WM.Add(6 + 1, 3 + 1, MTmp);

	/* node 2 force derivative */
	WM.Sub(18 + 1, 3 + 1, MTmp);

	/* node 1 moment derivative */
	WM.Sub(9 + 1, 1, MTmp);
	WM.Add(9 + 1, 12 + 1, MTmp);

	MTmp = Mat3x3(FPrimeTmp, f2Ref);

	/* node 1 moment derivative */
	WM.Sub(9 + 1, 15 + 1, MTmp);

	/* node 2 moment derivative */
	WM.Add(21 + 1, 15 + 1, MTmp);

	MTmp = Mat3x3(dRef, FPrimeTmp);

	/* node 1 moment derivative */
	WM.Add(9 + 1, 3 + 1, MTmp);

	MTmp = Mat3x3(f2Ref, FPrimeTmp);

	/* node 1 moment derivative */
	WM.Sub(21 + 1, 3 + 1, MTmp);

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
ImposedDisplacementJoint::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering ImposedDisplacementJoint::InitialAssRes()" << std::endl);

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
	integer iFirstReactionPrimeIndex = iFirstReactionIndex + 1;

	/* Setta gli indici del vettore */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode1FirstVelIndex + iCnt);
		WorkVec.PutRowIndex(12 + iCnt, iNode2FirstPosIndex + iCnt);
		WorkVec.PutRowIndex(18 + iCnt, iNode2FirstVelIndex + iCnt);
	}
	WorkVec.PutRowIndex(24 + 1, iFirstReactionIndex + 1);
	WorkVec.PutRowIndex(24 + 2, iFirstReactionPrimeIndex + 1);

	F = XCurr(iFirstReactionIndex + 1);
	doublereal FPrime = XCurr(iFirstReactionPrimeIndex + 1);

	f2Ref = pNode2->GetRCurr()*f2;
	dRef = pNode2->GetXCurr() + f2Ref - pNode1->GetXCurr();
	e1Ref = pNode1->GetRCurr()*e1;

	Vec3 FTmp(e1Ref*F);
	Vec3 FPrimeTmp(e1Ref*FPrime);

	WorkVec.Add(1, FTmp);
	WorkVec.Add(3 + 1, dRef.Cross(FTmp));
	WorkVec.Add(6 + 1, FPrimeTmp);
	WorkVec.Add(9 + 1, dRef.Cross(FPrimeTmp));

	WorkVec.Sub(12 + 1, FTmp);
	WorkVec.Sub(15 + 1, f2Ref.Cross(FTmp));
	WorkVec.Sub(18 + 1, FPrimeTmp);
	WorkVec.Sub(21 + 1, f2Ref.Cross(FPrimeTmp));

	WorkVec.Sub(24 + 1, e1Ref*dRef - e1xf1 - dGet());

	Vec3 dPrimeTmp(pNode2->GetVCurr()
		+ pNode2->GetWCurr().Cross(f2Ref)
		- pNode1->GetVCurr());
	doublereal d = e1Ref*dPrimeTmp + dRef*(pNode1->GetWCurr().Cross(e1Ref));
	if (bIsDifferentiable()) {
		d -= dGetP();
	}
	WorkVec.Sub(24 + 2, d);

	return WorkVec;
}

/* ImposedDisplacementJoint - end */

/* ImposedDisplacementPinJoint - begin */

/* Costruttore non banale */
ImposedDisplacementPinJoint::ImposedDisplacementPinJoint(unsigned int uL,
	const DofOwner* pDO,
	const DriveCaller* pDC,
	const StructNode* pN,
	const Vec3& f,
	const Vec3& x,
	const Vec3& e,
	flag fOut)
: Elem(uL, fOut),
Joint(uL, pDO, fOut),
DriveOwner(pDC),
pNode(pN), f(f), x(x), e(e),
fRef(0.),
dRef(0.),
F(0.)
{
	ASSERT(pNode != NULL);
	ASSERT(pNode->GetNodeType() == Node::STRUCTURAL);
}


/* Distruttore */
ImposedDisplacementPinJoint::~ImposedDisplacementPinJoint(void)
{
	NO_OP;
}


/* Contributo al file di restart */
std::ostream&
ImposedDisplacementPinJoint::Restart(std::ostream& out) const
{
	Joint::Restart(out) << ", imposed displacement pin, "
		<< pNode->GetLabel() << ", "
		"reference, node, ", f.Write(out, ", ") << ", ",
		x.Write(out, ", ") << ", ",
		e.Write(out, ", ") << ", ",
		pGetDriveCaller()->Restart(out) << ';' << std::endl;
	return out;
}


void
ImposedDisplacementPinJoint::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		Vec3 d(pNode->GetXCurr() + pNode->GetRCurr()*f - x);
		Vec3 FTmp(e*F);
		Joint::Output(OH.Joints(), "ImposedDisplacementPinJoint", GetLabel(),
				FTmp, Zero3, FTmp, Zero3)
			<< " " << d << std::endl;
	}
}

void
ImposedDisplacementPinJoint::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
	if (ph) {
		for (unsigned i = 0; i < ph->size(); i++) {
			Joint::JointHint *pjh = dynamic_cast<Joint::JointHint *>((*ph)[i]);

			if (pjh) {
				if (dynamic_cast<Joint::OffsetHint<1> *>(pjh)) {
					Mat3x3 Rt(pNode->GetRCurr().Transpose());

					f = Rt*(x + e*dGet() - pNode->GetXCurr());

				} else if (dynamic_cast<Joint::OffsetHint<0> *>(pjh)) {
					Vec3 fTmp(pNode->GetRCurr()*f);

					x = pNode->GetXCurr() + fTmp - e*dGet();

				} else if (dynamic_cast<Joint::ReactionsHint *>(pjh)) {
					/* TODO */
				}
				continue;
			}

			DriveHint *pdh = dynamic_cast<DriveHint *>((*ph)[i]);

			if (pdh) {
				pedantic_cout("ImposedDisplacementPinJoint(" << uLabel << "): "
					"creating drive from hint[" << i << "]" << std::endl);

				DriveCaller *pDC = pdh->pCreateDrive(pDM);
				if (pDC == 0) {
					silent_cerr("ImposedDisplacementPinJoint(" << uLabel << "): "
						"unable to create drive "
						"after hint[" << i << "]" << std::endl);
					throw ErrGeneric();
				}

				DriveOwner::Set(pDC);
				continue;
			}
		}
	}
}

Hint *
ImposedDisplacementPinJoint::ParseHint(DataManager *pDM, const char *s) const
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
ImposedDisplacementPinJoint::DescribeDof(std::ostream& out,
	const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << ": "
			"reaction force" << std::endl;

	if (bInitial) {
		iIndex += 1;
		out
			<< prefix << iIndex + 1 << ": "
				"reaction force derivative" << std::endl;
	}

	return out;
}

void
ImposedDisplacementPinJoint::DescribeDof(std::vector<std::string>& desc,
	bool bInitial, int i) const
{
	int nself = 1;
	
	if (i == -1) {
		if (bInitial) {
			nself *= 2;
		}
	}
	desc.resize(nself);

	std::ostringstream os;
	os << "ImposedDisplacementPinJoint(" << GetLabel() << ")";

	if (i == -1) {
		std::string name(os.str());

		for (int i = 0; i < nself; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": " << dof[i];
			desc[i] = os.str();
		}

	} else {
		if (i < 0 || i >= nself) {
			// error
			throw ErrGeneric();
		}

		os << ": " << dof[i];
		desc[0] = os.str();
	}
}

std::ostream&
ImposedDisplacementPinJoint::DescribeEq(std::ostream& out,
		const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << ": "
			"position constraint" << std::endl;

	if (bInitial) {
		iIndex += 1;
		out
			<< prefix << iIndex + 1 << ": "
				"velocity constraint" << std::endl;
	}

	return out;

}

void
ImposedDisplacementPinJoint::DescribeEq(std::vector<std::string>& desc,
	bool bInitial, int i) const
{
	int nself = 1;
	
	if (i == -1) {
		if (bInitial) {
			nself *= 2;
		}
	}
	desc.resize(nself);

	std::ostringstream os;
	os << "ImposedDisplacementPinJoint(" << GetLabel() << ")";

	if (i == -1) {
		std::string name(os.str());

		for (int i = 0; i < nself; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": " << eq[i];
			desc[i] = os.str();
		}

	} else {
		if (i < 0 || i >= nself) {
			// error
			throw ErrGeneric();
		}

		os << ": " << eq[i];
		desc[0] = os.str();
	}
}

/* Dati privati (aggiungere magari le reazioni vincolari) */
unsigned int
ImposedDisplacementPinJoint::iGetNumPrivData(void) const
{
	return 2;
};

unsigned int
ImposedDisplacementPinJoint::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);
	ASSERT(s[0] != '\0');

	if (s[1] != '\0') {
		return 0;
	}

	switch (s[0]) {
	case 'd':
		return 1;

	case 'f':
		return 2;
	}

	return 0;
}

doublereal
ImposedDisplacementPinJoint::dGetPrivData(unsigned int i) const
{
	ASSERT(i >= 1 && i <= 2);

	switch (i) {
	case 1:
		return dGet();

	case 2:
		return F;
	}

	throw ErrGeneric();
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
ImposedDisplacementPinJoint::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering ImposedDisplacementPinJoint::AssJac()" << std::endl);

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
	WM.PutRowIndex(6 + 1, iFirstReactionIndex + 1);
	WM.PutColIndex(6 + 1, iFirstReactionIndex + 1);

	AssMat(WM, dCoef);

	return WorkMat;
}


void
ImposedDisplacementPinJoint::AfterPredict(VectorHandler& /* X */ ,
	VectorHandler& /* XP */ )
{
	/* Recupera i dati */
	fRef = pNode->GetRRef()*f;
	dRef = pNode->GetXCurr() + fRef - x;
}


void
ImposedDisplacementPinJoint::AssMat(FullSubMatrixHandler& WM, doublereal dCoef)
{
	Vec3 fxe(fRef.Cross(e));

	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		doublereal d = e(iCnt);

		/* node force */
		WM.IncCoef(iCnt, 6 + 1, d);

		/* node constraint */
		WM.IncCoef(6 + 1, iCnt, d);

		d = fxe(iCnt);

		/* node moment */
		WM.IncCoef(3 + iCnt, 6 + 1, d);

		/* node constraint */
		WM.IncCoef(6 + 1, 3 + iCnt, d);
	}

	Mat3x3 MTmp(e*(F*dCoef), fRef);

	/* node moment */
	WM.Add(3 + 1, 3 + 1, MTmp);
}


/* assemblaggio residuo */
SubVectorHandler&
ImposedDisplacementPinJoint::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering ImposedDisplacementPinJoint::AssRes()" << std::endl);

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

	WorkVec.PutRowIndex(6 + 1, iFirstReactionIndex + 1);

	F = XCurr(iFirstReactionIndex + 1);

	AssVec(WorkVec, dCoef);

	return WorkVec;
}


void
ImposedDisplacementPinJoint::AssVec(SubVectorHandler& WorkVec, doublereal dCoef)
{
	Vec3 fTmp = pNode->GetRCurr()*f;
	Vec3 d = pNode->GetXCurr() + fTmp - x;

	Vec3 FTmp(e*F);

	WorkVec.Sub(1, FTmp);
	WorkVec.Sub(3 + 1, fTmp.Cross(FTmp));

	if (dCoef != 0.) {
		WorkVec.Sub(6 + 1, (e*d - dGet())/dCoef);
	}
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
ImposedDisplacementPinJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering ImposedDisplacementPinJoint::InitialAssJac()" << std::endl);

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
	integer iFirstReactionPrimeIndex = iFirstReactionIndex + 1;

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNodeFirstPosIndex + iCnt);
		WM.PutColIndex(iCnt, iNodeFirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNodeFirstVelIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNodeFirstVelIndex + iCnt);
	}
	WM.PutRowIndex(12 + 1, iFirstReactionIndex + 1);
	WM.PutColIndex(12 + 1, iFirstReactionIndex + 1);
	WM.PutRowIndex(12 + 2, iFirstReactionPrimeIndex + 1);
	WM.PutColIndex(12 + 2, iFirstReactionPrimeIndex + 1);

	doublereal FPrime = XCurr(iFirstReactionPrimeIndex + 1);

	Vec3 fxe(fRef.Cross(e));
	Vec3 fxomegaxe(fRef.Cross(pNode->GetWRef().Cross(e)));

	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		doublereal d = e(iCnt);

		/* node force */
		WM.IncCoef(iCnt, 12 + 1, d);

		/* node force derivative */
		WM.IncCoef(6 + iCnt, 12 + 2, d);

		/* node constraint */
		WM.IncCoef(12 + 1, iCnt, d);

		/* node constraint derivative */
		WM.IncCoef(12 + 2, 6 + iCnt, d);

		d = fxe(iCnt);

		/* node couple */
		WM.IncCoef(3 + iCnt, 12 + 1, d);

		/* node couple derivative */
		WM.IncCoef(9 + iCnt, 12 + 2, d);

		/* node constraint */
		WM.IncCoef(12 + 1, 3 + iCnt, d);

		/* node constraint derivative */
		WM.IncCoef(12 + 2, 9 + iCnt, d);

		d = fxomegaxe(iCnt);

		/* node constraint derivative */
		WM.DecCoef(12 + 2, 3 + iCnt, d);
	}

	Mat3x3 MTmp(e, fRef*F);

	/* node moment */
	WM.Add(3 + 1, 3 + 1, MTmp);

	MTmp = Mat3x3(e, fRef*FPrime);

	/* node moment derivative */
	WM.Add(9 + 1, 3 + 1, MTmp);

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
ImposedDisplacementPinJoint::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering ImposedDisplacementPinJoint::InitialAssRes()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera gli indici */
	integer iNodeFirstPosIndex = pNode->iGetFirstPositionIndex();
	integer iFirstReactionIndex = iGetFirstIndex();
	integer iNodeFirstVelIndex = iNodeFirstPosIndex + 6;
	integer iFirstReactionPrimeIndex = iFirstReactionIndex + 1;

	/* Setta gli indici del vettore */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNodeFirstPosIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNodeFirstVelIndex + iCnt);
	}
	WorkVec.PutRowIndex(12 + 1, iFirstReactionIndex + 1);
	WorkVec.PutRowIndex(12 + 2, iFirstReactionPrimeIndex + 1);

	F = XCurr(iFirstReactionIndex + 1);
	doublereal FPrime = XCurr(iFirstReactionPrimeIndex + 1);

	fRef = pNode->GetRCurr()*f;
	dRef = pNode->GetXCurr() + fRef - x;

	Vec3 FTmp = e*F;
	Vec3 FPrimeTmp = e*FPrime;

	WorkVec.Sub(1, FTmp);
	WorkVec.Sub(3 + 1, fRef.Cross(FTmp));
	WorkVec.Sub(6 + 1, FPrimeTmp);
	WorkVec.Sub(9 + 1, fRef.Cross(FPrimeTmp));

	WorkVec.Sub(12 + 1, e*dRef - dGet());

	doublereal d = e*(pNode->GetVCurr() + pNode->GetWCurr().Cross(fRef));
	if (bIsDifferentiable()) {
		d -= dGetP();
	}
	WorkVec.Sub(12 + 2, d);

	return WorkVec;
}

/* ImposedDisplacementPinJoint - end */
