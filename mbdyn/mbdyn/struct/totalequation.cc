/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <iostream>
#include <fstream>

#include "totalequation.h"
#include "Rot.hh"
#include "hint_impl.h"

static const char idx2xyz[] = { 'x', 'y', 'z' };

/* TotalEquation - begin */

TotalEquation::TotalEquation(unsigned int uL, const DofOwner *pDO,
	bool bPos[3], bool bVel[3],
	TplDriveCaller<Vec3> *const pDCPos[3],
	bool bRot[3], bool bAgv[3],	/* Agv stands for AnGular Velocity */
	TplDriveCaller<Vec3> *const pDCRot[3],
	const StructNode *pN1,
	const Vec3& f1Tmp, const Mat3x3& R1hTmp, const Mat3x3& R1hrTmp,
	const StructNode *pN2,
	const Vec3& f2Tmp, const Mat3x3& R2hTmp, const Mat3x3& R2hrTmp,
	flag fOut)
: Elem(uL, fOut),
Joint(uL, pDO, fOut),
pNode1(pN1), pNode2(pN2),
f1(f1Tmp), R1h(R1hTmp), R1hr(R1hrTmp),
f2(f2Tmp), R2h(R2hTmp), R2hr(R2hrTmp),
XDrv(pDCPos[0]), XPDrv(pDCPos[1]), XPPDrv(pDCPos[2]),
ThetaDrv(pDCRot[0]), OmegaDrv(pDCRot[1]), OmegaPDrv(pDCRot[2]),
nConstraints(0), nPosConstraints(0), nRotConstraints(0),
nVelConstraints(0), nAgvConstraints(0),
tilde_f1(R1h.MulTV(f1)),
M(Zero3), F(Zero3), ThetaDelta(Zero3), ThetaDeltaPrev(Zero3)
{
	/* Equations 1->3: Positions
	 * Equations 4->6: Rotations */
	
	unsigned int index = 0;

	for (unsigned int i = 0; i < 3; i++) {
		ASSERT(bPos[i] == false || bVel[i] == false);

		bPosActive[i] = bPos[i];
		bVelActive[i] = bVel[i];
		
		if (bPosActive[i]) {
			iPosIncid[nPosConstraints] = i + 1;
			iPosEqIndex[nPosConstraints] = index;
			nPosConstraints++;
			index++;
		}
		
		if (bVelActive[i]) {
			iVelIncid[nVelConstraints] = i + 1;
			iVelEqIndex[nVelConstraints] = index;
			nVelConstraints++;
			index++;
		}
	}	
	
	index = 0;
	for (unsigned int i = 0; i < 3; i++) {
		ASSERT(bRot[i] == false || bAgv[i] == false);

		bRotActive[i] = bRot[i];
		bAgvActive[i] = bAgv[i];
		
		if (bRotActive[i]) {
			iRotIncid[nRotConstraints] = i + 1;
			iRotEqIndex[nRotConstraints] = nPosConstraints + nVelConstraints + index;
			nRotConstraints++;
			index++;
		}
		
		if (bAgvActive[i]) {
			iAgvIncid[nAgvConstraints] = i + 1;
			iAgvEqIndex[nAgvConstraints] = nPosConstraints + nVelConstraints + index;
			nAgvConstraints++;
			index++;
		}
	}
	nConstraints = nPosConstraints + nRotConstraints + nVelConstraints + nAgvConstraints;

}

TotalEquation::~TotalEquation(void)
{
	NO_OP;
};

/* FIXME: velocity stuff not implemented yet */
std::ostream&
TotalEquation::DescribeDof(std::ostream& out,
	const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

	if (nPosConstraints > 0 || nVelConstraints > 0) {
		out << prefix << iIndex + 1;
		if (nPosConstraints + nVelConstraints > 1) {
			out << "->" << iIndex + nPosConstraints + nVelConstraints;
		}
		out << ": ";
		out << "reaction force(s) [";

		for (unsigned int i = 0, cnt = 0; i < 3; i++) {
			if (bPosActive[i] || bVelActive[i]) {
				cnt++;
				if (cnt > 1) {
					out << ",";
				}
				out << "F" << idx2xyz[i];
			}
		}
		out << "]" << std::endl;
	}

	if (nRotConstraints > 0 || nAgvConstraints > 0) {
		out << prefix << iIndex + nPosConstraints + nVelConstraints + 1;
		if (nRotConstraints + nAgvConstraints > 1) {
			out << "->" << iIndex + nConstraints ;
		}
		out << ": ";
		out << "reaction couple(s) [";

		for (unsigned int i = 0, cnt = 0; i < 3; i++) {
			if (bRotActive[i] || bAgvActive[i]) {
				cnt++;
				if (cnt > 1) {
					out << ",";
				}
				out << "m" << idx2xyz[i];
			}
		}
		out << "]" << std::endl;
	}

	if (bInitial) {
		iIndex += nConstraints;

		if (nPosConstraints > 0) {
			out << prefix << iIndex + 1;
			if (nPosConstraints > 1) {
				out << "->" << iIndex + nPosConstraints;
			}
			out << ": ";
			out << "reaction force(s) derivative(s) [";

			for (unsigned int i = 0, cnt = 0; i < 3; i++) {
				if (bPosActive[i]) {
					cnt++;
					if (cnt > 1) {
						out << ",";
					}
					out << "FP" << idx2xyz[i];
				}
			}
			out << "]" << std::endl;
		}

		if (nRotConstraints > 0) {
			out << prefix << iIndex + nPosConstraints + 1;
			if (nRotConstraints > 1) {
				out << "->" << iIndex + nConstraints;
			}
			out << ": ";
			out << "reaction couple(s) derivative(s) [";

			for (unsigned int i = 0, cnt = 0; i < 3; i++) {
				if (bRotActive[i]) {
					cnt++;
					if (cnt > 1) {
						out << ",";
					}
					out << "mP" << idx2xyz[i];
				}
			}
			out << "]" << std::endl;
		}
	}

	return out;
}

void
TotalEquation::DescribeDof(std::vector<std::string>& desc,
	bool bInitial, int i) const
{
	// FIXME: todo
	int ndof = 1;
	if (i == -1) {
		if (bInitial) {
			ndof = iGetNumDof();

		} else {
			ndof = iGetInitialNumDof();
		}
	}
	desc.resize(ndof);

	std::ostringstream os;
	os << "TotalEquation(" << GetLabel() << ")";

	if (i == -1) {
		std::string name(os.str());

		for (i = 0; i < ndof; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": dof(" << i + 1 << ")";
			desc[i] = os.str();
		}

	} else {
		os << ": dof(" << i + 1 << ")";
		desc[0] = os.str();
	}
}

/* FIXME: velocity stuff not implemented yet */
std::ostream&
TotalEquation::DescribeEq(std::ostream& out,
	const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

	if (nPosConstraints > 0 || nVelConstraints > 0) {
		out << prefix << iIndex + 1;
		if (nPosConstraints + nVelConstraints > 1) {
			out << "->" << iIndex + nPosConstraints + nVelConstraints;
		}
		out << ": ";
		out << "position/velocity constraint(s) [";
	
		for (unsigned int i = 0, cnt = 0; i < 3; i++) {
			if (bPosActive[i]) {
				cnt++;
				if (cnt > 1) {
					out << ",";
				}
				out << "P" << idx2xyz[i] << "1=P" << idx2xyz[i] << "2";
			} else if (bVelActive[i]) {
				cnt++;
				if (cnt > 1) {
					out << ",";
				}
				out << "V" << idx2xyz[i] << "1=V" << idx2xyz[i] << "2";
			}
		}
	
		out << "]" << std::endl;
	}

	if (nRotConstraints > 0 || nAgvConstraints > 0) {
		out << prefix << iIndex + nPosConstraints + nVelConstraints + 1;
		if (nRotConstraints + nAgvConstraints > 1) {
			out << "->" << iIndex + nConstraints ;
		}
		out << ": ";
		out << "orientation/angular velocity constraint(s) [";
	
		for (unsigned int i = 0, cnt = 0; i < 3; i++) {
			if (bRotActive[i]) {
				cnt++;
				if (cnt > 1) {
					out << ",";
				}
				out << "g" << idx2xyz[i] << "1=g" << idx2xyz[i] << "2";
			} else if (bAgvActive[i]) {
				cnt++;
				if (cnt > 1) {
					out << ",";
				}
				out << "W" << idx2xyz[i] << "1=W" << idx2xyz[i] << "2";
			}
		}

		out << "]" << std::endl;
	}

	if (bInitial) {
		iIndex += nConstraints;

		if (nPosConstraints > 0) {
			out << prefix << iIndex + 1;
			if (nPosConstraints > 1) {
				out << "->" << iIndex + nPosConstraints;
			}
			out << ": ";
			out << "velocity constraint(s) [";
		
			for (unsigned int i = 0, cnt = 0; i < 3; i++) {
				if (bPosActive[i]) {
					cnt++;
					if (cnt > 1) {
						out << ",";
					}
					out << "v" << idx2xyz[i] << "1=v" << idx2xyz[i] << "2";
				}
			}

			out << "]" << std::endl;
		}
		
		if (nRotConstraints > 0) {
			out << prefix << iIndex + nPosConstraints + 1;
			if (nRotConstraints > 1) {
				out << "->" << iIndex + nConstraints ;
			}
			out << ": ";
			out << "angular velocity constraint(s) [";
		
			for (unsigned int i = 0, cnt = 0; i < 3; i++) {
				if (bRotActive[i]) {
					cnt++;
					if (cnt > 1) {
						out << ",";
					}
					out << "w" << idx2xyz[i] << "1=w" << idx2xyz[i] << "2";
				}
			}

			out << "]" << std::endl;
		}	
	}

	return out;
}

void
TotalEquation::DescribeEq(std::vector<std::string>& desc,
	bool bInitial, int i) const
{
	// FIXME: todo
	int ndof = 1;
	if (i == -1) {
		if (bInitial) {
			ndof = iGetNumDof();

		} else {
			ndof = iGetInitialNumDof();
		}
	}
	desc.resize(ndof);

	std::ostringstream os;
	os << "TotalEquation(" << GetLabel() << ")";

	if (i == -1) {
		std::string name(os.str());

		for (i = 0; i < ndof; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": equation(" << i + 1 << ")";
			desc[i] = os.str();
		}

	} else {
		os << ": equation(" << i + 1 << ")";
		desc[0] = os.str();
	}
}

void
TotalEquation::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	if (ph) {
		for (unsigned int i = 0; i < ph->size(); i++) {
			Joint::JointHint *pjh = dynamic_cast<Joint::JointHint *>((*ph)[i]);

			if (pjh) {

				if (dynamic_cast<Joint::OffsetHint<1> *>(pjh)) {
					const Mat3x3& R1(pNode1->GetRCurr());
					Vec3 fTmp2(pNode2->GetRCurr()*f2);

					f1 = R1.MulTV(pNode2->GetXCurr() + fTmp2 - pNode1->GetXCurr());
					tilde_f1 = R1h.MulTV(f1) - XDrv.Get();

				} else if (dynamic_cast<Joint::OffsetHint<2> *>(pjh)) {
					const Mat3x3& R2(pNode2->GetRCurr());
					Mat3x3 R1(pNode1->GetRCurr()*R1h);
					Vec3 fTmp1(pNode1->GetRCurr()*f1);

					f2 = R2.MulTV(pNode1->GetXCurr() + fTmp1 - pNode2->GetXCurr() + R1*XDrv.Get());

				} else if (dynamic_cast<Joint::HingeHint<1> *>(pjh)) {
					if (dynamic_cast<Joint::PositionHingeHint<1> *>(pjh)) {
						const Mat3x3& R1(pNode1->GetRCurr());
						R1h = R1.MulTM(pNode2->GetRCurr()*R2h);

					} else if (dynamic_cast<Joint::OrientationHingeHint<1> *>(pjh)) {
						const Mat3x3& R1(pNode1->GetRCurr());
						R1hr = R1.MulTM(pNode2->GetRCurr()*R2hr);
					}

				} else if (dynamic_cast<Joint::HingeHint<2> *>(pjh)) {
					if (dynamic_cast<Joint::PositionHingeHint<2> *>(pjh)) {
						const Mat3x3& R2(pNode2->GetRCurr());
						R2h = R2.MulTM(pNode1->GetRCurr()*R1h);

					} else if (dynamic_cast<Joint::OrientationHingeHint<2> *>(pjh)) {
						const Mat3x3& R2(pNode2->GetRCurr());
						R2hr = R2.MulTM(pNode1->GetRCurr()*R1hr);
					}

				} else if (dynamic_cast<Joint::JointDriveHint<Vec3> *>(pjh)) {
					Joint::JointDriveHint<Vec3> *pjdh
						= dynamic_cast<Joint::JointDriveHint<Vec3> *>(pjh);
					pedantic_cout("TotalEquation(" << uLabel << "): "
						"creating drive from hint[" << i << "]..." << std::endl);

					TplDriveCaller<Vec3> *pDC = pjdh->pTDH->pCreateDrive(pDM);
					if (pDC == 0) {
						silent_cerr("TotalEquation(" << uLabel << "): "
							"unable to create drive "
							"after hint #" << i << std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}

					if (dynamic_cast<Joint::PositionDriveHint<Vec3> *>(pjdh)) {
						XDrv.Set(pDC);

					} else if (dynamic_cast<Joint::VelocityDriveHint<Vec3> *>(pjdh)) {
						XPDrv.Set(pDC);

					} else if (dynamic_cast<Joint::AccelerationDriveHint<Vec3> *>(pjdh)) {
						XPPDrv.Set(pDC);

					} else if (dynamic_cast<Joint::OrientationDriveHint<Vec3> *>(pjdh)) {
						ThetaDrv.Set(pDC);

					} else if (dynamic_cast<Joint::AngularVelocityDriveHint<Vec3> *>(pjdh)) {
						OmegaDrv.Set(pDC);

					} else if (dynamic_cast<Joint::AngularAccelerationDriveHint<Vec3> *>(pjdh)) {
						OmegaPDrv.Set(pDC);

					} else {
						delete pDC;
					}

				} else if (dynamic_cast<Joint::ReactionsHint *>(pjh)) {
					/* TODO */
				}
				continue;
			}
		}
	}
}

Hint *
TotalEquation::ParseHint(DataManager *pDM, const char *s) const
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

	} else if (strncasecmp(s, "position-hinge{" /*}*/, STRLENOF("position-hinge{" /*}*/)) == 0) {
		s += STRLENOF("position-hinge{" /*}*/);

		if (strcmp(&s[1], /*{*/ "}") != 0) {
			return 0;
		}

		switch (s[0]) {
		case '1':
			return new Joint::HingeHint<1>;

		case '2':
			return new Joint::HingeHint<2>;
		}

	} else if (strncasecmp(s, "position-drive3{" /*}*/, STRLENOF("position-drive3{" /*}*/)) == 0) {
		s += STRLENOF("position-");

		Hint *pH = ::ParseHint(pDM, s);
		if (pH) {
			TplDriveHint<Vec3> *pTDH = dynamic_cast<TplDriveHint<Vec3> *>(pH);
			if (pTDH) {
				return new PositionDriveHint<Vec3>(pTDH);
			}
		}
		return 0;

	} else if (strncasecmp(s, "orientation-hinge{" /*}*/, STRLENOF("orientation-hinge{" /*}*/)) == 0) {
		s += STRLENOF("orientation-hinge{" /*}*/);

		if (strcmp(&s[1], /*{*/ "}") != 0) {
			return 0;
		}

		switch (s[0]) {
		case '1':
			return new Joint::HingeHint<1>;

		case '2':
			return new Joint::HingeHint<2>;
		}

	} else if (strncasecmp(s, "orientation-drive3{" /*}*/, STRLENOF("orientation-drive3{" /*}*/)) == 0) {
		s += STRLENOF("orientation-");

		Hint *pH = ::ParseHint(pDM, s);
		if (pH) {
			TplDriveHint<Vec3> *pTDH = dynamic_cast<TplDriveHint<Vec3> *>(pH);
			if (pTDH) {
				return new OrientationDriveHint<Vec3>(pTDH);
			}
		}
		return 0;
	}

	return 0;
}

void
TotalEquation::AfterConvergence(const VectorHandler& /* X */ ,
		const VectorHandler& /* XP */ )
{
	ThetaDeltaPrev = Unwrap(ThetaDeltaPrev, ThetaDelta);
}

void
TotalEquation::AfterConvergence(const VectorHandler& /* X */ ,
		const VectorHandler& /* XP */, const VectorHandler& /* XP */)
{
	ThetaDeltaPrev = Unwrap(ThetaDeltaPrev, ThetaDelta);
}

/* Contributo al file di restart */
/* FIXME: velocity stuffs not implemented yet */
std::ostream&
TotalEquation::Restart(std::ostream& out) const
{
	Joint::Restart(out) << ", total equation, "
		<< pNode1->GetLabel() << ", "
			<< "position, " << f1.Write(out, ", ") << ", "
			<< "position orientation, "
				"1 , ", R1h.GetVec(1).Write(out, ", ") << ", "
				"2 , ", R1h.GetVec(2).Write(out, ", ") << ", "
			<< "rotation orientation, "
				"1 , ", R1hr.GetVec(1).Write(out, ", ") << ", "
				"2 , ", R1hr.GetVec(2).Write(out, ", ") << ", "
		<< pNode2->GetLabel() << ", "
			<< "position, " << f2.Write(out, ", ") << ", "
			<< "position orientation, "
				"1 , ", R2h.GetVec(1).Write(out, ", ") << ", "
				"2 , ", R2h.GetVec(2).Write(out, ", ") << ", "
			<< "rotation orientation, "
				"1 , ", R2hr.GetVec(1).Write(out, ", ") << ", "
				"2 , ", R2hr.GetVec(2).Write(out, ", ");

	if (bPosActive[0] || bPosActive[1] || bPosActive[2]) {
		out << ", position constraint";
		for (unsigned i = 0; i < 3; i++) {
			if (bPosActive[i]) {
				out << ", active";
			} else {
				out << ", inactive";
			}
		}
	}

	if (bRotActive[0] || bRotActive[1] || bRotActive[2]) {
		out << ", orientation constraint";
		for (unsigned i = 0; i < 3; i++) {
			if (bRotActive[i]) {
				out << ", active";
			} else {
				out << ", inactive";
			}
		}
	}

	return out << std::endl;
}

/* Assemblaggio jacobiano */
VariableSubMatrixHandler&
TotalEquation::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	/*
	 See tecman.pdf for details
	 */
	 
	DEBUGCOUT("Entering TotalEquation::AssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Ridimensiona la sottomatrice in base alle esigenze */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici delle varie incognite */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iFirstReactionIndex = iGetFirstIndex();

	/* Setta gli indici delle equazioni */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	for (unsigned int iCnt = 1; iCnt <= nConstraints; iCnt++) {
		WM.PutRowIndex(iCnt, iFirstReactionIndex + iCnt);
	}

	/* Recupera i dati che servono */
	Mat3x3 R1(pNode1->GetRRef()*R1h);
	Mat3x3 R1r(pNode1->GetRRef()*R1hr);
	Vec3 b2(pNode2->GetRCurr()*f2);
	Vec3 b1(pNode2->GetXCurr() + b2 - pNode1->GetXCurr());

/* Phi/q is the only contribution... */

	Mat3x3 b1Cross_R1(b1.Cross(R1)); // = [ b1 x ] * R1
	Mat3x3 b2Cross_R1(b2.Cross(R1)); // = [ b2 x ] * R1

	for (unsigned iCnt = 0 ; iCnt < nPosConstraints; iCnt++) {
		Vec3 vR1(R1.GetVec(iPosIncid[iCnt]));
		Vec3 vb1Cross_R1(b1Cross_R1.GetVec(iPosIncid[iCnt]));
		Vec3 vb2Cross_R1(b2Cross_R1.GetVec(iPosIncid[iCnt]));

		/* Constraint, node 1 */
      		WM.SubT(1 + iPosEqIndex[iCnt], 1, vR1);
      		WM.SubT(1 + iPosEqIndex[iCnt], 3 + 1, vb1Cross_R1);

		/* Constraint, node 2 */
      		WM.AddT(1 + iPosEqIndex[iCnt], 6 + 1, vR1);
      		WM.AddT(1 + iPosEqIndex[iCnt], 9 + 1, vb2Cross_R1);
	}
	
	for (unsigned iCnt = 0 ; iCnt < nRotConstraints; iCnt++) {
		Vec3 vR1(R1r.GetVec(iRotIncid[iCnt]));

		/* Constraint, node 1 */
      		WM.SubT(1 + iRotEqIndex[iCnt], 3 + 1, vR1);

		/* Constraint, node 2 */
      		WM.AddT(1 + iRotEqIndex[iCnt], 9 + 1, vR1);
	}

	if (nVelConstraints > 0) {
		Mat3x3 Omega1Cross_R1(pNode1->GetWCurr().Cross(R1));
		Mat3x3 	Tmp12 =	(
				Mat3x3(MatCrossCross, pNode1->GetWCurr(), -b1)
				+ Mat3x3(MatCrossCross, b2, pNode1->GetWCurr())
				- Mat3x3(MatCross, pNode2->GetVCurr())
				- Mat3x3(MatCrossCross, pNode2->GetWCurr(), b2)
				+ Mat3x3(MatCross, pNode1->GetVCurr())
				) * R1;
		Mat3x3 Tmp22(MatCrossCross, pNode2->GetWCurr(), b2);

		for (unsigned iCnt = 0 ; iCnt < nVelConstraints; iCnt++) {
			Vec3 vR1(R1.GetVec(iVelIncid[iCnt]));
			Vec3 vb1Cross_R1(b1Cross_R1.GetVec(iVelIncid[iCnt]));
			Vec3 vb2Cross_R1(b2Cross_R1.GetVec(iVelIncid[iCnt]));

			Vec3 vOmega1Cross_R1(Omega1Cross_R1.GetVec(iVelIncid[iCnt]));
			Vec3 vTmp12(Tmp12.GetVec(iVelIncid[iCnt]));	
			Vec3 vTmp22(Tmp22.GetVec(iVelIncid[iCnt]));	
			
			/* Constraint, node 1 */
			/* The same as in position constraint*/
      			WM.SubT(1 + iVelEqIndex[iCnt], 1, vR1);				// delta_v1
      			WM.SubT(1 + iVelEqIndex[iCnt], 3 + 1, vb1Cross_R1);		// delta_W1
			/* New contributions, related to delta_x1 and delta_g1 */
     			WM.SubT(1 + iVelEqIndex[iCnt], 1, vOmega1Cross_R1 * dCoef);	// delta_x1
     			WM.AddT(1 + iVelEqIndex[iCnt], 3 + 1, vTmp12 * dCoef);		// delta_g1
			
			/* Constraint, node 2 */
			/* The same as in position constraint*/
      			WM.AddT(1 + iVelEqIndex[iCnt], 6 + 1, vR1);			// delta_v2
      			WM.AddT(1 + iVelEqIndex[iCnt], 9 + 1, vb2Cross_R1);		// delta_W2
			/* New contributions, related to delta_x1 and delta_g1 */
     			WM.AddT(1 + iVelEqIndex[iCnt], 6 + 1, vOmega1Cross_R1 * dCoef);	// delta_x2
     			WM.AddT(1 + iVelEqIndex[iCnt], 9 + 1, vTmp22 * dCoef);		// delta_g2
		}
	}
	
	if(nAgvConstraints > 0)	{
		Mat3x3 DeltaWCross_R1((pNode1->GetWCurr() - pNode2->GetWCurr()).Cross(R1r));
		for (unsigned iCnt = 0 ; iCnt < nAgvConstraints; iCnt++) {
			Vec3 vR1(R1r.GetVec(iAgvIncid[iCnt]));
			Vec3 vDeltaWCross_R1(DeltaWCross_R1.GetVec(iAgvIncid[iCnt]));

			/* Constraint, node 1 */
      			WM.SubT(1 + iAgvEqIndex[iCnt], 3 + 1, vR1);	// delta_W1
			/* New contribution, related to delta_g1 */
      			WM.SubT(1 + iAgvEqIndex[iCnt], 3 + 1, vDeltaWCross_R1 * dCoef);	// delta_g1
	
			/* Constraint, node 2 */
      			WM.AddT(1 + iAgvEqIndex[iCnt], 9 + 1, vR1);// delta_W2
		}
	}

	return WorkMat;
}

/* Assemblaggio residuo */
SubVectorHandler&
TotalEquation::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering TotalEquation::AssRes()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici */
	integer iFirstReactionIndex = iGetFirstIndex();

	/* Indici del vincolo */
	for (unsigned int iCnt = 1; iCnt <= nConstraints; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstReactionIndex + iCnt);
	}

	Vec3 b2(pNode2->GetRCurr()*f2);
	Vec3 b1(pNode2->GetXCurr() + b2 - pNode1->GetXCurr());

	Mat3x3 R1 = pNode1->GetRCurr()*R1h;
	Mat3x3 R1r = pNode1->GetRCurr()*R1hr;
	Mat3x3 R2r = pNode2->GetRCurr()*R2hr;

	Vec3 XDelta = R1.MulTV(b1) - tilde_f1 - XDrv.Get();
	Vec3 VDelta = R1.MulTV(
			pNode2->GetVCurr()
			- b2.Cross(pNode2->GetWCurr())
			- pNode1->GetVCurr()
			+ b1.Cross(pNode1->GetWCurr())
			) 
			- XDrv.Get(); 	
			
	Mat3x3 R0T = RotManip::Rot(-ThetaDrv.Get());	// -Theta0 to get R0 transposed
	Mat3x3 RDelta = R1r.MulTM(R2r*R0T);
	ThetaDelta = RotManip::VecRot(RDelta);
	
	Vec3 WDelta = R1r.MulTV(pNode2->GetWCurr() - pNode1->GetWCurr()) - ThetaDrv.Get();

	/* Constraint equations are divided by dCoef
	 * ONLY if the constraint is on position */
	ASSERT(dCoef != 0.);

	/* Position constraint:  */
	for (unsigned iCnt = 0; iCnt < nPosConstraints; iCnt++) {
		WorkVec.PutCoef(1 + iPosEqIndex[iCnt],
			-(XDelta(iPosIncid[iCnt])/dCoef));
	}

	/* Rotation constraints: */
	for (unsigned iCnt = 0; iCnt < nRotConstraints; iCnt++) {
		WorkVec.PutCoef(1 + iRotEqIndex[iCnt],
			-(ThetaDelta(iRotIncid[iCnt])/dCoef));
	}

	/* Linear Velocity Constraint */
	for (unsigned iCnt = 0; iCnt < nVelConstraints; iCnt++) {
		WorkVec.PutCoef(1 + iVelEqIndex[iCnt],
			-(VDelta(iVelIncid[iCnt])));
	}
		
	/* Angular Velocity Constraint */
	for (unsigned iCnt = 0; iCnt < nAgvConstraints; iCnt++) {
		WorkVec.PutCoef(1 + iAgvEqIndex[iCnt],
			-(WDelta(iAgvIncid[iCnt])));
	}

	return WorkVec;
}

void
TotalEquation::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		const Vec3& X1(pNode1->GetXCurr());
		const Vec3& X2(pNode2->GetXCurr());
		const Mat3x3& R1(pNode1->GetRCurr());
		const Mat3x3& R2(pNode2->GetRCurr());
		const Vec3& V1(pNode1->GetVCurr());
		const Vec3& V2(pNode2->GetVCurr());
		const Vec3& Omega1(pNode1->GetWCurr());
		const Vec3& Omega2(pNode2->GetWCurr());

		Mat3x3 R1Tmp(R1*R1h);
		Mat3x3 R1rTmp(R1*R1hr);
		Mat3x3 R2rTmp(R2*R2hr);
		Mat3x3 RTmp(R1rTmp.MulTM(R2rTmp));
		Vec3 b2(R2*f2);
		Vec3 b1(X2 + b2 - X1);
		Vec3 X(R1Tmp.MulTV(b1) - tilde_f1);

		Joint::Output(OH.Joints(), "TotalEquation", GetLabel(),
			Zero3, Zero3, Zero3, Zero3)
			<< " " << X
			<< " " << RotManip::VecRot(RTmp)
			<< " " << R1Tmp.MulTV(V2 + Omega2.Cross(b2) - V1 - Omega1.Cross(b1))
			<< " " << R1rTmp.MulTV(Omega2 - Omega1)
			<< std::endl;
			// accelerations?
	}
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
TotalEquation::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{

	/* Per ora usa la matrice piena; eventualmente si puo'
	 * passare a quella sparsa quando si ottimizza */
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Equazioni: vedi joints.dvi */

	/*	 equazioni ed incognite
	 * F1					  Delta_x1	1	   
	 * M1					  Delta_g1	3 + 1 
	 * FP1  				  Delta_xP1	6 + 1
	 * MP1  				  Delta_w1	9 + 1 
	 * F2					  Delta_x2	12 + 1
	 * M2					  Delta_g2	15 + 1
	 * FP2  				  Delta_xP2	18 + 1  
	 * MP2  				  Delta_w2	21 + 1 
	 * vincolo spostamento  		  Delta_F	24 + 1 
	 * vincolo rotazione			  Delta_M	24 + nPosConstraints  
	 * derivata vincolo spostamento 	  Delta_FP	24 + nConstraints  
	 * derivata vincolo rotazione		  Delta_MP	24 + nConstraints + nPosConstraints  
	 */
	
	
	/* Indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstVelIndex = iNode1FirstPosIndex + 6;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstVelIndex = iNode2FirstPosIndex + 6;
	integer iFirstReactionIndex = iGetFirstIndex();
	integer iReactionPrimeIndex = iFirstReactionIndex + nConstraints;


	/* Setta gli indici dei nodi */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode1FirstVelIndex + iCnt);
		WM.PutColIndex(12 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutColIndex(18 + iCnt, iNode2FirstVelIndex + iCnt);
	}

	/* Setta gli indici delle reazioni */
	for (unsigned int iCnt = 1; iCnt <=  nConstraints; iCnt++) {
		WM.PutRowIndex(24 + iCnt, iFirstReactionIndex + iCnt);
		WM.PutRowIndex(24 + nConstraints + iCnt, iReactionPrimeIndex + iCnt);
	}
	
	/* Recupera i dati che servono */
	Mat3x3 R1(pNode1->GetRRef()*R1h);
	Mat3x3 R1r(pNode1->GetRRef()*R1hr);
	
	Vec3 b2(pNode2->GetRCurr()*f2);
	Vec3 b1(pNode2->GetXCurr() + b2 - pNode1->GetXCurr());
	
	Vec3 Omega1(pNode1->GetWCurr());
	Vec3 Omega2(pNode2->GetWCurr());
	
	Vec3 b1Prime(pNode2->GetVCurr() + Omega2.Cross(b2) - pNode1->GetVCurr());
	
	/* Constraints: Add only active rows/columns*/	
	
	/* Positions contribution:*/

	Mat3x3 b1Cross_R1(b1.Cross(R1)); // = [ b1 x ] * R1
	Mat3x3 b2Cross_R1(b2.Cross(R1)); // = [ b2 x ] * R1

	for (unsigned iCnt = 0 ; iCnt < nPosConstraints; iCnt++) {
		Vec3 vR1(R1.GetVec(iPosIncid[iCnt]));
		Vec3 vb1Cross_R1(b1Cross_R1.GetVec(iPosIncid[iCnt]));
		Vec3 vb2Cross_R1(b2Cross_R1.GetVec(iPosIncid[iCnt]));

		/* Constraint, node 1 */
      		WM.SubT(24 + 1 + iCnt, 1, vR1);	// * Delta_x1
      		WM.SubT(24 + 1 + iCnt, 3 + 1, vb1Cross_R1);	// * Delta_g1

		/* d/dt(Constraint), node 1 */
      		WM.SubT(24 + 1 + nConstraints + iCnt, 6 + 1, vR1);	// * Delta_xP1
      		WM.SubT(24 + 1 + nConstraints + iCnt, 9 + 1, vb1Cross_R1);	// * Delta_W1
		
		/* Constraint, node 2 */
      		WM.AddT(24 + 1 + iCnt, 12 + 1, vR1);	// * Delta_x2
      		WM.AddT(24 + 1 + iCnt, 15 + 1, vb2Cross_R1);	// * Delta_g2
	
		/* d/dt(Constraint), node 2 */
      		WM.AddT(24 + 1 + nConstraints +  iCnt, 18 + 1, vR1);	// * Delta_xP2
      		WM.AddT(24 + 1 + nConstraints +  iCnt, 21 + 1, vb2Cross_R1);	// * Delta_W2
	}

	for (unsigned iCnt = 0 ; iCnt < nRotConstraints; iCnt++) {
		Vec3 vR1(R1r.GetVec(iRotIncid[iCnt]));

		/* Constraint, node 1 */
      		WM.SubT(24 + 1 + nPosConstraints + iCnt, 3 + 1, vR1);	// * Delta_g1

		/* d/dt(Constraint), node 1 */
      		WM.SubT(24 + 1 + nConstraints + nPosConstraints + iCnt, 9 + 1, vR1);	// * Delta_W1

		/* Constraint, node 2 */
      		WM.AddT(24 + 1 + nPosConstraints +  iCnt, 15 + 1, vR1);	// * Delta_g2

		/* d/dt(Constraint), node 2 */
      		WM.AddT(24 + 1 + nConstraints + nPosConstraints + iCnt, 21 + 1, vR1);	// * Delta_W2
	}

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
TotalEquation::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{

	DEBUGCOUT("Entering TotalEquation::InitialAssRes()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstVelIndex = iNode1FirstPosIndex + 6;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstVelIndex = iNode2FirstPosIndex + 6;
	integer iFirstReactionIndex = iGetFirstIndex();
	integer iReactionPrimeIndex = iFirstReactionIndex + nConstraints;

	/* Setta gli indici */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode1FirstVelIndex + iCnt);
		WorkVec.PutRowIndex(12 + iCnt, iNode2FirstPosIndex + iCnt);
		WorkVec.PutRowIndex(18 + iCnt, iNode2FirstVelIndex + iCnt);
	}

	for (unsigned int iCnt = 1; iCnt <= nConstraints; iCnt++)	{
		WorkVec.PutRowIndex(24 + iCnt, iFirstReactionIndex + iCnt);
		WorkVec.PutRowIndex(24 + nConstraints + iCnt, iReactionPrimeIndex + iCnt);
	}

	/* Recupera i dati */
	/* Recupera i dati che servono */
	Mat3x3 R1(pNode1->GetRRef()*R1h);
	Mat3x3 R1r(pNode1->GetRRef()*R1hr);
	Mat3x3 R2r(pNode2->GetRCurr()*R2hr);
	
	Vec3 b2(pNode2->GetRCurr()*f2);
	Vec3 b1(pNode2->GetXCurr() + b2 - pNode1->GetXCurr());
	
	Vec3 Omega1(pNode1->GetWCurr());
	Vec3 Omega2(pNode2->GetWCurr());
	
	Vec3 b1Prime(pNode2->GetVCurr() + Omega2.Cross(b2) - pNode1->GetVCurr());
	

	/* Constraint Equations */
	
	Vec3 XDelta = R1.MulTV(b1) - tilde_f1 - XDrv.Get();
	Vec3 XDeltaPrime = R1.MulTV(b1Prime + b1.Cross(Omega1));
	
	if(XDrv.bIsDifferentiable())	{
		XDeltaPrime -= XDrv.GetP();
	}
	
	Mat3x3 R0T = RotManip::Rot(-ThetaDrv.Get());	// -Theta0 to get R0 transposed
	Mat3x3 RDelta = R1r.MulTM(R2r*R0T);
	ThetaDelta = RotManip::VecRot(RDelta);
	Vec3 ThetaDeltaPrime = R1r.MulTV(Omega2 - Omega1);

	if(ThetaDrv.bIsDifferentiable())	{
		ThetaDeltaPrime -= RDelta*ThetaDrv.GetP();
	}


	/* Position constraint:  */
	for (unsigned iCnt = 0; iCnt < nPosConstraints; iCnt++) {
		WorkVec.PutCoef(24 + 1 + iCnt, -XDelta(iPosIncid[iCnt]));
		WorkVec.PutCoef(24 + 1 + nConstraints + iCnt, -XDeltaPrime(iPosIncid[iCnt]));
	}

	/* Rotation constraints: */
	for (unsigned iCnt = 0; iCnt < nRotConstraints; iCnt++) {
		WorkVec.PutCoef(24 + 1 + nPosConstraints + iCnt, -ThetaDelta(iRotIncid[iCnt]));
		WorkVec.PutCoef(24 + 1 + nPosConstraints + nConstraints +  iCnt, -ThetaDeltaPrime(iRotIncid[iCnt]));
	}

return WorkVec;
}


unsigned int
TotalEquation::iGetNumPrivData(void) const
{
	return 24;
}

unsigned int
TotalEquation::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);

	unsigned int off = 0;

	switch (s[0]) {
	case 'p':
		/* relative position */
		break;

	case 'r':
		/* relative orientation */
		off += 3;
		break;

	case 'F':
		/* force */
		off += 6;
		break;

	case 'M':
		/* moment */
		off += 9;
		break;

	case 'd':
		/* imposed relative position */
		off += 12;
		break;

	case 't':
		/* imposed relative orientation */
		off += 15;
		break;
	
	case 'v':
		/* relative linear velocity */
		off += 18;
		break;

	case 'w':
		/* relative angular velocity */
		off += 21;
		break;


	default:
		return 0;
	}

	switch (s[1]) {
	case 'x':
		return off + 1;

	case 'y':
		return off + 2;

	case 'z':
		return off + 3;
	}

	return 0;
}

doublereal
TotalEquation::dGetPrivData(unsigned int i) const
{
	switch (i) {
	case 1:
	case 2:
	case 3: {
		Vec3 x(pNode1->GetRCurr().MulTV(
			pNode2->GetXCurr() + pNode2->GetRCurr()*f2
				- pNode1->GetXCurr()) - f1);
			return R1h.GetVec(i)*x;
		}

	case 4:
	case 5:
	case 6: {
		Vec3 Theta(Unwrap(ThetaDeltaPrev, ThetaDelta) + ThetaDrv.Get());
		return Theta(i - 3);
		}

	case 7:
	case 8:
	case 9:
		return 0.;

	case 10:
	case 11:
	case 12:
		return 0.;

	case 13:
	case 14:
	case 15:
		if (!bPosActive[i - 13]) {
			return 0.;
		}
		return XDrv.Get()(i - 12);

	case 16:
	case 17:
	case 18:
		if (!bRotActive[i - 16]) {
			return 0.;
		}
		return ThetaDrv.Get()(i - 15);
	case 19:
	case 20:
	case 21:
		// FIXME: blind fix
		{
		Vec3 v(pNode1->GetRCurr().MulTV(
			(pNode2->GetVCurr() + pNode2->GetWCurr().Cross(pNode2->GetRCurr()*f2)
				- pNode1->GetVCurr()) - 
			pNode1->GetWCurr().Cross(pNode2->GetXCurr() + 
				pNode2->GetRCurr()*f2
				- pNode1->GetXCurr() - f1) 
			)
		);
		
			return R1h.GetVec(i-18)*v;
		}
	case 22:
	case 23:
	case 24:
		{
		Vec3 W(pNode1->GetRCurr().MulTV(pNode2->GetWCurr() - pNode1->GetWCurr()))
		;
		
			return R1hr.GetVec(i-21)*W;
		}



	default:
		ASSERT(0);
	}

	return 0.;
}


/* TotalEquation - end */


/* TotalReaction - begin */

TotalReaction::TotalReaction(unsigned int uL, const DofOwner *pDO,
	bool bPos[3],
	bool bRot[3],
	const StructNode *pN1,
	const Vec3& f1Tmp, const Mat3x3& R1hTmp, const Mat3x3& R1hrTmp,
	const StructNode *pN2,
	const Vec3& f2Tmp, const Mat3x3& R2hTmp, const Mat3x3& R2hrTmp,
	TotalEquation* t_elm, 
	flag fOut)
: Elem(uL, fOut),
Joint(uL, pDO, fOut),
pNode1(pN1), pNode2(pN2),
f1(f1Tmp), R1h(R1hTmp), R1hr(R1hrTmp),
f2(f2Tmp), R2h(R2hTmp), R2hr(R2hrTmp),
nConstraints(0), nPosConstraints(0), nRotConstraints(0),
tilde_f1(R1h.MulTV(f1)),
M(Zero3), F(Zero3), ThetaDelta(Zero3), ThetaDeltaPrev(Zero3),
total_equation_element(t_elm)
{
	/* Equations 1->3: Positions
	 * Equations 4->6: Rotations */

	for (unsigned int i = 0; i < 3; i++) {
		bPosActive[i] = bPos[i];
		bRotActive[i] = bRot[i];
		if (bPosActive[i]) {
			iPosIncid[nPosConstraints] = i + 1;
			nPosConstraints++;
		}
		if (bRotActive[i]) {
			iRotIncid[nRotConstraints] = i + 1;
			nRotConstraints++;
		}
	}
	nConstraints = nPosConstraints + nRotConstraints;
	if (nConstraints != total_equation_element->nConstraints) {
		silent_cerr("Error: TotalReaction element " <<
			GetLabel() << " has " << nConstraints <<
			" reactions, while the corresponding TotalEquation joint "
			<< total_equation_element->GetLabel() <<
			" defines only " << total_equation_element->nConstraints <<
			" contraints" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

TotalReaction::~TotalReaction(void)
{
	NO_OP;
};

std::ostream&
TotalReaction::DescribeDof(std::ostream& out,
	const char *prefix, bool bInitial) const
{
	integer iIndex = total_equation_element->iGetFirstIndex();

	if (nPosConstraints > 0) {
		out << prefix << iIndex + 1;
		if (nPosConstraints > 1) {
			out << "->" << iIndex + nPosConstraints;
		}
		out << ": ";
		out << "reaction force(s) [";

		for (unsigned int i = 0, cnt = 0; i < 3; i++) {
			if (bPosActive[i]) {
				cnt++;
				if (cnt > 1) {
					out << ",";
				}
				out << "F" << idx2xyz[i];
			}
		}
		out << "]" << std::endl;
	}

	if (nRotConstraints > 0) {
		out << prefix << iIndex + nPosConstraints + 1;
		if (nRotConstraints > 1) {
			out << "->" << iIndex + nConstraints ;
		}
		out << ": ";
		out << "reaction couple(s) [";

		for (unsigned int i = 0, cnt = 0; i < 3; i++) {
			if (bRotActive[i]) {
				cnt++;
				if (cnt > 1) {
					out << ",";
				}
				out << "m" << idx2xyz[i];
			}
		}
		out << "]" << std::endl;
	}

	if (bInitial) {
		iIndex += nConstraints;

		if (nPosConstraints > 0) {
			out << prefix << iIndex + 1;
			if (nPosConstraints > 1) {
				out << "->" << iIndex + nPosConstraints;
			}
			out << ": ";
			out << "reaction force(s) derivative(s) [";

			for (unsigned int i = 0, cnt = 0; i < 3; i++) {
				if (bPosActive[i]) {
					cnt++;
					if (cnt > 1) {
						out << ",";
					}
					out << "FP" << idx2xyz[i];
				}
			}
			out << "]" << std::endl;
		}

		if (nRotConstraints > 0) {
			out << prefix << iIndex + nPosConstraints + 1;
			if (nRotConstraints > 1) {
				out << "->" << iIndex + nConstraints;
			}
			out << ": ";
			out << "reaction couple(s) derivative(s) [";

			for (unsigned int i = 0, cnt = 0; i < 3; i++) {
				if (bRotActive[i]) {
					cnt++;
					if (cnt > 1) {
						out << ",";
					}
					out << "mP" << idx2xyz[i];
				}
			}
			out << "]" << std::endl;
		}
	}
	return out;
}

void
TotalReaction::DescribeDof(std::vector<std::string>& desc,
	bool bInitial, int i) const
{
	desc.resize(0);
}

void
TotalReaction::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	if (ph) {
		for (unsigned int i = 0; i < ph->size(); i++) {
			Joint::JointHint *pjh = dynamic_cast<Joint::JointHint *>((*ph)[i]);

			if (pjh) {

				if (dynamic_cast<Joint::OffsetHint<1> *>(pjh)) {
					const Mat3x3& R1(pNode1->GetRCurr());
					Vec3 fTmp2(pNode2->GetRCurr()*f2);

					f1 = R1.MulTV(pNode2->GetXCurr() + fTmp2 - pNode1->GetXCurr());
					tilde_f1 = R1h.MulTV(f1);
					

				} else if (dynamic_cast<Joint::OffsetHint<2> *>(pjh)) {
					const Mat3x3& R2(pNode2->GetRCurr());
					Vec3 fTmp1(pNode1->GetRCurr()*f1);

					f2 = R2.MulTV(pNode1->GetXCurr() + fTmp1 - pNode2->GetXCurr());

				} else if (dynamic_cast<Joint::HingeHint<1> *>(pjh)) {
					if (dynamic_cast<Joint::PositionHingeHint<1> *>(pjh)) {
						const Mat3x3& R1(pNode1->GetRCurr());
						R1h = R1.MulTM(pNode2->GetRCurr()*R2h);

					} else if (dynamic_cast<Joint::OrientationHingeHint<1> *>(pjh)) {
						const Mat3x3& R1(pNode1->GetRCurr());
						R1hr = R1.MulTM(pNode2->GetRCurr()*R2hr);
					}

				} else if (dynamic_cast<Joint::HingeHint<2> *>(pjh)) {
					if (dynamic_cast<Joint::PositionHingeHint<2> *>(pjh)) {
						const Mat3x3& R2(pNode2->GetRCurr());
						R2h = R2.MulTM(pNode1->GetRCurr()*R1h);

					} else if (dynamic_cast<Joint::OrientationHingeHint<2> *>(pjh)) {
						const Mat3x3& R2(pNode2->GetRCurr());
						R2hr = R2.MulTM(pNode1->GetRCurr()*R1hr);
					}

#if 0
 				} else if (dynamic_cast<Joint::JointDriveHint<Vec3> *>(pjh)) {
 					Joint::JointDriveHint<Vec3> *pjdh
 						= dynamic_cast<Joint::JointDriveHint<Vec3> *>(pjh);
 					pedantic_cout("TotalReaction(" << uLabel << "): "
 						"creating drive from hint[" << i << "]..." << std::endl);
 
 					TplDriveCaller<Vec3> *pDC = pjdh->pTDH->pCreateDrive(pDM);
 					if (pDC == 0) {
 						silent_cerr("TotalReaction(" << uLabel << "): "
 							"unable to create drive "
 							"after hint #" << i << std::endl);
 						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
 					}
 
 					if (dynamic_cast<Joint::PositionDriveHint<Vec3> *>(pjdh)) {
 						XDrv.Set(pDC);
 
 					} else if (dynamic_cast<Joint::VelocityDriveHint<Vec3> *>(pjdh)) {
 						XPDrv.Set(pDC);
 
 					} else if (dynamic_cast<Joint::AccelerationDriveHint<Vec3> *>(pjdh)) {
 						XPPDrv.Set(pDC);
 
 					} else if (dynamic_cast<Joint::OrientationDriveHint<Vec3> *>(pjdh)) {
 						ThetaDrv.Set(pDC);
 
 					} else if (dynamic_cast<Joint::AngularVelocityDriveHint<Vec3> *>(pjdh)) {
 						OmegaDrv.Set(pDC);
 
 					} else if (dynamic_cast<Joint::AngularAccelerationDriveHint<Vec3> *>(pjdh)) {
 						OmegaPDrv.Set(pDC);
 
 					} else {
 						delete pDC;
 					}
#endif
 
				} else if (dynamic_cast<Joint::ReactionsHint *>(pjh)) {
					/* TODO */
				}

				continue;
			}
		}
	}
}

Hint *
TotalReaction::ParseHint(DataManager *pDM, const char *s) const
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

	} else if (strncasecmp(s, "position-hinge{" /*}*/, STRLENOF("position-hinge{" /*}*/)) == 0) {
		s += STRLENOF("position-hinge{" /*}*/);

		if (strcmp(&s[1], /*{*/ "}") != 0) {
			return 0;
		}

		switch (s[0]) {
		case '1':
			return new Joint::HingeHint<1>;

		case '2':
			return new Joint::HingeHint<2>;
		}

	} else if (strncasecmp(s, "position-drive3{" /*}*/, STRLENOF("position-drive3{" /*}*/)) == 0) {
		s += STRLENOF("position-");

		Hint *pH = ::ParseHint(pDM, s);
		if (pH) {
			TplDriveHint<Vec3> *pTDH = dynamic_cast<TplDriveHint<Vec3> *>(pH);
			if (pTDH) {
				return new PositionDriveHint<Vec3>(pTDH);
			}
		}
		return 0;

	} else if (strncasecmp(s, "orientation-hinge{" /*}*/, STRLENOF("orientation-hinge{" /*}*/)) == 0) {
		s += STRLENOF("orientation-hinge{" /*}*/);

		if (strcmp(&s[1], /*{*/ "}") != 0) {
			return 0;
		}

		switch (s[0]) {
		case '1':
			return new Joint::HingeHint<1>;

		case '2':
			return new Joint::HingeHint<2>;
		}

	} else if (strncasecmp(s, "orientation-drive3{" /*}*/, STRLENOF("orientation-drive3{" /*}*/)) == 0) {
		s += STRLENOF("orientation-");

		Hint *pH = ::ParseHint(pDM, s);
		if (pH) {
			TplDriveHint<Vec3> *pTDH = dynamic_cast<TplDriveHint<Vec3> *>(pH);
			if (pTDH) {
				return new OrientationDriveHint<Vec3>(pTDH);
			}
		}
		return 0;
	}

	return 0;
}

void
TotalReaction::AfterConvergence(const VectorHandler& /* X */ ,
		const VectorHandler& /* XP */ )
{
	ThetaDeltaPrev = Unwrap(ThetaDeltaPrev, ThetaDelta);
}

void
TotalReaction::AfterConvergence(const VectorHandler& /* X */ ,
		const VectorHandler& /* XP */, const VectorHandler& /* XP */)
{
	ThetaDeltaPrev = Unwrap(ThetaDeltaPrev, ThetaDelta);
}

/* Contributo al file di restart */
std::ostream&
TotalReaction::Restart(std::ostream& out) const
{
	Joint::Restart(out) << ", total equation, "
		<< pNode1->GetLabel() << ", "
			<< "position, " << f1.Write(out, ", ") << ", "
			<< "position orientation, "
				"1 , ", R1h.GetVec(1).Write(out, ", ") << ", "
				"2 , ", R1h.GetVec(2).Write(out, ", ") << ", "
			<< "rotation orientation, "
				"1 , ", R1hr.GetVec(1).Write(out, ", ") << ", "
				"2 , ", R1hr.GetVec(2).Write(out, ", ") << ", "
		<< pNode2->GetLabel() << ", "
			<< "position, " << f2.Write(out, ", ") << ", "
			<< "position orientation, "
				"1 , ", R2h.GetVec(1).Write(out, ", ") << ", "
				"2 , ", R2h.GetVec(2).Write(out, ", ") << ", "
			<< "rotation orientation, "
				"1 , ", R2hr.GetVec(1).Write(out, ", ") << ", "
				"2 , ", R2hr.GetVec(2).Write(out, ", ");

	if (bPosActive[0] || bPosActive[1] || bPosActive[2]) {
		out << ", position constraint";
		for (unsigned i = 0; i < 3; i++) {
			if (bPosActive[i]) {
				out << ", active";
			} else {
				out << ", inactive";
			}
		}
	}

	if (bRotActive[0] || bRotActive[1] || bRotActive[2]) {
		out << ", orientation constraint";
		for (unsigned i = 0; i < 3; i++) {
			if (bRotActive[i]) {
				out << ", active";
			} else {
				out << ", inactive";
			}
		}
	}

	return out << std::endl;
}

/* Assemblaggio jacobiano */
VariableSubMatrixHandler&
TotalReaction::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	/*
	 * Constraint Equations:
	 * Position: 	R1^T(x2 + R2*f2 -x1 - R1*f1) - d = x^delta
	 * 		==> d(Vec3) = imposed displacement in node 1 local R.F.
	 * 		==> x^delta is used to activate/deactivate the constraint
	 * 			equation along the corresponding direction. If each
	 * 			component is set to 0, all relative displacement
	 * 			are forbidden (or imposed by the drive). If a component
	 * 			of x^delta is left free, the corrsponding equation is
	 * 			dropped
	 *
	 * Orientation:	Theta - Theta0 = ax(exp^-1(R1^T * R2 * R0^T)) = Theta^delta
	 * 		==> Theta = ax(exp^-1(R1^T * R2)) = Relative orientation in node1 R.F.
	 * 		==> Theta0 = Imposed relative orientation = ax(exp^-1(R0))
	 * 		==> Theta^delta is used to activate/deactivate the constraint
	 * 			equation along the corresponding direction. If each
	 * 			component is set to 0, all relative rotation
	 * 			are forbidden (or imposed by the drive). If a component
	 * 			of Theta^delta is left free, the corrsponding equation is
	 * 			dropped
	 *Jacobian Matrix:
	 *       x1  	     g1       	   x2    	g2       	 F	      M
	 * Q1 |  0   	     F1X           0            0              -R1            0	 | | x1 |
	 * G1 |-(F1)X  (b1)X(F1)X+(M1)X  (F1)X     -(F1)X(b2)X       (b1)X(R1)      -R1r | | g1 |
	 * Q2 |  0          -F1X           0    	0   	         R1	      0  | | x2 |
	 * G2 |  0    -(b2)X(F1)X-(M1)X    0        (F1)X(b2)X       (b2)X(R1)       R1r | | g2 |
	 * F  |-c*R1^T  c*R1^T*[(b1)X]   c*R1^T   -c*R1^T*[(b2)X]        0            0	 | | F  |if(bPos)
	 * M  |  0        -c*R1r^T         0         c* R1r^T            0   	      0	 | | M  |if(bRot)
	 *           	                                               if(bPos)    if(bRot)
	 *with: _ b1 = (x2 + R2*f2 - x1)
	 *      _ b2 = (R2*f2)
	 *      _ R1 = R1*R1h
	 *      _ R2 = R2*R2h
	 *      _ R1r = R1*R1hr
	 *      _ F1 = R1*F
	 *      _ M1 = R1*M
	 *      _ X = "Cross" operator
	 *
	 *     */

	DEBUGCOUT("Entering TotalReaction::AssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Ridimensiona la sottomatrice in base alle esigenze */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici delle varie incognite */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
	integer iFirstReactionIndex = total_equation_element->iGetFirstIndex();

	/* Setta gli indici delle equazioni */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	for (unsigned int iCnt = 1; iCnt <= nConstraints; iCnt++) {
		WM.PutColIndex(12 + iCnt, iFirstReactionIndex + iCnt);
	}

	/* Recupera i dati che servono */
	Mat3x3 R1(pNode1->GetRRef()*R1h);
	Mat3x3 R1r(pNode1->GetRRef()*R1hr);
	Vec3 b2(pNode2->GetRCurr()*f2);
	Vec3 b1(pNode2->GetXCurr() + b2 - pNode1->GetXCurr());

 	/* Moltiplica il momento e la forza per il coefficiente del metodo */
 	Vec3 FTmp(R1*(F*dCoef));
 	Vec3 MTmp(R1r*(M*dCoef));
 
 	/* Equilibrium: ((Phi/q)^T*Lambda)/q */
 
 	Mat3x3 Tmp;
 
 	/* [ F x ] */
 	Tmp = Mat3x3(MatCross, FTmp);
 
 	/* Lines 1->3: */
 	WM.Add(1, 3 + 1, Tmp);
 
 	/* Lines 4->6: */
 	WM.Sub(3 + 1, 1, Tmp);
 
 	WM.Add(3 + 1, 6 + 1, Tmp);
 
 	/* Lines 7->9: */
 	WM.Sub(6 + 1, 3 + 1, Tmp);
 
 	/* [ F x ] [ b2 x ] */
 	Tmp = Mat3x3(MatCrossCross, FTmp, b2);
 
 	/* Lines 4->6: */
 	WM.Sub(3 + 1, 9 + 1, Tmp);
 
 	/* Lines 10->12: */
 	WM.Add(9 + 1, 9 + 1, Tmp);
 
 	/* [ b1 x ] [ F x ] + [ M x ] */
 
 	/* Lines 4->6: */
 	WM.Add(3 + 1, 3 + 1, Mat3x3(MatCrossCross, b1, FTmp) + Mat3x3(MatCross, MTmp));
 
 	/* [ b2 x ] [ F x ] + [ M x ] */
 
 	/* Lines 10->12: */
 	WM.Sub(9 + 1, 3 + 1, Mat3x3(MatCrossCross, b2, FTmp) + Mat3x3(MatCross, MTmp));

/* Phi/q and (Phi/q)^T */

	Mat3x3 b1Cross_R1(b1.Cross(R1)); // = [ b1 x ] * R1
	Mat3x3 b2Cross_R1(b2.Cross(R1)); // = [ b2 x ] * R1

	for (unsigned iCnt = 0 ; iCnt < nPosConstraints; iCnt++) {
		Vec3 vR1(R1.GetVec(iPosIncid[iCnt]));
		Vec3 vb1Cross_R1(b1Cross_R1.GetVec(iPosIncid[iCnt]));
		Vec3 vb2Cross_R1(b2Cross_R1.GetVec(iPosIncid[iCnt]));

		/* Equilibrium, node 1 */
      		WM.Sub(1, 12 + 1 + iCnt, vR1);
      		WM.Sub(3 + 1, 12 + 1 + iCnt, vb1Cross_R1);

		/* Equilibrium, node 2 */
      		WM.Add(6 + 1, 12 + 1 + iCnt, vR1);
      		WM.Add(9 + 1, 12 + 1 + iCnt, vb2Cross_R1);
	}

	for (unsigned iCnt = 0 ; iCnt < nRotConstraints; iCnt++) {
		Vec3 vR1(R1r.GetVec(iRotIncid[iCnt]));

		/* Equilibrium, node 1 */
      		WM.Sub(3 + 1, 12 + 1 + nPosConstraints +  iCnt, vR1);

		/* Equilibrium, node 2 */
      		WM.Add(9 + 1, 12 + 1 + nPosConstraints + iCnt, vR1);
	}

	return WorkMat;
}

/* Assemblaggio residuo */
SubVectorHandler&
TotalReaction::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering TotalReaction::AssRes()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici */
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
	integer iFirstReactionIndex = total_equation_element->iGetFirstIndex();

	/* Indici dei nodi */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WorkVec.PutRowIndex(6+iCnt, iNode2FirstMomIndex + iCnt);
	}

	/* Get constraint reactions */

	for (unsigned iCnt = 0; iCnt < nPosConstraints; iCnt++) {
		F(iPosIncid[iCnt]) = XCurr(iFirstReactionIndex + 1 + iCnt);
	}

	for (unsigned iCnt = 0; iCnt < nRotConstraints; iCnt++) {
		M(iRotIncid[iCnt]) = XCurr(iFirstReactionIndex + 1 + nPosConstraints + iCnt);
	}


	Vec3 b2(pNode2->GetRCurr()*f2);
	Vec3 b1(pNode2->GetXCurr() + b2 - pNode1->GetXCurr());

	Mat3x3 R1 = pNode1->GetRCurr()*R1h;
	Mat3x3 R1r = pNode1->GetRCurr()*R1hr;
	Mat3x3 R2r = pNode2->GetRCurr()*R2hr;

	Vec3 FTmp(R1*F);
	Vec3 MTmp(R1r*M);

	/* Equilibrium, node 1 */
	WorkVec.Add(1, FTmp);
	WorkVec.Add(3 + 1, MTmp + b1.Cross(FTmp));

	/* Equilibrium, node 2 */
	WorkVec.Sub(6 + 1, FTmp);
	WorkVec.Sub(9 + 1, MTmp + b2.Cross(FTmp));

	return WorkVec;
}

DofOrder::Order
TotalReaction::GetEqType(unsigned int i) const
{
	ASSERTMSGBREAK(i < iGetNumDof(),
		"INDEX ERROR in TotalReaction::GetEqType");

	return DofOrder::ALGEBRAIC;
}

/* Output (da mettere a punto), per ora solo reazioni */
void
TotalReaction::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		const Vec3& X1(pNode1->GetXCurr());
		const Vec3& X2(pNode2->GetXCurr());
		const Mat3x3& R1(pNode1->GetRCurr());
		const Mat3x3& R2(pNode2->GetRCurr());
		const Vec3& V1(pNode1->GetVCurr());
		const Vec3& V2(pNode2->GetVCurr());
		const Vec3& Omega1(pNode1->GetWCurr());
		const Vec3& Omega2(pNode2->GetWCurr());

		Mat3x3 R1Tmp(R1*R1h);
		Mat3x3 R1rTmp(R1*R1hr);
		Mat3x3 R2rTmp(R2*R2hr);
		Mat3x3 RTmp(R1rTmp.MulTM(R2rTmp));
		Vec3 b2(R2*f2);
		Vec3 b1(X2 + b2 - X1);
		Vec3 X(R1Tmp.MulTV(b1) - tilde_f1);

		Joint::Output(OH.Joints(), "TotalReaction", GetLabel(),
			F, M, R1Tmp*F, R1rTmp*M)
			<< " " << X
			<< " " << RotManip::VecRot(RTmp)
			<< " " << R1Tmp.MulTV(V2 + Omega2.Cross(b2) - V1 - Omega1.Cross(b1))
			<< " " << R1rTmp.MulTV(Omega2 - Omega1)
			<< std::endl;
			// accelerations?
	}
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
TotalReaction::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{

	/* Per ora usa la matrice piena; eventualmente si puo'
	 * passare a quella sparsa quando si ottimizza */
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Equazioni: vedi joints.dvi */

	/*	 equazioni ed incognite
	 * F1					  Delta_x1	1	   
	 * M1					  Delta_g1	3 + 1 
	 * FP1  				  Delta_xP1	6 + 1
	 * MP1  				  Delta_w1	9 + 1 
	 * F2					  Delta_x2	12 + 1
	 * M2					  Delta_g2	15 + 1
	 * FP2  				  Delta_xP2	18 + 1  
	 * MP2  				  Delta_w2	21 + 1 
	 * vincolo spostamento  		  Delta_F	24 + 1 
	 * vincolo rotazione			  Delta_M	24 + nPosConstraints  
	 * derivata vincolo spostamento 	  Delta_FP	24 + nConstraints  
	 * derivata vincolo rotazione		  Delta_MP	24 + nConstraints + nPosConstraints  
	 */
	
	
	/* Indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstVelIndex = iNode1FirstPosIndex + 6;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstVelIndex = iNode2FirstPosIndex + 6;
	integer iFirstReactionIndex = total_equation_element->iGetFirstIndex();
	integer iReactionPrimeIndex = iFirstReactionIndex + nConstraints;


	/* Setta gli indici dei nodi */
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

	/* Setta gli indici delle reazioni */
	for (unsigned int iCnt = 1; iCnt <=  nConstraints; iCnt++) {
		WM.PutRowIndex(24 + iCnt, iFirstReactionIndex + iCnt);
		WM.PutColIndex(24 + iCnt, iFirstReactionIndex + iCnt);
		WM.PutRowIndex(24 + nConstraints + iCnt, iReactionPrimeIndex + iCnt);
		WM.PutColIndex(24 + nConstraints + iCnt, iReactionPrimeIndex + iCnt);
	}
	
	/* Recupera i dati che servono */
	Mat3x3 R1(pNode1->GetRRef()*R1h);
	Mat3x3 R1r(pNode1->GetRRef()*R1hr);
	
	Vec3 b2(pNode2->GetRCurr()*f2);
	Vec3 b1(pNode2->GetXCurr() + b2 - pNode1->GetXCurr());
	
	Vec3 Omega1(pNode1->GetWCurr());
	Vec3 Omega2(pNode2->GetWCurr());
	
	Vec3 b1Prime(pNode2->GetVCurr() + Omega2.Cross(b2) - pNode1->GetVCurr());
	
	/* F ed M sono gia' state aggiornate da InitialAssRes;
	 * Recupero FPrime e MPrime*/
	Vec3 MPrime(Zero3);
	Vec3 FPrime(Zero3);

	for (unsigned iCnt = 0; iCnt < nPosConstraints; iCnt++) {
		FPrime(iPosIncid[iCnt]) = XCurr(iReactionPrimeIndex + 1 + iCnt);
	}
	
	for (unsigned iCnt = 0; iCnt < nRotConstraints; iCnt++) {
		MPrime(iRotIncid[iCnt]) = XCurr(iReactionPrimeIndex + 1 + nPosConstraints + iCnt);
	}

	Vec3 FTmp(R1 * F);
	Vec3 MTmp(R1r * M);
	Vec3 FPrimeTmp(R1 * FPrime);
	Vec3 MPrimeTmp(R1r * MPrime);

	Mat3x3 Tmp;

	/* [ F x ] */
	Tmp = Mat3x3(MatCross, FTmp);

	/* Force, Node 1, Lines 1->3: */
	WM.Add(1, 3 + 1, Tmp);	// * Delta_g1

	/* Moment, Node 1, Lines 4->6: */
	WM.Sub(3 + 1, 1, Tmp);	// * Delta_x1

	WM.Add(3 + 1, 12 + 1, Tmp);	// * Delta_x2

	/* Force, Node 2, Lines 13->15: */
	WM.Sub(12 + 1, 3 + 1, Tmp);	// * Delta_g1
	
	/* [ FP x ] */
	Tmp = Mat3x3(MatCross, FPrimeTmp);

	/* d/dt(Force), Node 1, Lines 7->9: */
	WM.Add(6 + 1, 9 + 1, Tmp);	// * Delta_W1

	/* d/dt(Moment), Node 1, Lines 10->12: */
	WM.Sub(9 + 1, 6 + 1, Tmp);	// * Delta_x1P

	WM.Add(9 + 1, 18 + 1, Tmp);	// * Delta_x2P

	/* d/dt(Force), Node 2, Lines 19->21: */
	WM.Sub(18 + 1 , 9 + 1, Tmp);	// * Delta_W1
	
	/* [ F x ] [ b2 x ] */
	Tmp = Mat3x3(MatCrossCross, FTmp, b2);

	/* Moment, Node1, Lines 4->6: */
	WM.Sub(3 + 1, 15 + 1, Tmp);	// * Delta_g2

	/* Moment, Node2, Lines 16->18: */
	WM.Add(15 + 1, 15 + 1, Tmp);	// * Delta_g2
	
	/* [ FP x ] [ b2 x ] */
	Tmp = Mat3x3(MatCrossCross, FPrimeTmp, b2);

	/* d/dt(Moment), Node1, Lines 10->12: */
	WM.Sub(9 + 1, 21 + 1, Tmp);	// * Delta_W2

	/* d/dt(Moment), Node2, Lines 22->24: */
	WM.Add(21 + 1, 21 + 1, Tmp);	// * Delta_W2

	/* [ b1 x ] [ F x ] + [ M x ] */

	/* Moment, Node1, Lines 4->6: */
	WM.Add(3 + 1, 3 + 1, Mat3x3(MatCrossCross, b1, FTmp) + Mat3x3(MatCross, MTmp));	// *Delta_g1

	/* d/dt(Moment), Node1, Lines 10->12: */
	WM.Add(9 + 1, 9 + 1, Mat3x3(MatCrossCross, b1, FPrimeTmp) + Mat3x3(MatCross, MPrimeTmp));// *Delta_W1
	
	/* [ b2 x ] [ F x ] + [ M x ] */

	/* Moment, Node2, Lines 16->18: */
	WM.Sub(15 + 1, 3 + 1, Mat3x3(MatCrossCross, b2, FTmp) + Mat3x3(MatCross, MTmp));	// * Delta_g1
	
	/* d/dt(Moment), Node2, Lines 22->24: */
	WM.Sub(21 + 1, 9 + 1, Mat3x3(MatCrossCross, b2, FTmp) + Mat3x3(MatCross, MTmp));	// * Delta_W1
	
	/* Constraints: Add only active rows/columns*/	
	
	/* Positions contribution:*/

	Mat3x3 b1Cross_R1(b1.Cross(R1)); // = [ b1 x ] * R1
	Mat3x3 b2Cross_R1(b2.Cross(R1)); // = [ b2 x ] * R1

	for (unsigned iCnt = 0 ; iCnt < nPosConstraints; iCnt++) {
		Vec3 vR1(R1.GetVec(iPosIncid[iCnt]));
		Vec3 vb1Cross_R1(b1Cross_R1.GetVec(iPosIncid[iCnt]));
		Vec3 vb2Cross_R1(b2Cross_R1.GetVec(iPosIncid[iCnt]));

		/* Equilibrium, node 1 */
      		WM.Sub(1, 24 + 1 + iCnt, vR1);	// * Delta_F
      		WM.Sub(3 + 1, 24 + 1 + iCnt, vb1Cross_R1);	// * Delta_F

		/* Constraint, node 1 */
      		WM.SubT(24 + 1 + iCnt, 1, vR1);	// * Delta_x1
      		WM.SubT(24 + 1 + iCnt, 3 + 1, vb1Cross_R1);	// * Delta_g1

		/* d/dt(Equilibrium), node 1 */
      		WM.Sub(6 + 1, 24 + 1 + nConstraints + iCnt, vR1);	// * Delta_FP
      		WM.Sub(9 + 1, 24 + 1 + nConstraints + iCnt, vb1Cross_R1);	// * Delta_FP

		/* d/dt(Constraint), node 1 */
      		WM.SubT(24 + 1 + nConstraints + iCnt, 6 + 1, vR1);	// * Delta_xP1
      		WM.SubT(24 + 1 + nConstraints + iCnt, 9 + 1, vb1Cross_R1);	// * Delta_W1
		
		/* Equilibrium, node 2 */
      		WM.Add(12 + 1, 24 + 1 + iCnt, vR1);	// * Delta_F
      		WM.Add(15 + 1, 24 + 1 + iCnt, vb2Cross_R1);	// * Delta_F

		/* Constraint, node 2 */
      		WM.AddT(24 + 1 + iCnt, 12 + 1, vR1);	// * Delta_x2
      		WM.AddT(24 + 1 + iCnt, 15 + 1, vb2Cross_R1);	// * Delta_g2
	
		/* d/dt(Equilibrium), node 2 */
      		WM.Add(18 + 1, 24 + 1 + nConstraints + iCnt, vR1);	// * Delta_FP
      		WM.Add(21 + 1, 24 + 1 + nConstraints + iCnt, vb2Cross_R1);	// * Delta_FP

		/* d/dt(Constraint), node 2 */
      		WM.AddT(24 + 1 + nConstraints +  iCnt, 18 + 1, vR1);	// * Delta_xP2
      		WM.AddT(24 + 1 + nConstraints +  iCnt, 21 + 1, vb2Cross_R1);	// * Delta_W2
	}

	for (unsigned iCnt = 0 ; iCnt < nRotConstraints; iCnt++) {
		Vec3 vR1(R1r.GetVec(iRotIncid[iCnt]));

		/* Equilibrium, node 1 */
      		WM.Sub(3 + 1, 24 + 1 + nPosConstraints +  iCnt, vR1);	// * Delta_M

		/* Constraint, node 1 */
      		WM.SubT(24 + 1 + nPosConstraints + iCnt, 3 + 1, vR1);	// * Delta_g1

		/* d/dt(Equilibrium), node 1 */
      		WM.Sub(9 + 1, 24 + 1 + nConstraints + nPosConstraints +  iCnt, vR1);	// * Delta_MP

		/* d/dt(Constraint), node 1 */
      		WM.SubT(24 + 1 + nConstraints + nPosConstraints + iCnt, 9 + 1, vR1);	// * Delta_W1

		/* Equilibrium, node 2 */
      		WM.Add(15 + 1, 24 + 1 + nPosConstraints + iCnt, vR1);	// * Delta_M

		/* Constraint, node 2 */
      		WM.AddT(24 + 1 + nPosConstraints +  iCnt, 15 + 1, vR1);	// * Delta_g2

		/* d/dt(Equilibrium), node 2 */
      		WM.Add(21 + 1, 24 + 1 + nConstraints + nPosConstraints + iCnt, vR1);	// * Delta_MP

		/* d/dt(Constraint), node 2 */
      		WM.AddT(24 + 1 + nConstraints + nPosConstraints + iCnt, 21 + 1, vR1);	// * Delta_W2
	}

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
TotalReaction::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{

	DEBUGCOUT("Entering TotalReaction::InitialAssRes()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstVelIndex = iNode1FirstPosIndex + 6;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstVelIndex = iNode2FirstPosIndex + 6;
	integer iFirstReactionIndex = total_equation_element->iGetFirstIndex();
	integer iReactionPrimeIndex = iFirstReactionIndex + nConstraints;

	/* Setta gli indici */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode1FirstVelIndex + iCnt);
		WorkVec.PutRowIndex(12 + iCnt, iNode2FirstPosIndex + iCnt);
		WorkVec.PutRowIndex(18 + iCnt, iNode2FirstVelIndex + iCnt);
	}

	for (unsigned int iCnt = 1; iCnt <= nConstraints; iCnt++)	{
		WorkVec.PutRowIndex(24 + iCnt, iFirstReactionIndex + iCnt);
		WorkVec.PutRowIndex(24 + nConstraints + iCnt, iReactionPrimeIndex + iCnt);
	}

	/* Recupera i dati */
	/* Recupera i dati che servono */
	Mat3x3 R1(pNode1->GetRRef()*R1h);
	Mat3x3 R1r(pNode1->GetRRef()*R1hr);
	Mat3x3 R2r(pNode2->GetRCurr()*R2hr);
	
	Vec3 b2(pNode2->GetRCurr()*f2);
	Vec3 b1(pNode2->GetXCurr() + b2 - pNode1->GetXCurr());
	
	Vec3 Omega1(pNode1->GetWCurr());
	Vec3 Omega2(pNode2->GetWCurr());
	
	Vec3 b1Prime(pNode2->GetVCurr() + Omega2.Cross(b2) - pNode1->GetVCurr());
	
	Vec3 FPrime(Zero3);
	Vec3 MPrime(Zero3);

	/* Aggiorna F ed M, che restano anche per InitialAssJac */
	for (unsigned iCnt = 0; iCnt < nPosConstraints; iCnt++) {
		F(iPosIncid[iCnt]) = XCurr(iFirstReactionIndex + 1 + iCnt);
		FPrime(iPosIncid[iCnt]) = XCurr(iReactionPrimeIndex + 1 + iCnt);
	}

	for (unsigned iCnt = 0; iCnt < nRotConstraints; iCnt++) {
		M(iRotIncid[iCnt]) = XCurr(iFirstReactionIndex + 1 + nPosConstraints + iCnt);
		MPrime(iRotIncid[iCnt]) = XCurr(iReactionPrimeIndex + 1 + nPosConstraints + iCnt);
	}

	Vec3 FTmp(R1 * F);
	Vec3 MTmp(R1r * M);
	Vec3 FPrimeTmp(R1 * FPrime);
	Vec3 MPrimeTmp(R1r * MPrime);

	/* Equilibrium, node 1 */
	WorkVec.Add(1, FTmp);
	WorkVec.Add(3 + 1, b1.Cross(FTmp) + MTmp);

	/* d/dt(Equilibrium), node 1 */
	WorkVec.Add(6 + 1, FPrimeTmp);
	WorkVec.Add(9 + 1, b1.Cross(FPrimeTmp) + MPrimeTmp);

	/* Equilibrium, node 2 */
	WorkVec.Sub(12 + 1, FTmp);
	WorkVec.Sub(15 + 1, b2.Cross(FTmp) + MTmp);

	/* d/dt( Equilibrium ) , node 2 */
	WorkVec.Sub(18 + 1, FPrimeTmp);
	WorkVec.Sub(21 + 1, b2.Cross(FPrimeTmp) + MPrimeTmp);

	return WorkVec;
}


unsigned int
TotalReaction::iGetNumPrivData(void) const
{
	return 21;
}

unsigned int
TotalReaction::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);

	unsigned int off = 0;

	switch (s[0]) {
	case 'p':
		/* relative position */
		break;

	case 'r':
		/* relative orientation */
		off += 3;
		break;

	case 'F':
		/* force */
		off += 6;
		break;

	case 'M':
		/* moment */
		off += 9;
		break;

	case 'd':
		/* imposed relative position */
		off += 12;
		break;

	case 't':
		/* imposed relative orientation */
		off += 15;
		break;
	
	case 'v':
		/* relative linear velocity */
		off += 18;
		break;

	default:
		return 0;
	}

	switch (s[1]) {
	case 'x':
		return off + 1;

	case 'y':
		return off + 2;

	case 'z':
		return off + 3;
	}

	return 0;
}

doublereal
TotalReaction::dGetPrivData(unsigned int i) const
{
	switch (i) {
	case 1:
	case 2:
	case 3: {
		Vec3 x(pNode1->GetRCurr().MulTV(pNode2->GetXCurr() + pNode2->GetRCurr()*f2 - pNode1->GetXCurr()) - f1);
			return R1h.GetVec(i)*x;
		}

	case 4:
	case 5:
	case 6: {
#if 0
 		Vec3 Theta(Unwrap(ThetaDeltaPrev, ThetaDelta) + ThetaDrv.Get());
 		return Theta(i - 3);
#endif
		return 0.;
		}

	case 7:
	case 8:
	case 9:
		return F(i - 6);

	case 10:
	case 11:
	case 12:
		return M(i - 9);

	case 13:
	case 14:
	case 15:
		if (!bPosActive[i - 13]) {
			return 0.;
		}
// 		return XDrv.Get()(i - 12);
		return 0.;

	case 16:
	case 17:
	case 18:
		if (!bRotActive[i - 16]) {
			return 0.;
		}
// 		return ThetaDrv.Get()(i - 15);
		return 0.;
	case 19:
	case 20:
	case 21:
		{
		// FIXME: blind fix
		Vec3 v(pNode1->GetRCurr().MulTV(
			(pNode2->GetVCurr() + pNode2->GetWCurr().Cross(pNode2->GetRCurr()*f2)
				- pNode1->GetVCurr()) -
			pNode1->GetWCurr().Cross( pNode2->GetXCurr() + 
				pNode2->GetRCurr()*f2
				- pNode1->GetXCurr() -f1) 
			)
		);
		
			return R1h.GetVec(i-18)*v;
		}


	default:
		ASSERT(0);
	}

	return 0.;
}

/* TotalReaction - end */
