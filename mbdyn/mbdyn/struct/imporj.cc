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

/* Continuano i vincoli di rotazione piani */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <fstream>

#include "imporj.h"
#include "Rot.hh"
#include "hint_impl.h"

ImposedOrientationJoint::ImposedOrientationJoint(unsigned int uL, const DofOwner* pDO,
	bool b[3],
	const TplDriveCaller<Vec3>* pDC,
	const StructNode* pN1, const StructNode* pN2,
	const Mat3x3& R1hTmp, const Mat3x3& R2hTmp,
	flag fOut)
: Elem(uL, fOut),
Joint(uL, pDO, fOut),
TplDriveOwner<Vec3>(pDC),
pNode1(pN1), pNode2(pN2),
R1h(R1hTmp), R2h(R2hTmp),
nConstraints(0),
M(0.), Theta(0.)
{
	for (unsigned int i = 0; i < 3; i++) {
		bActive[i] = b[i];
		if (bActive[i]) {
			nConstraints++;
		}
	}
}

ImposedOrientationJoint::~ImposedOrientationJoint(void)
{
	NO_OP;
};

static const char xyz[] = "xyz";

std::ostream&
ImposedOrientationJoint::DescribeDof(std::ostream& out,
	const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

	out << prefix << iIndex + 1;
	if (nConstraints > 1) {
		out << "->" << iIndex + nConstraints;
	}
	out << ": ";

	out << "reaction couple(s) [";
	for (unsigned int i = 0, cnt = 0; i < 3; i++) {
		if (bActive[i]) {
			cnt++;
			if (cnt > 1) {
				out << ",";
			}
			out << "m" << xyz[i];
		}
	}
	out << "]" << std::endl;

	if (bInitial) {
		iIndex += nConstraints;

		out << prefix << iIndex + 1;
		if (nConstraints > 1) {
			out << "->" << iIndex + nConstraints;
		}
		out << ": ";

		out << "reaction couple derivative(s) [";
	
		for (unsigned int i = 0, cnt = 0; i < 3; i++) {
			if (bActive[i]) {
				cnt++;
				if (cnt > 1) {
					out << ",";
				}
				out << "mP" << xyz[i];
			}
		}
		out << "]" << std::endl;
	}

	return out;
}

static const char *dof[] = {
	"reaction couple m",
	"reaction couple derivative mP"
};
static const char *eq[] = {
	"orientation constraint g",
	"orientation constraint derivative w"
};

void
ImposedOrientationJoint::DescribeDof(std::vector<std::string>& desc,
	bool bInitial, int i) const
{
	int nself = 1;
	
	if (i == -1) {
		nself = nConstraints;
		if (bInitial) {
			nself *= 2;
		}
	}
	desc.resize(nself);

	unsigned iend = 3;
	if (bInitial) {
		iend *= 2;
	}

	std::ostringstream os;
	os << "ImposedOrientationJoint(" << GetLabel() << ")";

	if (i == -1) {
		std::string name(os.str());

		int n = 0;
		for (unsigned i = 0; i < iend; i++) {
			if (bActive[i%3]) {
				os.str(name);
				os.seekp(0, std::ios_base::end);
				os << ": " << dof[i/3] << xyz[i%3];
				desc[n] = os.str();
				if (++n == nself) {
					break;
				}
			}
		}

	} else {
		int n = 0;
		for (unsigned j = 0; j < iend; j++) {
			if (bActive[j%3]) {
				if (n == i) {
					os << ": " << dof[i/3] << xyz[i%3];
					desc[0] = os.str();
					break;
				}
			}
		}
	}
}

std::ostream&
ImposedOrientationJoint::DescribeEq(std::ostream& out,
	const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

	out << prefix << iIndex + 1;
	if (nConstraints > 1) {
		out << "->" << iIndex + nConstraints;
	}
	out << ": ";

	out << "orientation constraint(s) [";
	
	for (unsigned int i = 0, cnt = 0; i < 3; i++) {
		if (bActive[i]) {
			cnt++;
			if (cnt > 1) {
				out << ",";
			}
			out << "g" << xyz[i] << "1=g" << xyz[i] << "2";
		}
	}
	out << "]" << std::endl;

	if (bInitial) {
		iIndex += nConstraints;

		out << prefix << iIndex + 1;
		if (nConstraints > 1) {
			out << "->" << iIndex + nConstraints;
		}
		out << ": ";

		out << "angular velocity constraint(s) [";
	
		for (unsigned int i = 0, cnt = 0; i < 3; i++) {
			if (bActive[i]) {
				cnt++;
				if (cnt > 1) {
					out << ",";
				}
				out << "w" << xyz[i] << "1=w" << xyz[i] << "2";
			}
		}
		out << "]" << std::endl;
	}

	return out;
}

void
ImposedOrientationJoint::DescribeEq(std::vector<std::string>& desc,
	bool bInitial, int i) const
{
	int nself = 1;
	
	if (i == -1) {
		nself = nConstraints;
		if (bInitial) {
			nself *= 2;
		}
	}
	desc.resize(nself);

	unsigned iend = 3;
	if (bInitial) {
		iend *= 2;
	}

	std::ostringstream os;
	os << "ImposedOrientationJoint(" << GetLabel() << ")";

	if (i == -1) {
		std::string name(os.str());

		int n = 0;
		for (unsigned i = 0; i < iend; i++) {
			if (bActive[i%3]) {
				os.str(name);
				os.seekp(0, std::ios_base::end);
				os << ": " << eq[i/3] << xyz[i%3];
				desc[n] = os.str();
				if (++n == nself) {
					break;
				}
			}
		}

	} else {
		int n = 0;
		for (unsigned j = 0; j < iend; j++) {
			if (bActive[j%3]) {
				if (n == i) {
					os << ": " << eq[i/3] << xyz[i%3];
					desc[0] = os.str();
					break;
				}
			}
		}
	}
}

void
ImposedOrientationJoint::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	if (ph) {
		for (unsigned int i = 0; i < ph->size(); i++) {
			Joint::JointHint *pjh = dynamic_cast<Joint::JointHint *>((*ph)[i]);

			if (pjh == 0) {
				continue;
			}

			if (dynamic_cast<Joint::HingeHint<1> *>(pjh)) {
				R1h = pNode1->GetRCurr().Transpose()*pNode2->GetRCurr()*R2h;

			} else if (dynamic_cast<Joint::HingeHint<2> *>(pjh)) {
				R2h = pNode2->GetRCurr().Transpose()*pNode1->GetRCurr()*R1h;

			} else if (dynamic_cast<Joint::ReactionsHint *>(pjh)) {
				/* TODO */
			}
		}
	}
}

Hint *
ImposedOrientationJoint::ParseHint(DataManager *pDM, const char *s) const
{
	if (strncasecmp(s, "hinge{" /*}*/, STRLENOF("hinge{" /*}*/)) == 0) {
		s += STRLENOF("hinge{" /*}*/);

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

	return 0;
}

void
ImposedOrientationJoint::AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP)
{
#if 0
	Mat3x3 RTmp(((pNode1->GetRCurr()*R1h).Transpose()
			*pNode1->GetRPrev()*R1h).Transpose()
			*((pNode2->GetRCurr()*R2h).Transpose()
			*pNode2->GetRPrev()*R2h));
	Vec3 v(MatR2EulerAngles(RTmp.Transpose()));

	dTheta += v.dGet(3);
#endif
}

/* Contributo al file di restart */
std::ostream&
ImposedOrientationJoint::Restart(std::ostream& out) const
{
#if 0
	Joint::Restart(out) << ", imposed orientation, "
		<< pNode1->GetLabel() << ", hinge, reference, node, "
			"1, ", (R1h.GetVec(1)).Write(out, ", ") << ", "
			"2, ", (R1h.GetVec(2)).Write(out, ", ") << ", "
		<< pNode2->GetLabel() << ", hinge, reference, node, "
			"1, ", (R2h.GetVec(1)).Write(out, ", ") << ", "
			"2, ", (R2h.GetVec(2)).Write(out, ", ") << ";"
		<< std::endl;
#endif

	return out;
}

/* Assemblaggio jacobiano */
VariableSubMatrixHandler&
ImposedOrientationJoint::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering ImposedOrientationJoint::AssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Ridimensiona la sottomatrice in base alle esigenze */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici delle varie incognite */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex() + 3;
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex() + 3;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex() + 3;
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex() + 3;
	integer iFirstReactionIndex = iGetFirstIndex();

	/* Setta gli indici delle equazioni */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(3 + iCnt, iNode2FirstMomIndex + iCnt);
		WM.PutColIndex(3 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	for (unsigned int iCnt = 0; iCnt <= nConstraints; iCnt++) {
		WM.PutRowIndex(6 + iCnt, iFirstReactionIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iFirstReactionIndex + iCnt);
	}

	/* Recupera i dati che servono */
	Mat3x3 R1(pNode1->GetRRef()*R1h);

	/* Moltiplica il momento per il coefficiente del metodo */
	Mat3x3 MTmp(M*dCoef);

	WM.Add(1, 1, MTmp);
	WM.Sub(3 + 1, 1, MTmp);

	for (int iCnt = 0, iCurr = 0; iCnt < 3; iCnt++) {
		if (bActive[iCnt]) {
			iCurr++;
			for (int iIdx = 1; iIdx <= 3; iIdx++) {
				doublereal d = R1(iIdx, iCnt + 1);

				/* Equilibrium, node 1 */
      				WM.PutCoef(iIdx, 6 + iCurr, -d);

				/* Equilibrium, node 2 */
      				WM.PutCoef(3 + iIdx, 6 + iCurr, d);

				/* Constraint, node 1 */
      				WM.PutCoef(6 + iCurr, iIdx, -d);

				/* Constraint, node 2 */
      				WM.PutCoef(6 + iCurr, 3 + iIdx, d);
			}
		}
	}

	return WorkMat;
}

/* Assemblaggio residuo */
SubVectorHandler&
ImposedOrientationJoint::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering ImposedOrientationJoint::AssRes()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici */
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex() + 3;
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex() + 3;
	integer iFirstReactionIndex = iGetFirstIndex();

	/* Indici dei nodi */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WorkVec.PutRowIndex(3 + iCnt, iNode2FirstMomIndex + iCnt);
	}

	/* Indici del vincolo */
	for (unsigned int iCnt = 1; iCnt <= nConstraints; iCnt++) {
		WorkVec.PutRowIndex(6 + iCnt, iFirstReactionIndex + iCnt);
	}

	/* Get constraint reactions */
	for (int iCnt = 0, iCurr = 0; iCnt < 3; iCnt++) {
		if (bActive[iCnt]) {
			iCurr++;
			M(iCnt + 1) = XCurr(iFirstReactionIndex + iCurr);

		} else {
			M(iCnt + 1) = 0.;
		}
	}

	Mat3x3 R1 = pNode1->GetRCurr()*R1h;
	Mat3x3 R2 = pNode2->GetRCurr()*R2h;
	Vec3 Theta0 = Get();
	Mat3x3 R0T = RotManip::Rot(-Theta0);	// -Theta0 to get R0 transposed
	Mat3x3 Rdelta = R1.Transpose()*R2*R0T;
	Vec3 Thetadelta = RotManip::VecRot(Rdelta);

	Vec3 MTmp = R1*M;

	/* Equilibrium, node 1 */
	WorkVec.Add(1, MTmp);

	/* Equilibrium, node 2 */
	WorkVec.Sub(3 + 1, MTmp);

	/* Constraint equations are divided by dCoef */
	if (dCoef != 0.) {

		/* Equazioni di vincolo di rotazione */
		for (int iCnt = 0, iCurr = 0; iCnt < 3; iCnt++) {
			if (bActive[iCnt]) {
				iCurr++;
				WorkVec.PutCoef(6 + iCurr, -Thetadelta(iCnt + 1)/dCoef);
			}
		}
	}

	return WorkVec;
}

DofOrder::Order
ImposedOrientationJoint::GetEqType(unsigned int i) const
{
#if 0
	ASSERTMSGBREAK(i < iGetNumDof(),
		"INDEX ERROR in ImposedOrientationJoint::GetEqType");
#endif
	return DofOrder::ALGEBRAIC;
}

/* Output (da mettere a punto) */
void
ImposedOrientationJoint::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		Mat3x3 R2Tmp(pNode2->GetRCurr()*R2h);
		Mat3x3 RTmp((pNode1->GetRCurr()*R1h).Transpose()*R2Tmp);
		Mat3x3 R2TmpT(R2Tmp.Transpose());

		Joint::Output(OH.Joints(), "ImposedOrientation", GetLabel(),
			Zero3, M, Zero3, R2Tmp*M)
			<< " " << MatR2EulerAngles(RTmp)*dRaDegr
			<< " " << R2TmpT*(pNode2->GetWCurr()-pNode1->GetWCurr()) << std::endl;
	}
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
ImposedOrientationJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
#if 0
   /* Per ora usa la matrice piena; eventualmente si puo'
    * passare a quella sparsa quando si ottimizza */
   FullSubMatrixHandler& WM = WorkMat.SetFull();

   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeReset(iNumRows, iNumCols);

   /* Equazioni: vedi joints.dvi */

   /*       equazioni ed incognite
    * F1                                     Delta_x1         0+1 =  1
    * M1                                     Delta_g1         3+1 =  4
    * FP1                                    Delta_xP1        6+1 =  7
    * MP1                                    Delta_w1         9+1 = 10
    * F2                                     Delta_x2        12+1 = 13
    * M2                                     Delta_g2        15+1 = 16
    * FP2                                    Delta_xP2       18+1 = 19
    * MP2                                    Delta_w2        21+1 = 22
    * vincolo spostamento                    Delta_F         24+1 = 25
    * vincolo rotazione                      Delta_M         27+1 = 28
    * derivata vincolo spostamento           Delta_FP        29+1 = 30
    * derivata vincolo rotazione             Delta_MP        32+1 = 33
    */


   /* Indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex()+3;
   integer iNode1FirstVelIndex = iNode1FirstPosIndex+6+3;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex()+3;
   integer iNode2FirstVelIndex = iNode2FirstPosIndex+6+3;
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+2;

   /* Nota: le reazioni vincolari sono:
    * Forza,       3 incognite, riferimento globale,
    * Momento,     2 incognite, riferimento locale
    */

   /* Setta gli indici dei nodi */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutRowIndex(3+iCnt, iNode1FirstVelIndex+iCnt);
      WM.PutColIndex(3+iCnt, iNode1FirstVelIndex+iCnt);
      WM.PutRowIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
      WM.PutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
      WM.PutRowIndex(9+iCnt, iNode2FirstVelIndex+iCnt);
      WM.PutColIndex(9+iCnt, iNode2FirstVelIndex+iCnt);
   }

   /* Setta gli indici delle reazioni */
   for (int iCnt = 1; iCnt <= 4; iCnt++) {
      WM.PutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
      WM.PutColIndex(12+iCnt, iFirstReactionIndex+iCnt);
   }


   /* Recupera i dati */
   const Mat3x3& R1(pNode1->GetRRef());
   const Mat3x3& R2(pNode2->GetRRef());
   const Vec3& Omega1(pNode1->GetWRef());
   const Vec3& Omega2(pNode2->GetWRef());

   /* F ed M sono gia' state aggiornate da InitialAssRes */
   Vec3 MPrime(XCurr.dGetCoef(iReactionPrimeIndex+1),
	       XCurr.dGetCoef(iReactionPrimeIndex+2),
	       0.);

   /* Matrici identita' */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      /* Contributo di forza all'equazione della forza, nodo 1 */
      WM.PutCoef(iCnt, 12+iCnt, 1.);

      /* Contrib. di der. di forza all'eq. della der. della forza, nodo 1 */
      WM.PutCoef(3+iCnt, 14+iCnt, 1.);

      /* Contributo di forza all'equazione della forza, nodo 2 */
      WM.PutCoef(6+iCnt, 12+iCnt, -1.);

      /* Contrib. di der. di forza all'eq. della der. della forza, nodo 2 */
      WM.PutCoef(9+iCnt, 14+iCnt, -1.);
   }

   /* Matrici di rotazione dai nodi alla cerniera nel sistema globale */
   Mat3x3 R1hTmp(R1*R1h);
   Mat3x3 R2hTmp(R2*R2h);

   Vec3 e3a(R1hTmp.GetVec(3));
   Vec3 e1b(R2hTmp.GetVec(1));
   Vec3 e2b(R2hTmp.GetVec(2));

   /* Ruota il momento e la sua derivata con le matrici della cerniera
    * rispetto ai nodi */
   Vec3 MTmp(e2b*M.dGet(1)-e1b*M.dGet(2));
   Vec3 MPrimeTmp(e2b*MPrime.dGet(1)-e1b*MPrime.dGet(2));

   Mat3x3 MDeltag1((Mat3x3(Omega2.Cross(MTmp)+MPrimeTmp)+
		    Mat3x3(MTmp, Omega1))*Mat3x3(e3a));
   Mat3x3 MDeltag2(Mat3x3(Omega1.Cross(e3a), MTmp)+
		   Mat3x3(e3a, MPrimeTmp)+
		   Mat3x3(e3a)*Mat3x3(Omega2, MTmp));

   /* Vettori temporanei */
   Vec3 Tmp1(e2b.Cross(e3a));
   Vec3 Tmp2(e3a.Cross(e1b));

   /* Prodotto vettore tra il versore 3 della cerniera secondo il nodo 1
    * ed il versore 1 della cerniera secondo il nodo 2. A convergenza
    * devono essere ortogonali, quindi il loro prodotto vettore deve essere
    * unitario */

   /* Error handling: il programma si ferma, pero' segnala dov'e' l'errore */
   if (Tmp1.Dot() <= DBL_EPSILON || Tmp2.Dot() <= DBL_EPSILON) {
      silent_cerr("ImposedOrientationJoint(" << GetLabel() << "): "
	      "first and second node hinge axes are (nearly) orthogonal"
	      << std::endl);
      throw Joint::ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

   Vec3 TmpPrime1(e2b.Cross(Omega1.Cross(e3a))-e3a.Cross(Omega2.Cross(e2b)));
   Vec3 TmpPrime2(e3a.Cross(Omega2.Cross(e1b))-e1b.Cross(Omega1.Cross(e3a)));

   /* Equazione di momento, nodo 1 */
   WM.Sub(4, 4, Mat3x3(MTmp, e3a));
   WM.Add(4, 16, Mat3x3(e3a, MTmp));
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutCoef(iCnt, 13, Tmp1.dGet(iCnt));
      WM.PutCoef(iCnt, 14, Tmp2.dGet(iCnt));
   }

   /* Equazione di momento, nodo 2 */
   WM.Add(7, 1, Mat3x3(MTmp, e3a));
   WM.Sub(7, 7, Mat3x3(e3a, MTmp));
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutCoef(6+iCnt, 13, -Tmp1.dGet(iCnt));
      WM.PutCoef(6+iCnt, 14, -Tmp2.dGet(iCnt));
   }

   /* Derivata dell'equazione di momento, nodo 1 */
   WM.Sub(4, 1, MDeltag1);
   WM.Sub(4, 4, Mat3x3(MTmp, e3a));
   WM.Add(4, 7, MDeltag2);
   WM.Add(4, 10, Mat3x3(e3a, MTmp));

   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutCoef(3+iCnt, 13, TmpPrime1.dGet(iCnt));
      WM.PutCoef(3+iCnt, 14, TmpPrime2.dGet(iCnt));
      WM.PutCoef(3+iCnt, 15, Tmp1.dGet(iCnt));
      WM.PutCoef(3+iCnt, 16, Tmp2.dGet(iCnt));
   }

   /* Derivata dell'equazione di momento, nodo 2 */
   WM.Add(10, 1, MDeltag1);
   WM.Add(10, 4, Mat3x3(MTmp, e3a));
   WM.Sub(10, 7, MDeltag2);
   WM.Sub(10, 10, Mat3x3(e3a, MTmp));

   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutCoef(9+iCnt, 13, -TmpPrime1.dGet(iCnt));
      WM.PutCoef(9+iCnt, 14, -TmpPrime2.dGet(iCnt));
      WM.PutCoef(9+iCnt, 15, -Tmp1.dGet(iCnt));
      WM.PutCoef(9+iCnt, 16, -Tmp2.dGet(iCnt));
   }

   /* Equazioni di vincolo di rotazione: e1b~e3a, e2b~e3a */

   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = Tmp1.dGet(iCnt);
      WM.PutCoef(13, iCnt, d);
      WM.PutCoef(13, 6+iCnt, -d);

      /* Queste sono per la derivata dell'equazione, sono qui solo per
       * ottimizzazione */
      WM.PutCoef(15, 3+iCnt, d);
      WM.PutCoef(15, 9+iCnt, -d);

      d = Tmp2.dGet(iCnt);
      WM.PutCoef(14, iCnt, -d);
      WM.PutCoef(14, 6+iCnt, d);

      /* Queste sono per la derivata dell'equazione, sono qui solo per
       * ottimizzazione */
      WM.PutCoef(16, 3+iCnt, -d);
      WM.PutCoef(16, 9+iCnt, d);
   }

   /* Derivate delle equazioni di vincolo di rotazione: e1b~e3a, e2b~e3a */
   Vec3 O1mO2(Omega1-Omega2);
   TmpPrime1 = e3a.Cross(O1mO2.Cross(e2b));
   TmpPrime2 = e2b.Cross(e3a.Cross(O1mO2));
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutCoef(15, iCnt, TmpPrime1.dGet(iCnt));
      WM.PutCoef(15, 6+iCnt, TmpPrime2.dGet(iCnt));
   }

   TmpPrime1 = e3a.Cross(O1mO2.Cross(e1b));
   TmpPrime2 = e1b.Cross(e3a.Cross(O1mO2));
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutCoef(16, iCnt, TmpPrime1.dGet(iCnt));
      WM.PutCoef(16, 6+iCnt, TmpPrime2.dGet(iCnt));
   }
#endif
   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
ImposedOrientationJoint::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
#if 0
   DEBUGCOUT("Entering ImposedOrientationJoint::InitialAssRes()" << std::endl);

   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);

   /* Indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex()+3;
   integer iNode1FirstVelIndex = iNode1FirstPosIndex+6+3;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex()+3;
   integer iNode2FirstVelIndex = iNode2FirstPosIndex+6+3;
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+2;

   /* Setta gli indici */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WorkVec.PutRowIndex(3+iCnt, iNode1FirstVelIndex+iCnt);
      WorkVec.PutRowIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
      WorkVec.PutRowIndex(9+iCnt, iNode2FirstVelIndex+iCnt);
   }

   for (int iCnt = 1; iCnt <= 4; iCnt++) {
      WorkVec.PutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
   }

   /* Recupera i dati */
   const Mat3x3& R1(pNode1->GetRCurr());
   const Mat3x3& R2(pNode2->GetRCurr());
   const Vec3& Omega1(pNode1->GetWCurr());
   const Vec3& Omega2(pNode2->GetWCurr());

   /* Aggiorna F ed M, che restano anche per InitialAssJac */
   M = Vec3(XCurr.dGetCoef(iFirstReactionIndex+1),
	    XCurr.dGetCoef(iFirstReactionIndex+2),
	    0.);
   Vec3 MPrime(XCurr.dGetCoef(iReactionPrimeIndex+1),
	       XCurr.dGetCoef(iReactionPrimeIndex+2),
	       0.);

   /* Distanza nel sistema globale */
   Mat3x3 R1hTmp(R1*R1h);
   Mat3x3 R2hTmp(R2*R2h);

   Vec3 e3a(R1hTmp.GetVec(3));
   Vec3 e1b(R2hTmp.GetVec(1));
   Vec3 e2b(R2hTmp.GetVec(2));

   /* Ruota il momento e la sua derivata con le matrici della cerniera
    * rispetto ai nodi */
   Vec3 MTmp(e2b*M.dGet(1)-e1b*M.dGet(2));
   Vec3 MPrimeTmp(e3a.Cross(MTmp.Cross(Omega2))+MTmp.Cross(Omega1.Cross(e3a))+
		  e2b.Cross(e3a)*MPrime.dGet(1)+e3a.Cross(e1b)*MPrime.dGet(2));

   /* Equazioni di equilibrio, nodo 1 */
   WorkVec.Sub(1, MTmp.Cross(e3a));

   /* Derivate delle equazioni di equilibrio, nodo 1 */
   WorkVec.Sub(4, MPrimeTmp);

   /* Equazioni di equilibrio, nodo 2 */
   WorkVec.Add(7, MTmp.Cross(e3a));

   /* Derivate delle equazioni di equilibrio, nodo 2 */
   WorkVec.Add(10, MPrimeTmp);

   /* Equazioni di vincolo di rotazione */
   WorkVec.PutCoef(13, e2b.Dot(e3a));
   WorkVec.PutCoef(14, e1b.Dot(e3a));

   /* Derivate delle equazioni di vincolo di rotazione: e1b~e3a, e2b~e3a */
   Vec3 Tmp((Omega1-Omega2).Cross(e3a));
   WorkVec.PutCoef(15, e2b.Dot(Tmp));
   WorkVec.PutCoef(16, e1b.Dot(Tmp));
#endif

	return WorkVec;
}


unsigned int
ImposedOrientationJoint::iGetNumPrivData(void) const
{
	return 5;
}

unsigned int
ImposedOrientationJoint::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);

	unsigned int idx = 0;

	switch (s[0]) {
	case 'w':
		idx++;
	case 'r':
		idx++;
		if (s[1] == 'z') {
			return idx;
		}
		break;

	case 'M':
		idx += 2;

		switch (s[1]) {
		case 'x':
			return idx + 1;

		case 'y':
			return idx + 2;

		case 'z':
			return idx + 3;
		}
	}

	return 0;
}

doublereal
ImposedOrientationJoint::dGetPrivData(unsigned int i) const
{
#if 0
   ASSERT(i >= 1 && i <= iGetNumPrivData());

   switch (i) {
    case 1: {
       Mat3x3 RTmp(((pNode1->GetRCurr()*R1h).Transpose()
			*pNode1->GetRPrev()*R1h).Transpose()
			*((pNode2->GetRCurr()*R2h).Transpose()
			*pNode2->GetRPrev()*R2h));
       Vec3 v(MatR2EulerAngles(RTmp.Transpose()));

       return dTheta + v(3);
    }

    case 2: {
       Mat3x3 R2TmpT((pNode2->GetRCurr()*R2h).Transpose());
       Vec3 v(R2TmpT*(pNode2->GetWCurr()-pNode1->GetWCurr()));

       return v(3);
    }

    case 3:
    case 4:
    case 5:
	    return M(i - 2);
   }

   silent_cerr("ImposedOrientationJoint(" << GetLabel() << "): "
	   "illegal private data " << i << std::endl);
   throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif
	return 0.;
}

/* ImposedOrientationJoint - end */

