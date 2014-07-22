/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <iostream>
#include <limits>

#include "planej.h"
#include "Rot.hh"
#include "hint_impl.h"

/* PlaneHingeJoint - begin */

const unsigned int PlaneHingeJoint::NumSelfDof(5);
const unsigned int PlaneHingeJoint::NumDof(17);

/* Costruttore non banale */
PlaneHingeJoint::PlaneHingeJoint(unsigned int uL, const DofOwner* pDO,
		const StructNode* pN1, const StructNode* pN2,
		const Vec3& dTmp1, const Vec3& dTmp2,
		const Mat3x3& R1hTmp, const Mat3x3& R2hTmp,
		const OrientationDescription& od,
		flag fOut, 
		const bool _calcInitdTheta,
		const doublereal initDTheta,
		const doublereal rr,
		const doublereal pref,
		BasicShapeCoefficient *const sh,
		BasicFriction *const f)
: Elem(uL, fOut), 
Joint(uL, pDO, fOut), 
pNode1(pN1), pNode2(pN2),
d1(dTmp1), R1h(R1hTmp), d2(dTmp2), R2h(R2hTmp), F(Zero3), M(Zero3),
#ifdef USE_NETCDF
Var_Phi(0),
Var_Omega(0),
//Var_MFR(0),
//Var_MU(0),
#endif // USE_NETCDF
calcInitdTheta(_calcInitdTheta), NTheta(0), dTheta(initDTheta), dThetaWrapped(initDTheta),
Sh_c(sh), fc(f), preF(pref), r(rr),
od(od)
{
	NO_OP;
	char * fname = NULL;
	int n = (uL > 0 ? 1 + (int)log10(uL) : 1);
	int len = STRLENOF("hinge") + n + STRLENOF(".out") + 1;
	SAFENEWARR(fname, char, len);
	snprintf(fname, len, "hinge%.*d.out", n, uL);
	SAFEDELETEARR(fname);
}


/* Distruttore banale */
PlaneHingeJoint::~PlaneHingeJoint(void)
{
	if (Sh_c) {
		delete Sh_c;
	}

	if (fc) {
		delete fc;
	}
}

std::ostream&
PlaneHingeJoint::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
			"reaction forces [Fx,Fy,Fz]" << std::endl
		<< prefix << iIndex + 4 << "->" << iIndex + 5 << ": "
			"reaction couples [mx,my]" << std::endl;

	if (bInitial) {
		iIndex += NumSelfDof;
		out
			<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
				"reaction force derivatives [FPx,FPy,FPz]" << std::endl
			<< prefix << iIndex + 4 << "->" << iIndex + 5 << ": "
				"reaction couple derivatives [mPx,mPy]" << std::endl;
	}

	iIndex += NumSelfDof;
	if (fc) {
		integer iFCDofs = fc->iGetNumDof();
		if (iFCDofs > 0) {
			out << prefix << iIndex + 1;
			if (iFCDofs > 1) {
				out << "->" << iIndex + iFCDofs;
			}
			out << ": friction dof(s)" << std::endl
				<< "        ", fc->DescribeDof(out, prefix, bInitial);
		}
	}

	return out;
}

static const char xyz[] = "xyz";

void
PlaneHingeJoint::DescribeDof(std::vector<std::string>& desc, bool bInitial, int i) const
{
	std::ostringstream os;
	os << "PlaneHingeJoint(" << GetLabel() << ")";

	unsigned short nself = NumSelfDof;
	if (bInitial) {
		nself *= 2;
	}
	if (fc && (i == -1 || i >= nself)) {
		fc->DescribeDof(desc, bInitial, i - nself);
		if (i != -1) {
			desc[0] = os.str() + ": " + desc[0];
			return;
		}
	}

	if (i == -1) {
		// move fc desc to the end
		unsigned short nfc = 0;
		if (fc) {
			nfc = desc.size();
		}
		desc.resize(nfc + nself);
		for (unsigned i = nfc; i-- > 0; ) {
			desc[nself + i] = os.str() + ": " + desc[nfc];
		}

		std::string name = os.str();

		for (unsigned i = 0; i < 3; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": reaction force f" << xyz[i];
			desc[i] = os.str();
		}

		for (unsigned i = 0; i < 2; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": reaction couple m" << xyz[i];
			desc[3 + i] = os.str();
		}

		if (bInitial) {
			for (unsigned i = 0; i < 3; i++) {
				os.str(name);
				os.seekp(0, std::ios_base::end);
				os << ": reaction force derivative fP" << xyz[i];
				desc[3 + 2 + i] = os.str();
			}
	
			for (unsigned i = 0; i < 2; i++) {
				os.str(name);
				os.seekp(0, std::ios_base::end);
				os << ": reaction couple derivative mP" << xyz[i];
				desc[3 + 2 + 3 + i] = os.str();
			}
		}

	} else {
		if (i < -1) {
			// error
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (i >= nself) {
			// error
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		desc.resize(1);

		switch (i) {
		case 0:
		case 1:
		case 2:
			os << ": reaction force f" << xyz[i];
			break;

		case 3:
		case 4:
			os << ": reaction couple m" << xyz[i - 3];
			break;

		case 5:
		case 6:
		case 7:
			os << ": reaction force derivative fP" << xyz[i - 3 - 2];
			break;

		case 8:
		case 9:
			os << ": reaction couple derivative mP" << xyz[i - 3 - 2 - 3];
			break;
		}
		desc[0] = os.str();
	}
}

std::ostream&
PlaneHingeJoint::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
			"position constraints [Px1=Px2,Py1=Py2,Pz1=Pz2]" << std::endl
		<< prefix << iIndex + 4 << "->" << iIndex + 5 << ": "
			"orientation constraints [gx1=gx2,gy1=gy2]" << std::endl;

	if (bInitial) {
		iIndex += NumSelfDof;
		out
			<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
				"velocity constraints [vx1=vx2,vy1=vy2,vz1=vz2]" << std::endl
			<< prefix << iIndex + 4 << "->" << iIndex + 5 << ": "
				"angular velocity constraints [wx1=wx2,wy1=wy2]" << std::endl;
	}

	iIndex += NumSelfDof;
	if (fc) {
		integer iFCDofs = fc->iGetNumDof();
		if (iFCDofs > 0) {
			out << prefix << iIndex + 1;
			if (iFCDofs > 1) {
				out << "->" << iIndex + iFCDofs;
			}
			out << ": friction equation(s)" << std::endl
				<< "        ", fc->DescribeEq(out, prefix, bInitial);
		}
	}

	return out;
}

void
PlaneHingeJoint::DescribeEq(std::vector<std::string>& desc, bool bInitial, int i) const
{
	std::ostringstream os;
	os << "PlaneHingeJoint(" << GetLabel() << ")";

	unsigned short nself = NumSelfDof;
	if (bInitial) {
		nself *= 2;
	}
	if (fc && (i == -1 || i >= nself)) {
		fc->DescribeEq(desc, bInitial, i - nself);
		if (i != -1) {
			desc[0] = os.str() + ": " + desc[0];
			return;
		}
	}

	if (i == -1) {
		// move fc desc to the end
		unsigned short nfc = 0;
		if (fc) {
			nfc = desc.size();
		}
		desc.resize(nfc + nself);
		for (unsigned i = nfc; i-- > 0; ) {
			desc[nself + i] = os.str() + ": " + desc[nfc];
		}

		std::string name = os.str();

		for (unsigned i = 0; i < 3; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": position constraint P" << xyz[i];
			desc[i] = os.str();
		}

		for (unsigned i = 0; i < 2; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": orientation constraint g" << xyz[i];
			desc[3 + i] = os.str();
		}

		if (bInitial) {
			for (unsigned i = 0; i < 3; i++) {
				os.str(name);
				os.seekp(0, std::ios_base::end);
				os << ": position constraint derivative v" << xyz[i];
				desc[3 + 2 + i] = os.str();
			}
	
			for (unsigned i = 0; i < 2; i++) {
				os.str(name);
				os.seekp(0, std::ios_base::end);
				os << ": orientation constraint derivative w" << xyz[i];
				desc[3 + 2 + 3 + i] = os.str();
			}
		}

	} else {
		if (i < -1) {
			// error
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (i >= nself) {
			// error
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		desc.resize(1);

		switch (i) {
		case 0:
		case 1:
		case 2:
			os << ": position constraint P" << xyz[i];
			break;

		case 3:
		case 4:
			os << ": orientation constraint g" << xyz[i - 3];
			break;

		case 5:
		case 6:
		case 7:
			os << ": position constraint derivative v" << xyz[i - 3 - 2];
			break;

		case 8:
		case 9:
			os << ": orientation constraint derivative w" << xyz[i - 3 - 2 - 3];
			break;
		}
		desc[0] = os.str();
	}
}

void
PlaneHingeJoint::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
	if (ph) {
		for (unsigned i = 0; i < ph->size(); i++) {
			Joint::JointHint *pjh = dynamic_cast<Joint::JointHint *>((*ph)[i]);

			if (pjh == 0) {
				continue;
			}

			if (dynamic_cast<Joint::OffsetHint<1> *>(pjh)) {
				const Mat3x3& R1(pNode1->GetRCurr());
				Vec3 dTmp2(pNode2->GetRCurr()*d2);
   
				d1 = R1.MulTV(pNode2->GetXCurr() + dTmp2 - pNode1->GetXCurr());

			} else if (dynamic_cast<Joint::OffsetHint<2> *>(pjh)) {
				const Mat3x3& R2(pNode2->GetRCurr());
				Vec3 dTmp1(pNode1->GetRCurr()*d1);
   
				d2 = R2.MulTV(pNode1->GetXCurr() + dTmp1 - pNode2->GetXCurr());

			} else if (dynamic_cast<Joint::HingeHint<1> *>(pjh)) {
				R1h = pNode1->GetRCurr().MulTM(pNode2->GetRCurr()*R2h);

			} else if (dynamic_cast<Joint::HingeHint<2> *>(pjh)) {
				R2h = pNode2->GetRCurr().MulTM(pNode1->GetRCurr()*R1h);

			} else if (dynamic_cast<Joint::ReactionsHint *>(pjh)) {
				/* TODO */
			}
		}
	}

	if (calcInitdTheta) {
		Vec3 v(RotManip::VecRot((pNode1->GetRCurr()*R1h).MulTM(pNode2->GetRCurr()*R2h)));
		
		dThetaWrapped = dTheta = v.dGet(3);
	}
	
#if 0
	std::cerr << "F: " << F << std::endl;
	std::cerr << "M: " << M << std::endl;
#endif
	
	integer iFirstReactionIndex = iGetFirstIndex();
	X.Put(iFirstReactionIndex + 1, F);
	X.PutCoef(iFirstReactionIndex + 4, M.dGet(1));
	X.PutCoef(iFirstReactionIndex + 5, M.dGet(2));

	if (fc) {
		fc->SetValue(pDM, X, XP, ph, iGetFirstIndex() + NumSelfDof);
	}
}

Hint *
PlaneHingeJoint::ParseHint(DataManager *pDM, const char *s) const
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

	} else if (strncasecmp(s, "hinge{" /*}*/, STRLENOF("hinge{" /*}*/)) == 0) {
		s += STRLENOF("hinge{" /*}*/);

		if (strcmp(&s[1], /* { */ "}") != 0) {
			return 0;
		}

		switch (s[0]) {
		case '1':
			return new Joint::HingeHint<1>;

		case '2':
			return new Joint::HingeHint<2>;
		}

	} else if (fc) {
		return fc->ParseHint(pDM, s);
	}

	return 0;
}

void
PlaneHingeJoint::AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP)
{
	Vec3 v(RotManip::VecRot((pNode1->GetRCurr()*R1h).MulTM(pNode2->GetRCurr()*R2h)));
	doublereal dThetaTmp(v(3));

	// unwrap
	if (dThetaTmp - dThetaWrapped < -M_PI) {
		NTheta++;
	}

	if (dThetaTmp - dThetaWrapped > M_PI) {
		NTheta--;
	}

	// save new wrapped angle
	dThetaWrapped = dThetaTmp;

	// compute new unwrapped angle
	dTheta = 2*M_PI*NTheta + dThetaWrapped;

	if (fc) {
		Mat3x3 R1(pNode1->GetRCurr());
		Mat3x3 R1hTmp(R1*R1h);
		Vec3 e3a(R1hTmp.GetVec(3));
		Vec3 Omega1(pNode1->GetWCurr());
		Vec3 Omega2(pNode2->GetWCurr());
		//relative velocity
		doublereal v = (Omega1-Omega2).Dot(e3a)*r;
		//reaction norm
		doublereal modF = std::max(F.Norm(), preF);;
		fc->AfterConvergence(modF,v,X,XP,iGetFirstIndex()+NumSelfDof);
	}
}

/* Funzione che legge lo stato iniziale dal file di input */
void
PlaneHingeJoint::ReadInitialState(MBDynParser& HP)
{
	F = HP.GetVec3();
	M = Vec3(HP.GetReal(), HP.GetReal(), 0.);
}


/* Contributo al file di restart */
std::ostream&
PlaneHingeJoint::Restart(std::ostream& out) const
{
   Joint::Restart(out) << ", revolute hinge, "
     << pNode1->GetLabel() << ", reference, node, ",
     d1.Write(out, ", ")
     << ", hinge, reference, node, 1, ", (R1h.GetVec(1)).Write(out, ", ")
     << ", 2, ", (R1h.GetVec(2)).Write(out, ", ") << ", "
     << pNode2->GetLabel() << ", reference, node, ",
     d2.Write(out, ", ")
     << ", hinge, reference, node, 1, ", (R2h.GetVec(1)).Write(out, ", ")
     << ", 2, ", (R2h.GetVec(2)).Write(out, ", ") << ", " 
     << "initial theta, " << dTheta << ", "
     << "initial state, ", F.Write(out, ", ") 
     << ", " << M.dGet(1) << ", " << M.dGet(2) << ';' << std::endl;
   
   return out;
}


/* Assemblaggio jacobiano */
VariableSubMatrixHandler& 
PlaneHingeJoint::AssJac(VariableSubMatrixHandler& WorkMat,
			    doublereal dCoef,
			    const VectorHandler& XCurr,
			    const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering PlaneHingeJoint::AssJac()" << std::endl);
   
   /* Setta la sottomatrice come piena (e' un po' dispersivo, ma lo jacobiano 
    * e' complicato */					
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
   integer iFirstReactionIndex = iGetFirstIndex();

   /* Setta gli indici delle equazioni */
   for (unsigned int iCnt = 1; iCnt <= 6; iCnt++) {	
      WM.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
      WM.PutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
   }
   
   for (unsigned int iCnt = 1; iCnt <= iGetNumDof(); iCnt++) {	
      WM.PutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
      WM.PutColIndex(12+iCnt, iFirstReactionIndex+iCnt);
   }
   
   /* Recupera i dati che servono */
   const Mat3x3& R1(pNode1->GetRRef());
   const Mat3x3& R2(pNode2->GetRRef());   
   Vec3 d1Tmp(R1*d1);
   Vec3 d2Tmp(R2*d2);
   Mat3x3 R1hTmp(R1*R1h);
   Mat3x3 R2hTmp(R2*R2h);
   
   /* Suppongo che le reazioni F, M siano gia' state aggiornate da AssRes;
    * ricordo che la forza F e' nel sistema globale, mentre la coppia M
    * e' nel sistema locale ed il terzo termine, M(3), e' nullo in quanto
    * diretto come l'asse attorno al quale la rotazione e' consentita */
   
      
   /* 
    * La cerniera piana ha le prime 3 equazioni uguali alla cerniera sferica;
    * inoltre ha due equazioni che affermano la coincidenza dell'asse 3 del
    * riferimento solidale con la cerniera visto dai due nodi.
    * 
    *      (R1 * R1h * e1)^T * (R2 * R2h * e3) = 0
    *      (R1 * R1h * e2)^T * (R2 * R2h * e3) = 0
    * 
    * A queste equazioni corrisponde una reazione di coppia agente attorno 
    * agli assi 1 e 2 del riferimento della cerniera. La coppia attorno 
    * all'asse 3 e' nulla per definizione. Quindi la coppia nel sistema 
    * globale e':
    *      -R1 * R1h * M       per il nodo 1,
    *       R2 * R2h * M       per il nodo 2
    * 
    * 
    *       xa   ga                   xb   gb                     F     M 
    * Qa |  0    0                     0    0                     I     0  | | xa |   | -F           |
    * Ga |  0    c*(F/\da/\-(Sa*M)/\)  0    0                     da/\  Sa | | ga |   | -da/\F-Sa*M |
    * Qb |  0    0                     0    0                    -I     0  | | xb | = |  F           |
    * Gb |  0    0                     0   -c*(F/\db/\-(Sb*M)/\) -db/\ -Sb | | gb |   |  db/\F+Sb*M |
    * F  | -c*I  c*da/\                c*I -c*db/\                0     0  | | F  |   |  xa+da-xb-db |
    * M1 |  0    c*Tb1^T*Ta3/\         0    c*Ta3^T*Tb1/\         0     0  | | M  |   |  Sb^T*Ta3    |
    * M2 |  0    c*Tb2^T*Ta3/\         0    c*Ta3^T*Tb2/\         0     0  | 
    * 
    * con da = R1*d01, db = R2*d02, c = dCoef,
    * Sa = R1*R1h*[e1,e2], Sb = R2*R2h*[e1, e2],
    * Ta3 = R1*R1h*e3, Tb1 = R2*R2h*e1, Tb2 = R2*R2h*e2
    */

   /* Contributo della forza alle equazioni di equilibrio dei due nodi */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WM.PutCoef(iCnt, 12+iCnt, 1.);
      WM.PutCoef(6+iCnt, 12+iCnt, -1.);
   }
   
   WM.Add(4, 13, Mat3x3(MatCross, d1Tmp));
   WM.Sub(10, 13, Mat3x3(MatCross, d2Tmp));   
   
   /* Moltiplica la forza ed il momento per il coefficiente
    * del metodo */
   Vec3 FTmp = F*dCoef;
   Vec3 MTmp = M*dCoef;

   Vec3 e3a(R1hTmp.GetVec(3));
   Vec3 e1b(R2hTmp.GetVec(1));
   Vec3 e2b(R2hTmp.GetVec(2));
   MTmp = e2b*MTmp.dGet(1)-e1b*MTmp.dGet(2);
   
   Mat3x3 MWedgee3aWedge(MatCrossCross, MTmp, e3a);
   Mat3x3 e3aWedgeMWedge(MatCrossCross, e3a, MTmp);
   
   WM.Add(4, 4, Mat3x3(MatCrossCross, FTmp, d1Tmp) - MWedgee3aWedge);
   WM.Add(4, 10, e3aWedgeMWedge);
   
   WM.Add(10, 4, MWedgee3aWedge);   
   WM.Sub(10, 10, Mat3x3(MatCrossCross, FTmp, d2Tmp) + e3aWedgeMWedge);      
   
   /* Contributo del momento alle equazioni di equilibrio dei nodi */
   Vec3 Tmp1(e2b.Cross(e3a));
   Vec3 Tmp2(e3a.Cross(e1b));
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = Tmp1.dGet(iCnt);
      WM.PutCoef(3+iCnt, 16, d);
      WM.PutCoef(9+iCnt, 16, -d);
      d = Tmp2.dGet(iCnt);
      WM.PutCoef(3+iCnt, 17, d);
      WM.PutCoef(9+iCnt, 17, -d);
   }         
   
   /* Modifica: divido le equazioni di vincolo per dCoef */
   
   /* Equazioni di vincolo degli spostamenti */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutCoef(12+iCnt, iCnt, -1.);
      WM.PutCoef(12+iCnt, 6+iCnt, 1.);
   }
   
   WM.Add(13, 4, Mat3x3(MatCross, d1Tmp));
   WM.Sub(13, 10, Mat3x3(MatCross, d2Tmp));
   
   /* Equazione di vincolo del momento
    * 
    * Attenzione: bisogna scrivere il vettore trasposto
    *   (Sb[1]^T*(Sa[3]/\))*dCoef
    * Questo pero' e' uguale a:
    *   (-Sa[3]/\*Sb[1])^T*dCoef,
    * che puo' essere ulteriormente semplificato:
    *   (Sa[3].Cross(Sb[1])*(-dCoef))^T;
    */
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = Tmp1.dGet(iCnt);
      WM.PutCoef(16, 3+iCnt, d);
      WM.PutCoef(16, 9+iCnt, -d);
      d = Tmp2.dGet(iCnt);
      WM.PutCoef(17, 3+iCnt, -d);
      WM.PutCoef(17, 9+iCnt, d);
   }   

   if (fc) {
      //retrive
          //friction coef
      doublereal f = fc->fc();
          //shape function
      doublereal shc = Sh_c->Sh_c();
          //omega and omega rif
      const Vec3& Omega1(pNode1->GetWCurr());
      const Vec3& Omega2(pNode2->GetWCurr());
      // const Vec3& Omega1r(pNode1->GetWRef());
      // const Vec3& Omega2r(pNode2->GetWRef());   
      //compute 
          //relative velocity
      doublereal v = (Omega1-Omega2).Dot(e3a)*r;
          //reaction norm
      doublereal modF = std::max(F.Norm(), preF);
          //reaction moment
      //doublereal M3 = shc*modF*r;
      
      ExpandableRowVector dfc;
      ExpandableRowVector dF;
      ExpandableRowVector dv;
          //variation of reaction force
      dF.ReDim(3);
      if ((modF == 0.) or (F.Norm() < preF)) {
          dF.Set(Vec3(Zero3),1,12+1);
      } else {
          dF.Set(F/modF,1,12+1);
      }
          //variation of relative velocity
      dv.ReDim(6);
      
/* old (wrong?) relative velocity linearization */

//       dv.Set((e3a.dGet(1)*1.-( e3a.dGet(2)*Omega1r.dGet(3)-e3a.dGet(3)*Omega1r.dGet(2))*dCoef)*r,1,0+4);
//       dv.Set((e3a.dGet(2)*1.-(-e3a.dGet(1)*Omega1r.dGet(3)+e3a.dGet(3)*Omega1r.dGet(1))*dCoef)*r,2,0+5);
//       dv.Set((e3a.dGet(3)*1.-( e3a.dGet(1)*Omega1r.dGet(2)-e3a.dGet(2)*Omega1r.dGet(1))*dCoef)*r,3,0+6);
//       
//       dv.Set(-(e3a.dGet(1)*1.-( e3a.dGet(2)*Omega2r.dGet(3)-e3a.dGet(3)*Omega2r.dGet(2))*dCoef)*r,4,6+4);
//       dv.Set(-(e3a.dGet(2)*1.-(-e3a.dGet(1)*Omega2r.dGet(3)+e3a.dGet(3)*Omega2r.dGet(1))*dCoef)*r,5,6+5);
//       dv.Set(-(e3a.dGet(3)*1.-( e3a.dGet(1)*Omega2r.dGet(2)-e3a.dGet(2)*Omega2r.dGet(1))*dCoef)*r,6,6+6);

/* new (exact?) relative velocity linearization */
// 
//       ExpandableRowVector domega11, domega12, domega13;
//       ExpandableRowVector domega21, domega22, domega23;
//       domega11.ReDim(3); domega12.ReDim(3); domega13.ReDim(3);
//       domega21.ReDim(3); domega22.ReDim(3); domega23.ReDim(3);
//       
//       domega11.Set(1., 1, 0+4);
//           domega11.Set( Omega1r.dGet(3)*dCoef, 2, 0+5);
//           domega11.Set(-Omega1r.dGet(2)*dCoef, 3, 0+6);
//       domega21.Set(1., 1, 6+4);
//           domega21.Set( Omega2r.dGet(3)*dCoef, 2, 6+5);
//           domega21.Set(-Omega2r.dGet(2)*dCoef, 3, 6+6);
//       domega12.Set(1., 1, 0+5);
//           domega12.Set(-Omega1r.dGet(3)*dCoef, 2, 0+4);
//           domega12.Set( Omega1r.dGet(1)*dCoef, 3, 0+6);
//       domega22.Set(1., 1, 6+5);
//           domega22.Set(-Omega2r.dGet(3)*dCoef, 2, 6+4);
//           domega22.Set( Omega2r.dGet(1)*dCoef, 3, 6+6);
//       domega13.Set(1., 1, 0+6);
//           domega13.Set( Omega1r.dGet(2)*dCoef, 2, 0+4);
//           domega13.Set(-Omega1r.dGet(1)*dCoef, 3, 0+5);
//       domega23.Set(1., 1, 6+6);
//           domega23.Set( Omega2r.dGet(2)*dCoef, 2, 6+4);
//           domega23.Set(-Omega2r.dGet(1)*dCoef, 3, 6+5);
// 
//       Vec3 domega = Omega1-Omega2;
//       dv.Set((e3a.dGet(1)*1.-( 
//       		e3a.dGet(2)*(Omega1.dGet(3)-Omega2.dGet(3))-
// 		e3a.dGet(3)*(domega.dGet(2)))*dCoef)*r,1);
// 		dv.Link(1, &domega11);
//       dv.Set((e3a.dGet(2)*1.-(
//       		-e3a.dGet(1)*(Omega1.dGet(3)-Omega2.dGet(3))+
// 		e3a.dGet(3)*(domega.dGet(1)))*dCoef)*r,2);
// 		dv.Link(2, &domega12);
//       dv.Set((e3a.dGet(3)*1.-( 
//       		e3a.dGet(1)*(Omega1.dGet(2)-Omega2.dGet(2))-
// 		e3a.dGet(2)*(domega.dGet(1)))*dCoef)*r,3);
// 		dv.Link(3, &domega13);
// 
//       dv.Set(-(e3a.dGet(1)*1.)*r,4,6+4); dv.Link(4, &domega21);
//       dv.Set(-(e3a.dGet(2)*1.)*r,5,6+5); dv.Link(5, &domega22);
//       dv.Set(-(e3a.dGet(3)*1.)*r,6,6+6); dv.Link(6, &domega23);


/* new (approximate: assume constant triads orientations) 
 * relative velocity linearization 
*/

      dv.Set(e3a*r,1, 0+4);
      dv.Set(-e3a*r,4, 6+4);

      //assemble friction states
      fc->AssJac(WM,dfc,12+NumSelfDof,iFirstReactionIndex+NumSelfDof,dCoef,modF,v,
      		XCurr,XPrimeCurr,dF,dv);
      ExpandableMatrix dM3;
      ExpandableRowVector dShc;
      //compute 
          //variation of shape function
      Sh_c->dSh_c(dShc,f,modF,v,dfc,dF,dv);
          //variation of moment component
      dM3.ReDim(3,2);
      dM3.SetBlockDim(1,1);
      dM3.SetBlockDim(2,1);
      dM3.Set(e3a*shc*r,1,1); dM3.Link(1,&dF);
      dM3.Set(e3a*modF*r,1,2); dM3.Link(2,&dShc);
      //assemble first node
          //variation of moment component
      dM3.Add(WM, 4, 1.);
      //assemble second node
          //variation of moment component
      dM3.Sub(WM, 6+4, 1.);
   }
   
   return WorkMat;
}


/* Assemblaggio residuo */
SubVectorHandler& PlaneHingeJoint::AssRes(SubVectorHandler& WorkVec,
					  doublereal dCoef,
					  const VectorHandler& XCurr, 
					  const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering PlaneHingeJoint::AssRes()" << std::endl);
      
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);
 
   /* Indici */
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   
   /* Indici dei nodi */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {	
      WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.PutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
   }   
   
   /* Indici del vincolo */
   for (unsigned int iCnt = 1; iCnt <= iGetNumDof(); iCnt++) {
      WorkVec.PutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
   }

   /* Aggiorna i dati propri */
   F = Vec3(XCurr, iFirstReactionIndex+1);
   M = Vec3(XCurr(iFirstReactionIndex+4),
	    XCurr(iFirstReactionIndex+5),
	    0.);

   /*
    * FIXME: provare a mettere "modificatori" di forza/momento sui gdl
    * residui: attrito, rigidezze e smorzamenti, ecc.
    */
   
   /* Recupera i dati */
   const Vec3& x1(pNode1->GetXCurr());
   const Vec3& x2(pNode2->GetXCurr());
   const Mat3x3& R1(pNode1->GetRCurr());
   const Mat3x3& R2(pNode2->GetRCurr());
   
   /* Costruisce i dati propri nella configurazione corrente */
   Vec3 dTmp1(R1*d1);
   Vec3 dTmp2(R2*d2);
   Mat3x3 R1hTmp(R1*R1h);
   Mat3x3 R2hTmp(R2*R2h);
   
   Vec3 e3a(R1hTmp.GetVec(3));
   Vec3 e1b(R2hTmp.GetVec(1));
   Vec3 e2b(R2hTmp.GetVec(2));
   
   Vec3 MTmp(e2b.Cross(e3a)*M.dGet(1)+e3a.Cross(e1b)*M.dGet(2));
   
   /* Equazioni di equilibrio, nodo 1 */
   WorkVec.Sub(1, F);
   WorkVec.Add(4, F.Cross(dTmp1)-MTmp); /* Sfrutto  F/\d = -d/\F */
   
   /* Equazioni di equilibrio, nodo 2 */
   WorkVec.Add(7, F);
   WorkVec.Add(10, dTmp2.Cross(F)+MTmp);

   /* Modifica: divido le equazioni di vincolo per dCoef */
   ASSERT(dCoef != 0.);
	
   /* Equazione di vincolo di posizione */
   WorkVec.Add(13, (x1+dTmp1-x2-dTmp2)/dCoef);
      
   /* Equazioni di vincolo di rotazione */
   Vec3 Tmp = R1hTmp.GetVec(3);
   WorkVec.PutCoef(16, Tmp.Dot(R2hTmp.GetVec(2))/dCoef);
   WorkVec.PutCoef(17, Tmp.Dot(R2hTmp.GetVec(1))/dCoef);

   if (fc) {
      bool ChangeJac(false);
      const Vec3& Omega1(pNode1->GetWCurr());
      const Vec3& Omega2(pNode2->GetWCurr());
      doublereal v = (Omega1-Omega2).Dot(e3a)*r;
      doublereal modF = std::max(F.Norm(), preF);
      try {
          fc->AssRes(WorkVec,12+NumSelfDof,iFirstReactionIndex+NumSelfDof,modF,v,XCurr,XPrimeCurr);
      }
      catch (Elem::ChangedEquationStructure) {
          ChangeJac = true;
      }
      doublereal f = fc->fc();
      doublereal shc = Sh_c->Sh_c(f,modF,v);
      M3 = shc*modF*r;
      WorkVec.Sub(4,e3a*M3);
      WorkVec.Add(10,e3a*M3);
//!!!!!!!!!!!!!!
//      M += e3a*M3;
      if (ChangeJac) {
          throw Elem::ChangedEquationStructure(MBDYN_EXCEPT_ARGS);
      }
   }
   
   return WorkVec;
}

unsigned int PlaneHingeJoint::iGetNumDof(void) const {
   unsigned int i = NumSelfDof;
   if (fc) {
       i+=fc->iGetNumDof();
   } 
   return i;
};


DofOrder::Order
PlaneHingeJoint::GetDofType(unsigned int i) const {
   ASSERT(i >= 0 && i < iGetNumDof());
   if (i<NumSelfDof) {
       return DofOrder::ALGEBRAIC; 
   } else {
       return fc->GetDofType(i-NumSelfDof);
   }
};

DofOrder::Order
PlaneHingeJoint::GetEqType(unsigned int i) const
{
	ASSERTMSGBREAK(i < iGetNumDof(), 
		"INDEX ERROR in PlaneHingeJoint::GetEqType");
   if (i<NumSelfDof) {
       return DofOrder::ALGEBRAIC; 
   } else {
       return fc->GetEqType(i-NumSelfDof);
   }
}

void
PlaneHingeJoint::OutputPrepare(OutputHandler& OH)
{
	if (fToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			std::string name;
			OutputPrepare_int("revolute hinge", OH, name);

			Var_Phi = OH.CreateRotationVar(name, "", od, "global");

			Var_Omega = OH.CreateVar<Vec3>(name + "Omega", "radian/s",
				"local relative angular velocity (x, y, z)");

/* TODO
			Var_MFR = OH.CreateVar<doublereal>(name + "MFR", "Nm",
				"friciton moment ");

			Var_MU = OH.CreateVar<doublereal>(name + "MU", "--",
					"friction model specific data: friction coefficient?");
*/
		}
#endif // USE_NETCDF
	}
}

/* Output (da mettere a punto) */
void PlaneHingeJoint::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) {
      Mat3x3 R2Tmp(pNode2->GetRCurr()*R2h);
      Mat3x3 RTmp((pNode1->GetRCurr()*R1h).MulTM(R2Tmp));
      Vec3 OmegaTmp(R2Tmp.MulTV(pNode2->GetWCurr()-pNode1->GetWCurr()));
		Vec3 E;
		switch (od) {
		case EULER_123:
			E = MatR2EulerAngles123(RTmp)*dRaDegr;
			break;

		case EULER_313:
			E = MatR2EulerAngles313(RTmp)*dRaDegr;
			break;

		case EULER_321:
			E = MatR2EulerAngles321(RTmp)*dRaDegr;
			break;

		case ORIENTATION_VECTOR:
			E = RotManip::VecRot(RTmp);
			break;

		case ORIENTATION_MATRIX:
			break;

		default:
			/* impossible */
			break;
		}

#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			Var_F_local->put_rec((R2Tmp.MulTV(F)).pGetVec(), OH.GetCurrentStep());
			Var_M_local->put_rec(M.pGetVec(), OH.GetCurrentStep());
			Var_F_global->put_rec(F.pGetVec(), OH.GetCurrentStep());
			Var_M_global->put_rec((R2Tmp*M).pGetVec(), OH.GetCurrentStep());

			switch (od) {
			case EULER_123:
			case EULER_313:
			case EULER_321:
			case ORIENTATION_VECTOR:
				Var_Phi->put_rec(E.pGetVec(), OH.GetCurrentStep());
				break;

			case ORIENTATION_MATRIX:
				Var_Phi->put_rec(RTmp.pGetMat(), OH.GetCurrentStep());
				break;

			default:
				/* impossible */
				break;
			}

			Var_Omega->put_rec(OmegaTmp.pGetVec(), OH.GetCurrentStep());
/*
			if (fc) {
					Var_MFR->put_rec(&M3, OH.GetCurrentStep());
					Var_MU->put_rec(fc->fc(), OH.GetCurrentStep());
			}
			else
			{
				Var_MFR->put_rec(0, OH.GetCurrentStep());
				Var_MU->put_rec(0, OH.GetCurrentStep());
			}
*/
		}
#endif // USE_NETCDF
		if (OH.UseText(OutputHandler::JOINTS)) {
			  std::ostream &of = Joint::Output(OH.Joints(), "PlaneHinge", GetLabel(),
					R2Tmp.MulTV(F), M, F, R2Tmp*M)
			<< " ";

			switch (od) {
			case EULER_123:
			case EULER_313:
			case EULER_321:
			case ORIENTATION_VECTOR:
				of << E;
				break;

			case ORIENTATION_MATRIX:
				of << RTmp;
				break;

			default:
				/* impossible */
				break;
			}

			of << " " << R2Tmp.MulTV(pNode2->GetWCurr()-pNode1->GetWCurr());
			  if (fc) {
				  of << " " << M3 << " " << fc->fc();
			  }
			  of << std::endl;
		}
   }
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
PlaneHingeJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
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
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstVelIndex = iNode1FirstPosIndex+6;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstVelIndex = iNode2FirstPosIndex+6;
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+5;
   
   /* Nota: le reazioni vincolari sono: 
    * Forza,       3 incognite, riferimento globale, 
    * Momento,     2 incognite, riferimento locale
    */

   /* Setta gli indici dei nodi */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutRowIndex(6+iCnt, iNode1FirstVelIndex+iCnt);
      WM.PutColIndex(6+iCnt, iNode1FirstVelIndex+iCnt);
      WM.PutRowIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
      WM.PutColIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
      WM.PutRowIndex(18+iCnt, iNode2FirstVelIndex+iCnt);
      WM.PutColIndex(18+iCnt, iNode2FirstVelIndex+iCnt);
   }
   
   /* Setta gli indici delle reazioni */
   for (int iCnt = 1; iCnt <= 10; iCnt++) {
      WM.PutRowIndex(24+iCnt, iFirstReactionIndex+iCnt);
      WM.PutColIndex(24+iCnt, iFirstReactionIndex+iCnt);	
   }   

   
   /* Recupera i dati */
   const Mat3x3& R1(pNode1->GetRRef());
   const Mat3x3& R2(pNode2->GetRRef());
   const Vec3& Omega1(pNode1->GetWRef());
   const Vec3& Omega2(pNode2->GetWRef());
   
   /* F ed M sono gia' state aggiornate da InitialAssRes */
   Vec3 FPrime(XCurr, iReactionPrimeIndex+1);
   Vec3 MPrime(XCurr(iReactionPrimeIndex+4),
	       XCurr(iReactionPrimeIndex+5),
	       0.);
   
   /* Matrici identita' */   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      /* Contributo di forza all'equazione della forza, nodo 1 */
      WM.PutCoef(iCnt, 24+iCnt, 1.);
      
      /* Contrib. di der. di forza all'eq. della der. della forza, nodo 1 */
      WM.PutCoef(6+iCnt, 29+iCnt, 1.);
      
      /* Contributo di forza all'equazione della forza, nodo 2 */
      WM.PutCoef(12+iCnt, 24+iCnt, -1.);
      
      /* Contrib. di der. di forza all'eq. della der. della forza, nodo 2 */
      WM.PutCoef(18+iCnt, 29+iCnt, -1.);
      
      /* Equazione di vincolo, nodo 1 */
      WM.PutCoef(24+iCnt, iCnt, -1.);
      
      /* Derivata dell'equazione di vincolo, nodo 1 */
      WM.PutCoef(29+iCnt, 6+iCnt, -1.);
      
      /* Equazione di vincolo, nodo 2 */
      WM.PutCoef(24+iCnt, 12+iCnt, 1.);
      
      /* Derivata dell'equazione di vincolo, nodo 2 */
      WM.PutCoef(29+iCnt, 18+iCnt, 1.);
   }
      
   /* Distanze e matrici di rotazione dai nodi alla cerniera 
    * nel sistema globale */
   Vec3 d1Tmp(R1*d1);
   Vec3 d2Tmp(R2*d2);
   Mat3x3 R1hTmp(R1*R1h);
   Mat3x3 R2hTmp(R2*R2h);
   
   Vec3 e3a(R1hTmp.GetVec(3));
   Vec3 e1b(R2hTmp.GetVec(1));
   Vec3 e2b(R2hTmp.GetVec(2));
   
   /* Ruota il momento e la sua derivata con le matrici della cerniera 
    * rispetto ai nodi */
   Vec3 MTmp(e2b*M.dGet(1)-e1b*M.dGet(2));
   Vec3 MPrimeTmp(e2b*MPrime.dGet(1)-e1b*MPrime.dGet(2));

   /* Matrici F/\d1/\, -F/\d2/\ */
   Mat3x3 FWedged1Wedge(MatCrossCross, F, d1Tmp);
   Mat3x3 FWedged2Wedge(MatCrossCross, F, -d2Tmp);
   
   /* Matrici (omega1/\d1)/\, -(omega2/\d2)/\ */
   Mat3x3 O1Wedged1Wedge(MatCross, Omega1.Cross(d1Tmp));
   Mat3x3 O2Wedged2Wedge(MatCross, d2Tmp.Cross(Omega2));
   
   Mat3x3 MDeltag1((Mat3x3(MatCross, Omega2.Cross(MTmp) + MPrimeTmp)
		   + Mat3x3(MatCrossCross, MTmp, Omega1))*Mat3x3(MatCross, e3a));
   Mat3x3 MDeltag2(Mat3x3(MatCrossCross, Omega1.Cross(e3a), MTmp)
		   + Mat3x3(MatCrossCross, e3a, MPrimeTmp)
		   + e3a.Cross(Mat3x3(MatCrossCross, Omega2, MTmp)));

   /* Vettori temporanei */
   Vec3 Tmp1(e2b.Cross(e3a));
   Vec3 Tmp2(e3a.Cross(e1b));
   
   /* Prodotto vettore tra il versore 3 della cerniera secondo il nodo 1
    * ed il versore 1 della cerniera secondo il nodo 2. A convergenza
    * devono essere ortogonali, quindi il loro prodotto vettore deve essere 
    * unitario */

   /* Error handling: il programma si ferma, pero' segnala dov'e' l'errore */
   if (Tmp1.Dot() <= std::numeric_limits<doublereal>::epsilon() || Tmp2.Dot() <= std::numeric_limits<doublereal>::epsilon()) {
      silent_cerr("PlaneHingeJoint(" << GetLabel() << "): "
	      "first and second node hinge axes are (nearly) orthogonal"
	      << std::endl);
      throw ErrNullNorm(MBDYN_EXCEPT_ARGS);
   }   
   
   Vec3 TmpPrime1(e2b.Cross(Omega1.Cross(e3a))-e3a.Cross(Omega2.Cross(e2b)));
   Vec3 TmpPrime2(e3a.Cross(Omega2.Cross(e1b))-e1b.Cross(Omega1.Cross(e3a)));
   
   /* Equazione di momento, nodo 1 */
   WM.Add(4, 4, FWedged1Wedge - Mat3x3(MatCrossCross, MTmp, e3a));
   WM.Add(4, 16, Mat3x3(MatCrossCross, e3a, MTmp));
   WM.Add(4, 25, Mat3x3(MatCross, d1Tmp));
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutCoef(3+iCnt, 28, Tmp1(iCnt));
      WM.PutCoef(3+iCnt, 29, Tmp2(iCnt));	
   }
   
   /* Equazione di momento, nodo 2 */
   WM.Add(16, 4, Mat3x3(MatCrossCross, MTmp, e3a));
   WM.Add(16, 16, FWedged2Wedge - Mat3x3(MatCrossCross, e3a, MTmp));
   WM.Sub(16, 25, Mat3x3(MatCross, d2Tmp));
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutCoef(15+iCnt, 28, -Tmp1(iCnt));
      WM.PutCoef(15+iCnt, 29, -Tmp2(iCnt));	
   }
   
   /* Derivata dell'equazione di momento, nodo 1 */
   WM.Add(10, 4, (Mat3x3(MatCross, FPrime) + Mat3x3(MatCrossCross, F, Omega1))*Mat3x3(MatCross, d1Tmp) - MDeltag1);
   WM.Add(10, 10, FWedged1Wedge - Mat3x3(MatCrossCross, MTmp, e3a));
   WM.Add(10, 16, MDeltag2);
   WM.Add(10, 22, Mat3x3(MatCrossCross, e3a, MTmp));
   WM.Add(10, 25, O1Wedged1Wedge);
   WM.Add(10, 30, Mat3x3(MatCross, d1Tmp));
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutCoef(9+iCnt, 28, TmpPrime1(iCnt));
      WM.PutCoef(9+iCnt, 29, TmpPrime2(iCnt));
      WM.PutCoef(9+iCnt, 33, Tmp1(iCnt));
      WM.PutCoef(9+iCnt, 34, Tmp2(iCnt));	
   }
   
   /* Derivata dell'equazione di momento, nodo 2 */
   WM.Add(22, 4, MDeltag1);
   WM.Add(22, 10, Mat3x3(MatCrossCross, MTmp, e3a));
   WM.Sub(22, 16, (Mat3x3(MatCross, FPrime) + Mat3x3(MatCrossCross, F, Omega2))*Mat3x3(MatCross, d2Tmp) + MDeltag2);
   WM.Add(22, 22, FWedged2Wedge - Mat3x3(MatCrossCross, e3a, MTmp));
   WM.Add(22, 25, O2Wedged2Wedge);
   WM.Sub(22, 30, Mat3x3(MatCross, d2Tmp));

   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutCoef(21+iCnt, 28, -TmpPrime1(iCnt));
      WM.PutCoef(21+iCnt, 29, -TmpPrime2(iCnt));	
      WM.PutCoef(21+iCnt, 33, -Tmp1(iCnt));
      WM.PutCoef(21+iCnt, 34, -Tmp2(iCnt));
   }
   
   /* Equazione di vincolo di posizione */
   WM.Add(25, 4, Mat3x3(MatCross, d1Tmp));
   WM.Sub(25, 16, Mat3x3(MatCross, d2Tmp));
   
   /* Derivata dell'equazione di vincolo di posizione */
   WM.Add(30, 4, O1Wedged1Wedge);
   WM.Add(30, 10, Mat3x3(MatCross, d1Tmp));
   WM.Add(30, 16, O2Wedged2Wedge);
   WM.Sub(30, 22, Mat3x3(MatCross, d2Tmp));
   
   /* Equazioni di vincolo di rotazione: e1b~e3a, e2b~e3a */
            
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = Tmp1(iCnt);
      WM.PutCoef(28, 3+iCnt, d);
      WM.PutCoef(28, 15+iCnt, -d);
      
      /* Queste sono per la derivata dell'equazione, sono qui solo per 
       * ottimizzazione */
      WM.PutCoef(33, 9+iCnt, d);
      WM.PutCoef(33, 21+iCnt, -d);
      
      d = Tmp2.dGet(iCnt);
      WM.PutCoef(29, 3+iCnt, -d);
      WM.PutCoef(29, 15+iCnt, d);
      
      /* Queste sono per la derivata dell'equazione, sono qui solo per 
       * ottimizzazione */
      WM.PutCoef(34, 9+iCnt, -d);
      WM.PutCoef(34, 21+iCnt, d);
   }   
   
   /* Derivate delle equazioni di vincolo di rotazione: e1b~e3a, e2b~e3a */
   Vec3 O1mO2(Omega1-Omega2);
   TmpPrime1 = e3a.Cross(O1mO2.Cross(e2b));   
   TmpPrime2 = e2b.Cross(e3a.Cross(O1mO2));
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WM.PutCoef(33, 3+iCnt, TmpPrime1.dGet(iCnt));
      WM.PutCoef(33, 15+iCnt, TmpPrime2.dGet(iCnt));
   }
   
   TmpPrime1 = e3a.Cross(O1mO2.Cross(e1b));
   TmpPrime2 = e1b.Cross(e3a.Cross(O1mO2));
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WM.PutCoef(34, 3+iCnt, TmpPrime1.dGet(iCnt));
      WM.PutCoef(34, 15+iCnt, TmpPrime2.dGet(iCnt));
   }   
   
   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
PlaneHingeJoint::InitialAssRes(SubVectorHandler& WorkVec,
			       const VectorHandler& XCurr)
{   
   DEBUGCOUT("Entering PlaneHingeJoint::InitialAssRes()" << std::endl);
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);

   /* Indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstVelIndex = iNode1FirstPosIndex+6;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstVelIndex = iNode2FirstPosIndex+6;
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+5;
   
   /* Setta gli indici */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {	
      WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WorkVec.PutRowIndex(6+iCnt, iNode1FirstVelIndex+iCnt);
      WorkVec.PutRowIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
      WorkVec.PutRowIndex(18+iCnt, iNode2FirstVelIndex+iCnt);
   }
   
   for (int iCnt = 1; iCnt <= 10; iCnt++) {
      WorkVec.PutRowIndex(24+iCnt, iFirstReactionIndex+iCnt);
   }

   /* Recupera i dati */
   const Vec3& x1(pNode1->GetXCurr());
   const Vec3& x2(pNode2->GetXCurr());
   const Vec3& v1(pNode1->GetVCurr());
   const Vec3& v2(pNode2->GetVCurr());
   const Mat3x3& R1(pNode1->GetRCurr());
   const Mat3x3& R2(pNode2->GetRCurr());
   const Vec3& Omega1(pNode1->GetWCurr());
   const Vec3& Omega2(pNode2->GetWCurr());
   
   /* Aggiorna F ed M, che restano anche per InitialAssJac */
   F = Vec3(XCurr, iFirstReactionIndex+1);
   M = Vec3(XCurr(iFirstReactionIndex+4),
	    XCurr(iFirstReactionIndex+5),
	    0.);
   Vec3 FPrime(XCurr, iReactionPrimeIndex+1);
   Vec3 MPrime(XCurr(iReactionPrimeIndex+4),
	       XCurr(iReactionPrimeIndex+5),
	       0.);
   
   /* Distanza nel sistema globale */
   Vec3 d1Tmp(R1*d1);
   Vec3 d2Tmp(R2*d2);
   Mat3x3 R1hTmp(R1*R1h);
   Mat3x3 R2hTmp(R2*R2h);

   /* Vettori omega1/\d1, -omega2/\d2 */
   Vec3 O1Wedged1(Omega1.Cross(d1Tmp));
   Vec3 O2Wedged2(Omega2.Cross(d2Tmp));
   
   Vec3 e3a(R1hTmp.GetVec(3));
   Vec3 e1b(R2hTmp.GetVec(1));
   Vec3 e2b(R2hTmp.GetVec(2));  
   
   /* Ruota il momento e la sua derivata con le matrici della cerniera 
    * rispetto ai nodi */
   Vec3 MTmp(e2b*M.dGet(1)-e1b*M.dGet(2));       
   Vec3 MPrimeTmp(e3a.Cross(MTmp.Cross(Omega2))+MTmp.Cross(Omega1.Cross(e3a))+
		  e2b.Cross(e3a)*MPrime.dGet(1)+e3a.Cross(e1b)*MPrime.dGet(2)); 
   
   /* Equazioni di equilibrio, nodo 1 */
   WorkVec.Sub(1, F);
   WorkVec.Add(4, F.Cross(d1Tmp)-MTmp.Cross(e3a)); /* Sfrutto il fatto che F/\d = -d/\F */
   
   /* Derivate delle equazioni di equilibrio, nodo 1 */
   WorkVec.Sub(7, FPrime);
   WorkVec.Add(10, FPrime.Cross(d1Tmp)-O1Wedged1.Cross(F)-MPrimeTmp);
   
   /* Equazioni di equilibrio, nodo 2 */
   WorkVec.Add(13, F);
   WorkVec.Add(16, d2Tmp.Cross(F)+MTmp.Cross(e3a)); 
   
   /* Derivate delle equazioni di equilibrio, nodo 2 */
   WorkVec.Add(19, FPrime);
   WorkVec.Add(22, d2Tmp.Cross(FPrime)+O2Wedged2.Cross(F)+MPrimeTmp);
   
   /* Equazione di vincolo di posizione */
   WorkVec.Add(25, x1+d1Tmp-x2-d2Tmp);
   
   /* Derivata dell'equazione di vincolo di posizione */
   WorkVec.Add(30, v1+O1Wedged1-v2-O2Wedged2);

   /* Equazioni di vincolo di rotazione */
   WorkVec.PutCoef(28, e2b.Dot(e3a));
   WorkVec.PutCoef(29, e1b.Dot(e3a));
   
   /* Derivate delle equazioni di vincolo di rotazione: e1b~e3a, e2b~e3a */
   Vec3 Tmp((Omega1-Omega2).Cross(e3a));
   WorkVec.PutCoef(33, e2b.Dot(Tmp));
   WorkVec.PutCoef(34, e1b.Dot(Tmp));

   return WorkVec;
}


unsigned int
PlaneHingeJoint::iGetNumPrivData(void) const
{
	return 8;
}

unsigned int
PlaneHingeJoint::iGetPrivDataIdx(const char *s) const
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
		idx += 3;
	case 'F':
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

doublereal PlaneHingeJoint::dGetPrivData(unsigned int i) const
{
   ASSERT(i >= 1 && i <= iGetNumPrivData());
   
   switch (i) {
    case 1: {
	Vec3 v(RotManip::VecRot((pNode1->GetRCurr()*R1h).MulTM(pNode2->GetRCurr()*R2h)));
	doublereal dThetaTmp(v(3));

	int n = 0;

	if (dThetaTmp - dThetaWrapped < -M_PI) {
		n++;
	}

	if (dThetaTmp - dThetaWrapped > M_PI) {
		n--;
	}

	return 2*M_PI*(NTheta + n) + dThetaTmp;
    }
      
    case 2: {
       Mat3x3 R2Tmp((pNode2->GetRCurr()*R2h));
       Vec3 v(R2Tmp.MulTV(pNode2->GetWCurr()-pNode1->GetWCurr()));
       
       return v(3);
    }

    case 3:
    case 4:
    case 5:
	    return F(i - 2);

    case 6:
    case 7:
    case 8:
	    return M(i - 5);
   }
      
   silent_cerr("PlaneHingeJoint(" << GetLabel() << "): "
	      "illegal private data " << i << std::endl);
   throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

/* PlaneHingeJoint - end */


/* PlaneRotationJoint - begin */

/* Costruttore non banale */
PlaneRotationJoint::PlaneRotationJoint(unsigned int uL, const DofOwner* pDO,
				 const StructNode* pN1, const StructNode* pN2,
				 const Mat3x3& R1hTmp, const Mat3x3& R2hTmp,
				 const OrientationDescription& od,
				 flag fOut)
: Elem(uL, fOut), 
Joint(uL, pDO, fOut), 
pNode1(pN1), pNode2(pN2),
R1h(R1hTmp), R2h(R2hTmp), M(Zero3),
NTheta(0), dTheta(0.),
dThetaWrapped(0.),
#ifdef USE_NETCDF
Var_Phi(0),
Var_Omega(0),
//Var_MFR(0),
//Var_MU(0),
#endif // USE_NETCDF
od(od)
{
   NO_OP;
}


/* Distruttore banale */
PlaneRotationJoint::~PlaneRotationJoint(void)
{
   NO_OP;
};


std::ostream&
PlaneRotationJoint::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << "->" << iIndex + 2 << ": "
			"reaction couples [mx,my]" << std::endl;

	if (bInitial) {
		iIndex += 2;
		out
			<< prefix << iIndex + 1 << "->" << iIndex + 2 << ": "
				"reaction couple derivatives [mPx,mPy]" << std::endl;
	}

	return out;
}

void
PlaneRotationJoint::DescribeDof(std::vector<std::string>& desc, bool bInitial, int i) const
{
	std::ostringstream os;
	os << "PlaneRotationJoint(" << GetLabel() << ")";

	unsigned short nself = 2;
	if (bInitial) {
		nself *= 2;
	}

	if (i == -1) {
		desc.resize(nself);
		std::string name = os.str();

		for (unsigned i = 0; i < 2; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": reaction couple m" << xyz[i];
			desc[i] = os.str();
		}

		if (bInitial) {
			for (unsigned i = 0; i < 2; i++) {
				os.str(name);
				os.seekp(0, std::ios_base::end);
				os << ": reaction couple derivative mP" << xyz[i];
				desc[2 + i] = os.str();
			}
		}

	} else {
		if (i < -1) {
			// error
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (i >= nself) {
			// error
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		desc.resize(1);

		switch (i) {
		case 0:
		case 1:
			os << ": reaction couple m" << xyz[i];
			break;

		case 2:
		case 3:
			os << ": reaction couple derivative mP" << xyz[i - 2];
			break;
		}
		desc[0] = os.str();
	}
}

std::ostream&
PlaneRotationJoint::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << "->" << iIndex + 2 << ": "
			"orientation constraints [gx1=gx2,gy1=gy2]" << std::endl;

	if (bInitial) {
		iIndex += 2;
		out
			<< prefix << iIndex + 1 << "->" << iIndex + 2 << ": "
				"angular velocity constraints [wx1=wx2,wy1=wy2]" << std::endl;
	}

	return out;
}

void
PlaneRotationJoint::DescribeEq(std::vector<std::string>& desc, bool bInitial, int i) const
{
	std::ostringstream os;
	os << "PlaneRotationJoint(" << GetLabel() << ")";

	unsigned short nself = 2;
	if (bInitial) {
		nself *= 2;
	}

	if (i == -1) {
		desc.resize(nself);
		std::string name = os.str();

		for (unsigned i = 0; i < 2; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": orientation constraint g" << xyz[i];
			desc[i] = os.str();
		}

		if (bInitial) {
			for (unsigned i = 0; i < 2; i++) {
				os.str(name);
				os.seekp(0, std::ios_base::end);
				os << ": orientation constraint derivative w" << xyz[i];
				desc[2 + i] = os.str();
			}
		}

	} else {
		if (i < -1) {
			// error
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (i >= nself) {
			// error
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		desc.resize(1);

		switch (i) {
		case 0:
		case 1:
			os << ": orientation constraint g" << xyz[i];
			break;

		case 2:
		case 3:
			os << ": orientation constraint derivative w" << xyz[i - 2];
			break;
		}
		desc[0] = os.str();
	}
}

void
PlaneRotationJoint::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
	if (ph) {
		for (unsigned i = 0; i < ph->size(); i++) {
			Joint::JointHint *pjh = dynamic_cast<Joint::JointHint *>((*ph)[i]);

			if (pjh == 0) {
				continue;
			}

			if (dynamic_cast<Joint::HingeHint<1> *>(pjh)) {
				R1h = pNode1->GetRCurr().MulTM(pNode2->GetRCurr()*R2h);

			} else if (dynamic_cast<Joint::HingeHint<2> *>(pjh)) {
				R2h = pNode2->GetRCurr().MulTM(pNode1->GetRCurr()*R1h);

			} else if (dynamic_cast<Joint::ReactionsHint *>(pjh)) {
				/* TODO */
			}
		}
	}

	Vec3 v(RotManip::VecRot((pNode1->GetRCurr()*R1h).MulTM(pNode2->GetRCurr()*R2h)));

	dThetaWrapped = dTheta = v.dGet(3);
}

Hint *
PlaneRotationJoint::ParseHint(DataManager *pDM, const char *s) const
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
PlaneRotationJoint::AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP)
{
	Vec3 v(RotManip::VecRot((pNode1->GetRCurr()*R1h).MulTM(pNode2->GetRCurr()*R2h)));
	doublereal dThetaTmp(v(3));

	// unwrap
	if (dThetaTmp - dThetaWrapped < -M_PI) {
		NTheta++;
	}

	if (dThetaTmp - dThetaWrapped > M_PI) {
		NTheta--;
	}

	// save new wrapped angle
	dThetaWrapped = dThetaTmp;

	// compute new unwrapped angle
	dTheta = 2*M_PI*NTheta + dThetaWrapped;
}


/* Contributo al file di restart */
std::ostream& PlaneRotationJoint::Restart(std::ostream& out) const
{
   Joint::Restart(out) << ", revolute hinge, "
     << pNode1->GetLabel() 
     << ", hinge, reference, node, 1, ", (R1h.GetVec(1)).Write(out, ", ")
     << ", 2, ", (R1h.GetVec(2)).Write(out, ", ") << ", "
     << pNode2->GetLabel() 
     << ", hinge, reference, node, 1, ", (R2h.GetVec(1)).Write(out, ", ")
     << ", 2, ", (R2h.GetVec(2)).Write(out, ", ") << ';' << std::endl;
   
   return out;
}


/* Assemblaggio jacobiano */
VariableSubMatrixHandler& 
PlaneRotationJoint::AssJac(VariableSubMatrixHandler& WorkMat,
			    doublereal dCoef,
			    const VectorHandler& /* XCurr */ ,
			    const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering PlaneRotationJoint::AssJac()" << std::endl);
   
   /* Setta la sottomatrice come piena (e' un po' dispersivo, ma lo jacobiano 
    * e' complicato */					
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Ridimensiona la sottomatrice in base alle esigenze */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeReset(iNumRows, iNumCols);
   
   /* Recupera gli indici delle varie incognite */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex()+3;
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex()+3;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex()+3;
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex()+3;
   integer iFirstReactionIndex = iGetFirstIndex();

   /* Setta gli indici delle equazioni */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WM.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutRowIndex(3+iCnt, iNode2FirstMomIndex+iCnt);
      WM.PutColIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
   }
   
   for (int iCnt = 1; iCnt <= 2; iCnt++) {	
      WM.PutRowIndex(6+iCnt, iFirstReactionIndex+iCnt);
      WM.PutColIndex(6+iCnt, iFirstReactionIndex+iCnt);
   }
   
   /* Recupera i dati che servono */
   const Mat3x3& R1(pNode1->GetRRef());
   const Mat3x3& R2(pNode2->GetRRef());   
   Mat3x3 R1hTmp(R1*R1h);
   Mat3x3 R2hTmp(R2*R2h);
   
   /* Suppongo che le reazioni F, M siano gia' state aggiornate da AssRes;
    * ricordo che la forza F e' nel sistema globale, mentre la coppia M
    * e' nel sistema locale ed il terzo termine, M(3), e' nullo in quanto
    * diretto come l'asse attorno al quale la rotazione e' consentita */
   
      
   /* 
    * La cerniera piana ha le prime 3 equazioni uguali alla cerniera sferica;
    * inoltre ha due equazioni che affermano la coincidenza dell'asse 3 del
    * riferimento solidale con la cerniera visto dai due nodi.
    * 
    *      (R1 * R1h * e1)^T * (R2 * R2h * e3) = 0
    *      (R1 * R1h * e2)^T * (R2 * R2h * e3) = 0
    * 
    * A queste equazioni corrisponde una reazione di coppia agente attorno 
    * agli assi 1 e 2 del riferimento della cerniera. La coppia attorno 
    * all'asse 3 e' nulla per definizione. Quindi la coppia nel sistema 
    * globale e':
    *      -R1 * R1h * M       per il nodo 1,
    *       R2 * R2h * M       per il nodo 2
    * 
    * 
    *       xa   ga                   xb   gb                     F     M 
    * Qa |  0    0                     0    0                     I     0  | | xa |   | -F           |
    * Ga |  0    c*(F/\da/\-(Sa*M)/\)  0    0                     da/\  Sa | | ga |   | -da/\F-Sa*M |
    * Qb |  0    0                     0    0                    -I     0  | | xb | = |  F           |
    * Gb |  0    0                     0   -c*(F/\db/\-(Sb*M)/\) -db/\ -Sb | | gb |   |  db/\F+Sb*M |
    * F  | -c*I  c*da/\                c*I -c*db/\                0     0  | | F  |   |  xa+da-xb-db |
    * M1 |  0    c*Tb1^T*Ta3/\         0    c*Ta3^T*Tb1/\         0     0  | | M  |   |  Sb^T*Ta3    |
    * M2 |  0    c*Tb2^T*Ta3/\         0    c*Ta3^T*Tb2/\         0     0  | 
    * 
    * con da = R1*d01, db = R2*d02, c = dCoef,
    * Sa = R1*R1h*[e1,e2], Sb = R2*R2h*[e1, e2],
    * Ta3 = R1*R1h*e3, Tb1 = R2*R2h*e1, Tb2 = R2*R2h*e2
    */

   /* Moltiplica il momento per il coefficiente del metodo */
   Vec3 MTmp = M*dCoef;

   Vec3 e3a(R1hTmp.GetVec(3));
   Vec3 e1b(R2hTmp.GetVec(1));
   Vec3 e2b(R2hTmp.GetVec(2));
   MTmp = e2b*MTmp.dGet(1)-e1b*MTmp.dGet(2);
   
   Mat3x3 MWedgee3aWedge(MatCrossCross, MTmp, e3a);
   Mat3x3 e3aWedgeMWedge(MatCrossCross, e3a, MTmp);
   
   WM.Sub(1, 1, MWedgee3aWedge);
   WM.Add(1, 4, e3aWedgeMWedge);
   
   WM.Add(4, 1, MWedgee3aWedge);   
   WM.Sub(4, 4, e3aWedgeMWedge);      
   
   /* Contributo del momento alle equazioni di equilibrio dei nodi */
   Vec3 Tmp1(e2b.Cross(e3a));
   Vec3 Tmp2(e3a.Cross(e1b));
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = Tmp1(iCnt);
      WM.PutCoef(iCnt, 7, d);
      WM.PutCoef(3+iCnt, 7, -d);
      d = Tmp2.dGet(iCnt);
      WM.PutCoef(iCnt, 8, d);
      WM.PutCoef(3+iCnt, 8, -d);
   }         
   
   /* Modifica: divido le equazioni di vincolo per dCoef */
   
   /* Equazione di vincolo del momento
    * 
    * Attenzione: bisogna scrivere il vettore trasposto
    *   (Sb[1]^T*(Sa[3]/\))*dCoef
    * Questo pero' e' uguale a:
    *   (-Sa[3]/\*Sb[1])^T*dCoef,
    * che puo' essere ulteriormente semplificato:
    *   (Sa[3].Cross(Sb[1])*(-dCoef))^T;
    */
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = Tmp1.dGet(iCnt);
      WM.PutCoef(7, iCnt, d);
      WM.PutCoef(7, 3+iCnt, -d);
      d = Tmp2.dGet(iCnt);
      WM.PutCoef(8, iCnt, -d);
      WM.PutCoef(8, 3+iCnt, d);
   }   
   
   return WorkMat;
}


/* Assemblaggio residuo */
SubVectorHandler& PlaneRotationJoint::AssRes(SubVectorHandler& WorkVec,
					  doublereal dCoef,
					  const VectorHandler& XCurr, 
					  const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering PlaneRotationJoint::AssRes()" << std::endl);
      
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);
 
   /* Indici */
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex()+3;
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex()+3;
   integer iFirstReactionIndex = iGetFirstIndex();
   
   /* Indici dei nodi */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.PutRowIndex(3+iCnt, iNode2FirstMomIndex+iCnt);
   }   
   
   /* Indici del vincolo */
   for (int iCnt = 1; iCnt <= 2; iCnt++) {
      WorkVec.PutRowIndex(6+iCnt, iFirstReactionIndex+iCnt);
   }

   /* Aggiorna i dati propri */
   M = Vec3(XCurr(iFirstReactionIndex+1),
	    XCurr(iFirstReactionIndex+2),
	    0.);

   /*
    * FIXME: provare a mettere "modificatori" di forza/momento sui gdl
    * residui: attrito, rigidezze e smorzamenti, ecc.
    */
   
   /* Recupera i dati */
   const Mat3x3& R1(pNode1->GetRCurr());
   const Mat3x3& R2(pNode2->GetRCurr());
   
   /* Costruisce i dati propri nella configurazione corrente */
   Mat3x3 R1hTmp(R1*R1h);
   Mat3x3 R2hTmp(R2*R2h);
   
   Vec3 e3a(R1hTmp.GetVec(3));
   Vec3 e1b(R2hTmp.GetVec(1));
   Vec3 e2b(R2hTmp.GetVec(2));
   
   Vec3 MTmp(e2b.Cross(e3a)*M.dGet(1)+e3a.Cross(e1b)*M.dGet(2));
   
   /* Equazioni di equilibrio, nodo 1 */
   WorkVec.Sub(1, MTmp); /* Sfrutto  F/\d = -d/\F */
   
   /* Equazioni di equilibrio, nodo 2 */
   WorkVec.Add(4, MTmp);

   /* Modifica: divido le equazioni di vincolo per dCoef */
   ASSERT(dCoef != 0.);
      
   /* Equazioni di vincolo di rotazione */
   Vec3 Tmp = R1hTmp.GetVec(3);
   WorkVec.PutCoef(7, Tmp.Dot(R2hTmp.GetVec(2))/dCoef);
   WorkVec.PutCoef(8, Tmp.Dot(R2hTmp.GetVec(1))/dCoef);
   
   return WorkVec;
}

DofOrder::Order
PlaneRotationJoint::GetEqType(unsigned int i) const
{
	ASSERTMSGBREAK(i < iGetNumDof(), 
		"INDEX ERROR in PlaneRotationJoint::GetEqType");
	return DofOrder::ALGEBRAIC;
}


void
PlaneRotationJoint::OutputPrepare(OutputHandler& OH)
{
	if (fToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			std::string name;
			OutputPrepare_int("revolute rotation", OH, name);

			Var_Phi = OH.CreateRotationVar(name, "", od, "global");

			Var_Omega = OH.CreateVar<Vec3>(name + "Omega", "radian/s",
				"local relative angular velocity (x, y, z)");
		}
#endif // USE_NETCDF
	}
}

/* Output (da mettere a punto) */
void PlaneRotationJoint::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) {
      Mat3x3 R2Tmp(pNode2->GetRCurr()*R2h);
      Mat3x3 RTmp((pNode1->GetRCurr()*R1h).MulTM(R2Tmp));
		Vec3 E;
		switch (od) {
		case EULER_123:
			E = MatR2EulerAngles123(RTmp)*dRaDegr;
			break;

		case EULER_313:
			E = MatR2EulerAngles313(RTmp)*dRaDegr;
			break;

		case EULER_321:
			E = MatR2EulerAngles321(RTmp)*dRaDegr;
			break;

		case ORIENTATION_VECTOR:
			E = RotManip::VecRot(RTmp);
			break;

		case ORIENTATION_MATRIX:
			break;

		default:
			/* impossible */
			break;
		}
      Vec3 OmegaTmp(R2Tmp.MulTV(pNode2->GetWCurr()-pNode1->GetWCurr()));
      
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			Var_F_local->put_rec(Zero3.pGetVec(), OH.GetCurrentStep());
			Var_M_local->put_rec(M.pGetVec(), OH.GetCurrentStep());
			Var_F_global->put_rec(Zero3.pGetVec(), OH.GetCurrentStep());
			Var_M_global->put_rec((R2Tmp*M).pGetVec(), OH.GetCurrentStep());
			switch (od) {
			case EULER_123:
			case EULER_313:
			case EULER_321:
			case ORIENTATION_VECTOR:
				Var_Phi->put_rec(E.pGetVec(), OH.GetCurrentStep());
				break;

			case ORIENTATION_MATRIX:
				Var_Phi->put_rec(RTmp.pGetMat(), OH.GetCurrentStep());
				break;

			default:
				/* impossible */
				break;
			}

			Var_Omega->put_rec(OmegaTmp.pGetVec(), OH.GetCurrentStep());
		}
#endif // USE_NETCDF
		if (OH.UseText(OutputHandler::JOINTS)) {
			  Joint::Output(OH.Joints(), "PlaneRotation", GetLabel(),
					Zero3, M, Zero3, R2Tmp*M)
			<< " ";

				switch (od) {
				case EULER_123:
				case EULER_313:
				case EULER_321:
				case ORIENTATION_VECTOR:
					OH.Joints() << E;
					break;

				case ORIENTATION_MATRIX:
					OH.Joints() << RTmp;
					break;

				default:
					/* impossible */
					break;
				}

				OH.Joints() << " " << OmegaTmp << std::endl;
		}
   }
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
PlaneRotationJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
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
   Vec3 MPrime(XCurr(iReactionPrimeIndex+1),
	       XCurr(iReactionPrimeIndex+2),
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

   Mat3x3 MDeltag1((Mat3x3(MatCross, Omega2.Cross(MTmp) + MPrimeTmp)
			   + Mat3x3(MatCrossCross, MTmp, Omega1))*Mat3x3(MatCross, e3a));
   Mat3x3 MDeltag2(Mat3x3(MatCrossCross, Omega1.Cross(e3a), MTmp)
		   + Mat3x3(MatCrossCross, e3a, MPrimeTmp)
		   + e3a.Cross(Mat3x3(MatCrossCross, Omega2, MTmp)));

   /* Vettori temporanei */
   Vec3 Tmp1(e2b.Cross(e3a));
   Vec3 Tmp2(e3a.Cross(e1b));
 
   /* Prodotto vettore tra il versore 3 della cerniera secondo il nodo 1
    * ed il versore 1 della cerniera secondo il nodo 2. A convergenza
    * devono essere ortogonali, quindi il loro prodotto vettore deve essere 
    * unitario */

   /* Error handling: il programma si ferma, pero' segnala dov'e' l'errore */
   if (Tmp1.Dot() <= std::numeric_limits<doublereal>::epsilon() || Tmp2.Dot() <= std::numeric_limits<doublereal>::epsilon()) {
      silent_cerr("PlaneRotationJoint(" << GetLabel() << "): "
	      "first and second node hinge axes are (nearly) orthogonal"
	      << std::endl);
      throw ErrNullNorm(MBDYN_EXCEPT_ARGS);
   }   
   
   Vec3 TmpPrime1(e2b.Cross(Omega1.Cross(e3a))-e3a.Cross(Omega2.Cross(e2b)));
   Vec3 TmpPrime2(e3a.Cross(Omega2.Cross(e1b))-e1b.Cross(Omega1.Cross(e3a)));
   
   /* Equazione di momento, nodo 1 */
   WM.Sub(4, 4, Mat3x3(MatCrossCross, MTmp, e3a));
   WM.Add(4, 16, Mat3x3(MatCrossCross, e3a, MTmp));
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutCoef(iCnt, 13, Tmp1.dGet(iCnt));
      WM.PutCoef(iCnt, 14, Tmp2.dGet(iCnt));	
   }
   
   /* Equazione di momento, nodo 2 */
   WM.Add(7, 1, Mat3x3(MatCrossCross, MTmp, e3a));
   WM.Sub(7, 7, Mat3x3(MatCrossCross, e3a, MTmp));
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutCoef(6+iCnt, 13, -Tmp1(iCnt));
      WM.PutCoef(6+iCnt, 14, -Tmp2(iCnt));	
   }
   
   /* Derivata dell'equazione di momento, nodo 1 */
   WM.Sub(4, 1, MDeltag1);
   WM.Sub(4, 4, Mat3x3(MatCrossCross, MTmp, e3a));
   WM.Add(4, 7, MDeltag2);
   WM.Add(4, 10, Mat3x3(MatCrossCross, e3a, MTmp));
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutCoef(3+iCnt, 13, TmpPrime1(iCnt));
      WM.PutCoef(3+iCnt, 14, TmpPrime2(iCnt));
      WM.PutCoef(3+iCnt, 15, Tmp1(iCnt));
      WM.PutCoef(3+iCnt, 16, Tmp2(iCnt));	
   }
   
   /* Derivata dell'equazione di momento, nodo 2 */
   WM.Add(10, 1, MDeltag1);
   WM.Add(10, 4, Mat3x3(MatCrossCross, MTmp, e3a));
   WM.Sub(10, 7, MDeltag2);
   WM.Sub(10, 10, Mat3x3(MatCrossCross, e3a, MTmp));

   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutCoef(9+iCnt, 13, -TmpPrime1(iCnt));
      WM.PutCoef(9+iCnt, 14, -TmpPrime2(iCnt));	
      WM.PutCoef(9+iCnt, 15, -Tmp1(iCnt));
      WM.PutCoef(9+iCnt, 16, -Tmp2(iCnt));	
   }
   
   /* Equazioni di vincolo di rotazione: e1b~e3a, e2b~e3a */
            
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      doublereal d = Tmp1(iCnt);
      WM.PutCoef(13, iCnt, d);
      WM.PutCoef(13, 6+iCnt, -d);
      
      /* Queste sono per la derivata dell'equazione, sono qui solo per 
       * ottimizzazione */
      WM.PutCoef(15, 3+iCnt, d);
      WM.PutCoef(15, 9+iCnt, -d);
      
      d = Tmp2(iCnt);
      WM.PutCoef(14, iCnt, -d);
      WM.PutCoef(14, 6+iCnt, d);
      
      /* Queste sono per la derivata dell'equazione, sono qui solo per 
       * ottimizzazione */
      WM.PutCoef(16, 3+iCnt, -d);
      WM.PutCoef(16, 9+iCnt, d);
   }   
   
   /* Derivate delle equazioni di vincolo di rotazione: e1b~e3a, e2b~e3a */
   Vec3 O1mO2(Omega1 - Omega2);
   TmpPrime1 = e3a.Cross(O1mO2.Cross(e2b));   
   TmpPrime2 = e2b.Cross(e3a.Cross(O1mO2));
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WM.PutCoef(15, iCnt, TmpPrime1(iCnt));
      WM.PutCoef(15, 6+iCnt, TmpPrime2(iCnt));
   }
   
   TmpPrime1 = e3a.Cross(O1mO2.Cross(e1b));
   TmpPrime2 = e1b.Cross(e3a.Cross(O1mO2));
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WM.PutCoef(16, iCnt, TmpPrime1(iCnt));
      WM.PutCoef(16, 6+iCnt, TmpPrime2(iCnt));
   }   
   
   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
PlaneRotationJoint::InitialAssRes(SubVectorHandler& WorkVec,
			       const VectorHandler& XCurr)
{   
   DEBUGCOUT("Entering PlaneRotationJoint::InitialAssRes()" << std::endl);
   
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
   M = Vec3(XCurr(iFirstReactionIndex+1),
	    XCurr(iFirstReactionIndex+2),
	    0.);
   Vec3 MPrime(XCurr(iReactionPrimeIndex+1),
	       XCurr(iReactionPrimeIndex+2),
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

   return WorkVec;
}


unsigned int
PlaneRotationJoint::iGetNumPrivData(void) const
{
	return 5;
}

unsigned int
PlaneRotationJoint::iGetPrivDataIdx(const char *s) const
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

doublereal PlaneRotationJoint::dGetPrivData(unsigned int i) const
{
   ASSERT(i >= 1 && i <= iGetNumPrivData());
   
   switch (i) {
    case 1: {
	Vec3 v(RotManip::VecRot((pNode1->GetRCurr()*R1h).MulTM(pNode2->GetRCurr()*R2h)));
	doublereal dThetaTmp(v(3));

	int n = 0;

	if (dThetaTmp - dThetaWrapped < -M_PI) {
		n++;
	}

	if (dThetaTmp - dThetaWrapped > M_PI) {
		n--;
	}

	return 2*M_PI*(NTheta + n) + dThetaTmp;
    }
      
    case 2: {
       Mat3x3 R2Tmp((pNode2->GetRCurr()*R2h));
       Vec3 v(R2Tmp.MulTV(pNode2->GetWCurr()-pNode1->GetWCurr()));
       
       return v(3);
    }

    case 3:
    case 4:
    case 5:
	    return M(i - 2);
   }
      
   silent_cerr("PlaneRotationJoint(" << GetLabel() << "): "
	   "illegal private data " << i << std::endl);
   throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

/* PlaneRotationJoint - end */


/* AxialRotationJoint - begin */

const unsigned int AxialRotationJoint::NumSelfDof(6);
const unsigned int AxialRotationJoint::NumDof(18);

/* Costruttore non banale */
AxialRotationJoint::AxialRotationJoint(unsigned int uL, const DofOwner* pDO,
		const StructNode* pN1, 
		const StructNode* pN2,
		const Vec3& dTmp1, const Vec3& dTmp2,
		const Mat3x3& R1hTmp, 
		const Mat3x3& R2hTmp,
		const DriveCaller* pDC,
		const OrientationDescription& od,
		flag fOut,
		const doublereal rr,
		const doublereal pref,
		BasicShapeCoefficient *const sh,
		BasicFriction *const f)
: Elem(uL, fOut), 
Joint(uL, pDO, fOut), 
DriveOwner(pDC), 
pNode1(pN1), pNode2(pN2), 
d1(dTmp1), R1h(R1hTmp), d2(dTmp2), R2h(R2hTmp), F(Zero3), M(Zero3),
NTheta(0), dTheta(0.), dThetaWrapped(0.),
#ifdef USE_NETCDF
Var_Phi(0),
Var_Omega(0),
//Var_MFR(0),
//Var_MU(0),
#endif // USE_NETCDF
Sh_c(sh), fc(f), preF(pref), r(rr),
od(od)
{
	NO_OP;
}


/* Distruttore banale */
AxialRotationJoint::~AxialRotationJoint(void)
{
	if (Sh_c) {
		delete Sh_c;
	}

	if (fc) {
		delete fc;
	}
}


std::ostream&
AxialRotationJoint::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
			"reaction forces [Fx,Fy,Fz]" << std::endl
		<< prefix << iIndex + 4 << "->" << iIndex + 6 << ": "
			"reaction couples [mx,my,mz]" << std::endl;

	if (bInitial) {
		iIndex += NumSelfDof;
		out
			<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
				"reaction force derivatives [FPx,FPy,FPz]" << std::endl
			<< prefix << iIndex + 4 << "->" << iIndex + 5 << ": "
				"reaction couple derivatives [mPx,mPy]" << std::endl;
	}
	
	iIndex += NumSelfDof;
	if (fc) {
		integer iFCDofs = fc->iGetNumDof();
		if (iFCDofs > 0) {
			out << prefix << iIndex + 1;
			if (iFCDofs > 1) {
				out << "->" << iIndex + iFCDofs;
			}
			out << ": friction dof(s)" << std::endl
				<< "        ", fc->DescribeDof(out, prefix, bInitial);
		}
	}

	return out;
}

void
AxialRotationJoint::DescribeDof(std::vector<std::string>& desc, bool bInitial, int i) const
{
	std::ostringstream os;
	os << "AxialRotationJoint(" << GetLabel() << ")";

	unsigned short nself = NumSelfDof;
	if (bInitial) {
		nself += NumSelfDof - 1;
	}
	if (fc && (i == -1 || i >= nself)) {
		fc->DescribeDof(desc, bInitial, i - nself);
		if (i != -1) {
			desc[0] = os.str() + ": " + desc[0];
			return;
		}
	}

	if (i == -1) {
		// move fc desc to the end
		unsigned short nfc = 0;
		if (fc) {
			nfc = desc.size();
		}
		desc.resize(nfc + nself);
		for (unsigned i = nfc; i-- > 0; ) {
			desc[nself + i] = os.str() + ": " + desc[nfc];
		}

		std::string name = os.str();

		for (unsigned i = 0; i < 3; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": reaction force f" << xyz[i];
			desc[i] = os.str();
		}

		for (unsigned i = 0; i < 3; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": reaction couple m" << xyz[i];
			desc[3 + i] = os.str();
		}

		if (bInitial) {
			for (unsigned i = 0; i < 3; i++) {
				os.str(name);
				os.seekp(0, std::ios_base::end);
				os << ": reaction force derivative fP" << xyz[i];
				desc[6 + i] = os.str();
			}
	
			for (unsigned i = 0; i < 2; i++) {
				os.str(name);
				os.seekp(0, std::ios_base::end);
				os << ": reaction couple derivative mP" << xyz[i];
				desc[9 + i] = os.str();
			}
		}

	} else {
		if (i < -1) {
			// error
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (i >= nself) {
			// error
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		desc.resize(1);

		switch (i) {
		case 0:
		case 1:
		case 2:
			os << ": reaction force f" << xyz[i];
			break;

		case 3:
		case 4:
		case 5:
			os << ": reaction couple m" << xyz[i - 3];
			break;

		case 6:
		case 7:
		case 8:
			os << ": reaction force derivative fP" << xyz[i - 6];
			break;

		case 9:
		case 10:
			os << ": reaction couple derivative mP" << xyz[i - 9];
			break;
		}
		desc[0] = os.str();
	}
}

std::ostream&
AxialRotationJoint::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
			"position constraints [Px1=Px2,Py1=Py2,Pz1=Pz2]" << std::endl
		<< prefix << iIndex + 4 << "->" << iIndex + 5 << ": "
			"orientation constraints [gx1=gx2,gy1=gy2]" << std::endl
		<< prefix << iIndex + 6 << ": "
			"angular velocity constraint wz2-wz1=Omega]" << std::endl;

	if (bInitial) {
		iIndex += NumSelfDof;
		out
			<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
				"velocity constraints [vx1=vx2,vy1=vy2,vz1=vz2]" << std::endl
			<< prefix << iIndex + 4 << "->" << iIndex + 5 << ": "
				"reaction couple derivatives [wx1=wx2,wy1=wy2]" << std::endl;
	}
	
	iIndex += NumSelfDof;
	if (fc) {
		integer iFCDofs = fc->iGetNumDof();
		if (iFCDofs > 0) {
			out << prefix << iIndex + 1;
			if (iFCDofs > 1) {
				out << "->" << iIndex + iFCDofs;
			}
			out << ": friction equation(s)" << std::endl
				<< "        ", fc->DescribeEq(out, prefix, bInitial);
		}
	}

	return out;
}

void
AxialRotationJoint::DescribeEq(std::vector<std::string>& desc, bool bInitial, int i) const
{
	std::ostringstream os;
	os << "AxialRotationJoint(" << GetLabel() << ")";

	unsigned short nself = NumSelfDof;
	if (bInitial) {
		nself += NumSelfDof - 1;
	}
	if (fc && (i == -1 || i >= nself)) {
		fc->DescribeEq(desc, bInitial, i - nself);
		if (i != -1) {
			desc[0] = os.str() + ": " + desc[0];
			return;
		}
	}

	if (i == -1) {
		// move fc desc to the end
		unsigned short nfc = 0;
		if (fc) {
			nfc = desc.size();
		}
		desc.resize(nfc + nself);
		for (unsigned i = nfc; i-- > 0; ) {
			desc[nself + i] = os.str() + ": " + desc[nfc];
		}

		std::string name = os.str();

		for (unsigned i = 0; i < 3; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": position constraint P" << xyz[i];
			desc[i] = os.str();
		}

		for (unsigned i = 0; i < 2; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": orientation constraint g" << xyz[i];
			desc[3 + i] = os.str();
		}

		os.str(name);
		os.seekp(0, std::ios_base::end);
		os << ": angular velocity constraint wz";
		desc[5] = os.str();

		if (bInitial) {
			for (unsigned i = 0; i < 3; i++) {
				os.str(name);
				os.seekp(0, std::ios_base::end);
				os << ": position constraint derivative v" << xyz[i];
				desc[6 + i] = os.str();
			}
	
			for (unsigned i = 0; i < 2; i++) {
				os.str(name);
				os.seekp(0, std::ios_base::end);
				os << ": orientation constraint derivative w" << xyz[i];
				desc[9 + i] = os.str();
			}
		}

	} else {
		if (i < -1) {
			// error
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (i >= nself) {
			// error
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		desc.resize(1);

		switch (i) {
		case 0:
		case 1:
		case 2:
			os << ": position constraint P" << xyz[i];
			break;

		case 3:
		case 4:
			os << ": orientation constraint g" << xyz[i - 3];
			break;

		case 5:
			os << ": angular velocity constraint wz";
			break;

		case 6:
		case 7:
		case 8:
			os << ": position constraint derivative v" << xyz[i - 6];
			break;

		case 9:
		case 10:
			os << ": orientation constraint derivative w" << xyz[i - 9];
			break;
		}
		desc[0] = os.str();
	}
}

void
AxialRotationJoint::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
	if (ph) {
		for (unsigned i = 0; i < ph->size(); i++) {
			Joint::JointHint *pjh = dynamic_cast<Joint::JointHint *>((*ph)[i]);

			if (pjh == 0) {
				DriveHint *pdh = dynamic_cast<DriveHint *>((*ph)[i]);
				if (pdh) {
					DriveCaller *pDC = pdh->pCreateDrive(pDM);
					if (pDC == 0) {
						silent_cerr("AxialRotationJoint::SetValue: "
							"unable to create drive after hint "
							"#" << i << std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}

					DriveOwner::Set(pDC);
				}
				continue;
			}

			if (dynamic_cast<Joint::OffsetHint<1> *>(pjh)) {
				const Mat3x3& R1(pNode1->GetRCurr());
				Vec3 dTmp2(pNode2->GetRCurr()*d2);
   
				d1 = R1.MulTV(pNode2->GetXCurr() + dTmp2 - pNode1->GetXCurr());

			} else if (dynamic_cast<Joint::OffsetHint<2> *>(pjh)) {
				const Mat3x3& R2(pNode2->GetRCurr());
				Vec3 dTmp1(pNode1->GetRCurr()*d1);
   
				d2 = R2.MulTV(pNode1->GetXCurr() + dTmp1 - pNode2->GetXCurr());

			} else if (dynamic_cast<Joint::HingeHint<1> *>(pjh)) {
				R1h = pNode1->GetRCurr().MulTM(pNode2->GetRCurr()*R2h);

			} else if (dynamic_cast<Joint::HingeHint<2> *>(pjh)) {
				R2h = pNode2->GetRCurr().MulTM(pNode1->GetRCurr()*R1h);

			} else if (dynamic_cast<Joint::ReactionsHint *>(pjh)) {
				/* TODO */
			}
		}
	}

	Mat3x3 RTmp((pNode1->GetRCurr()*R1h).MulTM(pNode2->GetRCurr()*R2h));
	Vec3 v(MatR2EulerAngles(RTmp));

	dThetaWrapped = dTheta = v.dGet(3);
	
	if (fc) {
		fc->SetValue(pDM, X, XP, ph, iGetFirstIndex() + NumSelfDof);
	}
}

Hint *
AxialRotationJoint::ParseHint(DataManager *pDM, const char *s) const
{
	if (strncasecmp(s, "offset{" /*}*/, STRLENOF("offset{" /*}*/)) == 0) {
		s += STRLENOF("offset{" /*}*/);

		if (strcmp(&s[1], /*{*/ "}") != 0) {
			return 0;
		}

		switch (s[0]) {
		case '1':
			return new Joint::OffsetHint<1>;

		case '2':
			return new Joint::OffsetHint<2>;
		}

	} else if (strncasecmp(s, "hinge{" /*}*/, STRLENOF("hinge{" /*}*/)) == 0) {
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

	} else if (fc) {
		return fc->ParseHint(pDM, s);
	}

	return 0;
}

void
AxialRotationJoint::AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP)
{
	Vec3 v(RotManip::VecRot((pNode1->GetRCurr()*R1h).MulTM(pNode2->GetRCurr()*R2h)));
	doublereal dThetaTmp(v(3));

	// unwrap
	if (dThetaTmp - dThetaWrapped < -M_PI) {
		NTheta++;
	}

	if (dThetaTmp - dThetaWrapped > M_PI) {
		NTheta--;
	}

	// save new wrapped angle
	dThetaWrapped = dThetaTmp;

	// compute new unwrapped angle
	dTheta = 2*M_PI*NTheta + dThetaWrapped;

	if (fc) {
		const Mat3x3& R1(pNode1->GetRCurr());
		Mat3x3 R1hTmp(R1*R1h);
		Vec3 e3a(R1hTmp.GetVec(3));
		const Vec3& Omega1(pNode1->GetWCurr());
		const Vec3& Omega2(pNode2->GetWCurr());
		//relative velocity
		doublereal v = (Omega1-Omega2).Dot(e3a)*r;
		//reaction norm
		doublereal modF = std::max(F.Norm(), preF);;
		fc->AfterConvergence(modF,v,X,XP,iGetFirstIndex()+NumSelfDof);
	}
}


/* Contributo al file di restart */
std::ostream& AxialRotationJoint::Restart(std::ostream& out) const
{
   Joint::Restart(out) << ", axial rotation, "
     << pNode1->GetLabel() 
     << ", reference, node, ", d1.Write(out, ", ")
     << ", hinge, reference, node, 1, ", (R1h.GetVec(1)).Write(out, ", ")
     << ", 2, ", (R1h.GetVec(2)).Write(out, ", ") << ", "
     << pNode2->GetLabel() 
     << ", reference, node, ", d2.Write(out, ", ") 
     << ", hinge, reference, node, 1, ", (R2h.GetVec(1)).Write(out, ", ")
     << ", 2, ", (R2h.GetVec(2)).Write(out, ", ") << ", ";
   
   return pGetDriveCaller()->Restart(out) << ';' << std::endl;
}


/* Assemblaggio jacobiano */
VariableSubMatrixHandler& 
AxialRotationJoint::AssJac(VariableSubMatrixHandler& WorkMat,
			   doublereal dCoef,
			   const VectorHandler& XCurr,
			   const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering AxialRotationJoint::AssJac()" << std::endl);
   
   /* Setta la sottomatrice come piena (e' un po' dispersivo, ma lo jacobiano 
    * e' complicato */					
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
   integer iFirstReactionIndex = iGetFirstIndex();

   /* Recupera i dati che servono */
   const Mat3x3& R1(pNode1->GetRRef());
   const Mat3x3& R2(pNode2->GetRRef());   
   const Vec3& Omega2(pNode2->GetWRef()); /* Serve dopo */
   Vec3 d1Tmp(R1*d1);
   Vec3 d2Tmp(R2*d2);
   Mat3x3 R1hTmp(R1*R1h);
   Mat3x3 R2hTmp(R2*R2h);
   /* Suppongo che le reazioni F, M siano gia' state aggiornate da AssRes;
    * ricordo che la forza F e' nel sistema globale, mentre la coppia M
    * e' nel sistema locale ed il terzo termine, M(3), e' nullo in quanto
    * diretto come l'asse attorno al quale la rotazione e' consentita */
   
      
   /* 
    * La cerniera piana ha le prime 3 equazioni uguali alla cerniera sferica;
    * inoltre ha due equazioni che affermano la coincidenza dell'asse 3 del
    * riferimento solidale con la cerniera visto dai due nodi.
    * 
    *      (R1 * R1h * e1)^T * (R2 * R2h * e3) = 0
    *      (R1 * R1h * e2)^T * (R2 * R2h * e3) = 0
    * 
    * A queste equazioni corrisponde una reazione di coppia agente attorno 
    * agli assi 1 e 2 del riferimento della cerniera. La coppia attorno 
    * all'asse 3 e' nulla per definizione. Quindi la coppia nel sistema 
    * globale e':
    *      -R1 * R1h * M       per il nodo 1,
    *       R2 * R2h * M       per il nodo 2
    * 
    * 
    *       xa   ga                   xb   gb                     F     M 
    * Qa |  0    0                     0    0                     I     0  | | xa |   | -F           |
    * Ga |  0    c*(F/\da/\-(Sa*M)/\)  0    0                     da/\  Sa | | ga |   | -da/\F-Sa*M |
    * Qb |  0    0                     0    0                    -I     0  | | xb | = |  F           |
    * Gb |  0    0                     0   -c*(F/\db/\-(Sb*M)/\) -db/\ -Sb | | gb |   |  db/\F+Sb*M |
    * F  | -c*I  c*da/\                c*I -c*db/\                0     0  | | F  |   |  xa+da-xb-db |
    * M1 |  0    c*Tb1^T*Ta3/\         0    c*Ta3^T*Tb1/\         0     0  | | M  |   |  Sb^T*Ta3    |
    * M2 |  0    c*Tb2^T*Ta3/\         0    c*Ta3^T*Tb2/\         0     0  | 
    * 
    * con da = R1*d01, db = R2*d02, c = dCoef,
    * Sa = R1*R1h*[e1,e2], Sb = R2*R2h*[e1, e2],
    * Ta3 = R1*R1h*e3, Tb1 = R2*R2h*e1, Tb2 = R2*R2h*e2
    */

   /* Setta gli indici delle equazioni */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
      WM.PutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
   }
   
   for (unsigned int iCnt = 1; iCnt <= iGetNumDof(); iCnt++) {
      WM.PutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
      WM.PutColIndex(12+iCnt, iFirstReactionIndex+iCnt);
   }
   
   /* Contributo della forza alle equazioni di equilibrio dei due nodi */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WM.PutCoef(iCnt, 12+iCnt, 1.);
      WM.PutCoef(6+iCnt, 12+iCnt, -1.);
   }
   
   WM.Add(4, 13, Mat3x3(MatCross, d1Tmp));
   WM.Sub(10, 13, Mat3x3(MatCross, d2Tmp));   
   
   /* Moltiplica la forza ed il momento per il coefficiente
    * del metodo */
   Vec3 FTmp = F*dCoef;
   Vec3 MTmp = M*dCoef;
   
   Vec3 e3a(R1hTmp.GetVec(3));
   Vec3 e1b(R2hTmp.GetVec(1));
   Vec3 e2b(R2hTmp.GetVec(2));
   MTmp = e2b*MTmp.dGet(1) - e1b*MTmp.dGet(2);
   
   Mat3x3 MWedgee3aWedge(MatCrossCross, MTmp, e3a);
   Mat3x3 e3aWedgeMWedge(MatCrossCross, e3a, MTmp);
   
   WM.Add(4, 4, Mat3x3(MatCrossCross, FTmp, d1Tmp) - MWedgee3aWedge - Mat3x3(MatCross, e3a*M(3)));
   WM.Add(4, 10, e3aWedgeMWedge);
   
   WM.Add(10, 4, MWedgee3aWedge);
   WM.Sub(10, 10, Mat3x3(MatCrossCross, FTmp, d2Tmp) + e3aWedgeMWedge
		   - Mat3x3(MatCross, e3a*M(3)));
 
   /* Contributo del momento alle equazioni di equilibrio dei nodi */
   Vec3 Tmp1(e2b.Cross(e3a));
   Vec3 Tmp2(e3a.Cross(e1b));
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = Tmp1(iCnt);
      WM.PutCoef(3+iCnt, 16, d);
      WM.PutCoef(9+iCnt, 16, -d);
      
      d = Tmp2(iCnt);
      WM.PutCoef(3+iCnt, 17, d);
      WM.PutCoef(9+iCnt, 17, -d);
      
      d = e3a(iCnt);
      WM.PutCoef(3+iCnt, 18, d);	
      WM.PutCoef(9+iCnt, 18, -d);	
   }         
   
   /* Modifica: divido le equazioni di vincolo per dCoef */
   
   /* Equazioni di vincolo degli spostamenti */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutCoef(12+iCnt, iCnt, -1.);
      WM.PutCoef(12+iCnt, 6+iCnt, 1.);
   }
   
   WM.Add(13, 4, Mat3x3(MatCross, d1Tmp));
   WM.Sub(13, 10, Mat3x3(MatCross, d2Tmp));
   
   /* Equazione di vincolo del momento
    * 
    * Attenzione: bisogna scrivere il vettore trasposto
    *   (Sb[1]^T*(Sa[3]/\))*dCoef
    * Questo pero' e' uguale a:
    *   (-Sa[3]/\*Sb[1])^T*dCoef,
    * che puo' essere ulteriormente semplificato:
    *   (Sa[3].Cross(Sb[1])*(-dCoef))^T;
    */
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = Tmp1(iCnt);
      WM.PutCoef(16, 3+iCnt, d);
      WM.PutCoef(16, 9+iCnt, -d);
      
      d = Tmp2.dGet(iCnt);
      WM.PutCoef(17, 3+iCnt, -d);
      WM.PutCoef(17, 9+iCnt, d);
   }
   
   /* Questa equazione non viene divisa per dCoef */
   
   /* Equazione di imposizione della velocita' di rotazione: 
    * (e3a+c(wb/\e3a))^T*(Delta_gPb-Delta_gPa) */
   Vec3 Tmp = e3a + Omega2.Cross(e3a*dCoef);
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = Tmp(iCnt);
      WM.PutCoef(18, 3+iCnt, -d);
      WM.PutCoef(18, 9+iCnt, d);
   }   

   if (fc) {
      //retrive
          //friction coef
      doublereal f = fc->fc();
          //shape function
      doublereal shc = Sh_c->Sh_c();
          //omega and omega rif
      const Vec3& Omega1(pNode1->GetWCurr());
      const Vec3& Omega2(pNode2->GetWCurr());
      // const Vec3& Omega1r(pNode1->GetWRef());
      // const Vec3& Omega2r(pNode2->GetWRef());   
      //compute 
          //relative velocity
      doublereal v = (Omega1-Omega2).Dot(e3a)*r;
          //reaction norm
      doublereal modF = std::max(F.Norm(), preF);
          //reaction moment
      //doublereal M3 = shc*modF*f*r;
      
      ExpandableRowVector dfc;
      ExpandableRowVector dF;
      ExpandableRowVector dv;
          //variation of reaction force
      dF.ReDim(3);
      if ((modF == 0.) or (F.Norm() > preF)) {
          dF.Set(0.,1,12+1);
          dF.Set(0.,2,12+2);
          dF.Set(0.,3,12+3);
      } else {
          dF.Set(F.dGet(1)/modF,1,12+1);
          dF.Set(F.dGet(2)/modF,2,12+2);
          dF.Set(F.dGet(3)/modF,3,12+3);
      }
          //variation of relative velocity
      dv.ReDim(6);
      
/* (approximate: assume constant triads orientations) 
 * relative velocity linearization 
*/
      dv.Set(e3a*r,1, 0+4);
      dv.Set(-e3a*r,4, 6+4);

      //assemble friction states
      fc->AssJac(WM,dfc,12+NumSelfDof,iFirstReactionIndex+NumSelfDof,dCoef,modF,v,
      		XCurr,XPrimeCurr,dF,dv);
      ExpandableRowVector dM3;
      ExpandableRowVector dShc;
      //compute 
          //variation of shape function
      Sh_c->dSh_c(dShc,f,modF,v,dfc,dF,dv);
          //variation of moment component
      dM3.ReDim(2);
      dM3.Set(shc * r,1); dM3.Link(1,&dF);
      dM3.Set(modF * r,2); dM3.Link(2,&dShc);
      //assemble first node
          //variation of moment component
      dM3.Add(WM,0+4,e3a.dGet(1));
      dM3.Add(WM,0+5,e3a.dGet(2));
      dM3.Add(WM,0+6,e3a.dGet(3));
      //assemble second node
          //variation of moment component
      dM3.Sub(WM,6+4,e3a.dGet(1));
      dM3.Sub(WM,6+5,e3a.dGet(2));
      dM3.Sub(WM,6+6,e3a.dGet(3));
   }
   
   return WorkMat;
}


/* Assemblaggio residuo */
SubVectorHandler& AxialRotationJoint::AssRes(SubVectorHandler& WorkVec,
					     doublereal dCoef,
					     const VectorHandler& XCurr,
					     const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering AxialRotationJoint::AssRes()" << std::endl);

   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);
 
   /* Indici */
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   
   /* Indici dei nodi */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {	
      WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.PutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
   }
   /* Indici del vincolo */
   for (unsigned int iCnt = 1; iCnt <= iGetNumDof(); iCnt++) {
      WorkVec.PutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
   }
   
   /* Aggiorna i dati propri */
   F = Vec3(XCurr, iFirstReactionIndex+1);
   M = Vec3(XCurr, iFirstReactionIndex+4);
   
   /* Recupera i dati */
   const Vec3& x1(pNode1->GetXCurr());
   const Vec3& x2(pNode2->GetXCurr());
   const Mat3x3& R1(pNode1->GetRCurr());
   const Mat3x3& R2(pNode2->GetRCurr());
   
   /* Costruisce i dati propri nella configurazione corrente */
   Vec3 dTmp1(R1*d1);
   Vec3 dTmp2(R2*d2);
   Mat3x3 R1hTmp(R1*R1h);
   Mat3x3 R2hTmp(R2*R2h);
   
   Vec3 e3a(R1hTmp.GetVec(3));
   Vec3 e1b(R2hTmp.GetVec(1));
   Vec3 e2b(R2hTmp.GetVec(2));
   
   Vec3 MTmp(e2b.Cross(e3a)*M.dGet(1)+e3a.Cross(e1b)*M.dGet(2));
   
   /* Equazioni di equilibrio, nodo 1 */
   WorkVec.Sub(1, F);
   WorkVec.Add(4, F.Cross(dTmp1)-MTmp-e3a*M.dGet(3)); /* Sfrutto  F/\d = -d/\F */
   
   /* Equazioni di equilibrio, nodo 2 */
   WorkVec.Add(7, F);
   WorkVec.Add(10, dTmp2.Cross(F)+MTmp+e3a*M.dGet(3));

   /* Modifica: divido le equazioni di vincolo per dCoef */
   ASSERT(dCoef != 0.);
      
   /* Equazione di vincolo di posizione */
   WorkVec.Add(13, (x1+dTmp1-x2-dTmp2)/dCoef);
      
   /* Equazioni di vincolo di rotazione */
   Vec3 Tmp = Vec3(R1hTmp.GetVec(3));
   WorkVec.PutCoef(16, Tmp.Dot(R2hTmp.GetVec(2))/dCoef);
   WorkVec.PutCoef(17, Tmp.Dot(R2hTmp.GetVec(1))/dCoef);
   
   /* Questa equazione non viene divisa per dCoef */
   
   /* Equazione di vincolo di velocita' di rotazione */
   const Vec3& Omega1(pNode1->GetWCurr());
   const Vec3& Omega2(pNode2->GetWCurr());
   doublereal dOmega0 = pGetDriveCaller()->dGet();
   WorkVec.PutCoef(18, dOmega0-e3a.Dot(Omega2-Omega1));

   if (fc) {
      bool ChangeJac(false);
      doublereal v = (Omega1-Omega2).Dot(e3a)*r;
      doublereal modF = std::max(F.Norm(), preF);
      try {
          fc->AssRes(WorkVec,12+NumSelfDof,iFirstReactionIndex+NumSelfDof,modF,v,XCurr,XPrimeCurr);
      }
      catch (Elem::ChangedEquationStructure) {
          ChangeJac = true;
      }
      doublereal f = fc->fc();
      doublereal shc = Sh_c->Sh_c(f,modF,v);
      M3 = shc*modF*r;
      WorkVec.Sub(4,e3a*M3);
      WorkVec.Add(10,e3a*M3);
//!!!!!!!!!!!!!!
//      M += e3a*M3;
      if (ChangeJac) {
          throw Elem::ChangedEquationStructure(MBDYN_EXCEPT_ARGS);
      }
   }
   
   return WorkVec;
}

DofOrder::Order
AxialRotationJoint::GetEqType(unsigned int i) const
{
	ASSERTMSGBREAK(i < iGetNumDof(),
		"INDEX ERROR in AxialRotationJoint::GetEqType");
	if (i == 5) {
		return DofOrder::DIFFERENTIAL;
	} else if (i < NumSelfDof) {
		return DofOrder::ALGEBRAIC;
	} else {
		return fc->GetEqType(i-NumSelfDof);
	}
}

void
AxialRotationJoint::OutputPrepare(OutputHandler& OH)
{
	if (fToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			std::string name;
			OutputPrepare_int("axial rotation", OH, name);

			Var_Phi = OH.CreateRotationVar(name, "", od, "global");

			Var_Omega = OH.CreateVar<Vec3>(name + "Omega", "radian/s",
				"local relative angular velocity (x, y, z)");

/* TODO
			Var_MFR = OH.CreateVar<doublereal>(name + "MFR", "Nm",
				"friciton moment ");

			Var_MU = OH.CreateVar<doublereal>(name + "MU", "--",
					"friction model specific data: friction coefficient?");
*/
		}
#endif // USE_NETCDF
	}
}

/* Output (da mettere a punto) */
void AxialRotationJoint::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) {
      Mat3x3 R2Tmp(pNode2->GetRCurr()*R2h);
      Mat3x3 RTmp((pNode1->GetRCurr()*R1h).MulTM(R2Tmp));
		Vec3 E;
		switch (od) {
		case EULER_123:
			E = MatR2EulerAngles123(RTmp)*dRaDegr;
			break;

		case EULER_313:
			E = MatR2EulerAngles313(RTmp)*dRaDegr;
			break;

		case EULER_321:
			E = MatR2EulerAngles321(RTmp)*dRaDegr;
			break;

		case ORIENTATION_VECTOR:
			E = RotManip::VecRot(RTmp);
			break;

		case ORIENTATION_MATRIX:
			break;

		default:
			/* impossible */
			break;
		}
      Vec3 OmegaTmp(R2Tmp.MulTV(pNode2->GetWCurr()-pNode1->GetWCurr()));
      
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			Var_F_local->put_rec((R2Tmp.MulTV(F)).pGetVec(), OH.GetCurrentStep());
			Var_M_local->put_rec(M.pGetVec(), OH.GetCurrentStep());
			Var_F_global->put_rec(F.pGetVec(), OH.GetCurrentStep());
			Var_M_global->put_rec((R2Tmp*M).pGetVec(), OH.GetCurrentStep());
			switch (od) {
			case EULER_123:
			case EULER_313:
			case EULER_321:
			case ORIENTATION_VECTOR:
				Var_Phi->put_rec(E.pGetVec(), OH.GetCurrentStep());
				break;

			case ORIENTATION_MATRIX:
				Var_Phi->put_rec(RTmp.pGetMat(), OH.GetCurrentStep());
				break;

			default:
				/* impossible */
				break;
			}
			Var_Omega->put_rec(OmegaTmp.pGetVec(), OH.GetCurrentStep());
/*
			if (fc) {
					Var_MFR->put_rec(&M3, OH.GetCurrentStep());
					Var_MU->put_rec(fc->fc(), OH.GetCurrentStep());
			}
			else
			{
				Var_MFR->put_rec(0, OH.GetCurrentStep());
				Var_MU->put_rec(0, OH.GetCurrentStep());
			}
*/
		}
#endif // USE_NETCDF
		if (OH.UseText(OutputHandler::JOINTS)) {
		  std::ostream &of = Joint::Output(OH.Joints(), "AxialRotation", GetLabel(),
				R2Tmp.MulTV(F), M, F, R2Tmp*M)
		  << " ";

			switch (od) {
			case EULER_123:
			case EULER_313:
			case EULER_321:
			case ORIENTATION_VECTOR:
				OH.Joints() << E;
				break;

			case ORIENTATION_MATRIX:
				OH.Joints() << RTmp;
				break;

			default:
				/* impossible */
				break;
			}

			OH.Joints() << " " << dGet()
		  << " " << OmegaTmp;
		  if (fc) {
			  of << " " << M3 << " " << fc->fc();
		  }
		of << std::endl;
		}
   }
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
AxialRotationJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
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
    * vincolo di velocita' di rotazione
    */
        
    
   /* Indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstVelIndex = iNode1FirstPosIndex+6;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstVelIndex = iNode2FirstPosIndex+6;
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+5;
   
   /* Nota: le reazioni vincolari sono: 
    * Forza,       3 incognite, riferimento globale, 
    * Momento,     2 incognite, riferimento locale
    */

   /* Setta gli indici dei nodi */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutRowIndex(6+iCnt, iNode1FirstVelIndex+iCnt);
      WM.PutColIndex(6+iCnt, iNode1FirstVelIndex+iCnt);
      WM.PutRowIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
      WM.PutColIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
      WM.PutRowIndex(18+iCnt, iNode2FirstVelIndex+iCnt);
      WM.PutColIndex(18+iCnt, iNode2FirstVelIndex+iCnt);
   }
   
   /* Setta gli indici delle reazioni */
   for (int iCnt = 1; iCnt <= 5; iCnt++) {
      WM.PutRowIndex(24+iCnt, iFirstReactionIndex+iCnt);
      WM.PutColIndex(24+iCnt, iFirstReactionIndex+iCnt);	
   }   
   
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.PutRowIndex(29+iCnt, iReactionPrimeIndex+iCnt);
      WM.PutColIndex(29+iCnt, iReactionPrimeIndex+iCnt);	
   }   
   
   /* Recupera i dati */
   const Mat3x3& R1(pNode1->GetRRef());
   const Mat3x3& R2(pNode2->GetRRef());
   const Vec3& Omega1(pNode1->GetWRef());
   const Vec3& Omega2(pNode2->GetWRef());   
   /* F ed M sono gia' state aggiornate da InitialAssRes */
   Vec3 FPrime(XCurr, iReactionPrimeIndex+1);
   Vec3 MPrime(XCurr, iReactionPrimeIndex+4);

   /* Matrici identita' */   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      /* Contributo di forza all'equazione della forza, nodo 1 */
      WM.PutCoef(iCnt, 24+iCnt, 1.);
      
      /* Contrib. di der. di forza all'eq. della der. della forza, nodo 1 */
      WM.PutCoef(6+iCnt, 29+iCnt, 1.);
      
      /* Contributo di forza all'equazione della forza, nodo 2 */
      WM.PutCoef(12+iCnt, 24+iCnt, -1.);
      
      /* Contrib. di der. di forza all'eq. della der. della forza, nodo 2 */
      WM.PutCoef(18+iCnt, 29+iCnt, -1.);
      
      /* Equazione di vincolo, nodo 1 */
      WM.PutCoef(24+iCnt, iCnt, -1.);
      
      /* Derivata dell'equazione di vincolo, nodo 1 */
      WM.PutCoef(29+iCnt, 6+iCnt, -1.);
      
      /* Equazione di vincolo, nodo 2 */
      WM.PutCoef(24+iCnt, 12+iCnt, 1.);
      
      /* Derivata dell'equazione di vincolo, nodo 2 */
      WM.PutCoef(29+iCnt, 18+iCnt, 1.);
   }
      
   /* Distanze e matrici di rotazione dai nodi alla cerniera 
    * nel sistema globale */
   Vec3 d1Tmp(R1*d1);
   Vec3 d2Tmp(R2*d2);
   Mat3x3 R1hTmp(R1*R1h);
   Mat3x3 R2hTmp(R2*R2h);
   
   Vec3 e3a(R1hTmp.GetVec(3));
   Vec3 e1b(R2hTmp.GetVec(1));
   Vec3 e2b(R2hTmp.GetVec(2));
   
   /* Ruota il momento e la sua derivata con le matrici della cerniera 
    * rispetto ai nodi */
   Vec3 MTmp(e2b*M.dGet(1)-e1b*M.dGet(2));
   Vec3 MPrimeTmp(e2b*MPrime.dGet(1)-e1b*MPrime.dGet(2));

   /* Matrici F/\d1/\, -F/\d2/\ */
   Mat3x3 FWedged1Wedge(MatCrossCross, F, d1Tmp);
   Mat3x3 FWedged2Wedge(MatCrossCross, F, -d2Tmp);
   
   /* Matrici (omega1/\d1)/\, -(omega2/\d2)/\ */
   Mat3x3 O1Wedged1Wedge(MatCross, Omega1.Cross(d1Tmp));
   Mat3x3 O2Wedged2Wedge(MatCross, d2Tmp.Cross(Omega2));
   
   Mat3x3 MDeltag1((Mat3x3(MatCross, Omega2.Cross(MTmp) + MPrimeTmp)
			   + Mat3x3(MatCrossCross, MTmp, Omega1))*Mat3x3(MatCross, e3a));
   Mat3x3 MDeltag2(Mat3x3(MatCrossCross, Omega1.Cross(e3a), MTmp)
		   + Mat3x3(MatCrossCross, e3a, MPrimeTmp)
		   + e3a.Cross(Mat3x3(MatCrossCross, Omega2, MTmp)));

   /* Vettori temporanei */
   Vec3 Tmp1(e2b.Cross(e3a));   
   Vec3 Tmp2(e3a.Cross(e1b));
   
   /* Prodotto vettore tra il versore 3 della cerniera secondo il nodo 1
    * ed il versore 1 della cerniera secondo il nodo 2. A convergenza
    * devono essere ortogonali, quindi il loro prodotto vettore deve essere 
    * unitario */

   /* Error handling: il programma si ferma, pero' segnala dov'e' l'errore */
   if (Tmp1.Dot() <= std::numeric_limits<doublereal>::epsilon() || Tmp2.Dot() <= std::numeric_limits<doublereal>::epsilon()) {
      silent_cerr("AxialRotationJoint(" << GetLabel() << "): "
	      "first and second node hinge axes are (nearly) orthogonal" 
	      << std::endl);
      throw ErrNullNorm(MBDYN_EXCEPT_ARGS);
   }   
   
   Vec3 TmpPrime1(e2b.Cross(Omega1.Cross(e3a))-e3a.Cross(Omega2.Cross(e2b)));
   Vec3 TmpPrime2(e3a.Cross(Omega2.Cross(e1b))-e1b.Cross(Omega1.Cross(e3a)));
   
   /* Equazione di momento, nodo 1 */
   WM.Add(4, 4, FWedged1Wedge - Mat3x3(MatCrossCross, MTmp, e3a));
   WM.Add(4, 16, Mat3x3(MatCrossCross, e3a, MTmp));
   WM.Add(4, 25, Mat3x3(MatCross, d1Tmp));
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutCoef(3+iCnt, 28, Tmp1(iCnt));
      WM.PutCoef(3+iCnt, 29, Tmp2(iCnt));
   }
   
   /* Equazione di momento, nodo 2 */
   WM.Add(4, 16, Mat3x3(MatCrossCross, MTmp, e3a));
   WM.Add(16, 16, FWedged2Wedge - Mat3x3(MatCrossCross, e3a, MTmp));
   WM.Sub(16, 25, Mat3x3(MatCross, d2Tmp));
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutCoef(15+iCnt, 28, -Tmp1(iCnt));
      WM.PutCoef(15+iCnt, 29, -Tmp2(iCnt));	
   }
 
   /* Derivata dell'equazione di momento, nodo 1 */
   WM.Add(10, 4, (Mat3x3(MatCross, FPrime) + Mat3x3(MatCrossCross, F, Omega1))*Mat3x3(MatCross, d1Tmp)
		   - MDeltag1
		   - Mat3x3(MatCross, e3a*MPrime(3)));
   WM.Add(10, 10, FWedged1Wedge
		   - Mat3x3(MatCrossCross, MTmp, e3a)
		   - Mat3x3(MatCross, e3a*MPrime(3)));
   WM.Add(10, 16, MDeltag2);
   WM.Add(10, 22, Mat3x3(MatCrossCross, e3a, MTmp));
   WM.Add(10, 25, O1Wedged1Wedge);
   WM.Add(10, 30, Mat3x3(MatCross, d1Tmp));
 
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutCoef(9+iCnt, 28, TmpPrime1(iCnt));
      WM.PutCoef(9+iCnt, 29, TmpPrime2(iCnt));
      WM.PutCoef(9+iCnt, 33, Tmp1(iCnt));
      WM.PutCoef(9+iCnt, 34, Tmp2(iCnt));	
      
      /* Contributo del momento attorno all'asse 3, dovuto alla velocita' 
       * imposta */
      WM.PutCoef(9+iCnt, 35, e3a(iCnt));	
   }
   
   /* Derivata dell'equazione di momento, nodo 2 */
   WM.Add(22, 4, MDeltag1 + Mat3x3(MatCross, e3a*MPrime(3)));
   WM.Add(22, 10, Mat3x3(MatCrossCross, MTmp, e3a) + Mat3x3(MatCross, e3a*MPrime(3)));
   WM.Sub(22, 16, (Mat3x3(MatCross, FPrime) + Mat3x3(MatCrossCross, F, Omega2))*Mat3x3(MatCross, d2Tmp) + MDeltag2);
   WM.Add(22, 22, FWedged2Wedge - Mat3x3(MatCrossCross, e3a, MTmp));
   WM.Add(22, 25, O2Wedged2Wedge);
   WM.Sub(22, 30, Mat3x3(MatCross, d2Tmp));

   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutCoef(21+iCnt, 28, -TmpPrime1(iCnt));
      WM.PutCoef(21+iCnt, 29, -TmpPrime2(iCnt));
      WM.PutCoef(21+iCnt, 33, -Tmp1(iCnt));
      WM.PutCoef(21+iCnt, 34, -Tmp2(iCnt));
      
      /* Contributo del momento attorno all'asse 3, dovuto alla velocita' 
       * imposta */
      WM.PutCoef(21+iCnt, 35, -e3a(iCnt));
   }
   
   /* Equazione di vincolo di posizione */
   WM.Add(25, 4, Mat3x3(MatCross, d1Tmp));
   WM.Sub(25, 16, Mat3x3(MatCross, d2Tmp));
   
   /* Derivata dell'equazione di vincolo di posizione */
   WM.Add(30, 4, O1Wedged1Wedge);
   WM.Add(30, 10, Mat3x3(MatCross, d1Tmp));
   WM.Add(30, 16, O2Wedged2Wedge);
   WM.Sub(30, 22, Mat3x3(MatCross, d2Tmp));
   
   /* Equazioni di vincolo di rotazione: e1b~e3a, e2b~e3a */
            
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = Tmp1(iCnt);
      WM.PutCoef(28, 3+iCnt, d);
      WM.PutCoef(28, 15+iCnt, -d);
      
      /* Queste sono per la derivata dell'equazione, sono qui solo per 
       * ottimizzazione */
      WM.PutCoef(33, 9+iCnt, d);
      WM.PutCoef(33, 21+iCnt, -d);
      
      d = Tmp2.dGet(iCnt);
      WM.PutCoef(29, 3+iCnt, -d);
      WM.PutCoef(29, 15+iCnt, d);
      
      /* Queste sono per la derivata dell'equazione, sono qui solo per 
       * ottimizzazione */
      WM.PutCoef(34, 9+iCnt, -d);
      WM.PutCoef(34, 21+iCnt, d);
   }   
   
   /* Derivate delle equazioni di vincolo di rotazione: e1b~e3a, e2b~e3a */
   Vec3 O1mO2(Omega1 - Omega2);
   TmpPrime1 = e3a.Cross(O1mO2.Cross(e2b));   
   TmpPrime2 = e2b.Cross(e3a.Cross(O1mO2));
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WM.PutCoef(33, 3+iCnt, TmpPrime1(iCnt));
      WM.PutCoef(33, 15+iCnt, TmpPrime2(iCnt));
   }
   
   TmpPrime1 = e3a.Cross(O1mO2.Cross(e1b));
   TmpPrime2 = e1b.Cross(e3a.Cross(O1mO2));
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WM.PutCoef(34, 3+iCnt, TmpPrime1(iCnt));
      WM.PutCoef(34, 15+iCnt, TmpPrime2(iCnt));
   }
   
   /* Equazione di vincolo di velocita' di rotazione; viene posta qui perche'
    * a questo numero di equazione corrisponde il numero della 
    * relativa incognita */
   Vec3 Tmp = O1mO2.Cross(e3a);
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutCoef(35, 3+iCnt, Tmp(iCnt));
      doublereal d = e3a(iCnt);
      WM.PutCoef(35, 9+iCnt, -d);
      WM.PutCoef(35, 21+iCnt, d);
   }
   
   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
AxialRotationJoint::InitialAssRes(SubVectorHandler& WorkVec,
				  const VectorHandler& XCurr)
{   
   DEBUGCOUT("Entering AxialRotationJoint::InitialAssRes()" << std::endl);
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);

   /* Indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstVelIndex = iNode1FirstPosIndex+6;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstVelIndex = iNode2FirstPosIndex+6;
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+5;
   
   /* Setta gli indici */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {	
      WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WorkVec.PutRowIndex(6+iCnt, iNode1FirstVelIndex+iCnt);
      WorkVec.PutRowIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
      WorkVec.PutRowIndex(18+iCnt, iNode2FirstVelIndex+iCnt);
   }
   
   for (int iCnt = 1; iCnt <= 5; iCnt++) {
      WorkVec.PutRowIndex(24+iCnt, iFirstReactionIndex+iCnt);
   }

   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.PutRowIndex(29+iCnt, iReactionPrimeIndex+iCnt);
   }

   /* Recupera i dati */
   const Vec3& x1(pNode1->GetXCurr());
   const Vec3& x2(pNode2->GetXCurr());
   const Vec3& v1(pNode1->GetVCurr());
   const Vec3& v2(pNode2->GetVCurr());
   const Mat3x3& R1(pNode1->GetRCurr());
   const Mat3x3& R2(pNode2->GetRCurr());
   const Vec3& Omega1(pNode1->GetWCurr());
   const Vec3& Omega2(pNode2->GetWCurr());
   
   /* Aggiorna F ed M, che restano anche per InitialAssJac */
   F = Vec3(XCurr, iFirstReactionIndex+1);
   M = Vec3(XCurr(iFirstReactionIndex+4),
	    XCurr(iFirstReactionIndex+5),
	    0.);
   Vec3 FPrime(XCurr, iReactionPrimeIndex+1);
   Vec3 MPrime(XCurr, iReactionPrimeIndex+4);
   
   /* Distanza nel sistema globale */
   Vec3 d1Tmp(R1*d1);
   Vec3 d2Tmp(R2*d2);
   Mat3x3 R1hTmp(R1*R1h);
   Mat3x3 R2hTmp(R2*R2h);

   /* Vettori omega1/\d1, -omega2/\d2 */
   Vec3 O1Wedged1(Omega1.Cross(d1Tmp));
   Vec3 O2Wedged2(Omega2.Cross(d2Tmp));
   
   Vec3 e3a(R1hTmp.GetVec(3));
   Vec3 e1b(R2hTmp.GetVec(1));
   Vec3 e2b(R2hTmp.GetVec(2));  
   
   /* Ruota il momento e la sua derivata con le matrici della cerniera 
    * rispetto ai nodi */
   Vec3 MTmp(e2b*M.dGet(1)-e1b*M.dGet(2));       
   Vec3 MPrimeTmp(e3a.Cross(MTmp.Cross(Omega2))+MTmp.Cross(Omega1.Cross(e3a))+
		  e2b.Cross(e3a)*MPrime.dGet(1)+e3a.Cross(e1b)*MPrime.dGet(2)); 
   
   /* Equazioni di equilibrio, nodo 1 */
   WorkVec.Sub(1, F);
   WorkVec.Add(4, F.Cross(d1Tmp)-MTmp.Cross(e3a)); /* Sfrutto il fatto che F/\d = -d/\F */
   
   /* Derivate delle equazioni di equilibrio, nodo 1 */
   WorkVec.Sub(7, FPrime);
   WorkVec.Add(10, FPrime.Cross(d1Tmp)-O1Wedged1.Cross(F)-MPrimeTmp-
	       e3a*MPrime.dGet(3));
   
   /* Equazioni di equilibrio, nodo 2 */
   WorkVec.Add(13, F);
   WorkVec.Add(16, d2Tmp.Cross(F)+MTmp.Cross(e3a));
   
   /* Derivate delle equazioni di equilibrio, nodo 2 */
   WorkVec.Add(19, FPrime);
   WorkVec.Add(22, d2Tmp.Cross(FPrime)+O2Wedged2.Cross(F)+MPrimeTmp+
	       e3a*MPrime.dGet(3));
   
   /* Equazione di vincolo di posizione */
   WorkVec.Add(25, x1+d1Tmp-x2-d2Tmp);
   
   /* Derivata dell'equazione di vincolo di posizione */
   WorkVec.Add(30, v1+O1Wedged1-v2-O2Wedged2);

   /* Equazioni di vincolo di rotazione */
   WorkVec.PutCoef(28, e2b.Dot(e3a));
   WorkVec.PutCoef(29, e1b.Dot(e3a));
   
   /* Derivate delle equazioni di vincolo di rotazione: e1b~e3a, e2b~e3a */
   Vec3 Tmp((Omega1-Omega2).Cross(e3a));
   WorkVec.PutCoef(33, e2b.Dot(Tmp));
   WorkVec.PutCoef(34, e1b.Dot(Tmp));

   /* Equazione di vincolo di velocita' di rotazione */
   doublereal Omega0 = pGetDriveCaller()->dGet();
   WorkVec.PutCoef(35, Omega0-e3a.Dot(Omega2-Omega1));

   return WorkVec;
}

unsigned int
AxialRotationJoint::iGetNumPrivData(void) const
{
	return 8;
}

unsigned int
AxialRotationJoint::iGetPrivDataIdx(const char *s) const
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
		idx += 3;
	case 'F':
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
AxialRotationJoint::dGetPrivData(unsigned int i) const
{
	ASSERT(i >= 1 && i <= iGetNumPrivData());
   
	switch (i) {
	case 1: {
		Vec3 v(RotManip::VecRot((pNode1->GetRCurr()*R1h).MulTM(pNode2->GetRCurr()*R2h)));
		doublereal dThetaTmp(v(3));

		int n = 0;

		if (dThetaTmp - dThetaWrapped < -M_PI) {
			n++;
		}

		if (dThetaTmp - dThetaWrapped > M_PI) {
			n--;
		}

		return 2*M_PI*(NTheta + n) + dThetaTmp;
	}
      
	case 2: 
		return dGet();
      
	case 3:
	case 4:
	case 5:
		return F(i - 2);
      
	case 6:
	case 7:
	case 8:
		return M(i - 5);
	}

	silent_cerr("AxialRotationJoint(" << GetLabel() << "): "
		"illegal private data " << i << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

/* AxialRotationJoint - end */


/* PlanePinJoint - begin */

/* Costruttore non banale */
PlanePinJoint::PlanePinJoint(unsigned int uL, const DofOwner* pDO,	       
			     const StructNode* pN,
			     const Vec3& X0Tmp, const Mat3x3& R0Tmp, 
			     const Vec3& dTmp, const Mat3x3& RhTmp,
			     flag fOut, const bool _calcInitdTheta,
			     const doublereal initDTheta)
: Elem(uL, fOut), 
Joint(uL, pDO, fOut), 
pNode(pN), 
X0(X0Tmp), R0(R0Tmp), d(dTmp), Rh(RhTmp),
F(Zero3), M(Zero3),
calcInitdTheta(_calcInitdTheta),
NTheta(0), dTheta(initDTheta), dThetaWrapped(0.)
{
   NO_OP;
}


/* Distruttore banale */
PlanePinJoint::~PlanePinJoint(void)
{
   NO_OP;
};


std::ostream&
PlanePinJoint::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
			"reaction forces [Fx,Fy,Fz]" << std::endl
		<< prefix << iIndex + 4 << "->" << iIndex + 5 << ": "
			"reaction couples [mx,my]" << std::endl;

	if (bInitial) {
		iIndex += 5;
		out
			<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
				"reaction force derivatives [FPx,FPy,FPz]" << std::endl
			<< prefix << iIndex + 4 << "->" << iIndex + 5 << ": "
				"reaction couple derivatives [mPx,mPy]" << std::endl;
	}

	return out;
}

void
PlanePinJoint::DescribeDof(std::vector<std::string>& desc, bool bInitial, int i) const
{
	std::ostringstream os;
	os << "PlanePinJoint(" << GetLabel() << ")";

	unsigned short nself = 5;
	if (bInitial) {
		nself *= 2;
	}

	if (i == -1) {
		desc.resize(nself);
		std::string name = os.str();

		for (unsigned i = 0; i < 3; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": reaction force f" << xyz[i];
			desc[i] = os.str();
		}

		for (unsigned i = 0; i < 2; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": reaction couple m" << xyz[i];
			desc[3 + i] = os.str();
		}

		if (bInitial) {
			for (unsigned i = 0; i < 3; i++) {
				os.str(name);
				os.seekp(0, std::ios_base::end);
				os << ": reaction force derivative fP" << xyz[i];
				desc[5 + i] = os.str();
			}
	
			for (unsigned i = 0; i < 2; i++) {
				os.str(name);
				os.seekp(0, std::ios_base::end);
				os << ": reaction couple derivative mP" << xyz[i];
				desc[8 + i] = os.str();
			}
		}

	} else {
		if (i < -1) {
			// error
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (i >= nself) {
			// error
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		desc.resize(1);

		switch (i) {
		case 0:
		case 1:
		case 2:
			os << ": reaction force f" << xyz[i];
			break;

		case 3:
		case 4:
			os << ": reaction couple m" << xyz[i - 3];
			break;

		case 5:
		case 6:
		case 7:
			os << ": reaction force derivative fP" << xyz[i - 5];
			break;

		case 8:
		case 9:
			os << ": reaction couple derivative mP" << xyz[i - 8];
			break;
		}
		desc[0] = os.str();
	}
}

std::ostream&
PlanePinJoint::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
			"position constraints [Px=Px0,Py=Py0,Pz=Pz0]" << std::endl
		<< prefix << iIndex + 4 << "->" << iIndex + 5 << ": "
			"orientation constraints [gx=gx0,gy=gy0]" << std::endl;

	if (bInitial) {
		iIndex += 5;
		out
			<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
				"velocity constraints [vx=0,vy=0,vz=0]" << std::endl
			<< prefix << iIndex + 4 << "->" << iIndex + 5 << ": "
				"angular velocity constraints [wx=0,wy=0]" << std::endl;
	}

	return out;
}

void
PlanePinJoint::DescribeEq(std::vector<std::string>& desc, bool bInitial, int i) const
{
	std::ostringstream os;
	os << "PlanePinJoint(" << GetLabel() << ")";

	unsigned short nself = 5;
	if (bInitial) {
		nself *= 2;
	}

	if (i == -1) {
		desc.resize(nself);
		std::string name = os.str();

		for (unsigned i = 0; i < 3; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": position constraint P" << xyz[i];
			desc[i] = os.str();
		}

		for (unsigned i = 0; i < 2; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": orientation constraint g" << xyz[i];
			desc[3 + i] = os.str();
		}

		if (bInitial) {
			for (unsigned i = 0; i < 3; i++) {
				os.str(name);
				os.seekp(0, std::ios_base::end);
				os << ": position constraint derivative v" << xyz[i];
				desc[5 + i] = os.str();
			}
	
			for (unsigned i = 0; i < 2; i++) {
				os.str(name);
				os.seekp(0, std::ios_base::end);
				os << ": orientation constraint derivative w" << xyz[i];
				desc[8 + i] = os.str();
			}
		}

	} else {
		if (i < -1) {
			// error
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (i >= nself) {
			// error
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		desc.resize(1);

		switch (i) {
		case 0:
		case 1:
		case 2:
			os << ": position constraint P" << xyz[i];
			break;

		case 3:
		case 4:
			os << ": orientation constraint g" << xyz[i - 3];
			break;

		case 5:
		case 6:
		case 7:
			os << ": position constraint derivative v" << xyz[i - 5];
			break;

		case 8:
		case 9:
			os << ": orientation constraint derivative w" << xyz[i - 8];
			break;
		}
		desc[0] = os.str();
	}
}

void
PlanePinJoint::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
	if (ph) {
		for (unsigned i = 0; i < ph->size(); i++) {
			Joint::JointHint *pjh = dynamic_cast<Joint::JointHint *>((*ph)[i]);

			if (pjh == 0) {
				continue;
			}

			if (dynamic_cast<Joint::OffsetHint<1> *>(pjh)) {
				const Mat3x3& R(pNode->GetRCurr());
   
				d = R.MulTV(X0 - pNode->GetXCurr());

			} else if (dynamic_cast<Joint::OffsetHint<0> *>(pjh)) {
				Vec3 dTmp(pNode->GetRCurr()*d);
   
				X0 = R0.MulTV(pNode->GetXCurr() + dTmp);

			} else if (dynamic_cast<Joint::HingeHint<1> *>(pjh)) {
				Rh = pNode->GetRCurr().MulTM(R0);

			} else if (dynamic_cast<Joint::HingeHint<2> *>(pjh)) {
				R0 = pNode->GetRCurr()*Rh;

			} else if (dynamic_cast<Joint::ReactionsHint *>(pjh)) {
				/* TODO */
			}
		}
	}

	if (calcInitdTheta) {
		Vec3 v(RotManip::VecRot(R0.MulTM(pNode->GetRCurr()*Rh)));
		dThetaWrapped = dTheta = v.dGet(3);
	}
}

Hint *
PlanePinJoint::ParseHint(DataManager *pDM, const char *s) const
{
	if (strncasecmp(s, "offset{" /*}*/, STRLENOF("offset{" /*}*/)) == 0) {
		s += STRLENOF("offset{" /*}*/);

		if (strcmp(&s[1], /* { */ "}") != 0) {
			return 0;
		}

		switch (s[0]) {
		case '1':
			return new Joint::OffsetHint<1>;

		case '0':
			return new Joint::OffsetHint<0>;
		}

	} else if (strncasecmp(s, "hinge{" /*}*/, STRLENOF("hinge{" /*}*/)) == 0) {
		s += STRLENOF("hinge{" /*}*/);

		if (strcmp(&s[1], /*{*/ "}") != 0) {
			return 0;
		}

		switch (s[0]) {
		case '1':
			return new Joint::HingeHint<1>;

		case '0':
			return new Joint::HingeHint<0>;
		}
	}

	return 0;
}

DofOrder::Order
PlanePinJoint::GetEqType(unsigned int i) const
{
	ASSERTMSGBREAK(i < iGetNumDof(), 
		"INDEX ERROR in PlanePinJoint::GetEqType");
   return DofOrder::ALGEBRAIC; 
}

void
PlanePinJoint::AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP)
{
	Vec3 v(RotManip::VecRot(R0.MulTM(pNode->GetRCurr()*Rh)));
	doublereal dThetaTmp(v(3));

	// unwrap
	if (dThetaTmp - dThetaWrapped < -M_PI) {
		NTheta++;
	}

	if (dThetaTmp - dThetaWrapped > M_PI) {
		NTheta--;
	}

	// save new wrapped angle
	dThetaWrapped = dThetaTmp;

	// compute new unwrapped angle
	dTheta = 2*M_PI*NTheta + dThetaWrapped;

}

void
PlanePinJoint::ReadInitialState(MBDynParser& HP)
{
	F = HP.GetVec3();
	M = Vec3(HP.GetReal(), HP.GetReal(), 0.);
}

/* Contributo al file di restart */
std::ostream& PlanePinJoint::Restart(std::ostream& out) const
{
   Joint::Restart(out) << ", revolute pin, "
     << pNode->GetLabel() 
     << ", reference, node, ", d.Write(out, ", ") 
     << ", hinge, reference, node, 1, ", 
     (Rh.GetVec(1)).Write(out, ", ") << ", 2, ", 
     (Rh.GetVec(2)).Write(out, ", ") 
     << ", reference, global, ", X0.Write(out, ", ") 
     << ",hinge, reference, global, 1, ",
     (R0.GetVec(1)).Write(out, ", ") << ", 2, ", 
     (R0.GetVec(2)).Write(out, ", ")  << ", " 
     << "initial theta, " << dTheta << ", "
     << "initial state, ", F.Write(out, ", ") 
     << ", " << M.dGet(1) << ", " << M.dGet(2) << ';' << std::endl;
   
   return out;
}


/* Assemblaggio jacobiano */
VariableSubMatrixHandler& 
PlanePinJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		      doublereal dCoef,
		      const VectorHandler& /* XCurr */ ,
		      const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering PlanePinJoint::AssJac()" << std::endl);
      
   SparseSubMatrixHandler& WM = WorkMat.SetSparse();
   WM.ResizeReset(39, 0);
   
   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();

   const Mat3x3& R(pNode->GetRRef());
   Vec3 dTmp(R*d);
   Mat3x3 RhTmp(R*Rh);
   
   
   /* 
    * L'equazione di vincolo afferma che il punto in cui si trova la
    * cerniera deve essere fissato:
    *      x + d = x0
    *      e3_0^Te1 = 0
    *      e3_0^Te2 = 0
    * 
    * con: d = R * d_0
    * 
    * La forza e' data dalla reazione vincolare F, nel sistema globale
    * La coppia dovuta all'eccentricita' e' data rispettivamente da:
    *     d /\ F
    *
    * 
    *       x      g         F
    * Q1 |  0      0             I   0    | | x |   | -F          |
    * G1 |  0      cF/\d1/\-M/\  d/\ e1e2 | | g |   | -d/\F-M     |
    * F  |  I      d/\           0   0    | | F |   |  (x+d-x0)/c |
    * M  |  0      e_0/\e1,e2    0   0    | | M |   |  e_0^Te1,e2 |
    * 
    * con d = R*d_0, c = dCoef
    */


   
   /* Moltiplica la forza ed il momento per il coefficiente
    * del metodo */

   Vec3 e3(R0.GetVec(3));
   Vec3 e1(RhTmp.GetVec(1));
   Vec3 e2(RhTmp.GetVec(2));
   Vec3 MTmp(e2*M.dGet(1)-e1*M.dGet(2));
            
   Vec3 Tmp1((e2).Cross(e3));
   Vec3 Tmp2((e3).Cross(e1));
   
   /* termini di reazione sul nodo (forza e momento) */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutItem(iCnt, iFirstMomentumIndex+iCnt, 
		  iFirstReactionIndex+iCnt, 1.);
      WM.PutItem(3+iCnt, 3+iFirstMomentumIndex+iCnt,
		  iFirstReactionIndex+4, Tmp1.dGet(iCnt));
      WM.PutItem(6+iCnt, 3+iFirstMomentumIndex+iCnt, 
		  iFirstReactionIndex+5, Tmp2.dGet(iCnt));
   }   
   
   WM.PutCross(10, iFirstMomentumIndex+3,
		iFirstReactionIndex, dTmp);
      
   
   /* Nota: F ed M, le reazioni vincolari, sono state aggiornate da AssRes */
   
   /* Termini diagonali del tipo: c*F/\d/\Delta_g 
    * nota: la forza e' gia' moltiplicata per dCoef */      
   WM.PutMat3x3(16, iFirstMomentumIndex+3, iFirstPositionIndex+3, 
		 Mat3x3(MatCrossCross, F*dCoef, dTmp) + Mat3x3(MatCrossCross, e3, MTmp*dCoef));

   /* Modifica: divido le equazioni di vincolo per dCoef */
   
   /* termini di vincolo dovuti al nodo 1 */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutItem(24+iCnt, iFirstReactionIndex+iCnt, 
		  iFirstPositionIndex+iCnt, -1.);
   }
   
   WM.PutCross(28, iFirstReactionIndex,
		iFirstPositionIndex+3, dTmp);
   
   for (int iCnt = 1; iCnt <= 3; iCnt ++) {
      WM.PutItem(33+iCnt, iFirstReactionIndex+4, 
		  iFirstPositionIndex+3+iCnt, -Tmp1.dGet(iCnt));	
      WM.PutItem(36+iCnt, iFirstReactionIndex+5, 
		  iFirstPositionIndex+3+iCnt, Tmp2.dGet(iCnt));	
   }
   
   return WorkMat;
}


/* Assemblaggio residuo */
SubVectorHandler& PlanePinJoint::AssRes(SubVectorHandler& WorkVec,
					      doublereal dCoef,
					      const VectorHandler& XCurr,
					      const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering PlanePinJoint::AssRes()" << std::endl);
      
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);
 
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   
   /* Indici dei nodi */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex+iCnt);
   }
     
   
   /* Indici del vincolo */
   for (int iCnt = 1; iCnt <= 5; iCnt++) {
      WorkVec.PutRowIndex(6+iCnt, iFirstReactionIndex+iCnt);   
   }
   
   F = Vec3(XCurr, iFirstReactionIndex+1);
   M = Vec3(XCurr(iFirstReactionIndex+4),
	    XCurr(iFirstReactionIndex+5),
	    0.);
   
   const Vec3& x(pNode->GetXCurr());
   const Mat3x3& R(pNode->GetRCurr());
   
   Vec3 dTmp(R*d);
   Mat3x3 RhTmp(R*Rh);
   
   Vec3 e3(R0.GetVec(3));
   Vec3 e1(RhTmp.GetVec(1));
   Vec3 e2(RhTmp.GetVec(2));
   
   WorkVec.Sub(1, F);
   WorkVec.Add(4, F.Cross(dTmp)-(e2*M.dGet(1)-e1*M.dGet(2)).Cross(e3)); /* Sfrutto il fatto che F/\d = -d/\F */
   
   /* Modifica: divido le equazioni di vincolo per dCoef */
   ASSERT(dCoef != 0.);
   WorkVec.Add(7, (x+dTmp-X0)/dCoef);

   WorkVec.PutCoef(10, e3.Dot(e2)/dCoef);
   WorkVec.PutCoef(11, e3.Dot(e1)/dCoef);

   return WorkVec;
}

/* Output (da mettere a punto) */
void PlanePinJoint::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) {
      Mat3x3 RTmp(pNode->GetRCurr()*Rh);
      Mat3x3 R0Tmp(R0.MulTM(RTmp));
      
      Joint::Output(OH.Joints(), "PlanePin", GetLabel(),
		    RTmp.MulTV(F), M, F, RTmp*M) 
	<< " " << MatR2EulerAngles(R0Tmp)*dRaDegr
	<< " " << RTmp.MulTV(pNode->GetWCurr()) << std::endl;      
   }
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
PlanePinJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
				   const VectorHandler& XCurr)
{
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeReset(iNumRows, iNumCols);

   /* Equazioni: vedi joints.dvi */
    
   /* Indici */
   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
   integer iFirstVelocityIndex = iFirstPositionIndex+6;
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+5;

   /* Setto gli indici */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.PutRowIndex(iCnt, iFirstPositionIndex+iCnt);
      WM.PutColIndex(iCnt, iFirstPositionIndex+iCnt);
      WM.PutRowIndex(6+iCnt, iFirstVelocityIndex+iCnt);
      WM.PutColIndex(6+iCnt, iFirstVelocityIndex+iCnt);
   }
   
   for (int iCnt = 1; iCnt <= 10; iCnt++) {
      WM.PutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
      WM.PutColIndex(12+iCnt, iFirstReactionIndex+iCnt);
   }   
   
   /* Matrici identita' */
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      /* Contributo di forza all'equazione della forza */
      WM.PutCoef(iCnt, 12+iCnt, 1.);
      
      /* Contrib. di der. di forza all'eq. della der. della forza */
      WM.PutCoef(6+iCnt, 17+iCnt, 1.);
      
      /* Equazione di vincolo */
      WM.PutCoef(12+iCnt, iCnt, -1.);
      
      /* Derivata dell'equazione di vincolo */
      WM.PutCoef(17+iCnt, 6+iCnt, -1.);
   }
   
   /* Recupera i dati */
   const Mat3x3& R(pNode->GetRRef());
   const Vec3& Omega(pNode->GetWRef());
   /* F, M sono state aggiornate da InitialAssRes */
   Vec3 FPrime(XCurr, iReactionPrimeIndex+1);
   Vec3 MPrime(XCurr(iReactionPrimeIndex+4),
	       XCurr(iReactionPrimeIndex+5),
	       0.);
   
   /* Distanza nel sistema globale */
   Vec3 dTmp(R*d);
   Mat3x3 RhTmp(R*Rh);

   Vec3 e3(R0.GetVec(3));
   Vec3 e1(RhTmp.GetVec(1));
   Vec3 e2(RhTmp.GetVec(2));

   /* Vettori temporanei */
   Vec3 Tmp1(e2.Cross(e3));   
   Vec3 Tmp2(e3.Cross(e1));
   
   /* Prodotto vettore tra il versore 3 della cerniera secondo il nodo 1
    * ed il versore 1 della cerniera secondo il nodo 2. A convergenza
    * devono essere ortogonali, quindi il loro prodotto vettore deve essere 
    * unitario */

   /* Error handling: il programma si ferma, pero' segnala dov'e' l'errore */
   if (Tmp1.Dot() <= std::numeric_limits<doublereal>::epsilon() || Tmp2.Dot() <= std::numeric_limits<doublereal>::epsilon()) {
      silent_cerr("PlanePinJoint(" << GetLabel() << "): "
	      "node and fixed point hinge axes are (nearly) orthogonal" 
	      << std::endl);
      throw ErrNullNorm(MBDYN_EXCEPT_ARGS);
   }   
   
   Vec3 TmpPrime1(e3.Cross(e2.Cross(Omega)));
   Vec3 TmpPrime2(e3.Cross(Omega.Cross(e1)));   
   
   /* Ruota il momento e la sua derivata con le matrici della cerniera 
    * rispetto ai nodi */
   Vec3 MTmp(e2*M.dGet(1)-e1*M.dGet(2));
   Vec3 MPrimeTmp(e2*MPrime.dGet(1)-e1*MPrime.dGet(2));
         
   Mat3x3 MDeltag(Mat3x3(MatCrossCross, e3, MPrimeTmp) + e3.Cross(Mat3x3(MatCrossCross, Omega, MTmp)));
   
   /* Matrici F/\d/\ */
   Mat3x3 FWedgedWedge(MatCrossCross, F, dTmp);
   
   /* Matrici (omega/\d)/\ */
   Mat3x3 OWedgedWedge(MatCross, Omega.Cross(dTmp));
   
   /* Equazione di momento */
   WM.Add(4, 4, FWedgedWedge + Mat3x3(MatCrossCross, e3, MTmp));
   WM.Add(4, 13, Mat3x3(MatCross, dTmp));
   
   /* Derivata dell'equazione di momento */
   WM.Add(10, 4, (Mat3x3(MatCross, FPrime) + Mat3x3(MatCrossCross, F, Omega))*Mat3x3(MatCross, dTmp) + MDeltag);
   WM.Add(10, 10, FWedgedWedge + Mat3x3(MatCrossCross, e3, MTmp));
   WM.Add(10, 13, OWedgedWedge);
   WM.Add(10, 18, Mat3x3(MatCross, dTmp));
 
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = Tmp1(iCnt);
      WM.PutCoef(3+iCnt, 16, d);
      WM.PutCoef(9+iCnt, 21, d);
      
      d = Tmp2(iCnt);
      WM.PutCoef(3+iCnt, 17, d);
      WM.PutCoef(9+iCnt, 22, d);
 
      WM.PutCoef(9+iCnt, 16, TmpPrime1(iCnt));
      WM.PutCoef(9+iCnt, 17, TmpPrime2(iCnt));
   }
 
   /* Equazione di vincolo */
   WM.Add(13, 4, Mat3x3(MatCross, dTmp));
 
   /* Derivata dell'equazione di vincolo */
   WM.Add(18, 4, OWedgedWedge);
   WM.Add(18, 10, Mat3x3(MatCross, dTmp));

   /* Equazioni di vincolo di rotazione: e1b~e3a, e2b~e3a */            
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      doublereal d = -Tmp1(iCnt);
      WM.PutCoef(16, 3+iCnt, d);
      
      /* Queste sono per la derivata dell'equazione, sono qui solo per 
       * ottimizzazione */
      WM.PutCoef(21, 9+iCnt, d);
      
      d = Tmp2(iCnt);
      WM.PutCoef(17, 3+iCnt, d);
      
      /* Queste sono per la derivata dell'equazione, sono qui solo per 
       * ottimizzazione */
      WM.PutCoef(22, 9+iCnt, d);
   }   
   
   /* Derivate delle equazioni di vincolo di rotazione: e1b~e3a, e2b~e3a */
   TmpPrime2 = e2.Cross(Omega.Cross(e3));
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WM.PutCoef(21, 3+iCnt, TmpPrime2(iCnt));
   }
   
   TmpPrime2 = e1.Cross(Omega.Cross(e3));
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutCoef(22, 3+iCnt, TmpPrime2(iCnt));
   }
 
   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
PlanePinJoint::InitialAssRes(SubVectorHandler& WorkVec,
			const VectorHandler& XCurr)
{   
   DEBUGCOUT("Entering PlanePinJoint::InitialAssRes()" << std::endl);
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);

   /* Indici */
   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
   integer iFirstVelocityIndex = iFirstPositionIndex+6;
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+5;
   
   /* Setta gli indici */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {	
      WorkVec.PutRowIndex(iCnt, iFirstPositionIndex+iCnt);
      WorkVec.PutRowIndex(6+iCnt, iFirstVelocityIndex+iCnt);
   }
   
   for (int iCnt = 1; iCnt <= 10; iCnt++) {	
      WorkVec.PutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
   }
   
   /* Recupera i dati */
   const Vec3& x(pNode->GetXCurr());
   const Vec3& v(pNode->GetVCurr());
   const Mat3x3& R(pNode->GetRCurr());
   const Vec3& Omega(pNode->GetWCurr());
   
   Mat3x3 RhTmp(R*Rh);
   
   F = Vec3(XCurr, iFirstReactionIndex+1);
   M = Vec3(XCurr(iFirstReactionIndex+4),
	    XCurr(iFirstReactionIndex+5),
	    0.);
   Vec3 FPrime(XCurr, iReactionPrimeIndex+1);
   Vec3 MPrime(XCurr(iReactionPrimeIndex+4),
	       XCurr(iReactionPrimeIndex+5),
	       0.);
   
   /* Versori delle cerniere */
   Vec3 e3(R0.GetVec(3));
   Vec3 e1(RhTmp.GetVec(1));
   Vec3 e2(RhTmp.GetVec(2));

   /* Vettori temporanei */
   Vec3 Tmp1(e2.Cross(e3));   
   Vec3 Tmp2(e3.Cross(e1));
      
   Vec3 TmpPrime1(e3.Cross(e2.Cross(Omega)));
   Vec3 TmpPrime2(e3.Cross(Omega.Cross(e1)));   
   
   /* Distanza nel sistema globale */
   Vec3 dTmp(R*d);

   /* Vettori omega/\d */
   Vec3 OWedged(Omega.Cross(dTmp));
   
   /* Ruota il momento e la sua derivata con le matrici della cerniera 
    * rispetto ai nodi */
   Vec3 MTmp(e2*M.dGet(1)-e1*M.dGet(2));       
   Vec3 MPrimeTmp(e3.Cross(MTmp.Cross(Omega))+
		  e2.Cross(e3)*MPrime.dGet(1)+e3.Cross(e1)*MPrime.dGet(2)); 
   
   /* Equazioni di equilibrio */
   WorkVec.Sub(1, F);
   WorkVec.Add(4, F.Cross(dTmp)-MTmp.Cross(e3)); /* Sfrutto il fatto che F/\d = -d/\F */
   
   /* Derivate delle equazioni di equilibrio, nodo 1 */
   WorkVec.Sub(7, FPrime);
   WorkVec.Add(10, FPrime.Cross(dTmp)-OWedged.Cross(F)-MPrimeTmp);
   
   /* Equazione di vincolo di posizione */
   WorkVec.Add(13, x+dTmp-X0);
   
   /* Equazioni di vincolo di rotazione */
   WorkVec.PutCoef(16, e2.Dot(e3));
   WorkVec.PutCoef(17, e1.Dot(e3));

   /* Derivata dell'equazione di vincolo di posizione */
   WorkVec.Add(18, v+OWedged);
   
   /* Derivate delle equazioni di vincolo di rotazione: e1b~e3a, e2b~e3a */
   Vec3 Tmp(e3.Cross(Omega));
   WorkVec.PutCoef(21, e2.Dot(Tmp));
   WorkVec.PutCoef(22, e1.Dot(Tmp));
      
   return WorkVec;
}

unsigned int
PlanePinJoint::iGetNumPrivData(void) const
{
	return 2;
}

unsigned int
PlanePinJoint::iGetPrivDataIdx(const char *s) const
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
		idx += 3;
	case 'F':
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
PlanePinJoint::dGetPrivData(unsigned int i) const
{
   ASSERT(i >= 1 && i <= iGetNumPrivData());
   
   switch (i) {
    case 1: {
	Vec3 v(RotManip::VecRot(R0.MulTM(pNode->GetRCurr()*Rh)));
	doublereal dThetaTmp(v(3));

	int n = 0;

	if (dThetaTmp - dThetaWrapped < -M_PI) {
		n++;
	}

	if (dThetaTmp - dThetaWrapped > M_PI) {
		n--;
	}

	return 2*M_PI*(NTheta + n) + dThetaTmp;
    }
      
    case 2: {
       Mat3x3 RTmp(pNode->GetRCurr()*Rh);
       Vec3 v(RTmp.MulTV(pNode->GetWCurr()));
       
       return v(3);
    }

    case 3:
    case 4:
    case 5:
	    return F(i - 2);

    case 6:
    case 7:
    case 8:
	    return M(i - 5);
   }
      
   silent_cerr("PlanePinJoint(" << GetLabel() << "): "
	   "illegal private data " << i << std::endl);
   throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

/* PlanePinJoint - end */
