/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2005
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

/* elementi di massa */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <limits>

#include "body.h"
#include "dataman.h"

/* Body - begin */

Body::Body(unsigned int uL,
	const StructNode *pNode,
	doublereal dMass,
	const Vec3& Xgc,
	const Mat3x3& J,
	flag fOut)
: Elem(uL, fOut),
ElemGravityOwner(uL, fOut),
InitialAssemblyElem(uL, fOut),
pNode(pNode),
dMass(dMass),
Xgc(Xgc),
S0(Xgc*dMass),
J0(J)
{
	ASSERT(pNode != NULL);
	ASSERT(pNode->GetNodeType() == Node::STRUCTURAL);
	ASSERT(dMass > 0.);
}


/* distruttore */
Body::~Body(void)
{
	NO_OP;
}


/* momento statico */
Vec3
Body::GetS_int(void) const
{
	return pNode->GetXCurr()*dMass + pNode->GetRCurr()*S0;
}


/* momento d'inerzia */
Mat3x3
Body::GetJ_int(void) const
{
	Vec3 s = pNode->GetRCurr()*S0;
	const Vec3& x = pNode->GetXCurr();

	return pNode->GetRCurr()*J0.MulMT(pNode->GetRCurr())
		- Mat3x3(MatCrossCross, x, x*dMass)
		- Mat3x3(MatCrossCross, s, x)
		- Mat3x3(MatCrossCross, x, s);
}


/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
Body::Restart_int(std::ostream& out) const
{
	out << "  body: " << GetLabel() << ", "
		<< pNode->GetLabel() << ", " << dMass << ", "
		<< "reference, node, ", Xgc.Write(out, ", ") << ", "
		<< "reference, node, ", (J0 + Mat3x3(MatCrossCross, S0, Xgc)).Write(out, ", ");

	return out;
}


void
Body::AfterPredict(VectorHandler& /* X */ , VectorHandler& /* XP */ )
{
	const Mat3x3& R = pNode->GetRRef();

	STmp = R*S0;
	JTmp = R*J0.MulMT(R);
}

/* massa totale */
doublereal
Body::dGetM(void) const
{
	return dMass;
}

/* momento statico */
Vec3
Body::GetS(void) const
{
	return GetS_int();
}

/* momento d'inerzia */
Mat3x3
Body::GetJ(void) const
{
	return GetJ_int();
}

/* nodo */
const StructNode *
Body::pGetNode(void) const
{
	return pNode;
}

/* Accesso ai dati privati */
unsigned int
Body::iGetNumPrivData(void) const
{
	return 1;
}

unsigned int
Body::iGetPrivDataIdx(const char *s) const
{
	if (strcmp(s, "E") == 0) {
		// kinetic energy
		return 1;
	}

	if (strcmp(s, "V") == 0) {
		// potential energy
		return 2;
	}

	return 0;
}

doublereal
Body::dGetPrivData(unsigned int i) const
{
	if (i == 1) {
		// kinetic energy
		const Mat3x3& Rn = pNode->GetRCurr();
		const Vec3& Vn = pNode->GetVCurr();
		const Vec3& Wn = pNode->GetWCurr();

		Vec3 X = Rn*Xgc;
		Vec3 V = Vn + Wn.Cross(X);
		Vec3 W = Rn*Wn;

		Mat3x3 Jgc = J0 + Mat3x3(MatCrossCross, Xgc, Xgc*dMass);

		return ((V*V)*dMass + W*(Jgc*W))/2.;
	}

	if (i == 2) {
		// potential energy
		const Vec3& Xn = pNode->GetXCurr();

		Vec3 GravityAcceleration;
		if (GravityOwner::bGetGravity(pNode->GetXCurr(), GravityAcceleration)) {
			return -Xn.Dot(GravityAcceleration)*dMass;
		}
	}

	return 0.;
}

void
Body::AssVecRBK_int(SubVectorHandler& WorkVec)
{
	const RigidBodyKinematics *pRBK = pNode->pGetRBK();

	Vec3 s0;

	integer iIdx = 0;
	if (dynamic_cast<DynamicBody *>(this)) {
		iIdx = 6;
	}

	s0 = pNode->GetXCurr()*dMass + STmp;

	// force
	Vec3 f;
	f = pRBK->GetXPP()*dMass;
	f += pRBK->GetWP().Cross(s0);
	f += pRBK->GetW().Cross(pRBK->GetW().Cross(s0));

	WorkVec.Sub(iIdx + 1, f);

	// moment
	Vec3 a;
	a = pRBK->GetXPP();
	a += pRBK->GetWP().Cross(pNode->GetXCurr());
	a += pRBK->GetW().Cross(pRBK->GetW().Cross(pNode->GetXCurr()));
	a += pRBK->GetW().Cross(pNode->GetVCurr());

	Vec3 m;
	m = STmp.Cross(a);
	m += pRBK->GetW().Cross(JTmp*pRBK->GetW());
	m += JTmp*pRBK->GetWP();
	m += pNode->GetWCurr().Cross(JTmp*pRBK->GetW());
	m -= JTmp*(pNode->GetWCurr().Cross(pRBK->GetW()));
	m += pNode->GetVCurr().Cross(pRBK->GetW().Cross(STmp));

	WorkVec.Sub(iIdx + 3 + 1, m);
}

void
Body::AssMatsRBK_int(
	FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB,
	const doublereal& dCoef,
	const Vec3& Sc)
{
	const RigidBodyKinematics *pRBK = pNode->pGetRBK();

	Mat3x3 MTmp;
	Vec3 VTmp;

	integer iIdx = 0;
	if (dynamic_cast<DynamicBody *>(this)) {
		iIdx = 6;
	}

	// f: delta x
	MTmp = Mat3x3(MatCross, pRBK->GetWP());
	MTmp += Mat3x3(MatCrossCross, pRBK->GetW(), pRBK->GetW());

	WMA.Add(iIdx + 1, 1, MTmp*(dMass*dCoef));


	// f: theta delta

	WMA.Sub(iIdx + 1, 3 + 1, MTmp*Mat3x3(MatCross, Sc));


	// m: delta x
	MTmp = Mat3x3(MatCrossCross, Sc, pRBK->GetWP());
	MTmp += Sc.Cross(Mat3x3(MatCrossCross, pRBK->GetW(), pRBK->GetW()));

	WMA.Add(iIdx + 3 + 1, 1, MTmp);


	// m: theta delta

	VTmp = pRBK->GetXPP();
	VTmp += pRBK->GetWP().Cross(pNode->GetXCurr());
	VTmp += pRBK->GetW().Cross(pRBK->GetW().Cross(pNode->GetXCurr()));
	VTmp += pRBK->GetW().Cross(pNode->GetVCurr());

	MTmp = Mat3x3(MatCrossCross, VTmp, Sc);

	VTmp = (pRBK->GetW() + pNode->GetWCurr())*dCoef;

	Mat3x3 MTmp2(JTmp*Mat3x3(MatCross, pRBK->GetW()) - Mat3x3(MatCross, JTmp*pRBK->GetW()));
	MTmp += VTmp.Cross(MTmp2);

	VTmp = (pRBK->GetWP() + pRBK->GetW().Cross(pNode->GetWCurr()))*dCoef;

	MTmp += JTmp*Mat3x3(MatCross, VTmp);
	MTmp -= Mat3x3(MatCross, JTmp*VTmp);

	MTmp -= pNode->GetVCurr().Cross(Mat3x3(MatCrossCross, pRBK->GetW(), Sc));

	WMA.Add(iIdx + 3 + 1, 3 + 1, MTmp);


	// m: delta dot x
	MTmp = Mat3x3(MatCrossCross, STmp, pRBK->GetW());
	MTmp -= Mat3x3(MatCross, pRBK->GetW().Cross(STmp));

	WMB.Add(iIdx + 3 + 1, 1, MTmp);

	// m: delta omega
	WMB.Add(iIdx + 3 + 1, 3 + 1, MTmp2);
}

/* Body - end */


/* DynamicBody - begin */

DynamicBody::DynamicBody(unsigned int uL,
	const DynamicStructNode* pNode,
	doublereal dMass,
	const Vec3& Xgc,
	const Mat3x3& J,
	flag fOut)
: Elem(uL, fOut),
Body(uL, pNode, dMass, Xgc, J, fOut)
{
	NO_OP;
}


/* distruttore */
DynamicBody::~DynamicBody(void)
{
	NO_OP;
}


/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
DynamicBody::Restart(std::ostream& out) const
{
	Body::Restart_int(out) << ';' << std::endl;

	return out;
}


VariableSubMatrixHandler&
DynamicBody::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("DynamicBody::AssJac");

	/* Casting di WorkMat */
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	Vec3 GravityAcceleration;
	bool g = GravityOwner::bGetGravity(pNode->GetXCurr(),
		GravityAcceleration);

	integer iNumRows = 6;
	if (g || pNode->pGetRBK()) {
		iNumRows = 12;
	}

	/* Dimensiona e resetta la matrice di lavoro */
	WM.ResizeReset(iNumRows, 6);

	/* Setta gli indici della matrice - le incognite sono ordinate come:
	 *   - posizione (3)
	 *   - parametri di rotazione (3)
	 *   - quantita' di moto (3)
	 *   - momento della quantita' di moto
	 * e gli indici sono consecutivi. La funzione pGetFirstPositionIndex()
	 * ritorna il valore del primo indice -1, in modo che l'indice i-esimo
	 * e' dato da iGetFirstPositionIndex()+i */
	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
		WM.PutColIndex(iCnt, iFirstPositionIndex + iCnt);
	}

	if (iNumRows == 12) {
		integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
		for (integer iCnt = 1; iCnt <= 6; iCnt++) {
			WM.PutRowIndex(6 + iCnt, iFirstMomentumIndex + iCnt);
		}
	}

	AssMats(WM, WM, dCoef, g, GravityAcceleration);

	return WorkMat;
}


void
DynamicBody::AssMats(VariableSubMatrixHandler& WorkMatA,
	VariableSubMatrixHandler& WorkMatB,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("DynamicBody::AssMats");

	/* Casting di WorkMat */
	FullSubMatrixHandler& WMA = WorkMatA.SetFull();
	FullSubMatrixHandler& WMB = WorkMatB.SetFull();

	Vec3 GravityAcceleration;
	bool g = GravityOwner::bGetGravity(pNode->GetXCurr(),
		GravityAcceleration);

	integer iNumRows = 6;
	if (g) {
		iNumRows = 12;
	}

	/* Dimensiona e resetta la matrice di lavoro */
	WMA.ResizeReset(iNumRows, 6);
	WMB.ResizeReset(6, 6);

	/* Setta gli indici della matrice - le incognite sono ordinate come:
	 *   - posizione (3)
	 *   - parametri di rotazione (3)
	 *   - quantita' di moto (3)
	 *   - momento della quantita' di moto
	 * e gli indici sono consecutivi. La funzione pGetFirstPositionIndex()
	 * ritorna il valore del primo indice -1, in modo che l'indice i-esimo
	 * e' dato da iGetFirstPositionIndex()+i */
	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WMA.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
		WMA.PutColIndex(iCnt, iFirstPositionIndex + iCnt);

		WMB.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
		WMB.PutColIndex(iCnt, iFirstPositionIndex + iCnt);
	}

	if (g) {
		integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
		for (integer iCnt = 1; iCnt <= 6; iCnt++) {
			WMA.PutRowIndex(6 + iCnt, iFirstMomentumIndex + iCnt);
		}
	}

	AssMats(WMA, WMB, 1., g, GravityAcceleration);
}


void
DynamicBody::AssMats(FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB,
	doublereal dCoef,
	bool bGravity,
	const Vec3& GravityAcceleration)
{
	DEBUGCOUTFNAME("DynamicBody::AssMats");

	const Vec3& V(pNode->GetVCurr());
	const Vec3& W(pNode->GetWCurr());

	// STmp, JTmp computed by AssRes()
	// const Mat3x3& R(pNode->GetRCurr());
	// STmp = R*S0;
	// JTmp = R*J0.MulMT(R);

	Mat3x3 SWedge(MatCross, STmp);			/* S /\ */
	Vec3 Sc(STmp*dCoef);

	/*
	 * momentum:
	 *
	 * m * I DeltaV - S /\ DeltagP + ( S /\ W ) /\ Deltag
	 */
	WMB.IncCoef(1, 1, dMass);
	WMB.IncCoef(2, 2, dMass);
	WMB.IncCoef(3, 3, dMass);

	WMB.Sub(1, 3 + 1, SWedge);
	WMA.Add(1, 3 + 1, Mat3x3(MatCross, Sc.Cross(W)));

	/*
	 * momenta moment:
	 *
	 * S /\ DeltaV + J DeltagP + ( V /\ S /\ - ( J * W ) /\ ) Deltag
	 */
	WMB.Add(3 + 1, 1, SWedge);

	WMB.Add(3 + 1, 3 + 1, JTmp);
	WMA.Add(3 + 1, 3 + 1, Mat3x3(MatCrossCross, V, Sc) - Mat3x3(MatCross, JTmp*(W*dCoef)));

	if (bGravity) {
		WMA.Sub(9 + 1, 3 + 1, Mat3x3(MatCrossCross, GravityAcceleration, Sc));
	}

	const RigidBodyKinematics *pRBK = pNode->pGetRBK();
	if (pRBK) {
		AssMatsRBK_int(WMA, WMB, dCoef, Sc);
	}
}


SubVectorHandler&
DynamicBody::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUTFNAME("Body::DynamicAssRes");

	/* Se e' definita l'accelerazione di gravita',
	 * la aggiunge (solo al residuo) */
	Vec3 GravityAcceleration;
	bool g = GravityOwner::bGetGravity(pNode->GetXCurr(),
		GravityAcceleration);

	const RigidBodyKinematics *pRBK = pNode->pGetRBK();

	integer iNumRows = 6;
	if (g || pRBK) {
		iNumRows = 12;
	}

	WorkVec.ResizeReset(iNumRows);

	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	for (integer iCnt = 1; iCnt <= iNumRows; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
	}

	const Vec3& V(pNode->GetVCurr());
	const Vec3& W(pNode->GetWCurr());

	/* Aggiorna i suoi dati (saranno pronti anche per AssJac) */
	const Mat3x3& R(pNode->GetRCurr());
	STmp = R*S0;
	JTmp = R*J0.MulMT(R);

	/* Quantita' di moto: R[1] = Q - M * V - W /\ S */
	WorkVec.Sub(1, V*dMass + W.Cross(STmp));

	/* Momento della quantita' di moto: R[2] = G - S /\ V - J * W */
	WorkVec.Sub(3 + 1, JTmp*W + STmp.Cross(V));

	if (g) {
		WorkVec.Add(6 + 1, GravityAcceleration*dMass);
		/* FIXME: this should go into Jacobian matrix
		 * as Gravity /\ S /\ Delta g */
		WorkVec.Add(9 + 1, STmp.Cross(GravityAcceleration));
	}

	if (pRBK) {
		AssVecRBK_int(WorkVec);
	}

	const DynamicStructNode *pDN = dynamic_cast<const DynamicStructNode *>(pNode);
	ASSERT(pDN != 0);

	pDN->AddInertia(dMass, STmp, JTmp);

	return WorkVec;
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
DynamicBody::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUTFNAME("DynamicBody::InitialAssJac");

	/* Casting di WorkMat */
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Setta gli indici della matrice - le incognite sono ordinate come:
	 *   - posizione (3)
	 *   - parametri di rotazione (3)
	 *   - quantita' di moto (3)
	 *   - momento della quantita' di moto
	 * e gli indici sono consecutivi. La funzione pGetFirstPositionIndex()
	 * ritorna il valore del primo indice -1, in modo che l'indice i-esimo
	 * e' dato da iGetFirstPositionIndex()+i */
	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	integer iFirstVelocityIndex = iFirstPositionIndex+6;
	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iFirstPositionIndex+iCnt);
		WM.PutRowIndex(6+iCnt, iFirstVelocityIndex+iCnt);
		WM.PutColIndex(iCnt, iFirstPositionIndex+iCnt);
	}

	/* Prepara matrici e vettori */

	/* Velocita' angolare corrente */
	const Vec3& W(pNode->GetWRef());

	Vec3 SWedgeW(STmp.Cross(W));
	Mat3x3 WWedgeSWedge(MatCrossCross, -W, STmp);
	Mat3x3 WWedge(MatCross, W);
	Mat3x3 WWedgeWWedgeSWedge(W.Cross(WWedgeSWedge));
	Mat3x3 FDeltaW(Mat3x3(MatCross, SWedgeW) + WWedgeSWedge);

	// STmp, JTmp computed by InitialAssRes()
	Vec3 JW(JTmp*W);
	Mat3x3 JWWedge(MatCross, JW);
	Mat3x3 MDeltag(W.Cross(JTmp*WWedge - JWWedge));
	Mat3x3 MDeltaW(W.Cross(JTmp) - JWWedge);

	/* Forza */
	WM.Add(1, 1, WWedgeWWedgeSWedge);
	WM.Add(1, 4, FDeltaW);

	/* Momento */
	WM.Add(4, 1, MDeltag);
	WM.Add(4, 4, MDeltaW);

	/* Derivata forza */
	WM.Add(7, 1, Mat3x3(MatCross, W.Cross(SWedgeW)) + W.Cross(FDeltaW));
	WM.Add(7, 4, W.Cross(WWedgeWWedgeSWedge));

	/* Derivata Momento */
	WM.Add(4, 1, W.Cross(MDeltag));
	WM.Add(4, 4, W.Cross(MDeltaW) - Mat3x3(MatCross, W.Cross(JW)));

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
DynamicBody::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUTFNAME("DynamicBody::InitialAssRes");

	integer iNumRows;
	integer iNumCols;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	for (integer iCnt = 1; iCnt <= 12; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstPositionIndex+iCnt);
	}

	const Vec3& X(pNode->GetXCurr());
	const Vec3& W(pNode->GetWCurr());

	// Aggiorna i suoi dati (saranno pronti anche per InitialAssJac)
	const Mat3x3& R(pNode->GetRCurr());
	STmp = R*S0;
	JTmp = R*J0.MulMT(R);

	Vec3 FC(-W.Cross(W.Cross(STmp)));
	Vec3 MC(-W.Cross(JTmp*W));

	/* Forza */
	WorkVec.Add(1, FC);

	/* Momento */
	WorkVec.Add(4, MC);

	/* Derivata forza */
	WorkVec.Add(7, W.Cross(FC));

	/* Derivata momento */
	WorkVec.Add(10, W.Cross(MC));

	/* Se e' definita l'accelerazione di gravita',
	 * la aggiunge (solo al residuo) */
	Vec3 GravityAcceleration;
	if (GravityOwner::bGetGravity(X, GravityAcceleration)) {
		WorkVec.Add(1, GravityAcceleration*dMass);
		WorkVec.Add(3 + 1, STmp.Cross(GravityAcceleration));
		WorkVec.Add(9 + 1, (W.Cross(STmp)).Cross(GravityAcceleration));
	}

	return WorkVec;
}


/* Usata per inizializzare la quantita' di moto */
void
DynamicBody::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& /* XP */ ,
	SimulationEntity::Hints *ph)
{
	integer iFirstIndex = pNode->iGetFirstMomentumIndex();

	const Vec3& V(pNode->GetVCurr());
	const Vec3& W(pNode->GetWCurr());
	const Mat3x3& R(pNode->GetRCurr());
	STmp = R*S0;
	JTmp = R*J0.MulMT(R);
	X.Add(iFirstIndex + 1, V*dMass + W.Cross(STmp));
	X.Add(iFirstIndex + 4, STmp.Cross(V) + JTmp*W);
}

/* momentum */
Vec3
DynamicBody::GetB_int(void) const
{
	const Vec3& V(pNode->GetVCurr());
	const Vec3& W(pNode->GetWCurr());
	const Mat3x3& R(pNode->GetRCurr());

	return V*dMass + W.Cross(R*S0);
}


/* momenta moment */
Vec3
DynamicBody::GetG_int(void) const
{
	const Vec3& X(pNode->GetXCurr());
	const Mat3x3& R(pNode->GetRCurr());
	const Vec3& V(pNode->GetVCurr());
	const Vec3& W(pNode->GetWCurr());

	Vec3 STmp(R*S0);

	// NOTE: with respect to the origin of the global reference frame!
	return (STmp + X*dMass).Cross(V) + R*(J0*(R.MulTV(W)))
		- X.Cross(STmp.Cross(W));
}


/* DynamicBody - end */


/* StaticBody - begin */

StaticBody::StaticBody(unsigned int uL,
	const StaticStructNode* pNode,
	doublereal dMass,
	const Vec3& Xgc,
	const Mat3x3& J,
	flag fOut)
: Elem(uL, fOut),
Body(uL, pNode, dMass, Xgc, J, fOut)
{
	NO_OP;
}


/* distruttore */
StaticBody::~StaticBody(void)
{
	NO_OP;
}


/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
StaticBody::Restart(std::ostream& out) const
{
	Body::Restart_int(out);

	out << ';' << std::endl;

	return out;
}


VariableSubMatrixHandler&
StaticBody::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("StaticBody::AssJac");

	/* Casting di WorkMat */
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	WM.ResizeReset(6, 6);

	/* Setta gli indici della matrice - le incognite sono ordinate come:
	 *   - posizione (3)
	 *   - parametri di rotazione (3)
	 *   - quantita' di moto (3)
	 *   - momento della quantita' di moto
	 * e gli indici sono consecutivi. La funzione pGetFirstPositionIndex()
	 * ritorna il valore del primo indice -1, in modo che l'indice i-esimo
	 * e' dato da iGetFirstPositionIndex() + i */
	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iFirstMomentumIndex + iCnt);
		WM.PutColIndex(iCnt, iFirstPositionIndex + iCnt);
	}

	if (AssMats(WM, WM, dCoef)) {
		WorkMat.SetNullMatrix();
	}

	return WorkMat;
}


void
StaticBody::AssMats(VariableSubMatrixHandler& WorkMatA,
	VariableSubMatrixHandler& WorkMatB,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("StaticBody::AssMats");

	/* Casting di WorkMat */
	FullSubMatrixHandler& WMA = WorkMatA.SetFull();
	FullSubMatrixHandler& WMB = WorkMatB.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	WMA.ResizeReset(6, 6);
	WMB.ResizeReset(6, 6);

	/* Setta gli indici della matrice - le incognite sono ordinate come:
	 *   - posizione (3)
	 *   - parametri di rotazione (3)
	 *   - quantita' di moto (3)
	 *   - momento della quantita' di moto
	 * e gli indici sono consecutivi. La funzione pGetFirstPositionIndex()
	 * ritorna il valore del primo indice -1, in modo che l'indice i-esimo
	 * e' dato da iGetFirstPositionIndex() + i */
	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WMA.PutRowIndex(iCnt, iFirstMomentumIndex + iCnt);
		WMA.PutColIndex(iCnt, iFirstPositionIndex + iCnt);

		WMB.PutRowIndex(iCnt, iFirstMomentumIndex + iCnt);
		WMB.PutColIndex(iCnt, iFirstPositionIndex + iCnt);
	}

	if (AssMats(WMA, WMB, 1.)) {
		WorkMatA.SetNullMatrix();
		WorkMatB.SetNullMatrix();
	}
}


bool
StaticBody::AssMats(FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB,
	doublereal dCoef)
{
	DEBUGCOUTFNAME("StaticBody::AssMats");

	/* Se e' definita l'accelerazione di gravita',
	 * la aggiunge (solo al residuo) */
	Vec3 Acceleration(Zero3);
	bool g = GravityOwner::bGetGravity(pNode->GetXCurr(), Acceleration);

	/* TODO: reference */
	Vec3 W(Zero3);

	const RigidBodyKinematics *pRBK = pNode->pGetRBK();

	if (!g && !pRBK) {
		/* Caller will set WMA & WMB to null matrix */
		return true;
	}

	Vec3 Sc(STmp*dCoef);

	if (g) {
		WMA.Add(3 + 1, 3 + 1, Mat3x3(MatCrossCross, Acceleration, Sc));
	}

	if (pRBK) {
		AssMatsRBK_int(WMA, WMB, dCoef, Sc);
	}

	return false;
}


SubVectorHandler&
StaticBody::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUTFNAME("StaticBody::AssRes");

	/* Se e' definita l'accelerazione di gravita',
	 * la aggiunge (solo al residuo) */
	Vec3 Acceleration(Zero3);
	bool g = GravityOwner::bGetGravity(pNode->GetXCurr(), Acceleration);

	/* W is uninitialized because its use is conditioned by w */
	const RigidBodyKinematics *pRBK = pNode->pGetRBK();

	if (!g && !pRBK) {
		WorkVec.Resize(0);
		return WorkVec;
	}

	WorkVec.ResizeReset(6);

	integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex + iCnt);
	}

	/* Aggiorna i suoi dati (saranno pronti anche per AssJac) */
	const Mat3x3& R(pNode->GetRCurr());
	STmp = R*S0;
	JTmp = R*J0.MulMT(R);

	if (g) {
		WorkVec.Add(1, Acceleration*dMass);
		WorkVec.Add(3 + 1, STmp.Cross(Acceleration));
	}

	if (pRBK) {
		AssVecRBK_int(WorkVec);
	}

	return WorkVec;
}


/* inverse dynamics capable element */
bool
StaticBody::bInverseDynamics(void) const
{
	return true;
}


SubVectorHandler&
StaticBody::AssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ ,
	const VectorHandler& /* XPrimePrimeCurr */ ,
	InverseDynamics::Order iOrder)
{
	DEBUGCOUTFNAME("DynamicBody::AssRes");
	
	ASSERT(iOrder == InverseDynamics::INVERSE_DYNAMICS);

	/* Se e' definita l'accelerazione di gravita', la aggiunge */
	Vec3 GravityAcceleration;
	bool g = GravityOwner::bGetGravity(pNode->GetXCurr(),
		GravityAcceleration);

	WorkVec.ResizeReset(6);

	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
	}

	const Mat3x3& R(pNode->GetRCurr());
	Vec3 XgcTmp = R*Xgc;
	STmp = R*S0;
	JTmp = R*J0.MulMT(R);

	Vec3 Acceleration = pNode->GetXPPCurr()
		+ pNode->GetWPCurr().Cross(XgcTmp)
		+ pNode->GetWCurr().Cross(pNode->GetWCurr().Cross(XgcTmp));
	if (g) {
		Acceleration -= GravityAcceleration;
	}

	WorkVec.Sub(1, Acceleration*dMass);

	Vec3 M = JTmp*pNode->GetWPCurr()
		+ STmp.Cross(pNode->GetXPPCurr())
		+ STmp.Cross(pNode->GetWCurr().Cross(pNode->GetWCurr().Cross(XgcTmp)));
	if (g) {
		M -= GravityAcceleration*dMass;
	}

	WorkVec.Sub(4, M);

	return WorkVec;
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
StaticBody::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUTFNAME("StaticBody::InitialAssJac");

	WorkMat.SetNullMatrix();

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
StaticBody::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUTFNAME("StaticBody::InitialAssRes");

	WorkVec.ResizeReset(0);

	return WorkVec;
}


/* Usata per inizializzare la quantita' di moto */
void
StaticBody::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& /* XP */ ,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

/* StaticBody - end */


/* Legge un corpo rigido */
Elem*
ReadBody(DataManager* pDM, MBDynParser& HP, unsigned int uLabel)
{
	DEBUGCOUTFNAME("ReadBody");

	const char* sKeyWords[] = {
		NULL
	};

	/* enum delle parole chiave */
	enum KeyWords {
		UNKNOWN = -1,
		LASTKEYWORD = 0
	};

	/* tabella delle parole chiave */
	KeyTable K(HP, sKeyWords);

	/* nodo collegato */
	Node* pStrNode = pDM->ReadNode(HP, Node::STRUCTURAL);

	/* may be determined by a special DataManager parameter... */
	bool bStaticModel = pDM->bIsStaticModel();
	bool bInverseDynamics = pDM->bIsInverseDynamics();

	integer iNumMasses = 1;
	if (HP.IsKeyWord("condense")) {
		iNumMasses = HP.GetInt();
		if (iNumMasses < 1) {
			silent_cerr("Body(" << uLabel << "): "
				"at least one mass is required in \"condense\" "
				"at line " << HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		DEBUGLCOUT(MYDEBUG_INPUT,
			iNumMasses << " masses will be condensed" << std::endl);

		/* The inertia is calculated as follows:
		 *
		 * dm = Sum(dm_i)
		 *
		 * Xgc = Sum(Xgc_i*dm_i)/Sum(dm_i)
		 *
		 * J = Sum(J_i)-Sum(dm_i*(Xgc_i-Xgc)/\(Xgc_i-Xgc)/\)
		 *
		 * and it can be accomplished by accumulating:
		 *
		 * dm = Sum(dm_i)
		 *
		 * ~S = Sum(Xgc_i*dm_i)
		 *
		 * ~J = Sum(J_i)-Sum(dm_i*Xgc_i/\*Xgc_i/\)
		 *
		 * then calculating
		 *
		 * Xgc = S/dm
		 *
		 * and finally:
		 *
		 * J = ~J-Xgc/\(dm*Xgc/\-2*~S)
		 *
		 */
	}

	doublereal dm = 0.;
	ReferenceFrame RF(dynamic_cast<StructNode *>(pStrNode));
	Vec3 Xgc(Zero3);
	Vec3 STmp(Zero3);
	Mat3x3 J(Zero3x3);
	bool bNegative(false);

	if (HP.IsKeyWord("allow" "negative" "mass")) {
		bNegative = true;
	}

	for (int iCnt = 1; iCnt <= iNumMasses; iCnt++) {
		/* massa */
		doublereal dMTmp = HP.GetReal();
		if (!bNegative && dMTmp < 0.) {
			silent_cerr("Body(" << uLabel << "): "
				"negative mass is not allowed at line "
				<< HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		DEBUGLCOUT(MYDEBUG_INPUT, "Mass(" << iCnt << ") = " << dMTmp << std::endl);
		dm += dMTmp;

		/* posiz. c.g. */
		Vec3 XgcTmp(HP.GetPosRel(RF));
		if (iNumMasses == 1) {
			Xgc = XgcTmp;

		} else {
			STmp += XgcTmp*dMTmp;
		}

		DEBUGLCOUT(MYDEBUG_INPUT, "position of mass(" << iCnt
			<< ") center of gravity = " << XgcTmp << std::endl);

		/*
		 * matrice del mom. d'inerzia
		 *
		 * Usa la funzione che legge una matrice qualsiasi con parole chiave
		 * per forme abbreviate:
		 *   - null: matrice vuota
		 *   - eye:  matrice identita'
		 *   - diag: matrice diagonale, seguita da 3 reali
		 *   - sym:  matrice simmetrica, seguita da 6 reali,
		 *           letta come triangolare superiore ordinata per righe:
		 *           m11, m12, m13,    m22, m23,    m33
		 *   - matrice generica, costituita da 9 reali, letta per righe:
		 *           m11, m12, m13,    m21, m22, m23,   m31, m32, m33
		 *
		 * Si assume inoltre che la matrice dei momenti di inerzia
		 * sia espressa nel centro di massa del corpo, quindi viene
		 * corretta per l'eventuale offset rispetto al nodo
		 */
		Mat3x3 JTmp(HP.GetMatRel(RF));
		DEBUGLCOUT(MYDEBUG_INPUT, "Inertia matrix of mass(" << iCnt
			<< ") =" << std::endl << JTmp << std::endl);

		if (HP.IsKeyWord("inertial")) {
			DEBUGLCOUT(MYDEBUG_INPUT,
				"supplied in inertial reference frame" << std::endl);
			if (HP.IsKeyWord("node")) {
				NO_OP;
			} else {
				Mat3x3 RTmp(HP.GetRotRel(RF));
				JTmp = RTmp*JTmp.MulMT(RTmp);
			}
			DEBUGLCOUT(MYDEBUG_INPUT,
				"Inertia matrix of mass(" << iCnt << ") "
				"in current frame =" << JTmp << std::endl);
		}

		J += JTmp - Mat3x3(MatCrossCross, XgcTmp, XgcTmp*dMTmp);
	}

	if (!bNegative && dm < 0.) {
		silent_cerr("Body(" << uLabel << "): "
			"negative mass is not allowed at line "
			<< HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (iNumMasses > 1) {
		if (dm < std::numeric_limits<doublereal>::epsilon()) {
			silent_cerr("Body(" << uLabel << "): "
				"mass value " << dm << " is too small at line "
				<< HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		Xgc = STmp/dm;
	}

	DEBUGLCOUT(MYDEBUG_INPUT, "Total mass: " << dm << std::endl
		<< "Center of mass: " << Xgc << std::endl
		<< "Inertia matrix:" << std::endl << J << std::endl);

	DynamicStructNode* pDynamicNode = 0;
	StaticStructNode* pStaticNode = 0;

	if (bStaticModel || bInverseDynamics) {
		/* static */
		pStaticNode = dynamic_cast<StaticStructNode*>(pStrNode);
		if (pStaticNode == 0 || pStaticNode->GetStructNodeType() != StructNode::STATIC) {
			silent_cerr("Body(" << uLabel << "): "
				"illegal structural node type "
				"for StructNode(" << pStrNode->GetLabel() << ") "
				"at line " << HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else {
		pDynamicNode = dynamic_cast<DynamicStructNode*>(pStrNode);
		if (pDynamicNode == 0 || pDynamicNode->GetStructNodeType() != StructNode::DYNAMIC) {
			silent_cerr("Body(" << uLabel << "): "
				"illegal structural node type "
				"for StructNode(" << pStrNode->GetLabel() << ") "
				"at line " << HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	flag fOut = pDM->fReadOutput(HP, Elem::BODY);

	/* Allocazione e costruzione */
	Elem* pEl = NULL;
	if (bStaticModel || bInverseDynamics) {
		/* static */
		SAFENEWWITHCONSTRUCTOR(pEl, StaticBody,
			StaticBody(uLabel, pStaticNode,
				dm, Xgc, J, fOut));

	} else {
		SAFENEWWITHCONSTRUCTOR(pEl, DynamicBody,
			DynamicBody(uLabel, pDynamicNode,
				dm, Xgc, J, fOut));
	}

	pDM->GetLogFile()
		<< "body: " << uLabel
		<< ' ' << pStrNode->GetLabel()
		<< ' ' << dm
		<< ' ' << Xgc
		<< ' ' << J
		<< std::endl;

	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		silent_cerr("semicolon expected "
			"at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pEl;
} /* End of ReadBody() */

