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
#include "body_vm.h"
#include "dataman.h"

/* Mass - begin */

Mass::Mass(unsigned int uL,
	const StructDispNode *pNode,
	doublereal dMass,
	flag fOut)
: Elem(uL, fOut),
ElemGravityOwner(uL, fOut),
InitialAssemblyElem(uL, fOut),
pNode(pNode),
dMass(dMass)
{
	ASSERT(pNode != NULL);
	ASSERT(pNode->GetNodeType() == Node::STRUCTURAL);
	ASSERT(dMass > 0.);
}


/* distruttore */
Mass::~Mass(void)
{
	NO_OP;
}


/* momento statico */
Vec3
Mass::GetS_int(void) const
{
	return pNode->GetXCurr()*dMass;
}


/* momento d'inerzia */
Mat3x3
Mass::GetJ_int(void) const
{
	const Vec3& x = pNode->GetXCurr();

	return Mat3x3(MatCrossCross, x, x*(-dMass));
}


/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
Mass::Restart(std::ostream& out) const
{
	out << "  body: " << GetLabel() << ", "
		<< pNode->GetLabel() << ", " << dMass << ';' << std::endl;

	return out;
}


/* massa totale */
doublereal
Mass::dGetM(void) const
{
	return dMass;
}

/* momento statico */
Vec3
Mass::GetS(void) const
{
	return GetS_int();
}

/* momento d'inerzia */
Mat3x3
Mass::GetJ(void) const
{
	return GetJ_int();
}

/* nodo */
const StructDispNode *
Mass::pGetNode(void) const
{
	return pNode;
}

/* Accesso ai dati privati */
unsigned int
Mass::iGetNumPrivData(void) const
{
	return 3;
}

unsigned int
Mass::iGetPrivDataIdx(const char *s) const
{
	if (s[1] == '\0') {
		switch (s[0]) {
		case 'E':
			// kinetic energy
			return 1;

		case 'V':
			// potential energy
			return 2;

		case 'm':
			return 3;
		}
	}

	return 0;
}

doublereal
Mass::dGetPrivData(unsigned int i) const
{
	switch (i) {
	case 1: {
		// kinetic energy
		const Vec3& Vn = pNode->GetVCurr();

		return Vn.Dot()*dMass;
		}

	case 2: {
		// potential energy
		const Vec3& Xn = pNode->GetXCurr();

		Vec3 GravityAcceleration;
		if (GravityOwner::bGetGravity(Xn, GravityAcceleration)) {
			return -Xn.Dot(GravityAcceleration)*dMass;
		}
		break;
		}

	case 3:
 		// mass
 		return dMass;
	}

	return 0.;
}

void
Mass::AssVecRBK_int(SubVectorHandler& WorkVec)
{
	const RigidBodyKinematics *pRBK = pNode->pGetRBK();

	Vec3 s0;

	integer iIdx = 0;
	if (dynamic_cast<DynamicMass *>(this)) {
		iIdx = 3;
	}

	s0 = pNode->GetXCurr()*dMass;

	// force
	Vec3 f;
	f = pRBK->GetXPP()*dMass;

	WorkVec.Sub(iIdx + 1, f);
}

void
Mass::AssMatsRBK_int(
	FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB,
	const doublereal& dCoef)
{
	const RigidBodyKinematics *pRBK = pNode->pGetRBK();

	integer iIdx = 0;
	if (dynamic_cast<DynamicMass *>(this)) {
		iIdx = 3;
	}

	// f: delta x
	Mat3x3 MTmp(MatCross, pRBK->GetWP());
	MTmp += Mat3x3(MatCrossCross, pRBK->GetW(), pRBK->GetW());

	WMA.Add(iIdx + 1, 1, MTmp*(dMass*dCoef));
}

/* Mass - end */


/* DynamicMass - begin */

DynamicMass::DynamicMass(unsigned int uL,
	const DynamicStructDispNode* pNode,
	doublereal dMass,
	flag fOut)
: Elem(uL, fOut),
Mass(uL, pNode, dMass, fOut)
{
	NO_OP;
}


/* distruttore */
DynamicMass::~DynamicMass(void)
{
	NO_OP;
}


VariableSubMatrixHandler&
DynamicMass::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("DynamicMass::AssJac");

	/* Casting di WorkMat */
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	Vec3 GravityAcceleration;
	bool g = GravityOwner::bGetGravity(pNode->GetXCurr(),
		GravityAcceleration);

	integer iNumRows = 3;
	if (pNode->pGetRBK()) {
		iNumRows = 6;
	}

	/* Dimensiona e resetta la matrice di lavoro */
	WM.ResizeReset(iNumRows, 3);

	/* Setta gli indici della matrice - le incognite sono ordinate come:
	 *   - posizione (3)
	 *   - parametri di rotazione (3)
	 *   - quantita' di moto (3)
	 *   - momento della quantita' di moto
	 * e gli indici sono consecutivi. La funzione pGetFirstPositionIndex()
	 * ritorna il valore del primo indice -1, in modo che l'indice i-esimo
	 * e' dato da iGetFirstPositionIndex()+i */
	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
		WM.PutColIndex(iCnt, iFirstPositionIndex + iCnt);
	}

	if (iNumRows == 6) {
		integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
		for (integer iCnt = 1; iCnt <= 3; iCnt++) {
			WM.PutRowIndex(3 + iCnt, iFirstMomentumIndex + iCnt);
		}
	}

	AssMats(WM, WM, dCoef, g, GravityAcceleration);

	return WorkMat;
}


void
DynamicMass::AssMats(VariableSubMatrixHandler& WorkMatA,
	VariableSubMatrixHandler& WorkMatB,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("DynamicMass::AssMats");

	/* Casting di WorkMat */
	FullSubMatrixHandler& WMA = WorkMatA.SetFull();
	FullSubMatrixHandler& WMB = WorkMatB.SetFull();

	Vec3 GravityAcceleration;
	bool g = GravityOwner::bGetGravity(pNode->GetXCurr(),
		GravityAcceleration);

	integer iNumRows = 3;
	if (g) {
		iNumRows = 6;
	}

	/* Dimensiona e resetta la matrice di lavoro */
	WMA.ResizeReset(iNumRows, 3);
	WMB.ResizeReset(3, 3);

	/* Setta gli indici della matrice - le incognite sono ordinate come:
	 *   - posizione (3)
	 *   - parametri di rotazione (3)
	 *   - quantita' di moto (3)
	 *   - momento della quantita' di moto
	 * e gli indici sono consecutivi. La funzione pGetFirstPositionIndex()
	 * ritorna il valore del primo indice -1, in modo che l'indice i-esimo
	 * e' dato da iGetFirstPositionIndex()+i */
	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WMA.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
		WMA.PutColIndex(iCnt, iFirstPositionIndex + iCnt);

		WMB.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
		WMB.PutColIndex(iCnt, iFirstPositionIndex + iCnt);
	}

	if (g) {
		integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
		for (integer iCnt = 1; iCnt <= 3; iCnt++) {
			WMA.PutRowIndex(3 + iCnt, iFirstMomentumIndex + iCnt);
		}
	}

	AssMats(WMA, WMB, 1., g, GravityAcceleration);
}


void
DynamicMass::AssMats(FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB,
	doublereal dCoef,
	bool bGravity,
	const Vec3& GravityAcceleration)
{
	DEBUGCOUTFNAME("DynamicMass::AssMats");

	/*
	 * momentum:
	 *
	 * m * I DeltaV - S /\ DeltagP + ( S /\ W ) /\ Deltag
	 */
	WMB.IncCoef(1, 1, dMass);
	WMB.IncCoef(2, 2, dMass);
	WMB.IncCoef(3, 3, dMass);

	const RigidBodyKinematics *pRBK = pNode->pGetRBK();
	if (pRBK) {
		AssMatsRBK_int(WMA, WMB, dCoef);
	}
}


SubVectorHandler&
DynamicMass::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUTFNAME("DynamicMass::AssRes");

	/* Se e' definita l'accelerazione di gravita',
	 * la aggiunge (solo al residuo) */
	Vec3 GravityAcceleration;
	bool g = GravityOwner::bGetGravity(pNode->GetXCurr(),
		GravityAcceleration);

	const RigidBodyKinematics *pRBK = pNode->pGetRBK();

	integer iNumRows = 3;
	if (g || pRBK) {
		iNumRows = 6;
	}

	WorkVec.ResizeReset(iNumRows);

	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	for (integer iCnt = 1; iCnt <= iNumRows; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
	}

	const Vec3& V(pNode->GetVCurr());

	/* Quantita' di moto: R[1] = Q */
	WorkVec.Sub(1, V*dMass);

	if (g) {
		WorkVec.Add(3 + 1, GravityAcceleration*dMass);
	}

	if (pRBK) {
		AssVecRBK_int(WorkVec);
	}

	const DynamicStructDispNode *pDN = dynamic_cast<const DynamicStructDispNode *>(pNode);
	ASSERT(pDN != 0);

	pDN->AddInertia(dMass);

	return WorkVec;
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
DynamicMass::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUTFNAME("DynamicMass::InitialAssJac");

	/* Casting di WorkMat */
	WorkMat.SetNullMatrix();

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
DynamicMass::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUTFNAME("DynamicMass::InitialAssRes");

	WorkVec.Resize(0);

	return WorkVec;
}


/* Usata per inizializzare la quantita' di moto */
void
DynamicMass::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& /* XP */ ,
	SimulationEntity::Hints *ph)
{
	integer iFirstIndex = pNode->iGetFirstMomentumIndex();

	const Vec3& V(pNode->GetVCurr());
	X.Add(iFirstIndex + 1, V*dMass);
}

/* momentum */
Vec3
DynamicMass::GetB_int(void) const
{
	const Vec3& V(pNode->GetVCurr());

	return V*dMass;
}

/* DynamicMass - end */


/* StaticMass - begin */

StaticMass::StaticMass(unsigned int uL,
	const StaticStructDispNode* pNode,
	doublereal dMass,
	flag fOut)
: Elem(uL, fOut),
Mass(uL, pNode, dMass, fOut)
{
	NO_OP;
}


/* distruttore */
StaticMass::~StaticMass(void)
{
	NO_OP;
}


VariableSubMatrixHandler&
StaticMass::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("StaticMass::AssJac");

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
StaticMass::AssMats(VariableSubMatrixHandler& WorkMatA,
	VariableSubMatrixHandler& WorkMatB,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("StaticMass::AssMats");

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
StaticMass::AssMats(FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB,
	doublereal dCoef)
{
	DEBUGCOUTFNAME("StaticMass::AssMats");

	/* Se e' definita l'accelerazione di gravita',
	 * la aggiunge (solo al residuo) */
	Vec3 Acceleration(Zero3);
	bool g = GravityOwner::bGetGravity(pNode->GetXCurr(), Acceleration);

	const RigidBodyKinematics *pRBK = pNode->pGetRBK();

	if (!g && !pRBK) {
		/* Caller will set WMA & WMB to null matrix */
		return true;
	}

	if (pRBK) {
		AssMatsRBK_int(WMA, WMB, dCoef);
	}

	return false;
}


SubVectorHandler&
StaticMass::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUTFNAME("StaticMass::AssRes");

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

	WorkVec.ResizeReset(3);

	integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex + iCnt);
	}

	if (g) {
		WorkVec.Add(1, Acceleration*dMass);
	}

	if (pRBK) {
		AssVecRBK_int(WorkVec);
	}

	return WorkVec;
}


/* inverse dynamics capable element */
bool
StaticMass::bInverseDynamics(void) const
{
	return true;
}


SubVectorHandler&
StaticMass::AssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ ,
	const VectorHandler& /* XPrimePrimeCurr */ ,
	InverseDynamics::Order iOrder)
{
	DEBUGCOUTFNAME("DynamicMass::AssRes");
	
	ASSERT(iOrder == InverseDynamics::INVERSE_DYNAMICS);

	/* Se e' definita l'accelerazione di gravita', la aggiunge */
	Vec3 GravityAcceleration;
	bool g = GravityOwner::bGetGravity(pNode->GetXCurr(),
		GravityAcceleration);

	WorkVec.ResizeReset(3);

	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
	}

	Vec3 Acceleration = pNode->GetXPPCurr();
	if (g) {
		Acceleration -= GravityAcceleration;
	}

	WorkVec.Sub(1, Acceleration*dMass);

	return WorkVec;
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
StaticMass::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUTFNAME("StaticMass::InitialAssJac");

	WorkMat.SetNullMatrix();

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
StaticMass::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUTFNAME("StaticMass::InitialAssRes");

	WorkVec.Resize(0);

	return WorkVec;
}


/* Usata per inizializzare la quantita' di moto */
void
StaticMass::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& /* XP */ ,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

/* StaticMass - end */


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
Body::Restart(std::ostream& out) const
{
	out << "  body: " << GetLabel() << ", "
		<< pNode->GetLabel() << ", " << dMass << ", "
		<< "reference, node, ", Xgc.Write(out, ", ") << ", "
		<< "reference, node, ", (J0 + Mat3x3(MatCrossCross, S0, Xgc)).Write(out, ", ")
		<< ";" << std::endl;

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
	return 3;
}

unsigned int
Body::iGetPrivDataIdx(const char *s) const
{
	if (s[1] == '\0') {
		switch (s[0]) {
		case 'E':
			// kinetic energy
			return 1;

		case 'V':
			// potential energy
			return 2;

		case 'm':
			return 3;
		}
	}

	return 0;
}

doublereal
Body::dGetPrivData(unsigned int i) const
{
	switch (i) {
	case 1: {
		// kinetic energy
		const Mat3x3& Rn = pNode->GetRCurr();
		const Vec3& Vn = pNode->GetVCurr();
		const Vec3& Wn = pNode->GetWCurr();

		Vec3 X = Rn*Xgc;
		Vec3 V = Vn + Wn.Cross(X);
		Vec3 W = Rn.MulTV(Wn);

		Mat3x3 Jgc = J0 + Mat3x3(MatCrossCross, Xgc, Xgc*dMass);

		return ((V*V)*dMass + W*(Jgc*W))/2.;
	}

	case 2: {
		// potential energy
		Vec3 X(pNode->GetXCurr() + pNode->GetRCurr()*Xgc);
		Vec3 GravityAcceleration;
		if (GravityOwner::bGetGravity(X, GravityAcceleration)) {
			return -X.Dot(GravityAcceleration)*dMass;
		}
		break;
	}

	case 3:
		return dMass;
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

#ifdef USE_SPARSE_AUTODIFF
template <typename T>
void
Body::AssVecRBK_int(const sp_grad::SpColVector<T, 3>& STmp,
                    const sp_grad::SpMatrix<T, 3, 3>& JTmp,
                    sp_grad::SpGradientAssVec<T>& WorkVec,
                    doublereal dCoef,
                    sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     const RigidBodyKinematics *pRBK = pNode->pGetRBK();

     SpColVector<T, 3> X(3, 1), V(3, 1), W(3, 0);

     pNode->GetXCurr(X, dCoef, func);
     pNode->GetVCurr(V, dCoef, func);
     pNode->GetWCurr(W, dCoef, func);

     SpColVector<T, 3> s0 = X * dMass + STmp;

     // force
     SpColVector<T, 3> F = pRBK->GetXPP() * -dMass
          - Cross(pRBK->GetWP(), s0)
          - Cross(pRBK->GetW(), Cross(pRBK->GetW(), s0));

     // moment
     SpColVector<T, 3> a = pRBK->GetXPP()
          + Cross(pRBK->GetWP(), X)
          + Cross(pRBK->GetW(), Cross(pRBK->GetW(), X))
          + Cross(pRBK->GetW(), V);

     const SpColVector<T, 3> JTmpWRBK = JTmp * pRBK->GetW();

     SpColVector<T, 3> M = -Cross(STmp, a)
          - Cross(pRBK->GetW(), JTmpWRBK)
          - JTmp * pRBK->GetWP()
          - Cross(W, JTmpWRBK)
          + JTmp * Cross(W, pRBK->GetW())
          - Cross(V, Cross(pRBK->GetW(), STmp));

     const integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();

     WorkVec.AddItem(iFirstMomentumIndex + 1, F);
     WorkVec.AddItem(iFirstMomentumIndex + 4, M);
}
#endif

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

void
DynamicBody::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{ 
	*piNumRows = 12; 
	*piNumCols = 6; 
}

void 
DynamicBody::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{ 
	*piNumRows = 12; 
	*piNumCols = 6; 
}

VariableSubMatrixHandler&
DynamicBody::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("DynamicBody::AssJac");

#ifndef USE_SPARSE_AUTODIFF
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

#else
        using namespace sp_grad;

        SpGradientAssVec<SpGradient>::AssJac(this,
                                             WorkMat.SetSparseGradient(),
                                             dCoef,
                                             XCurr,
                                             XPrimeCurr,
                                             SpFunctionCall::REGULAR_JAC);
#endif

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
                    doublereal dCoef,
                    const VectorHandler& XCurr ,
                    const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("DynamicBody::AssRes");

#ifndef USE_SPARSE_AUTODIFF
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
#else
        using namespace sp_grad;

        SpGradientAssVec<doublereal>::AssRes(this,
                                             WorkVec,
                                             dCoef,
                                             XCurr,
                                             XPrimeCurr,
                                             SpFunctionCall::REGULAR_RES);
#endif

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

	// TODO: make configurable
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

#ifdef USE_SPARSE_AUTODIFF
template <typename T>
void
DynamicBody::AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                    doublereal dCoef,
                    const sp_grad::SpGradientVectorHandler<T>& XCurr,
                    const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
                    sp_grad::SpFunctionCall func)
{
        using namespace sp_grad;

        Vec3 GravityAcceleration;
        bool g = GravityOwner::bGetGravity(pNode->GetXCurr(),
                                           GravityAcceleration);

        const RigidBodyKinematics *pRBK = pNode->pGetRBK();

        const integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();

        SpColVectorA<T, 3> V;
        SpColVectorA<T, 3> W;
        SpMatrixA<T, 3, 3> R;

        pNode->GetVCurr(V, dCoef, func);
        pNode->GetWCurr(W, dCoef, func);
        pNode->GetRCurr(R, dCoef, func);

        SpColVector<T, 3> STmp = R * S0;
        SpMatrix<T, 3, 3> JTmp = R * J0 * Transpose(R);
        SpColVector<T, 3> f1 = V * -dMass - Cross(W, STmp);
        SpColVector<T, 3> f2 = -Cross(STmp, V) - JTmp * W;

        WorkVec.AddItem(iFirstPositionIndex + 1, f1);
        WorkVec.AddItem(iFirstPositionIndex + 4, f2);

        if (g) {
                SpColVector<T, 3> f3 = GravityAcceleration * dMass;
                SpColVector<T, 3> f4 = Cross(STmp, GravityAcceleration);

                WorkVec.AddItem(iFirstPositionIndex + 7, f3);
                WorkVec.AddItem(iFirstPositionIndex + 10, f4);
        }

        if (pRBK) {
                AssVecRBK_int(STmp, JTmp, WorkVec, dCoef, func);
        }

        UpdateInertia(STmp, JTmp);
}

void DynamicBody::UpdateInertia(const sp_grad::SpColVector<doublereal, 3>& S,
                                const sp_grad::SpMatrix<doublereal, 3, 3>& J) const
{
        const DynamicStructNode *pDN = dynamic_cast<const DynamicStructNode*>(pNode);

        ASSERT(pDN != 0);

        for (integer i = 1; i <= 3; ++i) {
             STmp(i) = S(i);
        }

        for (integer j = 1; j <= 3; ++j) {
             for (integer i = 1; i <= 3; ++i) {
                  JTmp(i, j) = J(i, j);
             }
        }

        pDN->AddInertia(dMass, STmp, JTmp);
}

void DynamicBody::UpdateInertia(const sp_grad::SpColVector<sp_grad::SpGradient, 3>& STmp,
                                const sp_grad::SpMatrix<sp_grad::SpGradient, 3, 3>& JTmp) const
{
}
#endif

/* DynamicBody - end */


/* ModalBody - begin */

ModalBody::ModalBody(unsigned int uL,
	const ModalNode* pNode,
	doublereal dMass,
	const Vec3& Xgc,
	const Mat3x3& J,
	flag fOut)
: Elem(uL, fOut),
DynamicBody(uL, pNode, dMass, Xgc, J, fOut),
XPP(::Zero3), WP(::Zero3)
{
	NO_OP;
}

/* distruttore */
ModalBody::~ModalBody(void)
{
	NO_OP;
}

void
ModalBody::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{ 
	*piNumRows = 12; 
	*piNumCols = 12; 
}

VariableSubMatrixHandler&
ModalBody::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("ModalBody::AssJac");

	/* Casting di WorkMat */
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	Vec3 GravityAcceleration;
	bool g = GravityOwner::bGetGravity(pNode->GetXCurr(),
		GravityAcceleration);

	const integer iNumRows = 12;

	/* Dimensiona e resetta la matrice di lavoro */
	WM.ResizeReset(iNumRows, 12);

	/* Setta gli indici della matrice - le incognite sono ordinate come:
	 *   - posizione (3)
	 *   - parametri di rotazione (3)
	 *   - quantita' di moto (3)
	 *   - momento della quantita' di moto
	 * e gli indici sono consecutivi. La funzione pGetFirstPositionIndex()
	 * ritorna il valore del primo indice -1, in modo che l'indice i-esimo
	 * e' dato da iGetFirstPositionIndex()+i */
	const integer iFirstIndexModal = pNode->iGetFirstIndex();
	for (integer iCnt = 1; iCnt <= 12; iCnt++) {
		WM.PutRowIndex(iCnt, iFirstIndexModal + iCnt);
	}

	for (integer iCnt = 1; iCnt <= 12; ++iCnt) {
		WM.PutColIndex(iCnt, iFirstIndexModal + iCnt);
	}

	AssMats(WM, WM, dCoef, XCurr, XPrimeCurr, g, GravityAcceleration);

	return WorkMat;
}


void
ModalBody::AssMats(VariableSubMatrixHandler& WorkMatA,
	VariableSubMatrixHandler& WorkMatB,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("ModalBody::AssMats");

	/* Casting di WorkMat */
	FullSubMatrixHandler& WMA = WorkMatA.SetFull();
	FullSubMatrixHandler& WMB = WorkMatB.SetFull();

	Vec3 GravityAcceleration;
	bool g = GravityOwner::bGetGravity(pNode->GetXCurr(),
		GravityAcceleration);

	const integer iNumRows = 12;

	/* Dimensiona e resetta la matrice di lavoro */
	WMA.ResizeReset(iNumRows, 12);
	WMB.ResizeReset(iNumRows, 12);

	/* Setta gli indici della matrice - le incognite sono ordinate come:
	 *   - posizione (3)
	 *   - parametri di rotazione (3)
	 *   - quantita' di moto (3)
	 *   - momento della quantita' di moto
	 * e gli indici sono consecutivi. La funzione pGetFirstPositionIndex()
	 * ritorna il valore del primo indice -1, in modo che l'indice i-esimo
	 * e' dato da iGetFirstPositionIndex()+i */
	integer iFirstIndexModal = pNode->iGetFirstIndex();
	for (integer iCnt = 1; iCnt <= 12; iCnt++) {
		WMA.PutRowIndex(iCnt, iFirstIndexModal + iCnt);
		WMA.PutColIndex(iCnt, iFirstIndexModal + iCnt);

		WMB.PutRowIndex(iCnt, iFirstIndexModal + iCnt);
		WMB.PutColIndex(iCnt, iFirstIndexModal + iCnt);
	}

	AssMats(WMA, WMB, 1., XCurr, XPrimeCurr, g, GravityAcceleration);
}


void
ModalBody::AssMats(FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr,
	bool bGravity,
	const Vec3& GravityAcceleration)
{
	DEBUGCOUTFNAME("ModalBody::AssMats");

	const Vec3& W(pNode->GetWCurr());

	Vec3 Sc(STmp*dCoef);
        
	const Mat3x3& RRef = pNode->GetRRef();
	const Mat3x3& RCurr = pNode->GetRCurr();

	const Mat3x3 J12A = (Mat3x3(MatCross, WP) + Mat3x3(MatCrossCross, W, W)).MulVCross(RRef * S0 * (-dCoef));
	const Mat3x3 J13B(Mat3x3DEye, dMass);
	const Mat3x3 J14A = (Mat3x3(MatCrossCross, W, STmp) + Mat3x3(MatCross, W.Cross(STmp))) * (-dCoef);
	const Mat3x3 J14B(MatCross, -STmp);
	const Mat3x3 J22A = (Mat3x3(MatCrossCross, W, RRef * J0 * RCurr.MulTV(W))
		- Mat3x3(MatCrossCross, XPP, RRef * S0)
		- W.Cross(RCurr * J0 * RRef.MulTVCross(W))
		+ Mat3x3(MatCross, RRef * J0 * RCurr.MulTV(WP))
		- RCurr * J0 * RRef.MulTVCross(WP))* (-dCoef);
	const Mat3x3 J23B(MatCross, STmp);
	const Mat3x3 J24A = (Mat3x3(MatCross, JTmp * W) - W.Cross(JTmp)) * (-dCoef);
	const Mat3x3 J24B = JTmp;

	WMA.Add(6 + 1, 3 + 1, J12A);
	WMB.Add(6 + 1, 6 + 1, J13B);
	WMA.Add(6 + 1, 9 + 1, J14A);
	WMB.Add(6 + 1, 9 + 1, J14B);
	WMA.Add(9 + 1, 3 + 1, J22A);
	WMB.Add(9 + 1, 6 + 1, J23B);
	WMA.Add(9 + 1, 9 + 1, J24A);
	WMB.Add(9 + 1, 9 + 1, J24B);

	if (bGravity) {
		WMA.Sub(9 + 1, 3 + 1, Mat3x3(MatCrossCross, GravityAcceleration, Sc));
	}

	const RigidBodyKinematics *pRBK = pNode->pGetRBK();
	if (pRBK) {
		AssMatsRBK_int(WMA, WMB, dCoef, Sc);
	}
}


SubVectorHandler&
ModalBody::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("ModalBody::AssRes");

	/* Se e' definita l'accelerazione di gravita',
	 * la aggiunge (solo al residuo) */
	Vec3 GravityAcceleration;
	bool g = GravityOwner::bGetGravity(pNode->GetXCurr(),
		GravityAcceleration);

	const RigidBodyKinematics *pRBK = pNode->pGetRBK();

	const integer iNumRows = 12;

	WorkVec.ResizeReset(iNumRows);

	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	for (integer iCnt = 1; iCnt <= iNumRows; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
	}

	const Vec3& W(pNode->GetWCurr());

	/* Aggiorna i suoi dati (saranno pronti anche per AssJac) */
	const Mat3x3& R(pNode->GetRCurr());
	STmp = R*S0;
	JTmp = R*J0.MulMT(R);

	const integer iFirstIndexModal = pNode->iGetFirstIndex();

	for (integer i = 1; i <= 3; ++i) {
		XPP(i) = XPrimeCurr.dGetCoef(iFirstIndexModal + i + 6);
		WP(i) = XPrimeCurr.dGetCoef(iFirstIndexModal + i + 9);
	}

	const Vec3 F = XPP * -dMass - WP.Cross(STmp) - W.Cross(W.Cross(STmp));
	const Vec3 M = -STmp.Cross(XPP) - W.Cross(JTmp * W) - JTmp * WP;

	WorkVec.Add(6 + 1, F);
	WorkVec.Add(9 + 1, M);
        
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

void
ModalBody::SetValue(DataManager *pDM,
                    VectorHandler& X, VectorHandler& XP,
                    SimulationEntity::Hints *ph)
{
        // Attention! We must overwrite DynamicBody::SetValue here!
        // Depending on the order of our elements,
        // it is possible to get an incorrect
        // initial velocity and initial angular velocity!
}

/* ModalBody - end */


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
	DEBUGCOUTFNAME("StaticBody::AssRes");
	
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
		M -= STmp.Cross(GravityAcceleration);
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

	/* associated node */
	const StructDispNode *pStrDispNode = pDM->ReadNode<const StructDispNode, Node::STRUCTURAL>(HP);
	const StructNode *pStrNode = dynamic_cast<const StructNode *>(pStrDispNode);

	/* may be determined by a special DataManager parameter... */
	bool bStaticModel = pDM->bIsStaticModel();
	bool bInverseDynamics = pDM->bIsInverseDynamics();

	if (HP.IsKeyWord("variable" "mass")) {
		return ReadVariableBody(pDM, HP, uLabel, pStrNode);
	}

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
		 * J = Sum(J_i) - Sum(dm_i*(Xgc_i-Xgc)/\(Xgc_i-Xgc)/\)
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
	ReferenceFrame RF(pStrDispNode);
	Vec3 Xgc(::Zero3);
	Vec3 STmp(::Zero3);
	Mat3x3 J(::Zero3x3);
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

		if (pStrNode) {
			/* Center of gravity position */
			Vec3 XgcTmp(HP.GetPosRel(RF));
			if (iNumMasses == 1) {
				Xgc = XgcTmp;

			} else {
				STmp += XgcTmp*dMTmp;
			}

			DEBUGLCOUT(MYDEBUG_INPUT, "position of mass(" << iCnt
				<< ") center of gravity = " << XgcTmp << std::endl);
			/*
			 * inertia matrix
			 *
			 * Uses the method that reads a generic matrix with keywords for
			 * abbreviated forms:
			 *   - null: empty matrix
			 *   - eye:  identity matrix
			 *   - diag: diagonal matrix, 3 real numbers expected
			 *   - sym:  symmetric matrix, 6 real numbers expected,
			 *           read as row-oriented upper triangular:
			 *           m11, m12, m13,    m22, m23,    m33
			 *   - generic matrix, 9 real numbers expected, read row-oriented:
			 *           m11, m12, m13,    m21, m22, m23,    m31, m32, m33
			 *
			 * It is also assumed that the inertia matrix is expressed with 
			 * respect to the center of gravity of the body, therefore the 
			 * transport term due to the possible offset with respect to the 
			 * node is removed.
			 *
			 */

			Mat3x3 JTmp(HP.GetMatRel(RF));
			DEBUGLCOUT(MYDEBUG_INPUT, "Inertia matrix of mass(" << iCnt
				<< ") = " << JTmp << std::endl);
			if (!JTmp.IsSymmetric()) {
				silent_cerr("Body(" << uLabel << "): "
					"warning, non-symmetric inertia tensor at line " << HP.GetLineData() << std::endl);
			}

			if (HP.IsKeyWord("orientation")) {
				DEBUGLCOUT(MYDEBUG_INPUT,
					"supplied in relative frame" << std::endl);
				Mat3x3 RTmp(HP.GetRotRel(RF));
				JTmp = RTmp*JTmp.MulMT(RTmp);
			} else if (HP.IsKeyWord("inertial")) {
				silent_cerr("Body(" << uLabel << "): "
					"warning, using deprecated keyword \"inertial\" at line " << HP.GetLineData()
					<< " just use \"orientation\" instead" << std::endl);
				DEBUGLCOUT(MYDEBUG_INPUT,
					"supplied in inertial reference frame" << std::endl);
				if (HP.IsKeyWord("node")) {
					NO_OP;
				} else {
					Mat3x3 RTmp(HP.GetRotRel(RF));
					JTmp = RTmp*JTmp.MulMT(RTmp);
				}
			}
			DEBUGLCOUT(MYDEBUG_INPUT,
				"Inertia matrix of mass(" << iCnt << ") "
				"in current frame = " << JTmp << std::endl);

			J += JTmp - Mat3x3(MatCrossCross, XgcTmp, XgcTmp*dMTmp);
		}
	}

	if (!bNegative && dm < 0.) {
		silent_cerr("Body(" << uLabel << "): "
			"negative mass is not allowed at line "
			<< HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (iNumMasses > 1 && pStrNode) {
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

	const DynamicStructDispNode* pDynamicDispNode = 0;
	const StaticStructDispNode* pStaticDispNode = 0;
	const char *sElemName;
	const char *sNodeName;
	if (pStrNode) {
		sElemName = "Body";
		sNodeName = "StructNode";
	} else {
		sElemName = "Mass";
		sNodeName = "StructDispNode";
	}

	if (bStaticModel || bInverseDynamics) {
		/* static */
		pStaticDispNode = dynamic_cast<const StaticStructDispNode *>(pStrDispNode);
		if (pStaticDispNode == 0 || pStaticDispNode->GetStructDispNodeType() != StructDispNode::STATIC) {
			silent_cerr(sElemName << "(" << uLabel << "): "
				"illegal structural node type "
				"for " << sNodeName << "(" << pStrDispNode->GetLabel() << ") "
				"in " << (bStaticModel ? "static model" : "inverse dynamics") << " analysis "
				"at line " << HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else {
		pDynamicDispNode = dynamic_cast<const DynamicStructDispNode*>(pStrDispNode);
		if (pDynamicDispNode == 0 || pDynamicDispNode->GetStructDispNodeType() != StructDispNode::DYNAMIC) {
			silent_cerr(sElemName << "(" << uLabel << "): "
				"illegal structural node type "
				"for " << sNodeName << "(" << pStrDispNode->GetLabel() << ") "
				"at line " << HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	flag fOut = pDM->fReadOutput(HP, Elem::BODY);

	/* Allocazione e costruzione */
	Elem* pEl = NULL;
	if (bStaticModel || bInverseDynamics) {
		StaticBody *pSB = 0;
		StaticMass *pSM = 0;

		/* static */
		if (pStrNode) {
			SAFENEWWITHCONSTRUCTOR(pSB, StaticBody,
				StaticBody(uLabel, dynamic_cast<const StaticStructNode *>(pStaticDispNode),
					dm, Xgc, J, fOut));
			pEl = pSB;

		} else {
			SAFENEWWITHCONSTRUCTOR(pSM, StaticMass,
				StaticMass(uLabel, pStaticDispNode,
					dm, fOut));
			pEl = pSM;
		}

		if (bInverseDynamics) {
			bool bIsRightHandSide(true);
			bool bIsErgonomy(true);

			if (HP.IsKeyWord("inverse" "dynamics")) {
				bIsRightHandSide = false;
				if (HP.IsKeyWord("right" "hand" "side")) {
					bIsRightHandSide = HP.GetYesNoOrBool(bIsRightHandSide);
				}

				bIsErgonomy = false;
				if (HP.IsKeyWord("ergonomy")) {
					bIsErgonomy = HP.GetYesNoOrBool(bIsErgonomy);
				}
			}

			unsigned flags = 0;

			if (bIsRightHandSide) {
				flags |= InverseDynamics::RIGHT_HAND_SIDE;
			}

			if (bIsErgonomy) {
				flags |= InverseDynamics::ERGONOMY;
			}

			if (pSB) {
				pSB->SetInverseDynamicsFlags(flags);

			} else {
				pSM->SetInverseDynamicsFlags(flags);
			}
		}

	} else {
		if (pStrNode) {
			const DynamicStructNode* pDynamicStructNode = dynamic_cast<const DynamicStructNode*>(pDynamicDispNode);
			const ModalNode* pModalNode = dynamic_cast<const ModalNode*>(pDynamicDispNode);
			const RigidBodyKinematics* pRBK = pDynamicStructNode->pGetRBK();

			if (pModalNode && pRBK) {
				silent_cerr("Body(" << uLabel << ") "
					"is connected to ModalNode(" << pModalNode->GetLabel() << ") "
					"which uses rigid body kinematics "
					"at line " << HP.GetLineData() << std::endl);
				throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
			}

			if (pModalNode) {
				SAFENEWWITHCONSTRUCTOR(pEl, ModalBody,
					ModalBody(uLabel, pModalNode, dm, Xgc, J, fOut));

			} else {
				SAFENEWWITHCONSTRUCTOR(pEl, DynamicBody,
					DynamicBody(uLabel, pDynamicStructNode, dm, Xgc, J, fOut));
			}
		} else {
			SAFENEWWITHCONSTRUCTOR(pEl, DynamicMass,
				DynamicMass(uLabel, pDynamicDispNode,
					dm, fOut));
		}
	}

	pDM->GetLogFile()
		<< "body: " << uLabel
		<< ' ' << pStrDispNode->GetLabel()
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

