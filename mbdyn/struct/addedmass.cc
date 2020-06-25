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

/* added mass element */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <limits>

#include "addedmass.h"
//#include "body_vm.h"
#include "dataman.h"

/* AddedMass - begin */

AddedMass::AddedMass(unsigned int uL,
	const StructDispNode *pNode,
	Vec3 AddedMassValues,
	flag fOut)
: Elem(uL, fOut),
InitialAssemblyElem(uL, fOut),
pNode(pNode),
AddedMassValues(AddedMassValues)
{
	ASSERT(pNode != NULL);
	ASSERT(pNode->GetNodeType() == Node::STRUCTURAL);

	for (int i=0 ; i < 3 ; i++)
	{
        ASSERT(AddedMassValues[i] > 0.);
	}

	AddedMassDiagMat = Mat3x3(Mat3x3DEye, AddedMassValues);

}


/* distruttore */
AddedMass::~AddedMass(void)
{
	NO_OP;
}


/* momento statico */
Vec3
AddedMass::GetS_int(void) const
{
	return pNode->GetXCurr().EBEMult(AddedMassValues);
}


/* momento d'inerzia */
Mat3x3
AddedMass::GetJ_int(void) const
{
	const Vec3& x = pNode->GetXCurr();

	return Mat3x3(MatCrossCross, x, x.EBEMult(-AddedMassValues));
}


/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
AddedMass::Restart(std::ostream& out) const
{
	out << "  added mass: " << GetLabel() << ", "
		<< pNode->GetLabel() << ", " << AddedMassValues << ';' << std::endl;

	return out;
}


/* total added mass */
Vec3
AddedMass::GetM(void) const
{
	return AddedMassValues;
}

///* momento statico */
//Vec3
//AddedMass::GetS(void) const
//{
//	return GetS_int();
//}

/* momento d'inerzia */
Mat3x3
AddedMass::GetJ(void) const
{
	return GetJ_int();
}

/* nodo */
const StructDispNode *
AddedMass::pGetNode(void) const
{
	return pNode;
}

/* Accesso ai dati privati */
unsigned int
AddedMass::iGetNumPrivData(void) const
{
	return 3;
}

unsigned int
AddedMass::iGetPrivDataIdx(const char *s) const
{
  	ASSERT(s != NULL);

	if (strlen(s) != 2) {
		return 0;
	}

	unsigned int off = 0;

	switch (s[0]) {
	case 'm':
		/* relative position */
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
AddedMass::dGetPrivData(unsigned int i) const
{
	switch (i) {

	case 1:
 		// AddedMass
 		return AddedMassValues[0];

    case 2:
 		// AddedMass
 		return AddedMassValues[1];

    case 3:
 		// AddedMass
 		return AddedMassValues[2];
	}

	return 0.;
}

void
AddedMass::AssVecRBK_int(SubVectorHandler& WorkVec)
{
	const RigidBodyKinematics *pRBK = pNode->pGetRBK();

	Vec3 s0;

	integer iIdx = 0;
	if (dynamic_cast<DynamicAddedMass *>(this)) {
		iIdx = 3;
	}

	// elementwise multiplication of AddedMass and positions
	s0 = AddedMassValues.EBEMult (pNode->GetXCurr());

	// force in each direction
	Vec3 f;
	f = AddedMassValues.EBEMult (pRBK->GetXPP());

	WorkVec.Sub(iIdx + 1, f);
}

void
AddedMass::AssMatsRBK_int(
	FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB,
	const doublereal& dCoef)
{
	const RigidBodyKinematics *pRBK = pNode->pGetRBK();

	integer iIdx = 0;
	if (dynamic_cast<DynamicAddedMass *>(this)) {
		iIdx = 3;
	}

	// f: delta x
	Mat3x3 MTmp(MatCross, pRBK->GetWP());
	MTmp += Mat3x3(MatCrossCross, pRBK->GetW(), pRBK->GetW());

//    Mat3x3 AddedMassDiag = Mat3x3( AddedMassValues[0],         0.0,           0.0,
//                                            0.0, AddedMassValues[1],          0.0,
//                                            0.0,          0.0, AddedMassValues[1]  );

	WMA.Add(iIdx + 1, 1, MTmp*(AddedMassDiagMat*dCoef));
}

/* AddedMass - end */


/* DynamicAddedMass - begin */

DynamicAddedMass::DynamicAddedMass(unsigned int uL,
	const DynamicStructDispNode* pNode,
	Vec3 AddedMassValues,
	flag fOut)
: Elem(uL, fOut),
AddedMass(uL, pNode, AddedMassValues, fOut)
{
	NO_OP;
}


/* distruttore */
DynamicAddedMass::~DynamicAddedMass(void)
{
	NO_OP;
}


VariableSubMatrixHandler&
DynamicAddedMass::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("DynamicAddedMass::AssJac");

	/* Casting di WorkMat */
	FullSubMatrixHandler& WM = WorkMat.SetFull();

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

	AssMats(WM, WM, dCoef);

	return WorkMat;
}


void
DynamicAddedMass::AssMats(VariableSubMatrixHandler& WorkMatA,
	VariableSubMatrixHandler& WorkMatB,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("DynamicAddedMass::AssMats");

	/* Casting di WorkMat */
	FullSubMatrixHandler& WMA = WorkMatA.SetFull();
	FullSubMatrixHandler& WMB = WorkMatB.SetFull();

	integer iNumRows = 3;

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

	AssMats(WMA, WMB, 1.);
}


void
DynamicAddedMass::AssMats(FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB,
	doublereal dCoef)
{
	DEBUGCOUTFNAME("DynamicAddedMass::AssMats");

	/*
	 * momentum:
	 *
	 * m * I DeltaV - S /\ DeltagP + ( S /\ W ) /\ Deltag
	 */
	WMB.IncCoef(1, 1, AddedMassValues[0]);
	WMB.IncCoef(2, 2, AddedMassValues[1]);
	WMB.IncCoef(3, 3, AddedMassValues[2]);

	const RigidBodyKinematics *pRBK = pNode->pGetRBK();
	if (pRBK) {
		AssMatsRBK_int(WMA, WMB, dCoef);
	}
}


SubVectorHandler&
DynamicAddedMass::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUTFNAME("DynamicAddedMass::AssRes");


	const RigidBodyKinematics *pRBK = pNode->pGetRBK();

	integer iNumRows = 3;

	WorkVec.ResizeReset(iNumRows);

	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	for (integer iCnt = 1; iCnt <= iNumRows; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
	}

	const Vec3& V(pNode->GetVCurr());

	/* Quantita' di moto: R[1] = Q */
	WorkVec.Sub(1, V.EBEMult(AddedMassValues));

	if (pRBK) {
		AssVecRBK_int(WorkVec);
	}

	const DynamicStructDispNode *pDN = dynamic_cast<const DynamicStructDispNode *>(pNode);
	ASSERT(pDN != 0);

	pDN->AddInertia(AddedMassValues); // TODO: what do we do here?

	return WorkVec;
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
DynamicAddedMass::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUTFNAME("DynamicAddedMass::InitialAssJac");

	/* Casting di WorkMat */
	WorkMat.SetNullMatrix();

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
DynamicAddedMass::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUTFNAME("DynamicAddedMass::InitialAssRes");

	WorkVec.Resize(0);

	return WorkVec;
}


/* Usata per inizializzare la quantita' di moto */
void
DynamicAddedMass::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& /* XP */ ,
	SimulationEntity::Hints *ph)
{
	integer iFirstIndex = pNode->iGetFirstMomentumIndex();

	const Vec3& V(pNode->GetVCurr());
	X.Add(iFirstIndex + 1, V.EBEMult(AddedMassValues));
}

/* momentum */
Vec3
DynamicAddedMass::GetB_int(void) const
{
	const Vec3& V(pNode->GetVCurr());

	return V.EBEMult(AddedMassValues);
}

/* DynamicAddedMass - end */


/* StaticAddedMass - begin */

StaticAddedMass::StaticAddedMass(unsigned int uL,
	const StaticStructDispNode* pNode,
	Vec3 AddedMassValues,
	flag fOut)
: Elem(uL, fOut),
AddedMass(uL, pNode, AddedMassValues, fOut)
{
	NO_OP;
}


/* distruttore */
StaticAddedMass::~StaticAddedMass(void)
{
	NO_OP;
}


VariableSubMatrixHandler&
StaticAddedMass::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("StaticAddedMass::AssJac");

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
StaticAddedMass::AssMats(VariableSubMatrixHandler& WorkMatA,
	VariableSubMatrixHandler& WorkMatB,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("StaticAddedMass::AssMats");

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
StaticAddedMass::AssMats(FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB,
	doublereal dCoef)
{
	DEBUGCOUTFNAME("StaticAddedMass::AssMats");

	/* Se e' definita l'accelerazione di gravita',
	 * la aggiunge (solo al residuo) */
	Vec3 Acceleration(Zero3);
	const RigidBodyKinematics *pRBK = pNode->pGetRBK();

	if (!pRBK) {
		/* Caller will set WMA & WMB to null matrix */
		return true;
	}

	if (pRBK) {
		AssMatsRBK_int(WMA, WMB, dCoef);
	}

	return false;
}


SubVectorHandler&
StaticAddedMass::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUTFNAME("StaticAddedMass::AssRes");

	/* Se e' definita l'accelerazione di gravita',
	 * la aggiunge (solo al residuo) */
	Vec3 Acceleration(Zero3);

	/* W is uninitialized because its use is conditioned by w */
	const RigidBodyKinematics *pRBK = pNode->pGetRBK();

	if (!pRBK) {
		WorkVec.Resize(0);
		return WorkVec;
	}

	WorkVec.ResizeReset(3);

	integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex + iCnt);
	}

	if (pRBK) {
		AssVecRBK_int(WorkVec);
	}

	return WorkVec;
}


/* inverse dynamics capable element */
bool
StaticAddedMass::bInverseDynamics(void) const
{
	return true;
}


SubVectorHandler&
StaticAddedMass::AssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ ,
	const VectorHandler& /* XPrimePrimeCurr */ ,
	InverseDynamics::Order iOrder)
{
	DEBUGCOUTFNAME("DynamicAddedMass::AssRes");

	ASSERT(iOrder == InverseDynamics::INVERSE_DYNAMICS);

	WorkVec.ResizeReset(3);

	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
	}

	Vec3 Acceleration = pNode->GetXPPCurr();

	WorkVec.Sub(1, Acceleration.EBEMult(AddedMassValues));

	return WorkVec;
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
StaticAddedMass::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUTFNAME("StaticAddedMass::InitialAssJac");

	WorkMat.SetNullMatrix();

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
StaticAddedMass::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUTFNAME("StaticAddedMass::InitialAssRes");

	WorkVec.Resize(0);

	return WorkVec;
}


/* Usata per inizializzare la quantita' di moto */
void
StaticAddedMass::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& /* XP */ ,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

/* StaticAddedMass - end */


/* AddedMassAndInertia - begin */

AddedMassAndInertia::AddedMassAndInertia(unsigned int uL,
	const StructNode *pNode,
	Vec3 AddedMassValues,
	const Vec3& Xgc,
	const Mat3x3& J,
	flag fOut)
: Elem(uL, fOut),
InitialAssemblyElem(uL, fOut),
pNode(pNode),
AddedMassValues(AddedMassValues),
Xgc(Xgc),
S0(Xgc.EBEMult(AddedMassValues)),
J0(J)
{
	ASSERT(pNode != NULL);
	ASSERT(pNode->GetNodeType() == Node::STRUCTURAL);
	for (int i=0 ; i < 3 ; i++)
	{
        ASSERT(AddedMassValues[i] > 0.);
	}

	AddedMassDiagMat = Mat3x3(Mat3x3DEye, AddedMassValues);
}


/* distruttore */
AddedMassAndInertia::~AddedMassAndInertia(void)
{
	NO_OP;
}


/* momento statico */
Vec3
AddedMassAndInertia::GetS_int(void) const
{
	return pNode->GetXCurr().EBEMult(AddedMassValues) + pNode->GetRCurr()*S0;
}


/* momento d'inerzia */
Mat3x3
AddedMassAndInertia::GetJ_int(void) const
{
	Vec3 s = pNode->GetRCurr()*S0;
	const Vec3& x = pNode->GetXCurr();

	return pNode->GetRCurr()*J0.MulMT(pNode->GetRCurr())
		- Mat3x3(MatCrossCross, x, x.EBEMult(AddedMassValues))
		- Mat3x3(MatCrossCross, s, x)
		- Mat3x3(MatCrossCross, x, s);
}


/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
AddedMassAndInertia::Restart(std::ostream& out) const
{
	out << "  added mass: " << GetLabel() << ", "
		<< pNode->GetLabel() << ", " << AddedMassValues << ", "
		<< "reference, node, ", Xgc.Write(out, ", ") << ", "
		<< "reference, node, ", (J0 + Mat3x3(MatCrossCross, S0, Xgc)).Write(out, ", ")
		<< ";" << std::endl;

	return out;
}


void
AddedMassAndInertia::AfterPredict(VectorHandler& /* X */ , VectorHandler& /* XP */ )
{
	const Mat3x3& R = pNode->GetRRef();

	STmp = R*S0;
	JTmp = R*J0.MulMT(R);
}

/* massa totale */
Vec3
AddedMassAndInertia::GetM(void) const
{
	return AddedMassValues;
}

/* momento statico */
Vec3
AddedMassAndInertia::GetS(void) const
{
	return GetS_int();
}

/* momento d'inerzia */
Mat3x3
AddedMassAndInertia::GetJ(void) const
{
	return GetJ_int();
}

/* nodo */
const StructNode *
AddedMassAndInertia::pGetNode(void) const
{
	return pNode;
}

/* Accesso ai dati privati */
unsigned int
AddedMassAndInertia::iGetNumPrivData(void) const
{
	return 3;
}

unsigned int
AddedMassAndInertia::iGetPrivDataIdx(const char *s) const
{
  	ASSERT(s != NULL);

	if (strlen(s) != 2) {
		return 0;
	}

	unsigned int off = 0;

	switch (s[0]) {
	case 'm':
		/* relative position */
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
AddedMassAndInertia::dGetPrivData(unsigned int i) const
{
	switch (i) {
	case 1:
		return AddedMassValues[0];
    case 2:
		return AddedMassValues[1];
	case 3:
		return AddedMassValues[2];
	}

	return 0.;
}

void
AddedMassAndInertia::AssVecRBK_int(SubVectorHandler& WorkVec)
{
	const RigidBodyKinematics *pRBK = pNode->pGetRBK();

	Vec3 s0;

	integer iIdx = 0;
	if (dynamic_cast<DynamicAddedMassAndInertia *>(this)) {
		iIdx = 6;
	}

	s0 = AddedMassValues.EBEMult(pNode->GetXCurr()) + STmp;

	// force
	Vec3 f;
	f = AddedMassValues.EBEMult(pRBK->GetXPP());
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
AddedMassAndInertia::AssMatsRBK_int(
	FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB,
	const doublereal& dCoef,
	const Vec3& Sc)
{
	const RigidBodyKinematics *pRBK = pNode->pGetRBK();

	Mat3x3 MTmp;
	Vec3 VTmp;

	integer iIdx = 0;
	if (dynamic_cast<DynamicAddedMassAndInertia *>(this)) {
		iIdx = 6;
	}

	// f: delta x
	MTmp = Mat3x3(MatCross, pRBK->GetWP());
	MTmp += Mat3x3(MatCrossCross, pRBK->GetW(), pRBK->GetW());

	WMA.Add(iIdx + 1, 1, MTmp*(AddedMassDiagMat*dCoef));


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

/* AddedMassAndInertia - end */


/* DynamicAddedMassAndInertia - begin */

DynamicAddedMassAndInertia::DynamicAddedMassAndInertia(unsigned int uL,
	const DynamicStructNode* pNode,
	Vec3 AddedMassValues,
	const Vec3& Xgc,
	const Mat3x3& J,
	flag fOut)
: Elem(uL, fOut),
AddedMassAndInertia(uL, pNode, AddedMassValues, Xgc, J, fOut)
{
	NO_OP;
}


/* distruttore */
DynamicAddedMassAndInertia::~DynamicAddedMassAndInertia(void)
{
	NO_OP;
}

void
DynamicAddedMassAndInertia::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 12;
	*piNumCols = 6;
}

void
DynamicAddedMassAndInertia::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 12;
	*piNumCols = 6;
}

VariableSubMatrixHandler&
DynamicAddedMassAndInertia::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("DynamicAddedMassAndInertia::AssJac");

	/* Casting di WorkMat */
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	integer iNumRows = 6;
	if (pNode->pGetRBK()) {
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

	AssMats(WM, WM, dCoef);

	return WorkMat;
}


void
DynamicAddedMassAndInertia::AssMats(VariableSubMatrixHandler& WorkMatA,
	VariableSubMatrixHandler& WorkMatB,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("DynamicAddedMassAndInertia::AssMats");

	/* Casting di WorkMat */
	FullSubMatrixHandler& WMA = WorkMatA.SetFull();
	FullSubMatrixHandler& WMB = WorkMatB.SetFull();

	integer iNumRows = 6;

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

	AssMats(WMA, WMB, 1.);
}


void
DynamicAddedMassAndInertia::AssMats(FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB,
	doublereal dCoef)
{
	DEBUGCOUTFNAME("DynamicAddedMassAndInertia::AssMats");

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
	WMB.IncCoef(1, 1, AddedMassValues[0]);
	WMB.IncCoef(2, 2, AddedMassValues[1]);
	WMB.IncCoef(3, 3, AddedMassValues[2]);

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

	const RigidBodyKinematics *pRBK = pNode->pGetRBK();
	if (pRBK) {
		AssMatsRBK_int(WMA, WMB, dCoef, Sc);
	}
}


SubVectorHandler&
DynamicAddedMassAndInertia::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUTFNAME("DynamicAddedMassAndInertia::AssRes");

	const RigidBodyKinematics *pRBK = pNode->pGetRBK();

	integer iNumRows = 6;
	if (pRBK) {
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
	WorkVec.Sub(1, V.EBEMult(AddedMassValues) + W.Cross(STmp));

	/* Momento della quantita' di moto: R[2] = G - S /\ V - J * W */
	WorkVec.Sub(3 + 1, JTmp*W + STmp.Cross(V));

	if (pRBK) {
		AssVecRBK_int(WorkVec);
	}

	const DynamicStructNode *pDN = dynamic_cast<const DynamicStructNode *>(pNode);
	ASSERT(pDN != 0);

	pDN->AddInertia(AddedMassValues, STmp, JTmp);

	return WorkVec;
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
DynamicAddedMassAndInertia::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUTFNAME("DynamicAddedMassAndInertia::InitialAssJac");

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
DynamicAddedMassAndInertia::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUTFNAME("DynamicAddedMassAndInertia::InitialAssRes");

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

	return WorkVec;
}


/* Usata per inizializzare la quantita' di moto */
void
DynamicAddedMassAndInertia::SetValue(DataManager *pDM,
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
	X.Add(iFirstIndex + 1, V.EBEMult(AddedMassValues) + W.Cross(STmp));
	X.Add(iFirstIndex + 4, STmp.Cross(V) + JTmp*W);
}

/* momentum */
Vec3
DynamicAddedMassAndInertia::GetB_int(void) const
{
	const Vec3& V(pNode->GetVCurr());
	const Vec3& W(pNode->GetWCurr());
	const Mat3x3& R(pNode->GetRCurr());

	return V.EBEMult(AddedMassValues) + W.Cross(R*S0);
}


/* momenta moment */
Vec3
DynamicAddedMassAndInertia::GetG_int(void) const
{
	const Vec3& X(pNode->GetXCurr());
	const Mat3x3& R(pNode->GetRCurr());
	const Vec3& V(pNode->GetVCurr());
	const Vec3& W(pNode->GetWCurr());

	Vec3 STmp(R*S0);

	// NOTE: with respect to the origin of the global reference frame!
	return (STmp + X.EBEMult(AddedMassValues)).Cross(V) + R*(J0*(R.MulTV(W)))
		- X.Cross(STmp.Cross(W));
}


/* DynamicAddedMassAndInertia - end */


/* ModalAddedMassAndInertia - begin */

ModalAddedMassAndInertia::ModalAddedMassAndInertia(unsigned int uL,
	const ModalNode* pNode,
	Vec3 AddedMassValues,
	const Vec3& Xgc,
	const Mat3x3& J,
	flag fOut)
: Elem(uL, fOut),
DynamicAddedMassAndInertia(uL, pNode, AddedMassValues, Xgc, J, fOut),
XPP(::Zero3), WP(::Zero3)
{
	NO_OP;
}

/* distruttore */
ModalAddedMassAndInertia::~ModalAddedMassAndInertia(void)
{
	NO_OP;
}

void
ModalAddedMassAndInertia::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 12;
	*piNumCols = 12;
}

VariableSubMatrixHandler&
ModalAddedMassAndInertia::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("ModalAddedMassAndInertia::AssJac");

	/* Casting di WorkMat */
	FullSubMatrixHandler& WM = WorkMat.SetFull();

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

	AssMats(WM, WM, dCoef, XCurr, XPrimeCurr);

	return WorkMat;
}


void
ModalAddedMassAndInertia::AssMats(VariableSubMatrixHandler& WorkMatA,
	VariableSubMatrixHandler& WorkMatB,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("ModalAddedMassAndInertia::AssMats");

	/* Casting di WorkMat */
	FullSubMatrixHandler& WMA = WorkMatA.SetFull();
	FullSubMatrixHandler& WMB = WorkMatB.SetFull();

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

	AssMats(WMA, WMB, 1., XCurr, XPrimeCurr);
}


void
ModalAddedMassAndInertia::AssMats(FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("ModalAddedMassAndInertia::AssMats");

	const Vec3& W(pNode->GetWCurr());

	Vec3 Sc(STmp*dCoef);

	const Mat3x3& RRef = pNode->GetRRef();
	const Mat3x3& RCurr = pNode->GetRCurr();

	const Mat3x3 J12A = (Mat3x3(MatCross, WP) + Mat3x3(MatCrossCross, W, W)).MulVCross(RRef * S0 * (-dCoef));
	const Mat3x3 J13B(Mat3x3DEye, AddedMassValues);
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

	const RigidBodyKinematics *pRBK = pNode->pGetRBK();
	if (pRBK) {
		AssMatsRBK_int(WMA, WMB, dCoef, Sc);
	}
}


SubVectorHandler&
ModalAddedMassAndInertia::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("ModalAddedMassAndInertia::AssRes");

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

	const Vec3 F = XPP.EBEMult(-AddedMassValues) - WP.Cross(STmp) - W.Cross(W.Cross(STmp));
	const Vec3 M = -STmp.Cross(XPP) - W.Cross(JTmp * W) - JTmp * WP;

	WorkVec.Add(6 + 1, F);
	WorkVec.Add(9 + 1, M);

	if (pRBK) {
		AssVecRBK_int(WorkVec);
	}

	const DynamicStructNode *pDN = dynamic_cast<const DynamicStructNode *>(pNode);
	ASSERT(pDN != 0);

	pDN->AddInertia(AddedMassValues, STmp, JTmp);

	return WorkVec;
}

/* ModalAddedMassAndInertia - end */


/* StaticAddedMassAndInertia - begin */

StaticAddedMassAndInertia::StaticAddedMassAndInertia(unsigned int uL,
	const StaticStructNode* pNode,
	Vec3 AddedMassValues,
	const Vec3& Xgc,
	const Mat3x3& J,
	flag fOut)
: Elem(uL, fOut),
AddedMassAndInertia(uL, pNode, AddedMassValues, Xgc, J, fOut)
{
	NO_OP;
}


/* distruttore */
StaticAddedMassAndInertia::~StaticAddedMassAndInertia(void)
{
	NO_OP;
}


VariableSubMatrixHandler&
StaticAddedMassAndInertia::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("StaticAddedMassAndInertia::AssJac");

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
StaticAddedMassAndInertia::AssMats(VariableSubMatrixHandler& WorkMatA,
	VariableSubMatrixHandler& WorkMatB,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("StaticAddedMassAndInertia::AssMats");

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
StaticAddedMassAndInertia::AssMats(FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB,
	doublereal dCoef)
{
	DEBUGCOUTFNAME("StaticAddedMassAndInertia::AssMats");

	/* TODO: reference */
	Vec3 W(Zero3);

	const RigidBodyKinematics *pRBK = pNode->pGetRBK();

	if (!pRBK) {
		/* Caller will set WMA & WMB to null matrix */
		return true;
	}

	Vec3 Sc(STmp*dCoef);

	if (pRBK) {
		AssMatsRBK_int(WMA, WMB, dCoef, Sc);
	}

	return false;
}


SubVectorHandler&
StaticAddedMassAndInertia::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUTFNAME("StaticAddedMassAndInertia::AssRes");

	/* W is uninitialized because its use is conditioned by w */
	const RigidBodyKinematics *pRBK = pNode->pGetRBK();

	if (!pRBK) {
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

	if (pRBK) {
		AssVecRBK_int(WorkVec);
	}

	return WorkVec;
}


/* inverse dynamics capable element */
bool
StaticAddedMassAndInertia::bInverseDynamics(void) const
{
	return true;
}


SubVectorHandler&
StaticAddedMassAndInertia::AssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ ,
	const VectorHandler& /* XPrimePrimeCurr */ ,
	InverseDynamics::Order iOrder)
{
	DEBUGCOUTFNAME("StaticAddedMassAndInertia::AssRes");

	ASSERT(iOrder == InverseDynamics::INVERSE_DYNAMICS);

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

	WorkVec.Sub(1, Acceleration.EBEMult(AddedMassValues));

	Vec3 M = JTmp*pNode->GetWPCurr()
		+ STmp.Cross(pNode->GetXPPCurr())
		+ STmp.Cross(pNode->GetWCurr().Cross(pNode->GetWCurr().Cross(XgcTmp)));

	WorkVec.Sub(4, M);

	return WorkVec;
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
StaticAddedMassAndInertia::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUTFNAME("StaticAddedMassAndInertia::InitialAssJac");

	WorkMat.SetNullMatrix();

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
StaticAddedMassAndInertia::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUTFNAME("StaticAddedMassAndInertia::InitialAssRes");

	WorkVec.ResizeReset(0);

	return WorkVec;
}


/* Usata per inizializzare la quantita' di moto */
void
StaticAddedMassAndInertia::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& /* XP */ ,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

/* StaticAddedMassAndInertia - end */


/* Legge un corpo rigido */
Elem*
ReadAddedMass(DataManager* pDM, MBDynParser& HP, unsigned int uLabel)
{
	DEBUGCOUTFNAME("ReadAddedMass");

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
	const StructDispNode *pStrDispNode = pDM->ReadNode<const StructDispNode, Node::STRUCTURAL>(HP);
	const StructNode *pStrNode = dynamic_cast<const StructNode *>(pStrDispNode);

	/* may be determined by a special DataManager parameter... */
	bool bStaticModel = pDM->bIsStaticModel();
	bool bInverseDynamics = pDM->bIsInverseDynamics();
//
//	if (HP.IsKeyWord("variable" "added" "mass")) {
//		return ReadVariableAddedMass(pDM, HP, uLabel, pStrNode);
//	}

//	integer iNumMasses = 1;
//	if (HP.IsKeyWord("condense")) {
//		iNumMasses = HP.GetInt();
//		if (iNumMasses < 1) {
//			silent_cerr("AddedMassAndInertia(" << uLabel << "): "
//				"at least one mass is required in \"condense\" "
//				"at line " << HP.GetLineData() << std::endl);
//			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
//		}
//		DEBUGLCOUT(MYDEBUG_INPUT,
//			iNumMasses << " masses will be condensed" << std::endl);
//
//		/* The inertia is calculated as follows:
//		 *
//		 * AM = Sum(dm_i)
//		 *
//		 * Xgc = Sum(Xgc_i*dm_i)/Sum(dm_i)
//		 *
//		 * J = Sum(J_i)-Sum(dm_i*(Xgc_i-Xgc)/\(Xgc_i-Xgc)/\)
//		 *
//		 * and it can be accomplished by accumulating:
//		 *
//		 * AM = Sum(dm_i)
//		 *
//		 * ~S = Sum(Xgc_i*dm_i)
//		 *
//		 * ~J = Sum(J_i)-Sum(dm_i*Xgc_i/\*Xgc_i/\)
//		 *
//		 * then calculating
//		 *
//		 * Xgc = S/AM
//		 *
//		 * and finally:
//		 *
//		 * J = ~J-Xgc/\(AM*Xgc/\-2*~S)
//		 *
//		 */
//	}

	ReferenceFrame RF(pStrDispNode);
	Vec3 Xgc(::Zero3);
	Vec3 STmp(::Zero3);
	Mat3x3 J(::Zero3x3);
	bool bNegative(false);

	if (HP.IsKeyWord("allow" "negative" "mass")) {
		bNegative = true;
	}

//	for (int iCnt = 1; iCnt <= iNumMasses; iCnt++) {
    /* added mass in x y and z */
    Vec3 AM = HP.GetVec3();

    if (!bNegative && (AM[0] < 0. || AM[1] < 0. || AM[2] < 0.)) {
        silent_cerr("AddedMass(" << uLabel << "): "
            "negative added mass is not allowed at line "
            << HP.GetLineData() << std::endl);
        throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
    }

    if (pStrNode) {
        /* posiz. c.g. */
        Vec3 XgcTmp(HP.GetPosRel(RF));
        Xgc = XgcTmp;

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
        // TODO: instead of "inertial" or "node", use "orientation [ , reference, <ref> ] , <mat>"...
        Mat3x3 JTmp(HP.GetMatRel(RF));
        DEBUGLCOUT(MYDEBUG_INPUT, "Added inertia matrix of added mass(" << iCnt
            << ") =" << std::endl << JTmp << std::endl);
        if (!JTmp.IsSymmetric()) {
            silent_cerr("AddedMass(" << uLabel << "): "
                "warning, non-symmetric inertia tensor at line " << HP.GetLineData() << std::endl);
        }

//        if (HP.IsKeyWord("inertial")) {
//            DEBUGLCOUT(MYDEBUG_INPUT,
//                "supplied in inertial reference frame" << std::endl);
//            if (HP.IsKeyWord("node")) {
//                NO_OP;
//            } else {
//                Mat3x3 RTmp(HP.GetRotRel(RF));
//                JTmp = RTmp*JTmp.MulMT(RTmp);
//            }
//            DEBUGLCOUT(MYDEBUG_INPUT,
//                "Inertia matrix of mass(" << iCnt << ") "
//                "in current frame =" << JTmp << std::endl);
//        }

        J += JTmp - Mat3x3(MatCrossCross, XgcTmp, XgcTmp.EBEMult(AM));
    }
//	}

//	if (!bNegative && AM < 0.) {
//		silent_cerr("AddedMassAndInertia(" << uLabel << "): "
//			"negative mass is not allowed at line "
//			<< HP.GetLineData() << std::endl);
//		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
//	}

	const DynamicStructDispNode* pDynamicDispNode = 0;
	const StaticStructDispNode* pStaticDispNode = 0;
	const char *sElemName;
	const char *sNodeName;
	sElemName = "AddedMass";
	if (pStrNode) {
        sElemName = "AddedMassAndInertia";
		sNodeName = "StructNode";
	} else {
        sElemName = "AddedMass";
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

	flag fOut = pDM->fReadOutput(HP, Elem::ADDEDMASS);

	/* Allocazione e costruzione */
	Elem* pEl = NULL;
	if (bStaticModel || bInverseDynamics) {
		StaticAddedMassAndInertia *pSB = 0;
		StaticAddedMass *pSM = 0;

		/* static */
		if (pStrNode) {
			SAFENEWWITHCONSTRUCTOR(pSB, StaticAddedMassAndInertia,
				StaticAddedMassAndInertia(uLabel, dynamic_cast<const StaticStructNode *>(pStaticDispNode),
					AM, Xgc, J, fOut));
			pEl = pSB;

		} else {
			SAFENEWWITHCONSTRUCTOR(pSM, StaticAddedMass,
				StaticAddedMass(uLabel, pStaticDispNode,
					AM, fOut));
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
				silent_cerr("AddedMass(" << uLabel << ") "
					"is connected to ModalNode(" << pModalNode->GetLabel() << ") "
					"which uses rigid body kinematics "
					"at line " << HP.GetLineData() << std::endl);
				throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
			}

			if (pModalNode) {
				SAFENEWWITHCONSTRUCTOR(pEl, ModalAddedMassAndInertia,
					ModalAddedMassAndInertia(uLabel, pModalNode, AM, Xgc, J, fOut));

			} else {
				SAFENEWWITHCONSTRUCTOR(pEl, DynamicAddedMassAndInertia,
					DynamicAddedMassAndInertia(uLabel, pDynamicStructNode, AM, Xgc, J, fOut));
			}
		} else {
			SAFENEWWITHCONSTRUCTOR(pEl, DynamicAddedMass,
				DynamicAddedMass(uLabel, pDynamicDispNode,
					AM, fOut));
		}
	}

	pDM->GetLogFile()
		<< "added mass: " << uLabel
		<< ' ' << pStrDispNode->GetLabel()
		<< ' ' << AM
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
} /* End of ReadAddedMass() */

