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

#include "body_vm.h"
#include "dataman.h"

/* VariableBody - begin */

VariableBody::VariableBody(unsigned int uL,
	const StructNode *pNode,
	const DriveCaller *pDCMass,
	const TplDriveCaller<Vec3> *pDCXgc,
	const TplDriveCaller<Mat3x3> *pDCJgc,
	flag fOut)
: Elem(uL, fOut),
ElemGravityOwner(uL, fOut),
InitialAssemblyElem(uL, fOut),
pNode(pNode),
m_Mass(pDCMass),
m_Xgc(pDCXgc),
m_Jgc(pDCJgc)
{
	ASSERT(pNode != 0);
	ASSERT(pNode->GetNodeType() == Node::STRUCTURAL);
}


/* distruttore */
VariableBody::~VariableBody(void)
{
	NO_OP;
}


/* momento statico */
Vec3
VariableBody::GetS_int(void) const
{
	return (pNode->GetXCurr() + pNode->GetRCurr()*m_Xgc.Get())*m_Mass.dGet();
}


/* momento d'inerzia */
Mat3x3
VariableBody::GetJ_int(void) const
{
	Vec3 x = pNode->GetXCurr() + pNode->GetRCurr()*m_Xgc.Get();

	return pNode->GetRCurr()*m_Jgc.Get().MulMT(pNode->GetRCurr())
		- Mat3x3(MatCrossCross, x, x*m_Mass.dGet());
}


/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
VariableBody::Restart(std::ostream& out) const
{
	out << "  body: " << GetLabel() << ", "
		<< pNode->GetLabel() << ", variable mass, ", m_Mass.pGetDriveCaller()->Restart(out) << ", "
		<< "reference, node, ", m_Xgc.pGetDriveCaller()->Restart(out) << ", "
		<< "reference, node, ", m_Jgc.pGetDriveCaller()->Restart(out) << ";" << std::endl;

	return out;
}


void
VariableBody::AfterPredict(VectorHandler& /* X */ , VectorHandler& /* XP */ )
{
	const Mat3x3& R = pNode->GetRRef();

	dMTmp = m_Mass.dGet();
	Vec3 x(m_Xgc.Get());
	STmp = (R*x)*dMTmp;
	JTmp = R*(m_Jgc.Get() + Mat3x3(MatCrossCross, x, x*dMTmp)).MulMT(R);
}

/* massa totale */
doublereal
VariableBody::dGetM(void) const
{
	return m_Mass.dGet();
}

/* momento statico */
Vec3
VariableBody::GetS(void) const
{
	return GetS_int();
}

/* momento d'inerzia */
Mat3x3
VariableBody::GetJ(void) const
{
	return GetJ_int();
}

/* nodo */
const StructNode *
VariableBody::pGetNode(void) const
{
	return pNode;
}

/* Accesso ai dati privati */
unsigned int
VariableBody::iGetNumPrivData(void) const
{
	return 1;
}

unsigned int
VariableBody::iGetPrivDataIdx(const char *s) const
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
VariableBody::dGetPrivData(unsigned int i) const
{
	if (i == 1) {
		// kinetic energy
		const Mat3x3& Rn = pNode->GetRCurr();
		const Vec3& Vn = pNode->GetVCurr();
		const Vec3& Wn = pNode->GetWCurr();

		Vec3 DXgc = Rn*m_Xgc.Get();
		Vec3 V = Vn + Wn.Cross(DXgc);
		Vec3 W = Rn*Wn;

		return ((V*V)*m_Mass.dGet() + W*(m_Jgc.Get()*W))/2.;
	}

	if (i == 2) {
		// potential energy
		// NOTE: it is only valid for uniform gravity field
		Vec3 Xgc = pNode->GetXCurr() + pNode->GetRCurr()*m_Xgc.Get();

		Vec3 GravityAcceleration;
		if (GravityOwner::bGetGravity(Xgc, GravityAcceleration)) {
			return -Xgc.Dot(GravityAcceleration)*m_Mass.dGet();
		}
	}

	return 0.;
}

void
VariableBody::AssVecRBK_int(SubVectorHandler& WorkVec)
{
	const RigidBodyKinematics *pRBK = pNode->pGetRBK();

	// NOTE: dMTmp, STmp, JTmp updated earlier

	Vec3 s0;

	integer iIdx = 0;
	if (dynamic_cast<DynamicVariableBody *>(this)) {
		iIdx = 6;
	}

	s0 = pNode->GetXCurr()*m_Mass.dGet() + STmp;

	// force
	Vec3 f;
	f = pRBK->GetXPP()*dMTmp;
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
VariableBody::AssMatsRBK_int(
	FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB,
	const doublereal& dCoef,
	const Vec3& Sc)
{
	const RigidBodyKinematics *pRBK = pNode->pGetRBK();

	Mat3x3 MTmp;
	Vec3 VTmp;

	integer iIdx = 0;
	if (dynamic_cast<DynamicVariableBody *>(this)) {
		iIdx = 6;
	}

	// f: delta x
	MTmp = Mat3x3(MatCross, pRBK->GetWP());
	MTmp += Mat3x3(MatCrossCross, pRBK->GetW(), pRBK->GetW());

	WMA.Add(iIdx + 1, 1, MTmp*(dMTmp*dCoef));


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

/* VariableBody - end */


/* DynamicVariableBody - begin */

DynamicVariableBody::DynamicVariableBody(unsigned int uL,
	const DynamicStructNode* pNode,
	const DriveCaller *pDCMass,
	const TplDriveCaller<Vec3> *pDCXgc,
	const TplDriveCaller<Mat3x3> *pDCJgc,
	flag fOut)
: Elem(uL, fOut),
VariableBody(uL, pNode, pDCMass, pDCXgc, pDCJgc, fOut)
{
	NO_OP;
}


/* distruttore */
DynamicVariableBody::~DynamicVariableBody(void)
{
	NO_OP;
}


VariableSubMatrixHandler&
DynamicVariableBody::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("DynamicVariableBody::AssJac");

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
DynamicVariableBody::AssMats(VariableSubMatrixHandler& WorkMatA,
	VariableSubMatrixHandler& WorkMatB,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("DynamicVariableBody::AssMats");

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
DynamicVariableBody::AssMats(FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB,
	doublereal dCoef,
	bool bGravity,
	const Vec3& GravityAcceleration)
{
	DEBUGCOUTFNAME("DynamicVariableBody::AssMats");

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
	WMB.IncCoef(1, 1, dMTmp);
	WMB.IncCoef(2, 2, dMTmp);
	WMB.IncCoef(3, 3, dMTmp);

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
DynamicVariableBody::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUTFNAME("VariableBody::DynamicAssRes");

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

	dMTmp = m_Mass.dGet();
	Vec3 DXgc(R*m_Xgc.Get());
	STmp = DXgc*dMTmp;
	JTmp = R*(m_Jgc.Get() - Mat3x3(MatCrossCross, STmp, DXgc)).MulMT(R);

	/* Quantita' di moto: R[1] = Q - M * V - W /\ S */
	WorkVec.Sub(1, V*dMTmp + W.Cross(STmp));

	/* Momento della quantita' di moto: R[2] = G - S /\ V - J * W */
	WorkVec.Sub(3 + 1, JTmp*W + STmp.Cross(V));

	if (g) {
		WorkVec.Add(6 + 1, GravityAcceleration*dMTmp);
		/* FIXME: this should go into Jacobian matrix
		 * as Gravity /\ S /\ Delta g */
		WorkVec.Add(9 + 1, STmp.Cross(GravityAcceleration));
	}

	if (pRBK) {
		AssVecRBK_int(WorkVec);
	}

	const DynamicStructNode *pDN = dynamic_cast<const DynamicStructNode *>(pNode);
	ASSERT(pDN != 0);

	pDN->AddInertia(dMTmp, STmp, JTmp);

	return WorkVec;
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
DynamicVariableBody::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUTFNAME("DynamicVariableBody::InitialAssJac");

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
	integer iFirstVelocityIndex = iFirstPositionIndex + 6;
	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iFirstVelocityIndex + iCnt);
		WM.PutColIndex(iCnt, iFirstPositionIndex + iCnt);
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
DynamicVariableBody::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUTFNAME("DynamicVariableBody::InitialAssRes");

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

	dMTmp = m_Mass.dGet();
	Vec3 DXgc(R*m_Xgc.Get());
	STmp = DXgc*dMTmp;
	JTmp = R*(m_Jgc.Get() - Mat3x3(MatCrossCross, DXgc, STmp)).MulMT(R);

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
	if (GravityOwner::bGetGravity(X + DXgc, GravityAcceleration)) {
		WorkVec.Add(1, GravityAcceleration*dMTmp);
		WorkVec.Add(3 + 1, STmp.Cross(GravityAcceleration));
		WorkVec.Add(9 + 1, (W.Cross(STmp)).Cross(GravityAcceleration));
	}

	return WorkVec;
}


/* Usata per inizializzare la quantita' di moto */
void
DynamicVariableBody::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& /* XP */ ,
	SimulationEntity::Hints *ph)
{
	integer iFirstIndex = pNode->iGetFirstMomentumIndex();

	const Vec3& V(pNode->GetVCurr());
	const Vec3& W(pNode->GetWCurr());
	const Mat3x3& R(pNode->GetRCurr());

	// FIXME: HACK!
	dMTmp = m_Mass.dGet();
	Vec3 DXgc(R*m_Xgc.Get());
	STmp = DXgc*dMTmp;
	JTmp = R*(m_Jgc.Get() - Mat3x3(MatCrossCross, DXgc, STmp)).MulMT(R);
	X.Add(iFirstIndex + 1, V*dMTmp + W.Cross(STmp));
	X.Add(iFirstIndex + 4, STmp.Cross(V) + JTmp*W);
}

/* momentum */
Vec3
DynamicVariableBody::GetB_int(void) const
{
	const Vec3& V(pNode->GetVCurr());
	const Vec3& W(pNode->GetWCurr());
	const Mat3x3& R(pNode->GetRCurr());

	return (V + W.Cross(R*m_Xgc.Get()))*m_Mass.dGet();
}


/* momenta moment */
Vec3
DynamicVariableBody::GetG_int(void) const
{
	const Vec3& X(pNode->GetXCurr());
	const Mat3x3& R(pNode->GetRCurr());
	const Vec3& V(pNode->GetVCurr());
	const Vec3& W(pNode->GetWCurr());

	Vec3 DXgc(R*m_Xgc.Get());

	// NOTE: with respect to the origin of the global reference frame!
	return (X + DXgc).Cross((V + W.Cross(R*DXgc))*m_Mass.dGet()) + m_Jgc.Get()*W;
}

/* DynamicVariableBody - end */


/* StaticVariableBody - begin */

StaticVariableBody::StaticVariableBody(unsigned int uL,
	const StaticStructNode* pNode,
	const DriveCaller *pDCMass,
	const TplDriveCaller<Vec3> *pDCXgc,
	const TplDriveCaller<Mat3x3> *pDCJgc,
	flag fOut)
: Elem(uL, fOut),
VariableBody(uL, pNode, pDCMass, pDCXgc, pDCJgc, fOut)
{
	NO_OP;
}


/* distruttore */
StaticVariableBody::~StaticVariableBody(void)
{
	NO_OP;
}


VariableSubMatrixHandler&
StaticVariableBody::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("StaticVariableBody::AssJac");

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
StaticVariableBody::AssMats(VariableSubMatrixHandler& WorkMatA,
	VariableSubMatrixHandler& WorkMatB,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("StaticVariableBody::AssMats");

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
StaticVariableBody::AssMats(FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB,
	doublereal dCoef)
{
	DEBUGCOUTFNAME("StaticVariableBody::AssMats");

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
StaticVariableBody::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUTFNAME("StaticVariableBody::AssRes");

	const Vec3& X(pNode->GetXCurr());
	const Mat3x3& R(pNode->GetRCurr());
	Vec3 DXgc(R*m_Xgc.Get());

	/* Se e' definita l'accelerazione di gravita',
	 * la aggiunge (solo al residuo) */
	Vec3 Acceleration(Zero3);
	bool g = GravityOwner::bGetGravity(X + DXgc, Acceleration);

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

	dMTmp = m_Mass.dGet();
	STmp = DXgc*dMTmp;
	JTmp = R*(m_Jgc.Get() - Mat3x3(MatCrossCross, STmp, DXgc)).MulMT(R);

	if (g) {
		WorkVec.Add(1, Acceleration*dMTmp);
		WorkVec.Add(3 + 1, STmp.Cross(Acceleration));
	}

	if (pRBK) {
		AssVecRBK_int(WorkVec);
	}

	return WorkVec;
}


/* inverse dynamics capable element */
bool
StaticVariableBody::bInverseDynamics(void) const
{
	return true;
}


SubVectorHandler&
StaticVariableBody::AssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ ,
	const VectorHandler& /* XPrimePrimeCurr */ ,
	InverseDynamics::Order iOrder)
{
	DEBUGCOUTFNAME("DynamicVariableBody::AssRes");
	
	ASSERT(iOrder == InverseDynamics::INVERSE_DYNAMICS);

	const Vec3& X(pNode->GetXCurr());
	const Mat3x3& R(pNode->GetRCurr());
	Vec3 DXgc(R*m_Xgc.Get());

	/* Se e' definita l'accelerazione di gravita', la aggiunge */
	Vec3 GravityAcceleration;
	bool g = GravityOwner::bGetGravity(X + DXgc, GravityAcceleration);

	WorkVec.ResizeReset(6);

	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
	}

	dMTmp = m_Mass.dGet();
	STmp = DXgc*dMTmp;
	JTmp = R*(m_Jgc.Get() - Mat3x3(MatCrossCross, STmp, DXgc)).MulMT(R);

	Vec3 Acceleration = pNode->GetXPPCurr()
		+ pNode->GetWPCurr().Cross(DXgc)
		+ pNode->GetWCurr().Cross(pNode->GetWCurr().Cross(DXgc));
	if (g) {
		Acceleration -= GravityAcceleration;
	}

	WorkVec.Sub(1, Acceleration*dMTmp);

	Vec3 M = JTmp*pNode->GetWPCurr()
		+ STmp.Cross(pNode->GetXPPCurr())
		+ STmp.Cross(pNode->GetWCurr().Cross(pNode->GetWCurr().Cross(DXgc)));
	if (g) {
		M -= STmp.Cross(GravityAcceleration);
	}

	WorkVec.Sub(4, M);

	return WorkVec;
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
StaticVariableBody::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUTFNAME("StaticVariableBody::InitialAssJac");

	WorkMat.SetNullMatrix();

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
StaticVariableBody::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUTFNAME("StaticVariableBody::InitialAssRes");

	WorkVec.ResizeReset(0);

	return WorkVec;
}


/* Usata per inizializzare la quantita' di moto */
void
StaticVariableBody::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& /* XP */ ,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

/* StaticVariableBody - end */


/* Legge un corpo rigido */
Elem*
ReadVariableBody(DataManager* pDM, MBDynParser& HP, unsigned int uLabel,
	const StructNode* pStrNode)
{
	DEBUGCOUTFNAME("ReadVariableBody");

	/* may be determined by a special DataManager parameter... */
	bool bStaticModel = pDM->bIsStaticModel();
	bool bInverseDynamics = pDM->bIsInverseDynamics();

	ReferenceFrame RF(pStrNode);

	const DriveCaller *pDCMass = HP.GetDriveCaller();
	const TplDriveCaller<Vec3> *pDCXgc = HP.GetTplDriveCaller<Vec3>();
	const TplDriveCaller<Mat3x3> *pDCJgc = HP.GetTplDriveCaller<Mat3x3>();

	const DynamicStructNode* pDynamicNode = 0;
	const StaticStructNode* pStaticNode = 0;

	if (bStaticModel || bInverseDynamics) {
		/* static */
		pStaticNode = dynamic_cast<const StaticStructNode *>(pStrNode);
		if (pStaticNode == 0 || pStaticNode->GetStructNodeType() != StructNode::STATIC) {
			silent_cerr("VariableBody(" << uLabel << "): "
				"illegal structural node type "
				"for StructNode(" << pStrNode->GetLabel() << ") "
				"in " << (bStaticModel ? "static model" : "inverse dynamics") << " analysis "
				"at line " << HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else {
		pDynamicNode = dynamic_cast<const DynamicStructNode*>(pStrNode);
		if (pDynamicNode == 0 || pDynamicNode->GetStructNodeType() != StructNode::DYNAMIC) {
			silent_cerr("VariableBody(" << uLabel << "): "
				"illegal structural node type "
				"for StructNode(" << pStrNode->GetLabel() << ") "
				"at line " << HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	flag fOut = pDM->fReadOutput(HP, Elem::BODY);

	/* Allocazione e costruzione */
	Elem *pEl = 0;
	if (bStaticModel || bInverseDynamics) {
		/* static */
		SAFENEWWITHCONSTRUCTOR(pEl, StaticVariableBody,
			StaticVariableBody(uLabel, pStaticNode,
				pDCMass, pDCXgc, pDCJgc, fOut));

	} else {
		SAFENEWWITHCONSTRUCTOR(pEl, DynamicVariableBody,
			DynamicVariableBody(uLabel, pDynamicNode,
				pDCMass, pDCXgc, pDCJgc, fOut));
	}

	pDM->GetLogFile()
		<< "variable body: " << uLabel
		<< ' ' << pStrNode->GetLabel()
		<< ' ' << pDCMass->dGet()
		<< ' ' << pDCXgc->Get()
		<< ' ' << pDCJgc->Get()
		<< std::endl;

	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		silent_cerr("VariableBody(" << uLabel << "): semicolon expected "
			"at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pEl;
} /* End of ReadVariableBody() */

