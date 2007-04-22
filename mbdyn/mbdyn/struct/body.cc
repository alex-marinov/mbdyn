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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <body.h>
#include <dataman.h>

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

	return pNode->GetRCurr()*J0*pNode->GetRCurr().Transpose()
		- Mat3x3(x, x*dMass) - Mat3x3(s, x) - Mat3x3(x, s);
}


/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
Body::Restart_int(std::ostream& out) const
{
	out << "  body: " << GetLabel() << ", "
		<< pNode->GetLabel() << ", " << dMass << ", "
		<< "reference, node, ", Xgc.Write(out, ", ") << ", "
		<< "reference, node, ", (J0 + Mat3x3(S0, Xgc)).Write(out, ", ");

	return out;
}


void
Body::AfterPredict(VectorHandler& /* X */ , VectorHandler& /* XP */ )
{
	const Mat3x3& R = pNode->GetRRef();

	S = R*S0;
	J = R*(J0*R.Transpose());
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

	/* Dimensiona e resetta la matrice di lavoro */
	WM.ResizeReset(6, 6);

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
		WM.PutRowIndex(iCnt, iFirstPositionIndex+iCnt);
		WM.PutColIndex(iCnt, iFirstPositionIndex+iCnt);
	}

	AssMats(WM, WM, dCoef);

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
	 * e' dato da iGetFirstPositionIndex()+i */
	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WMA.PutRowIndex(iCnt, iFirstPositionIndex+iCnt);
		WMA.PutColIndex(iCnt, iFirstPositionIndex+iCnt);

		WMB.PutRowIndex(iCnt, iFirstPositionIndex+iCnt);
		WMB.PutColIndex(iCnt, iFirstPositionIndex+iCnt);
	}

	AssMats(WMA, WMB, 1.);
}


void
DynamicBody::AssMats(FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB,
	doublereal dCoef)
{
	DEBUGCOUTFNAME("DynamicBody::AssMats");

	const Vec3& V(pNode->GetVCurr());
	const Vec3& W(pNode->GetWRef());

#if 0
	/* ci pensa AfterPredict() */
	S = pNode->GetRRef()*S0;
	J = pNode->GetRRef()*(J0*(pNode->GetRRef()).Transpose());
#endif

	Mat3x3 SWedge(S);			/* S /\ */
	Vec3 Sc(S*dCoef);

	/*
	 * momentum:
	 *
	 * m * I DeltaV - S /\ DeltagP + ( S /\ W ) /\ Deltag
	 */
	WMB.IncCoef(1, 1, dMass);
	WMB.IncCoef(2, 2, dMass);
	WMB.IncCoef(3, 3, dMass);

	WMB.Sub(1, 4, SWedge);
	WMA.Add(1, 4, Mat3x3(Sc.Cross(W)));

	/*
	 * momenta moment:
	 *
	 * S /\ DeltaV + J DeltagP + ( V /\ S /\ - ( J * W ) /\ ) Deltag
	 */
	WMB.Add(4, 1, SWedge);

	WMB.Add(4, 4, J);
	WMA.Add(4, 4, Mat3x3(V, Sc)-Mat3x3(J*(W*dCoef)));
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
	Vec3 GravityAcceleration(0.);
	bool g = GravityOwner::bGetGravity(pNode->GetXCurr(),
		GravityAcceleration);

	integer iNumRows = 6;
	if (g) {
		iNumRows = 12;
	}

	WorkVec.ResizeReset(iNumRows);

	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	for (integer iCnt = 1; iCnt <= iNumRows; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstPositionIndex+iCnt);
	}

	const Vec3& V(pNode->GetVCurr());
	const Vec3& W(pNode->GetWCurr());

	/* Aggiorna i suoi dati (saranno pronti anche per AssJac) */
	const Mat3x3& R(pNode->GetRCurr());
	Vec3 STmp = R*S0;
	Mat3x3 JTmp = R*(J0*R.Transpose());

	/* Quantita' di moto: R[1] = Q - M * V - W /\ S */
	WorkVec.Sub(1, V*dMass + W.Cross(STmp));

	/* Momento della quantita' di moto: R[2] = G - S /\ V - J * W */
	WorkVec.Sub(4, JTmp*W + STmp.Cross(V));

	if (g) {
		WorkVec.Add(7, GravityAcceleration*dMass);
		/* FIXME: this should go into Jacobian matrix
		 * as Gravity /\ S /\ Delta g */
		WorkVec.Add(10, STmp.Cross(GravityAcceleration));
	}

	dynamic_cast<const DynamicStructNode *>(pNode)->AddInertia(dMass, STmp, JTmp);

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

	Vec3 SWedgeW(S.Cross(W));
	Mat3x3 WWedgeSWedge(-W, S);
	Mat3x3 WWedge(W);
	Mat3x3 WWedgeWWedgeSWedge(WWedge*WWedgeSWedge);
	Mat3x3 FDeltaW(Mat3x3(SWedgeW)+WWedgeSWedge);

	Vec3 JW(J*W);
	Mat3x3 JWWedge(JW);
	Mat3x3 MDeltag(WWedge*(J*WWedge-JWWedge));
	Mat3x3 MDeltaW(WWedge*J-JWWedge);

	/* Forza */
	WM.Add(1, 1, WWedgeWWedgeSWedge);
	WM.Add(1, 4, FDeltaW);

	/* Momento */
	WM.Add(4, 1, MDeltag);
	WM.Add(4, 4, MDeltaW);

	/* Derivata forza */
	WM.Add(7, 1, Mat3x3(W.Cross(SWedgeW))+WWedge*FDeltaW);
	WM.Add(7, 4, WWedge*WWedgeWWedgeSWedge);

	/* Derivata Momento */
	WM.Add(4, 1, WWedge*MDeltag);
	WM.Add(4, 4, WWedge*MDeltaW-Mat3x3(W.Cross(JW)));

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

	/* Aggiorna i suoi dati (saranno pronti anche per AssJac) */
	const Mat3x3& R(pNode->GetRCurr());
	Vec3 STmp = R*S0;
	Mat3x3 JTmp = R*J0*R.Transpose();

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
	Vec3 GravityAcceleration(0.);
	if (GravityOwner::bGetGravity(X, GravityAcceleration)) {
		WorkVec.Add(1, GravityAcceleration*dMass);
		WorkVec.Add(4, STmp.Cross(GravityAcceleration));
		WorkVec.Add(10, (W.Cross(STmp)).Cross(GravityAcceleration));
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
	S = R*S0;
	J = R*(J0*R.Transpose());
	X.Add(iFirstIndex + 1, V*dMass + W.Cross(S));
	X.Add(iFirstIndex + 4, S.Cross(V) + J*W);
}

/* DynamicBody - end */


/* StaticBody - begin */

StaticBody::StaticBody(unsigned int uL,
	const StaticStructNode* pNode,
	const StructNode* pRefNode,
	doublereal dMass,
	const Vec3& Xgc,
	const Mat3x3& J,
	flag fOut)
: Elem(uL, fOut),
Body(uL, pNode, dMass, Xgc, J, fOut),
pRefNode(pRefNode)
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

	if (pRefNode != 0) {
		out << ", reference node, "
			<< pRefNode->GetLabel();
	}
		
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

	AssMats(WM, WM, dCoef);

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

	AssMats(WMA, WMB, 1.);
}


void
StaticBody::AssMats(FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB,
	doublereal dCoef)
{
	DEBUGCOUTFNAME("StaticBody::AssMats");

	/* Se e' definita l'accelerazione di gravita',
	 * la aggiunge (solo al residuo) */
	Vec3 Acceleration(0.);
	bool g = GravityOwner::bGetGravity(pNode->GetXCurr(),
		Acceleration);

	/* TODO: reference */
	Vec3 W(0.);
	bool w(pRefNode != 0);

	if (!g && !w) {
		/* FIXME: should be (0, 0), but this requires
		 * setting WMA.SetNullMatrix().
		 * Need to change API */
		WMA.ResizeReset(1, 1);
		WMB.ResizeReset(1, 1);
		return;
	}

	if (w) {
		W = pRefNode->GetWCurr();
		Acceleration -= W.Cross(W.Cross(pNode->GetXCurr() - pRefNode->GetXCurr()));
	}

	if (g || w) {
		WMA.Add(4, 4, Mat3x3(Acceleration, S*dCoef));
	}

	if (w) {
		Mat3x3 wcwc(W, W);
		Mat3x3 Sc(S*dCoef);

		WMA.Sub(1, 1, wcwc*(dMass*dCoef));
		WMA.Sub(1, 4, wcwc*Sc);

		WMA.Sub(4, 1, Sc*wcwc);
		WMA.Sub(4, 4, Mat3x3(W*dCoef, J*W) + Mat3x3(W)*J*Mat3x3(W*dCoef));
	}
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
	Vec3 Acceleration(0.);
	bool g = GravityOwner::bGetGravity(pNode->GetXCurr(), Acceleration);

	/* W is uninitialized because its use is conditioned by w */
	Vec3 W;
	bool w(pRefNode != 0);

	if (!g && !w) {
		WorkVec.ResizeReset(0);
		return WorkVec;
	}

	if (w) {
		W = pRefNode->GetWCurr();
		Acceleration -= W.Cross(W.Cross(pNode->GetXCurr() - pRefNode->GetXCurr()));
	}

	WorkVec.ResizeReset(6);

	integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex + iCnt);
	}

	/* Aggiorna i suoi dati (saranno pronti anche per AssJac) */
	const Mat3x3& R(pNode->GetRCurr());
	Vec3 STmp = R*S0;
	Mat3x3 JTmp = R*(J0*R.Transpose());

	if (g || w) {
		WorkVec.Add(1, Acceleration*dMass);
		WorkVec.Add(4, STmp.Cross(Acceleration));
	}

	if (w) {
		WorkVec.Add(1, W.Cross(W.Cross(STmp)));
		WorkVec.Add(4, W.Cross(JTmp*W));
	}

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
Elem* ReadBody(DataManager* pDM, MBDynParser& HP, unsigned int uLabel)
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
	bool bStatic = pDM->bIsStaticModel();

	integer iNumMasses = 1;
	if (HP.IsKeyWord("condense")) {
		iNumMasses = HP.GetInt();
		if (iNumMasses < 1) {
			silent_cerr("Body(" << uLabel << "): "
				"at least one mass is required in \"condense\" "
				"at line " << HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric();
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
	Vec3 Xgc(0.);
	Vec3 STmp(0.);
	Mat3x3 J(0.);

	for (int iCnt = 1; iCnt <= iNumMasses; iCnt++) {
		/* massa */
		doublereal dMTmp = HP.GetReal();

		DEBUGLCOUT(MYDEBUG_INPUT, "Mass(" << iCnt << ") = " << dMTmp << std::endl);
		dm += dMTmp;

		/* posiz. c.g. */
		Vec3 XgcTmp(HP.GetPosRel(RF));
		STmp += XgcTmp*dMTmp;

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
				JTmp = RTmp*(JTmp*RTmp.Transpose());
			}
			DEBUGLCOUT(MYDEBUG_INPUT, "Inertia matrix of mass(" << iCnt
				<< ") in current frame =" << std::endl << JTmp << std::endl);
		}

		J += (JTmp-Mat3x3(XgcTmp, XgcTmp*dMTmp));
	}

	Xgc = STmp/dm;

	DEBUGLCOUT(MYDEBUG_INPUT, "Total mass: " << dm << std::endl
		<< "Center of mass: " << Xgc << std::endl
		<< "Inertia matrix:" << std::endl << J << std::endl);

	DynamicStructNode* pDynamicNode = 0;
	StaticStructNode* pStaticNode = 0;
	StructNode *pRefNode = 0;

	if (HP.IsKeyWord("reference" "node")) {
		Node *pRN = pDM->ReadNode(HP, Node::STRUCTURAL);
		pRefNode = dynamic_cast<StructNode *>(pRN);
		ASSERT(pRefNode != 0);
		bStatic = true;
	}

	if (bStatic) {
		/* static */
		pStaticNode = dynamic_cast<StaticStructNode*>(pStrNode);
		if (pStaticNode == 0 || pStaticNode->GetStructNodeType() != StructNode::STATIC) {
			silent_cerr("Body(" << uLabel << "): "
				"illegal structural node type "
				"for StructNode(" << pStrNode->GetLabel() << ") "
				"at line " << HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric();
		}

	} else {
		pDynamicNode = dynamic_cast<DynamicStructNode*>(pStrNode);
		if (pDynamicNode == 0 || pDynamicNode->GetStructNodeType() != StructNode::DYNAMIC) {
			silent_cerr("Body(" << uLabel << "): "
				"illegal structural node type "
				"for StructNode(" << pStrNode->GetLabel() << ") "
				"at line " << HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric();
		}
	}

	flag fOut = pDM->fReadOutput(HP, Elem::BODY);

	/* Allocazione e costruzione */
	Elem* pEl = NULL;
	if (bStatic) {
		/* static */
		SAFENEWWITHCONSTRUCTOR(pEl, StaticBody, StaticBody(uLabel, pStaticNode, pRefNode, dm, Xgc, J, fOut));

	} else {
		SAFENEWWITHCONSTRUCTOR(pEl, DynamicBody, DynamicBody(uLabel, pDynamicNode, dm, Xgc, J, fOut));
	}

	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		silent_cerr("semicolon expected "
			"at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric();
	}

	return pEl;
} /* End of ReadBody() */

