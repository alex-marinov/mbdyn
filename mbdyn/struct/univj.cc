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

/* Giunti universali */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <limits>

#include "univj.h"
#include "Rot.hh"

/* UniversalHingeJoint - begin */

/* Costruttore non banale */
UniversalHingeJoint::UniversalHingeJoint(unsigned int uL,
	const DofOwner* pDO,
	const StructNode* pN1,
	const StructNode* pN2,
	const Vec3& dTmp1,
	const Vec3& dTmp2,
	const Mat3x3& R1hTmp,
	const Mat3x3& R2hTmp,
	flag fOut)
: Elem(uL, fOut),
Joint(uL, pDO, fOut),
pNode1(pN1), pNode2(pN2),
d1(dTmp1), R1h(R1hTmp), d2(dTmp2), R2h(R2hTmp), F(Zero3), dM(0.)
{
	NO_OP;
}


/* Distruttore banale */
UniversalHingeJoint::~UniversalHingeJoint(void)
{
	NO_OP;
}


/* Contributo al file di restart */
std::ostream&
UniversalHingeJoint::Restart(std::ostream& out) const
{
	Joint::Restart(out) << ", universal hinge, "
		<< pNode1->GetLabel() << ", "
		"reference, node, ", d1.Write(out, ", ") << ", "
		"hinge, reference, node, "
			"1, ", (R1h.GetVec(1)).Write(out, ", ") << ", "
			"2, ", (R1h.GetVec(2)).Write(out, ", ") << ", "
		<< pNode2->GetLabel() << ", "
		"reference, node, ", d2.Write(out, ", ") << ", "
		"hinge, reference, node, "
			"1, ", (R2h.GetVec(1)).Write(out, ", ") << ", "
			"2, ", (R2h.GetVec(2)).Write(out, ", ") << ";"
		<< std::endl;

	return out;
}


/* Assemblaggio jacobiano */
VariableSubMatrixHandler&
UniversalHingeJoint::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering UniversalHingeJoint::AssJac()" << std::endl);

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
	Vec3 d1Tmp(R1*d1);
	Vec3 d2Tmp(R2*d2);

	/* Suppongo che le reazioni F, M siano gia' state aggiornate da AssRes;
	 * ricordo che la forza F e' nel sistema globale, mentre la coppia M
	 * e' nel sistema locale ed il terzo termine, M(3), e' nullo in quanto
	 * diretto come l'asse attorno al quale la rotazione e' consentita */

	/* Setta gli indici delle equazioni */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	for (int iCnt = 1; iCnt <= 4; iCnt++) {
		WM.PutRowIndex(12 + iCnt, iFirstReactionIndex + iCnt);
		WM.PutColIndex(12 + iCnt, iFirstReactionIndex + iCnt);
	}

	/* Contributo della forza alle equazioni di equilibrio dei due nodi */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutCoef(iCnt, 12 + iCnt, 1.);
		WM.PutCoef(6 + iCnt, 12 + iCnt, -1.);
	}

	WM.Add(4, 13, Mat3x3(MatCross, d1Tmp));
	WM.Sub(10, 13, Mat3x3(MatCross, d2Tmp));

	/* Moltiplica la forza ed il momento per il coefficiente del metodo */
	Vec3 FTmp(F*dCoef);

	Vec3 e3a(R1*R1h.GetVec(3));
	Vec3 e2b(R2*R2h.GetVec(2));
	Vec3 MTmp = e2b*(dM*dCoef);

	Mat3x3 MWedgee3aWedge(MatCrossCross, MTmp, e3a);
	Mat3x3 e3aWedgeMWedge(MatCrossCross, e3a, MTmp);

	WM.Add(4, 4, Mat3x3(MatCrossCross, FTmp, d1Tmp) - MWedgee3aWedge);
	WM.Add(4, 10, e3aWedgeMWedge);

	WM.Add(10, 4, MWedgee3aWedge);
	WM.Sub(10, 10, Mat3x3(MatCrossCross, FTmp, d2Tmp) + e3aWedgeMWedge);

	/* Contributo del momento alle equazioni di equilibrio dei nodi */
	Vec3 Tmp(e2b.Cross(e3a));

	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		doublereal d = Tmp(iCnt);
		WM.PutCoef(3 + iCnt, 16, d);
		WM.PutCoef(9 + iCnt, 16, -d);
	}

	/* Modifica: divido le equazioni di vincolo per dCoef */

	/* Equazioni di vincolo degli spostamenti */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutCoef(12 + iCnt, iCnt, 1.);
		WM.PutCoef(12 + iCnt, 6 + iCnt, -1.);
	}

	WM.Sub(13, 4, Mat3x3(MatCross, d1Tmp));
	WM.Add(13, 10, Mat3x3(MatCross, d2Tmp));

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
		doublereal d = Tmp(iCnt);
		WM.PutCoef(16, 3 + iCnt, d);
		WM.PutCoef(16, 9 + iCnt, -d);
	}

	return WorkMat;
}


/* Assemblaggio residuo */
SubVectorHandler&
UniversalHingeJoint::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering UniversalHingeJoint::AssRes()" << std::endl);

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
		WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
	}

	/* Indici del vincolo */
	for (int iCnt = 1; iCnt <= 4; iCnt++) {
		WorkVec.PutRowIndex(12 + iCnt, iFirstReactionIndex + iCnt);
	}

	/* Aggiorna i dati propri */
	F = Vec3(XCurr, iFirstReactionIndex + 1);
	dM = XCurr(iFirstReactionIndex + 4);

	/* Recupera i dati */
	const Vec3& x1(pNode1->GetXCurr());
	const Vec3& x2(pNode2->GetXCurr());
	const Mat3x3& R1(pNode1->GetRCurr());
	const Mat3x3& R2(pNode2->GetRCurr());

	/* Costruisce i dati propri nella configurazione corrente */
	Vec3 dTmp1(R1*d1);
	Vec3 dTmp2(R2*d2);

	Vec3 e3a(R1*R1h.GetVec(3));
	Vec3 e2b(R2*R2h.GetVec(2));

	Vec3 MTmp(e2b.Cross(e3a)*dM);

	/* Equazioni di equilibrio, nodo 1 */
	WorkVec.Sub(1, F);
	WorkVec.Sub(4, dTmp1.Cross(F) + MTmp);

	/* Equazioni di equilibrio, nodo 2 */
	WorkVec.Add(7, F);
	WorkVec.Add(10, dTmp2.Cross(F) + MTmp);

	/* Modifica: divido le equazioni di vincolo per dCoef */
	ASSERT(dCoef != 0.);

	/* Equazione di vincolo di posizione */
	WorkVec.Add(13, (x2 + dTmp2 - x1 - dTmp1)/dCoef);

	/* Equazione di vincolo di rotazione */
	WorkVec.PutCoef(16, (e3a*e2b)/dCoef);

	return WorkVec;
}

/* Output (da mettere a punto) */
void
UniversalHingeJoint::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {
		Mat3x3 R1Tmp(pNode1->GetRCurr()*R1h);
		Mat3x3 R2Tmp(pNode2->GetRCurr()*R2h);

		Vec3 vTmp(R2Tmp.GetVec(2).Cross(R1Tmp.GetVec(3)));

		Joint::Output(OH.Joints(), "CardanoHinge", GetLabel(),
			R1Tmp.Transpose()*F, Vec3(dM, 0., 0.), F, vTmp*dM)
			<< " " << MatR2EulerAngles(R2Tmp.MulTM(R1Tmp))*dRaDegr
			<< std::endl;
	}
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
UniversalHingeJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering UniversalHingeJoint::InitialAssJac()" << std::endl);

	/* Per ora usa la matrice piena; eventualmente si puo'
	 * passare a quella sparsa quando si ottimizza */
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Equazioni: vedi joints.dvi */

	/* Indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstVelIndex = iNode1FirstPosIndex + 6;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstVelIndex = iNode2FirstPosIndex + 6;
	integer iFirstReactionIndex = iGetFirstIndex();
	integer iReactionPrimeIndex = iFirstReactionIndex + 4;

	/* Nota: le reazioni vincolari sono:
	 * Forza,       3 incognite, riferimento globale,
	 * Momento,     1 incognita, riferimento locale
	 */

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
	for (int iCnt = 1; iCnt <= 8; iCnt++) {
		WM.PutRowIndex(24 + iCnt, iFirstReactionIndex + iCnt);
		WM.PutColIndex(24 + iCnt, iFirstReactionIndex + iCnt);
	}

	/* Matrici identita' */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		/* Contributo di forza all'equazione della forza, nodo 1 */
		WM.IncCoef(iCnt, 24 + iCnt, 1.);

		/* Contrib. di der. di forza all'eq. della der. della forza, nodo 1 */
		WM.IncCoef(6 + iCnt, 28 + iCnt, 1.);

		/* Contributo di forza all'equazione della forza, nodo 2 */
		WM.DecCoef(12 + iCnt, 24 + iCnt, 1.);

		/* Contrib. di der. di forza all'eq. della der. della forza, nodo 2 */
		WM.DecCoef(18 + iCnt, 28 + iCnt, 1.);

		/* Equazione di vincolo, nodo 1 */
		WM.DecCoef(24 + iCnt, iCnt, 1.);

		/* Derivata dell'equazione di vincolo, nodo 1 */
		WM.DecCoef(28 + iCnt, 6 + iCnt, 1.);

		/* Equazione di vincolo, nodo 2 */
		WM.IncCoef(24 + iCnt, 12 + iCnt, 1.);

		/* Derivata dell'equazione di vincolo, nodo 2 */
		WM.IncCoef(28 + iCnt, 18 + iCnt, 1.);
	}

	/* Recupera i dati */
	const Mat3x3& R1(pNode1->GetRRef());
	const Mat3x3& R2(pNode2->GetRRef());
	const Vec3& Omega1(pNode1->GetWRef());
	const Vec3& Omega2(pNode2->GetWRef());
	/* F ed M sono gia' state aggiornate da InitialAssRes */
	Vec3 FPrime(XCurr, iReactionPrimeIndex+1);
	doublereal dMPrime(XCurr(iReactionPrimeIndex+4));

	/* Distanze e matrici di rotazione dai nodi alla cerniera
	 * nel sistema globale */
	Vec3 d1Tmp(R1*d1);
	Vec3 d2Tmp(R2*d2);

	Vec3 e3a(R1*R1h.GetVec(3));
	Vec3 e2b(R2*R2h.GetVec(2));

	/* Ruota il momento e la sua derivata con le matrici della cerniera
	 * rispetto ai nodi */
	Vec3 MTmp(e2b*dM);
	Vec3 MPrimeTmp(e2b*dMPrime);

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
	Vec3 Tmp(e2b.Cross(e3a));

	/* Prodotto vettore tra il versore 3 della cerniera secondo il nodo 1
	 * ed il versore 1 della cerniera secondo il nodo 2. A convergenza
	 * devono essere ortogonali, quindi il loro prodotto vettore deve essere
	 * unitario */

	/* Error handling: il programma si ferma, pero' segnala dov'e' l'errore */
	if (Tmp.Dot() < std::numeric_limits<doublereal>::epsilon()) {
		silent_cerr("CardanoHingeJoint(" << GetLabel() << "): "
			"first and second node hinge axes are (nearly) orthogonal"
			<< std::endl);
		throw ErrNullNorm(MBDYN_EXCEPT_ARGS);
	}

	Vec3 TmpPrime(e2b.Cross(Omega1.Cross(e3a)) - e3a.Cross(Omega2.Cross(e2b)));

	/* Equazione di momento, nodo 1 */
	WM.Add(4, 4, FWedged1Wedge - Mat3x3(MatCrossCross, MTmp, e3a));
	WM.Add(4, 16, Mat3x3(MatCrossCross, e3a, MTmp));
	WM.Add(4, 25, Mat3x3(MatCross, d1Tmp));
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutCoef(3 + iCnt, 28, Tmp(iCnt));
	}

	/* Equazione di momento, nodo 2 */
	WM.Add(16, 4, Mat3x3(MatCrossCross, MTmp, e3a));
	WM.Add(16, 16, FWedged2Wedge - Mat3x3(MatCrossCross, e3a, MTmp));
	WM.Sub(16, 25, Mat3x3(MatCross, d2Tmp));
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.DecCoef(15 + iCnt, 28, Tmp(iCnt));
	}

	/* Derivata dell'equazione di momento, nodo 1 */
	WM.Add(10, 4, (Mat3x3(MatCross, FPrime) + Mat3x3(MatCrossCross, F, Omega1))*Mat3x3(MatCross, d1Tmp) - MDeltag1);
	WM.Add(10, 10, FWedged1Wedge - Mat3x3(MatCrossCross, MTmp, e3a));
	WM.Add(10, 16, MDeltag2);
	WM.Add(10, 22, Mat3x3(MatCrossCross, e3a, MTmp));
	WM.Add(10, 25, O1Wedged1Wedge);
	WM.Add(10, 29, Mat3x3(MatCross, d1Tmp));

	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutCoef(9 + iCnt, 28, TmpPrime(iCnt));
		WM.PutCoef(9 + iCnt, 32, Tmp(iCnt));
	}

	/* Derivata dell'equazione di momento, nodo 2 */
	WM.Add(22, 4, MDeltag1);
	WM.Add(22, 10, Mat3x3(MatCrossCross, MTmp, e3a));
	WM.Sub(22, 16, (Mat3x3(MatCross, FPrime) + Mat3x3(MatCrossCross, F, Omega2))*Mat3x3(MatCross, d2Tmp) + MDeltag2);
	WM.Add(22, 22, FWedged2Wedge - Mat3x3(MatCrossCross, e3a, MTmp));
	WM.Add(22, 25, O2Wedged2Wedge);
	WM.Sub(22, 29, Mat3x3(MatCross, d2Tmp));

	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.DecCoef(21 + iCnt, 28, TmpPrime(iCnt));
		WM.DecCoef(21 + iCnt, 32, Tmp(iCnt));
	}

	/* Equazione di vincolo di posizione */
	WM.Add(25, 4, Mat3x3(MatCross, d1Tmp));
	WM.Sub(25, 16, Mat3x3(MatCross, d2Tmp));

	/* Derivata dell'equazione di vincolo di posizione */
	WM.Add(29, 4, O1Wedged1Wedge);
	WM.Add(29, 10, Mat3x3(MatCross, d1Tmp));
	WM.Add(29, 16, O2Wedged2Wedge);
	WM.Sub(29, 22, Mat3x3(MatCross, d2Tmp));

	/* Equazioni di vincolo di rotazione: e2b~e3a */

	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		doublereal d = Tmp(iCnt);
		WM.IncCoef(28, 3 + iCnt, d);
		WM.DecCoef(28, 15 + iCnt, d);

		/* Queste sono per la derivata dell'equazione, sono qui solo per
		 * ottimizzazione */
		WM.IncCoef(32, 9 + iCnt, d);
		WM.DecCoef(32, 21 + iCnt, d);
	}

	/* Derivate delle equazioni di vincolo di rotazione: e2b~e3a */
	Vec3 O1mO2(Omega1 - Omega2);
	TmpPrime = e3a.Cross(O1mO2.Cross(e2b));
	Vec3 TmpPrime2 = e2b.Cross(e3a.Cross(O1mO2));
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutCoef(32, 3 + iCnt, TmpPrime(iCnt));
		WM.PutCoef(32, 15 + iCnt, TmpPrime2(iCnt));
	}

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
UniversalHingeJoint::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering UniversalHingeJoint::InitialAssRes()" << std::endl);

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
	integer iReactionPrimeIndex = iFirstReactionIndex + 4;

	/* Setta gli indici */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode1FirstVelIndex + iCnt);
		WorkVec.PutRowIndex(12 + iCnt, iNode2FirstPosIndex + iCnt);
		WorkVec.PutRowIndex(18 + iCnt, iNode2FirstVelIndex + iCnt);
	}

	for (int iCnt = 1; iCnt <= 8; iCnt++) {
		WorkVec.PutRowIndex(24 + iCnt, iFirstReactionIndex+iCnt);
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
	F = Vec3(XCurr, iFirstReactionIndex + 1);
	dM = XCurr(iFirstReactionIndex + 4);
	Vec3 FPrime(XCurr, iReactionPrimeIndex + 1);
	doublereal dMPrime(XCurr(iReactionPrimeIndex + 4));

	/* Distanza nel sistema globale */
	Vec3 d1Tmp(R1*d1);
	Vec3 d2Tmp(R2*d2);

	/* Vettori omega1/\d1, -omega2/\d2 */
	Vec3 O1Wedged1(Omega1.Cross(d1Tmp));
	Vec3 O2Wedged2(Omega2.Cross(d2Tmp));

	Vec3 e3a(R1*R1h.GetVec(3));
	Vec3 e2b(R2*R2h.GetVec(2));

	/* Ruota il momento e la sua derivata con le matrici della cerniera
	 * rispetto ai nodi */
	Vec3 MTmp(e2b*dM);
	Vec3 MPrimeTmp(e3a.Cross(MTmp.Cross(Omega2)) + e2b.Cross(e3a)*dMPrime);

	/* Equazioni di equilibrio, nodo 1 */
	WorkVec.Sub(1, F);
	WorkVec.Sub(4, d1Tmp.Cross(F) + MTmp.Cross(e3a));

	/* Derivate delle equazioni di equilibrio, nodo 1 */
	WorkVec.Sub(7, FPrime);
	WorkVec.Sub(10, d1Tmp.Cross(FPrime) - O1Wedged1.Cross(F) - MPrimeTmp);

	/* Equazioni di equilibrio, nodo 2 */
	WorkVec.Add(13, F);
	WorkVec.Add(16, d2Tmp.Cross(F) + MTmp.Cross(e3a));

	/* Derivate delle equazioni di equilibrio, nodo 2 */
	WorkVec.Add(19, FPrime);
	WorkVec.Add(22, d2Tmp.Cross(FPrime) + O2Wedged2.Cross(F) + MPrimeTmp);

	/* Equazione di vincolo di posizione */
	WorkVec.Add(25, x1 + d1Tmp - x2 - d2Tmp);

	/* Derivata dell'equazione di vincolo di posizione */
	WorkVec.Add(29, v1 + O1Wedged1 - v2 - O2Wedged2);

	/* Equazioni di vincolo di rotazione */
	WorkVec.PutCoef(28, e2b*e3a);

	/* Derivate delle equazioni di vincolo di rotazione: e2b~e3a */
	Vec3 Tmp((Omega1 - Omega2).Cross(e3a));
	WorkVec.PutCoef(32, e2b*Tmp);

	return WorkVec;
}

/* UniversalHingeJoint - end */


/* UniversalRotationJoint - begin */

/* Costruttore non banale */
UniversalRotationJoint::UniversalRotationJoint(unsigned int uL,
	const DofOwner* pDO,
	const StructNode* pN1,
	const StructNode* pN2,
	const Mat3x3& R1hTmp,
	const Mat3x3& R2hTmp,
	const OrientationDescription& od,
	flag fOut)
: Elem(uL, fOut),
Joint(uL, pDO, fOut),
pNode1(pN1), pNode2(pN2),
#ifdef USE_NETCDFC // netcdfcxx4 has non-pointer vars...
Var_Phi(0),
#endif // USE_NETCDFC
R1h(R1hTmp), R2h(R2hTmp), dM(0.), od(od)
{
	NO_OP;
}


/* Distruttore banale */
UniversalRotationJoint::~UniversalRotationJoint(void)
{
	NO_OP;
}


/* Contributo al file di restart */
std::ostream&
UniversalRotationJoint::Restart(std::ostream& out) const
{
	Joint::Restart(out) << ", universal rotation, "
		<< pNode1->GetLabel() << ", hinge, reference, node, "
			"1, ", (R1h.GetVec(1)).Write(out, ", ") << ", "
			"2, ", (R1h.GetVec(2)).Write(out, ", ") << ", "
		<< pNode2->GetLabel() << ", hinge, reference, node, "
			"1, ", (R2h.GetVec(1)).Write(out, ", ") << ", "
			"2, ", (R2h.GetVec(2)).Write(out, ", ") << ";"
		<< std::endl;

	return out;
}


/* Assemblaggio jacobiano */
VariableSubMatrixHandler&
UniversalRotationJoint::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering UniversalRotationJoint::AssJac()" << std::endl);

	/* Setta la sottomatrice come piena (e' un po' dispersivo, ma lo jacobiano
	 * e' complicato */
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

	WM.PutRowIndex(6 + 1, iFirstReactionIndex + 1);
	WM.PutColIndex(6 + 1, iFirstReactionIndex + 1);

	/* Recupera i dati che servono */
	const Mat3x3& R1(pNode1->GetRRef());
	const Mat3x3& R2(pNode2->GetRRef());

	/* Suppongo che le reazioni F, M siano gia' state aggiornate da AssRes;
	 * ricordo che la forza F e' nel sistema globale, mentre la coppia M
	 * e' nel sistema locale ed il terzo termine, M(3), e' nullo in quanto
	 * diretto come l'asse attorno al quale la rotazione e' consentita */

	Vec3 e3a(R1*R1h.GetVec(3));
	Vec3 e2b(R2*R2h.GetVec(2));
	Vec3 MTmp = e2b*(dM*dCoef);

	Mat3x3 MWedgee3aWedge(MatCrossCross, MTmp, e3a);
	Mat3x3 e3aWedgeMWedge(MatCrossCross, e3a, MTmp);

	WM.Sub(0 + 1, 0 + 1, MWedgee3aWedge);
	WM.Add(0 + 1, 3 + 1, e3aWedgeMWedge);

	WM.Add(3 + 1, 0 + 1, MWedgee3aWedge);
	WM.Sub(3 + 1, 3 + 1, e3aWedgeMWedge);

	/* Contributo del momento alle equazioni di equilibrio dei nodi */
	Vec3 Tmp(e2b.Cross(e3a));

	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		doublereal d = Tmp(iCnt);
		WM.IncCoef(iCnt, 6 + 1, d);
		WM.DecCoef(3 + iCnt, 6 + 1, d);

		/* Modifica: divido le equazioni di vincolo per dCoef */
		WM.IncCoef(6 + 1, iCnt, d);
		WM.DecCoef(6 + 1, 3 + iCnt, d);
	}

	return WorkMat;
}


/* Assemblaggio residuo */
SubVectorHandler&
UniversalRotationJoint::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering UniversalRotationJoint::AssRes()" << std::endl);

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
		WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WorkVec.PutRowIndex(3 + iCnt, iNode2FirstMomIndex + iCnt);
	}

	/* Indici del vincolo */
	WorkVec.PutRowIndex(6 + 1, iFirstReactionIndex + 1);

	/* Aggiorna i dati propri */
	dM = XCurr(iFirstReactionIndex+1);

	/* Recupera i dati */
	const Mat3x3& R1(pNode1->GetRCurr());
	const Mat3x3& R2(pNode2->GetRCurr());

	/* Costruisce i dati propri nella configurazione corrente */
	Vec3 e3a(R1*R1h.GetVec(3));
	Vec3 e2b(R2*R2h.GetVec(2));

	Vec3 MTmp(e2b.Cross(e3a)*dM);

	/* Equazioni di equilibrio, nodo 1 */
	WorkVec.Sub(1, MTmp); /* Sfrutto  F/\d = -d/\F */

	/* Equazioni di equilibrio, nodo 2 */
	WorkVec.Add(4, MTmp);

	/* Modifica: divido le equazioni di vincolo per dCoef */
	ASSERT(dCoef != 0.);
	/* Equazione di vincolo di rotazione */
	WorkVec.PutCoef(6 + 1, (e3a*e2b)/dCoef);

	return WorkVec;
}

void
UniversalRotationJoint::OutputPrepare(OutputHandler& OH)
{
	if (bToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			std::string name;
			OutputPrepare_int("cardano rotation", OH, name);

			Var_Phi = OH.CreateRotationVar(name, "", od, "global");

		}
#endif // USE_NETCDF
	}
}


/* Output (da mettere a punto) */
void
UniversalRotationJoint::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {
		Mat3x3 R1Tmp(pNode1->GetRCurr()*R1h);
		Mat3x3 R2Tmp(pNode2->GetRCurr()*R2h);
		Vec3 vTmp(R2Tmp.GetVec(2).Cross(R1Tmp.GetVec(3)));
		Mat3x3 RTmp(R2Tmp.Transpose()*R1Tmp);
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
			
			OH.WriteNcVar(Var_F_local, Zero3);
			OH.WriteNcVar(Var_M_local, Vec3(dM, 0., 0.));
			OH.WriteNcVar(Var_F_global, Zero3);
			OH.WriteNcVar(Var_M_global, (vTmp*dM));
			
			switch (od) {
			case EULER_123:
			case EULER_313:
			case EULER_321:
			case ORIENTATION_VECTOR:
				OH.WriteNcVar(Var_Phi, E);
				break;

			case ORIENTATION_MATRIX:
				OH.WriteNcVar(Var_Phi, RTmp);
				break;

			default:
				/* impossible */
				break;
			}
		}
#endif // USE_NETCDF
		if (OH.UseText(OutputHandler::JOINTS)) {
			Joint::Output(OH.Joints(), "CardanoHinge", GetLabel(),
				Zero3, Vec3(dM, 0., 0.), Zero3, vTmp*dM)
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

				OH.Joints() << std::endl;
		}
	}
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
UniversalRotationJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering UniversalRotationJoint::InitialAssJac()" << std::endl);

	/* Per ora usa la matrice piena; eventualmente si puo'
	 * passare a quella sparsa quando si ottimizza */
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Equazioni: vedi joints.dvi */

	/* Indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex() + 3;
	integer iNode1FirstVelIndex = iNode1FirstPosIndex + 6;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex() + 3;
	integer iNode2FirstVelIndex = iNode2FirstPosIndex + 6;
	integer iFirstReactionIndex = iGetFirstIndex();
	integer iReactionPrimeIndex = iFirstReactionIndex + 1;

	/* Nota: le reazioni vincolari sono:
	 * Forza,       3 incognite, riferimento globale,
	 * Momento,     1 incognita, riferimento locale
	 */

	/* Setta gli indici dei nodi */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(3 + iCnt, iNode1FirstVelIndex + iCnt);
		WM.PutColIndex(3 + iCnt, iNode1FirstVelIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutRowIndex(9 + iCnt, iNode2FirstVelIndex + iCnt);
		WM.PutColIndex(9 + iCnt, iNode2FirstVelIndex + iCnt);
	}

	/* Setta gli indici delle reazioni */
	for (int iCnt = 1; iCnt <= 2; iCnt++) {
		WM.PutRowIndex(12 + iCnt, iFirstReactionIndex + iCnt);
		WM.PutColIndex(12 + iCnt, iFirstReactionIndex + iCnt);
	}

	/* Recupera i dati */
	const Mat3x3& R1(pNode1->GetRRef());
	const Mat3x3& R2(pNode2->GetRRef());
	const Vec3& Omega1(pNode1->GetWRef());
	const Vec3& Omega2(pNode2->GetWRef());
	/* F ed M sono gia' state aggiornate da InitialAssRes */
	doublereal dMPrime(XCurr(iReactionPrimeIndex+1));

	/* Matrici di rotazione dai nodi alla cerniera
	 * nel sistema globale */
	Vec3 e3a(R1*R1h.GetVec(3));
	Vec3 e2b(R2*R2h.GetVec(2));

	/* Ruota il momento e la sua derivata con le matrici della cerniera
	 * rispetto ai nodi */
	Vec3 MTmp(e2b*dM);
	Vec3 MPrimeTmp(e2b*dMPrime);

	/* Matrici (omega1/\d1)/\, -(omega2/\d2)/\ */
	Mat3x3 MDeltag1((Mat3x3(MatCross, Omega2.Cross(MTmp) + MPrimeTmp)
		+ Mat3x3(MatCrossCross, MTmp, Omega1))*Mat3x3(MatCross, e3a));
	Mat3x3 MDeltag2(Mat3x3(MatCrossCross, Omega1.Cross(e3a), MTmp)
		+ Mat3x3(MatCrossCross, e3a, MPrimeTmp)
		+ e3a.Cross(Mat3x3(MatCrossCross, Omega2, MTmp)));

	/* Vettori temporanei */
	Vec3 Tmp(e2b.Cross(e3a));

	/* Prodotto vettore tra il versore 3 della cerniera secondo il nodo 1
	 * ed il versore 1 della cerniera secondo il nodo 2. A convergenza
	 * devono essere ortogonali, quindi il loro prodotto vettore deve essere
	 * unitario */

	/* Error handling: il programma si ferma, pero' segnala dov'e' l'errore */
	if (Tmp.Dot() < std::numeric_limits<doublereal>::epsilon()) {
		silent_cerr("CardanoRotationJoint(" << GetLabel() << "): "
			"first and second node hinge axes are (nearly) orthogonal"
			<< std::endl);
		throw ErrNullNorm(MBDYN_EXCEPT_ARGS);
	}

	Vec3 TmpPrime(e2b.Cross(Omega1.Cross(e3a)) - e3a.Cross(Omega2.Cross(e2b)));

	/* Equazione di momento, nodo 1 */
	WM.Sub(1, 1, Mat3x3(MatCrossCross, MTmp, e3a));
	WM.Add(4, 7, Mat3x3(MatCrossCross, e3a, MTmp));
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutCoef(iCnt, 13, Tmp.dGet(iCnt));
	}

	/* Equazione di momento, nodo 2 */
	WM.Add(7, 1, Mat3x3(MatCrossCross, MTmp, e3a));
	WM.Sub(7, 7, Mat3x3(MatCrossCross, e3a, MTmp));
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.DecCoef(6 + iCnt, 13, Tmp(iCnt));
	}

	/* Derivata dell'equazione di momento, nodo 1 */
	WM.Sub(4, 1, MDeltag1);
	WM.Sub(4, 4, Mat3x3(MatCrossCross, MTmp, e3a));
	WM.Add(4, 7, MDeltag2);
	WM.Add(4, 10, Mat3x3(MatCrossCross, e3a, MTmp));

	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutCoef(9 + iCnt, 13, TmpPrime(iCnt));
		WM.PutCoef(9 + iCnt, 14, Tmp(iCnt));
	}

	/* Derivata dell'equazione di momento, nodo 2 */
	WM.Add(4, 1, Mat3x3(MatCrossCross, MTmp, e3a));
	WM.Sub(4, 4, Mat3x3(MatCrossCross, e3a, MTmp));

	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.DecCoef(3 + iCnt, 13, TmpPrime(iCnt));
		WM.DecCoef(3 + iCnt, 14, Tmp(iCnt));
	}

	/* Equazioni di vincolo di rotazione: e2b~e3a */

	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		doublereal d = Tmp(iCnt);
		WM.IncCoef(13, iCnt, d);
		WM.DecCoef(13, 6 + iCnt, d);

		/* Queste sono per la derivata dell'equazione, sono qui solo per
		 * ottimizzazione */
		WM.IncCoef(14, 3 + iCnt, d);
		WM.DecCoef(14, 9 + iCnt, d);
	}

	/* Derivate delle equazioni di vincolo di rotazione: e2b~e3a */
	Vec3 O1mO2(Omega1 - Omega2);
	TmpPrime = e3a.Cross(O1mO2.Cross(e2b));
	Vec3 TmpPrime2 = e2b.Cross(e3a.Cross(O1mO2));
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutCoef(14, iCnt, TmpPrime(iCnt));
		WM.PutCoef(14, 6 + iCnt, TmpPrime2(iCnt));
	}

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
UniversalRotationJoint::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering UniversalRotationJoint::InitialAssRes()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex() + 3;
	integer iNode1FirstVelIndex = iNode1FirstPosIndex + 6;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex() + 3;
	integer iNode2FirstVelIndex = iNode2FirstPosIndex + 6;
	integer iFirstReactionIndex = iGetFirstIndex();
	integer iReactionPrimeIndex = iFirstReactionIndex + 1;

	/* Setta gli indici */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WorkVec.PutRowIndex(3 + iCnt, iNode1FirstVelIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
		WorkVec.PutRowIndex(9 + iCnt, iNode2FirstVelIndex + iCnt);
	}

	for (int iCnt = 1; iCnt <= 2; iCnt++) {
		WorkVec.PutRowIndex(12 + iCnt, iFirstReactionIndex + iCnt);
	}

	/* Recupera i dati */
	const Mat3x3& R1(pNode1->GetRCurr());
	const Mat3x3& R2(pNode2->GetRCurr());
	const Vec3& Omega1(pNode1->GetWCurr());
	const Vec3& Omega2(pNode2->GetWCurr());

	/* Aggiorna F ed M, che restano anche per InitialAssJac */
	dM = XCurr(iFirstReactionIndex + 1);
	doublereal dMPrime(XCurr(iReactionPrimeIndex + 1));

	/* orientazione nel sistema globale */
	Vec3 e3a(R1*R1h.GetVec(3));
	Vec3 e2b(R2*R2h.GetVec(2));

	/* Ruota il momento e la sua derivata con le matrici della cerniera
	 * rispetto ai nodi */
	Vec3 MTmp(e2b*dM);
	Vec3 MPrimeTmp(e3a.Cross(MTmp.Cross(Omega2)) + e2b.Cross(e3a)*dMPrime);

	/* Equazioni di equilibrio, nodo 1 */
	WorkVec.Sub(1, MTmp.Cross(e3a));

	/* Derivate delle equazioni di equilibrio, nodo 1 */
	WorkVec.Sub(4, MPrimeTmp);

	/* Equazioni di equilibrio, nodo 2 */
	WorkVec.Add(7, MTmp.Cross(e3a));

	/* Derivate delle equazioni di equilibrio, nodo 2 */
	WorkVec.Add(10, MPrimeTmp);

	/* Equazioni di vincolo di rotazione */
	WorkVec.PutCoef(13, e2b*e3a);

	/* Derivate delle equazioni di vincolo di rotazione: e2b~e3a */
	Vec3 Tmp((Omega1 - Omega2).Cross(e3a));
	WorkVec.PutCoef(14, e2b*Tmp);

	return WorkVec;
}

/* UniversalRotationJoint - end */


/* UniversalPinJoint - begin */

/* Costruttore non banale */
UniversalPinJoint::UniversalPinJoint(unsigned int uL, const DofOwner* pDO,
	const StructNode* pN,
	const Vec3& X0Tmp, const Mat3x3& R0Tmp,
	const Vec3& dTmp, const Mat3x3& RhTmp,
	flag fOut)
: Elem(uL, fOut),
Joint(uL, pDO, fOut),
pNode(pN),
X0(X0Tmp), R0(R0Tmp), d(dTmp), Rh(RhTmp), F(Zero3), dM(0.)
{
	NO_OP;
}


/* Distruttore banale */
UniversalPinJoint::~UniversalPinJoint(void)
{
	NO_OP;
}


/* Contributo al file di restart */
std::ostream&
UniversalPinJoint::Restart(std::ostream& out) const
{
	Joint::Restart(out) << ", universal pin, "
		<< pNode->GetLabel() << ", "
		"reference, node, ", d.Write(out, ", ") << ", "
		"hinge, reference, node, "
			"1, ", (Rh.GetVec(1)).Write(out, ", ") << ", "
			"2, ", (Rh.GetVec(2)).Write(out, ", ") << ", "
		"reference, global, ", X0.Write(out, ", ") << ", "
		"reference, global, "
			"1, ", (R0.GetVec(1)).Write(out, ", ") << ", "
			"2, ", (R0.GetVec(2)).Write(out, ", ") << ';'
		<< std::endl;

	return out;
}


/* Assemblaggio jacobiano */
VariableSubMatrixHandler&
UniversalPinJoint::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering UniversalPinJoint::AssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
	integer iFirstReactionIndex = iGetFirstIndex();

	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iFirstMomentumIndex + iCnt);
		WM.PutColIndex(iCnt, iFirstPositionIndex + iCnt);
	}

	for (int iCnt = 1; iCnt <= 4; iCnt++) {
		WM.PutRowIndex(6 + iCnt, iFirstReactionIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iFirstReactionIndex + iCnt);
	}

	const Mat3x3& R(pNode->GetRRef());
	Vec3 dTmp(R*d);

	/* termini di reazione sul nodo (forza e momento) */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutCoef(iCnt, 6 + iCnt, -1.);
	}

	WM.Sub(4, 7, Mat3x3(MatCross, dTmp));

	/* Note: F and dM updated by AssRes */
   	Vec3 e3(R0.GetVec(3));
	Vec3 e2(R*Rh.GetVec(2));
	Vec3 MTmp(e2*(dM*dCoef));
	Vec3 FTmp(F*dCoef);

	Mat3x3 e3aWedgeMWedge(MatCrossCross, e3, MTmp);

	WM.Sub(4, 4, Mat3x3(MatCrossCross, FTmp, dTmp) + e3aWedgeMWedge);

	Vec3 Tmp(e2.Cross(e3));

	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		doublereal d = Tmp(iCnt);
		WM.PutCoef(3 + iCnt, 10, -d);
	}

	/* Modifica: divido le equazioni di vincolo per dCoef */

	/* termini di vincolo dovuti al nodo 1 */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		// position, delta x
		WM.PutCoef(6 + iCnt, iCnt, -1.);

		// orientation, delta g
		doublereal d = Tmp(iCnt);
		WM.PutCoef(6 + 3 + 1, 3 + iCnt, -d);
	}

	// position, delta g
	WM.Add(6 + 1, 3 + 1, Mat3x3(MatCross, dTmp));

	return WorkMat;
}


/* Assemblaggio residuo */
SubVectorHandler& UniversalPinJoint::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering UniversalPinJoint::AssRes()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
	integer iFirstReactionIndex = iGetFirstIndex();

	/* Indici dei nodi */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex + iCnt);
	}

	/* Indici del vincolo */
	for (int iCnt = 1; iCnt <= 4; iCnt++) {
		WorkVec.PutRowIndex(6 + iCnt, iFirstReactionIndex + iCnt);
	}

	F = Vec3(XCurr, iFirstReactionIndex + 1);
	dM = XCurr(iFirstReactionIndex + 4);

	const Vec3& x(pNode->GetXCurr());
	const Mat3x3& R(pNode->GetRCurr());

	Vec3 dTmp(R*d);

	Vec3 e3(R0.GetVec(3));
	Vec3 e2(R*Rh.GetVec(2));

	WorkVec.Add(1, F);
	WorkVec.Add(4, dTmp.Cross(F) + e2.Cross(e3)*dM);

	/* Modifica: divido le equazioni di vincolo per dCoef */
	ASSERT(dCoef != 0.);
	WorkVec.Add(7, (x + dTmp - X0)/dCoef);

	WorkVec.PutCoef(10, e3.Dot(e2)/dCoef);

	return WorkVec;
}

/* Output (da mettere a punto) */
void
UniversalPinJoint::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {
		Mat3x3 RTmp(pNode->GetRCurr()*Rh);
		Vec3 vTmp(RTmp.GetVec(2).Cross(R0.GetVec(3)));

		Joint::Output(OH.Joints(), "CardanoPin", GetLabel(),
			RTmp.MulTV(F), Vec3(dM, 0., 0.), F, vTmp*dM)
			<< " " << MatR2EulerAngles(R0.Transpose()*RTmp)*dRaDegr
			<< std::endl;
	}
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
UniversalPinJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering UniversalPinJoint::InitialAssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Equazioni: vedi joints.dvi */

	/* Indici */
	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	integer iFirstVelocityIndex = iFirstPositionIndex + 6;
	integer iFirstReactionIndex = iGetFirstIndex();
	integer iReactionPrimeIndex = iFirstReactionIndex + 4;

	/* Setto gli indici */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
		WM.PutColIndex(iCnt, iFirstPositionIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iFirstVelocityIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iFirstVelocityIndex + iCnt);
	}

	for (int iCnt = 1; iCnt <= 8; iCnt++) {
		WM.PutRowIndex(12 + iCnt, iFirstReactionIndex + iCnt);
		WM.PutColIndex(12 + iCnt, iFirstReactionIndex + iCnt);
	}

	/* Recupera i dati */
	const Mat3x3& R(pNode->GetRRef());
	const Vec3& Omega(pNode->GetWRef());
	/* F, M sono state aggiornate da InitialAssRes */
	Vec3 FPrime(XCurr, iReactionPrimeIndex + 1);
	doublereal dMPrime(XCurr(iReactionPrimeIndex + 4));

	/* Matrici identita' */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		/* Contributo di forza all'equazione della forza */
		WM.PutCoef(iCnt, 12 + iCnt, 1.);

		/* Contrib. di der. di forza all'eq. della der. della forza */
		WM.PutCoef(6 + iCnt, 16 + iCnt, 1.);

		/* Equazione di vincolo */
		WM.PutCoef(12 + iCnt, iCnt, -1.);

		/* Derivata dell'equazione di vincolo */
		WM.PutCoef(16 + iCnt, 6 + iCnt, -1.);
	}

	/* Distanza nel sistema globale */
	Vec3 dTmp(R*d);

	Vec3 e3(R0.GetVec(3));
	Vec3 e2(R*Rh.GetVec(2));

	/* Vettori temporanei */
	Vec3 Tmp(e2.Cross(e3));

	/* Prodotto vettore tra il versore 3 della cerniera secondo il nodo 1
	 * ed il versore 1 della cerniera secondo il nodo 2. A convergenza
	 * devono essere ortogonali, quindi il loro prodotto vettore deve essere
	 * unitario */

	/* Error handling: il programma si ferma, pero' segnala dov'e' l'errore */
	if (Tmp.Dot() < std::numeric_limits<doublereal>::epsilon()) {
		silent_cerr("CardanoPinJoint(" << GetLabel() << "): "
			"node and fixed point hinge axes are (nearly) orthogonal"
			<< std::endl);
		throw ErrNullNorm(MBDYN_EXCEPT_ARGS);
	}

	Vec3 TmpPrime(e3.Cross(e2.Cross(Omega)));

	/* Ruota il momento e la sua derivata con le matrici della cerniera
	 * rispetto ai nodi */
	Vec3 MTmp(e2*dM);
	Vec3 MPrimeTmp(e2*dMPrime);

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
	WM.Add(10, 17, Mat3x3(MatCross, dTmp));

	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		doublereal d = Tmp(iCnt);
		WM.PutCoef(3 + iCnt, 16, d);
		WM.PutCoef(9 + iCnt, 20, d);

		WM.PutCoef(9 + iCnt, 16, TmpPrime(iCnt));
	}

	/* Equazione di vincolo */
	WM.Add(13, 4, Mat3x3(MatCross, dTmp));

	/* Derivata dell'equazione di vincolo */
	WM.Add(17, 4, OWedgedWedge);
	WM.Add(17, 10, Mat3x3(MatCross, dTmp));

	/* Equazioni di vincolo di rotazione: e2b~e3a */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		doublereal d = -Tmp(iCnt);
		WM.PutCoef(16, 3 + iCnt, d);

		/* Queste sono per la derivata dell'equazione, sono qui solo per
		 * ottimizzazione */
		WM.PutCoef(20, 9 + iCnt, d);
	}

	/* Derivate delle equazioni di vincolo di rotazione: e2b~e3a */
	TmpPrime = e2.Cross(Omega.Cross(e3));
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutCoef(20, 3 + iCnt, TmpPrime(iCnt));
	}

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
UniversalPinJoint::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering UniversalPinJoint::InitialAssRes()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici */
	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	integer iFirstVelocityIndex = iFirstPositionIndex + 6;
	integer iFirstReactionIndex = iGetFirstIndex();
	integer iReactionPrimeIndex = iFirstReactionIndex + 4;

	/* Setta gli indici */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iFirstVelocityIndex + iCnt);
	}

	for (int iCnt = 1; iCnt <= 8; iCnt++) {
		WorkVec.PutRowIndex(12 + iCnt, iFirstReactionIndex + iCnt);
	}

	/* Recupera i dati */
	const Vec3& x(pNode->GetXCurr());
	const Vec3& v(pNode->GetVCurr());
	const Mat3x3& R(pNode->GetRCurr());
	const Vec3& Omega(pNode->GetWCurr());

	F = Vec3(XCurr, iFirstReactionIndex + 1);
	dM = XCurr(iFirstReactionIndex + 4);
	Vec3 FPrime(XCurr, iReactionPrimeIndex + 1);
	doublereal dMPrime(XCurr(iReactionPrimeIndex + 4));

	/* Versori delle cerniere */
	Vec3 e3(R0.GetVec(3));
	Vec3 e2(R*Rh.GetVec(2));

	/* Vettori temporanei */
	Vec3 Tmp(e2.Cross(e3));

	Vec3 TmpPrime(e3.Cross(e2.Cross(Omega)));

	/* Distanza nel sistema globale */
	Vec3 dTmp(R*d);

	/* Vettori omega/\d */
	Vec3 OWedged(Omega.Cross(dTmp));

	/* Ruota il momento e la sua derivata con le matrici della cerniera
	 * rispetto ai nodi */
	Vec3 MTmp(e2*dM);
	Vec3 MPrimeTmp(e3.Cross(MTmp.Cross(Omega)) + e2.Cross(e3)*dMPrime);

	/* Equazioni di equilibrio */
	WorkVec.Sub(1, F);
	WorkVec.Sub(4, dTmp.Cross(F) - MTmp.Cross(e3));

	/* Derivate delle equazioni di equilibrio, nodo 1 */
	WorkVec.Sub(7, FPrime);
	WorkVec.Sub(10, dTmp.Cross(FPrime) + OWedged.Cross(F) + MPrimeTmp);

	/* Equazione di vincolo di posizione */
	WorkVec.Add(13, x + dTmp - X0);

	/* Equazioni di vincolo di rotazione */
	WorkVec.PutCoef(16, e2*e3);

	/* Derivata dell'equazione di vincolo di posizione */
	WorkVec.Add(17, v + OWedged);

	/* Derivate delle equazioni di vincolo di rotazione: e2b~e3a */
	Tmp = e3.Cross(Omega);
	WorkVec.PutCoef(20, e2*Tmp);

	return WorkVec;
}

/* UniversalPinJoint - end */
