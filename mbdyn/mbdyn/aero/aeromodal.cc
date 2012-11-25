/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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
 /* Elemento aerodinamico modale */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cerrno>

#include "aeromodal.h"
#include "dataman.h"

/* AerodynamicModal - begin */

AerodynamicModal::AerodynamicModal(unsigned int uLabel,
	const StructNode* pN,
	const Modal* pMJ,
	const Mat3x3& RaTmp,
	const DofOwner* pDO,
	doublereal Cd,
	const int NModal,
	const int NAero,
	RigidF_t rgF,
	const int Gust,
	const doublereal Vff,
	SpMapMatrixHandler* pAMat,
	FullMatrixHandler* pBMat,
	FullMatrixHandler* pCMat,
	FullMatrixHandler* pD0Mat,
	FullMatrixHandler* pD1Mat,
	FullMatrixHandler* pD2Mat,
	flag fout)
: Elem(uLabel, fout),
AerodynamicElem(uLabel, pDO, fout),
InitialAssemblyElem(uLabel, fout),
pModalNode(pN), pModalJoint(pMJ),
Ra(RaTmp),
Chord(Cd),
NStModes(NModal),
NAeroStates(NAero),
NGust(Gust),
pA(pAMat), pB(pBMat), pC(pCMat),
pD0(pD0Mat), pD1(pD1Mat), pD2(pD2Mat),
pq(0), pqPrime(0), pqSec(0),
pxa(0), pxaPrime(0),
pgs(0), pgsPrime(0),
gustVff(Vff), gustXi(0.707),
RigidF(rgF)
{
	DEBUGCOUTFNAME("AerodynamicModal::AerodynamicModal");

	R0 = pModalNode->GetRCurr()*Ra;
	P0 = R0.MulTV(pModalNode->GetXCurr());

	SAFENEWWITHCONSTRUCTOR(pq, MyVectorHandler,
		MyVectorHandler(NStModes+RigidF));
	SAFENEWWITHCONSTRUCTOR(pqPrime, MyVectorHandler,
		MyVectorHandler(NStModes+RigidF));
	SAFENEWWITHCONSTRUCTOR(pqSec, MyVectorHandler,
		MyVectorHandler(NStModes+RigidF));
	SAFENEWWITHCONSTRUCTOR(pxa, MyVectorHandler,
		MyVectorHandler(NAeroStates));
	SAFENEWWITHCONSTRUCTOR(pxaPrime, MyVectorHandler,
		MyVectorHandler(NAeroStates));
	SAFENEWWITHCONSTRUCTOR(pgs, MyVectorHandler,
		MyVectorHandler(2*NGust));
	SAFENEWWITHCONSTRUCTOR(pgsPrime, MyVectorHandler,
		MyVectorHandler(2*NGust));

	ASSERT(pModalNode != 0);
	ASSERT(pModalNode->GetNodeType() == Node::STRUCTURAL);
}

AerodynamicModal::~AerodynamicModal(void)
{
	DEBUGCOUTFNAME("AerodynamicModal::~AerodynamicModal");
	if (pq != 0) {
		SAFEDELETE(pq);
	}
	if (pqPrime != 0) {
		SAFEDELETE(pqPrime);
	}
	if (pqSec != 0) {
		SAFEDELETE(pqSec);
	}
	if (pxa != 0) {
		SAFEDELETE(pxa);
	}
	if (pxaPrime != 0) {
		SAFEDELETE(pxaPrime);
	}
	if (pA != 0) {
		SAFEDELETE(pA);
	}
	if (pB != 0) {
		SAFEDELETE(pB);
	}
	if (pC != 0) {
		SAFEDELETE(pC);
	}
	if (pD0 != 0) {
		SAFEDELETE(pD0);
	}
	if (pD1 != 0) {
		SAFEDELETE(pD1);
	}
	if (pD2 != 0) {
		SAFEDELETE(pD2);
	}
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
AerodynamicModal::Restart(std::ostream& out) const
{
	return out << "  /* aerodynamic modal: not implemented yet */" << std::endl;
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
AerodynamicModal::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr )
{
	DEBUGCOUT("Entering AerodynamicModal::AssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	integer iModalIndex = pModalJoint->iGetModalIndex();
	for (unsigned int iCnt = 1; iCnt <= NStModes; iCnt++) {
		WM.PutRowIndex(iCnt, iModalIndex+NStModes+iCnt);
		WM.PutColIndex(iCnt, iModalIndex+iCnt);
		WM.PutColIndex(NStModes+iCnt, iModalIndex+NStModes+iCnt);
	}

	integer iAeroIndex = iGetFirstIndex();
	for (unsigned int iCnt = 1; iCnt <= NAeroStates+2*NGust; iCnt++) {
		WM.PutRowIndex(NStModes+iCnt, iAeroIndex+iCnt);
		WM.PutColIndex(2*NStModes+iCnt, iAeroIndex+iCnt);
	}

	/* Dati del nodo rigido */
	const Vec3& X0(pModalNode->GetXCurr());
	const Vec3& V0(pModalNode->GetVCurr());
	const Mat3x3& Rn(pModalNode->GetRCurr());
	Mat3x3 RR(Rn*Ra);
	Vec3 Vr(V0);

	doublereal rho, vs, p, T;
	GetAirProps(X0, rho, vs, p, T);		/* p, T are not used yet */

	Vec3 VTmp(Zero3);
	if (fGetAirVelocity(VTmp, X0)) {
		Vr -= VTmp;
	}

	/* velocità nel riferimento nodale aerodinamico */
	VTmp = RR.MulTV(Vr);
	doublereal nV = std::abs(VTmp(1));
	/* doublereal CV=Chord/(2*nV); */
	doublereal qd  = 0.5*rho*nV*nV;
	doublereal qd1 = 0.25*rho*nV*Chord; /* qd*CV */
	doublereal qd2 = 0.125*rho*Chord*Chord; /* qd*CV*CV */

	/* parte deformabile :
	 *
	 * |                                ||aP|
	 * |                                ||  |
	 * |   cKae     Mae+cCae     -cqC   ||bP|
	 * |                     	    ||  |
	 * |-c(1/CV)B     0      I-c(1/CV)A ||xa|
	 *
	 * con
	 * Mae=-qd CV^2 D2
	 * Cae=-qd CV D1
	 * Kae=-qd D0
	 */
	for (unsigned int i = 1; i <= NStModes; i++) {
		for (unsigned int j = 1; j <= NStModes; j++) {
			WM.IncCoef(i, j,
				-dCoef*qd*pD0->operator()(RigidF + i, RigidF + j));
			WM.IncCoef(i, j + NStModes,
				-qd2*pD2->operator()(RigidF + i, RigidF + j)
				-dCoef*qd1*pD1->operator()(RigidF + i, RigidF + j));
		}
	}

	for (unsigned int j = 1; j <= NAeroStates; j++) {
		for (unsigned int i = 1; i <= NStModes; i++) {
			WM.IncCoef(i, j + 2*NStModes,
				-dCoef*qd*pC->operator()(RigidF + i, j));
			WM.IncCoef(j + NStModes, i,
				-dCoef*(2*nV/Chord)*pB->operator()(j, RigidF + i));
		}
	}

	for (unsigned int j = 1; j <= NAeroStates; j++) {
		for (unsigned int i = 1; i <= NAeroStates; i++) {
			WM.IncCoef(i + NStModes, j + 2*NStModes,
				1.*(i == j) - dCoef*(2*nV/Chord)*pA->operator()(i, j));
		}
	}

	if (NGust) {
		for (unsigned int i = 1; i <= NStModes; i++) {
			WM.IncCoef(i, 2*NStModes + NAeroStates + 1,
				-dCoef*qd*pD0->operator()(RigidF + i, RigidF + NStModes + 1));
			WM.IncCoef(i, 2*NStModes + NAeroStates + 3,
				-dCoef*qd*pD0->operator()(RigidF + i, RigidF + NStModes + 2));
			WM.IncCoef(i, 2*NStModes + NAeroStates + 2,
				-qd2*pD2->operator()(RigidF + i, RigidF + NStModes + 1)
				-dCoef*qd1*pD1->operator()(RigidF + i, RigidF + NStModes + 1));
			WM.IncCoef(i, 2*NStModes + NAeroStates + 4,
				-qd2*pD2->operator()(RigidF + i, RigidF + NStModes + 2)
				-dCoef*qd1*pD1->operator()(RigidF + i, RigidF + NStModes + 2));
		}

		for (unsigned int i = 1; i <= NAeroStates; i++) {
			WM.IncCoef(i + NStModes, 2*NStModes + NAeroStates + 1,
				-dCoef*(2*nV/Chord)*pB->operator()(i, RigidF + NStModes + 1));
			WM.IncCoef(i + NStModes, 2*NStModes + NAeroStates + 3,
				-dCoef*(2*nV/Chord)*pB->operator()(i, RigidF + NStModes + 2));
		}

		for (unsigned int i = 0; i < NGust; i++) {
			WM.IncCoef(i*2 + 1 + NStModes+NAeroStates,
				i*2 + 1 + 2*NStModes + NAeroStates, 1.);
			WM.IncCoef(i*2 + 1 + NStModes+NAeroStates,
				i*2 + 2 + 2*NStModes + NAeroStates, -dCoef);
			WM.IncCoef(i*2 + 2 + NStModes+NAeroStates,
				i*2 + 1 + 2*NStModes+NAeroStates,
				dCoef*gustVff*gustVff);
			WM.IncCoef(i*2 + 2 + NStModes + NAeroStates,
				i*2 + 2 + 2*NStModes + NAeroStates,
				1. + dCoef*2*gustXi*gustVff);
		}
	}

	return WorkMat;
}


SubVectorHandler&
AerodynamicModal::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("AerodynamicModal::AssRes");
	WorkVec.ResizeReset(RigidF + NStModes + NAeroStates + 2*NGust);

	if (RigidF) {
		integer iFirstIndex = pModalNode->iGetFirstMomentumIndex();
		for (int iCnt = 1; iCnt <= RigidF; iCnt++) {
			WorkVec.PutRowIndex(iCnt, iFirstIndex + iCnt);
		}

		const Vec3& X0(pModalNode->GetXCurr());
		const Mat3x3& Rn(pModalNode->GetRCurr());

		Mat3x3 RR(Rn*Ra);

		// q
		pq->Put(1, RR.MulTV(X0 - P0));

		// FIXME: use orientation vector?
		Vec3 g(CGR_Rot::Param, RR.MulMT(R0));
		doublereal d(g.Norm());
		if (d > std::numeric_limits<doublereal>::epsilon()) {
			pq->Put(4, g*(2./d*atan(d/2.)));

		} else {
			pq->Put(4, Zero3);
		}

		// dot{q}
		const Vec3& V0(pModalNode->GetVCurr());
		const Vec3& W0(pModalNode->GetWCurr());
		pqPrime->Put(1, RR.MulTV(V0 + X0.Cross(W0)));
		pqPrime->Put(4, RR.MulTV(W0));

		// ddot{q}
		const Vec3& XPP0(pModalNode->GetXPPCurr());
		const Vec3& WP0(pModalNode->GetWPCurr());
		pqSec->Put(1, RR.MulTV(XPP0));	// verificare?
		pqSec->Put(4, RR.MulTV(WP0));
	}

	integer iModalIndex = pModalJoint->iGetModalIndex();

	for (unsigned int iCnt = 1; iCnt <= NStModes; iCnt++) {
		WorkVec.PutRowIndex(RigidF + iCnt, iModalIndex + NStModes + iCnt);
	}

	integer iAeroIndex = iGetFirstIndex();
	for (unsigned int iCnt = 1; iCnt <= NAeroStates + 2*NGust; iCnt++) {
		WorkVec.PutRowIndex(RigidF + NStModes+iCnt, iAeroIndex + iCnt);
	}

	/* Recupera i vettori {a} e {aP} e {aS}(deformate modali) */
	// FIXME: get this info from modal joint?
	for (unsigned int  iCnt = 1; iCnt <= NStModes; iCnt++) {
		pq->PutCoef(iCnt + RigidF, XCurr(iModalIndex + iCnt));
		pqPrime->PutCoef(iCnt + RigidF, XPrimeCurr(iModalIndex + iCnt));
		pqSec->PutCoef(iCnt + RigidF, XPrimeCurr(iModalIndex + NStModes + iCnt));
	}

	/* Recupera i vettori {xa} e {xaP}  */
	for (unsigned int  iCnt = 1; iCnt <= NAeroStates; iCnt++) {
		pxa->PutCoef(iCnt, XCurr(iAeroIndex + iCnt));
		pxaPrime->PutCoef(iCnt, XPrimeCurr(iAeroIndex + iCnt));
	}

	for (unsigned int  iCnt = 1; iCnt <= 2*NGust; iCnt++) {
		pgs->PutCoef(iCnt, XCurr(iAeroIndex + NAeroStates + iCnt));
		pgsPrime->PutCoef(iCnt, XPrimeCurr(iAeroIndex + NAeroStates + iCnt));
	}

	AssVec(WorkVec);

	return WorkVec;
}

SubVectorHandler&
AerodynamicModal::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	DEBUGCOUTFNAME("AerodynamicModal::AssRes");
	const Vec3& X0(pModalNode->GetXCurr());
	const Mat3x3& Rn(pModalNode->GetRCurr());
	R0 = Rn*Ra;
	P0 = R0.MulTV(X0);

	WorkVec.ResizeReset(RigidF + NStModes + NAeroStates + 2*NGust);

	if (RigidF) {
		integer iFirstIndex = pModalNode->iGetFirstIndex();
		for (int iCnt = 1; iCnt <= RigidF; iCnt++) {
			WorkVec.PutRowIndex(iCnt, iFirstIndex + RigidF + iCnt);
		}

		const Vec3& X0(pModalNode->GetXCurr());
		const Mat3x3& Rn(pModalNode->GetRCurr());
		Mat3x3 RR(Rn*Ra);

		// q
		pq->Put(1, RR.MulTV(X0) - P0);

		Vec3 g(CGR_Rot::Param, RR.MulMT(R0));
		doublereal d(g.Norm());
		if (d > std::numeric_limits<doublereal>::epsilon()) {
			pq->Put(4, g*(2./d*atan(d/2.)));

		} else {
			pq->Put(4, Zero3);
		}

		// dot{q}
		const Vec3& V0(pModalNode->GetVCurr());
		const Vec3& W0(pModalNode->GetWCurr());
		pqPrime->Put(1, RR.MulTV(V0));
		pqPrime->Put(4, RR.MulTV(W0));

		// ddot{q}
		pqPrime->Put(1, Zero3);
		pqPrime->Put(4, Zero3);
	}

	integer iModalIndex = pModalJoint->iGetModalIndex();

	for (unsigned int iCnt = 1; iCnt <= NStModes; iCnt++) {
		WorkVec.PutRowIndex(RigidF + iCnt, iModalIndex + NStModes + iCnt);
	}

	integer iAeroIndex = iGetFirstIndex();
	for (unsigned int iCnt = 1; iCnt <= NAeroStates + 2*NGust; iCnt++) {
		WorkVec.PutRowIndex(RigidF + NStModes + iCnt, iAeroIndex + iCnt);
	}

	/* Recupera i vettori {a} e {aP} e {aS} (deformate modali) */
	for (unsigned int iCnt = 1; iCnt <= NStModes; iCnt++) {
		pq->PutCoef(iCnt + RigidF, XCurr(iModalIndex + iCnt));
		pqPrime->PutCoef(iCnt + RigidF, 0);
		pqSec->PutCoef(iCnt + RigidF, 0);
	}

	/* Recupera i vettori {xa} e {xaP}  */
	for (unsigned int iCnt = 1; iCnt <= NAeroStates; iCnt++) {
		pxa->PutCoef(iCnt, XCurr(iAeroIndex + iCnt));
		pxaPrime->PutCoef(iCnt, 0);
	}

	for (unsigned int iCnt = 1; iCnt <= 2*NGust; iCnt++) {
		pgs->PutCoef(iCnt, XCurr(iAeroIndex + NAeroStates + iCnt));
		pgsPrime->PutCoef(iCnt, 0.);
	}

	AssVec(WorkVec);

	return WorkVec;
}

/* assemblaggio residuo */
void
AerodynamicModal::AssVec(SubVectorHandler& WorkVec)
{
	DEBUGCOUTFNAME("AerodynamicModal::AssVec");

	/* Dati del nodo rigido */
	const Vec3& X0(pModalNode->GetXCurr());
	Vec3 V0(pModalNode->GetVCurr());
	const Mat3x3& Rn(pModalNode->GetRCurr());
	const Mat3x3& RR(Rn*Ra);
	Vec3 Vr(V0);

	doublereal rho, vs, p, T;
	GetAirProps(X0, rho, vs, p, T);	/* p, T are not used yet */

	Vec3 Vair(Zero3);
	if (fGetAirVelocity(Vair, X0)) {
		Vr -= Vair;
	}
	/* velocita' nel riferimento nodale aerodinamico */
	V0 = RR.MulTV(Vr);
	doublereal nV = std::abs(V0(1));
	MyVectorHandler TmpA(NAeroStates);
	TmpA.Reset();
	MyVectorHandler Tmp(RigidF + NStModes);
	Tmp.Reset();
	MyVectorHandler TmpP(RigidF + NStModes);
	TmpP.Reset();
	MyVectorHandler TmpS(RigidF + NStModes);
	TmpS.Reset();

	// FIXME: use matrix/vector product?
	for (unsigned int i = 1; i <= NAeroStates; i++) {
		for (unsigned int j = 1; j <= NStModes + RigidF; j++) {
			TmpA.IncCoef(i, pB->operator()(i, j)*(pq->operator()(j)));
		}
	}

	pA->MatVecIncMul(TmpA, *pxa);

	/* doublereal CV = Chord/(2*nV); */
	doublereal qd = 0.5*rho*nV*nV;
	doublereal qd1 = 0.25*rho*nV*Chord; /* qd*CV */
	doublereal qd2 = 0.125*rho*Chord*Chord; /* qd*CV*CV */

	for (unsigned int i = 1; i <= NAeroStates; i++) {
		WorkVec.IncCoef(RigidF + NStModes + i, -pxaPrime->operator()(i) + (2*nV/Chord)*(TmpA(i)));
	}

	for (unsigned int i = 1; i <= NStModes + RigidF; i++) {
		for (unsigned int j = 1; j <= NStModes + RigidF; j++) {
			Tmp.IncCoef(i, pD0->operator()(i, j)*(pq->operator()(j)));
			TmpP.IncCoef(i, pD1->operator()(i, j)*(pqPrime->operator()(j)));
			TmpS.IncCoef(i, pD2->operator()(i, j)*(pqSec->operator()(j)));
		}
	}

	for (unsigned int i = 1; i <= NStModes + RigidF; i++) {
		for (unsigned int j = 1; j <= NAeroStates; j++) {
			Tmp.IncCoef(i, pC->operator()(i, j)*(pxa->operator()(j)));
		}
	}

	if (RigidF) {
		Vec3 F(Zero3);
		Vec3 M(Zero3);

		for (unsigned int i = 1; i <= 3; i++) {
			F.Put(i, qd*Tmp(i) + qd1*TmpP(i) + qd2*TmpS(i));
			M.Put(i, qd*Tmp(i + 3) + qd1*TmpP(i + 3) + qd2*TmpS(i + 3));
		}

		// std::cout << F << std::endl;
		F = RR*(-F);
		M = RR*M;
		WorkVec.Add(1, F);
		WorkVec.Add(4, M);
	}

	for (unsigned int i = 1 + RigidF; i <= NStModes + RigidF; i++) {
		WorkVec.IncCoef(i, qd*Tmp(i) + qd1*TmpP(i) + qd2*TmpS(i));
	}

	if (NGust) {
#if 0 /* unused */
		doublereal Vyg = Vair(2)/nV;
		doublereal Vzg = Vair(3)/nV;
#endif
		for (unsigned int i = 1; i <= NAeroStates; i++) {
			WorkVec.IncCoef(RigidF + NStModes + i,
				(2*nV/Chord)*pB->operator()(i, RigidF + NStModes + 1)*pgs->operator()(1)
				+ (2*nV/Chord)*pB->operator()(i, RigidF + NStModes + 2)*pgs->operator()(3));
		}

		if (RigidF) {
			Vec3 F(Zero3);
			Vec3 M(Zero3);

			for (unsigned int i = 1; i <= 3; i++) {
				F.Put(i, qd*pD0->operator()(i, RigidF + NStModes + 1)*pgs->operator()(1)
					+ qd*pD0->operator()(i, RigidF + NStModes + 2)*pgs->operator()(3)
					+ qd1*pD1->operator()(i, RigidF + NStModes + 1)*pgs->operator()(2)
					+ qd1*pD1->operator()(i, RigidF + NStModes + 2)*pgs->operator()(4)
					+ qd2*pD2->operator()(i, RigidF + NStModes + 1)*pgsPrime->operator()(2)
					+ qd2*pD2->operator()(i, RigidF + NStModes + 2)*pgsPrime->operator()(4));
				M.Put(i, qd*pD0->operator()(i + 3, RigidF + NStModes + 1)*pgs->operator()(1)
					+ qd*pD0->operator()(i + 3, RigidF + NStModes + 2)*pgs->operator()(3)
					+ qd1*pD1->operator()(i + 3, RigidF + NStModes + 1)*pgs->operator()(2)
					+ qd1*pD1->operator()(i + 3, RigidF + NStModes + 2)*pgs->operator()(4)
					+ qd2*pD2->operator()(i + 3, RigidF + NStModes + 1)*pgsPrime->operator()(2)
					+ qd2*pD2->operator()(i + 3, RigidF + NStModes + 2)*pgsPrime->operator()(4));
			}

			// std::cout << F << std::endl;
			F = RR*(-F);
			M = RR*M;
			WorkVec.Add(1, F);
			WorkVec.Add(4, M);
		}

		for (unsigned int i = 1 + RigidF; i <= NStModes + RigidF; i++) {
			WorkVec.IncCoef(i, qd*pD0->operator()(i, RigidF + NStModes + 1)*pgs->operator()(1)
				+ qd*pD0->operator()(i, RigidF + NStModes + 2)*pgs->operator()(3)
				+ qd1*pD1->operator()(i, RigidF + NStModes + 1)*pgs->operator()(2)
				+ qd1*pD1->operator()(i, RigidF + NStModes + 2)*pgs->operator()(4)
				+ qd2*pD2->operator()(i, RigidF + NStModes + 1)*pgsPrime->operator()(2)
				+ qd2*pD2->operator()(i, RigidF + NStModes + 2)*pgsPrime->operator()(4));
		}

		for (unsigned int i = 0; i < NGust; i++) {
			WorkVec.IncCoef(RigidF + i*2 + 1 + NStModes + NAeroStates,
				-pgsPrime->operator()(1 + i*2) + pgs->operator()(2 + i*2));
			if (nV != 0) {
				WorkVec.IncCoef(RigidF + i*2 + 2 + NStModes + NAeroStates,
					-pgsPrime->operator()(2 + i*2)
					- gustVff*gustVff*pgs->operator()(1 + i*2)
					- 2*gustXi*gustVff*pgs->operator()(2 + i*2)
					+ gustVff*gustVff*Vair(2 + i)/nV);

			} else {
				WorkVec.IncCoef(RigidF + i*2 + 2 + NStModes + NAeroStates,
					-pgsPrime->operator()(2 + i*2)
					- gustVff*gustVff*pgs->operator()(1 + i*2)
					- 2*gustXi*gustVff*pgs->operator()(2 + i*2));
			}
		}

		// std::cout << X0 << std::endl;
		// std::cout << WorkVec(3) << std::endl;
	}
}

/* output; si assume che ogni tipo di elemento sappia, attraverso
 * l'OutputHandler, dove scrivere il proprio output */
void
AerodynamicModal::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		OH.AeroModals() << std::setw(8) << GetLabel() << " ";
		for (unsigned int iCnt = 1; iCnt <= NAeroStates; iCnt++) {
			OH.AeroModals() << " " << pxa->operator()(iCnt);
		}
		for (unsigned int iCnt = 1; iCnt <= NAeroStates + 2*NGust; iCnt++) {
			OH.AeroModals() << " " << pxaPrime->operator()(iCnt);
		}
		OH.AeroModals() << std::endl;
	}
}

/* AerodynamicModal - end */

Elem *
ReadAerodynamicModal(DataManager* pDM,
	MBDynParser& HP,
	const DofOwner* pDO,
	unsigned int uLabel)
{
	DEBUGCOUTFNAME("ReadAerodynamicModal(" << uLabel << ")" << std::endl);

	/* formato dell'input:
	 *
	 *  label,
	 *  modal node,
	 *  reference modal joint,
	 *  rotation of the aerodynamic reference in respect of nodal reference,
	 *  reference cord,
	 *  number of aerodynamic states,
	 *  {rigid, gust, gust filter cut-off frequency}
	 *  file name containing state space model matrices;
	 */

	/* giunto modale collegato */
	const Modal *pModalJoint = pDM->ReadElem<const Modal, const Joint, Elem::JOINT>(HP);
	ReferenceFrame RF;
	const StructNode *pModalNode = pModalJoint->pGetModalNode();
	if (pModalNode) {
		RF = ReferenceFrame(pModalNode);

	} else {
		RF = AbsRefFrame;
	}

	Mat3x3 Ra(HP.GetRotRel(RF));

	doublereal Chord = HP.GetReal();

	unsigned int AeroN = HP.GetInt();
	/* numero modi e FEM */
	unsigned int NModes = pModalJoint->uGetNModes();

	AerodynamicModal::RigidF_t rigidF = AerodynamicModal::NO_RIGID;
	if (HP.IsKeyWord("rigid")) {
		rigidF = AerodynamicModal::RIGID;
		NModes += rigidF;
	}

	/* Eventuale raffica */
	unsigned int GustN = 0;
	doublereal Vff = 0.;
	if (HP.IsKeyWord("gust")) {
		GustN = 2;
		Vff = HP.GetReal();
	}

	/* apre il file contenente le matrici A B C D0 D1 D2 */
	const char *sFileData = HP.GetFileName();
	std::ifstream fdat(sFileData);
	DEBUGCOUT("Reading Aerodynamic State Space Matrices from file '"
		<< sFileData << '\'' << std::endl);
	if (!fdat) {
		int save_errno = errno;
		silent_cerr(std::endl << "Unable to open file \"" << sFileData << "\" "
			"(" << save_errno << ": " << strerror(save_errno) << ")" << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	SpMapMatrixHandler* pAMat = 0;
	FullMatrixHandler* pBMat = 0;
	FullMatrixHandler* pCMat = 0;
	FullMatrixHandler* pD0Mat = 0;
	FullMatrixHandler* pD1Mat = 0;
	FullMatrixHandler* pD2Mat = 0;
	SAFENEWWITHCONSTRUCTOR(pAMat, SpMapMatrixHandler, SpMapMatrixHandler(AeroN, AeroN));
	pAMat->Reset();
	SAFENEWWITHCONSTRUCTOR(pBMat, FullMatrixHandler, FullMatrixHandler(AeroN, NModes + GustN));
	SAFENEWWITHCONSTRUCTOR(pCMat, FullMatrixHandler, FullMatrixHandler(NModes, AeroN));
	SAFENEWWITHCONSTRUCTOR(pD0Mat, FullMatrixHandler, FullMatrixHandler(NModes, NModes + GustN));
	SAFENEWWITHCONSTRUCTOR(pD1Mat, FullMatrixHandler, FullMatrixHandler(NModes, NModes + GustN));
	SAFENEWWITHCONSTRUCTOR(pD2Mat, FullMatrixHandler, FullMatrixHandler(NModes, NModes + GustN));

	doublereal d;
	char str[BUFSIZ];

	/* parsing del file */
	while (!fdat.eof()) {
		fdat.getline(str, sizeof(str));

		/* legge il primo blocco (HEADER) */
		if (!strncmp("*** MATRIX A", str, STRLENOF("*** MATRIX A"))) {
			for (unsigned int iCnt = 1; iCnt <= AeroN; iCnt++)  {
				for (unsigned int jCnt = 1; jCnt <= AeroN; jCnt++)  {
					fdat >> d;
					pAMat->PutCoef(iCnt, jCnt, d);
				}
			}

		/* legge il primo blocco (HEADER) */
		} else if (!strncmp("*** MATRIX B", str, STRLENOF("*** MATRIX B"))) {
			for (unsigned int iCnt = 1; iCnt <= AeroN; iCnt++)  {
				for (unsigned int jCnt = 1; jCnt <= NModes + GustN; jCnt++)  {
					fdat >> d;
					pBMat->PutCoef(iCnt, jCnt, d);
				}
			}

		/* legge il primo blocco (HEADER) */
		} else if (!strncmp("*** MATRIX C", str, STRLENOF("*** MATRIX C"))) {
			for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++)  {
				for (unsigned int jCnt = 1; jCnt <= AeroN; jCnt++)  {
					fdat >> d;
					pCMat->PutCoef(iCnt, jCnt, d);
				}
			}

		/* legge il primo blocco (HEADER) */
		} else if (!strncmp("*** MATRIX D0", str, STRLENOF("*** MATRIX D0"))) {
			for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++)  {
				for (unsigned int jCnt = 1; jCnt <= NModes + GustN; jCnt++)  {
					fdat >> d;
					pD0Mat->PutCoef(iCnt, jCnt, d);
				}
			}

		/* legge il primo blocco (HEADER) */
		} else if (!strncmp("*** MATRIX D1", str, STRLENOF("*** MATRIX D1"))) {
			for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++)  {
				for (unsigned int jCnt = 1; jCnt <= NModes + GustN; jCnt++)  {
					fdat >> d;
					pD1Mat->PutCoef(iCnt, jCnt, d);
				}
			}

		/* legge il primo blocco (HEADER) */
		} else if (!strncmp("*** MATRIX D2", str, STRLENOF("*** MATRIX D2"))) {
			for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++)  {
				for (unsigned int jCnt = 1; jCnt <= NModes + GustN; jCnt++)  {
					fdat >> d;
					pD2Mat->PutCoef(iCnt, jCnt, d);
				}
			}
		}
	}
	fdat.close();

	flag fOut = pDM->fReadOutput(HP, Elem::AEROMODAL);

	Elem* pEl = 0;
	SAFENEWWITHCONSTRUCTOR(pEl,
		AerodynamicModal,
		AerodynamicModal(uLabel, pModalNode, pModalJoint,
			Ra, pDO, Chord, NModes - rigidF,
			AeroN, rigidF, GustN, Vff,
			pAMat, pBMat, pCMat, pD0Mat, pD1Mat, pD2Mat, fOut));

	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		silent_cerr("semicolon expected at line "
			<< HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pEl;
} /* End of DataManager::ReadAerodynamicModal() */

