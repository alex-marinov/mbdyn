/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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

/*
   classe per l'introduzione della flessibilita' modale nel codice multi-corpo

   29/7/99:  implementata la parte di corpo rigido (solo massa e momenti
             d'inerzia) verificata con massa isolata soggetta a forze
	     e momenti esterni

   13/9/99:  corpo rigido OK (anche momenti statici ed eq. di vincolo)
             corpo deformabile isolato OK (verificato con osc. armonico)

   23/9/99:  corpo deformabile + eq. di vincolo (no moto rigido)  OK
   	     (verificato con trave incastrata)
             corpo rigido + deformabile OK (verificato con trave libera)

     10/99:  validazioni con NASTRAN: trave incastrata, trave libera
             con forza in mezzeria, trave libera con forza all'estremita'

   4/11/99:  completate le funzioni InitialAssJac e InitialAssRes
             per l'assemblaggio iniziale del 'vincolo'

  22/11/99:  modificata la funzione ReadModal per leggere il file
             generato da DADS

  26/11/99:  validazione con piastra vincolata con elementi elastici

     12/99:  validazione con piastra e bauchau

   01/2000:  modi rotanti

  30/11/99:  aggiunto lo smorzamento strutturale tra i dati d'ingresso
             aggiunte le inerzie concentrate

   1/12/99:  modifiche alla lettura dati d'ingresso (l'estensione *.fem
             al file con i dati del modello ad elementi viene aggiunta
	     automaticamente dal programma, che crea un file di
             output con l'estensione *.mod)
             corretto un bug nella scrittura delle equazioni di vincolo

  17/12/99:  aggiunta la possibilita' di definire uno smorzamento strutturale
             variabile con le frequenze

  23/12/99:  nuova modifica alla lettura dati di ingresso

10/01/2000:  introdotta la possibilita' di definire un fattore di scala
             per i dati del file d'ingresso

22/01/2000:  tolto il fattore di scala perche' non funziona

   03/2000:  aggiunte funzioni che restituiscono dati dell'elemento
             (autovettori, modi ecc.) aggiunta la possibilita'
	     di imporre delle deformate modali iniziali

   Modifiche fatte al resto del programma:

   Nella classe joint    : aggiunto il vincolo modale
   Nella classe strnode  : aggiunto il nodo modale
   Nella classe MatVec3n : aggiunte classi e funzioni per gestire matrici
   			   Nx3 e NxN
   Nella classe SubMat   : corretto un bug nelle funzioni Add Mat3xN ecc.,
                           aggiunte funzioni per trattare sottomatrici Nx3
*/

/*
 * Copyright 1999-2004 Felice Felippone <ffelipp@tin.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */

/*
 * Copyright 1999-2004 Pierangelo Masarati  <masarati@aero.polimi.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 *
 * Modified by Pierangelo Masarati
 */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

/* FIXME: gravity in modal elements is eXperimental; undefine to disable */
#define MODAL_USE_GRAVITY

#include <modal.h>
#include <dataman.h>

Modal::Modal(unsigned int uL,
		const ModalNode* pR,
		const Vec3& x0,
		const Mat3x3& R0,
		const DofOwner* pDO,
		unsigned int NM,         /* numero modi */
		unsigned int NI,         /* numero nodi d'interfaccia */
		unsigned int NF,         /* numero nodi FEM */
		doublereal dMassTmp,     /* inv. inerzia (m, m.stat., d'in.) */
		const Vec3& STmp,
		const Mat3x3& JTmp,
		MatNxN *pGenMass,
		MatNxN *pGenStiff,
		MatNxN *pGenDamp,
		unsigned int *IdFemNodes, /* label nodi FEM */
		unsigned int *intFEMNodes,/* label nodi FEM d'interfaccia */
		unsigned int *intMBNodes, /* label nodi MB d'interfaccia */
		Mat3xN *pN,               /* posizione dei nodi FEM */
		Mat3xN *pOffsetfemNodes,
		Mat3xN *pOffsetmbNodes,
		Mat3xN *pRotmbNodes,
		const StructNode** pIN,   /* nodi d'interfaccia */
		Mat3xN *pPHItStrNode,     /* forme modali nodi d'interfaccia */
		Mat3xN *pPHIrStrNode,
		Mat3xN *pModeShapest,     /* autovettori: servono a aeromodal */
		Mat3xN *pModeShapesr,
		Mat3xN *pInv3,            /* invarianti d'inerzia I3...I11 */
		Mat3xN *pInv4,
		Mat3xN *pInv5,
		Mat3xN *pInv8,
		Mat3xN *pInv9,
		Mat3xN *pInv10,
		Mat3xN *pInv11,
		VecN *aa,
		VecN *bb,
		const char *sFileMod,
		DataManager* pDM,  /* non serve */
		MBDynParser& HP,   /* non serve */
		flag fOut)
: Elem(uL, Elem::JOINT, fOut),
Joint(uL, Joint::MODAL, pDO, fOut),
pModalNode(pR),
iRigidOffset(pModalNode ? 12 : 0),
x(x0),
R(R0),
RT(R0.Transpose()),
NModes(NM),
NStrNodes(NI),
NFemNodes(NF),
dMass(dMassTmp),
Inv2(STmp),
Inv7(JTmp),
pModalMass(pGenMass),
pModalStiff(pGenStiff),
pModalDamp(pGenDamp),
IdFemNodes(IdFemNodes),
IntFEMNodes(intFEMNodes),
IntMBNodes(intMBNodes),
pXYZFemNodes(pN),
pOffsetFEMNodes(pOffsetfemNodes),
pOffsetMBNodes(pOffsetmbNodes),
pRotMBNodes(pRotmbNodes),
pInterfaceNodes(pIN),
pPHIt(pPHItStrNode),
pPHIr(pPHIrStrNode),
pModeShapest(pModeShapest),
pModeShapesr(pModeShapesr),
pCurrXYZ(NULL),
pCurrXYZVel(NULL),
pInv3(pInv3),
pInv4(pInv4),
pInv5(pInv5),
pInv8(pInv8),
pInv9(pInv9),
pInv10(pInv10),
pInv11(pInv11),
Inv3jaj(0.),
Inv3jaPj(0.),
Inv8jaj(0.),
Inv8jaPj(0.),
Inv5jaj(NModes, 0.),
Inv5jaPj(NModes, 0.),
Inv4j(0.),
VInv5jaj(0.),
VInv5jaPj(0.),
Inv8jTranspose(0.),
Inv9jkak(0.),
#ifdef MODAL_USE_INV9
Inv9jkajaPk(0.),
#endif /* MODAL_USE_INV9 */
a(*aa),
aPrime(NModes, 0.),
b(*bb),
bPrime(NModes, 0.),
pd1tot(NULL),
pd2(NULL),
pR1tot(NULL),
pR2(NULL),
pF(NULL),
pM(NULL),
fOutFlex(sFileMod, std::ios::out /* | std::ios::noreplace */ )
{
	ASSERT(pModalNode->GetStructNodeType() == StructNode::MODAL);

	SAFENEWARR(pd1tot, Vec3, NStrNodes);
	SAFENEWARR(pd2, Vec3, NStrNodes);
	SAFENEWARR(pR1tot, Mat3x3, NStrNodes);
	SAFENEWARR(pR2, Mat3x3, NStrNodes);
	SAFENEWARR(pF, Vec3, NStrNodes);
	SAFENEWARR(pM, Vec3, NStrNodes);

	if (!fOutFlex) {
		silent_cerr("Modal(" << GetLabel() << "): "
			"unable to open output file \"" << sFileMod << "\""
			<< std::endl);
		throw ErrGeneric();
	}

	for (unsigned int i = 0; i < NStrNodes; i++) {
		pd1tot[i] = Zero3;
		pd2[i] = Zero3;
		pR1tot[i] = Eye3;
		pR2[i] = Eye3;
		pF[i] = Zero3;
		pM[i] = Zero3;
	}
}

Modal::~Modal(void)
{
	if (pd1tot != NULL) {
		SAFEDELETEARR(pd1tot);
	}
	if (pd2 != NULL) {
		SAFEDELETEARR(pd2);
	}
	if (pF != NULL) {
		SAFEDELETEARR(pF);
	}
	if (pM != NULL) {
		SAFEDELETEARR(pM);
	}
	if (pR1tot != NULL) {
		SAFEDELETEARR(pR1tot);
	}
	if (pR2 != NULL) {
		SAFEDELETEARR(pR2);
	}

	/* FIXME: destroy all the other data ... */
}

Joint::Type
Modal::GetJointType(void) const
{
	return Joint::MODAL;
}

std::ostream&
Modal::Restart(std::ostream& out) const
{
	return out << "modal; # not implemented yet" << std::endl;
}

unsigned int
Modal::iGetNumDof(void) const
{
	/* i gradi di liberta' propri sono:
	 *   2*M per i modi
	 *   6*N per le reazioni vincolari dei nodi d'interfaccia */
	return 2*NModes + 6*NStrNodes;
}

std::ostream&
Modal::DescribeDof(std::ostream& out, char *prefix, bool bInitial, int i) const
{
	integer iModalIndex = iGetFirstIndex();

	if (i >= 0) {
		silent_cerr("Modal(" << GetLabel() << "): "
			"DescribeDof(" << i << ") "
			"not implemented yet" << std::endl);
		throw ErrGeneric();
	}

	out 
		<< prefix << iModalIndex + 1 << "->" << iModalIndex + NModes
		<< ": modal displacements [q(1.." << NModes << ")]" << std::endl
		<< prefix << iModalIndex + NModes + 1 << "->" << iModalIndex + 2*NModes
		<< ": modal velocities [qP(1.." << NModes << ")]" << std::endl;
	iModalIndex += 2*NModes;
	for (unsigned iStrNodem1 = 0; iStrNodem1 < NStrNodes; iStrNodem1++, iModalIndex += 6) {
		out
			<< prefix << iModalIndex + 1 << "->" << iModalIndex + 3 << ": "
				"StructNode(" << pInterfaceNodes[iStrNodem1]->GetLabel() << ") reaction forces [Fx,Fy,Fz]" << std::endl
			<< prefix << iModalIndex + 4 << "->" << iModalIndex + 6 << ": "
				"StructNode(" << pInterfaceNodes[iStrNodem1]->GetLabel() << ") reaction couples [mx,my,mz]" << std::endl;
		if (bInitial) {
			iModalIndex += 6;
			out
				<< prefix << iModalIndex + 1 << "->" << iModalIndex + 3 << ": "
					"StructNode(" << pInterfaceNodes[iStrNodem1]->GetLabel() << ") reaction force derivatives [FPx,FPy,FPz]" << std::endl
				<< prefix << iModalIndex + 4 << "->" << iModalIndex + 6 << ": "
					"StructNode(" << pInterfaceNodes[iStrNodem1]->GetLabel() << ") reaction couple derivatives [mPx,mPy,mPz]" << std::endl;
		}
	}

	return out;
}

DofOrder::Order
Modal::GetDofType(unsigned int i) const
{
	ASSERT(i < iGetNumDof());

	if (i < 2*NModes) {
		/* gradi di liberta' differenziali (eq. modali) */
		return DofOrder::DIFFERENTIAL;

	} /* else if (i >= 2*NModes && i < iGetNumDof()) { */
	/* gradi di liberta' algebrici (eq. di vincolo) */
	return DofOrder::ALGEBRAIC;
	/* } */
}

DofOrder::Order
Modal::GetEqType(unsigned int i) const
{
	ASSERT(i < iGetNumDof());

	if (i < 2*NModes) {
		/* gradi di liberta' differenziali (eq. modali) */
		return DofOrder::DIFFERENTIAL;

	} /* else if (i >= 2*NModes && i < iGetNumDof()) { */
	/* gradi di liberta' algebrici (eq. di vincolo) */
	return DofOrder::ALGEBRAIC;
	/* } */
}

void
Modal::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	/* la matrice e' gestita come piena (c'e' un po' di spreco...) */
	*piNumCols = *piNumRows = iRigidOffset + iGetNumDof() + 6*NStrNodes;
}

VariableSubMatrixHandler&
Modal::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering Modal::AssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* gli indici sono ordinati cosi': i primi 6 sono le equazioni
	 * per abbassare di grado il sistema,
	 * quindi le 6 equazioni del moto rigido, quindi le 2*M modali,
	 * quindi le eq. vincolari   */

	/* FIXME: by now, I add a test everywhere it's needed;
	 * later, I'll try to group conditional parts in separate tests */

	/* indici della parte rigida */

	if (pModalNode) {
		integer iRigidIndex = pModalNode->iGetFirstIndex();

		for (unsigned int iCnt = 1; iCnt <= iRigidOffset; iCnt++) {
			WM.PutRowIndex(iCnt, iRigidIndex+iCnt);
			WM.PutColIndex(iCnt, iRigidIndex+iCnt);
		}
	}

	/* indici della parte deformabile e delle reazioni vincolari */
	integer iFlexIndex = iGetFirstIndex();
	unsigned int iNumDof = iGetNumDof();

	for (unsigned int iCnt = 1; iCnt <= iNumDof; iCnt++) {
		WM.PutRowIndex(iRigidOffset+iCnt, iFlexIndex+iCnt);
		WM.PutColIndex(iRigidOffset+iCnt, iFlexIndex+iCnt);
	}

	/* indici delle equazioni vincolari (solo per il nodo 2) */
	for (unsigned int iStrNodem1 = 0; iStrNodem1 < NStrNodes; iStrNodem1++) {
		integer iNodeFirstMomIndex =
			pInterfaceNodes[iStrNodem1]->iGetFirstMomentumIndex();
		integer iNodeFirstPosIndex =
			pInterfaceNodes[iStrNodem1]->iGetFirstPositionIndex();

		integer iOffset = iRigidOffset+iNumDof+6*iStrNodem1;
		for (integer iCnt = 1; iCnt <= 6; iCnt++) {
			WM.PutRowIndex(iOffset+iCnt, iNodeFirstMomIndex+iCnt);
			WM.PutColIndex(iOffset+iCnt, iNodeFirstPosIndex+iCnt);
		}
	}

	/* Assemblaggio dello Jacobiano */

	Vec3 wr(Zero3);
	Mat3x3 wrWedge(Zero3x3);
	Mat3x3 J(Zero3x3);
	Vec3 S(Zero3);
	if (pModalNode) {
		/* recupera e aggiorna i dati necessari */
		/* i prodotti Inv3j*aj ecc. sono gia' stati fatti da AssRes() */

		wr = pModalNode->GetWRef();
		wrWedge = Mat3x3(wr);
#if 0	/* updated by AssRes() */
		R = pModalNode->GetRRef();
		RT = R.Transpose();
#endif
		J = R*(Inv7+Inv8jaj.Symm2())*RT;
		S = R*(Inv2+Inv3jaj);

		/* matrice di massa:        J[1,1] = Mtot  */
		for (int iCnt = 6+1; iCnt <= 6+3; iCnt++) {
			WM.PutCoef(iCnt, iCnt, dMass);
		}

		/* momenti statici J[1,2] = -[S/\]+c[-2w/\S/\+S/\w/\] */
		Mat3x3 SWedge(S);
		WM.Add(6+1, 9+1, ((SWedge*wrWedge-wrWedge*(SWedge*(2.)))*dCoef)-SWedge);

		/* J[2,1] = [S/\] */
		WM.Add(9+1, 6+1, SWedge);

		/* momenti d'inerzia:       J[2,2] = J+c[-(Jw)/\+(w/\)J];    */
		WM.Add(9+1, 9+1, J+((wrWedge*J-Mat3x3(J*wr))*dCoef));

		/* completa lo Jacobiano con l'aggiunta delle equazioni {xP} = {v}
		 {gP} - [wr/\]{g} = {w} */
		for (int iCnt = 1; iCnt <= 6; iCnt++) {
			WM.IncCoef(iCnt, iCnt, 1.);
			WM.DecCoef(iCnt, 6+iCnt, dCoef);
		}
		WM.Sub(3+1, 3+1, wrWedge*dCoef);
	}

	/* parte deformabile :
	 *
	 * | I  -cI  ||aP|
	 * |         ||  |
	 * |cK  M+cC ||bP|
	 */

	for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
		unsigned int iiCnt = iRigidOffset+iCnt;

		WM.PutCoef(iiCnt, iiCnt, 1.);
		WM.PutCoef(iiCnt, iiCnt+NModes, -dCoef);
		for (unsigned int jCnt = 1; jCnt <= NModes; jCnt++) {
			unsigned int jjCnt = iRigidOffset+jCnt;

			WM.PutCoef(iiCnt+NModes, jjCnt,
					dCoef*pModalStiff->dGet(iCnt, jCnt));
			WM.PutCoef(iiCnt+NModes, jjCnt+NModes,
					pModalMass->dGet(iCnt, jCnt)
					+dCoef*pModalDamp->dGet(iCnt, jCnt));
		}
	}

	if (pModalNode) {

		/* termini di accoppiamento moto rigido-deformabile;
		 * eventualmente l'utente potra' scegliere 
		 * se trascurarli tutti, una parte o considerarli tutti */

		/* linearizzazione delle OmegaPrime:
		 * J13 = R*Inv3jaj
		 * J23 = R*Inv4+Inv5jaj
		 * (questi termini ci vogliono sempre)
		 */

		Mat3xN Jac13(NModes, 0.), Jac23(NModes, 0.),
			Inv5jajRef(NModes, 0.);
		MatNx3 Jac13T(NModes, 0.), Jac23T(NModes, 0.);

		Jac13.LeftMult(R, *pInv3);
		Jac23.LeftMult(R, *pInv4);
		Jac23 += Inv5jajRef.LeftMult(R, Inv5jaj);

		WM.Add(6+1, iRigidOffset+NModes+1, Jac13);
		WM.Add(9+1, iRigidOffset+NModes+1, Jac23);
		WM.Add(iRigidOffset+NModes+1, 6+1, Jac13T.Transpose(Jac13));
		WM.Add(iRigidOffset+NModes+1, 9+1, Jac23T.Transpose(Jac23));

		/* termini di Coriolis: linearizzazione delle Omega
		 * (si puo' evitare se non ho grosse vel. angolari):
		 * J13 = -2*R*[Inv3jaPj/\]*RT
		 * J23 = 2*R*[Inv8jaPj-Inv9jkajaPk]*RT */

		Mat3x3 Inv3jaPjWedge(Inv3jaPj);
		WM.Sub(6+1, 9+1, R*Inv3jaPjWedge*(RT*(2.*dCoef)));
#ifdef MODAL_USE_INV9
		WM.Add(9+1, 9+1, R*((Inv8jaPj-Inv9jkajaPk)*(RT*(2.*dCoef))));
#else /* !MODAL_USE_INV9 */
		WM.Add(9+1, 9+1, R*(Inv8jaPj*(RT*(2.*dCoef))));
#endif /* !MODAL_USE_INV9 */

		/* termini di Coriolis: linearizzazione delle b;
		 * si puo' evitare 'quasi' sempre: le velocita' 
		 * di deformazione dovrebbero essere sempre piccole
		 * rispetto ai termini rigidi
		 * Jac13 = 2*[Omega/\]*R*PHI
		 */
#if 0
		Jac13.LeftMult(wrWedge*R*2*dCoef, *pInv3);
		WM.Add(6+1, iRigidOffset+NModes+1, Jac13);
#endif
		/* nota: non riesco a tirar fuori la Omega dall'equazione
		 * dei momenti:
		 * 2*[ri/\]*[Omega/\]*R*PHIi*{DeltaaP},  i = 1,...nnodi
		 * quindi questa equazione non si puo' linearizzare
		 * a meno di ripetere la sommatoria estesa a tutti i nodi
		 * a ogni passo del Newton-Rapson... */

		/* linearizzazione delle forze centrifughe e di Coriolis
		 * in base modale (anche questi termini dovrebbero essere
		 * trascurabili) */
#if 0
		for (unsigned int iMode = 1; iMode <= NModes; iMode++) {
			Inv4j = pInv4->GetVec(iMode);
			VInv5jaj = Inv5jaj.GetVec(iMode);
			VInv5jaPj = Inv5jaPj.GetVec(iMode);
			unsigned int jOffset = (iMode-1)*3+1;
			Inv8jTranspose = (pInv8->GetMat3x3(jOffset)).Transpose();
#ifdef MODAL_USE_INV9
			Inv9jkak = 0.;
			for (unsigned int kModem1 = 0; kModem1 < NModes; kModem1++)  {
		 		doublereal a_kMode = a.dGet(kModem1+1);
		 		integer iOffset = (iMode-1)*3*NModes+kModem1*3+1;
				Inv9jkak += pInv9->GetMat3x3ScalarMult(iOffset, a_kMode);
			}
#endif /* !MODAL_USE_INV9 */
			for (int iCnt = 1; iCnt <= 3; iCnt++) {
				doublereal temp1 = 0., temp2 = 0.;
				for (int jCnt = 1; jCnt<=3; jCnt++) {
#ifdef MODAL_USE_INV9
		 			temp1 += -2*wr.dGet(jCnt)*(R*(Inv8jTranspose-Inv9jkak)*RT).dGet(iCnt, jCnt);
#else /* !MODAL_USE_INV9 */
		 			temp1 += -2*wr.dGet(jCnt)*(R*(Inv8jTranspose)*RT).dGet(iCnt, jCnt);
#endif /* !MODAL_USE_INV9 */
		 			temp2 += -(R*(Inv4j+VInv5jaj)).dGet(jCnt)*wrWedge.dGet(iCnt, jCnt);
				}
				WM.IncCoef(iRigidOffset+NModes+iMode, 9+iCnt,
						dCoef*(temp1+temp2));
			}
		 	for (int iCnt = 1; iCnt<=3; iCnt++) {
				doublereal temp1 = 0.;
				for (int jCnt = 1; jCnt<=3; jCnt++) {
		 			temp1 += (R*VInv5jaPj*2).dGet(jCnt);
				}
				WM.IncCoef(iRigidOffset+NModes+iMode, 9+iCnt,
						dCoef*temp1);
			}
		}
#endif
	}

	/* ciclo esteso a tutti i nodi d'interfaccia */
	for (unsigned int iStrNode = 1; iStrNode <= NStrNodes; iStrNode++) {
		unsigned int iStrNodem1 = iStrNode - 1;

		/* recupero le forme modali del nodo vincolato */
		Mat3xN PHIt(NModes), PHIr(NModes);
		for (unsigned int iMode = 1; iMode <= NModes; iMode++) {
			integer iOffset = (iMode-1)*NStrNodes+iStrNode;

			PHIt.PutVec(iMode, pPHIt->GetVec(iOffset));
			PHIr.PutVec(iMode, pPHIr->GetVec(iOffset));
		}

		MatNx3 PHItT(NModes), PHIrT(NModes);
		PHItT.Transpose(PHIt);
		PHIrT.Transpose(PHIr);

		/* nota: le operazioni
		 * d1tot = d1+PHIt*a, R1tot = R*[I+(PHIr*a)/\]
		 * sono gia' state fatte da AssRes */

		Mat3xN SubMat1(NModes), SubMat2(NModes);
		MatNx3 SubMat3(NModes);
		MatNxN SubMat4(NModes);

		/* cerniera sferica */
		/* F e' aggiornata da AssRes */

		integer iReactionIndex = iRigidOffset+2*NModes+6*iStrNodem1;
		integer iStrNodeIndex = iRigidOffset+iNumDof+6*iStrNodem1;

		Mat3x3 FTmpWedge(pF[iStrNodem1]*dCoef);

		Mat3xN PHItR(NModes);
		PHItR.LeftMult(R, PHIt);

		for (int iCnt = 1; iCnt <= 3; iCnt++) {
			/* termini di reazione sui nodi */
			WM.DecCoef(iStrNodeIndex+iCnt, iReactionIndex+iCnt, 1.);

			/* termini di vincolo dovuti ai nodi */
			WM.IncCoef(iReactionIndex+iCnt, iStrNodeIndex+iCnt, 1.);
		}

		if (pModalNode) {
			for (int iCnt = 1; iCnt <= 3; iCnt++) {
				/* termini di reazione sui nodi */
				WM.IncCoef(6+iCnt, iReactionIndex+iCnt, 1.);
		
				/* termini di vincolo dovuti ai nodi */
				WM.DecCoef(iReactionIndex+iCnt, iCnt, 1.);
			}

			/* pd1Tot e' il puntatore all'array
			 * che contiene le posizioni del nodi FEM */
			Mat3x3 dTmp1Wedge(R*pd1tot[iStrNodem1]);

			WM.Add(9+1, iReactionIndex+1, dTmp1Wedge);
			
			/* termini del tipo: c*F/\*d/\*Deltag */
			WM.Add(9+1, 3+1, FTmpWedge*dTmp1Wedge);
			
			/* termini di vincolo dovuti ai nodi */
			WM.Add(iReactionIndex+1, 3+1, dTmp1Wedge);
			
			/* termini aggiuntivi dovuti alla deformabilita' */
			SubMat3.RightMult(PHItT, RT*FTmpWedge);
			WM.Add(iRigidOffset+NModes+1, 3+1, SubMat3);
			
			SubMat1.LeftMult(FTmpWedge, PHItR);
			WM.Sub(9+1, iRigidOffset+1, SubMat1);
		}
			
		/* contributo dovuto alla flessibilita' */
		WM.Sub(iReactionIndex+1, iRigidOffset+1, PHItR);

		/* idem per pd2 e pR2 */
		Mat3x3 dTmp2Wedge(pR2[iStrNodem1]*pd2[iStrNodem1]);
		WM.Sub(iStrNodeIndex+3+1, iReactionIndex+1, dTmp2Wedge);

		/* termini del tipo: c*F/\*d/\*Deltag */
		WM.Sub(iStrNodeIndex+3+1, iStrNodeIndex+3+1,
				FTmpWedge*dTmp2Wedge);

		/* termini aggiuntivi dovuti alla deformabilita' */
		SubMat3.RightMult(PHItT, RT);
		WM.Add(iRigidOffset+NModes+1, iReactionIndex+1, SubMat3);

		/* modifica: divido le equazioni di vincolo per dCoef */

		/* termini di vincolo dovuti al nodo 2 */
		WM.Sub(iReactionIndex+1, iStrNodeIndex+3+1, dTmp2Wedge);

		/* fine cerniera sferica */

		/* equazioni di vincolo : giunto prismatico */

		/* Vec3 M(XCurr, iModalIndex+2*NModes+6*iStrNodem1+3+1); */
		Vec3 MTmp = pM[iStrNodem1]*dCoef;

		Mat3x3 R1totTranspose = pR1tot[iStrNodem1].Transpose();

		Vec3 e1tota(pR1tot[iStrNodem1].GetVec(1));
		Vec3 e2tota(pR1tot[iStrNodem1].GetVec(2));
		Vec3 e3tota(pR1tot[iStrNodem1].GetVec(3));
		Vec3 e1b(pR2[iStrNodem1].GetVec(1));
		Vec3 e2b(pR2[iStrNodem1].GetVec(2));
		Vec3 e3b(pR2[iStrNodem1].GetVec(3));

		Mat3x3 MWedge(Mat3x3(e3b, e2tota*MTmp.dGet(1))
				+ Mat3x3(e1b, e3tota*MTmp.dGet(2))
				+ Mat3x3(e2b, e1tota*MTmp.dGet(3)));
		Mat3x3 MWedgeT(MWedge.Transpose());

		/* Eq. dei momenti (termini del tipo [e3b/\][e2a/\]M) */
		if (pModalNode) {
			WM.Add(9+1, 3+1, MWedge);
			WM.Sub(9+1, iStrNodeIndex+3+1, MWedgeT);

			WM.Sub(iStrNodeIndex+3+1, 3+1, MWedge);
		}
		WM.Add(iStrNodeIndex+3+1, iStrNodeIndex+3+1, MWedgeT);

		/* Eq. dei momenti, contributo modale */
		Vec3 e1ra(R.GetVec(1));
		Vec3 e2ra(R.GetVec(2));
		Vec3 e3ra(R.GetVec(3));
		Mat3x3 M1Wedge(Mat3x3(e3b, e2ra*MTmp.dGet(1))
				+ Mat3x3(e1b, e3ra*MTmp.dGet(2))
				+ Mat3x3(e2b, e1ra*MTmp.dGet(3)));

		SubMat1.LeftMult(M1Wedge, PHIr);
		if (pModalNode) {
			WM.Add(9+1, iRigidOffset+1, SubMat1);
		}
		WM.Sub(iStrNodeIndex+3+1, iRigidOffset+1, SubMat1);

		/* Eq. d'equilibrio ai modi */
		if (pModalNode) {
			SubMat3.RightMult(PHIrT, R1totTranspose*MWedge);
			WM.Add(iRigidOffset+NModes+1, 3+1, SubMat3);
		}
		SubMat3.RightMult(PHIrT, R1totTranspose*MWedgeT);
		WM.Sub(iRigidOffset+NModes+1, iStrNodeIndex+1, SubMat3);

		Vec3 MA(Mat3x3(e2tota.Cross(e3b), e3tota.Cross(e1b),
				e1tota.Cross(e2b))*(pM[iStrNodem1]*dCoef));
		Mat3x3 MAWedge(MA);
		if (pModalNode) {
			SubMat3.RightMult(PHIrT, R1totTranspose*MAWedge);
			WM.Add(iRigidOffset+NModes+1, 3+1, SubMat3);
		}

		Mat3x3 R1TMAWedge(RT*MAWedge);
		SubMat1.LeftMult(R1TMAWedge, PHIr);
		SubMat2.LeftMult(R1totTranspose*M1Wedge, PHIr);
		SubMat1 += SubMat2;
		SubMat4.Mult(PHIrT, SubMat1);
		for (unsigned int iMode = 1; iMode <= NModes; iMode++) {
			for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
				WM.IncCoef(iRigidOffset+NModes+iMode,
						iRigidOffset+jMode,
						SubMat4.dGet(iMode, jMode));
			}
		}

		/* Termine e2a/\e3b + e3a/\e1b + e1a/\e2b */
		Vec3 v1(e2tota.Cross(e3b));
		Vec3 v2(e3tota.Cross(e1b));
		Vec3 v3(e1tota.Cross(e2b));

		MWedge = Mat3x3(v1, v2, v3);

		if (pModalNode) {
			WM.Add(9+1, iReactionIndex+3+1, MWedge);
		}
		WM.Sub(iStrNodeIndex+3+1, iReactionIndex+3+1, MWedge);

		/* contributo modale:
		 * PHIrT*RT*(e2a/\e3b + e3a/\e1b + e1a/\e2b) */
		SubMat3.RightMult(PHIrT, R1totTranspose*MWedge);
		WM.Add(iRigidOffset+NModes+1, iReactionIndex+3+1, SubMat3);

		/* Modifica: divido le equazioni di vincolo per dCoef */
		MWedge = MWedge.Transpose();

		/* Eq. di vincolo */
		if (pModalNode) {
			WM.Add(iReactionIndex+3+1, 3+1, MWedge);
		}
		WM.Sub(iReactionIndex+3+1, iStrNodeIndex+3+1, MWedge);

		/* Eq. di vincolo, termine aggiuntivo modale */
		Vec3 u1(e2ra.Cross(e3b));
		Vec3 u2(e3ra.Cross(e1b));
		Vec3 u3(e1ra.Cross(e2b));

		M1Wedge = (Mat3x3(u1, u2, u3)).Transpose();
		SubMat1.LeftMult(M1Wedge, PHIr);
		WM.Add(iReactionIndex+3+1, iRigidOffset+1, SubMat1);
	}

	return WorkMat;
}

SubVectorHandler&
Modal::AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering Modal::AssRes()" << std::endl);

	integer iNumRows;
	integer iNumCols;

	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/*
	 * Equations:
	 *
	 * 1  		-> 12:		rigid body eq.
	 *
	 * 13 		-> 12 + 2*NM:	modes eq. (\dot{a}=b); m\dot{b}+ka=f)
	 *
	 * 12 + 2*NM 	-> 12 + 2*NM + 6*NN:		node constraints
	 *
	 * 12 + 2*NM + 6*NN	-> 12 + 2*NM + 6*NN + 6*NN: constr. nodes eq.
	 *
	 *
	 * Unknowns:
	 *
	 * rigid body:			from modal node
	 * modes:			iGetFirstIndex()
	 * node reactions:		iGetFirstIndex() + 2*NM
	 * nodes:			from constraint nodes
	 */

	/* modal dofs and node constraints indices */
	integer iModalIndex = iGetFirstIndex();
	unsigned int iNumDof = iGetNumDof();
	for (unsigned int iCnt = 1; iCnt <= iNumDof; iCnt++) {
		WorkVec.PutRowIndex(iRigidOffset+iCnt, iModalIndex+iCnt);
	}

	/* interface nodes equilibrium indices */
	integer iOffset = iRigidOffset+iNumDof;
	for (unsigned iStrNodem1 = 0; iStrNodem1 < NStrNodes; iStrNodem1++) {
		integer iNodeFirstMomIndex =
			pInterfaceNodes[iStrNodem1]->iGetFirstMomentumIndex();

		for (unsigned int iCnt = 1; iCnt <= 6; iCnt++) {
			WorkVec.PutRowIndex(iOffset+iCnt,
					iNodeFirstMomIndex+iCnt);
		}

		iOffset += 6;
	}

	a.Copy(XCurr, iModalIndex + 1);
	aPrime.Copy(XPrimeCurr, iModalIndex + 1);
	b.Copy(XCurr, iModalIndex + NModes + 1);
	bPrime.Copy(XPrimeCurr, iModalIndex + NModes + 1);

	VecN Ka(NModes), CaP(NModes), MaPP(NModes);

	/*
	 * aggiornamento invarianti
	 */
	Ka.Mult(*pModalStiff, a);
	CaP.Mult(*pModalDamp, b);
	MaPP.Mult(*pModalMass, bPrime);

#if 0
	std::cerr << "### Stiff" << std::endl;
	for (unsigned int i = 1; i <= NModes; i++) {
		for (unsigned int j = 1; j <= NModes; j++) {
			std::cerr << std::setw(2) << i << ", "
				<< std::setw(2) << j << " "
				<< std::setw(20) << (*pModalStiff)(i,j)
				<< std::endl;

		}
	}
	std::cerr << "### Damp" << std::endl;
	for (unsigned int i = 1; i <= NModes; i++) {
		for (unsigned int j = 1; j <= NModes; j++) {
			std::cerr << std::setw(2) << i << ", "
				<< std::setw(2) << j << " "
				<< std::setw(20) << (*pModalDamp)(i,j)
				<< std::endl;
		}
	}
	std::cerr << "### Mass" << std::endl;
	for (unsigned int i = 1; i <= NModes; i++) {
		for (unsigned int j = 1; j <= NModes; j++) {
			std::cerr << std::setw(2) << i << ", "
				<< std::setw(2) << j << " "
				<< std::setw(20) << (*pModalMass)(i,j)
				<< std::endl;
		}
	}
	std::cerr << "### a" << std::endl;
	for (unsigned int i = 1; i <= NModes; i++) {
		std::cerr << std::setw(2) << i <<  " "
			<< std::setw(20) << a(i) << std::endl;
	}
	std::cerr << "### b" << std::endl;
	for (unsigned int i = 1; i <= NModes; i++) {
		std::cerr << std::setw(2) << i <<  " "
			<< std::setw(20) << b(i) << std::endl;
	}
	std::cerr << "### bP" << std::endl;
	for (unsigned int i = 1; i <= NModes; i++) {
		std::cerr << std::setw(2) << i <<  " "
			<< std::setw(20) << bPrime(i) << std::endl;
	}
#endif

	Inv3jaj = *pInv3 * a;
	Inv3jaPj = *pInv3 * b;
	Vec3 Inv3jaPPj = *pInv3 * bPrime;

	/* invarianti rotazionali */
	Vec3 Inv11jaPj;
	Inv11jaPj = R*(*pInv11 * b);

	Inv8jaj = 0.;
	Inv8jaPj = 0.;
	Inv5jaj.Reset(0.);

	Inv5jaPj.Reset(0.);
	Mat3x3 MatTmp1(Zero3x3), MatTmp2(Zero3x3);
#ifdef MODAL_USE_INV9
	Mat3x3 Inv9jkajak(Zero3x3);
	Inv9jkajaPk = 0.;
#endif /* MODAL_USE_INV9 */
	Mat3x3 Inv10jaPj(Zero3x3);

	for (unsigned int iMode = 1; iMode <= NModes; iMode++)  {
		Mat3x3 Inv8jajTmp;
		Mat3x3 Inv10jaPjTmp;

		doublereal a_iMode = a.dGet(iMode);
		doublereal aP_iMode = b.dGet(iMode);

		for (unsigned int jCnt = 1; jCnt <= 3; jCnt++) {
			unsigned int iMjC = (iMode-1)*3+jCnt;

			Inv8jajTmp.PutVec(jCnt, pInv8->GetVec(iMjC));
			Inv10jaPjTmp.PutVec(jCnt, pInv10->GetVec(iMjC)*aP_iMode);
		}

		for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
			Vec3 v = pInv5->GetVec((iMode-1)*NModes+jMode);

			Inv5jaj.AddVec(jMode, v*a_iMode);
			Inv5jaPj.AddVec(jMode, v*aP_iMode);
		}
				
		Inv8jaj += Inv8jajTmp * a_iMode;
		Inv8jaPj += Inv8jajTmp * aP_iMode;
		Inv10jaPj += Inv10jaPjTmp;

#ifdef MODAL_USE_INV9
		/*
		 * questi termini si possono commentare perche' sono
		 * (sempre ?) piccoli (termini del tipo a*a o a*b)
		 * eventualmente dare all'utente la possibilita'
		 * di scegliere se trascurarli o no
		 */
		for (unsigned int kMode = 1; kMode <= NModes; kMode++) {
			doublereal a_kMode = a.dGet(kMode);
			doublereal aP_kMode = b.dGet(kMode);
			unsigned int iOffset = (iMode-1)*3*NModes+(kMode-1)*3+1;
			Inv9jkajak += pInv9->GetMat3x3ScalarMult(iOffset, a_iMode*a_kMode );
			Inv9jkajaPk += pInv9->GetMat3x3ScalarMult(iOffset, a_iMode*aP_kMode);
		}
#endif /* MODAL_USE_INV9 */
	} /* fine ciclo sui modi */

#ifdef MODAL_USE_GRAVITY
	/* forza di gravita' (decidere come inserire g) */
	/* FIXME: use a reasonable reference point where compute gravity */
	Vec3 GravityAcceleration(Zero3);
	bool bGravity = GravityOwner::fGetAcceleration(x, GravityAcceleration);
#endif /* MODAL_USE_GRAVITY */

	Vec3 vP(Zero3);
	Vec3 wP(Zero3);
	Vec3 w(Zero3);
	Vec3 RTw(Zero3);
	if (pModalNode) {
		/* rigid body indices */
		integer iRigidIndex = pModalNode->iGetFirstIndex();
		for (unsigned int iCnt = 1; iCnt <= iRigidOffset; iCnt++) {
			WorkVec.PutRowIndex(iCnt, iRigidIndex+iCnt);
		}

		/* recupera i dati necessari */
		x = pModalNode->GetXCurr();
		Vec3 xP = pModalNode->GetVCurr();
		Vec3 g = pModalNode->GetgCurr();
		Vec3 gP = pModalNode->GetgPCurr();
		vP = pModalNode->GetXPPCurr();
		Vec3 wr = pModalNode->GetWRef();
		wP = pModalNode->GetWPCurr();
		
		Vec3 v(XCurr, iRigidIndex+6+1);
		w = Vec3(XCurr, iRigidIndex+9+1);

		R = pModalNode->GetRCurr();
		RT = R.Transpose();
		RTw = RT*w;

		Mat3x3 J = R*(Inv7+Inv8jaj.Symm2()
#ifdef MODAL_USE_INV9
				-Inv9jkajak
#endif /* MODAL_USE_INV9 */
				)*RT;
		Vec3 S = R*(Inv2+Inv3jaj);

		Mat3xN Inv4Curr(NModes, 0);
		Mat3xN Inv5jajCurr(NModes, 0);
		Inv4Curr.LeftMult(R, *pInv4);
		Inv5jajCurr.LeftMult(R, Inv5jaj);

		/*
		 * fine aggiornamento invarianti
		 */

		/* forze d'inerzia */
		WorkVec.Sub(6+1, vP*dMass-S.Cross(wP)
				+R*Inv3jaPPj+w.Cross(w.Cross(S))
				+(w.Cross(R*Inv3jaPj))*2.);

#if 0
		std::cerr << "m=" << dMass << "; S=" << S
			<< "; a="  << a << "; aPrime =" << aPrime
			<< "; b=" << b <<  "; bPrime= " << bPrime
			<< "; tot=" << vP*dMass-S.Cross(wP)+w.Cross(w.Cross(S))
			+(w.Cross(R*Inv3jaPj))*2+R*Inv3jaPPj << std::endl;
#endif

		WorkVec.Sub(9+1, S.Cross(vP)+J*wP+w.Cross(J*w)
#ifdef MODAL_USE_INV9
				+(R*((Inv8jaPj-Inv9jkajaPk)*RTw))*2.
#else /* !MODAL_USE_INV9 */
				+(R*(Inv8jaPj*RTw))*2.
#endif /* !MODAL_USE_INV9 */
				+Inv4Curr*bPrime+Inv5jajCurr*bPrime);

		/* termini dovuti alle inerzie rotazionali */
		WorkVec.Sub(9+1, R*(Inv10jaPj.Symm2()*RTw)+w.Cross(Inv11jaPj));

#ifdef MODAL_USE_GRAVITY
		/* forza di gravita' (decidere come inserire g) */
		if (bGravity) {
			WorkVec.Add(6+1, GravityAcceleration*dMass);
			WorkVec.Add(9+1, S.Cross(GravityAcceleration));
		}
#endif /* MODAL_USE_GRAVITY */
	
		/* equazioni per abbassare di grado il sistema */
		WorkVec.Add(1, v-xP);
		WorkVec.Add(3+1, w-Mat3x3(MatG, g)*gP-Mat3x3(MatR, g)*wr);
	}

	/* forze modali */
	Mat3x3 Inv10j;
	for (unsigned int iMode = 1; iMode <= NModes; iMode++) {
		Vec3 Inv3j = pInv3->GetVec(iMode);

		Inv4j = pInv4->GetVec(iMode);
		VInv5jaj = Inv5jaj.GetVec(iMode);
		VInv5jaPj = Inv5jaPj.GetVec(iMode);

		unsigned int jOffset = (iMode-1)*3+1;
		Inv8jTranspose = (pInv8->GetMat3x3(jOffset)).Transpose();
		Inv10j = pInv10->GetMat3x3(jOffset);

#ifdef MODAL_USE_INV9
		Inv9jkak = 0.;

		for (unsigned int kModem1 = 0; kModem1 < NModes; kModem1++) {
			doublereal a_kMode = a.dGet(kModem1+1);
			integer iOffset = (iMode-1)*3*NModes+kModem1*3+1;

			Inv9jkak += pInv9->GetMat3x3ScalarMult(iOffset, a_kMode);
		}
#endif /* !MODAL_USE_INV9 */

		WorkVec.IncCoef(iRigidOffset+NModes+iMode,
				-(R*Inv3j).Dot(vP)-(R*(Inv4j+VInv5jaj)).Dot(wP)
#ifdef MODAL_USE_INV9
				+w.Dot(R*((Inv8jTranspose-Inv9jkak+Inv10j)*RTw))
#else /* !MODAL_USE_INV9 */
				+w.Dot(R*((Inv8jTranspose+Inv10j)*RTw))
#endif /* !MODAL_USE_INV9 */
				-(R*VInv5jaPj).Dot(w)*2.
				-MaPP.dGet(iMode)-CaP.dGet(iMode)-Ka.dGet(iMode));

#if 0
		std::cerr << "(R*Inv3j).Dot(vP)=" << (R*Inv3j).Dot(vP)
			<< std::endl
			<< "(R*(Inv4j+VInv5jaj)).Dot(wP)="
			<< (R*(Inv4j+VInv5jaj)).Dot(wP) << std::endl
			<< "w.Dot(R*((Inv8jTranspose+Inv10j)*RTw))"
			<< w.Dot(R*((Inv8jTranspose+Inv10j)*RTw)) << std::endl
			<< " -(R*VInv5jaPj).Dot(w)*2."
			<< -(R*VInv5jaPj).Dot(w)*2.<< std::endl
			<< " -CaP.dGet(iMode)" << -CaP.dGet(iMode)<< std::endl
			<< "-Ka.dGet(iMode) "
			<< -CaP.dGet(iMode)-Ka.dGet(iMode) << std::endl;
#endif
#ifdef MODAL_USE_GRAVITY
		/* forza di gravita': */
		if (bGravity) {
			WorkVec.IncCoef(iRigidOffset+NModes+iMode,
					(R*Inv3j).Dot(GravityAcceleration));
		}
#endif /* MODAL_USE_GRAVITY */
	}

	/* equazioni per abbassare di grado il sistema */
	for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
		WorkVec.PutCoef(iRigidOffset+iCnt, b.dGet(iCnt)-aPrime.dGet(iCnt));
	}

	/* equazioni di vincolo */
	for (unsigned int iStrNode = 1; iStrNode <= NStrNodes; iStrNode++) {
		unsigned int iStrNodem1 = iStrNode - 1;
		integer iReactionIndex = iRigidOffset+2*NModes+6*iStrNodem1;
		integer iStrNodeIndex = iReactionIndex + 6*NStrNodes;

		/* nodo 1 (FEM) e 2 (Multi - Body): recupera la posizione */
		Vec3 d1rig(pOffsetFEMNodes->GetVec(iStrNode));
		pd2[iStrNodem1] = pOffsetMBNodes->GetVec(iStrNode);

		/* recupero le forme modali del nodo vincolato */
		/* FIXME: what about using Blitz++ ? :) */
		Mat3xN PHIt(NModes), PHIr(NModes);
		for (unsigned int iMode = 1; iMode <= NModes; iMode++) {
			integer iOffset = (iMode-1)*NStrNodes+iStrNode;
			PHIt.PutVec(iMode, pPHIt->GetVec(iOffset));
			PHIr.PutVec(iMode, pPHIr->GetVec(iOffset));
		}

		/*
		 * aggiorno d1 e R con il contributo dovuto
		 * alla flessibilita':
		 * d1tot = d1+PHIt*a, R1tot = R*[I+(PHIr*a)/\]
		 */
		pd1tot[iStrNodem1] = d1rig + PHIt*a;
		pR1tot[iStrNodem1] = R*Mat3x3(1., PHIr*a);

		/* constraint reaction (force) */
		pF[iStrNodem1] = Vec3(XCurr,
				iModalIndex+2*NModes+6*iStrNodem1+1);
		Vec3 x2 = pInterfaceNodes[iStrNodem1]->GetXCurr();

		/* FIXME: R2??? */
		pR2[iStrNodem1] = pInterfaceNodes[iStrNodem1]->GetRCurr();

		/* cerniera sferica */
		Vec3 dTmp1(R*pd1tot[iStrNodem1]);
		Vec3 dTmp2(pR2[iStrNodem1]*pd2[iStrNodem1]);

		/* Eq. d'equilibrio, nodo 1 */
		if (pModalNode) {
			WorkVec.Sub(6+1, pF[iStrNodem1]);
			WorkVec.Sub(9+1, dTmp1.Cross(pF[iStrNodem1]));
		}

		/* termine aggiuntivo dovuto alla deformabilita':
		 * -PHItiT*RT*F */
		Vec3 vtemp = RT*pF[iStrNodem1];
		for (unsigned int iMode = 1; iMode <= NModes; iMode++) {
			/* PHItTranspose(iMode) * vtemp */
			WorkVec.DecCoef(iRigidOffset+NModes+iMode,
					PHIt.GetVec(iMode)*vtemp);
		}

		/* Eq. d'equilibrio, nodo 2 */
		WorkVec.Add(iStrNodeIndex+1, pF[iStrNodem1]);
		WorkVec.Add(iStrNodeIndex+3+1, dTmp2.Cross(pF[iStrNodem1]));

		/* Eq. di vincolo */
		if (dCoef != 0.) {
			WorkVec.Add(iReactionIndex+1, (x+dTmp1-x2-dTmp2)/dCoef);
		}

		/* constraint reaction (moment) */
		pM[iStrNodem1] = Vec3(XCurr,
				iModalIndex+2*NModes+6*iStrNodem1+3+1);

		/* giunto prismatico */
		Vec3 e1a(pR1tot[iStrNodem1].GetVec(1));
		Vec3 e2a(pR1tot[iStrNodem1].GetVec(2));
		Vec3 e3a(pR1tot[iStrNodem1].GetVec(3));
		Vec3 e1b(pR2[iStrNodem1].GetVec(1));
		Vec3 e2b(pR2[iStrNodem1].GetVec(2));
		Vec3 e3b(pR2[iStrNodem1].GetVec(3));

		Vec3 MTmp(Mat3x3(e2a.Cross(e3b), e3a.Cross(e1b),
					e1a.Cross(e2b))*pM[iStrNodem1]);

		/* Equazioni di equilibrio, nodo 1 */
		if (pModalNode) {
			WorkVec.Sub(9+1, MTmp);
		}

		/* Equazioni di equilibrio, nodo 2 */
		WorkVec.Add(iStrNodeIndex+3+1, MTmp);

		/* Contributo dovuto alla flessibilita' :-PHIrT*RtotT*M */
		vtemp = pR1tot[iStrNodem1].Transpose()*MTmp;
		for (unsigned int iMode = 1; iMode <= NModes; iMode++) {
			/* PHIrTranspose(iMode) * vtemp */
			WorkVec.DecCoef(iRigidOffset+NModes+iMode,
					PHIr.GetVec(iMode)*vtemp);
		}

		/* Modifica: divido le equazioni di vincolo per dCoef */
		if (dCoef != 0.) {
			/* Equazioni di vincolo di rotazione */
			WorkVec.DecCoef(iReactionIndex+3+1,
					e3b.Dot(e2a)/dCoef);
			WorkVec.DecCoef(iReactionIndex+3+2,
					e1b.Dot(e3a)/dCoef);
			WorkVec.DecCoef(iReactionIndex+3+3,
					e2b.Dot(e1a)/dCoef);
		}
	}

#if 0
	std::cerr << "###" << std::endl;
	for (int i = 1; i <= WorkVec.iGetSize(); i++) {
		std::cerr << std::setw(2) << i << "("
			<< std::setw(2) << WorkVec.iGetRowIndex(i) << ")"
			<< std::setw(20) << WorkVec.dGetCoef(i) << std::endl;
	}
#endif

	return WorkVec;
}

void
Modal::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		/* stampa sul file di output i modi */

		for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
			fOutFlex << std::setw(8) << iCnt
				<< " " << a.dGet(iCnt)
				<< " " << b.dGet(iCnt)
				<< " " << bPrime.dGet(iCnt) << std::endl;
		}
	}
}

unsigned int
Modal::iGetInitialNumDof(void) const
{
	return 2*NModes + 12*NStrNodes;
}

void
Modal::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = *piNumCols = iRigidOffset + iGetInitialNumDof() + 12*NStrNodes;
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
Modal::InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering Modal::InitialAssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	if (pModalNode) {
		integer iRigidIndex = pModalNode->iGetFirstIndex();

		for (unsigned int iCnt = 1; iCnt <= iRigidOffset; iCnt++) {
			WM.PutRowIndex(iCnt, iRigidIndex+iCnt);
			WM.PutColIndex(iCnt, iRigidIndex+iCnt);
		}
	}

	integer iFlexIndex = iGetFirstIndex();

	for (unsigned int iCnt = 1; iCnt <= iGetInitialNumDof(); iCnt++) {
		WM.PutRowIndex(iRigidOffset+iCnt, iFlexIndex+iCnt);
		WM.PutColIndex(iRigidOffset+iCnt, iFlexIndex+iCnt);
	}

	for (unsigned iStrNodem1 = 0; iStrNodem1 < NStrNodes; iStrNodem1++) {
		integer iNodeFirstPosIndex = 
			pInterfaceNodes[iStrNodem1]->iGetFirstPositionIndex();
		integer iNodeFirstVelIndex = iNodeFirstPosIndex+6;

		integer iOffset = iRigidOffset + iGetInitialNumDof() + 12*iStrNodem1;
		for (unsigned int iCnt = 1; iCnt <= 6; iCnt++) {
			WM.PutRowIndex(iOffset+iCnt,
					iNodeFirstPosIndex+iCnt);
			WM.PutColIndex(iOffset+iCnt,
					iNodeFirstPosIndex+iCnt);
			WM.PutRowIndex(iOffset+6+iCnt,
					iNodeFirstVelIndex+iCnt);
			WM.PutColIndex(iOffset+6+iCnt,
					iNodeFirstVelIndex+iCnt);
		}
	}

	/* comincia l'assemblaggio dello Jacobiano */

	/* nota: nelle prime iRigidOffset + 2*M equazioni 
	 * metto dei termini piccoli per evitare singolarita'
	 * quando non ho vincoli esterni o ho azzerato i modi */
	for (unsigned int iCnt = 1; iCnt <= iRigidOffset+2*NModes; iCnt++) {
		WM.PutCoef(iCnt, iCnt, 1.e-12);
	}

	/* forze elastiche e viscose */
	for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
		for (unsigned int jCnt = 1; jCnt <= NModes; jCnt++) {
			WM.PutCoef(iRigidOffset+iCnt, iRigidOffset+jCnt,
					pModalStiff->dGet(iCnt,jCnt));
			WM.PutCoef(iRigidOffset+iCnt, iRigidOffset+NModes+jCnt,
					pModalDamp->dGet(iCnt, jCnt));
		}
	}

	/* equazioni di vincolo */

	/* cerniera sferica */
	for (unsigned int iStrNode = 1; iStrNode <= NStrNodes; iStrNode++) {
		unsigned int iStrNodem1 = iStrNode - 1;

		/* recupera i dati */
		Mat3x3 R2(pInterfaceNodes[iStrNodem1]->GetRRef());
		Vec3 Omega1(Zero3);
		if (pModalNode) {
			Omega1 = pModalNode->GetWRef();
		}
		Vec3 Omega2(pInterfaceNodes[iStrNodem1]->GetWRef());
		Vec3 F(XCurr, iFlexIndex+2*NModes+12*iStrNodem1+1);
		Vec3 FPrime(XCurr, iFlexIndex+2*NModes+12*iStrNodem1+6+1);

		Mat3xN PHIt(NModes,0), PHIr(NModes, 0);
		MatNx3 PHItT(NModes, 0), PHIrT(NModes, 0);

		for (unsigned int iMode = 1; iMode <= NModes; iMode++) {
			PHIt.PutVec(iMode, pPHIt->GetVec((iMode-1)*NStrNodes+iStrNode));
			PHIr.PutVec(iMode, pPHIr->GetVec((iMode-1)*NStrNodes+iStrNode));
		}

		PHItT.Transpose(PHIt);
		PHIrT.Transpose(PHIr);

		Vec3 d1rig(pOffsetFEMNodes->GetVec(iStrNode));
		Vec3 d2(pOffsetMBNodes->GetVec(iStrNode));

		Vec3 d1tot = d1rig+PHIt*a;
		Mat3x3 R1tot = R*Mat3x3(1., PHIr*a);
		Mat3xN SubMat1(NModes, 0.);
		MatNx3 SubMat2(NModes, 0.);

		integer iOffset1 = iRigidOffset+2*NModes+12*iStrNodem1;
		integer iOffset2 = iRigidOffset+iGetInitialNumDof()+12*iStrNodem1;
		if (pModalNode) {
			for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
				/* Contributo di forza all'equazione
				 * della forza, nodo 1 */
				WM.IncCoef(iCnt, iOffset1+iCnt, 1.);

				/* Contrib. di der. di forza all'eq. della der. 
				 * della forza, nodo 1 */
				WM.IncCoef(6+iCnt, iOffset1+6+iCnt, 1.);

				/* Equazione di vincolo, nodo 1 */
				WM.DecCoef(iOffset1+iCnt, iCnt, 1.);

				/* Derivata dell'equazione di vincolo, nodo 1 */
				WM.DecCoef(iOffset1+6+iCnt, 6+iCnt, 1.);
			}
		}

		for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
			/* Equazione di vincolo, nodo 2 */
			WM.IncCoef(iOffset1+iCnt, iOffset2+iCnt, 1.);

			/* Derivata dell'equazione di vincolo, nodo 2 */
			WM.IncCoef(iOffset1+6+iCnt, iOffset2+6+iCnt, 1.);

			/* Contributo di forza all'equazione della forza,
			 * nodo 2 */
			WM.DecCoef(iOffset2+iCnt, iOffset1+iCnt, 1.);

			/* Contrib. di der. di forza all'eq. della der.
			 * della forza, nodo 2 */
			WM.DecCoef(iOffset2+6+iCnt, iOffset1+6+iCnt, 1.);
		}

		/* Distanza nel sistema globale */
		Vec3 d1Tmp(R*d1tot);
		Vec3 d2Tmp(R2*d2);

		/* Matrici F/\d1/\, -F/\d2/\ , F/\omega1/\ */
		Mat3x3 FWedged1Wedge(F, d1Tmp);
		Mat3x3 FWedged2Wedge(F, d2Tmp);
		Mat3x3 FWedgeO1Wedge(F, Omega1);

		/* Matrici (omega1/\d1)/\, -(omega2/\d2)/\ */
		Mat3x3 O1Wedged1Wedge(Omega1.Cross(d1Tmp));
		Mat3x3 O2Wedged2Wedge(d2Tmp.Cross(Omega2));

		/* d1Prime= w1/\d1 + R*PHIt*b */
		Mat3xN R1PHIt(NModes);
		R1PHIt.LeftMult(R, PHIt);
		Vec3 d1Prime(Omega1.Cross(d1Tmp)+R1PHIt*b);

		Mat3x3 FWedge(F);
		if (pModalNode) {
			/* Equazione di momento, nodo 1 */
			WM.Add(3+1, 3+1, FWedged1Wedge);
			WM.Add(3+1, iOffset1+1, Mat3x3(d1Tmp));

			/* Equazione di momento, contributo modale */
			SubMat1.LeftMult(FWedge*R, PHIt);
			WM.Sub(3+1, iRigidOffset+1, SubMat1);
			
			/* Eq. di equilibrio ai modi */
			SubMat2.RightMult(PHItT, RT*FWedge);
			WM.Add(iRigidOffset+1, 3+1, SubMat2);
		}

		/* Equazione di momento, nodo 2 */
		WM.Sub(iOffset2+3+1,iOffset2+3+1, FWedged2Wedge);
		WM.Sub(iOffset2+3+1,iOffset1+1, Mat3x3(d2Tmp));

		/* Eq. di equilibrio ai modi */
		SubMat2.RightMult(PHItT, RT);
		WM.Add(iRigidOffset+1, iRigidOffset+2*NModes+1, SubMat2);

		if (pModalNode) {
			/* derivata dell'equazione di momento, nodo 1 */
			WM.Add(9+1, 3+1,
					(Mat3x3(FPrime)+Mat3x3(F, Omega1))*Mat3x3(d1Tmp)
					+Mat3x3(F, R*(PHIt*b)));
			WM.Add(9+1, 9+1, FWedged1Wedge);
			WM.Add(9+1, iOffset1+1, O1Wedged1Wedge+Mat3x3(R*(PHIt*b)));
			WM.Add(9+1, iOffset1+6+1, Mat3x3(d1Tmp));

			/* derivata dell'equazione di momento, contributo modale */
			SubMat1.LeftMult((-FWedgeO1Wedge-Mat3x3(FPrime))*R, PHIt);
			WM.Add(9+1, iRigidOffset+1, SubMat1);
			SubMat1.LeftMult(-FWedge*R, PHIt);
			WM.Add(9+1, iRigidOffset+NModes+1, SubMat1);

			/* derivata dell'eq. di equilibrio ai modi */
			SubMat2.RightMult(PHItT, RT*FWedge);
			WM.Add(iRigidOffset+NModes+1, 9+1, SubMat2);
			SubMat2.RightMult(PHItT, RT*FWedgeO1Wedge);
			WM.Add(iRigidOffset+NModes+1, 3+1, SubMat2);
		}
		
		SubMat2.RightMult(PHItT, RT);
		WM.Add(iRigidOffset+NModes+1, iRigidOffset+2*NModes+6+1, SubMat2);

		/* Derivata dell'equazione di momento, nodo 2 */
		WM.Add(iOffset2+9+1, iOffset2+3+1,
				(Mat3x3(FPrime) +Mat3x3(F, Omega2))*Mat3x3(-d2Tmp));
		WM.Sub(iOffset2+9+1, iOffset2+9+1, FWedged2Wedge);
		WM.Add(iOffset2+9+1, iOffset1+1, O2Wedged2Wedge);
		WM.Sub(iOffset2+9+1, iOffset1+6+1, Mat3x3(d2Tmp));

		/* Equazione di vincolo */
		if (pModalNode) {
			WM.Add(iOffset1+1, 3+1, Mat3x3(d1Tmp));
		}
		WM.Sub(iOffset1+1, iOffset2+3+1, Mat3x3(d2Tmp));

		/* Equazione di vincolo, contributo modale */
		SubMat1.LeftMult(R, PHIt);
		WM.Sub(iOffset1+1, iRigidOffset+1, SubMat1);

		/* Derivata dell'equazione di vincolo */
		if (pModalNode) {
			WM.Add(iOffset1+6+1, 3+1, O1Wedged1Wedge+R*(PHIt*b));
			WM.Add(iOffset1+6+1, 9+1, Mat3x3(d1Tmp));
		}
		WM.Add(iOffset1+6+1, iOffset2+3+1, O2Wedged2Wedge);
		WM.Sub(iOffset1+6+1, iOffset2+9+1, Mat3x3(d2Tmp));

		/* Derivata dell'equazione di vincolo, contributo modale */
		SubMat1.LeftMult(Mat3x3(Omega1)*R, PHIt);
		WM.Sub(iOffset1+6+1, iRigidOffset+1, SubMat1);
		SubMat1.LeftMult(R, PHIt);
		WM.Sub(iOffset1+6+1, iRigidOffset+NModes+1, SubMat1);

		/* giunto prismatico */
		Vec3 M(XCurr, iFlexIndex+2*NModes+12*iStrNodem1+3+1);
		Vec3 MPrime(XCurr, iFlexIndex+2*NModes+12*iStrNodem1+9+1);
		Vec3 e1tota(R1tot.GetVec(1));
		Vec3 e2tota(R1tot.GetVec(2));
		Vec3 e3tota(R1tot.GetVec(3));
		Vec3 e1b(R2.GetVec(1));
		Vec3 e2b(R2.GetVec(2));
		Vec3 e3b(R2.GetVec(3));

		/* */
		Mat3x3 MWedge(Mat3x3(e3b, e2tota*M.dGet(1))
				+Mat3x3(e1b, e3tota*M.dGet(2))
				+Mat3x3(e2b, e1tota*M.dGet(3)));
		Mat3x3 MWedgeT(MWedge.Transpose());

		/* Equilibrio */
		if (pModalNode) {
			WM.Add(3+1, 3+1, MWedge);
			WM.Sub(3+1, iOffset2+3+1, MWedgeT);

			WM.Add(iOffset2+3+1, 3+1, MWedgeT);
		}
		WM.Sub(iOffset2+3+1, iOffset2+3+1, MWedge);

		/* Equilibrio dei momenti, termini aggiuntivi modali */
		Vec3 e1ra(R.GetVec(1));
		Vec3 e2ra(R.GetVec(2));
		Vec3 e3ra(R.GetVec(3));
		Mat3x3 M1Wedge(Mat3x3(e3b, e2ra*M.dGet(1))
				+Mat3x3(e1b, e3ra*M.dGet(2))
				+Mat3x3(e2b, e1ra*M.dGet(3)));

		SubMat1.LeftMult(M1Wedge, PHIr);
		if (pModalNode) {
			WM.Add(3+1, iRigidOffset+1, SubMat1);
		}
		WM.Sub(iOffset2+3+1, iRigidOffset+1, SubMat1);

		/* Equilibrio ai modi */
		Mat3x3 R1totT(R1tot.Transpose());
		if (pModalNode) {
			SubMat2.RightMult(PHIrT, R1totT*MWedge);
			WM.Add(iRigidOffset+1, 3+1, SubMat2);
		}
		SubMat2.RightMult(PHIrT, R1totT*MWedgeT);
		WM.Sub(iRigidOffset+1, iOffset2+3+1, SubMat2);

		Vec3 MA(Mat3x3(e2tota.Cross(e3b), e3tota.Cross(e1b),
					e1tota.Cross(e2b))*M);
		if (pModalNode) {
			Mat3x3 MAWedge(MA);
			SubMat2.RightMult(PHIrT, R1totT*MAWedge);
			WM.Add(iRigidOffset+1, 3+1, SubMat2);
		}

		Mat3x3 R1TMAWedge(RT*MA);
		SubMat1.LeftMult(R1TMAWedge, PHIr);
		Mat3xN SubMat3(NModes);
		MatNxN SubMat4(NModes);
		SubMat3.LeftMult(R1totT*M1Wedge, PHIr);
		SubMat1 += SubMat3;
		SubMat4.Mult(PHIrT, SubMat1);
		for (unsigned int iMode = 1; iMode <= NModes; iMode++) {
			for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
				WM.IncCoef(iRigidOffset+iMode, iRigidOffset+jMode,
						SubMat4.dGet(iMode,jMode));
			}
		}

		/* Derivate dell'equilibrio */
		if (pModalNode) {
			WM.Add(9+1, 9+1, MWedge);
			WM.Sub(9+1, iOffset2+9+1, MWedgeT);

			WM.Add(iOffset2+9+1, 9+1, MWedgeT);
		}
		WM.Sub(iOffset2+9+1, iOffset2+9+1, MWedge);

		MWedge = ((Mat3x3(e3b, Omega1)+Mat3x3(Omega2.Cross(e3b))*M.dGet(1))
				+Mat3x3(e3b)*MPrime.dGet(1))*Mat3x3(e2tota)
				+((Mat3x3(e1b, Omega1)+Mat3x3(Omega2.Cross(e1b))*M.dGet(2))
				+Mat3x3(e1b)*MPrime.dGet(2))*Mat3x3(e3tota)
				+((Mat3x3(e2b, Omega1)+Mat3x3(Omega2.Cross(e2b))*M.dGet(3))
				+Mat3x3(e2b)*MPrime.dGet(3))*Mat3x3(e1tota);
		
		if (pModalNode) {
			WM.Add(9+1, 3+1, MWedge);
			WM.Sub(iOffset2+9+1, 3+1, MWedge);
		}

		MWedge = ((Mat3x3(e2tota, Omega2)+Mat3x3(Omega1.Cross(e2tota))*M.dGet(1))
				+Mat3x3(e2tota)*MPrime.dGet(1))*Mat3x3(e3b)
				+((Mat3x3(e3tota, Omega2)+Mat3x3(Omega1.Cross(e3tota))*M.dGet(2))
				+Mat3x3(e3tota)*MPrime.dGet(2))*Mat3x3(e1b)
				+((Mat3x3(e1tota, Omega2)+Mat3x3(Omega1.Cross(e1tota))*M.dGet(3))
				+Mat3x3(e1tota)*MPrime.dGet(3))*Mat3x3(e2b);

		if (pModalNode) {
			WM.Sub(9+1, iOffset2+3+1, MWedge);
		}
		WM.Add(iOffset2+9+1, iOffset2+3+1, MWedge);

		/* Derivate dell'equilibrio, termini aggiuntivi modali */
		SubMat1.LeftMult(M1Wedge, PHIr);
		if (pModalNode) {
			WM.Add(9+1, iRigidOffset+NModes+1, SubMat1);
		}
		WM.Sub(iOffset2+9+1, iRigidOffset+NModes+1, SubMat1);

		Vec3 v1(e2tota.Cross(e3b));
		Vec3 v2(e3tota.Cross(e1b));
		Vec3 v3(e1tota.Cross(e2b));

		/* Error handling: il programma si ferma,
		 * pero' segnala dov'e' l'errore */
		if (v1.Dot() < DBL_EPSILON
				|| v2.Dot() < DBL_EPSILON
				|| v3.Dot() < DBL_EPSILON)
		{
			silent_cerr("Modal(" << GetLabel() << "):" 
				<< "warning, first node hinge axis "
				"and second node hinge axis "
				"are (nearly) orthogonal; aborting ..."
				<< std::endl);
			throw Joint::ErrGeneric();
		}
		
		MWedge = Mat3x3(v1, v2, v3);

		/* Equilibrio */
		if (pModalNode) {
			WM.Add(3+1, iOffset1+3+1, MWedge);
		}
		WM.Sub(iOffset2+3+1, iOffset1+3+1, MWedge);

		SubMat2.RightMult(PHIrT, R1totT*MWedge);
		WM.Add(iRigidOffset+1, iOffset1+3+1, SubMat2);

		/* Derivate dell'equilibrio */
		if (pModalNode) {
			WM.Add(9+1, iOffset1+9+1, MWedge);
		}
		WM.Sub(iOffset2+9+1, iOffset1+9+1, MWedge);

		MWedge = MWedge.Transpose();

		/* Equaz. di vincolo */
		if (pModalNode) {
			WM.Add(iOffset1+3+1, 3+1, MWedge);
		}
		WM.Sub(iOffset1+3+1, iOffset2+3+1, MWedge);

		/* Equaz. di vincolo: termine aggiuntivo
		 * dovuto alla flessibilita' */
		Vec3 u1(e2ra.Cross(e3b));
		Vec3 u2(e3ra.Cross(e1b));
		Vec3 u3(e1ra.Cross(e2b));

		M1Wedge = (Mat3x3(u1, u2, u3)).Transpose();
		SubMat1.LeftMult(M1Wedge, PHIr);
		WM.Add(iOffset1+3+1, iRigidOffset+1, SubMat1);

		/* Derivate delle equaz. di vincolo */
		if (pModalNode) {
			WM.Add(iOffset1+9+1, 9+1, MWedge);
		}
		WM.Sub(iOffset1+9+1, iOffset2+9+1, MWedge);

		v1 = e3b.Cross(e2tota.Cross(Omega1))
			+e2tota.Cross(Omega2.Cross(e3b));
		v2 = e1b.Cross(e3tota.Cross(Omega1))
			+e3tota.Cross(Omega2.Cross(e1b));
		v3 = e2b.Cross(e1tota.Cross(Omega1))
			+e1tota.Cross(Omega2.Cross(e2b));

		MWedge = Mat3x3(v1, v2, v3);

		/* Derivate dell'equilibrio */
		if (pModalNode) {
			WM.Add(9+1, iOffset1+3+1, MWedge);
		}
		WM.Sub(iOffset2+9+1, iOffset1+3+1, MWedge);

		/* Derivate delle equaz. di vincolo */
		Omega1 = Omega1-Omega2;

		v1 = e2tota.Cross(e3b.Cross(Omega1));
		Vec3 v1p(e3b.Cross(Omega1.Cross(e2tota)));
		v2 = e3tota.Cross(e1b.Cross(Omega1));
		Vec3 v2p(e1b.Cross(Omega1.Cross(e3tota)));
		v3 = e1tota.Cross(e2b.Cross(Omega1));
		Vec3 v3p(e2b.Cross(Omega1.Cross(e1tota)));

		if (pModalNode) {
			for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
				doublereal d;

				d = v1.dGet(iCnt);
				WM.IncCoef(iOffset1+9+1, 3+iCnt, d);

				d = v2.dGet(iCnt);
				WM.IncCoef(iOffset1+9+2, 3+iCnt, d);

				d = v3.dGet(iCnt);
				WM.IncCoef(iOffset1+9+3, 3+iCnt, d);
			}
		}

		for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
			doublereal d;

			d = v1p.dGet(iCnt);
			WM.IncCoef(iOffset1+9+1, iOffset2+3+iCnt, d);

			d = v2p.dGet(iCnt);
			WM.IncCoef(iOffset1+9+2, iOffset2+3+iCnt, d);

			d = v3p.dGet(iCnt);
			WM.IncCoef(iOffset1+9+3, iOffset2+3+iCnt, d);
		}

	} /* fine ciclo sui nodi vincolati */

	return WorkMat;
}

/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
Modal::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering Modal::InitialAssRes()" << std::endl);

	integer iNumRows;
	integer iNumCols;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	VecN Ka(NModes);

	if (pModalNode) {
		integer iRigidIndex = pModalNode->iGetFirstIndex();

		for (unsigned int iCnt = 1; iCnt <= iRigidOffset; iCnt++) {
			WorkVec.PutRowIndex(iCnt, iRigidIndex+iCnt);
		}
	}

	integer iFlexIndex = iGetFirstIndex();

	for (unsigned int iCnt = 1; iCnt <= iGetInitialNumDof(); iCnt++) {
		WorkVec.PutRowIndex(iRigidOffset+iCnt, iFlexIndex+iCnt);
	}

	for (unsigned iStrNodem1 = 0; iStrNodem1 < NStrNodes; iStrNodem1++) {
		integer iNodeFirstPosIndex = 
			pInterfaceNodes[iStrNodem1]->iGetFirstPositionIndex();
		integer iNodeFirstVelIndex = iNodeFirstPosIndex + 6;

		integer iOffset2 = iRigidOffset+iGetInitialNumDof()+12*iStrNodem1;
		for (unsigned int iCnt = 1; iCnt <= 6; iCnt++) {
			WorkVec.PutRowIndex(iOffset2+iCnt,
					iNodeFirstPosIndex+iCnt);
			WorkVec.PutRowIndex(iOffset2+6+iCnt,
					iNodeFirstVelIndex+iCnt);
		}
	}

	/* aggiorna le forze modali : K*a, C*aP */
	for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
		a.Put(iCnt, XCurr.dGetCoef(iFlexIndex+iCnt));
		b.Put(iCnt, XCurr.dGetCoef(iFlexIndex+NModes+iCnt));
	}

	for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
		doublereal temp = 0.;

		for (unsigned int jCnt = 1; jCnt <= NModes; jCnt++) {
			temp += pModalStiff->dGet(iCnt, jCnt)*a.dGet(jCnt)
				+ pModalDamp->dGet(iCnt, jCnt)*b.dGet(jCnt);
		}
		WorkVec.DecCoef(iRigidOffset+iCnt, temp);
	}

	/* equazioni di vincolo */

	/* Recupera i dati */
	Vec3 v1;
	Vec3 Omega1;

	if (pModalNode) {
		v1 = pModalNode->GetVCurr();
		R = pModalNode->GetRCurr();
		RT = R.Transpose();
		Omega1 = pModalNode->GetWCurr();

	} else {
		v1 = Zero3;
		Omega1 = Zero3;
	}

	for (unsigned int iStrNode = 1; iStrNode <= NStrNodes; iStrNode++) {
		unsigned int iStrNodem1 = iStrNode - 1;
		Mat3xN PHIt(NModes,0), PHIr(NModes,0);

		integer iOffset1 = iRigidOffset+2*NModes+12*iStrNodem1;
		integer iOffset2 = iRigidOffset+iGetInitialNumDof()+12*iStrNodem1;

		for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
			for (unsigned int iMode = 1; iMode <= NModes; iMode++) {
				PHIt.Put(iCnt, iMode, pPHIt->dGet(iCnt, (iMode-1)*NStrNodes+iStrNode));
				PHIr.Put(iCnt, iMode, pPHIr->dGet(iCnt, (iMode-1)*NStrNodes+iStrNode));
			}
		}

		MatNx3 PHItT(NModes), PHIrT(NModes);
		PHItT.Transpose(PHIt);
		PHIrT.Transpose(PHIr);

		Vec3 d1rig(pOffsetFEMNodes->GetVec(iStrNode));
		Vec3 d2(pOffsetMBNodes->GetVec(iStrNode));

		Vec3 d1tot = d1rig+PHIt*a;
		Mat3x3 R1tot = R*Mat3x3(1., PHIr*a);
		Mat3x3 R1totT(R1tot.Transpose());

		Vec3 x2(pInterfaceNodes[iStrNodem1]->GetXCurr());
		Vec3 v2(pInterfaceNodes[iStrNodem1]->GetVCurr());
		Mat3x3 R2(pInterfaceNodes[iStrNodem1]->GetRCurr());
		Vec3 Omega2(pInterfaceNodes[iStrNodem1]->GetWCurr());
		Vec3 F(XCurr, iFlexIndex+2*NModes+12*iStrNodem1+1);
		Vec3 FPrime(XCurr, iFlexIndex+2*NModes+12*iStrNodem1+6+1);

		/* cerniera sferica */

		/* Distanza nel sistema globale */
		Vec3 d1Tmp(R*d1tot);
		Vec3 d2Tmp(R2*d2);

		/* Vettori omega1/\d1, -omega2/\d2 */
		Vec3 O1Wedged1(Omega1.Cross(d1Tmp));
		Vec3 O2Wedged2(Omega2.Cross(d2Tmp));

		/* d1Prime= w1/\d1 + R*PHIt*b */
		Mat3xN R1PHIt(NModes);
		R1PHIt.LeftMult(R, PHIt);
		Vec3 d1Prime(O1Wedged1+R1PHIt*b);

		/* Equazioni di equilibrio, nodo 1 */
		if (pModalNode) {
			WorkVec.Sub(1, F);
			WorkVec.Add(3+1, F.Cross(d1Tmp)); /* F/\d = -d/\F */
		}

		/* Eq. d'equilibrio ai modi */
		for (unsigned int iMode = 1; iMode <= NModes; iMode++) {
			doublereal temp = 0.;
			for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
				temp += PHIt.dGet(iCnt, iMode)*(RT*F).dGet(iCnt);
			}
			WorkVec.DecCoef(iRigidOffset+iMode, temp);
		}

		/* Derivate delle equazioni di equilibrio, nodo 1 */
		if (pModalNode) {
			WorkVec.Sub(6+1, FPrime);
			WorkVec.Add(9+1, FPrime.Cross(d1Tmp)-d1Prime.Cross(F));
		}

		/* Derivata dell'eq. di equilibrio ai modi */
		MatNx3 PHItTR1T(NModes);
		PHItTR1T.RightMult(PHItT, RT);
		for (unsigned int iMode = 1; iMode <= NModes; iMode++) {
			doublereal temp = 0.;

			for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
				temp += PHItTR1T.dGet(iMode,iCnt)*(Omega1.Cross(F)-FPrime).dGet(iCnt);
			}
			WorkVec.IncCoef(iRigidOffset+NModes+iMode, temp);
		}

		/* Equazioni di equilibrio, nodo 2 */
		WorkVec.Add(iOffset2+1, F);
		WorkVec.Add(iOffset2+3+1, d2Tmp.Cross(F));

		/* Derivate delle equazioni di equilibrio, nodo 2 */
		WorkVec.Add(iOffset2+6+1, FPrime);
		WorkVec.Add(iOffset2+9+1, d2Tmp.Cross(FPrime)+O2Wedged2.Cross(F));

		/* Equazione di vincolo */
		WorkVec.Add(iOffset1+1, x+d1Tmp-x2-d2Tmp);

		/* Derivata dell'equazione di vincolo */
		WorkVec.Add(iOffset1+6+1, v1+d1Prime-v2-O2Wedged2);

		/* giunto prismatico */
		Vec3 M(XCurr, iFlexIndex+2*NModes+12*iStrNodem1+3+1);
		Vec3 MPrime(XCurr, iFlexIndex+2*NModes+12*iStrNodem1+9+1);

		Vec3 e1a(R1tot.GetVec(1));
		Vec3 e2a(R1tot.GetVec(2));
		Vec3 e3a(R1tot.GetVec(3));
		Vec3 e1b(R2.GetVec(1));
		Vec3 e2b(R2.GetVec(2));
		Vec3 e3b(R2.GetVec(3));

		Vec3 MTmp(e2a.Cross(e3b*M.dGet(1))
				+e3a.Cross(e1b*M.dGet(2))
				+e1a.Cross(e2b*M.dGet(3)));

		/* Equazioni di equilibrio, nodo 1 */
		if (pModalNode) {
			WorkVec.Sub(3+1, MTmp);
		}

		/* Equazioni di equilibrio, nodo 2 */
		WorkVec.Add(iOffset2+3+1, MTmp);

		/* Equazioni di equilibrio, contributo modale */
		for (unsigned int iMode = 1; iMode <= NModes; iMode++) {
			doublereal temp = 0;
			
			for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
				temp += PHIrT.dGet(iMode, iCnt)*(R1totT*MTmp).dGet(iCnt);
			}
			WorkVec.DecCoef(iRigidOffset+iMode, temp);
		}

		/* eaPrime = w/\ea + R*[(PHIr*b)/\]ia */
		Mat3x3 Tmp = R*Mat3x3(PHIr*b);
		Vec3 e1aPrime = Omega1.Cross(e2a)+Tmp.GetVec(1);
		Vec3 e2aPrime = Omega1.Cross(e2a)+Tmp.GetVec(2);
		Vec3 e3aPrime = Omega1.Cross(e2a)+Tmp.GetVec(3);
		Vec3 MTmpPrime(0.);
		MTmpPrime =
			 (e2a.Cross(Omega2.Cross(e3b))-e3b.Cross(e2aPrime))*M.dGet(1)
			+(e3a.Cross(Omega2.Cross(e1b))-e1b.Cross(e3aPrime))*M.dGet(1)
			+(e1a.Cross(Omega2.Cross(e2b))-e2b.Cross(e1aPrime))*M.dGet(1)
			+e2a.Cross(e3b*MPrime.dGet(1))
			+e3a.Cross(e1b*MPrime.dGet(2))
			+e1a.Cross(e2b*MPrime.dGet(3));

		/* Derivate delle equazioni di equilibrio, nodo 1 */
		if (pModalNode) {
			WorkVec.Sub(9+1, MTmpPrime);
		}

		/* Derivate delle equazioni di equilibrio, nodo 2 */
		WorkVec.Add(iOffset2+9+1, MTmpPrime);

		/* Derivate delle equazioni di equilibrio, contributo modale */
		MatNx3 SubMat1(NModes);
		SubMat1.RightMult(PHIrT, R1totT);

		/* FIXME: temporary ((PHIr*b).Cross(RT*MTmp)) */
		Vec3 T1 = MTmpPrime-Omega1.Cross(MTmp);
		Vec3 T2 = (PHIr*b).Cross(RT*MTmp);

		for (unsigned int iMode = 1; iMode <= NModes; iMode++) {
			doublereal temp = 0;
			
			for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
				temp += SubMat1.dGet(iMode, iCnt)*T1.dGet(iCnt)
					-PHIrT.dGet(iMode, iCnt)*T2.dGet(iCnt);
			}
			WorkVec.DecCoef(iRigidOffset+NModes+iMode, temp);
		}

		/* Equazioni di vincolo di rotazione */
		WorkVec.DecCoef(iOffset1+3+1, e3b.Dot(e2a));
		WorkVec.DecCoef(iOffset1+3+2, e1b.Dot(e3a));
		WorkVec.DecCoef(iOffset1+3+3, e2b.Dot(e1a));

		/* Derivate delle equazioni di vincolo di rotazione */
		WorkVec.IncCoef(iOffset1+9+1,
				e3b.Dot(Omega2.Cross(e2a)-e2aPrime));
		WorkVec.IncCoef(iOffset1+9+2,
				e1b.Dot(Omega2.Cross(e3a)-e3aPrime));
		WorkVec.IncCoef(iOffset1+9+3,
				e2b.Dot(Omega2.Cross(e1a)-e1aPrime));

	} /* fine equazioni di vincolo */

	return WorkVec;
}

void
Modal::SetValue(VectorHandler& X, VectorHandler& XP) const
{
	/* inizializza la soluzione e la sua derivata
	 * subito dopo l'assemblaggio iniziale
	 * e prima del calcolo delle derivate */

	int iFlexIndex = iGetFirstIndex();

	for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
		/* modal multipliers */
		X.Put(iFlexIndex+iCnt, a.dGet(iCnt));

		/* derivatives of modal multipliers */
		X.Put(iFlexIndex+NModes+iCnt, b.dGet(iCnt));
		XP.Put(iFlexIndex+iCnt, b.dGet(iCnt));
	}
}

unsigned int
Modal::iGetNumPrivData(void) const
{
	return 3*NModes;
}

unsigned int
Modal::iGetPrivDataIdx(const char *s) const
{
	/*
	 * x[n] | xP[n] | xPP[n]
	 *
	 * dove n e' il numero del modo (a base 1)
	 */

	unsigned int idx = 0;

	/* che cosa e' richiesto? */
	if (strncmp(s, "x[", sizeof("x[") - 1) == 0) {
		s += sizeof("x[") - 1;

	} else if (strncmp(s, "xP[", sizeof("xP[") - 1) == 0) {
		s += sizeof("xP[") - 1;
		idx += NModes;

	} else if (strncmp(s, "xPP[", sizeof("xPP[") - 1) == 0) {
		s += sizeof("xPP[") - 1;
		idx += 2*NModes;

	} else {
		return 0;
	}

	/* trova la parentesi chiusa (e il terminatore di stringa) */
	char *end = strchr(s, ']');
	if (end == NULL || end[1] != '\0') {
		return 0;
	}

	/* buffer per numero (dimensione massima: 32 bit) */
	char buf[] = "18446744073709551615UL";
	size_t		len = end - s;

	ASSERT(len < sizeof(buf));

	strncpy(buf, s, len);
	buf[len] = '\0';

	/* leggi il numero */
#ifdef HAVE_STRTOUL
	unsigned long n = strtoul(buf, &end, 10);
	if (end == buf || end[0] != '\0') {
		return 0;
	}
#else /* !HAVE_STRTOUL */
	long n = atol(buf);
	if (n <= 0) {
		return 0;
	}
#endif /* !HAVE_STRTOUL */

	if (n <= 0 || n > NModes) {
		return 0;
	}

	return idx + n;
}


doublereal
Modal::dGetPrivData(unsigned int i) const
{
	ASSERT(i > 0 && i < iGetNumPrivData());
	if (i <= NModes) {
		return a.dGet(i);

	} else {
		return b.dGet(i - NModes);
	}
}

Mat3xN *
Modal::GetCurrFemNodesPosition(void)
{
	if (pCurrXYZ == NULL) {
		   SAFENEWWITHCONSTRUCTOR(pCurrXYZ, Mat3xN, Mat3xN(NFemNodes, 0.));
	}

	if (pModalNode) {
		R = pModalNode->GetRCurr();
		x = pModalNode->GetXCurr();
	}

	for (unsigned int iMode = 1; iMode <= NModes; iMode++) {
		for (unsigned int iNode = 1; iNode <= NFemNodes; iNode++) {
      			for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
				pCurrXYZ->Add(iCnt, iNode,
					pXYZFemNodes->dGet(iCnt, iNode) + pModeShapest->dGet(iCnt, (iMode-1)*NFemNodes+iNode) * a.dGet(iMode) );
			}
		}
	}
	
	pCurrXYZ->LeftMult(R);
	for (unsigned int iNode = 1; iNode <= NFemNodes; iNode++) {
      		for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
	 		pCurrXYZ->Add(iCnt, iNode, x.dGet(iCnt));
      		}
   	}

      	return pCurrXYZ;
}

Mat3xN *
Modal::GetCurrFemNodesVelocity(void)
{
	if (pCurrXYZ == NULL) {
		   SAFENEWWITHCONSTRUCTOR(pCurrXYZVel, Mat3xN, Mat3xN(NFemNodes, 0.));
	}

	Vec3 w(Zero3);
	Vec3 v(Zero3);

	if (pModalNode) {
      		R = pModalNode->GetRCurr();
		w = pModalNode->GetWCurr();
		v = pModalNode->GetVCurr();
	}

	Mat3x3 wWedge(w);

      	Mat3xN CurrXYZTmp(NFemNodes, 0.);
	for (unsigned int iMode = 1; iMode <= NModes; iMode++) {
		for (unsigned int iNode = 1; iNode <= NFemNodes; iNode++) {
      			for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
				doublereal d = pModeShapest->dGet(iCnt,(iMode-1)*NFemNodes+iNode);
				pCurrXYZVel->Add(iCnt, iNode, pXYZFemNodes->dGet(iCnt,iNode) + d * a.dGet(iMode) );
				CurrXYZTmp.Add(iCnt, iNode,d * b.dGet(iMode) );
			}
		}
	}
	pCurrXYZVel->LeftMult(wWedge*R);
	CurrXYZTmp.LeftMult(R);
	for (unsigned int iNode = 1; iNode <= NFemNodes; iNode++) {
      		for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
	 		pCurrXYZVel->Add(iCnt, iNode, CurrXYZTmp.dGet(iCnt, iNode) + v.dGet(iCnt));
      		}
   	}

      	return pCurrXYZVel;
}

/* from gravity.h */
/* massa totale */
doublereal
Modal::dGetM(void) const
{
	return dMass;
}

/* momento statico */
Vec3
Modal::GetS_int(void) const
{
	Mat3x3 R = pModalNode->GetRCurr();

	Vec3 S = R*(Inv2+Inv3jaj);

	return S;
}

/* momento d'inerzia */
Mat3x3
Modal::GetJ_int(void) const
{
	Mat3x3 R = pModalNode->GetRCurr();
	Mat3x3 RT = R.Transpose();

	Mat3x3 J = R*(Inv7+Inv8jaj.Symm2()
#ifdef MODAL_USE_INV9
			-Inv9jkajak
#endif /* MODAL_USE_INV9 */
			)*RT;

	return J;
}


Joint *
ReadModal(DataManager* pDM,
		MBDynParser& HP,
		const DofOwner* pDO,
		unsigned int uLabel)
{
	/* legge i dati d'ingresso e li passa al costruttore dell'elemento */
	Joint* pEl = NULL;
	
	ModalNode* pModalNode = 0;
	Vec3 X0(Zero3);
	Mat3x3 R(Eye3);

	/* If the modal element is clamped, a fixed position 
	 * and orientation of the reference point can be added */
	if (HP.IsKeyWord("clamped")) {
		if (HP.IsKeyWord("position")) {
			X0 = HP.GetPosAbs(AbsRefFrame);
		}

		if (HP.IsKeyWord("orientation")) {
			R = HP.GetRotAbs(AbsRefFrame);
		}

	/* otherwise a structural node of type "modal" must be given;
	 * the "origin node" of the FEM model, if given, or the 0,0,0 
	 * coordinate point of the mesh will be put in the modal node */
	} else {
		unsigned int uNode = HP.GetInt();

		DEBUGCOUT("Linked to Modal Node: " << uNode << std::endl);

		/* verifica di esistenza del nodo */
		StructNode* pTmpNode = pDM->pFindStructNode(uNode);
		if (pTmpNode == NULL) {
			silent_cerr("Modal(" << uLabel << "): "
				"StructuralNode(" << uNode << ") "
				"at line " << HP.GetLineData()
				<< " not defined" << std::endl);
			throw DataManager::ErrGeneric();
		}

		if (pTmpNode->GetStructNodeType() != StructNode::MODAL) {
			silent_cerr("Modal(" << uLabel << "): "
				"illegal type for "
				"StructuralNode(" << uNode << ") "
				"at line " << HP.GetLineData() 
				<< std::endl);
			throw DataManager::ErrGeneric();
		}
		pModalNode = dynamic_cast<ModalNode *>(pTmpNode);

		/* la posizione del nodo modale e' quella dell'origine del SdR
		 * del corpo flessibile */
		X0 = pModalNode->GetXCurr();

		/* orientazione del corpo flessibile, data come orientazione
		 * del nodo modale collegato */
		R = pModalNode->GetRCurr();
	}

	/* Legge i dati relativi al corpo flessibile */
	int tmpNModes = HP.GetInt();     /* numero di modi */
	if (tmpNModes <= 0) {
		silent_cerr("Modal(" << uLabel << "): "
			"illegal number of modes " << tmpNModes << " at line "
			<< HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric();
	}
	unsigned int NModes = (unsigned int)tmpNModes;

	unsigned int *uModeNumber = NULL;
	SAFENEWARR(uModeNumber, unsigned int, NModes);
	if (HP.IsKeyWord("list")) {
		for (unsigned int iCnt = 0; iCnt < NModes; iCnt++) {
			int n = HP.GetInt();

			if (n <= 0) {
				char *th = NULL;

				switch ((iCnt+1)%10) {
				case 1:
					th = "st";
					break;

				case 2:
					th = "nd";
					break;

				case 3:
					th = "rd";
					break;

				default:
					th = "th";
					break;
				}

				silent_cerr("Modal(" << uLabel << "): "
					"illegal " << iCnt+1 << "'" << th 
					<< " mode number " << n 
					<< std::endl);
				throw ErrGeneric();
			}

			/* FIXME: check for duplicates? */
			uModeNumber[iCnt] = n;
		}

	} else {
		for (unsigned int iCnt = 0; iCnt < NModes; iCnt++) {
			uModeNumber[iCnt] = iCnt+1;
		}
	}

	/* numero di nodi FEM del modello */
	int tmpNFemNodes = HP.GetInt();
	if (tmpNFemNodes <= 0) {
		silent_cerr("Modal(" << uLabel << "): "
			"illegal number of FEM nodes " << tmpNFemNodes
			<< " at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric();
	}
	unsigned int NFemNodes = (unsigned int)tmpNFemNodes;

#ifdef MODAL_SCALE_DATA
	/* Legge gli eventuali fattori di scala per le masse nodali
	 * (scale masses) e per i modi (scale modes)   */

	/* NOTA: AL MOMENTO NON E' USATO */
	doublereal scalemasses = 1.;
	doublereal scaleinertia = 1.;
	doublereal scalemodes = 1.;

	if (HP.IsKeyWord("scale" "masses")) {
		scalemasses = HP.GetReal();
	}

	if (HP.IsKeyWord("scale" "modes")) {
		scalemodes = HP.GetReal();
	}

	scaleinertia = scalemasses*(scalemodes*scalemodes);
#endif /* MODAL_SCALE_DATA */

	/* Legge i coefficienti di smorzamento */
	doublereal cdamp = 0.;
	VecN DampRatios(NModes, 0.);
	integer iDampFlag = 0;

	if (HP.IsKeyWord("no" "damping")) {
		cdamp = 0.;
	} else if (HP.IsKeyWord("proportional" "damping")) {
		cdamp = HP.GetReal();
	} else if (HP.IsKeyWord("diag" "damping"))  {
		iDampFlag = 1;
		
		for (unsigned int iCnt = 1; iCnt <= NModes; iCnt ++) {
			integer iDampedMode =  HP.GetInt();
			cdamp = HP.GetReal();
			DampRatios.Put(iDampedMode, cdamp);
		}
	} else {
		silent_cout("Modal(" << uLabel << "): "
				"no damping is assumed at line "
				<< HP.GetLineData() << " (deprecated)"
				<< std::endl);
	}

	DEBUGCOUT("Number of Modes Imported : " << NModes << std::endl);
	DEBUGCOUT("Number of FEM Nodes Imported : " << NFemNodes << std::endl);
	DEBUGCOUT("Origin of FEM Model : " << X0 << std::endl);
	DEBUGCOUT("Orientation of FEM Model : " << R << std::endl);
	/* DEBUGCOUT("Damping coefficient: "<< cdamp << std::endl); */

	doublereal dMass = 0;              /* massa totale */
	Vec3 STmp(0.);                     /* momenti statici  */
	Mat3x3 JTmp(0.);                   /* inerzia totale  */
	VecN FemMass(NFemNodes, 0.);       /* masse nodali   */
	Mat3xN FemJ(NFemNodes, 0.);        /* inerzie nodali (sono diagonali) */

	Mat3xN *pModeShapest = NULL;       /* forme di traslaz. e rotaz. */
	Mat3xN *pModeShapesr = NULL;
	Mat3xN PHIti(NModes, 0.);          /* forme nodo i-esimo: 3*nmodi */
	Mat3xN PHIri(NModes, 0.);
	Mat3xN *pXYZFemNodes = NULL;       /* punt. alle coordinate nodali */
	Mat3xN *pOffsetFEMNodes = NULL;    /* punt. offset FEM (per vincoli) */
	Mat3xN *pOffsetMBNodes = NULL;     /* punt. offset MB (per vincoli) */
	Mat3xN *pRotMBNodes = NULL;        /* punt. orient. MB (per vincoli) */
	MatNxN *pGenMass = NULL;           /* punt. masse e rigidezze modali */
	MatNxN *pGenStiff = NULL;
	MatNxN *pGenDamp = NULL;

	Mat3xN *pInv3 = NULL;     	      /* invarianti d'inerzia */
	Mat3xN *pInv4 = NULL;
	Mat3xN *pInv5 = NULL;
	Mat3xN *pInv8 = NULL;
	Mat3xN *pInv9 = NULL;
	Mat3xN *pInv10 = NULL;
	Mat3xN *pInv11 = NULL;

	VecN *a = NULL;                    /* spostamenti e velocita' modali */
	VecN *aP = NULL;

	unsigned int iNode, iMode, jMode, iStrNode;

	/* input file */
	const char *s = HP.GetFileName();
	if (s == NULL) {
		silent_cerr("Modal(" << uLabel << "): unable to get "
			"modal data file name" << std::endl);
		throw ErrGeneric();
	}

	const char *sFileFem = NULL;
	SAFESTRDUP(sFileFem, s);

	/* apre il file con i dati del modello FEM */
	std::ifstream fdat(sFileFem);
	if (!fdat) {
		silent_cerr("Modal(" << uLabel << "): "
			"unable to open file \"" << sFileFem << "\""
			"at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric();
	}
	DEBUGCOUT("Reading Flexible Body Data from file "
			<< sFileFem << std::endl);

	/* carica i dati relativi a coordinate nodali, massa, momenti statici
	 * e d'inerzia, massa e rigidezza generalizzate dal file nomefile.
	 */

	doublereal d;
	unsigned int NFemNodesDADS = 0, NModesDADS = 0, NRejModes = 0;
	char str[BUFSIZ];

	/* alloca la memoria per le matrici necessarie a memorizzare i dati
	 * relativi al corpo flessibile
	 * nota: devo usare i puntatori perche' altrimenti non si riesce
	 * a passarli al costruttore       */
	SAFENEWWITHCONSTRUCTOR(pXYZFemNodes, Mat3xN, Mat3xN(NFemNodes, 0.));
	SAFENEWWITHCONSTRUCTOR(pGenMass,  MatNxN, MatNxN(NModes, 0.));
	SAFENEWWITHCONSTRUCTOR(pGenStiff, MatNxN, MatNxN(NModes, 0.));
	SAFENEWWITHCONSTRUCTOR(pGenDamp,  MatNxN, MatNxN(NModes, 0.));
	SAFENEWWITHCONSTRUCTOR(pModeShapest, Mat3xN,
			Mat3xN(NFemNodes*NModes, 0.));
	SAFENEWWITHCONSTRUCTOR(pModeShapesr, Mat3xN,
			Mat3xN(NFemNodes*NModes, 0.));
	SAFENEWWITHCONSTRUCTOR(pInv3, Mat3xN, Mat3xN(NModes, 0.));
	SAFENEWWITHCONSTRUCTOR(pInv4, Mat3xN, Mat3xN(NModes, 0.));
	/* Inv5 e' un 3xMxM */
	SAFENEWWITHCONSTRUCTOR(pInv5, Mat3xN, Mat3xN(NModes*NModes, 0.));
	/* Inv8 e' un 3x3xM */
	SAFENEWWITHCONSTRUCTOR(pInv8, Mat3xN, Mat3xN(3*NModes, 0.));
	/* Inv9 e' un 3x3xMxM */
	SAFENEWWITHCONSTRUCTOR(pInv9, Mat3xN, Mat3xN(3*NModes*NModes, 0.));
	/* Inv10 e' un 3x3xM */
	SAFENEWWITHCONSTRUCTOR(pInv10,Mat3xN, Mat3xN(3*NModes, 0.));
	SAFENEWWITHCONSTRUCTOR(pInv11,Mat3xN, Mat3xN(NModes, 0.));
	SAFENEWWITHCONSTRUCTOR(a,  VecN, VecN(NModes, 0.));
	SAFENEWWITHCONSTRUCTOR(aP, VecN, VecN(NModes, 0.));

	/* array contenente le label dei nodi FEM */
	unsigned int *IdFemNodes = NULL;
	SAFENEWARR(IdFemNodes, unsigned int, NFemNodes);

	bool *bActiveModes = NULL;

	while (!fdat.eof()) {        /* parsing del file */
		fdat.getline(str, sizeof(str));

		/* legge il primo blocco (HEADER) */
		if (!strncmp("** RECORD GROUP 1,", str,
					sizeof("** RECORD GROUP 1,") - 1)) {
	 		fdat.getline(str, sizeof(str));
			fdat >> str;

			/* FEM nodes number */
			fdat >> NFemNodesDADS;

			/* add to modes number */
			fdat >> NModesDADS;

			unsigned int i;
			fdat >> i;
			NModesDADS += i;
			fdat >> i;
			NModesDADS += i;

			/* "rejected" modes (subtract from modes number) */
			fdat >> NRejModes;
			NModesDADS -= NRejModes;
		
			/* consistency checks */
			if (NFemNodes != NFemNodesDADS) {
				silent_cerr("Modal(" << uLabel << "), "
					"file \"" << sFileFem << "\": "
					"FEM nodes number " << NFemNodes
					<< " does not match node number "
					<< NFemNodesDADS
					<< std::endl);
				throw DataManager::ErrGeneric();
			}

			if (NModes != NModesDADS) {
				silent_cout("Modal(" << uLabel
						<< "), file '" << sFileFem
						<< "': using " << NModes
						<< " of " << NModesDADS
						<< " modes" << std::endl);
			}

			if (bActiveModes != NULL) {
				throw ErrGeneric();
			}

			SAFENEWARR(bActiveModes, bool, NModesDADS+1);

			for (unsigned int iCnt = 1; iCnt <= NModesDADS; iCnt++) {
				bActiveModes[iCnt] = false;
			}

			for (unsigned int iCnt = 0; iCnt < NModes; iCnt++) {
				if (uModeNumber[iCnt] > NModesDADS) {
					silent_cerr("Modal(" << uLabel << "): "
						"mode " << uModeNumber[iCnt]
						<< " is not available (max = "
						<< NModesDADS << ")"
						<< std::endl);
					throw ErrGeneric();
				}
				bActiveModes[uModeNumber[iCnt]] = true;
			}
			SAFEDELETEARR(uModeNumber);
			uModeNumber = NULL;

		/* legge il secondo blocco (Id.nodi) */
		} else if (!strncmp("** RECORD GROUP 2,", str,
					sizeof("** RECORD GROUP 2,") - 1)) {
			for (iNode = 1; iNode <= NFemNodes; iNode++) {
				fdat >> IdFemNodes[iNode-1];
			}

		/* deformate iniziali dei modi */
		} else if (!strncmp("** RECORD GROUP 3,", str,
					sizeof("** RECORD GROUP 3,") - 1)) {
			unsigned int iCnt = 1;

			if (bActiveModes == NULL) {
				silent_cerr("Modal(" << uLabel << "): "
					"input file \"" << sFileFem << "\""
					"is bogus (RECORD GROUP 3)"
					<< std::endl);
				throw ErrGeneric();
			}

			for (iMode = 1; iMode <= NModesDADS; iMode++) {
				fdat >> d;
				if (!bActiveModes[iMode]) {
					continue;
				}
				a->Put(iCnt, d);
				iCnt++;
			}

		/* velocita' iniziali dei modi */
		} else if (!strncmp("** RECORD GROUP 4,", str,
					sizeof("** RECORD GROUP 4,") - 1)) {
			unsigned int iCnt = 1;

			if (bActiveModes == NULL) {
				silent_cerr("Modal(" << uLabel << "): "
					"input file \"" << sFileFem << "\""
					"is bogus (RECORD GROUP 4)"
					<< std::endl);
				throw ErrGeneric();
			}

			for (iMode = 1; iMode <= NModesDADS; iMode++) {
				fdat >> d;
				if (!bActiveModes[iMode]) {
					continue;
				}
				aP->Put(iCnt, d);
				iCnt++;
			}

		/* Coordinate X dei nodi*/
		} else if (!strncmp("** RECORD GROUP 5,", str,
					sizeof("** RECORD GROUP 5,") - 1)) {
			for (iNode = 1; iNode <= NFemNodes; iNode++) {
				fdat >> d;
#ifdef MODAL_SCALE_DATA
				d *= scalemodes;
#endif /* MODAL_SCALE_DATA */
				pXYZFemNodes->Put(1, iNode, d);
			}

		/* Coordinate Y dei nodi*/
		} else if (!strncmp("** RECORD GROUP 6,", str,
					sizeof("** RECORD GROUP 6,") - 1)) {
			for (iNode = 1; iNode <= NFemNodes; iNode++) {
				fdat >> d;
#ifdef MODAL_SCALE_DATA
				d *= scalemodes;
#endif /* MODAL_SCALE_DATA */
				pXYZFemNodes->Put(2, iNode, d);
			}

		/* Coordinate Z dei nodi*/
		} else if (!strncmp("** RECORD GROUP 7,", str,
					sizeof("** RECORD GROUP 7,") - 1)) {
			for (iNode = 1; iNode <= NFemNodes; iNode++) {
				fdat >> d;
#ifdef MODAL_SCALE_DATA
				d *= scalemodes;
#endif /* MODAL_SCALE_DATA */
				pXYZFemNodes->Put(3, iNode, d);
			}

		/* Forme modali */
		} else if (!strncmp("** RECORD GROUP 8,", str,
					sizeof("** RECORD GROUP 8,") - 1)) {
			for (iMode = 1; iMode <= NRejModes; iMode++) {
				/* FIXME: siamo sicuri di avere 
				 * raggiunto '\n'? */
				fdat.getline(str, sizeof(str));
				fdat.getline(str, sizeof(str));
			}
			
			if (bActiveModes == NULL) {
				silent_cerr("Modal(" << uLabel << "): "
					"input file \"" << sFileFem << "\""
					"is bogus (RECORD GROUP 8)"
					<< std::endl);
				throw ErrGeneric();
			}

			unsigned int iCnt = 1;
			for (iMode = 1; iMode <= NModesDADS; iMode++) {
				fdat.getline(str, sizeof(str));
				for (iNode = 1; iNode <= NFemNodes; iNode++) {
					doublereal t1, t2, t3, r1, r2, r3;

					fdat >> t1 >> t2 >> t3
						>> r1 >> r2 >> r3;

					if (!bActiveModes[iMode]) {
						continue;
					}

#ifdef MODAL_SCALE_DATA
					t1 *= scalemodes;
					t2 *= scalemodes;
					t3 *= scalemodes;
#endif /* MODAL_SCALE_DATA */
					pModeShapest->PutVec((iCnt-1)*NFemNodes+iNode, Vec3(t1, t2, t3));
					pModeShapesr->PutVec((iCnt-1)*NFemNodes+iNode, Vec3(r1, r2, r3));
				}

				if (bActiveModes[iMode]) {
					iCnt++;
				}
				fdat.getline(str, sizeof(str));
			}

		/* Matrice di massa  modale */
		} else if (!strncmp("** RECORD GROUP 9,", str,
					sizeof("** RECORD GROUP 9,") - 1)) {
			unsigned int iCnt = 1;

			if (bActiveModes == NULL) {
				silent_cerr("Modal(" << uLabel << "): "
					"input file \"" << sFileFem << "\""
					"is bogus (RECORD GROUP 9)"
					<< std::endl);
				throw ErrGeneric();
			}

			for (iMode = 1; iMode <= NModesDADS; iMode++) {
				unsigned int jCnt = 1;

				for (jMode = 1; jMode <= NModesDADS; jMode++) {
					fdat >> d;
					if (!bActiveModes[iMode] || !bActiveModes[jMode]) {
						continue;
					}
					pGenMass->Put(iCnt, jCnt, d);
					jCnt++;
				}

				if (bActiveModes[iMode]) {
					iCnt++;
				}
			}

			/* Matrice di rigidezza  modale */
		} else if (!strncmp("** RECORD GROUP 10,", str,
					sizeof("** RECORD GROUP 10,") - 1)) {
			unsigned int iCnt = 1;

			if (bActiveModes == NULL) {
				silent_cerr("Modal(" << uLabel << "): "
					"input file \"" << sFileFem << "\""
					"is bogus (RECORD GROUP 10)"
					<< std::endl);
				throw ErrGeneric();
			}

			for (iMode = 1; iMode <= NModesDADS; iMode++) {
				unsigned int jCnt = 1;

				for (jMode = 1; jMode <= NModesDADS; jMode++) {
					fdat >> d;
					if (!bActiveModes[iMode] || !bActiveModes[jMode]) {
						continue;
					}
					pGenStiff->Put(iCnt, jCnt, d);
					jCnt++;
				}

				if (bActiveModes[iMode]) {
					iCnt++;
				}
			}

			/* Lumped Masses */
		} else if (!strncmp("** RECORD GROUP 11,", str,
					sizeof("** RECORD GROUP 11,") - 1)) {
			for (iNode = 1; iNode <= NFemNodes; iNode++) {
				for (unsigned int jCnt = 1; jCnt <= 6; jCnt++) {
					fdat >> d;
					switch (jCnt) {
					case 1:
#ifdef MODAL_SCALE_DATA
						d *= scalemass;
#endif /* MODAL_SCALE_DATA */
						FemMass.Put(iNode, d);
						break;
					
					case 4:
#ifdef MODAL_SCALE_DATA
						d *= scaleinertia;
#endif /* MODAL_SCALE_DATA */
						FemJ.Put(1, iNode, d);
						break;

					case 5:
#ifdef MODAL_SCALE_DATA
						d *= scaleinertia;
#endif /* MODAL_SCALE_DATA */
						FemJ.Put(2, iNode, d);
						break;
					
					case 6:
#ifdef MODAL_SCALE_DATA
						d *= scaleinertia;
#endif /* MODAL_SCALE_DATA */
						FemJ.Put(3, iNode, d);
						break;
					}
				}
			}
		} /* fine parser del file */
	}
	SAFEDELETEARR(bActiveModes);
	bActiveModes = NULL;

	SAFEDELETEARR(sFileFem);
	sFileFem = NULL;

	fdat.close();

	/* lettura dati di vincolo:
	 * l'utente specifica la label del nodo FEM e del nodo rigido
	 * d'interfaccia.
	 * L'orientamento del nodo FEM e' quello del nodo modale, la
	 * posizione e' la somma di quella modale e di quella FEM   */

	/* puntatori ai nodi multibody */
	const StructNode** pInterfaceNodes = NULL;
	/* array contenente le label dei nodi d'interfaccia */
	unsigned int *IntFEMNodes = NULL;
	unsigned int *IntMBNodes = NULL;
	/* array contenente le forme modali dei nodi d'interfaccia */
	Mat3xN* pPHItStrNode = NULL;
	Mat3xN* pPHIrStrNode = NULL;

	/* traslazione origine delle coordinate */
	if (HP.IsKeyWord("origin" "node")) {
		/* numero di nodi FEM del modello */
		unsigned int NFemOriginNode = HP.GetInt();

		for (iNode = 0; iNode < NFemNodes; iNode++) {
			if (NFemOriginNode == IdFemNodes[iNode]) {
				break;
			}
		}

		if (iNode == NFemNodes) {
			silent_cerr("Modal(" << uLabel << "): "
				"FEM node " << NFemOriginNode
				<< " at line " << HP.GetLineData()
				<< " not defined " << std::endl);
			throw DataManager::ErrGeneric();
		}

		iNode++;
		Vec3 Origin(pXYZFemNodes->GetVec(iNode));

		pedantic_cout("Modal(" << uLabel << "): origin x={" << Origin << "}" << std::endl);

		for (iStrNode = 1; iStrNode <= NFemNodes; iStrNode++) {
			pXYZFemNodes->SubVec(iStrNode, Origin);
		}
	}

	/* numero di nodi d'interfaccia */
	unsigned int NStrNodes = HP.GetInt();
	DEBUGCOUT("Number of Interface Nodes : " << NStrNodes << std::endl);

	SAFENEWWITHCONSTRUCTOR(pOffsetFEMNodes, Mat3xN, Mat3xN(NStrNodes, 0.));
	SAFENEWWITHCONSTRUCTOR(pOffsetMBNodes, Mat3xN, Mat3xN(NStrNodes, 0.));
	SAFENEWWITHCONSTRUCTOR(pRotMBNodes, Mat3xN, Mat3xN(3*NStrNodes, 0.));
	SAFENEWWITHCONSTRUCTOR(pPHItStrNode, Mat3xN,
			Mat3xN(NStrNodes*NModes, 0.));
	SAFENEWWITHCONSTRUCTOR(pPHIrStrNode, Mat3xN,
			Mat3xN(NStrNodes*NModes, 0.));

	SAFENEWARR(pInterfaceNodes, const StructNode*, NStrNodes);
	SAFENEWARR(IntFEMNodes, unsigned int, NStrNodes);
	SAFENEWARR(IntMBNodes, unsigned int, NStrNodes);

	for (iStrNode = 1; iStrNode <= NStrNodes; iStrNode++) {
		/* nodo collegato 1 (è il nodo FEM) */
		unsigned int uNode1 = (unsigned int)HP.GetInt();
		DEBUGCOUT("Linked to FEM Node " << uNode1 << std::endl);

		/* verifica di esistenza del nodo 1*/
		for (iNode = 0; iNode < NFemNodes; iNode++) {
			if (uNode1 == IdFemNodes[iNode]) {
				break;
			}
		}

		if (iNode == NFemNodes) {
			silent_cerr("Modal(" << uLabel << "): "
				"FEM node " << uNode1
				<< " at line " << HP.GetLineData()
				<< " not defined " << std::endl);
			throw DataManager::ErrGeneric();
		}

		iNode++;
		
		int iNodeCurr = iNode;

		/* recupera la posizione del nodo FEM, somma di posizione
		 * e eventuale offset;
		 *
		 * HEADS UP: non piu' offset per i nodi FEM !!!!!!!!!
		 * 
		 * nota: iNodeCurr contiene la posizione a cui si trova
		 * il nodo FEM d'interfaccia nell'array pXYZNodes */
		pOffsetFEMNodes->PutVec(iStrNode,
				pXYZFemNodes->GetVec(iNodeCurr));

		/* salva le forme modali del nodo d'interfaccia
		 * nell'array pPHIStrNode */
		for (iMode = 0; iMode < NModes; iMode++) {
			pPHItStrNode->PutVec(iMode*NStrNodes+iStrNode,
					pModeShapest->GetVec(iMode*NFemNodes+iNodeCurr));
			pPHIrStrNode->PutVec(iMode*NStrNodes+iStrNode,
					pModeShapesr->GetVec(iMode*NFemNodes+iNodeCurr));
		}

		/* nodo collegato 2 (e' il nodo multibody) */
		unsigned int uNode2 = (unsigned int)HP.GetInt();
		DEBUGCOUT("Linked to Multi-Body Node " << uNode2 << std::endl);

		/* verifica di esistenza del nodo 2 */
		pInterfaceNodes[iStrNode-1] = pDM->pFindStructNode(uNode2);
		if (pInterfaceNodes[iStrNode-1] == NULL) {
			silent_cerr("Modal(" << uLabel << "): "
				"StructuralNode(" << uNode2 << ") "
				"at line " << HP.GetLineData()
				<< " not defined" << std::endl);
			throw DataManager::ErrGeneric();
		}

		/* offset del nodo Multi-Body */
		ReferenceFrame RF = ReferenceFrame(pInterfaceNodes[iStrNode-1]);
		Vec3 d2(HP.GetPosRel(RF));
		Mat3x3 R2(Eye3);
		if (HP.IsKeyWord("hinge") || HP.IsKeyWord("orientation")) {
			R2 = HP.GetRotRel(RF);
		}

		pOffsetMBNodes->PutVec(iStrNode, d2);
		pRotMBNodes->PutMat3x3(3*(iStrNode-1)+1, R2);

		DEBUGCOUT("Multibody Node reference frame d2:" << std::endl
				<< d2 << std::endl);

		/* salva le label dei nodi vincolati nell'array IntNodes;
		 * puo' servire per il restart? */
		IntFEMNodes[iStrNode-1] = uNode1;
		IntMBNodes[iStrNode-1] = uNode2;

		pedantic_cout("Interface node " << iStrNode << ":" << std::endl
				<< "        MB node " << uNode2 << " x={" << pInterfaceNodes[iStrNode-1]->GetXCurr() << "}" << std::endl);
		if (pModalNode) {
			Vec3 u = pXYZFemNodes->GetVec(iNodeCurr);
			pedantic_cout("        FEM node " << uNode1 << " x={" << pModalNode->GetXCurr() + pModalNode->GetRCurr()*u
					<< "} xrel={" << u << "}" << std::endl);
		} else {
			pedantic_cout("        FEM node " << uNode1 << " x={" << pXYZFemNodes->GetVec(iNodeCurr) << "}" << std::endl);
		}
	}  /* fine ciclo sui nodi d'interfaccia */

	/* fine ciclo caricamento dati */

	/*
	 * calcola gli invarianti d'inerzia (massa, momenti statici
	 * e d'inerzia, termini di accoppiamento nel SdR locale)
	 */

	/* inizio ciclo scansione nodi */
	for (iNode = 1; iNode <= NFemNodes; iNode++) {
		doublereal mi = FemMass.dGet(iNode);

		/* massa totale (Inv 1) */
		dMass += mi;

		/* posizione nodi FEM */
		Vec3 ui = pXYZFemNodes->GetVec(iNode);

		Mat3x3 uivett(ui);
		Mat3x3 JiNodeTmp(0.);

		JiNodeTmp.Put(1, 1, FemJ.dGet(1, iNode));
		JiNodeTmp.Put(2, 2, FemJ.dGet(2, iNode));
		JiNodeTmp.Put(3, 3, FemJ.dGet(3, iNode));

		JTmp += JiNodeTmp-Mat3x3(ui, ui*mi);
		STmp += ui*mi;

		/* estrae le forme modali del nodo i-esimo */
		for (iMode = 1; iMode <= NModes; iMode++) {
			unsigned int iOffset = (iMode-1)*NFemNodes+iNode;

			PHIti.PutVec(iMode, pModeShapest->GetVec(iOffset));
			PHIri.PutVec(iMode, pModeShapesr->GetVec(iOffset));
		}

		Mat3xN Inv3Tmp(NModes, 0.);
		Mat3xN Inv4Tmp(NModes, 0.);
		Mat3xN Inv4JTmp(NModes, 0.);
		Inv3Tmp.Copy(PHIti);

		/* Inv3 = mi*PHIti,      i = 1,...nnodi */
		Inv3Tmp *= mi;

		/* Inv4 = mi*ui/\*PHIti+Ji*PHIri, i = 1,...nnodi */
		Inv4Tmp.LeftMult(uivett*mi, PHIti);
		Inv4JTmp.LeftMult(JiNodeTmp, PHIri);
		Inv4Tmp += Inv4JTmp;
		*pInv3 += Inv3Tmp;
		*pInv4 += Inv4Tmp;
		*pInv11 += Inv4JTmp;

		/* inizio ciclo scansione modi */
		for (iMode = 1; iMode <= NModes; iMode++) {
			Vec3 PHItij = PHIti.GetVec(iMode);
			Vec3 PHIrij = PHIri.GetVec(iMode);

			Mat3x3 PHItijvett_mi(PHItij*mi);
			Mat3xN Inv5jTmp(NModes, 0);

			/* Inv5 = mi*PHItij/\*PHIti,
			 * i = 1,...nnodi, j = 1,...nmodi */
			Inv5jTmp.LeftMult(PHItijvett_mi, PHIti);
			for (jMode = 1; jMode <= NModes; jMode++)  {
				pInv5->AddVec((iMode-1)*NModes+jMode,
						Inv5jTmp.GetVec(jMode));
			}

			/* Inv8 = -mi*ui/\*PHItij/\,
			 * i = 1,...nnodi, j = 1,...nmodi */
			Mat3x3 Inv8jTmp = -uivett*PHItijvett_mi;
			pInv8->AddMat3x3((iMode-1)*3+1, Inv8jTmp);

			/* Inv9 = mi*PHItij/\*PHItik/\,
			 * i = 1,...nnodi, j, k = 1...nmodi */
			for (unsigned int kMode = 1; kMode <= NModes; kMode++) {
				Mat3x3 PHItikvett(PHIti.GetVec(kMode));
				Mat3x3 Inv9jkTmp = PHItijvett_mi*PHItikvett;

				pInv9->AddMat3x3((iMode-1)*3*NModes+(kMode-1)*3+1, Inv9jkTmp);
			}

			/* Inv10 = [PHIrij/\][J0i],
			 * i = 1,...nnodi, j = 1,...nmodi */
			Mat3x3 Inv10jTmp = Mat3x3(PHIrij)*JiNodeTmp;
			pInv10->AddMat3x3((iMode-1)*3+1, Inv10jTmp);
		} /*  fine ciclo scansione modi */
	} /* fine ciclo scansione nodi */

	/*
	 * costruisce la matrice di smorzamento:
	 * il termine diagonale i-esimo e' pari a
	 * cii = 2*cdampi*(ki*mi)^.5
	 */
	for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
		doublereal d = sqrt(pGenStiff->dGet(iCnt, iCnt)
				*pGenMass->dGet(iCnt, iCnt ));

		if (!iDampFlag) {
			pGenDamp->Put(iCnt, iCnt, 2.*cdamp*d);
		} else {
			pGenDamp->Put(iCnt, iCnt, 2.*DampRatios.dGet(iCnt)*d);
		}
	}

#ifdef DEBUG
	DEBUGCOUT("Total Mass : " << dMass << std::endl);
	DEBUGCOUT("Inertia Matrix : " << std::endl << JTmp << std::endl);
	DEBUGCOUT("Static Moment Vector : " << STmp << std::endl);

	DEBUGCOUT("Generalized Stiffness: " << std::endl);
	for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
		for (unsigned int jCnt = 1; jCnt <= NModes; jCnt++) {
			std::cout << " " << pGenStiff->dGet(iCnt, jCnt);
		}
		std::cout << std::endl;
	}

	DEBUGCOUT("Generalized Mass: " << std::endl);
	for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
		for (unsigned int jCnt = 1; jCnt <= NModes; jCnt++) {
			std::cout << " " << pGenMass->dGet(iCnt,jCnt);
		}
		std::cout << std::endl;
	}

	DEBUGCOUT("Generalized Damping: " << std::endl);
	for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
		for (unsigned int jCnt = 1; jCnt <= NModes; jCnt++) {
			std::cout << " " << pGenDamp->dGet(iCnt,jCnt);
		}
		std::cout << std::endl;
	}

	DEBUGCOUT("Inv3 : " << std::endl);
	for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
		for (unsigned int jCnt = 1; jCnt <= NModes; jCnt++) {
			std::cout << " " << pInv3->dGet(iCnt,jCnt);
		}
		std::cout << std::endl;
	}

	DEBUGCOUT("Inv4 : " << std::endl);
	for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
		for (unsigned int jCnt = 1; jCnt <= NModes; jCnt++) {
			std::cout << " " << pInv4->dGet(iCnt,jCnt);
		}
		std::cout << std::endl;
	}
	
	for (iMode = 1; iMode <= NModes; iMode++) {
		DEBUGCOUT("Inv5j : " << " j = " << iMode << std::endl);
		for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
			for (jMode = 1; jMode <= NModes; jMode++) {
				std::cout << " " << pInv5->dGet(iCnt, (iMode-1)*NModes+jMode);
			}
			std::cout << std::endl;
		}
	}
	
	for (iMode = 1; iMode <= NModes; iMode++) {
		DEBUGCOUT("Inv8j : " << " j = " << iMode << std::endl);
		for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
			for (unsigned int jCnt = 1; jCnt <= 3; jCnt++) {
				std::cout << " " << pInv8->dGet(iCnt, (iMode-1)*3+jCnt);
			}
			std::cout << std::endl;
		}
	}

	for (iMode = 1; iMode <= NModes; iMode++) {
		for (jMode = 1; jMode <= NModes; jMode++) {
			DEBUGCOUT("Inv9jk : " << " j = " << iMode << " k = " << jMode << std::endl);
			for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
				for (unsigned int jCnt = 1; jCnt <= 3; jCnt++) {
					std::cout << " " << pInv9->dGet(iCnt, (iMode-1)*3*NModes+(jMode-1)*3+jCnt);
				}
				std::cout << std::endl;
			}
		}
	}
#endif /* DEBUG */

	const char *sFileMod = HP.GetFileName();
	flag fOut = pDM->fReadOutput(HP, Elem::JOINT);

	SAFENEWWITHCONSTRUCTOR(pEl,
			Modal,
			Modal(uLabel,
				pModalNode,
				X0,
				R,
				pDO,
				NModes,
				NStrNodes,
				NFemNodes,
				dMass,
				STmp,
				JTmp,
				pGenMass,
				pGenStiff,
				pGenDamp,
				IdFemNodes,
				IntFEMNodes,
				IntMBNodes,
				pXYZFemNodes,
				pOffsetFEMNodes,
				pOffsetMBNodes,
				pRotMBNodes,
				pInterfaceNodes,
				pPHItStrNode,
				pPHIrStrNode,
				pModeShapest,
				pModeShapesr,
				pInv3,
				pInv4,
				pInv5,
				pInv8,
				pInv9,
				pInv10,
				pInv11,
				a,
				aP,
				sFileMod,
				pDM,
				HP,
				fOut));

	return pEl;
}

