/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2011
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
 * Copyright 1999-2011 Felice Felippone <ffelipp@tin.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */

/*
 * Copyright 1999-2011 Pierangelo Masarati  <masarati@aero.polimi.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 *
 * Modified by Pierangelo Masarati
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

/* FIXME: gravity in modal elements is eXperimental; undefine to disable */
#define MODAL_USE_GRAVITY

#include <cerrno>
#include <cstring>
#include <sys/stat.h>
#include <limits>
#include <algorithm>

#include "modal.h"
#include "dataman.h"
#include "Rot.hh"
#include "hint_impl.h"

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
	const std::vector<unsigned int>& uModeNumber,
	MatNxN *pGenMass,
	MatNxN *pGenStiff,
	MatNxN *pGenDamp,
	std::vector<std::string>& IdFEMNodes,	/* label nodi FEM */
	std::vector<std::string>& IntFEMNodes,	/* label nodi FEM d'interfaccia */
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
	DataManager* pDM,  /* non serve */
	MBDynParser& HP,   /* non serve */
	flag fOut)
: Elem(uL, fOut),
Joint(uL, pDO, fOut),
pModalNode(pR),
iRigidOffset(pModalNode ? 12 : 0),
x(x0),
R(R0),
RT(R0.Transpose()),
NModes(NM),
NStrNodes(NI),
NFEMNodes(NF),
dMass(dMassTmp),
Inv2(STmp),
Inv7(JTmp),
uModeNumber(uModeNumber),
pModalMass(pGenMass),
pModalStiff(pGenStiff),
pModalDamp(pGenDamp),
IdFEMNodes(IdFEMNodes),
IntFEMNodes(IntFEMNodes),
pXYZFEMNodes(pN),
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
Inv3jaj(Zero3),
Inv3jaPj(Zero3),
Inv8jaj(Zero3x3),
Inv8jaPj(Zero3x3),
Inv5jaj(NModes, 0.),
Inv5jaPj(NModes, 0.),
VInv5jaj(Zero3),
VInv5jaPj(Zero3),
Inv8jTranspose(Eye3),
Inv9jkajak(Eye3),
Inv9jkajaPk(Eye3),
a(*aa),
aPrime(NModes, 0.),
b(*bb),
bPrime(NModes, 0.),
pd1tot(NULL),
pd2(NULL),
pR1tot(NULL),
pR2(NULL),
pF(NULL),
pM(NULL)
{
	ASSERT(pModalNode->GetStructNodeType() == StructNode::MODAL);

	SAFENEWARRNOFILL(pd1tot, Vec3, NStrNodes);
	SAFENEWARRNOFILL(pd2, Vec3, NStrNodes);
	SAFENEWARRNOFILL(pR1tot, Mat3x3, NStrNodes);
	SAFENEWARRNOFILL(pR2, Mat3x3, NStrNodes);
	SAFENEWARRNOFILL(pF, Vec3, NStrNodes);
	SAFENEWARRNOFILL(pM, Vec3, NStrNodes);

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
	if (pInterfaceNodes) {
		SAFEDELETEARR(pInterfaceNodes);
	}
	if (pXYZFEMNodes) {
		SAFEDELETE(pXYZFEMNodes);
	}
	if (pModalMass) {
		SAFEDELETE(pModalMass);
	}
	if (pModalStiff) {
		SAFEDELETE(pModalStiff);
	}
	if (pModalDamp) {
		SAFEDELETE(pModalDamp);
	}
	if (pModeShapest) {
		SAFEDELETE(pModeShapest);
	}
	if (pModeShapesr) {
		SAFEDELETE(pModeShapesr);
	}
	if (pOffsetFEMNodes) {
		SAFEDELETE(pOffsetFEMNodes);
	}
	if (pOffsetMBNodes) {
		SAFEDELETE(pOffsetMBNodes);
	}
	if (pRotMBNodes) {
		SAFEDELETE(pRotMBNodes);
	}
	if (pPHIt) {
		SAFEDELETE(pPHIt);
	}
	if (pPHIr) {
		SAFEDELETE(pPHIr);
	}
	if (pInv3) {
		SAFEDELETE(pInv3);
	}
	if (pInv4) {
		SAFEDELETE(pInv4);
	}
	if (pInv5) {
		SAFEDELETE(pInv5);
	}
	if (pInv8) {
		SAFEDELETE(pInv8);
	}
	if (pInv9) {
		SAFEDELETE(pInv9);
	}
	if (pInv10) {
		SAFEDELETE(pInv10);
	}
	if (pInv11) {
		SAFEDELETE(pInv11);
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
Modal::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
	integer iModalIndex = iGetFirstIndex();

	out 
		<< prefix << iModalIndex + 1 << "->" << iModalIndex + NModes
		<< ": modal displacements [q(1.." << NModes << ")]" << std::endl
		<< prefix << iModalIndex + NModes + 1 << "->" << iModalIndex + 2*NModes
		<< ": modal velocities [qP(1.." << NModes << ")]" << std::endl;
	iModalIndex += 2*NModes;
	for (unsigned iStrNodem1 = 0; iStrNodem1 < NStrNodes; iStrNodem1++, iModalIndex += 6) {
		out
			<< prefix << iModalIndex + 1 << "->" << iModalIndex + 3 << ": "
				"StructNode(" << pInterfaceNodes[iStrNodem1]->GetLabel() << ") "
				"reaction forces [Fx,Fy,Fz]" << std::endl
			<< prefix << iModalIndex + 4 << "->" << iModalIndex + 6 << ": "
				"StructNode(" << pInterfaceNodes[iStrNodem1]->GetLabel() << ") "
				"reaction couples [Mx,My,Mz]" << std::endl;
		if (bInitial) {
			iModalIndex += 6;
			out
				<< prefix << iModalIndex + 1 << "->" << iModalIndex + 3 << ": "
					"StructNode(" << pInterfaceNodes[iStrNodem1]->GetLabel() << ") "
					"reaction force derivatives [FPx,FPy,FPz]" << std::endl
				<< prefix << iModalIndex + 4 << "->" << iModalIndex + 6 << ": "
					"StructNode(" << pInterfaceNodes[iStrNodem1]->GetLabel() << ") "
					"reaction couple derivatives [MPx,MPy,MPz]" << std::endl;
		}
	}

	return out;
}

static const char xyz[] = "xyz";
static const char *mdof[] = {
	"modal displacement q",
	"modal velocity qP"
};
static const char *rdof[] = {
	"reaction force F",
	"reaction couple M",
	"reaction force derivative FP",
	"reaction couple derivative MP"
};
static const char *meq[] = {
	"modal velocity definition q",
	"modal equilibrium equation f"
};
static const char *req[] = {
	"position constraint P",
	"orientation constraint g",
	"velocity constraint v",
	"angular velocity constraint w"
};

void
Modal::DescribeDof(std::vector<std::string>& desc, bool bInitial, int i) const
{
	int iend = 1;
	if (i == -1) {
		iend = 2*NModes + 6*NStrNodes;
		if (bInitial) {
			iend += 6*NStrNodes;
		}
	}
	desc.resize(iend);

	unsigned modulo = 6;
	if (bInitial) {
		modulo = 12;
	}

	std::ostringstream os;
	os << "Modal(" << GetLabel() << ")";

	if (i < -1 ) {
		// error
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);

	} else if (i == -1) {
		std::string name(os.str());

		for (unsigned i = 0; i < 2*NModes; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": " << mdof[i/NModes] << "(" << i%NModes + 1 << ")";
			desc[i] = os.str();
		}

		for (unsigned i = 0; i < modulo*NStrNodes; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": StructNode(" << pInterfaceNodes[i/modulo]->GetLabel() << ") "
				<< rdof[(i/3)%(modulo/3)] << xyz[i%3];
			desc[2*NModes + i] = os.str();
		}

	} else {
		if (unsigned(i) < 2*NModes) {
			// modes
			os << ": " << mdof[i/NModes] << "(" << i%NModes + 1 << ")";

		} else {
			// reactions
			i -= 2*NModes;

			if (unsigned(i) >= modulo*NStrNodes) {
				// error
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			os << ": StructNode(" << pInterfaceNodes[i/modulo]->GetLabel() << ") "
				<< rdof[(i/3)%(modulo/3)] << xyz[i%3];
		}
		desc[0] = os.str();
	}
}

std::ostream&
Modal::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
	integer iModalIndex = iGetFirstIndex();

	out 
		<< prefix << iModalIndex + 1 << "->" << iModalIndex + NModes
		<< ": modal velocity definitions [q(1.." << NModes << ")=qP(1.." << NModes << ")]" << std::endl
		<< prefix << iModalIndex + NModes + 1 << "->" << iModalIndex + 2*NModes
		<< ": modal equilibrium equations [1.." << NModes << "]" << std::endl;
	iModalIndex += 2*NModes;
	for (unsigned iStrNodem1 = 0; iStrNodem1 < NStrNodes; iStrNodem1++, iModalIndex += 6) {
		out
			<< prefix << iModalIndex + 1 << "->" << iModalIndex + 3 << ": "
				"StructNode(" << pInterfaceNodes[iStrNodem1]->GetLabel() << ") "
				"position constraints [Px=PxN,Py=PyN,Pz=PzN]" << std::endl
			<< prefix << iModalIndex + 4 << "->" << iModalIndex + 6 << ": "
				"StructNode(" << pInterfaceNodes[iStrNodem1]->GetLabel() << ") "
				"orientation constraints [gx=gxN,gy=gyN,gz=gzN]" << std::endl;
		if (bInitial) {
			iModalIndex += 6;
			out
				<< prefix << iModalIndex + 1 << "->" << iModalIndex + 3 << ": "
					"StructNode(" << pInterfaceNodes[iStrNodem1]->GetLabel() << ") "
					"velocity constraints [vx=vxN,vy=vyN,vz=vzN]" << std::endl
				<< prefix << iModalIndex + 4 << "->" << iModalIndex + 6 << ": "
					"StructNode(" << pInterfaceNodes[iStrNodem1]->GetLabel() << ") "
					"angular velocity constraints [wx=wxN,wy=wyN,wz=wzN]" << std::endl;
		}
	}

	return out;
}

void
Modal::DescribeEq(std::vector<std::string>& desc, bool bInitial, int i) const
{
	int iend = 1;
	if (i == -1) {
		iend = 2*NModes + 6*NStrNodes;
		if (bInitial) {
			iend += 6*NStrNodes;
		}
	}
	desc.resize(iend);

	unsigned modulo = 6;
	if (bInitial) {
		modulo = 12;
	}

	std::ostringstream os;
	os << "Modal(" << GetLabel() << ")";

	if (i < -1 ) {
		// error
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);

	} else if (i == -1) {
		std::string name(os.str());

		for (unsigned i = 0; i < 2*NModes; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": " << meq[i/NModes] << "(" << i%NModes + 1 << ")";
			desc[i] = os.str();
		}

		for (unsigned i = 0; i < modulo*NStrNodes; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": StructNode(" << pInterfaceNodes[i/modulo]->GetLabel() << ") "
				<< req[(i/3)%(modulo/3)] << xyz[i%3];
			desc[2*NModes + i] = os.str();
		}

	} else {
		if (unsigned(i) < 2*NModes) {
			// modes
			os << ": " << meq[i/NModes] << "(" << i%NModes + 1 << ")";

		} else {
			// reactions
			i -= 2*NModes;

			if (unsigned(i) >= modulo*NStrNodes) {
				// error
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			os << ": StructNode(" << pInterfaceNodes[i/modulo]->GetLabel() << ") "
				<< req[(i/3)%(modulo/3)] << xyz[i%3];
		}
		desc[0] = os.str();
	}
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
			WM.PutRowIndex(iCnt, iRigidIndex + iCnt);
			WM.PutColIndex(iCnt, iRigidIndex + iCnt);
		}
	}

	/* indici della parte deformabile e delle reazioni vincolari */
	integer iFlexIndex = iGetFirstIndex();
	unsigned int iNumDof = iGetNumDof();

	for (unsigned int iCnt = 1; iCnt <= iNumDof; iCnt++) {
		WM.PutRowIndex(iRigidOffset + iCnt, iFlexIndex + iCnt);
		WM.PutColIndex(iRigidOffset + iCnt, iFlexIndex + iCnt);
	}

	/* indici delle equazioni vincolari (solo per il nodo 2) */
	for (unsigned int iStrNodem1 = 0; iStrNodem1 < NStrNodes; iStrNodem1++) {
		integer iNodeFirstMomIndex =
			pInterfaceNodes[iStrNodem1]->iGetFirstMomentumIndex();
		integer iNodeFirstPosIndex =
			pInterfaceNodes[iStrNodem1]->iGetFirstPositionIndex();

		integer iOffset = iRigidOffset + iNumDof + 6*iStrNodem1;
		for (integer iCnt = 1; iCnt <= 6; iCnt++) {
			WM.PutRowIndex(iOffset + iCnt, iNodeFirstMomIndex + iCnt);
			WM.PutColIndex(iOffset + iCnt, iNodeFirstPosIndex + iCnt);
		}
	}

	/* Assemblaggio dello Jacobiano */

	Vec3 wr(Zero3);
	Mat3x3 J(Zero3x3);
	Vec3 S(Zero3);
	if (pModalNode) {
		/* recupera e aggiorna i dati necessari */
		/* i prodotti Inv3j*aj ecc. sono gia' stati fatti da AssRes() */

		wr = pModalNode->GetWRef();
		/* R and RT are updated by AssRes() */
		Mat3x3 JTmp = Inv7;
		if (pInv8 != 0) {
			JTmp += Inv8jaj.Symm2();
		}
		J = R*JTmp*RT;
		Vec3 STmp = Inv2;
		if (pInv3 != 0) {
			STmp += Inv3jaj;
		}
		S = R*STmp;

		/* matrice di massa:        J[1,1] = Mtot  */
		for (int iCnt = 6 + 1; iCnt <= 6 + 3; iCnt++) {
			WM.PutCoef(iCnt, iCnt, dMass);
		}

		/* momenti statici J[1,2] = -[S/\] + c[-2w/\S/\ + S/\w/\] */
		WM.Add(6 + 1, 9 + 1, Mat3x3(MatCrossCross, S, wr*dCoef)
			- Mat3x3(MatCrossCross, wr, S*(2.*dCoef))
			- Mat3x3(MatCross, S));

		/* J[2,1] = [S/\] */
		WM.Add(9 + 1, 6 + 1, Mat3x3(MatCross, S));

		/* momenti d'inerzia:       J[2,2] = J + c[ - (Jw)/\ + (w/\)J];    */
		WM.Add(9 + 1, 9 + 1, J + ((wr.Cross(J) - Mat3x3(MatCross, J*wr))*dCoef));

		/* completa lo Jacobiano con l'aggiunta delle equazioni {xP} = {v}
		 {gP} - [wr/\]{g} = {w} */
		for (int iCnt = 1; iCnt <= 6; iCnt++) {
			WM.IncCoef(iCnt, iCnt, 1.);
			WM.DecCoef(iCnt, 6 + iCnt, dCoef);
		}
		WM.Sub(3 + 1, 3 + 1, Mat3x3(MatCross, wr*dCoef));
	}

	/* parte deformabile :
	 *
	 * | I  -cI  ||aP|
	 * |         ||  |
	 * |cK  M + cC ||bP|
	 */

	for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
		unsigned int iiCnt = iRigidOffset + iCnt;

		WM.PutCoef(iiCnt, iiCnt, 1.);
		WM.PutCoef(iiCnt, iiCnt + NModes, -dCoef);
		for (unsigned int jCnt = 1; jCnt <= NModes; jCnt++) {
			unsigned int jjCnt = iRigidOffset + jCnt;

			WM.PutCoef(iiCnt + NModes, jjCnt,
					dCoef*pModalStiff->dGet(iCnt, jCnt));
			WM.PutCoef(iiCnt + NModes, jjCnt + NModes,
					pModalMass->dGet(iCnt, jCnt)
					+ dCoef*pModalDamp->dGet(iCnt, jCnt));
		}
	}

	if (pModalNode) {

		/* termini di accoppiamento moto rigido-deformabile;
		 * eventualmente l'utente potra' scegliere 
		 * se trascurarli tutti, una parte o considerarli tutti */

		/* linearizzazione delle OmegaPrime:
		 * J13 = R*Inv3jaj
		 * J23 = R*Inv4 + Inv5jaj
		 * (questi termini ci vogliono sempre)
		 */

		if (pInv3 != 0) {
			Mat3xN Jac13(NModes, 0.);
			WM.Add(6 + 1, iRigidOffset + NModes + 1, Jac13.LeftMult(R, *pInv3));

			MatNx3 Jac13T(NModes, 0.);
			WM.Add(iRigidOffset + NModes + 1, 6 + 1, Jac13T.Transpose(Jac13));

			Mat3x3 Inv3jaPjWedge(MatCross, Inv3jaPj);
			WM.Sub(6 + 1, 9 + 1, R*Inv3jaPjWedge*(RT*(2.*dCoef)));
		}

		if (pInv4 != 0) {
			Mat3xN Jac23(NModes, 0.);
			Jac23.LeftMult(R, *pInv4);
			if (pInv5 != 0) {
				Mat3xN Inv5jajRef(NModes, 0.);

				Jac23 += Inv5jajRef.LeftMult(R, Inv5jaj);
			}
			WM.Add(9 + 1, iRigidOffset + NModes + 1, Jac23);

			MatNx3 Jac23T(NModes, 0.);
			WM.Add(iRigidOffset + NModes + 1, 9 + 1, Jac23T.Transpose(Jac23));
		}


		/* termini di Coriolis: linearizzazione delle Omega
		 * (si puo' evitare se non ho grosse vel. angolari):
		 * J13 = -2*R*[Inv3jaPj/\]*RT
		 * J23 = 2*R*[Inv8jaPj - Inv9jkajaPk]*RT */

		if (pInv8 != 0 ) {
			Mat3x3 JTmp = Inv8jaPj;
			if (pInv9 != 0) {
				JTmp -= Inv9jkajaPk;
			}
			WM.Add(9 + 1, 9 + 1, R*(JTmp*(RT*(2.*dCoef))));
		}

		/* termini di Coriolis: linearizzazione delle b;
		 * si puo' evitare 'quasi' sempre: le velocita' 
		 * di deformazione dovrebbero essere sempre piccole
		 * rispetto ai termini rigidi
		 * Jac13 = 2*[Omega/\]*R*PHI
		 */
#if 0
		Jac13.LeftMult(wrWedge*R*2*dCoef, *pInv3);
		WM.Add(6 + 1, iRigidOffset + NModes + 1, Jac13);
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
		for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
			Inv4j = pInv4->GetVec(jMode);
			VInv5jaj = Inv5jaj.GetVec(jMode);
			VInv5jaPj = Inv5jaPj.GetVec(jMode);
			unsigned int jOffset = (jMode - 1)*3 + 1;
			Inv8jTranspose = (pInv8->GetMat3x3(jOffset)).Transpose();
			if (pInv9 != 0 ) {
				Inv9jkak = 0.;
				for (unsigned int kModem1 = 0; kModem1 < NModes; kModem1++)  {
		 			doublereal a_kMode = a.dGet(kModem1 + 1);
		 			integer iOffset = (jMode - 1)*3*NModes + kModem1*3 + 1;
					Inv9jkak += pInv9->GetMat3x3ScalarMult(iOffset, a_kMode);
				}
			}
			for (int iCnt = 1; iCnt <= 3; iCnt++) {
				doublereal temp1 = 0., temp2 = 0.;
				for (int jCnt = 1; jCnt<=3; jCnt++) {
					if (pInv9 != 0) {
		 				temp1 += -2*wr.dGet(jCnt)*(R*(Inv8jTranspose - Inv9jkak)*RT).dGet(iCnt, jCnt);
					} else {
		 				temp1 += -2*wr.dGet(jCnt)*(R*(Inv8jTranspose)*RT).dGet(iCnt, jCnt);
					}
		 			temp2 += -(R*(Inv4j + VInv5jaj)).dGet(jCnt)*wrWedge.dGet(iCnt, jCnt);
				}
				WM.IncCoef(iRigidOffset + NModes + jMode, 9 + iCnt,
						dCoef*(temp1 + temp2));
			}
		 	for (int iCnt = 1; iCnt<=3; iCnt++) {
				doublereal temp1 = 0.;
				for (int jCnt = 1; jCnt<=3; jCnt++) {
		 			temp1 += (R*VInv5jaPj*2).dGet(jCnt);
				}
				WM.IncCoef(iRigidOffset + NModes + jMode, 9 + iCnt,
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
		for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
			integer iOffset = (jMode - 1)*NStrNodes + iStrNode;

			PHIt.PutVec(jMode, pPHIt->GetVec(iOffset));
			PHIr.PutVec(jMode, pPHIr->GetVec(iOffset));
		}

		MatNx3 PHItT(NModes), PHIrT(NModes);
		PHItT.Transpose(PHIt);
		PHIrT.Transpose(PHIr);

		/* nota: le operazioni
		 * d1tot = d1 + PHIt*a, R1tot = R*[I + (PHIr*a)/\]
		 * sono gia' state fatte da AssRes */

		Mat3xN SubMat1(NModes), SubMat2(NModes);
		MatNx3 SubMat3(NModes);
		MatNxN SubMat4(NModes);

		/* cerniera sferica */
		/* F e' aggiornata da AssRes */

		integer iReactionIndex = iRigidOffset + 2*NModes + 6*iStrNodem1;
		integer iStrNodeIndex = iRigidOffset + iNumDof + 6*iStrNodem1;

		Mat3x3 FTmpWedge(MatCross, pF[iStrNodem1]*dCoef);

		Mat3xN PHItR(NModes);
		PHItR.LeftMult(R, PHIt);

		for (int iCnt = 1; iCnt <= 3; iCnt++) {
			/* termini di reazione sui nodi */
			WM.DecCoef(iStrNodeIndex + iCnt, iReactionIndex + iCnt, 1.);

			/* termini di vincolo dovuti ai nodi */
			WM.DecCoef(iReactionIndex + iCnt, iStrNodeIndex + iCnt, 1.);
		}

		if (pModalNode) {
			for (int iCnt = 1; iCnt <= 3; iCnt++) {
				/* termini di reazione sul nodo modale */
				WM.IncCoef(6 + iCnt, iReactionIndex + iCnt, 1.);
		
				/* termini di vincolo dovuti al nodo modale */
				WM.IncCoef(iReactionIndex + iCnt, iCnt, 1.);
			}

			/* pd1Tot e' il puntatore all'array
			 * che contiene le posizioni del nodi FEM */
			Mat3x3 dTmp1Wedge(MatCross, R*pd1tot[iStrNodem1]);

			WM.Add(9 + 1, iReactionIndex + 1, dTmp1Wedge);
			
			/* termini del tipo: c*F/\*d/\*Deltag */
			WM.Add(9 + 1, 3 + 1, FTmpWedge*dTmp1Wedge);
			
			/* termini di vincolo dovuti ai nodi */
			WM.Sub(iReactionIndex + 1, 3 + 1, dTmp1Wedge);
			
			/* termini aggiuntivi dovuti alla deformabilita' */
			SubMat3.RightMult(PHItT, RT*FTmpWedge);
			WM.Add(iRigidOffset + NModes + 1, 3 + 1, SubMat3);
			
			SubMat1.LeftMult(FTmpWedge, PHItR);
			WM.Sub(9 + 1, iRigidOffset + 1, SubMat1);
		}
			
		/* contributo dovuto alla flessibilita' */
		WM.Add(iReactionIndex + 1, iRigidOffset + 1, PHItR);

		/* idem per pd2 e pR2 */
		Mat3x3 dTmp2Wedge(MatCross, pR2[iStrNodem1]*pd2[iStrNodem1]);
		WM.Sub(iStrNodeIndex + 3 + 1, iReactionIndex + 1, dTmp2Wedge);

		/* termini del tipo: c*F/\*d/\*Deltag */
		WM.Sub(iStrNodeIndex + 3 + 1, iStrNodeIndex + 3 + 1,
				FTmpWedge*dTmp2Wedge);

		/* termini aggiuntivi dovuti alla deformabilita' */
		SubMat3.RightMult(PHItT, RT);
		WM.Add(iRigidOffset + NModes + 1, iReactionIndex + 1, SubMat3);

		/* modifica: divido le equazioni di vincolo per dCoef */

		/* termini di vincolo dovuti al nodo 2 */
		WM.Add(iReactionIndex + 1, iStrNodeIndex + 3 + 1, dTmp2Wedge);

		/* fine cerniera sferica */

		/* equazioni di vincolo : giunto prismatico */

		/* Vec3 M(XCurr, iModalIndex + 2*NModes + 6*iStrNodem1 + 3 + 1); */
		Vec3 MTmp = pR2[iStrNodem1]*(pM[iStrNodem1]*dCoef);
		Mat3x3 MTmpWedge(MatCross, MTmp);

		/* Eq. dei momenti */
		if (pModalNode) {
			WM.Sub(9 + 1, iStrNodeIndex + 3 + 1, MTmpWedge);
		}
		WM.Add(iStrNodeIndex + 3 + 1, iStrNodeIndex + 3 + 1, MTmpWedge);

		/* Eq. d'equilibrio ai modi */
		SubMat3.RightMult(PHIrT, R*MTmpWedge);
		if (pModalNode) {
			WM.Add(iRigidOffset + NModes + 1, 3 + 1, SubMat3);
		}
		WM.Sub(iRigidOffset + NModes + 1, iStrNodeIndex + 3 + 1, SubMat3);

		/* */
		Mat3x3 R2T = pR2[iStrNodem1].Transpose();
		if (pModalNode) {
			WM.Add(iReactionIndex + 3 + 1, 3 + 1, R2T);
			WM.Add(9 + 1, iReactionIndex + 3 + 1, pR2[iStrNodem1]);
		}
		WM.Sub(iReactionIndex + 3 + 1, iStrNodeIndex + 3 + 1, R2T);
		WM.Sub(iStrNodeIndex + 3 + 1, iReactionIndex + 3 + 1, pR2[iStrNodem1]);

		/* */
		SubMat3.RightMult(PHIrT, RT*pR2[iStrNodem1]);
		WM.Add(iRigidOffset + NModes + 1, iReactionIndex + 3 + 1, SubMat3);
		for (unsigned iMode = 1; iMode <= NModes; iMode++) {
			for (unsigned jIdx = 1; jIdx <= 3; jIdx++) {
				WM.IncCoef(iReactionIndex + 3 + jIdx, iRigidOffset + iMode,
					SubMat3(iMode, jIdx));
			}
		}
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
	 * 13 		-> 12 + 2*NM:	modes eq. (\dot{a}=b); m\dot{b} + ka=f)
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
		WorkVec.PutRowIndex(iRigidOffset + iCnt, iModalIndex + iCnt);
	}

	/* interface nodes equilibrium indices */
	integer iOffset = iRigidOffset + iNumDof;
	for (unsigned iStrNodem1 = 0; iStrNodem1 < NStrNodes; iStrNodem1++) {
		integer iNodeFirstMomIndex =
			pInterfaceNodes[iStrNodem1]->iGetFirstMomentumIndex();

		for (unsigned int iCnt = 1; iCnt <= 6; iCnt++) {
			WorkVec.PutRowIndex(iOffset + iCnt,
					iNodeFirstMomIndex + iCnt);
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

	Vec3 Inv3jaPPj(Zero3);
	if (pInv3 != 0) {
		Inv3jaj = *pInv3 * a;
		Inv3jaPj = *pInv3 * b;
		Inv3jaPPj = *pInv3 * bPrime;
	}

	/* invarianti rotazionali */
	if (pInv8 != 0) {
		Inv8jaj.Reset();
		Inv8jaPj.Reset();
	}
	if (pInv5 != 0) {
		Inv5jaj.Reset(0.);
		Inv5jaPj.Reset(0.);
	}
	Mat3x3 MatTmp1(Zero3x3), MatTmp2(Zero3x3);
	if (pInv9) {
		Inv9jkajak.Reset();
		Inv9jkajaPk.Reset();
	}
	Mat3x3 Inv10jaPj(Zero3x3);

	if (pInv5 != 0 || pInv8 != 0 || pInv9 != 0 || pInv10 != 0) {
		for (unsigned int jMode = 1; jMode <= NModes; jMode++)  {
			doublereal a_jMode = a.dGet(jMode);
			doublereal aP_jMode = b.dGet(jMode);

			if (pInv5 != 0) {
				for (unsigned int kMode = 1; kMode <= NModes; kMode++) {
					Vec3 v = pInv5->GetVec((jMode - 1)*NModes + kMode);

					Inv5jaj.AddVec(kMode, v*a_jMode);
					Inv5jaPj.AddVec(kMode, v*aP_jMode);
				}
			}

			if (pInv8 != 0) {
				Mat3x3 Inv8jajTmp;
				
				for (unsigned int jCnt = 1; jCnt <= 3; jCnt++) {
					unsigned int iMjC = (jMode - 1)*3 + jCnt;

					Inv8jajTmp.PutVec(jCnt, pInv8->GetVec(iMjC));
				}

				Inv8jaj += Inv8jajTmp * a_jMode;
				Inv8jaPj += Inv8jajTmp * aP_jMode;

				if (pInv9 != 0) {
					/*
					 * questi termini si possono commentare perche' sono
					 * (sempre ?) piccoli (termini del tipo a*a o a*b)
					 * eventualmente dare all'utente la possibilita'
					 * di scegliere se trascurarli o no
					 */
					for (unsigned int kMode = 1; kMode <= NModes; kMode++) {
						doublereal a_kMode = a.dGet(kMode);
						doublereal aP_kMode = b.dGet(kMode);
						unsigned int iOffset = (jMode - 1)*3*NModes + (kMode - 1)*3 + 1;
						Inv9jkajak += pInv9->GetMat3x3ScalarMult(iOffset, a_jMode*a_kMode);
						Inv9jkajaPk += pInv9->GetMat3x3ScalarMult(iOffset, a_jMode*aP_kMode);
					}
				}
			}

			if (pInv10 != 0) {
				Mat3x3 Inv10jaPjTmp;

				for (unsigned int jCnt = 1; jCnt <= 3; jCnt++) {
					unsigned int iMjC = (jMode - 1)*3 + jCnt;

					Inv10jaPjTmp.PutVec(jCnt, pInv10->GetVec(iMjC)*aP_jMode);
				}

			
				Inv10jaPj += Inv10jaPjTmp;
			}
		}
	} /* fine ciclo sui modi */

#ifdef MODAL_USE_GRAVITY
	/* forza di gravita' (decidere come inserire g) */
	/* FIXME: use a reasonable reference point where compute gravity */
	Vec3 GravityAcceleration(Zero3);
	bool bGravity = GravityOwner::bGetGravity(x, GravityAcceleration);
#endif /* MODAL_USE_GRAVITY */

	Vec3 vP(Zero3);
	Vec3 wP(Zero3);
	Vec3 w(Zero3);
	Vec3 RTw(Zero3);
	if (pModalNode) {
		/* rigid body indices */
		integer iRigidIndex = pModalNode->iGetFirstIndex();
		for (unsigned int iCnt = 1; iCnt <= iRigidOffset; iCnt++) {
			WorkVec.PutRowIndex(iCnt, iRigidIndex + iCnt);
		}

		/* recupera i dati necessari */
		x = pModalNode->GetXCurr();
		Vec3 xP = pModalNode->GetVCurr();
		Vec3 g = pModalNode->GetgCurr();
		Vec3 gP = pModalNode->GetgPCurr();
		vP = pModalNode->GetXPPCurr();
		Vec3 wr = pModalNode->GetWRef();
		wP = pModalNode->GetWPCurr();
		
		Vec3 v(XCurr, iRigidIndex + 6 + 1);
		w = Vec3(XCurr, iRigidIndex + 9 + 1);

		R = pModalNode->GetRCurr();
		RT = R.Transpose();
		RTw = RT*w;

		Mat3x3 JTmp = Inv7;
		if (pInv8 != 0) {
			JTmp += Inv8jaj.Symm2();
			if (pInv9 != 0) {
				JTmp -= Inv9jkajak;
			}
		}
		Mat3x3 J = R*JTmp*RT;
		Vec3 STmp = Inv2;
		if (pInv3 != 0) {
			STmp += Inv3jaj;
		}
		Vec3 S = R*STmp;

		/*
		 * fine aggiornamento invarianti
		 */

		/* forze d'inerzia */
		Vec3 FTmp = vP*dMass - S.Cross(wP) + w.Cross(w.Cross(S));
		if (pInv3 != 0) {
			FTmp += R*Inv3jaPPj + (w.Cross(R*Inv3jaPj))*2.;
		}
		WorkVec.Sub(6 + 1, FTmp);

#if 0
		std::cerr << "m=" << dMass << "; S=" << S
			<< "; a="  << a << "; aPrime =" << aPrime
			<< "; b=" << b <<  "; bPrime= " << bPrime
			<< "; tot=" << vP*dMass - S.Cross(wP) + w.Cross(w.Cross(S))
			+ (w.Cross(R*Inv3jaPj))*2 + R*Inv3jaPPj << std::endl;
#endif

		FTmp = S.Cross(vP) + J*wP + w.Cross(J*w);
		if (pInv4 != 0) {
			Mat3xN Inv4Curr(NModes, 0);
			Inv4Curr.LeftMult(R, *pInv4);
			FTmp += Inv4Curr*bPrime;
		}
		if (pInv5 != 0) {
			Mat3xN Inv5jajCurr(NModes, 0);
			Inv5jajCurr.LeftMult(R, Inv5jaj);
			FTmp += Inv5jajCurr*bPrime;
		}
		if (pInv8 != 0) {
			Mat3x3 MTmp = Inv8jaPj;
			if (pInv9 != 0) {
				MTmp -= Inv9jkajaPk;
			}
			FTmp += R*MTmp*(RTw*2.);
		}
		/* termini dovuti alle inerzie rotazionali */
		if (pInv10 != 0) {
			Vec3 VTmp = Inv10jaPj.Symm2()*RTw;
			if (pInv11 != 0) {
				VTmp += w.Cross(R*(*pInv11*b));
			}
			FTmp += R*VTmp;
		}
		WorkVec.Sub(9 + 1, FTmp);


#ifdef MODAL_USE_GRAVITY
		/* forza di gravita' (decidere come inserire g) */
		if (bGravity) {
			WorkVec.Add(6 + 1, GravityAcceleration*dMass);
			WorkVec.Add(9 + 1, S.Cross(GravityAcceleration));
		}
#endif /* MODAL_USE_GRAVITY */
	
		/* equazioni per abbassare di grado il sistema */
		WorkVec.Add(1, v - xP);
		WorkVec.Add(3 + 1, w - Mat3x3(CGR_Rot::MatG, g)*gP
			- Mat3x3(CGR_Rot::MatR, g)*wr);
	}

	/* forze modali */
	for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
		unsigned int jOffset = (jMode - 1)*3 + 1;
		doublereal d = - MaPP.dGet(jMode) - CaP.dGet(jMode) - Ka.dGet(jMode);

		if (pInv3 != 0) {
			Vec3 Inv3j = pInv3->GetVec(jMode);

			d -= (R*Inv3j).Dot(vP);
			
#ifdef MODAL_USE_GRAVITY
			/* forza di gravita': */
			if (bGravity) {
				WorkVec.IncCoef(iRigidOffset + NModes + jMode,
						(R*Inv3j).Dot(GravityAcceleration));
			}
#endif /* MODAL_USE_GRAVITY */
		}

		if (pInv4 != 0) {
			Vec3 Inv4j = pInv4->GetVec(jMode);
			if (pInv5 != 0) {
				VInv5jaj = Inv5jaj.GetVec(jMode);
				Inv4j += VInv5jaj;

				VInv5jaPj = Inv5jaPj.GetVec(jMode);
				d -= ((R*VInv5jaPj).Dot(w))*2.;
			}

			d -= (R*Inv4j).Dot(wP);
		}

		if (pInv8 != 0 || pInv9 != 0 || pInv10 != 0) {
			Mat3x3 MTmp(Zero3x3);

			if (pInv8 != 0) {
				MTmp += (pInv8->GetMat3x3(jOffset)).Transpose();
				if (pInv9 != 0) {
					for (unsigned int kModem1 = 0; kModem1 < NModes; kModem1++) {
						doublereal a_kMode = a.dGet(kModem1 + 1);
						integer kOffset = (jMode - 1)*3*NModes + kModem1*3 + 1;
	
						MTmp -= pInv9->GetMat3x3ScalarMult(kOffset, a_kMode);
					}
				}
			}

			if (pInv10 != 0) {
				MTmp += pInv10->GetMat3x3(jOffset);
			}

			d += w.Dot(R*(MTmp*RTw));
		}

		WorkVec.IncCoef(iRigidOffset + NModes + jMode, d);

#if 0
		std::cerr << "(R*Inv3j).Dot(vP)=" << (R*Inv3j).Dot(vP)
			<< std::endl
			<< "(R*(Inv4j + VInv5jaj)).Dot(wP)="
			<< (R*(Inv4j + VInv5jaj)).Dot(wP) << std::endl
			<< "w.Dot(R*((Inv8jTranspose + Inv10j)*RTw))"
			<< w.Dot(R*((Inv8jTranspose + Inv10j)*RTw)) << std::endl
			<< " -(R*VInv5jaPj).Dot(w)*2."
			<< -(R*VInv5jaPj).Dot(w)*2.<< std::endl
			<< " -CaP.dGet(jMode)" << -CaP.dGet(jMode)<< std::endl
			<< "-Ka.dGet(jMode) "
			<< -CaP.dGet(jMode) - Ka.dGet(jMode) << std::endl;
#endif
	}

	/* equazioni per abbassare di grado il sistema */
	for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
		WorkVec.PutCoef(iRigidOffset + iCnt, b.dGet(iCnt) - aPrime.dGet(iCnt));
	}

	/* equazioni di vincolo */
	for (unsigned int iStrNode = 1; iStrNode <= NStrNodes; iStrNode++) {
		unsigned int iStrNodem1 = iStrNode - 1;
		integer iReactionIndex = iRigidOffset + 2*NModes + 6*iStrNodem1;
		integer iStrNodeIndex = iReactionIndex + 6*NStrNodes;

		/* nodo 1 (FEM) e 2 (Multi - Body): recupera la posizione */
		Vec3 d1rig(pOffsetFEMNodes->GetVec(iStrNode));
		pd2[iStrNodem1] = pOffsetMBNodes->GetVec(iStrNode);

		/* recupero le forme modali del nodo vincolato */
		/* FIXME: what about using Blitz++ ? :) */
		Mat3xN PHIt(NModes), PHIr(NModes);
		for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
			integer iOffset = (jMode - 1)*NStrNodes + iStrNode;
			PHIt.PutVec(jMode, pPHIt->GetVec(iOffset));
			PHIr.PutVec(jMode, pPHIr->GetVec(iOffset));
		}

		/*
		 * aggiorno d1 e R con il contributo dovuto
		 * alla flessibilita':
		 * d1tot = d1 + PHIt*a, R1tot = R*[I + (PHIr*a)/\]
		 */
		pd1tot[iStrNodem1] = d1rig + PHIt*a;
		pR1tot[iStrNodem1] = R*Mat3x3(1., PHIr*a);

		/* constraint reaction (force) */
		pF[iStrNodem1] = Vec3(XCurr,
				iModalIndex + 2*NModes + 6*iStrNodem1 + 1);
		Vec3 x2 = pInterfaceNodes[iStrNodem1]->GetXCurr();

		/* FIXME: R2??? */
		pR2[iStrNodem1] = pInterfaceNodes[iStrNodem1]->GetRCurr();

		/* cerniera sferica */
		Vec3 dTmp1(R*pd1tot[iStrNodem1]);
		Vec3 dTmp2(pR2[iStrNodem1]*pd2[iStrNodem1]);

		/* Eq. d'equilibrio, nodo 1 */
		if (pModalNode) {
			WorkVec.Sub(6 + 1, pF[iStrNodem1]);
			WorkVec.Sub(9 + 1, dTmp1.Cross(pF[iStrNodem1]));
		}

		/* termine aggiuntivo dovuto alla deformabilita':
		 * -PHItiT*RT*F */
		Vec3 vtemp = RT*pF[iStrNodem1];
		for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
			/* PHItTranspose(jMode) * vtemp */
			WorkVec.DecCoef(iRigidOffset + NModes + jMode,
					PHIt.GetVec(jMode)*vtemp);
		}

		/* Eq. d'equilibrio, nodo 2 */
		WorkVec.Add(iStrNodeIndex + 1, pF[iStrNodem1]);
		WorkVec.Add(iStrNodeIndex + 3 + 1, dTmp2.Cross(pF[iStrNodem1]));

		/* Eq. di vincolo */
		if (dCoef != 0.) {
			WorkVec.Add(iReactionIndex + 1, (x2 + dTmp2 - x - dTmp1)/dCoef);
		}

		/* constraint reaction (moment) */
		pM[iStrNodem1] = Vec3(XCurr,
				iModalIndex + 2*NModes + 6*iStrNodem1 + 3 + 1);

		/* giunto prismatico */
		Vec3 ThetaCurr = RotManip::VecRot(pR2[iStrNodem1].MulTM(pR1tot[iStrNodem1]));
		Vec3 MTmp = pR2[iStrNodem1]*pM[iStrNodem1];

		/* Equazioni di equilibrio, nodo 1 */
		if (pModalNode) {
			WorkVec.Sub(9 + 1, MTmp);
		}

		/* Equazioni di equilibrio, nodo 2 */
		WorkVec.Add(iStrNodeIndex + 3 + 1, MTmp);

		/* Contributo dovuto alla flessibilita' :-PHIrT*RtotT*M */
		vtemp = RT*MTmp;
		for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
			/* PHIrTranspose(jMode) * vtemp */
			WorkVec.DecCoef(iRigidOffset + NModes + jMode,
					PHIr.GetVec(jMode)*vtemp);
		}

		/* Modifica: divido le equazioni di vincolo per dCoef */
		if (dCoef != 0.) {
			/* Equazioni di vincolo di rotazione */
			WorkVec.Sub(iReactionIndex + 3 + 1, ThetaCurr/dCoef);
		}
	}

#if 0
	std::cerr << "###" << std::endl;
	for (int i = 1; i <= WorkVec.iGetSize(); i++) {
		std::cerr << std::setw(2) << i << "("
			<< std::setw(2) << WorkVec.iGetRowIndex(i) << ")"
			<< std::setw(20) << WorkVec(i) << std::endl;
	}
#endif

	return WorkVec;
}

void
Modal::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		/* stampa sul file di output i modi */
		std::ostream& out = OH.Modal();

		for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
			out << std::setw(8) << GetLabel() << '.' << uModeNumber[iCnt - 1]
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
			WM.PutRowIndex(iCnt, iRigidIndex + iCnt);
			WM.PutColIndex(iCnt, iRigidIndex + iCnt);
		}
	}

	integer iFlexIndex = iGetFirstIndex();

	for (unsigned int iCnt = 1; iCnt <= iGetInitialNumDof(); iCnt++) {
		WM.PutRowIndex(iRigidOffset + iCnt, iFlexIndex + iCnt);
		WM.PutColIndex(iRigidOffset + iCnt, iFlexIndex + iCnt);
	}

	for (unsigned iStrNodem1 = 0; iStrNodem1 < NStrNodes; iStrNodem1++) {
		integer iNodeFirstPosIndex = 
			pInterfaceNodes[iStrNodem1]->iGetFirstPositionIndex();
		integer iNodeFirstVelIndex = iNodeFirstPosIndex + 6;

		integer iOffset = iRigidOffset + iGetInitialNumDof() + 12*iStrNodem1;
		for (unsigned int iCnt = 1; iCnt <= 6; iCnt++) {
			WM.PutRowIndex(iOffset + iCnt,
					iNodeFirstPosIndex + iCnt);
			WM.PutColIndex(iOffset + iCnt,
					iNodeFirstPosIndex + iCnt);
			WM.PutRowIndex(iOffset + 6 + iCnt,
					iNodeFirstVelIndex + iCnt);
			WM.PutColIndex(iOffset + 6 + iCnt,
					iNodeFirstVelIndex + iCnt);
		}
	}

	/* comincia l'assemblaggio dello Jacobiano */

	/* nota: nelle prime iRigidOffset + 2*M equazioni 
	 * metto dei termini piccoli per evitare singolarita'
	 * quando non ho vincoli esterni o ho azzerato i modi */
	for (unsigned int iCnt = 1; iCnt <= iRigidOffset + 2*NModes; iCnt++) {
		WM.PutCoef(iCnt, iCnt, 1.e-12);
	}

	/* forze elastiche e viscose */
	for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
		for (unsigned int jCnt = 1; jCnt <= NModes; jCnt++) {
			WM.PutCoef(iRigidOffset + iCnt, iRigidOffset + jCnt,
					pModalStiff->dGet(iCnt,jCnt));
			WM.PutCoef(iRigidOffset + iCnt, iRigidOffset + NModes + jCnt,
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
		Vec3 F(XCurr, iFlexIndex + 2*NModes + 12*iStrNodem1 + 1);
		Vec3 FPrime(XCurr, iFlexIndex + 2*NModes + 12*iStrNodem1 + 6 + 1);

		Mat3xN PHIt(NModes, 0), PHIr(NModes, 0);
		MatNx3 PHItT(NModes, 0), PHIrT(NModes, 0);

		for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
			PHIt.PutVec(jMode, pPHIt->GetVec((jMode - 1)*NStrNodes + iStrNode));
			PHIr.PutVec(jMode, pPHIr->GetVec((jMode - 1)*NStrNodes + iStrNode));
		}

		PHItT.Transpose(PHIt);
		PHIrT.Transpose(PHIr);

		Vec3 d1rig(pOffsetFEMNodes->GetVec(iStrNode));
		Vec3 d2(pOffsetMBNodes->GetVec(iStrNode));

		Vec3 d1tot = d1rig + PHIt*a;
		Mat3x3 R1tot = R*Mat3x3(1., PHIr*a);
		Mat3xN SubMat1(NModes, 0.);
		MatNx3 SubMat2(NModes, 0.);

		integer iOffset1 = iRigidOffset + 2*NModes + 12*iStrNodem1;
		integer iOffset2 = iRigidOffset + iGetInitialNumDof() + 12*iStrNodem1;
		if (pModalNode) {
			for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
				/* Contributo di forza all'equazione
				 * della forza, nodo 1 */
				WM.IncCoef(iCnt, iOffset1 + iCnt, 1.);

				/* Contrib. di der. di forza all'eq. della der. 
				 * della forza, nodo 1 */
				WM.IncCoef(6 + iCnt, iOffset1 + 6 + iCnt, 1.);

				/* Equazione di vincolo, nodo 1 */
				WM.DecCoef(iOffset1 + iCnt, iCnt, 1.);

				/* Derivata dell'equazione di vincolo, nodo 1 */
				WM.DecCoef(iOffset1 + 6 + iCnt, 6 + iCnt, 1.);
			}
		}

		for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
			/* Equazione di vincolo, nodo 2 */
			WM.IncCoef(iOffset1 + iCnt, iOffset2 + iCnt, 1.);

			/* Derivata dell'equazione di vincolo, nodo 2 */
			WM.IncCoef(iOffset1 + 6 + iCnt, iOffset2 + 6 + iCnt, 1.);

			/* Contributo di forza all'equazione della forza,
			 * nodo 2 */
			WM.DecCoef(iOffset2 + iCnt, iOffset1 + iCnt, 1.);

			/* Contrib. di der. di forza all'eq. della der.
			 * della forza, nodo 2 */
			WM.DecCoef(iOffset2 + 6 + iCnt, iOffset1 + 6 + iCnt, 1.);
		}

		/* Distanza nel sistema globale */
		Vec3 d1Tmp(R*d1tot);
		Vec3 d2Tmp(R2*d2);

		/* Matrici F/\d1/\, -F/\d2/\ , F/\omega1/\ */
		Mat3x3 FWedged1Wedge(MatCrossCross, F, d1Tmp);
		Mat3x3 FWedged2Wedge(MatCrossCross, F, d2Tmp);
		Mat3x3 FWedgeO1Wedge(MatCrossCross, F, Omega1);

		/* Matrici (omega1/\d1)/\, -(omega2/\d2)/\ */
		Mat3x3 O1Wedged1Wedge(MatCross, Omega1.Cross(d1Tmp));
		Mat3x3 O2Wedged2Wedge(MatCross, d2Tmp.Cross(Omega2));

		/* d1Prime= w1/\d1 + R*PHIt*b */
		Mat3xN R1PHIt(NModes);
		R1PHIt.LeftMult(R, PHIt);
		Vec3 d1Prime(Omega1.Cross(d1Tmp) + R1PHIt*b);

		Mat3x3 FWedge(MatCross, F);
		if (pModalNode) {
			/* Equazione di momento, nodo 1 */
			WM.Add(3 + 1, 3 + 1, FWedged1Wedge);
			WM.Add(3 + 1, iOffset1 + 1, Mat3x3(MatCross, d1Tmp));

			/* Equazione di momento, contributo modale */
			SubMat1.LeftMult(FWedge*R, PHIt);
			WM.Sub(3 + 1, iRigidOffset + 1, SubMat1);
			
			/* Eq. di equilibrio ai modi */
			SubMat2.RightMult(PHItT, RT*FWedge);
			WM.Add(iRigidOffset + 1, 3 + 1, SubMat2);
		}

		/* Equazione di momento, nodo 2 */
		WM.Sub(iOffset2 + 3 + 1,iOffset2 + 3 + 1, FWedged2Wedge);
		WM.Sub(iOffset2 + 3 + 1,iOffset1 + 1, Mat3x3(MatCross, d2Tmp));

		/* Eq. di equilibrio ai modi */
		SubMat2.RightMult(PHItT, RT);
		WM.Add(iRigidOffset + 1, iRigidOffset + 2*NModes + 1, SubMat2);

		if (pModalNode) {
			/* derivata dell'equazione di momento, nodo 1 */
			WM.Add(9 + 1, 3 + 1,
				Mat3x3(MatCrossCross, FPrime, d1Tmp) + F.Cross(Mat3x3(MatCrossCross, Omega1, d1Tmp))
				+ Mat3x3(MatCrossCross, F, R*(PHIt*b)));
			WM.Add(9 + 1, 9 + 1, FWedged1Wedge);
			WM.Add(9 + 1, iOffset1 + 1, O1Wedged1Wedge + Mat3x3(MatCross, R*(PHIt*b)));
			WM.Add(9 + 1, iOffset1 + 6 + 1, Mat3x3(MatCross, d1Tmp));

			/* derivata dell'equazione di momento, contributo modale */
			SubMat1.LeftMult((-FWedgeO1Wedge - Mat3x3(MatCross, FPrime))*R, PHIt);
			WM.Add(9 + 1, iRigidOffset + 1, SubMat1);
			SubMat1.LeftMult(-FWedge*R, PHIt);
			WM.Add(9 + 1, iRigidOffset + NModes + 1, SubMat1);

			/* derivata dell'eq. di equilibrio ai modi */
			SubMat2.RightMult(PHItT, RT*FWedge);
			WM.Add(iRigidOffset + NModes + 1, 9 + 1, SubMat2);
			SubMat2.RightMult(PHItT, RT*FWedgeO1Wedge);
			WM.Add(iRigidOffset + NModes + 1, 3 + 1, SubMat2);
		}
		
		SubMat2.RightMult(PHItT, RT);
		WM.Add(iRigidOffset + NModes + 1, iRigidOffset + 2*NModes + 6 + 1, SubMat2);

		/* Derivata dell'equazione di momento, nodo 2 */
		WM.Add(iOffset2 + 9 + 1, iOffset2 + 3 + 1,
			Mat3x3(MatCrossCross, FPrime, -d2Tmp)
			- F.Cross(Mat3x3(MatCrossCross, Omega2, d2Tmp)));
		WM.Sub(iOffset2 + 9 + 1, iOffset2 + 9 + 1, FWedged2Wedge);
		WM.Add(iOffset2 + 9 + 1, iOffset1 + 1, O2Wedged2Wedge);
		WM.Sub(iOffset2 + 9 + 1, iOffset1 + 6 + 1, Mat3x3(MatCross, d2Tmp));

		/* Equazione di vincolo */
		if (pModalNode) {
			WM.Add(iOffset1 + 1, 3 + 1, Mat3x3(MatCross, d1Tmp));
		}
		WM.Sub(iOffset1 + 1, iOffset2 + 3 + 1, Mat3x3(MatCross, d2Tmp));

		/* Equazione di vincolo, contributo modale */
		SubMat1.LeftMult(R, PHIt);
		WM.Sub(iOffset1 + 1, iRigidOffset + 1, SubMat1);

		/* Derivata dell'equazione di vincolo */
		if (pModalNode) {
			WM.Add(iOffset1 + 6 + 1, 3 + 1, O1Wedged1Wedge + Mat3x3(MatCross, R*(PHIt*b)));
			WM.Add(iOffset1 + 6 + 1, 9 + 1, Mat3x3(MatCross, d1Tmp));
		}
		WM.Add(iOffset1 + 6 + 1, iOffset2 + 3 + 1, O2Wedged2Wedge);
		WM.Sub(iOffset1 + 6 + 1, iOffset2 + 9 + 1, Mat3x3(MatCross, d2Tmp));

		/* Derivata dell'equazione di vincolo, contributo modale */
		SubMat1.LeftMult(Omega1.Cross(R), PHIt);
		WM.Sub(iOffset1 + 6 + 1, iRigidOffset + 1, SubMat1);
		SubMat1.LeftMult(R, PHIt);
		WM.Sub(iOffset1 + 6 + 1, iRigidOffset + NModes + 1, SubMat1);

		/* giunto prismatico */
		Vec3 M(XCurr, iFlexIndex + 2*NModes + 12*iStrNodem1 + 3 + 1);
		Vec3 MPrime(XCurr, iFlexIndex + 2*NModes + 12*iStrNodem1 + 9 + 1);
		Vec3 e1tota(R1tot.GetVec(1));
		Vec3 e2tota(R1tot.GetVec(2));
		Vec3 e3tota(R1tot.GetVec(3));
		Vec3 e1b(R2.GetVec(1));
		Vec3 e2b(R2.GetVec(2));
		Vec3 e3b(R2.GetVec(3));

		/* */
		Mat3x3 MWedge(Mat3x3(MatCrossCross, e3b, e2tota*M.dGet(1))
			+ Mat3x3(MatCrossCross, e1b, e3tota*M.dGet(2))
			+ Mat3x3(MatCrossCross, e2b, e1tota*M.dGet(3)));
		Mat3x3 MWedgeT(MWedge.Transpose());

		/* Equilibrio */
		if (pModalNode) {
			WM.Add(3 + 1, 3 + 1, MWedge);
			WM.Sub(3 + 1, iOffset2 + 3 + 1, MWedgeT);

			WM.Add(iOffset2 + 3 + 1, 3 + 1, MWedgeT);
		}
		WM.Sub(iOffset2 + 3 + 1, iOffset2 + 3 + 1, MWedge);

		/* Equilibrio dei momenti, termini aggiuntivi modali */
		Vec3 e1ra(R.GetVec(1));
		Vec3 e2ra(R.GetVec(2));
		Vec3 e3ra(R.GetVec(3));
		Mat3x3 M1Wedge(Mat3x3(MatCrossCross, e3b, e2ra*M.dGet(1))
			+ Mat3x3(MatCrossCross, e1b, e3ra*M.dGet(2))
			+ Mat3x3(MatCrossCross, e2b, e1ra*M.dGet(3)));

		SubMat1.LeftMult(M1Wedge, PHIr);
		if (pModalNode) {
			WM.Add(3 + 1, iRigidOffset + 1, SubMat1);
		}
		WM.Sub(iOffset2 + 3 + 1, iRigidOffset + 1, SubMat1);

		/* Equilibrio ai modi */
		if (pModalNode) {
			SubMat2.RightMult(PHIrT, R1tot.MulTM(MWedge));
			WM.Add(iRigidOffset + 1, 3 + 1, SubMat2);
		}
		SubMat2.RightMult(PHIrT, R1tot.MulTM(MWedgeT));
		WM.Sub(iRigidOffset + 1, iOffset2 + 3 + 1, SubMat2);

		Vec3 MA(Mat3x3(e2tota.Cross(e3b), e3tota.Cross(e1b),
					e1tota.Cross(e2b))*M);
		if (pModalNode) {
			SubMat2.RightMult(PHIrT, R1tot.MulTM(Mat3x3(MatCross, MA)));
			WM.Add(iRigidOffset + 1, 3 + 1, SubMat2);
		}

		Mat3x3 R1TMAWedge(MatCross, RT*MA);
		SubMat1.LeftMult(R1TMAWedge, PHIr);
		Mat3xN SubMat3(NModes);
		MatNxN SubMat4(NModes);
		SubMat3.LeftMult(R1tot.MulTM(M1Wedge), PHIr);
		SubMat1 += SubMat3;
		SubMat4.Mult(PHIrT, SubMat1);
		for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
			for (unsigned int kMode = 1; kMode <= NModes; kMode++) {
				WM.IncCoef(iRigidOffset + jMode, iRigidOffset + kMode,
						SubMat4.dGet(jMode,kMode));
			}
		}

		/* Derivate dell'equilibrio */
		if (pModalNode) {
			WM.Add(9 + 1, 9 + 1, MWedge);
			WM.Sub(9 + 1, iOffset2 + 9 + 1, MWedgeT);

			WM.Add(iOffset2 + 9 + 1, 9 + 1, MWedgeT);
		}
		WM.Sub(iOffset2 + 9 + 1, iOffset2 + 9 + 1, MWedge);

		MWedge = ((Mat3x3(MatCrossCross, e3b, Omega1) + Mat3x3(MatCross, Omega2.Cross(e3b*M(1))))
			+ Mat3x3(MatCross, e3b*MPrime(1)))*Mat3x3(MatCross, e2tota)
			+ ((Mat3x3(MatCrossCross, e1b, Omega1) + Mat3x3(MatCross, Omega2.Cross(e1b*M(2))))
			+ Mat3x3(MatCross, e1b*MPrime(2)))*Mat3x3(MatCross, e3tota)
			+ ((Mat3x3(MatCrossCross, e2b, Omega1) + Mat3x3(MatCross, Omega2.Cross(e2b*M(3))))
			+ Mat3x3(MatCross, e2b*MPrime(3)))*Mat3x3(MatCross, e1tota);
		
		if (pModalNode) {
			WM.Add(9 + 1, 3 + 1, MWedge);
			WM.Sub(iOffset2 + 9 + 1, 3 + 1, MWedge);
		}

		MWedge = ((Mat3x3(MatCrossCross, e2tota, Omega2) + Mat3x3(MatCross, Omega1.Cross(e2tota*M(1))))
			+ Mat3x3(MatCross, e2tota*MPrime(1)))*Mat3x3(MatCross, e3b)
			+ ((Mat3x3(MatCrossCross, e3tota, Omega2) + Mat3x3(MatCross, Omega1.Cross(e3tota*M(2))))
			+ Mat3x3(MatCross, e3tota*MPrime(2)))*Mat3x3(MatCross, e1b)
			+ ((Mat3x3(MatCrossCross, e1tota, Omega2) + Mat3x3(MatCross, Omega1.Cross(e1tota*M(3))))
			+ Mat3x3(MatCross, e1tota*MPrime(3)))*Mat3x3(MatCross, e2b);

		if (pModalNode) {
			WM.Sub(9 + 1, iOffset2 + 3 + 1, MWedge);
		}
		WM.Add(iOffset2 + 9 + 1, iOffset2 + 3 + 1, MWedge);

		/* Derivate dell'equilibrio, termini aggiuntivi modali */
		SubMat1.LeftMult(M1Wedge, PHIr);
		if (pModalNode) {
			WM.Add(9 + 1, iRigidOffset + NModes + 1, SubMat1);
		}
		WM.Sub(iOffset2 + 9 + 1, iRigidOffset + NModes + 1, SubMat1);

		Vec3 v1(e2tota.Cross(e3b));
		Vec3 v2(e3tota.Cross(e1b));
		Vec3 v3(e1tota.Cross(e2b));

		/* Error handling: il programma si ferma,
		 * pero' segnala dov'e' l'errore */
		if (v1.Dot() < std::numeric_limits<doublereal>::epsilon()
				|| v2.Dot() < std::numeric_limits<doublereal>::epsilon()
				|| v3.Dot() < std::numeric_limits<doublereal>::epsilon())
		{
			silent_cerr("Modal(" << GetLabel() << "):" 
				<< "warning, first node hinge axis "
				"and second node hinge axis "
				"are (nearly) orthogonal; aborting..."
				<< std::endl);
			throw ErrNullNorm(MBDYN_EXCEPT_ARGS);
		}
		
		MWedge = Mat3x3(v1, v2, v3);

		/* Equilibrio */
		if (pModalNode) {
			WM.Add(3 + 1, iOffset1 + 3 + 1, MWedge);
		}
		WM.Sub(iOffset2 + 3 + 1, iOffset1 + 3 + 1, MWedge);

		SubMat2.RightMult(PHIrT, R1tot.MulTM(MWedge));
		WM.Add(iRigidOffset + 1, iOffset1 + 3 + 1, SubMat2);

		/* Derivate dell'equilibrio */
		if (pModalNode) {
			WM.Add(9 + 1, iOffset1 + 9 + 1, MWedge);
		}
		WM.Sub(iOffset2 + 9 + 1, iOffset1 + 9 + 1, MWedge);

		MWedge = MWedge.Transpose();

		/* Equaz. di vincolo */
		if (pModalNode) {
			WM.Add(iOffset1 + 3 + 1, 3 + 1, MWedge);
		}
		WM.Sub(iOffset1 + 3 + 1, iOffset2 + 3 + 1, MWedge);

		/* Equaz. di vincolo: termine aggiuntivo
		 * dovuto alla flessibilita' */
		Vec3 u1(e2ra.Cross(e3b));
		Vec3 u2(e3ra.Cross(e1b));
		Vec3 u3(e1ra.Cross(e2b));

		M1Wedge = (Mat3x3(u1, u2, u3)).Transpose();
		SubMat1.LeftMult(M1Wedge, PHIr);
		WM.Add(iOffset1 + 3 + 1, iRigidOffset + 1, SubMat1);

		/* Derivate delle equaz. di vincolo */
		if (pModalNode) {
			WM.Add(iOffset1 + 9 + 1, 9 + 1, MWedge);
		}
		WM.Sub(iOffset1 + 9 + 1, iOffset2 + 9 + 1, MWedge);

		v1 = e3b.Cross(e2tota.Cross(Omega1))
			+ e2tota.Cross(Omega2.Cross(e3b));
		v2 = e1b.Cross(e3tota.Cross(Omega1))
			+ e3tota.Cross(Omega2.Cross(e1b));
		v3 = e2b.Cross(e1tota.Cross(Omega1))
			+ e1tota.Cross(Omega2.Cross(e2b));

		MWedge = Mat3x3(v1, v2, v3);

		/* Derivate dell'equilibrio */
		if (pModalNode) {
			WM.Add(9 + 1, iOffset1 + 3 + 1, MWedge);
		}
		WM.Sub(iOffset2 + 9 + 1, iOffset1 + 3 + 1, MWedge);

		/* Derivate delle equaz. di vincolo */
		Omega1 = Omega1 - Omega2;

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
				WM.IncCoef(iOffset1 + 9 + 1, 3 + iCnt, d);

				d = v2.dGet(iCnt);
				WM.IncCoef(iOffset1 + 9 + 2, 3 + iCnt, d);

				d = v3.dGet(iCnt);
				WM.IncCoef(iOffset1 + 9 + 3, 3 + iCnt, d);
			}
		}

		for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
			doublereal d;

			d = v1p.dGet(iCnt);
			WM.IncCoef(iOffset1 + 9 + 1, iOffset2 + 3 + iCnt, d);

			d = v2p.dGet(iCnt);
			WM.IncCoef(iOffset1 + 9 + 2, iOffset2 + 3 + iCnt, d);

			d = v3p.dGet(iCnt);
			WM.IncCoef(iOffset1 + 9 + 3, iOffset2 + 3 + iCnt, d);
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
			WorkVec.PutRowIndex(iCnt, iRigidIndex + iCnt);
		}
	}

	integer iFlexIndex = iGetFirstIndex();

	for (unsigned int iCnt = 1; iCnt <= iGetInitialNumDof(); iCnt++) {
		WorkVec.PutRowIndex(iRigidOffset + iCnt, iFlexIndex + iCnt);
	}

	for (unsigned iStrNodem1 = 0; iStrNodem1 < NStrNodes; iStrNodem1++) {
		integer iNodeFirstPosIndex = 
			pInterfaceNodes[iStrNodem1]->iGetFirstPositionIndex();
		integer iNodeFirstVelIndex = iNodeFirstPosIndex + 6;

		integer iOffset2 = iRigidOffset + iGetInitialNumDof() + 12*iStrNodem1;
		for (unsigned int iCnt = 1; iCnt <= 6; iCnt++) {
			WorkVec.PutRowIndex(iOffset2 + iCnt,
					iNodeFirstPosIndex + iCnt);
			WorkVec.PutRowIndex(iOffset2 + 6 + iCnt,
					iNodeFirstVelIndex + iCnt);
		}
	}

	/* aggiorna le forze modali : K*a, C*aP */
	for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
		a.Put(iCnt, XCurr(iFlexIndex + iCnt));
		b.Put(iCnt, XCurr(iFlexIndex + NModes + iCnt));
	}

	for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
		doublereal temp = 0.;

		for (unsigned int jCnt = 1; jCnt <= NModes; jCnt++) {
			temp += pModalStiff->dGet(iCnt, jCnt)*a.dGet(jCnt)
				+ pModalDamp->dGet(iCnt, jCnt)*b.dGet(jCnt);
		}
		WorkVec.DecCoef(iRigidOffset + iCnt, temp);
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

		integer iOffset1 = iRigidOffset + 2*NModes + 12*iStrNodem1;
		integer iOffset2 = iRigidOffset + iGetInitialNumDof() + 12*iStrNodem1;

		for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
			for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
				PHIt.Put(iCnt, jMode, pPHIt->dGet(iCnt, (jMode - 1)*NStrNodes + iStrNode));
				PHIr.Put(iCnt, jMode, pPHIr->dGet(iCnt, (jMode - 1)*NStrNodes + iStrNode));
			}
		}

		MatNx3 PHItT(NModes), PHIrT(NModes);
		PHItT.Transpose(PHIt);
		PHIrT.Transpose(PHIr);

		Vec3 d1rig(pOffsetFEMNodes->GetVec(iStrNode));
		Vec3 d2(pOffsetMBNodes->GetVec(iStrNode));

		Vec3 d1tot = d1rig + PHIt*a;
		Mat3x3 R1tot = R*Mat3x3(1., PHIr*a);

		Vec3 x2(pInterfaceNodes[iStrNodem1]->GetXCurr());
		Vec3 v2(pInterfaceNodes[iStrNodem1]->GetVCurr());
		Mat3x3 R2(pInterfaceNodes[iStrNodem1]->GetRCurr());
		Vec3 Omega2(pInterfaceNodes[iStrNodem1]->GetWCurr());
		Vec3 F(XCurr, iFlexIndex + 2*NModes + 12*iStrNodem1 + 1);
		Vec3 FPrime(XCurr, iFlexIndex + 2*NModes + 12*iStrNodem1 + 6 + 1);

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
		Vec3 d1Prime(O1Wedged1 + R1PHIt*b);

		/* Equazioni di equilibrio, nodo 1 */
		if (pModalNode) {
			WorkVec.Sub(1, F);
			WorkVec.Add(3 + 1, F.Cross(d1Tmp)); /* F/\d = -d/\F */
		}

		/* Eq. d'equilibrio ai modi */
		for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
			doublereal temp = 0.;
			for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
				temp += PHIt.dGet(iCnt, jMode)*(RT*F).dGet(iCnt);
			}
			WorkVec.DecCoef(iRigidOffset + jMode, temp);
		}

		/* Derivate delle equazioni di equilibrio, nodo 1 */
		if (pModalNode) {
			WorkVec.Sub(6 + 1, FPrime);
			WorkVec.Add(9 + 1, FPrime.Cross(d1Tmp) - d1Prime.Cross(F));
		}

		/* Derivata dell'eq. di equilibrio ai modi */
		MatNx3 PHItTR1T(NModes);
		PHItTR1T.RightMult(PHItT, RT);
		for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
			doublereal temp = 0.;

			for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
				temp += PHItTR1T.dGet(jMode,iCnt)*(Omega1.Cross(F) - FPrime).dGet(iCnt);
			}
			WorkVec.IncCoef(iRigidOffset + NModes + jMode, temp);
		}

		/* Equazioni di equilibrio, nodo 2 */
		WorkVec.Add(iOffset2 + 1, F);
		WorkVec.Add(iOffset2 + 3 + 1, d2Tmp.Cross(F));

		/* Derivate delle equazioni di equilibrio, nodo 2 */
		WorkVec.Add(iOffset2 + 6 + 1, FPrime);
		WorkVec.Add(iOffset2 + 9 + 1, d2Tmp.Cross(FPrime) + O2Wedged2.Cross(F));

		/* Equazione di vincolo */
		WorkVec.Add(iOffset1 + 1, x + d1Tmp - x2 - d2Tmp);

		/* Derivata dell'equazione di vincolo */
		WorkVec.Add(iOffset1 + 6 + 1, v1 + d1Prime - v2 - O2Wedged2);

		/* giunto prismatico */
		Vec3 M(XCurr, iFlexIndex + 2*NModes + 12*iStrNodem1 + 3 + 1);
		Vec3 MPrime(XCurr, iFlexIndex + 2*NModes + 12*iStrNodem1 + 9 + 1);

		Vec3 e1a(R1tot.GetVec(1));
		Vec3 e2a(R1tot.GetVec(2));
		Vec3 e3a(R1tot.GetVec(3));
		Vec3 e1b(R2.GetVec(1));
		Vec3 e2b(R2.GetVec(2));
		Vec3 e3b(R2.GetVec(3));

		Vec3 MTmp(e2a.Cross(e3b*M.dGet(1))
				+ e3a.Cross(e1b*M.dGet(2))
				+ e1a.Cross(e2b*M.dGet(3)));

		/* Equazioni di equilibrio, nodo 1 */
		if (pModalNode) {
			WorkVec.Sub(3 + 1, MTmp);
		}

		/* Equazioni di equilibrio, nodo 2 */
		WorkVec.Add(iOffset2 + 3 + 1, MTmp);

		/* Equazioni di equilibrio, contributo modale */
		Vec3 MTmpRel(R1tot.MulTV(MTmp));
		for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
			WorkVec.DecCoef(iRigidOffset + jMode, PHIrT.GetVec(jMode)*MTmpRel);
		}

		/* eaPrime = w/\ea + R*[(PHIr*b)/\]ia */
		Mat3x3 Tmp = R*Mat3x3(MatCross, PHIr*b);
		Vec3 e1aPrime = Omega1.Cross(e1a) + Tmp.GetVec(1);
		Vec3 e2aPrime = Omega1.Cross(e2a) + Tmp.GetVec(2);
		Vec3 e3aPrime = Omega1.Cross(e3a) + Tmp.GetVec(3);
		Vec3 MTmpPrime(Zero3);
		MTmpPrime =
			(e2a.Cross(Omega2.Cross(e3b)) - e3b.Cross(e2aPrime))*M.dGet(1)
			+ (e3a.Cross(Omega2.Cross(e1b)) - e1b.Cross(e3aPrime))*M.dGet(1)
			+ (e1a.Cross(Omega2.Cross(e2b)) - e2b.Cross(e1aPrime))*M.dGet(1)
			+ e2a.Cross(e3b*MPrime.dGet(1))
			+ e3a.Cross(e1b*MPrime.dGet(2))
			+ e1a.Cross(e2b*MPrime.dGet(3));

		/* Derivate delle equazioni di equilibrio, nodo 1 */
		if (pModalNode) {
			WorkVec.Sub(9 + 1, MTmpPrime);
		}

		/* Derivate delle equazioni di equilibrio, nodo 2 */
		WorkVec.Add(iOffset2 + 9 + 1, MTmpPrime);

		/* Derivate delle equazioni di equilibrio, contributo modale */
		MatNx3 SubMat1(NModes);
		SubMat1.RightMult(PHIrT, R1tot.Transpose());

		/* FIXME: temporary ((PHIr*b).Cross(RT*MTmp)) */
		Vec3 T1 = MTmpPrime - Omega1.Cross(MTmp);
		Vec3 T2 = (PHIr*b).Cross(RT*MTmp);

		for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
			doublereal temp = 0;
			
			for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
				temp += SubMat1.dGet(jMode, iCnt)*T1.dGet(iCnt)
					-PHIrT.dGet(jMode, iCnt)*T2.dGet(iCnt);
			}
			WorkVec.DecCoef(iRigidOffset + NModes + jMode, temp);
		}

		/* Equazioni di vincolo di rotazione */
		WorkVec.DecCoef(iOffset1 + 3 + 1, e3b.Dot(e2a));
		WorkVec.DecCoef(iOffset1 + 3 + 2, e1b.Dot(e3a));
		WorkVec.DecCoef(iOffset1 + 3 + 3, e2b.Dot(e1a));

		/* Derivate delle equazioni di vincolo di rotazione */
		WorkVec.IncCoef(iOffset1 + 9 + 1,
				e3b.Dot(Omega2.Cross(e2a) - e2aPrime));
		WorkVec.IncCoef(iOffset1 + 9 + 2,
				e1b.Dot(Omega2.Cross(e3a) - e3aPrime));
		WorkVec.IncCoef(iOffset1 + 9 + 3,
				e2b.Dot(Omega2.Cross(e1a) - e1aPrime));

	} /* fine equazioni di vincolo */

	return WorkVec;
}

void
Modal::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
	/* inizializza la soluzione e la sua derivata
	 * subito dopo l'assemblaggio iniziale
	 * e prima del calcolo delle derivate */

	int iFlexIndex = iGetFirstIndex();

	for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
		/* modal multipliers */
		X.PutCoef(iFlexIndex + iCnt, a.dGet(iCnt));

		/* derivatives of modal multipliers */
		X.PutCoef(iFlexIndex + NModes + iCnt, b.dGet(iCnt));
		XP.PutCoef(iFlexIndex + iCnt, b.dGet(iCnt));
	}

	if (pModalNode) {
		integer iRigidIndex = pModalNode->iGetFirstIndex();

		X.Put(iRigidIndex + 6 + 1, pModalNode->GetVCurr());
		X.Put(iRigidIndex + 9 + 1, pModalNode->GetWCurr());
		XP.Put(iRigidIndex + 6 + 1, pModalNode->GetXPPCurr());
		XP.Put(iRigidIndex + 9 + 1, pModalNode->GetWPCurr());
	}
}

#if 0
/* Aggiorna dati durante l'iterazione fittizia iniziale */
void
Modal::DerivativesUpdate(const VectorHandler& X,
	const VectorHandler& XP)
{
	// NOTE: identical to SetValue()...
	int iFlexIndex = iGetFirstIndex();

	for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
		/* modal multipliers */
		const_cast<VectorHandler &>(X).PutCoef(iFlexIndex + iCnt, a.dGet(iCnt));

		/* derivatives of modal multipliers */
		const_cast<VectorHandler &>(X).PutCoef(iFlexIndex + NModes + iCnt, b.dGet(iCnt));
		const_cast<VectorHandler &>(XP).PutCoef(iFlexIndex + iCnt, b.dGet(iCnt));
	}
}
#endif

unsigned int
Modal::iGetNumPrivData(void) const
{
	return 3*NModes + 18*NFEMNodes;
}

unsigned int
Modal::iGetPrivDataIdx(const char *s) const
{
	/*
	 * q[n]   | qP[n]   | qPP[n]
	 *
	 * where n is the 1-based mode number
	 *
	 * x[n,i] | xP[n,i] | xPP[n,i]
	 *          w[n,i]  | wP[n,i]
	 *
	 * where n is the FEM node label and i is the index number (1-3)
	 */

	unsigned p = 0;
	char what = s[0];

	/* che cosa e' richiesto? */
	switch (what) {
	case 'q':
		/* modal dofs */
		break;

	case 'x':
		/* structural nodes positions */
		break;

	case 'w':
		/* structural nodes angular velocities */
		p = 1;
		break;

	default:
		/* unknown priv data */
		return 0;
	}

	while ((++s)[0] == 'P') {
		p++;
		if (p > 2) {
			/* no more than two derivative levels allowed */
			return 0;
		}
	}

	if (s[0] != '[') {
		return 0;
	}
	s++;

	/* trova la parentesi chiusa (e il terminatore di stringa) */
	const char *end = std::strchr(s, ']');
	if (end == 0 || end[1] != '\0') {
		return 0;
	}

	if (what == 'q') {
		bool bIdx(false);

		if (s[0] == '#') {
			bIdx = true;
			s++;
		}

		/* buffer per numero (dimensione massima: 32 bit) */
		char buf[sizeof("18446744073709551615,18446744073709551615")];
		size_t		len = end - s;

		ASSERT(len < sizeof(buf));

		strncpy(buf, s, len);
		buf[len] = '\0';

		/* leggi il numero */
		if (buf[0] == '-') {
			return 0;
		}

		errno = 0;
		char *next;
		unsigned long n = strtoul(buf, &next, 10);
		int save_errno = errno;
		end = next;
		if (end == buf) {
			return 0;

		} else if (end[0] != '\0') {
			return 0;

		} else if (save_errno == ERANGE) {
			silent_cerr("Modal(" << GetLabel() << "): "
				"warning, private data mode index "
				<< std::string(buf, end - buf)
				<< " overflows" << std::endl);
			return 0;
		}

		if (bIdx) {
			if (n > NModes) {
				return 0;
			}

		} else {
			std::vector<unsigned int>::const_iterator iv
				= std::find(uModeNumber.begin(), uModeNumber.end(), n);
			if (iv == uModeNumber.end()) {
				return 0;
			}

			n = iv - uModeNumber.begin() + 1;
		}

		return p*NModes + n;
	}

	end = std::strchr(s, ',');
	if (end == NULL) {
		return 0;
	}

	std::string FEMLabel(s, end - s);

	end++;
	if (end[0] == '-') {
		return 0;
	}

	char *next;
	errno = 0;
	unsigned int index = strtoul(end, &next, 10);
	int save_errno = errno;
	if (next == end || next[0] != ']') {
		return 0;

	} else if (save_errno == ERANGE) {
		silent_cerr("Modal(" << GetLabel() << "): "
			"warning, private data FEM node component "
			<< std::string(end, next - end)
			<< " overflows" << std::endl);
		return 0;
	}

	if (index == 0 || index > 3) {
		return 0;
	}

	unsigned int i;
	for (i = 0; i < NFEMNodes; i++) {
		if (IdFEMNodes[i] == FEMLabel) {
			break;
		}
	}

	if (i == NFEMNodes) {
		return 0;
	}

	return 3*NModes + 18*i + (what == 'w' ? 3 : 0) + 6*p + index;
}


doublereal
Modal::dGetPrivData(unsigned int i) const
{
	ASSERT(i > 0 && i < iGetNumPrivData());

	if (i <= NModes) {
		return a.dGet(i);
	}

	i -= NModes;
	if (i <= NModes) {
		return b.dGet(i);
	}

	i -= NModes;
	if (i <= NModes) {
		return bPrime.dGet(i);
	}

	i -= NModes;
	unsigned int n = (i - 1) / 18;
	unsigned int index = (i - 1) % 18 + 1;
	unsigned int p = (index - 1) / 6;
	unsigned int c = (index - 1) % 6 + 1;
	unsigned int w = (c - 1) / 3;
	if (w) {
		c -= 3;
	}

	switch (p) {
	case 0:
		if (w) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		} else {
			if (pModalNode) {
				R = pModalNode->GetRCurr();
				x = pModalNode->GetXCurr();
			}
			Vec3 Xn(pXYZFEMNodes->GetVec(n + 1));
			for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
				integer iOffset = (jMode - 1)*NFEMNodes + n + 1;
				Xn += pModeShapest->GetVec(iOffset)*a(jMode);
			}

			Vec3 X(x + R*Xn);
			return X(c);
		}
		break;

	case 1:
		if (w) {
			if (pModalNode) {
				R = pModalNode->GetRCurr();
			}
			Vec3 Wn(Zero3);
			for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
				integer iOffset = (jMode - 1)*NFEMNodes + n + 1;
				Wn += pModeShapesr->GetVec(iOffset)*b(jMode);
			}

			Vec3 W(R*Wn);
			if (pModalNode) {
				W += pModalNode->GetWCurr();
			}
			return W(c);

		} else {
			if (pModalNode) {
				R = pModalNode->GetRCurr();
			}
			Vec3 Vn(Zero3);
			Vec3 Xn(pXYZFEMNodes->GetVec(n + 1));
			for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
				integer iOffset = (jMode - 1)*NFEMNodes + n + 1;
				Vec3 v = pModeShapest->GetVec(iOffset);
				Vn += v*b(jMode);
				Xn += v*a(jMode);
			}
			Vec3 V(R*Vn);
			if (pModalNode) {
				V += pModalNode->GetWCurr().Cross(R*Xn);
				V += pModalNode->GetVCurr();
			}
			return V(c);
		}
		break;

	case 2:
		if (w) {
			if (pModalNode) {
				R = pModalNode->GetRCurr();
			}
			Vec3 WPn(Zero3);
			Vec3 Wn(Zero3);
			for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
				integer iOffset = (jMode - 1)*NFEMNodes + n + 1;
				Vec3 v = pModeShapesr->GetVec(iOffset);
				WPn += v*bPrime(jMode);
				Wn += v*b(jMode);
			}

			Vec3 WP(R*WPn);
			if (pModalNode) {
				WP += pModalNode->GetWCurr().Cross(R*Wn);
				WP += pModalNode->GetWPCurr();
			}
			return WP(c);

		} else {
			if (pModalNode) {
				R = pModalNode->GetRCurr();
			}
			Vec3 XPPn(Zero3);
			Vec3 XPn(Zero3);
			Vec3 Xn(pXYZFEMNodes->GetVec(n + 1));
			for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
				integer iOffset = (jMode - 1)*NFEMNodes + n + 1;
				Vec3 v = pModeShapest->GetVec(iOffset);
				XPPn += v*bPrime(jMode);
				XPn += v*b(jMode);
				Xn += v*a(jMode);
			}
			Vec3 XPP(R*XPPn);
			if (pModalNode) {
				Vec3 W0 = pModalNode->GetWCurr();
				XPP += W0.Cross(R*(XPn*2.));

				Vec3 X(R*Xn);
				XPP += W0.Cross(W0.Cross(X));
				XPP += pModalNode->GetWPCurr().Cross(X);
				XPP += pModalNode->GetXPPCurr();
			}
			return XPP(c);
		}
		break;
	}

	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

Mat3xN *
Modal::GetCurrFEMNodesPosition(void)
{
	if (pCurrXYZ == NULL) {
		   SAFENEWWITHCONSTRUCTOR(pCurrXYZ, Mat3xN, Mat3xN(NFEMNodes, 0.));
	}

	if (pModalNode) {
		R = pModalNode->GetRCurr();
		x = pModalNode->GetXCurr();
	}

	for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
		for (unsigned int iNode = 1; iNode <= NFEMNodes; iNode++) {
      			for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
				pCurrXYZ->Add(iCnt, iNode,
					pXYZFEMNodes->dGet(iCnt, iNode) + pModeShapest->dGet(iCnt, (jMode - 1)*NFEMNodes + iNode) * a.dGet(jMode) );
			}
		}
	}
	
	pCurrXYZ->LeftMult(R);
	for (unsigned int iNode = 1; iNode <= NFEMNodes; iNode++) {
      		for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
	 		pCurrXYZ->Add(iCnt, iNode, x.dGet(iCnt));
      		}
   	}

      	return pCurrXYZ;
}

Mat3xN *
Modal::GetCurrFEMNodesVelocity(void)
{
	if (pCurrXYZ == NULL) {
		   SAFENEWWITHCONSTRUCTOR(pCurrXYZVel, Mat3xN, Mat3xN(NFEMNodes, 0.));
	}

	Vec3 w(Zero3);
	Vec3 v(Zero3);

	if (pModalNode) {
      		R = pModalNode->GetRCurr();
		w = pModalNode->GetWCurr();
		v = pModalNode->GetVCurr();
	}

	Mat3x3 wWedge(MatCross, w);

      	Mat3xN CurrXYZTmp(NFEMNodes, 0.);
	for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
		for (unsigned int iNode = 1; iNode <= NFEMNodes; iNode++) {
      			for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
				doublereal d = pModeShapest->dGet(iCnt,(jMode - 1)*NFEMNodes + iNode);
				pCurrXYZVel->Add(iCnt, iNode, pXYZFEMNodes->dGet(iCnt,iNode) + d * a.dGet(jMode) );
				CurrXYZTmp.Add(iCnt, iNode,d * b.dGet(jMode) );
			}
		}
	}
	pCurrXYZVel->LeftMult(wWedge*R);
	CurrXYZTmp.LeftMult(R);
	for (unsigned int iNode = 1; iNode <= NFEMNodes; iNode++) {
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
	if (pModalNode != 0) {
		x = pModalNode->GetXCurr();
		R = pModalNode->GetRCurr();
	}

	Vec3 STmp = Inv2;
	if (pInv3 != 0) {
		STmp += Inv3jaj;
	}

	return R*STmp + x*dMass;
}

/* momento d'inerzia */
Mat3x3
Modal::GetJ_int(void) const
{
	if (pModalNode != 0) {
		x = pModalNode->GetXCurr();
		R = pModalNode->GetRCurr();
	}

	Vec3 STmp = Inv2;
	if (pInv3 != 0) {
		STmp += Inv3jaj;
	}
	Vec3 s(R*STmp);

	Mat3x3 JTmp = Inv7;
	if (pInv8 != 0) {
		JTmp += Inv8jaj.Symm2();
		if (pInv9 != 0) {
			JTmp -= Inv9jkajak;
		}
	}

	return R*JTmp.MulMT(R)
		- Mat3x3(MatCrossCross, x, x*dMass)
		- Mat3x3(MatCrossCross, x, s)
		- Mat3x3(MatCrossCross, s, x);
}

// momentum
Vec3
Modal::GetB_int(void) const
{
	// FIXME: gross estimate
	if (pModalNode != 0) {
		return pModalNode->GetVCurr()*dMass;
	}

	return Zero3;
}

// momenta moment
Vec3
Modal::GetG_int(void) const
{
	return Zero3;
}

Joint *
ReadModal(DataManager* pDM,
	MBDynParser& HP,
	const DofOwner* pDO,
	unsigned int uLabel)
{
	/* legge i dati d'ingresso e li passa al costruttore dell'elemento */
	Joint* pEl = NULL;
	
	const ModalNode *pModalNode = 0;
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
		const StructNode* pTmpNode = pDM->pFindNode<const StructNode, Node::STRUCTURAL>(uNode);
		if (pTmpNode == 0) {
			silent_cerr("Modal(" << uLabel << "): "
				"StructuralNode(" << uNode << ") "
				"at line " << HP.GetLineData()
				<< " not defined" << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (pTmpNode->GetStructNodeType() != StructNode::MODAL) {
			silent_cerr("Modal(" << uLabel << "): "
				"illegal type for "
				"StructuralNode(" << uNode << ") "
				"at line " << HP.GetLineData() 
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		pModalNode = dynamic_cast<const ModalNode *>(pTmpNode);

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
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	unsigned int NModes = (unsigned int)tmpNModes;
	unsigned int uLargestMode(0);

	std::vector<unsigned int> uModeNumber(NModes);
	if (HP.IsKeyWord("list")) {
		for (unsigned int iCnt = 0; iCnt < NModes; iCnt++) {
			int n = HP.GetInt();

			if (n <= 0) {
				silent_cerr("Modal(" << uLabel << "): "
					"illegal mode number "
					"for mode #" << iCnt + 1
					<< " at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			std::vector<unsigned int>::const_iterator
				iv = std::find(uModeNumber.begin(),
					uModeNumber.end(), (unsigned int)n);
			if (iv != uModeNumber.end()) {
				silent_cerr("Modal(" << uLabel << "): "
					"mode #" << iCnt + 1
					<< " already defined "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if ((unsigned int)n > uLargestMode) {
				uLargestMode = (unsigned int)n;
			}

			uModeNumber[iCnt] = n;
		}

	} else {
		for (unsigned int iCnt = 0; iCnt < NModes; iCnt++) {
			uModeNumber[iCnt] = iCnt + 1;
		}
		uLargestMode = NModes;
	}

	std::vector<doublereal> InitialValues(NModes, 0.);
	std::vector<doublereal> InitialValuesP(NModes, 0.);
	bool bInitialValues(false);
	if (HP.IsKeyWord("initial" "value")) {
		bInitialValues = true;

		if (HP.IsKeyWord("mode")) {
			int iCnt = 0;
			do {
				int n = HP.GetInt();

				if (n <= 0) {
					silent_cerr("Modal(" << uLabel << "): "
						"illegal mode number "
						"for initial value #" << iCnt + 1
						<< " at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				std::vector<unsigned int>::const_iterator
					iv = std::find(uModeNumber.begin(),
						uModeNumber.end(), (unsigned int)n);
				if (iv == uModeNumber.end()) {
					silent_cerr("Modal(" << uLabel << "): "
						"mode " << n << " not defined "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				unsigned iIdx = iv - uModeNumber.begin();
				InitialValues[iIdx] = HP.GetReal();
				InitialValuesP[iIdx] = HP.GetReal();

				iCnt++;
			} while (HP.IsKeyWord("mode"));

		} else {
			for (unsigned int iCnt = 0; iCnt < NModes; iCnt++) {
				InitialValues[0] = HP.GetReal();
				InitialValuesP[0] = HP.GetReal();
			}
		}
	}

	/* numero di nodi FEM del modello */
	unsigned int NFEMNodes = unsigned(-1);
	if (!HP.IsKeyWord("from" "file")) {
		int tmpNFEMNodes = HP.GetInt();
		if (tmpNFEMNodes <= 0) {
			silent_cerr("Modal(" << uLabel << "): "
				"illegal number of FEM nodes " << tmpNFEMNodes
				<< " at line " << HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		NFEMNodes = unsigned(tmpNFEMNodes);
	}

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
	bool bDampFlag = false;

	if (HP.IsKeyWord("no" "damping")) {
		cdamp = 0.;

	} else if (HP.IsKeyWord("proportional" "damping")) {
		cdamp = HP.GetReal();

	} else if (HP.IsKeyWord("diag" "damping"))  {
		bDampFlag = true;

		if (HP.IsKeyWord("all")) {
			for (unsigned int iCnt = 1; iCnt <= NModes; iCnt ++) {
				cdamp = HP.GetReal();
				DampRatios.Put(iCnt, cdamp);
			}

		} else {
			int iDM = HP.GetInt();
			if (iDM <= 0 || unsigned(iDM) > NModes) {
				silent_cerr("invalid number of damped modes "
					"at line " << HP.GetLineData() <<std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			std::vector<bool> gotit(NModes);
			fill(gotit.begin(), gotit.end(), false);

			for (unsigned int iCnt = 1; iCnt <= unsigned(iDM); iCnt ++) {
				integer iDampedMode =  HP.GetInt();
				if (iDampedMode <= 0 || unsigned(iDampedMode) > NModes) {
					silent_cerr("invalid index " << iDampedMode
						<< " for damped mode #" << iCnt
						<< " of " << iDM
						<< " at line " << HP.GetLineData() <<std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (gotit[iDampedMode - 1]) {
					silent_cerr("damping already provided for mode " << iDampedMode
						<< " (damped mode #" << iCnt
						<< " of " << iDM
						<< ") at line " << HP.GetLineData() <<std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);

				}
				gotit[iDampedMode - 1] = true;
				cdamp = HP.GetReal();
				DampRatios.Put(iDampedMode, cdamp);
			}
		}

	} else {
		silent_cout("Modal(" << uLabel << "): "
				"no damping is assumed at line "
				<< HP.GetLineData() << " (deprecated)"
				<< std::endl);
	}

	DEBUGCOUT("Number of Modes Imported : " << NModes << std::endl);
	DEBUGCOUT("Number of FEM Nodes Imported : " << NFEMNodes << std::endl);
	DEBUGCOUT("Origin of FEM Model : " << X0 << std::endl);
	DEBUGCOUT("Orientation of FEM Model : " << R << std::endl);
	/* DEBUGCOUT("Damping coefficient: "<< cdamp << std::endl); */

	doublereal dMass = 0;              /* massa totale */
	Vec3 STmp(Zero3);                  /* momenti statici  */
	Mat3x3 JTmp(Zero3x3);              /* inerzia totale  */
	std::vector<doublereal> FEMMass;   /* masse nodali   */
	std::vector<Vec3> FEMJ;            /* inerzie nodali (sono diagonali) */

	Mat3xN *pModeShapest = NULL;       /* forme di traslaz. e rotaz. */
	Mat3xN *pModeShapesr = NULL;
	Mat3xN PHIti(NModes, 0.);          /* forme nodo i-esimo: 3*nmodi */
	Mat3xN PHIri(NModes, 0.);
	Mat3xN *pXYZFEMNodes = NULL;       /* punt. alle coordinate nodali */
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

	/* input file */
	const char *s = HP.GetFileName();
	if (s == NULL) {
		silent_cerr("Modal(" << uLabel << "): unable to get "
			"modal data file name at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	std::string sFileFEM(s);

	doublereal dMassThreshold = 0.;
	doublereal dStiffnessThreshold = 0.;
	while (true) {
		if (HP.IsKeyWord("threshold")) {
			dStiffnessThreshold = HP.GetReal();
			if (dStiffnessThreshold < 0.) {
				silent_cerr("Modal(" << uLabel << "): "
					"invalid threshold " << dStiffnessThreshold
					<< " at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			dMassThreshold = dStiffnessThreshold;

		} else if (HP.IsKeyWord("mass" "threshold")) {
			dMassThreshold = HP.GetReal();
			if (dMassThreshold < 0.) {
				silent_cerr("Modal(" << uLabel << "): "
					"invalid mass threshold " << dMassThreshold
					<< " at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

		} else if (HP.IsKeyWord("stiffness" "threshold")) {
			dStiffnessThreshold = HP.GetReal();
			if (dStiffnessThreshold < 0.) {
				silent_cerr("Modal(" << uLabel << "): "
					"invalid stiffness threshold " << dStiffnessThreshold
					<< " at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

		} else {
			break;
		}
	}

	/* loop for binary keywords */
	bool bCreateBinary = false;
	bool bUseBinary = false;
	bool bUpdateBinary = false;
	while (true) {
		if (HP.IsKeyWord("create" "binary")) {
			bCreateBinary = true;

		} else if (HP.IsKeyWord("use" "binary")) {
			bUseBinary = true;

		} else if (HP.IsKeyWord("update" "binary")) {
			bUpdateBinary = true;

		} else {
			break;
		}
	}

	/* stuff for binary FEM file handling */
	std::string	sBinFileFEM;
	struct stat	stFEM, stBIN;
	bool		bReadFEM = true,
			bWriteBIN = false,
			bCheckBIN = false;

	if (stat(sFileFEM.c_str(), &stFEM) == -1) {
		int	save_errno = errno;
		char	*errmsg = strerror(save_errno);

		silent_cerr("Modal(" << uLabel << "): "
				"unable to stat(\"" << sFileFEM << "\") "
				"(" << save_errno << ": " << errmsg << ")"
				<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (bUseBinary || bCreateBinary || bUpdateBinary) {
		sBinFileFEM = sFileFEM + ".bin";

		if (stat(sBinFileFEM.c_str(), &stBIN) == -1) {
			int	save_errno = errno;
			char	*errmsg = strerror(save_errno);

			if (save_errno == ENOENT && (bCreateBinary || bUpdateBinary)) {
				silent_cerr("Modal(" << uLabel << "): "
						"creating binary file "
						"\"" << sBinFileFEM << "\" "
						"from file "
						"\"" << sFileFEM << "\""
						<< std::endl);

			} else {
				silent_cerr("Modal(" << uLabel << "): "
						"unable to stat(\"" << sBinFileFEM << "\") "
						"(" << save_errno << ": " << errmsg << ")"
						<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			
		} else {
			/* can check bin */
			bCheckBIN = true;
		}
	}

	if (bUseBinary && bCheckBIN) {
		/* if timestamp of binary is later than of ASCII use binary */
		if (stBIN.st_mtime > stFEM.st_mtime) {
			bReadFEM = false;

		/* otherwise, if requested, update binary */
		} else if (bUpdateBinary) {
			bWriteBIN = true;
		}

	} else if (bCreateBinary || bUpdateBinary) {
		bWriteBIN = true;
	}

	/* alloca la memoria per le matrici necessarie a memorizzare i dati
	 * relativi al corpo flessibile
	 * nota: devo usare i puntatori perche' altrimenti non si riesce
	 * a passarli al costruttore
	 */
	SAFENEWWITHCONSTRUCTOR(pGenMass,  MatNxN, MatNxN(NModes, 0.));
	SAFENEWWITHCONSTRUCTOR(pGenStiff, MatNxN, MatNxN(NModes, 0.));
	SAFENEWWITHCONSTRUCTOR(pGenDamp,  MatNxN, MatNxN(NModes, 0.));

	SAFENEWWITHCONSTRUCTOR(a,  VecN, VecN(NModes, 0.));
	SAFENEWWITHCONSTRUCTOR(aP, VecN, VecN(NModes, 0.));
	
	/* array contenente le label dei nodi FEM */
	std::vector<std::string> IdFEMNodes;
	std::vector<bool> bActiveModes;

	/* NOTE: increment this each time the binary format changes! */
	unsigned	BinVersion = 2;

	unsigned	currBinVersion;
	char		checkPoint;
	unsigned	NFEMNodesFEM = 0, NModesFEM = 0;

	bool		bBuildInvariants = false;

	/* Note: to keep it readable, we use a base 1 array */
	bool		bRecordGroup[] = {
		false,		// use 1-based array

		false,		// record 1: header
		false,		// record 2: FEM node labels
		false,		// record 3: initial modal amplitudes
		false,		// record 4: initial modal amplitude rates
		false,		// record 5: FEM node X
		false,		// record 6: FEM node Y
		false,		// record 7: FEM node Z
		false,		// record 8: shapes
		false,		// record 9: generalized mass matrix
		false,		// record 10: generalized stiffness matrix
		false,		// record 11: diagonal mass matrix
		false,		// record 12: rigid body inertia

		false		// record -1: end of file
	};

	if (bReadFEM) {
		/* apre il file con i dati del modello FEM */
		std::ifstream fdat(sFileFEM.c_str());
		if (!fdat) {
			silent_cerr("Modal(" << uLabel << "): "
				"unable to open file \"" << sFileFEM << "\""
				"at line " << HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		/* don't leave behind a corrupted .bin file */
		try {

		std::ofstream fbin;
		if (bWriteBIN) {
			fbin.open(sBinFileFEM.c_str());
			if (!fbin) {
				silent_cerr("Modal(" << uLabel << "): "
					"unable to open file \"" << sBinFileFEM << "\""
					"at line " << HP.GetLineData() << std::endl);
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			currBinVersion = BinVersion;
			fbin.write((char *)&currBinVersion, sizeof(currBinVersion));
		}

		/* carica i dati relativi a coordinate nodali, massa, momenti statici
		 * e d'inerzia, massa e rigidezza generalizzate dal file nomefile.
		 */
	
		doublereal	d;
		unsigned int	NRejModes = 0;
		char		str[BUFSIZ];

		while (fdat && !fdat.eof()) {        /* parsing del file */
			fdat.getline(str, sizeof(str));
	
			/* legge il primo blocco (HEADER) */
			if (!strncmp("** RECORD GROUP 1,", str,
				STRLENOF("** RECORD GROUP 1,")))
			{
				if (bRecordGroup[1]) {
					silent_cerr("file=\"" << sFileFEM << "\": "
						"\"RECORD GROUP 1\" already parsed"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

		 		fdat.getline(str, sizeof(str));

				/* ignore REV by now */
				fdat >> str;
	
				/* FEM nodes number */
				fdat >> NFEMNodesFEM;
	
				/* add to modes number */
				fdat >> NModesFEM;
	
				unsigned int i;
				fdat >> i;
				NModesFEM += i;
				fdat >> i;
				NModesFEM += i;

				/* "rejected" modes (subtract from modes number) */
				fdat >> NRejModes;
				NModesFEM -= NRejModes;
		
				/* consistency checks */
				if (NFEMNodes == unsigned(-1)) {
					NFEMNodes = NFEMNodesFEM;

				} else if (NFEMNodes != NFEMNodesFEM) {
					silent_cerr("Modal(" << uLabel << "), "
						"file \"" << sFileFEM << "\": "
						"FEM nodes number " << NFEMNodes
						<< " does not match node number "
						<< NFEMNodesFEM
						<< std::endl);
					throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (NModes > NModesFEM) {
					silent_cerr("Modal(" << uLabel << "), "
						"file '" << sFileFEM << "': "
						"number of requested modes "
						"(" << NModes << ") "
						"exceeds number of available "
						"modes (" << NModesFEM << ")"
						<< std::endl);
					throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);

				}

				if (NModes < NModesFEM) {
					silent_cout("Modal(" << uLabel << "), "
						"file '" << sFileFEM << "': "
						"using " << NModes
						<< " of " << NModesFEM
						<< " available modes"
						<< std::endl);
				}

				if (uLargestMode > NModesFEM) {
					silent_cerr("Modal(" << uLabel << "), "
						"file '" << sFileFEM << "': "
						"largest requested mode "
						"(" << uLargestMode << ") "
						"exceeds available modes "
						"(" << NModesFEM << ")"
						<< std::endl);
					throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (!bActiveModes.empty()) {
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				FEMMass.resize(NFEMNodes);
				FEMJ.resize(NFEMNodes);
				IdFEMNodes.resize(NFEMNodes);

				SAFENEWWITHCONSTRUCTOR(pXYZFEMNodes, Mat3xN, Mat3xN(NFEMNodes, 0.));
				SAFENEWWITHCONSTRUCTOR(pModeShapest, Mat3xN, Mat3xN(NFEMNodes*NModes, 0.));
				SAFENEWWITHCONSTRUCTOR(pModeShapesr, Mat3xN, Mat3xN(NFEMNodes*NModes, 0.));

				bActiveModes.resize(NModesFEM + 1);

				for (unsigned int iCnt = 1; iCnt <= NModesFEM; iCnt++) {
					bActiveModes[iCnt] = false;
				}

				for (unsigned int iCnt = 0; iCnt < NModes; iCnt++) {
					if (uModeNumber[iCnt] > NModesFEM) {
						silent_cerr("Modal(" << uLabel << "): "
							"mode " << uModeNumber[iCnt]
							<< " is not available (max = "
							<< NModesFEM << ")"
							<< std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
					bActiveModes[uModeNumber[iCnt]] = true;
				}

				if (bWriteBIN) {
					checkPoint = 1;
					fbin.write(&checkPoint, sizeof(checkPoint));
					fbin.write((char *)&NFEMNodesFEM, sizeof(NFEMNodesFEM));
					fbin.write((char *)&NModesFEM, sizeof(NModesFEM));
				}

				bRecordGroup[1] = true;

			/* legge il secondo blocco (Id.nodi) */
			} else if (!strncmp("** RECORD GROUP 2,", str,
				STRLENOF("** RECORD GROUP 2,")))
			{
				if (bRecordGroup[2]) {
					silent_cerr("file=\"" << sFileFEM << "\": "
						"\"RECORD GROUP 2\" already parsed"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (!bRecordGroup[1]) {
					silent_cerr("file=\"" << sFileFEM << "\": "
						"\"RECORD GROUP 1\" not parsed yet"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				for (unsigned int iNode = 0; iNode < NFEMNodes; iNode++) {
					// NOTE: this precludes the possibility
					// to have blanks in node labels
					fdat >> IdFEMNodes[iNode];
				}

				if (bWriteBIN) {
					checkPoint = 2;
					fbin.write(&checkPoint, sizeof(checkPoint));
					for (unsigned int iNode = 0; iNode < NFEMNodes; iNode++) {
						size_t len = IdFEMNodes[iNode].size();
						fbin.write((char *)&len, sizeof(size_t));
						fbin.write((char *)IdFEMNodes[iNode].c_str(), len);
					}
				}

				bRecordGroup[2] = true;

			/* deformate iniziali dei modi */
			} else if (!strncmp("** RECORD GROUP 3,", str,
				STRLENOF("** RECORD GROUP 3,")))
			{
				if (bRecordGroup[3]) {
					silent_cerr("file=\"" << sFileFEM << "\": "
						"\"RECORD GROUP 3\" already parsed"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (!bRecordGroup[1]) {
					silent_cerr("file=\"" << sFileFEM << "\": "
						"\"RECORD GROUP 1\" not parsed yet"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				unsigned int iCnt = 1;

				if (bActiveModes.empty()) {
					silent_cerr("Modal(" << uLabel << "): "
						"input file \"" << sFileFEM << "\" "
						"is bogus (RECORD GROUP 3)"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (bWriteBIN) {
					checkPoint = 3;
					fbin.write(&checkPoint, sizeof(checkPoint));
				}

				for (unsigned int jMode = 1; jMode <= NModesFEM; jMode++) {
					fdat >> d;

					if (bWriteBIN) {
						fbin.write((char *)&d, sizeof(d));
					}
					
					if (!bActiveModes[jMode]) {
						continue;
					}
					a->Put(iCnt, d);
					iCnt++;
				}

				bRecordGroup[3] = true;
	
			/* velocita' iniziali dei modi */
			} else if (!strncmp("** RECORD GROUP 4,", str,
				STRLENOF("** RECORD GROUP 4,")))
			{
				if (bRecordGroup[4]) {
					silent_cerr("file=\"" << sFileFEM << "\": "
						"\"RECORD GROUP 4\" already parsed"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (!bRecordGroup[1]) {
					silent_cerr("file=\"" << sFileFEM << "\": "
						"\"RECORD GROUP 1\" not parsed yet"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				unsigned int iCnt = 1;
	
				if (bActiveModes.empty()) {
					silent_cerr("Modal(" << uLabel << "): "
						"input file \"" << sFileFEM << "\" "
						"is bogus (RECORD GROUP 4)"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
	
				if (bWriteBIN) {
					checkPoint = 4;
					fbin.write(&checkPoint, sizeof(checkPoint));
				}

				for (unsigned int jMode = 1; jMode <= NModesFEM; jMode++) {
					fdat >> d;

					if (bWriteBIN) {
						fbin.write((char *)&d, sizeof(d));
					}
					
					if (!bActiveModes[jMode]) {
						continue;
					}
					aP->Put(iCnt, d);
					iCnt++;
				}

				bRecordGroup[4] = true;
	
			/* Coordinate X dei nodi */
			} else if (!strncmp("** RECORD GROUP 5,", str,
				STRLENOF("** RECORD GROUP 5,")))
			{
				if (bRecordGroup[5]) {
					silent_cerr("file=\"" << sFileFEM << "\": "
						"\"RECORD GROUP 5\" already parsed"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (!bRecordGroup[1]) {
					silent_cerr("file=\"" << sFileFEM << "\": "
						"\"RECORD GROUP 1\" not parsed yet"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (bWriteBIN) {
					checkPoint = 5;
					fbin.write(&checkPoint, sizeof(checkPoint));
				}

				for (unsigned int iNode = 1; iNode <= NFEMNodes; iNode++) {
					fdat >> d;
					
					if (bWriteBIN) {
						fbin.write((char *)&d, sizeof(d));
					}

#ifdef MODAL_SCALE_DATA
					d *= scalemodes;
#endif /* MODAL_SCALE_DATA */
					pXYZFEMNodes->Put(1, iNode, d);
				}

				bRecordGroup[5] = true;

			/* Coordinate Y dei nodi*/
			} else if (!strncmp("** RECORD GROUP 6,", str,
				STRLENOF("** RECORD GROUP 6,")))
			{
				if (bRecordGroup[6]) {
					silent_cerr("file=\"" << sFileFEM << "\": "
						"\"RECORD GROUP 6\" already parsed"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (!bRecordGroup[1]) {
					silent_cerr("file=\"" << sFileFEM << "\": "
						"\"RECORD GROUP 1\" not parsed yet"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (bWriteBIN) {
					checkPoint = 6;
					fbin.write(&checkPoint, sizeof(checkPoint));
				}

				for (unsigned int iNode = 1; iNode <= NFEMNodes; iNode++) {
					fdat >> d;

					if (bWriteBIN) {
						fbin.write((char *)&d, sizeof(d));
					}

#ifdef MODAL_SCALE_DATA
					d *= scalemodes;
#endif /* MODAL_SCALE_DATA */
					pXYZFEMNodes->Put(2, iNode, d);
				}

				bRecordGroup[6] = true;
	
			/* Coordinate Z dei nodi*/
			} else if (!strncmp("** RECORD GROUP 7,", str,
				STRLENOF("** RECORD GROUP 7,")))
			{
				if (bRecordGroup[7]) {
					silent_cerr("file=\"" << sFileFEM << "\": "
						"\"RECORD GROUP 7\" already parsed"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (!bRecordGroup[1]) {
					silent_cerr("file=\"" << sFileFEM << "\": "
						"\"RECORD GROUP 1\" not parsed yet"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (bWriteBIN) {
					checkPoint = 7;
					fbin.write(&checkPoint, sizeof(checkPoint));
				}

				for (unsigned int iNode = 1; iNode <= NFEMNodes; iNode++) {
					fdat >> d;

					if (bWriteBIN) {
						fbin.write((char *)&d, sizeof(d));
					}

#ifdef MODAL_SCALE_DATA
					d *= scalemodes;
#endif /* MODAL_SCALE_DATA */
					pXYZFEMNodes->Put(3, iNode, d);
				}

				bRecordGroup[7] = true;
	
			/* Forme modali */
			} else if (!strncmp("** RECORD GROUP 8,", str,
				STRLENOF("** RECORD GROUP 8,")))
			{
				if (bRecordGroup[8]) {
					silent_cerr("file=\"" << sFileFEM << "\": "
						"\"RECORD GROUP 8\" already parsed"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (!bRecordGroup[1]) {
					silent_cerr("file=\"" << sFileFEM << "\": "
						"\"RECORD GROUP 1\" not parsed yet"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (bWriteBIN) {
					checkPoint = 8;
					fbin.write(&checkPoint, sizeof(checkPoint));
				}

				// we assume rejected modes come first
				for (unsigned int jMode = 1; jMode <= NRejModes; jMode++) {
					// header
					fdat.getline(str, sizeof(str));

					silent_cout("Modal(" << uLabel << "): rejecting mode #" << jMode
						<< " (" << str << ")" << std::endl);

					// mode
					for (unsigned int iNode = 1; iNode <= NFEMNodes; iNode++) {
						doublereal n[6];
	
						fdat >> n[0] >> n[1] >> n[2]
							>> n[3] >> n[4] >> n[5];
					}

					/* FIXME: siamo sicuri di avere 
					 * raggiunto '\n'? */
					fdat.getline(str, sizeof(str));
				}
				
				if (bActiveModes.empty()) {
					silent_cerr("Modal(" << uLabel << "): "
						"input file \"" << sFileFEM << "\" "
						"is bogus (RECORD GROUP 8)"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
	
				unsigned int iCnt = 1;
				for (unsigned int jMode = 1; jMode <= NModesFEM; jMode++) {
					fdat.getline(str, sizeof(str));
					if (str[0] != '\0' && strncmp("**", str, STRLENOF("**")) != 0)
					{
						if (NRejModes) {
							silent_cerr("Modal(" << uLabel << "): "
								"input file \"" << sFileFEM << "\""
								"is bogus (\"**\" expected before mode #" << jMode
								<< " of " << NModesFEM << "; "
								"#" << jMode + NRejModes
								<< " of " << NModesFEM + NRejModes
								<< " including rejected)"
								<< std::endl);
						} else {
							silent_cerr("Modal(" << uLabel << "): "
								"input file \"" << sFileFEM << "\""
								"is bogus (\"**\" expected before mode #" << jMode
								<< " of " << NModesFEM << ")"
								<< std::endl);
						}
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
					
					for (unsigned int iNode = 1; iNode <= NFEMNodes; iNode++) {
						doublereal n[6];
	
						fdat >> n[0] >> n[1] >> n[2]
							>> n[3] >> n[4] >> n[5];

						if (bWriteBIN) {
							fbin.write((char *)n, sizeof(n));
						}
	
						if (!bActiveModes[jMode]) {
							continue;
						}

#ifdef MODAL_SCALE_DATA
						n[0] *= scalemodes;
						n[1] *= scalemodes;
						n[2] *= scalemodes;
#endif /* MODAL_SCALE_DATA */
						pModeShapest->PutVec((iCnt - 1)*NFEMNodes + iNode, Vec3(&n[0]));
						pModeShapesr->PutVec((iCnt - 1)*NFEMNodes + iNode, Vec3(&n[3]));
					}

					if (bActiveModes[jMode]) {
						iCnt++;
					}
					fdat.getline(str, sizeof(str));
				}

				bRecordGroup[8] = true;

			/* Matrice di massa  modale */
			} else if (!strncmp("** RECORD GROUP 9,", str,
				STRLENOF("** RECORD GROUP 9,")))
			{
				if (bRecordGroup[9]) {
					silent_cerr("file=\"" << sFileFEM << "\": "
						"\"RECORD GROUP 9\" already parsed"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (!bRecordGroup[1]) {
					silent_cerr("file=\"" << sFileFEM << "\": "
						"\"RECORD GROUP 1\" not parsed yet"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (bActiveModes.empty()) {
					silent_cerr("Modal(" << uLabel << "): "
						"input file \"" << sFileFEM << "\" "
						"is bogus (RECORD GROUP 9)"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (bWriteBIN) {
					checkPoint = 9;
					fbin.write(&checkPoint, sizeof(checkPoint));
				}

				unsigned int iCnt = 1;
				for (unsigned int jMode = 1; jMode <= NModesFEM; jMode++) {
					unsigned int jCnt = 1;
	
					for (unsigned int kMode = 1; kMode <= NModesFEM; kMode++) {
						fdat >> d;

						if (dMassThreshold > 0.0 && fabs(d) < dMassThreshold) {
							d = 0.0;
						}

						if (bWriteBIN) {
							fbin.write((char *)&d, sizeof(d));
						}
						
						if (!bActiveModes[jMode] || !bActiveModes[kMode]) {
							continue;
						}
						pGenMass->Put(iCnt, jCnt, d);
						jCnt++;
					}
	
					if (bActiveModes[jMode]) {
						iCnt++;
					}
				}

				bRecordGroup[9] = true;
	
			/* Matrice di rigidezza  modale */
			} else if (!strncmp("** RECORD GROUP 10,", str,
				STRLENOF("** RECORD GROUP 10,")))
			{
				if (bRecordGroup[10]) {
					silent_cerr("file=\"" << sFileFEM << "\": "
						"\"RECORD GROUP 10\" already parsed"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (!bRecordGroup[1]) {
					silent_cerr("file=\"" << sFileFEM << "\": "
						"\"RECORD GROUP 1\" not parsed yet"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (bActiveModes.empty()) {
					silent_cerr("Modal(" << uLabel << "): "
						"input file \"" << sFileFEM << "\" "
						"is bogus (RECORD GROUP 10)"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
	
				if (bWriteBIN) {
					checkPoint = 10;
					fbin.write(&checkPoint, sizeof(checkPoint));
				}

				unsigned int iCnt = 1;
				for (unsigned int jMode = 1; jMode <= NModesFEM; jMode++) {
					unsigned int jCnt = 1;
	
					for (unsigned int kMode = 1; kMode <= NModesFEM; kMode++) {
						fdat >> d;
						
						if (dStiffnessThreshold > 0.0 && fabs(d) < dStiffnessThreshold) {
							d = 0.0;
						}

						if (bWriteBIN) {
							fbin.write((char *)&d, sizeof(d));
						}
						
						if (!bActiveModes[jMode] || !bActiveModes[kMode]) {
							continue;
						}
						pGenStiff->Put(iCnt, jCnt, d);
						jCnt++;
					}

					if (bActiveModes[jMode]) {
						iCnt++;
					}
				}

				bRecordGroup[10] = true;

			/* Lumped Masses */
			} else if (!strncmp("** RECORD GROUP 11,", str,
				STRLENOF("** RECORD GROUP 11,")))
			{
				if (bRecordGroup[11]) {
					silent_cerr("file=\"" << sFileFEM << "\": "
						"\"RECORD GROUP 11\" already parsed"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (!bRecordGroup[1]) {
					silent_cerr("file=\"" << sFileFEM << "\": "
						"\"RECORD GROUP 1\" not parsed yet"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (bWriteBIN) {
					checkPoint = 11;
					fbin.write(&checkPoint, sizeof(checkPoint));
				}

				for (unsigned int iNode = 1; iNode <= NFEMNodes; iNode++) {
					doublereal tmpmass = 0.;
					for (unsigned int jCnt = 1; jCnt <= 6; jCnt++) {
						fdat >> d;

						if (bWriteBIN) {
							fbin.write((char *)&d, sizeof(d));
						}
						
						switch (jCnt) {
						case 1:
							tmpmass = d;
#ifdef MODAL_SCALE_DATA
							d *= scalemass;
#endif /* MODAL_SCALE_DATA */
							FEMMass[iNode - 1] = d;
							break;

						case 2:
						case 3:
							if (d != tmpmass) {
								silent_cout("Warning: component(" << jCnt << ") "
									"of FEM node(" << iNode << ") mass "
									"differs from component(1)=" << tmpmass << std::endl);
							}
							break;
					
						case 4:
						case 5:
						case 6:
#ifdef MODAL_SCALE_DATA
							d *= scaleinertia;
#endif /* MODAL_SCALE_DATA */
							FEMJ[iNode - 1](jCnt - 3) = d;
							break;
						}
					}
				}

				bBuildInvariants = true;

				bRecordGroup[11] = true;

			/* Rigid body inertia */
			} else if (!strncmp("** RECORD GROUP 12,", str,
				STRLENOF("** RECORD GROUP 12,")))
			{
				if (bRecordGroup[12]) {
					silent_cerr("file=\"" << sFileFEM << "\": "
						"\"RECORD GROUP 12\" already parsed"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (!bRecordGroup[1]) {
					silent_cerr("file=\"" << sFileFEM << "\": "
						"\"RECORD GROUP 1\" not parsed yet"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (bWriteBIN) {
					checkPoint = 12;
					fbin.write(&checkPoint, sizeof(checkPoint));
				}

				fdat >> d;

				if (bWriteBIN) {
					fbin.write((char *)&d, sizeof(d));
				}

				dMass = d;

				// NOTE: the CM location is read and temporarily stored in STmp
				// later, it will be multiplied by the mass
				for (int iRow = 1; iRow <= 3; iRow++) {
					fdat >> d;

					if (bWriteBIN) {
						fbin.write((char *)&d, sizeof(d));
					}

					STmp(iRow) = d;
				}

				for (int iRow = 1; iRow <= 3; iRow++) {
					for (int iCol = 1; iCol <= 3; iCol++) {
						fdat >> d;

						if (bWriteBIN) {
							fbin.write((char *)&d, sizeof(d));
						}

						JTmp(iRow, iCol) = d;
					}
				}

				JTmp -= Mat3x3(MatCrossCross, STmp, STmp*dMass);

				bRecordGroup[12] = true;
				
			} /* fine parser del file */
		}

		fdat.close();
		if (bWriteBIN) {
			checkPoint = -1;
			fbin.write(&checkPoint, sizeof(checkPoint));
			fbin.close();
		}

		/* unlink binary file if create/update failed */
		} catch (...) {
			if (bWriteBIN) {
				// NOTE: might be worth leaving in place
				// for debugging purposes
				(void)unlink(sBinFileFEM.c_str());
			}
			throw;
		}

		if (!bRecordGroup[1]) {
			silent_cerr("Modal(" << uLabel << "): "
				"incomplete input file \"" << sFileFEM << "\", "
				"record group 1 missing "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (!bRecordGroup[2]) {
			silent_cerr("Modal(" << uLabel << "): "
				"incomplete input file \"" << sFileFEM << "\", "
				"record group 2 missing "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (!bRecordGroup[3]) {
			silent_cerr("Modal(" << uLabel << "): "
				"incomplete input file \"" << sFileFEM << "\", "
				"record group 3 missing "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (!bRecordGroup[4]) {
			silent_cerr("Modal(" << uLabel << "): "
				"incomplete input file \"" << sFileFEM << "\", "
				"record group 4 missing "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (!bRecordGroup[5]) {
			silent_cerr("Modal(" << uLabel << "): "
				"incomplete input file \"" << sFileFEM << "\", "
				"record group 5 missing "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (!bRecordGroup[6]) {
			silent_cerr("Modal(" << uLabel << "): "
				"incomplete input file \"" << sFileFEM << "\", "
				"record group 6 missing "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (!bRecordGroup[7]) {
			silent_cerr("Modal(" << uLabel << "): "
				"incomplete input file \"" << sFileFEM << "\", "
				"record group 7 missing "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (!bRecordGroup[8]) {
			silent_cerr("Modal(" << uLabel << "): "
				"incomplete input file \"" << sFileFEM << "\", "
				"record group 8 missing "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (!bRecordGroup[9]) {
			silent_cerr("Modal(" << uLabel << "): "
				"incomplete input file \"" << sFileFEM << "\", "
				"record group 9 missing "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (!bRecordGroup[10]) {
			silent_cerr("Modal(" << uLabel << "): "
				"incomplete input file \"" << sFileFEM << "\", "
				"record group 10 missing "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else {
		std::ifstream fbin(sBinFileFEM.c_str());
		if (!fbin) {
			silent_cerr("Modal(" << uLabel << "): "
				"unable to open file \"" << sBinFileFEM << "\""
				"at line " << HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		silent_cout("Modal(" << uLabel << "): "
			"reading flexible body data from file "
			"\"" << sBinFileFEM << "\"" << std::endl);

		/* check version */
		// NOTE: any time the format of the binary file changes,
		// the BinVersion should be incremented.
		// When reading the binary file, the version of the file
		// is compared with that of the code.  If they do not match:
		// if the file is newer than the code, abort.
		// if the file is older than the code, try to deal with it
		// as much as possible
		//
		// History:
		//
		// 1: first implementation
		//
		// 2: July 2008, at LaMSID
		//	- string FEM labels
		//	- 3 & 4 optional
		//	- 11 & 12 optional and no longer mutually exclusive
		//	- -1 as the last checkpoint
		//	- version 1 tolerated
		fbin.read((char *)&currBinVersion, sizeof(currBinVersion));
		if (currBinVersion > BinVersion) {
			silent_cerr("Modal(" << uLabel << "): "
					"file \"" << sBinFileFEM << "\" "
					"version " << currBinVersion << " unsupported" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		} else if (currBinVersion < BinVersion) {
			silent_cout("Modal(" << uLabel << "): "
					"file \"" << sBinFileFEM << "\" "
					"version " << currBinVersion << " is outdated; "
					"trying to read anyway..." << std::endl);
		}

		pedantic_cout("Modal(" << uLabel << "): "
			"binary version " << currBinVersion << std::endl);

		// NOTE: the binary file is expected to be in order;
		// however, if generated from a non-ordered textual file
		// it will be out of order...  perhaps also the textual
		// file should always be ordered?

		/* legge il primo blocco (HEADER) */
		fbin.read(&checkPoint, sizeof(checkPoint));
		if (checkPoint != 1) {
			silent_cerr("Modal(" << uLabel << "): "
					"file \"" << sBinFileFEM << "\" "
					"looks broken (expecting checkpoint 1)"
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		pedantic_cout("Modal(" << uLabel << "): "
			"reading block " << int(checkPoint) << std::endl);

		fbin.read((char *)&NFEMNodesFEM, sizeof(NFEMNodesFEM));
		fbin.read((char *)&NModesFEM, sizeof(NModesFEM));

		/* consistency checks */
		if (NFEMNodes == unsigned(-1)) {
			NFEMNodes = NFEMNodesFEM;

		} else if (NFEMNodes != NFEMNodesFEM) {
			silent_cerr("Modal(" << uLabel << "), "
				"file \"" << sFileFEM << "\": "
				"FEM nodes number " << NFEMNodes
				<< " does not match node number "
				<< NFEMNodesFEM
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (NModes != NModesFEM) {
			silent_cout("Modal(" << uLabel
					<< "), file '" << sFileFEM
					<< "': using " << NModes
					<< " of " << NModesFEM
					<< " modes" << std::endl);
		}

		if (!bActiveModes.empty()) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		FEMMass.resize(NFEMNodes);
		FEMJ.resize(NFEMNodes);
		IdFEMNodes.resize(NFEMNodes);

		SAFENEWWITHCONSTRUCTOR(pXYZFEMNodes, Mat3xN, Mat3xN(NFEMNodes, 0.));
		SAFENEWWITHCONSTRUCTOR(pModeShapest, Mat3xN, Mat3xN(NFEMNodes*NModes, 0.));
		SAFENEWWITHCONSTRUCTOR(pModeShapesr, Mat3xN, Mat3xN(NFEMNodes*NModes, 0.));

		bActiveModes.resize(NModesFEM + 1);

		for (unsigned int iCnt = 1; iCnt <= NModesFEM; iCnt++) {
			bActiveModes[iCnt] = false;
		}

		for (unsigned int iCnt = 0; iCnt < NModes; iCnt++) {
			if (uModeNumber[iCnt] > NModesFEM) {
				silent_cerr("Modal(" << uLabel << "): "
					"mode " << uModeNumber[iCnt]
					<< " is not available (max = "
					<< NModesFEM << ")"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			bActiveModes[uModeNumber[iCnt]] = true;
		}

		bRecordGroup[1] = true;

		/* legge il secondo blocco (Id.nodi) */
		fbin.read(&checkPoint, sizeof(checkPoint));
		if (checkPoint != 2) {
			silent_cerr("Modal(" << uLabel << "): "
					"file \"" << sBinFileFEM << "\" "
					"looks broken (expecting checkpoint 2)"
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		pedantic_cout("Modal(" << uLabel << "): "
			"reading block " << int(checkPoint) << std::endl);

		if (currBinVersion == 1) {
			// NOTE: accept unsigned int FEM labels (old binary format)
			for (unsigned int iNode = 0; iNode < NFEMNodes; iNode++) {
				unsigned uFEMLabel;
				fbin.read((char *)&uFEMLabel, sizeof(uFEMLabel));
				std::ostringstream os;
				os << uFEMLabel;
				IdFEMNodes[iNode] = os.str();
			}

		} else {
			for (unsigned int iNode = 0; iNode < NFEMNodes; iNode++) {
				size_t len;
				fbin.read((char *)&len, sizeof(size_t));
				ASSERT(len > 0);
				IdFEMNodes[iNode].resize(len);
				fbin.read((char *)IdFEMNodes[iNode].c_str(), len);
			}
		}

		bRecordGroup[2] = true;

		/* deformate iniziali dei modi */
		fbin.read(&checkPoint, sizeof(checkPoint));
		if (checkPoint != 3) {
			if (currBinVersion == 1) {
				silent_cerr("Modal(" << uLabel << "): "
						"file \"" << sBinFileFEM << "\" "
						"looks broken (expecting checkPoint 3)"
						<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			silent_cerr("Modal(" << uLabel << "): "
					"warning: file \"" << sBinFileFEM << "\" "
					"does not contain modal coordinate initial values"
					<< std::endl);

		} else {
			pedantic_cout("Modal(" << uLabel << "): "
				"reading block " << int(checkPoint) << std::endl);

			for (unsigned int iCnt = 1, jMode = 1; jMode <= NModesFEM; jMode++) {
				doublereal	d;
	
				fbin.read((char *)&d, sizeof(d));
					
				if (!bActiveModes[jMode]) {
					continue;
				}
	
				a->Put(iCnt, d);
				iCnt++;
			}
	
			bRecordGroup[3] = true;

			fbin.read(&checkPoint, sizeof(checkPoint));
		}

		/* velocita' iniziali dei modi */
		if (checkPoint != 4) {
			if (currBinVersion == 1) {
				silent_cerr("Modal(" << uLabel << "): "
						"file \"" << sBinFileFEM << "\" "
						"looks broken (expecting checkPoint 4)"
						<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			silent_cerr("Modal(" << uLabel << "): "
					"warning: file \"" << sBinFileFEM << "\" "
					"does not contain modal coordinate derivative initial values"
					<< std::endl);

		} else {
			pedantic_cout("Modal(" << uLabel << "): "
				"reading block " << int(checkPoint) << std::endl);

			for (unsigned int iCnt = 1, jMode = 1; jMode <= NModesFEM; jMode++) {
				doublereal	d;

				fbin.read((char *)&d, sizeof(d));
				
				if (!bActiveModes[jMode]) {
					continue;
				}

				aP->Put(iCnt, d);
				iCnt++;
			}
	
			bRecordGroup[4] = true;

			fbin.read(&checkPoint, sizeof(checkPoint));
		}

		/* Coordinate X dei nodi */
		if (checkPoint != 5) {
			silent_cerr("Modal(" << uLabel << "): "
					"file \"" << sBinFileFEM << "\" "
					"looks broken (expecting checkpoint 5)"
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		pedantic_cout("Modal(" << uLabel << "): "
			"reading block " << int(checkPoint) << std::endl);

		for (unsigned int iNode = 1; iNode <= NFEMNodes; iNode++) {
			doublereal	d;

			fbin.read((char *)&d, sizeof(d));

#ifdef MODAL_SCALE_DATA
			d *= scalemodes;
#endif /* MODAL_SCALE_DATA */

			pXYZFEMNodes->Put(1, iNode, d);
		}

		bRecordGroup[5] = true;

		/* Coordinate Y dei nodi*/
		fbin.read(&checkPoint, sizeof(checkPoint));
		if (checkPoint != 6) {
			silent_cerr("Modal(" << uLabel << "): "
					"file \"" << sBinFileFEM << "\" "
					"looks broken (expecting checkpoint 6)"
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		pedantic_cout("Modal(" << uLabel << "): "
			"reading block " << int(checkPoint) << std::endl);

		for (unsigned int iNode = 1; iNode <= NFEMNodes; iNode++) {
			doublereal	d;

			fbin.read((char *)&d, sizeof(d));

#ifdef MODAL_SCALE_DATA
			d *= scalemodes;
#endif /* MODAL_SCALE_DATA */

			pXYZFEMNodes->Put(2, iNode, d);
		}

		bRecordGroup[6] = true;

		/* Coordinate Z dei nodi*/
		fbin.read(&checkPoint, sizeof(checkPoint));
		if (checkPoint != 7) {
			silent_cerr("Modal(" << uLabel << "): "
					"file \"" << sBinFileFEM << "\" "
					"looks broken (expecting checkpoint 7)"
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		pedantic_cout("Modal(" << uLabel << "): "
			"reading block " << int(checkPoint) << std::endl);

		for (unsigned int iNode = 1; iNode <= NFEMNodes; iNode++) {
			doublereal	d;

			fbin.read((char *)&d, sizeof(d));

#ifdef MODAL_SCALE_DATA
			d *= scalemodes;
#endif /* MODAL_SCALE_DATA */

			pXYZFEMNodes->Put(3, iNode, d);
		}

		bRecordGroup[7] = true;

		/* Forme modali */
		fbin.read(&checkPoint, sizeof(checkPoint));
		if (checkPoint != 8) {
			silent_cerr("Modal(" << uLabel << "): "
					"file \"" << sBinFileFEM << "\" "
					"looks broken (expecting checkpoint 8)"
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		pedantic_cout("Modal(" << uLabel << "): "
			"reading block " << int(checkPoint) << std::endl);

		for (unsigned int iCnt = 1, jMode = 1; jMode <= NModesFEM; jMode++) {
			for (unsigned int iNode = 1; iNode <= NFEMNodes; iNode++) {
				doublereal n[6];
	
				fbin.read((char *)n, sizeof(n));
	
				if (!bActiveModes[jMode]) {
					continue;
				}

#ifdef MODAL_SCALE_DATA
				n[0] *= scalemodes;
				n[1] *= scalemodes;
				n[2] *= scalemodes;
#endif /* MODAL_SCALE_DATA */
				pModeShapest->PutVec((iCnt - 1)*NFEMNodes + iNode, Vec3(&n[0]));
				pModeShapesr->PutVec((iCnt - 1)*NFEMNodes + iNode, Vec3(&n[3]));
			}

			if (bActiveModes[jMode]) {
				iCnt++;
			}
		}

		bRecordGroup[8] = true;

		/* Matrice di massa modale */
		fbin.read(&checkPoint, sizeof(checkPoint));
		if (checkPoint != 9) {
			silent_cerr("Modal(" << uLabel << "): "
					"file \"" << sBinFileFEM << "\" "
					"looks broken (expecting checkpoint 9)"
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		pedantic_cout("Modal(" << uLabel << "): "
			"reading block " << int(checkPoint) << std::endl);

		for (unsigned int iCnt = 1, jMode = 1; jMode <= NModesFEM; jMode++) {
			unsigned int jCnt = 1;
	
			for (unsigned int kMode = 1; kMode <= NModesFEM; kMode++) {
				doublereal	d;

				fbin.read((char *)&d, sizeof(d));

				if (!bActiveModes[jMode] || !bActiveModes[kMode]) {
					continue;
				}
				pGenMass->Put(iCnt, jCnt, d);
				jCnt++;
			}
	
			if (bActiveModes[jMode]) {
				iCnt++;
			}
		}

		bRecordGroup[9] = true;

		/* Matrice di rigidezza  modale */
		fbin.read(&checkPoint, sizeof(checkPoint));
		if (checkPoint != 10) {
			silent_cerr("Modal(" << uLabel << "): "
					"file \"" << sBinFileFEM << "\" "
					"looks broken (expecting checkpoint 10)"
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		pedantic_cout("Modal(" << uLabel << "): "
			"reading block " << int(checkPoint) << std::endl);

		for (unsigned int iCnt = 1, jMode = 1; jMode <= NModesFEM; jMode++) {
			unsigned int jCnt = 1;
	
			for (unsigned int kMode = 1; kMode <= NModesFEM; kMode++) {
				doublereal	d;

				fbin.read((char *)&d, sizeof(d));

				if (!bActiveModes[jMode] || !bActiveModes[kMode]) {
					continue;
				}
				pGenStiff->Put(iCnt, jCnt, d);
				jCnt++;
			}
	
			if (bActiveModes[jMode]) {
				iCnt++;
			}
		}

		bRecordGroup[10] = true;

		/* Lumped Masses */
		fbin.read(&checkPoint, sizeof(checkPoint));
		switch (checkPoint) {
		case 11:
			pedantic_cout("Modal(" << uLabel << "): "
				"reading block " << int(checkPoint) << std::endl);

			for (unsigned int iNode = 1; iNode <= NFEMNodes; iNode++) {
				for (unsigned int jCnt = 1; jCnt <= 6; jCnt++) {
					doublereal	d;

					fbin.read((char *)&d, sizeof(d));
						
					switch (jCnt) {
					case 1:
#ifdef MODAL_SCALE_DATA
						d *= scalemass;
#endif /* MODAL_SCALE_DATA */
						FEMMass[iNode - 1] = d;
						break;

					case 4:
					case 5:
					case 6:
#ifdef MODAL_SCALE_DATA
						d *= scaleinertia;
#endif /* MODAL_SCALE_DATA */
						FEMJ[iNode - 1](jCnt - 3) = d;
						break;
					}
				}
			}

			bBuildInvariants = true;

			bRecordGroup[11] = true;

			// read next checkpoint
			fbin.read(&checkPoint, sizeof(checkPoint));
			break;

		case 12: {
			pedantic_cout("Modal(" << uLabel << "): "
				"reading block " << int(checkPoint) << std::endl);

			doublereal      d;
			fbin.read((char *)&d, sizeof(d));
			dMass = d;

			// NOTE: the CM location is read and temporarily stored in STmp
			// later, it will be multiplied by the mass
			for (int iRow = 1; iRow <= 3; iRow++) {
				fbin.read((char *)&d, sizeof(d));

				STmp(iRow) = d;
			}

			for (int iRow = 1; iRow <= 3; iRow++) {
				for (int iCol = 1; iCol <= 3; iCol++) {
					fbin.read((char *)&d, sizeof(d));

					JTmp(iRow, iCol) = d;
				}
			}

			JTmp -= Mat3x3(MatCrossCross, STmp, STmp*dMass);

			bRecordGroup[12] = true;

			// read next checkpoint
			fbin.read(&checkPoint, sizeof(checkPoint));

			} break;
		}

		if (currBinVersion != 1) {
			if (checkPoint != -1) {
				silent_cerr("Modal(" << uLabel << "): "
					"file \"" << sBinFileFEM << "\" "
					"looks broken (expecting final checkpoint)"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		pedantic_cout("Modal(" << uLabel << "): "
			"reading block " << int(checkPoint) << std::endl);

		fbin.close();
	}

	/* lettura dati di vincolo:
	 * l'utente specifica la label del nodo FEM e del nodo rigido
	 * d'interfaccia.
	 * L'orientamento del nodo FEM e' quello del nodo modale, la
	 * posizione e' la somma di quella modale e di quella FEM   */

	/* puntatori ai nodi multibody */
	const StructNode** pInterfaceNodes = NULL;
	/* array contenente le label dei nodi d'interfaccia */
	std::vector<std::string> IntFEMNodes;
	/* array contenente le forme modali dei nodi d'interfaccia */
	Mat3xN* pPHItStrNode = NULL;
	Mat3xN* pPHIrStrNode = NULL;

	/* traslazione origine delle coordinate */
	Vec3 Origin(Zero3);
	bool bOrigin(false);

	if (HP.IsKeyWord("origin" "node")) {
		/* numero di nodi FEM del modello */
		std::string FEMOriginNode;
		if (HP.IsStringWithDelims()) {
			FEMOriginNode = HP.GetStringWithDelims();

		} else {
			pedantic_cerr("Modal(" << uLabel << "): "
				"origin node expected as string with delimiters"
				<< std::endl);
			std::ostringstream os;
			os << HP.GetInt();
			FEMOriginNode = os.str();
		}

		unsigned int iNode;
		for (iNode = 0; iNode < NFEMNodes; iNode++) {
			if (FEMOriginNode == IdFEMNodes[iNode]) {
				break;
			}
		}

		if (iNode == NFEMNodes) {
			silent_cerr("Modal(" << uLabel << "): "
				"FEM node \"" << FEMOriginNode << "\""
				<< " at line " << HP.GetLineData()
				<< " not defined " << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		iNode++;
		Origin = pXYZFEMNodes->GetVec(iNode);

		pedantic_cout("Modal(" << uLabel << "): origin x={" << Origin << "}" << std::endl);
		bOrigin = true;

	} else if (HP.IsKeyWord("origin" "position")) {
		Origin = HP.GetPosAbs(AbsRefFrame);
		bOrigin = true;
	}

	if (bOrigin) {
		for (unsigned int iStrNode = 1; iStrNode <= NFEMNodes; iStrNode++) {
			pXYZFEMNodes->SubVec(iStrNode, Origin);
		}

		if (!bBuildInvariants) {
			STmp -= Origin;
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
	IntFEMNodes.resize(NStrNodes);

	for (unsigned int iStrNode = 1; iStrNode <= NStrNodes; iStrNode++) {
		/* nodo collegato 1 (e' il nodo FEM) */
		std::string Node1;
		if (HP.IsStringWithDelims()) {
			Node1 = HP.GetStringWithDelims();
			// NOTE: we should check whether a label includes
			// whitespace.  In any case, it wouldn't match
			// any of the labels read, so it is pointless,
			// but could help debugging...
			if (Node1.find(' ') != std::string::npos ) {
				silent_cout("Modal(" << uLabel << "): "
					"FEM node \"" << Node1 << "\""
					<< " at line " << HP.GetLineData()
					<< " contains a blank" << std::endl);
			}

		} else {
			std::ostringstream os;
			os << HP.GetInt();
			Node1 = os.str();
		}

		DEBUGCOUT("Linked to FEM Node " << Node1 << std::endl);

		unsigned int iNode;
		/* verifica di esistenza del nodo 1*/
		for (iNode = 0; iNode < NFEMNodes; iNode++) {
			if (Node1 == IdFEMNodes[iNode]) {
				break;
			}
		}

		if (iNode == NFEMNodes) {
			silent_cerr("Modal(" << uLabel << "): "
				"FEM node \"" << Node1 << "\""
				<< " at line " << HP.GetLineData()
				<< " not defined " << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
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
				pXYZFEMNodes->GetVec(iNodeCurr));

		/* salva le forme modali del nodo d'interfaccia
		 * nell'array pPHIStrNode */
		for (unsigned int jMode = 0; jMode < NModes; jMode++) {
			pPHItStrNode->PutVec(jMode*NStrNodes + iStrNode,
					pModeShapest->GetVec(jMode*NFEMNodes + iNodeCurr));
			pPHIrStrNode->PutVec(jMode*NStrNodes + iStrNode,
					pModeShapesr->GetVec(jMode*NFEMNodes + iNodeCurr));
		}

		/* nodo collegato 2 (e' il nodo multibody) */
		unsigned int uNode2 = (unsigned int)HP.GetInt();
		DEBUGCOUT("Linked to Multi-Body Node " << uNode2 << std::endl);

		/* verifica di esistenza del nodo 2 */
		pInterfaceNodes[iStrNode - 1] = pDM->pFindNode<const StructNode, Node::STRUCTURAL>(uNode2);
		if (pInterfaceNodes[iStrNode - 1] == NULL) {
			silent_cerr("Modal(" << uLabel << "): "
				"StructuralNode(" << uNode2 << ") "
				"at line " << HP.GetLineData()
				<< " not defined" << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		/* offset del nodo Multi-Body */
		ReferenceFrame RF = ReferenceFrame(pInterfaceNodes[iStrNode - 1]);
		Vec3 d2(HP.GetPosRel(RF));
		Mat3x3 R2(Eye3);
		if (HP.IsKeyWord("hinge") || HP.IsKeyWord("orientation")) {
			R2 = HP.GetRotRel(RF);
		}

		pOffsetMBNodes->PutVec(iStrNode, d2);
		pRotMBNodes->PutMat3x3(3*(iStrNode - 1) + 1, R2);

		DEBUGCOUT("Multibody Node reference frame d2: " << std::endl
				<< d2 << std::endl);

		/* salva le label dei nodi vincolati nell'array IntNodes;
		 * puo' servire per il restart? */
		IntFEMNodes[iStrNode - 1] = Node1;

		Vec3 xMB(pInterfaceNodes[iStrNode - 1]->GetXCurr());
		pedantic_cout("Interface node " << iStrNode << ":" << std::endl
				<< "        MB node " << uNode2 << " x={" << xMB << "}" << std::endl);
		Vec3 xFEMRel(pXYZFEMNodes->GetVec(iNodeCurr));
		Vec3 xFEM(X0 + R*xFEMRel);
		pedantic_cout("        FEM node \"" << Node1 << "\" x={" << xFEM << "} "
			"xrel={" << xFEMRel << "}" << std::endl);
		pedantic_cout("        offset={" << xFEM - xMB << "}" << std::endl);
	}  /* fine ciclo sui nodi d'interfaccia */

	/* fine ciclo caricamento dati */

	/*
	 * calcola gli invarianti d'inerzia (massa, momenti statici
	 * e d'inerzia, termini di accoppiamento nel SdR locale)
	 */

	/* NOTE: right now, either they are all used (except for Inv9)
	 * or none is used, and the global inertia of the body is expected */
	if (bBuildInvariants) {
		doublereal dMassInv = 0.;
		Vec3 STmpInv(Zero3);
		Mat3x3 JTmpInv(Zero3x3);

		MatNxN GenMass(NModes, 0.);

		/* TODO: select what terms are used */

		SAFENEWWITHCONSTRUCTOR(pInv3, Mat3xN, Mat3xN(NModes, 0.));
		SAFENEWWITHCONSTRUCTOR(pInv4, Mat3xN, Mat3xN(NModes, 0.));

		/* NOTE: we assume that Inv5 is used only if Inv4 is used as well */
		/* Inv5 e' un 3xMxM */
		SAFENEWWITHCONSTRUCTOR(pInv5, Mat3xN, Mat3xN(NModes*NModes, 0.));

		/* Inv8 e' un 3x3xM */
		SAFENEWWITHCONSTRUCTOR(pInv8, Mat3xN, Mat3xN(3*NModes, 0.));

		/* NOTE: we assume that Inv9 is used only if Inv8 is used as well */
		/* by default: no */
		if (HP.IsKeyWord("use" "invariant" "9")) {
			/* Inv9 e' un 3x3xMxM */
			SAFENEWWITHCONSTRUCTOR(pInv9, Mat3xN, Mat3xN(3*NModes*NModes, 0.));
		}
		/* Inv10 e' un 3x3xM */
		SAFENEWWITHCONSTRUCTOR(pInv10,Mat3xN, Mat3xN(3*NModes, 0.));
		SAFENEWWITHCONSTRUCTOR(pInv11,Mat3xN, Mat3xN(NModes, 0.));

		/* inizio ciclo scansione nodi */
		for (unsigned int iNode = 1; iNode <= NFEMNodes; iNode++) {
			doublereal mi = FEMMass[iNode - 1];
	
			/* massa totale (Inv 1) */
			dMassInv += mi;
	
			/* posizione nodi FEM */
			Vec3 ui = pXYZFEMNodes->GetVec(iNode);
	
			Mat3x3 uiWedge(MatCross, ui);
			Mat3x3 JiNodeTmp(Zero3x3);
	
			JiNodeTmp(1, 1) = FEMJ[iNode - 1](1);
			JiNodeTmp(2, 2) = FEMJ[iNode - 1](2);
			JiNodeTmp(3, 3) = FEMJ[iNode - 1](3);
	
			JTmpInv += JiNodeTmp - Mat3x3(MatCrossCross, ui, ui*mi);
			STmpInv += ui*mi;
	
			/* estrae le forme modali del nodo i-esimo */
			for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
				unsigned int iOffset = (jMode - 1)*NFEMNodes + iNode;
	
				PHIti.PutVec(jMode, pModeShapest->GetVec(iOffset));
				PHIri.PutVec(jMode, pModeShapesr->GetVec(iOffset));
			}
	
			/* TODO: only build what is required */
	
			Mat3xN Inv3Tmp(NModes, 0.);
			Mat3xN Inv4Tmp(NModes, 0.);
			Mat3xN Inv4JTmp(NModes, 0.);
			Inv3Tmp.Copy(PHIti);
	
			/* Inv3 = mi*PHIti,      i = 1,...nnodi */
			Inv3Tmp *= mi;

			/* Inv4 = mi*ui/\*PHIti + Ji*PHIri, i = 1,...nnodi */
			Inv4Tmp.LeftMult(uiWedge*mi, PHIti);
			Inv4JTmp.LeftMult(JiNodeTmp, PHIri);
			Inv4Tmp += Inv4JTmp;
			*pInv3 += Inv3Tmp;
			*pInv4 += Inv4Tmp;
			*pInv11 += Inv4JTmp;
	
			/* inizio ciclo scansione modi */
			for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
				Vec3 PHItij = PHIti.GetVec(jMode);
				Vec3 PHIrij = PHIri.GetVec(jMode);
	
				Mat3x3 PHItijvett_mi(MatCross, PHItij*mi);
				Mat3xN Inv5jTmp(NModes, 0);
	
				/* Inv5 = mi*PHItij/\*PHIti,
				 * i = 1,...nnodi, j = 1,...nmodi */
				Inv5jTmp.LeftMult(PHItijvett_mi, PHIti);
				for (unsigned int kMode = 1; kMode <= NModes; kMode++)  {
					pInv5->AddVec((jMode - 1)*NModes + kMode,
							Inv5jTmp.GetVec(kMode));

					/* compute the modal mass matrix
					 * using the FEM inertia and the
					 * mode shapes */
					GenMass(jMode, kMode) += (PHItij*PHIti.GetVec(kMode))*mi
						+ PHIrij*(JiNodeTmp*PHIri.GetVec(kMode));
				}
	
				/* Inv8 = -mi*ui/\*PHItij/\,
				 * i = 1,...nnodi, j = 1,...nmodi */
				Mat3x3 Inv8jTmp = -uiWedge*PHItijvett_mi;
				pInv8->AddMat3x3((jMode - 1)*3 + 1, Inv8jTmp);
	
				/* Inv9 = mi*PHItij/\*PHItik/\,
				 * i = 1,...nnodi, j, k = 1...nmodi */
				if (pInv9 != 0) {
					for (unsigned int kMode = 1; kMode <= NModes; kMode++) {
						Mat3x3 PHItikvett(MatCross, PHIti.GetVec(kMode));
						Mat3x3 Inv9jkTmp = PHItijvett_mi*PHItikvett;
	
						pInv9->AddMat3x3((jMode - 1)*3*NModes + (kMode - 1)*3 + 1, Inv9jkTmp);
					}
				}
	
				/* Inv10 = [PHIrij/\][J0i],
				 * i = 1,...nnodi, j = 1,...nmodi */
				Mat3x3 Inv10jTmp = PHIrij.Cross(JiNodeTmp);
				pInv10->AddMat3x3((jMode - 1)*3 + 1, Inv10jTmp);
			} /*  fine ciclo scansione modi */
		} /* fine ciclo scansione nodi */

		if (bRecordGroup[12]) {
			Mat3x3 DJ = JTmp - JTmpInv;
			pedantic_cerr("  Rigid-body mass: "
				"input - computed" << std::endl
				<< "    " << dMass - dMassInv << std::endl);
			pedantic_cerr("  Rigid-body CM location: "
				"input - computed" << std::endl
				<< "    " << STmp - STmpInv/dMassInv << std::endl);
			pedantic_cerr("  Rigid-body inertia: "
				"input - computed" << std::endl
				<< "    " << DJ.GetVec(1) << std::endl
				<< "    " << DJ.GetVec(2) << std::endl
				<< "    " << DJ.GetVec(3) << std::endl);
		}

		/*
		 * TODO: check modal mass
		 */
		pedantic_cerr("  Generalized mass: input - computed" << std:: endl);
		for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
			pedantic_cerr("   ");
			for (unsigned int kMode = 1; kMode <= NModes; kMode++) {
				pedantic_cerr(" " << pGenMass->dGet(jMode, kMode) - GenMass(jMode, kMode));
			}
			pedantic_cerr(std::endl);
		}

		// if record 12 was read, leave its data in place
		if (!bRecordGroup[12]) {
			dMass = dMassInv;
			STmp = STmpInv;
			JTmp = JTmpInv;
		}

	}

	// if record 12 was read, fix static moment
	if (bRecordGroup[12]) {
		/* left over when reading XCG */
		STmp *= dMass;
	}

	/*
	 * TODO: Check rank of modal stiffness matrix
	 */
	
	/*
	 * costruisce la matrice di smorzamento:
	 * il termine diagonale i-esimo e' pari a
	 * cii = 2*cdampi*(ki*mi)^.5
	 */
	for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
		doublereal d = sqrt(pGenStiff->dGet(iCnt, iCnt)
				*pGenMass->dGet(iCnt, iCnt));

		if (bDampFlag) {
			pGenDamp->Put(iCnt, iCnt, 2.*DampRatios.dGet(iCnt)*d);

		} else {
			pGenDamp->Put(iCnt, iCnt, 2.*cdamp*d);
		}
	}

#if 1 // def DEBUG
	if (pedantic_output) {
		std::ostream &out = std::cout;

		out << "  Total Mass: " << dMass << std::endl;
		out << "  Inertia Matrix: " << std::endl
			<< "    " << JTmp << std::endl;
		out << "  Static Moment Vector: " << STmp << std::endl;

		out << "  Generalized Stiffness: " << std::endl;
		for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
			out << "   ";
			for (unsigned int jCnt = 1; jCnt <= NModes; jCnt++) {
				out << " " << pGenStiff->dGet(iCnt, jCnt);
			}
			out << std::endl;
		}

		out << "  Generalized Mass: " << std::endl;
		for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
			out << "   ";
			for (unsigned int jCnt = 1; jCnt <= NModes; jCnt++) {
				out << " " << pGenMass->dGet(iCnt, jCnt);
			}
			out << std::endl;
		}

		out << "  Generalized Damping: " << std::endl;
		for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
			out << "   ";
			for (unsigned int jCnt = 1; jCnt <= NModes; jCnt++) {
				out << " " << pGenDamp->dGet(iCnt, jCnt);
			}
			out << std::endl;
		}

		if (pInv3 != 0) {
			out << "  Inv3: " << std::endl;
			for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
				out << "   ";
				for (unsigned int jCnt = 1; jCnt <= NModes; jCnt++) {
					out << " " << pInv3->dGet(iCnt, jCnt);
				}
				out << std::endl;
			}
		} else {
			out << "  Inv3: unused" << std::endl;
		}

		if (pInv4 != 0) {
			out << "  Inv4: " << std::endl;
			for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
				out << "   ";
				for (unsigned int jCnt = 1; jCnt <= NModes; jCnt++) {
					out << " " << pInv4->dGet(iCnt, jCnt);
				}
				out << std::endl;
			}
		} else {
			out << "  Inv4: unused" << std::endl;
		}

		if (pInv5 != 0) {
			for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
				out << "  Inv5j(j=" << jMode << "): " << std::endl;
				for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
					out << "   ";
					for (unsigned int kMode = 1; kMode <= NModes; kMode++) {
						out << " " << pInv5->dGet(iCnt, (jMode - 1)*NModes + kMode);
					}
					out << std::endl;
				}
			}
		} else {
			out << "  Inv5: unused" << std::endl;
		}

		if (pInv8 != 0) {
			for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
				out << "  Inv8j(j=" << jMode << "): " << std::endl;
				for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
					out << "   ";
					for (unsigned int jCnt = 1; jCnt <= 3; jCnt++) {
						out << " " << pInv8->dGet(iCnt, (jMode - 1)*3 + jCnt);
					}
					out << std::endl;
				}
			}
		} else {
			out << "  Inv8: unused" << std::endl;
		}

		if (pInv9 != 0) {
			for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
				for (unsigned int kMode = 1; kMode <= NModes; kMode++) {
					out << "  Inv9jk(j=" << jMode << ",k=" << kMode << "): " << std::endl;
					for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
						out << "   ";
						for (unsigned int jCnt = 1; jCnt <= 3; jCnt++) {
							out << " " << pInv9->dGet(iCnt, (jMode - 1)*3*NModes + (kMode - 1)*3 + jCnt);
						}
						out << std::endl;
					}
				}
			}
		} else {
			out << "  Inv9: unused" << std::endl;
		}
	}
#endif /* DEBUG */

	if (HP.IsStringWithDelims()) {
		const char *sTmp = HP.GetFileName();
		silent_cout("Modal(" << uLabel << "): warning, the syntax changed "
			"since 1.2.7; the output now occurs to a common \".mod\" file, "
			"the per-element file \"" << sTmp << "\" is no longer required, "
			"and will actually be ignored." << std::endl);
	}
	flag fOut = pDM->fReadOutput(HP, Elem::JOINT);

	if (bInitialValues) {
		for (unsigned int iCnt = 0; iCnt < NModes; iCnt++) {
			a->Put(iCnt + 1, InitialValues[iCnt]);
			aP->Put(iCnt + 1, InitialValuesP[iCnt]);
		}
	}

	SAFENEWWITHCONSTRUCTOR(pEl,
			Modal,
			Modal(uLabel,
				pModalNode,
				X0,
				R,
				pDO,
				NModes,
				NStrNodes,
				NFEMNodes,
				dMass,
				STmp,
				JTmp,
				uModeNumber,
				pGenMass,
				pGenStiff,
				pGenDamp,
				IdFEMNodes,
				IntFEMNodes,
				pXYZFEMNodes,
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
				pDM,
				HP,
				fOut));

	if (fOut) {
		pDM->OutputOpen(OutputHandler::MODAL);
	}

	SAFEDELETE(a);
	SAFEDELETE(aP);

	return pEl;
}

