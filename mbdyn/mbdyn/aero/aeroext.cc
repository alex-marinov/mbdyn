/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#ifdef USE_EXTERNAL
#include "aeroext.h"
#include "myassert.h"
#include "except.h"
#include "dataman.h"
#include "external.h"

/* le label delle comunicazioni vengono costruite in questo modo:
 *
 *	blockID*10   +  1   label e numero nodi aggiuntivi
 *			2   posizioni
 *			3   velocita'
 *			4   forze
 *			5   momenti
 *
 */

/* Aerodynamic External - begin */
AerodynamicExternal::AerodynamicExternal(unsigned int uLabel,
	const DofOwner *pDO,
	const DriveCaller* pDC,
	int NN,
	const StructNode** ppN,
	const doublereal* pRefL,
	MPI::Intercomm* pIC,
	flag fOut,
	bool VF,
	bool MF)
: Elem(uLabel, fOut),
AerodynamicElem(uLabel, pDO, fOut),
DriveOwner(pDC),
pdBuffer(NULL),
pdBufferVel(NULL),
NodeN(NN),
ppNode(ppN),
pRefLength(pRefL),
OffN(3),
pOffsetVectors(NULL),
pInterfComm(pIC),
VelFlag(VF),
MomFlag(MF),
SentFlag(false),
pForce(NULL),
pMoms(NULL),
pLabList(NULL)
{
	/* costruisce la mat 3xN degli offset */
	SAFENEWWITHCONSTRUCTOR(pOffsetVectors, Mat3xN, Mat3xN(OffN, 0.));
	for (int i=1; i<=3; i++) pOffsetVectors->Put(i,i,1.);
	ConstructAndInitialize();
}

AerodynamicExternal::AerodynamicExternal(unsigned int uLabel,
	const DofOwner* pDO,
	const DriveCaller* pDC,
	int NN,
	const StructNode** ppN,
	const doublereal* pRefL,
	MPI::Intercomm* pIC,
	int	ON,
	Mat3xN* OV,
	flag fOut,
	bool VF,
	bool MF)
: Elem(uLabel, fOut),
AerodynamicElem(uLabel, pDO, fOut),
DriveOwner(pDC),
pdBuffer(NULL),
pdBufferVel(NULL),
NodeN(NN),
ppNode(ppN),
pRefLength(pRefL),
OffN(ON),
pOffsetVectors(OV),
pInterfComm(pIC),
VelFlag(VF),
MomFlag(MF),
SentFlag(false),
pForce(NULL),
pMoms(NULL),
pLabList(NULL)
{
	ConstructAndInitialize();
}

AerodynamicExternal::~AerodynamicExternal(void)
{
	if (pdBuffer != NULL) SAFEDELETE(pdBuffer);
	if (pdBufferVel != NULL) SAFEDELETE(pdBufferVel);
	if (ppNode != NULL) SAFEDELETEARR(ppNode);
	if (pRefLength != NULL) SAFEDELETEARR(pRefLength);
	if (pForce != NULL) SAFEDELETEARR(pForce);
	if (pMoms != NULL) SAFEDELETEARR(pMoms);
	if (pLabList != NULL) SAFEDELETEARR(pLabList);
}

void
AerodynamicExternal::ConstructAndInitialize(void)
{
	/* dimensiona l'array per il buffer */
	DEBUGCOUTFNAME("AerodynamicExternal::ConstructAndInitialize");
	 if (bToBeOutput()) {
		SAFENEWWITHCONSTRUCTOR(pForce, Mat3xN, Mat3xN(NodeN+1,0.));
		SAFENEWWITHCONSTRUCTOR(pMoms, Mat3xN, Mat3xN(NodeN+1,0.));
	}
	SAFENEWWITHCONSTRUCTOR(pdBuffer, MyVectorHandler,
			MyVectorHandler(8+3*NodeN*(OffN+1)));
	if (VelFlag || MomFlag) SAFENEWWITHCONSTRUCTOR(pdBufferVel, MyVectorHandler,
			MyVectorHandler(3*NodeN*(OffN+1)));

	/* invio delle posizioni iniziali dei nodi per il calcolo della matrice di interpolazione */
	SAFENEWARR(pLabList, unsigned int, 2*NodeN);
	Mat3xN MatTmp(OffN,0.);
	Vec3 X;
	for (int i=0; i < NodeN; i++) {
		pLabList[i] = ppNode[i]->GetLabel();
		pLabList[i+NodeN] = 1+OffN;
		X = ppNode[i]->GetXCurr();
		pdBuffer->Put(i*(OffN+1)*3+1, X);
		for (int j=1; j <= OffN; j++) {
			MatTmp.LeftMult(ppNode[i]->GetRCurr(), *pOffsetVectors);
			MatTmp *= pRefLength[i];
			pdBuffer->Put(i*(OffN+1)*3+j*3+1, X + MatTmp.GetVec(j));
		}
	}
	/* il rank del processo con cui si comunica e' zero perche' il comuncatire del
	codice di interfaccia e' sempre fatto da una macchina sola */
	pInterfComm->Send(pLabList, 2*NodeN, MPI::UNSIGNED, 0,(this->GetLabel())*10+1);
	pInterfComm->Send(pdBuffer->pdGetVec(), 3*NodeN*(OffN+1), MPI::DOUBLE, 0,(this->GetLabel())*10+2);
}


/* invia posizione e velocita' predetti */
void
AerodynamicExternal::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	DEBUGCOUTFNAME("AerodynamicExternal::AfterPredict");
	//std::cout << "AerodynamicExternal::AfterPredict" << std::endl;
	Send(X, XP);
}

/* invia posizione e velocita' predetti */
void
AerodynamicExternal::Update(const VectorHandler& X, const VectorHandler& XP)
{
	DEBUGCOUTFNAME("AerodynamicExternal::Update");
	//std::cout << "AerodynamicExternal::Update" << std::endl;
	Send(X, XP);
	//std::cout << "done" << std::endl;
}

void
AerodynamicExternal::Send(const VectorHandler& X, const VectorHandler& XP)
{
	if (!SentFlag) {
		//std::cout << "AerodynamicExternal::Send " <<  SentFlag << std::endl;
		Vec3 Vinf(0.);
		Vec3 Xc(0.), V, W, Xt;
		Mat3x3 R;
		doublereal rho, c, p, T;
		GetAirProps(Xc, rho, c, p, T);
		fGetAirVelocity(Vinf,Xc);
		/* current time */
		pdBuffer->PutCoef(1, DriveOwner::dGet());
		pdBuffer->Put(2, Vinf);
		pdBuffer->PutCoef(5, rho);
		pdBuffer->PutCoef(6,c);
		pdBuffer->PutCoef(7,p);
		pdBuffer->PutCoef(8,T);

		Mat3xN MatTmp(OffN,0.);
		for (int i=0; i < NodeN; i++) {
			Xc = ppNode[i]->GetXCurr();
			pdBuffer->Put(8+i*(OffN+1)*3+1, Xc);
                        //std::cout << "Block " << (this->GetLabel()) << " Node " << i << " Pos " << Xc << std::endl;
			if (VelFlag) {
				V = ppNode[i]->GetVCurr();
                                W = ppNode[i]->GetWCurr();
				pdBufferVel->Put(i*(OffN+1)*3+1, V);
			}
			R = ppNode[i]->GetRCurr();
                        //std::cout <<  "Vel " << V  << " R " << R << " W " << W << std::endl;
			for (int j=1; j <= OffN; j++) {
				MatTmp.LeftMult(R, *pOffsetVectors);
				MatTmp *= pRefLength[i];
				Xt = MatTmp.GetVec(j);
				pdBuffer->Put(8+i*(OffN+1)*3+j*3+1, Xc + Xt);
				if (VelFlag) pdBufferVel->Put(i*(OffN+1)*3+j*3+1, V + W.Cross(Xt));
				//std::cout << " OffVec " << j << ": Pos " << Xt << " Vel " << V + W.Cross(Xt) << std::endl;
			}
		}
		pInterfComm->Send(pdBuffer->pdGetVec(), 8+3*NodeN*(OffN+1), MPI::DOUBLE, 0,(this->GetLabel())*10+2);
		if (VelFlag) pInterfComm->Send(pdBufferVel->pdGetVec(), 3*NodeN*(OffN+1), MPI::DOUBLE, 0,(this->GetLabel())*10+3);
		SentFlag = true;
	}
}

SubVectorHandler&
AerodynamicExternal::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("AerodynamicExternal::AssRes");
	/* verifica il completamento dei send */
	//std::cout << "AerodynamicExternal::AssRes " <<  SentFlag << std::endl;

	Send(XCurr, XPrimeCurr);

	WorkVec.ResizeReset(6*NodeN);

	for (int i=0; i < NodeN; i++) {
		integer iFirstIndex = ppNode[i]->iGetFirstMomentumIndex();
		for (int iCnt = 1; iCnt <= 6; iCnt++) {
			WorkVec.PutRowIndex(i*6 + iCnt, iFirstIndex + iCnt);
		}
	}

	pInterfComm->Recv(pdBuffer->pdGetVec(), 3*NodeN*(OffN+1), MPI::DOUBLE, 0,(this->GetLabel())*10+4);
	if(MomFlag) pInterfComm->Recv(pdBufferVel->pdGetVec(), 3*NodeN*(OffN+1), MPI::DOUBLE, 0,(this->GetLabel())*10+5);
	/* calcola le forze su ciascun nodo */
	/* le forze e i momenti sono gia' espressi nel sistema di riferimento global */
	if (bToBeOutput()) {
		pForce->Reset();
		pMoms->Reset();
	}
	for (int i=0; i < NodeN; i++) {
		Vec3 F(0.);
		Vec3 Fo(0.);
		Vec3 M(0.);
		F[0] +=  pdBuffer->operator()(i*(OffN+1)*3+1);
		F[1] +=  pdBuffer->operator()(i*(OffN+1)*3+2);
		F[2] +=  pdBuffer->operator()(i*(OffN+1)*3+3);
		if (MomFlag) {
			M[0] +=  pdBufferVel->operator()(i*(OffN+1)*3+1);
			M[1] +=  pdBufferVel->operator()(i*(OffN+1)*3+2);
			M[2] +=  pdBufferVel->operator()(i*(OffN+1)*3+3);
		}
		for (int j=1; j <= OffN; j++) {
			Fo[0] =  pdBuffer->operator()(i*(OffN+1)*3+j*3+1);
			Fo[1] =  pdBuffer->operator()(i*(OffN+1)*3+j*3+2);
			Fo[2] =  pdBuffer->operator()(i*(OffN+1)*3+j*3+3);
			F += Fo;
			if (MomFlag) {
				M[0] +=  pdBufferVel->operator()(i*(OffN+1)*3+j*3+1);
				M[1] +=  pdBufferVel->operator()(i*(OffN+1)*3+j*3+2);
				M[2] +=  pdBufferVel->operator()(i*(OffN+1)*3+j*3+3);
			}
			M -= Fo.Cross( ppNode[i]->GetRCurr() * pOffsetVectors->GetVec(j) * pRefLength[i]);
		}
		if (bToBeOutput()) {
			pForce->AddVec(1,F);
			pMoms->AddVec(1,M);
			pForce->PutVec(i+2,F);
			pMoms->PutVec(i+2,M);
		}
		WorkVec.Add(i*6+1, F);
		WorkVec.Add(i*6+4, M);
	}
	//std::cout << "FT " <<  pForce->GetVec(1) << std::endl;
	//std::cout << "MT " <<  pMoms->GetVec(1) << std::endl;

	SentFlag = false;
	return WorkVec;
}

void
AerodynamicExternal::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {
		integer lab = -1*GetLabel();
		std::ostream& out = OH.Externals()
                        << std::setw(8) << lab
			<< " " << pForce->GetVec(1)
			<< " " << pMoms->GetVec(1) << std::endl;
		if (NodeN > 1) {
			for (int i=1; i <= NodeN; i++) {
				out << std::setw(8) << pLabList[i-1]
					<< " " << pForce->GetVec(i+1)
					<< " " << pMoms->GetVec(i+1) << std::endl;
			}
		}
	}
}

Elem *
ReadAerodynamicExternal(DataManager* pDM, MBDynParser& HP,
	const DofOwner *pDO, unsigned int uLabel)
{
	/* legge i dati d'ingresso e li passa al costruttore dell'elemento */
	AerodynamicExternal* pEl = NULL;

	unsigned int NodeN = HP.GetInt();

	std::cout << "Aerodynamic External Block made of  " << NodeN << " Nodes." << std::endl;

	const StructNode** ppNodeList = NULL;
	SAFENEWARR(ppNodeList, const StructNode*, NodeN);

	for (unsigned int iN = 0; iN < NodeN; iN++) {
		unsigned int NodeL = HP.GetInt();
		const StructNode* pTmpNode = pDM->pFindNode<const StructNode, Node::STRUCTURAL>(NodeL);
		if (pTmpNode == 0) {
			std::cerr << "Aerodynamic External(" << uLabel
				<< "): structural node " << NodeL
				<< " at line " << HP.GetLineData()
				<< " not defined" << std::endl;
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		} else {
			ppNodeList[iN] = pTmpNode;
		}
	}
	doublereal* pRefLenght = NULL;
	SAFENEWARR(pRefLenght, doublereal, NodeN);
	if (HP.IsKeyWord("list")) {
		for (unsigned int iN=0; iN < NodeN; iN++) {
			pRefLenght[iN] = HP.GetReal();
		}
	} else {
		pRefLenght[0] = HP.GetReal();
		for (unsigned int iN=1; iN < NodeN; iN++) {
			pRefLenght[iN] = pRefLenght[iN-1];
		}
	}

	int OffN = 0;
	Mat3xN* pOffVec = NULL;
	if (HP.IsKeyWord("offset")) {
		OffN = HP.GetInt();
		SAFENEWWITHCONSTRUCTOR(pOffVec, Mat3xN, Mat3xN(OffN));
		for (int iOf = 1; iOf <= OffN; iOf++) {
			for (int iCnt = 1; iCnt <= 3; iCnt++) {
				pOffVec->Put(iCnt, iOf, HP.GetReal());
			}
		}
	}


	int comm = 0;
	if (HP.IsKeyWord("comm")) {
		comm = HP.GetInt();
	}
	std::list<MPI::Intercomm>::iterator iComm =  InterfaceComms.begin();
	std::list<MPI::Intercomm>::iterator iComm_end =  InterfaceComms.end();
	if (iComm == iComm_end) {
                   std::cerr << "Aerodynamic External(" << uLabel
                                << "): there are no external communicators "
                                << "Aborting." << std::endl;
                   throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	for (int i = 0; i++ < comm; ++iComm) {
		if (iComm == iComm_end) {
                   std::cerr << "Aerodynamic External(" << uLabel
                                << "): external communicator " << comm
				<< "does not exist. Last communicator is " << i-1
                                << ". Aborting." << std::endl;
                   throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
	MPI::Intercomm* pInterC = &(*(iComm));

	bool VelFlag = false;
	bool MomFlag = false;

	if (HP.IsKeyWord("velocity")) VelFlag = true;
	if (HP.IsKeyWord("moment")) MomFlag = true;

	flag fOut = pDM->fReadOutput(HP, Elem::AERODYNAMIC);
	DriveCaller* pDC = NULL;
	SAFENEWWITHCONSTRUCTOR(pDC, TimeDriveCaller, TimeDriveCaller(pDM->pGetDrvHdl()));

	if (OffN) {
		SAFENEWWITHCONSTRUCTOR(pEl,
			AerodynamicExternal,
			AerodynamicExternal(uLabel,
				pDO,
				pDC,
				NodeN,
				ppNodeList,
				pRefLenght,
				pInterC,
				OffN,
				pOffVec,
				fOut,
				VelFlag,
				MomFlag));
	} else {
		SAFENEWWITHCONSTRUCTOR(pEl,
			AerodynamicExternal,
			AerodynamicExternal(uLabel,
				pDO,
				pDC,
				NodeN,
				ppNodeList,
				pRefLenght,
				pInterC,
				fOut,
				VelFlag,
				MomFlag));
	}

	return pEl;
}


/* Aerodynamic External - end */


/* Aerodynamic External Modal - begin */

AerodynamicExternalModal::AerodynamicExternalModal(unsigned int uLabel,
	const DofOwner *pDO,
	const DriveCaller* pDC,
	Modal* pM,
	MPI::Intercomm* IC,
	flag fOut,
	bool VF,
	bool MF)
:Elem(uLabel, fOut),
AerodynamicElem(uLabel, pDO, fOut),
DriveOwner(pDC),
pdBuffer(NULL),
pdBufferVel(NULL),
pModal(pM),
ModalNodes(pM->uGetNFemNodes()),
pInterfComm(IC),
VelFlag(VF),
MomFlag(MF),
SentFlag(false),
pForce(NULL)
{
	SAFENEWWITHCONSTRUCTOR(pdBuffer, MyVectorHandler,
		MyVectorHandler(8+3*ModalNodes));

	if (bToBeOutput()) {
		integer NModes = pModal->uGetNModes();
		SAFENEWWITHCONSTRUCTOR(pForce, MyVectorHandler,
			MyVectorHandler(NModes));
	}

	if (VelFlag || MomFlag) SAFENEWWITHCONSTRUCTOR(pdBufferVel, MyVectorHandler,
			MyVectorHandler(3*ModalNodes));

	unsigned int*  pLabList = NULL;
	SAFENEWARR(pLabList, unsigned int, 2);
	pLabList[0] = pModal->GetLabel();
	pLabList[1] = ModalNodes;
	/* il rank del processo con cui si comunica e' zero perche' il comuncatore del
	codice di interfaccia e' sempre fatto da una macchina sola */
	pInterfComm->Send(pLabList, 2, MPI::UNSIGNED, 0,(this->GetLabel())*10+1);

	Mat3xN* pNodePos = pModal->GetCurrFemNodesPosition();
	for (int i=0; i < ModalNodes; i++) {
		pdBuffer->Put(i*3, pNodePos->GetVec(i+1));
	}
	pInterfComm->Send(pdBuffer->pdGetVec(), 3*ModalNodes, MPI::DOUBLE, 0,(this->GetLabel())*10+2);
}

AerodynamicExternalModal::~AerodynamicExternalModal(void)
{
	if (pdBuffer != NULL) SAFEDELETE(pdBuffer);
	if (pdBufferVel != NULL) SAFEDELETE(pdBufferVel);
	if (pForce != NULL) SAFEDELETE(pForce);
}

/* invia posizione e velocita' predetti */
void
AerodynamicExternalModal::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	DEBUGCOUTFNAME("AerodynamicExternalModal::AfterPredict");
	Send(X, XP);
}

/* invia posizione e velocita' predetti */
void AerodynamicExternalModal::Update(const VectorHandler& X  ,
					const VectorHandler&  XP  )
{
	DEBUGCOUTFNAME("AerodynamicExternalModal::Update");
	Send(X, XP);
}

void AerodynamicExternalModal::Send(const VectorHandler& X  ,
					const VectorHandler&  XP  )
{

	Vec3 Vinf(0.);
	doublereal rho, c, p, T;
	GetAirProps(Vinf, rho, c, p, T);
	/* current time */
	pdBuffer->PutCoef(1, DriveOwner::dGet());
	pdBuffer->Put(2, Vinf);
	pdBuffer->PutCoef(5, rho);
	pdBuffer->PutCoef(6,c);
	pdBuffer->PutCoef(7,p);
	pdBuffer->PutCoef(8,T);

	Mat3xN* pNodePos = pModal->GetCurrFemNodesPosition();
	Mat3xN* pNodeVel = NULL;
	if (VelFlag) pNodeVel = pModal->GetCurrFemNodesVelocity();
	for (int i=0; i < ModalNodes; i++) {
		pdBuffer->Put(8+i*3, pNodePos->GetVec(i+1));
		if (VelFlag) pdBufferVel->Put(i*3, pNodeVel->GetVec(i+1));
	}
	pInterfComm->Send(pdBuffer->pdGetVec(), 8+3*ModalNodes, MPI::DOUBLE, 0,(this->GetLabel())*10+2);
	if (VelFlag) pInterfComm->Send(pdBufferVel->pdGetVec(), 3*ModalNodes, MPI::DOUBLE, 0,(this->GetLabel())*10+3);
	SentFlag = true;

}


SubVectorHandler&
AerodynamicExternalModal::AssRes(SubVectorHandler& WorkVec,
				doublereal dCoef,
				const VectorHandler& XCurr,
				const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("AerodynamicExternalModal::AssRes");
	if (!SentFlag) {
		Send(XCurr, XPrimeCurr);
	}
	integer NModes = pModal->uGetNModes();
	WorkVec.ResizeReset(NModes);


	integer iFirstIndex = pModal->iGetFirstIndex() + NModes;
	for (int iCnt = 0; iCnt < NModes; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstIndex + iCnt);
	}

	pInterfComm->Recv(pdBuffer->pdGetVec(), 3*ModalNodes, MPI::DOUBLE, 0, GetLabel()*10 + 4);
	if (MomFlag) {
		pInterfComm->Recv(pdBufferVel->pdGetVec(), 3*ModalNodes, MPI::DOUBLE, 0, GetLabel()*10 + 5);
	}

	/* calcola le forze su ciascun nodo
	 * le forze e i momenti sono espressi nel sistema di riferimento global
	 * vanno quindi riportati nel sistema di riferimento del nodo modal e
	 * poi proiettati sulla base modale
	 */

	Vec3 F(0.);
	Vec3 M(0.);
	Mat3x3 RT(pModal->pGetModalNode()->GetRCurr().Transpose());

#warning "FIXME: review"
	for (int iNode = 0; iNode < ModalNodes; iNode++) {
		F[0] +=  pdBuffer->operator()(iNode*ModalNodes*3 + 1);
		F[1] +=  pdBuffer->operator()(iNode*ModalNodes*3 + 2);
		F[2] +=  pdBuffer->operator()(iNode*ModalNodes*3 + 3);
		if (MomFlag) {
			M[0] +=  pdBufferVel->operator()(iNode*ModalNodes*3 + 1);
			M[1] +=  pdBufferVel->operator()(iNode*ModalNodes*3 + 2);
			M[2] +=  pdBufferVel->operator()(iNode*ModalNodes*3 + 3);
		}
		for (int iMode = 0; iMode < NModes; iMode++) {
			WorkVec.Add(iMode + 1, RT*F*(pModal->pGetPHIt().GetVec(iMode*ModalNodes + iNode + 1)));
			if (MomFlag) {
				WorkVec.Add(iMode + 1, RT*M*(pModal->pGetPHIr().GetVec(iMode*ModalNodes + iNode + 1)));
			}
		}
	}

	if (bToBeOutput()) {
		*pForce = WorkVec;
	}

	SentFlag = false;
	return WorkVec;
}


void
AerodynamicExternalModal::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {
		std::ostream& out = OH.Externals()
                        << std::setw(8) << -1*GetLabel();
		int size = pForce->iGetSize();
		for (int i = 0; i < size; i++) {
			out << std::setw(8)
				<< " " << pForce->operator()(i+1) << std::endl;
		}
	}
}

Elem *
ReadAerodynamicExternalModal(DataManager* pDM, MBDynParser& HP,
	const DofOwner *pDO, unsigned int uLabel)
{
	/* legge i dati d'ingresso e li passa al costruttore dell'elemento */
	AerodynamicExternalModal* pEl = NULL;

	/* giunto modale collegato */
	const Modal *pModalJoint = pDM->ReadElem<const Modal, const Joint, Elem::JOINT>(HP);
	int comm = 0;
	if (HP.IsKeyWord("comm")) {
		comm = HP.GetInt();
	}

	std::list<MPI::Intercomm>::iterator iComm =  InterfaceComms.begin();
	for (int i = 0; i++ < comm; ++iComm) {}
	MPI::Intercomm* pInterC = &(*(iComm));

	bool VelFlag = false;
	bool MomFlag = false;

	if (HP.IsKeyWord("velocity")) VelFlag = true;
	if (HP.IsKeyWord("moment")) MomFlag = true;

	flag fOut = pDM->fReadOutput(HP, Elem::AERODYNAMIC);
	DriveCaller* pDC = NULL;
	SAFENEWWITHCONSTRUCTOR(pDC, TimeDriveCaller, TimeDriveCaller(pDM->pGetDrvHdl()));

	SAFENEWWITHCONSTRUCTOR(pEl,
		AerodynamicExternalModal,
		AerodynamicExternalModal(uLabel,
			pDO,
			pDC,
			pModalJoint,
			pInterC,
			fOut,
			VelFlag,
			MomFlag));

	return pEl;
}

/* Aerodynamic External Modal - end */

#endif /* USE_EXTERNAL */

