/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2002
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
 
#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include"aeroext.h"
#include <myassert.h>
#include <except.h>

/* le label delle comuncazioni vengono costruite in questo modo:
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
					int NN,
			    		const StructNodes** ppN,
			    		const doublereal* pRefL,
			    		MPI::Intercomm* pIC,
					flag fOut,
					bool VF,
					bool MF)
:Elem(uLabel, Elem::AERODYNAMIC, fOut), 
AerodynamicElem(uLabel, AerodynamicElem::AERODYNAMICEXTERNAL, fOut),
pdBuffer(NULL),
pdBufferVel(NULL),
NodeN(NN),
ppNode(ppN),
pRefLength(pRefL),
OffN(3),
pOffsetVectors(NULL),
pInterfComm(pIC),
pSenReq(NULL),
pRecReq(NULL),
VelFlag(VF),
MomFlag(MF)
{
	/* costruisce la mat 3xN degli offset */
	SAFENEWWITHCONSTRUCTOR(pOffsetVectors, Mat3xN, Mat3xN(OffN, 0.));
	for (int i=0; i<3; i++) pOffsetVectors(i,i,1.);
	ConstructAndInitialize();
	
		
}

AerodynamicExternal::AerodynamicExternal(unsigned int uLabel,
					int NN,
			    		const StructNodes** ppN,
			    		const doublereal* pRefL,
			    		MPI::Intercomm* pIC,
					int	ON,
			    		Mat3xN* OV,
			    		flag fOut,
			    		bool VF,
					bool MF)
:Elem(uLabel, Elem::AERODYNAMIC, fOut), 
AerodynamicElem(uLabel, AerodynamicElem::AERODYNAMICEXTERENAL, fOut),
pdBuffer(NULL),
pdBufferVel(NULL),
NodeN(NN),
ppNode(ppN),
pRefLength(pRefL),
OffN(ON),
pOffsetVectors(OV),
pInterfComm(pIC),
pSenReq(NULL),
pRecReq(NULL),
VelFlag(VF),
MomFlag(MF)
{	
	ConstructAndInitialize();
}	

AerodynamicExternal::~AerodynamicExternal(void)
{
	if (pdBuffer != NULL) SAFEDELETE(pdBuffer);
	if (pdBufferVel != NULL) SAFEDELETE(pdBufferVel);
	if (ppNode != NULL) SAFEDELETEARR(ppNode);
	if (pRefLength != NULL) SAFEDELETEARR(pRefLength);
	if (pSenReq != NULL) SAFEDELETEARR(pSenReq);
	if (pRecReq != NULL) SAFEDELETEARR(pRecReq);
}
		
void AerodynamicExternal::ConstructAndInitialize(void) 
{
	/* dimensiona l'array per il buffer */
	DEBUGCOUTFNAME("AerodynamicExternal::ConstructAndInitialize");
	SAFENEWWITHCONSTRUCTOR(pdBuffer, MyVectorHandler,
			MyVectorHandler(3*NodeN*(OffN+1)));
	
	if (VelFlag || MomFlag) SAFENEWWITHCONSTRUCTOR(pdBufferVel, MyVectorHandler,
			MyVectorHandler(3*NodeN*(OffN+1)));
	
	/* invio delle posizioni iniziali dei nodi per il calcolo della matrice di interpolazione */
	unsigned int*  pLabList = NULL; 
	SAFENEWARR(pLabList, unsigned int, 2*NodeN);
	Mat3xN MatTmp(OffN,0.);
	for (int i=0; i < NodeN; i++) {
		pLabList[i] = ppNode[i]->GetLabel();
		pLabList[i+NodeN] = 1+OffN;
		pdBuffer->Put(i*(OffN+1)*3+1, ppNode[i]->GetXCurr());
		for (int j=0; j < OffN; j++) {
			MatTmp.LeftMult(ppNode[i]->GetRCurr(), *pOffsetVectors);
			MatTmp *= pRefLength[i];
			pdBuffer->Put(i*(OffN+1)*3+(j+1)*3+1, ppNode[i]->GetXCurr() + MatTmp.GetVec(j));
		}
	}
	
	/* il rank del processo con cui si comunica e' zero perche' il comuncatire del 
	codice di interfaccia e' sempre fatto da una macchina sola */
	pInterfComm->Send(pLabList, 2*NodeN, MPI::UNSIGNED, 0,(this->GetLabel())*10+1); 
	pInterfComm->Send(pdBuffer->pdGetVec(), 3*NodeN*(OffN+1), MPI::DOUBLE, 0,(this->GetLabel())*10+2);
	
	/* creo le request per la trasmissione e ricezione dei dati */
	 
	if (VelFlag) {
		SAFENEWARR(pSenReq, MPI::Prequest, 2);
		pSenReq[0] = pInterfComm->Send_init(pdBuffer->pdGetVec(), 3*NodeN*(OffN+1), MPI::DOUBLE, 0,(this->GetLabel())*10+2);
		pSenReq[1] = pInterfComm->Send_init(pdBufferVel->pdGetVec(), 3*NodeN*(OffN+1), MPI::DOUBLE, 0,(this->GetLabel())*10+3);
	} else { 
		SAFENEWARR(pSenReq, MPI::Prequest, 1);
		pSenReq[0] = pInterfComm->Send_init(pdBuffer->pdGetVec(), 3*NodeN*(OffN+1), MPI::DOUBLE, 0,(this->GetLabel())*10+2);
	}			
	if (MomFlag) {
		SAFENEWARR(pRecReq, MPI::Prequest, 2);
		RecReq[0] = pInterfComm->Recv_init(pdBuffer->pdGetVec(), 3*NodeN*(OffN+1), MPI::DOUBLE, 0,(this->GetLabel())*10+4);
		RecReq[1] = pInterfComm->Recv_init(pdBufferVel->pdGetVec(), 3*NodeN*(OffN+1), MPI::DOUBLE, 0,(this->GetLabel())*10+5);
	}else {
		SAFENEWARR(pRecReq, MPI::Prequest, 1);
		RecReq[0] = pInterfComm->Recv_init(pdBuffer->pdGetVec(), 3*NodeN*(OffN+1), MPI::DOUBLE, 0,(this->GetLabel())*10+4);
	}
	SAFEDELETEARR(pLabList);	
}


/* invia posizione e velocita' predetti */ 
void AerodynamicExternal::AfterPredict(VectorHandler& X  , 
					VectorHandler&  XP  )
{

	DEBUGCOUTFNAME("AerodynamicExternal::AfterPredict");
	Mat3xN MatTmp(OffN,0.);
	for (int i=0; i < NodeN; i++) {
		pdBuffer->Put(i*(OffN+1)*3+1, ppNode[i]->GetXCurr());
		if (FlagVel) pdBufferVel->Put(i*(OffN+1)*3+1, ppNode[i]->GetVCurr());
		for (int j=0; j < OffN; j++) {
			MatTmp.LeftMult(ppNode[i]->GetRCurr(), *pOffsetVectors);
			MatTmp *= pRefLength[i];
			pdBuffer->Put(i*(OffN+1)*3+(j+1)*3+1, ppNode[i]->GetXCurr() + MatTmp.GetVec(j));
			if (FlagVel) pdBufferVel->Put(i*(OffN+1)*3+(j+1)*3+1, ppNode[i]->GetVCurr() + ( Mat3x3(ppNode[i]->GetWCurr()) * ppNode[i]->GetRCurr() )* MatTmp.GetVec(j));
		}
	}
		
	MPI::Prequest::Startall(1 + VelFlag, pSenReq);
	/* attiva la ricezione delle forze */
	MPI::Prequest::Startall(1 , *RecReq);
}

/* invia posizione e velocita' predetti */ 
void AerodynamicExternal::Update(VectorHandler& X  , 
					VectorHandler&  XP  )
{
	
	DEBUGCOUTFNAME("AerodynamicExternal::Update");
	Mat3xN MatTmp(OffN,0.);
	for (int i=0; i < NodeN; i++) {
		pdBuffer->Put(i*(OffN+1)*3+1, ppNode[i]->GetXCurr());
		if (FlagVel) pdBufferVel->Put(i*(OffN+1)*3+1, ppNode[i]->GetVCurr());
		for (int j=0; j < OffN; j++) {
			MatTmp.LeftMult(ppNode[i]->GetRCurr(), *pOffsetVectors);
			MatTmp *= pRefLength[i];
			pdBuffer->Put(i*(OffN+1)*3+(j+1)*3+1, ppNode[i]->GetXCurr() + MatTmp.GetVec(j));
			(ppNode[i]->GetXCurr() + MatTmp.GetVec(j)).PutTo(pBuffer[i*(OffN+1)*3+(j+1)*3]);
			if (FlagVel) pdBufferVel->Put(i*(OffN+1)*3+(j+1)*3+1, ppNode[i]->GetVCurr() + ( Mat3x3(ppNode[i]->GetWCurr()) * ppNode[i]->GetRCurr() )* MatTmp.GetVec(j));
		}
	}
	MPI::Prequest::Startall(1 + VelFlag, pSenReq);
	/* attiva la ricezione delle forze */
	MPI::Prequest::Startall(1 , *RecReq);
 
}


SubVectorHandler&
	AerodynamicExternal::AssRes(SubVectorHandler& WorkVec,
	       doublereal dCoef,
	       const VectorHandler& XCurr,
	       const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("AerodynamicExternal::AssRes");
	/* verifica il completamento dei send */
	MPI::Request::Waitall(1 + VelFlag, pSenReq);
   	
	WorkVec.Resize(6*NodeN);
   	WorkVec.Reset(0.);

	
	for (int i=0; i < NodeN; i++) {
		integer iFirstIndex = ppNode[i]->iGetFirstMomentumIndex();
   		for (int iCnt = 1; iCnt <= 6; iCnt++) {
      			WorkVec.fPutRowIndex(iCnt, iFirstIndex+iCnt);
   		}
	}
		
	/* attende il completamento dei receive */
	MPI::Request::Waitall(1 + MomFlag, *RecReq);
	
	/* calcola le forze su ciascun nodo */
	/* le forze e i momeni sono gia' espressi nel sistema di riferimeto global */
	Vec3 F(0.);
	Vec3 M(0.);
	for (int i=0; i < NodeN; i++) {
		F[0] +=  pdBuffer->dGet(i*(OffN+1)*3+1);
		F[1] +=  pdBuffer->dGet(i*(OffN+1)*3+2);
		F[2] +=  pdBuffer->dGet(i*(OffN+1)*3+3);
		if (MomFlag) {
			M[0] +=  pdBufferVel->dGet(i*(OffN+1)*3+1);
			M[1] +=  pdBufferVel->dGet(i*(OffN+1)*3+2);
			M[2] +=  pdBufferVel->dGet(i*(OffN+1)*3+3);
		}
		for (int j=0; j < OffN; j++) {
			F[0] +=  pdBuffer->dGet(i*(OffN+1)*3+(j+1)*3+1);
			F[1] +=  pdBuffer->dGet(i*(OffN+1)*3+(j+1)*3+2);
			F[2] +=  pdBuffer->dGet(i*(OffN+1)*3+(j+1)*3+3);			
			if (MomFlag) {
				M[0] +=  pdBufferVel->dGet(i*(OffN+1)*3+(j+1)*3+1);
				M[1] +=  pdBufferVel->dGet(i*(OffN+1)*3+(j+1)*3+2);
				M[2] +=  pdBufferVel->dGet(i*(OffN+1)*3+(j+1)*3+3);
			}
			Vec3 TmpF(*pdBuffer, i*(OffN+1)*3+(j+1)*3+1);
			M += TmpF.Cross( ppNode[i]->GetRCurr() * MatTmp.GetVec(j) );
		}
		WorkVec.Add(i*6+1, F);
   		WorkVec.Add(i*6+4, M);
	}
	return WorkVec;			
}

Elem *
ReadAerodynamicExternal(DataManager* pDM, MBDynParser& HP, unsigned int uLabel)
{
	
}
/* Aerodynamic External - end */


AerodynamicExternalModal::AerodynamicExternalModal(unsigned int uLabel,
			    				const Modal* pM,
			    				MPI::Intercomm* IC,
			    				flag fOut,
			    				bool VelFlag,
							bool MomFlag)
:Elem(uLabel, Elem::AERODYNAMIC, fOut), 
AerodynamicElem(uLabel, AerodynamicElem::AERODYNAMICEXTERNALMODAL, fOut),
pdBuffer(NULL),
pdBufferVel(NULL),
pModal(pM),
ModalNodes(pM->uGetNFemNodes()),
pModalPos(pM->pGetFemNodesPosition),
pInterfComm(IC),
pSenReq(NULL),
pRecReq(NULL),
VelFlag(VF),
MomFlag(MF)
{
	SAFENEWWITHCONSTRUCTOR(pdBuffer, MyVectorHandler,
			MyVectorHandler(3*ModalNode));
	
	if (VelFlag || MomFlag) SAFENEWWITHCONSTRUCTOR(pdBufferVel, MyVectorHandler,
			MyVectorHandler(3*ModalNode));
	
	unsigned int*  pLabList = NULL; 
	SAFENEWARR(pLabList, unsigned int, 2);
	pLabList[0] = pModal->GetLabel();
	pLabList[1] = ModalNodes;
	/* il rank del processo con cui si comunica e' zero perche' il comuncatore del 
	codice di interfaccia e' sempre fatto da una macchina sola */
	pInterfComm->Send(pLabList, 2, MPI::UNSIGNED, 0,(this->GetLabel())*10+1); 
	
	for (int i=0; i < ModalNodes; i++) {
		pdBuffer->Put(i*3, (pModalPos->GetVec(i)));
		
}

/* Aerodynamic External Modal - begin */

	

Elem *
ReadAerodynamicExternalModal(DataManager* pDM, MBDynParser& HP, unsigned int uLabel)
{
	
}

/* Aerodynamic External Modal - end */
