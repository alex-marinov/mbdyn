/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <aeromodal.h>
#include <dataman.h>

/* AerodynamicModal - begin */

AerodynamicModal::AerodynamicModal(unsigned int uLabel, 
		   		const StructNode* pN, 
		   		const Modal* pMJ,
				const Mat3x3& RaTmp,
		   		const DofOwner* pDO,
				doublereal Cd,
				const int NModal,
		   		const int NAero,
		   		const int Gust,
		   		std::vector<double>* pAMat,
		   		FullMatrixHandler* pBMat,
		   		FullMatrixHandler* pCMat,
		   		FullMatrixHandler* pD0Mat,
		   		FullMatrixHandler* pD1Mat,
		   		FullMatrixHandler* pD2Mat,
		   		const DriveCaller* pDC, flag fout)
: Elem(uLabel, Elem::AERODYNAMIC, fout), 
AerodynamicElem(uLabel, AerodynamicElem::AEROMODAL, fout), 
InitialAssemblyElem(uLabel, Elem::AEROMODAL, fout),
ElemWithDofs(uLabel, Elem::AEROMODAL, pDO, fout),
DriveOwner(pDC),
pModalNode(pN), pModalJoint(pMJ),
Ra(RaTmp),
Chord(Cd),
NStModes(NModal),
NAeroStates(NAero), NGust(Gust),
pA(pAMat), pB(pBMat), pC(pCMat),
pD0(pD0Mat), pD1(pD1Mat), pD2(pD2Mat)
{
   DEBUGCOUTFNAME("AerodynamicModal::AerodynamicModal");
   SAFENEWWITHCONSTRUCTOR(pq, MyVectorHandler,
   				 MyVectorHandler(NStModes));
   SAFENEWWITHCONSTRUCTOR(pqPrime, MyVectorHandler,
   				 MyVectorHandler(NStModes));
   SAFENEWWITHCONSTRUCTOR(pqSec, MyVectorHandler,
   				 MyVectorHandler(NStModes));
   SAFENEWWITHCONSTRUCTOR(pxa, MyVectorHandler,
   				 MyVectorHandler(NAeroStates));
   SAFENEWWITHCONSTRUCTOR(pxaPrime, MyVectorHandler,
   				 MyVectorHandler(NAeroStates));
				 
   ASSERT(pModalNode != NULL);
   ASSERT(pModalNode->GetNodeType() == Node::STRUCTURAL);
   
}


AerodynamicModal::~AerodynamicModal(void)
{
   	DEBUGCOUTFNAME("AerodynamicModal::~AerodynamicModal");
  	if (pq != NULL) {
      		SAFEDELETE(pq);
   	}
  	if (pqPrime != NULL) {
      		SAFEDELETE(pqPrime);
   	}
	if (pqSec != NULL) {
      		SAFEDELETE(pqSec);
   	}
	if (pxa != NULL) {
      		SAFEDELETE(pxa);
   	}
	if (pxaPrime != NULL) {
      		SAFEDELETE(pxaPrime);
   	}
	if (pB != NULL) {
      		SAFEDELETE(pB);
   	}
	if (pC != NULL) {
      		SAFEDELETE(pC);
   	}
	if (pD0 != NULL) {
      		SAFEDELETE(pD0);
   	}
	if (pD1 != NULL) {
      		SAFEDELETE(pD1);
   	}
	if (pD2 != NULL) {
      		SAFEDELETE(pD2);
   	}

}


/* Scrive il contributo dell'elemento al file di restart */
std::ostream& AerodynamicModal::Restart(std::ostream& out) const
{

   return out << "  /* aerodynamic modal: not implemented yet */" << std::endl;
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler& 
AerodynamicModal::AssJac(VariableSubMatrixHandler& WorkMat,
	    		doublereal dCoef ,
	    		const VectorHandler& XCurr ,
	    		const VectorHandler&  XPrimeCurr )
{
	DEBUGCOUT("Entering AerodynamicModal::AssJac()" << std::endl);  
   	FullSubMatrixHandler& WM = WorkMat.SetFull();
   	integer iNumRows = 0;
   	integer iNumCols = 0;
   	WorkSpaceDim(&iNumRows, &iNumCols);
   	WM.ResizeInit(iNumRows, iNumCols, 0.);

#if 0
   	integer iFirstIndex = pModalNode->iGetFirstIndex();
#endif
   	integer iModalIndex = pModalJoint->iGetModalIndex();

   
    	for (unsigned int iCnt = 1; iCnt <= 2*NStModes; iCnt++) {
      		WM.fPutRowIndex(iCnt, iModalIndex+iCnt);
      		WM.fPutColIndex(iCnt, iModalIndex+iCnt);
    	} 
   
    	integer iAeroIndex = iGetFirstIndex();		   
    	for (unsigned int iCnt = 1; iCnt <= NAeroStates; iCnt++) {
      		WM.fPutRowIndex(2*NStModes+iCnt, iAeroIndex+iCnt);
      		WM.fPutColIndex(2*NStModes+iCnt, iAeroIndex+iCnt);
    	}

   	/* Dati del nodo rigido */
   	Vec3 X0(pModalNode->GetXCurr());
   	Vec3 V0(pModalNode->GetVCurr());
   	Mat3x3 Rn(pModalNode->GetRCurr());
   	Mat3x3 RR(Rn*Ra);
   	Mat3x3 RRT(RR.Transpose());
   	Vec3 Vr(V0);
   	doublereal rho = dGetAirDensity(X0);
#if 0
   	doublereal vs  = dGetSoundSpeed(X0);
#endif
   	Vec3 VTmp(0.);
   	if(fGetAirVelocity(VTmp, X0)) {
   		Vr -= VTmp;
   	}
   	/* velocità nel riferimento nodale aerodinamico */
   	VTmp=RRT*Vr;
   	doublereal nV = VTmp.Norm();
   	doublereal qd =0.5*rho*nV*nV;	
        doublereal CV=Chord/(2*nV);
   
	/* parte deformabile :
    	 * 
         * |     0          0          0    ||aP|
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
        for (unsigned int i=1; i <= NStModes; i++) {
        	for (unsigned int j=1; j <= NStModes; j++) {
			WM.fIncCoef(i+NStModes,j,-dCoef*qd*pD0->dGetCoef(i,j));
	        	WM.fIncCoef(i+NStModes,j+NStModes,-qd*CV*CV*pD2->dGetCoef(i,j)-dCoef*qd*CV*pD1->dGetCoef(i,j));
		}
	}
        for (unsigned int j=1; j <= NAeroStates; j++) {
	        WM.fIncCoef(j+2*NStModes,j+2*NStModes, 1-dCoef*(1/CV)*(*pA)[j-1]);
        	for (unsigned int i=1; i <= NStModes; i++) {
			WM.fIncCoef(i+NStModes,j+2*NStModes,-dCoef*qd*pC->dGetCoef(i,j));
	        	WM.fIncCoef(j+2*NStModes,i,-dCoef*(1/CV)*pB->dGetCoef(j,i));
		}
	}
	
    	return WorkMat;   
}


SubVectorHandler& AerodynamicModal::AssRes(SubVectorHandler& WorkVec,
					  doublereal /* dCoef */ ,
					  const VectorHandler& XCurr,
					  const VectorHandler&  XPrimeCurr )
{
   DEBUGCOUTFNAME("AerodynamicModal::AssRes");
   WorkVec.Resize(NStModes+NAeroStates);
   WorkVec.Reset(0.);

#if 0 
   integer iFirstIndex = pModalNode->iGetFirstIndex();
#endif
   integer iModalIndex = pModalJoint->iGetModalIndex();

   
    for (unsigned int iCnt = 1; iCnt <= NStModes; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iModalIndex+NStModes+iCnt);
    } 
   
    integer iAeroIndex = iGetFirstIndex();		   
    for (unsigned int iCnt = 1; iCnt <= NAeroStates; iCnt++) {
      WorkVec.fPutRowIndex(NStModes+iCnt, iAeroIndex+iCnt);
    } 


   /* Recupera i vettori {a} e {aP} e {aS}(deformate modali) */  
   for (unsigned int  iCnt=1; iCnt<=NStModes; iCnt++) {
     pq->fPutCoef(iCnt, XCurr.dGetCoef(iModalIndex+iCnt));
     pqPrime->fPutCoef(iCnt, XPrimeCurr.dGetCoef(iModalIndex+iCnt));
     pqSec->fPutCoef(iCnt, XPrimeCurr.dGetCoef(iModalIndex+NStModes+iCnt));
   }
  
   /* Recupera i vettori {xa} e {xaP}  */  
   for (unsigned int  iCnt=1; iCnt<=NAeroStates; iCnt++) {
     pxa->fPutCoef(iCnt, XCurr.dGetCoef(iAeroIndex+iCnt));
     pxaPrime->fPutCoef(iCnt, XPrimeCurr.dGetCoef(iAeroIndex+iCnt));
   }
   AssVec(WorkVec);
   
   return WorkVec;
}


SubVectorHandler& AerodynamicModal::InitialAssRes(SubVectorHandler& WorkVec,
						 const VectorHandler&  XCurr )
{
   DEBUGCOUTFNAME("AerodynamicModal::AssRes");
   WorkVec.Resize(NStModes+NAeroStates);
   WorkVec.Reset(0.);

#if 0 
   integer iFirstIndex = pModalNode->iGetFirstIndex();
#endif
   integer iModalIndex = pModalJoint->iGetModalIndex();

   
    for (unsigned int iCnt = 1; iCnt <= NStModes; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iModalIndex+NStModes+iCnt);
    } 
   
    integer iAeroIndex = iGetFirstIndex();		   
    for (unsigned int iCnt = 1; iCnt <= NAeroStates; iCnt++) {
      WorkVec.fPutRowIndex(NStModes+iCnt, iAeroIndex+iCnt);
    } 


   /* Recupera i vettori {a} e {aP} e {aS}(deformate modali) */  
   for (unsigned int  iCnt=1; iCnt<=NStModes; iCnt++) {
     pq->fPutCoef(iCnt, XCurr.dGetCoef(iModalIndex+iCnt));
     pqPrime->fPutCoef(iCnt, 0);
     pqSec->fPutCoef(iCnt, 0);
   }
  
   /* Recupera i vettori {xa} e {xaP}  */  
   for (unsigned int  iCnt=1; iCnt<=NAeroStates; iCnt++) {
     pxa->fPutCoef(iCnt, XCurr.dGetCoef(iAeroIndex+iCnt));
     pxaPrime->fPutCoef(iCnt, 0);
   }
   AssVec(WorkVec);
   
   return WorkVec;
}


/* assemblaggio residuo */
void AerodynamicModal::AssVec(SubVectorHandler& WorkVec)
{
   DEBUGCOUTFNAME("AerodynamicModal::AssVec");
      
   /* Dati del nodo rigido */
   Vec3 X0(pModalNode->GetXCurr());
   Vec3 V0(pModalNode->GetVCurr());
   Mat3x3 Rn(pModalNode->GetRCurr());
   Mat3x3 RR(Rn*Ra);
   Mat3x3 RRT(RR.Transpose());
   Vec3 Vr(V0);
   doublereal rho = dGetAirDensity(X0);
#if 0
   doublereal vs  = dGetSoundSpeed(X0);
#endif
   Vec3 VTmp(0.);
   if(fGetAirVelocity(VTmp, X0)) {
   	Vr -= VTmp;
   }
   /* velocita' nel riferimento nodale aerodinamico */
   VTmp=RRT*Vr;
   doublereal nV = VTmp.Norm();
   doublereal qd =0.5*rho*nV*nV;	
   MyVectorHandler Tmp(NAeroStates);
   MyVectorHandler TmpP(NAeroStates);
   MyVectorHandler TmpS(NAeroStates);
   pB->MatVecMul(Tmp,*pq);
   doublereal CV=Chord/(2*nV);
   for (unsigned int i=1; i <= NAeroStates; i++) {
	WorkVec.Add(NStModes+i, -pxaPrime->dGetCoef(i)+(1/CV)*(*pA)[i-1]*(pxa->dGetCoef(i))+(1/CV)*Tmp.dGetCoef(i));
    }
   pD0->MatVecMul(Tmp,*pq);
   pC->MatVecIncMul(Tmp,*pxa);   
   pD1->MatVecMul(TmpP,*pqPrime);
   pD2->MatVecMul(TmpS,*pqSec);
   for (unsigned int i=1; i <= NStModes; i++) {
	WorkVec.Add(i,+qd*Tmp.dGetCoef(i)+qd*CV*TmpP.dGetCoef(i)+qd*CV*CV*TmpS.dGetCoef(i));
    }


}


/* output; si assume che ogni tipo di elemento sappia, attraverso
 * l'OutputHandler, dove scrivere il proprio output */
void AerodynamicModal::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		OH.AeroModals() << std::setw(8) << GetLabel() << " ";
		for (unsigned int iCnt=0; iCnt < NAeroStates; iCnt++){
			OH.AeroModals() << " " << pxa->dGetCoef(iCnt);
		}
		for (unsigned int iCnt=0; iCnt < NAeroStates; iCnt++){
			OH.AeroModals() << " " << pxaPrime->dGetCoef(iCnt);
		}
		OH.AeroModals() << std::endl;
	}
}

/* AerodynamicModal - end */




Elem* ReadAerodynamicModal(DataManager* pDM,
			  MBDynParser& HP,
			  const DofOwner* pDO,
			  unsigned int uLabel)
{
    DEBUGCOUTFNAME("ReadAerodynamicModal");

   /* formato dell'input: 
    *
    *  label, 
    *  modal node,
    *  reference modal joint,   
    *  rotation of the aerodynamic reference in respect of nodal reference,
    *  reference cord,
    *  number of aerodynamic states,
    *  {number of gusts,
    *  gust driver,}
    *  file name containing state space model matrices;
    */             
 
    	  
   /* Nodo */
   StructNode* pModalNode = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);

   /* giunto modale collegato */		     
   Elem* pM = pDM->ReadElem(HP, Elem::JOINT);
   Modal* pModalJoint = (Modal*)pM->pGet();
   if (pModalJoint->GetJointType() != Joint::MODAL) {
      std::cerr << "Joint " << pModalJoint->GetLabel()
	      << " is required to be a modal joint" << std::endl;
      THROW(DataManager::ErrGeneric());
   }
   
   ReferenceFrame RF(pModalNode);

   Mat3x3 Ra(HP.GetRotRel(RF));

   doublereal Chord = HP.GetReal();

   unsigned int AeroN = HP.GetInt();
   /* numero modi e FEM */
   unsigned int NModes = pModalJoint->uGetNModes();


   unsigned int GustN = 0;
   DriveCaller* pDC = NULL;
   /* Eventuale raffica */
   if (HP.IsKeyWord("gust")) {
	GustN = HP.GetInt();
      	pDC = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());   
   }   

    /* apre il file contenente le matrici A B C D0 D1 D2 */
   const char *sFileData = HP.GetFileName();
   std::ifstream fdat(sFileData);       
   DEBUGCOUT("Reading Aerodynamic State Space Matrices from file " << sFileData << std::endl);
   if (!fdat) {
     std::cerr << std::endl << "Unable to open file " << sFileData << std::endl;
     THROW(DataManager::ErrGeneric());
   }	
   std::vector<double>* pAMat = new  std::vector<double>(AeroN);

   FullMatrixHandler* pBMat = NULL;	   
   FullMatrixHandler* pCMat = NULL;		   
   FullMatrixHandler* pD0Mat = NULL;
   FullMatrixHandler* pD1Mat = NULL;
   FullMatrixHandler* pD2Mat = NULL;
   SAFENEWWITHCONSTRUCTOR(pBMat, FullMatrixHandler, FullMatrixHandler(AeroN, NModes+GustN, 0.));
   SAFENEWWITHCONSTRUCTOR(pCMat, FullMatrixHandler, FullMatrixHandler(NModes, AeroN, 0.));
   SAFENEWWITHCONSTRUCTOR(pD0Mat, FullMatrixHandler, FullMatrixHandler(NModes, NModes+GustN, 0.));
   SAFENEWWITHCONSTRUCTOR(pD1Mat, FullMatrixHandler, FullMatrixHandler(NModes, NModes+GustN, 0.));
   SAFENEWWITHCONSTRUCTOR(pD2Mat, FullMatrixHandler, FullMatrixHandler(NModes, NModes+GustN, 0.));
   
   doublereal d;
   char str[256];

   while (!fdat.eof()) {        /* parsing del file */ 
      fdat.getline(str, sizeof(str));

	/* legge il primo blocco (HEADER) */
      if (!strncmp("*** MATRIX A", str, 12)) {  
	for (unsigned int iCnt = 0; iCnt < AeroN; iCnt++)  { 
		fdat >> d;
		(*pAMat)[iCnt] = d; 
	}
	
	/* legge il primo blocco (HEADER) */
      } else if (!strncmp("*** MATRIX B", str, 12)) {  
	for (unsigned int iCnt = 0; iCnt < AeroN; iCnt++)  { 
		for (unsigned int jCnt = 0; jCnt < NModes+GustN; jCnt++)  { 
			fdat >> d;
			pBMat->fPutCoef(iCnt+1,jCnt+1,d);
		} 
	}
	
	/* legge il primo blocco (HEADER) */
      } else if (!strncmp("*** MATRIX C", str, 12)) {  
	for (unsigned int iCnt = 0; iCnt < NModes; iCnt++)  { 
		for (unsigned int jCnt = 0; jCnt < AeroN; jCnt++)  { 
			fdat >> d;
			pCMat->fPutCoef(iCnt+1,jCnt+1,d);
		} 
	}

	/* legge il primo blocco (HEADER) */
      } else if (!strncmp("*** MATRIX D0", str, 13)) {  
	for (unsigned int iCnt = 0; iCnt < NModes; iCnt++)  { 
		for (unsigned int jCnt = 0; jCnt < NModes+GustN; jCnt++)  { 
			fdat >> d;
			pD0Mat->fPutCoef(iCnt+1,jCnt+1,d);
		} 
	}

	/* legge il primo blocco (HEADER) */
      } else if (!strncmp("*** MATRIX D1", str, 13)) {
	for (unsigned int iCnt = 0; iCnt < NModes; iCnt++)  { 
		for (unsigned int jCnt = 0; jCnt < NModes+GustN; jCnt++)  { 
			fdat >> d;
			pD1Mat->fPutCoef(iCnt+1,jCnt+1,d);
		} 
	}
	
	/* legge il primo blocco (HEADER) */
      } else if (!strncmp("*** MATRIX D2", str, 13)) {  
	for (unsigned int iCnt = 0; iCnt < NModes; iCnt++)  { 
		for (unsigned int jCnt = 0; jCnt < NModes+GustN; jCnt++)  { 
			fdat >> d;
			pD2Mat->fPutCoef(iCnt+1,jCnt+1,d);
		} 
	}   
      }
   }
   fdat.close();  

   flag fOut = pDM->fReadOutput(HP, Elem::AEROMODAL);
   
   Elem* pEl = NULL;
   SAFENEWWITHCONSTRUCTOR(pEl, 
			  AerodynamicModal,
			  AerodynamicModal(uLabel, pModalNode, pModalJoint, 
			                   Ra, pDO, Chord, NModes, 
					   AeroN, GustN, 
					   pAMat, pBMat, pCMat, pD0Mat, pD1Mat,
					   pD2Mat, pDC, fOut));
   
   /* Se non c'e' il punto e virgola finale */
   if (HP.fIsArg()) {
      std::cerr << "semicolon expected at line " 
	      << HP.GetLineData() << std::endl;      
      THROW(DataManager::ErrGeneric());
   }   
   
   return pEl;
} /* End of DataManager::ReadAerodynamicModal() */

