/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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
//#include <cmath>
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
: Elem(uLabel, Elem::AERODYNAMIC, fout), 
AerodynamicElem(uLabel, AerodynamicElem::AEROMODAL, fout), 
InitialAssemblyElem(uLabel, Elem::AEROMODAL, fout),
ElemWithDofs(uLabel, Elem::AEROMODAL, pDO, fout),
pModalNode(pN), pModalJoint(pMJ),
Ra(RaTmp), 
Chord(Cd),
NStModes(NModal),
NAeroStates(NAero),
NGust(Gust),
pA(pAMat), pB(pBMat), pC(pCMat),
pD0(pD0Mat), pD1(pD1Mat), pD2(pD2Mat),
pq(NULL), pqPrime(NULL), pqSec(NULL), 
pxa(NULL), pxaPrime(NULL), 
pgs(NULL), pgsPrime(NULL),
gustVff(Vff), gustXi(0.707),
RigidF(rgF)
{
   
   R0 = pModalNode->GetRCurr()*Ra;
   P0 = R0.Transpose()*pModalNode->GetXCurr();   
   DEBUGCOUTFNAME("AerodynamicModal::AerodynamicModal");
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
	if (pA != NULL) {
      		SAFEDELETE(pA);
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
	
   	integer iModalIndex = pModalJoint->iGetModalIndex();

   
    	for (unsigned int iCnt = 1; iCnt <= NStModes; iCnt++) {
      		WM.PutRowIndex(iCnt, iModalIndex+NStModes+iCnt);
      		WM.PutColIndex(iCnt, iModalIndex+iCnt);
      		WM.PutColIndex(iCnt+NStModes, iModalIndex+NStModes+iCnt);
    	} 
   
    	integer iAeroIndex = iGetFirstIndex();		   
    	for (unsigned int iCnt = 1; iCnt <= NAeroStates+2*NGust; iCnt++) {
      		WM.PutRowIndex(NStModes+iCnt, iAeroIndex+iCnt);
      		WM.PutColIndex(2*NStModes+iCnt, iAeroIndex+iCnt);
    	}

   	/* Dati del nodo rigido */
   	Vec3 X0(pModalNode->GetXCurr());
   	Vec3 V0(pModalNode->GetVCurr());
   	Mat3x3 Rn(pModalNode->GetRCurr());
   	Mat3x3 RR(Rn*Ra);
   	Mat3x3 RRT(RR.Transpose());
   	Vec3 Vr(V0);

   	doublereal rho, vs, p, T;
	GetAirProps(X0, rho, vs, p, T);		/* p, T are not used yet */

   	Vec3 VTmp(0.);
   	if(fGetAirVelocity(VTmp, X0)) {
   		Vr -= VTmp;
   	}
   	/* velocità nel riferimento nodale aerodinamico */
   	VTmp=RRT*Vr;
   	doublereal nV = fabs(VTmp.dGet(1));
   	doublereal qd =0.5*rho*nV*nV;	
        doublereal CV=Chord/(2*nV);
   
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
    	for (unsigned int i=1; i <= NStModes; i++) {
        	for (unsigned int j=1; j <= NStModes; j++) {
			WM.IncCoef(i, j, -dCoef*qd*pD0->dGetCoef(RigidF+i,
						RigidF+j));
	        	WM.IncCoef(i, j+NStModes, 
					-qd*CV*CV*pD2->dGetCoef(RigidF+i,
						RigidF+j)
					-dCoef*qd*CV*pD1->dGetCoef(RigidF+i,
						RigidF+j));
		}
	}
        for (unsigned int j=1; j <= NAeroStates; j++) {
        	for (unsigned int i=1; i <= NStModes; i++) {
			WM.IncCoef(i, j+2*NStModes,
					-dCoef*qd*pC->dGetCoef(RigidF+i, j));
	        	WM.IncCoef(j+NStModes, i, 
					-dCoef*(1/CV)*pB->dGetCoef(j, 
						RigidF+i));
		}
	}
        for (unsigned int j=1; j <= NAeroStates; j++) {
        	for (unsigned int i=1; i <= NAeroStates; i++) {
		        WM.IncCoef(i+NStModes,j+2*NStModes, 1*(i==j)-dCoef*(1/CV)*pA->dGetCoef(i,j));
		}
	}

	if (NGust) {
        	for (unsigned int i=1; i <= NStModes; i++) {
			WM.IncCoef(i, 2*NStModes+NAeroStates+1,
					-dCoef*qd*pD0->dGetCoef(RigidF+i,
						RigidF+NStModes+1));
			WM.IncCoef(i, 2*NStModes+NAeroStates+3,
					-dCoef*qd*pD0->dGetCoef(RigidF+i,
						RigidF+NStModes+2));
			WM.IncCoef(i, 2*NStModes+NAeroStates+2,
					-CV*CV*qd*pD2->dGetCoef(RigidF+i,
						RigidF+NStModes+1)
					-dCoef*CV*qd*pD1->dGetCoef(RigidF+i,
						RigidF+NStModes+1));
			WM.IncCoef(i, 2*NStModes+NAeroStates+4,
					-CV*CV*qd*pD2->dGetCoef(RigidF+i,
						RigidF+NStModes+2)
					-dCoef*CV*qd*pD1->dGetCoef(RigidF+i,
						RigidF+NStModes+2));
		}
        	for (unsigned int i=1; i <= NAeroStates; i++) {
			WM.IncCoef(i+NStModes, 2*NStModes+NAeroStates+1,
					-dCoef*(1/CV)*pB->dGetCoef(i, RigidF+NStModes+1));
			WM.IncCoef(i+NStModes, 2*NStModes+NAeroStates+3,
					-dCoef*(1/CV)*pB->dGetCoef(i, RigidF+NStModes+2));
		}
        	for (unsigned int i=0; i < NGust; i++) {
			WM.IncCoef(i*2+1+NStModes+NAeroStates,i*2+1+2*NStModes+NAeroStates,1);
			WM.IncCoef(i*2+1+NStModes+NAeroStates,i*2+2+2*NStModes+NAeroStates,-dCoef);
			WM.IncCoef(i*2+2+NStModes+NAeroStates,i*2+1+2*NStModes+NAeroStates,dCoef*gustVff*gustVff);
			WM.IncCoef(i*2+2+NStModes+NAeroStates,i*2+2+2*NStModes+NAeroStates,1+dCoef*2*gustXi*gustVff);
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
   WorkVec.Resize(RigidF+NStModes+NAeroStates+2*NGust);
   WorkVec.Reset(0.);


   if (RigidF) {
   	integer iFirstIndex = pModalNode->iGetFirstIndex();
    	for (int iCnt = 1; iCnt <= RigidF; iCnt++) {
      		WorkVec.PutRowIndex(iCnt, iFirstIndex+RigidF+iCnt);
	}
   	Vec3 X0(pModalNode->GetXCurr());
   	Vec3 V0(pModalNode->GetVCurr());
   	Mat3x3 Rn(pModalNode->GetRCurr());
   	Mat3x3 RR(Rn*Ra);
   	Mat3x3 RRT(RR.Transpose());
      	Mat3x3 R_i = RR*R0.Transpose();
	Vec3 g(gparam(R_i));
	doublereal d(g.Norm());
	Vec3 Phi_i(0.);
	if (d != 0) {
		Phi_i = g*(2./d*atan(d/2.));
	}     
	Vec3 X_i  = RRT*X0 - P0;
	Vec3 V_i  = RRT*V0;
   	Vec3 W_i  = RRT*pModalNode->GetWCurr();
   	Vec3 A_i  = RRT*(Vec3(XPrimeCurr, iFirstIndex+7));
	Vec3 WP_i = RRT*(Vec3(XPrimeCurr, iFirstIndex+10));
	for (unsigned int  iCnt=1; iCnt<=3; iCnt++) {
     		pq->PutCoef(iCnt, X_i.dGet(iCnt));
     		pq->PutCoef(iCnt+3, Phi_i.dGet(iCnt));
     		pqPrime->PutCoef(iCnt, V_i.dGet(iCnt));
     		pqPrime->PutCoef(iCnt+3, W_i.dGet(iCnt));
     		pqSec->PutCoef(iCnt, A_i.dGet(iCnt));
     		pqSec->PutCoef(iCnt+3, WP_i.dGet(iCnt));
	}		
   }
   	
   integer iModalIndex = pModalJoint->iGetModalIndex();

   
    for (unsigned int iCnt = 1; iCnt <= NStModes; iCnt++) {
      WorkVec.PutRowIndex(RigidF+iCnt, iModalIndex+NStModes+iCnt);
    } 
   
    integer iAeroIndex = iGetFirstIndex();		   
    for (unsigned int iCnt = 1; iCnt <= NAeroStates+2*NGust; iCnt++) {
      WorkVec.PutRowIndex(RigidF+NStModes+iCnt, iAeroIndex+iCnt);
    } 
    	


   /* Recupera i vettori {a} e {aP} e {aS}(deformate modali) */  
   for (unsigned int  iCnt=1; iCnt<=NStModes; iCnt++) {
     pq->PutCoef(iCnt+RigidF, XCurr.dGetCoef(iModalIndex+iCnt));
     pqPrime->PutCoef(iCnt+RigidF, XPrimeCurr.dGetCoef(iModalIndex+iCnt));
     pqSec->PutCoef(iCnt+RigidF, XPrimeCurr.dGetCoef(iModalIndex+NStModes+iCnt));
   }
  
   /* Recupera i vettori {xa} e {xaP}  */  
   for (unsigned int  iCnt=1; iCnt<=NAeroStates; iCnt++) {
     pxa->PutCoef(iCnt, XCurr.dGetCoef(iAeroIndex+iCnt));
     pxaPrime->PutCoef(iCnt, XPrimeCurr.dGetCoef(iAeroIndex+iCnt));
   }
   for (unsigned int  iCnt=1; iCnt<=2*NGust; iCnt++) {
     pgs->PutCoef(iCnt, XCurr.dGetCoef(iAeroIndex+NAeroStates+iCnt));
     pgsPrime->PutCoef(iCnt, XPrimeCurr.dGetCoef(iAeroIndex+NAeroStates+iCnt));

   }
   
   AssVec(WorkVec);
   
   return WorkVec;
}


SubVectorHandler& AerodynamicModal::InitialAssRes(SubVectorHandler& WorkVec,
						 const VectorHandler&  XCurr )
{
   DEBUGCOUTFNAME("AerodynamicModal::AssRes");
   Vec3 X0(pModalNode->GetXCurr());
   Mat3x3 Rn(pModalNode->GetRCurr());
   R0 = Rn*Ra;
   Mat3x3 RRT(R0.Transpose());
   P0 = RRT*X0;
    
   WorkVec.Resize(RigidF+NStModes+NAeroStates+2*NGust);
   WorkVec.Reset(0.);

   if (RigidF) {
   	integer iFirstIndex = pModalNode->iGetFirstIndex();
    	for (int iCnt = 1; iCnt <= RigidF; iCnt++) {
      		WorkVec.PutRowIndex(iCnt, iFirstIndex+RigidF+iCnt);
	}
   	Vec3 X0(pModalNode->GetXCurr());
   	Vec3 V0(pModalNode->GetVCurr());
   	Mat3x3 Rn(pModalNode->GetRCurr());
   	Mat3x3 RR(Rn*Ra);
   	Mat3x3 RRT(RR.Transpose());
      	Mat3x3 R_i = RR*R0.Transpose();
	Vec3 g(gparam(R_i));
	doublereal d(g.Norm());
	Vec3 Phi_i(0.);
	if (d != 0) {
		Phi_i = g*(2./d*atan(d/2.));
	}     
	Vec3 X_i  = RRT*X0 - P0;
	Vec3 V_i  = RRT*V0;
   	Vec3 W_i  = RRT*pModalNode->GetWCurr();
	for (unsigned int  iCnt=1; iCnt<=3; iCnt++) {
     		pq->PutCoef(iCnt, X_i.dGet(iCnt));
     		pq->PutCoef(iCnt+3, Phi_i.dGet(iCnt));
     		pqPrime->PutCoef(iCnt, V_i.dGet(iCnt));
     		pqPrime->PutCoef(iCnt+3, W_i.dGet(iCnt));
     		pqSec->PutCoef(iCnt, 0.);
     		pqSec->PutCoef(iCnt+3, 0.);
	}		
   }
   
   integer iModalIndex = pModalJoint->iGetModalIndex();

   
    for (unsigned int iCnt = 1; iCnt <= NStModes; iCnt++) {
      WorkVec.PutRowIndex(RigidF+iCnt, iModalIndex+NStModes+iCnt);
    } 
   
    integer iAeroIndex = iGetFirstIndex();		   
    for (unsigned int iCnt = 1; iCnt <= NAeroStates+2*NGust; iCnt++) {
      WorkVec.PutRowIndex(RigidF+NStModes+iCnt, iAeroIndex+iCnt);
    } 


   /* Recupera i vettori {a} e {aP} e {aS}(deformate modali) */  
   for (unsigned int  iCnt=1; iCnt<=NStModes; iCnt++) {
     pq->PutCoef(iCnt+RigidF, XCurr.dGetCoef(iModalIndex+iCnt));
     pqPrime->PutCoef(iCnt+RigidF, 0);
     pqSec->PutCoef(iCnt+RigidF, 0);
   }
  
   /* Recupera i vettori {xa} e {xaP}  */  
   for (unsigned int  iCnt=1; iCnt<=NAeroStates; iCnt++) {
     pxa->PutCoef(iCnt, XCurr.dGetCoef(iAeroIndex+iCnt));
     pxaPrime->PutCoef(iCnt, 0);
   }
   for (unsigned int  iCnt=1; iCnt<=2*NGust; iCnt++) {
     pgs->PutCoef(iCnt, XCurr.dGetCoef(iAeroIndex+NAeroStates+iCnt));
     pgsPrime->PutCoef(iCnt, 0.);
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

   doublereal rho, vs, p, T;
   GetAirProps(X0, rho, vs, p, T);	/* p, T are not used yet */

   Vec3 Vair(0.);
   if(fGetAirVelocity(Vair, X0)) {
   	Vr -= Vair;
   }
    /* velocita' nel riferimento nodale aerodinamico */
   V0=RRT*Vr;
   doublereal nV = fabs(V0.dGet(1));
   doublereal qd =0.5*rho*nV*nV;	
   MyVectorHandler TmpA(NAeroStates);
   TmpA.Reset();
   MyVectorHandler Tmp(RigidF+NStModes);
   Tmp.Reset();
   MyVectorHandler TmpP(RigidF+NStModes);
   TmpP.Reset();
   MyVectorHandler TmpS(RigidF+NStModes);
   TmpS.Reset();
   for (unsigned int i=1; i <= NAeroStates; i++) {
   	for (unsigned int j=1; j <= NStModes+RigidF; j++) {
		TmpA.IncCoef(i, pB->dGetCoef(i,j)*(pq->dGetCoef(j)));
	}
   }   
   pA->MatVecIncMul(TmpA,*pxa);
   doublereal CV=Chord/(2*nV);
   for (unsigned int i=1; i <= NAeroStates; i++) {
 	WorkVec.Add(RigidF+NStModes+i, -pxaPrime->dGetCoef(i)+(1/CV)*(TmpA.dGetCoef(i)));
   }
   for (unsigned int i=1; i <= NStModes+RigidF; i++) {
   	for (unsigned int j=1; j <= NStModes+RigidF; j++) {
		Tmp.IncCoef(i, pD0->dGetCoef(i,j)*(pq->dGetCoef(j)));
		TmpP.IncCoef(i, pD1->dGetCoef(i,j)*(pqPrime->dGetCoef(j)));
		TmpS.IncCoef(i, pD2->dGetCoef(i,j)*(pqSec->dGetCoef(j)));
	}
   }   
   for (unsigned int i=1; i <= NStModes+RigidF; i++) {
   	for (unsigned int j=1; j <= NAeroStates; j++) {
		Tmp.IncCoef(i, pC->dGetCoef(i,j)*(pxa->dGetCoef(j)));
        }
   }
   if (RigidF) {
   	Vec3 F(0.);
   	Vec3 M(0.);	
   	for (unsigned int i=1; i <= 3; i++) {
		F.Put(i,qd*Tmp.dGetCoef(i)+qd*CV*TmpP.dGetCoef(i)+qd*CV*CV*TmpS.dGetCoef(i));
		M.Put(i,qd*Tmp.dGetCoef(i+3)+qd*CV*TmpP.dGetCoef(i+3)+qd*CV*CV*TmpS.dGetCoef(i+3));
   	}
//	std::cout << F << std::endl;
	F = RR*(-F);
	M = RR*M;
	WorkVec.Add(1,F);
	WorkVec.Add(4,M);
   }	
   for (unsigned int i=1+RigidF; i <= NStModes+RigidF; i++) {
	WorkVec.Add(i,qd*Tmp.dGetCoef(i)+qd*CV*TmpP.dGetCoef(i)+qd*CV*CV*TmpS.dGetCoef(i));
   }
   
   if (NGust) {
#if 0 /* unused */
   	doublereal Vyg = Vair.dGet(2)/nV; 
   	doublereal Vzg = Vair.dGet(3)/nV; 
#endif
   	for (unsigned int i=1; i <= NAeroStates; i++) {
 		WorkVec.Add(RigidF+NStModes+i, (1/CV)*pB->dGetCoef(i,RigidF+NStModes+1)*pgs->dGetCoef(1)
		+(1/CV)*pB->dGetCoef(i,RigidF+NStModes+2)*pgs->dGetCoef(3));
   	}
   	if (RigidF) {
   		Vec3 F(0.);
   		Vec3 M(0.);	
   		for (unsigned int i=1; i <= 3; i++) {
			F.Put(i,qd*pD0->dGetCoef(i,RigidF+NStModes+1)*pgs->dGetCoef(1)+qd*pD0->dGetCoef(i,RigidF+NStModes+2)*pgs->dGetCoef(3)
		              +qd*CV*pD1->dGetCoef(i,RigidF+NStModes+1)*pgs->dGetCoef(2)
			      +qd*CV*pD1->dGetCoef(i,RigidF+NStModes+2)*pgs->dGetCoef(4)
			      +qd*CV*CV*pD2->dGetCoef(i,RigidF+NStModes+1)*pgsPrime->dGetCoef(2)
			      +qd*CV*CV*pD2->dGetCoef(i,RigidF+NStModes+2)*pgsPrime->dGetCoef(4));
			M.Put(i,qd*pD0->dGetCoef(i+3,RigidF+NStModes+1)*pgs->dGetCoef(1)+qd*pD0->dGetCoef(i+3,RigidF+NStModes+2)*pgs->dGetCoef(3)
		              +qd*CV*pD1->dGetCoef(i+3,RigidF+NStModes+1)*pgs->dGetCoef(2)
			      +qd*CV*pD1->dGetCoef(i+3,RigidF+NStModes+2)*pgs->dGetCoef(4)
			      +qd*CV*CV*pD2->dGetCoef(i+3,RigidF+NStModes+1)*pgsPrime->dGetCoef(2)
			      +qd*CV*CV*pD2->dGetCoef(i+3,RigidF+NStModes+2)*pgsPrime->dGetCoef(4));
   		}
	//	std::cout << F << std::endl;
		F = RR*(-F);
		M = RR*M;
		WorkVec.Add(1,F);
		WorkVec.Add(4,M);
   	}	
   	for (unsigned int i=1+RigidF; i <= NStModes+RigidF; i++) {
		WorkVec.Add(i,qd*pD0->dGetCoef(i,RigidF+NStModes+1)*pgs->dGetCoef(1)+qd*pD0->dGetCoef(i,RigidF+NStModes+2)*pgs->dGetCoef(3)
		              +qd*CV*pD1->dGetCoef(i,RigidF+NStModes+1)*pgs->dGetCoef(2)
			      +qd*CV*pD1->dGetCoef(i,RigidF+NStModes+2)*pgs->dGetCoef(4)
			      +qd*CV*CV*pD2->dGetCoef(i,RigidF+NStModes+1)*pgsPrime->dGetCoef(2)
			      +qd*CV*CV*pD2->dGetCoef(i,RigidF+NStModes+2)*pgsPrime->dGetCoef(4));
	}
	for (unsigned int i=0; i < NGust; i++) {
		WorkVec.Add(RigidF+i*2+1+NStModes+NAeroStates,-pgsPrime->dGetCoef(1+i*2)+pgs->dGetCoef(2+i*2));
		WorkVec.Add(RigidF+i*2+2+NStModes+NAeroStates,-pgsPrime->dGetCoef(2+i*2)
			-gustVff*gustVff*pgs->dGetCoef(1+i*2)-2*gustXi*gustVff*pgs->dGetCoef(2+i*2)
			+gustVff*gustVff*Vair.dGet(2+i)/nV);
	}
 //  std::cout << X0 << std::endl;
 //  std::cout << WorkVec.dGetCoef(3) << std::endl;
   }	
}


/* output; si assume che ogni tipo di elemento sappia, attraverso
 * l'OutputHandler, dove scrivere il proprio output */
void AerodynamicModal::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		OH.AeroModals() << std::setw(8) << GetLabel() << " ";
		for (unsigned int iCnt=1; iCnt <= NAeroStates; iCnt++){
			OH.AeroModals() << " " << pxa->dGetCoef(iCnt);
		}
		for (unsigned int iCnt=1; iCnt <= NAeroStates+2*NGust; iCnt++){
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
    *  {rigid, gust, gust filter cut-off frequency}
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

   AerodynamicModal::RigidF_t rigidF = AerodynamicModal::NO_RIGID;
   if (HP.IsKeyWord("rigid")) {
	rigidF = AerodynamicModal::RIGID;
	NModes += rigidF;
   }   
   unsigned int GustN = 0;
   doublereal Vff = 0.;
   /* Eventuale raffica */
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
     std::cerr << std::endl << "Unable to open file '" << sFileData << '\''
	     << std::endl;
     THROW(DataManager::ErrGeneric());
   }	
   SpMapMatrixHandler* pAMat = NULL;
   FullMatrixHandler* pBMat = NULL;	   
   FullMatrixHandler* pCMat = NULL;		   
   FullMatrixHandler* pD0Mat = NULL;
   FullMatrixHandler* pD1Mat = NULL;
   FullMatrixHandler* pD2Mat = NULL;
   SAFENEWWITHCONSTRUCTOR(pAMat, SpMapMatrixHandler, SpMapMatrixHandler(AeroN, AeroN));
   pAMat->Reset();
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
		for (unsigned int jCnt = 0; jCnt < AeroN; jCnt++)  { 
			fdat >> d;
			if (d != 0) { pAMat->PutCoef(iCnt+1,jCnt+1,d);} 
		}
	}
	/* legge il primo blocco (HEADER) */
      } else if (!strncmp("*** MATRIX B", str, 12)) {  
	for (unsigned int iCnt = 0; iCnt < AeroN; iCnt++)  { 
		for (unsigned int jCnt = 0; jCnt < NModes+GustN; jCnt++)  { 
			fdat >> d;
			pBMat->PutCoef(iCnt+1,jCnt+1,d);
		} 
	}
	
	/* legge il primo blocco (HEADER) */
      } else if (!strncmp("*** MATRIX C", str, 12)) {  
	for (unsigned int iCnt = 0; iCnt < NModes; iCnt++)  { 
		for (unsigned int jCnt = 0; jCnt < AeroN; jCnt++)  { 
			fdat >> d;
			pCMat->PutCoef(iCnt+1,jCnt+1,d);
		} 
	}

	/* legge il primo blocco (HEADER) */
      } else if (!strncmp("*** MATRIX D0", str, 13)) {  
	for (unsigned int iCnt = 0; iCnt < NModes; iCnt++)  { 
		for (unsigned int jCnt = 0; jCnt < NModes+GustN; jCnt++)  { 
			fdat >> d;
			pD0Mat->PutCoef(iCnt+1,jCnt+1,d);
		} 
	}

	/* legge il primo blocco (HEADER) */
      } else if (!strncmp("*** MATRIX D1", str, 13)) {
	for (unsigned int iCnt = 0; iCnt < NModes; iCnt++)  { 
		for (unsigned int jCnt = 0; jCnt < NModes+GustN; jCnt++)  { 
			fdat >> d;
			pD1Mat->PutCoef(iCnt+1,jCnt+1,d);
		} 
	}
	
	/* legge il primo blocco (HEADER) */
      } else if (!strncmp("*** MATRIX D2", str, 13)) {  
	for (unsigned int iCnt = 0; iCnt < NModes; iCnt++)  { 
		for (unsigned int jCnt = 0; jCnt < NModes+GustN; jCnt++)  { 
			fdat >> d;
			pD2Mat->PutCoef(iCnt+1,jCnt+1,d);
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
			                   Ra, pDO, Chord, NModes-rigidF, 
					   AeroN, rigidF, GustN, Vff,
					   pAMat, pBMat, pCMat, pD0Mat, pD1Mat,
					   pD2Mat, fOut));
   
   /* Se non c'e' il punto e virgola finale */
   if (HP.IsArg()) {
      std::cerr << "semicolon expected at line " 
	      << HP.GetLineData() << std::endl;      
      THROW(DataManager::ErrGeneric());
   }   
   
   return pEl;
} /* End of DataManager::ReadAerodynamicModal() */

