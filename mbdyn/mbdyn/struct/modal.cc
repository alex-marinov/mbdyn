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

/* 
   classe per l'introduzione della flessibilita' modale nel codice multi-corpo

   29/7/99:  implementata la parte di corpo rigido (solo massa e momenti d'inerzia)
             verificata con massa isolata soggetta a forze e momenti esterni
           
   13/9/99:  corpo rigido OK (anche momenti statici ed eq. di vincolo)
             corpo deformabile isolato OK (verificato con osc. armonico)

   23/9/99:  corpo deformabile + eq. di vincolo (no moto rigido)  OK (verificato con trave incastrata)   
             corpo rigido + deformabile OK (verificato con trave libera)   
 
     10/99:  validazioni con NASTRAN: trave incastrata, trave libera con forza in mezzeria,
                                      trave libera con forza all'estremita'
  
   4/11/99:  completate le funzioni InitialAssJac e InitialAssRes per l'assemblaggio iniziale del 'vincolo' 

  22/11/99:  modificata la funzione ReadModal per leggere il file generato da DADS

  26/11/99:  validazione con piastra vincolata con elementi elastici 

     12/99:  validazione con piastra e bauchau

   01/2000:  modi rotanti

  30/11/99:  aggiunto lo smorzamento strutturale tra i dati d'ingresso
             aggiunte le inerzie concentrate 

   1/12/99:  modifiche alla lettura dati d'ingresso (l'estensione *.fem al file con i dati del
             modello ad elementi viene aggiunta automaticamente dal programma, che crea un file di
             output con l'estensione *.mod)
             corretto un bug nella scrittura delle equazioni di vincolo 

  17/12/99:  aggiunta la possibilita' di definire uno smorzamento strutturale variabile con le frequenze

  23/12/99:  nuova modifica alla lettura dati di ingresso

10/01/2000:  introdotta la possibilita' di definire un fattore di scala per i dati del file d'ingresso  

22/01/2000:  tolto il fattore di scala perche' non funziona

 03/2000:  aggiunte funzioni che restituiscono dati dell'elemento (autovettori, modi ecc.)
             aggiunta la possibilita' di imporre delle deformate modali iniziali
 
   Modifiche fatte al resto del programma:

   Nella classe joint    : aggiunto il vincolo modale
   Nella classe strnode  : aggiunto il nodo modale
   Nella classe MatVec3n : aggiunte classi e funzioni per gestire matrici Nx3 e NxN
   Nella classe SubMat   : corretto un bug nelle funzioni Add Mat3xN ecc.,
                           aggiunte funzioni per trattare sottomatrici Nx3
*/

/* 
 * Copyright 1999-2000 Felice Felippone <ffelipp@tin.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <modal.h>
#include <dataman.h>

Modal::Modal(unsigned int uL,
             const ModalNode* pR,
             const DofOwner* pDO,
	     unsigned int N,        /* numero modi */
	     unsigned int NS,       /* numero nodi d'interfaccia */
             unsigned int NFN,      /* numero nodi FEM */
             doublereal dMassTmp,   /* invarianti d'inerzia (massa, mom. statici e d'inerzia) */ 
             const Vec3& STmp,
             const Mat3x3& JTmp,
             MatNxN *pGenMass,
             MatNxN *pGenStiff,
             MatNxN *pGenDamp,
             unsigned int* IdFemNodes, /* array contenente le label dei nodi FEM */
             unsigned int* IntNodes,   /* array contenente le label dei nodi d'interfaccia */
	     Mat3xN *pN,               /* posizione dei nodi FEM */
             Mat3xN *pXYZOffsetNodes,
             const StructNode** pIN,  /* array di puntatori ai nodi d'interfaccia */
             Mat3xN *pPHItStrNode,    /* forme modali dei nodi d'interfaccia */
             Mat3xN *pPHIrStrNode,
	     Mat3xN *pModeShapest,    /* autovettori: non servono a modal ma a altre classi (es. aeromodal) */
             Mat3xN *pModeShapesr,
             Mat3xN *pInv3,           /* invarianti d'inerzia I3...I11 */  
             Mat3xN *pInv4,
             Mat3xN *pInv5,
             Mat3xN *pInv8,
             Mat3xN *pInv9,
             Mat3xN *pInv10,
             Mat3xN *pInv11,
	     VecN *a,
	     VecN *b,
	     const char *sFileMod,
             DataManager* pDM,  /* non serve */
             MBDynParser& HP,   /* non serve */
	     flag fOut)

: Elem(uL, Elem::JOINT, fOut),
Joint(uL, Joint::MODAL, pDO, fOut),
pModalNode(pR),
NModes(N), 
NStrNodes(NS), 
NFemNodes(NFN), 
dMass(dMassTmp), 
Inv2(STmp), 
Inv7(JTmp),
pModalMass(pGenMass), 
pModalStiff(pGenStiff), 
pModalDamp(pGenDamp),
IdFemNodes(IdFemNodes),
IntNodes(IntNodes),
pXYZFemNodes(pN),
pOffsetNodes(pXYZOffsetNodes), 
pInterfaceNodes(pIN), 
pPHIt(pPHItStrNode), 
pPHIr(pPHIrStrNode),
pModeShapest(pModeShapest),
pModeShapesr(pModeShapesr),
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
Inv5jaj(NModes,0),
Inv4j(0.), 
VInv5jaj(0.), 
VInv5jaPj(0.), 
Inv8j(0.), 
Inv9jkak(0.), 
Inv9jkajaPk(0.),
a(*a), 
b(*b),
bPrime(NModes,0.),
pd1tot(NULL), 
pd2(NULL), 
pR1tot(NULL), 
pR2(NULL), 
pF(NULL), 
pM(NULL),
fOutFlex(sFileMod, std::ios::out /* | std::ios::noreplace */ )
{
   SAFENEWARR(pd1tot, Vec3, NStrNodes);
   SAFENEWARR(pd2, Vec3, NStrNodes);
   SAFENEWARR(pR1tot, Mat3x3, NStrNodes);
   SAFENEWARR(pR2, Mat3x3, NStrNodes);
   SAFENEWARR(pF, Vec3, NStrNodes);
   SAFENEWARR(pM, Vec3, NStrNodes);
   
   if (!fOutFlex) {
      std::cerr << "Modal(" << GetLabel() 
	<< "): unable to open output file \"" << sFileMod 
	<< "\"" << std::endl;
      THROW(ErrGeneric());
   }
   
   for (unsigned int i = 0; i < NStrNodes; i++) {
      pd1tot[i] = Vec3(0.);
      pd2[i] = Vec3(0.);
      pR1tot[i] = Mat3x3(0.);
      pR2[i] = Mat3x3(0.);
      pF[i] = Vec3(0.);
      pM[i] = Vec3(0.);
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
}


Joint::Type 
Modal::GetJointType(void) const
{
   return Joint::MODAL;
}


std::ostream& 
Modal::Restart(std::ostream& out) const      
{
   return out << "not implemented yet" << std::endl;
}


unsigned int 
Modal::iGetNumDof(void) const
{
   /* i gradi di liberta' propri sono: 2*M per i modi
    * 6*N nodi d'interfaccia per le reazioni vincolari */
   return 2*NModes+6*NStrNodes;    
}


DofOrder::Order 
Modal::SetDof(unsigned int i) const
{
   /* gradi di liberta' differenziali (eq. modali) */   
   ASSERT(i < 2*NModes+6*NStrNodes);
   if (i < 2*NModes) {
      return DofOrder::DIFFERENTIAL;

   /* gradi di liberta' algebrici (eq. di vincolo) */
   } /* else if (i >= 2*NModes && i < iGetNumDof()) { */
      return DofOrder::ALGEBRAIC;
   /* } */
   THROW(ErrGeneric());
}


void 
Modal::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
   /* la matrice e' gestita come piena (c'e' un po' di spreco...) */
   *piNumCols = *piNumRows = 12+iGetNumDof()+6*NStrNodes;
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
   WM.ResizeInit(iNumRows, iNumCols, 0.);
   
   /* gli indici sono ordinati cosi': i primi 6 sono le equazioni 
    * per abbassare di grado il sistema,
    * quindi le 6 equazioni del moto rigido, quindi le 2*M modali, 
    * quindi le eq. vincolari   */
   
   /* indici della parte rigida */   
   integer iRigidIndex = pModalNode->iGetFirstIndex();

   for (integer iCnt = 1; iCnt <= 12; iCnt++) {
      WM.fPutRowIndex(iCnt, iRigidIndex+iCnt);
      WM.fPutColIndex(iCnt, iRigidIndex+iCnt);
   }
   
   /* indici della parte deformabile e delle reazioni vincolari */
   integer iFlexIndex = iGetFirstIndex();
   
   for (unsigned int iCnt = 1; iCnt <= iGetNumDof(); iCnt++) {
      WM.fPutRowIndex(12+iCnt, iFlexIndex+iCnt);
      WM.fPutColIndex(12+iCnt, iFlexIndex+iCnt);
   }
   
   /* indici delle equazioni vincolari (solo per il nodo 2) */   
   for (unsigned int iStrNodem1 = 0; iStrNodem1 < NStrNodes; iStrNodem1++) {
      integer iNodeFirstMomIndex = 
	      pInterfaceNodes[iStrNodem1]->iGetFirstMomentumIndex();
      integer iNodeFirstPosIndex = 
	      pInterfaceNodes[iStrNodem1]->iGetFirstPositionIndex();

      integer iOffset = 12+iGetNumDof()+6*iStrNodem1;
      for (integer iCnt = 1; iCnt <= 6; iCnt++) {   
	 WM.fPutRowIndex(iOffset+iCnt, iNodeFirstMomIndex+iCnt);
	 WM.fPutColIndex(iOffset+iCnt, iNodeFirstPosIndex+iCnt);
      }
   }

   /* Assemblaggio dello Jacobiano */
  
   /* recupera e aggiorna i dati necessari */
   /* i prodotti Inv3j*aj ecc. sono gia' stati fatti da AssRes() */
   
   Vec3 wr = pModalNode->GetWRef();
   Mat3x3 wrWedge(wr); 
   Mat3x3 R(pModalNode->GetRRef());      
   Mat3x3 RT = R.Transpose();
   Mat3x3 J = R*(Inv7+Inv8jaj.Symm2())*RT;  
   Vec3 S = R*(Inv2+Inv3jaj);  
   
   /* matrice di massa:        J[1,1] = Mtot  */
   for (int iCnt = 7; iCnt <= 9; iCnt++) {      
      WM.fPutCoef(iCnt, iCnt, dMass);
   }   
   
   /* momenti statici J[1,2] = -[S/\]+c[-2w/\S/\+S/\w/\] */
   Mat3x3 SWedge(S);
   WM.Add(7, 10, ((SWedge*wrWedge-wrWedge*(SWedge*(2.)))*dCoef)-SWedge);
   
   /* J[2,1] = [S/\] */   
   WM.Add(10, 7, SWedge);
   
   /* momenti d'inerzia:       J[2,2] = J+c[-(Jw)/\+(w/\)J];    */   
   Mat3x3 JwWedge(J*wr);
   WM.Add(10, 10, J+((wrWedge*J-JwWedge)*dCoef));
   
   /* completa lo Jacobiano con l'aggiunta delle equazioni {xP} = {v}
    {gP} - [wr/\]{g} = {w} */   
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.fPutCoef(iCnt, iCnt, 1);
      WM.fPutCoef(iCnt, iCnt+6, -dCoef);
   }   
   WM.Add(4, 4, wrWedge*(-dCoef));
   
   /* parte deformabile :
    * 
    * | I  -cI  ||aP|
    * |         ||  |
    * |cK  M+cC ||bP|
    */
   
   for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
      WM.fPutCoef(12+iCnt, 12+iCnt, 1.);
      WM.fPutCoef(12+iCnt, 12+iCnt+NModes, -dCoef);
      for (unsigned int jCnt = 1; jCnt <= NModes; jCnt++) {
	 WM.fPutCoef(12+NModes+iCnt, 12+jCnt, dCoef*pModalStiff->dGet(iCnt, jCnt));
	 WM.fPutCoef(12+NModes+iCnt, 12+NModes+jCnt, pModalMass->dGet(iCnt, jCnt)+dCoef*pModalDamp->dGet(iCnt, jCnt));   
      }
   }
   
 /* termini di accoppiamento moto rigido-deformabile; eventualmente 
  * l'utente potra' scegliere se trascurarli tutti, una parte o considerarli 
  * tutti */
   
  /* linearizzazione delle OmegaPrime: 
   * J13 = R*Inv3jaj 
   * J23 = R*Inv4+Inv5jaj
   * (questi termini ci vogliono sempre)
   */
   
   Mat3xN Jac13(NModes,0.), Jac23(NModes,0.), Inv5jajRef(NModes, 0.);
   MatNx3 Jac13T(NModes, 0.), Jac23T(NModes, 0.);
   Jac13.LeftMult(R, *pInv3);
   Jac23.LeftMult(R, *pInv4);
   Jac23 += Inv5jajRef.LeftMult(R, Inv5jaj);
   
   WM.Add(7, 12+NModes+1, Jac13); 
   WM.Add(10, 12+NModes+1, Jac23);
   WM.Add(12+NModes+1, 7, Jac13T.Transpose(Jac13)); 
   WM.Add(12+NModes+1, 10, Jac23T.Transpose(Jac23));

  /* termini di Coriolis: linearizzazione delle Omega 
   * (si puo' evitare se non ho grosse vel. angolari):
   * J13 = -2*R*[Inv3jaPj/\]*RT
   * J23 = 2*R*[Inv8jaPj-Inv9jkajaPk]*RT */

   Mat3x3 Inv3jaPjWedge(Inv3jaPj);
   WM.Add(7, 10, R*Inv3jaPjWedge*RT*dCoef*-2);
   WM.Add(10,10, R*( /* ( */ Inv8jaPj /* -Inv9jkajaPk ) */ *(RT*(2.*dCoef))));
   
   /* termini di Coriolis: linearizzazione delle b;
    * si puo' evitare 'quasi' sempre: le vel. di deformazione dovrebbero
    * essere sempre piccole rispetto ai termini rigidi
    * Jac13 = 2*[Omega/\]*R*PHI
    */
#if 0
   Jac13.LeftMult(wrWedge*R*2*dCoef, *pInv3);
   WM.Add(7, 12+NModes+1, Jac13);
#endif
   /* nota: non risco a tirar fuori la Omega dall'eq. dei momenti:
    * 2*[ri/\]*[Omega/\]*R*PHIi*{DeltaaP},  i = 1,...nnodi 
    * quindi questa equazione non si puo' linearizzare
    * a meno di ripetere la sommatoria estesa a tutti i nodi a ogni passo 
    * del Newton-Rapson... */
  
   /* linearizzazione delle forze centrifughe e di Coriolis in base modale 
    * (anche questi termini dovrebbero essere trascurabili) 
    */
#if 0
   for (int iCnt = 1; iCnt<=3; iCnt++) { 
      double temp1 = 0., temp2 = 0.;
      for (int jCnt = 1; jCnt<=3; jCnt++) {
	 temp1 += -2*wr.dGet(jCnt)*(R*(Inv8j.Transpose()-Inv9jkak)*RT).dGet(iCnt, jCnt);
	 temp2 += -(R*(Inv4j+VInv5jaj)).dGet(jCnt)*wrWedge.dGet(iCnt, jCnt);
      }
      WM.fIncCoef(12+NModes+1, 9+iCnt, dCoef*(temp1+temp2));
   }
   
   for (int iCnt = 1; iCnt<=3; iCnt++) { 
      double temp1 = 0.;
      for (int jCnt = 1; jCnt<=3; jCnt++) {
	 temp1 += (R*VInv5jaPj*2).dGet(jCnt);
      }
      WM.fIncCoef(12+NModes+1, 9+iCnt, dCoef*temp1);
   } 
#endif
   
   /* ciclo esteso a tutti i nodi d'interfaccia */
   for (unsigned int iStrNode = 1; iStrNode <= NStrNodes; iStrNode++) {   
      unsigned int iStrNodem1 = iStrNode - 1;
      
      /* recupero le forme modali del nodo vincolato */
      Mat3xN PHIt(NModes), PHIr(NModes);
      for (int iCnt = 1; iCnt <= 3; iCnt++) {   
	 for (unsigned int iMode = 1; iMode <= NModes; iMode++) {
	    PHIt.Put(iCnt, iMode, pPHIt->dGet(iCnt, (iMode-1)*NStrNodes+iStrNode));
	    PHIr.Put(iCnt, iMode, pPHIr->dGet(iCnt, (iMode-1)*NStrNodes+iStrNode));
	 }
      }
      
      MatNx3 PHItT(NModes), PHIrT(NModes);
      PHItT.Transpose(PHIt);
      PHIrT.Transpose(PHIr);
      
      /* nota: le operazioni d1tot = d1+PHIt*a, R1tot = R*[I+(PHIr*a)/\] 
       * sono gia' state fatte da AssRes */
      
      Mat3xN SubMat1(NModes), SubMat2(NModes);
      MatNx3 SubMat3(NModes); 
      MatNxN SubMat4(NModes);
      
      /* cerniera sferica */
      /* F e' aggiornata da AssRes */
      
      integer iReactionIndex = 12+2*NModes+6*iStrNodem1;
      integer iStrNodeIndex = 12+iGetNumDof()+6*iStrNodem1;

      /* termini di reazione sui nodi */      
      for (int iCnt = 1; iCnt <= 3; iCnt++) { 
	 WM.fIncCoef(iCnt+6, iReactionIndex+iCnt, 1.);
	 WM.fIncCoef(iStrNodeIndex+iCnt, iReactionIndex+iCnt, -1.);
      }

      Vec3 dTmp1(R*pd1tot[iStrNodem1]); 
      /* ppd1Tot e' il puntatore all'array 
       * che contiene le posizioni del nodi FEM */ 

      Mat3x3 dTmp1Wedge(dTmp1); 
      WM.Add(10, iReactionIndex+1, dTmp1Wedge);
      
      /* idem per pd2 e pR2 */
      Vec3 dTmp2(pR2[iStrNodem1]*pd2[iStrNodem1]);
      
      Mat3x3 dTmp2Wedge(dTmp2);
      WM.Add(iStrNodeIndex+4, iReactionIndex+1, -dTmp2Wedge);
      
      Mat3x3 FTmpWedge(pF[iStrNodem1]*dCoef);
      
      /* termini del tipo: c*F/\*d/\*Deltag */
      
      WM.Add(10, 4, FTmpWedge*dTmp1Wedge);
      WM.Add(iStrNodeIndex+4, iStrNodeIndex+4, -FTmpWedge*dTmp2Wedge);
      
      /* termini aggiuntivi dovuti alla deformabilita' */
      
      SubMat3.RightMult(PHItT, RT*FTmpWedge);
      WM.Add(12+NModes+1, 4, SubMat3);
      SubMat3.RightMult(PHItT, RT);
      WM.Add(12+NModes+1, iReactionIndex+1, SubMat3);
      SubMat1.LeftMult(-FTmpWedge, PHIt);
      WM.Add(10, 13, SubMat1);
      
      /* modifica: divido le equazioni di vincolo per dCoef */
      
      /* termini di vincolo dovuti ai nodi */
      for (int iCnt = 1; iCnt <= 3; iCnt++) {
	 WM.fIncCoef(iReactionIndex+iCnt, iCnt, -1.);
	 WM.fIncCoef(iReactionIndex+iCnt, iStrNodeIndex+iCnt, 1);
      }
      WM.Add(iReactionIndex+1, 4, dTmp1Wedge);
      
      /* contributo dovuto alla flessibilita' */
      SubMat1.LeftMult(-R, PHIt);
      WM.Add(iReactionIndex+1, 13, SubMat1); 
      
      /*termini di vincolo dovuti al nodo 2 */
      WM.Add(iReactionIndex+1, iStrNodeIndex+4, -dTmp2Wedge);
      
      /* fine cerniera sferica */
      
      /* equazioni di vincolo : giunto prismatico */
      
      /* Vec3 M(XCurr, iModalIndex+2*NModes+6*iStrNodem1+4); */
      Vec3 MTmp = pM[iStrNodem1]*dCoef; 

      Mat3x3 R1totTranspose = pR1tot[iStrNodem1].Transpose();
      
      Vec3 e1tota(pR1tot[iStrNodem1].GetVec(1));
      Vec3 e2tota(pR1tot[iStrNodem1].GetVec(2));
      Vec3 e3tota(pR1tot[iStrNodem1].GetVec(3));
      Vec3 e1b(pR2[iStrNodem1].GetVec(1));
      Vec3 e2b(pR2[iStrNodem1].GetVec(2));
      Vec3 e3b(pR2[iStrNodem1].GetVec(3));
      
      Mat3x3 MWedge(Mat3x3(e3b, e2tota*MTmp.dGet(1))
		    +Mat3x3(e1b, e3tota*MTmp.dGet(2))
		    +Mat3x3(e2b, e1tota*MTmp.dGet(3)));
      Mat3x3 MWedgeT(MWedge.Transpose());
      
      /* Eq. dei momenti (termini del tipo [e3b/\][e2a/\]M) */      
      WM.Add(10, 4, MWedge);
      WM.Sub(10, iStrNodeIndex+4, MWedgeT);
      
      WM.Sub(iStrNodeIndex+4, 4, MWedge);
      WM.Add(iStrNodeIndex+4, iStrNodeIndex+4, MWedgeT);
      
      /* Eq. dei momenti, contributo modale */      
      Vec3 e1ra(R.GetVec(1));
      Vec3 e2ra(R.GetVec(2));
      Vec3 e3ra(R.GetVec(3));
      Mat3x3 M1Wedge(Mat3x3(e3b, e2ra*MTmp.dGet(1))
		     +Mat3x3(e1b, e3ra*MTmp.dGet(2))
		     +Mat3x3(e2b, e1ra*MTmp.dGet(3)));
      
      SubMat1.LeftMult(M1Wedge, PHIr);
      WM.Add(10, 13, SubMat1);
      WM.Sub(iStrNodeIndex+4, 13, SubMat1);
      
      /* Eq. d'equilibrio ai modi */      
      SubMat3.RightMult(PHIrT, R1totTranspose*MWedge);
      WM.Add(12+NModes+1, 4, SubMat3);
      SubMat3.RightMult(PHIrT, R1totTranspose*-MWedgeT);	// FIXME ?!?
      WM.Add(12+NModes+1, iStrNodeIndex+1, SubMat3);
      
      Vec3 MA(Mat3x3(e2tota.Cross(e3b), e3tota.Cross(e1b), e1tota.Cross(e2b))*(pM[iStrNodem1]*dCoef));
      Mat3x3 MAWedge(MA);
      SubMat3.RightMult(PHIrT, R1totTranspose*MAWedge);
      WM.Add(12+NModes+1, 4, SubMat3);
      
      Mat3x3 R1TMAWedge(RT*MA);
      SubMat1.LeftMult(R1TMAWedge, PHIr);   
      SubMat2.LeftMult(R1totTranspose*M1Wedge, PHIr);
      SubMat1 += SubMat2;
      SubMat4.Mult(PHIrT, SubMat1);
      for (unsigned int iMode = 1; iMode <= NModes; iMode++) {
	 for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
	    WM.fIncCoef(12+NModes+iMode, 12+jMode, SubMat4.dGet(iMode, jMode));
	 }
      }
      
      /* Termine e2a/\e3b + e3a/\e1b + e1a/\e2b */      
      Vec3 v1(e2tota.Cross(e3b));
      Vec3 v2(e3tota.Cross(e1b));
      Vec3 v3(e1tota.Cross(e2b));
      
      MWedge = Mat3x3(v1, v2, v3);
      
      WM.Add(10, iReactionIndex+4, MWedge);
      WM.Sub(iStrNodeIndex+4, iReactionIndex+4, MWedge);
      
      /* contributo modale: PHIrT*R1T*(e2a/\e3b + e3a/\e1b + e1a/\e2b)  */
      SubMat3.RightMult(PHIrT, R1totTranspose*MWedge);
      WM.Add(12+NModes+1, iReactionIndex+4, SubMat3);
      
      /* Modifica: divido le equazioni di vincolo per dCoef */
      MWedge = MWedge.Transpose();
      
      /* Eq. di vincolo */
      WM.Add(iReactionIndex+4, 4, MWedge);
      WM.Sub(iReactionIndex+4, iStrNodeIndex+4, MWedge);
      
      /* Eq. di vincolo, termine aggiuntivo modale */
      Vec3 u1(e2ra.Cross(e3b));
      Vec3 u2(e3ra.Cross(e1b));
      Vec3 u3(e1ra.Cross(e2b));
      
      M1Wedge = (Mat3x3(u1, u2, u3)).Transpose();
      SubMat1.LeftMult(M1Wedge, PHIr);
      WM.Add(iReactionIndex+4, 13, SubMat1); 
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
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);

   /*
    * 1  		-> 12:				rigid body
    *
    * 13 		-> 12 + 2*NM:			modes
    *
    * 12 + 2*NM 	-> 12 + 2*NM + 6*NN:		node reactions
    *
    * 12 + 2*NM + 6*NN	-> 12 + 2*NM + 6*NN + 6*NN:	nodes
    */
   
   /* rigid body indices */
   integer iRigidIndex = pModalNode->iGetFirstIndex();
   for (unsigned int iCnt = 1; iCnt <= 12; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iRigidIndex+iCnt);
   }
  
   /* modal dofs indices */
   integer iModalIndex = iGetFirstIndex();
   for (unsigned int iCnt = 1; iCnt <= iGetNumDof(); iCnt++) {
      WorkVec.fPutRowIndex(12+iCnt, iModalIndex+iCnt);
   }
  
   /* interface nodes indices */
   for (unsigned int iStrNodem1 = 0; iStrNodem1 < NStrNodes; iStrNodem1++) {
      integer iNodeFirstMomIndex = 
	      pInterfaceNodes[iStrNodem1]->iGetFirstMomentumIndex();
      
      for (unsigned int iCnt = 1; iCnt <= 6; iCnt++) { 
	 WorkVec.fPutRowIndex(12+iGetNumDof()+6*iStrNodem1+iCnt, 
			 iNodeFirstMomIndex+iCnt);
      }
   }
   
   /* recupera i dati necessari */   
   Vec3 x = (pModalNode->GetXCurr());
   Vec3 xP = pModalNode->GetVCurr();
   Vec3 g = pModalNode->GetgCurr();
   Vec3 gP = pModalNode->GetgPCurr();
   Vec3 v(XCurr, iRigidIndex+7);
   Vec3 vP = pModalNode->GetXPPCurr(); 
   Vec3 w(XCurr, iRigidIndex+10); 
   Vec3 wr = pModalNode->GetWRef();     
   Vec3 wP = pModalNode->GetWPCurr();

   a.Copy(XCurr, iModalIndex + 1);
   VecN aPrime(XPrimeCurr, NModes, iModalIndex + 1);
   b.Copy(XCurr, iModalIndex+NModes+1);
   bPrime.Copy(XPrimeCurr, iModalIndex + NModes + 1);
   
   Mat3x3 R(pModalNode->GetRCurr());
   Mat3x3 RT = R.Transpose();
   
   VecN Ka(NModes), CaP(NModes), MaPP(NModes);
   
   /* aggiorna gli invarianti */   
   Ka.Mult(*pModalStiff, a);
   CaP.Mult(*pModalDamp, b);
   MaPP.Mult(*pModalMass, bPrime);
   
   Inv3jaj = *pInv3*a;
   Inv3jaPj = *pInv3*b;
   Vec3 Inv3jaPPj = *pInv3*bPrime;
   
   /* Aggiorna gli invarianti rotazionali */
   Vec3 Inv11jaPj;
   Inv11jaPj = R*(*pInv11*b);
   
   Inv8jaj = 0.; 
   Inv8jaPj = 0.;  
   Inv5jaj.Reset();
   
   Mat3xN Inv5jaPj(NModes, 0.);
   Mat3x3 MatTmp1(0.), MatTmp2(0.);
#ifdef MODAL_USE_INV9
   Mat3x3 Inv9jkajak(0.);
#endif /* MODAL_USE_INV9 */
   Mat3x3 Inv10jaPj(0.);
   
   for (unsigned int iMode = 1; iMode <= NModes; iMode++)  {
      Mat3x3 Inv8jajTmp;
      Mat3x3 Inv10jaPjTmp;

      doublereal a_iMode = a.dGet(iMode);
      doublereal aP_iMode = b.dGet(iMode);
      
      for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
	 for (unsigned int jCnt = 1; jCnt <= 3; jCnt++) { 
	    Inv8jajTmp.Put(iCnt, jCnt, pInv8->dGet(iCnt, (iMode-1)*3+jCnt));
	    Inv10jaPjTmp.Put(iCnt, jCnt, 
			    pInv10->dGet(iCnt, (iMode-1)*3+jCnt)*aP_iMode);
	 }

	 for (unsigned int jMode = 1; jMode <= NModes; jMode++) {
	    doublereal d = pInv5->dGet(iCnt, (iMode-1)*NModes+jMode);

	    Inv5jaj.Add(iCnt, jMode, d*a_iMode);
	    Inv5jaPj.Add(iCnt, jMode, d*aP_iMode);
	 }
      }


      Inv8jaj += Inv8jajTmp*a_iMode;
      Inv8jaPj += Inv8jajTmp*aP_iMode; 
      Inv10jaPj += Inv10jaPjTmp;
      
#ifdef MODAL_USE_INV9
      /* questi termini si possono commentare perchè sono (sempre ?) piccoli 
       * (termini del tipo a*a o a*b)
       * eventualmente dare all'utente la possibilita' di scegliere 
       * se trascurarli o no */
      for (int kMode = 1; kMode <= NModes; kMode++)  {
	 doublereal a_kMode = a.dGet(kMode);
	 doublereal aP_kMode = b.dGet(kMode);
	 unsigned int iOffset = (iMode-1)*3*NModes+(kMode-1)*3;

	 for (int iCnt = 1; iCnt <= 3; iCnt++) {
	    for (int jCnt = 1; jCnt <= 3; jCnt++) {
	       MatTmp1.Put(iCnt, jCnt, pInv9->dGet(iCnt, iOffset+jCnt)*a_iMode*a_kMode);
	       MatTmp2.Put(iCnt, jCnt, pInv9->dGet(iCnt, iOffset+jCnt)*a_iMode*aP_kMode);
	    }
	 }
	 Inv9jkajak  += MatTmp1;
	 Inv9jkajaPk += MatTmp2; 
      }
#endif /* MODAL_USE_INV9 */
      
   } /* fine ciclo sui modi */
   
   Mat3x3 J = R*(Inv7+Inv8jaj.Symm2()
#ifdef MODAL_USE_INV9
		   -Inv9jkajak
#endif /* MODAL_USE_INV9 */ 
		   )*RT;
   Vec3 S = R*(Inv2+Inv3jaj);
   
   Mat3xN Inv4Curr(NModes, 0), Inv5jajCurr(NModes, 0);
   Inv4Curr.LeftMult(R, *pInv4);
   Inv5jajCurr.LeftMult(R, Inv5jaj);
   
   /* fine aggiornamento invarianti */
   
   /* forze d'inerzia */
   WorkVec.Sub(7, vP*dMass-S.Cross(wP)+w.Cross(w.Cross(S))+(w.Cross(R*Inv3jaPj))*2.+R*Inv3jaPPj);

#if 0
   std::cerr << "m*vP=" << vP*dMass << "; R*I3*aPP=" << R*Inv3jaPPj << "; tot=" 
   	<< vP*dMass-S.Cross(wP)+w.Cross(w.Cross(S))+(w.Cross(R*Inv3jaPj))*2+R*Inv3jaPPj << std::endl;
#endif

   WorkVec.Sub(10, S.Cross(vP)+J*wP+w.Cross(J*w)
		   +(R*( /* ( */ Inv8jaPj /* -Inv9jkajaPk) */ *(RT*w)))*2.
		   +Inv4Curr*bPrime+Inv5jajCurr*bPrime);
   
   /* termini dovuti alle inerzie rotazionali */
   WorkVec.Sub(10, R*(Inv10jaPj.Symm2()*(RT*w))+w.Cross(Inv11jaPj));
   
   /* forza di gravita' (decidere come inserire g) */
#ifdef MODAL_USE_GRAVITY
   WorkVec.Add(7, g*dMass);
   WorkVec.Add(10, S.Cross(g));
#endif /* MODAL_USE_GRAVITY */
    
   /* forze modali */
   Vec3 Inv3j; 
   Mat3x3 Inv10j;
   for (unsigned int iMode = 1; iMode <= NModes; iMode++) { 
      for (int iCnt = 1; iCnt <= 3; iCnt++) {
	 Inv3j.Put(iCnt, pInv3->dGet(iCnt, iMode));
	 Inv4j.Put(iCnt, pInv4->dGet(iCnt, iMode));
	 VInv5jaj.Put(iCnt, Inv5jaj.dGet(iCnt, iMode));
	 VInv5jaPj.Put(iCnt, Inv5jaPj.dGet(iCnt, iMode));
	 for (unsigned int jCnt = 1; jCnt <= 3; jCnt++) {
	    Inv8j.Put(iCnt, jCnt, pInv8->dGet(iCnt, (iMode-1)*3+jCnt));
	    Inv10j.Put(iCnt, jCnt, pInv10->dGet(iCnt, (iMode-1)*3+jCnt));
	 }
      }

      Inv9jkak = 0.;

      for (unsigned int kModem1 = 0; kModem1 < NModes; kModem1++)  {
	 doublereal a_kMode = a.dGet(kModem1+1);

	 for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
	    for (unsigned int jCnt = 1; jCnt <= 3; jCnt++) 
	      MatTmp1.Put(iCnt, jCnt, pInv9->dGet(iCnt,(iMode-1)*3*NModes+kModem1*3+jCnt)*a_kMode);
	 }
	 Inv9jkak += MatTmp1;
      }

      WorkVec.fIncCoef(12+NModes+iMode, 
      		-(R*Inv3j).Dot(vP)-(R*(Inv4j+VInv5jaj)).Dot(wP)
		+w.Dot(R*((Inv8j.Transpose()-Inv9jkak+Inv10j)*(RT*w)))
		-(R*VInv5jaPj).Dot(w)*2.
		-MaPP.dGet(iMode)-CaP.dGet(iMode)-Ka.dGet(iMode));
      
      /* forza di gravita': */
#ifdef MODAL_USE_GRAVITY
      WorkVec.fIncCoef(12+NModes+iMode, (R*Inv3j).Dot(g));
#endif /* MODAL_USE_GRAVITY */
   }
   
   /* equazioni per abbassare di grado il sistema */
   WorkVec.Add(1, v-xP);
   WorkVec.Add(4, w-Mat3x3(MatG, g)*gP-Mat3x3(MatR, g)*wr);
   for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
      WorkVec.fPutCoef(12+iCnt, b.dGet(iCnt)-aPrime.dGet(iCnt));
   }
   
   /* equazioni di vincolo */
   for (unsigned int iStrNode = 1; iStrNode <= NStrNodes; iStrNode++) {
      unsigned int iStrNodem1 = iStrNode - 1;
      integer iReactionIndex = 12+2*NModes+6*iStrNodem1;
      integer iStrNodeIndex = iReactionIndex + 6*NStrNodes;
      
      /* nodo 1 (FEM) e 2 (Multi - Body): recupera la posizione */
      Vec3 d1rig(pOffsetNodes->dGet(1, iStrNode), 
		      pOffsetNodes->dGet(2, iStrNode), 
		      pOffsetNodes->dGet(3, iStrNode));
      
      pd2[iStrNodem1] = Vec3(pOffsetNodes->dGet(1, NStrNodes+iStrNode),
		      pOffsetNodes->dGet(2, NStrNodes+iStrNode), 
		      pOffsetNodes->dGet(3, NStrNodes+iStrNode));
      
      /* recupero le forme modali del nodo vincolato */
      Mat3xN PHIt(NModes), PHIr(NModes);
      for (int iCnt = 1; iCnt <= 3; iCnt++) {   
	 for (unsigned int iMode = 1; iMode <= NModes; iMode++) {
	    PHIt.Put(iCnt, iMode, pPHIt->dGet(iCnt, (iMode-1)*NStrNodes+iStrNode));
	    PHIr.Put(iCnt, iMode, pPHIr->dGet(iCnt, (iMode-1)*NStrNodes+iStrNode));
	 }
      }
      
      MatNx3 PHItT(NModes), PHIrT(NModes);
      PHItT.Transpose(PHIt);
      PHIrT.Transpose(PHIr);
      
      /*
       * aggiorno d1 e R1 con il contributo dovuto alla flessibilita': 
       * d1tot = d1+PHIt*a, R1tot = R*[I+(PHIr*a)/\]
       */
      pd1tot[iStrNodem1] = d1rig+PHIt*a;
      pR1tot[iStrNodem1] = R*Mat3x3(1., PHIr*a);
     
      /* constraint reaction (force) */
      pF[iStrNodem1] = Vec3(XCurr, iModalIndex+2*NModes+6*iStrNodem1+1);
      Vec3 x2 = pInterfaceNodes[iStrNodem1]->GetXCurr();
      pR2[iStrNodem1] = pInterfaceNodes[iStrNodem1]->GetRCurr();
      
      /* cerniera sferica */
      Vec3 dTmp1(R*pd1tot[iStrNodem1]);
      Vec3 dTmp2(pR2[iStrNodem1]*pd2[iStrNodem1]);
      
      /* Eq. d'equilibrio, nodo 1 */
      WorkVec.Sub(7, pF[iStrNodem1]);
      WorkVec.Sub(10, dTmp1.Cross(pF[iStrNodem1]));
      
      /*termine aggiuntivo dovuto alla deformabilita': -PHItiT*RT*F */
      Vec3 vtemp = RT*pF[iStrNodem1];
      for (unsigned int iMode = 1; iMode <= NModes; iMode++) {
	 doublereal temp = 0.;

	 for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
	    temp += PHItT.dGet(iMode, iCnt)*vtemp.dGet(iCnt); 
	 }

	 WorkVec.fIncCoef(12+NModes+iMode, -temp);
      }
      /* Eq. d'equilibrio, nodo 2 */
      WorkVec.Add(iStrNodeIndex+1, pF[iStrNodem1]);
      WorkVec.Add(iStrNodeIndex+4, dTmp2.Cross(pF[iStrNodem1]));
      
      /* Eq. di vincolo */
      if (dCoef != 0.) {
	 WorkVec.Add(iReactionIndex+1, (x+dTmp1-x2-dTmp2)/dCoef);
      }
      
      /* giunto prismatico */
      pM[iStrNodem1] = Vec3(XCurr, iModalIndex+2*NModes+6*iStrNodem1+4);
      
      Vec3 e1a(pR1tot[iStrNodem1].GetVec(1));
      Vec3 e2a(pR1tot[iStrNodem1].GetVec(2));
      Vec3 e3a(pR1tot[iStrNodem1].GetVec(3));
      Vec3 e1b(pR2[iStrNodem1].GetVec(1));
      Vec3 e2b(pR2[iStrNodem1].GetVec(2));
      Vec3 e3b(pR2[iStrNodem1].GetVec(3));
      
      Vec3 MTmp(Mat3x3(e2a.Cross(e3b), e3a.Cross(e1b), e1a.Cross(e2b))*pM[iStrNodem1]);

      /* Equazioni di equilibrio, nodo 1 */
      WorkVec.Sub(10, MTmp); 
      
      /* Equazioni di equilibrio, nodo 2 */
      WorkVec.Add(iStrNodeIndex+4, MTmp);
      
      /* Contributo dovuto alla flessibilita' :-PHIrT*RtotT*M */
      vtemp = pR1tot[iStrNodem1].Transpose()*MTmp;
      for (unsigned int iMode = 1; iMode <= NModes; iMode++) {
	 doublereal temp = 0;

	 for (int iCnt = 1; iCnt <= 3; iCnt++) {
	    temp += PHIrT.dGet(iMode, iCnt)*vtemp.dGet(iCnt);
	 } 

	 WorkVec.fIncCoef(12+NModes+iMode, -temp);
      }
      
      /* Modifica: divido le equazioni di vincolo per dCoef */
      if (dCoef != 0.) {
	 /* Equazioni di vincolo di rotazione */
	 WorkVec.fIncCoef(iReactionIndex+4, -(e3b.Dot(e2a)/dCoef));
	 WorkVec.fIncCoef(iReactionIndex+5, -(e1b.Dot(e3a)/dCoef));
	 WorkVec.fIncCoef(iReactionIndex+6, -(e2b.Dot(e1a)/dCoef));
      }   
   }
   
   return WorkVec;
}


void 
Modal::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) {
      /* stampa sul file di output i modi */       
      
      for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
	 fOutFlex << " " << iCnt
	   << " " << a.dGet(iCnt)
	     << " " << b.dGet(iCnt)
	       << " " << bPrime.dGet(iCnt) << std::endl;
      }
   }
}


unsigned int 
Modal::iGetInitialNumDof(void) const
{
   return 2*NModes+12*NStrNodes;
}


void 
Modal::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
   *piNumRows = *piNumCols = 12+iGetInitialNumDof()+12*NStrNodes;
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
   WM.ResizeInit(iNumRows, iNumCols, 0.);
   
   integer iRigidIndex = pModalNode->iGetFirstIndex();
   
   for (unsigned int iCnt = 1; iCnt <= 12; iCnt++) {
      WM.fPutRowIndex(iCnt, iRigidIndex+iCnt);
      WM.fPutColIndex(iCnt, iRigidIndex+iCnt);
   }
   
   integer iFlexIndex = iGetFirstIndex();
   
   for (unsigned int iCnt = 1; iCnt <= iGetInitialNumDof(); iCnt++) {
      WM.fPutRowIndex(12+iCnt, iFlexIndex+iCnt);
      WM.fPutColIndex(12+iCnt, iFlexIndex+iCnt);
   } 
   
   for (unsigned int iStrNodem1 = 0; iStrNodem1 < NStrNodes; iStrNodem1++) {
      integer iNodeFirstPosIndex = pInterfaceNodes[iStrNodem1]->iGetFirstPositionIndex();
      integer iNodeFirstVelIndex = iNodeFirstPosIndex+6;
      
      for (unsigned int iCnt = 1; iCnt <= 6; iCnt++) {   
	 WM.fPutRowIndex(12+iGetInitialNumDof()+12*iStrNodem1+iCnt, 
			 iNodeFirstPosIndex+iCnt);
	 WM.fPutColIndex(12+iGetInitialNumDof()+12*iStrNodem1+iCnt, 
			 iNodeFirstPosIndex+iCnt);
	 WM.fPutRowIndex(12+iGetInitialNumDof()+12*iStrNodem1+6+iCnt, 
			 iNodeFirstVelIndex+iCnt);
	 WM.fPutColIndex(12+iGetInitialNumDof()+12*iStrNodem1+6+iCnt, 
			 iNodeFirstVelIndex+iCnt);
      }
   }
   
   /* comincia l'assemblaggio dello Jacobiano */
   
   /* nota: nelle prime 12+2M equazioni metto dei termini piccoli per evitare singolarita' 
    * quando non ho vincoli esterni o ho azzerato i modi */
   for (unsigned int iCnt = 1; iCnt <= 12+2*NModes; iCnt++) {
      WM.fPutCoef(iCnt, iCnt, 1.e-12);
   }
   
   /* forze elastiche e viscose */
   for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
      for (unsigned int jCnt = 1; jCnt <= NModes; jCnt++) {
	 WM.fPutCoef(12+iCnt, 12+jCnt, pModalStiff->dGet(iCnt,jCnt));
	 WM.fPutCoef(12+iCnt, 12+NModes+jCnt, pModalDamp->dGet(iCnt, jCnt));
      }
   }
   
   /* equazioni di vincolo */
   
   /* cerniera sferica */
   for (unsigned int iStrNode = 1; iStrNode <= NStrNodes; iStrNode++) {
      unsigned int iStrNodem1 = iStrNode - 1;
      
      /* recupera i dati */
      Mat3x3 R1(pModalNode->GetRRef());
      Mat3x3 R2(pInterfaceNodes[iStrNodem1]->GetRRef());
      Vec3 Omega1(pModalNode->GetWRef());
      Vec3 Omega2(pInterfaceNodes[iStrNodem1]->GetWRef());
      Vec3 F(XCurr, iFlexIndex+2*NModes+12*iStrNodem1+1);
      Vec3 FPrime(XCurr, iFlexIndex+2*NModes+12*iStrNodem1+7);
      
      Mat3xN PHIt(NModes,0), PHIr(NModes, 0);
      MatNx3 PHItT(NModes, 0), PHIrT(NModes, 0);
      
      for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
	 for (unsigned int iMode = 1; iMode <= NModes; iMode++) {
	    PHIt.Put(iCnt, iMode, pPHIt->dGet(iCnt, (iMode-1)*NStrNodes+iStrNode));
	    PHIr.Put(iCnt, iMode, pPHIr->dGet(iCnt, (iMode-1)*NStrNodes+iStrNode));
	 }
      }
      PHItT.Transpose(PHIt); 
      PHIrT.Transpose(PHIr);
      
      Vec3 d1rig(pOffsetNodes->dGet(1, iStrNode),
		      pOffsetNodes->dGet(2, iStrNode),
		      pOffsetNodes->dGet(3, iStrNode));
      Vec3 d2(pOffsetNodes->dGet(1, NStrNodes+iStrNode),
		      pOffsetNodes->dGet(2, NStrNodes+iStrNode),
		      pOffsetNodes->dGet(3, NStrNodes+iStrNode)); 
      
      Vec3   d1tot = d1rig+PHIt*a;
      Mat3x3 R1tot = R1*Mat3x3(1., PHIr*a);   
      Mat3xN SubMat1(NModes, 0.);
      MatNx3 SubMat2(NModes, 0.);
      
      for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
	 /* Contributo di forza all'equazione della forza, nodo 1 */
	 WM.fIncCoef(iCnt, 12+2*NModes+12*iStrNodem1+iCnt, 1.);
	 
	 /* Contrib. di der. di forza all'eq. della der. della forza, nodo 1 */
	 WM.fIncCoef(6+iCnt, 12+2*NModes+12*iStrNodem1+6+iCnt, 1.);
	 
	 /* Contributo di forza all'equazione della forza, nodo 2 */
	 WM.fIncCoef(12+iGetInitialNumDof()+12*iStrNodem1+iCnt, 12+2*NModes+12*iStrNodem1+iCnt, -1.);
	 
	 /* Contrib. di der. di forza all'eq. della der. della forza, nodo 2 */
	 WM.fIncCoef(12+iGetInitialNumDof()+12*iStrNodem1+6+iCnt, 12+2*NModes+12*iStrNodem1+6+iCnt, -1.);
	 
	 /* Equazione di vincolo, nodo 1 */
	 WM.fIncCoef(12+2*NModes+12*iStrNodem1+iCnt, iCnt, -1.);
	 
	 /* Derivata dell'equazione di vincolo, nodo 1 */
	 WM.fIncCoef(12+2*NModes+12*iStrNodem1+6+iCnt, 6+iCnt, -1.);
	 
	 /* Equazione di vincolo, nodo 2 */
	 WM.fIncCoef(12+2*NModes+12*iStrNodem1+iCnt, 12+iGetInitialNumDof()+12*iStrNodem1+iCnt, 1.);
	 
	 /* Derivata dell'equazione di vincolo, nodo 2 */
	 WM.fIncCoef(12+2*NModes+12*iStrNodem1+6+iCnt, 12+iGetInitialNumDof()+12*iStrNodem1+6+iCnt, 1.);
      }
      
      /* Distanza nel sistema globale */
      Vec3 d1Tmp(R1*d1tot);
      Vec3 d2Tmp(R2*d2);
      
      /* Matrici F/\d1/\, -F/\d2/\ , F/\omega1/\ */
      Mat3x3 FWedged1Wedge(F, d1Tmp);
      Mat3x3 FWedged2Wedge(F, -d2Tmp);
      Mat3x3 FWedgeO1Wedge(F, Omega1);
      
      /* Matrici (omega1/\d1)/\, -(omega2/\d2)/\ */
      Mat3x3 O1Wedged1Wedge(Omega1.Cross(d1Tmp));
      Mat3x3 O2Wedged2Wedge(d2Tmp.Cross(Omega2));
      
      /* d1Prime= w1/\d1 + R*PHIt*b */
      Mat3xN R1PHIt(NModes);
      R1PHIt.LeftMult(R1, PHIt);
      Vec3 d1Prime(Omega1.Cross(d1Tmp)+R1PHIt*b);
      
      /* Equazione di momento, nodo 1 */
      WM.Add(4, 4, FWedged1Wedge);
      WM.Add(4, 12+2*NModes+12*iStrNodem1+1, Mat3x3(d1Tmp));
      
      /* Equazione di momento, nodo 2 */
      WM.Add(12+iGetInitialNumDof()+12*iStrNodem1+4,12+iGetInitialNumDof()+12*iStrNodem1+4,  FWedged2Wedge);
      WM.Add(12+iGetInitialNumDof()+12*iStrNodem1+4,12+2*NModes+12*iStrNodem1+1, Mat3x3(-d2Tmp));
      
      /* Equazione di momento, contributo modale */
      Mat3x3 FWedge(F);
      SubMat1.LeftMult(-FWedge*R1, PHIt);
      WM.Add(4, 13, SubMat1);
      
      /* Eq. di equilibrio ai modi */
      SubMat2.RightMult(PHItT, R1.Transpose()*FWedge);
      WM.Add(13, 4, SubMat2);
      SubMat2.RightMult(PHItT, R1.Transpose());
      WM.Add(13, 12+2*NModes+1, SubMat2);
      
      /* derivata dell'equazione di momento, nodo 1 */
      WM.Add(10, 4, (Mat3x3(FPrime)+Mat3x3(F, Omega1))*Mat3x3(d1Tmp)+Mat3x3(F, R1*(PHIt*b)));
      WM.Add(10, 10, FWedged1Wedge);
      WM.Add(10, 12+2*NModes+12*iStrNodem1+1, O1Wedged1Wedge+Mat3x3(R1*(PHIt*b)));
      WM.Add(10, 12+2*NModes+12*iStrNodem1+7, Mat3x3(d1Tmp));
      
      /* derivata dell'equazione di momento, contributo modale */
      SubMat1.LeftMult((-FWedgeO1Wedge-Mat3x3(FPrime))*R1, PHIt);
      WM.Add(10, 13, SubMat1);
      SubMat1.LeftMult(-FWedge*R1, PHIt);   
      WM.Add(10, 12+NModes+1, SubMat1);
      
      /* derivata dell'eq. di equilibrio ai modi */      
      SubMat2.RightMult(PHItT, R1.Transpose()*FWedge);
      WM.Add(12+NModes+1, 10, SubMat2);
      SubMat2.RightMult(PHItT, R1.Transpose()*FWedgeO1Wedge);
      WM.Add(12+NModes+1, 4, SubMat2);
      SubMat2.RightMult(PHItT, R1.Transpose());
      WM.Add(12+NModes+1, 12+2*NModes+7, SubMat2); 
      
      /* Derivata dell'equazione di momento, nodo 2 */
      WM.Add(12+iGetInitialNumDof()+12*iStrNodem1+10, 12+iGetInitialNumDof()+12*iStrNodem1+4, 
	     (Mat3x3(FPrime)+Mat3x3(F, Omega2))*Mat3x3(-d2Tmp));
      WM.Add(12+iGetInitialNumDof()+12*iStrNodem1+10, 12+iGetInitialNumDof()+12*iStrNodem1+10, FWedged2Wedge);
      WM.Add(12+iGetInitialNumDof()+12*iStrNodem1+10, 12+2*NModes+12*iStrNodem1+1, O2Wedged2Wedge);
      WM.Add(12+iGetInitialNumDof()+12*iStrNodem1+10, 12+2*NModes+12*iStrNodem1+7, Mat3x3(-d2Tmp));
      
      /* Equazione di vincolo */
      WM.Add(12+2*NModes+12*iStrNodem1+1, 4, Mat3x3(d1Tmp));
      WM.Add(12+2*NModes+12*iStrNodem1+1, 12+iGetInitialNumDof()+12*iStrNodem1+4, Mat3x3(-d2Tmp));
      
      /* Equazione di vincolo, contributo modale */
      SubMat1.LeftMult(-R1, PHIt);
      WM.Add(12+2*NModes+12*iStrNodem1+1, 13, SubMat1);
      
      /* Derivata dell'equazione di vincolo */
      WM.Add(12+2*NModes+12*iStrNodem1+7, 4, O1Wedged1Wedge+R1*(PHIt*b));
      WM.Add(12+2*NModes+12*iStrNodem1+7, 10, Mat3x3(d1Tmp));
      WM.Add(12+2*NModes+12*iStrNodem1+7, 12+iGetInitialNumDof()+12*iStrNodem1+4, O2Wedged2Wedge);
      WM.Add(12+2*NModes+12*iStrNodem1+7, 12+iGetInitialNumDof()+12*iStrNodem1+10, Mat3x3(-d2Tmp));
      
      /* Derivata dell'equazione di vincolo, contributo modale */
      SubMat1.LeftMult(-Mat3x3(Omega1)*R1, PHIt);
      WM.Add(12+2*NModes+12*iStrNodem1+7, 13, SubMat1);
      SubMat1.LeftMult(-R1, PHIt);
      WM.Add(12+2*NModes+12*iStrNodem1+7, 12+NModes+1, SubMat1);
      
      /* giunto prismatico */      
      Vec3 M(XCurr, iFlexIndex+2*NModes+12*iStrNodem1+4);
      Vec3 MPrime(XCurr, iFlexIndex+2*NModes+12*iStrNodem1+10);
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
      WM.Add(4, 4, MWedge);
      WM.Add(4, 12+iGetInitialNumDof()+12*iStrNodem1+4, -MWedgeT);
      
      WM.Add(12+iGetInitialNumDof()+12*iStrNodem1+4, 4, MWedgeT);   
      WM.Add(12+iGetInitialNumDof()+12*iStrNodem1+4, 12+iGetInitialNumDof()+12*iStrNodem1+4, -MWedge);
      
      /* Equilibrio dei momenti, termini aggiuntivi modali */      
      Vec3 e1ra(R1.GetVec(1));
      Vec3 e2ra(R1.GetVec(2));
      Vec3 e3ra(R1.GetVec(3));
      Mat3x3 M1Wedge(Mat3x3(e3b, e2ra*M.dGet(1))
		     +Mat3x3(e1b, e3ra*M.dGet(2))
		     +Mat3x3(e2b, e1ra*M.dGet(3)));
      
      SubMat1.LeftMult(M1Wedge, PHIr);
      WM.Add(4, 13, SubMat1);
      WM.Sub(12+iGetInitialNumDof()+12*iStrNodem1+4, 13, SubMat1);
      
      /* Equilibrio ai modi */
      SubMat2.RightMult(PHIrT, R1tot.Transpose()*MWedge);
      WM.Add(13, 4, SubMat2);
      SubMat2.RightMult(PHIrT, R1tot.Transpose()*-MWedgeT);
      WM.Add(13, 12+iGetInitialNumDof()+12*iStrNodem1+4, SubMat2);
      
      Vec3 MA(Mat3x3(e2tota.Cross(e3b), e3tota.Cross(e1b), e1tota.Cross(e2b))*M);
      Mat3x3 MAWedge(MA);
      SubMat2.RightMult(PHIrT, R1tot.Transpose()*MAWedge);
      WM.Add(13, 4, SubMat2);
      
      Mat3x3 R1TMAWedge(R1.Transpose()*MA);
      SubMat1.LeftMult(R1TMAWedge, PHIr);
      Mat3xN SubMat3(NModes);
      MatNxN SubMat4(NModes);   
      SubMat3.LeftMult(R1tot.Transpose()*M1Wedge, PHIr);
      SubMat1+=SubMat3;
      SubMat4.Mult(PHIrT, SubMat1);
      for (unsigned int iMode = 1; iMode <= NModes; iMode++) {
	 for (unsigned int jMode = 1; jMode <= NModes; jMode++) { 
	    WM.fIncCoef(12+iMode, 12+jMode, SubMat4.dGet(iMode,jMode));
	 }
      } 
      
      /* Derivate dell'equilibrio */
      WM.Add(10, 10, MWedge);
      WM.Add(10, 12+iGetInitialNumDof()+12*iStrNodem1+10, -MWedgeT);
      
      WM.Add(12+iGetInitialNumDof()+12*iStrNodem1+10, 10, MWedgeT);   
      WM.Add(12+iGetInitialNumDof()+12*iStrNodem1+10, 12+iGetInitialNumDof()+12*iStrNodem1+10, -MWedge);
      
      MWedge = 
	( (Mat3x3(e3b, Omega1)+Mat3x3(Omega2.Cross(e3b))*M.dGet(1))
	 +Mat3x3(e3b)*MPrime.dGet(1) )*Mat3x3(e2tota)
	 +( (Mat3x3(e1b, Omega1)+Mat3x3(Omega2.Cross(e1b))*M.dGet(2))
	   +Mat3x3(e1b)*MPrime.dGet(2) )*Mat3x3(e3tota)
	   +( (Mat3x3(e2b, Omega1)+Mat3x3(Omega2.Cross(e2b))*M.dGet(3))
	     +Mat3x3(e2b)*MPrime.dGet(3) )*Mat3x3(e1tota);
      
      WM.Add(10, 4, MWedge);
      WM.Add(12+iGetInitialNumDof()+12*iStrNodem1+10, 4, -MWedge);
      
      MWedge =
	( (Mat3x3(e2tota, Omega2)+Mat3x3(Omega1.Cross(e2tota))*M.dGet(1))
	 +Mat3x3(e2tota)*MPrime.dGet(1) )*Mat3x3(e3b)
	 +( (Mat3x3(e3tota, Omega2)+Mat3x3(Omega1.Cross(e3tota))*M.dGet(2))
	   +Mat3x3(e3tota)*MPrime.dGet(2) )*Mat3x3(e1b)
	   +( (Mat3x3(e1tota, Omega2)+Mat3x3(Omega1.Cross(e1tota))*M.dGet(3))
	     +Mat3x3(e1tota)*MPrime.dGet(3) )*Mat3x3(e2b);
      
      WM.Add(10, 12+iGetInitialNumDof()+12*iStrNodem1+4, -MWedge);
      WM.Add(12+iGetInitialNumDof()+12*iStrNodem1+10, 12+iGetInitialNumDof()+12*iStrNodem1+4, MWedge);
      
      /* Derivate dell'equilibrio, termini aggiuntivi modali */
      SubMat1.LeftMult(M1Wedge, PHIr);
      WM.Add(10, 12+NModes+1, SubMat1);
      WM.Sub(12+iGetInitialNumDof()+12*iStrNodem1+10, 12+NModes+1, SubMat1); 
      
      Vec3 v1(e2tota.Cross(e3b)); 
      Vec3 v2(e3tota.Cross(e1b)); 
      Vec3 v3(e1tota.Cross(e2b));
      
      /* Error handling: il programma si ferma, pero' segnala dov'e' l'errore */
      if (v1.Dot() < DBL_EPSILON || v2.Dot() < DBL_EPSILON || v3.Dot() < DBL_EPSILON) {
	 std::cerr << "joint " << GetLabel() << ':' << std::endl
	   << "warning, first node hinge axis and second node hinge axis are (nearly) orthogonal;" << std::endl
	   << "aborting ..." << std::endl;
	 THROW(Joint::ErrGeneric());
      }      
      MWedge = Mat3x3(v1, v2, v3);
      
      /* Equilibrio */
      WM.Add(4, 12+2*NModes+12*iStrNodem1+4, MWedge);
      WM.Add(12+iGetInitialNumDof()+12*iStrNodem1+4, 12+2*NModes+12*iStrNodem1+4, -MWedge);
      
      SubMat2.RightMult(PHIrT, R1tot.Transpose()*MWedge);
      WM.Add(13, 12+2*NModes+12*iStrNodem1+4, SubMat2);   
      
      /* Derivate dell'equilibrio */
      WM.Add(10, 12+2*NModes+12*iStrNodem1+10, MWedge);
      WM.Add(12+iGetInitialNumDof()+12*iStrNodem1+10, 12+2*NModes+12*iStrNodem1+10, -MWedge);
      
      MWedge = MWedge.Transpose();
      
      /* Equaz. di vincolo */
      WM.Add(12+2*NModes+12*iStrNodem1+4, 4, MWedge);
      WM.Add(12+2*NModes+12*iStrNodem1+4, 12+iGetInitialNumDof()+12*iStrNodem1+4, -MWedge);
      
      /* Equaz. di vincolo: termine aggiuntivo dovuto alla flessibilita' */
      Vec3 u1(e2ra.Cross(e3b));
      Vec3 u2(e3ra.Cross(e1b));
      Vec3 u3(e1ra.Cross(e2b));
      
      M1Wedge = (Mat3x3(u1, u2, u3)).Transpose();
      SubMat1.LeftMult(M1Wedge, PHIr);
      WM.Add(12+2*NModes+12*iStrNodem1+4, 13, SubMat1); 
      
      /* Derivate delle equaz. di vincolo */
      WM.Add(12+2*NModes+12*iStrNodem1+10, 10, MWedge);
      WM.Add(12+2*NModes+12*iStrNodem1+10, 12+iGetInitialNumDof()+12*iStrNodem1+10, -MWedge);
      
      v1 = e3b.Cross(e2tota.Cross(Omega1))+e2tota.Cross(Omega2.Cross(e3b));
      v2 = e1b.Cross(e3tota.Cross(Omega1))+e3tota.Cross(Omega2.Cross(e1b));
      v3 = e2b.Cross(e1tota.Cross(Omega1))+e1tota.Cross(Omega2.Cross(e2b));
      
      MWedge = Mat3x3(v1, v2, v3);
      
      /* Derivate dell'equilibrio */
      WM.Add(10, 12+2*NModes+12*iStrNodem1+4, MWedge);
      WM.Add(12+iGetInitialNumDof()+12*iStrNodem1+10, 12+2*NModes+12*iStrNodem1+4, -MWedge);
      
      /* Derivate delle equaz. di vincolo */
      Omega1 = Omega1-Omega2;
      
      v1 = e2tota.Cross(e3b.Cross(Omega1));
      Vec3 v1p(e3b.Cross(Omega1.Cross(e2tota)));
      v2 = e3tota.Cross(e1b.Cross(Omega1));
      Vec3 v2p(e1b.Cross(Omega1.Cross(e3tota)));
      v3 = e1tota.Cross(e2b.Cross(Omega1));
      Vec3 v3p(e2b.Cross(Omega1.Cross(e1tota)));
      
      for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
	 doublereal d = v1.dGet(iCnt);
	 WM.fIncCoef(12+2*NModes+12*iStrNodem1+10, 3+iCnt, d);
	 d = v1p.dGet(iCnt);
	 WM.fIncCoef(12+2*NModes+12*iStrNodem1+10, 12+iGetInitialNumDof()+12*iStrNodem1+3+iCnt, d);
	 
	 d = v2.dGet(iCnt);
	 WM.fIncCoef(12+2*NModes+12*iStrNodem1+11, 3+iCnt, d);
	 d = v2p.dGet(iCnt);
	 WM.fIncCoef(12+2*NModes+12*iStrNodem1+11, 12+iGetInitialNumDof()+12*iStrNodem1+3+iCnt, d);
	 
	 d = v3.dGet(iCnt);
	 WM.fIncCoef(12+2*NModes+12*iStrNodem1+12, 3+iCnt, d);
	 d = v3p.dGet(iCnt);
	 WM.fIncCoef(12+2*NModes+12*iStrNodem1+12, 12+iGetInitialNumDof()+12*iStrNodem1+3+iCnt, d);
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
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
   
   VecN Ka(NModes);
   
   integer iRigidIndex = pModalNode->iGetFirstIndex();
   
   for (integer iCnt = 1; iCnt <= 12; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iRigidIndex+iCnt);
   }
   
   integer iFlexIndex = iGetFirstIndex();
   
   for (unsigned int iCnt = 1; iCnt <= iGetInitialNumDof(); iCnt++) {
      WorkVec.fPutRowIndex(12+iCnt, iFlexIndex+iCnt);
   }
   
   for (unsigned int iStrNodem1 = 0; iStrNodem1 < NStrNodes; iStrNodem1++) {
      integer iNodeFirstPosIndex = pInterfaceNodes[iStrNodem1]->iGetFirstPositionIndex();
      integer iNodeFirstVelIndex = iNodeFirstPosIndex + 6;
      for (unsigned int iCnt = 1; iCnt <= 6; iCnt++) { 
	 WorkVec.fPutRowIndex(12+iGetInitialNumDof()+12*iStrNodem1+iCnt, 
			 iNodeFirstPosIndex+iCnt);
	 WorkVec.fPutRowIndex(12+iGetInitialNumDof()+12*iStrNodem1+6+iCnt, 
			 iNodeFirstVelIndex+iCnt);
      }
   }
   
   /* aggiorna le forze modali : K*a, C*aP */
   for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
      a.Put(iCnt, XCurr.dGetCoef(iFlexIndex+iCnt));
      b.Put(iCnt, XCurr.dGetCoef(iFlexIndex+NModes+iCnt));
   }
   
   for (unsigned int iCnt = 1; iCnt <= NModes; iCnt++) {
      double temp1 = 0, temp2 = 0; 
      for (unsigned int jCnt = 1; jCnt <= NModes; jCnt++) { 
	 temp1 += pModalStiff->dGet(iCnt,jCnt)*a.dGet(jCnt);
	 temp2 += pModalDamp->dGet(iCnt,jCnt)*b.dGet(jCnt); 
      }
      WorkVec.fIncCoef(12+iCnt, -temp1-temp2); 
   }
   
   /* equazioni di vincolo */
   
   /* Recupera i dati */
   Vec3 x1(pModalNode->GetXCurr());
   Vec3 v1(pModalNode->GetVCurr());
   Mat3x3 R1(pModalNode->GetRCurr());
   Vec3 Omega1(pModalNode->GetWCurr());
   
   for (unsigned int iStrNode = 1; iStrNode <= NStrNodes; iStrNode++) {
      unsigned int iStrNodem1 = iStrNode - 1;
      Mat3xN PHIt(NModes,0), PHIr(NModes,0);
      
      for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
	 for (unsigned int iMode = 1; iMode <= NModes; iMode++) {
	    PHIt.Put(iCnt, iMode, pPHIt->dGet(iCnt, (iMode-1)*NStrNodes+iStrNode));
	    PHIr.Put(iCnt, iMode, pPHIr->dGet(iCnt, (iMode-1)*NStrNodes+iStrNode));
	 }
      }
      
      MatNx3 PHItT(NModes), PHIrT(NModes);
      PHItT.Transpose(PHIt);
      PHIrT.Transpose(PHIr);
      
      Vec3 d1rig, d2;  
      d1rig.Put(1, pOffsetNodes->dGet(1, iStrNode)); 
      d1rig.Put(2, pOffsetNodes->dGet(2, iStrNode)); 
      d1rig.Put(3, pOffsetNodes->dGet(3, iStrNode));
      d2.Put(1, pOffsetNodes->dGet(1, NStrNodes+iStrNode)); 
      d2.Put(2, pOffsetNodes->dGet(2, NStrNodes+iStrNode)); 
      d2.Put(3, pOffsetNodes->dGet(3, NStrNodes+iStrNode)); 
      
      Vec3 d1tot = d1rig+PHIt*a;
      Mat3x3 R1tot = R1*Mat3x3(1., PHIr*a);
      
      Vec3 x2(pInterfaceNodes[iStrNodem1]->GetXCurr());
      Vec3 v2(pInterfaceNodes[iStrNodem1]->GetVCurr());
      Mat3x3 R2(pInterfaceNodes[iStrNodem1]->GetRCurr());
      Vec3 Omega2(pInterfaceNodes[iStrNodem1]->GetWCurr());
      Vec3 F(XCurr, iFlexIndex+2*NModes+12*iStrNodem1+1);
      Vec3 FPrime(XCurr, iFlexIndex+2*NModes+12*iStrNodem1+7);
      
      /* cerniera sferica */
      
      /* Distanza nel sistema globale */
      Vec3 d1Tmp(R1*d1tot);
      Vec3 d2Tmp(R2*d2);
      
      /* Vettori omega1/\d1, -omega2/\d2 */
      Vec3 O1Wedged1(Omega1.Cross(d1Tmp));
      Vec3 O2Wedged2(Omega2.Cross(d2Tmp));
      
      /* d1Prime= w1/\d1 + R*PHIt*b */
      Mat3xN R1PHIt(NModes);
      R1PHIt.LeftMult(R1, PHIt);
      Vec3 d1Prime(O1Wedged1+R1PHIt*b);
      
      /* Equazioni di equilibrio, nodo 1 */
      WorkVec.Add(1, -F);
      WorkVec.Add(4, F.Cross(d1Tmp)); /* Sfrutto il fatto che F/\d = -d/\F */
      
      /* Eq. d'equilibrio ai modi */
      for (unsigned int iMode = 1; iMode <= NModes; iMode++) {
	 doublereal temp = 0.;
	 for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
	    temp += PHIt.dGet(iCnt, iMode)*(R1.Transpose()*F).dGet(iCnt); 
	 }
	 WorkVec.fIncCoef(12+iMode, -temp);
      }
      
      /* Derivate delle equazioni di equilibrio, nodo 1 */
      WorkVec.Add(7, -FPrime);
      WorkVec.Add(10, FPrime.Cross(d1Tmp)-d1Prime.Cross(F));
      
      /* Derivata dell'eq. di equilibrio ai modi */
      MatNx3 PHItTR1T(NModes);
      PHItTR1T.RightMult(PHItT, R1.Transpose());
      for (unsigned int iMode = 1; iMode <= NModes; iMode++) {
	 doublereal temp = 0.;
	 for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
	    temp += PHItTR1T.dGet(iMode,iCnt)*(Omega1.Cross(F)-FPrime).dGet(iCnt); 
	 }
	 WorkVec.fIncCoef(12+NModes+iMode, temp);
      }
      
      /* Equazioni di equilibrio, nodo 2 */
      WorkVec.Add(12+iGetInitialNumDof()+12*iStrNodem1+1, F);
      WorkVec.Add(12+iGetInitialNumDof()+12*iStrNodem1+4, d2Tmp.Cross(F)); 
      
      /* Derivate delle equazioni di equilibrio, nodo 2 */
      WorkVec.Add(12+iGetInitialNumDof()+12*iStrNodem1+7, FPrime);
      WorkVec.Add(12+iGetInitialNumDof()+12*iStrNodem1+10, d2Tmp.Cross(FPrime)+O2Wedged2.Cross(F));
      
      /* Equazione di vincolo */
      WorkVec.Add(12+2*NModes+12*iStrNodem1+1, x1+d1Tmp-x2-d2Tmp);
      
      /* Derivata dell'equazione di vincolo */
      WorkVec.Add(12+2*NModes+12*iStrNodem1+7, v1+d1Prime-v2-O2Wedged2);
      
      /* giunto prismatico */
      Vec3 M(XCurr, iFlexIndex+2*NModes+12*iStrNodem1+4);
      Vec3 MPrime(XCurr, iFlexIndex+2*NModes+12*iStrNodem1+10); 
      
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
      WorkVec.Add(4, -MTmp);
      
      /* Equazioni di equilibrio, nodo 2 */
      WorkVec.Add(12+iGetInitialNumDof()+12*iStrNodem1+4, MTmp);
      
      /* Equazioni di equilibrio, contributo modale */
      for (unsigned int iMode = 1; iMode <= NModes; iMode++) {
	 doublereal temp = 0;
	 for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
	    temp += PHIrT.dGet(iMode, iCnt)*(R1tot.Transpose()*MTmp).dGet(iCnt);
	 } 
	 WorkVec.fIncCoef(12+iMode, -temp);
      }
      
      /* eaPrime = w/\ea + R1*[(PHIr*b)/\]ia */
      Vec3 i1(1.,0.,0.);
      Vec3 i2(0.,1.,0.);
      Vec3 i3(0.,0.,1.);

      Mat3x3 Tmp = R1*Mat3x3(PHIr*b);
      Vec3 e1aPrime = Omega1.Cross(e2a)+Tmp*i1; 
      Vec3 e2aPrime = Omega1.Cross(e2a)+Tmp*i2; 
      Vec3 e3aPrime = Omega1.Cross(e2a)+Tmp*i3;
      Vec3 MTmpPrime(0.); 
      MTmpPrime = 
	(e2a.Cross(Omega2.Cross(e3b))-e3b.Cross(e2aPrime))*M.dGet(1)
	  +(e3a.Cross(Omega2.Cross(e1b))-e1b.Cross(e3aPrime))*M.dGet(1)
	    +(e1a.Cross(Omega2.Cross(e2b))-e2b.Cross(e1aPrime))*M.dGet(1)
	      +e2a.Cross(e3b*MPrime.dGet(1))
		+e3a.Cross(e1b*MPrime.dGet(2))
		  +e1a.Cross(e2b*MPrime.dGet(3));   
      
      /* Derivate delle equazioni di equilibrio, nodo 1 */
      WorkVec.Add(10, -MTmpPrime);
      
      /* Derivate delle equazioni di equilibrio, nodo 2 */
      WorkVec.Add(12+iGetInitialNumDof()+12*iStrNodem1+10, MTmpPrime);
      
      /* Derivate delle equazioni di equilibrio, contributo modale */
      MatNx3 SubMat1(NModes);
      SubMat1.RightMult(PHIrT, R1tot.Transpose());

      // FIXME: temporary ((PHIr*b).Cross(R1.Transpose()*MTmp))
      Vec3 T1 = MTmpPrime-Omega1.Cross(MTmp);
      Vec3 T2 = (PHIr*b).Cross(R1.Transpose()*MTmp);
      
      for (unsigned int iMode = 1; iMode <= NModes; iMode++) {
	 doublereal temp = 0;
	 for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
	    temp += SubMat1.dGet(iMode, iCnt)*T1.dGet(iCnt)
		    - PHIrT.dGet(iMode, iCnt)*T2.dGet(iCnt);
	 } 
	 WorkVec.fDecCoef(12+NModes+iMode, temp);
      }
      
      /* Equazioni di vincolo di rotazione */
      WorkVec.fIncCoef(12+2*NModes+12*iStrNodem1+4, -(e3b.Dot(e2a)));
      WorkVec.fIncCoef(12+2*NModes+12*iStrNodem1+5, -(e1b.Dot(e3a)));
      WorkVec.fIncCoef(12+2*NModes+12*iStrNodem1+6, -(e2b.Dot(e1a)));
      
      /* Derivate delle equazioni di vincolo di rotazione */
      WorkVec.fIncCoef(12+2*NModes+12*iStrNodem1+10, e3b.Dot(Omega2.Cross(e2a)-e2aPrime));
      WorkVec.fIncCoef(12+2*NModes+12*iStrNodem1+11, e1b.Dot(Omega2.Cross(e3a)-e3aPrime));
      WorkVec.fIncCoef(12+2*NModes+12*iStrNodem1+12, e2b.Dot(Omega2.Cross(e1a)-e1aPrime));
      
   } /* fine equazioni di vincolo */
   
   return WorkVec;
}


void 
Modal::SetValue(VectorHandler& X, VectorHandler& XP) const
{
   /* inizializza la soluzione e la sua derivata
    * subito dopo l'assemblaggio iniziale 
    * e prima del calcolo delle derivate */
   
   /* integer iRigidIndex = pModalNode->iGetFirstIndex(); */
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
   return 2*NModes;
}


doublereal 
Modal::dGetPrivData(unsigned int i) const
{
   ASSERT(i > 0 && i < iGetNumPrivData());
   if (i <= NModes) {
      return a.dGet(i);
   } else {
      return b.dGet(i-NModes);
   }
}


Joint* 
ReadModal(DataManager* pDM, 
	  MBDynParser& HP, 
	  const DofOwner* pDO, 
	  unsigned int uLabel)
{
   /* legge i dati d'ingresso e li passa al costruttore dell'elemento */
   Joint* pEl = NULL;
   unsigned int uNode = HP.GetInt(); 
   
   DEBUGCOUT("Linked to Modal Node: " << uNode << std::endl);    
   
   /* verifica di esistenza del nodo */  
   StructNode* pTmpNode; 
   ModalNode* pModalNode;
   
   if ((pTmpNode = pDM->pFindStructNode(uNode)) == NULL) {
      std::cerr << "structural node " << uNode
	<< " at line " << HP.GetLineData()
	<< " not defined" << std::endl;
      THROW(DataManager::ErrGeneric());
   }
   
   if (pTmpNode->GetStructNodeType() != StructNode::MODAL) {
      std::cerr << "Illegal structural node type for body " << uLabel 
	      << " at line " << HP.GetLineData() << std::endl;
      THROW(DataManager::ErrGeneric());
   }
   pModalNode = (ModalNode *)pTmpNode;
   
   /* la posizione del nodo modale e' quella dell'origine del SdR 
    * del corpo flessibile */
   Vec3 X0(pModalNode->GetXCurr());
   
   /* orientazione del corpo flessibile, data come orientazione 
    * del nodo modale collegato */
   Mat3x3 R(pModalNode->GetRCurr());
   
   /* Legge i dati relativi al corpo flessibile */      
   unsigned int NModes = HP.GetInt();     /* numero di modi */
   if (NModes == 0) {
      std::cerr << "Illegal number of modes at line " << HP.GetLineData()
	      << std::endl;
      THROW(DataManager::ErrGeneric());
   }

   unsigned int NFemNodes = HP.GetInt();  /* numero di nodi FEM del modello */
   if (NModes == 0) {
      std::cerr << "Illegal number of FEM nodes at line " << HP.GetLineData()
	      << std::endl;
      THROW(DataManager::ErrGeneric());
   }
   
#ifdef MODAL_SCALE_DATA
   /* Legge gli eventuali fattori di scala per le masse nodali (scale masses)
    * e per i modi (scale modes)   */
   
   /* NOTA: E' COMMENTATO PERCHE' AL MOMENTO NON SERVE */
   doublereal scalemasses = 1.;
   doublereal scalemodes = 1.;
   
   if (HP.IsKeyWord("scalemasses")) {
      scalemasses = HP.GetReal();
   }
   
   if (HP.IsKeyWord("scalemodes")) {
      scalemodes = HP.GetReal();
   }
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
      for (unsigned int iCnt = 1; iCnt <= NModes; iCnt ++) {
         iDampFlag = 1;
         integer iDampedMode =  HP.GetInt();
         cdamp = HP.GetReal();
         DampRatios.Put(iDampedMode, cdamp);    
      }
   } else {
      silent_cout("no damping is assumed at line " 
		  << HP.GetLineData() << " (deprecated)" << std::endl);
   }    
   
   DEBUGCOUT("Number of Modes Imported : " << NModes << std::endl);    
   DEBUGCOUT("Number of FEM Nodes Imported : " << NFemNodes << std::endl);
   DEBUGCOUT("Origin of FEM Model : " << X0 << std::endl);
   DEBUGCOUT("Orientation of FEM Model : " << R << std::endl);
   /* DEBUGCOUT("Damping coefficient: "<< cdamp << std::endl); */
   
   doublereal dMass = 0;              /* massa totale */
   Vec3 STmp(0.);                     /* momenti statici  */
   Mat3x3 JTmp(0.);                   /* inerzia totale  */
   VecN   FemMass(NFemNodes, 0.);     /* masse nodali   */
   Mat3xN FemJ(NFemNodes, 0.);        /* inerzie nodali (sono diagonali) */
   
   Mat3xN *pModeShapest = NULL;       /* forme modali di traslazione e rotaz. */
   Mat3xN *pModeShapesr = NULL;
   Mat3xN PHIti(NModes, 0.);          /* forme modali nodo i-esimo: 3*nmodi */
   Mat3xN PHIri(NModes, 0.);
   Vec3   PHItij(0.);                 /* j-esima forma del nodo i-esimo */
   Vec3   PHIrij(0.);
   Mat3xN *pXYZFemNodes = NULL;       /* puntatore alle coordinate nodali */
   Mat3xN *pXYZOffsetNodes = NULL;    /* punt. offset nodali (per vincoli) */
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
   const char *sFileFem = HP.GetFileName();
   
   /* apre il file con i dati del modello FEM */ 
   std::ifstream fdat(sFileFem);
   if (!fdat) {
      std::cerr << std::endl << "Unable to open file " << sFileFem 
	      << " at line " << HP.GetLineData() << std::endl;
      THROW(DataManager::ErrGeneric());
   }
   DEBUGCOUT("Reading Flexible Body Data from file " << sFileFem << std::endl);
   
   /* carica i dati relativi a coordinate nodali, massa, momenti statici 
    * e d'inerzia, massa e rigidezza generalizzate dal file nomefile.flex.
    * Questo file e' quello generato da DADS */
   
   doublereal d;
   unsigned int i, NFemNodesDADS = 0, NModesDADS = 0, NRejModes = 0;
   char str[BUFSIZ];
   
   /* alloca la memoria per le matrici necessarie a memorizzare i dati 
    * relativi al corpo flessibile
    * nota: devo usare i puntatori perche' altrimenti non si riesce 
    * a passarli al costruttore       */
   SAFENEWWITHCONSTRUCTOR(pXYZFemNodes, Mat3xN, Mat3xN(NFemNodes, 0.));
   SAFENEWWITHCONSTRUCTOR(pGenMass,  MatNxN, MatNxN(NModes, 0.));
   SAFENEWWITHCONSTRUCTOR(pGenStiff, MatNxN, MatNxN(NModes, 0.));
   SAFENEWWITHCONSTRUCTOR(pGenDamp,  MatNxN, MatNxN(NModes, 0.));
   SAFENEWWITHCONSTRUCTOR(pModeShapest, Mat3xN, Mat3xN(NFemNodes*NModes, 0.));
   SAFENEWWITHCONSTRUCTOR(pModeShapesr, Mat3xN, Mat3xN(NFemNodes*NModes, 0.));
   SAFENEWWITHCONSTRUCTOR(pInv3, Mat3xN, Mat3xN(NModes, 0.));
   SAFENEWWITHCONSTRUCTOR(pInv4, Mat3xN, Mat3xN(NModes, 0.));
   SAFENEWWITHCONSTRUCTOR(pInv5, Mat3xN, Mat3xN(NModes*NModes, 0.));   /* Inv5 e' un 3xMxM */
   SAFENEWWITHCONSTRUCTOR(pInv8, Mat3xN, Mat3xN(3*NModes, 0.));        /* Inv8 e' un 3x3xM */
   SAFENEWWITHCONSTRUCTOR(pInv9, Mat3xN, Mat3xN(3*NModes*NModes, 0.)); /* Inv9 e' un 3x3xMxM */
   SAFENEWWITHCONSTRUCTOR(pInv10,Mat3xN, Mat3xN(3*NModes, 0.));        /* Inv10 e' un 3x3xM */
   SAFENEWWITHCONSTRUCTOR(pInv11,Mat3xN, Mat3xN(NModes, 0.)); 
   SAFENEWWITHCONSTRUCTOR(a,  VecN, VecN(NModes, 0.));
   SAFENEWWITHCONSTRUCTOR(aP, VecN, VecN(NModes, 0.));

   unsigned int *IdFemNodes;        /* array contenente le label dei nodi FEM */
   IdFemNodes = new unsigned int[NFemNodes];
   if (!IdFemNodes) {    
      THROW(DataManager::ErrGeneric());
   }
  
#if 0 
   doublereal scalfact = scalemodes*scalemodes;
#endif

   while (!fdat.eof()) {        /* parsing del file */ 
      fdat.getline(str, sizeof(str));

      /* legge il primo blocco (HEADER) */
      if (!strncmp("** RECORD GROUP 1,", str, sizeof("** RECORD GROUP 1,") - 1)) {
	 fdat.getline(str, sizeof(str));
	 fdat >> str;
	 for (unsigned int iCnt = 0; iCnt < 5; iCnt++)  { 
	    fdat >> i;  
	    if (iCnt == 0) {
	       NFemNodesDADS = i;
	    } else if (iCnt >= 1 && iCnt <=3) {
	       NModesDADS += i;
	    } else if (iCnt == 4) {
	       NRejModes = i;
	       NModesDADS -= NRejModes;
	    }
	 }
	 if (NFemNodes != NFemNodesDADS) { 
	    std::cerr << "file '" << sFileFem << "': FEM nodes " << NFemNodes
		    << " do not match node number " << NFemNodesDADS
		    << "of modal joint " << uLabel << std::endl;
	    THROW(DataManager::ErrGeneric());
	 }
	 if (NModes != NModesDADS) { 
	    std::cerr << "file '" << sFileFem << "': using " << NModes
		    << " of " << NModesDADS
		    << " modes for modal joint " << uLabel << std::endl;
	 }

      /* legge il secondo blocco (Id.nodi) */
      } else if (!strncmp("** RECORD GROUP 2,", str, sizeof("** RECORD GROUP 2,") - 1)) {
	 unsigned int ui;
	 for (iNode = 1; iNode <= NFemNodes; iNode++) {
	    fdat >> ui;
	    IdFemNodes[iNode-1] = ui;
	 } 
      
      /* deformate iniziali dei modi */
      } else if (!strncmp("** RECORD GROUP 3,", str, sizeof("** RECORD GROUP 3,") - 1)) {
	 for (iMode = 1; iMode <= NModesDADS; iMode++) {
	    fdat >> d;
	    if (iMode > NModes) {
	       continue;
	    }
	    a->Put(iMode, d);
	 } 
	
      /* velocita' iniziali dei modi */   
      } else if (!strncmp("** RECORD GROUP 4,", str, sizeof("** RECORD GROUP 4,") - 1)) {
	 for (iMode = 1; iMode <= NModesDADS; iMode++) {
	    fdat >> d;
	    if (iMode > NModes) {
	       continue;
	    }
	    aP->Put(iMode, d);
	 } 

      /* Coordinate X dei nodi*/
      } else if (!strncmp("** RECORD GROUP 5,", str, sizeof("** RECORD GROUP 5,") - 1)) {
	 for (iNode = 1; iNode <= NFemNodes; iNode++) {
	    fdat >> d;  
	    pXYZFemNodes->Put(1, iNode, d /**scalemodes*/ );
	 }

      /* Coordinate Y dei nodi*/    
      } else if (!strncmp("** RECORD GROUP 6,", str, sizeof("** RECORD GROUP 6,") - 1)) {
	 for (iNode = 1; iNode <= NFemNodes; iNode++) {
	    fdat >> d;  
	    pXYZFemNodes->Put(2, iNode, d /**scalemodes*/ );
	 }

      /* Coordinate Z dei nodi*/   
      } else if (!strncmp("** RECORD GROUP 7,", str, sizeof("** RECORD GROUP 7,") - 1)) {
	 for (iNode = 1; iNode <= NFemNodes; iNode++) {
	    fdat >> d;  
	    pXYZFemNodes->Put(3, iNode, d /**scalemodes*/ );
	 }

      /* Forme modali */
      } else if (!strncmp("** RECORD GROUP 8,", str, sizeof("** RECORD GROUP 8,") - 1)) {
	 for (iMode = 1; iMode <= NRejModes; iMode++) {
     		fdat.getline(str,sizeof(str));
     		fdat.getline(str,sizeof(str));
	 }
	 for (iMode = 1; iMode <= NModesDADS; iMode++) {
     	    fdat.getline(str,sizeof(str));
	    for (iNode = 1; iNode <= NFemNodes; iNode++) {
	       for (int iCnt = 1; iCnt <= 3; iCnt++) {
		  fdat >> d;
		  if (iMode > NModes) {
		     continue;
		  }
		  pModeShapest->Put(iCnt, (iMode-1)*NFemNodes+iNode, 
				  d /* *scalemodes */ );
	       }
	       for (int iCnt = 1; iCnt <= 3; iCnt++) {
	          fdat >> d;
		  if (iMode > NModes) {
		     continue;
		  }
		  pModeShapesr->Put(iCnt, (iMode-1)*NFemNodes+iNode, d);
	       }
	    }
	    fdat.getline(str,sizeof(str));
	 }   
   
      /* Matrice di massa  modale */
      }	else if (!strncmp("** RECORD GROUP 9,", str, sizeof("** RECORD GROUP 9,") - 1)) {
         for (iMode = 1; iMode <= NModes; iMode++) {
	    for (jMode = 1; jMode <= NModes; jMode++) {
	       fdat >> d;
	       if (iMode > NModes || jMode > NModes) {
		  continue;
	       }
	       pGenMass->Put(iMode, jMode, d /* *scalfact */ );
	    }
         }

      /* Matrice di rigidezza  modale */
      } else if (!strncmp("** RECORD GROUP 10,", str, sizeof("** RECORD GROUP 10,") - 1)) {
         for (iMode = 1; iMode <= NModes; iMode++) {
	    for (jMode = 1; jMode <= NModes; jMode++) {
	       fdat >> d;
	       if (iMode > NModes || jMode > NModes) {
		  continue;
	       }
	       pGenStiff->Put(iMode, jMode, d/* *scalfact */);
	    }
         }

      /* Lumped Masses */
      } else if (!strncmp("** RECORD GROUP 11,", str, sizeof("** RECORD GROUP 11,") - 1)) {
         for (iNode = 1; iNode <= NFemNodes; iNode++) {
	    for (unsigned int jCnt = 1; jCnt <= 6; jCnt++) {
	       fdat >> d;
	       switch (jCnt) {
	       case 1:
	          FemMass.Put(iNode, d/**scalemasses*/);
	          break;
	       case 4:
	          FemJ.Put(1, iNode, d/**scalemasses*/);
	          break;
	       case 5:
	          FemJ.Put(2, iNode, d/**scalemasses*/);
		  
		  
	          break;
	       case 6:
	          FemJ.Put(3, iNode, d/**scalemasses*/);
	          break;
	       }
	    }
         }
      } /* fine parser del file */
   }
   fdat.close();   
   
   /* lettura dati di vincolo:
    * l'utente specifica la label del nodo FEM e del nodo rigido 
    * d'interfaccia.
    * L'orientamento del nodo FEM e' quello del nodo modale, la
    * posizione e' la somma di quella modale e di quella FEM   */    

   const StructNode** pInterfaceNodes; /* puntatori ai nodi multibody */
   unsigned int *IntNodes;       /* array contenente le label dei nodi d'interfaccia */
   Mat3xN* pPHItStrNode = NULL;  /* array contenente le forme modali dei nodi d'interfaccia */
   Mat3xN* pPHIrStrNode = NULL;
   
   unsigned int NStrNodes = HP.GetInt();  /* numero di nodi d'interfaccia */
   DEBUGCOUT("Number of Interface Nodes : " << NStrNodes << std::endl);
   
   SAFENEWWITHCONSTRUCTOR(pXYZOffsetNodes, Mat3xN, Mat3xN(2*NStrNodes+1, 0.));
   SAFENEWWITHCONSTRUCTOR(pPHItStrNode, Mat3xN, Mat3xN(NStrNodes*NModes, 0.));
   SAFENEWWITHCONSTRUCTOR(pPHIrStrNode, Mat3xN, Mat3xN(NStrNodes*NModes, 0.));
   
   SAFENEWARR(pInterfaceNodes, const StructNode*, NStrNodes);
   if (pInterfaceNodes == NULL) {    
      THROW(DataManager::ErrGeneric());
   }
   IntNodes = new unsigned int[2*NStrNodes];
   if (!IntNodes) {    
      THROW(DataManager::ErrGeneric());
   }
   
   for (iStrNode = 1; iStrNode <= NStrNodes; iStrNode++) {
      
      /* nodo collegato 1 (è il nodo FEM) */
      unsigned int uNode1 = (unsigned int)HP.GetInt();	    
      DEBUGCOUT("Linked to FEM Node " << uNode1 << std::endl);
      
      /* offset del nodo FEM */
      ReferenceFrame RF(pModalNode);
      Vec3 d1(HP.GetPosRel(RF));
      
      DEBUGCOUT("Fem Node reference frame d1:" << std::endl << d1 << std::endl);
      
      /* verifica di esistenza del nodo 1*/  
      for (iNode = 1; iNode <= NFemNodes; iNode++) {
	 if (uNode1 == IdFemNodes[iNode-1]) {
	    break;
	 }
	 if (iNode == NFemNodes) { 
	    std::cerr << "FEM node " << uNode1
	      << " at line " << HP.GetLineData() 
		<< " not defined " << std::endl;
	    THROW(DataManager::ErrGeneric());
	 }
      }
      
      int iNodeCurr = iNode;
      
      /* recupera la posizione del nodo FEM, somma di posizione 
       * e eventuale offset;
       * nota: iNodeCurr contiene la posizione a cui si trova
       * il nodo FEM d'interfaccia nell'array pXYZNodes */
      for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
	 pXYZOffsetNodes->Put(iCnt, iStrNode, pXYZFemNodes->dGet(iCnt, iNodeCurr)+d1.dGet(iCnt));
      }
      /* salva le forme modali del nodo d'interfaccia nell'array pPHIStrNode */
      
      for (iMode = 1; iMode <= NModes; iMode++) {
         for (unsigned int iCnt = 1; iCnt <= 6; iCnt++) {
	    if (iCnt <= 3) {
	      pPHItStrNode->Put(iCnt, (iMode-1)*NStrNodes+iStrNode, 
				pModeShapest->dGet(iCnt, (iMode-1)*NFemNodes+iNodeCurr));
	    } else {
	      pPHIrStrNode->Put(iCnt-3, (iMode-1)*NStrNodes+iStrNode, 
				pModeShapesr->dGet(iCnt-3, (iMode-1)*NFemNodes+iNodeCurr));
	    }
	 }
      }      
      
      /* nodo collegato 2 (e' il nodo multibody) */
      unsigned int uNode2 = (unsigned int)HP.GetInt();	     
      DEBUGCOUT("Linked to Multi-Body Node " << uNode2 << std::endl);
      
      /* verifica di esistenza del nodo 2 */
      pInterfaceNodes[iStrNode-1] = pDM->pFindStructNode(uNode2);
      if (pInterfaceNodes[iStrNode-1] == NULL) {
         std::cerr << "structural node " << uNode2
           << " at line " << HP.GetLineData() 
	     << " not defined" << std::endl;	  
         THROW(DataManager::ErrGeneric());
      }
      
      /* offset del nodo Multi-Body */
      RF = ReferenceFrame(pInterfaceNodes[iStrNode-1]);
      Vec3 d2(HP.GetPosRel(RF));
      
      for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
	 pXYZOffsetNodes->Put(iCnt, NStrNodes+iStrNode, d2.dGet(iCnt));
      }
      
      DEBUGCOUT("Multibody Node reference frame d2:" << std::endl << d2 << std::endl);
      
      /* salva le label dei nodi vincolati nell'array IntNodes
       * non serve piu' ricordarsi di toglierlo */
      IntNodes[iStrNode-1] = uNode1;
      IntNodes[NStrNodes+iStrNode-1] = uNode2;
      
   }  /* fine ciclo sui nodi d'interfaccia */
   
   
   /* fine ciclo caricamento dati */ 
   
   doublereal mi;
   Vec3 ui; 
   
   /* 
    * calcola gli invarianti d'inerzia (massa, momenti statici e d'inerzia, 
    * termini di accoppiamento nel SdR locale) 
    */
   
   /* inizio ciclo scansione nodi */
   for (iNode = 1; iNode <= NFemNodes; iNode++) {
      mi = FemMass.dGet(iNode);
      /* massa totale (Inv 1) */
      dMass += mi;                                   
     
      for (unsigned int iCnt = 1; iCnt <= 3;iCnt++) {
      	 /* vettore posizione indeformata del nodo {ui} */
	 ui.Put(iCnt, pXYZFemNodes->dGet(iCnt, iNode));
      }
      
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
	 for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
	    PHIti.Put(iCnt, iMode, pModeShapest->dGet(iCnt, iOffset));          
	    PHIri.Put(iCnt, iMode, pModeShapesr->dGet(iCnt, iOffset));
	 }
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
	 for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
	    /* estrae la j-esima funzione di forma del nodo i-esimo */
	    PHItij.Put(iCnt, PHIti.dGet(iCnt, iMode));
	    PHIrij.Put(iCnt, PHIri.dGet(iCnt, iMode));  
	 }
	 
	 Mat3x3 PHItijvett_mi(PHItij*mi);
	 Mat3xN Inv5jTmp(NModes, 0);

    	 /* Inv5 = mi*PHItij/\*PHIti, i = 1,...nnodi, j = 1,...nmodi */
	 Inv5jTmp.LeftMult(PHItijvett_mi, PHIti);
	 for (jMode = 1; jMode <= NModes; jMode++)  {
	    for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
	       pInv5->Add(iCnt, (iMode-1)*NModes+jMode, 
			       Inv5jTmp.dGet(iCnt, jMode));
	    }
	 }

	 /* Inv8 = -mi*ui/\*PHItij/\, i = 1,...nnodi, j = 1,...nmodi */
	 Mat3x3 Inv8jTmp = -uivett*PHItijvett_mi;
	 for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
	    for (unsigned int jCnt = 1; jCnt <= 3; jCnt++) {
	       pInv8->Add(iCnt, (iMode-1)*3+jCnt, Inv8jTmp.dGet(iCnt,jCnt));
	    }
	 }

	 /* Inv9 = mi*PHItij/\*PHItik/\, i = 1,...nnodi, j, k = 1...nmodi */
	 Vec3 PHItik;
	 for (unsigned int kMode = 1; kMode <= NModes; kMode++) { 
	    for (unsigned int kCnt = 1; kCnt <= 3; kCnt++) {
	       PHItik.Put(kCnt, PHIti.dGet(kCnt, kMode));
	    }
	    Mat3x3 PHItikvett(PHItik); 
	    Mat3x3 Inv9jkTmp = PHItijvett_mi*PHItikvett;
	    for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
	       for (unsigned int jCnt = 1; jCnt <= 3; jCnt++) {
		  pInv9->Add(iCnt, (iMode-1)*3*NModes+(kMode-1)*3+jCnt, 
				  Inv9jkTmp.dGet(iCnt,jCnt));
	       }
	    }
	 }    
	 
         /* Inv10 = [PHIrij/\][J0i], i = 1,...nnodi, j = 1,...nmodi */
	 Mat3x3 Inv10jTmp = Mat3x3(PHIrij)*JiNodeTmp;
	 for (unsigned int iCnt = 1; iCnt <= 3; iCnt++) {
	    for (unsigned int jCnt = 1; jCnt <= 3; jCnt++) {
	       pInv10->Add(iCnt, (iMode-1)*3+jCnt, Inv10jTmp.dGet(iCnt,jCnt));
	    }
	 }
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
				IntNodes, 
				pXYZFemNodes, 
				pXYZOffsetNodes,
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
