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
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
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

/* Atti di moto piani */

#include <mbconfig.h>

#include <planedj.h>

/* PlaneDispJoint - begin */

/* Costruttore non banale */
PlaneDispJoint::PlaneDispJoint(unsigned int uL, 
					 const DofOwner* pDO,
					 const StructNode* pN1, 
					 const StructNode* pN2,
					 const Vec3& dTmp1, 
					 const Vec3& dTmp2,
					 const Mat3x3& R1hTmp, 
					 const Mat3x3& R2hTmp,
					 flag fOut)
: Elem(uL, ElemType::JOINT, fOut),
Joint(uL, JointType::PLANEDISP, pDO, fOut),
pNode1(pN1), pNode2(pN2),
d1(dTmp1), R1h(R1hTmp), d2(dTmp2), R2h(R2hTmp), F(0.), M(0.)
{
   NO_OP;
}


/* Distruttore banale */
PlaneDispJoint::~PlaneDispJoint(void)
{
   NO_OP;
};


/* Contributo al file di restart */
ostream& 
PlaneDispJoint::Restart(ostream& out) const
{
   return out << "plane disp joint not implemented yet" << endl;
}


/* Assemblaggio jacobiano */
VariableSubMatrixHandler& 
PlaneDispJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		       doublereal dCoef,
		       const VectorHandler& /* XCurr */ ,
		       const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUTFNAME("PlaneDispJoint::AssJac");
   
   /* Setta la sottomatrice come piena (e' un po' dispersivo, ma lo jacobiano 
    * e' complicato)
    * Nota: i primi 6+6 indici si riferiscono ai nodi;
    * l'indice 13 si riferisce alla reazione normale al piano;
    * gli indici 14 e 15 si riferiscono alle coppie */
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Ridimensiona la sottomatrice in base alle esigenze */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->WorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeInit(iNumRows, iNumCols, 0.);
   
   /* Recupera gli indici delle varie incognite */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();

   /* Setta gli indici delle equazioni */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {	
      WM.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WM.fPutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.fPutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
      WM.fPutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
   }
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WM.fPutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
      WM.fPutColIndex(12+iCnt, iFirstReactionIndex+iCnt);
   }
   
   /* Recupera i dati che servono */
   Mat3x3 R1(pNode1->GetRRef());
   Mat3x3 R2(pNode2->GetRRef());
   Vec3 d2Tmp(R2*d2);
   Mat3x3 R1hTmp(R1*R1h);
   Mat3x3 R2hTmp(R2*R2h);
   
   Vec3 e3a(R1hTmp.GetVec(3));
   Vec3 e1b(R2hTmp.GetVec(1));
   Vec3 e2b(R2hTmp.GetVec(2));
   
   /* Suppongo che le reazioni F, M siano gia' state aggiornate da AssRes;
    * ricordo che la forza F, F(3), come la coppia M,
    * e' nel sistema locale ed il terzo termine, M(3), e' nullo in quanto
    * diretto come l'asse attorno al quale la rotazione e' consentita */
   
   Vec3 x2d2mx1(pNode2->GetXCurr()+d2Tmp-pNode1->GetXCurr());

   /* Contributo della forza alle equazioni di equilibrio dei due nodi */
   Vec3 Tmp1(e3a.Cross(x2d2mx1));
   Vec3 Tmp2(d2Tmp.Cross(e3a));
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = e3a.dGet(iCnt);
      WM.fPutCoef(iCnt, 13, -d);
      WM.fPutCoef(6+iCnt, 13, d);
      
      WM.fPutCoef(13, iCnt, -d);
      WM.fPutCoef(13, 6+iCnt, d);
      
      d = Tmp1.dGet(iCnt);
      WM.fPutCoef(3+iCnt, 13, d);
      
      WM.fPutCoef(13, 3+iCnt, d);
      
      d = Tmp2.dGet(iCnt);
      WM.fPutCoef(9+iCnt, 13, d);
      
      WM.fPutCoef(13, 9+iCnt, d);
   }
   
   Vec3 FTmp(F*dCoef);
   Mat3x3 FCross(FTmp);
   
   WM.Add(1, 4, FCross);
   WM.Sub(4, 1, FCross);
   
   WM.Add(4, 4, Mat3x3(x2d2mx1, FTmp));

#ifdef __GNUC__
#warning "rivedere!"
#endif
   WM.Add(4, 7, FCross);
   WM.Sub(7, 10, FCross);
   
   WM.Sub(10, 10, Mat3x3(d2Tmp, FTmp));
      
   /* Moltiplica la forza per il coefficiente del metodo */
   Vec3 MTmp(e2b*(M.dGet(1)*dCoef)-e1b*(M.dGet(2)*dCoef));
   
   Mat3x3 MWedgee3aWedge(MTmp, e3a);
   Mat3x3 e3aWedgeMWedge(e3a, MTmp);
   
   WM.Sub(4, 4, MWedgee3aWedge);
   WM.Add(4, 10, e3aWedgeMWedge);
   
   WM.Add(10, 4, MWedgee3aWedge);
   WM.Sub(10, 10, e3aWedgeMWedge);
   
   /* Contributo del momento alle equazioni di equilibrio dei nodi */
   Tmp1 = Vec3(e2b.Cross(e3a));
   Tmp2 = Vec3(e3a.Cross(e1b));
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = Tmp1.dGet(iCnt);
      WM.fPutCoef(3+iCnt, 14, d);
      WM.fPutCoef(9+iCnt, 14, -d);
      
      WM.fPutCoef(14, 3+iCnt, d);
      WM.fPutCoef(14, 9+iCnt, -d);
      
      d = Tmp2.dGet(iCnt);
      WM.fPutCoef(3+iCnt, 15, d);
      WM.fPutCoef(9+iCnt, 15, -d);
      
      WM.fPutCoef(15, 3+iCnt, -d);
      WM.fPutCoef(15, 9+iCnt, d);
   }
   
   return WorkMat;
}


/* Assemblaggio residuo */
SubVectorHandler& 
PlaneDispJoint::AssRes(SubVectorHandler& WorkVec,
		       doublereal dCoef,
		       const VectorHandler& XCurr, 
		       const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUTFNAME("PlaneDispJoint::AssRes");
      
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
      
   /* Indici */
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   
   /* Indici dei nodi */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {	
      WorkVec.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.fPutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
   }   
   
   /* Indici del vincolo */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WorkVec.fPutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
   }
   
   /* Recupera i dati */
   Vec3 x1(pNode1->GetXCurr());
   Vec3 x2(pNode2->GetXCurr());
   Mat3x3 R1(pNode1->GetRCurr());
   Mat3x3 R2(pNode2->GetRCurr());
   
   /* Costruisce i dati propri nella configurazione corrente */
   Vec3 d1Tmp(R1*d1);
   Vec3 d2Tmp(R2*d2);
   Mat3x3 R1hTmp(R1*R1h);
   Mat3x3 R2hTmp(R2*R2h);
   
   Vec3 e3a(R1hTmp.GetVec(3));
   Vec3 e1b(R2hTmp.GetVec(1));
   Vec3 e2b(R2hTmp.GetVec(2));
   
   /* Aggiorna i dati propri */
   F = e3a*XCurr.dGetCoef(iFirstReactionIndex+1);
   M = Vec3(XCurr.dGetCoef(iFirstReactionIndex+2),
	    XCurr.dGetCoef(iFirstReactionIndex+3),
	    0.);
   Vec3 MTmp(e2b.Cross(e3a)*M.dGet(1)+e3a.Cross(e1b)*M.dGet(2));
   Vec3 x2d2mx1(x2+d2Tmp-x1);
      
   /* Equazioni di equilibrio, nodo 1 */
   WorkVec.Add(1, F);
   WorkVec.Add(4, x2d2mx1.Cross(F)-MTmp);
   
   /* Equazioni di equilibrio, nodo 2 */
   WorkVec.Sub(7, F);
   WorkVec.Sub(10, d2Tmp.Cross(F)-MTmp);

   /* Modifica: divido le equazioni di vincolo per dCoef */
   if (dCoef != 0.) {

      /* Equazione di vincolo di posizione */
      WorkVec.fPutCoef(13, (e3a.Dot(x2d2mx1-d1Tmp))/dCoef);
      
      /* Equazioni di vincolo di rotazione */
      WorkVec.fPutCoef(14, e3a.Dot(e2b)/dCoef);
      WorkVec.fPutCoef(15, e3a.Dot(e1b)/dCoef);
   }
   
   return WorkVec;
}

/* Output (da mettere a punto) */
void 
PlaneDispJoint::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) {
      Mat3x3 R2Tmp(pNode2->GetRCurr()*R2h);
      Mat3x3 RTmp((pNode1->GetRCurr()*R1h).Transpose()*R2Tmp);
      Mat3x3 R2TmpT(R2Tmp.Transpose());
      
      Joint::Output(OH.Joints(), "PlaneDisp", GetLabel(),
		    R2TmpT*F, M, F, R2Tmp*M)
	<< " " << EulerAngles(RTmp) 
	  << " " << R2TmpT*(pNode1->GetWCurr()-pNode2->GetWCurr()) << endl;
   }   
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
PlaneDispJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
				   const VectorHandler& XCurr)
{
   DEBUGCOUTFNAME("PlaneDispJoint::InitialAssJac");
   
   /* Per ora usa la matrice piena; eventualmente si puo' 
    * passare a quella sparsa quando si ottimizza */
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeInit(iNumRows, iNumCols, 0.);

   /* Equazioni: vedi joints.dvi */
   
   /* Indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstVelIndex = iNode1FirstPosIndex+6;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstVelIndex = iNode2FirstPosIndex+6;
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+3;
   
   /* Nota: le reazioni vincolari sono:
    * Forza,       1 incognita, riferimento locale,
    * Momento,     2 incognite, riferimento locale
    */

   /* Setta gli indici dei nodi */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.fPutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.fPutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.fPutRowIndex(6+iCnt, iNode1FirstVelIndex+iCnt);
      WM.fPutColIndex(6+iCnt, iNode1FirstVelIndex+iCnt);
      WM.fPutRowIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
      WM.fPutColIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
      WM.fPutRowIndex(18+iCnt, iNode2FirstVelIndex+iCnt);
      WM.fPutColIndex(18+iCnt, iNode2FirstVelIndex+iCnt);
   }
   
   /* Setta gli indici delle reazioni */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.fPutRowIndex(24+iCnt, iFirstReactionIndex+iCnt);
      WM.fPutColIndex(24+iCnt, iFirstReactionIndex+iCnt);	
   }   
   
   /* Recupera i dati */
   Mat3x3 R1(pNode1->GetRRef());
   Mat3x3 R2(pNode2->GetRRef());
   Vec3 Omega1(pNode1->GetWRef());
   Vec3 Omega2(pNode2->GetWRef());
   
   /* Distanze e matrici di rotazione dai nodi alla cerniera 
    * nel sistema globale */
   Vec3 d1Tmp(R1*d1);
   Vec3 d2Tmp(R2*d2);
   Mat3x3 R1hTmp(R1*R1h);
   Mat3x3 R2hTmp(R2*R2h);
   
   Vec3 e3a(R1hTmp.GetVec(3));
   Vec3 e1b(R2hTmp.GetVec(1));
   Vec3 e2b(R2hTmp.GetVec(2));
   
   Vec3 MPrime(XCurr.dGetCoef(iReactionPrimeIndex+4),
	       XCurr.dGetCoef(iReactionPrimeIndex+5),
	       0.);
   
   Vec3 x2pqmx1(pNode2->GetXCurr()+d2Tmp-pNode1->GetXCurr());
   Vec3 xp2pqpmxp1(pNode2->GetVCurr()+Omega2.Cross(d2Tmp)-pNode1->GetVCurr());
   
   /* Aggiorna i dati propri */
   doublereal dFPrime = XCurr.dGetCoef(iReactionPrimeIndex+1);   
   Vec3 FPrime(e3a*dFPrime);

   Vec3 Tmp1(e3a.Cross(x2pqmx1));
   Vec3 Tmp2(Omega1.Cross(e3a));
   Vec3 Tmp3((Omega1.Cross(x2pqmx1)-xp2pqpmxp1).Cross(e3a));
   Vec3 Tmp4(-(xp2pqpmxp1.Cross(e3a)+x2pqmx1.Cross(e3a)));
   
   Vec3 Tmp5(d2Tmp.Cross(e3a));
   Vec3 Tmp6(d2Tmp.Cross(e3a.Cross(Omega2-Omega1)));
   Vec3 Tmp7(d2Tmp.Cross(Omega1.Cross(e3a))-e3a.Cross(Omega2.Cross(d2Tmp)));
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = e3a.dGet(iCnt);
      WM.fPutCoef(iCnt, 25, -d);    
      WM.fPutCoef(12+iCnt, 25, d);  
      
      WM.fPutCoef(25, iCnt, -d);   
      WM.fPutCoef(25, 12+iCnt, d);  

      WM.fPutCoef(6+iCnt, 28, -d); 
      WM.fPutCoef(18+iCnt, 28, d);  

      WM.fPutCoef(28, 6+iCnt, -d);  
      WM.fPutCoef(28, 18+iCnt, d);  

      d = Tmp1.dGet(iCnt);
      WM.fPutCoef(3+iCnt, 25, d);   
      WM.fPutCoef(25, 3+iCnt, d);   

      WM.fPutCoef(28, 9+iCnt, d);   
      WM.fPutCoef(9+iCnt, 28, d);   
	
      d = Tmp2.dGet(iCnt);
      WM.fPutCoef(28, iCnt, -d);    
      WM.fPutCoef(28, 12+iCnt, d);  
      
      WM.fPutCoef(6+iCnt, 25, -d);  
      WM.fPutCoef(18+iCnt, 25, d);  
      
      d = Tmp3.dGet(iCnt);
      WM.fPutCoef(28, 3+iCnt, d);   
      
      d = Tmp4.dGet(iCnt);
      WM.fPutCoef(9+iCnt, 25, d);   

      d = Tmp5.dGet(iCnt);
      WM.fPutCoef(15+iCnt, 25, d);   
      WM.fPutCoef(21+iCnt, 28, d);   

      WM.fPutCoef(25, 15+iCnt, d);   
      WM.fPutCoef(28, 21+iCnt, d);
      
      d = Tmp6.dGet(iCnt);
      WM.fPutCoef(28, 15+iCnt, d);   

      d = Tmp7.dGet(iCnt);
      WM.fPutCoef(21+iCnt, 25, d);   
   }   

   Mat3x3 MTmp(F);
   WM.Add(1, 4, MTmp);             
   WM.Add(4, 13, MTmp);            
   
   WM.Add(7, 10, MTmp);            
   WM.Add(10, 19, MTmp);           
   
   WM.Sub(4, 1, MTmp);
   WM.Sub(13, 4, MTmp);
   
   WM.Sub(19, 10, MTmp);           
   WM.Sub(10, 7, MTmp);            
   
   MTmp = Mat3x3(x2pqmx1, F);
   WM.Add(4, 4, MTmp);             

   WM.Add(10, 10, MTmp);
 
   MTmp = Mat3x3(Omega1, F)+Mat3x3(FPrime);
   WM.Add(7, 4, MTmp);   
   WM.Sub(19, 4, MTmp);
      
   MTmp = Mat3x3(Omega1.Cross(F)+FPrime);
   WM.Sub(10, 1, MTmp);           
   WM.Add(10, 13, MTmp);           
   
   MTmp = (Mat3x3(xp2pqpmxp1)+Mat3x3(x2pqmx1, Omega1))*Mat3x3(F)
     +Mat3x3(x2pqmx1, FPrime);
   WM.Add(10, 4, MTmp);    
   
   MTmp = Mat3x3(F, d2Tmp);
   WM.Add(16, 16, MTmp);
   WM.Add(22, 22, MTmp);
   
   WM.Sub(4, 16, MTmp);
   WM.Sub(10, 22, MTmp);
   
   MTmp = MTmp.Transpose();
   WM.Add(16, 4, MTmp);
   WM.Add(22, 10, MTmp);
   
   MTmp = (Mat3x3(F, Omega2)+Mat3x3(Omega1.Cross(F)+FPrime))*Mat3x3(d2Tmp);
   WM.Add(22, 16, MTmp);
   WM.Sub(10, 16, MTmp);
   
   MTmp = Mat3x3(d2Tmp.Cross(Omega2), F)
     -Mat3x3(d2Tmp)*(Mat3x3(Omega1, F)+Mat3x3(FPrime));
   WM.Add(22, 4, MTmp);   
      
   /* Ruota il momento e la sua derivata con le matrici della cerniera 
    * rispetto ai nodi */
   Vec3 MTmp2(e2b*M.dGet(1)-e1b*M.dGet(2));
   Vec3 MPrimeTmp(e2b*MPrime.dGet(1)-e1b*MPrime.dGet(2));

   Mat3x3 MDeltag1((Mat3x3(Omega2.Cross(MTmp2)+MPrimeTmp)+
		    Mat3x3(MTmp2, Omega1))*Mat3x3(e3a));
   Mat3x3 MDeltag2(Mat3x3(Omega1.Cross(e3a), MTmp2)+
		   Mat3x3(e3a, MPrimeTmp)+
		   Mat3x3(e3a)*Mat3x3(Omega2, MTmp2));

   /* Vettori temporanei */
   Tmp1 = e2b.Cross(e3a);
   Tmp2 = e3a.Cross(e1b);
   
   /* Prodotto vettore tra il versore 3 della cerniera secondo il nodo 1
    * ed il versore 1 della cerniera secondo il nodo 2. A convergenza
    * devono essere ortogonali, quindi il loro prodotto vettore deve essere 
    * unitario */

   /* Error handling: il programma si ferma, pero' segnala dov'e' l'errore */
   if (Tmp1.Dot() <= DBL_EPSILON || Tmp2.Dot() <= DBL_EPSILON) {
      cerr << "joint " << GetLabel() << ':' << endl
	<< "warning, first node hinge axis and second node hinge axis are (nearly) orthogonal;" << endl
	<< "aborting ..." << endl;
      THROW(Joint::ErrGeneric());
   }   
   
   Vec3 TmpPrime1(e2b.Cross(Omega1.Cross(e3a))-e3a.Cross(Omega2.Cross(e2b)));
   Vec3 TmpPrime2(e3a.Cross(Omega2.Cross(e1b))-e1b.Cross(Omega1.Cross(e3a)));
   
   /* Equazione di momento, nodo 1 */
   WM.Sub(4, 4, Mat3x3(MTmp2, e3a));
   WM.Add(4, 16, Mat3x3(e3a, MTmp2));

   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.fPutCoef(3+iCnt, 26, Tmp1.dGet(iCnt));
      WM.fPutCoef(3+iCnt, 27, Tmp2.dGet(iCnt));	
   }
   
   /* Equazione di momento, nodo 2 */
   WM.Add(16, 4, Mat3x3(MTmp2, e3a));
   WM.Sub(16, 16, Mat3x3(e3a, MTmp2));
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.fPutCoef(15+iCnt, 26, -Tmp1.dGet(iCnt));
      WM.fPutCoef(15+iCnt, 27, -Tmp2.dGet(iCnt));	
   }
   
   /* Derivata dell'equazione di momento, nodo 1 */
   WM.Sub(10, 4, MDeltag1);
   WM.Sub(10, 10, Mat3x3(MTmp2, e3a));
   WM.Add(10, 16, MDeltag2);
   WM.Add(10, 22, Mat3x3(e3a, MTmp2));
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.fPutCoef(9+iCnt, 26, TmpPrime1.dGet(iCnt));
      WM.fPutCoef(9+iCnt, 27, TmpPrime2.dGet(iCnt));
      WM.fPutCoef(9+iCnt, 29, Tmp1.dGet(iCnt));
      WM.fPutCoef(9+iCnt, 30, Tmp2.dGet(iCnt));	
   }
   
   /* Derivata dell'equazione di momento, nodo 2 */
   WM.Add(22, 4, MDeltag1);
   WM.Add(22, 10, Mat3x3(MTmp2, e3a));
   WM.Sub(22, 16, MDeltag2);
   WM.Sub(22, 22, Mat3x3(e3a, MTmp2));

   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.fPutCoef(21+iCnt, 26, -TmpPrime1.dGet(iCnt));
      WM.fPutCoef(21+iCnt, 27, -TmpPrime2.dGet(iCnt));	
      WM.fPutCoef(21+iCnt, 29, -Tmp1.dGet(iCnt));
      WM.fPutCoef(21+iCnt, 30, -Tmp2.dGet(iCnt));	
   }
   
   /* Equazioni di vincolo di rotazione: e1b~e3a, e2b~e3a */            
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      doublereal d = Tmp1.dGet(iCnt);
      WM.fPutCoef(26, 3+iCnt, d);
      WM.fPutCoef(26, 15+iCnt, -d);
      
      /* Queste sono per la derivata dell'equazione, sono qui solo per 
       * ottimizzazione */
      WM.fPutCoef(29, 9+iCnt, d);
      WM.fPutCoef(29, 21+iCnt, -d);
      
      d = Tmp2.dGet(iCnt);
      WM.fPutCoef(27, 3+iCnt, -d);
      WM.fPutCoef(27, 15+iCnt, d);
      
      /* Queste sono per la derivata dell'equazione, sono qui solo per 
       * ottimizzazione */
      WM.fPutCoef(30, 9+iCnt, -d);
      WM.fPutCoef(30, 21+iCnt, d);
   }   
   
   /* Derivate delle equazioni di vincolo di rotazione: e1b~e3a, e2b~e3a */
   Vec3 O1mO2(Omega1-Omega2);
   TmpPrime1 = e3a.Cross(O1mO2.Cross(e2b));   
   TmpPrime2 = e2b.Cross(e3a.Cross(O1mO2));
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WM.fPutCoef(29, 3+iCnt, TmpPrime1.dGet(iCnt));
      WM.fPutCoef(29, 15+iCnt, TmpPrime2.dGet(iCnt));
   }
   
   TmpPrime1 = e3a.Cross(O1mO2.Cross(e1b));
   TmpPrime2 = e1b.Cross(e3a.Cross(O1mO2));
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WM.fPutCoef(30, 3+iCnt, TmpPrime1.dGet(iCnt));
      WM.fPutCoef(30, 15+iCnt, TmpPrime2.dGet(iCnt));
   }   
   
   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
PlaneDispJoint::InitialAssRes(SubVectorHandler& WorkVec,
				   const VectorHandler& XCurr)
{   
   DEBUGCOUTFNAME("PlaneDispJoint::InitialAssRes");
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
   
   /* Indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstVelIndex = iNode1FirstPosIndex+6;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstVelIndex = iNode2FirstPosIndex+6;
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+3;
   
   /* Setta gli indici */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {	
      WorkVec.fPutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WorkVec.fPutRowIndex(6+iCnt, iNode1FirstVelIndex+iCnt);
      WorkVec.fPutRowIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
      WorkVec.fPutRowIndex(18+iCnt, iNode2FirstVelIndex+iCnt);
   }
   
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.fPutRowIndex(24+iCnt, iFirstReactionIndex+iCnt);
   }

   /* Dati */   
   Vec3 x1(pNode1->GetXCurr());
   Vec3 x2(pNode2->GetXCurr());
   Vec3 v1(pNode1->GetVCurr());
   Vec3 v2(pNode2->GetVCurr());
   Mat3x3 R1(pNode1->GetRCurr());
   Mat3x3 R2(pNode2->GetRCurr());
   Vec3 Omega1(pNode1->GetWCurr());
   Vec3 Omega2(pNode2->GetWCurr());
   
   /* Aggiorna F ed M, che restano anche per InitialAssJac */
   M = Vec3(XCurr.dGetCoef(iFirstReactionIndex+2),
	    XCurr.dGetCoef(iFirstReactionIndex+3),
	    0.);
   Vec3 MPrime(XCurr.dGetCoef(iReactionPrimeIndex+2),
	       XCurr.dGetCoef(iReactionPrimeIndex+3),
	       0.);
   
   /* Distanza nel sistema globale */
   Vec3 d1Tmp(R1*d1);
   Vec3 d2Tmp(R2*d2);
   Mat3x3 R1hTmp(R1*R1h);
   Mat3x3 R2hTmp(R2*R2h);

   Vec3 e3a(R1hTmp.GetVec(3));
   Vec3 e1b(R2hTmp.GetVec(1));
   Vec3 e2b(R2hTmp.GetVec(2));  
   
   Vec3 x2pqmx1(pNode2->GetXCurr()+d2Tmp-pNode1->GetXCurr());
   Vec3 xp2pqpmxp1(pNode2->GetVCurr()+Omega2.Cross(d2Tmp)-pNode1->GetVCurr());
   
   /* Aggiorna i dati propri */
   F = e3a*XCurr.dGetCoef(iFirstReactionIndex+1);
   doublereal dFPrime = XCurr.dGetCoef(iReactionPrimeIndex+1);
   Vec3 FPrime(e3a*dFPrime);
   Vec3 Tmp(Omega1.Cross(F)+FPrime);
   
   WorkVec.Add(1, F);
   WorkVec.Add(4, x2pqmx1.Cross(F));
   WorkVec.Add(7, Tmp);
   WorkVec.Add(10, xp2pqpmxp1.Cross(F)+x2pqmx1.Cross(Tmp));
   WorkVec.Sub(13, F);
   WorkVec.Add(16, F.Cross(d2Tmp));
   WorkVec.Sub(19, Tmp);
   WorkVec.Add(22, F.Cross(Omega2.Cross(d2Tmp))-d2Tmp.Cross(Tmp));   
   
   WorkVec.fPutCoef(25, e3a.Dot(d1Tmp-x2pqmx1));
   WorkVec.fPutCoef(28, x2pqmx1.Dot(e3a.Cross(Omega1))-e3a.Dot(xp2pqpmxp1));
   
   /* Ruota il momento e la sua derivata con le matrici della cerniera 
    * rispetto ai nodi */
   Vec3 MTmp(e2b*M.dGet(1)-e1b*M.dGet(2));       
   Vec3 MPrimeTmp(e3a.Cross(MTmp.Cross(Omega2))+MTmp.Cross(Omega1.Cross(e3a))+
		  e2b.Cross(e3a)*MPrime.dGet(1)+e3a.Cross(e1b)*MPrime.dGet(2)); 

   /* Equazioni di equilibrio, nodo 1 */
   WorkVec.Sub(4, MTmp.Cross(e3a));
   
   /* Derivate delle equazioni di equilibrio, nodo 1 */
   WorkVec.Sub(10, MPrimeTmp);
   
   /* Equazioni di equilibrio, nodo 2 */
   WorkVec.Add(16, MTmp.Cross(e3a)); 
   
   /* Derivate delle equazioni di equilibrio, nodo 2 */
   WorkVec.Add(22, MPrimeTmp);
   
   /* Equazioni di vincolo di rotazione */
   WorkVec.fPutCoef(26, e2b.Dot(e3a));
   WorkVec.fPutCoef(27, e1b.Dot(e3a));
   
   /* Derivate delle equazioni di vincolo di rotazione: e1b~e3a, e2b~e3a */
   Tmp = (Omega1-Omega2).Cross(e3a);
   WorkVec.fPutCoef(29, e2b.Dot(Tmp));
   WorkVec.fPutCoef(30, e1b.Dot(Tmp));

   return WorkVec;
}


doublereal 
PlaneDispJoint::dGetPrivData(unsigned int i) const
{
   ASSERT(i >= 1 && i <= 3);
   
   double d;
   
   Mat3x3 R2Tmp(pNode2->GetRCurr()*R2h);
   Mat3x3 RTmp((pNode1->GetRCurr()*R1h).Transpose()*R2Tmp);
   
   switch (i) {
    case 1:      
    case 2: {
       Vec3 x1 = pNode1->GetXCurr()+pNode1->GetRCurr()*d1;
       Vec3 x2 = pNode2->GetXCurr()+pNode2->GetRCurr()*d2;
       
       Vec3 v(R1h.Transpose()*(pNode1->GetRCurr().Transpose()*(x2-x1)));
       
       d = v.dGet(i);
       break;
    }
            
    case 3: {
       Vec3 v(EulerAngles(RTmp));
       
       d = v.dGet(3);
       break;
    }
      
    case 4:
    case 5: {
       Vec3 v1(pNode1->GetVCurr()
	       +pNode1->GetWCurr().Cross(pNode1->GetRCurr()*d1));
       Vec3 v2(pNode2->GetVCurr()
	       +pNode2->GetWCurr().Cross(pNode2->GetRCurr()*d2));
       
       Vec3 v(R1h.Transpose()*(pNode1->GetRCurr().Transpose()*(v2-v1)));
      
       d = v.dGet(i);
       break;       
    }
      
    case 6: {
       Mat3x3 R2TmpT(R2Tmp.Transpose());
       Vec3 v(R2TmpT*(pNode1->GetWCurr()-pNode2->GetWCurr()));
       
       d = v.dGet(i-3);
       break;
    }
      
    default:
      cerr << "Illegal private data" << endl;
      THROW(ErrGeneric());
   }
   
   return d;
}

/* PlaneDispJoint - end */


/* PlaneDispPinJoint - begin */

/* Costruttore non banale */
PlaneDispPinJoint::PlaneDispPinJoint(unsigned int uL, 
				     const DofOwner* pDO,	       
				     const StructNode* pN,
				     const Vec3& X0Tmp, 
				     const Mat3x3& R0Tmp, 
				     const Vec3& dTmp, 
				     const Mat3x3& RhTmp,
				     flag fOut)
: Elem(uL, ElemType::JOINT, fOut), 
Joint(uL, JointType::PLANEDISPPIN, pDO, fOut), 
pNode(pN), 
X0(X0Tmp), R0(R0Tmp), d(dTmp), Rh(RhTmp),
F(0.), M(0.)
{
   NO_OP;
}


/* Distruttore banale */
PlaneDispPinJoint::~PlaneDispPinJoint(void)
{
   NO_OP;
};


/* Contributo al file di restart */
ostream& PlaneDispPinJoint::Restart(ostream& out) const
{
   Joint::Restart(out) << ", plane pin, "
     << pNode->GetLabel() 
     << ", reference, node, ", d.Write(out, ", ") 
     << ", hinge, reference, node, 1, ", 
     (Rh.GetVec(1)).Write(out, ", ") << ", 2, ", 
     (Rh.GetVec(2)).Write(out, ", ") 
     << ", reference, global, ", X0.Write(out, ", ") 
     << ", reference, global, 1, ",
     (R0.GetVec(1)).Write(out, ", ") << ", 2, ", 
     (R0.GetVec(2)).Write(out, ", ") << ';' << endl;
   
   return out;
}


/* Assemblaggio jacobiano */
VariableSubMatrixHandler& 
PlaneDispPinJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		      doublereal dCoef,
		      const VectorHandler& /* XCurr */ ,
		      const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUTFNAME("PlaneDispPinJoint::AssJac");
      
   SparseSubMatrixHandler& WM = WorkMat.SetSparse();
   WM.ResizeInit(39, 0, 0.);
   
   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();

   Mat3x3 R(pNode->GetRRef());
   Vec3 dTmp(R*d);
   Mat3x3 RhTmp(R*Rh);
   
   
   /* 
    * L'equazione di vincolo afferma che il punto in cui si trova la
    * cerniera deve essere fissato:
    *      x + d = x0
    *      e3_0^Te1 = 0
    *      e3_0^Te2 = 0
    * 
    * con: d = R * d_0
    * 
    * La forza e' data dalla reazione vincolare F, nel sistema globale
    * La coppia dovuta all'eccentricita' e' data rispettivamente da:
    *     d /\ F
    *
    * 
    *       x      g         F
    * Q1 |  0      0             I   0    | | x |   | -F          |
    * G1 |  0      cF/\d1/\-M/\  d/\ e1e2 | | g |   | -d/\F-M     |
    * F  |  I      d/\           0   0    | | F |   |  (x+d-x0)/c |
    * M  |  0      e_0/\e1,e2    0   0    | | M |   |  e_0^Te1,e2 |
    * 
    * con d = R*d_0, c = dCoef
    */


   
   /* Moltiplica la forza ed il momento per il coefficiente
    * del metodo */

   Vec3 e3(R0.GetVec(3));
   Vec3 e1(RhTmp.GetVec(1));
   Vec3 e2(RhTmp.GetVec(2));
   Vec3 MTmp(e2*M.dGet(1)-e1*M.dGet(2));
            
   Vec3 Tmp1((e2).Cross(e3));
   Vec3 Tmp2((e3).Cross(e1));
   
   /* termini di reazione sul nodo (forza e momento) */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.fPutItem(iCnt, iFirstMomentumIndex+iCnt, 
		  iFirstReactionIndex+iCnt, 1.);
      WM.fPutItem(3+iCnt, 3+iFirstMomentumIndex+iCnt,
		  iFirstReactionIndex+4, Tmp1.dGet(iCnt));
      WM.fPutItem(6+iCnt, 3+iFirstMomentumIndex+iCnt, 
		  iFirstReactionIndex+5, Tmp2.dGet(iCnt));
   }   
   
   WM.fPutCross(10, iFirstMomentumIndex+3,
		iFirstReactionIndex, dTmp);
      
   
   /* Nota: F ed M, le reazioni vincolari, sono state aggiornate da AssRes */
   
   /* Termini diagonali del tipo: c*F/\d/\Delta_g 
    * nota: la forza e' gia' moltiplicata per dCoef */      
   WM.fPutMat3x3(16, iFirstMomentumIndex+3, iFirstPositionIndex+3, 
		 Mat3x3(F*dCoef, dTmp)+Mat3x3(e3, MTmp*dCoef));

   /* Modifica: divido le equazioni di vincolo per dCoef */
   
   /* termini di vincolo dovuti al nodo 1 */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.fPutItem(24+iCnt, iFirstReactionIndex+iCnt, 
		  iFirstPositionIndex+iCnt, -1.);
   }
   
   WM.fPutCross(28, iFirstReactionIndex,
		iFirstPositionIndex+3, dTmp);
   
   for (int iCnt = 1; iCnt <= 3; iCnt ++) {
      WM.fPutItem(33+iCnt, iFirstReactionIndex+4, 
		  iFirstPositionIndex+3+iCnt, Tmp1.dGet(iCnt));	
      WM.fPutItem(36+iCnt, iFirstReactionIndex+5, 
		  iFirstPositionIndex+3+iCnt, -Tmp2.dGet(iCnt));	
   }
   
   return WorkMat;
}


/* Assemblaggio residuo */
SubVectorHandler& PlaneDispPinJoint::AssRes(SubVectorHandler& WorkVec,
					      doublereal dCoef,
					      const VectorHandler& XCurr,
					      const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUTFNAME("PlaneDispPinJoint::AssRes");
      
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
        
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   
   /* Indici dei nodi */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iFirstMomentumIndex+iCnt);
   }
     
   
   /* Indici del vincolo */
   for (int iCnt = 1; iCnt <= 5; iCnt++) {
      WorkVec.fPutRowIndex(6+iCnt, iFirstReactionIndex+iCnt);   
   }
   
   F = Vec3(XCurr, iFirstReactionIndex+1);
   M = Vec3(XCurr.dGetCoef(iFirstReactionIndex+4),
	    XCurr.dGetCoef(iFirstReactionIndex+5),
	    0.);
   
   Vec3 x(pNode->GetXCurr());
   Mat3x3 R(pNode->GetRCurr());
   
   Vec3 dTmp(R*d);
   Mat3x3 RhTmp(R*Rh);
   
   Vec3 e3(R0.GetVec(3));
   Vec3 e1(RhTmp.GetVec(1));
   Vec3 e2(RhTmp.GetVec(2));
   
   WorkVec.Add(1, -F);
   WorkVec.Add(4, F.Cross(dTmp)-(e2*M.dGet(1)-e1*M.dGet(2)).Cross(e3)); /* Sfrutto il fatto che F/\d = -d/\F */
   
   /* Modifica: divido le equazioni di vincolo per dCoef */
   if (dCoef != 0.) {	
      WorkVec.Add(7, (x+dTmp-X0)/dCoef);
      
      WorkVec.fPutCoef(10, e3.Dot(e2)/dCoef);
      WorkVec.fPutCoef(11, e3.Dot(e1)/dCoef);
   }   
   
   return WorkVec;
}

/* Output (da mettere a punto) */
void PlaneDispPinJoint::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) {
#ifdef DEBUG
      OH.Output() << "Joint " << uLabel << ", type \""
	<< psJointNames[JointType::PLANEPIN]
	<< "\", linked to node " << pNode->GetLabel() << ':' << endl
	<< "Distance from node (node reference frame): " << endl << d << endl
	<< "Hinge rotation matrix (node reference frame): " << endl << Rh << endl
	<< "Current reaction force: " << endl << F << endl
	<< "Current reaction couple (hinge reference frame):" << endl 
	<< M << endl;
#endif
      
      Mat3x3 RTmp(pNode->GetRCurr()*Rh);
      Mat3x3 RTmpT(RTmp.Transpose());
      Mat3x3 R0Tmp(R0.Transpose()*RTmp);
      
      Joint::Output(OH.Joints(), "PlaneDispPin", GetLabel(),
		    RTmpT*F, M, F, RTmp*M) 
	<< " " << EulerAngles(R0Tmp)
	<< " " << RTmpT*(pNode->GetWCurr()) << endl;      
   }
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
PlaneDispPinJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
				   const VectorHandler& XCurr)
{
   DEBUGCOUTFNAME("PlaneDispPinJoint::InitialAssJac");

   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeInit(iNumRows, iNumCols, 0.);

   /* Equazioni: vedi joints.dvi */
    
   /* Indici */
   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
   integer iFirstVelocityIndex = iFirstPositionIndex+6;
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+5;

   /* Setto gli indici */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.fPutRowIndex(iCnt, iFirstPositionIndex+iCnt);
      WM.fPutColIndex(iCnt, iFirstPositionIndex+iCnt);
      WM.fPutRowIndex(6+iCnt, iFirstVelocityIndex+iCnt);
      WM.fPutColIndex(6+iCnt, iFirstVelocityIndex+iCnt);
   }
   
   for (int iCnt = 1; iCnt <= 10; iCnt++) {
      WM.fPutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
      WM.fPutColIndex(12+iCnt, iFirstReactionIndex+iCnt);
   }   
   
   /* Matrici identita' */
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      /* Contributo di forza all'equazione della forza */
      WM.fPutCoef(iCnt, 12+iCnt, 1.);
      
      /* Contrib. di der. di forza all'eq. della der. della forza */
      WM.fPutCoef(6+iCnt, 17+iCnt, 1.);
      
      /* Equazione di vincolo */
      WM.fPutCoef(12+iCnt, iCnt, -1.);
      
      /* Derivata dell'equazione di vincolo */
      WM.fPutCoef(17+iCnt, 6+iCnt, -1.);
   }
   
   /* Recupera i dati */
   Mat3x3 R(pNode->GetRRef());
   Vec3 Omega(pNode->GetWRef());
   /* F, M sono state aggiornate da InitialAssRes */
   Vec3 FPrime(XCurr, iReactionPrimeIndex+1);
   Vec3 MPrime(XCurr.dGetCoef(iReactionPrimeIndex+4),
	       XCurr.dGetCoef(iReactionPrimeIndex+5),
	       0.);
   
   /* Distanza nel sistema globale */
   Vec3 dTmp(R*d);
   Mat3x3 RhTmp(R*Rh);

   Vec3 e3(R0.GetVec(3));
   Vec3 e1(RhTmp.GetVec(1));
   Vec3 e2(RhTmp.GetVec(2));

   /* Vettori temporanei */
   Vec3 Tmp1(e2.Cross(e3));   
   Vec3 Tmp2(e3.Cross(e1));
   
   /* Prodotto vettore tra il versore 3 della cerniera secondo il nodo 1
    * ed il versore 1 della cerniera secondo il nodo 2. A convergenza
    * devono essere ortogonali, quindi il loro prodotto vettore deve essere 
    * unitario */

   /* Error handling: il programma si ferma, pero' segnala dov'e' l'errore */
   if (Tmp1.Dot() <= DBL_EPSILON || Tmp2.Dot() <= DBL_EPSILON) {
      cerr << "joint " << GetLabel() << ':' << endl
	<< "warning, node hinge axis and fixed point hinge axis are (nearly) orthogonal;" << endl
	<< "aborting ..." << endl;
      THROW(Joint::ErrGeneric());
   }   
   
   Vec3 TmpPrime1(e3.Cross(e2.Cross(Omega)));
   Vec3 TmpPrime2(e3.Cross(Omega.Cross(e1)));   
   
   /* Ruota il momento e la sua derivata con le matrici della cerniera 
    * rispetto ai nodi */
   Vec3 MTmp(e2*M.dGet(1)-e1*M.dGet(2));
   Vec3 MPrimeTmp(e2*MPrime.dGet(1)-e1*MPrime.dGet(2));
         
   Mat3x3 MDeltag(Mat3x3(e3, MPrimeTmp)+Mat3x3(e3)*Mat3x3(Omega, MTmp));
   
   /* Matrici F/\d/\ */
   Mat3x3 FWedgedWedge(F, dTmp);
   
   /* Matrici (omega/\d)/\ */
   Mat3x3 OWedgedWedge(Omega.Cross(dTmp));
   
   /* Equazione di momento */
   WM.Add(4, 4, FWedgedWedge+Mat3x3(e3, MTmp));
   WM.Add(4, 13, Mat3x3(dTmp));
   
   /* Derivata dell'equazione di momento */
   WM.Add(10, 4, (Mat3x3(FPrime)+Mat3x3(F, Omega))*Mat3x3(dTmp)+MDeltag);
   WM.Add(10, 10, FWedgedWedge+Mat3x3(e3, MTmp));
   WM.Add(10, 13, OWedgedWedge);
   WM.Add(10, 18, Mat3x3(dTmp));
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = Tmp1.dGet(iCnt);
      WM.fPutCoef(3+iCnt, 16, d);
      WM.fPutCoef(9+iCnt, 21, d);
      
      d = Tmp2.dGet(iCnt);
      WM.fPutCoef(3+iCnt, 17, d);
      WM.fPutCoef(9+iCnt, 22, d);
      
      WM.fPutCoef(9+iCnt, 16, TmpPrime1.dGet(iCnt));
      WM.fPutCoef(9+iCnt, 17, TmpPrime2.dGet(iCnt));
   }
   
   /* Equazione di vincolo */
   WM.Add(13, 4, Mat3x3(dTmp));
   
   /* Derivata dell'equazione di vincolo */
   WM.Add(18, 4, OWedgedWedge);
   WM.Add(18, 10, Mat3x3(dTmp));

   /* Equazioni di vincolo di rotazione: e1b~e3a, e2b~e3a */            
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      doublereal d = -Tmp1.dGet(iCnt);
      WM.fPutCoef(16, 3+iCnt, d);
      
      /* Queste sono per la derivata dell'equazione, sono qui solo per 
       * ottimizzazione */
      WM.fPutCoef(21, 9+iCnt, d);
      
      d = Tmp2.dGet(iCnt);
      WM.fPutCoef(17, 3+iCnt, d);
      
      /* Queste sono per la derivata dell'equazione, sono qui solo per 
       * ottimizzazione */
      WM.fPutCoef(22, 9+iCnt, d);
   }   
   
   /* Derivate delle equazioni di vincolo di rotazione: e1b~e3a, e2b~e3a */
   TmpPrime2 = e2.Cross(Omega.Cross(e3));
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WM.fPutCoef(21, 3+iCnt, TmpPrime2.dGet(iCnt));
   }
   
   TmpPrime2 = e1.Cross(Omega.Cross(e3));
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WM.fPutCoef(22, 3+iCnt, TmpPrime2.dGet(iCnt));
   }      
   
   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
PlaneDispPinJoint::InitialAssRes(SubVectorHandler& WorkVec,
			const VectorHandler& XCurr)
{   
   DEBUGCOUTFNAME("PlaneDispPinJoint::InitialAssRes");
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
   
   /* Indici */
   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
   integer iFirstVelocityIndex = iFirstPositionIndex+6;
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+5;
   
   /* Setta gli indici */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {	
      WorkVec.fPutRowIndex(iCnt, iFirstPositionIndex+iCnt);
      WorkVec.fPutRowIndex(6+iCnt, iFirstVelocityIndex+iCnt);
   }
   
   for (int iCnt = 1; iCnt <= 10; iCnt++) {	
      WorkVec.fPutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
   }
   
   /* Recupera i dati */
   Vec3 x(pNode->GetXCurr());
   Vec3 v(pNode->GetVCurr());
   Mat3x3 R(pNode->GetRCurr());
   Vec3 Omega(pNode->GetWCurr());
   
   Mat3x3 RhTmp(R*Rh);
   
   F = Vec3(XCurr, iFirstReactionIndex+1);
   M = Vec3(XCurr.dGetCoef(iFirstReactionIndex+4),
	    XCurr.dGetCoef(iFirstReactionIndex+5),
	    0.);
   Vec3 FPrime(XCurr, iReactionPrimeIndex+1);
   Vec3 MPrime(XCurr.dGetCoef(iReactionPrimeIndex+4),
	       XCurr.dGetCoef(iReactionPrimeIndex+5),
	       0.);
   
   /* Versori delle cerniere */
   Vec3 e3(R0.GetVec(3));
   Vec3 e1(RhTmp.GetVec(1));
   Vec3 e2(RhTmp.GetVec(2));

   /* Vettori temporanei */
   Vec3 Tmp1(e2.Cross(e3));   
   Vec3 Tmp2(e3.Cross(e1));
      
   Vec3 TmpPrime1(e3.Cross(e2.Cross(Omega)));
   Vec3 TmpPrime2(e3.Cross(Omega.Cross(e1)));   
   
   /* Distanza nel sistema globale */
   Vec3 dTmp(R*d);

   /* Vettori omega/\d */
   Vec3 OWedged(Omega.Cross(dTmp));
   
   /* Ruota il momento e la sua derivata con le matrici della cerniera 
    * rispetto ai nodi */
   Vec3 MTmp(e2*M.dGet(1)-e1*M.dGet(2));       
   Vec3 MPrimeTmp(e3.Cross(MTmp.Cross(Omega))+
		  e2.Cross(e3)*MPrime.dGet(1)+e3.Cross(e1)*MPrime.dGet(2)); 
   
   /* Equazioni di equilibrio */
   WorkVec.Add(1, -F);
   WorkVec.Add(4, F.Cross(dTmp)-MTmp.Cross(e3)); /* Sfrutto il fatto che F/\d = -d/\F */
   
   /* Derivate delle equazioni di equilibrio, nodo 1 */
   WorkVec.Add(7, -FPrime);
   WorkVec.Add(10, FPrime.Cross(dTmp)-OWedged.Cross(F)-MPrimeTmp);
   
   /* Equazione di vincolo di posizione */
   WorkVec.Add(13, x+dTmp-X0);
   
   /* Equazioni di vincolo di rotazione */
   WorkVec.fPutCoef(16, e2.Dot(e3));
   WorkVec.fPutCoef(17, e1.Dot(e3));

   /* Derivata dell'equazione di vincolo di posizione */
   WorkVec.Add(18, v+OWedged);
   
   /* Derivate delle equazioni di vincolo di rotazione: e1b~e3a, e2b~e3a */
   Vec3 Tmp(e3.Cross(Omega));
   WorkVec.fPutCoef(21, e2.Dot(Tmp));
   WorkVec.fPutCoef(22, e1.Dot(Tmp));
      
   return WorkVec;
}

/* PlaneDispPinJoint - end */
