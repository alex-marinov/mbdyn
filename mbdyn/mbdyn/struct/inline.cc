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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <inline.h>

/* InLineJoint - begin */

/* Costruttore non banale */
InLineJoint::InLineJoint(unsigned int uL, const DofOwner* pDO,
			 const StructNode* pN1, const StructNode* pN2, 
			 const Mat3x3& RvTmp, const Vec3& pTmp, flag fOut)
: Elem(uL, Elem::JOINT, fOut), Joint(uL, Joint::INLINE, pDO, fOut),
pNode1(pN1), pNode2(pN2), Rv(RvTmp), p(pTmp), F(0.)
{
   NO_OP;
};

   
InLineJoint::~InLineJoint(void)
{
   NO_OP;
}


/* Contributo al file di restart */
ostream& InLineJoint::Restart(ostream& out) const
{
   Joint::Restart(out) << ", in line, not implemented yet" << endl;
   return out;
}

   
VariableSubMatrixHandler& 
InLineJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		    doublereal dCoef,
		    const VectorHandler& /* XCurr */ ,
		    const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUTFNAME("InLineJoint::AssJac");
   SparseSubMatrixHandler& WM = WorkMat.SetSparse();
   
   WM.ResizeInit(69, 0, 0.);
   
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   
   Vec3 x2mx1(pNode2->GetXCurr()-pNode1->GetXCurr());
   Mat3x3 RvTmp(pNode1->GetRRef()*Rv);
   Vec3 FTmp(RvTmp*(F*dCoef));
   
   
   Vec3 Tmp1_1(RvTmp.GetVec(1).Cross(x2mx1));
   Vec3 Tmp1_2(RvTmp.GetVec(2).Cross(x2mx1));
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = RvTmp.dGet(iCnt, 1);
      WM.fPutItem(iCnt, iNode1FirstMomIndex+iCnt,
		  iFirstReactionIndex+1, -d);
      WM.fPutItem(3+iCnt, iNode2FirstMomIndex+iCnt,
		  iFirstReactionIndex+1, d);      
      
      WM.fPutItem(6+iCnt, iFirstReactionIndex+1,
		  iNode1FirstPosIndex+iCnt, -d);
      WM.fPutItem(9+iCnt, iFirstReactionIndex+1,
		  iNode2FirstPosIndex+iCnt, d);
      
      d = RvTmp.dGet(iCnt, 2);
      WM.fPutItem(12+iCnt, iNode1FirstMomIndex+iCnt,
		  iFirstReactionIndex+2, -d);
      WM.fPutItem(15+iCnt, iNode2FirstMomIndex+iCnt,
		  iFirstReactionIndex+2, d);      
      
      WM.fPutItem(18+iCnt, iFirstReactionIndex+2,
		  iNode1FirstPosIndex+iCnt, -d);
      WM.fPutItem(21+iCnt, iFirstReactionIndex+2,
		  iNode2FirstPosIndex+iCnt, d);

      d = Tmp1_1.dGet(iCnt);
      WM.fPutItem(24+iCnt, iNode1FirstMomIndex+3+iCnt,
		  iFirstReactionIndex+1, d);
      
      WM.fPutItem(27+iCnt, iFirstReactionIndex+1,
		  iNode1FirstPosIndex+3+iCnt, d);
      
      d = Tmp1_2.dGet(iCnt);
      WM.fPutItem(30+iCnt, iNode1FirstMomIndex+3+iCnt,
		  iFirstReactionIndex+2, d);
      
      WM.fPutItem(33+iCnt, iFirstReactionIndex+2,
		  iNode1FirstPosIndex+3+iCnt, d);
   }
   
   WM.fPutCross(36+1, iNode1FirstMomIndex,
		iNode1FirstPosIndex+3, FTmp);
   WM.fPutCross(42+1, iNode1FirstMomIndex+3,
		iNode1FirstPosIndex, -FTmp);
   WM.fPutMat3x3(48+1, iNode1FirstMomIndex+3,		
		iNode1FirstPosIndex+3, Mat3x3(x2mx1, FTmp));
   WM.fPutCross(57+1, iNode1FirstMomIndex+3,
		iNode2FirstPosIndex, FTmp);
   WM.fPutCross(63+1, iNode2FirstMomIndex,
		iNode2FirstPosIndex+3, -FTmp);
   
   return WorkMat;
}


SubVectorHandler& 
InLineJoint::AssRes(SubVectorHandler& WorkVec,
		    doublereal dCoef,
		    const VectorHandler& XCurr, 
		    const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUTFNAME("InLineJoint::AssRes");
   WorkVec.Resize(14);
   WorkVec.Reset(0.);
 
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();

   /* Indici equazioni */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.fPutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
   }
   
   /* Indice equazione vincolo */
   WorkVec.fPutRowIndex(13, iFirstReactionIndex+1);
   WorkVec.fPutRowIndex(14, iFirstReactionIndex+2);

   Vec3 x2mx1(pNode2->GetXCurr()-pNode1->GetXCurr());
   Mat3x3 RvTmp = pNode1->GetRCurr()*Rv;

   /* Aggiorna i dati propri */
   F.Put(1, XCurr.dGetCoef(iFirstReactionIndex+1));
   F.Put(2, XCurr.dGetCoef(iFirstReactionIndex+2));
   Vec3 FTmp(RvTmp*F);
   
   WorkVec.Add(1, FTmp);
   WorkVec.Add(4, x2mx1.Cross(FTmp)); /* ( = -p/\F) */
   WorkVec.Sub(7, FTmp);
   
   if (dCoef != 0.) {
      WorkVec.fPutCoef(13, (Rv.GetVec(1).Dot(p)-RvTmp.GetVec(1).Dot(x2mx1))/dCoef);
      WorkVec.fPutCoef(14, (Rv.GetVec(2).Dot(p)-RvTmp.GetVec(2).Dot(x2mx1))/dCoef);
   }
      
   return WorkVec;
}


void InLineJoint::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) {      
      Mat3x3 RvTmp(pNode1->GetRCurr()*Rv);
      Joint::Output(OH.Joints(), "inline", GetLabel(),
		    F, Zero3, RvTmp*F, Zero3) << endl;
   }   
}
 

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
InLineJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
			    const VectorHandler& XCurr)
{
   DEBUGCOUTFNAME("InLineJoint::InitialAssJac");
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.ResizeInit(28, 28, 0.);
   
   /* Indici gdl */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+2;
   
   /* Indici equazioni nodi */
   for(int iCnt = 1; iCnt <= 12; iCnt++) {
      WM.fPutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.fPutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.fPutRowIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
      WM.fPutColIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
   }
   
   /* Indici vincoli */
   for(int iCnt = 1; iCnt <= 4; iCnt++) {
      WM.fPutRowIndex(24+iCnt, iFirstReactionIndex+iCnt);
      WM.fPutColIndex(24+iCnt, iFirstReactionIndex+iCnt);
   }
   
   /* Dati */
   Vec3 x2mx1(pNode2->GetXCurr()-pNode1->GetXCurr());
   Vec3 xp2mxp1(pNode2->GetVCurr()-pNode1->GetVCurr());
   Vec3 Omega(pNode1->GetWRef());
   
   /* Aggiorna i dati propri */
   Vec3 FPrime(XCurr.dGetCoef(iReactionPrimeIndex+1),
	       XCurr.dGetCoef(iReactionPrimeIndex+2),
	       0.);
   Mat3x3 RvTmp(pNode1->GetRRef()*Rv);
   Vec3 FTmp(RvTmp*F);
   Vec3 FPrimeTmp(RvTmp*FPrime);

   Vec3 v1(RvTmp.GetVec(1));
   Vec3 Tmp1_1(v1.Cross(x2mx1));
   Vec3 Tmp2_1(Omega.Cross(v1));
   Vec3 Tmp3_1((Omega.Cross(x2mx1)-xp2mxp1).Cross(v1));
   Vec3 Tmp4_1(-(xp2mxp1.Cross(v1)+x2mx1.Cross(Tmp2_1)));
   Vec3 v2(RvTmp.GetVec(2));
   Vec3 Tmp1_2(v2.Cross(x2mx1));
   Vec3 Tmp2_2(Omega.Cross(v2));
   Vec3 Tmp3_2((Omega.Cross(x2mx1)-xp2mxp1).Cross(v2));
   Vec3 Tmp4_2(-(xp2mxp1.Cross(v2)+x2mx1.Cross(Tmp2_2)));
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = v1.dGet(iCnt);
      WM.fPutCoef(iCnt, 25, -d);
      WM.fPutCoef(12+iCnt, 25, d);
      
      WM.fPutCoef(25, iCnt, -d);
      WM.fPutCoef(25, 12+iCnt, d);

      WM.fPutCoef(6+iCnt, 27, -d);  
      WM.fPutCoef(18+iCnt, 27, d);

      WM.fPutCoef(27, 6+iCnt, -d);  
      WM.fPutCoef(27, 18+iCnt, d);
      
      d = v2.dGet(iCnt);
      WM.fPutCoef(iCnt, 26, -d);
      WM.fPutCoef(12+iCnt, 26, d);
      
      WM.fPutCoef(26, iCnt, -d);
      WM.fPutCoef(26, 12+iCnt, d);

      WM.fPutCoef(6+iCnt, 28, -d);  
      WM.fPutCoef(18+iCnt, 28, d);

      WM.fPutCoef(28, 6+iCnt, -d);  
      WM.fPutCoef(28, 18+iCnt, d);

      d = Tmp1_1.dGet(iCnt);
      WM.fPutCoef(3+iCnt, 25, d);   
      WM.fPutCoef(25, 3+iCnt, d);   

      WM.fPutCoef(27, 9+iCnt, d);   
      WM.fPutCoef(9+iCnt, 27, d);
      
      d = Tmp1_2.dGet(iCnt);
      WM.fPutCoef(3+iCnt, 26, d);   
      WM.fPutCoef(26, 3+iCnt, d);   

      WM.fPutCoef(28, 9+iCnt, d);   
      WM.fPutCoef(9+iCnt, 28, d);
      
      d = Tmp2_1.dGet(iCnt);
      WM.fPutCoef(27, iCnt, -d);    
      WM.fPutCoef(27, 12+iCnt, d);  

      WM.fPutCoef(6+iCnt, 25, -d);  
      WM.fPutCoef(18+iCnt, 25, d);
      
      d = Tmp2_2.dGet(iCnt);
      WM.fPutCoef(28, iCnt, -d);    
      WM.fPutCoef(28, 12+iCnt, d);  

      WM.fPutCoef(6+iCnt, 26, -d);  
      WM.fPutCoef(18+iCnt, 26, d);

      d = Tmp3_1.dGet(iCnt);
      WM.fPutCoef(27, 3+iCnt, d);   

      d = Tmp3_2.dGet(iCnt);
      WM.fPutCoef(28, 3+iCnt, d);   

      d = Tmp4_1.dGet(iCnt);
      WM.fPutCoef(9+iCnt, 25, d);   

      d = Tmp4_2.dGet(iCnt);
      WM.fPutCoef(9+iCnt, 27, d);   
   }   

   Mat3x3 MTmp(FTmp);
   WM.Add(1, 4, MTmp);             
   WM.Add(4, 13, MTmp);            
   
   WM.Add(7, 10, MTmp);            
   WM.Add(10, 19, MTmp);           
   
   WM.Sub(4, 1, MTmp);
   WM.Sub(13, 4, MTmp);
   
   WM.Sub(19, 10, MTmp);           
   WM.Sub(10, 7, MTmp);            
   
   MTmp = Mat3x3(x2mx1, FTmp);
   WM.Add(4, 4, MTmp);             

   WM.Add(10, 10, MTmp);           
 
   MTmp = Mat3x3(Omega, FTmp)+Mat3x3(FPrimeTmp);
   WM.Add(7, 4, MTmp);             
   WM.Add(19, 4, -MTmp);           
   
   MTmp = Mat3x3(Omega.Cross(FTmp)+FPrimeTmp);
   WM.Add(10, 1, -MTmp);
   WM.Add(10, 13, MTmp);
   
   MTmp = Mat3x3((Mat3x3(xp2mxp1)+Mat3x3(x2mx1, Omega))*Mat3x3(FTmp)
		 +Mat3x3(x2mx1, FPrimeTmp));
   WM.Add(10, 4, MTmp);

   return WorkMat;
}

   
/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
InLineJoint::InitialAssRes(SubVectorHandler& WorkVec,
			   const VectorHandler& XCurr)
{
   DEBUGCOUTFNAME("InLineJoint::InitialAssRes");
   WorkVec.Resize(28);
   WorkVec.Reset(0.);
   
   /* Indici gdl */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+2;
   
   /* Indici equazioni nodi */
   for(int iCnt = 1; iCnt <= 12; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WorkVec.fPutRowIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
   }
   
   /* Indici equazioni vincoli */
   WorkVec.fPutRowIndex(25, iFirstReactionIndex+1);
   WorkVec.fPutRowIndex(26, iFirstReactionIndex+2);
   WorkVec.fPutRowIndex(27, iReactionPrimeIndex+1);
   WorkVec.fPutRowIndex(28, iReactionPrimeIndex+2);
      
   /* Dati */
   Vec3 x2mx1(pNode2->GetXCurr()-pNode1->GetXCurr());
   Vec3 xp2mxp1(pNode2->GetVCurr()-pNode1->GetVCurr());
   Mat3x3 RvTmp(pNode1->GetRCurr()*Rv);
   Vec3 Omega(pNode1->GetWCurr());
   
   /* Aggiorna i dati propri */
   F.Put(1, XCurr.dGetCoef(iFirstReactionIndex+1));
   F.Put(2, XCurr.dGetCoef(iFirstReactionIndex+2));
   Vec3 FPrime(XCurr.dGetCoef(iReactionPrimeIndex+1),
	       XCurr.dGetCoef(iReactionPrimeIndex+2),
	       0.);
   Vec3 FTmp(RvTmp*F);
   Vec3 FPrimeTmp(RvTmp*FPrime);
   Vec3 Tmp(Omega.Cross(FTmp)+FPrimeTmp);
   
   WorkVec.Add(1, FTmp);
   WorkVec.Add(4, x2mx1.Cross(FTmp)); /* ( = -p/\F) */
   WorkVec.Add(7, Tmp);
   WorkVec.Add(10, xp2mxp1.Cross(FTmp)+x2mx1.Cross(Tmp));
   WorkVec.Sub(13, FTmp);
   WorkVec.Sub(19, Tmp);   
   
   WorkVec.fPutCoef(25, Rv.GetVec(1).Dot(p)-RvTmp.GetVec(1).Dot(x2mx1));
   WorkVec.fPutCoef(26, Rv.GetVec(2).Dot(p)-RvTmp.GetVec(2).Dot(x2mx1));
   WorkVec.fPutCoef(27, x2mx1.Dot(RvTmp.GetVec(1).Cross(Omega))-RvTmp.GetVec(1).Dot(xp2mxp1));
   WorkVec.fPutCoef(28, x2mx1.Dot(RvTmp.GetVec(2).Cross(Omega))-RvTmp.GetVec(2).Dot(xp2mxp1));
   
   return WorkVec;
}

   
/* Setta il valore iniziale delle proprie variabili */
void InLineJoint::SetInitialValue(VectorHandler& /* X */ ) const
{ 
   NO_OP;
}


void InLineJoint::SetValue(VectorHandler& /* X */ , VectorHandler& /* XP */ ) const
{ 
   NO_OP;
}

/* InLIneJoint - end */


/* InLineWithOffsetJoint - begin */

/* Costruttore non banale */
InLineWithOffsetJoint::InLineWithOffsetJoint(unsigned int uL, 
					     const DofOwner* pDO,
					     const StructNode* pN1, 
					     const StructNode* pN2, 
					     const Mat3x3& RvTmp, 
					     const Vec3& pTmp,
					     const Vec3& qTmp, 
					     flag fOut)
: Elem(uL, Elem::JOINT, fOut), Joint(uL, Joint::INLINE, pDO, fOut),
pNode1(pN1), pNode2(pN2), Rv(RvTmp), p(pTmp), q(qTmp), F(0.)
{
   NO_OP;
};

   
InLineWithOffsetJoint::~InLineWithOffsetJoint(void)
{
   NO_OP;
}


/* Contributo al file di restart */
ostream& InLineWithOffsetJoint::Restart(ostream& out) const
{
   Joint::Restart(out) << ", in line, not implemented yet" << endl;
   return out;
}

   
VariableSubMatrixHandler& 
InLineWithOffsetJoint::AssJac(VariableSubMatrixHandler& WorkMat,
			      doublereal dCoef,
			      const VectorHandler& /* XCurr */ ,
			      const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUTFNAME("InLineWithOffsetJoint::AssJac");
   SparseSubMatrixHandler& WM = WorkMat.SetSparse();
   
   WM.ResizeInit(108, 0, 0.);
   
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();

   Vec3 qTmp(pNode2->GetRCurr()*q);
   Vec3 x2qmx1(pNode2->GetXCurr()+qTmp-pNode1->GetXCurr());
   Mat3x3 RvTmp(pNode1->GetRRef()*Rv);
   Vec3 FTmp(RvTmp*(F*dCoef));
   
   
   Vec3 Tmp1_1(RvTmp.GetVec(1).Cross(x2qmx1));
   Vec3 Tmp1_2(RvTmp.GetVec(2).Cross(x2qmx1));
   Vec3 Tmp2_1(qTmp.Cross(RvTmp.GetVec(1)));
   Vec3 Tmp2_2(qTmp.Cross(RvTmp.GetVec(2)));
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = RvTmp.dGet(iCnt, 1);
      WM.fPutItem(iCnt, iNode1FirstMomIndex+iCnt,
		  iFirstReactionIndex+1, -d);
      WM.fPutItem(3+iCnt, iNode2FirstMomIndex+iCnt,
		  iFirstReactionIndex+1, d);      
      
      WM.fPutItem(6+iCnt, iFirstReactionIndex+1,
		  iNode1FirstPosIndex+iCnt, -d);
      WM.fPutItem(9+iCnt, iFirstReactionIndex+1,
		  iNode2FirstPosIndex+iCnt, d);
      
      d = RvTmp.dGet(iCnt, 2);
      WM.fPutItem(12+iCnt, iNode1FirstMomIndex+iCnt,
		  iFirstReactionIndex+2, -d);
      WM.fPutItem(15+iCnt, iNode2FirstMomIndex+iCnt,
		  iFirstReactionIndex+2, d);      
      
      WM.fPutItem(18+iCnt, iFirstReactionIndex+2,
		  iNode1FirstPosIndex+iCnt, -d);
      WM.fPutItem(21+iCnt, iFirstReactionIndex+2,
		  iNode2FirstPosIndex+iCnt, d);

      d = Tmp1_1.dGet(iCnt);
      WM.fPutItem(24+iCnt, iNode1FirstMomIndex+3+iCnt,
		  iFirstReactionIndex+1, d);
      
      WM.fPutItem(27+iCnt, iFirstReactionIndex+1,
		  iNode1FirstPosIndex+3+iCnt, d);
      
      d = Tmp1_2.dGet(iCnt);
      WM.fPutItem(30+iCnt, iNode1FirstMomIndex+3+iCnt,
		  iFirstReactionIndex+2, d);
      
      WM.fPutItem(33+iCnt, iFirstReactionIndex+2,
		  iNode1FirstPosIndex+3+iCnt, d);

      d = Tmp2_1.dGet(iCnt);
      WM.fPutItem(36+iCnt, iNode2FirstMomIndex+3+iCnt,
		  iFirstReactionIndex+1, d);
      
      WM.fPutItem(39+iCnt, iFirstReactionIndex+1,
		  iNode2FirstPosIndex+3+iCnt, d);

      d = Tmp2_2.dGet(iCnt);
      WM.fPutItem(42+iCnt, iNode2FirstMomIndex+3+iCnt,
		  iFirstReactionIndex+2, d);
      
      WM.fPutItem(45+iCnt, iFirstReactionIndex+2,
		  iNode2FirstPosIndex+3+iCnt, d);
   }
   
   WM.fPutCross(48+1, iNode1FirstMomIndex,
		iNode1FirstPosIndex+3, FTmp);
   WM.fPutCross(54+1, iNode1FirstMomIndex+3,
		iNode1FirstPosIndex, -FTmp);
   WM.fPutMat3x3(60+1, iNode1FirstMomIndex+3,		
		iNode1FirstPosIndex+3, Mat3x3(x2qmx1, FTmp));
   WM.fPutCross(69+1, iNode1FirstMomIndex+3,
		iNode2FirstPosIndex, FTmp);
   WM.fPutCross(75+1, iNode2FirstMomIndex,
		iNode2FirstPosIndex+3, -FTmp);

   Mat3x3 MTmp(FTmp, qTmp);
   WM.fPutMat3x3(81+1, iNode1FirstMomIndex+3,
		 iNode2FirstPosIndex+3, -MTmp);
   WM.fPutMat3x3(90+1, iNode2FirstMomIndex+3,
		 iNode2FirstPosIndex+3, -MTmp);

   WM.fPutMat3x3(99+1, iNode2FirstMomIndex+3,		
		 iNode1FirstPosIndex+3, Mat3x3(-qTmp, FTmp));
   
   return WorkMat;
}


SubVectorHandler& 
InLineWithOffsetJoint::AssRes(SubVectorHandler& WorkVec,
			      doublereal dCoef,
			      const VectorHandler& XCurr, 
			      const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUTFNAME("InLineWithOffsetJoint::AssRes");
   WorkVec.Resize(14);
   WorkVec.Reset(0.);
 
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();

   /* Indici equazioni */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.fPutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
   }
   
   /* Indice equazione vincolo */
   WorkVec.fPutRowIndex(13, iFirstReactionIndex+1);
   WorkVec.fPutRowIndex(14, iFirstReactionIndex+2);

   Vec3 qTmp(pNode2->GetRCurr()*q);
   Vec3 x2qmx1(pNode2->GetXCurr()+qTmp-pNode1->GetXCurr());
   Mat3x3 RvTmp = pNode1->GetRCurr()*Rv;

   /* Aggiorna i dati propri */
   F.Put(1, XCurr.dGetCoef(iFirstReactionIndex+1));
   F.Put(2, XCurr.dGetCoef(iFirstReactionIndex+2));
   Vec3 FTmp(RvTmp*F);
   
   WorkVec.Add(1, FTmp);
   WorkVec.Add(4, x2qmx1.Cross(FTmp));
   WorkVec.Sub(7, FTmp);
   WorkVec.Sub(10, qTmp.Cross(FTmp));
   
   if (dCoef != 0.) {
      WorkVec.fPutCoef(13, (Rv.GetVec(1).Dot(p)-RvTmp.GetVec(1).Dot(x2qmx1))/dCoef);
      WorkVec.fPutCoef(14, (Rv.GetVec(2).Dot(p)-RvTmp.GetVec(2).Dot(x2qmx1))/dCoef);
   }

   return WorkVec;
}


void InLineWithOffsetJoint::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) {      
      Mat3x3 RvTmp(pNode1->GetRCurr()*Rv);
      Joint::Output(OH.Joints(), "inline", GetLabel(),
		    F, Zero3, RvTmp*F, Zero3) << endl;
   }
}
 

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
InLineWithOffsetJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
				     const VectorHandler& XCurr)
{
   DEBUGCOUTFNAME("InLineWithOffsetJoint::InitialAssJac");
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.ResizeInit(28, 28, 0.);
   
   /* Indici gdl */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+2;
   
   /* Indici equazioni nodi */
   for(int iCnt = 1; iCnt <= 12; iCnt++) {
      WM.fPutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.fPutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.fPutRowIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
      WM.fPutColIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
   }
   
   /* Indici vincoli */
   for(int iCnt = 1; iCnt <= 4; iCnt++) {
      WM.fPutRowIndex(24+iCnt, iFirstReactionIndex+iCnt);
      WM.fPutColIndex(24+iCnt, iFirstReactionIndex+iCnt);
   }
   
   /* Dati */
   Vec3 Omega1(pNode1->GetWRef());
   Vec3 Omega2(pNode2->GetWRef());
   Mat3x3 RvTmp(pNode1->GetRRef()*Rv);
   Vec3 qTmp(pNode2->GetRRef()*q);
   Vec3 x2qmx1(pNode2->GetXCurr()+qTmp-pNode1->GetXCurr());
   Vec3 xp2qmxp1(pNode2->GetVCurr()+Omega2.Cross(qTmp)-pNode1->GetVCurr());
   
   /* Aggiorna i dati propri */
   Vec3 FPrime(XCurr.dGetCoef(iReactionPrimeIndex+1),
	       XCurr.dGetCoef(iReactionPrimeIndex+2),
	       0.);
   Vec3 FTmp(RvTmp*F);
   Vec3 FPrimeTmp(RvTmp*FPrime);

   Vec3 v1(RvTmp.GetVec(1));
   Vec3 Tmp1_1(v1.Cross(x2qmx1));
   Vec3 Tmp2_1(Omega1.Cross(v1));
   Vec3 Tmp3_1((Omega1.Cross(x2qmx1)-xp2qmxp1).Cross(v1));
   Vec3 Tmp4_1(-(xp2qmxp1.Cross(v1)+x2qmx1.Cross(Tmp2_1)));
   
   Vec3 Tmp5_1(qTmp.Cross(v1));
   Vec3 Tmp6_1(qTmp.Cross(v1.Cross(Omega2-Omega1)));
   Vec3 Tmp7_1(qTmp.Cross(Omega1.Cross(v1))-v1.Cross(Omega2.Cross(qTmp)));
   
   Vec3 v2(RvTmp.GetVec(2));
   Vec3 Tmp1_2(v2.Cross(x2qmx1));
   Vec3 Tmp2_2(Omega1.Cross(v2));
   Vec3 Tmp3_2((Omega1.Cross(x2qmx1)-xp2qmxp1).Cross(v2));
   Vec3 Tmp4_2(-(xp2qmxp1.Cross(v2)+x2qmx1.Cross(Tmp2_2)));
   
   Vec3 Tmp5_2(qTmp.Cross(v2));
   Vec3 Tmp6_2(qTmp.Cross(v2.Cross(Omega2-Omega1)));
   Vec3 Tmp7_2(qTmp.Cross(Omega1.Cross(v2))-v2.Cross(Omega2.Cross(qTmp)));
   
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = v1.dGet(iCnt);
      WM.fPutCoef(iCnt, 25, -d);
      WM.fPutCoef(12+iCnt, 25, d);
      
      WM.fPutCoef(25, iCnt, -d);
      WM.fPutCoef(25, 12+iCnt, d);

      WM.fPutCoef(6+iCnt, 27, -d);  
      WM.fPutCoef(18+iCnt, 27, d);

      WM.fPutCoef(27, 6+iCnt, -d);  
      WM.fPutCoef(27, 18+iCnt, d);
      
      d = v2.dGet(iCnt);
      WM.fPutCoef(iCnt, 26, -d);
      WM.fPutCoef(12+iCnt, 26, d);
      
      WM.fPutCoef(26, iCnt, -d);
      WM.fPutCoef(26, 12+iCnt, d);

      WM.fPutCoef(6+iCnt, 28, -d);  
      WM.fPutCoef(18+iCnt, 28, d);

      WM.fPutCoef(28, 6+iCnt, -d);  
      WM.fPutCoef(28, 18+iCnt, d);

      d = Tmp1_1.dGet(iCnt);
      WM.fPutCoef(3+iCnt, 25, d);   
      WM.fPutCoef(25, 3+iCnt, d);   

      WM.fPutCoef(27, 9+iCnt, d);   
      WM.fPutCoef(9+iCnt, 27, d);
      
      d = Tmp1_2.dGet(iCnt);
      WM.fPutCoef(3+iCnt, 26, d);   
      WM.fPutCoef(26, 3+iCnt, d);   

      WM.fPutCoef(28, 9+iCnt, d);   
      WM.fPutCoef(9+iCnt, 28, d);
      
      d = Tmp2_1.dGet(iCnt);
      WM.fPutCoef(27, iCnt, -d);    
      WM.fPutCoef(27, 12+iCnt, d);  

      WM.fPutCoef(6+iCnt, 25, -d);  
      WM.fPutCoef(18+iCnt, 25, d);
      
      d = Tmp2_2.dGet(iCnt);
      WM.fPutCoef(28, iCnt, -d);    
      WM.fPutCoef(28, 12+iCnt, d);  

      WM.fPutCoef(6+iCnt, 26, -d);  
      WM.fPutCoef(18+iCnt, 26, d);

      d = Tmp3_1.dGet(iCnt);
      WM.fPutCoef(27, 3+iCnt, d);   

      d = Tmp3_2.dGet(iCnt);
      WM.fPutCoef(28, 3+iCnt, d);   

      d = Tmp4_1.dGet(iCnt);
      WM.fPutCoef(9+iCnt, 25, d);   

      d = Tmp4_2.dGet(iCnt);
      WM.fPutCoef(9+iCnt, 27, d);   
      
      d = Tmp5_1.dGet(iCnt);
      WM.fPutCoef(15+iCnt, 25, d);
      WM.fPutCoef(21+iCnt, 27, d);

      WM.fPutCoef(25, 15+iCnt, d);
      WM.fPutCoef(27, 21+iCnt, d);

      d = Tmp5_2.dGet(iCnt);
      WM.fPutCoef(15+iCnt, 26, d);
      WM.fPutCoef(21+iCnt, 28, d);

      WM.fPutCoef(26, 15+iCnt, d);
      WM.fPutCoef(28, 21+iCnt, d);

      d = Tmp6_1.dGet(iCnt);
      WM.fPutCoef(27, 15+iCnt, d);

      d = Tmp6_2.dGet(iCnt);
      WM.fPutCoef(28, 15+iCnt, d);

      d = Tmp7_1.dGet(iCnt);
      WM.fPutCoef(21+iCnt, 25, d);
      
      d = Tmp7_2.dGet(iCnt);
      WM.fPutCoef(21+iCnt, 26, d);
   }   

   Mat3x3 MTmp(FTmp);
   WM.Add(1, 4, MTmp);             
   WM.Add(4, 13, MTmp);            
   
   WM.Add(7, 10, MTmp);            
   WM.Add(10, 19, MTmp);           
   
   WM.Sub(4, 1, MTmp);             
   WM.Sub(13, 4, MTmp);            
   
   WM.Sub(19, 10, MTmp);           
   WM.Sub(10, 7, MTmp);            
   
   MTmp = Mat3x3(x2qmx1, FTmp);
   WM.Add(4, 4, MTmp);             

   WM.Add(10, 10, MTmp);           
 
   MTmp = Mat3x3(Omega1, FTmp)+Mat3x3(FPrimeTmp);
   WM.Add(7, 4, MTmp);             
   WM.Add(19, 4, -MTmp);           
   
   MTmp = Mat3x3(Omega1.Cross(FTmp)+FPrimeTmp);
   WM.Add(10, 1, -MTmp);
   WM.Add(10, 13, MTmp);
   
   MTmp = Mat3x3((Mat3x3(xp2qmxp1)+Mat3x3(x2qmx1, Omega1))*Mat3x3(FTmp)
		 +Mat3x3(x2qmx1, FPrimeTmp));
   WM.Add(10, 4, MTmp);

   MTmp = Mat3x3(FTmp, qTmp);
   WM.Add(16, 16, MTmp);
   WM.Add(22, 22, MTmp);
   
   WM.Sub(4, 16, MTmp);
   WM.Sub(10, 22, MTmp);
   
   MTmp = MTmp.Transpose();
   WM.Add(16, 4, MTmp);
   WM.Add(22, 10, MTmp);
   
   MTmp = (Mat3x3(FTmp, Omega2)+Mat3x3(Omega1.Cross(FTmp)+FPrimeTmp))*Mat3x3(qTmp);
   WM.Add(22, 16, MTmp);
   WM.Add(10, 16, -MTmp);
   
   MTmp = Mat3x3(qTmp.Cross(Omega2), FTmp)
     -Mat3x3(qTmp)*(Mat3x3(Omega1, FTmp)+Mat3x3(FPrimeTmp));
   WM.Add(22, 4, MTmp);

   return WorkMat;
}

   
/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
InLineWithOffsetJoint::InitialAssRes(SubVectorHandler& WorkVec,
				     const VectorHandler& XCurr)
{
   DEBUGCOUTFNAME("InLineWithOffsetJoint::InitialAssRes");
   WorkVec.Resize(28);
   WorkVec.Reset(0.);
   
   /* Indici gdl */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+2;
   
   /* Indici equazioni nodi */
   for(int iCnt = 1; iCnt <= 12; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WorkVec.fPutRowIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
   }
   
   /* Indici equazioni vincoli */
   WorkVec.fPutRowIndex(25, iFirstReactionIndex+1);
   WorkVec.fPutRowIndex(26, iFirstReactionIndex+2);
   WorkVec.fPutRowIndex(27, iReactionPrimeIndex+1);
   WorkVec.fPutRowIndex(28, iReactionPrimeIndex+2);
      
   /* Dati */
   Mat3x3 RvTmp(pNode1->GetRCurr()*Rv);
   Vec3 qTmp(pNode2->GetRCurr()*q);
   Vec3 Omega1(pNode1->GetWCurr());
   Vec3 Omega2(pNode2->GetWCurr());
   Vec3 x2qmx1(pNode2->GetXCurr()+qTmp-pNode1->GetXCurr());
   Vec3 xp2qmxp1(pNode2->GetVCurr()+Omega2.Cross(qTmp)-pNode1->GetVCurr());
   
   /* Aggiorna i dati propri */
   F.Put(1, XCurr.dGetCoef(iFirstReactionIndex+1));
   F.Put(2, XCurr.dGetCoef(iFirstReactionIndex+2));
   Vec3 FPrime(XCurr.dGetCoef(iReactionPrimeIndex+1),
	       XCurr.dGetCoef(iReactionPrimeIndex+2),
	       0.);
   Vec3 FTmp(RvTmp*F);
   Vec3 FPrimeTmp(RvTmp*FPrime);
   Vec3 Tmp(Omega1.Cross(FTmp)+FPrimeTmp);
   
   WorkVec.Add(1, FTmp);
   WorkVec.Add(4, x2qmx1.Cross(FTmp)); /* ( = -p/\F) */
   WorkVec.Add(7, Tmp);
   WorkVec.Add(10, xp2qmxp1.Cross(FTmp)+x2qmx1.Cross(Tmp));
   WorkVec.Sub(13, FTmp);
   WorkVec.Sub(16, qTmp.Cross(FTmp));
   WorkVec.Sub(19, Tmp);
   WorkVec.Sub(22, qTmp.Cross(Tmp)+(Omega2.Cross(qTmp)).Cross(FTmp));
   
   WorkVec.fPutCoef(25, Rv.GetVec(1).Dot(p)-RvTmp.GetVec(1).Dot(x2qmx1));
   WorkVec.fPutCoef(26, Rv.GetVec(2).Dot(p)-RvTmp.GetVec(2).Dot(x2qmx1));
   WorkVec.fPutCoef(27, x2qmx1.Dot(RvTmp.GetVec(1).Cross(Omega1))-RvTmp.GetVec(1).Dot(xp2qmxp1));
   WorkVec.fPutCoef(28, x2qmx1.Dot(RvTmp.GetVec(2).Cross(Omega1))-RvTmp.GetVec(2).Dot(xp2qmxp1));
   
   return WorkVec;
}

   
/* Setta il valore iniziale delle proprie variabili */
void InLineWithOffsetJoint::SetInitialValue(VectorHandler& /* X */ ) const
{ 
   NO_OP;
}


void InLineWithOffsetJoint::SetValue(VectorHandler& /* X */ , VectorHandler& /* XP */ ) const
{ 
   NO_OP;
}

/* InLineWithOffsetJoint - end */
