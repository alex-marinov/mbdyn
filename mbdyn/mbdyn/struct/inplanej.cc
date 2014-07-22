/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

/* Vincoli relativi a movimenti lineari */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "inplanej.h"

/* InPlaneJoint - begin */

/* Costruttore non banale */
InPlaneJoint::InPlaneJoint(unsigned int uL, const DofOwner* pDO,
			   const StructNode* pN1, const StructNode* pN2, 
			   const Vec3& vTmp, const Vec3& pTmp, flag fOut)
: Elem(uL, fOut), Joint(uL, pDO, fOut),
pNode1(pN1), pNode2(pN2), v(vTmp), p(pTmp), dF(0.)
{
   NO_OP;
};

   
InPlaneJoint::~InPlaneJoint(void)
{
   NO_OP;
}

/* Contributo al file di restart */
std::ostream& InPlaneJoint::Restart(std::ostream& out) const
{
   Joint::Restart(out) << ", in plane, "
     << pNode1->GetLabel() 
     << ", reference, node, ",
     p.Write(out, ", ") 
     << ", reference, node, ";
   v.Write(out, ", ") << ", "
     << pNode2->GetLabel() << ';' << std::endl;
   return out;
}

   
VariableSubMatrixHandler& 
InPlaneJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		     doublereal dCoef,
		     const VectorHandler& /* XCurr */ ,
		     const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering InPlaneJoint::AssJac()" << std::endl);
   SparseSubMatrixHandler& WM = WorkMat.SetSparse();
   WM.ResizeReset(51, 0);
   
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   
   Vec3 x2mx1(pNode2->GetXCurr()-pNode1->GetXCurr());
   Vec3 vTmp(pNode1->GetRRef()*v);
   Vec3 F(vTmp*(dF*dCoef));
   
   
   Vec3 Tmp(vTmp.Cross(x2mx1));
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = vTmp.dGet(iCnt);
      WM.PutItem(iCnt, iNode1FirstMomIndex+iCnt,
		  iFirstReactionIndex+1, -d);
      WM.PutItem(3+iCnt, iNode2FirstMomIndex+iCnt,
		  iFirstReactionIndex+1, d);
      
      WM.PutItem(6+iCnt, iFirstReactionIndex+1,
		  iNode1FirstPosIndex+iCnt, -d);
      WM.PutItem(9+iCnt, iFirstReactionIndex+1,
		  iNode2FirstPosIndex+iCnt, d);

      d = Tmp.dGet(iCnt);
      WM.PutItem(12+iCnt, iNode1FirstMomIndex+3+iCnt,
		  iFirstReactionIndex+1, d);
      
      WM.PutItem(15+iCnt, iFirstReactionIndex+1,
		  iNode1FirstPosIndex+3+iCnt, d);
   }
   
   WM.PutCross(19, iNode1FirstMomIndex,
		iNode1FirstPosIndex+3, F);
   WM.PutCross(25, iNode1FirstMomIndex+3,
		iNode1FirstPosIndex, -F);
   WM.PutMat3x3(31, iNode1FirstMomIndex+3,		
		iNode1FirstPosIndex+3, Mat3x3(MatCrossCross, x2mx1, F));
   WM.PutCross(40, iNode1FirstMomIndex+3,
		iNode2FirstPosIndex, F);
   WM.PutCross(46, iNode2FirstMomIndex,
		iNode2FirstPosIndex+3, -F);
   
   return WorkMat;
}


SubVectorHandler& InPlaneJoint::AssRes(SubVectorHandler& WorkVec,
				       doublereal dCoef,
				       const VectorHandler& XCurr, 
				       const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering InPlaneJoint::AssRes()" << std::endl);
   WorkVec.ResizeReset(13);
 
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();

   /* Indici equazioni */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.PutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
   }
   
   /* Indice equazione vincolo */
   WorkVec.PutRowIndex(13, iFirstReactionIndex+1);   

   Vec3 x2mx1(pNode2->GetXCurr()-pNode1->GetXCurr());
   Mat3x3 R(pNode1->GetRCurr());
   Vec3 vTmp = R*v;
   // Vec3 pTmp = R*p;
   
   /* Aggiorna i dati propri */
   dF = XCurr(iFirstReactionIndex+1);
   Vec3 F(vTmp*dF);
   
   WorkVec.Add(1, F);
   WorkVec.Add(4, x2mx1.Cross(F)); /* ( = -p/\F) */
   WorkVec.Add(7, -F);
   ASSERT(dCoef != 0.);
   WorkVec.PutCoef(13, (v.Dot(p)-vTmp.Dot(x2mx1))/dCoef);

   return WorkVec;
}


void InPlaneJoint::Output(OutputHandler& OH) const
{
   if(fToBeOutput()) {      
      Vec3 vTmp(pNode1->GetRCurr()*v);
      Joint::Output(OH.Joints(), "InPlane", GetLabel(),
		    Vec3(dF, 0., 0.), Zero3, vTmp*dF, Zero3) << std::endl;      
   }   
}
 

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
InPlaneJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
			    const VectorHandler& XCurr)
{
   DEBUGCOUT("Entering InPlaneJoint::InitialAssJac()" << std::endl);
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.ResizeReset(26, 26);
   
   /* Indici gdl */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+1;
   
   /* Indici equazioni nodi */
   for(int iCnt = 1; iCnt <= 12; iCnt++) {
      WM.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutRowIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
      WM.PutColIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
   }
   
   /* Indici vincoli */
   for(int iCnt = 1; iCnt <= 2; iCnt++) {
      WM.PutRowIndex(24+iCnt, iFirstReactionIndex+iCnt);
      WM.PutColIndex(24+iCnt, iFirstReactionIndex+iCnt);
   }
   
   /* Dati */
   Vec3 x2mx1(pNode2->GetXCurr()-pNode1->GetXCurr());
   Vec3 xp2mxp1(pNode2->GetVCurr()-pNode1->GetVCurr());
   Vec3 Omega(pNode1->GetWRef());
   
   /* Aggiorna i dati propri */
   doublereal dFPrime = XCurr(iReactionPrimeIndex+1);
   Vec3 vTmp(pNode1->GetRRef()*v);
   Vec3 F(vTmp*dF);
   Vec3 FPrime(vTmp*dFPrime);

   Vec3 Tmp1(vTmp.Cross(x2mx1));
   Vec3 Tmp2(Omega.Cross(vTmp));
   Vec3 Tmp3((Omega.Cross(x2mx1)-xp2mxp1).Cross(vTmp));
   Vec3 Tmp4(-(xp2mxp1.Cross(vTmp)+x2mx1.Cross(Tmp2)));
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = vTmp.dGet(iCnt);
      WM.PutCoef(iCnt, 25, -d);    
      WM.PutCoef(12+iCnt, 25, d);  
      
      WM.PutCoef(25, iCnt, -d);    
      WM.PutCoef(25, 12+iCnt, d);  

      WM.PutCoef(6+iCnt, 26, -d);  
      WM.PutCoef(18+iCnt, 26, d);  

      WM.PutCoef(26, 6+iCnt, -d);  
      WM.PutCoef(26, 18+iCnt, d);  

      d = Tmp1.dGet(iCnt);
      WM.PutCoef(3+iCnt, 25, d);   
      WM.PutCoef(25, 3+iCnt, d);   

      WM.PutCoef(26, 9+iCnt, d);   
      WM.PutCoef(9+iCnt, 26, d);   
      
      d = Tmp2.dGet(iCnt);
      WM.PutCoef(26, iCnt, -d);    
      WM.PutCoef(26, 12+iCnt, d);  

      WM.PutCoef(6+iCnt, 25, -d);  
      WM.PutCoef(18+iCnt, 25, d);  

      d = Tmp3.dGet(iCnt);
      WM.PutCoef(26, 3+iCnt, d);   

      d = Tmp4.dGet(iCnt);
      WM.PutCoef(9+iCnt, 25, d);   
   }   

   Mat3x3 MTmp(MatCross, F);
   WM.Add(1, 4, MTmp);             
   WM.Add(4, 13, MTmp);            
   
   WM.Add(7, 10, MTmp);            
   WM.Add(10, 19, MTmp);           
   
   MTmp -= MTmp;
   WM.Add(4, 1, MTmp);             
   WM.Add(13, 4, MTmp);            
   
   WM.Add(19, 10, MTmp);           
   WM.Add(10, 7, MTmp);            
   
   MTmp = Mat3x3(MatCrossCross, x2mx1, F);
   WM.Add(4, 4, MTmp);             

   WM.Add(10, 10, MTmp);           
 
   MTmp = Mat3x3(MatCrossCross, Omega, F) + Mat3x3(MatCross, FPrime);
   WM.Add(7, 4, MTmp);             
   WM.Sub(19, 4, MTmp);
   
   MTmp = Mat3x3(MatCross, Omega.Cross(F) + FPrime);
   WM.Sub(10, 1, MTmp);           
   WM.Add(10, 13, MTmp);           
   
   MTmp = Mat3x3((Mat3x3(MatCross, xp2mxp1) + Mat3x3(MatCrossCross, x2mx1, Omega))*Mat3x3(MatCross, F)
		 + Mat3x3(MatCrossCross, x2mx1, FPrime));
   WM.Add(10, 4, MTmp); 
      
   return WorkMat;
}

   
/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& InPlaneJoint::InitialAssRes(SubVectorHandler& WorkVec,
					      const VectorHandler& XCurr)
{
   DEBUGCOUT("Entering InPlaneJoint::InitialAssRes()" << std::endl);
   WorkVec.ResizeReset(26);
   
   /* Indici gdl */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+1;
   
   /* Indici equazioni nodi */
   for(int iCnt = 1; iCnt <= 12; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WorkVec.PutRowIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
   }
   
   /* Indici equazioni vincoli */
   WorkVec.PutRowIndex(25, iFirstReactionIndex+1);
   WorkVec.PutRowIndex(26, iReactionPrimeIndex+1);
      
   /* Dati */
   Vec3 x2mx1(pNode2->GetXCurr()-pNode1->GetXCurr());
   Vec3 xp2mxp1(pNode2->GetVCurr()-pNode1->GetVCurr());
   Mat3x3 R(pNode1->GetRCurr());
   Vec3 Omega(pNode1->GetWCurr());
   Vec3 vTmp = R*v;
   // Vec3 pTmp = R*p;   
   
   /* Aggiorna i dati propri */
   dF = XCurr(iFirstReactionIndex+1);
   doublereal dFPrime = XCurr(iReactionPrimeIndex+1);
   Vec3 F(vTmp*dF);
   Vec3 FPrime(vTmp*dFPrime);
   Vec3 Tmp(Omega.Cross(F)+FPrime);
   
   WorkVec.Add(1, F);
   WorkVec.Add(4, x2mx1.Cross(F)); /* ( = -p/\F) */
   WorkVec.Add(7, Tmp);
   WorkVec.Add(10, xp2mxp1.Cross(F)+x2mx1.Cross(Tmp));
   WorkVec.Add(13, -F);
   WorkVec.Add(19, -Tmp);
   
   WorkVec.PutCoef(25, v.Dot(p)-vTmp.Dot(x2mx1));
   WorkVec.PutCoef(26, x2mx1.Dot(vTmp.Cross(Omega))-vTmp.Dot(xp2mxp1));
      
   return WorkVec;
}

   
/* Setta il valore iniziale delle proprie variabili */
void InPlaneJoint::SetInitialValue(VectorHandler& /* X */ )
{ 
   NO_OP;
}

/* InPlaneJoint - end */


/* InPlaneWithOffsetJoint - begin */

/* Costruttore non banale */
InPlaneWithOffsetJoint::InPlaneWithOffsetJoint(unsigned int uL, 
					       const DofOwner* pDO,
					       const StructNode* pN1, 
					       const StructNode* pN2, 
					       const Vec3& vT, 
					       const Vec3& pT, 
					       const Vec3& qT,
					       flag fOut)
: Elem(uL, fOut), Joint(uL, pDO, fOut),
pNode1(pN1), pNode2(pN2), v(vT), p(pT), q(qT), dF(0.)
{
   NO_OP;
};

   
InPlaneWithOffsetJoint::~InPlaneWithOffsetJoint(void)
{
   NO_OP;
}

/* Contributo al file di restart */
std::ostream& InPlaneWithOffsetJoint::Restart(std::ostream& out) const
{
   Joint::Restart(out) << ", in plane, "
     << pNode1->GetLabel() 
     << ", reference, node, ",
     p.Write(out, ", ") 
     << ", reference, node, ",
     v.Write(out, ", ") << ", "
     << pNode2->GetLabel() << ", offset, reference, node, ";
   return q.Write(out, ", ") << ';' << std::endl;
}

   
VariableSubMatrixHandler& 
InPlaneWithOffsetJoint::AssJac(VariableSubMatrixHandler& WorkMat,
			       doublereal dCoef,
			       const VectorHandler& /* XCurr */ ,
			       const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering InPlaneWithOffsetJoint::AssJac()" << std::endl);
   SparseSubMatrixHandler& WM = WorkMat.SetSparse();
   WM.ResizeReset(84, 0);
   
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   
   Vec3 vTmp(pNode1->GetRRef()*v);
   Vec3 qTmp(pNode2->GetRRef()*q);
   Vec3 x2pqmx1(pNode2->GetXCurr()+qTmp-pNode1->GetXCurr());
   Vec3 F(vTmp*(dF*dCoef));
   
   
   Vec3 Tmp1(vTmp.Cross(x2pqmx1));
   Vec3 Tmp2(qTmp.Cross(vTmp));
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = vTmp.dGet(iCnt);
      WM.PutItem(iCnt, iNode1FirstMomIndex+iCnt,
		  iFirstReactionIndex+1, -d);
      WM.PutItem(3+iCnt, iNode2FirstMomIndex+iCnt,
		  iFirstReactionIndex+1, d);
      
      WM.PutItem(6+iCnt, iFirstReactionIndex+1,
		  iNode1FirstPosIndex+iCnt, -d);
      WM.PutItem(9+iCnt, iFirstReactionIndex+1,
		  iNode2FirstPosIndex+iCnt, d);

      d = Tmp1.dGet(iCnt);
      WM.PutItem(12+iCnt, iNode1FirstMomIndex+3+iCnt,
		  iFirstReactionIndex+1, d);
      
      WM.PutItem(15+iCnt, iFirstReactionIndex+1,
		  iNode1FirstPosIndex+3+iCnt, d);

      d = Tmp2.dGet(iCnt);
      WM.PutItem(18+iCnt, iNode2FirstMomIndex+3+iCnt,
		  iFirstReactionIndex+1, d);
      
      WM.PutItem(21+iCnt, iFirstReactionIndex+1,
		  iNode2FirstPosIndex+3+iCnt, d);
   }
   
   WM.PutCross(25, iNode1FirstMomIndex,
		iNode1FirstPosIndex+3, F);
   WM.PutCross(31, iNode1FirstMomIndex+3,
		iNode1FirstPosIndex, -F);
   WM.PutMat3x3(37, iNode1FirstMomIndex+3,		
		 iNode1FirstPosIndex+3, Mat3x3(MatCrossCross, x2pqmx1, F));
   WM.PutCross(46, iNode1FirstMomIndex+3,
		iNode2FirstPosIndex, F);
   WM.PutCross(52, iNode2FirstMomIndex,
		iNode2FirstPosIndex+3, -F);
   
   Mat3x3 MTmp(MatCrossCross, F, qTmp);
   WM.PutMat3x3(58, iNode1FirstMomIndex+3,		
		 iNode2FirstPosIndex+3, -MTmp);
   WM.PutMat3x3(67, iNode2FirstMomIndex+3,		
		 iNode2FirstPosIndex+3, -MTmp);

   WM.PutMat3x3(76, iNode2FirstMomIndex+3,		
		 iNode1FirstPosIndex+3, Mat3x3(MatCrossCross, -qTmp, F));
      
   return WorkMat;
}


SubVectorHandler& 
InPlaneWithOffsetJoint::AssRes(SubVectorHandler& WorkVec,
			       doublereal dCoef,
			       const VectorHandler& XCurr,
			       const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering InPlaneWithOffsetJoint::AssRes()" << std::endl);
   WorkVec.ResizeReset(13);
 
   // integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   // integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();

   /* Indici equazioni */
   for(int iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.PutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
   }
   
   /* Indice equazione vincolo */
   WorkVec.PutRowIndex(13, iFirstReactionIndex+1);   
  
   Vec3 vTmp(pNode1->GetRCurr()*v);  
   Vec3 qTmp(pNode2->GetRCurr()*q);
   Vec3 x2pqmx1(pNode2->GetXCurr()+qTmp-pNode1->GetXCurr());
   
   /* Aggiorna i dati propri */
   dF = XCurr(iFirstReactionIndex+1);
   Vec3 F(vTmp*dF);
   
   WorkVec.Add(1, F);
   WorkVec.Add(4, x2pqmx1.Cross(F));
   WorkVec.Add(7, -F);
   WorkVec.Add(10, F.Cross(qTmp));
   if(dCoef != 0.) {
      WorkVec.PutCoef(13, (v.Dot(p)-vTmp.Dot(x2pqmx1))/dCoef);
   }
      
   return WorkVec;
}


void InPlaneWithOffsetJoint::Output(OutputHandler& OH) const
{
   if(fToBeOutput()) {            
      Vec3 vTmp(pNode1->GetRCurr()*v);
      Joint::Output(OH.Joints(), "InPlaneWithOffs", GetLabel(),
		    Vec3(dF, 0., 0.), Zero3, vTmp*dF, Zero3) << std::endl;      
   }   
}
 

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
InPlaneWithOffsetJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
				      const VectorHandler& XCurr)
{
   DEBUGCOUT("Entering InPlaneWithOffsetJoint::InitialAssJac()" << std::endl);
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.ResizeReset(26, 26);
   
   /* Indici gdl */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+1;
   
   /* Indici equazioni nodi */
   for(int iCnt = 1; iCnt <= 12; iCnt++) {
      WM.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutRowIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
      WM.PutColIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
   }
   
   /* Indici vincoli */
   for(int iCnt = 1; iCnt <= 2; iCnt++) {
      WM.PutRowIndex(24+iCnt, iFirstReactionIndex+iCnt);
      WM.PutColIndex(24+iCnt, iFirstReactionIndex+iCnt);
   }
   
   /* Dati */
   Vec3 Omega1(pNode1->GetWRef());
   Vec3 Omega2(pNode2->GetWRef());
   Vec3 vTmp(pNode1->GetRRef()*v);
   Vec3 qTmp(pNode2->GetRRef()*q);
   Vec3 x2pqmx1(pNode2->GetXCurr()+qTmp-pNode1->GetXCurr());
   Vec3 xp2pqpmxp1(pNode2->GetVCurr()+Omega2.Cross(qTmp)-pNode1->GetVCurr());
   
   /* Aggiorna i dati propri */
   doublereal dFPrime = XCurr(iReactionPrimeIndex+1);
   Vec3 F(vTmp*dF);
   Vec3 FPrime(vTmp*dFPrime);

   Vec3 Tmp1(vTmp.Cross(x2pqmx1));
   Vec3 Tmp2(Omega1.Cross(vTmp));
   Vec3 Tmp3((Omega1.Cross(x2pqmx1) - xp2pqpmxp1).Cross(vTmp));
   Vec3 Tmp4(-(xp2pqpmxp1.Cross(vTmp) + x2pqmx1.Cross(Tmp2)));
   
   Vec3 Tmp5(qTmp.Cross(vTmp));
   Vec3 Tmp6(qTmp.Cross(vTmp.Cross(Omega2-Omega1)));
   Vec3 Tmp7(qTmp.Cross(Omega1.Cross(vTmp)) - vTmp.Cross(Omega2.Cross(qTmp)));
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = vTmp.dGet(iCnt);
      WM.PutCoef(iCnt, 25, -d);    
      WM.PutCoef(12+iCnt, 25, d);  
      
      WM.PutCoef(25, iCnt, -d);   
      WM.PutCoef(25, 12+iCnt, d);  

      WM.PutCoef(6+iCnt, 26, -d); 
      WM.PutCoef(18+iCnt, 26, d);  

      WM.PutCoef(26, 6+iCnt, -d);  
      WM.PutCoef(26, 18+iCnt, d);  

      d = Tmp1.dGet(iCnt);
      WM.PutCoef(3+iCnt, 25, d);   
      WM.PutCoef(25, 3+iCnt, d);   

      WM.PutCoef(26, 9+iCnt, d);   
      WM.PutCoef(9+iCnt, 26, d);   
	
      d = Tmp2.dGet(iCnt);
      WM.PutCoef(26, iCnt, -d);    
      WM.PutCoef(26, 12+iCnt, d);  
      
      WM.PutCoef(6+iCnt, 25, -d);  
      WM.PutCoef(18+iCnt, 25, d);  
      
      d = Tmp3.dGet(iCnt);
      WM.PutCoef(26, 3+iCnt, d);   
      
      d = Tmp4.dGet(iCnt);
      WM.PutCoef(9+iCnt, 25, d);   

      d = Tmp5.dGet(iCnt);
      WM.PutCoef(15+iCnt, 25, d);   
      WM.PutCoef(21+iCnt, 26, d);   

      WM.PutCoef(25, 15+iCnt, d);   
      WM.PutCoef(26, 21+iCnt, d);
      
      d = Tmp6.dGet(iCnt);
      WM.PutCoef(26, 15+iCnt, d);   

      d = Tmp7.dGet(iCnt);
      WM.PutCoef(21+iCnt, 25, d);   
   }   

   Mat3x3 MTmp(MatCross, F);
   WM.Add(1, 4, MTmp);             
   WM.Add(4, 13, MTmp);            
   
   WM.Add(7, 10, MTmp);            
   WM.Add(10, 19, MTmp);           
   
   WM.Sub(4, 1, MTmp);
   WM.Sub(13, 4, MTmp);
   
   WM.Sub(19, 10, MTmp);           
   WM.Sub(10, 7, MTmp);            
   
   MTmp = Mat3x3(MatCrossCross, x2pqmx1, F);
   WM.Add(4, 4, MTmp);             

   WM.Add(10, 10, MTmp);
 
   MTmp = Mat3x3(MatCrossCross, Omega1, F) + Mat3x3(MatCross, FPrime);
   WM.Add(7, 4, MTmp);   
   WM.Sub(19, 4, MTmp);
      
   MTmp = Mat3x3(MatCross, Omega1.Cross(F) + FPrime);
   WM.Sub(10, 1, MTmp);
   WM.Add(10, 13, MTmp);
   
   MTmp = (Mat3x3(MatCross, xp2pqpmxp1) + Mat3x3(MatCrossCross, x2pqmx1, Omega1))*Mat3x3(MatCross, F)
     + Mat3x3(MatCrossCross, x2pqmx1, FPrime);
   WM.Add(10, 4, MTmp);
   
   MTmp = Mat3x3(MatCrossCross, F, qTmp);
   WM.Add(16, 16, MTmp);
   WM.Add(22, 22, MTmp);
   
   WM.Sub(4, 16, MTmp);
   WM.Sub(10, 22, MTmp);
   
   MTmp = MTmp.Transpose();
   WM.Add(16, 4, MTmp);
   WM.Add(22, 10, MTmp);
   
   MTmp = (Mat3x3(MatCrossCross, F, Omega2) + Mat3x3(MatCross, Omega1.Cross(F) + FPrime))*Mat3x3(MatCross, qTmp);
   WM.Add(22, 16, MTmp);
   WM.Sub(10, 16, MTmp);
   
   MTmp = Mat3x3(MatCrossCross, qTmp.Cross(Omega2), F)
     - qTmp.Cross(Mat3x3(MatCrossCross, Omega1, F) + Mat3x3(MatCross, FPrime));
   WM.Add(22, 4, MTmp);
   
   return WorkMat;
}

   
/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
InPlaneWithOffsetJoint::InitialAssRes(SubVectorHandler& WorkVec,
				      const VectorHandler& XCurr)
{
   DEBUGCOUT("Entering InPlaneWithOffsetJoint::InitialAssRes()" << std::endl);
   WorkVec.ResizeReset(26);
   
   /* Indici gdl */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+1;
   
   /* Indici equazioni nodi */
   for(int iCnt = 1; iCnt <= 12; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WorkVec.PutRowIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
   }
   
   /* Indici equazioni vincoli */
   WorkVec.PutRowIndex(25, iFirstReactionIndex+1);
   WorkVec.PutRowIndex(26, iReactionPrimeIndex+1);
      
   /* Dati */   
   Vec3 Omega1(pNode1->GetWCurr());
   Vec3 Omega2(pNode2->GetWCurr());
   Vec3 vTmp(pNode1->GetRCurr()*v);
   Vec3 qTmp(pNode2->GetRCurr()*q);

   Vec3 x2pqmx1(pNode2->GetXCurr()+qTmp-pNode1->GetXCurr());
   Vec3 xp2pqpmxp1(pNode2->GetVCurr()+Omega2.Cross(qTmp)-pNode1->GetVCurr());
   
   /* Aggiorna i dati propri */
   dF = XCurr(iFirstReactionIndex+1);
   doublereal dFPrime = XCurr(iReactionPrimeIndex+1);
   Vec3 F(vTmp*dF);
   Vec3 FPrime(vTmp*dFPrime);
   Vec3 Tmp(Omega1.Cross(F)+FPrime);
   
   WorkVec.Add(1, F);
   WorkVec.Add(4, x2pqmx1.Cross(F));
   WorkVec.Add(7, Tmp);
   WorkVec.Add(10, xp2pqpmxp1.Cross(F)+x2pqmx1.Cross(Tmp));
   WorkVec.Add(13, -F);
   WorkVec.Add(16, F.Cross(qTmp));
   WorkVec.Add(19, -Tmp);
   WorkVec.Add(22, F.Cross(Omega2.Cross(qTmp))-qTmp.Cross(Tmp));   
   
   WorkVec.PutCoef(25, v.Dot(p)-vTmp.Dot(x2pqmx1));
   WorkVec.PutCoef(26, x2pqmx1.Dot(vTmp.Cross(Omega1))
		    -vTmp.Dot(xp2pqpmxp1));
      
   return WorkVec;
}

   
/* Setta il valore iniziale delle proprie variabili */
void InPlaneWithOffsetJoint::SetInitialValue(VectorHandler& /* X */ )
{ 
   NO_OP;
}

/* InPlaneWithOffsetJoint - end */
