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

#include <mbconfig.h>

#include <kinj.h>

/* Costruttore */
KinJoint::KinJoint(unsigned int uL, 
		   const DofOwner* pDO,
		   const StructNode* pN,
		   const Kinematics* pK,
		   flag fOut)
: Elem(uL, ElemType::JOINT, fOut),
Joint(uL, JointType::IMPOSEDKINEMATICS, pDO, fOut),
pNode(pN), pKin(pK)
{
   ASSERT(pNode != NULL);
   ASSERT(pNode->GetNodeType() == NodeType::STRUCTURAL);
   ASSERT(pNode->GetStructNodeType() == StructNodeType::STATIC);
}


KinJoint::~KinJoint(void)
{
   SAFEDELETE(pKin, DMmm);
}

/* Contributo al file di restart */
ostream& 
KinJoint::Restart(ostream& out) const
{
   return out << "not implemented yet" << endl;
}


VariableSubMatrixHandler& 
KinJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		 doublereal dCoef,
		 const VectorHandler& XCurr, 
		 const VectorHandler& XPrimeCurr)
{
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.ResizeInit(24, 24, 0.);
   
   integer iNodeFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
   integer iNodeFirstPositionIndex = pNode->iGetFirstPositionIndex();
   integer iFirstIndex = iGetFirstIndex();
   
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.fPutRowIndex(iCnt, iNodeFirstMomentumIndex+iCnt);
      WM.fPutRowIndex(6+iCnt, iFirstIndex+iCnt);
      WM.fPutRowIndex(12+iCnt, iFirstIndex+6+iCnt);
      WM.fPutRowIndex(18+iCnt, iFirstIndex+12+iCnt);      
      
      WM.fPutColIndex(iCnt, iNodeFirstPositionIndex+iCnt);
      WM.fPutColIndex(6+iCnt, iFirstIndex+iCnt);
      WM.fPutColIndex(12+iCnt, iFirstIndex+6+iCnt);
      WM.fPutColIndex(18+iCnt, iFirstIndex+12+iCnt);      
   }
     
   /* suppongo dati: */   
   /* Vec3 X0(pKin->GetXCurr()); */
   Mat3x3 R0(pKin->GetRCurr());
   /* Vec3 V0(pKin->GetVCurr()); */
   /* Vec3 W0(pKin->GetWCurr()); */
   
   Mat3x3 R(pNode->GetRRef());
   Vec3 W(pNode->GetWRef());
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.fPutCoef(iCnt, 18+iCnt, -1.);
      
      WM.fPutCoef(6+iCnt, iCnt, 1.);
      WM.fPutCoef(6+iCnt, 6+iCnt, -1.);
      WM.fPutCoef(6+iCnt, 12+iCnt, 1.);      

      WM.fPutCoef(9+iCnt, 3+iCnt, 1.);
      WM.fPutCoef(9+iCnt, 9+iCnt, -1.);
      WM.fPutCoef(9+iCnt, 15+iCnt, 1.);
      
      WM.fPutCoef(12+iCnt, iCnt, dCoef);
      
      WM.fPutCoef(18+iCnt, 6+iCnt, 1.);
      WM.fPutCoef(21+iCnt, 9+iCnt, 1.);
   }
   
   Mat3x3 H(R.GetVec(2).Cross(R0.GetVec(3)),
	    R.GetVec(3).Cross(R0.GetVec(1)),
	    R.GetVec(1).Cross(R0.GetVec(2)));
   
   WM.Sub(4, 22, H);
   WM.Add(16, 4, H.Transpose()*dCoef);
   
   WM.Sub(4, 4,
	  Mat3x3(R0.GetVec(3), R.GetVec(2)*(dCoef*lambda_w.dGet(1)))
	  +Mat3x3(R0.GetVec(1), R.GetVec(3)*(dCoef*lambda_w.dGet(2)))
	  +Mat3x3(R0.GetVec(2), R.GetVec(1)*(dCoef*lambda_w.dGet(3))));
   
   WM.Sub(10, 4, Mat3x3(W*dCoef));
   
   return WorkMat;
}


SubVectorHandler& 
KinJoint::AssRes(SubVectorHandler& WorkVec,
		 doublereal dCoef,
		 const VectorHandler& XCurr, 
		 const VectorHandler& XPrimeCurr)
{
   WorkVec.Resize(24);
   WorkVec.Reset(0.);
   
   integer iNodeFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
   integer iFirstIndex = iGetFirstIndex();
   
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iNodeFirstMomentumIndex+iCnt);
      WorkVec.fPutRowIndex(6+iCnt, iFirstIndex+iCnt);
      WorkVec.fPutRowIndex(12+iCnt, iFirstIndex+6+iCnt);
      WorkVec.fPutRowIndex(18+iCnt, iFirstIndex+12+iCnt);      
   }
   
   v = Vec3(XCurr, iFirstIndex+1);
   w = Vec3(XCurr, iFirstIndex+4);
   mu_v = Vec3(XCurr, iFirstIndex+7);
   mu_w = Vec3(XCurr, iFirstIndex+10);
   lambda_v = Vec3(XCurr, iFirstIndex+13);
   lambda_w = Vec3(XCurr, iFirstIndex+16);
  
   /* suppongo dati: */
   Vec3 X0(pKin->GetXCurr());
   Mat3x3 R0(pKin->GetRCurr());
   Vec3 V0(pKin->GetVCurr());
   Vec3 W0(pKin->GetWCurr());
   
   Vec3 X(pNode->GetXCurr());
   Mat3x3 R(pNode->GetRCurr());
   Vec3 V(pNode->GetVCurr());
   Vec3 W(pNode->GetWCurr());
   
   DEBUGCOUT(endl
	     << "v = " << v << endl
	     << "w = " << w << endl
	     << "mu_v = " << mu_v << endl
	     << "mu_w = " << mu_w << endl
	     << "l_v = " << lambda_v << endl
	     << "l_w = " << lambda_w << endl);
   
   Mat3x3 H(R.GetVec(2).Cross(R0.GetVec(3)),
	    R.GetVec(3).Cross(R0.GetVec(1)),
	    R.GetVec(1).Cross(R0.GetVec(2)));
   
   WorkVec.Put(1, lambda_v);
   WorkVec.Put(4, H*lambda_w);
   
   WorkVec.Put(7, v-mu_v-V);
   WorkVec.Put(10, w-H*mu_w-W);
   
   WorkVec.Put(13, X0-X);
   WorkVec.Put(16, Vec3(-R.GetVec(2).Dot(R0.GetVec(3)),
			-R.GetVec(3).Dot(R0.GetVec(1)),
			-R.GetVec(1).Dot(R0.GetVec(2))));
   
   WorkVec.Put(19, V0-v);
   WorkVec.Put(22, W0-w);
   
   return WorkVec;
}


void 
KinJoint::Output(OutputHandler& OH) const
{
   if(fToBeOutput()) {
      Mat3x3 R(pNode->GetRCurr());
      Mat3x3 RT(R.Transpose());
      
      Joint::Output(OH.Joints(), "Kinematic", GetLabel(),
                    RT*lambda_v, lambda_w, lambda_v, R*lambda_w) << " " 
	<< mu_v << " " << mu_w << endl;
   }
}

   
/* funzioni usate nell'assemblaggio iniziale */

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
KinJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
			const VectorHandler& XCurr)
{
   WorkMat.SetNullMatrix();
   return WorkMat;
}
   
/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
KinJoint::InitialAssRes(SubVectorHandler& WorkVec,
			const VectorHandler& XCurr)
{
   WorkVec.Resize(0);
   return WorkVec;
}

/* Setta il valore iniziale delle proprie variabili */
void 
KinJoint::SetInitialValue(VectorHandler& /* X */ ) const
{
   NO_OP;
}

void 
KinJoint::SetValue(VectorHandler& X, VectorHandler& /* XP */ ) const
{
   integer iIndex = iGetFirstIndex();
   X.Put(iIndex+1, pNode->GetVCurr());
   X.Put(iIndex+4, pNode->GetWCurr());
}
