/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2007
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

#include <kinj.h>

/* Costruttore */
KinJoint::KinJoint(unsigned int uL, 
		   const DofOwner* pDO,
		   const StructNode* pN,
		   const Kinematics* pK,
		   flag fOut)
: Elem(uL, fOut),
Joint(uL, pDO, fOut),
pNode(pN), pKin(pK)
{
   ASSERT(pNode != NULL);
   ASSERT(pNode->GetNodeType() == Node::STRUCTURAL);
   ASSERT(pNode->GetStructNodeType() == StructNode::STATIC);
}


KinJoint::~KinJoint(void)
{
   SAFEDELETE(pKin);
}

/* Contributo al file di restart */
std::ostream& 
KinJoint::Restart(std::ostream& out) const
{
   return out << "not implemented yet" << std::endl;
}


VariableSubMatrixHandler& 
KinJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		 doublereal dCoef,
		 const VectorHandler& XCurr, 
		 const VectorHandler& XPrimeCurr)
{
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.ResizeReset(24, 24);
   
   integer iNodeFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
   integer iNodeFirstPositionIndex = pNode->iGetFirstPositionIndex();
   integer iFirstIndex = iGetFirstIndex();
   
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.PutRowIndex(iCnt, iNodeFirstMomentumIndex+iCnt);
      WM.PutRowIndex(6+iCnt, iFirstIndex+iCnt);
      WM.PutRowIndex(12+iCnt, iFirstIndex+6+iCnt);
      WM.PutRowIndex(18+iCnt, iFirstIndex+12+iCnt);      
      
      WM.PutColIndex(iCnt, iNodeFirstPositionIndex+iCnt);
      WM.PutColIndex(6+iCnt, iFirstIndex+iCnt);
      WM.PutColIndex(12+iCnt, iFirstIndex+6+iCnt);
      WM.PutColIndex(18+iCnt, iFirstIndex+12+iCnt);      
   }
     
   /* suppongo dati: */   
   /* Vec3 X0(pKin->GetXCurr()); */
   Mat3x3 R0(pKin->GetRCurr());
   /* Vec3 V0(pKin->GetVCurr()); */
   /* Vec3 W0(pKin->GetWCurr()); */
   
   Mat3x3 R(pNode->GetRRef());
   Vec3 W(pNode->GetWRef());
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutCoef(iCnt, 18+iCnt, -1.);
      
      WM.PutCoef(6+iCnt, iCnt, 1.);
      WM.PutCoef(6+iCnt, 6+iCnt, -1.);
      WM.PutCoef(6+iCnt, 12+iCnt, 1.);      

      WM.PutCoef(9+iCnt, 3+iCnt, 1.);
      WM.PutCoef(9+iCnt, 9+iCnt, -1.);
      WM.PutCoef(9+iCnt, 15+iCnt, 1.);
      
      WM.PutCoef(12+iCnt, iCnt, dCoef);
      
      WM.PutCoef(18+iCnt, 6+iCnt, 1.);
      WM.PutCoef(21+iCnt, 9+iCnt, 1.);
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
   WorkVec.ResizeReset(24);
 
   integer iNodeFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
   integer iFirstIndex = iGetFirstIndex();
   
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iNodeFirstMomentumIndex+iCnt);
      WorkVec.PutRowIndex(6+iCnt, iFirstIndex+iCnt);
      WorkVec.PutRowIndex(12+iCnt, iFirstIndex+6+iCnt);
      WorkVec.PutRowIndex(18+iCnt, iFirstIndex+12+iCnt);      
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
   
   DEBUGCOUT(std::endl
	     << "v = " << v << std::endl
	     << "w = " << w << std::endl
	     << "mu_v = " << mu_v << std::endl
	     << "mu_w = " << mu_w << std::endl
	     << "l_v = " << lambda_v << std::endl
	     << "l_w = " << lambda_w << std::endl);
   
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
	<< mu_v << " " << mu_w << std::endl;
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
KinJoint::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& /* XP */ ,
		SimulationEntity::Hints *ph)
{
   integer iIndex = iGetFirstIndex();
   X.Put(iIndex+1, pNode->GetVCurr());
   X.Put(iIndex+4, pNode->GetWCurr());
}
