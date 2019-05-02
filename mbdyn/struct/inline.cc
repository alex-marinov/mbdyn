/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

#include "inline.h"

/* InLineJoint - begin */

const unsigned int InLineJoint::NumSelfDof(2);
const unsigned int InLineJoint::NumDof(14);

/* Costruttore non banale */
InLineJoint::InLineJoint(unsigned int uL, const DofOwner* pDO,
			 const StructNode* pN1, const StructNode* pN2, 
			 const Mat3x3& RvTmp, const Vec3& pTmp, flag fOut,
          const doublereal pref,
          BasicShapeCoefficient *const sh,
          BasicFriction *const f)
: Elem(uL, fOut), Joint(uL, pDO, fOut),
pNode1(pN1), pNode2(pN2), Rv(RvTmp), p(pTmp), F(Zero3),
preF(pref), Sh_c(sh), 
#ifdef USE_NETCDFC // netcdfcxx4 has non-pointer vars...
Var_FF(0),
Var_fc(0),
#endif // USE_NETCDFC
fc(f)
{
   NO_OP;
};

   
InLineJoint::~InLineJoint(void)
{
	if (Sh_c) {
		delete Sh_c;
	}

	if (fc) {
		delete fc;
	}
}

DofOrder::Order
InLineJoint::GetDofType(unsigned int i) const {
   ASSERT(i >= 0 && i < iGetNumDof());
   if (i<NumSelfDof) {
       return DofOrder::ALGEBRAIC; 
   } else {
       return fc->GetDofType(i-NumSelfDof);
   }
};

DofOrder::Order
InLineJoint::GetEqType(unsigned int i) const
{
	ASSERTMSGBREAK(i < iGetNumDof(),
		"INDEX ERROR in InLineJoint::GetEqType");
   if (i<NumSelfDof) {
      return DofOrder::ALGEBRAIC;
   }
	return fc->GetEqType(i-NumSelfDof);
}


/* Contributo al file di restart */
std::ostream& InLineJoint::Restart(std::ostream& out) const
{
   Joint::Restart(out) << ", in line, not implemented yet" << std::endl;
   return out;
}

   
VariableSubMatrixHandler& 
InLineJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		    doublereal dCoef,
		    const VectorHandler& XCurr,
		    const VectorHandler& XPrimeCurr)
{
   DEBUGCOUTFNAME("InLineJoint::AssJac");

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

   /* Full matrix is required for the use of ExpandableRowVectors */
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeReset(iNumRows, iNumCols);

   /* Setta gli indici delle equazioni */
   for (unsigned int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
      WM.PutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
   }
   
   for (unsigned int iCnt = 1; iCnt <= iGetNumDof(); iCnt++) {        // iGet gives self DOFs
      WM.PutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
      WM.PutColIndex(12+iCnt, iFirstReactionIndex+iCnt);
   }
   
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = RvTmp.dGet(iCnt, 1);
      WM.PutCoef(iCnt, 12+1, -d);
      WM.PutCoef(6+iCnt, 12+1, d);	
      
      WM.PutCoef(12+1, iCnt, -d);
      WM.PutCoef(12+1, 6+iCnt, d);
      
      d = RvTmp.dGet(iCnt, 2);
      WM.PutCoef(iCnt, 12+2, -d);
      WM.PutCoef(6+iCnt, 12+2, d);	
      
      WM.PutCoef(12+2, iCnt, -d);
      WM.PutCoef(12+2, 6+iCnt, d);

      d = Tmp1_1.dGet(iCnt);
      WM.PutCoef(3+iCnt, 12+1, d);
      
      WM.PutCoef(12+1, 3+iCnt, d);
      
      d = Tmp1_2.dGet(iCnt);
      WM.PutCoef(3+iCnt, 12+2, d);
      
      WM.PutCoef(12+2, 3+iCnt, d);
      }
      
   WM.PutCross(0, 3, FTmp);  // TODO: change putcross FTmp signs back once FullMat putcross is fixed!
   WM.PutCross(3, 0, -FTmp);
   WM.Put(4, 4, Mat3x3(MatCrossCross, x2mx1, FTmp));
   WM.PutCross(3, 6, FTmp);
   WM.PutCross(6, 6+3, -FTmp);
   
   if (fc) {
       // friction specific contributions:     
       Vec3 e3a(RvTmp.GetVec(3));
       //retrieve friction coef
       doublereal f = fc->fc();
       //shape function
       doublereal shc = Sh_c->Sh_c();
       //compute relative velocity
       doublereal v = (pNode1->GetVCurr()-pNode2->GetVCurr()).Dot(e3a);
       //reaction norm
       doublereal modF = std::max(F.Norm(), preF); // F is not updated inside AssJac?
       //reaction moment
   //doublereal M3 = shc*modF*r;
   
       ExpandableRowVector dfc;
       ExpandableRowVector dF;
       ExpandableRowVector dv;
       //variation of reaction force
       dF.ReDim(3);
       if ((modF == 0.) or (F.Norm() < preF)) {
           dF.Set(Vec3(Zero3),1,12+1); // dF/d(13) ? = dF/d(1st reaction) ?
       } else {
           dF.Set(F/modF,1,12+1);
       }
       //variation of relative velocity
       dv.ReDim(6);
   
       dv.Set(e3a,1, 0+1);
       dv.Set(-e3a,4, 6+1);

       //assemble friction states
       fc->AssJac(WM,dfc,12+NumSelfDof,iFirstReactionIndex+NumSelfDof,dCoef,modF,v,
   	         XCurr,XPrimeCurr,dF,dv);
       ExpandableMatrix dF3;
       ExpandableRowVector dShc;
       //compute variation of shape function
       Sh_c->dSh_c(dShc,f,modF,v,dfc,dF,dv);
       //variation of force component
       dF3.ReDim(3,2);
       dF3.SetBlockDim(1,1);
       dF3.SetBlockDim(2,1);
       dF3.Set(e3a*shc,1,1); dF3.Link(1,&dF); // dF3/dF * dF/d(pos1?)
       dF3.Set(e3a*modF,1,2); dF3.Link(2,&dShc); // dF3/dShc * dShc/d(?)
       //assemble first node variation of force component
       dF3.Add(WM, 1, 1.);
       //assemble second node variation of force component
       dF3.Sub(WM, 6+1, 1.);
   }
   
   return WorkMat;
}


SubVectorHandler& 
InLineJoint::AssRes(SubVectorHandler& WorkVec,
		    doublereal dCoef,
		    const VectorHandler& XCurr, 
		    const VectorHandler& XPrimeCurr )
{
   DEBUGCOUTFNAME("InLineJoint::AssRes");
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);
 
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();

   /* Indici equazioni */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.PutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
   }
   
   /* Indice equazione vincolo */
   for (unsigned int iCnt = 1; iCnt <= iGetNumDof(); iCnt++) { // iGet gives self DOFs
      WorkVec.PutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
   }

   Vec3 x2mx1(pNode2->GetXCurr()-pNode1->GetXCurr());
   Mat3x3 RvTmp = pNode1->GetRCurr()*Rv;

   /* Aggiorna i dati propri */
   F.Put(1, XCurr(iFirstReactionIndex+1));
   F.Put(2, XCurr(iFirstReactionIndex+2));
   Vec3 FTmp(RvTmp*F);
   
   WorkVec.Add(1, FTmp);
   WorkVec.Add(4, x2mx1.Cross(FTmp)); /* ( = -p/\F) */
   WorkVec.Sub(7, FTmp);
   
   ASSERT(dCoef != 0.);
   WorkVec.PutCoef(13, (Rv.GetVec(1).Dot(p)-RvTmp.GetVec(1).Dot(x2mx1))/dCoef);
   WorkVec.PutCoef(14, (Rv.GetVec(2).Dot(p)-RvTmp.GetVec(2).Dot(x2mx1))/dCoef);
   
   if (fc) {
      bool ChangeJac(false);

      Vec3 e3a(RvTmp.GetVec(3));
      doublereal v = (pNode1->GetVCurr()-pNode2->GetVCurr()).Dot(e3a);
      doublereal modF = std::max(F.Norm(), preF);
      try {
          fc->AssRes(WorkVec,12+NumSelfDof,iFirstReactionIndex+NumSelfDof,modF,v,XCurr,XPrimeCurr);
      }
      catch (Elem::ChangedEquationStructure& e) {
          ChangeJac = true;
      }
      doublereal f = fc->fc();
      doublereal shc = Sh_c->Sh_c(f,modF,v);
      F3 = shc*modF;  // or M(3) with a Vec3 ?
      WorkVec.Sub(1,e3a*F3); // subtracting a vector
      WorkVec.Add(7,e3a*F3); // adding a vector
      if (ChangeJac) {
          throw Elem::ChangedEquationStructure(MBDYN_EXCEPT_ARGS);
      }
   }

   return WorkVec;
}

unsigned int InLineJoint::iGetNumDof(void) const {
   unsigned int i = NumSelfDof;
   if (fc) {
       i+=fc->iGetNumDof();
   } 
   return i;
};

void
InLineJoint::OutputPrepare(OutputHandler& OH)
{
	if (bToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			std::string name;
			OutputPrepare_int("inline", OH, name);

			Var_FF = OH.CreateVar<Vec3>(name + "FF", "N",
				"friction force (x, y, z)");

			Var_fc = OH.CreateVar<doublereal>(name + "fc", "-",
				"friction coefficient");
		}
#endif // USE_NETCDF
	}
}

void InLineJoint::Output(OutputHandler& OH) const
{
   if (bToBeOutput()) {
      Mat3x3 RvTmp(pNode1->GetRCurr()*Rv);
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			Joint::NetCDFOutput(OH, F, Zero3, RvTmp*F, Zero3);
			OH.WriteNcVar(Var_FF, F3);
			OH.WriteNcVar(Var_fc, fc->fc());
		}
#endif // USE_NETCDF
		if (OH.UseText(OutputHandler::JOINTS)) {

			std::ostream &of = Joint::Output(OH.Joints(), "inline", GetLabel(),
				F, Zero3, RvTmp*F, Zero3);
		  	// TODO: output relative position and orientation
			if (fc) {
				of << " " << F3 << " " << fc->fc();
			}
			of << std::endl;
		}
   }   
}
 

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
InLineJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
			    const VectorHandler& XCurr)
{
   DEBUGCOUTFNAME("InLineJoint::InitialAssJac");
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.ResizeReset(28, 28);
   
   /* Indici gdl */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+2;
   
   /* Indici equazioni nodi */
   for(int iCnt = 1; iCnt <= 12; iCnt++) {
      WM.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutRowIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
      WM.PutColIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
   }
   
   /* Indici vincoli */
   for(int iCnt = 1; iCnt <= 4; iCnt++) {
      WM.PutRowIndex(24+iCnt, iFirstReactionIndex+iCnt);
      WM.PutColIndex(24+iCnt, iFirstReactionIndex+iCnt);
   }
   
   /* Dati */
   Vec3 x2mx1(pNode2->GetXCurr()-pNode1->GetXCurr());
   Vec3 xp2mxp1(pNode2->GetVCurr()-pNode1->GetVCurr());
   Vec3 Omega(pNode1->GetWRef());
   
   /* Aggiorna i dati propri */
   Vec3 FPrime(XCurr(iReactionPrimeIndex+1),
	       XCurr(iReactionPrimeIndex+2),
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
      WM.PutCoef(iCnt, 25, -d);
      WM.PutCoef(12+iCnt, 25, d);
      
      WM.PutCoef(25, iCnt, -d);
      WM.PutCoef(25, 12+iCnt, d);

      WM.PutCoef(6+iCnt, 27, -d);  
      WM.PutCoef(18+iCnt, 27, d);

      WM.PutCoef(27, 6+iCnt, -d);  
      WM.PutCoef(27, 18+iCnt, d);
      
      d = v2.dGet(iCnt);
      WM.PutCoef(iCnt, 26, -d);
      WM.PutCoef(12+iCnt, 26, d);
      
      WM.PutCoef(26, iCnt, -d);
      WM.PutCoef(26, 12+iCnt, d);

      WM.PutCoef(6+iCnt, 28, -d);  
      WM.PutCoef(18+iCnt, 28, d);

      WM.PutCoef(28, 6+iCnt, -d);  
      WM.PutCoef(28, 18+iCnt, d);

      d = Tmp1_1.dGet(iCnt);
      WM.PutCoef(3+iCnt, 25, d);   
      WM.PutCoef(25, 3+iCnt, d);   

      WM.PutCoef(27, 9+iCnt, d);   
      WM.PutCoef(9+iCnt, 27, d);
      
      d = Tmp1_2.dGet(iCnt);
      WM.PutCoef(3+iCnt, 26, d);   
      WM.PutCoef(26, 3+iCnt, d);   

      WM.PutCoef(28, 9+iCnt, d);   
      WM.PutCoef(9+iCnt, 28, d);
      
      d = Tmp2_1.dGet(iCnt);
      WM.PutCoef(27, iCnt, -d);    
      WM.PutCoef(27, 12+iCnt, d);  

      WM.PutCoef(6+iCnt, 25, -d);  
      WM.PutCoef(18+iCnt, 25, d);
      
      d = Tmp2_2.dGet(iCnt);
      WM.PutCoef(28, iCnt, -d);    
      WM.PutCoef(28, 12+iCnt, d);  

      WM.PutCoef(6+iCnt, 26, -d);  
      WM.PutCoef(18+iCnt, 26, d);

      d = Tmp3_1.dGet(iCnt);
      WM.PutCoef(27, 3+iCnt, d);   

      d = Tmp3_2.dGet(iCnt);
      WM.PutCoef(28, 3+iCnt, d);   

      d = Tmp4_1.dGet(iCnt);
      WM.PutCoef(9+iCnt, 25, d);   

      d = Tmp4_2.dGet(iCnt);
      WM.PutCoef(9+iCnt, 27, d);   
   }   

   Mat3x3 MTmp(MatCross, FTmp);
   WM.Add(1, 4, MTmp);             
   WM.Add(4, 13, MTmp);            
   
   WM.Add(7, 10, MTmp);            
   WM.Add(10, 19, MTmp);           
   
   WM.Sub(4, 1, MTmp);
   WM.Sub(13, 4, MTmp);
   
   WM.Sub(19, 10, MTmp);           
   WM.Sub(10, 7, MTmp);            
   
   MTmp = Mat3x3(MatCrossCross, x2mx1, FTmp);
   WM.Add(4, 4, MTmp);             

   WM.Add(10, 10, MTmp);           
 
   MTmp = Mat3x3(MatCrossCross, Omega, FTmp) + Mat3x3(MatCross, FPrimeTmp);
   WM.Add(7, 4, MTmp);             
   WM.Sub(19, 4, MTmp);
   
   MTmp = Mat3x3(MatCross, Omega.Cross(FTmp) + FPrimeTmp);
   WM.Sub(10, 1, MTmp);
   WM.Add(10, 13, MTmp);
   
   MTmp = Mat3x3((Mat3x3(MatCross, xp2mxp1) + Mat3x3(MatCrossCross, x2mx1, Omega))*Mat3x3(MatCross, FTmp)
		 + Mat3x3(MatCrossCross, x2mx1, FPrimeTmp));
   WM.Add(10, 4, MTmp);

   return WorkMat;
}

   
/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
InLineJoint::InitialAssRes(SubVectorHandler& WorkVec,
			   const VectorHandler& XCurr)
{
   DEBUGCOUTFNAME("InLineJoint::InitialAssRes");
   WorkVec.ResizeReset(28);
   
   /* Indici gdl */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+2;
   
   /* Indici equazioni nodi */
   for(int iCnt = 1; iCnt <= 12; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WorkVec.PutRowIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
   }
   
   /* Indici equazioni vincoli */
   WorkVec.PutRowIndex(25, iFirstReactionIndex+1);
   WorkVec.PutRowIndex(26, iFirstReactionIndex+2);
   WorkVec.PutRowIndex(27, iReactionPrimeIndex+1);
   WorkVec.PutRowIndex(28, iReactionPrimeIndex+2);
      
   /* Dati */
   Vec3 x2mx1(pNode2->GetXCurr()-pNode1->GetXCurr());
   Vec3 xp2mxp1(pNode2->GetVCurr()-pNode1->GetVCurr());
   Mat3x3 RvTmp(pNode1->GetRCurr()*Rv);
   Vec3 Omega(pNode1->GetWCurr());
   
   /* Aggiorna i dati propri */
   F.Put(1, XCurr(iFirstReactionIndex+1));
   F.Put(2, XCurr(iFirstReactionIndex+2));
   Vec3 FPrime(XCurr(iReactionPrimeIndex+1),
	       XCurr(iReactionPrimeIndex+2),
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
   
   WorkVec.PutCoef(25, Rv.GetVec(1).Dot(p)-RvTmp.GetVec(1).Dot(x2mx1));
   WorkVec.PutCoef(26, Rv.GetVec(2).Dot(p)-RvTmp.GetVec(2).Dot(x2mx1));
   WorkVec.PutCoef(27, x2mx1.Dot(RvTmp.GetVec(1).Cross(Omega))-RvTmp.GetVec(1).Dot(xp2mxp1));
   WorkVec.PutCoef(28, x2mx1.Dot(RvTmp.GetVec(2).Cross(Omega))-RvTmp.GetVec(2).Dot(xp2mxp1));
   
   return WorkVec;
}

/* InLineJoint - end */


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
: Elem(uL, fOut), Joint(uL, pDO, fOut),
pNode1(pN1), pNode2(pN2), Rv(RvTmp), p(pTmp), q(qTmp), F(Zero3)
{
   NO_OP;
};

   
InLineWithOffsetJoint::~InLineWithOffsetJoint(void)
{
   NO_OP;
}


/* Contributo al file di restart */
std::ostream& InLineWithOffsetJoint::Restart(std::ostream& out) const
{
   Joint::Restart(out) << ", in line, not implemented yet" << std::endl;
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
   
   WM.ResizeReset(108, 0);
   
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
      WM.PutItem(iCnt, iNode1FirstMomIndex+iCnt,
		  iFirstReactionIndex+1, -d);
      WM.PutItem(3+iCnt, iNode2FirstMomIndex+iCnt,
		  iFirstReactionIndex+1, d);      
      
      WM.PutItem(6+iCnt, iFirstReactionIndex+1,
		  iNode1FirstPosIndex+iCnt, -d);
      WM.PutItem(9+iCnt, iFirstReactionIndex+1,
		  iNode2FirstPosIndex+iCnt, d);
      
      d = RvTmp.dGet(iCnt, 2);
      WM.PutItem(12+iCnt, iNode1FirstMomIndex+iCnt,
		  iFirstReactionIndex+2, -d);
      WM.PutItem(15+iCnt, iNode2FirstMomIndex+iCnt,
		  iFirstReactionIndex+2, d);      
      
      WM.PutItem(18+iCnt, iFirstReactionIndex+2,
		  iNode1FirstPosIndex+iCnt, -d);
      WM.PutItem(21+iCnt, iFirstReactionIndex+2,
		  iNode2FirstPosIndex+iCnt, d);

      d = Tmp1_1.dGet(iCnt);
      WM.PutItem(24+iCnt, iNode1FirstMomIndex+3+iCnt,
		  iFirstReactionIndex+1, d);
      
      WM.PutItem(27+iCnt, iFirstReactionIndex+1,
		  iNode1FirstPosIndex+3+iCnt, d);
      
      d = Tmp1_2.dGet(iCnt);
      WM.PutItem(30+iCnt, iNode1FirstMomIndex+3+iCnt,
		  iFirstReactionIndex+2, d);
      
      WM.PutItem(33+iCnt, iFirstReactionIndex+2,
		  iNode1FirstPosIndex+3+iCnt, d);

      d = Tmp2_1.dGet(iCnt);
      WM.PutItem(36+iCnt, iNode2FirstMomIndex+3+iCnt,
		  iFirstReactionIndex+1, d);
      
      WM.PutItem(39+iCnt, iFirstReactionIndex+1,
		  iNode2FirstPosIndex+3+iCnt, d);

      d = Tmp2_2.dGet(iCnt);
      WM.PutItem(42+iCnt, iNode2FirstMomIndex+3+iCnt,
		  iFirstReactionIndex+2, d);
      
      WM.PutItem(45+iCnt, iFirstReactionIndex+2,
		  iNode2FirstPosIndex+3+iCnt, d);
   }
   
   WM.PutCross(48+1, iNode1FirstMomIndex,
		iNode1FirstPosIndex+3, FTmp);
   WM.PutCross(54+1, iNode1FirstMomIndex+3,
		iNode1FirstPosIndex, -FTmp);
   WM.PutMat3x3(60+1, iNode1FirstMomIndex+3,		
		iNode1FirstPosIndex+3, Mat3x3(MatCrossCross, x2qmx1, FTmp));
   WM.PutCross(69+1, iNode1FirstMomIndex+3,
		iNode2FirstPosIndex, FTmp);
   WM.PutCross(75+1, iNode2FirstMomIndex,
		iNode2FirstPosIndex+3, -FTmp);

   Mat3x3 MTmp(MatCrossCross, FTmp, qTmp);
   WM.PutMat3x3(81+1, iNode1FirstMomIndex+3,
		 iNode2FirstPosIndex+3, -MTmp);
   WM.PutMat3x3(90+1, iNode2FirstMomIndex+3,
		 iNode2FirstPosIndex+3, -MTmp);

   WM.PutMat3x3(99+1, iNode2FirstMomIndex+3,		
		 iNode1FirstPosIndex+3, Mat3x3(MatCrossCross, -qTmp, FTmp));
   
   return WorkMat;
}


SubVectorHandler& 
InLineWithOffsetJoint::AssRes(SubVectorHandler& WorkVec,
			      doublereal dCoef,
			      const VectorHandler& XCurr, 
			      const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUTFNAME("InLineWithOffsetJoint::AssRes");
   WorkVec.ResizeReset(14);
 
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
   WorkVec.PutRowIndex(14, iFirstReactionIndex+2);

   Vec3 qTmp(pNode2->GetRCurr()*q);
   Vec3 x2qmx1(pNode2->GetXCurr()+qTmp-pNode1->GetXCurr());
   Mat3x3 RvTmp = pNode1->GetRCurr()*Rv;

   /* Aggiorna i dati propri */
   F.Put(1, XCurr(iFirstReactionIndex+1));
   F.Put(2, XCurr(iFirstReactionIndex+2));
   Vec3 FTmp(RvTmp*F);
   
   WorkVec.Add(1, FTmp);
   WorkVec.Add(4, x2qmx1.Cross(FTmp));
   WorkVec.Sub(7, FTmp);
   WorkVec.Sub(10, qTmp.Cross(FTmp));
   
   ASSERT(dCoef != 0.);
   WorkVec.PutCoef(13, (Rv.GetVec(1).Dot(p)-RvTmp.GetVec(1).Dot(x2qmx1))/dCoef);
   WorkVec.PutCoef(14, (Rv.GetVec(2).Dot(p)-RvTmp.GetVec(2).Dot(x2qmx1))/dCoef);

   return WorkVec;
}

DofOrder::Order
InLineWithOffsetJoint::GetEqType(unsigned int i) const
{
	ASSERTMSGBREAK((i>=0) and (i<15), 
		"INDEX ERROR in InLineWithOffsetJoint::GetEqType");
	if ((i==13) or (i==14)) {
		return DofOrder::ALGEBRAIC;
	}
	return DofOrder::DIFFERENTIAL;
}

void InLineWithOffsetJoint::Output(OutputHandler& OH) const
{
   if (bToBeOutput()) {      
      Mat3x3 RvTmp(pNode1->GetRCurr()*Rv);
      Joint::Output(OH.Joints(), "inline", GetLabel(),
		    F, Zero3, RvTmp*F, Zero3) << std::endl;
   }
}
 

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
InLineWithOffsetJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
				     const VectorHandler& XCurr)
{
   DEBUGCOUTFNAME("InLineWithOffsetJoint::InitialAssJac");
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.ResizeReset(28, 28);
   
   /* Indici gdl */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+2;
   
   /* Indici equazioni nodi */
   for(int iCnt = 1; iCnt <= 12; iCnt++) {
      WM.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutRowIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
      WM.PutColIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
   }
   
   /* Indici vincoli */
   for(int iCnt = 1; iCnt <= 4; iCnt++) {
      WM.PutRowIndex(24+iCnt, iFirstReactionIndex+iCnt);
      WM.PutColIndex(24+iCnt, iFirstReactionIndex+iCnt);
   }
   
   /* Dati */
   Vec3 Omega1(pNode1->GetWRef());
   Vec3 Omega2(pNode2->GetWRef());
   Mat3x3 RvTmp(pNode1->GetRRef()*Rv);
   Vec3 qTmp(pNode2->GetRRef()*q);
   Vec3 x2qmx1(pNode2->GetXCurr()+qTmp-pNode1->GetXCurr());
   Vec3 xp2qmxp1(pNode2->GetVCurr()+Omega2.Cross(qTmp)-pNode1->GetVCurr());
   
   /* Aggiorna i dati propri */
   Vec3 FPrime(XCurr(iReactionPrimeIndex+1),
	       XCurr(iReactionPrimeIndex+2),
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
      WM.PutCoef(iCnt, 25, -d);
      WM.PutCoef(12+iCnt, 25, d);
      
      WM.PutCoef(25, iCnt, -d);
      WM.PutCoef(25, 12+iCnt, d);

      WM.PutCoef(6+iCnt, 27, -d);  
      WM.PutCoef(18+iCnt, 27, d);

      WM.PutCoef(27, 6+iCnt, -d);  
      WM.PutCoef(27, 18+iCnt, d);
      
      d = v2.dGet(iCnt);
      WM.PutCoef(iCnt, 26, -d);
      WM.PutCoef(12+iCnt, 26, d);
      
      WM.PutCoef(26, iCnt, -d);
      WM.PutCoef(26, 12+iCnt, d);

      WM.PutCoef(6+iCnt, 28, -d);  
      WM.PutCoef(18+iCnt, 28, d);

      WM.PutCoef(28, 6+iCnt, -d);  
      WM.PutCoef(28, 18+iCnt, d);

      d = Tmp1_1.dGet(iCnt);
      WM.PutCoef(3+iCnt, 25, d);   
      WM.PutCoef(25, 3+iCnt, d);   

      WM.PutCoef(27, 9+iCnt, d);   
      WM.PutCoef(9+iCnt, 27, d);
      
      d = Tmp1_2.dGet(iCnt);
      WM.PutCoef(3+iCnt, 26, d);   
      WM.PutCoef(26, 3+iCnt, d);   

      WM.PutCoef(28, 9+iCnt, d);   
      WM.PutCoef(9+iCnt, 28, d);
      
      d = Tmp2_1.dGet(iCnt);
      WM.PutCoef(27, iCnt, -d);    
      WM.PutCoef(27, 12+iCnt, d);  

      WM.PutCoef(6+iCnt, 25, -d);  
      WM.PutCoef(18+iCnt, 25, d);
      
      d = Tmp2_2.dGet(iCnt);
      WM.PutCoef(28, iCnt, -d);    
      WM.PutCoef(28, 12+iCnt, d);  

      WM.PutCoef(6+iCnt, 26, -d);  
      WM.PutCoef(18+iCnt, 26, d);

      d = Tmp3_1.dGet(iCnt);
      WM.PutCoef(27, 3+iCnt, d);   

      d = Tmp3_2.dGet(iCnt);
      WM.PutCoef(28, 3+iCnt, d);   

      d = Tmp4_1.dGet(iCnt);
      WM.PutCoef(9+iCnt, 25, d);   

      d = Tmp4_2.dGet(iCnt);
      WM.PutCoef(9+iCnt, 27, d);   
      
      d = Tmp5_1.dGet(iCnt);
      WM.PutCoef(15+iCnt, 25, d);
      WM.PutCoef(21+iCnt, 27, d);

      WM.PutCoef(25, 15+iCnt, d);
      WM.PutCoef(27, 21+iCnt, d);

      d = Tmp5_2.dGet(iCnt);
      WM.PutCoef(15+iCnt, 26, d);
      WM.PutCoef(21+iCnt, 28, d);

      WM.PutCoef(26, 15+iCnt, d);
      WM.PutCoef(28, 21+iCnt, d);

      d = Tmp6_1.dGet(iCnt);
      WM.PutCoef(27, 15+iCnt, d);

      d = Tmp6_2.dGet(iCnt);
      WM.PutCoef(28, 15+iCnt, d);

      d = Tmp7_1.dGet(iCnt);
      WM.PutCoef(21+iCnt, 25, d);
      
      d = Tmp7_2.dGet(iCnt);
      WM.PutCoef(21+iCnt, 26, d);
   }   

   Mat3x3 MTmp(MatCross, FTmp);
   WM.Add(1, 4, MTmp);             
   WM.Add(4, 13, MTmp);            
   
   WM.Add(7, 10, MTmp);            
   WM.Add(10, 19, MTmp);           
   
   WM.Sub(4, 1, MTmp);             
   WM.Sub(13, 4, MTmp);            
   
   WM.Sub(19, 10, MTmp);           
   WM.Sub(10, 7, MTmp);            
   
   MTmp = Mat3x3(MatCrossCross, x2qmx1, FTmp);
   WM.Add(4, 4, MTmp);             

   WM.Add(10, 10, MTmp);           
 
   MTmp = Mat3x3(MatCrossCross, Omega1, FTmp) + Mat3x3(MatCross, FPrimeTmp);
   WM.Add(7, 4, MTmp);             
   WM.Sub(19, 4, MTmp);           
   
   MTmp = Mat3x3(MatCross, Omega1.Cross(FTmp) + FPrimeTmp);
   WM.Sub(10, 1, MTmp);
   WM.Add(10, 13, MTmp);
   
   MTmp = Mat3x3((Mat3x3(MatCross, xp2qmxp1) + Mat3x3(MatCrossCross, x2qmx1, Omega1))*Mat3x3(MatCross, FTmp)
		 + Mat3x3(MatCrossCross, x2qmx1, FPrimeTmp));
   WM.Add(10, 4, MTmp);

   MTmp = Mat3x3(MatCrossCross, FTmp, qTmp);
   WM.Add(16, 16, MTmp);
   WM.Add(22, 22, MTmp);
   
   WM.Sub(4, 16, MTmp);
   WM.Sub(10, 22, MTmp);
   
   MTmp = MTmp.Transpose();
   WM.Add(16, 4, MTmp);
   WM.Add(22, 10, MTmp);
   
   MTmp = (Mat3x3(MatCrossCross, FTmp, Omega2) + Mat3x3(MatCross, Omega1.Cross(FTmp) + FPrimeTmp))*Mat3x3(MatCross, qTmp);
   WM.Add(22, 16, MTmp);
   WM.Sub(10, 16, MTmp);
   
   MTmp = Mat3x3(MatCrossCross, qTmp.Cross(Omega2), FTmp)
     - Mat3x3(MatCross, qTmp)*(Mat3x3(MatCrossCross, Omega1, FTmp) + Mat3x3(MatCross, FPrimeTmp));
   WM.Add(22, 4, MTmp);

   return WorkMat;
}

   
/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
InLineWithOffsetJoint::InitialAssRes(SubVectorHandler& WorkVec,
				     const VectorHandler& XCurr)
{
   DEBUGCOUTFNAME("InLineWithOffsetJoint::InitialAssRes");
   WorkVec.ResizeReset(28);
   
   /* Indici gdl */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+2;
   
   /* Indici equazioni nodi */
   for(int iCnt = 1; iCnt <= 12; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WorkVec.PutRowIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
   }
   
   /* Indici equazioni vincoli */
   WorkVec.PutRowIndex(25, iFirstReactionIndex+1);
   WorkVec.PutRowIndex(26, iFirstReactionIndex+2);
   WorkVec.PutRowIndex(27, iReactionPrimeIndex+1);
   WorkVec.PutRowIndex(28, iReactionPrimeIndex+2);
      
   /* Dati */
   Mat3x3 RvTmp(pNode1->GetRCurr()*Rv);
   Vec3 qTmp(pNode2->GetRCurr()*q);
   Vec3 Omega1(pNode1->GetWCurr());
   Vec3 Omega2(pNode2->GetWCurr());
   Vec3 x2qmx1(pNode2->GetXCurr()+qTmp-pNode1->GetXCurr());
   Vec3 xp2qmxp1(pNode2->GetVCurr()+Omega2.Cross(qTmp)-pNode1->GetVCurr());
   
   /* Aggiorna i dati propri */
   F.Put(1, XCurr(iFirstReactionIndex+1));
   F.Put(2, XCurr(iFirstReactionIndex+2));
   Vec3 FPrime(XCurr(iReactionPrimeIndex+1),
	       XCurr(iReactionPrimeIndex+2),
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
   
   WorkVec.PutCoef(25, Rv.GetVec(1).Dot(p)-RvTmp.GetVec(1).Dot(x2qmx1));
   WorkVec.PutCoef(26, Rv.GetVec(2).Dot(p)-RvTmp.GetVec(2).Dot(x2qmx1));
   WorkVec.PutCoef(27, x2qmx1.Dot(RvTmp.GetVec(1).Cross(Omega1))-RvTmp.GetVec(1).Dot(xp2qmxp1));
   WorkVec.PutCoef(28, x2qmx1.Dot(RvTmp.GetVec(2).Cross(Omega1))-RvTmp.GetVec(2).Dot(xp2qmxp1));
   
   return WorkVec;
}

/* InLineWithOffsetJoint - end */

