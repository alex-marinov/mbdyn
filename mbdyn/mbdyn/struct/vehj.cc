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

/* Cerniera deformabile */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <vehj.h>

#include <matvecexp.h>
#include <Rot.hh>

/* DeformableHingeJoint - begin */

/* Costruttore non banale */
DeformableHingeJoint::DeformableHingeJoint(unsigned int uL,
					   DefHingeType::Type T,
					   const DofOwner* pDO, 
					   const ConstitutiveLaw3D* pCL,
					   const StructNode* pN1, 
					   const StructNode* pN2,
					   const Mat3x3& R1,
					   const Mat3x3& R2,
					   flag fOut)
: Elem(uL, Elem::JOINT, fOut), 
Joint(uL, Joint::DEFORMABLEHINGE, pDO, fOut), 
ConstitutiveLaw3DOwner(pCL), DefHingeT(T),
pNode1(pN1), pNode2(pN2), R1h(R1), R2h(R2), fFirstRes(1)
{
   ASSERT(pNode1 != NULL);
   ASSERT(pNode2 != NULL);
   ASSERT(pNode1->GetNodeType() == Node::STRUCTURAL);
   ASSERT(pNode2->GetNodeType() == Node::STRUCTURAL);
}

   
/* Distruttore */
DeformableHingeJoint::~DeformableHingeJoint(void)
{
   NO_OP;
}

   
/* Contributo al file di restart */
std::ostream& DeformableHingeJoint::Restart(std::ostream& out) const
{
   Joint::Restart(out) << ", deformable hinge, "
     << pNode1->GetLabel() << ", reference, node, 1, ",
     (R1h.GetVec(1)).Write(out, ", ")
     << ", 2, ", (R1h.GetVec(2)).Write(out, ", ") << ", "
     << pNode2->GetLabel() << ", reference, node, 1, ",
     (R2h.GetVec(1)).Write(out, ", ")
     << ", 2, ", (R2h.GetVec(2)).Write(out, ", ") << ", ";
   return pGetConstLaw()->Restart(out) << ';' << std::endl;
}


void DeformableHingeJoint::Output(OutputHandler& OH) const
{   
   if (fToBeOutput()) {
      Vec3 d(EulerAngles(pNode1->GetRCurr().Transpose()*pNode2->GetRCurr()));
      Vec3 v(GetF());
      Joint::Output(OH.Joints(), "DeformableHinge", GetLabel(),
		    Zero3, v, Zero3, pNode1->GetRCurr()*v) 
	<< " " << d << std::endl;    
   }     
}

/* DeformableHingeJoint - end */


/* ElasticHingeJoint - begin */

ElasticHingeJoint::ElasticHingeJoint(unsigned int uL, 
				     const DofOwner* pDO, 
				     const ConstitutiveLaw3D* pCL,
				     const StructNode* pN1, 
				     const StructNode* pN2,
				     const Mat3x3& R1,
				     const Mat3x3& R2,
				     flag fOut)
: Elem(uL, Elem::JOINT, fOut), 
DeformableHingeJoint(uL, DefHingeType::ELASTIC, 
		     pDO, pCL, pN1, pN2, R1, R2, fOut),
ThetaRef(0.), TaCurr(0.), TbCurr(0.), FDE(0.)
{
   /* 
    * Chiede la matrice tangente di riferimento e la porta nel sistema globale 
    */
   Mat3x3 R1(pNode1->GetRRef()*R1h);
   FDE = R1*GetFDE()*R1.Transpose();   
}


ElasticHingeJoint::~ElasticHingeJoint(void)
{
   NO_OP;
}

   
/* assemblaggio jacobiano */
VariableSubMatrixHandler& 
ElasticHingeJoint::AssJac(VariableSubMatrixHandler& WorkMat,
			  doublereal dCoef, 
			  const VectorHandler& /* XCurr */ ,
			  const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering ElasticHingeJoint::AssJac()" << std::endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeInit(iNumRows, iNumCols, 0.);

   /* Recupera gli indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex()+3;
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex()+3;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex()+3;
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex()+3;
   
   /* Setta gli indici della matrice */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WM.fPutColIndex(iCnt, iNode1FirstPosIndex+iCnt);	
      WM.fPutRowIndex(3+iCnt, iNode2FirstMomIndex+iCnt);
      WM.fPutColIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
   }  
   
   AssMat(WM, dCoef);

   return WorkMat;
}


void ElasticHingeJoint::AfterPredict(VectorHandler& /* X */ ,
				     VectorHandler& /* XP */ )
{
   /* Calcola le deformazioni, aggiorna il legame costitutivo e crea la FDE */

   /* Recupera i dati */
   Mat3x3 R1(pNode1->GetRRef()*R1h);
   Mat3x3 R2(pNode2->GetRRef()*R2h);
   Mat3x3 R1T(R1.Transpose());
   
   /* Calcola la deformazione corrente nel sistema locale (nodo a) */
   ThetaCurr = ThetaRef = RotManip::VecRot(R1T*R2);

   /* Calcola l'inversa di Gamma di ThetaRef */
   Mat3x3 GammaRefm1 = RotManip::DRot_I(ThetaRef);
   
   /* Aggiorna il legame costitutivo */
   ConstitutiveLaw3DOwner::Update(ThetaRef);
      
   /* Chiede la matrice tangente di riferimento e la porta 
    * nel sistema globale */
   FDE = R1*ConstitutiveLaw3DOwner::GetFDE()*GammaRefm1*R1T;
						     
   fFirstRes = flag(1);						     
}


void ElasticHingeJoint::AssMat(FullSubMatrixHandler& WM, doublereal dCoef)   
{   
   Mat3x3 R1(pNode1->GetRRef()*R1h);
  
   Mat3x3 FDETmp = FDE*dCoef;
   
   WM.Add(4, 4, FDETmp);
   WM.Sub(1, 4, FDETmp);

   FDETmp += Mat3x3(R1*(GetF()*dCoef));
   WM.Add(1, 1, FDETmp);   
   WM.Sub(4, 1, FDETmp);
}


/* assemblaggio residuo */
SubVectorHandler& 
ElasticHingeJoint::AssRes(SubVectorHandler& WorkVec,
			  doublereal /* dCoef */ ,
			  const VectorHandler& /* XCurr */ ,
			  const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering ElasticHingeJoint::AssRes()" << std::endl);   
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
   

   /* Recupera gli indici */
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex()+3;
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex()+3;
   
   /* Setta gli indici della matrice */
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.fPutRowIndex(3+iCnt, iNode2FirstMomIndex+iCnt);
   }  
   
   AssVec(WorkVec);
   
   return WorkVec;
}

   
void ElasticHingeJoint::AssVec(SubVectorHandler& WorkVec)
{
   Mat3x3 R1(pNode1->GetRCurr()*R1h);

   if (fFirstRes) {
      fFirstRes = 0;

   } else {
      Mat3x3 R2(pNode2->GetRCurr()*R2h);
      ThetaCurr = RotManip::VecRot(R1.Transpose()*R2);
      ConstitutiveLaw3DOwner::Update(ThetaCurr);
   }

   Vec3 F(R1*GetF());
   
   WorkVec.Add(1, F);
   WorkVec.Sub(4, F);
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
ElasticHingeJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
				 const VectorHandler& /* XCurr */ )
{
   DEBUGCOUT("Entering ElasticHingeJoint::InitialAssJac()" << std::endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeInit(iNumRows, iNumCols, 0.);

   /* Recupera gli indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex()+3;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex()+3;
   
   /* Setta gli indici della matrice */
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.fPutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.fPutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.fPutRowIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
      WM.fPutColIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
   }  
   
   AssMat(WM, 1.);

   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
ElasticHingeJoint::InitialAssRes(SubVectorHandler& WorkVec,
				 const VectorHandler& /* XCurr */ )
{
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
   

   /* Recupera gli indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex()+3;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex()+3;
   
   /* Setta gli indici della matrice */
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WorkVec.fPutRowIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
   }  
   
   AssVec(WorkVec);
   
   return WorkVec;
}

/* ElasticHingeJoint - end */


/* ViscousHingeJoint - begin */

ViscousHingeJoint::ViscousHingeJoint(unsigned int uL, 
				     const DofOwner* pDO, 
				     const ConstitutiveLaw3D* pCL,
				     const StructNode* pN1, 
				     const StructNode* pN2,
				     const Mat3x3& R1,
				     const Mat3x3& R2,
				     flag fOut)
: Elem(uL, Elem::JOINT, fOut), 
DeformableHingeJoint(uL, DefHingeType::VISCOUS, 
		     pDO, pCL, pN1, pN2, R1, R2, fOut),
ThetaRefPrime(0.), ThetaCurrPrime(0.), TaCurrPrime(0.), TbCurrPrime(0.)
{
   Mat3x3 Ra(pNode1->GetRRef());
   FDEPrime = Ra*ConstitutiveLaw3DOwner::GetFDEPrime()*Ra.Transpose();
}


ViscousHingeJoint::~ViscousHingeJoint(void)
{
   NO_OP;
}

   
void ViscousHingeJoint::AfterPredict(VectorHandler& /* X */ ,
				     VectorHandler& /* XP */ )
{
   // Calcola le deformazioni, aggiorna il legame costitutivo e crea la FDE
   
   // Recupera i dati
   Mat3x3 Ra(pNode1->GetRRef()*R1h);
   Mat3x3 RaT(Ra.Transpose());
   
   Vec3 ga(pNode1->GetgRef());
   Vec3 gb(pNode2->GetgRef());
   
   Vec3 Wa(pNode1->GetWRef());
   Vec3 Wb(pNode2->GetWRef());
   
   Vec3 gPa(pNode1->GetgPRef());
   Vec3 gPb(pNode2->GetgPRef());
   
   // Aggiorna le deformazioni di riferimento nel sistema globale
   // Nota: non occorre G*g in quanto g/\g e' zero.
   TaCurrPrime = Mat3x3(MatG, ga)*gPa-Wa.Cross(ga*(4./(4.+ga.Dot())));
   TbCurrPrime = Mat3x3(MatG, gb)*gPb-Wa.Cross(gb*(4./(4.+gb.Dot())));

   // Calcola la deformazione corrente nel sistema locale (nodo a)  
   ThetaCurrPrime = ThetaRefPrime = RaT*(TbCurrPrime-TaCurrPrime); // +ThetaCurrPrime;
   
   // Aggiorna il legame costitutivo
   ConstitutiveLaw3DOwner::Update(0., ThetaRefPrime);
      
   // Chiede la matrice tangente di riferimento e la porta nel sistema globale
   FDEPrime = Ra*GetFDEPrime()*RaT;

   fFirstRes = flag(1);						     
}
						    
						     
/* assemblaggio jacobiano */
VariableSubMatrixHandler& 
ViscousHingeJoint::AssJac(VariableSubMatrixHandler& WorkMat,
			  doublereal dCoef,
			  const VectorHandler& /* XCurr */ ,
			  const VectorHandler& /* XPrimeCurr */ )
{
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeInit(iNumRows, iNumCols, 0.);

   /* Recupera gli indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex()+3;
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex()+3;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex()+3;
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex()+3;
   
   /* Setta gli indici della matrice */
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WM.fPutColIndex(iCnt, iNode1FirstPosIndex+iCnt);	
      WM.fPutRowIndex(3+iCnt, iNode2FirstMomIndex+iCnt);
      WM.fPutColIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
   }     

   Mat3x3 Ra(pNode1->GetRRef()*R1h);
   Vec3 Wa(pNode1->GetWRef());
 
   Mat3x3 Tmp(FDEPrime-FDEPrime*Mat3x3(Wa*dCoef));
         
   WM.Put(4, 4, Tmp);
   WM.Put(1, 4, -Tmp);

   Tmp += Mat3x3(Ra*(ConstitutiveLaw3DOwner::GetF()*dCoef));
   WM.Put(1, 1, Tmp);   
   WM.Put(4, 1, -Tmp);

   return WorkMat;
}


/* assemblaggio residuo */
SubVectorHandler& 
ViscousHingeJoint::AssRes(SubVectorHandler& WorkVec,
			  doublereal /* dCoef */ ,
			  const VectorHandler& /* XCurr */ ,
			  const VectorHandler& /* XPrimeCurr */ )
{
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
   

   /* Recupera gli indici */
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex()+3;
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex()+3;
   
   /* Setta gli indici della matrice */
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.fPutRowIndex(3+iCnt, iNode2FirstMomIndex+iCnt);
   }  

   Mat3x3 Ra(pNode1->GetRCurr()*R1h);
   
   if (fFirstRes) {      
      fFirstRes = flag(0);
   } else {
      Mat3x3 RaT(Ra.Transpose());
      
      Vec3 Wa(pNode1->GetWCurr());
      Vec3 Wb(pNode2->GetWCurr());
      
      Vec3 ga(pNode1->GetgCurr());
      Vec3 gb(pNode2->GetgCurr());
      
      Vec3 gPa(pNode1->GetgPCurr());
      Vec3 gPb(pNode2->GetgPCurr());

      TaCurrPrime = Mat3x3(MatG, ga)*gPa-Wa.Cross(ga*(4./(4.+ga.Dot())));
      TbCurrPrime = Mat3x3(MatG, gb)*gPb-Wa.Cross(gb*(4./(4.+gb.Dot())));
        
      ThetaCurrPrime = RaT*(TbCurrPrime-TaCurrPrime)+ThetaRefPrime;   
   
      ConstitutiveLaw3DOwner::Update(0., ThetaCurrPrime);
   }   
   
   Vec3 F(Ra*ConstitutiveLaw3DOwner::GetF());
   
   WorkVec.Put(1, F);
   WorkVec.Put(4, -F);
   
   return WorkVec;
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
ViscousHingeJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat, 
				 const VectorHandler& /* XCurr */ )
{
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeInit(iNumRows, iNumCols, 0.);

   /* Recupera gli indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex()+3;
   integer iNode1FirstVelIndex = iNode1FirstPosIndex+6;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex()+3;
   integer iNode2FirstVelIndex = iNode2FirstPosIndex+6;
   
   /* Setta gli indici della matrice */
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.fPutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.fPutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.fPutColIndex(3+iCnt, iNode1FirstVelIndex+iCnt);
      
      WM.fPutRowIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
      WM.fPutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
      WM.fPutColIndex(9+iCnt, iNode2FirstVelIndex+iCnt);
   }  

   Mat3x3 Ra(pNode1->GetRCurr()*R1h);
   
   Vec3 Wa(pNode1->GetWCurr());
   Vec3 Wb(pNode2->GetWCurr());
  
   FDEPrime = Ra*ConstitutiveLaw3DOwner::GetFDEPrime()*Ra.Transpose();
   
   Mat3x3 Tmp(FDEPrime*Mat3x3(Wb)+Mat3x3(Ra*ConstitutiveLaw3DOwner::GetF()));
   WM.Put(1, 1, Tmp);
   WM.Put(4, 1, -Tmp);
   
   Tmp = FDEPrime*Mat3x3(Wb);
   WM.Put(4, 7, Tmp);
   WM.Put(1, 7, -Tmp);
   
   WM.Put(1, 4, FDEPrime);
   WM.Put(4, 10, FDEPrime);
   
   Tmp = -FDEPrime;
   WM.Put(1, 10, Tmp);   
   WM.Put(4, 4, Tmp);

   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
ViscousHingeJoint::InitialAssRes(SubVectorHandler& WorkVec,
				 const VectorHandler& /* XCurr */ )
{
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
   

   /* Recupera gli indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex()+3;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex()+3;
   
   /* Setta gli indici della matrice */
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WorkVec.fPutRowIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
   }  

   Mat3x3 Ra(pNode1->GetRCurr()*R1h);
   
   if (fFirstRes) {      
      fFirstRes = flag(0);
   } else {
      Mat3x3 RaT(Ra.Transpose());
      
      Vec3 Wa(pNode1->GetWCurr());
      Vec3 Wb(pNode2->GetWCurr());
      
      Vec3 ga(pNode1->GetgCurr());
      Vec3 gb(pNode2->GetgCurr());
      
      TaCurrPrime = Wa-Wa.Cross(ga*(4./(4.+ga.Dot())));
      TbCurrPrime = Wb-Wa.Cross(gb*(4./(4.+gb.Dot())));
         
      ThetaCurrPrime = RaT*(TbCurrPrime-TaCurrPrime);
   
      ConstitutiveLaw3DOwner::Update(0., ThetaCurrPrime);
   }   
   
   Vec3 F(Ra*ConstitutiveLaw3DOwner::GetF());
   
   WorkVec.Put(1, F);
   WorkVec.Put(4, -F);   
      
   return WorkVec;
}

/* ViscousHingeJoint - end */


/* ViscoElasticHingeJoint - begin */

ViscoElasticHingeJoint::ViscoElasticHingeJoint(unsigned int uL, 
					       const DofOwner* pDO, 
					       const ConstitutiveLaw3D* pCL,
					       const StructNode* pN1, 
					       const StructNode* pN2,
					       const Mat3x3& R1,
					       const Mat3x3& R2, 
					       flag fOut)
: Elem(uL, Elem::JOINT, fOut), 
DeformableHingeJoint(uL, DefHingeType::VISCOELASTIC, 
		     pDO, pCL, pN1, pN2, R1, R2, fOut),
ThetaRef(0.), ThetaCurr(0.), ThetaRefPrime(0.), ThetaCurrPrime(0.),
TaCurr(0.), TbCurr(0.), TaCurrPrime(0.), TbCurrPrime(0.)
{
   Mat3x3 Ra(pNode1->GetRRef());
   Mat3x3 RaT(Ra.Transpose());
   FDE = Ra*ConstitutiveLaw3DOwner::GetFDE()*RaT;
   FDEPrime = Ra*ConstitutiveLaw3DOwner::GetFDEPrime()*RaT;   
}


ViscoElasticHingeJoint::~ViscoElasticHingeJoint(void)
{
   NO_OP;
}

   
void ViscoElasticHingeJoint::AfterPredict(VectorHandler& /* X */ ,
					  VectorHandler& /* XP */ )
{
   /* Calcola le deformazioni, aggiorna il legame costitutivo e crea la FDE */
   
   /* Recupera i dati */
   Mat3x3 Ra(pNode1->GetRRef()*R1h);
   Mat3x3 RaT(Ra.Transpose());
   
   Vec3 ga(pNode1->GetgRef());
   Vec3 gb(pNode2->GetgRef());
   
   Vec3 Wa(pNode1->GetWRef());
   Vec3 Wb(pNode2->GetWRef());
   
   Vec3 gPa(pNode1->GetgPRef());
   Vec3 gPb(pNode2->GetgPRef());
   
   /* Aggiorna le deformazioni di riferimento nel sistema globale
    * Nota: non occorre G*g in quanto g x g e' zero. */
   TaCurr = ga*(4./(4.+ga.Dot()));
   TbCurr = gb*(4./(4.+gb.Dot()));

   TaCurrPrime = Mat3x3(MatG, ga)*gPa-Wa.Cross(TaCurr);
   TbCurrPrime = Mat3x3(MatG, gb)*gPb-Wa.Cross(TbCurr);

   /* Calcola la deformazione corrente nel sistema locale (nodo a) */
   ThetaCurr = ThetaRef = RaT*(TbCurr-TaCurr)+ThetaCurr;
   ThetaCurrPrime = ThetaRefPrime = RaT*(TbCurrPrime-TaCurrPrime);
   
   /* Aggiorna il legame costitutivo */
   ConstitutiveLaw3DOwner::Update(ThetaRef, ThetaRefPrime);
      
   /* Chiede la matrice tangente di riferimento e la porta 
    * nel sistema globale */
   FDE = Ra*GetFDE()*RaT;
   FDEPrime = Ra*GetFDEPrime()*RaT;
						     
   fFirstRes = flag(1);
}


/* assemblaggio jacobiano */
VariableSubMatrixHandler& 
ViscoElasticHingeJoint::AssJac(VariableSubMatrixHandler& WorkMat,
			       doublereal dCoef, 
			       const VectorHandler& /* XCurr */ ,
			       const VectorHandler& /* XPrimeCurr */ )
{
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeInit(iNumRows, iNumCols, 0.);

   /* Recupera gli indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex()+3;
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex()+3;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex()+3;
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex()+3;
   
   /* Setta gli indici della matrice */
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WM.fPutColIndex(iCnt, iNode1FirstPosIndex+iCnt);	
      WM.fPutRowIndex(3+iCnt, iNode2FirstMomIndex+iCnt);
      WM.fPutColIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
   }  

   Mat3x3 Ra(pNode1->GetRRef()*R1h);
   Vec3 Wa(pNode1->GetWRef());
  
   Mat3x3 Tmp(FDE*dCoef+FDEPrime-FDEPrime*Mat3x3(Wa*dCoef));
         
   WM.Put(4, 4, Tmp);
   WM.Put(1, 4, -Tmp);
   
   Tmp += Mat3x3(Ra*(ConstitutiveLaw3DOwner::GetF()*dCoef));
   WM.Put(1, 1, Tmp);   
   WM.Put(4, 1, -Tmp);

   return WorkMat;
}


/* assemblaggio residuo */
SubVectorHandler& 
ViscoElasticHingeJoint::AssRes(SubVectorHandler& WorkVec,
			       doublereal /* dCoef */ ,
			       const VectorHandler& /* XCurr */ ,
			       const VectorHandler& /* XPrimeCurr */ )
{
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
   

   /* Recupera gli indici */
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex()+3;
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex()+3;
   
   /* Setta gli indici della matrice */
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.fPutRowIndex(3+iCnt, iNode2FirstMomIndex+iCnt);
   }  
   
   Mat3x3 Ra(pNode1->GetRCurr()*R1h);
   
   if (fFirstRes) {      
      fFirstRes = flag(0);
   } else {
      Mat3x3 RaT(Ra.Transpose());
      
      Vec3 Wa(pNode1->GetWCurr());
      Vec3 Wb(pNode2->GetWCurr());
      
      Vec3 ga(pNode1->GetgCurr());
      Vec3 gb(pNode2->GetgCurr());
      
      Vec3 gPa(pNode1->GetgPCurr());
      Vec3 gPb(pNode2->GetgPCurr());

      TaCurr = ga*(4./(4.+ga.Dot()));
      TbCurr = gb*(4./(4.+gb.Dot()));
      
      TaCurrPrime = Mat3x3(MatG, ga)*gPa-Wa.Cross(TaCurr);
      TbCurrPrime = Mat3x3(MatG, gb)*gPb-Wa.Cross(TbCurr);
      
      ThetaCurr = RaT*(TbCurr-TaCurr)+ThetaRef;
      ThetaCurrPrime = RaT*(TbCurrPrime-TaCurrPrime)+ThetaRefPrime;   
   
      ConstitutiveLaw3DOwner::Update(ThetaCurr, ThetaCurrPrime);
   }   
   
   Vec3 F(Ra*ConstitutiveLaw3DOwner::GetF());
   
   WorkVec.Put(1, F);
   WorkVec.Put(4, -F);   
   
   return WorkVec;
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
ViscoElasticHingeJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat, 
				      const VectorHandler& /* XCurr */ )
{
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeInit(iNumRows, iNumCols, 0.);

   /* Recupera gli indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex()+3;
   integer iNode1FirstVelIndex = iNode1FirstPosIndex+6;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex()+3;
   integer iNode2FirstVelIndex = iNode2FirstPosIndex+6;
   
   /* Setta gli indici della matrice */
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.fPutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.fPutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.fPutColIndex(3+iCnt, iNode1FirstVelIndex+iCnt);
      
      WM.fPutRowIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
      WM.fPutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
      WM.fPutColIndex(9+iCnt, iNode2FirstVelIndex+iCnt);
   }  

   Mat3x3 Ra(pNode1->GetRCurr()*R1h);
   Mat3x3 RaT(Ra.Transpose());
   Vec3 Wa(pNode1->GetWCurr());
   Vec3 Wb(pNode2->GetWCurr());
  
   FDE = Ra*ConstitutiveLaw3DOwner::GetFDE()*RaT;
   FDEPrime = Ra*ConstitutiveLaw3DOwner::GetFDEPrime()*RaT;      
   
   Mat3x3 Tmp(FDE-FDEPrime*Mat3x3(Wb)+Mat3x3(Ra*ConstitutiveLaw3DOwner::GetF()));
   WM.Put(1, 1, Tmp);
   WM.Put(4, 1, -Tmp);
   
   Tmp = FDE-FDEPrime*Mat3x3(Wb);
   WM.Put(4, 7, Tmp);
   WM.Put(1, 7, -Tmp);
   
   WM.Put(1, 4, FDEPrime);
   WM.Put(4, 10, FDEPrime);
   
   Tmp = -FDEPrime;
   WM.Put(1, 10, Tmp);   
   WM.Put(4, 4, Tmp);

   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
ViscoElasticHingeJoint::InitialAssRes(SubVectorHandler& WorkVec,
				      const VectorHandler& /* XCurr */ )
{
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
   

   /* Recupera gli indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex()+3;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex()+3;
   
   /* Setta gli indici della matrice */
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WorkVec.fPutRowIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
   }  
   
   Mat3x3 Ra(pNode1->GetRCurr()*R1h);
   
   if (fFirstRes) {      
      fFirstRes = flag(0);
   } else {
      Mat3x3 RaT(Ra.Transpose());
      
      Vec3 Wa(pNode1->GetWCurr());
      Vec3 Wb(pNode2->GetWCurr());
      
      Vec3 ga(pNode1->GetgCurr());
      Vec3 gb(pNode2->GetgCurr());
      
      TaCurr = ga*(4./(4.+ga.Dot()));
      TbCurr = gb*(4./(4.+gb.Dot()));
      
      TaCurrPrime = Wa-Wa.Cross(TaCurr);
      TbCurrPrime = Wb-Wa.Cross(TbCurr);
      
      ThetaCurr = RaT*(TbCurr-TaCurr)+ThetaRef;
      ThetaCurrPrime = RaT*(TbCurrPrime-TaCurrPrime);
   
      ConstitutiveLaw3DOwner::Update(ThetaCurr, ThetaCurrPrime);
   }   
   
   Vec3 F(Ra*ConstitutiveLaw3DOwner::GetF());
   
   WorkVec.Put(1, F);
   WorkVec.Put(4, -F);   
   
   return WorkVec;
}

/* ViscoElasticHingeJoint - end */
