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

/* Cerniera deformabile */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <vehj2.h>

/* DeformableDispHingeJoint - begin */

/* Costruttore non banale */
DeformableDispHingeJoint::DeformableDispHingeJoint(unsigned int uL,
						   DefHingeType::Type T,
						   const DofOwner* pDO, 
						   const ConstitutiveLaw3D* pCL,
						   const StructNode* pN1,
						   const StructNode* pN2,
						   const Vec3& f1Tmp,
						   const Vec3& f2Tmp,
						   const Mat3x3& R1,
						   const Mat3x3& R2,
						   flag fOut)
: Elem(uL, ElemType::JOINT, fOut), 
Joint(uL, JointType::DEFORMABLEHINGE, pDO, fOut), 
ConstitutiveLaw3DOwner(pCL), DefHingeT(T),
pNode1(pN1), pNode2(pN2), f1(f1Tmp), f2(f2Tmp), R1h(R1), R2h(R2)
{
   ASSERT(pNode1 != NULL);
   ASSERT(pNode2 != NULL);
   ASSERT(pNode1->GetNodeType() == NodeType::STRUCTURAL);
   ASSERT(pNode2->GetNodeType() == NodeType::STRUCTURAL);
}

   
/* Distruttore */
DeformableDispHingeJoint::~DeformableDispHingeJoint(void)
{
   NO_OP;
}

   
/* Contributo al file di restart */
ostream& DeformableDispHingeJoint::Restart(ostream& out) const
{
   Joint::Restart(out) << ", deformable displacement hinge, " 
     << pNode1->GetLabel() << ", reference, node, ",
     f1.Write(out, ", ") << ", hinge, reference, node, 1, ", 
     (R1h.GetVec(1)).Write(out, ", ")
     << ", 2, ", (R1h.GetVec(2)).Write(out, ", ") << ", "
     << pNode2->GetLabel() << ", reference, node, ",
     f2.Write(out, ", ") << ", hinge, reference, node, 1, ",
     (R2h.GetVec(1)).Write(out, ", ")
       << ", 2, ", (R2h.GetVec(2)).Write(out, ", ") << ", ";
   return pGetConstLaw()->Restart(out) << ';' << endl;
}


void DeformableDispHingeJoint::Output(OutputHandler& OH) const
{
   
   if(fToBeOutput()) {
      // Vec3 d(EulerAngles(pNode1->GetRCurr().Transpose()*pNode2->GetRCurr()));
      Vec3 v(GetF());
      Joint::Output(OH.Joints(), "DeformableHinge", GetLabel(),
		    v, Zero3, pNode1->GetRCurr()*v, Zero3) << endl;	
   }     
}

/* DeformableDispHingeJoint - end */


/* ElasticDispHingeJoint - begin */

ElasticDispHingeJoint::ElasticDispHingeJoint(unsigned int uL, 
					     const DofOwner* pDO, 
					     const ConstitutiveLaw3D* pCL,
					     const StructNode* pN1, 
					     const StructNode* pN2,
					     const Vec3& f1Tmp,
					     const Vec3& f2Tmp,
					     const Mat3x3& R1,
					     const Mat3x3& R2,
					     flag fOut)
: Elem(uL, ElemType::JOINT, fOut), 
DeformableDispHingeJoint(uL, DefHingeType::ELASTIC, 
			 pDO, pCL, pN1, pN2, 
			 f1Tmp, f2Tmp, R1, R2, fOut)
{
   NO_OP;
}


ElasticDispHingeJoint::~ElasticDispHingeJoint(void)
{
   NO_OP;
}

   
/* assemblaggio jacobiano */
VariableSubMatrixHandler& 
ElasticDispHingeJoint::AssJac(VariableSubMatrixHandler& WorkMat,
			      doublereal dCoef, 
			      const VectorHandler& /* XCurr */ ,
			      const VectorHandler& /* XPrimeCurr */ )
{
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->WorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeInit(iNumRows, iNumCols, 0.);

   /* Recupera gli indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   
   /* Setta gli indici della matrice */
   for(int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WM.fPutColIndex(iCnt, iNode1FirstPosIndex+iCnt);	
      WM.fPutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
      WM.fPutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
   }  
   
   AssMat(WM, dCoef);

   return WorkMat;
}


void ElasticDispHingeJoint::AssMat(FullSubMatrixHandler& WM, doublereal dCoef)   
{      
   Vec3 f1Tmp(pNode1->GetRRef()*f1);
   Vec3 f2Tmp(pNode2->GetRRef()*f2);
   Vec3 d1(pNode2->GetXCurr()+f2Tmp-pNode1->GetXCurr());

   Mat3x3 R1(pNode1->GetRRef()*R1h);
   Vec3 F(R1*(GetF()*dCoef));
   Mat3x3 MatF(F);
   Mat3x3 FDE(R1*GetFDE()*(R1.Transpose()*dCoef));
   
   Mat3x3 MTmp(FDE); /* F/e */
   WM.Put(1, 1, MTmp);
   WM.Put(7, 7, MTmp);
   MTmp = -MTmp; /* -F/e */
   WM.Put(1, 7, MTmp);
   WM.Put(7, 1, MTmp);
   
   MTmp = FDE*Mat3x3(f2Tmp); /* F/e * f2/\ */
   WM.Put(1, 10, MTmp);   
   WM.Put(7, 10, -MTmp);
   
   WM.Put(4, 10, Mat3x3(f1Tmp)*MTmp); /* f1/\ F/e f2/\ */
   
   MTmp = Mat3x3(f2Tmp)*FDE; /* f2/\ F/e */
   WM.Put(10, 7, MTmp);
   WM.Put(10, 1, -MTmp);
   
   WM.Put(10, 10, (MatF-MTmp)*Mat3x3(f2Tmp)); /* (F/\ - f2/\ F/e) f2/\ */
   
   MTmp = Mat3x3(f1Tmp)*FDE; /* f1/\ F/e */
   WM.Put(4, 1, MTmp);
   WM.Put(4, 7, -MTmp);
   
   MTmp = FDE*Mat3x3(d1)-MatF; /* F/e * d1/\ - F/\ */
   WM.Put(7, 4, MTmp);
   WM.Put(1, 4, -MTmp);
   WM.Put(10, 4, Mat3x3(f2Tmp)*MTmp); /* f2/\ (F/e * d1/\ - F/\) */
   
   WM.Put(4, 4, Mat3x3(f1Tmp, F)-Mat3x3(f1Tmp)*FDE*Mat3x3(d1)); /* ... */
}


/* assemblaggio residuo */
SubVectorHandler& 
ElasticDispHingeJoint::AssRes(SubVectorHandler& WorkVec,
			      doublereal /* dCoef */ ,
			      const VectorHandler& /* XCurr */ ,
			      const VectorHandler& /* XPrimeCurr */ )
{
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
   
   /* Recupera gli indici */
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   
   /* Setta gli indici della matrice */
   for(int iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.fPutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
   }  
   
   AssVec(WorkVec);
   
   return WorkVec;
}

   
void ElasticDispHingeJoint::AssVec(SubVectorHandler& WorkVec)
{   
   Mat3x3 R1(pNode1->GetRCurr()*R1h);
   Vec3 f1Tmp(pNode1->GetRCurr()*f1);
   Vec3 f2Tmp(pNode2->GetRCurr()*f2);
   Vec3 d1(pNode2->GetXCurr()+f2Tmp-pNode1->GetXCurr());   

   /* k = R1^T*(d1-f1) */   
   Vec3 k(R1.Transpose()*(d1-f1Tmp));
   ConstitutiveLaw3DOwner::Update(k);
   
   Vec3 F(R1*GetF());
   
   WorkVec.Put(1, F);
   WorkVec.Put(4, f1Tmp.Cross(F));
   WorkVec.Put(7, -F);   
   WorkVec.Put(10, F.Cross(f2Tmp));
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
ElasticDispHingeJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
				     const VectorHandler& /* XCurr */ )
{
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeInit(iNumRows, iNumCols, 0.);

   /* Recupera gli indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   
   /* Setta gli indici della matrice */
   for(int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.fPutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.fPutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.fPutRowIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
      WM.fPutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
   }  
   
   AssMat(WM, 1.);

   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
ElasticDispHingeJoint::InitialAssRes(SubVectorHandler& WorkVec,
				     const VectorHandler& /* XCurr */ )
{
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);   

   /* Recupera gli indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   
   /* Setta gli indici della matrice */
   for(int iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WorkVec.fPutRowIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
   }  
 
   AssVec(WorkVec);
   
   return WorkVec;
}

/* ElasticDispHingeJoint - end */


/* ViscousDispHingeJoint - begin */

ViscousDispHingeJoint::ViscousDispHingeJoint(unsigned int uL, 
					     const DofOwner* pDO, 
					     const ConstitutiveLaw3D* pCL,
					     const StructNode* pN1, 
					     const StructNode* pN2,
					     const Vec3& f1Tmp,
					     const Vec3& f2Tmp,
					     const Mat3x3& R1,
					     const Mat3x3& R2,
					     flag fOut)
: Elem(uL, ElemType::JOINT, fOut), 
DeformableDispHingeJoint(uL, DefHingeType::VISCOUS, 
			 pDO, pCL, pN1, pN2, f1Tmp, f2Tmp, R1, R2, fOut)
{
   NO_OP;
   
   /* Temporary */
   cerr << "DeformableHingeJoint(): warning, this element is not implemented yet" << endl;
   THROW(ErrNotImplementedYet());
}


ViscousDispHingeJoint::~ViscousDispHingeJoint(void)
{
   NO_OP;
}

   
/* assemblaggio jacobiano */
VariableSubMatrixHandler& 
ViscousDispHingeJoint::AssJac(VariableSubMatrixHandler& WorkMat,
			      doublereal dCoef,
			      const VectorHandler& /* XCurr */ ,
			      const VectorHandler& /* XPrimeCurr */ )
{
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->WorkSpaceDim(&iNumRows, &iNumCols);
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

   Mat3x3 R1(pNode1->GetRRef()*R1h);
   Vec3 Omega2(pNode2->GetWRef());

   Vec3 F(R1*GetF());
   Mat3x3 FDEPrime(R1*GetFDEPrime()*R1.Transpose());
   Mat3x3 Tmp(FDEPrime-FDEPrime*Mat3x3(Omega2*dCoef));
   
   WM.Put(4, 4, Tmp);
   WM.Put(1, 4, -Tmp);
   
   Tmp += Mat3x3(F*dCoef);
   WM.Put(1, 1, Tmp);   
   WM.Put(4, 1, -Tmp);

   return WorkMat;
}


/* assemblaggio residuo */
SubVectorHandler& 
ViscousDispHingeJoint::AssRes(SubVectorHandler& WorkVec,
			      doublereal /* dCoef */ ,
			      const VectorHandler& /* XCurr */ ,
			      const VectorHandler& /* XPrimeCurr */ )
{
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->WorkSpaceDim(&iNumRows, &iNumCols);
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
 
   Mat3x3 R1(pNode1->GetRCurr()*R1h);
   Vec3 Omega1(pNode1->GetWCurr());
   Vec3 Omega2(pNode2->GetWCurr());
   
   Vec3 g1(pNode1->GetgCurr());
   Vec3 g2(pNode2->GetgCurr());
   
   Vec3 kPrime(R1.Transpose()*(Omega2-Omega1));
   IncrementalUpdate(Vec3(0.), kPrime);
   
   Vec3 F(R1*GetF());
   
   WorkVec.Put(1, F);
   WorkVec.Put(4, -F);   
   
   return WorkVec;
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
ViscousDispHingeJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat, 
				     const VectorHandler& /* XCurr */ )
{
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->InitialWorkSpaceDim(&iNumRows, &iNumCols);
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

   Mat3x3 R1(pNode1->GetRRef()*R1h);
   Vec3 Omega1(pNode1->GetWRef());
   Vec3 Omega2(pNode2->GetWRef());

   Vec3 F(R1*GetF());
   Mat3x3 FDEPrime(R1*GetFDEPrime()*R1.Transpose());      
   Mat3x3 Tmp(Mat3x3(F)-FDEPrime*Mat3x3(Omega2-Omega1));
   
   WM.Put(1, 1, Tmp);
   WM.Put(4, 1, -Tmp);
   
   WM.Put(1, 4, FDEPrime);
   WM.Put(4, 10, FDEPrime);
   
   FDEPrime = -FDEPrime;
   WM.Put(1, 10, FDEPrime);   
   WM.Put(4, 4, FDEPrime);

   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
ViscousDispHingeJoint::InitialAssRes(SubVectorHandler& WorkVec,
				     const VectorHandler& /* XCurr */ )
{
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->InitialWorkSpaceDim(&iNumRows, &iNumCols);
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
   
   Mat3x3 R1(pNode1->GetRCurr()*R1h);
   Vec3 Omega1(pNode1->GetWCurr());
   Vec3 Omega2(pNode2->GetWCurr());
   
   Vec3 g1(pNode1->GetgCurr());
   Vec3 g2(pNode2->GetgCurr());
   
   /* Aggiornamento: k += R1^T(G(g2)*g2-G(g1)*g1) */
   Vec3 k(R1.Transpose()*(g2*(4./(4.+g2.Dot()))
			  -g1*(4./(4.+g1.Dot()))));
   Vec3 kPrime(R1.Transpose()*(Omega2-Omega1));
   IncrementalUpdate(k, kPrime);
   
   Vec3 F(R1*GetF());
   
   WorkVec.Put(1, F);
   WorkVec.Put(4, -F);   
   
   return WorkVec;
}

/* ViscousDispHingeJoint - end */


/* ViscoElasticDispHingeJoint - begin */

ViscoElasticDispHingeJoint::ViscoElasticDispHingeJoint(unsigned int uL, 
						       const DofOwner* pDO, 
						       const ConstitutiveLaw3D* pCL,
						       const StructNode* pN1, 
						       const StructNode* pN2,
						       const Vec3& f1Tmp,
						       const Vec3& f2Tmp,
						       const Mat3x3& R1,
						       const Mat3x3& R2, 
						       flag fOut)
: Elem(uL, ElemType::JOINT, fOut), 
DeformableDispHingeJoint(uL, DefHingeType::VISCOELASTIC, 
			 pDO, pCL, pN1, pN2, f1Tmp, f2Tmp, R1, R2, fOut)
{
   /* Temporary */
   cerr << "DeformableHingeJoint(): warning, this element is not implemented yet" << endl;
   THROW(ErrNotImplementedYet());
}


ViscoElasticDispHingeJoint::~ViscoElasticDispHingeJoint(void)
{
   NO_OP;
}

   
/* assemblaggio jacobiano */
VariableSubMatrixHandler& 
ViscoElasticDispHingeJoint::AssJac(VariableSubMatrixHandler& WorkMat,
				   doublereal dCoef, 
				   const VectorHandler& /* XCurr */ ,
				   const VectorHandler& /* XPrimeCurr */ )
{
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->WorkSpaceDim(&iNumRows, &iNumCols);
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

   Mat3x3 R1(pNode1->GetRRef()*R1h);
   Mat3x3 R1t(R1.Transpose());
   Vec3 Omega2(pNode2->GetWRef());

   Vec3 F(R1*GetF());
   Mat3x3 FDE(R1*GetFDE()*(R1t*dCoef));
   Mat3x3 FDEPrime(R1*GetFDEPrime()*R1t);
   
   Mat3x3 Tmp(FDEPrime-FDEPrime*Mat3x3(Omega2*dCoef)+FDE);
         
   WM.Put(4, 4, Tmp);
   WM.Put(1, 4, -Tmp);
   
   Tmp += Mat3x3(F*dCoef);
   WM.Put(1, 1, Tmp);   
   WM.Put(4, 1, -Tmp);

   return WorkMat;
}


/* assemblaggio residuo */
SubVectorHandler& 
ViscoElasticDispHingeJoint::AssRes(SubVectorHandler& WorkVec,
				   doublereal /* dCoef */ ,
				   const VectorHandler& /* XCurr */ ,
				   const VectorHandler& /* XPrimeCurr */ )
{
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->WorkSpaceDim(&iNumRows, &iNumCols);
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
   
   Mat3x3 R1(pNode1->GetRCurr()*R1h);
   Mat3x3 R1t(R1.Transpose());
   Vec3 Omega1(pNode1->GetWCurr());
   Vec3 Omega2(pNode2->GetWCurr());
   
   Vec3 g1(pNode1->GetgCurr());
   Vec3 g2(pNode2->GetgCurr());
   
   /* Aggiornamento: k += R1^T(G(g2)*g2-G(g1)*g1) */
   Vec3 k(R1t*(g2*(4./(4.+g2.Dot()))
	       -g1*(4./(4.+g1.Dot()))));
   Vec3 kPrime(R1t*(Omega2-Omega1));
   IncrementalUpdate(k, kPrime);
   
   Vec3 F(R1*GetF());
   
   WorkVec.Put(1, F);
   WorkVec.Put(4, -F);   
   
   return WorkVec;
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
ViscoElasticDispHingeJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat, 
					  const VectorHandler& /* XCurr */ )
{
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->InitialWorkSpaceDim(&iNumRows, &iNumCols);
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
   
   Mat3x3 R1(pNode1->GetRRef()*R1h);
   Mat3x3 R1t(R1.Transpose());
   Vec3 Omega1(pNode1->GetWRef());
   Vec3 Omega2(pNode2->GetWRef());

   Vec3 F(R1*GetF());
   Mat3x3 FDE(R1*GetFDE()*R1t);
   Mat3x3 FDEPrime(R1*GetFDEPrime()*R1t);
   
   Mat3x3 Tmp(Mat3x3(F)-FDEPrime*Mat3x3(Omega2-Omega1)+FDE);
   
   WM.Put(1, 1, Tmp);
   WM.Put(4, 1, -Tmp);
   
   WM.Put(4, 7, FDE);
   WM.Put(1, 7, -FDE);
   
   WM.Put(1, 4, FDEPrime);
   WM.Put(4, 10, FDEPrime);
   
   FDEPrime = -FDEPrime;
   WM.Put(1, 10, FDEPrime);   
   WM.Put(4, 4, FDEPrime);

   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
ViscoElasticDispHingeJoint::InitialAssRes(SubVectorHandler& WorkVec,
					  const VectorHandler& /* XCurr */ )
{
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->InitialWorkSpaceDim(&iNumRows, &iNumCols);
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
   
   Mat3x3 R1(pNode1->GetRCurr()*R1h);
   Mat3x3 R1t(R1.Transpose());
   Vec3 Omega1(pNode1->GetWCurr());
   Vec3 Omega2(pNode2->GetWCurr());
   
   Vec3 g1(pNode1->GetgCurr());
   Vec3 g2(pNode2->GetgCurr());
   
   /* Aggiornamento: k += R1^T(G(g2)*g2-G(g1)*g1) */
   Vec3 k(R1t*(g2*(4./(4.+g2.Dot()))
	       -g1*(4./(4.+g1.Dot()))));
   Vec3 kPrime(R1t*(Omega2-Omega1));
   IncrementalUpdate(k, kPrime);
   
   Vec3 F(R1*GetF());
   
   WorkVec.Put(1, F);
   WorkVec.Put(4, -F);   
   
   return WorkVec;
}

/* ViscoElasticDispHingeJoint - end */
