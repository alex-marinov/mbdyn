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

/* Rods */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <rodj.h>

extern "C" {
#include <mymath.h>
}

/* RodJoint - begin */

/* Costruttore non banale */
RodJoint::RodJoint(unsigned int uL, const DofOwner* pDO, 
		   const ConstitutiveLaw1D* pCL,
		   const StructNode* pN1, const StructNode* pN2,
		   doublereal dLength, flag fOut, flag fHasOffsets)
: Elem(uL, ElemType::JOINT, fOut), 
Joint(uL, JointType::ROD, pDO, fOut),
ConstitutiveLaw1DOwner(pCL), RodT(RodType::ELASTIC),
pNode1(pN1), pNode2(pN2), dL0(dLength), v(0.), dElle(0.)
{
#ifdef DEBUG   
   /* Verifica di consistenza dei dati iniziali */   
   ASSERT(pN1 != NULL);
   ASSERT(pN1->GetNodeType() == NodeType::STRUCTURAL);
   ASSERT(pN2 != NULL);
   ASSERT(pN2->GetNodeType() == NodeType::STRUCTURAL);
   
   if (!fHasOffsets) {
      v = pN2->GetXCurr()-pN1->GetXCurr();
   
      ASSERT(v.Dot() > DBL_EPSILON);
   }
   
   ASSERT(dLength > DBL_EPSILON);
#endif	
}


/* Distruttore */
RodJoint::~RodJoint(void) 
{
   NO_OP;
}


/* Contributo al file di restart */
ostream& RodJoint::Restart(ostream& out) const
{
   Joint::Restart(out) << ", rod, "
     << pNode1->GetLabel() << ", "
     << pNode2->GetLabel() << ", " 
     << dL0 << ", ";
   return pGetConstLaw()->Restart(out) << ';' << endl;
}


void RodJoint::AssMat(FullSubMatrixHandler& WorkMat, doublereal dCoef)
{
   /* v = x2-x1 */
   /* v = pNode2->GetXCurr()-pNode1->GetXCurr(); */
   doublereal dCross = v.Dot();
   
   /* Verifica che la distanza non sia nulla */
   if (dCross <= DBL_EPSILON) {
      cerr << endl << "Null distance between nodes " << pNode1->GetLabel() 
	<< " and " << pNode2->GetLabel() 
	<< " in Rod Joint " << uLabel << ';' << endl
	<< "aborting ..." << endl;
      THROW(Joint::ErrGeneric());
   }   

   /* Lunghezza corrente */
   dElle = sqrt(dCross);

   /* Forza e slope */
   doublereal dF = GetF();
   doublereal dFDE = GetFDE();

   Mat3x3 K(Mat3x3(v, v*((-dF*dCoef)/(dElle*dCross)))
    	    +v.Tens(v*((dFDE*dCoef)/(dL0*dCross))));
         
   /* Termini diagonali */
   WorkMat.Add(1, 1, K);
   WorkMat.Add(4, 4, K);
   
   /* termini extradiagonali */
   K = -K;
   WorkMat.Add(1, 4, K);
   WorkMat.Add(4, 1, K);
}


void RodJoint::AssVec(SubVectorHandler& WorkVec)
{
   /* v = x2-x1 */
   v = pNode2->GetXCurr()-pNode1->GetXCurr();
   doublereal dCross = v.Dot();
   
   /* Verifica che la distanza non sia nulla */
   if (dCross <= DBL_EPSILON) {
      cerr << endl << "Null distance between nodes " << pNode1->GetLabel() 
	<< " and " << pNode2->GetLabel() 
	<< " in Rod Joint " << uLabel << ';' << endl
	<< "aborting ..." << endl;
      THROW(Joint::ErrGeneric());
   }   
   
   /* Deformazione */
   dElle = sqrt(dCross);
   doublereal dEpsilon = dElle/dL0-1.;
   
   /* Ampiezza della forza */
   ConstitutiveLaw1DOwner::Update(dEpsilon);
   doublereal dF = GetF();
   
   /* Vettore forza */
   Vec3 F = v*(dF/dElle);
   
   WorkVec.Add(1, F);
   WorkVec.Add(4, -F);
}


VariableSubMatrixHandler& RodJoint::AssJac(VariableSubMatrixHandler& WorkMat,
					   doublereal dCoef,
					   const VectorHandler& /* XCurr */ ,
					   const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering RodJoint::AssJac()" << endl);
   
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
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WM.fPutColIndex(iCnt, iNode1FirstPosIndex+iCnt);	
      WM.fPutRowIndex(3+iCnt, iNode2FirstMomIndex+iCnt);
      WM.fPutColIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
   }  
   
   /* Genera la matrice */
   AssMat(WM, dCoef);
   
   return WorkMat;
}


void RodJoint::AssEig(VariableSubMatrixHandler& WorkMatA,
		      VariableSubMatrixHandler& WorkMatB,
		      const VectorHandler& /* XCurr */ ,
		      const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering RodJoint::AssEig()" << endl);
   
   WorkMatB.SetNullMatrix();   
   FullSubMatrixHandler& WMA = WorkMatA.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->WorkSpaceDim(&iNumRows, &iNumCols);
   WMA.ResizeInit(iNumRows, iNumCols, 0.);

   /* Recupera gli indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   
   /* Setta gli indici della matrice */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WMA.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WMA.fPutColIndex(iCnt, iNode1FirstPosIndex+iCnt);	
      WMA.fPutRowIndex(3+iCnt, iNode2FirstMomIndex+iCnt);
      WMA.fPutColIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
   }  
   
   /* Genera la matrice */
   AssMat(WMA, 1.);   
}


SubVectorHandler& RodJoint::AssRes(SubVectorHandler& WorkVec,
				   doublereal /* dCoef */ ,
				   const VectorHandler& /* XCurr */ ,
				   const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering RodJoint::AssRes()" << endl);
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
   
   /* Indici */
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   
   /* Setta gli indici */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WorkVec.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.fPutRowIndex(3+iCnt, iNode2FirstMomIndex+iCnt);
   }
   
   /* Costruisce il vettore */
   AssVec(WorkVec);
   
   return WorkVec;   
}
   

void RodJoint::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) {      
#ifdef DEBUG   
      OH.Output() << "Joint " << uLabel << ", type \""
	<< psJointNames[JointType::ROD]
	<< "\", linked to nodes " 
	<< pNode1->GetLabel() << " and " << pNode2->GetLabel() << ':' << endl
	<< "Initial length = " << dL0 << endl;
#endif   
      
      ASSERT(dElle > DBL_EPSILON);	
      Vec3 vTmp(v/dElle);      
      doublereal d = GetF();
      
      Joint::Output(OH.Joints(), "Rod", GetLabel(),
		    Vec3(d, 0., 0.), Zero3, vTmp*d, Zero3)
	<< " " << dElle << " " << vTmp << endl;
   }
}
 

VariableSubMatrixHandler& 
RodJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
			const VectorHandler& /* XCurr */ )
{
   DEBUGCOUT("Entering RodJoint::InitialAssJac()" << endl);
   
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
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.fPutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.fPutColIndex(iCnt, iNode1FirstPosIndex+iCnt);	
      WM.fPutRowIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
      WM.fPutColIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
   }  
   
   /* Genera la matrice */
   AssMat(WM);
   
   return WorkMat;
}


SubVectorHandler& RodJoint::InitialAssRes(SubVectorHandler& WorkVec,
					  const VectorHandler& /* XCurr */ )
{
   DEBUGCOUT("Entering RodJoint::InitialAssRes()" << endl);
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
   
   /* Indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   
   /* Setta gli indici */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WorkVec.fPutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WorkVec.fPutRowIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
   }
   
   /* Costruisce il vettore */
   AssVec(WorkVec);
   
   return WorkVec;   
}


void 
RodJoint::GetAdamsDummyPart(unsigned int part,
			    Vec3& x, 
			    Mat3x3& R) const 
{
   ASSERT(part == 1);
   x = pNode1->GetXCurr();
   R = pNode1->GetRCurr();
}


ostream& 
RodJoint::WriteAdamsDummyPartCmd(ostream& out,
				 unsigned int part, 
				 unsigned int firstId) const
{
   Vec3 x1 = pNode1->GetXCurr();
   Vec3 x2 = pNode2->GetXCurr();
     
   Vec3 v1 = x2-x1;
   doublereal l = v1.Norm();
   v1 /= l;
   
   Mat3x3 Rx(Eye3-v1.Tens(v1));
   int index = 1;
   if (fabs(v1.dGet(2)) < fabs(v1.dGet(index))) {
      index = 2;
   }
   if (fabs(v1.dGet(3)) < fabs(v1.dGet(index))) {
      index = 3;
   }
   
   Vec3 v2(Rx.GetVec(index));
   v2 /= v2.Norm();
   
   Vec3 e(EulerAngles(MatR2vec(1, v1, 2, v2)));
   
   return out 
     << psAdamsElemCode[GetElemType()] << "_" << GetLabel() << "_" << part << endl
     << firstId << " "
     << x1 << " "
     << EulerAngles(pNode1->GetRCurr()) << " "
     << x1 << " "
     << e << " "
     << l << " " << 0. << " " << 0. << " "
     << Zero3 << endl;
}

/* RodJoint - end */


/* ViscoElasticRodJoint - begin */

/* Costruttore non banale */
ViscoElasticRodJoint::ViscoElasticRodJoint(unsigned int uL, 
					   const DofOwner* pDO,
					   const ConstitutiveLaw1D* pCL,
					   const StructNode* pN1, 
					   const StructNode* pN2,
					   doublereal dLength, flag fOut)
: Elem(uL, ElemType::JOINT, fOut), 
RodJoint(uL, pDO, pCL, pN1, pN2, dLength, fOut)
{
   SetRodType(RodType::VISCOELASTIC);
}


/* Distruttore */
ViscoElasticRodJoint::~ViscoElasticRodJoint(void) 
{ 
   NO_OP; 
}


VariableSubMatrixHandler& 
ViscoElasticRodJoint::AssJac(VariableSubMatrixHandler& WorkMat,
			     doublereal dCoef,
			     const VectorHandler& /* XCurr */ ,
			     const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering ViscoElasticRodJoint::AssJac()" << endl);
   
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
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WM.fPutColIndex(iCnt, iNode1FirstPosIndex+iCnt);	
      WM.fPutRowIndex(3+iCnt, iNode2FirstMomIndex+iCnt);
      WM.fPutColIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
   }           
   
   /* v = x2-x1 */
   /* v(pNode2->GetXCurr()-pNode1->GetXCurr()); */
   doublereal dCross = v.Dot();
   
   /* Verifica che la distanza non sia nulla */
   if (dCross <= DBL_EPSILON) {
      cerr << endl << "Null distance between nodes " << pNode1->GetLabel() 
	<< " and " << pNode2->GetLabel() 
	<< " in Rod Joint " << uLabel << ';' << endl
	<< "aborting ..." << endl;
      THROW(Joint::ErrGeneric());
   }   

   /* Lunghezza corrente */
   dElle = sqrt(dCross);

   /* Velocita' di deformazione */
   Vec3 vPrime(pNode2->GetVCurr()-pNode1->GetVCurr());

   /* Forza e slopes */
   doublereal dF = GetF();
   doublereal dFDE = GetFDE();
   doublereal dFDEPrime = GetFDEPrime();
   
   Mat3x3 K(Mat3x3( v, v*((-dF*dCoef)/(dElle*dCross)) )
	    + v.Tens( v*((dFDE*dCoef+dFDEPrime)/(dL0*dCross)) )
	    + v.Tens( v.Cross( vPrime.Cross( v*((dFDEPrime*dCoef)/(dL0*dCross*dCross)) ) ) ));
         
   /* Termini diagonali */
   WM.Add(1, 1, K);
   WM.Add(4, 4, K);
   
   /* termini extradiagonali */
   K = -K;
   WM.Add(1, 4, K);
   WM.Add(4, 1, K);
   
   return WorkMat;
}


SubVectorHandler& 
ViscoElasticRodJoint::AssRes(SubVectorHandler& WorkVec,
			     doublereal /* dCoef */ ,
			     const VectorHandler& /* XCurr */ ,
			     const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering ViscoElasticRodJoint::AssRes()" << endl);
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
   
   /* Indici */
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   
   /* Setta gli indici */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WorkVec.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.fPutRowIndex(3+iCnt, iNode2FirstMomIndex+iCnt);
   }
         
   /* v = x2-x1 */
   v = pNode2->GetXCurr()-pNode1->GetXCurr();
   doublereal dCross = v.Dot();
   
   /* Verifica che la distanza non sia nulla */
   if (dCross <= DBL_EPSILON) {
      cerr << endl << "Null distance between nodes " << pNode1->GetLabel() 
	<< " and " << pNode2->GetLabel() 
	<< " in Rod Joint " << uLabel << ';' << endl
	<< "aborting ..." << endl;
      THROW(Joint::ErrGeneric());
   }   

   /* Lunghezza corrente */
   dElle = sqrt(dCross);

   /* Deformazione */
   doublereal dEpsilon = dElle/dL0-1.;
   
   /* Velocita' di deformazione */
   Vec3 vPrime(pNode2->GetVCurr()-pNode1->GetVCurr());
   doublereal dEpsilonPrime  = (v.Dot(vPrime))/(dElle*dL0);
   
   /* Ampiezza della forza */
   ConstitutiveLaw1DOwner::Update(dEpsilon, dEpsilonPrime);
   doublereal dF = GetF();
   
   /* Vettore forza */
   Vec3 F(v*(dF/dElle));
   
   WorkVec.Add(1, F);
   WorkVec.Add(4, -F);
   
   return WorkVec;
}


VariableSubMatrixHandler& 
ViscoElasticRodJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
				    const VectorHandler& /* XCurr */ )
{
   DEBUGCOUT("Entering ViscoElasticRodJoint::InitialAssJac()" << endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeInit(iNumRows, iNumCols, 0.);

   /* Recupera gli indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstVelIndex = iNode1FirstPosIndex+6;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstVelIndex = iNode2FirstPosIndex+6;
   
   /* Setta gli indici della matrice */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.fPutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.fPutColIndex(iCnt, iNode1FirstPosIndex+iCnt);	
      WM.fPutColIndex(3+iCnt, iNode1FirstVelIndex+iCnt);
      
      WM.fPutRowIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
      WM.fPutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
      WM.fPutColIndex(9+iCnt, iNode2FirstVelIndex+iCnt);
   }           
   
   /* v = x2-x1 */
   /* v(pNode2->GetXCurr()-pNode1->GetXCurr()); */
   doublereal dCross = v.Dot();
   
   /* Verifica che la distanza non sia nulla */
   if (dCross <= DBL_EPSILON) {
      cerr << endl << "Null distance between nodes " << pNode1->GetLabel() 
	<< " and " << pNode2->GetLabel() 
	<< " in Rod Joint " << uLabel << ';' << endl
	<< "aborting ..." << endl;
      THROW(Joint::ErrGeneric());
   }   

   /* Lunghezza corrente */
   dElle = sqrt(dCross);

   /* Velocita' di deformazione */
   Vec3 vPrime(pNode2->GetVCurr()-pNode1->GetVCurr());

   /* Forza e slopes */
   doublereal dF = GetF();
   doublereal dFDE = GetFDE();
   doublereal dFDEPrime = GetFDEPrime();

   Mat3x3 K( Mat3x3( v, v*((-dF)/(dElle*dCross)) )
	    + v.Tens( v*((dFDE)/(dL0*dCross)) )
	    + v.Tens( v.Cross( vPrime.Cross( v*((dFDEPrime)/(dL0*dCross*dCross)) ) ) ));
   Mat3x3 KPrime( v.Tens( v*((dFDEPrime)/(dL0*dCross)) ) );
         
   /* Termini diagonali */
   WM.Add(1, 1, K);
   WM.Add(4, 7, K);
   
   WM.Add(1, 4, KPrime);
   WM.Add(4, 10, KPrime);
   
   /* termini extradiagonali */
   K = -K;
   WM.Add(1, 7, K);
   WM.Add(4, 1, K);

   KPrime = -KPrime;
   WM.Add(1, 10, KPrime);
   WM.Add(4, 4, KPrime);
   
   return WorkMat;
}


SubVectorHandler& 
ViscoElasticRodJoint::InitialAssRes(SubVectorHandler& WorkVec,
				    const VectorHandler& /* XCurr */ )
{
   DEBUGCOUT("Entering ViscoElasticRodJoint::InitialAssRes()" << endl);
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
   
   /* Indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   
   /* Setta gli indici */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WorkVec.fPutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WorkVec.fPutRowIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
   }
         
   /* v = x2-x1 */
   v = pNode2->GetXCurr()-pNode1->GetXCurr();
   doublereal dCross = v.Dot();
   
   /* Verifica che la distanza non sia nulla */
   if (dCross <= DBL_EPSILON) {
      cerr << endl << "Null distance between nodes " << pNode1->GetLabel() 
	<< " and " << pNode2->GetLabel() 
	<< " in Rod Joint " << uLabel << ';' << endl
	<< "aborting ..." << endl;
      THROW(Joint::ErrGeneric());
   }   

   /* Lunghezza corrente */
   dElle = sqrt(dCross);

   /* Deformazione */
   doublereal dEpsilon = dElle/dL0-1.;
   
   /* Velocita' di deformazione */
   Vec3 vPrime(pNode2->GetVCurr()-pNode1->GetVCurr());
   doublereal dEpsilonPrime  = (v.Dot(vPrime))/(dElle*dL0);
   
   /* Ampiezza della forza */
   ConstitutiveLaw1DOwner::Update(dEpsilon, dEpsilonPrime);
   doublereal dF = GetF();
   
   /* Vettore forza */
   Vec3 F(v*(dF/dElle));
   
   WorkVec.Add(1, F);
   WorkVec.Add(4, -F);
   
   return WorkVec;
}

/* ViscoElasticRodJoint - end */


/* RodWithOffsetJoint - begin */

/* Costruttore non banale */
RodWithOffsetJoint::RodWithOffsetJoint(unsigned int uL, 
				       const DofOwner* pDO,
				       const ConstitutiveLaw1D* pCL,
				       const StructNode* pN1, 
				       const StructNode* pN2,
				       const Vec3& f1Tmp, 
				       const Vec3& f2Tmp,
				       doublereal dLength, 
				       flag fOut)
: Elem(uL, ElemType::JOINT, fOut), 
RodJoint(uL, pDO, pCL, pN1, pN2, dLength, fOut, 1),
f1(f1Tmp), f2(f2Tmp)
{
#ifdef DEBUG   
   /* Verifica di consistenza dei dati iniziali */   
   ASSERT(pN1 != NULL);
   ASSERT(pN1->GetNodeType() == NodeType::STRUCTURAL);
   ASSERT(pN2 != NULL);
   ASSERT(pN2->GetNodeType() == NodeType::STRUCTURAL);
   
   v = pN2->GetXCurr()+(pN2->GetRCurr()*f2Tmp)
     -pN1->GetXCurr()-(pN1->GetRCurr()*f1Tmp);
   
   ASSERT(v.Dot() > DBL_EPSILON);
   ASSERT(dLength > DBL_EPSILON);
#endif	   
   
   SetRodType(RodType::VISCOELASTICWITHOFFSET);
}


/* Distruttore */
RodWithOffsetJoint::~RodWithOffsetJoint(void)
{
   NO_OP;
}
      
   
/* Contributo al file di restart */
ostream& RodWithOffsetJoint::Restart(ostream& out) const
{
   Joint::Restart(out) << ", rod, "
     << pNode1->GetLabel() << ", "
     << pNode2->GetLabel() << ", " 
     << dL0 << ", offset, reference, node, ",
     f1.Write(out, ", ") << ", reference, node, ",
     f2.Write(out, ", ") << ", ";
   return pGetConstLaw()->Restart(out) << ';' << endl;
}

         
VariableSubMatrixHandler& 
RodWithOffsetJoint::AssJac(VariableSubMatrixHandler& WorkMat, 
			   doublereal dCoef,
			   const VectorHandler& /* XCurr */ ,
			   const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering RodWithOffsetJoint::AssJac()" << endl);
   
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
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WM.fPutColIndex(iCnt, iNode1FirstPosIndex+iCnt);	
      WM.fPutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
      WM.fPutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
   }           
   
   Mat3x3 R1(pNode1->GetRRef());
   Mat3x3 R2(pNode2->GetRRef());
   Vec3 f1Tmp(R1*f1);
   Vec3 f2Tmp(R2*f2);
   Vec3 x1(pNode1->GetXCurr());
   Vec3 x2(pNode2->GetXCurr());
   
   Vec3 v1(pNode1->GetVCurr());
   Vec3 v2(pNode2->GetVCurr());
   Vec3 Omega1(pNode1->GetWRef());
   Vec3 Omega2(pNode2->GetWRef());
   
   /* v = p2-p1 */
   /* v = x2+f2Tmp-x1-f1Tmp; */
   doublereal dCross = v.Dot();
   
   /* Verifica che la distanza non sia nulla */
   if (dCross <= DBL_EPSILON) {
      cerr << endl << "Null distance between nodes " << pNode1->GetLabel() 
	<< " and " << pNode2->GetLabel() 
	<< " in Rod Joint " << uLabel << ';' << endl
	<< "aborting ..." << endl;
      THROW(Joint::ErrGeneric());
   }   

   /* Lunghezza corrente */
   dElle = sqrt(dCross);

   /* Velocita' di deformazione */
   Vec3 vPrime(v2+Omega2.Cross(f2Tmp)-v1-Omega1.Cross(f1Tmp));

   /* Forza e slopes */        
   doublereal dF = GetF();
   doublereal dFDE = GetFDE();
   doublereal dFDEPrime = GetFDEPrime();

   /* Vettore forza */
   Vec3 F = v*(dF/dElle);
   
   Mat3x3 K(Mat3x3(v,v*((-dF*dCoef)/(dElle*dCross)))
	    +v.Tens(v*((dFDE*dCoef)/(dL0*dCross)))
	    +v.Tens(v.Cross(vPrime.Cross(v*((dFDEPrime*dCoef)/(dL0*dCross*dCross))))));

   Mat3x3 KPrime(v.Tens(v*((dFDEPrime)/(dL0*dCross))));
   
   /* Termini di forza diagonali */
   Mat3x3 Tmp(K+KPrime);   
   WM.Add(1, 1, Tmp);
   WM.Add(7, 7, Tmp);
   
   /* Termini di coppia, nodo 1 */
   Mat3x3 Tmp2 = Mat3x3(f1Tmp)*Tmp;
   WM.Add(4, 1, Tmp2);
   WM.Add(4, 7, -Tmp2);
   
   /* Termini di coppia, nodo 2 */
   Tmp2 = Mat3x3(f2Tmp)*Tmp;
   WM.Add(10, 7, Tmp2);
   WM.Add(10, 1, -Tmp2);
   
   /* termini di forza extradiagonali */
   Tmp = -Tmp;
   WM.Add(1, 7, Tmp);
   WM.Add(7, 1, Tmp);
   
   
   /* Termini di rotazione, Delta g1 */
   Tmp = K*Mat3x3(-f1Tmp)
     +KPrime*(Mat3x3(Omega1.Cross(f1Tmp*dCoef))-Mat3x3(f1Tmp));
   WM.Add(1, 4, Tmp);
   WM.Add(7, 4, -Tmp);
   
   /* Termini di coppia, Delta g1 */
   Tmp2 = Mat3x3(f1Tmp)*Tmp;
   WM.Add(10, 4, -Tmp2);
   Tmp2 += Mat3x3(f1Tmp*dCoef, F);
   WM.Add(4, 4, Tmp2);     
           
   /* Termini di rotazione, Delta g2 */
   Tmp = K*Mat3x3(-f2Tmp)
     +KPrime*(Mat3x3(Omega2.Cross(f2Tmp*dCoef))-Mat3x3(f2Tmp));   
   WM.Add(7, 10, Tmp);
   WM.Add(1, 10, -Tmp);

   /* Termini di coppia, Delta g2 */
   Tmp2 = Mat3x3(f2Tmp)*Tmp;
   WM.Add(4, 10, -Tmp2);
   Tmp2 += Mat3x3(f2Tmp*dCoef, F);
   WM.Add(10, 10, Tmp2);     
   
   return WorkMat;   
}

	   
SubVectorHandler& 
RodWithOffsetJoint::AssRes(SubVectorHandler& WorkVec,
			   doublereal /* dCoef */ ,
			   const VectorHandler& /* XCurr */ ,
			   const VectorHandler& /* XPrimeCurr */ )
{   
   DEBUGCOUT("RodWithOffsetJoint::AssRes()" << endl);
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
   
   /* Indici */
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   
   /* Setta gli indici */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {	
      WorkVec.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.fPutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
   }

   /* Dati */
   Mat3x3 R1(pNode1->GetRCurr());
   Mat3x3 R2(pNode2->GetRCurr());
   Vec3 f1Tmp(R1*f1);
   Vec3 f2Tmp(R2*f2);
   Vec3 x1(pNode1->GetXCurr());
   Vec3 x2(pNode2->GetXCurr());
   
   Vec3 v1(pNode1->GetVCurr());
   Vec3 v2(pNode2->GetVCurr());
   Vec3 Omega1(pNode1->GetWCurr());
   Vec3 Omega2(pNode2->GetWCurr());
   
   
   /* v = x2-x1 */
   v = x2+f2Tmp-x1-f1Tmp;
   doublereal dCross = v.Dot();
   
   /* Verifica che la distanza non sia nulla */
   if (dCross <= DBL_EPSILON) {
      cerr << endl << "Null distance between nodes " << pNode1->GetLabel() 
	<< " and " << pNode2->GetLabel() 
	<< " in Rod Joint " << uLabel << ';' << endl
	<< "aborting ..." << endl;
      THROW(Joint::ErrGeneric());
   }   

   /* Lunghezza corrente */
   dElle = sqrt(dCross);

   /* Deformazione */
   doublereal dEpsilon = dElle/dL0-1.;
   
   /* Velocita' di deformazione */
   Vec3 vPrime(v2+Omega2.Cross(f2Tmp)-v1-Omega1.Cross(f1Tmp));
   doublereal dEpsilonPrime  = (v.Dot(vPrime))/(dElle*dL0);
   
   /* Ampiezza della forza */
   ConstitutiveLaw1DOwner::Update(dEpsilon, dEpsilonPrime);
   doublereal dF = GetF();
   
   /* Vettore forza */
   Vec3 F(v*(dF/dElle));
   
   WorkVec.Add(1, F);
   WorkVec.Add(4, f1Tmp.Cross(F));
   WorkVec.Add(7, -F);
   WorkVec.Add(10, F.Cross(f2Tmp));
   
   return WorkVec;
}

   
void RodWithOffsetJoint::Output(OutputHandler& OH) const
{
   /* Mettere magari l'output della forza, 
    * della deformazione e della velocita' di deformazione ? */
   
   if (fToBeOutput()) {
      ASSERT(dElle > DBL_EPSILON);
      Vec3 vTmp(v/dElle);
      doublereal d = GetF();
      Joint::Output(OH.Joints(), "RodWithOffs", GetLabel(),
		    Vec3(d, 0., 0.), Zero3, vTmp*d, Zero3) 
	<< " " << dElle << " " << vTmp << endl;      
   }
} 

   
/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
RodWithOffsetJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat, 
				  const VectorHandler& /* XCurr */ )
{
   DEBUGCOUT("Entering RodWithOffsetJoint::InitialAssJac()" << endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeInit(iNumRows, iNumCols, 0.);

   /* Recupera gli indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstVelIndex = iNode1FirstPosIndex+6;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstVelIndex = iNode2FirstPosIndex+6;
   
   /* Setta gli indici della matrice */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.fPutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.fPutColIndex(iCnt, iNode1FirstPosIndex+iCnt);	
      WM.fPutColIndex(6+iCnt, iNode1FirstVelIndex+iCnt);
      
      WM.fPutRowIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
      WM.fPutColIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
      WM.fPutColIndex(18+iCnt, iNode2FirstVelIndex+iCnt);
   }           
   
   Mat3x3 R1(pNode1->GetRRef());
   Mat3x3 R2(pNode2->GetRRef());
   Vec3 f1Tmp(R1*f1);
   Vec3 f2Tmp(R2*f2);
   Vec3 x1(pNode1->GetXCurr());
   Vec3 x2(pNode2->GetXCurr());
   
   Vec3 v1(pNode1->GetVCurr());
   Vec3 v2(pNode2->GetVCurr());
   Vec3 Omega1(pNode1->GetWRef());
   Vec3 Omega2(pNode2->GetWRef());
   
   /* v = p2-p1 */
   v = x2+f2Tmp-x1-f1Tmp;
   doublereal dCross = v.Dot();
   
   /* Verifica che la distanza non sia nulla */
   if (dCross <= DBL_EPSILON) {
      cerr << endl << "Null distance between nodes " << pNode1->GetLabel() 
	<< " and " << pNode2->GetLabel() 
	<< " in Rod Joint " << uLabel << ';' << endl
	<< "aborting ..." << endl;
      THROW(Joint::ErrGeneric());
   }   

   /* Lunghezza corrente */
   dElle = sqrt(dCross);

   /* Velocita' di deformazione */
   Vec3 vPrime(v2+Omega2.Cross(f2Tmp)-v1-Omega1.Cross(f1Tmp));

   /* Forza e slopes */
   doublereal dF = GetF();
   doublereal dFDE = GetFDE();
   doublereal dFDEPrime = GetFDEPrime();

   /* Vettore forza */
   Vec3 F = v*(dF/dElle);
   
   Mat3x3 K( Mat3x3(v,v*((-dF)/(dElle*dCross)))
	    +v.Tens(v*((dFDE)/(dL0*dCross)))
	    +v.Tens(v.Cross(vPrime.Cross(v*((dFDEPrime)/(dL0*dCross*dCross))))));

   Mat3x3 KPrime(v.Tens(v*((dFDEPrime)/(dL0*dCross))));
   
   /* Termini di forza diagonali */
   Mat3x3 Tmp = K;
   WM.Add(1, 1, Tmp);
   WM.Add(7, 13, Tmp);
   
   /* Termini di coppia, nodo 1 */
   Mat3x3 Tmp2 = Mat3x3(f1Tmp)*Tmp;
   WM.Add(4, 1, Tmp2);
   WM.Add(4, 13, -Tmp2);
   
   /* Termini di coppia, nodo 2 */
   Tmp2 = Mat3x3(f2Tmp)*Tmp;
   WM.Add(10, 13, Tmp2);
   WM.Add(10, 1, -Tmp2);
   
   /* termini di forza extradiagonali */
   Tmp = -Tmp;
   WM.Add(1, 13, Tmp);
   WM.Add(7, 1, Tmp);

   
   /* Termini di forza, velocita' */
   Tmp = KPrime;
   WM.Add(1, 7, Tmp);
   WM.Add(7, 19, Tmp);
   
   /* Termini di coppia, nodo 1 */
   Tmp2 = Mat3x3(f1Tmp)*Tmp;
   WM.Add(4, 7, Tmp2);
   WM.Add(4, 19, -Tmp2);
   
   /* Termini di coppia, nodo 2 */
   Tmp2 = Mat3x3(f2Tmp)*Tmp;
   WM.Add(10, 19, Tmp2);
   WM.Add(10, 7, -Tmp2);
   
   /* termini di forza extradiagonali */
   Tmp = -Tmp;
   WM.Add(1, 19, Tmp);
   WM.Add(7, 7, Tmp);
      
   
   /* Termini di rotazione, Delta g1 */
   Tmp = K*Mat3x3(-f1Tmp)+KPrime*Mat3x3(Omega1, -f1Tmp);	    
   WM.Add(1, 4, Tmp);
   WM.Add(7, 4, -Tmp);
   
   /* Termini di coppia, Delta g1 */
   Tmp2 = Mat3x3(f1Tmp)*Tmp;
   WM.Add(10, 4, -Tmp2);
   Tmp2 += Mat3x3(f1Tmp, F);
   WM.Add(4, 4, Tmp2);     
           
   /* Termini di rotazione, Delta g2 */
   Tmp = K*Mat3x3(-f2Tmp)+KPrime*Mat3x3(Omega2, -f2Tmp);
   WM.Add(7, 16, Tmp);
   WM.Add(1, 16, -Tmp);

   /* Termini di coppia, Delta g2 */
   Tmp2 = Mat3x3(f2Tmp)*Tmp;
   WM.Add(4, 16, -Tmp2);
   Tmp2 += Mat3x3(f2Tmp, F);
   WM.Add(10, 16, Tmp2);     
   
   
   /* Termini di velocita' angolare, Delta Omega1 */
   Tmp = KPrime*Mat3x3(-f1Tmp);
   WM.Add(1, 10, Tmp);
   WM.Add(7, 10, -Tmp);
   
   /* Termini di coppia, Delta Omega1 */
   Tmp2 = Mat3x3(f1Tmp)*Tmp;
   WM.Add(4, 10, Tmp2);
   WM.Add(10, 10, -Tmp2);     
   
   /* Termini di velocita' angolare, Delta Omega2 */
   Tmp = KPrime*Mat3x3(-f2Tmp);
   WM.Add(7, 22, Tmp);
   WM.Add(1, 22, -Tmp);
   
   /* Termini di coppia, Delta Omega2 */
   Tmp2 = Mat3x3(f2Tmp)*Tmp;
   WM.Add(10, 22, Tmp2);
   WM.Add(4, 22, -Tmp2);     
   
   return WorkMat;   
}

   
/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
RodWithOffsetJoint::InitialAssRes(SubVectorHandler& WorkVec,
				  const VectorHandler& /* XCurr */ )
{
   DEBUGCOUT("RodWithOffsetJoint::InitialAssRes()" << endl);
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
   
   /* Indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   
   /* Setta gli indici */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {	
      WorkVec.fPutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WorkVec.fPutRowIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
   }

   /* Dati */
   Mat3x3 R1(pNode1->GetRCurr());
   Mat3x3 R2(pNode2->GetRCurr());
   Vec3 f1Tmp(R1*f1);
   Vec3 f2Tmp(R2*f2);
   Vec3 x1(pNode1->GetXCurr());
   Vec3 x2(pNode2->GetXCurr());
   
   Vec3 v1(pNode1->GetVCurr());
   Vec3 v2(pNode2->GetVCurr());
   Vec3 Omega1(pNode1->GetWCurr());
   Vec3 Omega2(pNode2->GetWCurr());
   
   
   /* v = x2-x1 */
   v = x2+f2Tmp-x1-f1Tmp;
   doublereal dCross = v.Dot();
   
   /* Verifica che la distanza non sia nulla */
   if (dCross <= DBL_EPSILON) {
      cerr << endl << "Null distance between nodes " << pNode1->GetLabel() 
	<< " and " << pNode2->GetLabel() 
	<< " in Rod Joint " << uLabel << ';' << endl
	<< "aborting ..." << endl;
      THROW(Joint::ErrGeneric());
   }   

   /* Lunghezza corrente */
   dElle = sqrt(dCross);

   /* Deformazione */
   doublereal dEpsilon = dElle/dL0-1.;
   
   /* Velocita' di deformazione */
   Vec3 vPrime(v2+Omega2.Cross(f2Tmp)-v1-Omega1.Cross(f1Tmp));
   doublereal dEpsilonPrime  = (v.Dot(vPrime))/(dElle*dL0);
   
   /* Ampiezza della forza */
   ConstitutiveLaw1DOwner::Update(dEpsilon, dEpsilonPrime);
   doublereal dF = GetF();
   
   /* Vettore forza */
   Vec3 F(v*(dF/dElle));
   
   WorkVec.Add(1, F);
   WorkVec.Add(4, f1Tmp.Cross(F));
   WorkVec.Add(7, -F);
   WorkVec.Add(10, F.Cross(f2Tmp));
   
   return WorkVec;
}


void 
RodWithOffsetJoint::GetAdamsDummyPart(unsigned int part,
			    Vec3& x, 
			    Mat3x3& R) const 
{
   ASSERT(part == 1);
   x = pNode1->GetXCurr()+pNode1->GetRCurr()*f1;
   R = pNode1->GetRCurr();
}


ostream& 
RodWithOffsetJoint::WriteAdamsDummyPartCmd(ostream& out,
					   unsigned int part, 
					   unsigned int firstId) const
{
   Vec3 x1 = pNode1->GetXCurr()+pNode1->GetRCurr()*f1;
   Vec3 x2 = pNode2->GetXCurr()+pNode2->GetRCurr()*f2;
     
   Vec3 v1 = x2-x1;
   doublereal l = v1.Norm();
   v1 /= l;
   
   Mat3x3 Rx(Eye3-v1.Tens(v1));
   int index = 1;
   if (fabs(v1.dGet(2)) < fabs(v1.dGet(index))) {
      index = 2;
   }
   if (fabs(v1.dGet(3)) < fabs(v1.dGet(index))) {
      index = 3;
   }
   
   Vec3 v2(Rx.GetVec(index));
   v2 /= v2.Norm();
   
   Vec3 e(EulerAngles(MatR2vec(1, v1, 2, v2)));
   
   return out 
     << psAdamsElemCode[GetElemType()] << "_" << GetLabel() << "_" << part << endl
     << firstId << " "
     << x1 << " "
     << EulerAngles(pNode1->GetRCurr()) << " "
     << x1 << " "
     << e << " "
     << l << " " << 0. << " " << 0. << " "
     << Zero3 << endl;   
}

/* RodWithOffsetJoint - end */
