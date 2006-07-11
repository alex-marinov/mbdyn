/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2006
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

/* Rods */

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "ac/math.h"

#include "dataman.h"
#include "rodj.h"


/* Rod - begin */

/* Costruttore non banale */
Rod::Rod(unsigned int uL, const DofOwner* pDO, 
		   const ConstitutiveLaw1D* pCL,
		   const StructNode* pN1, const StructNode* pN2,
		   doublereal dLength, flag fOut, flag fHasOffsets)
: Elem(uL, fOut), 
Joint(uL, pDO, fOut),
ConstitutiveLaw1DOwner(pCL),
pNode1(pN1),
pNode2(pN2),
dL0(dLength),
v(0.),
dElle(0.),
dEpsilon(0.),
dEpsilonPrime(0.)
{
   /* Verifica di consistenza dei dati iniziali */   
   ASSERT(pN1 != NULL);
   ASSERT(pN1->GetNodeType() == Node::STRUCTURAL);
   ASSERT(pN2 != NULL);
   ASSERT(pN2->GetNodeType() == Node::STRUCTURAL);
   
   if (!fHasOffsets) {
      v = pN2->GetXCurr()-pN1->GetXCurr();
   
      doublereal dDot = v.Dot();
      if (dDot <= DBL_EPSILON) {
	 silent_cerr("Rod(" << GetLabel() << "): "
		 "initial length must be non-null" << std::endl);
	 throw ErrGeneric();
      }

      dElle = sqrt(dDot);
   }
   
   ASSERT(dLength > DBL_EPSILON);
}


/* Distruttore */
Rod::~Rod(void) 
{
   NO_OP;
}

void
Rod::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
	ConstitutiveLaw1DOwner::AfterConvergence(dEpsilon);
}

/* Contributo al file di restart */
std::ostream& Rod::Restart(std::ostream& out) const
{
   Joint::Restart(out) << ", rod, "
     << pNode1->GetLabel() << ", "
     << pNode2->GetLabel() << ", " 
     << dL0 << ", ";
   return pGetConstLaw()->Restart(out) << ';' << std::endl;
}


void Rod::AssMat(FullSubMatrixHandler& WorkMat, doublereal dCoef)
{
   /* v = x2-x1 */
   /* v = pNode2->GetXCurr()-pNode1->GetXCurr(); */
   doublereal dCross = v.Dot();
   
   /* Verifica che la distanza non sia nulla */
   if (dCross <= DBL_EPSILON) {
      silent_cerr("Rod(" << GetLabel() << "): "
	      "null distance between nodes " << pNode1->GetLabel() 
      	      << " and " << pNode2->GetLabel() << std::endl);
      throw Joint::ErrGeneric();
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
   WorkMat.Sub(1, 4, K);
   WorkMat.Sub(4, 1, K);
}


void Rod::AssVec(SubVectorHandler& WorkVec)
{
   /* v = x2-x1 */
   v = pNode2->GetXCurr()-pNode1->GetXCurr();
   doublereal dCross = v.Dot();
   
   /* Verifica che la distanza non sia nulla */
   if (dCross <= DBL_EPSILON) {
      silent_cerr("Rod(" << GetLabel() << "): "
	      "null distance between nodes " << pNode1->GetLabel() 
      	      << " and " << pNode2->GetLabel() << std::endl);
      throw Joint::ErrGeneric();
   }   
   
   /* Deformazione */
   dElle = sqrt(dCross);
   dEpsilon = dElle/dL0-1.;
   
   Vec3 vPrime(pNode2->GetVCurr()-pNode1->GetVCurr());
   dEpsilonPrime = (v.Dot(vPrime))/(dElle*dL0);

   /* Ampiezza della forza */
   bool ChangeJac(false);
   try {
      ConstitutiveLaw1DOwner::Update(dEpsilon);
   }
   catch (Elem::ChangedEquationStructure) {
      ChangeJac = true;
   }

   doublereal dF = GetF();
   
   /* Vettore forza */
   Vec3 F = v*(dF/dElle);
   
   WorkVec.Add(1, F);
   WorkVec.Sub(4, F);

   if (ChangeJac) {
      throw Elem::ChangedEquationStructure();
   }
}


VariableSubMatrixHandler& Rod::AssJac(VariableSubMatrixHandler& WorkMat,
					   doublereal dCoef,
					   const VectorHandler& /* XCurr */ ,
					   const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering Rod::AssJac()" << std::endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeReset(iNumRows, iNumCols);

   /* Recupera gli indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   
   /* Setta gli indici della matrice */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);	
      WM.PutRowIndex(3+iCnt, iNode2FirstMomIndex+iCnt);
      WM.PutColIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
   }  
   
   /* Genera la matrice */
   AssMat(WM, dCoef);
   
   return WorkMat;
}


void Rod::AssMats(VariableSubMatrixHandler& WorkMatA,
		      VariableSubMatrixHandler& WorkMatB,
		      const VectorHandler& /* XCurr */ ,
		      const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering Rod::AssMats()" << std::endl);
   
   WorkMatB.SetNullMatrix();   
   FullSubMatrixHandler& WMA = WorkMatA.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WMA.ResizeReset(iNumRows, iNumCols);

   /* Recupera gli indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   
   /* Setta gli indici della matrice */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WMA.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WMA.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);	
      WMA.PutRowIndex(3+iCnt, iNode2FirstMomIndex+iCnt);
      WMA.PutColIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
   }  
   
   /* Genera la matrice */
   AssMat(WMA, 1.);   
}


SubVectorHandler& Rod::AssRes(SubVectorHandler& WorkVec,
				   doublereal /* dCoef */ ,
				   const VectorHandler& /* XCurr */ ,
				   const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering Rod::AssRes()" << std::endl);
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);

   /* Indici */
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   
   /* Setta gli indici */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.PutRowIndex(3+iCnt, iNode2FirstMomIndex+iCnt);
   }
   
   /* Costruisce il vettore */
   AssVec(WorkVec);
   
   return WorkVec;   
}
   

void Rod::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) {      
      ASSERT(dElle > DBL_EPSILON);	
      Vec3 vTmp(v/dElle);      
      doublereal d = GetF();

      std::ostream& out = OH.Joints();
      
      Joint::Output(out, "Rod", GetLabel(),
		    Vec3(d, 0., 0.), Zero3, vTmp*d, Zero3)
	<< " " << dElle << " " << vTmp << " " << dEpsilonPrime*dL0, 
        ConstitutiveLaw1DOwner::OutputAppend(out) << std::endl;
   }
}

#if 0
/* Output di un modello NASTRAN equivalente nella configurazione corrente */
void
Rod::Output_pch(std::ostream& out) const
{
#if defined(__HACK_NASTRAN_MODES__)
	if (fToBeOutput()) {
		unsigned int label = GetLabel();
		if (label > 9999999) {
			silent_cerr("Rod(" << label <<"): label is too large" 
					<< std::endl);
			throw ErrGeneric();
		}

		const char *name = GetName();
		out << "$ Rod " << GetLabel();
		if (name) {
			out << " (" << name << ")";
		}

#define __NASTRAN_FORMAT__ __HACK_NASTRAN_MODES__

#if __NASTRAN_FORMAT__ == __NASTRAN_FORMAT_FIXED__
		out << std::endl
			/* PBEAM */
			<< "PBEAM   "
			<< std::setw(8) << 20000000+label	/* label */
			<< std::setw(8) << 1			/* material */
			<< std::setw(8) << 1.		/* area */
			<< std::setw(8) << 1.		/* J1 */
			<< std::setw(8) << 1.		/* J2 */
			<< std::setw(8) << ""		/* J12 */
			<< std::setw(8) << 1.		/* Jp */
			<< std::endl

			/* CBEAM */
			<< "CBEAM   "
			<< std::setw(8) << 20000000+label	/* label */
			<< std::setw(8) << 20000000+label	/* prop */
			<< std::setw(8) << pNode1->GetLabel()	/* node 1 */
			<< std::setw(8) << pNode2->GetLabel()	/* node 2 */
			<< enld;
#elif __NASTRAN_FORMAT__ == __NASTRAN_FORMAT_FIXED16__
		out << std::endl
			/* PBEAM */
			<< "PBEAM*  "
			<< std::setw(16) << 20000000+label	/* label */
			<< std::setw(16) << 1		/* material */
			<< std::setw(16) << 1.		/* area */
			<< std::setw(16) << 1.		/* J1 */
			<< "*" << std::setw(7) << 1
			<< std::endl
			<< "*" << std::setw(7) << 1
			<< std::setw(16) << 1.		/* J2 */
			<< std::setw(16) << " "		/* J12 */
			<< std::setw(16) << 1.		/* Jp */
			<< std::endl

			/* CBEAM */
			<< "CBEAM*  "
			<< std::setw(16) << 20000000+label 	/* label */
			<< std::setw(16) << 20000000+label	/* prop */
			<< std::setw(16) << pNode1->GetLabel()	/* node 1 */
			<< std::setw(16) << pNode2->GetLabel()	/* node 2 */
			<< std::endl;
#elif __NASTRAN_FORMAT__ == __NASTRAN_FORMAT_FREE__
		out << std::endl
			/* PBEAM */
			<< "PBEAM," 
			<< 20000000+label << ","
			<< 1 << ","
			<< 1. << ","
			<< 1. << ","
			<< 1. << ","
			<< ","
			<< 1. << std::endl

			/* CBEAM */
			<< "CBEAM,"
			<< 20000000+label << ","
			<< 20000000+label << ","
			<< pNode1->GetLabel() << ","
			<< pNode2->GetLabel() << std::endl;
#else
#error "unknown NASTRAN format"
#endif
	}
#endif /* __HACK_NASTRAN_MODES__ */
}
#endif 

VariableSubMatrixHandler& 
Rod::InitialAssJac(VariableSubMatrixHandler& WorkMat,
			const VectorHandler& /* XCurr */ )
{
   DEBUGCOUT("Entering Rod::InitialAssJac()" << std::endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeReset(iNumRows, iNumCols);

   /* Recupera gli indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   
   /* Setta gli indici della matrice */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);	
      WM.PutRowIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
      WM.PutColIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
   }  
   
   /* Genera la matrice */
   AssMat(WM);
   
   return WorkMat;
}


SubVectorHandler& Rod::InitialAssRes(SubVectorHandler& WorkVec,
					  const VectorHandler& /* XCurr */ )
{
   DEBUGCOUT("Entering Rod::InitialAssRes()" << std::endl);
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);

   /* Indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   
   /* Setta gli indici */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WorkVec.PutRowIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
   }
   
   /* Costruisce il vettore */
   AssVec(WorkVec);
   
   return WorkVec;   
}


void 
Rod::GetDummyPartPos(unsigned int part,
			    Vec3& x, 
			    Mat3x3& R) const 
{
   ASSERT(part == 1);
   x = pNode1->GetXCurr();
   R = pNode1->GetRCurr();
}

void 
Rod::GetDummyPartVel(unsigned int part,
			    Vec3& v, 
			    Vec3& w) const 
{
   ASSERT(part == 1);
   v = pNode1->GetVCurr();
   w = pNode1->GetWCurr();
}


#ifdef USE_ADAMS
std::ostream& 
Rod::WriteAdamsDummyPartCmd(std::ostream& out,
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
   
   Vec3 e(MatR2EulerAngles(MatR2vec(1, v1, 2, v2))*dRaDegr);
   
   return out 
     << psAdamsElemCode[GetElemType()] << "_" << GetLabel() << "_" << part << std::endl
     << firstId << " "
     << x1 << " "
     << MatR2EulerAngles(pNode1->GetRCurr())*dRaDegr << " "
     << x1 << " "
     << e << " "
     << l << " " << 0. << " " << 0. << " "
     << Zero3 << std::endl;
}
#endif /* USE_ADAMS */


unsigned int
Rod::iGetNumPrivData(void) const 
{
	return 3 + ConstitutiveLaw1DOwner::iGetNumPrivData();
}

unsigned int
Rod::iGetPrivDataIdx(const char *s) const 
{
	ASSERT(s != NULL);

	if (strcmp(s, "F") == 0) {
		return 1;
	}

	if (strcmp(s, "L") == 0) {
		return 2;
	}

	if (strcmp(s, "LPrime") == 0) {
		return 3;
	}

	size_t l = sizeof("constitutiveLaw.") - 1;
	if (strncmp(s, "constitutiveLaw.", l) == 0) {
		return 3 + ConstitutiveLaw1DOwner::iGetPrivDataIdx(s + l);
	}

	/* error; handle later */
	return 0;
}

doublereal
Rod::dGetPrivData(unsigned int i) const
{
	ASSERT(i > 0);

	switch (i) {
	case 1:
		return GetF();

	case 2:
		return dElle;

	case 3:
		return dL0*dEpsilonPrime;
	}

	i -= 3;
	ASSERT(i <= ConstitutiveLaw1DOwner::iGetNumPrivData());

	return ConstitutiveLaw1DOwner::dGetPrivData(i);
}

/* Rod - end */


/* ViscoElasticRod - begin */

/* Costruttore non banale */
ViscoElasticRod::ViscoElasticRod(unsigned int uL, 
					   const DofOwner* pDO,
					   const ConstitutiveLaw1D* pCL,
					   const StructNode* pN1, 
					   const StructNode* pN2,
					   doublereal dLength, flag fOut)
: Elem(uL, fOut), 
Rod(uL, pDO, pCL, pN1, pN2, dLength, fOut)
{
   NO_OP;
}


/* Distruttore */
ViscoElasticRod::~ViscoElasticRod(void) 
{ 
   NO_OP; 
}

void
ViscoElasticRod::AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP)
{
	ConstitutiveLaw1DOwner::AfterConvergence(dEpsilon, dEpsilonPrime);
}

VariableSubMatrixHandler& 
ViscoElasticRod::AssJac(VariableSubMatrixHandler& WorkMat,
			     doublereal dCoef,
			     const VectorHandler& /* XCurr */ ,
			     const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering ViscoElasticRod::AssJac()" << std::endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeReset(iNumRows, iNumCols);

   /* Recupera gli indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   
   /* Setta gli indici della matrice */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);	
      WM.PutRowIndex(3+iCnt, iNode2FirstMomIndex+iCnt);
      WM.PutColIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
   }           
   
   /* v = x2-x1 */
   /* v(pNode2->GetXCurr()-pNode1->GetXCurr()); */
   doublereal dCross = v.Dot();
   
   /* Verifica che la distanza non sia nulla */
   if (dCross <= DBL_EPSILON) {
      silent_cerr("ViscoElasticRod(" << GetLabel() << "): "
	      "null distance between nodes " << pNode1->GetLabel() 
      	      << " and " << pNode2->GetLabel() << std::endl);
      throw Joint::ErrGeneric();
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
   WM.Sub(1, 4, K);
   WM.Sub(4, 1, K);
   
   return WorkMat;
}


SubVectorHandler& 
ViscoElasticRod::AssRes(SubVectorHandler& WorkVec,
			     doublereal /* dCoef */ ,
			     const VectorHandler& /* XCurr */ ,
			     const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering ViscoElasticRod::AssRes()" << std::endl);
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);

   /* Indici */
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   
   /* Setta gli indici */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.PutRowIndex(3+iCnt, iNode2FirstMomIndex+iCnt);
   }
         
   /* v = x2-x1 */
   v = pNode2->GetXCurr()-pNode1->GetXCurr();
   doublereal dCross = v.Dot();
   
   /* Verifica che la distanza non sia nulla */
   if (dCross <= DBL_EPSILON) {
      silent_cerr("ViscoElasticRod(" << GetLabel() << "): "
	      "null distance between nodes " << pNode1->GetLabel() 
      	      << " and " << pNode2->GetLabel() << std::endl);
      throw Joint::ErrGeneric();
   }   

   /* Lunghezza corrente */
   dElle = sqrt(dCross);

   /* Deformazione */
   dEpsilon = dElle/dL0-1.;
   
   /* Velocita' di deformazione */
   Vec3 vPrime(pNode2->GetVCurr()-pNode1->GetVCurr());
   dEpsilonPrime = (v.Dot(vPrime))/(dElle*dL0);
   
   /* Ampiezza della forza */
   bool ChangeJac(false);
   try {
      ConstitutiveLaw1DOwner::Update(dEpsilon, dEpsilonPrime);
   }
   catch (Elem::ChangedEquationStructure) {
      ChangeJac = true;
   }
   doublereal dF = GetF();
   
   /* Vettore forza */
   Vec3 F(v*(dF/dElle));
   
   WorkVec.Add(1, F);
   WorkVec.Sub(4, F);

   if (ChangeJac) {
      throw Elem::ChangedEquationStructure();
   }
   
   return WorkVec;
}


VariableSubMatrixHandler& 
ViscoElasticRod::InitialAssJac(VariableSubMatrixHandler& WorkMat,
				    const VectorHandler& /* XCurr */ )
{
   DEBUGCOUT("Entering ViscoElasticRod::InitialAssJac()" << std::endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeReset(iNumRows, iNumCols);

   /* Recupera gli indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstVelIndex = iNode1FirstPosIndex+6;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstVelIndex = iNode2FirstPosIndex+6;
   
   /* Setta gli indici della matrice */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);	
      WM.PutColIndex(3+iCnt, iNode1FirstVelIndex+iCnt);
      
      WM.PutRowIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
      WM.PutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
      WM.PutColIndex(9+iCnt, iNode2FirstVelIndex+iCnt);
   }           
   
   /* v = x2-x1 */
   /* v(pNode2->GetXCurr()-pNode1->GetXCurr()); */
   doublereal dCross = v.Dot();
   
   /* Verifica che la distanza non sia nulla */
   if (dCross <= DBL_EPSILON) {
      silent_cerr("ViscoElasticRod(" << GetLabel() << "): "
	      "null distance between nodes " << pNode1->GetLabel() 
      	      << " and " << pNode2->GetLabel() << std::endl);
      throw Joint::ErrGeneric();
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
   WM.Sub(1, 7, K);
   WM.Sub(4, 1, K);

   WM.Sub(1, 10, KPrime);
   WM.Sub(4, 4, KPrime);
   
   return WorkMat;
}


SubVectorHandler& 
ViscoElasticRod::InitialAssRes(SubVectorHandler& WorkVec,
				    const VectorHandler& /* XCurr */ )
{
   DEBUGCOUT("Entering ViscoElasticRod::InitialAssRes()" << std::endl);
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);

   /* Indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   
   /* Setta gli indici */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WorkVec.PutRowIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
   }
         
   /* v = x2-x1 */
   v = pNode2->GetXCurr()-pNode1->GetXCurr();
   doublereal dCross = v.Dot();
   
   /* Verifica che la distanza non sia nulla */
   if (dCross <= DBL_EPSILON) {
      silent_cerr("ViscoElasticRod(" << GetLabel() << "): "
	      "null distance between nodes " << pNode1->GetLabel() 
      	      << " and " << pNode2->GetLabel() << std::endl);
      throw Joint::ErrGeneric();
   }   

   /* Lunghezza corrente */
   dElle = sqrt(dCross);

   /* Deformazione */
   dEpsilon = dElle/dL0-1.;
   
   /* Velocita' di deformazione */
   Vec3 vPrime(pNode2->GetVCurr()-pNode1->GetVCurr());
   dEpsilonPrime = (v.Dot(vPrime))/(dElle*dL0);
   
   /* Ampiezza della forza */
   ConstitutiveLaw1DOwner::Update(dEpsilon, dEpsilonPrime);
   doublereal dF = GetF();
   
   /* Vettore forza */
   Vec3 F(v*(dF/dElle));
   
   WorkVec.Add(1, F);
   WorkVec.Sub(4, F);
   
   return WorkVec;
}

/* ViscoElasticRod - end */


/* RodWithOffset - begin */

/* Costruttore non banale */
RodWithOffset::RodWithOffset(unsigned int uL, 
				       const DofOwner* pDO,
				       const ConstitutiveLaw1D* pCL,
				       const StructNode* pN1, 
				       const StructNode* pN2,
				       const Vec3& f1Tmp, 
				       const Vec3& f2Tmp,
				       doublereal dLength, 
				       flag fOut)
: Elem(uL, fOut), 
Rod(uL, pDO, pCL, pN1, pN2, dLength, fOut, 1),
f1(f1Tmp),
f2(f2Tmp)
{
   /* Verifica di consistenza dei dati iniziali */   
   ASSERT(pN1 != NULL);
   ASSERT(pN1->GetNodeType() == Node::STRUCTURAL);
   ASSERT(pN2 != NULL);
   ASSERT(pN2->GetNodeType() == Node::STRUCTURAL);
   
   v = pN2->GetXCurr()+(pN2->GetRCurr()*f2Tmp)
     -pN1->GetXCurr()-(pN1->GetRCurr()*f1Tmp);
   
   doublereal dDot = v.Dot();
   if (dDot <= DBL_EPSILON) {
      silent_cerr("RodWithOffset(" << GetLabel() << "): "
	      "inital length must be non-null" << std::endl);
      throw ErrGeneric();
   }

   dElle = sqrt(dDot);
   
   ASSERT(dLength > DBL_EPSILON);
}


/* Distruttore */
RodWithOffset::~RodWithOffset(void)
{
   NO_OP;
}
      
   
/* Contributo al file di restart */
std::ostream& RodWithOffset::Restart(std::ostream& out) const
{
   Joint::Restart(out) << ", rod, "
     << pNode1->GetLabel() << ", "
     << pNode2->GetLabel() << ", " 
     << dL0 << ", offset, reference, node, ",
     f1.Write(out, ", ") << ", reference, node, ",
     f2.Write(out, ", ") << ", ";
   return pGetConstLaw()->Restart(out) << ';' << std::endl;
}

void
RodWithOffset::AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP)
{
	ConstitutiveLaw1DOwner::AfterConvergence(dEpsilon, dEpsilonPrime);
}
         
VariableSubMatrixHandler& 
RodWithOffset::AssJac(VariableSubMatrixHandler& WorkMat, 
			   doublereal dCoef,
			   const VectorHandler& /* XCurr */ ,
			   const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering RodWithOffset::AssJac()" << std::endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeReset(iNumRows, iNumCols);

   /* Recupera gli indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   
   /* Setta gli indici della matrice */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);	
      WM.PutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
      WM.PutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
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
      silent_cerr("RodWithOffset(" << GetLabel() << "): "
	      "null distance between nodes " << pNode1->GetLabel() 
      	      << " and " << pNode2->GetLabel() << std::endl);
      throw Joint::ErrGeneric();
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
   WM.Sub(4, 7, Tmp2);
   
   /* Termini di coppia, nodo 2 */
   Tmp2 = Mat3x3(f2Tmp)*Tmp;
   WM.Add(10, 7, Tmp2);
   WM.Sub(10, 1, Tmp2);
   
   /* termini di forza extradiagonali */
   WM.Sub(1, 7, Tmp);
   WM.Sub(7, 1, Tmp);
   
   
   /* Termini di rotazione, Delta g1 */
   Tmp = K*Mat3x3(-f1Tmp)
     +KPrime*(Mat3x3(Omega1.Cross(f1Tmp*dCoef))-Mat3x3(f1Tmp));
   WM.Add(1, 4, Tmp);
   WM.Sub(7, 4, Tmp);
   
   /* Termini di coppia, Delta g1 */
   Tmp2 = Mat3x3(f1Tmp)*Tmp;
   WM.Sub(10, 4, Tmp2);
   Tmp2 += Mat3x3(f1Tmp*dCoef, F);
   WM.Add(4, 4, Tmp2);     
           
   /* Termini di rotazione, Delta g2 */
   Tmp = K*Mat3x3(-f2Tmp)
     +KPrime*(Mat3x3(Omega2.Cross(f2Tmp*dCoef))-Mat3x3(f2Tmp));   
   WM.Add(7, 10, Tmp);
   WM.Sub(1, 10, Tmp);

   /* Termini di coppia, Delta g2 */
   Tmp2 = Mat3x3(f2Tmp)*Tmp;
   WM.Sub(4, 10, Tmp2);
   Tmp2 += Mat3x3(f2Tmp*dCoef, F);
   WM.Add(10, 10, Tmp2);     
   
   return WorkMat;   
}

	   
SubVectorHandler& 
RodWithOffset::AssRes(SubVectorHandler& WorkVec,
			   doublereal /* dCoef */ ,
			   const VectorHandler& /* XCurr */ ,
			   const VectorHandler& /* XPrimeCurr */ )
{   
   DEBUGCOUT("RodWithOffset::AssRes()" << std::endl);
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);

   /* Indici */
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   
   /* Setta gli indici */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {	
      WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.PutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
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
      silent_cerr("RodWithOffset(" << GetLabel() << "): "
	      "null distance between nodes " << pNode1->GetLabel() 
      	      << " and " << pNode2->GetLabel() << std::endl);
      throw Joint::ErrGeneric();
   }   

   /* Lunghezza corrente */
   dElle = sqrt(dCross);

   /* Deformazione */
   dEpsilon = dElle/dL0-1.;
   
   /* Velocita' di deformazione */
   Vec3 vPrime(v2+Omega2.Cross(f2Tmp)-v1-Omega1.Cross(f1Tmp));
   dEpsilonPrime  = (v.Dot(vPrime))/(dElle*dL0);
   
   /* Ampiezza della forza */
   bool ChangeJac(false);
   try {
      ConstitutiveLaw1DOwner::Update(dEpsilon, dEpsilonPrime);
   }
   catch (Elem::ChangedEquationStructure) {
      ChangeJac = true;
   }
   
   doublereal dF = GetF();
   
   /* Vettore forza */
   Vec3 F(v*(dF/dElle));
   
   WorkVec.Add(1, F);
   WorkVec.Add(4, f1Tmp.Cross(F));
   WorkVec.Sub(7, F);
   WorkVec.Sub(10, f2Tmp.Cross(F));

   if (ChangeJac) {
      throw Elem::ChangedEquationStructure();
   }
   
   return WorkVec;
}

#if 0
/* Output di un modello NASTRAN equivalente nella configurazione corrente */
void
RodWithOffset::Output_pch(std::ostream& out) const
{
#if defined(__HACK_NASTRAN_MODES__)
	if (fToBeOutput()) {
		unsigned int label = GetLabel();
		if (label > 9999999) {
			silent_cerr("Rod(" << label <<"): label is too large"
					<< std::endl);
			throw ErrGeneric();
		}

		const char *name = GetName();
		out << "$ Rod " << GetLabel();
		if (name) {
			out << " (" << name << ")";
		}

#define __NASTRAN_FORMAT__ __HACK_NASTRAN_MODES__
		Vec3 F1(pNode1->GetRCurr()*f1);
		Vec3 F2(pNode2->GetRCurr()*f2);

#if __NASTRAN_FORMAT__ == __NASTRAN_FORMAT_FIXED__
		out << std::endl
			/* PBEAM */
			<< "PBEAM   "
			<< std::setw(8) << 20000000+label	/* label */
			<< std::setw(8) << 1			/* material */
			<< std::setw(8) << 1.		/* area */
			<< std::setw(8) << 1.		/* J1 */
			<< std::setw(8) << 1.		/* J2 */
			<< std::setw(8) << " "		/* J12 */
			<< std::setw(8) << 1.		/* Jp */
			<< std::endl

			/* CBEAM */
			<< "CBEAM   "
			<< std::setw(8) << 20000000+label	/* label */
			<< std::setw(8) << 20000000+label	/* prop */
			<< std::setw(8) << pNode1->GetLabel()	/* node 1 */
			<< std::setw(8) << pNode2->GetLabel()	/* node 2 */
			<< std::setw(32) << " " 
			<< "+" << std::setw(7) << 1
			<< std::endl
			<< "+" << std::setw(7) << 1
			<< std::setw(16) << " "
			<< std::setw(8) << F1.dGet(1)
			<< std::setw(8) << F1.dGet(2)
			<< std::setw(8) << F1.dGet(3)
			<< std::setw(8) << F2.dGet(1)
			<< std::setw(8) << F2.dGet(2)
			<< std::setw(8) << F2.dGet(3)
			<< enld;
#elif __NASTRAN_FORMAT__ == __NASTRAN_FORMAT_FIXED16__
		out << std::endl
			/* PBEAM */
			<< "PBEAM*  "
			<< std::setw(16) << 20000000+label	/* label */
			<< std::setw(16) << 1		/* material */
			<< std::setw(16) << 1.		/* area */
			<< std::setw(16) << 1.		/* J1 */
			<< "*" << std::setw(7) << 1
			<< std::endl
			<< "*" << std::setw(7) << 1
			<< std::setw(16) << 1.		/* J2 */
			<< std::setw(16) << " "		/* J12 */
			<< std::setw(16) << 1.		/* Jp */
			<< std::endl

			/* CBEAM */
			<< "CBEAM*  "
			<< std::setw(16) << 20000000+label 	/* label */
			<< std::setw(16) << 20000000+label	/* prop */
			<< std::setw(16) << pNode1->GetLabel()	/* node 1 */
			<< std::setw(16) << pNode2->GetLabel()	/* node 2 */
			<< "*" << std::setw(7) << 1
			<< std::endl
			<< "*" << std::setw(7) << 1
			<< std::setw(64) << " "
			<< "*" << std::setw(7) << 2
			<< std::endl
			<< "*" << std::setw(7) << 2
			<< std::setw(32) << " "
			<< std::setw(16) << F1.dGet(1)
			<< std::setw(16) << F1.dGet(2)
			<< "*" << std::setw(7) << 3
			<< std::endl
			<< "*" << std::setw(7) << 3
			<< std::setw(16) << F1.dGet(3)
			<< std::setw(16) << F2.dGet(1)
			<< std::setw(16) << F2.dGet(2)
			<< std::setw(16) << F2.dGet(3)
			<< std::endl;
#elif __NASTRAN_FORMAT__ == __NASTRAN_FORMAT_FREE__
		out << std::endl
			/* PBEAM */
			<< "PBEAM," 
			<< 20000000+label << ","
			<< 1 << ","
			<< 1. << ","
			<< 1. << ","
			<< 1. << ","
			<< ","
			<< 1. << std::endl

			/* CBEAM */
			<< "CBEAM,"
			<< 20000000+label << ","
			<< 20000000+label << ","
			<< pNode1->GetLabel() << ","
			<< pNode2->GetLabel() << ",,,,"
#if 0
			<< "," 
#endif
			<< std::endl
#if 1
			<< "," 
#endif
			<< " ,,", F1.Write(out, ",") << ",", F2.Write(out, ",")
			<< std::endl;
#else
#error "unknown NASTRAN format"
#endif
	}
#endif /* __HACK_NASTRAN_MODES__ */
}
#endif

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
RodWithOffset::InitialAssJac(VariableSubMatrixHandler& WorkMat, 
				  const VectorHandler& /* XCurr */ )
{
   DEBUGCOUT("Entering RodWithOffset::InitialAssJac()" << std::endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeReset(iNumRows, iNumCols);

   /* Recupera gli indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstVelIndex = iNode1FirstPosIndex+6;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstVelIndex = iNode2FirstPosIndex+6;
   
   /* Setta gli indici della matrice */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);	
      WM.PutColIndex(6+iCnt, iNode1FirstVelIndex+iCnt);
      
      WM.PutRowIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
      WM.PutColIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
      WM.PutColIndex(18+iCnt, iNode2FirstVelIndex+iCnt);
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
      silent_cerr("RodWithOffset(" << GetLabel() << "): "
	      "null distance between nodes " << pNode1->GetLabel() 
      	      << " and " << pNode2->GetLabel() << std::endl);
      throw Joint::ErrGeneric();
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
   WM.Sub(4, 13, Tmp2);
   
   /* Termini di coppia, nodo 2 */
   Tmp2 = Mat3x3(f2Tmp)*Tmp;
   WM.Add(10, 13, Tmp2);
   WM.Sub(10, 1, Tmp2);
   
   /* termini di forza extradiagonali */
   WM.Sub(1, 13, Tmp);
   WM.Sub(7, 1, Tmp);

   
   /* Termini di forza, velocita' */
   Tmp = KPrime;
   WM.Add(1, 7, Tmp);
   WM.Add(7, 19, Tmp);
   
   /* Termini di coppia, nodo 1 */
   Tmp2 = Mat3x3(f1Tmp)*Tmp;
   WM.Add(4, 7, Tmp2);
   WM.Sub(4, 19, Tmp2);
   
   /* Termini di coppia, nodo 2 */
   Tmp2 = Mat3x3(f2Tmp)*Tmp;
   WM.Add(10, 19, Tmp2);
   WM.Sub(10, 7, Tmp2);
   
   /* termini di forza extradiagonali */
   WM.Sub(1, 19, Tmp);
   WM.Sub(7, 7, Tmp);
      
   
   /* Termini di rotazione, Delta g1 */
   Tmp = K*Mat3x3(-f1Tmp)+KPrime*Mat3x3(Omega1, -f1Tmp);	    
   WM.Add(1, 4, Tmp);
   WM.Sub(7, 4, Tmp);
   
   /* Termini di coppia, Delta g1 */
   Tmp2 = Mat3x3(f1Tmp)*Tmp;
   WM.Sub(10, 4, Tmp2);
   Tmp2 += Mat3x3(f1Tmp, F);
   WM.Add(4, 4, Tmp2);     
           
   /* Termini di rotazione, Delta g2 */
   Tmp = K*Mat3x3(-f2Tmp)+KPrime*Mat3x3(Omega2, -f2Tmp);
   WM.Add(7, 16, Tmp);
   WM.Sub(1, 16, Tmp);

   /* Termini di coppia, Delta g2 */
   Tmp2 = Mat3x3(f2Tmp)*Tmp;
   WM.Sub(4, 16, Tmp2);
   Tmp2 += Mat3x3(f2Tmp, F);
   WM.Add(10, 16, Tmp2);     
   
   
   /* Termini di velocita' angolare, Delta Omega1 */
   Tmp = KPrime*Mat3x3(-f1Tmp);
   WM.Add(1, 10, Tmp);
   WM.Sub(7, 10, Tmp);
   
   /* Termini di coppia, Delta Omega1 */
   Tmp2 = Mat3x3(f1Tmp)*Tmp;
   WM.Add(4, 10, Tmp2);
   WM.Sub(10, 10, Tmp2);     
   
   /* Termini di velocita' angolare, Delta Omega2 */
   Tmp = KPrime*Mat3x3(-f2Tmp);
   WM.Add(7, 22, Tmp);
   WM.Sub(1, 22, Tmp);
   
   /* Termini di coppia, Delta Omega2 */
   Tmp2 = Mat3x3(f2Tmp)*Tmp;
   WM.Add(10, 22, Tmp2);
   WM.Sub(4, 22, Tmp2);     
   
   return WorkMat;   
}

   
/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
RodWithOffset::InitialAssRes(SubVectorHandler& WorkVec,
				  const VectorHandler& /* XCurr */ )
{
   DEBUGCOUT("RodWithOffset::InitialAssRes()" << std::endl);
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);

   /* Indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   
   /* Setta gli indici */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {	
      WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WorkVec.PutRowIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
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
      silent_cerr("RodWithOffset(" << GetLabel() << "): "
	      "null distance between nodes " << pNode1->GetLabel() 
      	      << " and " << pNode2->GetLabel() << std::endl);
      throw Joint::ErrGeneric();
   }   

   /* Lunghezza corrente */
   dElle = sqrt(dCross);

   /* Deformazione */
   dEpsilon = dElle/dL0-1.;
   
   /* Velocita' di deformazione */
   Vec3 vPrime(v2+Omega2.Cross(f2Tmp)-v1-Omega1.Cross(f1Tmp));
   dEpsilonPrime  = (v.Dot(vPrime))/(dElle*dL0);
   
   /* Ampiezza della forza */
   ConstitutiveLaw1DOwner::Update(dEpsilon, dEpsilonPrime);
   doublereal dF = GetF();
   
   /* Vettore forza */
   Vec3 F(v*(dF/dElle));
   
   WorkVec.Add(1, F);
   WorkVec.Add(4, f1Tmp.Cross(F));
   WorkVec.Sub(7, F);
   WorkVec.Sub(10, f2Tmp.Cross(F));
   
   return WorkVec;
}


void 
RodWithOffset::GetDummyPartPos(unsigned int part,
			    Vec3& x, 
			    Mat3x3& R) const 
{
   ASSERT(part == 1);
   x = pNode1->GetXCurr()+pNode1->GetRCurr()*f1;
   R = pNode1->GetRCurr();
}

void 
RodWithOffset::GetDummyPartVel(unsigned int part,
			    Vec3& v, 
			    Vec3& w) const 
{
   ASSERT(part == 1);
   w = pNode1->GetWCurr();
   v = pNode1->GetVCurr()+w.Cross(pNode1->GetRCurr()*f1);
}

#ifdef USE_ADAMS
std::ostream& 
RodWithOffset::WriteAdamsDummyPartCmd(std::ostream& out,
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
   
   Vec3 e(MatR2EulerAngles(MatR2vec(1, v1, 2, v2))*dRaDegr);
   
   return out 
     << psAdamsElemCode[GetElemType()] << "_" << GetLabel() << "_" << part << std::endl
     << firstId << " "
     << x1 << " "
     << MatR2EulerAngles(pNode1->GetRCurr())*dRaDegr << " "
     << x1 << " "
     << e << " "
     << l << " " << 0. << " " << 0. << " "
     << Zero3 << std::endl;   
}
#endif /* USE_ADAMS */

/* RodWithOffset - end */

