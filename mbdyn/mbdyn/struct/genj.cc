/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2011
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

/* Vincoli generali */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <limits>

#include "dataman.h"
#include "genj.h"
#include "hint.h"
#include "hint_impl.h"
#include "Rot.hh"

#ifndef MBDYN_X_DISTANCE_JOINT

/* DistanceJoint - begin */

/* Costruttore non banale */
DistanceJoint::DistanceJoint(unsigned int uL, const DofOwner* pDO,
			     const StructNode* pN1, const StructNode* pN2,
			     const DriveCaller* pDC, flag fOut)
: Elem(uL, fOut), 
Joint(uL, pDO, fOut),
DriveOwner(pDC),
pNode1(pN1), pNode2(pN2), v(Zero3), dAlpha(0.)
{
   NO_OP;
}


/* Distruttore banale - ci pensa il DriveOwner a distruggere il DriveCaller */
DistanceJoint::~DistanceJoint(void) 
{ 
   NO_OP; 
}


/* Dati privati */
unsigned int
DistanceJoint::iGetNumPrivData(void) const
{
	return 1;
}

unsigned int
DistanceJoint::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);

	if (strcmp(s, "d") == 0) {
		return 1;
	}

	return 0;
}

doublereal
DistanceJoint::dGetPrivData(unsigned int i) const
{
	ASSERT(i == 1);

	if (i == 1) {
		return dGet();
	}

	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

/* Contributo al file di restart */
std::ostream& DistanceJoint::Restart(std::ostream& out) const
{
   Joint::Restart(out) << ", distance, " 
     << pNode1->GetLabel() << ", "
     << pNode2->GetLabel() << ", ";
   return pGetDriveCaller()->Restart(out) << ';' << std::endl;
}


/* Assemblaggio jacobiano */
VariableSubMatrixHandler& 
DistanceJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		      doublereal dCoef,
		      const VectorHandler& /* XCurr */ ,
		      const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering DistanceJoint::AssJac()" << std::endl);
      
   SparseSubMatrixHandler& WM = WorkMat.SetSparse();
   
   /* Dimensiona e resetta la matrice di lavoro */
   doublereal dDistance = pGetDriveCaller()->dGet();
   
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   
      /* 
       *       forza:     dAlpha*v
       *       distanza:  x2 - x1 = d*v
       *       equazione di normalizzazione di v:
       *                  (v^T*v)^1/2 = 1
       *
       * 
       *        x1   Q1   x2   Q2   v    a 
       * 
       * x1  |  0    0    0    0    0    0 |
       * Q1  |  0    0    0    0   -a*I -v |
       * x2  |  0    0    0    0    0    0 |
       * Q2  |  0    0    0    0    a*I  v |
       * v   | -c*I  0    c*I  0   -d*I  0 |
       * a   |  0    0    0    0    v^T  0 |
       *                        ----------
       *                        (vT*v)^1/2
       * 
       * con c = dCoef, a = dAlpha
       */
      
      WM.ResizeReset(24, 0);
      
      /* Nota: la direzione della reazione, v,
       * ed il coefficiente di amplificazione dAlpha
       * sono stati aggiornati durante il calcolo del residuo */
      
      for (int iCnt = 1; iCnt <= 3; iCnt++) {	   
	 /* termini di Delta_x1 */
	 WM.PutItem(iCnt, iFirstReactionIndex+iCnt,
		     iNode1FirstPosIndex+iCnt, -dCoef);
	 
	 /* termini di Delta_x2 */
	 WM.PutItem(3+iCnt, iFirstReactionIndex+iCnt,
		     iNode2FirstPosIndex+iCnt, dCoef);
	 
	 /* termini di Delta_F sul nodo 1 */
	 WM.PutItem(6+iCnt, iNode1FirstMomIndex+iCnt,
		     iFirstReactionIndex+iCnt, -dAlpha);
	 
	 /* termini di Delta_F sul nodo 2 */
	 WM.PutItem(9+iCnt, iNode2FirstMomIndex+iCnt,
		     iFirstReactionIndex+iCnt, dAlpha);
	 
	 /* termini diagonali di Delta_v */
	 WM.PutItem(12+iCnt, iFirstReactionIndex+iCnt,
		     iFirstReactionIndex+iCnt, -dDistance);
	 
	 /* termini di Delta_alpha sul nodo 1 */
	 WM.PutItem(15+iCnt, iNode1FirstMomIndex+iCnt,
		     iFirstReactionIndex+4, -v.dGet(iCnt));
	 
	 /* termini di Delta_alpha sul nodo 2 */
	 WM.PutItem(18+iCnt, iNode2FirstMomIndex+iCnt,
		     iFirstReactionIndex+4, v.dGet(iCnt));
      }	
      
      doublereal d = v.Dot();
      ASSERT(d > std::numeric_limits<doublereal>::epsilon());
      if (d > std::numeric_limits<doublereal>::epsilon()) {	   
	 d = std::sqrt(d);
	 /* termini di Delta_v su alpha */
	 for (int iCnt = 3; iCnt > 0; iCnt--) {	      
	    WM.PutItem(21+iCnt, iFirstReactionIndex+4,
			iFirstReactionIndex+iCnt, v.dGet(iCnt)/d);
	 }	   				
      }
   
   return WorkMat;
}
   

/* Assemblaggio residuo */
SubVectorHandler& DistanceJoint::AssRes(SubVectorHandler& WorkVec,
					doublereal dCoef,
					const VectorHandler& XCurr, 
					const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering DistanceJoint::AssRes()" << std::endl);
      
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);
 
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();  
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {      
      /* Indici del nodo 1 */
      WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
   
      /* Indici del nodo 2 */
      WorkVec.PutRowIndex(3+iCnt, iNode2FirstMomIndex+iCnt);
   }
   
   /* Indici del vincolo */
   for (int iCnt = 1; iCnt <= 4; iCnt++) {      
      WorkVec.PutRowIndex(6+iCnt, iFirstReactionIndex+iCnt);
   }   

   Vec3 x1(pNode1->GetXCurr());
   Vec3 x2(pNode2->GetXCurr());
   
   /* Aggiorna i dati propri */
   v = Vec3(XCurr, iFirstReactionIndex+1);
   dAlpha = XCurr(iFirstReactionIndex+4);   

   doublereal dDistance = pGetDriveCaller()->dGet();

   /* Distanza nulla */
   if (fabs(dDistance) <= std::numeric_limits<doublereal>::epsilon()) {	
      silent_cerr("DistanceJoint(" << GetLabel() << "): "
	      "near-zero distance" << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

   Vec3 TmpVec(v*dAlpha);
   WorkVec.Add(1, TmpVec);
   WorkVec.Sub(4, TmpVec);	
      
   WorkVec.Add(7, x1-x2+v*dDistance);
      
   doublereal d = v.Dot();
   ASSERT(d > std::numeric_limits<doublereal>::epsilon());
   if (d > std::numeric_limits<doublereal>::epsilon()) {	   
      d = std::sqrt(d);
   } else {	   
      d = 0.;
   }
      
   WorkVec.PutCoef(10, 1.-d);	  
   
   return WorkVec;
}


void DistanceJoint::Output(OutputHandler& OH) const
{
   if(fToBeOutput()) {      
      doublereal d = dGet();
      Vec3 vTmp;      
      if (fabs(d) > std::numeric_limits<doublereal>::epsilon()) {
	 vTmp = Vec3(dAlpha, 0., 0.);
      } else {      
	 vTmp = v;
      }
      Joint::Output(OH.Joints(), "Distance", GetLabel(),
		    vTmp, Zero3, v*dAlpha, Zero3)
	<< " " << v << " " << d << std::endl;
   }
}



/* Nota: vanno modificati in analogia al DistanceWithOffsetJoint */

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
DistanceJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
			     const VectorHandler& XCurr)
{
   DEBUGCOUT("Entering DistanceJoint::InitialAssJac()" << std::endl);
   
   SparseSubMatrixHandler& WM = WorkMat.SetSparse();

   /* Dimensiona e resetta la matrice di lavoro */
   doublereal dDistance = pGetDriveCaller()->dGet();
   
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstVelIndex = iNode1FirstPosIndex+6;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstVelIndex = iNode2FirstPosIndex+6;
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+4;
   
   /* Distanza nulla */
      /* 
       *       forza:     dAlpha*v
       *       distanza:  x2 - x1 = d*v
       *       equazione di normalizzazione di v:
       *                  (v^T*v)^1/2 = 1
       *
       * 
       *        x1   x1P  x2   x2P  v    a    w    b
       * 
       * x1  |  0    0    0    0   -a*I -v    0    0 |
       * x1P |  0    0    0    0    0    0   -b*I -w |
       * x2  |  0    0    0    0    a*I  v    0    0 |
       * x2P |  0    0    0    0    0    0    b*I  w |
       * v   | -I    0    I    0   -d*I  0    0    0 |
       * a   |  0    0    0    0    r    0    0    0 |
       * w   |  0    0    0    0   -I    0    I    0 |
       * b   |  0   -w^T  0    w^T  0    0    s    0 |
       * 
       *                           v^T
       *                 r =    ----------
       *                        (vT*v)^1/2
       * 
       *                 s = (x2P - x1P)^T
       * 
       * a = dAlpha
       */
      
      WM.ResizeReset(51, 1);
      
      /* Nota: la direzione della reazione, v,
       * ed il coefficiente di amplificazione dAlpha
       * sono stati aggiornati durante il calcolo del residuo */
      
      dAlpha = XCurr(iFirstReactionIndex+4);
      doublereal dBeta = XCurr(iReactionPrimeIndex+4);
      v = Vec3(XCurr, iFirstReactionIndex+1);
      Vec3 w(XCurr, iReactionPrimeIndex+1);
      Vec3 deltaXP(pNode2->GetVCurr()-pNode1->GetVCurr());
      
      for (int iCnt = 1; iCnt <= 3; iCnt++) {
	 doublereal d;

	 d = dAlpha;
	 /* termini di Delta_v sul nodo 1 */
	 WM.PutItem(6+iCnt, iNode1FirstPosIndex+iCnt,
		     iFirstReactionIndex+iCnt, -d);
	 
	 /* termini di Delta_v sul nodo 2 */
	 WM.PutItem(9+iCnt, iNode2FirstPosIndex+iCnt,
		     iFirstReactionIndex+iCnt, d);
	 
	 d = v.dGet(iCnt);
	 /* termini di Delta_alpha sul nodo 1 */
	 WM.PutItem(12+iCnt, iNode1FirstPosIndex+iCnt,
		     iFirstReactionIndex+4, -d);
	 
	 /* termini di Delta_alpha sul nodo 2 */
	 WM.PutItem(15+iCnt, iNode2FirstPosIndex+iCnt,
		     iFirstReactionIndex+4, d);
	 
	 d = dBeta;
	 /* termini di Delta_w sul nodo 1 */
	 WM.PutItem(18+iCnt, iNode1FirstVelIndex+iCnt,
		     iReactionPrimeIndex+iCnt, -d);
	 
	 /* termini di Delta_w sul nodo 2 */
	 WM.PutItem(21+iCnt, iNode2FirstVelIndex+iCnt,
		     iReactionPrimeIndex+iCnt, d);
	 
	 d = w.dGet(iCnt);
	 /* termini di Delta_beta sul nodo 1 */
	 WM.PutItem(24+iCnt, iNode1FirstVelIndex+iCnt,
		     iReactionPrimeIndex+4, -d);
	 
	 /* termini di Delta_beta sul nodo 2 */
	 WM.PutItem(27+iCnt, iNode2FirstVelIndex+iCnt,
		     iReactionPrimeIndex+4, d);
	 
	 d = dDistance;
	 /* termini diagonali di Delta_v */
	 WM.PutItem(30+iCnt, iFirstReactionIndex+iCnt,
		     iFirstReactionIndex+iCnt, -d);
	 
	 /* termini diagonali di Delta_w */
	 WM.PutItem(33+iCnt, iReactionPrimeIndex+iCnt,
		     iReactionPrimeIndex+iCnt, 1.);
	 
	 /* termini incrociati di Delta_w Delta_v */
	 WM.PutItem(36+iCnt, iReactionPrimeIndex+iCnt,
		     iFirstReactionIndex+iCnt, -1.);
	 
	 d = w.dGet(iCnt);
	 /* termini di beta per Delta_x1P: w^T */
	 WM.PutItem(39+iCnt, iReactionPrimeIndex+4,
		     iNode1FirstVelIndex+iCnt, -d);
	 
	 /* termini di beta per Delta_x2P: w^T */
	 WM.PutItem(42+iCnt, iReactionPrimeIndex+4,
		     iNode2FirstVelIndex+iCnt, d);
	 
	 d = deltaXP.dGet(iCnt);
	 /* termini di beta per Delta_w: (x2P-x1P)^T */
	 WM.PutItem(45+iCnt, iReactionPrimeIndex+4,
		     iReactionPrimeIndex+iCnt, d);
      }
      
      
      /* termini di Delta_v su alpha */
      doublereal d = v.Dot();
      ASSERT(d > std::numeric_limits<doublereal>::epsilon());
      if (d > std::numeric_limits<doublereal>::epsilon()) {
	 d = std::sqrt(d);	     
	 for (int iCnt = 1; iCnt <= 3; iCnt++) {		
	    WM.PutItem(48+iCnt, iFirstReactionIndex+4,
			iFirstReactionIndex+iCnt, v.dGet(iCnt)/d);
	 }	     
      } else {	     
	 for (int iCnt = 1; iCnt <= 3; iCnt++) {	      
	    WM.PutItem(48+iCnt, iFirstReactionIndex+4,
			iFirstReactionIndex+iCnt, 0.);
	 }	   
      }		  				
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {      
      /* termini di Delta_x1 */
      WM.PutItem(iCnt, iFirstReactionIndex+iCnt,
		  iNode1FirstPosIndex+iCnt, -1.);
      
      /* termini di Delta_x2 */
      WM.PutItem(iCnt+3, iFirstReactionIndex+iCnt,
		  iNode2FirstPosIndex+iCnt, 1.);
   }   
   
   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
DistanceJoint::InitialAssRes(SubVectorHandler& WorkVec,
			     const VectorHandler& XCurr)
{   
   DEBUGCOUT("Entering DistanceJoint::InitialAssRes()" << std::endl);
   
   /* Dimensiona e resetta la matrice di lavoro */
   WorkVec.ResizeReset(20);
      
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstVelIndex = iNode1FirstPosIndex+6;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstVelIndex = iNode2FirstPosIndex+6;
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+4;
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {      
      /* Indici del nodo 1, posizione */
      WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      
      /* Indici del nodo 1, velocita' */
      WorkVec.PutRowIndex(3+iCnt, iNode1FirstVelIndex+iCnt);
      
      /* Indici del nodo 2, posizione */
      WorkVec.PutRowIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
      
      /* Indici del nodo 2, velocita' */
      WorkVec.PutRowIndex(9+iCnt, iNode2FirstVelIndex+iCnt);
   }
   
   /* Indici del vincolo */
   for (int iCnt = 1; iCnt <= 8; iCnt++) {      
      WorkVec.PutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
   }   
   
   Vec3 x1(pNode1->GetXCurr());
   Vec3 x2(pNode2->GetXCurr());
   Vec3 x1P(pNode1->GetVCurr());
   Vec3 x2P(pNode2->GetVCurr());
   v = Vec3(XCurr, iFirstReactionIndex+1);
   Vec3 w(XCurr, iReactionPrimeIndex+1);
   dAlpha = XCurr(iFirstReactionIndex+4);   
   doublereal dBeta = XCurr(iReactionPrimeIndex+4);

   doublereal dDistance = pGetDriveCaller()->dGet();
   
   /* Distanza nulla */
   if (fabs(dDistance) <= std::numeric_limits<doublereal>::epsilon()) {
      silent_cerr("DistanceJoint(" << GetLabel() << "): "
	      "near-zero distance" << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }
   WorkVec.Add(0+1, v*dAlpha);
   WorkVec.Add(3+1, w*dBeta);
   WorkVec.Add(6+1, -v*dAlpha);
   WorkVec.Add(9+1, -w*dBeta);
   WorkVec.Add(12+1, x1-x2+v*dDistance);
      
   doublereal d = v.Dot();
   ASSERT(d > std::numeric_limits<doublereal>::epsilon());
   if(d > std::numeric_limits<doublereal>::epsilon()) {	 
      d = std::sqrt(d);
   } else {	 
      d = 0.;
   }
      
   WorkVec.PutCoef(15+1, 1.-d);
      
   WorkVec.Add(16+1, v-w);
   WorkVec.PutCoef(19+1, (x1P-x2P).Dot(w));

   return WorkVec;
}


void DistanceJoint::SetInitialValue(VectorHandler& X)
{
   integer iFirstIndex = iGetFirstIndex();   
   
   doublereal dDistance = pGetDriveCaller()->dGet();
   if (fabs(dDistance) <= std::numeric_limits<doublereal>::epsilon()) {
      silent_cerr("DistanceJoint(" << GetLabel() << "): "
	      "near-zero distance" << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }
     
   v = ((pNode2->GetXCurr())-(pNode1->GetXCurr()));
   doublereal d = v.Dot();
   if (d < std::numeric_limits<doublereal>::epsilon()) {
      silent_cerr("DistanceJoint(" << GetLabel() << "): "
	      "initial length is null" << std::endl);
      throw ErrNullNorm(MBDYN_EXCEPT_ARGS);
   }
   v /= std::sqrt(d);
   
   X.Put(iFirstIndex+1, v);
   X.Put(iFirstIndex+4+1, v);   
}


void
DistanceJoint::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& /* XP */ ,
		SimulationEntity::Hints *ph)
{
	if (ph) {
		for (unsigned i = 0; i < ph->size(); i++) {
			DriveHint *pdh = dynamic_cast<DriveHint *>((*ph)[i]);

			if (pdh) {
				pedantic_cout("DistanceJoint(" << uLabel << "): "
					"creating drive from hint[" << i << "]..." << std::endl);

				DriveCaller *pDC = pdh->pCreateDrive(pDM);
				if (pDC == 0) {
					silent_cerr("DistanceJoint(" << uLabel << "): "
						"unable to create drive after hint "
						"#" << i << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				DriveOwner::Set(pDC);
				continue;
			}
		}
	}
	
	doublereal dDistance = pGetDriveCaller()->dGet();
   
	/* Setta a 1 dAlpha, che e' indipendente dalle altre variabili
	 * in caso di distanza nulla */
	if (fabs(dDistance) <= std::numeric_limits<doublereal>::epsilon()) {	
		silent_cerr("DistanceJoint(" << uLabel << "): "
			"near-zero distance" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* Scrive la direzione della distanza. Se e' stata ottenuta con 
	 * l'assemblaggio iniziale bene, se no' la calcola */
	v = pNode2->GetXCurr() - pNode1->GetXCurr();
	doublereal d = v.Dot();
	if (d <= std::numeric_limits<doublereal>::epsilon()) {
    		silent_cerr("DistanceJoint(" << uLabel << ") "
			"linked to nodes " << pNode1->GetLabel()
			<< " and " << pNode2->GetLabel() << ": "
			"nodes are coincident;" << std::endl
	  		<< "initial joint assembly is recommended"
			<< std::endl);
		throw ErrNullNorm(MBDYN_EXCEPT_ARGS);
	}
	v /= std::sqrt(d);	     

	/* Scrittura sul vettore soluzione */
	X.Put(iGetFirstIndex() + 1, v);	 
}


void 
DistanceJoint::GetDummyPartPos(unsigned int part, 
				 Vec3& x, 
				 Mat3x3& R) const 
{
   ASSERT(part == 1);
   x = pNode1->GetXCurr();
   R = pNode1->GetRCurr();
}

void 
DistanceJoint::GetDummyPartVel(unsigned int part, 
				 Vec3& v, 
				 Vec3& w) const 
{
   ASSERT(part == 1);
   v = pNode1->GetVCurr();
   w = pNode1->GetWCurr();
}

#ifdef USE_ADAMS
std::ostream& 
DistanceJoint::WriteAdamsDummyPartCmd(std::ostream& out,
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

/* DistanceJoint - end */


/* DistanceJointWithOffset - begin */

/* Costruttore non banale */
DistanceJointWithOffset::DistanceJointWithOffset(unsigned int uL, 
						 const DofOwner* pDO,
						 const StructNode* pN1, 
						 const StructNode* pN2,
						 const Vec3& f1Tmp, 
						 const Vec3& f2Tmp,
						 const DriveCaller* pDC,
						 flag fOut)
: Elem(uL, fOut), 
Joint(uL, pDO, fOut), 
DriveOwner(pDC),
pNode1(pN1), pNode2(pN2), f1(f1Tmp), f2(f2Tmp), v(Zero3), dAlpha(0.)
{
   ASSERT(pDO != NULL);
   ASSERT(pDC != NULL);
   ASSERT(pNode1 != NULL);
   ASSERT(pNode1->GetNodeType() == Node::STRUCTURAL);
   ASSERT(pNode2 != NULL);
   ASSERT(pNode2->GetNodeType() == Node::STRUCTURAL);
}


/* Distruttore banale - ci pensa il DriveOwner a distruggere il DriveCaller */
DistanceJointWithOffset::~DistanceJointWithOffset(void) 
{ 
   NO_OP;
}

/* Dati privati */
unsigned int
DistanceJointWithOffset::iGetNumPrivData(void) const
{
	return 1;
}

unsigned int
DistanceJointWithOffset::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);

	if (strcmp(s, "d") == 0) {
		return 1;
	}

	return 0;
}

doublereal
DistanceJointWithOffset::dGetPrivData(unsigned int i) const
{
	ASSERT(i == 1);

	if (i == 1) {
		return dGet();
	}

	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}


/* Contributo al file di restart */
std::ostream& DistanceJointWithOffset::Restart(std::ostream& out) const
{
   Joint::Restart(out) << ", distance with offset, " 
     << pNode1->GetLabel() 
     << ", reference, node, ",
     f1.Write(out, ", ") << ", "
     << pNode2->GetLabel() 
     << ", reference, node, ",
     f2.Write(out, ", ") << ", ";
   return pGetDriveCaller()->Restart(out) << ';' << std::endl;
}


/* Assemblaggio jacobiano */
VariableSubMatrixHandler& 
DistanceJointWithOffset::AssJac(VariableSubMatrixHandler& WorkMat,
				doublereal dCoef,
				const VectorHandler& /* XCurr */ ,
				const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering DistanceJointWithOffset::AssJac()" << std::endl);
      
   SparseSubMatrixHandler& WM = WorkMat.SetSparse();
   
   /* Dimensiona e resetta la matrice di lavoro */
   doublereal dDistance = pGetDriveCaller()->dGet();
   
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   
   Vec3 f1Tmp(pNode1->GetRRef()*f1);
   Vec3 f2Tmp(pNode2->GetRRef()*f2);
   
   /* Distanza nulla */
      WM.ResizeReset(72, 0);
            
      for (int iCnt = 1; iCnt <= 3; iCnt++) {
	 
	 /* termini di Delta_x1 */
	 WM.PutItem(iCnt, iFirstReactionIndex+iCnt,
		     iNode1FirstPosIndex+iCnt, -dCoef);
      
	 /* termini di Delta_x2 */
	 WM.PutItem(3+iCnt, iFirstReactionIndex+iCnt,
		     iNode2FirstPosIndex+iCnt, dCoef);
      
	 /* termini di Delta_F sul nodo 1 */
	 WM.PutItem(6+iCnt, iNode1FirstMomIndex+iCnt,
		     iFirstReactionIndex+iCnt, -dAlpha);
      
	 /* termini di Delta_F sul nodo 2 */
	 WM.PutItem(9+iCnt, iNode2FirstMomIndex+iCnt,
		     iFirstReactionIndex+iCnt, dAlpha);
      
	 /* termini diagonali di Delta_v */
	 WM.PutItem(12+iCnt, iFirstReactionIndex+iCnt,
		     iFirstReactionIndex+iCnt, -dDistance);
      
	 doublereal d = v.dGet(iCnt);
	 
	 /* termini di Delta_alpha sul nodo 1 */
	 WM.PutItem(15+iCnt, iNode1FirstMomIndex+iCnt,
		     iFirstReactionIndex+4, -d);
      
	 /* termini di Delta_alpha sul nodo 2 */
	 WM.PutItem(18+iCnt, iNode2FirstMomIndex+iCnt,
		     iFirstReactionIndex+4, d);
      }


      /* Termini di offset nell'equazione di vincolo */
      WM.PutCross(22, iFirstReactionIndex, 
		   iNode1FirstPosIndex+3, f1Tmp*dCoef);

      WM.PutCross(28, iFirstReactionIndex, 
		   iNode2FirstPosIndex+3, f2Tmp*(-dCoef));

      
      Vec3 Tmp(v*(dAlpha*dCoef));
      
      WM.PutMat3x3(34, iNode1FirstMomIndex+3,
		    iNode1FirstPosIndex+3, Mat3x3(MatCrossCross, Tmp, -f1Tmp));      
      WM.PutMat3x3(43, iNode2FirstMomIndex+3,
		    iNode2FirstPosIndex+3, Mat3x3(MatCrossCross, Tmp, f2Tmp));

      
      WM.PutCross(52, iNode1FirstMomIndex+3,
		   iFirstReactionIndex, f1Tmp*(-dAlpha));
      WM.PutCross(58, iNode2FirstMomIndex+3,
		   iFirstReactionIndex, f2Tmp*dAlpha);
      
      Tmp = v.Cross(f1Tmp);
      for (int iCnt = 1; iCnt <= 3; iCnt++) {
	 WM.PutItem(63+iCnt, iNode1FirstMomIndex+3+iCnt,
		     iFirstReactionIndex+4, Tmp.dGet(iCnt));
      }
      
      Tmp = f2Tmp.Cross(v);
      for (int iCnt = 1; iCnt <= 3; iCnt++) {
	 WM.PutItem(66+iCnt, iNode2FirstMomIndex+3+iCnt,
		     iFirstReactionIndex+4, Tmp.dGet(iCnt));
      }
      
      
      doublereal d = v.Dot();
      ASSERT(d > std::numeric_limits<doublereal>::epsilon());
      if (d > std::numeric_limits<doublereal>::epsilon()) {
	 d = std::sqrt(d);
      } else {
	 d = 1.;
      }      
      
      /* termini di Delta_v su alpha */
      for (int iCnt = 1; iCnt <= 3; iCnt++) {	 
	 WM.PutItem(69+iCnt, iFirstReactionIndex+4,
		     iFirstReactionIndex+iCnt, v.dGet(iCnt)/d);      
      }      
   
   return WorkMat;
}
   

/* Assemblaggio residuo */
SubVectorHandler& 
DistanceJointWithOffset::AssRes(SubVectorHandler& WorkVec,
				doublereal dCoef,
				const VectorHandler& XCurr, 
				const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering DistanceJointWithOffset::AssRes()" << std::endl);
      
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);
      
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();   
   
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      /* Indici del nodo 1 */
      WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.PutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
   }
   
   /* Indici del vincolo */
   for(int iCnt = 1; iCnt <= 4; iCnt++) {
      WorkVec.PutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
   }

   Vec3 x1(pNode1->GetXCurr());
   Vec3 x2(pNode2->GetXCurr());
   Vec3 f1Tmp(pNode1->GetRCurr()*f1);
   Vec3 f2Tmp(pNode2->GetRCurr()*f2);

   /* Aggiorna i dati propri */
   v = Vec3(XCurr, iFirstReactionIndex+1);
   dAlpha = XCurr(iFirstReactionIndex+4);   

   doublereal dDistance = pGetDriveCaller()->dGet();
   
   /* Distanza nulla */
   if (fabs(dDistance) <= std::numeric_limits<doublereal>::epsilon()) {	
      silent_cerr("DistanceJoint(" << GetLabel() << "): "
		      "near-zero distance" << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }
   Vec3 TmpVec(v*dAlpha);
   WorkVec.Add(1, TmpVec);
   WorkVec.Add(4, f1Tmp.Cross(TmpVec));
   WorkVec.Sub(7, TmpVec);	
   WorkVec.Sub(10, f2Tmp.Cross(TmpVec));
      
   WorkVec.Add(13, x1+f1Tmp-x2-f2Tmp+v*dDistance);
      
   doublereal dTmp = v.Dot();
   ASSERT(dTmp >= 0.0);
   dTmp = std::sqrt(dTmp);

   WorkVec.PutCoef(16, 1.-dTmp);

   return WorkVec;
}


void DistanceJointWithOffset::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) {
      doublereal d = dGet();
      Vec3 vTmp;      
      if (fabs(d) > std::numeric_limits<doublereal>::epsilon()) {
	 vTmp = Vec3(dAlpha, 0., 0.);
      } else {      
	 vTmp = v;
      }
      Joint::Output(OH.Joints(), "DistanceWithOffs", GetLabel(),
		    vTmp, Zero3, v*dAlpha, Zero3)
	<< " " << v << " " << d << std::endl;
   }   
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
DistanceJointWithOffset::InitialAssJac(VariableSubMatrixHandler& WorkMat,
				       const VectorHandler& XCurr)
{
   DEBUGCOUT("Entering DistanceJointWithOffset::InitialAssJac()" << std::endl);

   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.ResizeReset(32, 32);

   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+4;
   
   for (int iCnt = 1; iCnt <= 12; iCnt++) {
      WM.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutRowIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
      WM.PutColIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
   }
   
   for (int iCnt = 1; iCnt <= 8; iCnt++) {
      WM.PutRowIndex(24+iCnt, iFirstReactionIndex+iCnt);
      WM.PutColIndex(24+iCnt, iFirstReactionIndex+iCnt);
   }
        
   doublereal dAlphaP = XCurr(iReactionPrimeIndex+4);  
   Vec3 vP(XCurr, iReactionPrimeIndex+1);   
   
   Vec3 f1Tmp(pNode1->GetRRef()*f1);
   Vec3 f2Tmp(pNode2->GetRRef()*f2);
   Vec3 Omega1(pNode1->GetWRef());
   Vec3 Omega2(pNode2->GetWRef());   
   
   doublereal dDistance = pGetDriveCaller()->dGet();
   
      /* Equazioni di equilibrio */
      for (int iCnt = 1; iCnt <= 3; iCnt++) {
	 WM.PutCoef(iCnt, 24+iCnt, -dAlpha);
	 WM.PutCoef(12+iCnt, 24+iCnt, dAlpha);

	 WM.PutCoef(6+iCnt, 28+iCnt, -dAlpha);
	 WM.PutCoef(18+iCnt, 28+iCnt, dAlpha);
	 
	 WM.PutCoef(6+iCnt, 24+iCnt, -dAlphaP);
	 WM.PutCoef(18+iCnt, 24+iCnt, dAlphaP);
	 
	 doublereal d = v.dGet(iCnt);
	 WM.PutCoef(iCnt, 28, -d);
	 WM.PutCoef(12+iCnt, 28, d);

	 WM.PutCoef(6+iCnt, 32, -d);
	 WM.PutCoef(18+iCnt, 32, d);
	 
	 d = vP.dGet(iCnt);
	 WM.PutCoef(6+iCnt, 28, -d);
	 WM.PutCoef(18+iCnt, 28, d);
	 
      }
      
      Vec3 Tmp(v.Cross(f1Tmp));
      for (int iCnt = 1; iCnt <= 3; iCnt++) {
	 doublereal d = Tmp.dGet(iCnt);
	 WM.PutCoef(3+iCnt, 28, d);
	 WM.PutCoef(9+iCnt, 32, d);
      }
      
      Tmp = f2Tmp.Cross(v);
      for (int iCnt = 1; iCnt <= 3; iCnt++) {
	 doublereal d = Tmp.dGet(iCnt);
	 WM.PutCoef(15+iCnt, 28, d);
	 WM.PutCoef(21+iCnt, 32, d);
      }
      
      Tmp = v.Cross(Omega1.Cross(f1Tmp))+vP.Cross(f1Tmp);
      for (int iCnt = 1; iCnt <= 3; iCnt++) {
	 doublereal d = Tmp.dGet(iCnt);
	 WM.PutCoef(9+iCnt, 28, d);
      }
      
      Tmp = v.Cross(f2Tmp.Cross(Omega2))-vP.Cross(f2Tmp);
      for (int iCnt = 1; iCnt <= 3; iCnt++) {
	 doublereal d = Tmp.dGet(iCnt);
	 WM.PutCoef(21+iCnt, 28, d);
      }
      
      
      Tmp = f1Tmp*(-dAlpha);
      
      Mat3x3 MTmp(MatCrossCross, v, Tmp);
      WM.Add(4, 4, MTmp);
      WM.Add(10, 10, MTmp);
      
      MTmp = Mat3x3(MatCross, Tmp);
      WM.Add(4, 25, MTmp);
      WM.Add(10, 29, MTmp);
      
      Tmp = f2Tmp*dAlpha;
      MTmp = Mat3x3(MatCrossCross, v, Tmp);
      WM.Add(16, 16, MTmp);
      WM.Add(22, 22, MTmp);
      
      MTmp = Mat3x3(MatCross, Tmp);
      WM.Add(16, 25, MTmp);
      WM.Add(22, 29, MTmp);
      
      MTmp = (Mat3x3(MatCrossCross, v*dAlpha, Omega1) + Mat3x3(MatCross, vP*dAlpha + v*dAlphaP))*Mat3x3(MatCross, f1Tmp);
      WM.Sub(10, 4, MTmp);
      MTmp = Mat3x3(MatCross, f1Tmp.Cross(Omega1*dAlpha) - f1Tmp*dAlphaP);
      WM.Add(10, 25, MTmp);
      
      MTmp = (Mat3x3(MatCrossCross, v*dAlpha, Omega2) + Mat3x3(MatCross, vP*dAlpha + v*dAlphaP))*Mat3x3(MatCross, f2Tmp);
      WM.Add(22, 16, MTmp);
      MTmp = Mat3x3(MatCross, f2Tmp.Cross(Omega2*(-dAlpha)) + f2Tmp*dAlphaP);
      WM.Add(22, 25, MTmp);            

      /* Equazioni di vincolo */
      for (int iCnt = 1; iCnt <= 3; iCnt++) {
	 WM.PutCoef(24+iCnt, iCnt, -1.);
	 WM.PutCoef(24+iCnt, 12+iCnt, 1.);
	 
	 WM.PutCoef(24+iCnt, 24+iCnt, -dDistance);

	 WM.PutCoef(28+iCnt, 6+iCnt, -1.);
	 WM.PutCoef(28+iCnt, 18+iCnt, 1.);
      
	 WM.PutCoef(28+iCnt, 28+iCnt, -dDistance);
      }
      
      MTmp = Mat3x3(MatCross, f1Tmp);
      WM.Add(25, 4, MTmp);
      WM.Add(29, 10, MTmp);
      
      MTmp = Mat3x3(MatCross, f2Tmp);
      WM.Sub(25, 16, MTmp);
      WM.Sub(29, 22, MTmp);
      
      MTmp = Mat3x3(MatCrossCross, Omega1, f1Tmp);
      WM.Add(29, 4, MTmp);
      
      MTmp = Mat3x3(MatCrossCross, Omega2, f2Tmp);
      WM.Sub(29, 16, MTmp);
      
      doublereal d = v.Dot();
      ASSERT(d > std::numeric_limits<doublereal>::epsilon());
      if (d > std::numeric_limits<doublereal>::epsilon()) {
	 d = std::sqrt(d);
	 
	 Tmp = v/d;
	 for (int iCnt = 1; iCnt <= 3; iCnt++) {
	    doublereal d = Tmp.dGet(iCnt);
	    WM.PutCoef(28, 24+iCnt, d);
	    WM.PutCoef(32, 28+iCnt, d);
	 }
	 
	 Tmp = vP/d-v*((vP.Dot(v))/pow(d, 3));
	 for (int iCnt = 1; iCnt <= 3; iCnt++) {
	    doublereal d = Tmp.dGet(iCnt);
	    WM.PutCoef(32, 24+iCnt, d);
	 }	 
      }                        
   
   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
DistanceJointWithOffset::InitialAssRes(SubVectorHandler& WorkVec,
				       const VectorHandler& XCurr)
{   
   DEBUGCOUT("Entering DistanceJointWithOffset::InitialAssRes()" << std::endl);
   
   /* Dimensiona e resetta la matrice di lavoro */
   WorkVec.ResizeReset(32);
      
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+4;

   for (int iCnt = 1; iCnt <= 12; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WorkVec.PutRowIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
   }
   
   for (int iCnt = 1; iCnt <= 8; iCnt++) {
      WorkVec.PutRowIndex(24+iCnt, iFirstReactionIndex+iCnt);
   }

   v = Vec3(XCurr, iFirstReactionIndex+1);
   Vec3 vP(XCurr, iReactionPrimeIndex+1);
   dAlpha = XCurr(iFirstReactionIndex+4);   
   doublereal dAlphaP = XCurr(iReactionPrimeIndex+4);   

   Vec3 x1(pNode1->GetXCurr());
   Vec3 x2(pNode2->GetXCurr());
   Vec3 v1(pNode1->GetVCurr());
   Vec3 v2(pNode2->GetVCurr());
   Vec3 f1Tmp(pNode1->GetRCurr()*f1);
   Vec3 f2Tmp(pNode2->GetRCurr()*f2);
   Vec3 Omega1(pNode1->GetWCurr());
   Vec3 Omega2(pNode2->GetWCurr());   
   
   doublereal dDistance = pGetDriveCaller()->dGet();
   
   if (fabs(dDistance) <= std::numeric_limits<doublereal>::epsilon()) {
      silent_cerr("DistanceJoint(" << GetLabel() << "): "
	      "near-zero distance" << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }
   Vec3 Tmp(v*dAlpha);
   WorkVec.Put(1, Tmp);
   WorkVec.Put(4, f1Tmp.Cross(Tmp));
   WorkVec.Put(13, -Tmp);
   WorkVec.Put(16, Tmp.Cross(f2Tmp));
      
   Tmp = vP*dAlpha+v*dAlphaP;
   WorkVec.Put(7, Tmp);
   WorkVec.Put(10, (Omega1.Cross(f1Tmp)).Cross(v*dAlpha)+f1Tmp.Cross(Tmp));
   WorkVec.Put(19, -Tmp);
   WorkVec.Put(22, (f2Tmp.Cross(Omega2)).Cross(v*dAlpha)+Tmp.Cross(f2Tmp));
      
   WorkVec.Put(25, v*dDistance-x2-f2Tmp+x1+f1Tmp);
   WorkVec.Put(29, vP*dDistance
		  -v2-Omega2.Cross(f2Tmp)+v1+Omega1.Cross(f1Tmp));
      
   doublereal d = v.Dot();
      
   ASSERT(d > std::numeric_limits<doublereal>::epsilon());
   d = std::sqrt(d);
   WorkVec.PutCoef(28, 1.-d);
   WorkVec.PutCoef(32, -v.Dot(vP)/d);
   
   return WorkVec;
}


void DistanceJointWithOffset::SetInitialValue(VectorHandler& X)
{
   integer iFirstIndex = iGetFirstIndex();   
   
   doublereal dDistance = pGetDriveCaller()->dGet();
   if (fabs(dDistance) <= std::numeric_limits<doublereal>::epsilon()) {
      silent_cerr("DistanceJoint(" << GetLabel() << "): "
	      "near-zero distance" << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }
   
   Vec3 x1(pNode1->GetXCurr());
   Vec3 x2(pNode2->GetXCurr());
   Vec3 v1(pNode1->GetVCurr());
   Vec3 v2(pNode2->GetVCurr());
   Vec3 f1Tmp((pNode1->GetRCurr())*f1);
   Vec3 f2Tmp((pNode2->GetRCurr())*f2);
   Vec3 Omega1(pNode1->GetWCurr());
   Vec3 Omega2(pNode2->GetWCurr());   
   
   v = x2+f2Tmp-x1-f1Tmp;
   doublereal d = v.Dot();
   if (d < std::numeric_limits<doublereal>::epsilon()) {
      silent_cerr("DistanceJoint(" << GetLabel() << "): "
	      "initial length is null" << std::endl);
      throw ErrNullNorm(MBDYN_EXCEPT_ARGS);
   }
   v /= std::sqrt(d);
   
   X.Put(iFirstIndex+1, v);
   X.Put(iFirstIndex+5, (v2+Omega2.Cross(f2Tmp)-v1-Omega1.Cross(f1Tmp))/d);
}


void
DistanceJointWithOffset::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& /* XP */ ,
		SimulationEntity::Hints *ph)
{
	if (ph) {
		for (unsigned i = 0; i < ph->size(); i++) {
			pedantic_cout("DistanceJointWithOffset(" << uLabel << "): "
				"creating drive from hint..." << std::endl);

			DriveHint *pdh = dynamic_cast<DriveHint *>((*ph)[i]);

			if (pdh) {
				DriveCaller *pDC = pdh->pCreateDrive(pDM);
				if (pDC == 0) {
					silent_cerr("DistanceJointWithOffset(" << uLabel << "): "
						"unable to create drive after hint "
						"#" << i << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				DriveOwner::Set(pDC);
				continue;
			}
		}
	}
	
	doublereal dDistance = pGetDriveCaller()->dGet();

	/* Setta a 1 dAlpha, che e' indipendente dalle altre variabili
	 * in caso di distanza nulla */
	if (fabs(dDistance) <= std::numeric_limits<doublereal>::epsilon()) {	
		silent_cerr("DistanceJoint(" << GetLabel() << "):"
			"near-zero distance" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	v = pNode2->GetXCurr() + pNode2->GetRCurr()*f2 - pNode1->GetXCurr() - pNode1->GetRCurr()*f1;
	doublereal d = v.Dot();
	if (d <= std::numeric_limits<doublereal>::epsilon()) {
		silent_cerr("DistanceJoint(" << GetLabel() << ") "
			"linked to nodes " << pNode1->GetLabel()
			<< " and " << pNode2->GetLabel() << ": "
			"nodes are coincident;" << std::endl
			<< "this is no longer supported" << std::endl);
		throw ErrNullNorm(MBDYN_EXCEPT_ARGS);
	}
	v /= std::sqrt(d);

	/* Scrittura sul vettore soluzione */
	X.Put(iGetFirstIndex() + 1, v);	 
}

void 
DistanceJointWithOffset::GetDummyPartPos(unsigned int part,
					   Vec3& x, 
					   Mat3x3& R) const 
{
   ASSERT(part == 1);
   x = pNode1->GetXCurr()+pNode1->GetRCurr()*f1;
   
   Vec3 x2 = pNode2->GetXCurr()+pNode2->GetRCurr()*f2;
   
   Vec3 v1 = x2-x;
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
   
   R = MatR2vec(1, v1, 2, v2);
   
   // R = pNode1->GetRCurr();
}

void 
DistanceJointWithOffset::GetDummyPartVel(unsigned int part,
					   Vec3& v, 
					   Vec3& w) const 
{
   ASSERT(part == 1);

#if 0
   x = pNode1->GetXCurr()+pNode1->GetRCurr()*f1;
   
   Vec3 x2 = pNode2->GetXCurr()+pNode2->GetRCurr()*f2;
   
   Vec3 v1 = x2-x;
   doublereal l = v1.Norm();
   v1 /= l;
   
   Mat3x3 Rx(Eye3-v1.Tens(v1));
   int index = 1;
   if (fabs(v1(2)) < fabs(v1(index))) {
      index = 2;
   }
   if (fabs(v1(3)) < fabs(v1(index))) {
      index = 3;
   }
   
   Vec3 v2(Rx.GetVec(index));
   v2 /= v2.Norm();
   
   R = MatR2vec(1, v1, 2, v2);
#endif

   w = pNode1->GetWCurr();
   v = pNode1->GetVCurr()+w.Cross(pNode1->GetRCurr()*f1);
}


#ifdef USE_ADAMS
std::ostream& 
DistanceJointWithOffset::WriteAdamsDummyPartCmd(std::ostream& out,
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
     << e /* MatR2EulerAngles(pNode1->GetRCurr())*dRaDegr */ << " "
     << x1 << " "
     << e << " "
     << l << " " << 0. << " " << 0. << " "
     << Zero3 << std::endl;
}
#endif /* USE_ADAMS */

/* DistanceJointWithOffset - end */

#endif


/* ClampJoint - begin */

/* Costruttore definitivo (da mettere a punto) */
ClampJoint::ClampJoint(unsigned int uL, const DofOwner*pD, 
		       const StructNode* pN, 
		       const Vec3& X0, const Mat3x3& R0, 
		       flag fOut)
: Elem(uL, fOut), 
Joint(uL, pD, fOut), 
pNode(pN), XClamp(X0), RClamp(R0), F(Zero3), M(Zero3)
{ 
   NO_OP; 
}


ClampJoint::~ClampJoint(void) 
{ 
   NO_OP; 
}


std::ostream&
ClampJoint::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
			"reaction forces [Fx,Fy,Fz]" << std::endl
		<< prefix << iIndex + 4 << "->" << iIndex + 6 << ": "
			"reaction couples [mx,my,mz]" << std::endl;
	
	if (bInitial) {
		iIndex += 6;
		out
			<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
				"reaction force derivatives [FPx,FPy,FPz]" << std::endl
			<< prefix << iIndex + 4 << "->" << iIndex + 6 << ": "
				"reaction couple derivatives [mPx,mPy,mPz]" << std::endl;
	}

	return out;
}

static const char xyz[] = "xyz";
static const char *dof[] = {
	"reaction force f",
	"reaction couple m",
	"reaction force derivative fP",
	"reaction couple derivative mP"
};
static const char *eq[] = {
	"position constraint P",
	"orientation constraint theta",
	"position constraint derivative v",
	"orientation constraint derivative w"
};

void
ClampJoint::DescribeDof(std::vector<std::string>& desc,
	bool bInitial, int i) const
{
	int iend = 1;
	if (i == -1) {
		if (bInitial) {
			iend = 12;

		} else {
			iend = 6;
		}
	}
	desc.resize(iend);

	std::ostringstream os;
	os << "ClampJoint(" << GetLabel() << ")";

	if (i == -1) {
		std::string name = os.str();
		for (i = 0; i < iend; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": " << dof[i/3] << xyz[i%3];

			desc[i] = os.str();
		}

	} else {
		os << ": " << dof[i/3] << xyz[i%3];
		desc[0] = os.str();
	}
}

std::ostream&
ClampJoint::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
			"position constraints [Px=Px0,Py=Py0,Pz=Pz0]" << std::endl
		<< prefix << iIndex + 4 << "->" << iIndex + 6 << ": "
			"orientation constraints [thetax=thetax0,thetay=thetay0,thetaz=thetaz0]" << std::endl;
	
	if (bInitial) {
		iIndex += 6;
		out
			<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": "
				"velocity constraints [vx=0,vy=0,vz=0]" << std::endl
			<< prefix << iIndex + 4 << "->" << iIndex + 6 << ": "
				"angular velocity constraints [wx=0,wy=0,wz=0]" << std::endl;
	}

	return out;
}

void
ClampJoint::DescribeEq(std::vector<std::string>& desc,
	bool bInitial, int i) const
{
	int iend = 1;
	if (i == -1) {
		if (bInitial) {
			iend = 12;

		} else {
			iend = 6;
		}
	}
	desc.resize(iend);

	std::ostringstream os;
	os << "ClampJoint(" << GetLabel() << ")";

	if (i == -1) {
		std::string name = os.str();
		for (i = 0; i < iend; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": " << eq[i/3] << xyz[i%3];

			desc[i] = os.str();
		}

	} else {
		os << ": " << eq[i/3] << xyz[i%3];
		desc[0] = os.str();
	}
}

/*Funzione che legge lo stato iniziale dal file di input*/
void
ClampJoint::ReadInitialState(MBDynParser& HP)
{
	F = HP.GetVec3();
	M = HP.GetVec3();
}

/* Contributo al file di restart */
std::ostream&
ClampJoint::Restart(std::ostream& out) const
{
	return Joint::Restart(out) << ", clamp, "
		<< pNode->GetLabel() << ", node, node, "
		<< "initial state, ", F.Write(out, ", ")
		<< ", ", M.Write(out, ", ") << ';' << std::endl;
}


VariableSubMatrixHandler& 
ClampJoint::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering ClampJoint::AssJac()" << std::endl);
      
	SparseSubMatrixHandler& WM = WorkMat.SetSparse();
   
	/* Dimensiona e resetta la matrice di lavoro */
	WM.ResizeReset(12, 1);

	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
	integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
	integer iFirstReactionIndex = iGetFirstIndex();
   
	/* Attenzione: modifico dividendo le equazioni di vincolo per dCoef */
	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutItem(iCnt, iFirstReactionIndex + iCnt, 
			iFirstPositionIndex + iCnt, 1.);    
		WM.PutItem(6 + iCnt, iFirstMomentumIndex + iCnt,
			iFirstReactionIndex + iCnt, 1.);
	}
 
	/* Con l'aggiunta dei nodi statici non occorre piu' evitare la
	 * singolarita' della matrice */

	return WorkMat;
}

/* Inverse Dynamics */
VariableSubMatrixHandler& 
ClampJoint::AssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUT("Entering ClampJoint::AssJac()" << std::endl);
   
  	SparseSubMatrixHandler& WM = WorkMat.SetSparse();
   
   	/* Dimensiona e resetta la matrice di lavoro */
   	WM.ResizeReset(6, 1);
   
   	integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
   	integer iFirstReactionIndex = iGetFirstIndex();
   	
	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
      		WM.PutItem(iCnt, iFirstReactionIndex + iCnt, 
			iFirstPositionIndex + iCnt, 1.);    
   	}
        
	return WorkMat;
}

/* assemblaggio matrici per autovalori */
void
ClampJoint::AssMats(VariableSubMatrixHandler& WorkMatA,
	VariableSubMatrixHandler& WorkMatB,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering ClampJoint::AssMats(); will result in call to AssJac()"
		<< std::endl);
   
	WorkMatA.SetNullMatrix();
	AssJac(WorkMatB, 1., XCurr, XPrimeCurr);
}


SubVectorHandler&
ClampJoint::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering ClampJoint::AssRes()" << std::endl);

	WorkVec.ResizeReset(12);

	/* Indici delle incognite del nodo e delle reazioni */
	integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
	integer iFirstReactionIndex = iGetFirstIndex();
	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iFirstReactionIndex + iCnt);
	}   

	/* Aggiorna le reazioni vincolari */
	F = Vec3(XCurr, iFirstReactionIndex + 1);
	M = Vec3(XCurr, iFirstReactionIndex + 3 + 1);

	/* Calcola posizione e parametri di rotazione */
	const Vec3& x(pNode->GetXCurr());
	const Mat3x3& R(pNode->GetRCurr());

	Vec3 theta_c(RotManip::VecRot(R.MulMT(RClamp)));

	/* Residuo della riga di equilibrio */
	WorkVec.Sub(1, F);
	WorkVec.Sub(3 + 1, M);
 
	/* Modifica: divido le equazioni di vincolo per dCoef */
	if (dCoef != 0.) {	
		/* Residuo dell'equazione di vincolo */
		WorkVec.Sub(6 + 1, (x - XClamp)/dCoef);
		WorkVec.Sub(9 + 1, theta_c/dCoef);   
	}

	return WorkVec;
}

/* Inverse Dynamics: AssRes()*/
SubVectorHandler&
ClampJoint::AssRes(SubVectorHandler& WorkVec,
	const VectorHandler& XCurr, 
	const VectorHandler& /* XPrimeCurr */,
	const VectorHandler& /* XPrimePrimeCurr */,
	InverseDynamics::Order iOrder)
{
   	DEBUGCOUT("Entering ClampJoint::AssRes()" << std::endl);

	/* The residual is != 0 only for position */
   	if (iOrder == InverseDynamics::POSITION) {
   		WorkVec.ResizeReset(6);

		/* FIXME: Indici delle incognite */
		integer iFirstReactionIndex = iGetFirstIndex();
		for (integer iCnt = 1; iCnt <= 6; iCnt++) {
			WorkVec.PutRowIndex(iCnt, iFirstReactionIndex+iCnt);
		}   
 
		/* Calcola posizione e parametri di rotazione */
		const Vec3& x(pNode->GetXCurr());
		const Mat3x3& R(pNode->GetRCurr());

		Vec3 theta_c(RotManip::VecRot(R.MulMT(RClamp)));

		/* Residuo dell'equazione di vincolo */
		WorkVec.Sub(1, x - XClamp);
		WorkVec.Sub(3 + 1, theta_c);   

	} else {
   		WorkVec.Resize(0);
	}

	return WorkVec;
}

/* Inverse Dynamics update */
void 
ClampJoint::Update(const VectorHandler& XCurr, InverseDynamics::Order iOrder)
{
	ASSERT(iOrder == InverseDynamics::INVERSE_DYNAMICS);

	integer iFirstReactionIndex = iGetFirstIndex();

   	/* Aggiorna le reazioni vincolari */
   	F = Vec3(XCurr, iFirstReactionIndex + 1);
   	M = Vec3(XCurr, iFirstReactionIndex + 4);
}

void
ClampJoint::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		const Mat3x3& R(pNode->GetRCurr());

		Joint::Output(OH.Joints(), "Clamp", GetLabel(),
			R.MulTV(F), R.MulTV(M), F, M) << std::endl;
	}
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
ClampJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& /* XCurr */ )
{
   DEBUGCOUT("Entering ClampJoint::InitialAssJac()" << std::endl);
   
   SparseSubMatrixHandler& WM = WorkMat.SetSparse();
   WM.ResizeReset(24, 1);
   
   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
   integer iFirstVelocityIndex = iFirstPositionIndex+6;
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+6;
   
   /* 
    * x  |  0   0   0   0  -I   0   0   0 || Delta_x  |   |  f    |
    * g  |  0   0   0   0   0  -I   0   0 || Delta_g  |   |  m    |
    * xP |  0   0   0   0   0   0  -I   0 || Delta_xP |   |  fP   |
    * w  |  0   0   0   0   0   0   0  -I || Delta_w  |   |  mP   |
    * f  |  I   0   0   0   0   0   0   0 || Delta_f  | = |  X0-x |
    * m  |  0   I   0   0   0   0   0   0 || Delta_m  |   | -xP   |
    * fP |  0   0   I   0   0   0   0   0 || Delta_fP |   | -g    |
    * mP |  0   0   0   I   0   0   0   0 || Delta_mP |   | -w    |
    *    
    * 
    * 
    */
   
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.PutItem(iCnt, iFirstReactionIndex+iCnt, 
		  iFirstPositionIndex+iCnt, 1.);
   
      WM.PutItem(6+iCnt, iReactionPrimeIndex+iCnt, 
		  iFirstVelocityIndex+iCnt, 1.);
   
      WM.PutItem(12+iCnt, iFirstPositionIndex+iCnt, 
		  iFirstReactionIndex+iCnt, -1.);
   
      WM.PutItem(18+iCnt, iFirstVelocityIndex+iCnt, 
		  iReactionPrimeIndex+iCnt, -1.);
   }
   
   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
ClampJoint::InitialAssRes(SubVectorHandler& WorkVec,
			  const VectorHandler& XCurr)
{   
   DEBUGCOUT("Entering ClampJoint::InitialAssRes()" << std::endl);
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);
      
   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();  
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+6;

   for (int iCnt = 1; iCnt <= 12; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iFirstPositionIndex+iCnt);
      WorkVec.PutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
   }   
   
   
   /* Forza */
   WorkVec.Put(1, Vec3(XCurr, iFirstReactionIndex+1));
   
   /* Coppia */
   WorkVec.Put(4, Vec3(XCurr, iFirstReactionIndex+4));
   
   /* Derivata della Forza */
   WorkVec.Put(7, Vec3(XCurr, iReactionPrimeIndex+1));
   
   /* Derivata della Coppia */
   WorkVec.Put(10, Vec3(XCurr, iReactionPrimeIndex+4));
   
   /* Posizione */
   WorkVec.Put(13, XClamp-pNode->GetXCurr());
   
   /* Parametri di rotazione; 
    * si sfrutta il fatto che g(R_Delta^T) = -g(R_Delta) */
   Mat3x3 R(pNode->GetRCurr());
   WorkVec.Put(16, Vec3(CGR_Rot::Param, RClamp.MulMT(R)));
   
   /* Velocita' */
   WorkVec.Put(19, -pNode->GetVCurr());
   
   /* Velocita' angolare */
   WorkVec.Put(22, -pNode->GetWCurr());
         
   return WorkVec;
}

void
ClampJoint::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
	integer iFirstReactionIndex = iGetFirstIndex();
	X.Put(iFirstReactionIndex+1,F);
	X.Put(iFirstReactionIndex+4,M);
}

/* Metodi per l'estrazione di dati "privati".
 * Si suppone che l'estrattore li sappia interpretare.
 * Come default non ci sono dati privati estraibili */
unsigned int
ClampJoint::iGetNumPrivData(void) const
{
	return 6;
}

unsigned int
ClampJoint::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);

	unsigned int idx = 0;

	switch (s[0]) {
	case 'F':
		break;

	case 'M':
		idx += 3;
		break;

	default:
		return 0;
	}

	switch (s[1]) {
	case 'x':
		idx += 1;
		break;

	case 'y':
		idx += 2;
		break;

	case 'z':
		idx += 3;
		break;

	default:
		return 0;
	}

	if (s[2] != '\0') {
		return 0;
	}

	return idx;
}

doublereal
ClampJoint::dGetPrivData(unsigned int i) const
{
	if (i >= 1 && i <= 3) {
		return F(i);
	}

	if (i >= 4 && i <= 6) {
		return M(i - 3);
	}

	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

/* ClampJoint - end */
