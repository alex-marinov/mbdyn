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

/* Vincoli generali */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <genj.h>

/* DistanceJoint - begin */

/* Costruttore non banale */
DistanceJoint::DistanceJoint(unsigned int uL, const DofOwner* pDO,
			     const StructNode* pN1, const StructNode* pN2,
			     const DriveCaller* pDC, flag fOut)
: Elem(uL, Elem::JOINT, fOut), 
Joint(uL, Joint::DISTANCE, pDO, fOut),
DriveOwner(pDC),
pNode1(pN1), pNode2(pN2), v(0.), dAlpha(0.)
{
   NO_OP;
}


/* Distruttore banale - ci pensa il DriveOwner a distruggere il DriveCaller */
DistanceJoint::~DistanceJoint(void) 
{ 
   NO_OP; 
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
   
   /* Distanza nulla */
   if (fabs(dDistance) <= DBL_EPSILON) {
      /* 
       *      forza:             v
       *      distanza (nulla):  x2 - x1 = 0
       *      equazione banale per dAlpha:
       *                         dAlpha = 0
       * 
       * 
       *        x1   Q1   x2   Q2   v    a
       * 
       * x1  |  0    0    0    0    0    0 |
       * Q1  |  0    0    0    0   -I    0 |
       * x2  |  0    0    0    0    0    0 |
       * Q2  |  0    0    0    0    I    0 |
       * v   | -c*I  0    c*I  0    0    0 |
       * a   |  0    0    0    0    0    1 |
       *
       * con c = dCoef
       * 
       * In realta' modifico lo jacobiano dividendo le equazioni 
       * di vincolo per dCoef dove possibile, in modo da migliorare
       * il condizionamento della matrice. */
      
      
      
      WM.ResizeInit(13, 1, 0.);
      
	
      /* Attenzione: modifico jacobiano e residuo dividendo per dCoef le
       * equazioni di vincolo (solo per d == 0) */
      
      for (int iCnt = 1; iCnt <= 3; iCnt++) {	   
	 /* termini di Delta_x1 */
	 WM.fPutItem(iCnt, iFirstReactionIndex+iCnt,
		     iNode1FirstPosIndex+iCnt, -1.);
	 
	 /* termini di Delta_x2 */
	 WM.fPutItem(3+iCnt, iFirstReactionIndex+iCnt,
		     iNode2FirstPosIndex+iCnt, 1.);	 
	 
	 /* termini di Delta_F sul nodo 1 */
	 WM.fPutItem(6+iCnt, iNode1FirstMomIndex+iCnt,
		     iFirstReactionIndex+iCnt, -1.);
	 
	 /* termini di Delta_F sul nodo 2 */
	 WM.fPutItem(9+iCnt, iNode2FirstMomIndex+iCnt,
		     iFirstReactionIndex+iCnt, 1.);
      }
      
      /* Termine diagonale per alpha */
      WM.fPutItem(13, iFirstReactionIndex+4,
		  iFirstReactionIndex+4, 1.);
      
   } else {
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
      
      WM.ResizeInit(24, 0, 0.);
      
      /* Nota: la direzione della reazione, v,
       * ed il coefficiente di amplificazione dAlpha
       * sono stati aggiornati durante il calcolo del residuo */
      
      for (int iCnt = 1; iCnt <= 3; iCnt++) {	   
	 /* termini di Delta_x1 */
	 WM.fPutItem(iCnt, iFirstReactionIndex+iCnt,
		     iNode1FirstPosIndex+iCnt, -dCoef);
	 
	 /* termini di Delta_x2 */
	 WM.fPutItem(3+iCnt, iFirstReactionIndex+iCnt,
		     iNode2FirstPosIndex+iCnt, dCoef);
	 
	 /* termini di Delta_F sul nodo 1 */
	 WM.fPutItem(6+iCnt, iNode1FirstMomIndex+iCnt,
		     iFirstReactionIndex+iCnt, -dAlpha);
	 
	 /* termini di Delta_F sul nodo 2 */
	 WM.fPutItem(9+iCnt, iNode2FirstMomIndex+iCnt,
		     iFirstReactionIndex+iCnt, dAlpha);
	 
	 /* termini diagonali di Delta_v */
	 WM.fPutItem(12+iCnt, iFirstReactionIndex+iCnt,
		     iFirstReactionIndex+iCnt, -dDistance);
	 
	 /* termini di Delta_alpha sul nodo 1 */
	 WM.fPutItem(15+iCnt, iNode1FirstMomIndex+iCnt,
		     iFirstReactionIndex+4, -v.dGet(iCnt));
	 
	 /* termini di Delta_alpha sul nodo 2 */
	 WM.fPutItem(18+iCnt, iNode2FirstMomIndex+iCnt,
		     iFirstReactionIndex+4, v.dGet(iCnt));
      }	
      
      doublereal d = v.Dot();
      ASSERT(d > DBL_EPSILON);
      if (d > DBL_EPSILON) {	   
	 d = sqrt(d);
	 /* termini di Delta_v su alpha */
	 for (int iCnt = 3; iCnt > 0; iCnt--) {	      
	    WM.fPutItem(21+iCnt, iFirstReactionIndex+4,
			iFirstReactionIndex+iCnt, v.dGet(iCnt)/d);
	 }	   				
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
   this->WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
        
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();  
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {      
      /* Indici del nodo 1 */
      WorkVec.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
   
      /* Indici del nodo 2 */
      WorkVec.fPutRowIndex(3+iCnt, iNode2FirstMomIndex+iCnt);
   }
   
   /* Indici del vincolo */
   for (int iCnt = 1; iCnt <= 4; iCnt++) {      
      WorkVec.fPutRowIndex(6+iCnt, iFirstReactionIndex+iCnt);
   }   

   Vec3 x1(pNode1->GetXCurr());
   Vec3 x2(pNode2->GetXCurr());
   
   /* Aggiorna i dati propri */
   v = Vec3(XCurr, iFirstReactionIndex+1);
   dAlpha = XCurr.dGetCoef(iFirstReactionIndex+4);   

   doublereal dDistance = pGetDriveCaller()->dGet();
   
   /* Distanza nulla */
   if (fabs(dDistance) <= DBL_EPSILON) {	
      WorkVec.Add(1, v);
      WorkVec.Sub(4, v);
      
      /* Modifica: se dCoef non e' nullo (caso normale), divido il residuo 
       * delle equazioni di vincolo (solo lo prime 3) per dCoef, altrimenti
       * lascio nullo il vincolo. Infatti dCoef = 0 significa che sta 
       * facendo il passo fittizio iniziale, nel quale si assume che 
       * le equazioni di vincolo siano soddisfatte */
      if (dCoef >= DBL_EPSILON) {	     
	 WorkVec.Add(7, (x1-x2)/dCoef);
      }
      
      WorkVec.fPutCoef(10, 1.-dAlpha);
      
   } else {	
      Vec3 TmpVec(v*dAlpha);
      WorkVec.Add(1, TmpVec);
      WorkVec.Sub(4, TmpVec);	
      
      WorkVec.Add(7, x1-x2+v*dDistance);
      
      doublereal d = v.Dot();
      ASSERT(d > DBL_EPSILON);
      if (d > DBL_EPSILON) {	   
	 d = sqrt(d);
      } else {	   
	 d = 0.;
      }
      
      WorkVec.fPutCoef(10, 1.-d);	  
   }
   
   return WorkVec;
}


void DistanceJoint::Output(OutputHandler& OH) const
{
   if(fToBeOutput()) {      
      doublereal d = dGet();
      Vec3 vTmp;      
      if (fabs(d) > DBL_EPSILON) {
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
   if (fabs(dDistance) <= DBL_EPSILON) {
      /* 
       *      forza:             v
       *      distanza (nulla):  x2 - x1 = 0
       *      equazione banale per dAlpha:
       *                         dAlpha = 1
       * 
       *      derivata della forza:    w
       *      derivata della distanza: xP2 - xP1 = 0
       *      derivata di alpha:       dAlphaP = 0
       * 
       * 
       *        x1   xP1  x2   xP2  v    w   a   aP 
       * 
       * x1  |  0    0    0    0   -I    0    0    0 || Delta_x1  |
       * xP1 |  0    0    0    0    0    0   -I    0 || Delta_xP1 |
       * x2  |  0    0    0    0    I    0    0    0 || Delta_x2  |
       * xP2 |  0    0    0    0    0    0    I    0 || Delta_xP2 |
       * v   | -I    0    I    0    0    0    0    0 || Delta_v   |
       * a   |  0    0    0    0    0    1    0    0 || Delta_a   |
       * w   |  0   -I    0    I    0    0    0    0 || Delta_w   |
       * aP  |  0    0    0    0    0    0    0    1 || Delta_aP  |
       * 
       * 
       */
	
      
      WM.ResizeInit(26, 1, 0.);
      
      for (int iCnt = 1; iCnt <= 3; iCnt++) {	   
	 /* derivata dell'equazione di vincolo, nodo 1 */
	 WM.fPutItem(6+iCnt, iReactionPrimeIndex+iCnt,
		     iNode1FirstVelIndex+iCnt, -1.);
	 
	 /* derivata dell'equazione di vincolo, nodo 2 */
	 WM.fPutItem(9+iCnt, iReactionPrimeIndex+iCnt,
		     iNode2FirstVelIndex+iCnt, 1.);
	 
	 /* termini di forza sul nodo 1 */
	 WM.fPutItem(12+iCnt, iNode1FirstPosIndex+iCnt,
		     iFirstReactionIndex+iCnt, -1.);
	 
	 /* termini di forza sul nodo 2 */
	 WM.fPutItem(15+iCnt, iNode2FirstPosIndex+iCnt,
		     iFirstReactionIndex+iCnt, 1.);
	
	 /* termini di derivata della forza sul nodo 1 */
	 WM.fPutItem(18+iCnt, iNode1FirstVelIndex+iCnt,
		     iReactionPrimeIndex+iCnt, -1.);
	 
	 /* termini di derivata della forza sul nodo 2 */
	 WM.fPutItem(21+iCnt, iNode2FirstVelIndex+iCnt,
		     iReactionPrimeIndex+iCnt, 1.);
      }
      
      /* Termine diagonale per alpha */
      WM.fPutItem(25, iFirstReactionIndex+4,
		  iFirstReactionIndex+4, 1.);
      
      /* Termine diagonale per alpha derivato */
      WM.fPutItem(26, iReactionPrimeIndex+4,
		  iReactionPrimeIndex+4, 1.);
      
   } else {
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
      
      WM.ResizeInit(51, 1, 0.);
      
      /* Nota: la direzione della reazione, v,
       * ed il coefficiente di amplificazione dAlpha
       * sono stati aggiornati durante il calcolo del residuo */
      
      dAlpha = XCurr.dGetCoef(iFirstReactionIndex+4);
      doublereal dBeta = XCurr.dGetCoef(iReactionPrimeIndex+4);
      v = Vec3(XCurr, iFirstReactionIndex+1);
      Vec3 w(XCurr, iReactionPrimeIndex+1);
      Vec3 deltaXP(pNode2->GetVCurr()-pNode1->GetVCurr());
      
      for (int iCnt = 1; iCnt <= 3; iCnt++) {	   
	 /* termini di Delta_v sul nodo 1 */
	 WM.fPutItem(6+iCnt, iNode1FirstPosIndex+iCnt,
		     iFirstReactionIndex+iCnt, -dAlpha);
	 
	 /* termini di Delta_v sul nodo 2 */
	 WM.fPutItem(9+iCnt, iNode2FirstPosIndex+iCnt,
		     iFirstReactionIndex+iCnt, dAlpha);
	 
	 /* termini di Delta_alpha sul nodo 1 */
	 WM.fPutItem(12+iCnt, iNode1FirstPosIndex+iCnt,
		     iFirstReactionIndex+4, -v.dGet(iCnt));
	 
	 /* termini di Delta_alpha sul nodo 2 */
	 WM.fPutItem(15+iCnt, iNode2FirstPosIndex+iCnt,
		     iFirstReactionIndex+4, v.dGet(iCnt));
	 
	 /* termini di Delta_w sul nodo 1 */
	 WM.fPutItem(18+iCnt, iNode1FirstVelIndex+iCnt,
		     iReactionPrimeIndex+iCnt, -dBeta);
	 
	 /* termini di Delta_w sul nodo 2 */
	 WM.fPutItem(21+iCnt, iNode2FirstVelIndex+iCnt,
		     iReactionPrimeIndex+iCnt, dBeta);
	 
	 /* termini di Delta_beta sul nodo 1 */
	 WM.fPutItem(24+iCnt, iNode1FirstVelIndex+iCnt,
		     iReactionPrimeIndex+4, -w.dGet(iCnt));
	 
	 /* termini di Delta_beta sul nodo 2 */
	 WM.fPutItem(27+iCnt, iNode2FirstVelIndex+iCnt,
		     iReactionPrimeIndex+4, w.dGet(iCnt));
	 
	 /* termini diagonali di Delta_v */
	 WM.fPutItem(30+iCnt, iFirstReactionIndex+iCnt,
		     iFirstReactionIndex+iCnt, -dDistance);
	 
	 /* termini diagonali di Delta_w */
	 WM.fPutItem(33+iCnt, iReactionPrimeIndex+iCnt,
		     iReactionPrimeIndex+iCnt, 1.);
	 
	 /* termini incrociati di Delta_w Delta_v */
	 WM.fPutItem(36+iCnt, iReactionPrimeIndex+iCnt,
		     iFirstReactionIndex+iCnt, -1.);
	 
	 /* termini di beta per Delta_x1P: w^T */
	 WM.fPutItem(39+iCnt, iReactionPrimeIndex+4,
		     iNode1FirstVelIndex+iCnt, -w.dGet(iCnt));
	 
	 /* termini di beta per Delta_x2P: w^T */
	 WM.fPutItem(42+iCnt, iReactionPrimeIndex+4,
		     iNode2FirstVelIndex+iCnt, w.dGet(iCnt));
	 
	 /* termini di beta per Delta_w: (x2P-x1P)^T */
	 WM.fPutItem(45+iCnt, iReactionPrimeIndex+4,
		     iReactionPrimeIndex+iCnt, deltaXP.dGet(iCnt));
      }
      
      
      /* termini di Delta_v su alpha */
      doublereal d = v.Dot();
      ASSERT(d > DBL_EPSILON);
      if (d > DBL_EPSILON) {
	 d = sqrt(d);	     
	 for (int iCnt = 1; iCnt <= 3; iCnt++) {		
	    WM.fPutItem(48+iCnt, iFirstReactionIndex+4,
			iFirstReactionIndex+iCnt, v.dGet(iCnt)/d);
	 }	     
      } else {	     
	 for (int iCnt = 1; iCnt <= 3; iCnt++) {	      
	    WM.fPutItem(48+iCnt, iFirstReactionIndex+4,
			iFirstReactionIndex+iCnt, 0.);
	 }	   
      }		  				
   }
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {      
      /* termini di Delta_x1 */
      WM.fPutItem(iCnt, iFirstReactionIndex+iCnt,
		  iNode1FirstPosIndex+iCnt, -1.);
      
      /* termini di Delta_x2 */
      WM.fPutItem(iCnt+3, iFirstReactionIndex+iCnt,
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
   WorkVec.Resize(20);
   WorkVec.Reset(0.);
      
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstVelIndex = iNode1FirstPosIndex+6;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstVelIndex = iNode2FirstPosIndex+6;
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+4;
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {      
      /* Indici del nodo 1, posizione */
      WorkVec.fPutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      
      /* Indici del nodo 1, velocita' */
      WorkVec.fPutRowIndex(3+iCnt, iNode1FirstVelIndex+iCnt);
      
      /* Indici del nodo 2, posizione */
      WorkVec.fPutRowIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
      
      /* Indici del nodo 2, velocita' */
      WorkVec.fPutRowIndex(9+iCnt, iNode2FirstVelIndex+iCnt);
   }
   
   /* Indici del vincolo */
   for (int iCnt = 1; iCnt <= 8; iCnt++) {      
      WorkVec.fPutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
   }   
   
   Vec3 x1(pNode1->GetXCurr());
   Vec3 x2(pNode2->GetXCurr());
   Vec3 x1P(pNode1->GetVCurr());
   Vec3 x2P(pNode2->GetVCurr());
   v = Vec3(XCurr, iFirstReactionIndex+1);
   Vec3 w(XCurr, iReactionPrimeIndex+1);
   dAlpha = XCurr.dGetCoef(iFirstReactionIndex+4);   
   doublereal dBeta = XCurr.dGetCoef(iReactionPrimeIndex+4);

   doublereal dDistance = pGetDriveCaller()->dGet();
   
   /* Distanza nulla */
   if (fabs(dDistance) <= DBL_EPSILON) {
      WorkVec.Add(1, v);
      WorkVec.Add(4, w);
      WorkVec.Add(7, -v);
      WorkVec.Add(10, -w);
      WorkVec.Add(13, x1-x2);
      WorkVec.fPutCoef(16, 1.-dAlpha);
      WorkVec.Add(17, x1P-x2P);
      WorkVec.fPutCoef(20, -dBeta);
   } else {			
      WorkVec.Add(1, v*dAlpha);
      WorkVec.Add(4, w*dBeta);
      WorkVec.Add(7, -v*dAlpha);
      WorkVec.Add(10, -w*dBeta);
      WorkVec.Add(13, x1-x2+v*dDistance);
      
      doublereal d = v.Dot();
      ASSERT(d > DBL_EPSILON);
      if(d > DBL_EPSILON) {	 
	d = sqrt(d);
      } else {	 
	d = 0.;
      }
      
      WorkVec.fPutCoef(16, 1.-d);
      
      WorkVec.Add(17, v-w);
      WorkVec.fPutCoef(20, (x1P-x2P).Dot(w));
   }

   return WorkVec;
}


void DistanceJoint::SetInitialValue(VectorHandler& X) const
{
   integer iFirstIndex = iGetFirstIndex();   
   
   doublereal dDistance = pGetDriveCaller()->dGet();
   if (fabs(dDistance) <= DBL_EPSILON) {
      X.Put(iFirstIndex+4, 1.);
      return;
   }
     
   (Vec3&)v = ((pNode2->GetXCurr())-(pNode1->GetXCurr()));
   doublereal d = v.Dot();
   ASSERT(d > DBL_EPSILON);
   if (d > DBL_EPSILON) {
      (Vec3&)v /= d;
   }     
   
   X.Put(iFirstIndex+1, v);
   X.Put(iFirstIndex+5, v);   
}


void DistanceJoint::SetValue(VectorHandler& X, VectorHandler& /* XP */ ) const
{
   doublereal dDistance = pGetDriveCaller()->dGet();
   
   /* Setta a 1 dAlpha, che e' indipendente dalle altre variabili
    * in caso di distanza nulla */
   if (fabs(dDistance) <= DBL_EPSILON) {	
      X.Put(iGetFirstIndex()+4, 1.);
   } else {      
   
   /* Scrive la direzione della distanza. Se e' stata ottenuta con 
    * l'assemblaggio iniziale bene, se no' la calcola */
      doublereal d = v.Dot();
      if (d <= DBL_EPSILON) {
	 (Vec3&)v = (pNode2->GetXCurr()-pNode1->GetXCurr());
	 d = v.Dot();
	 if(d <= DBL_EPSILON) {
	    std::cerr << "Joint " << uLabel << ", linked to nodes " 
	      << pNode1->GetLabel() << " and " << pNode2->GetLabel() 
	      << ": nodes are coincident." << std::endl
	      << "Initial joint assembly is recommended; aborting ... " 
	      << std::endl;
	    
	    THROW(ErrGeneric());
	 }
	 (Vec3&)v /= sqrt(d);	     
      }

      /* Scrittura sul vettore soluzione */
      X.Put(iGetFirstIndex()+1, v);	 
   }   
}


void 
DistanceJoint::GetAdamsDummyPart(unsigned int part, 
				 Vec3& x, 
				 Mat3x3& R) const 
{
   ASSERT(part == 1);
   x = pNode1->GetXCurr();
   R = pNode1->GetRCurr();
}


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
   
   Vec3 e(MatR2EulerAngles(MatR2vec(1, v1, 2, v2)));
   
   return out 
     << psAdamsElemCode[GetElemType()] << "_" << GetLabel() << "_" << part << std::endl
     << firstId << " "
     << x1 << " "
     << MatR2EulerAngles(pNode1->GetRCurr()) << " "
     << x1 << " "
     << e << " "
     << l << " " << 0. << " " << 0. << " "
     << Zero3 << std::endl;
}

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
: Elem(uL, Elem::JOINT, fOut), 
Joint(uL, Joint::DISTANCEWITHOFFSET, pDO, fOut), 
DriveOwner(pDC),
pNode1(pN1), pNode2(pN2), f1(f1Tmp), f2(f2Tmp), v(0.), dAlpha(0.)
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
   if (fabs(dDistance) <= DBL_EPSILON) {            
      
      WM.ResizeInit(55, 0, 0.);
      
      /* Attenzione: per d == 0. divido per dCoef */
            
      for (int iCnt = 1; iCnt <= 3; iCnt++) {	 

	 /* termini di Delta_x1 */
	 WM.fPutItem(iCnt, iFirstReactionIndex+iCnt,
		     iNode1FirstPosIndex+iCnt, -1.);
      
	 /* termini di Delta_x2 */
	 WM.fPutItem(3+iCnt, iFirstReactionIndex+iCnt,
		     iNode2FirstPosIndex+iCnt, 1.);
      
      
	 /* termini di Delta_F sul nodo 1 */
	 WM.fPutItem(6+iCnt, iNode1FirstMomIndex+iCnt,
		     iFirstReactionIndex+iCnt, -1.);
      
	 /* termini di Delta_F sul nodo 2 */
	 WM.fPutItem(9+iCnt, iNode2FirstMomIndex+iCnt,
		     iFirstReactionIndex+iCnt, 1.);
      }
      
      /* Termini di offset nell'equazione di vincolo */
      WM.fPutCross(13, iFirstReactionIndex, 
		   iNode1FirstPosIndex+3, f1Tmp);

      WM.fPutCross(19, iFirstReactionIndex, 
		   iNode2FirstPosIndex+3, -f2Tmp);
            
      /* Termini di momento nell'equazione di equilibrio nodo 1 */
      WM.fPutCross(25, iNode1FirstMomIndex+3,
		   iFirstReactionIndex, -f1Tmp);
      WM.fPutMat3x3(31, iNode1FirstMomIndex+3,
		    iNode1FirstPosIndex+3, Mat3x3(v, f1Tmp*(-dCoef)));

      
      /* Termini di momento nell'equazione di equilibrio nodo 2 */
      WM.fPutCross(40, iNode2FirstMomIndex+3,
		   iFirstReactionIndex, f2Tmp);
      WM.fPutMat3x3(46, iNode2FirstMomIndex+3,
		    iNode2FirstPosIndex+3, Mat3x3(v, f2Tmp*dCoef));
      

      /* Termine diagonale per alpha */
      WM.fPutItem(55, iFirstReactionIndex+4,
		  iFirstReactionIndex+4, 1.);
      
   } else {
      
      WM.ResizeInit(72, 0, 0.);
            
      for (int iCnt = 1; iCnt <= 3; iCnt++) {
	 
	 /* termini di Delta_x1 */
	 WM.fPutItem(iCnt, iFirstReactionIndex+iCnt,
		     iNode1FirstPosIndex+iCnt, -dCoef);
      
	 /* termini di Delta_x2 */
	 WM.fPutItem(3+iCnt, iFirstReactionIndex+iCnt,
		     iNode2FirstPosIndex+iCnt, dCoef);
      
	 /* termini di Delta_F sul nodo 1 */
	 WM.fPutItem(6+iCnt, iNode1FirstMomIndex+iCnt,
		     iFirstReactionIndex+iCnt, -dAlpha);
      
	 /* termini di Delta_F sul nodo 2 */
	 WM.fPutItem(9+iCnt, iNode2FirstMomIndex+iCnt,
		     iFirstReactionIndex+iCnt, dAlpha);
      
	 /* termini diagonali di Delta_v */
	 WM.fPutItem(12+iCnt, iFirstReactionIndex+iCnt,
		     iFirstReactionIndex+iCnt, -dDistance);
      
	 doublereal d = v.dGet(iCnt);
	 
	 /* termini di Delta_alpha sul nodo 1 */
	 WM.fPutItem(15+iCnt, iNode1FirstMomIndex+iCnt,
		     iFirstReactionIndex+4, -d);
      
	 /* termini di Delta_alpha sul nodo 2 */
	 WM.fPutItem(18+iCnt, iNode2FirstMomIndex+iCnt,
		     iFirstReactionIndex+4, d);
      }


      /* Termini di offset nell'equazione di vincolo */
      WM.fPutCross(22, iFirstReactionIndex, 
		   iNode1FirstPosIndex+3, f1Tmp*dCoef);

      WM.fPutCross(28, iFirstReactionIndex, 
		   iNode2FirstPosIndex+3, f2Tmp*(-dCoef));

      
      Vec3 Tmp(v*(dAlpha*dCoef));
      
      WM.fPutMat3x3(34, iNode1FirstMomIndex+3,
		    iNode1FirstPosIndex+3, Mat3x3(Tmp, -f1Tmp));      
      WM.fPutMat3x3(43, iNode2FirstMomIndex+3,
		    iNode2FirstPosIndex+3, Mat3x3(Tmp, f2Tmp));

      
      WM.fPutCross(52, iNode1FirstMomIndex+3,
		   iFirstReactionIndex, f1Tmp*(-dAlpha));
      WM.fPutCross(58, iNode2FirstMomIndex+3,
		   iFirstReactionIndex, f2Tmp*dAlpha);
      
      Tmp = v.Cross(f1Tmp);
      for (int iCnt = 1; iCnt <= 3; iCnt++) {
	 WM.fPutItem(63+iCnt, iNode1FirstMomIndex+3+iCnt,
		     iFirstReactionIndex+4, Tmp.dGet(iCnt));
      }
      
      Tmp = f2Tmp.Cross(v);
      for (int iCnt = 1; iCnt <= 3; iCnt++) {
	 WM.fPutItem(66+iCnt, iNode2FirstMomIndex+3+iCnt,
		     iFirstReactionIndex+4, Tmp.dGet(iCnt));
      }
      
      
      doublereal d = v.Dot();
      ASSERT(d > DBL_EPSILON);
      if (d > DBL_EPSILON) {
	 d = sqrt(d);
      } else {
	 d = 1.;
      }      
      
      /* termini di Delta_v su alpha */
      for (int iCnt = 1; iCnt <= 3; iCnt++) {	 
	 WM.fPutItem(69+iCnt, iFirstReactionIndex+4,
		     iFirstReactionIndex+iCnt, v.dGet(iCnt)/d);      
      }      
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
   this->WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
      
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();   
   
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      /* Indici del nodo 1 */
      WorkVec.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.fPutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
   }
   
   /* Indici del vincolo */
   for(int iCnt = 1; iCnt <= 4; iCnt++) {
      WorkVec.fPutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
   }

   Vec3 x1(pNode1->GetXCurr());
   Vec3 x2(pNode2->GetXCurr());
   Vec3 f1Tmp(pNode1->GetRCurr()*f1);
   Vec3 f2Tmp(pNode2->GetRCurr()*f2);

   /* Aggiorna i dati propri */
   v = Vec3(XCurr, iFirstReactionIndex+1);
   dAlpha = XCurr.dGetCoef(iFirstReactionIndex+4);   

   doublereal dDistance = pGetDriveCaller()->dGet();
   
   /* Distanza nulla */
   if (fabs(dDistance) <= DBL_EPSILON) {	
      WorkVec.Add(1, v);
      WorkVec.Add(4, f1Tmp.Cross(v));
      WorkVec.Sub(7, v);
      WorkVec.Sub(10, f2Tmp.Cross(v));
      
      /* Modifica: se dCoef non e' nullo (caso normale), divido il residuo 
       * delle equazioni di vincolo (solo lo prime 3) per dCoef, altrimenti
       * lascio nullo il vincolo. Infatti dCoef = 0 significa che sta 
       * facendo il passo fittizio iniziale, nel quale si assume che 
       * le equazioni di vincolo siano soddisfatte */
      if (fabs(dCoef) > DBL_EPSILON) {	     
	 WorkVec.Add(13, (x1+f1Tmp-x2-f2Tmp)/dCoef);
      }
	
      WorkVec.fPutCoef(16, 1.-dAlpha);
      
   } else {	
      Vec3 TmpVec(v*dAlpha);
      WorkVec.Add(1, TmpVec);
      WorkVec.Add(4, f1Tmp.Cross(TmpVec));
      WorkVec.Sub(7, TmpVec);	
      WorkVec.Sub(10, f2Tmp.Cross(TmpVec));
      
      WorkVec.Add(13, x1+f1Tmp-x2-f2Tmp+v*dDistance);
      
      doublereal dTmp = v.Dot();
      ASSERT(dTmp >= 0.0);
      if (dTmp > 0.) {	   
	 dTmp = sqrt(dTmp);
      } else {	   
	 dTmp = 0.;
      }	
      WorkVec.fPutCoef(16, 1.-dTmp);
      
   }
   
   return WorkVec;
}


void DistanceJointWithOffset::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) {
      doublereal d = dGet();
      Vec3 vTmp;      
      if (fabs(d) > DBL_EPSILON) {
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
   WM.ResizeInit(32, 32, 0.);

   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+4;
   
   for (int iCnt = 1; iCnt <= 12; iCnt++) {
      WM.fPutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.fPutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.fPutRowIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
      WM.fPutColIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
   }
   
   for (int iCnt = 1; iCnt <= 8; iCnt++) {
      WM.fPutRowIndex(24+iCnt, iFirstReactionIndex+iCnt);
      WM.fPutColIndex(24+iCnt, iFirstReactionIndex+iCnt);
   }
        
   doublereal dAlphaP = XCurr.dGetCoef(iReactionPrimeIndex+4);  
   Vec3 vP(XCurr, iReactionPrimeIndex+1);   
   
   Vec3 f1Tmp(pNode1->GetRRef()*f1);
   Vec3 f2Tmp(pNode2->GetRRef()*f2);
   Vec3 Omega1(pNode1->GetWRef());
   Vec3 Omega2(pNode2->GetWRef());   
   
   doublereal dDistance = pGetDriveCaller()->dGet();
   
   if (fabs(dDistance) <= DBL_EPSILON) {
      for (int iCnt = 1; iCnt <= 3; iCnt++) {	 
	 WM.fPutCoef(iCnt, 24+iCnt, -1.);
	 WM.fPutCoef(12+iCnt, 24+iCnt, 1.);
	 
	 WM.fPutCoef(6+iCnt, 28+iCnt, -1.);
	 WM.fPutCoef(18+iCnt, 28+iCnt, 1.);
	 
	 WM.fPutCoef(24+iCnt, iCnt, -1.);
	 WM.fPutCoef(24+iCnt, 12+iCnt, 1.);

	 WM.fPutCoef(28+iCnt, 6+iCnt, -1.);
	 WM.fPutCoef(28+iCnt, 18+iCnt, 1.);
      }
      
      Mat3x3 MTmp(-f1Tmp);
      WM.Add(4, 25, MTmp);
      WM.Add(10, 29, MTmp);

      MTmp = Mat3x3(f1Tmp);
      WM.Add(25, 4, MTmp);
      WM.Add(29, 10, MTmp);
      
      MTmp = Mat3x3(f2Tmp);
      WM.Add(16, 25, MTmp);
      WM.Add(22, 29, MTmp);

      MTmp = Mat3x3(-f2Tmp);
      WM.Add(25, 16, MTmp);
      WM.Add(29, 22, MTmp);
      
      MTmp = Mat3x3(-v, f1Tmp);
      WM.Add(4, 4, MTmp);
      WM.Add(10, 10, MTmp);
      
      MTmp = Mat3x3(v, f2Tmp);
      WM.Add(16, 16, MTmp);
      WM.Add(22, 22, MTmp);
      
      MTmp = Mat3x3(Omega1, f1Tmp);
      WM.Add(29, 4, MTmp);
      
      MTmp = Mat3x3(-Omega2, f2Tmp);
      WM.Add(29, 16, MTmp);
      
      MTmp = Mat3x3(f1Tmp.Cross(Omega1));
      WM.Add(10, 25, MTmp);
      
      MTmp = Mat3x3(Omega2.Cross(f2Tmp));
      WM.Add(22, 25, MTmp);
      
      MTmp = (Mat3x3(v, Omega1)+Mat3x3(vP))*Mat3x3(-f1Tmp);
      WM.Add(10, 4, MTmp);
      
      MTmp = (Mat3x3(v, Omega2)+Mat3x3(vP))*Mat3x3(f2Tmp);
      WM.Add(22, 16, MTmp);
      
      WM.fPutCoef(28, 28, 1.);
      WM.fPutCoef(32, 32, 1.);      
      
   } else {
      
      /* Equazioni di equilibrio */
      for (int iCnt = 1; iCnt <= 3; iCnt++) {
	 WM.fPutCoef(iCnt, 24+iCnt, -dAlpha);
	 WM.fPutCoef(12+iCnt, 24+iCnt, dAlpha);

	 WM.fPutCoef(6+iCnt, 28+iCnt, -dAlpha);
	 WM.fPutCoef(18+iCnt, 28+iCnt, dAlpha);
	 
	 WM.fPutCoef(6+iCnt, 24+iCnt, -dAlphaP);
	 WM.fPutCoef(18+iCnt, 24+iCnt, dAlphaP);
	 
	 doublereal d = v.dGet(iCnt);
	 WM.fPutCoef(iCnt, 28, -d);
	 WM.fPutCoef(12+iCnt, 28, d);

	 WM.fPutCoef(6+iCnt, 32, -d);
	 WM.fPutCoef(18+iCnt, 32, d);
	 
	 d = vP.dGet(iCnt);
	 WM.fPutCoef(6+iCnt, 28, -d);
	 WM.fPutCoef(18+iCnt, 28, d);
	 
      }
      
      Vec3 Tmp(v.Cross(f1Tmp));
      for (int iCnt = 1; iCnt <= 3; iCnt++) {
	 doublereal d = Tmp.dGet(iCnt);
	 WM.fPutCoef(3+iCnt, 28, d);
	 WM.fPutCoef(9+iCnt, 32, d);
      }
      
      Tmp = f2Tmp.Cross(v);
      for (int iCnt = 1; iCnt <= 3; iCnt++) {
	 doublereal d = Tmp.dGet(iCnt);
	 WM.fPutCoef(15+iCnt, 28, d);
	 WM.fPutCoef(21+iCnt, 32, d);
      }
      
      Tmp = v.Cross(Omega1.Cross(f1Tmp))+vP.Cross(f1Tmp);
      for (int iCnt = 1; iCnt <= 3; iCnt++) {
	 doublereal d = Tmp.dGet(iCnt);
	 WM.fPutCoef(9+iCnt, 28, d);
      }
      
      Tmp = v.Cross(f2Tmp.Cross(Omega2))-vP.Cross(f2Tmp);
      for (int iCnt = 1; iCnt <= 3; iCnt++) {
	 doublereal d = Tmp.dGet(iCnt);
	 WM.fPutCoef(21+iCnt, 28, d);
      }
      
      
      Tmp = f1Tmp*(-dAlpha);
      
      Mat3x3 MTmp(v, Tmp);
      WM.Add(4, 4, MTmp);
      WM.Add(10, 10, MTmp);
      
      MTmp = Mat3x3(Tmp);
      WM.Add(4, 25, MTmp);
      WM.Add(10, 29, MTmp);
      
      Tmp = f2Tmp*dAlpha;
      MTmp = Mat3x3(v, Tmp);
      WM.Add(16, 16, MTmp);
      WM.Add(22, 22, MTmp);
      
      MTmp = Mat3x3(Tmp);
      WM.Add(16, 25, MTmp);
      WM.Add(22, 29, MTmp);
      
      MTmp = (Mat3x3(v*dAlpha, Omega1)+
	      Mat3x3(vP*dAlpha+v*dAlphaP))*Mat3x3(-f1Tmp);
      WM.Add(10, 4, MTmp);
      MTmp = Mat3x3(f1Tmp.Cross(Omega1*dAlpha)-f1Tmp*dAlphaP);
      WM.Add(10, 25, MTmp);
      
      MTmp = (Mat3x3(v*dAlpha, Omega2)+
	      Mat3x3(vP*dAlpha+v*dAlphaP))*Mat3x3(f2Tmp);
      WM.Add(22, 16, MTmp);
      MTmp = Mat3x3(f2Tmp.Cross(Omega2*(-dAlpha))+f2Tmp*dAlphaP);
      WM.Add(22, 25, MTmp);            

      /* Equazioni di vincolo */
      for (int iCnt = 1; iCnt <= 3; iCnt++) {
	 WM.fPutCoef(24+iCnt, iCnt, -1.);
	 WM.fPutCoef(24+iCnt, 12+iCnt, 1.);
	 
	 WM.fPutCoef(24+iCnt, 24+iCnt, -dDistance);

	 WM.fPutCoef(28+iCnt, 6+iCnt, -1.);
	 WM.fPutCoef(28+iCnt, 18+iCnt, 1.);
      
	 WM.fPutCoef(28+iCnt, 28+iCnt, -dDistance);
      }
      
      MTmp = Mat3x3(f1Tmp);
      WM.Add(25, 4, MTmp);
      WM.Add(29, 10, MTmp);
      
      MTmp = Mat3x3(-f2Tmp);
      WM.Add(25, 16, MTmp);
      WM.Add(29, 22, MTmp);
      
      MTmp = Mat3x3(Omega1, f1Tmp);
      WM.Add(29, 4, MTmp);
      
      MTmp = Mat3x3(Omega2, -f2Tmp);
      WM.Add(29, 16, MTmp);
      
      doublereal d = v.Dot();
      ASSERT(d > DBL_EPSILON);
      if (d > DBL_EPSILON) {
	 d = sqrt(d);
	 
	 Tmp = v/d;
	 for (int iCnt = 1; iCnt <= 3; iCnt++) {
	    doublereal d = Tmp.dGet(iCnt);
	    WM.fPutCoef(28, 24+iCnt, d);
	    WM.fPutCoef(32, 28+iCnt, d);
	 }
	 
	 Tmp = vP/d-v*((vP.Dot(v))/pow(d, 3));
	 for (int iCnt = 1; iCnt <= 3; iCnt++) {
	    doublereal d = Tmp.dGet(iCnt);
	    WM.fPutCoef(32, 24+iCnt, d);
	 }	 
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
   WorkVec.Resize(32);
   WorkVec.Reset(0.);
      
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+4;

   for (int iCnt = 1; iCnt <= 12; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WorkVec.fPutRowIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
   }
   
   for (int iCnt = 1; iCnt <= 8; iCnt++) {
      WorkVec.fPutRowIndex(24+iCnt, iFirstReactionIndex+iCnt);
   }

   v = Vec3(XCurr, iFirstReactionIndex+1);
   Vec3 vP(XCurr, iReactionPrimeIndex+1);
   dAlpha = XCurr.dGetCoef(iFirstReactionIndex+4);   
   doublereal dAlphaP = XCurr.dGetCoef(iReactionPrimeIndex+4);   

   Vec3 x1(pNode1->GetXCurr());
   Vec3 x2(pNode2->GetXCurr());
   Vec3 v1(pNode1->GetVCurr());
   Vec3 v2(pNode2->GetVCurr());
   Vec3 f1Tmp(pNode1->GetRCurr()*f1);
   Vec3 f2Tmp(pNode2->GetRCurr()*f2);
   Vec3 Omega1(pNode1->GetWCurr());
   Vec3 Omega2(pNode2->GetWCurr());   
   
   doublereal dDistance = pGetDriveCaller()->dGet();
   
   if (fabs(dDistance) <= DBL_EPSILON) {
      WorkVec.Put(1, v);
      WorkVec.Put(4, f1Tmp.Cross(v));
      WorkVec.Put(7, vP);
      WorkVec.Put(10, (Omega1.Cross(f1Tmp)).Cross(v)+f1Tmp.Cross(vP));
      
      WorkVec.Put(13, -v);
      WorkVec.Put(16, v.Cross(f2Tmp));
      WorkVec.Put(19, -vP);
      WorkVec.Put(22, v.Cross(Omega2.Cross(f2Tmp))-f2Tmp.Cross(vP));
      
      WorkVec.Put(25, x1+f1Tmp-x2-f2Tmp);
      WorkVec.fPutCoef(28, 1.-dAlpha);
      
      WorkVec.Put(29, v1+Omega1.Cross(f1Tmp)-v2-Omega2.Cross(f2Tmp));
      WorkVec.fPutCoef(32, -dAlphaP);
      
   } else {
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
      ASSERT(d > DBL_EPSILON);
      if (d > DBL_EPSILON) {
	 d = sqrt(d);
	 WorkVec.fPutCoef(28, 1.-d);
	 WorkVec.fPutCoef(32, -v.Dot(vP)/d);
      }                  
   }
   
   return WorkVec;
}


void DistanceJointWithOffset::SetInitialValue(VectorHandler& X) const
{
   integer iFirstIndex = iGetFirstIndex();   
   
   doublereal dDistance = pGetDriveCaller()->dGet();
   if (fabs(dDistance) <= DBL_EPSILON) {
      X.Put(iFirstIndex+4, 1.);
      return;
   }
   
   Vec3 x1(pNode1->GetXCurr());
   Vec3 x2(pNode2->GetXCurr());
   Vec3 v1(pNode1->GetVCurr());
   Vec3 v2(pNode2->GetVCurr());
   Vec3 f1Tmp((pNode1->GetRCurr())*f1);
   Vec3 f2Tmp((pNode2->GetRCurr())*f2);
   Vec3 Omega1(pNode1->GetWCurr());
   Vec3 Omega2(pNode2->GetWCurr());   
   
   (Vec3&)v = x2+f2Tmp-x1-f1Tmp;
   doublereal d = v.Dot();
   ASSERT(d > DBL_EPSILON);
   if (d > DBL_EPSILON) {
      d = sqrt(d);
      (Vec3&)v = v/d;
   } else {
      (Vec3&)v = Vec3(1., 0., 0.);
   }
   
   X.Put(iFirstIndex+1, v);
   if (d > DBL_EPSILON) {      
      X.Put(iFirstIndex+5, (v2+Omega2.Cross(f2Tmp)-v1-Omega1.Cross(f1Tmp))/d);
   }   
}


void DistanceJointWithOffset::SetValue(VectorHandler& X, 
				       VectorHandler& /* XP */ ) const
{
   doublereal dDistance = pGetDriveCaller()->dGet();
   
   /* Setta a 1 dAlpha, che e' indipendente dalle altre variabili
    * in caso di distanza nulla */
   if (fabs(dDistance) <= DBL_EPSILON) {	
      X.Put(iGetFirstIndex()+4, 1.);
   } else {
      doublereal d = v.Dot();
      if (d <= DBL_EPSILON) {
	 (Vec3&)v = pNode2->GetXCurr()+pNode2->GetRCurr()*f2
	   -pNode1->GetXCurr()-pNode1->GetRCurr()*f1;
	 d = v.Dot();
	 if (d <= DBL_EPSILON) {
	    std::cerr << std::endl << "Joint " << uLabel << ", linked to nodes " 
	      << pNode1->GetLabel() << " and " << pNode2->GetLabel() 
	      << ": nodes are coincident." << std::endl
	      << "Initial joint assembly is recommended; aborting ... " 
	      << std::endl;
	    
	    THROW(ErrMemory());
	 }
	 (Vec3&)v /= sqrt(d);
      }
      
	/* Scrittura sul vettore soluzione */
      X.Put(iGetFirstIndex()+1, v);	 
   }   
}


void 
DistanceJointWithOffset::GetAdamsDummyPart(unsigned int part,
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
   
   Vec3 e(MatR2EulerAngles(MatR2vec(1, v1, 2, v2)));
   
   return out 
     << psAdamsElemCode[GetElemType()] << "_" << GetLabel() << "_" << part << std::endl
     << firstId << " "
     << x1 << " "
     << e /* MatR2EulerAngles(pNode1->GetRCurr()) */ << " "
     << x1 << " "
     << e << " "
     << l << " " << 0. << " " << 0. << " "
     << Zero3 << std::endl;
}

/* DistanceJointWithOffset - end */


/* ClampJoint - begin */

/* Costruttore definitivo (da mettere a punto) */
ClampJoint::ClampJoint(unsigned int uL, const DofOwner*pD, 
		       const StructNode* pN, 
		       const Vec3& X0, const Mat3x3& R0, 
		       flag fOut)
: Elem(uL, Elem::JOINT, fOut), 
Joint(uL, Joint::CLAMP, pD, fOut), 
pNode(pN), XClamp(X0), RClamp(R0), F(0.), M(0.)
{ 
   NO_OP; 
}


ClampJoint::~ClampJoint(void) 
{ 
   NO_OP; 
}


/* Contributo al file di restart */
std::ostream& ClampJoint::Restart(std::ostream& out) const
{
   return Joint::Restart(out) << ", clamp, "
     << pNode->GetLabel() << ", node, node;" << std::endl;
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
   WM.ResizeInit(12, 1, 0.);
   
   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   
   /* Matrice jacobiana del vincolo di incastro:
    * 
    *     posizione:  x = XClamp
    *     assetto:    R = RClamp
    * 
    *  -F          - Delta_F = 0
    *  -M          - Delta_M = 0
    *   x          + Delta_x = XClamp
    *  -g(R_Delta) + Delta_g = 0
    * 
    * con R_Delta = R*RClamp^T
    * si ottiene quindi:
    * 
    *  -Delta_F = F
    *  -Delta_M = M
    *   Delta_x = XClamp - x
    *   Delta_g = -g(R_Delta)
    * 
    * Nota: con la definizione dei nodi statici, questa parte non viene piu'
    * usata.
    * 
    * A questo si aggiunga la riga di definizione della quantita' di moto,
    * che deve essere identicamente nulla (si presume che per un nodo 
    * incastrato non sia definito elemento di massa)
    * 
    * 
    * Organizzazione del contributo allo jacobiano
    * 
    *     x, g    Q, Gamma   F, M     dot
    *  |   0      -dCoef*I    0  | |x, g    |   | Q, Gamma              |
    *  |   0         0       -I  | |Q, Gamma| = | F, M                  |
    *  | dCoef*I     0        0  | |F, M    |   | XClamp-x, -g(R_Delta) |
    * 
    * Il termine -dCoef*I viene messo in quanto si suppone che un nodo
    * vincolato non abbia un elemento corpo rigido definito.
    * In questo modo si impone che la quantita' di moto sia nulla
    * e la matrice jacobiana non e' singolare.
    * 
    * Viene usata la sottomatrice sparsa, in quanto sono richiesti solo
    * diciotto coefficienti non nulli */
   
   /* Attenzione: modifico dividendo le equazioni di vincolo per dCoef */
   
   for (integer iCnt = 1; iCnt <= 6; iCnt++) {
      WM.fPutItem(iCnt, iFirstReactionIndex+iCnt, 
		  iFirstPositionIndex+iCnt, 1.);    
      WM.fPutItem(6+iCnt, iFirstMomentumIndex+iCnt,
		  iFirstReactionIndex+iCnt, -1.);
   }
        
   /* Con l'aggiunta dei nodi statici non occorre piu' evitare la
    * singolarita' della matrice */
   
   return WorkMat;
}


/* assemblaggio matrici per autovalori */
void ClampJoint::AssMats(VariableSubMatrixHandler& WorkMatA,
			VariableSubMatrixHandler& WorkMatB,
			const VectorHandler& XCurr,
			const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering ClampJoint::AssMats(); will result in call to AssJac()"
		   << std::endl);
   
   WorkMatA.SetNullMatrix();
   AssJac(WorkMatB, 1., XCurr, XPrimeCurr);
}



SubVectorHandler& ClampJoint::AssRes(SubVectorHandler& WorkVec,
				     doublereal dCoef,
				     const VectorHandler& XCurr, 
				     const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering ClampJoint::AssRes()" << std::endl);

   integer iNumRows;
   integer iNumCols;
   this->WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);

   /* Per i commenti, vedi AssJac */
   
   /* Con l'aggiunta dei nodi statici non occorre piu' evitare 
    * la singolarita' della matrice */
   
   /* Indici delle incognite del nodo e delle reazioni */
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   for (integer iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iFirstMomentumIndex+iCnt);
      WorkVec.fPutRowIndex(6+iCnt, iFirstReactionIndex+iCnt);
   }   
   
   /* Aggiorna le reazioni vincolari */
   F = Vec3(XCurr, iFirstReactionIndex+1);
   M = Vec3(XCurr, iFirstReactionIndex+4);
   
   /* Calcola posizione e parametri di rotazione */
   Vec3 x(pNode->GetXCurr());
   Mat3x3 R(pNode->GetRCurr());
   /* Nota: si sfrutta il fatto che g(R^T) = -g(R) per avere -g(R_Delta) */
   Vec3 g(gparam(RClamp*R.Transpose()));
   

   /* Residuo della riga di equilibrio */
   WorkVec.Add(1, F);
   WorkVec.Add(4, M);
   
   /* Modifica: divido le equazioni di vincolo per dCoef */
   if (dCoef != 0.) {	
      /* Residuo dell'equazione di vincolo */
      WorkVec.Add(7, (XClamp-x)/dCoef);
      WorkVec.Add(10, g/dCoef);   
   }
   
   return WorkVec;
}


void ClampJoint::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) {
      Mat3x3 RT(pNode->GetRCurr().Transpose());

      Joint::Output(OH.Joints(), "Clamp", GetLabel(),
		    RT*F, RT*M, F, M) << std::endl;
   }
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
ClampJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
			  const VectorHandler& /* XCurr */ )
{
   DEBUGCOUT("Entering ClampJoint::InitialAssJac()" << std::endl);
   
   SparseSubMatrixHandler& WM = WorkMat.SetSparse();
   WM.ResizeInit(24, 1, 0.);
   
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
      WM.fPutItem(iCnt, iFirstReactionIndex+iCnt, 
		  iFirstPositionIndex+iCnt, 1.);
   
      WM.fPutItem(6+iCnt, iReactionPrimeIndex+iCnt, 
		  iFirstVelocityIndex+iCnt, 1.);
   
      WM.fPutItem(12+iCnt, iFirstPositionIndex+iCnt, 
		  iFirstReactionIndex+iCnt, -1.);
   
      WM.fPutItem(18+iCnt, iFirstVelocityIndex+iCnt, 
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
   this->InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
      
   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();  
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+6;

   for (int iCnt = 1; iCnt <= 12; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iFirstPositionIndex+iCnt);
      WorkVec.fPutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
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
   WorkVec.Put(16, gparam(RClamp*R.Transpose()));
   
   /* Velocita' */
   WorkVec.Put(19, -pNode->GetVCurr());
   
   /* Velocita' angolare */
   WorkVec.Put(22, -pNode->GetWCurr());
         
   return WorkVec;
}


/* Metodi per l'estrazione di dati "privati".
 * Si suppone che l'estrattore li sappia interpretare.
 * Come default non ci sono dati privati estraibili */
unsigned int ClampJoint::iGetNumPrivData(void) const
{
   return 6;
}


doublereal ClampJoint::dGetPrivData(unsigned int i) const
{
   if (i >= 1 && i <= 3) {
      return F.dGet(i);
   } else if (i >= 4 && i <= 6) {
      return M.dGet(i-3);
   } else {
      THROW(ErrGeneric());
   }
#ifndef USE_EXCEPTIONS
   return 0.;
#endif /* USE_EXCEPTIONS */
}

/* ClampJoint - end */
