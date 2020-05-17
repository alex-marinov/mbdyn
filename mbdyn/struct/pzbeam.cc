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

#include "dataman.h"
#include "pzbeam.h"
#include "shapefnc.h"

/* PiezoActuatorBeam - begin */

/* Funzioni di calcolo delle matrici */
void PiezoActuatorBeam::AssStiffnessMat(FullSubMatrixHandler& WMA,
					FullSubMatrixHandler& WMB,
					doublereal dCoef,
					const VectorHandler& XCurr,
					const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering PiezoActuatorBeam::AssStiffnessMat()" << std::endl);
   
   Beam::AssStiffnessMat(WMA, WMB, dCoef, XCurr, XPrimeCurr);
   
   Mat3xN tmp1(iNumElec);
   Mat3xN tmp2(iNumElec);

   tmp1.LeftMult(R[S_I]*dCoef, PiezoMat[STRAIN][S_I]);
   WMA.Sub(1, 19, tmp1);
   WMA.Add(7, 19, tmp1);
   tmp2.LeftMult(Mat3x3(MatCross, p[S_I] - pNode[NODE1]->GetXCurr()), tmp1);
   WMA.Add(4, 19, tmp2); 
   tmp2.LeftMult(Mat3x3(MatCross, p[S_I] - pNode[NODE2]->GetXCurr()), tmp1);
   WMA.Sub(10, 19, tmp2);
   
   tmp1.LeftMult(R[S_I]*dCoef, PiezoMat[CURVAT][S_I]);
   WMA.Sub(4, 19, tmp1);
   WMA.Add(10, 19, tmp1);

   tmp1.LeftMult(R[SII]*dCoef, PiezoMat[STRAIN][SII]);
   WMA.Sub(7, 19, tmp1);
   WMA.Add(13, 19, tmp1);
   tmp2.LeftMult(Mat3x3(MatCross, p[SII] - pNode[NODE2]->GetXCurr()), tmp1);
   WMA.Add(10, 19, tmp2);
   tmp2.LeftMult(Mat3x3(MatCross, p[SII] - pNode[NODE3]->GetXCurr()), tmp1);
   WMA.Sub(16, 19, tmp2);

   tmp1.LeftMult(R[SII]*dCoef, PiezoMat[CURVAT][SII]);
   WMA.Sub(10, 19, tmp1);
   WMA.Add(16, 19, tmp1);
}


void PiezoActuatorBeam::AssStiffnessVec(SubVectorHandler& WorkVec,
					doublereal dCoef,
					const VectorHandler& XCurr,
					const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering PiezoActuatorBeam::AssStiffnessVec()" << std::endl);
   
   /* Riceve il vettore gia' dimensionato e con gli indici a posto 
    * per scrivere il residuo delle equazioni di equilibrio dei tre nodi */
   
   /* Per la trattazione teorica, il riferimento e' il file ul-travi.tex 
    * (ora e' superato) */
   
   if (bFirstRes) {
      /* fFirstRes = flag(0);  AfterPredict ha gia' calcolato tutto */
   } else {    
      for (integer iCnt = 1; iCnt <= iNumElec; iCnt++) {
	 V.Put(iCnt, pvElecDofs[iCnt-1]->dGetX());
      }
   }
   
   Beam::AssStiffnessVec(WorkVec, dCoef, XCurr, XPrimeCurr);
}


void PiezoActuatorBeam::AddInternalForces(Vec6& AzLoc, unsigned int iSez)
{
   AzLoc += Vec6(PiezoMat[STRAIN][iSez]*V, PiezoMat[CURVAT][iSez]*V);
}
      

/* Costruttore normale */
PiezoActuatorBeam::PiezoActuatorBeam(
		unsigned int uL,
		const StructNode* pN1,
		const StructNode* pN2,
		const StructNode* pN3,
		const Vec3& F1,
		const Vec3& F2,
		const Vec3& F3,
		const Mat3x3& R1,
		const Mat3x3& R2,
		const Mat3x3& R3,
		const Mat3x3& r_I, const Mat3x3& rII,
		const ConstitutiveLaw6D* pD_I,
		const ConstitutiveLaw6D* pDII,
		int iEl,
		const ScalarDifferentialNode **pEDof,
		const Mat3xN& T_Ie, const Mat3xN& T_Ik,
		const Mat3xN& TIIe, const Mat3xN& TIIk,
		OrientationDescription ood,
		flag fOut
) : Elem(uL, fOut),
Beam(uL, pN1, pN2, pN3, F1, F2, F3, R1, R2, R3, r_I, rII, pD_I, pDII, ood, fOut),
iNumElec(iEl), pvElecDofs(pEDof), V(iEl)
{
#ifdef DEBUG
   ASSERT(iNumElec > 0);
   ASSERT(pvElecDofs != NULL);
   for (int i = iNumElec; i-- > 0; ) {
      ASSERT(pvElecDofs[i] != NULL);
   }   
#endif // DEBUG

   PiezoMat[STRAIN][S_I].Copy(T_Ie);
   PiezoMat[STRAIN][SII].Copy(TIIe);
   PiezoMat[CURVAT][S_I].Copy(T_Ik);
   PiezoMat[CURVAT][SII].Copy(TIIk);  
}

   
/* Distruttore banale */
PiezoActuatorBeam::~PiezoActuatorBeam(void)
{
   SAFEDELETEARR(pvElecDofs);
}

    
/* Contributo al file di restart */
std::ostream&
PiezoActuatorBeam::Restart(std::ostream& out) const
{
   Restart_(out);
   return out << "/* piezoelectric actuator NOT IMPLEMENTED YET */ ;" << std::endl;
}

    
/* Dimensioni del workspace; sono 36 righe perche' se genera anche le
 *     * forze d'inerzia consistenti deve avere accesso alle righe di definizione
 *     * della quantita' di moto */
void PiezoActuatorBeam::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
   if (bConsistentInertia) {
      *piNumRows = 36;
   } else {
      *piNumRows = 18;
   }
   *piNumCols = 18+iNumElec;
}


/* Settings iniziali, prima della prima soluzione */
void PiezoActuatorBeam::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
   /* se proprio non serve, si puo' eliminare */
   Beam::SetValue(pDM, X, XP, ph);
}


/* Prepara i parametri di riferimento dopo la predizione */
void PiezoActuatorBeam::AfterPredict(VectorHandler& X,
				     VectorHandler& XP)
{
   // Calcola le deformazioni, aggiorna il legame costitutivo e crea la FDE
   
   for (integer iCnt = 1; iCnt <= iNumElec; iCnt++) {
	 V.Put(iCnt, pvElecDofs[iCnt-1]->dGetX());
   }
   
   Beam::AfterPredict(X, XP);
}


/* assemblaggio jacobiano */
VariableSubMatrixHandler&
PiezoActuatorBeam::AssJac(VariableSubMatrixHandler& WorkMat,
			  doublereal dCoef,
			  const VectorHandler& XCurr,
			  const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering PiezoActuatorBeam::AssJac(); will result in call to AssStiffnessMat()" << std::endl);
   
   integer iNode1FirstMomIndex = pNode[NODE1]->iGetFirstMomentumIndex();
   integer iNode1FirstPosIndex = pNode[NODE1]->iGetFirstPositionIndex();
   integer iNode2FirstMomIndex = pNode[NODE2]->iGetFirstMomentumIndex();
   integer iNode2FirstPosIndex = pNode[NODE2]->iGetFirstPositionIndex();
   integer iNode3FirstMomIndex = pNode[NODE3]->iGetFirstMomentumIndex();
   integer iNode3FirstPosIndex = pNode[NODE3]->iGetFirstPositionIndex();
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();

   /* Dimensiona la matrice, la azzera e pone gli indici corretti */
   if(bConsistentInertia) {
      WM.ResizeReset(36, 18+iNumElec);
   } else {
      WM.ResizeReset(18, 18+iNumElec);
   }
   
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
      WM.PutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
      WM.PutRowIndex(12+iCnt, iNode3FirstMomIndex+iCnt);
      WM.PutColIndex(12+iCnt, iNode3FirstPosIndex+iCnt);
   }
 
   for (int iCnt = 1; iCnt <= iNumElec; iCnt++) {
      WM.PutColIndex(18+iCnt, pvElecDofs[iCnt-1]->iGetFirstColIndex()+1);
   }   
   
   AssStiffnessMat(WM, WM, dCoef, XCurr, XPrimeCurr);
   
   if (bConsistentInertia) {
      for (int iCnt = 1; iCnt <= 6; iCnt++) {
	 WM.PutRowIndex(18+iCnt, iNode1FirstPosIndex+iCnt);
	 WM.PutRowIndex(24+iCnt, iNode2FirstPosIndex+iCnt);
	 WM.PutRowIndex(30+iCnt, iNode3FirstPosIndex+iCnt);
      }
      Beam::AssInertiaMat(WM, WM, dCoef, XCurr, XPrimeCurr);
   }
  
   return WorkMat;        
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
PiezoActuatorBeam::InitialAssJac(VariableSubMatrixHandler& WorkMat,
				 const VectorHandler& XCurr)
{
   return Beam::InitialAssJac(WorkMat, XCurr);
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
PiezoActuatorBeam::InitialAssRes(SubVectorHandler& WorkVec,
				 const VectorHandler& XCurr)
{
   return Beam::InitialAssRes(WorkVec, XCurr);
}

/* PiezoActuatorBeam - end */


/* PiezoActuatorVEBeam - begin */

/* Funzioni di calcolo delle matrici */
void PiezoActuatorVEBeam::AssStiffnessMat(FullSubMatrixHandler& WMA,
					  FullSubMatrixHandler& WMB,
					  doublereal dCoef,
					  const VectorHandler& XCurr,
					  const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering PiezoActuatorVEBeam::AssStiffnessMat()" << std::endl);
   
   ViscoElasticBeam::AssStiffnessMat(WMA, WMB, dCoef, XCurr, XPrimeCurr);
   
   Mat3xN tmp1(iNumElec);
   Mat3xN tmp2(iNumElec);

   tmp1.LeftMult(R[S_I]*dCoef, PiezoMat[STRAIN][S_I]);
   WMA.Sub(1, 19, tmp1);
   WMA.Add(7, 19, tmp1);
   tmp2.LeftMult(Mat3x3(MatCross, p[S_I] - pNode[NODE1]->GetXCurr()), tmp1);
   WMA.Add(4, 19, tmp2); 
   tmp2.LeftMult(Mat3x3(MatCross, p[S_I] - pNode[NODE2]->GetXCurr()), tmp1);
   WMA.Sub(10, 19, tmp2);
   
   tmp1.LeftMult(R[S_I]*dCoef, PiezoMat[CURVAT][S_I]);
   WMA.Sub(4, 19, tmp1);
   WMA.Add(10, 19, tmp1);   

   tmp1.LeftMult(R[SII]*dCoef, PiezoMat[STRAIN][SII]);
   WMA.Sub(7, 19, tmp1);
   WMA.Add(13, 19, tmp1);
   tmp2.LeftMult(Mat3x3(MatCross, p[SII] - pNode[NODE2]->GetXCurr()), tmp1);
   WMA.Add(10, 19, tmp2);
   tmp2.LeftMult(Mat3x3(MatCross, p[SII] - pNode[NODE3]->GetXCurr()), tmp1);
   WMA.Sub(16, 19, tmp2);

   tmp1.LeftMult(R[SII]*dCoef, PiezoMat[CURVAT][SII]);
   WMA.Sub(10, 19, tmp1);
   WMA.Add(16, 19, tmp1);
}


void PiezoActuatorVEBeam::AssStiffnessVec(SubVectorHandler& WorkVec,
					doublereal dCoef,
					const VectorHandler& XCurr,
					const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering PiezoActuatorVEBeam::AssStiffnessVec()" << std::endl);
   
   /* Riceve il vettore gia' dimensionato e con gli indici a posto 
    * per scrivere il residuo delle equazioni di equilibrio dei tre nodi */
   
   /* Per la trattazione teorica, il riferimento e' il file ul-travi.tex 
    * (ora e' superato) */
   
   if (bFirstRes) {
      /* fFirstRes = flag(0);  AfterPredict ha gia' calcolato tutto */
   } else {    
      for (integer iCnt = 1; iCnt <= iNumElec; iCnt++) {
         V.Put(iCnt, pvElecDofs[iCnt-1]->dGetX());
      }
   }
   
   ViscoElasticBeam::AssStiffnessVec(WorkVec, dCoef, XCurr, XPrimeCurr);
}


void PiezoActuatorVEBeam::AddInternalForces(Vec6& AzLoc, unsigned int iSez)
{
   AzLoc += Vec6(PiezoMat[STRAIN][iSez]*V, PiezoMat[CURVAT][iSez]*V);
}
      

/* Costruttore normale */
PiezoActuatorVEBeam::PiezoActuatorVEBeam(
		unsigned int uL,
		const StructNode* pN1,
		const StructNode* pN2,
		const StructNode* pN3,
		const Vec3& F1,
		const Vec3& F2,
		const Vec3& F3,
		const Mat3x3& R1,
		const Mat3x3& R2,
		const Mat3x3& R3,
		const Mat3x3& r_I,
		const Mat3x3& rII,
		const ConstitutiveLaw6D* pD_I,
		const ConstitutiveLaw6D* pDII,
		int iEl,
		const ScalarDifferentialNode **pEDof,
		const Mat3xN& T_Ie,
		const Mat3xN& T_Ik,
		const Mat3xN& TIIe,
		const Mat3xN& TIIk,
		OrientationDescription ood,
		flag fOut
) : Elem(uL, fOut),
ViscoElasticBeam(uL, pN1, pN2, pN3, F1, F2, F3, R1, R2, R3,
		r_I, rII, pD_I, pDII, ood, fOut),
iNumElec(iEl), pvElecDofs(pEDof), V(iEl)
{
#ifdef DEBUG
   ASSERT(iNumElec > 0);
   ASSERT(pvElecDofs != NULL);
   for (int i = iNumElec; i-- > 0; ) {
      ASSERT(pvElecDofs[i] != NULL);
   }   
#endif /* DEBUG */

   PiezoMat[STRAIN][S_I].Copy(T_Ie);
   PiezoMat[STRAIN][SII].Copy(TIIe);
   PiezoMat[CURVAT][S_I].Copy(T_Ik);
   PiezoMat[CURVAT][SII].Copy(TIIk);  
}

   
/* Distruttore banale */
PiezoActuatorVEBeam::~PiezoActuatorVEBeam(void)
{
   SAFEDELETEARR(pvElecDofs);
}

    
/* Contributo al file di restart */
std::ostream&
PiezoActuatorVEBeam::Restart(std::ostream& out) const
{
   Restart_(out);
   return out << "/* piezoelectric actuator NOT IMPLEMENTED YET */ ;" << std::endl;
}


/* Dimensioni del workspace; sono 36 righe perche' se genera anche le
 * forze d'inerzia consistenti deve avere accesso alle righe di definizione
 * della quantita' di moto */
void
PiezoActuatorVEBeam::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
   if (bConsistentInertia) {
      *piNumRows = 36;
   } else {
      *piNumRows = 18;
   }
   *piNumCols = 18+iNumElec;
}


/* Settings iniziali, prima della prima soluzione */
void PiezoActuatorVEBeam::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
   /* se proprio non serve, si puo' eliminare */
   ViscoElasticBeam::SetValue(pDM, X, XP, ph);
}


/* Prepara i parametri di riferimento dopo la predizione */
void PiezoActuatorVEBeam::AfterPredict(VectorHandler& X,
				     VectorHandler& XP)
{
   // Calcola le deformazioni, aggiorna il legame costitutivo e crea la FDE
   
   for (integer iCnt = 1; iCnt <= iNumElec; iCnt++) {
	 V.Put(iCnt, pvElecDofs[iCnt-1]->dGetX());
   }
   
   ViscoElasticBeam::AfterPredict(X, XP);
}


/* assemblaggio jacobiano */
VariableSubMatrixHandler&
PiezoActuatorVEBeam::AssJac(VariableSubMatrixHandler& WorkMat,
			  doublereal dCoef,
			  const VectorHandler& XCurr,
			  const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering PiezoActuatorVEBeam::AssJac(); will result in call to AssStiffnessMat()" << std::endl);
   
   integer iNode1FirstMomIndex = pNode[NODE1]->iGetFirstMomentumIndex();
   integer iNode1FirstPosIndex = pNode[NODE1]->iGetFirstPositionIndex();
   integer iNode2FirstMomIndex = pNode[NODE2]->iGetFirstMomentumIndex();
   integer iNode2FirstPosIndex = pNode[NODE2]->iGetFirstPositionIndex();
   integer iNode3FirstMomIndex = pNode[NODE3]->iGetFirstMomentumIndex();
   integer iNode3FirstPosIndex = pNode[NODE3]->iGetFirstPositionIndex();
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();

   /* Dimensiona la matrice, la azzera e pone gli indici corretti */
   if(bConsistentInertia) {
      WM.ResizeReset(36, 18+iNumElec);
   } else {
      WM.ResizeReset(18, 18+iNumElec);
   }
   
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
      WM.PutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
      WM.PutRowIndex(12+iCnt, iNode3FirstMomIndex+iCnt);
      WM.PutColIndex(12+iCnt, iNode3FirstPosIndex+iCnt);
   }
 
   for (int iCnt = 1; iCnt <= iNumElec; iCnt++) {
      WM.PutColIndex(18+iCnt, pvElecDofs[iCnt-1]->iGetFirstColIndex()+1);
   }   
   
   AssStiffnessMat(WM, WM, dCoef, XCurr, XPrimeCurr);
   
   if (bConsistentInertia) {
      for (int iCnt = 1; iCnt <= 6; iCnt++) {
	 WM.PutRowIndex(18+iCnt, iNode1FirstPosIndex+iCnt);
	 WM.PutRowIndex(24+iCnt, iNode2FirstPosIndex+iCnt);
	 WM.PutRowIndex(30+iCnt, iNode3FirstPosIndex+iCnt);
      }
      ViscoElasticBeam::AssInertiaMat(WM, WM, dCoef, XCurr, XPrimeCurr);
   }
  
   return WorkMat;        
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
PiezoActuatorVEBeam::InitialAssJac(VariableSubMatrixHandler& WorkMat,
				 const VectorHandler& XCurr)
{
   return ViscoElasticBeam::InitialAssJac(WorkMat, XCurr);
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
PiezoActuatorVEBeam::InitialAssRes(SubVectorHandler& WorkVec,
				 const VectorHandler& XCurr)
{
   return ViscoElasticBeam::InitialAssRes(WorkVec, XCurr);
}

/* PiezoActuatorVEBeam - end */

/* Full Piezo beam */

/* Funzioni di calcolo delle matrici */
void PiezoBeam::AssStiffnessMat(FullSubMatrixHandler& WMA,
					  FullSubMatrixHandler& WMB,
					  doublereal dCoef,
					  const VectorHandler& XCurr,
					  const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering PiezoActuatorVEBeam::AssStiffnessMat()" << std::endl);
   
   ViscoElasticBeam::AssStiffnessMat(WMA, WMB, dCoef, XCurr, XPrimeCurr);
   
   Mat3xN tmp1(iNumElec);
   Mat3xN tmp2(iNumElec);

   tmp1.LeftMult(R[S_I]*dCoef, PiezoMat[STRAIN][S_I]);
   WMA.Sub(1, 19, tmp1);
   WMA.Add(7, 19, tmp1);
   tmp2.LeftMult(Mat3x3(MatCross, p[S_I] - pNode[NODE1]->GetXCurr()), tmp1);
   WMA.Add(4, 19, tmp2); 
   tmp2.LeftMult(Mat3x3(MatCross, p[S_I] - pNode[NODE2]->GetXCurr()), tmp1);
   WMA.Sub(10, 19, tmp2);
   
   tmp1.LeftMult(R[S_I]*dCoef, PiezoMat[CURVAT][S_I]);
   WMA.Sub(4, 19, tmp1);
   WMA.Add(10, 19, tmp1);   

   tmp1.LeftMult(R[SII]*dCoef, PiezoMat[STRAIN][SII]);
   WMA.Sub(7, 19, tmp1);
   WMA.Add(13, 19, tmp1);
   tmp2.LeftMult(Mat3x3(MatCross, p[SII] - pNode[NODE2]->GetXCurr()), tmp1);
   WMA.Add(10, 19, tmp2);
   tmp2.LeftMult(Mat3x3(MatCross, p[SII] - pNode[NODE3]->GetXCurr()), tmp1);
   WMA.Sub(16, 19, tmp2);

   tmp1.LeftMult(R[SII]*dCoef, PiezoMat[CURVAT][SII]);
   WMA.Sub(10, 19, tmp1);
   WMA.Add(16, 19, tmp1);
}

/* Costruttore normale */
PiezoBeam::PiezoBeam(
		unsigned int uL,
		const StructNode* pN1,
		const StructNode* pN2,
		const StructNode* pN3,
		const Vec3& F1,
		const Vec3& F2,
		const Vec3& F3,
		const Mat3x3& R1,
		const Mat3x3& R2,
		const Mat3x3& R3,
		const Mat3x3& r_I,
		const Mat3x3& rII,
		const ConstitutiveLaw6D* pD_I,
		const ConstitutiveLaw6D* pDII,
		int iEl,
		const ScalarDifferentialNode **pEDof,
		const Mat3xN& T_Ie,
		const Mat3xN& T_Ik,
		const Mat3xN& TIIe,
		const Mat3xN& TIIk,
		const MatNx3& Q_Ie,
		const MatNx3& Q_Ik,
		const MatNx3& QIIe,
		const MatNx3& QIIk,
      const MatNxN& QV_I,
      const MatNxN& QVII,
		OrientationDescription ood,
		flag fOut
) : Elem(uL, fOut),
   PiezoActuatorVEBeam(uL, pN1, pN2, pN3, F1, F2, F3, R1, R2, R3,
		r_I, rII, pD_I, pDII, iEl, pEDof, T_Ie, T_Ik, TIIe, TIIk, ood, fOut),
   Vdot(iEl),
   Qdot(iEl)
{
   PiezoMatQ[STRAIN][S_I].Copy(Q_Ie);
   PiezoMatQ[STRAIN][SII].Copy(QIIe);
   PiezoMatQ[CURVAT][S_I].Copy(Q_Ik);
   PiezoMatQ[CURVAT][SII].Copy(QIIk);
   PiezoMatQV[S_I].Copy(QV_I);
   PiezoMatQV[SII].Copy(QVII);
}

   
/* Distruttore banale */
PiezoBeam::~PiezoBeam(void)
{
   SAFEDELETEARR(pvElecDofs);
}

    
/* Contributo al file di restart */
std::ostream&
PiezoBeam::Restart(std::ostream& out) const
{
   Restart_(out);
   return out << "/* piezoelectric beam restart NOT IMPLEMENTED YET */ ;" << std::endl;
}


/* Dimensioni del workspace; sono 36 righe perche' se genera anche le
 * forze d'inerzia consistenti deve avere accesso alle righe di definizione
 * della quantita' di moto */
void
PiezoBeam::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
   if (bConsistentInertia) {
      *piNumRows = 36;
   } else {
      *piNumRows = 18;
   }
   *piNumRows += iNumElec;
   *piNumCols = 18+iNumElec;
}

SubVectorHandler&
PiezoBeam::AssRes(SubVectorHandler& WorkVec,
				    doublereal dCoef,
				    const VectorHandler& XCurr,
				    const VectorHandler& XPrimeCurr)
{
   DEBUGCOUTFNAME("PiezoBeam::AssRes => AssStiffnessVec");

	integer iNode1FirstMomIndex = pNode[NODE1]->iGetFirstMomentumIndex();
	integer iNode2FirstMomIndex = pNode[NODE2]->iGetFirstMomentumIndex();
	integer iNode3FirstMomIndex = pNode[NODE3]->iGetFirstMomentumIndex();

   int iNumInerzia = 18;
   if(bConsistentInertia) {
      iNumInerzia = 36;
   }

	/* Dimensiona il vettore, lo azzera e pone gli indici corretti */
	if (bConsistentInertia) {
		WorkVec.ResizeReset(iNumInerzia+iNumElec);
	} else {
		WorkVec.ResizeReset(iNumInerzia+iNumElec);
	}

	for (unsigned int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
		WorkVec.PutRowIndex(12 + iCnt, iNode3FirstMomIndex + iCnt);
	}

	if (bConsistentInertia) {
		integer iNode1FirstPosIndex = pNode[NODE1]->iGetFirstPositionIndex();
		integer iNode2FirstPosIndex = pNode[NODE2]->iGetFirstPositionIndex();
		integer iNode3FirstPosIndex = pNode[NODE3]->iGetFirstPositionIndex();

		for (unsigned int iCnt = 1; iCnt <= 6; iCnt++) {
			WorkVec.PutRowIndex(18 + iCnt, iNode1FirstPosIndex + iCnt);
			WorkVec.PutRowIndex(24 + iCnt, iNode2FirstPosIndex + iCnt);
			WorkVec.PutRowIndex(30 + iCnt, iNode3FirstPosIndex + iCnt);
		}
   }

   for (int iCnt = 1; iCnt <= iNumElec; iCnt++) {
      WorkVec.PutRowIndex(iNumInerzia+iCnt, pvElecDofs[iCnt-1]->iGetFirstRowIndex()+1);
   }

   PiezoActuatorVEBeam::AssStiffnessVec(WorkVec, dCoef, XCurr, XPrimeCurr);
   AssPiezoRes(WorkVec, dCoef, XCurr, XPrimeCurr);

	return WorkVec;

}

void PiezoBeam::AssPiezoJac(FullSubMatrixHandler& WM,
				doublereal dCoef,
				const VectorHandler& XCurr,
				const VectorHandler& XPrimeCurr) {
   int iNumInerzia = 18;
   if(bConsistentInertia) {
      iNumInerzia = 36;
   }

    Vec3 fTmp[NUMNODES];
    doublereal omegafTmp[NUMNODES];
    Mat3x3 omega_o_fTmp[NUMNODES];
    for (unsigned int i = 0; i < NUMNODES; i++) {
        fTmp[i] = pNode[i]->GetRCurr()*f[i];
        omegafTmp[i] = pNode[i]->GetWRef() * fTmp[i];
        omega_o_fTmp[i] = pNode[i]->GetWRef().Tens(fTmp[i]);
    }

    Mat6x6 DefPrimeTmp[NUMSEZ][NUMNODES];
    /* per ogni punto di valutazione: */
    for (unsigned int iSez = 0; iSez < NUMSEZ; iSez++) {
        for (unsigned int node = 0; node < NUMNODES; node++) {
            /* Delta - deformazioni */
            DefPrimeTmp[iSez][node]
                = Mat6x6(mb_deye<Mat3x3>(dN3P[iSez][node]*dsdxi[iSez]),
                    Zero3x3,
                    Mat3x3(MatCross, L[iSez]*dN3[iSez][node] - fTmp[node]*(dN3P[iSez][node]*dsdxi[iSez])),
                    mb_deye<Mat3x3>(dN3P[iSez][node]*dsdxi[iSez]));
            DefPrimeTmp[iSez][node] +=
                Mat6x6(Mat3x3(MatCross, Omega[iSez]*(-dN3P[iSez][node]*dsdxi[iSez]*dCoef)),
                    Zero3x3,
//                     (Mat3x3(MatCross, LPrime[iSez]) - 
//                     Mat3x3(MatCrossCross, Omega[iSez], L[iSez]))*(dN3[iSez][node]*dCoef) //TEMO SIA SBAGLIATO
                    Mat3x3(MatCross, LPrime[iSez] - Omega[iSez].Cross(L[iSez]))*(dN3[iSez][node]*dCoef) //VERSIONE MIA
                        + Mat3x3(MatCrossCross, Omega[iSez], fTmp[node]*(dN3P[iSez][node]*dsdxi[iSez]*dCoef))
                        + Mat3x3(MatCross, fTmp[node].Cross(pNode[node]->GetWCurr()*(dN3P[iSez][node]*dsdxi[iSez]*dCoef)))
                        - mb_deye<Mat3x3>(omegafTmp[node]*dN3P[iSez][node]*dsdxi[iSez]*dCoef)
                        + omega_o_fTmp[node] * (dN3P[iSez][node]*dsdxi[iSez]*dCoef)
                        + mb_deye<Mat3x3>(OmegaRef[iSez]*L[iSez]*dN3[iSez][node]*dCoef)
                        - OmegaRef[iSez].Tens(L[iSez])*(dN3[iSez][node]*dCoef),
                    Mat3x3(MatCross, Omega[iSez]*(-dN3P[iSez][node]*dsdxi[iSez]*dCoef))
                );
            DefPrimeTmp[iSez][node].GetMat11() = RRef[iSez].MulTM(DefPrimeTmp[iSez][node].GetMat11());
            DefPrimeTmp[iSez][node].GetMat12() = RRef[iSez].MulTM(DefPrimeTmp[iSez][node].GetMat12());
            DefPrimeTmp[iSez][node].GetMat21() = RRef[iSez].MulTM(DefPrimeTmp[iSez][node].GetMat21());
            DefPrimeTmp[iSez][node].GetMat22() = RRef[iSez].MulTM(DefPrimeTmp[iSez][node].GetMat22());
        }
	
    } /* end ciclo sui punti di valutazione */
   for (unsigned int iSez = 0; iSez < NUMSEZ; iSez++) {
      for (int r_el = 1; r_el <= iNumElec; r_el++) {
         for (int c_el = 1; c_el <= iNumElec; c_el++) {
            WM.DecCoef(iNumInerzia+r_el, 18+c_el, PiezoMatQV[iSez](r_el, c_el) / dsdxi[iSez]); //FIXME
         }
         for (int def_comp = 1; def_comp <= 3; def_comp++) {
            for (int node = 0; node < NUMNODES; node++) {
               for (int dof = 1; dof <= 3; dof++) {
                  WM.DecCoef(iNumInerzia+r_el, 6*node + dof, PiezoMatQ[0][iSez](r_el, def_comp) * DefPrimeTmp[iSez][node].GetMat11()(def_comp, dof) / dsdxi[iSez]);
                  WM.DecCoef(iNumInerzia+r_el, 6*node + 3 + dof, PiezoMatQ[0][iSez](r_el, def_comp) * DefPrimeTmp[iSez][node].GetMat12()(def_comp, dof) / dsdxi[iSez]);
                  WM.DecCoef(iNumInerzia+r_el, 6*node + dof, PiezoMatQ[1][iSez](r_el, def_comp) * DefPrimeTmp[iSez][node].GetMat21()(def_comp, dof) / dsdxi[iSez]);
                  WM.DecCoef(iNumInerzia+r_el, 6*node + 3 + dof, PiezoMatQ[1][iSez](r_el, def_comp) * DefPrimeTmp[iSez][node].GetMat22()(def_comp, dof) / dsdxi[iSez]);
               }
            }
         }
      }
   }
}

void PiezoBeam::AssPiezoRes(SubVectorHandler& WorkVec,
				doublereal dCoef,
				const VectorHandler& XCurr,
				const VectorHandler& XPrimeCurr) {

   int iNumInerzia = 18;
   if(bConsistentInertia) {
      iNumInerzia = 36;
   }

   for (integer iCnt = 1; iCnt <= iNumElec; iCnt++) {
	   Vdot.Put(iCnt, pvElecDofs[iCnt-1]->dGetXPrime());
   }

   Qdot.Reset();
   for (int iSez = 0; iSez < NUMSEZ; iSez++) {
      for (int r_el = 1; r_el <= iNumElec; r_el++) {
         doublereal Qt = 0.;
         for (int c_el = 1; c_el <= iNumElec; c_el++) {
            Qt += PiezoMatQV[iSez](r_el, c_el) * Vdot(c_el);
         }
         for (int def_component = 1; def_component <= 3; def_component++) {
            Qt += PiezoMatQ[0][iSez](r_el, def_component) * DefPrimeLoc[iSez].GetVec1()(def_component);
            Qt += PiezoMatQ[1][iSez](r_el, def_component) * DefPrimeLoc[iSez].GetVec2()(def_component);
         }
         Qt /= dsdxi[iSez]; //Gauss weight == 1
         Qdot(r_el) += Qt;
      }
   }
   
   for (int r = 1; r <= iNumElec; r++) {
      // positive current: flowing out
      // but we need to compute desidual = - sum(flowing out current)
      WorkVec.IncCoef(iNumInerzia+r, Qdot(r));
   }
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
PiezoBeam::AssJac(VariableSubMatrixHandler& WorkMat,
			  doublereal dCoef,
			  const VectorHandler& XCurr,
			  const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering PiezoActuatorVEBeam::AssJac(); will result in call to AssStiffnessMat()" << std::endl);
   
   integer iNode1FirstMomIndex = pNode[NODE1]->iGetFirstMomentumIndex();
   integer iNode1FirstPosIndex = pNode[NODE1]->iGetFirstPositionIndex();
   integer iNode2FirstMomIndex = pNode[NODE2]->iGetFirstMomentumIndex();
   integer iNode2FirstPosIndex = pNode[NODE2]->iGetFirstPositionIndex();
   integer iNode3FirstMomIndex = pNode[NODE3]->iGetFirstMomentumIndex();
   integer iNode3FirstPosIndex = pNode[NODE3]->iGetFirstPositionIndex();
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();

   /* Dimensiona la matrice, la azzera e pone gli indici corretti */
   int iNumInerzia = 18;
   if(bConsistentInertia) {
      iNumInerzia = 36;
   }
   WM.ResizeReset(iNumInerzia+iNumElec, 18+iNumElec);
   
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
      WM.PutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
      WM.PutRowIndex(12+iCnt, iNode3FirstMomIndex+iCnt);
      WM.PutColIndex(12+iCnt, iNode3FirstPosIndex+iCnt);
   }
 
   for (int iCnt = 1; iCnt <= iNumElec; iCnt++) {
      WM.PutRowIndex(iNumInerzia+iCnt, pvElecDofs[iCnt-1]->iGetFirstRowIndex()+1);
      WM.PutColIndex(18+iCnt, pvElecDofs[iCnt-1]->iGetFirstColIndex()+1);
   }
   
   PiezoActuatorVEBeam::AssStiffnessMat(WM, WM, dCoef, XCurr, XPrimeCurr);
   
   if (bConsistentInertia) {
      for (int iCnt = 1; iCnt <= 6; iCnt++) {
	 WM.PutRowIndex(18+iCnt, iNode1FirstPosIndex+iCnt);
	 WM.PutRowIndex(24+iCnt, iNode2FirstPosIndex+iCnt);
	 WM.PutRowIndex(30+iCnt, iNode3FirstPosIndex+iCnt);
      }
      ViscoElasticBeam::AssInertiaMat(WM, WM, dCoef, XCurr, XPrimeCurr);
   }

   AssPiezoJac(WM, dCoef, XCurr, XPrimeCurr);
  
   return WorkMat;        
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
PiezoBeam::InitialAssJac(VariableSubMatrixHandler& WorkMat,
				 const VectorHandler& XCurr)
{
   return ViscoElasticBeam::InitialAssJac(WorkMat, XCurr);
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
PiezoBeam::InitialAssRes(SubVectorHandler& WorkVec,
				 const VectorHandler& XCurr)
{
   return ViscoElasticBeam::InitialAssRes(WorkVec, XCurr);
}

/* Full Piezo beam  - end*/
