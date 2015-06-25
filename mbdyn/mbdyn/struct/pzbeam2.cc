/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

#include "mynewmem.h"
#include "pzbeam2.h"

/* PiezoActuatorBeam2 - begin */

/* Funzioni di calcolo delle matrici */
void
PiezoActuatorBeam2::AssStiffnessMat(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering PiezoActuatorBeam2::AssStiffnessMat()" << std::endl);
	
   Beam2::AssStiffnessMat(WMA, WMB, dCoef, XCurr, XPrimeCurr);
   
   Mat3xN tmp1(iNumElec);
   Mat3xN tmp2(iNumElec);

   tmp1.LeftMult(R*dCoef, PiezoMat[STRAIN]);
   WMA.Sub(1, 13, tmp1);
   WMA.Add(7, 13, tmp1);
   tmp2.LeftMult(Mat3x3(MatCross, p - pNode[NODE1]->GetXCurr()), tmp1);
   WMA.Add(4, 13, tmp2); 
   tmp2.LeftMult(Mat3x3(MatCross, p - pNode[NODE2]->GetXCurr()), tmp1);
   WMA.Sub(10, 13, tmp2);
   
   tmp1.LeftMult(R*dCoef, PiezoMat[CURVAT]);
   WMA.Sub(4, 13, tmp1);
   WMA.Add(10, 13, tmp1);   
}


void
PiezoActuatorBeam2::AssStiffnessVec(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering PiezoActuatorBeam2::AssStiffnessVec()" << std::endl);
	
	/*
	 * Riceve il vettore gia' dimensionato e con gli indici a posto 
	 * per scrivere il residuo delle equazioni di equilibrio dei tre nodi
	 */
	
	/* 
	 * Per la trattazione teorica, il riferimento e' il file ul-travi.tex 
	 * (ora e' superato)
	 */
	
	if (bFirstRes) {
		/*
		 * AfterPredict ha gia' calcolato tutto;
		 * bFirstRes viene resettato poi da Beam2::AssStiffnessVec
		 */

	} else {    
		for (integer iCnt = 1; iCnt <= iNumElec; iCnt++) {
			V.Put(iCnt, pvElecDofs[iCnt-1]->dGetX());
		}
	}
	
	Beam2::AssStiffnessVec(WorkVec, dCoef, XCurr, XPrimeCurr);
}


void
PiezoActuatorBeam2::AddInternalForces(Vec6& AzLoc)
{
	AzLoc += Vec6(PiezoMat[STRAIN]*V, PiezoMat[CURVAT]*V);
}
      

/* Costruttore normale */
PiezoActuatorBeam2::PiezoActuatorBeam2(unsigned int uL,
		const StructNode* pN1, const StructNode* pN2,
		const Vec3& F1, const Vec3& F2,
		const Mat3x3& R1, const Mat3x3& R2,
		const Mat3x3& r,
		const ConstitutiveLaw6D* pd,
		int iEl,
		const ScalarDifferentialNode **pEDof,
		const Mat3xN& Te, const Mat3xN& Tk,
		OrientationDescription ood,
		flag fOut)
: Elem(uL, fOut),
Beam2(uL, pN1, pN2, F1, F2, R1, R2, r, pd, ood, fOut),
iNumElec(iEl), pvElecDofs(pEDof), V(iEl)
{
	ASSERT(iNumElec > 0);
	ASSERT(pvElecDofs != NULL);

#ifdef DEBUG
	for (int i = iNumElec; i-- > 0; ) {
		ASSERT(pvElecDofs[i] != NULL);
	}   
#endif /* DEBUG */
	
	PiezoMat[STRAIN].Copy(Te);
	PiezoMat[CURVAT].Copy(Tk);
}

   
/* Distruttore banale */
PiezoActuatorBeam2::~PiezoActuatorBeam2(void)
{
	ASSERT(pvElecDofs != NULL);
	SAFEDELETEARR(pvElecDofs);
}

    
/* Contributo al file di restart */
std::ostream&
PiezoActuatorBeam2::Restart(std::ostream& out) const
{
	Restart_(out);
	return out << "/* piezoelectric actuator NOT IMPLEMENTED YET */ ;"
		<< std::endl;
}


/*
 * Dimensioni del workspace; sono 36 righe perche' se genera anche le
 * forze d'inerzia consistenti deve avere accesso alle righe di definizione
 * della quantita' di moto 
 */
void
PiezoActuatorBeam2::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 12;
	*piNumCols = 18+iNumElec;
}


/* Settings iniziali, prima della prima soluzione */
void
PiezoActuatorBeam2::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
	/* se proprio non serve, si puo' eliminare */
	Beam2::SetValue(pDM, X, XP, ph);
}


/* Prepara i parametri di riferimento dopo la predizione */
void
PiezoActuatorBeam2::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	/*
	 * Calcola le deformazioni, aggiorna il legame costitutivo
	 * e crea la FDE
	 */
	for (integer iCnt = 1; iCnt <= iNumElec; iCnt++) {
		V.Put(iCnt, pvElecDofs[iCnt-1]->dGetX());
	}
	
	Beam2::AfterPredict(X, XP);
}


/* assemblaggio jacobiano */
VariableSubMatrixHandler&
PiezoActuatorBeam2::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering PiezoActuatorBeam2::AssJac();"
			" will result in call to AssStiffnessMat()" << std::endl);
	
	integer iNode1FirstMomIndex = pNode[NODE1]->iGetFirstMomentumIndex();
	integer iNode1FirstPosIndex = pNode[NODE1]->iGetFirstPositionIndex();
	integer iNode2FirstMomIndex = pNode[NODE2]->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode[NODE2]->iGetFirstPositionIndex();
	
	FullSubMatrixHandler& WM = WorkMat.SetFull();
	
	/* Dimensiona la matrice, la azzera e pone gli indici corretti */
	WM.ResizeReset(12, 12+iNumElec);
	
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
		WM.PutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
		WM.PutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
	}
	
	for (int iCnt = 1; iCnt <= iNumElec; iCnt++) {
		WM.PutColIndex(12+iCnt,
				pvElecDofs[iCnt-1]->iGetFirstColIndex()+1);
	}   
	
	AssStiffnessMat(WM, WM, dCoef, XCurr, XPrimeCurr);
	
	return WorkMat;        
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
PiezoActuatorBeam2::InitialAssJac(VariableSubMatrixHandler& WorkMat,
				 const VectorHandler& XCurr)
{
	return Beam2::InitialAssJac(WorkMat, XCurr);
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
PiezoActuatorBeam2::InitialAssRes(SubVectorHandler& WorkVec,
				 const VectorHandler& XCurr)
{
	return Beam2::InitialAssRes(WorkVec, XCurr);
}

/* PiezoActuatorBeam2 - end */


/* PiezoActuatorVEBeam2 - begin */

/* Funzioni di calcolo delle matrici */
void
PiezoActuatorVEBeam2::AssStiffnessMat(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering PiezoActuatorVEBeam2::AssStiffnessMat()" << std::endl);
	
	ViscoElasticBeam2::AssStiffnessMat(WMA, WMB, dCoef, XCurr, XPrimeCurr);
	
	Mat3xN tmp1(iNumElec);
	Mat3xN tmp2(iNumElec);
	
	tmp1.LeftMult(R*dCoef, PiezoMat[STRAIN]);
	WMA.Sub(1, 13, tmp1);
	WMA.Add(7, 13, tmp1);
	
	tmp2.LeftMult(Mat3x3(MatCross, p - pNode[NODE1]->GetXCurr()), tmp1);
	WMA.Add(4, 13, tmp2); 

	tmp2.LeftMult(Mat3x3(MatCross, p - pNode[NODE2]->GetXCurr()), tmp1);
	WMA.Sub(10, 13, tmp2);
   
	tmp1.LeftMult(R*dCoef, PiezoMat[CURVAT]);
	WMA.Sub(4, 13, tmp1);
	WMA.Add(10, 13, tmp1);   
}


void
PiezoActuatorVEBeam2::AssStiffnessVec(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering PiezoActuatorVEBeam2::AssStiffnessVec()" << std::endl);
	
	/*
	 * Riceve il vettore gia' dimensionato e con gli indici a posto 
	 * per scrivere il residuo delle equazioni di equilibrio dei due nodi
	 */
	
	/*
	 * Per la trattazione teorica, il riferimento e' il file ul-travi.tex 
	 * (ora e' superato)
	 */
	
	if (bFirstRes) {
		/*
		 * AfterPredict ha gia' calcolato tutto;
		 * bFirstRes viene resettato da
		 * ViscoElasticBeam2::AssStiffnessVec
		 */

	} else {
		for (integer iCnt = 1; iCnt <= iNumElec; iCnt++) {
			V.Put(iCnt, pvElecDofs[iCnt-1]->dGetX());
		}
	}
	
	ViscoElasticBeam2::AssStiffnessVec(WorkVec, dCoef, XCurr, XPrimeCurr);
}


void
PiezoActuatorVEBeam2::AddInternalForces(Vec6& AzLoc)
{
	AzLoc += Vec6(PiezoMat[STRAIN]*V, PiezoMat[CURVAT]*V);
}


/* Costruttore normale */
PiezoActuatorVEBeam2::PiezoActuatorVEBeam2(unsigned int uL,
		const StructNode* pN1, const StructNode* pN2,
		const Vec3& F1, const Vec3& F2,
		const Mat3x3& R1, const Mat3x3& R2,
		const Mat3x3& r,
		const ConstitutiveLaw6D* pD,
		int iEl,
		const ScalarDifferentialNode **pEDof,
		const Mat3xN& Te, const Mat3xN& Tk,
		OrientationDescription ood,
		flag fOut)
: Elem(uL, fOut),
ViscoElasticBeam2(uL, pN1, pN2, F1, F2, R1, R2, r, pD, ood, fOut),
iNumElec(iEl), pvElecDofs(pEDof), V(iEl)
{
	ASSERT(iNumElec > 0);
	ASSERT(pvElecDofs != NULL);

#ifdef DEBUG
	for (int i = iNumElec; i-- > 0; ) {
		ASSERT(pvElecDofs[i] != NULL);
	}
#endif /* DEBUG */
	
	PiezoMat[STRAIN].Copy(Te);
	PiezoMat[CURVAT].Copy(Tk);
}


/* Distruttore banale */
PiezoActuatorVEBeam2::~PiezoActuatorVEBeam2(void)
{
	ASSERT(pvElecDofs != NULL);
	SAFEDELETEARR(pvElecDofs);
}


/* Contributo al file di restart */
std::ostream&
PiezoActuatorVEBeam2::Restart(std::ostream& out) const
{
	Restart_(out);
	return out << "/* piezoelectric actuator NOT IMPLEMENTED YET */ ;"
		<< std::endl;
}


/*
 * Dimensioni del workspace
 */
void
PiezoActuatorVEBeam2::WorkSpaceDim(integer* piNumRows,
		integer* piNumCols) const
{
	*piNumRows = 12;
	*piNumCols = 12+iNumElec;
}


/* Settings iniziali, prima della prima soluzione */
void
PiezoActuatorVEBeam2::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
	/* se proprio non serve, si puo' eliminare */
	ViscoElasticBeam2::SetValue(pDM, X, XP, ph);
}


/* Prepara i parametri di riferimento dopo la predizione */
void
PiezoActuatorVEBeam2::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	/*
	 * Calcola le deformazioni, aggiorna il legame costitutivo
	 * e crea la FDE
	 */
	for (integer iCnt = 1; iCnt <= iNumElec; iCnt++) {
		V.Put(iCnt, pvElecDofs[iCnt-1]->dGetX());
	}
	
	ViscoElasticBeam2::AfterPredict(X, XP);
}


/* assemblaggio jacobiano */
VariableSubMatrixHandler&
PiezoActuatorVEBeam2::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering PiezoActuatorVEBeam2::AssJac();"
			" will result in call to AssStiffnessMat()" << std::endl);
	
	integer iNode1FirstMomIndex = pNode[NODE1]->iGetFirstMomentumIndex();
	integer iNode1FirstPosIndex = pNode[NODE1]->iGetFirstPositionIndex();
	integer iNode2FirstMomIndex = pNode[NODE2]->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode[NODE2]->iGetFirstPositionIndex();
	
	FullSubMatrixHandler& WM = WorkMat.SetFull();
	
	/* Dimensiona la matrice, la azzera e pone gli indici corretti */
	WM.ResizeReset(12, 12+iNumElec);
	
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
		WM.PutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
		WM.PutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
	}
	
	for (int iCnt = 1; iCnt <= iNumElec; iCnt++) {
		WM.PutColIndex(12+iCnt,
				pvElecDofs[iCnt-1]->iGetFirstColIndex()+1);
	}
	
	AssStiffnessMat(WM, WM, dCoef, XCurr, XPrimeCurr);
	
	return WorkMat;        
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
PiezoActuatorVEBeam2::InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr)
{
	return ViscoElasticBeam2::InitialAssJac(WorkMat, XCurr);
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
PiezoActuatorVEBeam2::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr)
{
	return ViscoElasticBeam2::InitialAssRes(WorkVec, XCurr);
}

/* PiezoActuatorVEBeam2 - end */

