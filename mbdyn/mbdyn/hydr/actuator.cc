/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

/* 
 * Copyright 1999-2000 Lamberto Puggelli <puggelli@tiscalinet.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */

/*
 * Attuatore idraulico.
 * 
 * L'attuatore e' un elemento che applica le pressioni di due nodi idraulici
 * a due nodi strutturali; allo stesso tempo usa la cinematica dei due nodi
 * strutturali per determinare il bilancio di portata dei due nodi idraulici
 * che rappresentano le due camere.
 * 
 * L'elemento attuatore assume che il movimento avvenga lungo la direzione
 * identificata da "axis" nel sistema di riferimento del nodo strutturale 1;
 * quindi sono richiesti due vincoli strutturali di spostamento "inline" e
 * "prismatic" per assicurare che l'unico movimento relativo tra i due nodi
 * strutturali sia di scorrimento lungo un asse, avente la giacitura di "axis"
 * 
 * Non sono previsti trafilamenti tra le due camere, in quanto possono essere 
 * riprodotti mettendo in parallelo un elemento di perdita di carico 
 * concentrata tra i due nodi idraulici.
 * 
 * L'eventuale presenza di attrito ecc. e' delegata ai vincoli strutturali di
 * contorno.
 * 
 * La topologia dell'attuatore e' definita dagli offset "f1" e "f2", 
 * che rappresentano le posizioni, nei sistemi di riferimento 
 * dei rispettivi nodi, dei baricentri delle due superfici che si affacciano 
 * sulla camera 1; quindi i due punti coincidono quando la camera 1 e' a pacco.
 * Tra i dati in ingresso figura anche una dimensione "dl", che rappresenta
 * la massima estensione della camera 1, ovvero la lunghezza del cilindro
 * a meno dello spessore del pistone.
 * la dimensione della camera 2 e' data dall'area 2 per la differenza tra
 * la dimensione massima dell'area 1 e la sua dimensione corrente.
 * 
 * La forza scambiata e' calcolata come (p2*A2-p1*A1) e viene proiettata
 * nella direzione "axis".
 * Se gli offset non sono allineati con "axis", il loro prodotto vettore per 
 * la forza determina la coppia che e' applicata ai nodi strutturali.
 * La portata e' calcolata come la derivata della massa presente in ogni
 * camera, ed e' considerata positiva quando entra nella camera (nodo), quindi
 * 
 *     q = d/dt( rho * V )
 * 
 *       = d( rho )/dp * p' * V + rho * V'
 * 
 * dove i volumi V1 e V2 sono, rispettivamente,
 * 
 *     V1 = A1 * d
 * 
 *     V2 = A2 * ( l - d )
 * 
 * Per disporre delle derivate delle pressioni si sono usati due stati 
 * interni.
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cfloat>
#include <limits>

#include "actuator.h"

#define HF1 HF /* Hack to use the base member HF as HF1 in actuator */

/* Actuator - begin */

Actuator::Actuator(unsigned int uL, const DofOwner* pDO,
		   const PressureNode* p1, const PressureNode* p2,
		   const StructNode* pN1, 
		   const StructNode* pN2,
		   const Vec3& f1Tmp, 
		   const Vec3& f2Tmp,
		   const Vec3& axisTmp,
		   HydraulicFluid* hf1,
		   HydraulicFluid* hf2,
		   doublereal A_1, doublereal A_2, 
		   doublereal l,
		   flag fOut)
: Elem(uL, fOut),
HydraulicElem(uL, pDO, hf1, fOut),
pNodeHyd1(p1), pNodeHyd2(p2),
HF2(hf2),
area1(A_1), area2(A_2),
dl(l), axis(axisTmp),
dp1(0.),
dp2(0.),
dpP1(0.),
dpP2(0.),
flow1(0.),
flow2(0.),
Vol1(0.),
Vol2(0.),
density1(0.),
density2(0.),
pNodeStr1(pN1), pNodeStr2(pN2), f1(f1Tmp), f2(f2Tmp)
{
   ASSERT(pNodeHyd1 != NULL);
   ASSERT(pNodeHyd1->GetNodeType() == Node::HYDRAULIC);
   ASSERT(pNodeHyd2 != NULL);
   ASSERT(pNodeHyd2->GetNodeType() == Node::HYDRAULIC);
   ASSERT(area1 > std::numeric_limits<doublereal>::epsilon()); 
   ASSERT(area2 > std::numeric_limits<doublereal>::epsilon());
   ASSERT(dl > std::numeric_limits<doublereal>::epsilon());
   ASSERT(pNodeStr1 != NULL);
   ASSERT(pNodeStr1->GetNodeType() == Node::STRUCTURAL);
   ASSERT(pNodeStr2 != NULL);
   ASSERT(pNodeStr2->GetNodeType() == Node::STRUCTURAL);
   ASSERT(HF2 != NULL);
}

Actuator::~Actuator(void)
{
   if (HF2 != NULL) {
      SAFEDELETE(HF2);
   }
}
   
/* Tipo di elemento idraulico (usato solo per debug ecc.) */
HydraulicElem::Type Actuator::GetHydraulicType(void) const 
{
   return HydraulicElem::ACTUATOR;
}


/* Contributo al file di restart */
std::ostream& Actuator::Restart(std::ostream& out) const
{
   return out << "Actuator not implemented yet!" << std::endl;
}


unsigned int Actuator::iGetNumDof(void) const 
{
  return 2;
}


DofOrder::Order Actuator::GetDofType(unsigned int i) const 
{   
   return DofOrder::DIFFERENTIAL;
}  


DofOrder::Order Actuator::GetEqType(unsigned int i) const 
{   
   return DofOrder::DIFFERENTIAL;
}  


void 
Actuator::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const 
{   
   *piNumRows = 16;
   *piNumCols = 16;   
    
}


VariableSubMatrixHandler& 
Actuator::AssJac(VariableSubMatrixHandler& WorkMat,
		       doublereal dCoef,
		       const VectorHandler& XCurr, 
		       const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering Actuator::AssJac()" << std::endl);
      
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.ResizeReset(16, 16);   
   
   /* Recupera gli indici */
   integer iNode1FirstPosIndex = pNodeStr1->iGetFirstPositionIndex();
   integer iNode1FirstMomIndex = pNodeStr1->iGetFirstMomentumIndex();
   integer iNode2FirstPosIndex = pNodeStr2->iGetFirstPositionIndex();
   integer iNode2FirstMomIndex = pNodeStr2->iGetFirstMomentumIndex();
   integer iIndex = iGetFirstIndex();
   
   /* Setta gli indici della matrice */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);	
      WM.PutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
      WM.PutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
   }           
    
   integer iNode1RowIndex = pNodeHyd1->iGetFirstRowIndex()+1;
   integer iNode1ColIndex = pNodeHyd1->iGetFirstColIndex()+1;
   integer iNode2RowIndex = pNodeHyd2->iGetFirstRowIndex()+1;
   integer iNode2ColIndex = pNodeHyd2->iGetFirstColIndex()+1;
 
   WM.PutRowIndex(13, iNode1RowIndex);
   WM.PutColIndex(13, iNode1ColIndex);
   WM.PutRowIndex(14, iNode2RowIndex);
   WM.PutColIndex(14, iNode2ColIndex);
   
   WM.PutRowIndex(15, iIndex+1);
   WM.PutColIndex(15, iIndex+1);
   WM.PutRowIndex(16, iIndex+2);
   WM.PutColIndex(16, iIndex+2);

   doublereal p1 = pNodeHyd1->dGetX();
   doublereal p2 = pNodeHyd2->dGetX();
   
   doublereal density1 = HF1->dGetDensity(p1);
   doublereal density2 = HF2->dGetDensity(p2);
   doublereal densityDP1 = HF1->dGetDensityDPres(p1);
   doublereal densityDP2 = HF2->dGetDensityDPres(p2);
   
   Vec3 x1(pNodeStr1->GetXCurr());
   Vec3 x2(pNodeStr2->GetXCurr());
   Mat3x3 R1(pNodeStr1->GetRRef());
   Mat3x3 R2(pNodeStr2->GetRRef());
   
   Vec3 f1Tmp(R1*f1);
   Vec3 f2Tmp(R2*f2);
     
   Vec3 v1(pNodeStr1->GetVCurr());
   Vec3 v2(pNodeStr2->GetVCurr());
   Vec3 Omega1(pNodeStr1->GetWRef());
   Vec3 Omega2(pNodeStr2->GetWRef());
   
   Vec3 TmpAxis(R1*axis);
   
   Vec3 FHyd = TmpAxis*(area1*p1-area2*p2); /* forza idraulica (scalare) */
   
   doublereal dDist = TmpAxis.Dot(x2+f2Tmp-x1-f1Tmp);
   doublereal dVel = TmpAxis.Dot(v2+Omega2.Cross(f2Tmp)-v1-Omega1.Cross(f1Tmp));
   
   /* prima equazione & terza equazione
    * moltiplica deltagpuntonodo1 */
   Mat3x3 JacMatWedge1(MatCross, FHyd*dCoef);
   WM.Sub(1, 4, JacMatWedge1);
   WM.Add(7, 4, JacMatWedge1);
   
   /* moltiplica deltaP1 e deltaP2 */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = TmpAxis.dGet(iCnt);
      doublereal d1 = d*area1;
      doublereal d2 = d*area2;
      
      WM.PutCoef(iCnt, 13, d1);
      WM.PutCoef(6+iCnt, 13, -d1);
      

      WM.PutCoef(iCnt, 14, -d2);
      WM.PutCoef(6+iCnt, 14, d2);
   }

   /* seconda equazione moltiplica deltagpuntonodo1 */
   WM.Sub(4, 4, Mat3x3(MatCross, f1Tmp.Cross(FHyd*dCoef)));
   
   /* moltiplica deltaP1 e deltaP2 */
   Vec3 JacVec = f1Tmp.Cross(TmpAxis);
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = JacVec.dGet(iCnt);
      
      WM.PutCoef(3+iCnt, 13, d*area1);
      WM.PutCoef(3+iCnt, 14, -d*area2);
   }
 
   /* quarta equazione moltiplica deltagpuntonodo1 */
   Mat3x3 TmpMat(MatCross, f2Tmp.Cross(FHyd*dCoef));
   WM.Add(4, 4, TmpMat);
   WM.Sub(4, 10, TmpMat);

   /* moltiplica deltaP1 e deltaP2 */
   JacVec = f2Tmp.Cross(TmpAxis);
   for (int iCnt = 1; iCnt <=3; iCnt++) {
      doublereal d = JacVec.dGet(iCnt);
      
      WM.PutCoef(9+iCnt, 13, -d*area1);
      WM.PutCoef(9+iCnt, 14, d*area2);
   }

   /* inerzia del fluido */
   WM.PutCoef(13, 15, -densityDP1*dDist*area1);
   WM.PutCoef(13, 13, -densityDP1*dVel*area1);
   WM.PutCoef(14, 16, -densityDP2*(dl-dDist)*area2);
   WM.PutCoef(14, 14, densityDP2*dVel*area2);
   
   /* termini in spostamento e velocita' */
   doublereal rho1A1 = (density1+dCoef*densityDP1*dpP1)*area1;
   doublereal rho2A2 = (density2+dCoef*densityDP2*dpP2)*area2;

   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = TmpAxis.dGet(iCnt);
      WM.PutCoef(13, iCnt, d*rho1A1);
      WM.PutCoef(13, 6+iCnt, -d*rho1A1);
      
      WM.PutCoef(14, iCnt, -d*rho2A2);
      WM.PutCoef(14, 6+iCnt, d*rho2A2);
   }

   /* prima eq., nodo 1 */
   JacVec = (TmpAxis.Cross(x2+f2Tmp-x1))*(densityDP1*dpP1*dCoef)
     +(TmpAxis.Cross(f1Tmp+(x2+Omega2.Cross(f2Tmp)-x1)*dCoef))*density1;
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutCoef(13, 3+iCnt, -JacVec.dGet(iCnt)*area1);
   }
   
   /* prima eq., nodo 2 */
   JacVec = (TmpAxis.Cross(f2Tmp))*(densityDP1*dpP1*dCoef)
     +(TmpAxis.Cross(f2Tmp+(Omega2.Cross(f2Tmp))*dCoef))*density1;
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutCoef(13, 9+iCnt, JacVec.dGet(iCnt)*area1);
   }
   
   /* seconda eq., nodo 1 */
   JacVec = (TmpAxis.Cross(x2+f2Tmp-x1))*(densityDP2*dpP2*dCoef)
     +(TmpAxis.Cross(f1Tmp+(x2+Omega2.Cross(f2Tmp)-x1)*dCoef))*density2;
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutCoef(14, 3+iCnt, JacVec.dGet(iCnt)*area2);
   }
   
   /* seconda eq., nodo 2 */
   JacVec = (TmpAxis.Cross(f2Tmp))*(densityDP2*dpP2*dCoef)
     +(TmpAxis.Cross(f2Tmp+(Omega2.Cross(f2Tmp))*dCoef))*density2;
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutCoef(14, 9+iCnt, JacVec.dGet(iCnt)*area2);
   }
   
   /* definizione della pressione interna */
   
   WM.PutCoef(15, 15, dCoef);
   WM.PutCoef(15, 13, -1.);
    
   WM.PutCoef(16, 16, dCoef);
   WM.PutCoef(16, 14, -1.);

   return WorkMat;
}


SubVectorHandler& 
Actuator::AssRes(SubVectorHandler& WorkVec,
		 doublereal dCoef,
		 const VectorHandler& XCurr, 
		 const VectorHandler& XPrimeCurr)

{      
   WorkVec.ResizeReset(16);
   
   integer iNode1RowIndex = pNodeHyd1->iGetFirstRowIndex()+1;
   integer iNode2RowIndex = pNodeHyd2->iGetFirstRowIndex()+1;
   integer iNode1FirstMomIndex = pNodeStr1->iGetFirstMomentumIndex();
   integer iNode2FirstMomIndex = pNodeStr2->iGetFirstMomentumIndex();
   integer iIndex = iGetFirstIndex();
     
   /* Setta gli indici */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.PutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
   }
   
   WorkVec.PutRowIndex(13, iNode1RowIndex);
   WorkVec.PutRowIndex(14, iNode2RowIndex);
   
   WorkVec.PutRowIndex(15, iIndex+1);
   WorkVec.PutRowIndex(16, iIndex+2);
   
   /* Dati */
   doublereal p1 = pNodeHyd1->dGetX();
   doublereal p2 = pNodeHyd2->dGetX();
   
   dp1 = XCurr(iIndex+1);
   dp2 = XCurr(iIndex+2);
   
   dpP1 = XPrimeCurr(iIndex+1);
   dpP2 = XPrimeCurr(iIndex+2);
   
   Vec3 x1(pNodeStr1->GetXCurr());
   Vec3 x2(pNodeStr2->GetXCurr());
   Mat3x3 R1(pNodeStr1->GetRCurr());
   Mat3x3 R2(pNodeStr2->GetRCurr());
     
   Vec3 v1(pNodeStr1->GetVCurr());
   Vec3 v2(pNodeStr2->GetVCurr());
   Vec3 Omega1(pNodeStr1->GetWCurr());
   Vec3 Omega2(pNodeStr2->GetWCurr());
   
   Vec3 f1Tmp(R1*f1);
   Vec3 f2Tmp(R2*f2);
  
   density1 = HF1->dGetDensity(p1);
   density2 = HF2->dGetDensity(p2);
 
   doublereal densityDP1 = HF1->dGetDensityDPres(p1);
   doublereal densityDP2 = HF2->dGetDensityDPres(p2);
   
   /* asse nel sistema globale */
   Vec3 TmpAxis(R1*axis);
   
   /* Calcolo della forza */
   Vec3 FHyd = TmpAxis*(area1*p1-area2*p2);
   
   /* Calcolo delle portate */
   /*
    * FIXME: manca nello Jacobiano un pezzo della velocita' 
    * dell'attuatore dipendente dalla velocita' angolare dell'asse
    */
   Vec3 Dist = x2+f2Tmp-x1-f1Tmp;
   doublereal dDist = TmpAxis*Dist;
   doublereal dVel = TmpAxis.Dot(v2+Omega2.Cross(f2Tmp)-v1-Omega1.Cross(f1Tmp))
	   +Dist*(Omega1.Cross(TmpAxis));
   flow1 = area1*(densityDP1*dpP1*dDist+density1*dVel);
   flow2 = area2*(densityDP2*dpP2*(dl-dDist)-density2*dVel);
   Vol1 = area1*dDist;
   Vol2 = area2*(dl-dDist);

   /* Forze e coppia sul nodo strutturale 1 */
   WorkVec.Sub(1, FHyd);
   WorkVec.Sub(4, f1Tmp.Cross(FHyd));
   
   /* Forze e coppia sul nodo strutturale 2 */
   WorkVec.Add(7, FHyd);
   WorkVec.Add(10, f2Tmp.Cross(FHyd));

   /* Portata sul nodo idraulico 1 */
   WorkVec.PutCoef(13, flow1);
   
   /* Portata sul nodo idraulico 2 */
   WorkVec.PutCoef(14, flow2);
   
   /* pressioni */
   WorkVec.PutCoef(15, p1-dp1);
   WorkVec.PutCoef(16, p2-dp2);
   
   return WorkVec;
}


void Actuator::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) { 
      std::ostream& out = OH.Hydraulic();
      out << std::setw(8) << GetLabel()
	<< " " << flow1  << " " << flow2 
	<< " " << dp1 << " " << dp2
	<< " " << dpP1 << " " << dpP2
	<< " " << Vol1 << " " << Vol2
	<< " " << density1 << " " << density2
	<< std::endl;
   }  
}


void Actuator::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
   integer iIndex = iGetFirstIndex();
   (doublereal&)dp1 = pNodeHyd1->dGetX();
   (doublereal&)dp2 = pNodeHyd2->dGetX();
   (doublereal&)dpP1 = 0.;
   (doublereal&)dpP2 = 0.;
   
   X.PutCoef(iIndex+1, dp1);
   X.PutCoef(iIndex+2, dp2);
   XP.PutCoef(iIndex+1, dpP1);
   XP.PutCoef(iIndex+2, dpP2);
}

/* Actuator - end */

