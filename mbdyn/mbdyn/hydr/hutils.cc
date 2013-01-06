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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cfloat>
#include <limits>

#include "hutils.h"

/* Accumulator - begin */

Accumulator::Accumulator(unsigned int uL, const DofOwner* pDO, 
			 HydraulicFluid* hf,
			 const PressureNode* p1, 
			 doublereal St, doublereal start, doublereal As, 
			 doublereal A_pipe, doublereal ms,
			 doublereal h_in, doublereal h_out,
			 doublereal P0, doublereal Pmax,doublereal k,
			 doublereal Wg, doublereal Kspr,doublereal F0,
			 doublereal cs, doublereal cv, doublereal ca,
			 flag fOut)
: Elem(uL, fOut),
HydraulicElem(uL, pDO, hf, fOut),
pNode1(p1),
stroke(St), start(start), area(As), area_pipe(A_pipe), mass(ms),
press0(P0), press_max(Pmax),
Kappa(k), weight(Wg), spring(Kspr),force0(F0), h_in(h_in), h_out(h_out),
c_spost(cs), c_vel(cv), c_acc(ca)
{
   ASSERT(pNode1 != NULL);
   ASSERT(pNode1->GetNodeType() == Node::HYDRAULIC);
   ASSERT(St > std::numeric_limits<doublereal>::epsilon());
   ASSERT(As > std::numeric_limits<doublereal>::epsilon()); 
   ASSERT(A_pipe > std::numeric_limits<doublereal>::epsilon()); 
   ASSERT(ms > std::numeric_limits<doublereal>::epsilon());
   ASSERT(P0 >= 0.);
   ASSERT(Pmax >= 0.);
   ASSERT(k >= 0.);
   ASSERT(Wg >= 0.);
   ASSERT(Kspr >= 0.);
   ASSERT(F0 >= 0.);
  
   s_min_gas = 0.;
   
   if (Kappa > 0. && press_max > 0.) {
      s_min_gas = stroke*pow(press0/press_max, 1./Kappa);
   }
   
   s_max = .999*stroke-s_min_gas;  /* ho messo 0.999 perchè se uso un accumulatore senza gas avrei un termine che va a +infinito */    
   ratio2 = area*area/(area_pipe*area_pipe); /* rapporto (area accumulatore/area tubo)^2 */
}


Accumulator::~Accumulator(void)
{
   NO_OP;
}
   
/* Tipo di elemento idraulico (usato solo per debug ecc.) */
HydraulicElem::Type  
Accumulator::GetHydraulicType(void) const 
{
   return HydraulicElem::ACCUMULATOR;
}

/* Contributo al file di restart */
std::ostream&  
Accumulator::Restart(std::ostream& out) const
{
   return out << " Accumulator not implemented yet!" << std::endl;
}


unsigned int  
Accumulator::iGetNumDof(void) const 
{ 
   return 2;
}
   
DofOrder::Order  
Accumulator::GetDofType(unsigned int i) const 
{
   ASSERT(i >= 0 && i <= 1);
   return DofOrder::DIFFERENTIAL;
}

void  
Accumulator::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const 
{
   *piNumRows = 3; 
   *piNumCols = 3; 
}
      
VariableSubMatrixHandler& 
Accumulator::AssJac(VariableSubMatrixHandler& WorkMat,
		    doublereal dCoef,
		    const VectorHandler& XCurr, 
		    const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering  Accumulator::AssJac()" << std::endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.Resize(3, 3);
   
   integer iNode1RowIndex = pNode1->iGetFirstRowIndex()+1;
   integer iNode1ColIndex = pNode1->iGetFirstColIndex()+1;
   integer iFirstIndex = iGetFirstIndex();
   
   WM.PutRowIndex(1, iNode1RowIndex);
   WM.PutColIndex(1, iNode1ColIndex);
   WM.PutRowIndex(2, iFirstIndex+1);
   WM.PutColIndex(2, iFirstIndex+1);
   WM.PutRowIndex(3, iFirstIndex+2);
   WM.PutColIndex(3, iFirstIndex+2);
   doublereal density = HF->dGetDensity();
  
   /* unused? doublereal p1 = pNode1->dGetX(); */
   s = XCurr(iFirstIndex+1);  /* spostamento */
   v = XCurr(iFirstIndex+2);  /* velocita' */
   
   doublereal Jac11 = 0;
   doublereal Jac12 = 0;
   doublereal Jac13 = 0;
   doublereal Jac21 = 0;
   doublereal Jac22 = 0;
   doublereal Jac23 = 0;
   doublereal Jac31 = 0;
   doublereal Jac32 = 0;
   doublereal Jac33 = 0;
   
   if (s > s_max) {
#ifdef HYDR_DEVEL
      DEBUGCOUT("AssJac(): ho superato la s_max: s" << s <<std::endl);
#endif
      s = s_max;
   } else if(s < 0.) {
#ifdef HYDR_DEVEL
      DEBUGCOUT("AssJac(): sono negativo: s" << s <<std::endl);
#endif
      s = 0.;
   }

   Jac11 = 0.;
   Jac12 = 0.;
   Jac13 = -density*area*dCoef;

   Jac21 = area;
   Jac22 = -press0*area*Kappa*dCoef*pow(stroke/(stroke-s), Kappa+1.)/stroke
     -c1*dCoef-cf1*dCoef+dCoef*spring;
   Jac23 = -mass-weight-h*density*area*ratio2*fabs(v)*dCoef
     -c2*dCoef-c3-cf2*dCoef-cf3;
     
   Jac31 = 0.; 
   Jac32 = -1.;
   Jac33 = dCoef;
   
#ifdef HYDR_DEVEL
   DEBUGCOUT("Jac11: " << Jac11 << std::endl);
   DEBUGCOUT("Jac12: " << Jac12 << std::endl);
   DEBUGCOUT("Jac13: " << Jac13 << std::endl);
   DEBUGCOUT("Jac21: " << Jac21 << std::endl);
   DEBUGCOUT("Jac22: " << Jac22 << std::endl);
   DEBUGCOUT("Jac23: " << Jac23 << std::endl);
   DEBUGCOUT("Jac31: " << Jac31 << std::endl);
   DEBUGCOUT("Jac32: " << Jac32 << std::endl);
   DEBUGCOUT("Jac33: " << Jac33 << std::endl);
#endif

   WM.PutCoef(1, 1, Jac11);
   WM.PutCoef(1, 2, Jac12);
   WM.PutCoef(1, 3, Jac13);
   WM.PutCoef(2, 1, Jac21);
   WM.PutCoef(2, 2, Jac22);
   WM.PutCoef(2, 3, Jac23);
   WM.PutCoef(3, 1, Jac31);
   WM.PutCoef(3, 2, Jac32);
   WM.PutCoef(3, 3, Jac33);
   
   return WorkMat;
}


SubVectorHandler& 
Accumulator::AssRes(SubVectorHandler& WorkVec,
		    doublereal dCoef,
		    const VectorHandler& XCurr, 
		    const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering Accumulator::AssRes()" << std::endl);
   
   WorkVec.Resize(3);
  
   integer iNode1RowIndex = pNode1->iGetFirstRowIndex()+1;     
   integer iFirstIndex = iGetFirstIndex();
   
   doublereal p1 = pNode1->dGetX();
   s = XCurr(iFirstIndex+1);        /* spostamento */
   v = XCurr(iFirstIndex+2);        /* velocita' */
   sp = XPrimeCurr(iFirstIndex+1);  /* velocita' */
   vp = XPrimeCurr(iFirstIndex+2);  /* accelerazione */

   doublereal Res_1 = 0.;
   doublereal Res_2 = 0.;
   doublereal Res_3 = 0.;
   pgas = 0.;
   doublereal x0;
   doublereal x0spring;
   doublereal density = HF->dGetDensity();
  
   x0 = (stroke*(pow(p1/press0, 1./Kappa)-1.))/pow(p1/press0, 1./Kappa);
   x0spring = (p1*area-force0)/spring;
      
   if (s < 0.) { 
      c1 = c_spost;
      if (sp < 0.) {
	 c2 = c_vel; /* ho v negativa devo smorzare */
      } else {
	 c2 = 0.;
      }
      
      if (vp < 0.) {
	 c3 = c_acc; /* ho acc negativa devo smorzare */
      } else {
	 c3 = 0.;
      }

      c4 = 0.;  
   } else {
      c1 = 0.;
      c2 = 0.;
      c3 = 0.;
      c4 = 0.;
   }

   if (s > s_max) {
      cf1 = c_spost;
      
      if (sp > 0.) {
	 cf2 = c_vel; /* ho v positiva devo smorzare */
      } else {
	 cf2 = 0.;
      }

      if (vp > 0.) {
	 cf3 = c_acc; /* ho acc positiva devo smorzare */
      } else {
	 cf3 = 0.;
      }
      
      cf4 = 0.;  
   } else {
      cf1 = 0.;
      cf2 = 0.;
      cf3 = 0.;
      cf4 = 0.;
   }

   if (sp > 0.) {
      h = h_in;  /* perdita concentrata fluido entra nell'accumulatore */
   } else {
      h = h_out;  /* perdita concentrata fluido esce dall'accumulatore */
   }

   if (s > s_max) {
      pgas = press_max;
#ifdef HYDR_DEVEL
      DEBUGCOUT("AssJac(): ho superato la s_max: s " << s <<std::endl);
#endif
   } else if (s < 0.) {
      pgas = press0;
#ifdef HYDR_DEVEL
      DEBUGCOUT("AssJac(): sono negativo: s " << s <<std::endl);
#endif
   } else {
      pgas = press0*pow(stroke/(stroke-s), Kappa);
   }

   Res_1 = density*v*area;    
   Res_2 = (mass+weight)*vp-p1*area+pgas*area+force0+spring*s
     +copysign(h*.5*density*area*ratio2*pow(v, 2), v)+c1*s
     +c2*sp+c3*vp+c4+cf1*(s-s_max)+cf2*sp+cf3*vp+cf4;
   Res_3 = sp-v;

   flow = -Res_1;  /* portata nodo 1 (per l'output) */

#ifdef HYDR_DEVEL
   DEBUGCOUT("x0:         " << x0 << std::endl); 
   DEBUGCOUT("x0spring:   " << x0spring << std::endl); 
   DEBUGCOUT("s:          " << s << std::endl); 
   DEBUGCOUT("sp:         " << sp << std::endl); 
   DEBUGCOUT("v:          " << v << std::endl); 
   DEBUGCOUT("vp:         " << vp << std::endl); 
   DEBUGCOUT("p1:         " << p1 << std::endl); 
   DEBUGCOUT("pgas:       " << pgas << std::endl); 
   DEBUGCOUT("smorzatore: " 
	     << copysign(h*.5*density*area*ratio2*pow(v, 2), v) << std::endl); 
   DEBUGCOUT("stroke:     " << stroke << std::endl); 
   DEBUGCOUT("area:       " << area << std::endl); 
   DEBUGCOUT("area_pipe:  " << area_pipe << std::endl); 
   DEBUGCOUT("mass:       " << mass << std::endl); 
   DEBUGCOUT("press0:     " << press0 << std::endl); 
   DEBUGCOUT("press_max:  " << press_max << std::endl); 
   DEBUGCOUT("Kappa:      " << Kappa << std::endl); 
   DEBUGCOUT("weight:     " << weight << std::endl); 
   DEBUGCOUT("spring:     " << spring <<  std::endl); 
   DEBUGCOUT("force0:     " << force0 <<  std::endl); 
   DEBUGCOUT("s_min_gas:  " << s_min_gas << std::endl); 
   DEBUGCOUT("s_max:      " << s_max << std::endl); 
   DEBUGCOUT("c1:         " << c1 << std::endl);
   DEBUGCOUT("c2:         " << c2 << std::endl);
   DEBUGCOUT("c3:         " << c3 << std::endl);
   DEBUGCOUT("c4:         " << c4 << std::endl);
   DEBUGCOUT("cf1:        " << cf1 << std::endl);
   DEBUGCOUT("cf2:        " << cf2 << std::endl);
   DEBUGCOUT("cf3:        " << cf3 << std::endl);
   DEBUGCOUT("cf4:        " << cf4 << std::endl);
   DEBUGCOUT("-Res_1 (portata nodo1): " << -Res_1 << std::endl); 
   DEBUGCOUT("Res_2:                  " << Res_2 << std::endl); 
   DEBUGCOUT("Res_3:                  " << Res_3 << std::endl);
#endif
   
   WorkVec.PutItem(1, iNode1RowIndex, Res_1);
   WorkVec.PutItem(2, iFirstIndex+1, Res_2);         
   WorkVec.PutItem(3, iFirstIndex+2, Res_3);
       
   return WorkVec;
}


void 
Accumulator::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) { 
      std::ostream& out = OH.Hydraulic();
      out << std::setw(8) << GetLabel()
	<< " " << s << " " << v << " " << vp << " " << pgas  
	<< " " << flow << " " << density << std::endl;
   }
}

void 
Accumulator::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{  
   integer i = iGetFirstIndex();
   
   X.PutCoef(i+1, start);
   X.PutCoef(i+2, 0.);
   XP.PutCoef(i+1, 0.);
   XP.PutCoef(i+2, 0.);  
}

/* Accumulator - end */


/* Tank - begin */

Tank::Tank(unsigned int uL, const DofOwner* pDO, HydraulicFluid* hf,
	   const PressureNode* p1, const PressureNode* p2, 
	   doublereal Ps,doublereal A_pipe,
	   doublereal A_serb, doublereal lev, doublereal s_mx, 
	   doublereal s_mn, doublereal c_s, flag fOut) 
: Elem(uL, fOut),
HydraulicElem(uL, pDO, hf, fOut),
pNode1(p1), pNode2(p2),
press(Ps), area_pipe(A_pipe), area_serb(A_serb),level(lev),
s_max(s_mx), s_min(s_mn), c_spost(c_s)
{
   ASSERT(pNode1 != NULL);
   ASSERT(pNode1->GetNodeType() == Node::HYDRAULIC);
   ASSERT(pNode2 != NULL);
   ASSERT(pNode2->GetNodeType() == Node::HYDRAULIC);
   ASSERT(Ps > std::numeric_limits<doublereal>::epsilon());
   ASSERT(A_pipe > std::numeric_limits<doublereal>::epsilon());
   ASSERT(A_serb > std::numeric_limits<doublereal>::epsilon());
   ASSERT(lev >= 0.);
   ASSERT(s_mx >= 0.);
   ASSERT(s_mn >= 0.);
   
   Kappa1 = 1.;
   Kappa2 = .5;  
}


Tank::~Tank(void)
{
   NO_OP;
}


/* Tipo di elemento idraulico (usato solo per debug ecc.) */
HydraulicElem::Type 
Tank::GetHydraulicType(void) const 
{
   return HydraulicElem::TANK;
}


/* Contributo al file di restart */
std::ostream& 
Tank::Restart(std::ostream& out) const
{
   return out << "Tank not implemented yet!" << std::endl;
}


unsigned int 
Tank::iGetNumDof(void) const 
{
   return 1;
}


DofOrder::Order 
Tank::GetDofType(unsigned int i) const 
{
   
   ASSERT(i >= 0 && i <= 1);
   return DofOrder::DIFFERENTIAL;
}


void 
Tank::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
   *piNumRows = 3; 
   *piNumCols = 3; 
}


VariableSubMatrixHandler& 
Tank::AssJac(VariableSubMatrixHandler& WorkMat,
		  doublereal dCoef,
		  const VectorHandler& XCurr, 
		  const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering Tank::AssJac()" << std::endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.Resize(3, 3);
      
   integer iNode1RowIndex = pNode1->iGetFirstRowIndex()+1;
   integer iNode2RowIndex = pNode2->iGetFirstRowIndex()+1;
   integer iNode1ColIndex = pNode1->iGetFirstColIndex()+1;
   integer iNode2ColIndex = pNode2->iGetFirstColIndex()+1;
   integer iFirstIndex = iGetFirstIndex();
 
   WM.PutRowIndex(1, iNode1RowIndex);
   WM.PutRowIndex(2, iNode2RowIndex);
   WM.PutColIndex(1, iNode1ColIndex);
   WM.PutColIndex(2, iNode2ColIndex);

   WM.PutRowIndex(3, iFirstIndex+1);
   WM.PutColIndex(3, iFirstIndex+1);
   
   doublereal p1 = pNode1->dGetX();
   doublereal p2 = pNode2->dGetX(); 
   /* p1 = XCurr(pNode1->iGetFirstRowIndex()+1); */

   s = XCurr(iFirstIndex+1); /* livello */
     
   /* salto di pressione nodo1-Tank */
   doublereal jumpPres1S = fabs(p1-press);  
   /* evito di dividere per un numero troppo piccolo */
   if (jumpPres1S < 1.e8*std::numeric_limits<doublereal>::epsilon()) {
      jumpPres1S = 1.e8*std::numeric_limits<doublereal>::epsilon();
   }
   
   /* salto di pressione Tank-nodo2 */
   doublereal jumpPresS2 = fabs(press-p2);
   /* evito di dividere per un numero troppo piccolo */
   if (jumpPresS2 < 1.e8*std::numeric_limits<doublereal>::epsilon()) {
      jumpPresS2 = 1.e8*std::numeric_limits<doublereal>::epsilon();
   }
   doublereal density = HF->dGetDensity((p1+p2)/2.);
   
   doublereal Jac11 = -.5*density*area_pipe*sqrt(2./(Kappa1*density*jumpPres1S));
   doublereal Jac12 = 0.;
   doublereal Jac13 = 0.;
  
   doublereal Jac21 = 0.;
   doublereal Jac22 = -.5*density*area_pipe*sqrt(2./(Kappa2*density*jumpPresS2));
   doublereal Jac23 = 0.;
  
   doublereal Jac31 = .5*area_pipe/area_serb*sqrt(2./(Kappa1*density*jumpPres1S));
   doublereal Jac32 = .5*area_pipe/area_serb*sqrt(2./(Kappa2*density*jumpPresS2));
   doublereal Jac33 = -1.-c_spost*dCoef;
  
   if(s < s_min) {
      /* Livello dell'olio sotto la soglia minima: 
       * esco dalla presa di emergenza */
#ifdef HYDR_DEVEL
      DEBUGCOUT("Jac Esco dalla presa di emergenza: " << std::endl);
#endif
   }
   if (s > s_max && -flow1 > flow2) {
      /* se e' pieno non puo' entrare di piu' di quella che esce */
      Jac31 = Jac32 = 0.;
#ifdef HYDR_DEVEL
      DEBUGCOUT("Jac Serbatoio  pieno : " << std::endl);
#endif
   }

   if (s < 0. && flow2 > -flow1) {
      /* se e' vuoto non puo' uscire di piu' di quella che entra */
      Jac31 = Jac32 = 0.;
#ifdef HYDR_DEVEL
      DEBUGCOUT("Jac Serbatoio vuoto : "<< std::endl);
#endif
   }

   WM.PutCoef(1, 1, Jac11);
   WM.PutCoef(1, 2, Jac12);
   WM.PutCoef(1, 3, Jac13);
   WM.PutCoef(2, 1, Jac21);
   WM.PutCoef(2, 2, Jac22);
   WM.PutCoef(2, 3, Jac23);
   WM.PutCoef(3, 1, Jac31);
   WM.PutCoef(3, 2, Jac32);
   WM.PutCoef(3, 3, Jac33);
     
   return WorkMat;
}


SubVectorHandler& 
Tank::AssRes(SubVectorHandler& WorkVec,
	     doublereal dCoef,
	     const VectorHandler& XCurr, 
	     const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering Tank::AssRes()" << std::endl);
   
   WorkVec.Resize(3);
   
   integer iNode1RowIndex = pNode1->iGetFirstRowIndex()+1;
   integer iNode2RowIndex = pNode2->iGetFirstRowIndex()+1;
   
   doublereal p1 = pNode1->dGetX();
   doublereal p2 = pNode2->dGetX();
   
   integer iFirstIndex = iGetFirstIndex();
   
   s = XCurr(iFirstIndex+1);       /* livello */
   sp = XPrimeCurr(iFirstIndex+1); /* velocita' del livello */
  
   if (s < 0.) {
      silent_cerr("Tank(" << GetLabel() << ": negative fluid level "
		      << s << " impossible" << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }
     
   doublereal jumpPres1S = fabs(p1-press);
   doublereal jumpPresS2 = fabs(press-p2);
   doublereal Res_1 = 0.;
   doublereal Res_2 = 0.;
   doublereal Res_3 = 0.;
   /* unused? doublereal Res_4 = 0.; */
   c_spost = 0.;
   
   doublereal density = HF->dGetDensity((p1+p2/2.));

   Res_1 =  density*area_pipe*sqrt(2./(Kappa1*density))*copysign(sqrt(jumpPres1S), p1-press);
   Res_2 = -density*area_pipe*sqrt(2./(Kappa2*density))*copysign(sqrt(jumpPresS2), press-p2);

   if (s > s_max && Res_1 > -Res_2) {
      /* se e' pieno non puo' entrare di piu' di quella che esce */
      Res_1 = -Res_2;
      c_spost = 1.;
#ifdef HYDR_DEVEL
      DEBUGCOUT("Res Serbatoio pieno: " << std::endl);
#endif
   }

   if (s <s_min) {
      /* livello dell'olio sotto al livello minimo: 
       * esco dalla presa di emergenza */
#ifdef HYDR_DEVEL
      DEBUGCOUT("Esco dalla presa di emergenza: " << std::endl);
#endif
   }
 
   
   if (s < 0. && Res_2 < -Res_1) {
      c_spost = 0.;
      /* se e' vuoto non puo' uscire di piu' di quella che entra */
      Res_2 = -Res_1;
#ifdef HYDR_DEVEL
      DEBUGCOUT("Res Serbatoio  vuoto: " << std::endl);
#endif
   }
     
   Res_3 = sp+(-Res_1-Res_2)/(area_serb*density)+c_spost*(s-s_max);

   flow1 = -Res_1; /* per l'output */
   flow2 = -Res_2; /* per l'output */
   
#ifdef HYDR_DEVEL
   DEBUGCOUT("Kappa1:     " << Kappa1 << std::endl);
   DEBUGCOUT("Kappa2:     " << Kappa2 << std::endl);
   DEBUGCOUT("jumpPres1S: " << jumpPres1S  << std::endl);
   DEBUGCOUT("jumpPresS2: " << jumpPresS2  << std::endl);
   DEBUGCOUT("density:    " << density << std::endl);
   DEBUGCOUT("p1:         " << p1 << std::endl);
   DEBUGCOUT("p2:         " << p2 << std::endl);
   DEBUGCOUT("press:      " << press << std::endl);
   DEBUGCOUT("s:          " << s << std::endl);
   DEBUGCOUT("level:      " << level << std::endl);
   DEBUGCOUT("s_max:      " << s_max << std::endl);
   DEBUGCOUT("s_min:      " << s_min << std::endl);
   DEBUGCOUT("sp:         " << sp << std::endl);
   DEBUGCOUT("area_pipe:  " << area_pipe << std::endl);
   DEBUGCOUT("area_serb:  " << area_serb << std::endl);
   DEBUGCOUT("PORTATE AI VARI NODI (positive se entranti)" << std::endl);
   DEBUGCOUT("-Res_1:     " << -Res_1 << " (portata nodo1) " << std::endl); 
   DEBUGCOUT("-Res_2:     " << -Res_2 << " (portata nodo2) " << std::endl);   
   DEBUGCOUT("Res_3:      " << Res_3 << std::endl); 
#endif
   
   WorkVec.PutItem(1, iNode1RowIndex, Res_1);
   WorkVec.PutItem(2, iNode2RowIndex, Res_2);         
   WorkVec.PutItem(3, iFirstIndex+1, Res_3);  

   return WorkVec;
}
   
void 
Tank::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) {
      std::ostream& out = OH.Hydraulic();
      out << std::setw(8) << GetLabel();
      out << " " << s << " " << sp << " " << flow1 << " " << flow2 << std::endl;
   }
}

void 
Tank::SetValue(DataManager *pDM,
		VectorHandler& X , VectorHandler& /* XP */ ,
		SimulationEntity::Hints *ph)
{
   integer i = iGetFirstIndex();
   X.PutCoef(i+1, level);
}

/* Tank - end */
