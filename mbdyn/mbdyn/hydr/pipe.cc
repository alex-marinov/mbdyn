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

/* 
 * Copyright 1999-2000 Lamberto Puggelli <puggelli@tiscalinet.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cfloat>
#include <limits>

#include "pipe.h"

/* Pipe - begin */

Pipe::Pipe(unsigned int uL, const DofOwner* pDO, HydraulicFluid* hf,
		   const PressureNode* p1, const PressureNode* p2,
		   doublereal Dh, doublereal A, doublereal L, flag transition, 
		   doublereal q0, flag fOut)
: Elem(uL, fOut), 
HydraulicElem(uL, pDO, hf, fOut),
pNode1(p1), pNode2(p2),
diameter(Dh), area(A),
length(L), turbulent(transition), q0(q0)
{
   ASSERT(pNode1 != NULL);
   ASSERT(pNode1->GetNodeType() == Node::HYDRAULIC);
   ASSERT(pNode2 != NULL);
   ASSERT(pNode2->GetNodeType() == Node::HYDRAULIC);
   ASSERT(Dh > std::numeric_limits<doublereal>::epsilon());
   ASSERT(A > std::numeric_limits<doublereal>::epsilon());
   ASSERT(L> std::numeric_limits<doublereal>::epsilon());
   viscosity = HF->dGetViscosity();
   
   doublereal density = HF->dGetDensity((pNode1->dGetX()+pNode2->dGetX())/2.);
   
   klam = 32.*length*viscosity/((diameter*diameter)*area*density);
   ktra = .5*length/(diameter*density*area*area); 
   doublereal den = pow(diameter, 1.25)*pow(area, 1.75)*density;
   doublereal num = .1582*pow(viscosity, .25)*length;
   ktrb = num/den;
	
#ifdef HYDR_DEVEL
   DEBUGCOUT("Costruttore Laminare: klam:    " << klam <<std::endl); 
   DEBUGCOUT("Costruttore Transizione: ktra: " << ktra <<std::endl); 
   DEBUGCOUT("Costruttore Turbolento: ktrb:  " << ktrb <<std::endl);
#endif /* HYDR_DEVEL */
}

Pipe::~Pipe(void)
{
   NO_OP;
}
   
/* Tipo di elemento idraulico (usato solo per debug ecc.) */
HydraulicElem::Type Pipe::GetHydraulicType(void) const {
   return HydraulicElem::PIPE;
}

/* Contributo al file di restart */
std::ostream& Pipe::Restart(std::ostream& out) const
{
   return out << "Pipe not implemented yet!" << std::endl;
}
   
unsigned int Pipe::iGetNumDof(void) const { 
   return 1;
}
   
DofOrder::Order Pipe::GetDofType(unsigned int i) const {
   ASSERT(i == 0);
   return DofOrder::ALGEBRAIC;  
}  


void Pipe::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
   *piNumRows = 3; 
   *piNumCols = 3; 
}
     
VariableSubMatrixHandler& 
Pipe::AssJac(VariableSubMatrixHandler& WorkMat,
		  doublereal dCoef,
		  const VectorHandler& XCurr, 
		  const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering Pipe::AssJac()" << std::endl);
#ifdef HYDR_DEVEL
   DEBUGCOUT("dblepsilon INIZIO: "  << std::numeric_limits<doublereal>::epsilon() << std::endl);
#endif /* HYDR_DEVEL */
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.ResizeReset(3, 3);
   
   integer iNode1RowIndex = pNode1->iGetFirstRowIndex()+1;
   integer iNode2RowIndex = pNode2->iGetFirstRowIndex()+1;
   integer iNode1ColIndex = pNode1->iGetFirstColIndex()+1;
   integer iNode2ColIndex = pNode2->iGetFirstColIndex()+1;
   integer iFirstIndex = iGetFirstIndex()+1;
   
   WM.PutRowIndex(1, iNode1RowIndex);
   WM.PutRowIndex(2, iNode2RowIndex);
   WM.PutColIndex(1, iNode1ColIndex);
   WM.PutColIndex(2, iNode2ColIndex);
   WM.PutRowIndex(3, iFirstIndex);
   WM.PutColIndex(3, iFirstIndex);

   doublereal p1 = pNode1->dGetX();
   doublereal p2 = pNode2->dGetX();
   doublereal jumpPres = fabs(p1-p2);
   doublereal q = XCurr(iFirstIndex);  /* portata */
   doublereal Jac13 = -1.;
   doublereal Jac23 = 1.;
   doublereal Jac31 = 1.;
   doublereal Jac32 = -1.;
   doublereal Jac33;
   
   if (Re < HF->dGetRe(HydraulicFluid::LOWER)) {
      /****************************************
       * moto sicuramente laminare (jacobiano)
       ***************************************/
#ifdef HYDR_DEVEL
      DEBUGCOUT("Entering Pipe::AssJac() sono laminare" << std::endl);
#endif /* HYDR_DEVEL */
      Jac33 = -klam;
   }  else if (Re > HF->dGetRe(HydraulicFluid::UPPER)) {
      /*******************************************
       * moto sicuramente turbolento  (jacobiano)
       ******************************************/
#ifdef HYDR_DEVEL
      DEBUGCOUT("Entering Pipe::AssJac() sono turbolento" << std::endl);   
#endif /* HYDR_DEVEL */
      /* evito di dividere per un numero troppo piccolo */
      if (jumpPres < 1.e8*std::numeric_limits<doublereal>::epsilon()) {
	 jumpPres = 1.e8*std::numeric_limits<doublereal>::epsilon();
      }
#ifdef HYDR_DEVEL
      DEBUGCOUT("AssJac() JUMPPRES dopo: " << jumpPres << std::endl);
#endif /* HYDR_DEVEL */
          
      Jac33 = -7./4.*ktrb*pow(fabs(q),3./4.);
      
       } else {
      /***********************************
       * moto di transizione  (jacobiano)
       **********************************/
#ifdef HYDR_DEVEL
      DEBUGCOUT("Entering Pipe::AssJac() sono in transizione" << std::endl);
#endif /* HYDR_DEVEL */
      if (turbulent == 0) {
	 /*******************************************************
	  * moto di transizione laminare-turbolento  (jacobiano)
	  ******************************************************/
#ifdef HYDR_DEVEL
	 DEBUGCOUT("Sono in transizione lam->turb" << std::endl);
#endif /* HYDR_DEVEL */
	 if (Re < HF->dGetRe(HydraulicFluid::LOWER)*1.25) {
	    /* uso lo jacobiano laminare */
	    Jac33 = -klam;
	 } else {
	    doublereal dva = diameter/(viscosity*area);
	    doublereal c = -1.8e-5*dva;
	    doublereal dva2 = dva*dva;
	    doublereal b = 8.e-10*dva2;
	    doublereal dva3 = dva2*dva;
	    doublereal a = 7.e-13*dva3;
	    doublereal d = .0542;

	    doublereal aq = fabs(q);
	    doublereal aq2 = aq*aq;
	    doublereal aq3 = aq2*aq;
	    doublereal aq4 = aq3*aq;
	    
	    Jac33 = -5.*ktra*a*aq4-4.*ktra*b*aq3-3.*ktra*c*aq2-2.*ktra*d*aq;
	 } 
      } else {
	 /*******************************************************
	  * moto di transizione turbolento-laminare  (jacobiano)
	  ******************************************************/
#ifdef HYDR_DEVEL
	 DEBUGCOUT("Sono in transizione turb->lam" << std::endl);    
#endif /* HYDR_DEVEL */
	 if (Re > HF->dGetRe(HydraulicFluid::UPPER)*.775) {
	    /* uso lo jacobiano turbolento per la parte finale */
	    /* evito di dividere per un numero troppo piccolo */
	    if (jumpPres < 1.e8*std::numeric_limits<doublereal>::epsilon()) {
	       jumpPres = 1.e8*std::numeric_limits<doublereal>::epsilon();
	    }
#ifdef HYDR_DEVEL	    
	    DEBUGCOUT("AssJac() JUMPPRES dopo: " << jumpPres << std::endl);
#endif /* HYDR_DEVEL */
	    
	    Jac33 = -7./4.*ktrb*pow(fabs(q),3./4.);
           
      	 } else {
#ifdef HYDR_DEVEL
	    DEBUGCOUT("Jac  turb->lam: TRATTO DI INTERPOLAZIONE" << std::endl);
#endif /* HYDR_DEVEL */
	    doublereal dva = diameter/(viscosity*area);
	    doublereal c = -1.8e-5*dva;
	    doublereal dva2 = dva*dva;
	    doublereal b = 2.e-9*dva2;
	    doublereal dva3 = dva2*dva;
	    doublereal a = 9.e-13*dva3;
	    doublereal d = .0528;
	    
	    doublereal aq = fabs(q);
	    doublereal aq2 = aq*aq;
	    doublereal aq3 = aq2*aq;
	    doublereal aq4 = aq3*aq;
	    
	    Jac33 = -5.*ktra*a*aq4-4.*ktra*b*aq3-3.*ktra*c*aq2-2.*ktra*d*aq;
	 }
      }
   }

#ifdef HYDR_DEVEL
   DEBUGCOUT("JAC Re:        " << Re << std::endl);
   DEBUGCOUT("JAC density:   " << HF->dGetDensity((p1+p2)/2.) << std::endl);
   DEBUGCOUT("JAC p1:        " << p1 << std::endl);
   DEBUGCOUT("JAC p2:        " << p2 << std::endl);
   DEBUGCOUT("JAC jumpPres:  " << jumpPres << std::endl);
   DEBUGCOUT("JAC length:    " << length << std::endl);
   DEBUGCOUT("JAC diameter:  " << diameter << std::endl);
   DEBUGCOUT("JAC viscosity: " << viscosity << std::endl);
   DEBUGCOUT("JAC turbulent: " << turbulent << std::endl);
   DEBUGCOUT("JAC q:         " << q << std::endl);
   DEBUGCOUT("JAC Jac13:     " << Jac13 << std::endl);
   DEBUGCOUT("JAC Jac23:     " << Jac23 << std::endl);
   DEBUGCOUT("JAC Jac31:     " << Jac31 << std::endl);
   DEBUGCOUT("JAC Jac32:     " << Jac32 << std::endl);
   DEBUGCOUT("JAC Jac33:     " << Jac33 << std::endl);
#endif /* HYDR_DEVEL */
   
   WM.PutCoef(1, 3, Jac13);
   WM.PutCoef(2, 3, Jac23);
   WM.PutCoef(3, 1, Jac31); 
   WM.PutCoef(3, 2, Jac32);
   WM.PutCoef(3, 3, Jac33);
  
   return WorkMat;
}

SubVectorHandler& 
Pipe::AssRes(SubVectorHandler& WorkVec,
		 doublereal dCoef,
		 const VectorHandler& XCurr, 
		 const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering Pipe::AssRes()" << std::endl);
   WorkVec.Resize(3);
   
   integer iNode1RowIndex = pNode1->iGetFirstRowIndex()+1;
   integer iNode2RowIndex = pNode2->iGetFirstRowIndex()+1;
   integer iFirstIndex = iGetFirstIndex()+1;
   
   doublereal q = XCurr(iFirstIndex);  /* portata q */
   doublereal p1 = pNode1->dGetX();
   doublereal p2 = pNode2->dGetX();
   doublereal density = HF->dGetDensity((p1+p2)/2.);
   
   doublereal Res_1 = q;
   doublereal Res_2 = -q;
   doublereal Res_3 = 0.;
   
   flow = q;
   vel = flow/(density*area);
   Re = density*fabs(vel)*diameter/viscosity;
   if (Re < HF->dGetRe(HydraulicFluid::LOWER)) {
      /**************************************
       * moto sicuramente laminare (residuo)
       *************************************/
#ifdef HYDR_DEVEL      
      DEBUGCOUT("Entering Pipe::AssRes() SONO LAMINARE" << std::endl);
#endif /* HYDR_DEVEL */
      Res_3 = klam*q-p1+p2;
   } else if (Re > HF->dGetRe(HydraulicFluid::UPPER)) {
      /*****************************************
       * moto sicuramente turbolento  (residuo)
       ****************************************/
#ifdef HYDR_DEVEL
      DEBUGCOUT("Entering Pipe::AssRes() SONO TURBOLENTO" << std::endl);
#endif /* HYDR_DEVEL */
      Res_3= ktrb*copysign(pow(fabs(q),7./4.),q)-p1+p2;  
   } else {
      /*********************************
       * moto di transizione  (residuo)
       ********************************/
#ifdef HYDR_DEVEL
      DEBUGCOUT("SONO IN TRANSIZIONE " << std::endl);
      DEBUGCOUT("Re:  " << Re << std::endl);
      DEBUGCOUT("q:   " << q << std::endl);
      DEBUGCOUT("vel: " << vel << std::endl);
#endif /* HYDR_DEVEL */
      if (turbulent == 0) {
	 /*****************************************************
	  * moto di transizione laminare-turbolento  (residuo)
	  ****************************************************/
#ifdef HYDR_DEVEL
	 DEBUGCOUT("SONO IN TRANSIZIONE lam->turb" << std::endl);
#endif /* HYDR_DEVEL */
	 if (Re < HF->dGetRe(HydraulicFluid::LOWER)*1.25) {
	    /* uso il residuo laminare */
	    Res_3 = klam*q-p1+p2;
#ifdef HYDR_DEVEL
	    DEBUGCOUT("RES lam->turb: TRATTO 64/RE" << std::endl);
	    DEBUGCOUT("Re: " << Re << std::endl);
#endif /* HYDR_DEVEL */
	 } else {
#ifdef HYDR_DEVEL
	    DEBUGCOUT("RES lam->turb:TRATTO DI INTERPOLAZIONE"<< std::endl);
#endif /* HYDR_DEVEL */
	    doublereal dva = diameter/(viscosity*area);
	    doublereal c = -1.8e-5*dva;
	    doublereal dva2 = dva*dva;
	    doublereal b = 8.e-10*dva2;
	    doublereal dva3 = dva2*dva;
	    doublereal a = 7.e-13*dva3;
	    doublereal d = .0542;
	   
#ifdef HYDR_DEVEL	    
	    doublereal fa = ((a*q+b)*q+c)*q+d;
	    DEBUGCOUT("fa lam->turb: " << fa << std::endl);
#endif /* HYDR_DEVEL */
	    doublereal aq = fabs(q);
	    doublereal aq2 = aq*aq;
	    doublereal aq3 = aq2*aq;
	    doublereal aq4 = aq3*aq;
	    
	    Res_3 = -p1+p2+q*(ktra*a*aq4+ktra*b*aq3+ktra*c*aq2+ktra*d*aq);
	 }
      } else {
	 /****************************************************
	  * moto di transizione turbolento-laminare (residuo)
	  ***************************************************/
#ifdef HYDR_DEVEL	 
	 DEBUGCOUT("SONO IN TRANSIZIONE turb->lam" << std::endl);
#endif /* HYDR_DEVEL */
	 if (Re > HF->dGetRe(HydraulicFluid::UPPER)*.775) {
	     /* utilizzo il residuo turbolento */
	     Res_3= ktrb*copysign(pow(fabs(q),7./4.),q)-p1+p2;  
   
#ifdef HYDR_DEVEL
	    DEBUGCOUT("RES turb->lam:TRATTO 0.3164/Re^0.25" << std::endl);
#endif /* HYDR_DEVEL */
	 } else {
#ifdef HYDR_DEVEL
	    DEBUGCOUT("RES turb->lam:TRATTO DI INTERPOLAZIONE"<< std::endl);
#endif /* HYDR_DEVEL */
	    doublereal dva = diameter/(viscosity*area);
	    doublereal c = -1.8e-5*dva;
	    doublereal dva2 = dva*dva;
	    doublereal b = 2.e-9*dva2;
	    doublereal dva3 = dva2*dva;
	    doublereal a = 9.e-13*dva3;
	    doublereal d = .0528;
	    
	    doublereal aq = fabs(q);
	    doublereal aq2 = aq*aq;
	    doublereal aq3 = aq2*aq;
	    doublereal aq4 = aq3*aq;
	    
	    Res_3 = -p1+p2+q*(ktra*a*aq4+ktra*b*aq3+ktra*c*aq2+ktra*d*aq);
	 }
      }
   }

#ifdef HYDR_DEVEL
   DEBUGCOUT("RES density:   " << density << std::endl);
   DEBUGCOUT("RES p1:        " << p1 << std::endl);
   DEBUGCOUT("RES p2:        " << p2 << std::endl);
   DEBUGCOUT("RES jumpPres:  " << fabs(p1-p2) << std::endl);
   DEBUGCOUT("RES length:    " << length << std::endl);
   DEBUGCOUT("RES diameter:  " << diameter << std::endl);
   DEBUGCOUT("RES viscosity: " << viscosity << std::endl);
   DEBUGCOUT("RES area :     " << area << std::endl);
   DEBUGCOUT("RES q:         " << q << std::endl);
   DEBUGCOUT("***********************************************" << std::endl);
   DEBUGCOUT("RES velocita': " << vel << std::endl);
   DEBUGCOUT("    se positiva il fluido va dal nodo 1 al nodo 2" << std::endl);
   DEBUGCOUT("***********************************************" << std::endl);
   DEBUGCOUT("RES Re:        " << Re << std::endl);
   DEBUGCOUT("RES Res_1:     " << Res_1 << std::endl);
   DEBUGCOUT("RES Res_2:     " << Res_2 << std::endl);
   DEBUGCOUT("RES Res_3:     " << Res_3 << std::endl);
#endif /* HYDR_DEVEL */
   
   WorkVec.PutItem(1, iNode1RowIndex, Res_1);	
   WorkVec.PutItem(2, iNode2RowIndex, Res_2);
   WorkVec.PutItem(3, iFirstIndex, Res_3);

   return WorkVec;
}

void
Pipe::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
   if (Re < HF->dGetRe(HydraulicFluid::LOWER)) {
#ifdef HYDR_DEVEL
      DEBUGCOUT("Pipe(" << GetLabel() << "): laminar" << std::endl);
#endif /* HYDR_DEVEL */
      turbulent = 0;
   } else if (Re > HF->dGetRe(HydraulicFluid::UPPER)) {
#ifdef HYDR_DEVEL
      DEBUGCOUT("Pipe(" << GetLabel() << "): turbulent" << std::endl);
#endif /* HYDR_DEVEL */
      turbulent = 1;
#ifdef HYDR_DEVEL
   } else {
      DEBUGCOUT("Pipe(" << GetLabel() << "): transition (" 
		      << turbulent << ")" << std::endl);
#endif /* HYDR_DEVEL */
   }
}
   
void Pipe::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) { 
      OH.Hydraulic()
	<< std::setw(8) << GetLabel()
	<< " " <<  vel << " " << flow  << " " << Re << std::endl;
   }
}


void Pipe::SetValue(DataManager *pDM,
		VectorHandler& X , VectorHandler& /* XP */ ,
		SimulationEntity::Hints *ph)
{
   integer i = iGetFirstIndex();
   X.PutCoef(i+1, q0);  /* portata iniziale nodo 2 */
}

/* Pipe - end */


/* Dynamic_pipe - begin */

Dynamic_pipe::Dynamic_pipe(unsigned int uL, const DofOwner* pDO, HydraulicFluid* hf,
			   const PressureNode* p1, const PressureNode* p2,
			   doublereal Dh, 
			   doublereal A, doublereal L, flag transition, 
			   doublereal q0, flag fOut)
: Elem(uL, fOut), 
HydraulicElem(uL, pDO, hf, fOut),
pNode1(p1), pNode2(p2),
diameter(Dh), area(A),
length(L), turbulent(transition), q0(q0)
{
   ASSERT(pNode1 != NULL);
   ASSERT(pNode1->GetNodeType() == Node::HYDRAULIC);
   ASSERT(pNode2 != NULL);
   ASSERT(pNode2->GetNodeType() == Node::HYDRAULIC);
   ASSERT(Dh > std::numeric_limits<doublereal>::epsilon());
   ASSERT(A > std::numeric_limits<doublereal>::epsilon());
   ASSERT(L > std::numeric_limits<doublereal>::epsilon());
   
   doublereal viscosity = HF->dGetViscosity();
   doublereal density = HF->dGetDensity((pNode1->dGetX()+pNode2->dGetX())/2.);
   
   klam = 8.*length*viscosity/(diameter*diameter);
   ktra = .5/(density*area*diameter);
   ktrb = .5*length*.1582*pow(viscosity, .25)/(pow(diameter, 1.25)*pow(area, .75));
}

Dynamic_pipe::~Dynamic_pipe(void)
{
   NO_OP;
}
   
/* Tipo di elemento idraulico (usato solo per debug ecc.) */
HydraulicElem::Type Dynamic_pipe::GetHydraulicType(void) const {
   return HydraulicElem::DYNAMIC_PIPE;
}

/* Contributo al file di restart */
std::ostream& Dynamic_pipe::Restart(std::ostream& out) const
{
   return out << "Pipe not implemented yet!" << std::endl;
}
   
unsigned int Dynamic_pipe::iGetNumDof(void) const 
{
   return 3;
}
   
DofOrder::Order 
Dynamic_pipe::GetDofType(unsigned int i) const 
{
   ASSERT(i >= 0 && i <= 2);
   return DofOrder::DIFFERENTIAL;
}

void 
Dynamic_pipe::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const 
{
   *piNumRows = 5; 
   *piNumCols = 5; 
}
     
VariableSubMatrixHandler& 
Dynamic_pipe::AssJac(VariableSubMatrixHandler& WorkMat,
	     doublereal dCoef,
	     const VectorHandler& XCurr, 
	     const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering Pipe::AssJac()" << std::endl);
   DEBUGCOUT("Valore di dblepsilon INIZIO:"  << std::numeric_limits<doublereal>::epsilon() << std::endl); 
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.ResizeReset(5, 5);
  
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
   WM.PutRowIndex(4, iFirstIndex+2);
   WM.PutColIndex(4, iFirstIndex+2);
   WM.PutRowIndex(5, iFirstIndex+3);
   WM.PutColIndex(5, iFirstIndex+3);
   
   doublereal p1 = pNode1->dGetX();
   doublereal p2 = pNode2->dGetX();
   
   doublereal q1 = XCurr(iFirstIndex+2);        /* portata nodo 1 */
   doublereal q2 = XCurr(iFirstIndex+3);        /* portata nodo 2 */
   
   doublereal densityDPres = HF->dGetDensityDPres();
   doublereal densityS = HF->dGetDensity(p1);            /* densita' all'inizio del tubo */
   doublereal densityE = HF->dGetDensity(p2);            /* densita' alla fine del nodo */
   doublereal kappa1 = length*area*densityDPres;

   doublereal Jac14 = dCoef;
   doublereal Jac25 = dCoef;
   doublereal Jac31 = -.5;
   doublereal Jac32 = -.5;
   doublereal Jac33 = dCoef;
   doublereal Jac43 = kappa1;
   doublereal Jac44 = dCoef;
   doublereal Jac45 = dCoef;
   doublereal Jac51 = -area+q1*q1/(densityS*densityS*area)*densityDPres;
   doublereal Jac52 = area-q2*q2/(densityE*densityE*area)*densityDPres;
   doublereal Jac54 = -length/2.-2.*dCoef*(q1/(area*densityS));
   doublereal Jac55 = length/2.+2.*dCoef*(q2/(area*densityE));
  
   doublereal lack54;
   doublereal lack55;

   /* usata nel caso turbolento & transizione */
   doublereal kappa2 = length/(6.*diameter*area);
  
   if (Re < HF->dGetRe(HydraulicFluid::LOWER)) {
       /*******************************************
       * moto sicuramente laminare  (jacobiano)
       ******************************************/
      lack54 = dCoef*(-2.*klam/densityS);
      lack55 = dCoef*(2.*klam/densityE);
   } else if (Re > HF->dGetRe(HydraulicFluid::UPPER)) {
      /*******************************************
       * moto sicuramente turbolento  (jacobiano)
       ******************************************/
#ifdef HYDR_DEVEL
      DEBUGCOUT("AssJac() sono turbolento" << std::endl);   
#endif /* HYDR_DEVEL */
      fa = .3164/pow(Re, .25);
       
      lack54 = -dCoef*(fa*kappa2/densityS)*fabs(-q2+2.*q1);
      lack55 = dCoef*(fa*kappa2/densityE)*fabs(2.*q2-q1);    
#ifdef HYDR_DEVEL 
      DEBUGCOUT("fa " << fa << std::endl);
#endif /* HYDR_DEVEL */
   } else {
      /***********************************
       * moto di transizione  (jacobiano)
       **********************************/
#ifdef HYDR_DEVEL
      DEBUGCOUT("AssJac() sono in transizione" << std::endl);
#endif /* HYDR_DEVEL */
      
      if (turbulent == 0) {
	 /*******************************************************
	  * moto di transizione laminare-turbolento  (jacobiano)
	  ******************************************************/
#ifdef HYDR_DEVEL
	 DEBUGCOUT("Sono in transizione lam->turb" << std::endl);
#endif /* HYDR_DEVEL */
	 if (Re < HF->dGetRe(HydraulicFluid::LOWER)*1.25) {
	    /* uso lo jacobiano laminare */	
	    lack54 = dCoef*(-2.*klam/densityS);
	    lack55 = dCoef*(+2.*klam/densityE);
	 } else {
	    doublereal dva = diameter/(viscosity*area);
	    doublereal c = -1.8e-5*dva;
	    doublereal dva2 = dva*dva;
	    doublereal b = 8.e-10*dva2;
	    doublereal dva3 = dva2*dva;
	    doublereal a = 7.e-13*dva3;
	    doublereal d = .0542;
	    
	    doublereal qm1 = fabs(q1);
	    doublereal qm2 = fabs(q2);
	    
	    doublereal fa1 = ((a*qm1+b)*qm1+c)*qm1+d;
	    doublereal fa2 = ((a*qm2+b)*qm2+c)*qm2+d;
#ifdef HYDR_DEVEL	    
	    DEBUGCOUT("JAC fa1: " << fa1 << std::endl);
	    DEBUGCOUT("JAC fa2: " << fa2 << std::endl);
#endif /* HYDR_DEVEL */
	    lack54 = -dCoef*(fa1*length/(6.*diameter*area*densityS))*fabs(-q2+2.*q1);
	    lack55 = dCoef*(fa2*length/(6.*diameter*area*densityE))*fabs(2.*q2-q1);
	 } 
      } else {
	 /*******************************************************
	  * moto di transizione turbolento-laminare  (jacobiano)
	  ******************************************************/
#ifdef HYDR_DEVEL
	 DEBUGCOUT("Sono in transizione turb->lam" << std::endl);
#endif /* HYDR_DEVEL */
	 
	 if (Re > HF->dGetRe(HydraulicFluid::UPPER)*.775) {
	    /* uso lo jacobiano turbolento per la parte finale */
	    fa = .3164/pow(Re, .25);
	    
	    lack54 = -dCoef*(fa*kappa2/densityS)*fabs(-q2+2.*q1);
	    lack55 = dCoef*(fa*kappa2/densityE)*fabs(2.*q2-q1);    
#ifdef HYDR_DEVEL
	    DEBUGCOUT("fa " << fa << std::endl);   
#endif /* HYDR_DEVEL */
	 } else {
#ifdef HYDR_DEVEL
	    DEBUGCOUT("Jac turb->lam: TRATTO DI INTERPOLAZIONE"<< std::endl);
#endif /* HYDR_DEVEL */
	    doublereal dva = diameter/(viscosity*area);
	    doublereal c = -1.8e-5*dva;
	    doublereal dva2 = dva*dva;
	    doublereal b = 2.e-9*dva2;
	    doublereal dva3 = dva2*dva;
	    doublereal a = 9.e-13*dva3;
	    doublereal d = .0528;
	    
	    doublereal qm1 = fabs(q1);
	    doublereal qm2 = fabs(q2);
	    
	    doublereal fa1 = ((a*qm1+b)*qm1+c)*qm1+d;
	    doublereal fa2 = ((a*qm2+b)*qm2+c)*qm2+d;
#ifdef HYDR_DEVEL	 
	    DEBUGCOUT("JAC fa1: " << fa1 << std::endl);
	    DEBUGCOUT("JAC fa2: " << fa2 << std::endl);
#endif /* HYDR_DEVEL */
	    lack54 = -dCoef*(fa1*length/(6.*diameter*area*densityS))*fabs(-q2+2.*q1);
	    lack55 = dCoef*(fa2*length/(6.*diameter*area*densityE))*fabs(2.*q2-q1);
	 }
      }
   } 
   
   Jac54 += lack54;
   Jac55 += lack55;
   
#ifdef HYDR_DEVEL
   DEBUGCOUT("JAC Re:        " << Re << std::endl);
   DEBUGCOUT("JAC density:   " << HF->dGetDensity() << std::endl);
   DEBUGCOUT("JAC p1:        " << p1 << std::endl);
   DEBUGCOUT("JAC p2:        " << p2 << std::endl);
   DEBUGCOUT("JAC q1:        " << q1 << std::endl);
   DEBUGCOUT("JAC q2:        " << q2 << std::endl);
   doublereal q1p = XPrimeCurr(iFirstIndex+2); /* derivata q nodo 1 */
   doublereal q2p = XPrimeCurr(iFirstIndex+3); /* derivata q nodo 2 */
   DEBUGCOUT("JAC q1p:       " << q1p << std::endl);
   DEBUGCOUT("JAC q2p:       " << q2p << std::endl);
   DEBUGCOUT("JAC length:    " << length << std::endl);
   DEBUGCOUT("JAC diameter:  " << diameter << std::endl);
   DEBUGCOUT("JAC viscosity: " << viscosity << std::endl);
   DEBUGCOUT("JAC turbulent: " << turbulent << std::endl);
   DEBUGCOUT("JAC Jac14:     " << Jac14 << std::endl);
   DEBUGCOUT("JAC Jac25:     " << Jac25 << std::endl);
   DEBUGCOUT("JAC Jac31:     " << Jac31 << std::endl);
   DEBUGCOUT("JAC Jac32:     " << Jac32 << std::endl);
   DEBUGCOUT("JAC Jac33:     " << Jac33 << std::endl);
   DEBUGCOUT("JAC Jac43:     " << Jac43 << std::endl);
   DEBUGCOUT("JAC Jac44:     " << Jac44 << std::endl);
   DEBUGCOUT("JAC Jac45:     " << Jac45 << std::endl);
   DEBUGCOUT("JAC Jac51:     " << Jac51 << std::endl);
   DEBUGCOUT("JAC Jac52:     " << Jac52 << std::endl);
   DEBUGCOUT("JAC Jac54:     " << Jac54 << std::endl);
   DEBUGCOUT("JAC Jac55:     " << Jac55 << std::endl);
#endif /* HYDR_DEVEL */
   
   WM.PutCoef(1, 4, Jac14);
   WM.PutCoef(2, 5, Jac25);
   WM.PutCoef(3, 1, Jac31);
   WM.PutCoef(3, 2, Jac32);
   WM.PutCoef(3, 3, Jac33);
   WM.PutCoef(4, 3, Jac43);
   WM.PutCoef(4, 4, Jac44);
   WM.PutCoef(4, 5, Jac45);
   WM.PutCoef(5, 1, Jac51);   
   WM.PutCoef(5, 2, Jac52);
   WM.PutCoef(5, 4, Jac54);
   WM.PutCoef(5, 5, Jac55);
 
   return WorkMat;
}


SubVectorHandler& 
Dynamic_pipe::AssRes(SubVectorHandler& WorkVec,
	     doublereal dCoef,
	     const VectorHandler& XCurr, 
	     const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering Dynamic_pipe::AssRes()" << std::endl);
   WorkVec.Resize(5);
   
   integer iNode1RowIndex = pNode1->iGetFirstRowIndex()+1;
   integer iNode2RowIndex = pNode2->iGetFirstRowIndex()+1;
   
   doublereal p1 = pNode1->dGetX();
   doublereal p2 = pNode2->dGetX();
   integer iFirstIndex = iGetFirstIndex();
   
   doublereal pr = XCurr(iFirstIndex+1);        /* pressione */
   doublereal prp = XPrimeCurr(iFirstIndex+1);  /* derivata pressione */
   
   doublereal q1 = XCurr(iFirstIndex+2);        /* portata nodo 1 */
   doublereal q2 = XCurr(iFirstIndex+3);        /* portata nodo 2 */
   doublereal q1p = XPrimeCurr(iFirstIndex+2);  /* derivata portata nodo 1 */
   doublereal q2p = XPrimeCurr(iFirstIndex+3);  /* derivata portata nodo 2 */
   
   flow1 = q1;  /* per l'output */
   flow2 = q2;  /* per l'output */
   pp = prp;    /* per l'output */
  
   doublereal densityS = HF->dGetDensity(p1);           /* densita' all'inizio del tubo */
   doublereal densityE = HF->dGetDensity(p2);           /* densita' alla fine del nodo */
   doublereal densityM = HF->dGetDensity(pr);           /* densita' a meta' tubo */
   
   
   densitySt = densityS;  /* per l'output */
   densityEn = densityE;  /* per l'output */
   densityMe = densityM;  /* per l'output */
   
   viscosity = HF->dGetViscosity();
   densityDPres = HF->dGetDensityDPres(); 
   
   doublereal kappa1 = area*length*densityDPres;
   
   const doublereal x1 = -1./(sqrt(3.));
   const doublereal x2 = -x1;
   doublereal Qx1 = (q1+q2)*x1-q1+q2;
   doublereal Qx2 = (q1+q2)*x2-q1+q2;
   
   doublereal densityx1 = HF->dGetDensity(.5*(p2-p1)*x1+pr);
   doublereal densityx2 = HF->dGetDensity(.5*(p2-p1)*x2+pr);
		 		 
   doublereal Res_1 = -q1;
   doublereal Res_2 = -q2;
   doublereal Res_3 = .5*(p1+p2)-pr;
   doublereal Res_4 = -q2-q1-kappa1*prp;
   doublereal Res_5 = (length/2.)*(q1p-q2p)+q1*q1/(area*densityS)-q2*q2/(area*densityE)+area*(p1-p2);
   doublereal lack;
   
   if (Re < HF->dGetRe(HydraulicFluid::LOWER)) {
#ifdef HYDR_DEVEL
     DEBUGCOUT("SONO IN LAMINARE" << std::endl);
#endif /* HYDR_DEVEL */
     lack = -klam*(Qx1/densityx1+Qx2/densityx2);
   
   } else if (Re > HF->dGetRe(HydraulicFluid::UPPER)) {
#ifdef HYDR_DEVEL
      DEBUGCOUT("SONO IN TURBOLENTO" << std::endl);
#endif /* HYDR_DEVEL */
      const doublereal espo = 7./4.;

      lack = -ktrb*(copysign(pow(.5*fabs(Qx1), espo)/densityx1, Qx1)
		    +copysign(pow(.5*fabs(Qx2), espo)/densityx2, Qx2));
#ifdef HYDR_DEVEL
      DEBUGCOUT("SONO IN TURBOLENTO lack: " << lack << std::endl);
      DEBUGCOUT("SONO IN TURBOLENTO Qx1:  " << Qx1 << std::endl);
      DEBUGCOUT("SONO IN TURBOLENTO Qx2:  " << Qx2 << std::endl);
#endif /* HYDR_DEVEL */
   } else { 
      /*********************************
       * moto di transizione  (residuo)
       ********************************/
#ifdef HYDR_DEVEL
      DEBUGCOUT("SONO IN TRANSIZIONE" << std::endl);
      DEBUGCOUT("Re: " << Re << std::endl);
#endif /* HYDR_DEVEL */
      if (turbulent == 0) {  
	 /****************************************************
	  * moto di transizione laminare-turbolento (residuo)
	  ***************************************************/
#ifdef HYDR_DEVEL
	 DEBUGCOUT("SONO IN TRANSIZIONE lam->turb" << std::endl);
#endif /* HYDR_DEVEL */
	 if (Re < HF->dGetRe(HydraulicFluid::LOWER)*1.25) {
	    /* uso il residuo laminare */
	    lack = -klam*(Qx1/densityx1+Qx2/densityx2);
#ifdef HYDR_DEVEL
	    DEBUGCOUT("RES lam->turb: TRATTO 64/RE" << std::endl);
#endif /* HYDR_DEVEL */
	 } else {
#ifdef HYDR_DEVEL
	    DEBUGCOUT("RES lam->turb: TRATTO DI INTERPOLAZIONE"<< std::endl);
#endif /* HYDR_DEVEL */
	    doublereal dva = diameter/(viscosity*area);
	    doublereal c = -1.8e-5*dva;
	    doublereal dva2 = dva*dva;
	    doublereal b = 8.e-10*dva2;
	    doublereal dva3 = dva2*dva;
	    doublereal a = 7.e-13*dva3;
	    doublereal d = .0542;
	     
	    doublereal qx1 = fabs(-q1+((q2+q1)/length)*x1); /* usata solo per calcolare il giusto fa in x1 */
	    doublereal qx2 = fabs(-q1+((q2+q1)/length)*x2); /* usata solo per calcolare il giusto fa in x2 */
	    doublereal fax1 = ((a*qx1+b)*qx1+c)*qx1+d;
	    doublereal fax2 = ((a*qx2+b)*qx2+c)*qx2+d;

	    doublereal kappa3 = length/(16.*diameter*area);
#ifdef HYDR_DEVEL
	    DEBUGCOUT("RES lam->turb fax1: " << fax1 << std::endl);
	    DEBUGCOUT("RES lam->turb fax2: " << fax2 << std::endl);
#endif /* HYDR_DEVEL */
	    lack = -kappa3*(fax1*copysign(Qx1*Qx1, Qx1)/densityx1
			    +fax2*copysign(Qx2*Qx2, Qx2)/densityx2);
	 } 
      } else {
	 /****************************************************
	  * moto di transizione turbolento-laminare (residuo)
	  ***************************************************/
#ifdef HYDR_DEVEL
	 DEBUGCOUT("SONO IN TRANSIZIONE turb->lam" << std::endl);
#endif /* HYDR_DEVEL */
	 if (Re > HF->dGetRe(HydraulicFluid::UPPER)*.775) {
	    /* utilizzo il residuo turbolento */
	    const doublereal espo = 7./4.;

	    lack = -ktrb*(copysign(pow(.5*fabs(Qx1), espo)/densityx1, Qx1)
			  +copysign(pow(.5*fabs(Qx2), espo)/densityx2, Qx2));
  
	    /* DEBUGCOUT("RES turb->lam: TRATTO 0.3164/Re^0.25" << std::endl); */
	 } else {
#ifdef HYDR_DEVEL
	    DEBUGCOUT("RES turb->lam:TRATTO DI INTERPOLAZIONE" << std::endl);
#endif /* HYDR_DEVEL */
	    doublereal dva = diameter/(viscosity*area);
	    doublereal c = -1.8e-5*dva;
	    doublereal dva2 = dva*dva;
	    doublereal b = 2.e-9*dva2;
	    doublereal dva3 = dva2*dva;
	    doublereal a = 9.e-13*dva3;
	    doublereal d = .0528;

	    doublereal qx1 = fabs(-q1+((q2+q1)/length)*x1); /* usata solo per calcolare il giusto fa in x1 */
	    doublereal qx2 = fabs(-q1+((q2+q1)/length)*x2); /* usata solo per calcolare il giusto fa in x2 */
	    doublereal fax1 = ((a*qx1+b)*qx1+c)*qx1+d;
	    doublereal fax2 = ((a*qx2+b)*qx2+c)*qx2+d;
#ifdef HYDR_DEVEL
	    DEBUGCOUT("RES turb->lam fax1: " << fax1 << std::endl);
	    DEBUGCOUT("RES turb->lam fax2: " << fax2 << std::endl);
#endif /* HYDR_DEVEL */
	    doublereal kappa3 = length/(16.*diameter*area);

	    lack = -kappa3*(fax1*copysign(Qx1*Qx1, Qx1)/densityx1
			    +fax2*copysign(Qx2*Qx2, Qx2)/densityx2);
	 }
      }
   }
   
   Res_5 += lack;                                 /* aggiungo la parte dovuta all'attrito */
   VelS = -q1/(area*densityS);                    /* velocita' inizio */
   VelM = ((-q1+q2)/2.)/(area*densityM);          /* velocita' media */
   VelE = q2/(area*densityE);                     /* velocita' fine */
   Re = densityM*fabs(VelM)*diameter/viscosity;   /* numero di Reynolds medio */
   
#ifdef HYDR_DEVEL   
   DEBUGCOUT("RES density:          " << HF->dGetDensity() << std::endl);
   DEBUGCOUT("RES densityS:         " << densityS << std::endl);
   DEBUGCOUT("RES densityM:         " << densityM << std::endl);
   DEBUGCOUT("RES densityE:         " << densityE << std::endl);
   DEBUGCOUT("RES p1 :              " << p1 << std::endl);
   DEBUGCOUT("RES p2 :              " << p2 << std::endl);
   DEBUGCOUT("RES pr :              " << pr << std::endl);
   DEBUGCOUT("RES prp:              " << prp << std::endl);
   DEBUGCOUT("RES q1:               " << q1 << std::endl);
   DEBUGCOUT("RES q2:               " << q2 << std::endl);
   DEBUGCOUT("RES q1p:              " << q1p << std::endl);
   DEBUGCOUT("RES q2p:              " << q2p << std::endl);
   DEBUGCOUT("RES length:           " << length << std::endl);
   DEBUGCOUT("RES diameter:         " << diameter << std::endl);
   DEBUGCOUT("RES viscosity:        " << viscosity << std::endl);
   DEBUGCOUT("RES area:             " << area << std::endl);
   DEBUGCOUT("RES Re:               " << Re << std::endl);
   DEBUGCOUT("***********************************************"<< std::endl);
   DEBUGCOUT("RES velocita' Start:  " << VelS << std::endl);
   DEBUGCOUT("RES velocita' Middle: " << VelM << std::endl);
   DEBUGCOUT("RES velocita' End:    " << VelE << std::endl);
   DEBUGCOUT("    se positiva il fluido va dal nodo 1 al nodo 2 " << std::endl);
   DEBUGCOUT("***********************************************"<< std::endl);
   DEBUGCOUT("RES Res_1:            " << Res_1 << std::endl);
   DEBUGCOUT("RES Res_2:            " << Res_2 << std::endl);
   DEBUGCOUT("RES Res_3:            " << Res_3 << std::endl);
   DEBUGCOUT("RES Res_4:            " << Res_4 << std::endl);
   DEBUGCOUT("RES Res_5:            " << Res_5 << std::endl);
#endif   /* HYDR_DEVEL */ 
   
   WorkVec.PutItem(1, iNode1RowIndex, Res_1);	
   WorkVec.PutItem(2, iNode2RowIndex, Res_2);
   WorkVec.PutItem(3, iFirstIndex+1, Res_3);
   WorkVec.PutItem(4, iFirstIndex+2, Res_4);    
   WorkVec.PutItem(5, iFirstIndex+3, Res_5);
   
   return WorkVec;
}

void
Dynamic_pipe::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
   if (Re < HF->dGetRe(HydraulicFluid::LOWER)) {
      turbulent = 0;
   } else if (Re > HF->dGetRe(HydraulicFluid::UPPER)) {
      turbulent = 1;
   }
   
#ifdef HYDR_DEVEL
   DEBUGCOUT("Dynamic_Pipe(" << GetLabel() << "): turbulent mode =  " 
		   << turbulent << std::endl);
#endif /* HYDR_DEVEL */
}

void Dynamic_pipe::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) { 
      std::ostream& out = OH.Hydraulic();
      out 
	<< std::setw(8) << GetLabel()
        << " " << Re << " " << -flow1 << " " << flow2
	<< " " << densitySt << " " << densityMe << " " << densityEn
	<< " " << VelS << " " << VelM << " " << VelE
	<< " " << pp
	<< std::endl;
   }
}

void 
Dynamic_pipe::SetValue(DataManager *pDM,
		VectorHandler& X , VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
   integer i = iGetFirstIndex();
   
   doublereal p1 = pNode1->dGetX();
   doublereal p2 = pNode2->dGetX();
   
   X.PutCoef(i+1, .5*(p1+p2));
   X.PutCoef(i+2, q0);
   X.PutCoef(i+3, -q0);
   
   XP.PutCoef(i+1, 0.);
   XP.PutCoef(i+2, 0.);
   XP.PutCoef(i+3, 0.);
}

/* Dynamic_pipe - end */


/* DynamicPipe - begin */

DynamicPipe::DynamicPipe(unsigned int uL, 
			 const DofOwner* pDO, 
			 HydraulicFluid* hf,
			 const PressureNode* p1, 
			 const PressureNode* p2,
			 doublereal Dh, 
			 doublereal A, 
			 doublereal L, 
			 flag transition, 
			 doublereal q0, 
			 flag fOut)
: Elem(uL, fOut), 
HydraulicElem(uL, pDO, hf, fOut),
pNode1(p1), 
pNode2(p2),
diameter(Dh), 
area(A),
length(L), 
turbulent(transition), 
q0(q0)
{
   ASSERT(pNode1 != NULL);
   ASSERT(pNode1->GetNodeType() == Node::HYDRAULIC);
   ASSERT(pNode2 != NULL);
   ASSERT(pNode2->GetNodeType() == Node::HYDRAULIC);
   
   ASSERT(Dh > std::numeric_limits<doublereal>::epsilon());
   ASSERT(A > std::numeric_limits<doublereal>::epsilon());
   ASSERT(L > std::numeric_limits<doublereal>::epsilon());
   
   dKlam = 16.*length/(diameter*diameter);
   dKtra = .25*length/(diameter*area);
   dKtrb = .5*length*.1582/(pow(diameter, 1.25)*pow(area, .75));
}

DynamicPipe::~DynamicPipe(void)
{
   NO_OP;
}
   
/* Tipo di elemento idraulico (usato solo per debug ecc.) */
HydraulicElem::Type DynamicPipe::GetHydraulicType(void) const {
   return HydraulicElem::DYNAMIC_PIPE;
}

/* Contributo al file di restart */
std::ostream& DynamicPipe::Restart(std::ostream& out) const
{
   return out << "Pipe not implemented yet!" << std::endl;
}
   
unsigned int DynamicPipe::iGetNumDof(void) const 
{
   return 4;
}
   
DofOrder::Order 
DynamicPipe::GetDofType(unsigned int i) const 
{
   ASSERT(i >= 0 && i < 4);
   return DofOrder::DIFFERENTIAL;
}

DofOrder::Order 
DynamicPipe::GetEqType(unsigned int i) const 
{
   ASSERT(i >= 0 && i < 4);
   return DofOrder::DIFFERENTIAL;
}

void 
DynamicPipe::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const 
{
   *piNumRows = 6;
   *piNumCols = 6; 
}
     
VariableSubMatrixHandler& 
DynamicPipe::AssJac(VariableSubMatrixHandler& WorkMat,
	     doublereal dCoef,
	     const VectorHandler& XCurr, 
	     const VectorHandler& XPrimeCurr)
{
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.ResizeReset(6, 6);
  
   integer iNode1RowIndex = pNode1->iGetFirstRowIndex()+1;
   integer iNode2RowIndex = pNode2->iGetFirstRowIndex()+1;
   integer iNode1ColIndex = pNode1->iGetFirstColIndex()+1;
   integer iNode2ColIndex = pNode2->iGetFirstColIndex()+1;
   
   integer iFirstIndex = iGetFirstIndex();  
   
   WM.PutRowIndex(1, iNode1RowIndex);
   WM.PutRowIndex(2, iNode2RowIndex);
   WM.PutRowIndex(3, iFirstIndex+1);
   WM.PutRowIndex(4, iFirstIndex+2);
   WM.PutRowIndex(5, iFirstIndex+3);
   WM.PutRowIndex(6, iFirstIndex+4);
   
   WM.PutColIndex(1, iNode1ColIndex);
   WM.PutColIndex(2, iNode2ColIndex);
   WM.PutColIndex(3, iFirstIndex+1);
   WM.PutColIndex(4, iFirstIndex+2);
   WM.PutColIndex(5, iFirstIndex+3);
   WM.PutColIndex(6, iFirstIndex+4);

   doublereal dRDP1 = HF->dGetDensityDPres(p1);
   doublereal dRDP2 = HF->dGetDensityDPres(p2);
   doublereal dRDP0 = HF->dGetDensityDPres(.5*(p1+p2));

   doublereal qq1 = .75*q1+.25*q2;
   doublereal qq2 = .25*q1+.75*q2;
   
   doublereal dd1 = .5*(density1+density0);
   doublereal dd2 = .5*(density0+density2);      
   
   doublereal dLoss11 = 0.;
   doublereal dLoss12 = 0.;
   doublereal dLoss21 = 0.;
   doublereal dLoss22 = 0.;
   
   if ((Re < HF->dGetRe(HydraulicFluid::LOWER)) || ((turbulent == 0) 
			   && (Re < 1.25*HF->dGetRe(HydraulicFluid::LOWER)))) {
      /* laminar, or ascending transition ascendente (continuation) */
      doublereal dk = dKlam*viscosity*dCoef;

      dLoss11 = .75*dk/dd1;
      dLoss12 = .25*dk/dd1;
      dLoss21 = .25*dk/dd2;
      dLoss22 = .75*dk/dd2;
   } else if ((Re > HF->dGetRe(HydraulicFluid::UPPER)) || ((turbulent == 1) 
			   && (Re > .775*HF->dGetRe(HydraulicFluid::UPPER)))) {
      /* turbulent, or descending transition (continuation) */
      
      doublereal dk = dKtrb*1.75*pow(viscosity, .25)*dCoef;
      doublereal dl1 = dk*pow(fabs(qq1), .75)/dd1;
      doublereal dl2 = dk*pow(fabs(qq2), .75)/dd2;
      
      dLoss11 = .75*dl1;
      dLoss12 = .25*dl1;
      dLoss21 = .25*dl2;
      dLoss22 = .75*dl2;
   } else {
      /* transition */
      doublereal dva = diameter/(viscosity*area);
      doublereal dva2 = dva*dva;
      doublereal dva3 = dva2*dva;
      
      doublereal a, b, c, d;
      
      if (turbulent == 0) {
	 /* ascending */
	 a = 7.e-13*dva3;
	 b = 8.e-10*dva2;
	 c = -1.8e-5*dva;
	 d = .0542;
      } else /* if (turbulent == 1) */ {
	 /* descending */
	 a = 9.e-13*dva3;
	 b = 2.e-9*dva2;
	 c = -1.8e-5*dva;
	 d = .0528;
      }
	 
      doublereal fqq1 = fabs(qq1);
      doublereal fqq2 = fabs(qq2);
      
      doublereal fax1 = (((5.*a*fqq1+4.*b)*fqq1+3.*c)*fqq1+2.*d)*fqq1;
      doublereal fax2 = (((5.*a*fqq2+4.*b)*fqq2+3.*c)*fqq2+2.*d)*fqq2;
      
      doublereal dk = dKtrb*dCoef;
      
      doublereal dl1 = dk*fax1/dd1;
      doublereal dl2 = dk*fax2/dd2;
      
      dLoss11 = .75*dl1;
      dLoss12 = .25*dl1;
      dLoss21 = .25*dl2;
      dLoss22 = .75*dl2;
   }
   
   /* primo blocco: conservazione massa */
   doublereal dr = .5*dCoef;
   WM.PutCoef(1, 5, -dr);
   WM.PutCoef(1, 6, -dr);
   WM.PutCoef(2, 5, dr);
   WM.PutCoef(2, 6, dr);
 
   dr = area*length/8.;
   WM.PutCoef(1, 3, -dr*3.*densityDPres1);
   WM.PutCoef(1, 4, -dr*densityDPres1);
   WM.PutCoef(2, 3, -dr*densityDPres2);
   WM.PutCoef(2, 4, -dr*3.*densityDPres2);
   
   
   /* secondo blocco: equazione quantita' di moto */
   dr = .5*dCoef*area;
   doublereal q12 = .5*(q1+q2);
   doublereal ddq12 = .5*(q12*q12)/(density0*density0*area)*dRDP0*dCoef;
   doublereal ddq1 = q1*q1/(density1*density1*area)*dRDP1*dCoef;
   doublereal ddq2 = q2*q2/(density2*density2*area)*dRDP2*dCoef;
   WM.PutCoef(3, 3, -dr + ddq1 - ddq12);
   WM.PutCoef(3, 4, dr - ddq12);
   WM.PutCoef(4, 3, -dr + ddq12);
   WM.PutCoef(4, 4, dr + ddq12 - ddq2);
   
   /* manca la viscosita' */
   dr = length/8.;
   doublereal dq12 = .5*(q1+q2)/density0;   
   WM.PutCoef(3, 5, 3.*dr+(dq12-2.*q1/density1)/area*dCoef  + dLoss11);
   WM.PutCoef(3, 6, dr+dq12/area*dCoef + dLoss12);
   WM.PutCoef(4, 5, dr-dq12/area*dCoef + dLoss21);
   WM.PutCoef(4, 6, 3.*dr+(2.*q2/density2-dq12)/area*dCoef + dLoss22);
   
   /* terzo blocco: definizione delle pressioni nodali interne */   
   WM.PutCoef(5, 1, -1.);
   WM.PutCoef(6, 2, -1.);
   
   WM.PutCoef(5, 3, dCoef);
   WM.PutCoef(6, 4, dCoef);
   
   return WorkMat;
}


SubVectorHandler& 
DynamicPipe::AssRes(SubVectorHandler& WorkVec,
		    doublereal dCoef,
		    const VectorHandler& XCurr, 
		    const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering DynamicPipe::AssRes()" << std::endl);
   WorkVec.Resize(6);
   
   integer iNode1RowIndex = pNode1->iGetFirstRowIndex()+1;
   integer iNode2RowIndex = pNode2->iGetFirstRowIndex()+1;
   
   doublereal pn1 = pNode1->dGetX();
   doublereal pn2 = pNode2->dGetX();
   
   integer iFirstIndex = iGetFirstIndex();
   
   p1 = XCurr(iFirstIndex+1);        /* pressione */
   p2 = XCurr(iFirstIndex+2);        /* pressione */
   p1p = XPrimeCurr(iFirstIndex+1);  /* derivata pressione */
   p2p = XPrimeCurr(iFirstIndex+2);  /* derivata pressione */
   
   q1 = XCurr(iFirstIndex+3);        /* portata nodo 1 */
   q2 = XCurr(iFirstIndex+4);        /* portata nodo 2 */
   q1p = XPrimeCurr(iFirstIndex+3);  /* derivata portata nodo 1 */
   q2p = XPrimeCurr(iFirstIndex+4);  /* derivata portata nodo 2 */
   
   doublereal p0 = .5*(p1+p2);
   
   density1 = HF->dGetDensity(p1);            /* densita' all'inizio del tubo */
   density2 = HF->dGetDensity(p2);            /* densita' alla fine del nodo */
   density0 = HF->dGetDensity(p0);            /* densita' a meta' tubo */
      
   /* viscosity = HF->dGetViscosity(); */
   densityDPres1 = HF->dGetDensityDPres(.75*p1+.25*p2);
   densityDPres2 = HF->dGetDensityDPres(.25*p1+.75*p2);
   
   viscosity = HF->dGetViscosity(p0);
   
   doublereal qq1 = .75*q1+.25*q2;
   doublereal qq2 = .25*q1+.75*q2;
   
   doublereal dd1 = .5*(density1+density0);
   doublereal dd2 = .5*(density0+density2);
   
   doublereal q12 = .5*(q1+q2);

   /* determines whether the flow will be turbulent or laminar */
   Re = fabs(q12)/area*diameter/viscosity;

   doublereal dLoss1 = 0.;
   doublereal dLoss2 = 0.;

   if ((Re < HF->dGetRe(HydraulicFluid::LOWER)) || ((turbulent == 0) 
			   && (Re < 1.25*HF->dGetRe(HydraulicFluid::LOWER)))) {
      /* laminar, or ascending transition (continuation) */
      doublereal dk = dKlam*viscosity;
      
      dLoss1 = dk*qq1/dd1;
      dLoss2 = dk*qq2/dd2;
   } else if ((Re > HF->dGetRe(HydraulicFluid::UPPER)) || ((turbulent == 1) 
			   && (Re > .775*HF->dGetRe(HydraulicFluid::UPPER)))) {
      /* turbulent, or descending transition (continuation) */
      doublereal dk = dKtrb*pow(viscosity, .25);
      
      dLoss1 = dk*qq1*pow(fabs(qq1), .75)/dd1;
      dLoss2 = dk*qq2*pow(fabs(qq2), .75)/dd2;
   } else {
      /* transition */
      doublereal dva = diameter/(viscosity*area);
      doublereal dva2 = dva*dva;
      doublereal dva3 = dva2*dva;
      
      doublereal a, b, c, d;
      
      if (turbulent == 0) {
	 /* ascending */
	 a = 7.e-13*dva3;
	 b = 8.e-10*dva2;
	 c = -1.8e-5*dva;
	 d = .0542;
      } else /* if (turbulent == 1) */ {
	 /* descending */
	 a = 9.e-13*dva3;
	 b = 2.e-9*dva2;
	 c = -1.8e-5*dva;
	 d = .0528;
      }
	 
      doublereal fqq1 = fabs(qq1);
      doublereal fqq2 = fabs(qq2);
      
      doublereal fax1 = qq1*(((a*fqq1+b)*fqq1+c)*fqq1+d)*fqq1;
      doublereal fax2 = qq2*(((a*fqq2+b)*fqq2+c)*fqq2+d)*fqq2;
      
      dLoss1 = dKtra*fax1/dd1;
      dLoss2 = dKtra*fax2/dd2;
   }

   /* mass conservation */
   doublereal dr = area*length/8.;
   WorkVec.PutItem(1, iNode1RowIndex , q12+dr*densityDPres1*(3.*p1p+p2p));
   WorkVec.PutItem(2, iNode2RowIndex , -q12+dr*densityDPres2*(p1p+3.*p2p));
   
   /* momentum balance */
   doublereal dp = .5*area*(p2-p1);
   doublereal dq12 = q12*q12/(density0*area);
   dr = length/8.;
   WorkVec.PutItem(3, iFirstIndex+1,
		    -dr*(3.*q1p+q2p)-dq12+q1*q1/(density1*area)-dp-dLoss1);
   WorkVec.PutItem(4, iFirstIndex+2,
		    -dr*(q1p+3.*q2p)-q2*q2/(density2*area)+dq12-dp-dLoss2);
   
   /* differential pressure definition */
   WorkVec.PutItem(5, iFirstIndex+3, pn1-p1);
   WorkVec.PutItem(6, iFirstIndex+4, pn2-p2);
      
   return WorkVec;
}

void
DynamicPipe::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
   if (Re < HF->dGetRe(HydraulicFluid::LOWER)) {
      turbulent = 0;
   } else if (Re > HF->dGetRe(HydraulicFluid::UPPER)) {
      turbulent = 1;
   }
   
#ifdef HYDR_DEVEL
   DEBUGCOUT("DynamicPipe(" << GetLabel() << "): turbulent mode =  " 
		   << turbulent << std::endl);
#endif /* HYDR_DEVEL */
}

void DynamicPipe::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) {
      std::ostream& out = OH.Hydraulic();
      out 
	<< std::setw(8) << GetLabel()	/*  1 */
	<< " " << p1			/*  2 */
	<< " " << p2			/*  3 */
	<< " " << p1p			/*  4 */
	<< " " << p2p			/*  5 */
	<< " " << q1			/*  6 */
	<< " " << q2			/*  7 */
	<< " " << q1p			/*  8 */
	<< " " << q2p			/*  9 */
	<< " " << density1		/* 10 */
	<< " " << density0		/* 11 */
	<< " " << density2		/* 12 */
	<< " " << densityDPres1		/* 13 */
	<< " " << densityDPres2		/* 14 */
	<< " " << Re			/* 15 */
	<< " " << turbulent		/* 16 */
	<< std::endl;
   }
}

void 
DynamicPipe::SetValue(DataManager *pDM,
		VectorHandler& X , VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
   integer i = iGetFirstIndex();
   
   doublereal p1 = pNode1->dGetX();
   doublereal p2 = pNode2->dGetX();
   
   X.PutCoef(i+1, p1);
   X.PutCoef(i+2, p2);
   X.PutCoef(i+3, q0);
   X.PutCoef(i+4, -q0);
   
   XP.PutCoef(i+1, 0.);
   XP.PutCoef(i+2, 0.);
   XP.PutCoef(i+3, 0.);
   XP.PutCoef(i+4, 0.);
}

/* DynamicPipe - end */
