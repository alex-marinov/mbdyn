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

/* 
 * Copyright 1999-2000 Lamberto Puggelli <puggelli@tiscalinet.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

extern "C" {
#include <float.h>
}

#include <valve.h>

/* Control_valve - begin */

Control_valve::Control_valve(unsigned int uL, const DofOwner* pDO,
			     HydraulicFluid* hf,
			     const PressureNode* p1, const PressureNode* p2, 
			     const PressureNode* p3, const PressureNode* p4,  
			     doublereal A_max, doublereal Loss_A, const DriveCaller* pDC,
			     flag fOut) 
: Elem(uL, ElemType::HYDRAULIC, fOut),
HydraulicElem(uL, pDO, hf, fOut),
DriveOwner(pDC),
pNode1(p1), pNode2(p2), pNode3(p3), pNode4(p4),
area_max(A_max), loss_area(Loss_A)
{
   ASSERT(pNode1 != NULL);
   ASSERT(pNode1->GetNodeType() == NodeType::HYDRAULIC);
   ASSERT(pNode2 != NULL);
   ASSERT(pNode2->GetNodeType() == NodeType::HYDRAULIC);
   ASSERT(pNode3 != NULL);
   ASSERT(pNode3->GetNodeType() == NodeType::HYDRAULIC);
   ASSERT(pNode4 != NULL);
   ASSERT(pNode4->GetNodeType() == NodeType::HYDRAULIC);
   ASSERT(A_max > DBL_EPSILON);
   ASSERT(loss_area >= 0.);
   
   Cd = .611; /* coefficiente di perdita */
 //   W = .005;  /* larghezza del condotto (m): A=x*W 0.005; */
}

Control_valve::~Control_valve(void)
{
   NO_OP;
}
   
/* Tipo di elemento idraulico (usato solo per debug ecc.) */
HydraulicType::Type Control_valve::GetHydraulicType(void) const 
{
   return HydraulicType::CONTROL_VALVE;
}

/* Contributo al file di restart */
ostream& Control_valve::Restart(ostream& out) const
{
   return out << "Control_valve not implemented yet!" << endl;
}
   
unsigned int Control_valve::iGetNumDof(void) const 
{
   return 0;
}
   
DofOrder::Order Control_valve::SetDof(unsigned int i) const 
{
   cerr << "Control valve has no dofs!" << endl;
   THROW(ErrGeneric());
}


void 
Control_valve::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const 
{
   *piNumRows = 4; 
   *piNumCols = 4; 
}


VariableSubMatrixHandler& 
Control_valve::AssJac(VariableSubMatrixHandler& WorkMat,
		      doublereal dCoef,
		      const VectorHandler& XCurr, 
		      const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering Control_valve::AssJac()" << endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.Resize(4, 4);
   
   integer iNode1RowIndex = pNode1->iGetFirstRowIndex()+1;
   integer iNode2RowIndex = pNode2->iGetFirstRowIndex()+1;
   integer iNode3RowIndex = pNode3->iGetFirstRowIndex()+1;
   integer iNode4RowIndex = pNode4->iGetFirstRowIndex()+1;
   integer iNode1ColIndex = pNode1->iGetFirstColIndex()+1;
   integer iNode2ColIndex = pNode2->iGetFirstColIndex()+1;
   integer iNode3ColIndex = pNode3->iGetFirstColIndex()+1;
   integer iNode4ColIndex = pNode4->iGetFirstColIndex()+1;
  
   WM.fPutRowIndex(1, iNode1RowIndex);
   WM.fPutRowIndex(2, iNode2RowIndex);
   WM.fPutRowIndex(3, iNode3RowIndex);
   WM.fPutRowIndex(4, iNode4RowIndex);
   WM.fPutColIndex(1, iNode1ColIndex);
   WM.fPutColIndex(2, iNode2ColIndex);
   WM.fPutColIndex(3, iNode3ColIndex);
   WM.fPutColIndex(4, iNode4ColIndex);   
   
   doublereal p1 = pNode1->dGetX();
   doublereal p2 = pNode2->dGetX();
   doublereal p3 = pNode3->dGetX();
   doublereal p4 = pNode4->dGetX();
   doublereal density = HF->dGetDensity();
 
   // doublereal A1min = 1.e-8;         /* valore minimo di A1 (m^2)1.e-8 */
   // doublereal A2min = 5.e-9;         /* valore minimo di A2 (m^2)5.e-9 */
   // doublereal A3min = 1.e-8;         /* valore minimo di A3 (m^2)1.e-8 */
   // doublereal A4min = 5.e-9;         /* valore minimo di A4 (m^2)5.e-9 */
  
   doublereal A1min = 2.*area_max*loss_area;  // valore minimo di A1 (e' il doppio di A2min)
   doublereal A2min = area_max*loss_area;     // valore minimo di A2
   doublereal A3min = 2.*area_max*loss_area;  // valore minimo di A3 (e' il doppio di A4min)
   doublereal A4min = area_max*loss_area;     // valore minimo di A4;
  
					 
   doublereal jumpPres12 = fabs(p1-p2); /* salto di pressione nodo1 & nodo3 */
   doublereal jumpPres13 = fabs(p1-p3); /* salto di pressione nodo1 & nodo4 */
   doublereal jumpPres24 = fabs(p2-p4); /* salto di pressione nodo2 & nodo3 */
   doublereal jumpPres34 = fabs(p3-p4); /* salto di pressione nodo2 & nodo4 */

   
   if (jumpPres12 == 0.) {
      jumpPres12 = DBL_EPSILON;
   }
   if (jumpPres13 == 0.) {
      jumpPres13 = DBL_EPSILON;
   }
   if (jumpPres24 == 0.) {
      jumpPres24 = DBL_EPSILON;
   }
   if (jumpPres34 == 0.) {
      jumpPres34 = DBL_EPSILON;
   }
   
   doublereal primo = 0.;
   if (jumpPres12 != 0.) {
      primo = Cd*A1/sqrt(2*jumpPres12/density);
   }

   doublereal secondo = 0.;
   if (jumpPres13 != 0.) {
      secondo = Cd*A2/sqrt(2*jumpPres13/density);
   }

   doublereal terzo = 0.;
   if (jumpPres24 != 0.) {
      terzo = Cd*A4/sqrt(2*jumpPres24/density);
   }
   
   doublereal quarto = 0.;
   if (jumpPres34 != 0.) {
      quarto = Cd*A3/sqrt(2*jumpPres34/density);
   }

   doublereal Jac11 = -primo-secondo;
   doublereal Jac12 = primo;
   doublereal Jac13 = secondo;
   doublereal Jac14 = 0.;
   doublereal Jac21 = primo;
   doublereal Jac22 = -primo-terzo;
   doublereal Jac23 = 0.;
   doublereal Jac24 = terzo;
   doublereal Jac31 = secondo;
   doublereal Jac32 = 0.;
   doublereal Jac33 = -secondo-quarto;
   doublereal Jac34 = quarto;
   doublereal Jac41 = 0.;
   doublereal Jac42 = terzo;
   doublereal Jac43 = quarto;
   doublereal Jac44 = -terzo-quarto;

   WM.fPutCoef(1, 1, Jac11);
   WM.fPutCoef(1, 2, Jac12);
   WM.fPutCoef(1, 3, Jac13);
   WM.fPutCoef(1, 4, Jac14);
   WM.fPutCoef(2, 1, Jac21);
   WM.fPutCoef(2, 2, Jac22);
   WM.fPutCoef(2, 3, Jac23);
   WM.fPutCoef(2, 4, Jac24);
   WM.fPutCoef(3, 1, Jac31);
   WM.fPutCoef(3, 2, Jac32);
   WM.fPutCoef(3, 3, Jac33);
   WM.fPutCoef(3, 4, Jac24);
   WM.fPutCoef(4, 1, Jac41);
   WM.fPutCoef(4, 2, Jac42);
   WM.fPutCoef(4, 3, Jac43);
   WM.fPutCoef(4, 4, Jac44);

   return WorkMat;
}


SubVectorHandler& 
Control_valve::AssRes(SubVectorHandler& WorkVec,
		      doublereal dCoef,
		      const VectorHandler& XCurr, 
		      const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering Control_valve::AssRes()" << endl);
   
   WorkVec.Resize(4);
   
   integer iNode1RowIndex = pNode1->iGetFirstRowIndex()+1;
   integer iNode2RowIndex = pNode2->iGetFirstRowIndex()+1;
   integer iNode3RowIndex = pNode3->iGetFirstRowIndex()+1;
   integer iNode4RowIndex = pNode4->iGetFirstRowIndex()+1;
  
   doublereal p1 = pNode1->dGetX();
   doublereal p2 = pNode2->dGetX();
   doublereal p3 = pNode3->dGetX();
   doublereal p4 = pNode4->dGetX();
   doublereal density = HF->dGetDensity();
 
   Stato = pGetDriveCaller()->dGet();

   doublereal jumpPres12 = fabs(p1-p2); /* salto di pressione nodo1 & nodo2 */
   doublereal jumpPres13 = fabs(p1-p3); /* salto di pressione nodo1 & nodo3 */
   doublereal jumpPres24 = fabs(p2-p4); /* salto di pressione nodo2 & nodo4 */
   doublereal jumpPres34 = fabs(p3-p4); /* salto di pressione nodo3 & nodo4 */

   
    
   doublereal A1min = 2.*area_max*loss_area;  // valore minimo di A1 (e' il doppio di A2min)
   doublereal A2min = area_max*loss_area;     // valore minimo di A2
   doublereal A3min = 2.*area_max*loss_area;  // valore minimo di A3 (e' il doppio di A4min)
   doublereal A4min = area_max*loss_area;     // valore minimo di A4;
    
   					 
   if (Stato > 1.) {
      Stato = 1.;
   } else if (Stato < -1.) {
      Stato = -1.;
   }

   if (Stato > 0.) { 
      A1 = Stato*area_max+A1min;
      A2 = A2min;
      A3 = Stato*area_max+A3min;
      A4 = A4min;
   } else {
      A1 = A1min;
      A2 = -Stato*area_max+A2min;
      A3 = A3min;
      A4 = -Stato*area_max+A4min;
   }

   doublereal Q12 = density*Cd*A1*copysign(sqrt(2.*jumpPres12/density), p1-p2);
   doublereal Q13 = density*Cd*A2*copysign(sqrt(2.*jumpPres13/density), p1-p3);
   doublereal Q24 = density*Cd*A4*copysign(sqrt(2.*jumpPres24/density), p2-p4);
   doublereal Q34 = density*Cd*A3*copysign(sqrt(2.*jumpPres34/density), p3-p4);

   doublereal Res_1 = Q12+Q13;
   doublereal Res_2 = -Q12+Q24;
   doublereal Res_3 = -Q13+Q34;
   doublereal Res_4 = -Q34-Q24;
        
   flow1 = -Res_1;
   flow2 = -Res_2;
   flow3 = -Res_3;
   flow4 = -Res_4;

#ifdef HYDR_DEVEL
   DEBUGCOUT("Stato:   " << Stato << endl);
   DEBUGCOUT("A1:      " << A1 << endl);
   DEBUGCOUT("A2:      " << A2  << endl);
   DEBUGCOUT("A3:      " << A3 << endl);
   DEBUGCOUT("A4:      " << A4 << endl);
   DEBUGCOUT("p1:      " << p1 << endl);
   DEBUGCOUT("p2:      " << p2 << endl);
   DEBUGCOUT("p3:      " << p3 << endl);
   DEBUGCOUT("p4:      " << p4 << endl);
   DEBUGCOUT("Cd:      " << Cd << endl);
   DEBUGCOUT("density: " << density << endl);
   DEBUGCOUT("Area_max:             " << area_max << endl);
   DEBUGCOUT("Loss_area in %:       " << loss_area << endl);
   DEBUGCOUT("Q12:     " << Q12 << endl);
   DEBUGCOUT("Q13:     " << Q13 << endl);
   DEBUGCOUT("Q24:     " << Q24 << endl);
   DEBUGCOUT("Q34:     " << Q34 << endl);
   DEBUGCOUT("PORTATE AI VARI NODI (positive se entranti)"<< endl);
   DEBUGCOUT("-Res_1 (portata nodo1): " << -Res_1 << endl); 
   DEBUGCOUT("-Res_2 (portata nodo2): " << -Res_2 << endl); 
   DEBUGCOUT("-Res_3 (portata nodo3): " << -Res_3 << endl); 
   DEBUGCOUT("-Res_4 (portata nodo4): " << -Res_4 << endl); 
#endif
   
   WorkVec.fPutItem(1, iNode1RowIndex, Res_1);
   WorkVec.fPutItem(2, iNode2RowIndex, Res_2);         
   WorkVec.fPutItem(3, iNode3RowIndex, Res_3);
   WorkVec.fPutItem(4, iNode4RowIndex, Res_4);    

   return WorkVec;
}
  
void Control_valve::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) { 
      ostream& out = OH.Hydraulic();
      out << setw(8) << GetLabel()
	<< " " << Stato 
	<< " " << flow1 << " " << flow2 
	<< " " << flow3 << " " << flow4 << endl;
   }   
}

void 
Control_valve::SetValue(VectorHandler& /* X */ , VectorHandler& /* XP */ ) const
{
   NO_OP;
}

/* Control_valve - end */


/* Dynamic_control_valve - begin */

Dynamic_control_valve::Dynamic_control_valve(unsigned int uL, const DofOwner* pDO, 
					     HydraulicFluid* hf,
					     const PressureNode* p1, 
					     const PressureNode* p2, 
					     const PressureNode* p3, 
					     const PressureNode* p4,
					     const DriveCaller* pDC,
					     doublereal s0, doublereal s_mx, doublereal W, doublereal Loss_A, 
					     doublereal Valve_d, doublereal Valve_rho,
					     doublereal cs, doublereal cv, 
					     doublereal ca, flag fOut) 
: Elem(uL, ElemType::HYDRAULIC, fOut),
HydraulicElem(uL, pDO, hf, fOut), DriveOwner(pDC),
pNode1(p1), pNode2(p2), pNode3(p3), pNode4(p4),
start(s0), c_spost(cs), c_vel(cv), c_acc(ca), 
width(W), loss_area(Loss_A), 
valve_diameter(Valve_d), valve_density(Valve_rho), s_max(s_mx)
{
   ASSERT(pNode1 != NULL);
   ASSERT(pNode1->GetNodeType() == NodeType::HYDRAULIC);
   ASSERT(pNode2 != NULL);
   ASSERT(pNode2->GetNodeType() == NodeType::HYDRAULIC);
   ASSERT(pNode3 != NULL);
   ASSERT(pNode3->GetNodeType() == NodeType::HYDRAULIC);
   ASSERT(pNode4 != NULL);
   ASSERT(pNode4->GetNodeType() == NodeType::HYDRAULIC);
   ASSERT(s0 >= 0.);
   ASSERT(Valve_rho > DBL_EPSILON);
   ASSERT(Valve_d > DBL_EPSILON);
   ASSERT(W > DBL_EPSILON);
   ASSERT(Loss_A >= 0.);
   ASSERT(s_mx >= 0.);
    
   Cd = .611;                /* coefficiente di perdita */
   Mass = 1.175*valve_diameter*valve_diameter*valve_diameter*valve_density*M_PI; // massa della valvola 
   
   // W = .005;                 /* larghezza del condotto (m): A=x*W 0.005; */
   // diameter = .01;           /* diametro della valvola (m) */
   // densityvalve = 7900.;     /* densita' del corpo della valvola (acciaio) (kg/m^3) */
   // Mass = 1.175*diameter*diameter*diameter*densityvalve*M_PI; /* massa della valvola (Kg) */
   // s_max = 2.;               /* corsa massima della valvola (m) */
}


Dynamic_control_valve::~Dynamic_control_valve(void)
{
   NO_OP;
}
   
/* Tipo di elemento idraulico (usato solo per debug ecc.) */
HydraulicType::Type Dynamic_control_valve::GetHydraulicType(void) const 
{
   return HydraulicType::DYNAMIC_CONTROL_VALVE;
}

/* Contributo al file di restart */
ostream& Dynamic_control_valve::Restart(ostream& out) const
{
   return out << "Dynamic_control_valve not implemented yet!" << endl;
}
   
unsigned int Dynamic_control_valve::iGetNumDof(void) const 
{
   return 2;
}
   
DofOrder::Order Dynamic_control_valve::SetDof(unsigned int i) const 
{
   ASSERT(i >= 0 && i <= 1);
   return DofOrder::DIFFERENTIAL;
}


void 
Dynamic_control_valve::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const 
{
   *piNumRows = 6; 
   *piNumCols = 6; 
}
      
VariableSubMatrixHandler& 
Dynamic_control_valve::AssJac(VariableSubMatrixHandler& WorkMat,
		  doublereal dCoef,
		  const VectorHandler& XCurr, 
		  const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering Control_valve::AssJac()" << endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.ResizeInit(6, 6, 0.);
   
   integer iNode1RowIndex = pNode1->iGetFirstRowIndex()+1;
   integer iNode2RowIndex = pNode2->iGetFirstRowIndex()+1;
   integer iNode3RowIndex = pNode3->iGetFirstRowIndex()+1;
   integer iNode4RowIndex = pNode4->iGetFirstRowIndex()+1;
   integer iNode1ColIndex = pNode1->iGetFirstColIndex()+1;
   integer iNode2ColIndex = pNode2->iGetFirstColIndex()+1;
   integer iNode3ColIndex = pNode3->iGetFirstColIndex()+1;
   integer iNode4ColIndex = pNode4->iGetFirstColIndex()+1;
   integer iFirstIndex = iGetFirstIndex();
   
   WM.fPutRowIndex(1, iNode1RowIndex);
   WM.fPutRowIndex(2, iNode2RowIndex);
   WM.fPutRowIndex(3, iNode3RowIndex);
   WM.fPutRowIndex(4, iNode4RowIndex);
   WM.fPutColIndex(1, iNode1ColIndex);
   WM.fPutColIndex(2, iNode2ColIndex);
   WM.fPutColIndex(3, iNode3ColIndex);
   WM.fPutColIndex(4, iNode4ColIndex);
   
   WM.fPutRowIndex(5, iFirstIndex+1);
   WM.fPutColIndex(5, iFirstIndex+1);
   WM.fPutRowIndex(6, iFirstIndex+2);
   WM.fPutColIndex(6, iFirstIndex+2);
   
   doublereal p1 = pNode1->dGetX();
   doublereal p2 = pNode2->dGetX();
   doublereal p3 = pNode3->dGetX();
   doublereal p4 = pNode4->dGetX();
   doublereal density = HF->dGetDensity();

   s = XCurr.dGetCoef(iFirstIndex+1); /* spostamento */
   v = XCurr.dGetCoef(iFirstIndex+2); /* velocita' */
   
   doublereal jumpPres12 = fabs(p1-p2); /* salto di pressione nodo1 & nodo2 */
   doublereal jumpPres13 = fabs(p1-p3); /* salto di pressione nodo1 & nodo3 */
   doublereal jumpPres24 = fabs(p2-p4); /* salto di pressione nodo2 & nodo4 */
   doublereal jumpPres34 = fabs(p3-p4); /* salto di pressione nodo3 & nodo4 */
   
    if (jumpPres12 == 0.) {
      jumpPres12 = DBL_EPSILON;
   }
   if (jumpPres13 == 0.) {
      jumpPres13 = DBL_EPSILON;
   }
   if (jumpPres24 == 0.) {
      jumpPres24 = DBL_EPSILON;
   }
   if (jumpPres34 == 0.) {
      jumpPres34 = DBL_EPSILON;
   }
   doublereal Jac15;
   doublereal Jac25;
   doublereal Jac35;
   doublereal Jac45;
    
   doublereal primo   = Cd*A1/sqrt(2.*jumpPres12/density);
   doublereal secondo = Cd*A2/sqrt(2.*jumpPres13/density);
   doublereal terzo   = Cd*A4/sqrt(2.*jumpPres24/density);
   doublereal quarto  = Cd*A3/sqrt(2.*jumpPres34/density);

   doublereal Jac11 = -primo-secondo;
   doublereal Jac12 = primo;
   doublereal Jac13 = secondo;
   doublereal Jac14 = 0.;
   doublereal Jac21 = primo;
   doublereal Jac22 = -primo-terzo;
   doublereal Jac23 = 0.;
   doublereal Jac24 = terzo;
   doublereal Jac31 = secondo;
   doublereal Jac32 = 0.;
   doublereal Jac33 = -secondo-quarto;
   doublereal Jac34 = quarto;
   doublereal Jac41 = 0.;
   doublereal Jac42 = terzo;
   doublereal Jac43 = quarto;
   doublereal Jac44 = -terzo-quarto;
      
   doublereal costante = density*Cd*width*valve_diameter;
   /* 
   doublereal Jac51  = -.2*costante*sp/sqrt(density*deltaP)-.43*width*s; 
   doublereal Jac52  = -Jac51;
   doublereal Jac53  = Jac51;
   doublereal Jac54  = -Jac51;
   doublereal Jac55  = -.4*valve_diameter*Cd*width*sqrt(density*deltaP)-dCoef*.43*width*deltaP-dCoef*c1-c2-dCoef*cf1-cf2;
   doublereal Jac56  = -Mass-c3-cf3;
  */
   
   doublereal Jac51  = -.2*costante*sp/sqrt(density*deltaP)-.43*width*s; 
   doublereal Jac52  = 0.;
   doublereal Jac53  = 0.;
   doublereal Jac54  = 0.;
   doublereal Jac55  = -.4*valve_diameter*Cd*width*sqrt(density*deltaP)-dCoef*.43*width*deltaP-dCoef*c1-c2-dCoef*cf1-cf2;
   doublereal Jac56  = -Mass-c3-cf3;
   
   if (s > 0.) {  /* collegamento diretto 1->2   3->4 */
      Jac15 = dCoef*density*Cd*width*copysign(sqrt(2.*jumpPres12/density), p1-p2);
      Jac25 = -Jac15;
      Jac35 = dCoef*density*Cd*width*copysign(sqrt(2.*jumpPres34/density), p3-p4);
      Jac45 = -Jac35;
   } else {       /* collegamento inverso 1->3     2->4 */
      Jac15 = -dCoef*density*Cd*width*copysign(sqrt(2.*jumpPres13/density), p1-p3);
      Jac25 = -dCoef*density*Cd*width*copysign(sqrt(2.*jumpPres24/density), p2-p4);
      Jac35 = -Jac15;
      Jac45 = -Jac25;
   }
   
   doublereal Jac65  = -1;
   doublereal Jac66  = dCoef;  

   WM.fPutCoef(1, 1, Jac11);
   WM.fPutCoef(1, 2, Jac12);
   WM.fPutCoef(1, 3, Jac13);
   WM.fPutCoef(1, 4, Jac14);
   WM.fPutCoef(1, 5, Jac15);
   WM.fPutCoef(2, 1, Jac21);
   WM.fPutCoef(2, 2, Jac22);
   WM.fPutCoef(2, 3, Jac23);
   WM.fPutCoef(2, 4, Jac24);
   WM.fPutCoef(2, 5, Jac25);
   WM.fPutCoef(3, 1, Jac31);
   WM.fPutCoef(3, 2, Jac32);
   WM.fPutCoef(3, 3, Jac33);
   WM.fPutCoef(3, 4, Jac34);
   WM.fPutCoef(3, 5, Jac35);
   WM.fPutCoef(4, 1, Jac41);
   WM.fPutCoef(4, 2, Jac42);
   WM.fPutCoef(4, 3, Jac43);
   WM.fPutCoef(4, 4, Jac44);
   WM.fPutCoef(4, 5, Jac45);
   WM.fPutCoef(5, 1, Jac51);
   WM.fPutCoef(5, 2, Jac52);
   WM.fPutCoef(5, 3, Jac53);
   WM.fPutCoef(5, 5, Jac55);
   WM.fPutCoef(5, 6, Jac56);
   WM.fPutCoef(6, 5, Jac65);
   WM.fPutCoef(6, 6, Jac66);
 
   return WorkMat;
}

SubVectorHandler& 
Dynamic_control_valve::AssRes(SubVectorHandler& WorkVec,
			  doublereal dCoef,
			  const VectorHandler& XCurr, 
			  const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering Dynamic_control_valve::AssRes()" << endl);
   
   WorkVec.Resize(6);
   
   integer iNode1RowIndex = pNode1->iGetFirstRowIndex()+1;
   integer iNode2RowIndex = pNode2->iGetFirstRowIndex()+1;
   integer iNode3RowIndex = pNode3->iGetFirstRowIndex()+1;
   integer iNode4RowIndex = pNode4->iGetFirstRowIndex()+1;
   integer iFirstIndex = iGetFirstIndex();
   
   doublereal p1 = pNode1->dGetX();
   doublereal p2 = pNode2->dGetX();
   doublereal p3 = pNode3->dGetX();
   doublereal p4 = pNode4->dGetX();
   doublereal density = HF->dGetDensity();

   Force = pGetDriveCaller()->dGet();
   
   s = XCurr.dGetCoef(iFirstIndex+1);       /* spostamento */
   v = XCurr.dGetCoef(iFirstIndex+2);       /* velocita' */
   sp = XPrimeCurr.dGetCoef(iFirstIndex+1); /* velocita' */
   vp = XPrimeCurr.dGetCoef(iFirstIndex+2); /* accelerazione */
   
   doublereal jumpPres12 = fabs(p1-p2); /* salto di pressione nodo1 & nodo2 */
   doublereal jumpPres13 = fabs(p1-p3); /* salto di pressione nodo1 & nodo3 */
   doublereal jumpPres24 = fabs(p2-p4); /* salto di pressione nodo2 & nodo4 */
   doublereal jumpPres34 = fabs(p3-p4); /* salto di pressione nodo3 & nodo4 */

   /* qui decido se sono a fondo od ad inizio corsa e di conseguenza calcolo
    * i giusti coefficienti */

  
   doublereal area_max = width*s_max;         // Area massima
   doublereal area_min = area_max*loss_area;  // Area minimo 
   doublereal deltaA = s*width;

   if (s > 0.) { 
      A1 = area_min+deltaA;
      A2 = area_min;
      A3 = area_min+deltaA;
      A4 = area_min;
   } else { /* ho deltaA negativo */
      A1 = area_min;
      A2 = area_min-deltaA;
      A3 = area_min;
      A4 = area_min-deltaA;
   }

    if (s <= -s_max) { 
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

   } else {
       c1 = 0.;
       c2 = 0.;
       c3 = 0.;
      }
   
   if (s >= s_max) { 
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
      
   } else {
      cf1 = 0.;
      cf2 = 0.;
      cf3 = 0.;
   }
      
   
   doublereal Q12 = density*Cd*A1*copysign(sqrt(2.*jumpPres12/density), p1-p2);
   doublereal Q13 = density*Cd*A2*copysign(sqrt(2.*jumpPres13/density), p1-p3);
   doublereal Q24 = density*Cd*A4*copysign(sqrt(2.*jumpPres24/density), p2-p4);
   doublereal Q34 = density*Cd*A3*copysign(sqrt(2.*jumpPres34/density), p3-p4);

   doublereal Res_1 = Q12+Q13;
   doublereal Res_2 = -Q12+Q24;
   doublereal Res_3 = -Q13+Q34;
   doublereal Res_4 = -Q34-Q24;
   
   deltaP =p1; 
 
   if (deltaP == 0.) {
      deltaP = 10.; /* evito di dividere per zero nello jacobiano */
   }
  
   doublereal C = .4*valve_diameter*Cd*width*sqrt(density*deltaP);
   doublereal K = .43*width*deltaP;
   doublereal Res_5 = -Force+Mass*vp+C*sp+K*s+c1*(s+s_max)+c2*sp+c3*vp+cf1*(s-s_max)+cf2*sp+cf3*vp;
   doublereal Res_6 = sp-v;

   flow1 = -Res_1;
   flow2 = -Res_2;
   flow3 = -Res_3;
   flow4 = -Res_4;

#ifdef HYDR_DEVEL
   DEBUGCOUT("Force:         " << Force << endl);
   DEBUGCOUT("X equilibrio:  " << Force/K << endl);
   DEBUGCOUT("A1:            " << A1 << endl);
   DEBUGCOUT("A2:            " << A2  << endl);
   DEBUGCOUT("A3:            " << A3 << endl);
   DEBUGCOUT("A4:            " << A4 << endl);
   DEBUGCOUT("p1:            " << p1 << endl);
   DEBUGCOUT("p2:            " << p2 << endl);
   DEBUGCOUT("p3:            " << p3 << endl);
   DEBUGCOUT("p4:            " << p4 << endl);
   DEBUGCOUT("Cd:            " << Cd << endl);
   DEBUGCOUT("s_max :        " << s_max << endl);
   DEBUGCOUT("s :            " << s << endl);
   DEBUGCOUT("sp:            " << sp << endl);
   DEBUGCOUT("v :            " << v << endl);
   DEBUGCOUT("vp:            " << vp << endl);
   DEBUGCOUT("Valve_diameter:" << valve_diameter << endl);
   DEBUGCOUT("massa:         " << Mass << endl);
   DEBUGCOUT("smorzatore:    " << C << endl);
   DEBUGCOUT("molla:         " << K << endl);
   DEBUGCOUT("density:       " << density << endl);
   DEBUGCOUT("Valve_density: " << valve_density << endl);
   DEBUGCOUT("Area_max:      " << area_max << endl);
   DEBUGCOUT("Width:         " << width << endl);
   DEBUGCOUT("Loss_area:     " << loss_area << endl);
   DEBUGCOUT("c1:            " << c1 << endl);
   DEBUGCOUT("c2:            " << c2 << endl);
   DEBUGCOUT("c3:            " << c3 << endl);
   DEBUGCOUT("cf1:           " << cf1 << endl);
   DEBUGCOUT("cf2:           " << cf2 << endl);
   DEBUGCOUT("cf3:           " << cf3 << endl);
   DEBUGCOUT("Q12:           " << Q12 << endl);
   DEBUGCOUT("Q13:           " << Q13 << endl);
   DEBUGCOUT("Q24:           " << Q24 << endl);
   DEBUGCOUT("Q34:           " << Q34 << endl);
   DEBUGCOUT("PORTATE AI VARI NODI (positive se entranti)"<< endl);
   DEBUGCOUT("-Res_1 (portata nodo1): " << -Res_1 << endl); 
   DEBUGCOUT("-Res_2 (portata nodo2): " << -Res_2 << endl); 
   DEBUGCOUT("-Res_3 (portata nodo3): " << -Res_3 << endl); 
   DEBUGCOUT("-Res_4 (portata nodo4): " << -Res_4 << endl); 
   DEBUGCOUT("-Res_5  eq.dinamica   : " << -Res_5 << endl); 
   DEBUGCOUT("-Res_6 sp-v           : " << -Res_6 << endl); 
#endif
   
   WorkVec.fPutItem(1, iNode1RowIndex, Res_1);
   WorkVec.fPutItem(2, iNode2RowIndex, Res_2);         
   WorkVec.fPutItem(3, iNode3RowIndex, Res_3);
   WorkVec.fPutItem(4, iNode4RowIndex, Res_4);    
   WorkVec.fPutItem(5, iFirstIndex+1, Res_5);
   WorkVec.fPutItem(6, iFirstIndex+2, Res_6);    
    
   return WorkVec;
}


void Dynamic_control_valve::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) { 
      ostream& out = OH.Hydraulic();
      out << setw(8) << GetLabel()
	<< " " << s << " "  << sp << " "  << vp 
	<< " " << flow1 << " " << flow2 << " " << flow3 << " " << flow4
	<< " " << A1 << " "  << A2 << " "  << A3 << " " << A4 << " " << endl;
   }
}


void 
Dynamic_control_valve::SetValue(VectorHandler& X , VectorHandler& XP ) const 
{
   integer i = iGetFirstIndex();
   
   X.fPutCoef(i+1, start);
   X.fPutCoef(i+2, 0.);
   XP.fPutCoef(i+1, 0.);
   XP.fPutCoef(i+2, 0.);
}
 
/* Dynamic_control_valve - end */

					
/* Pressure_flow_control_valve - begin */

Pressure_flow_control_valve::Pressure_flow_control_valve(unsigned int uL, const DofOwner* pDO, 
				     HydraulicFluid* hf,
				     const PressureNode* p1, 
				     const PressureNode* p2, 
				     const PressureNode* p3, 
				     const PressureNode* p4,
				     const PressureNode* p5, 
				     const PressureNode* p6,
				     const DriveCaller* pDC,
				     doublereal s0, doublereal s_mx, doublereal W, doublereal Loss_A, 
				     doublereal Valve_d, doublereal Valve_rho,
				     doublereal cs, doublereal cv, 
				     doublereal ca, flag fOut) 
: Elem(uL, ElemType::HYDRAULIC, fOut),
HydraulicElem(uL, pDO, hf, fOut), DriveOwner(pDC),
pNode1(p1), pNode2(p2), pNode3(p3), pNode4(p4), pNode5(p5), pNode6(p6),
start(s0), s_max(s_mx), width(W), loss_area(Loss_A), 
valve_diameter(Valve_d), valve_density(Valve_rho),
c_spost(cs), c_vel(cv), c_acc(ca)
{
   ASSERT(pNode1 != NULL);
   ASSERT(pNode1->GetNodeType() == NodeType::HYDRAULIC);
   ASSERT(pNode2 != NULL);
   ASSERT(pNode2->GetNodeType() == NodeType::HYDRAULIC);
   ASSERT(pNode3 != NULL);
   ASSERT(pNode3->GetNodeType() == NodeType::HYDRAULIC);
   ASSERT(pNode4 != NULL);
   ASSERT(pNode4->GetNodeType() == NodeType::HYDRAULIC);
   ASSERT(pNode5 != NULL);
   ASSERT(pNode5->GetNodeType() == NodeType::HYDRAULIC);
   ASSERT(pNode6 != NULL);
   ASSERT(pNode6->GetNodeType() == NodeType::HYDRAULIC);
   ASSERT(s0 >= 0.);
   ASSERT(Valve_rho > DBL_EPSILON);
   ASSERT(Valve_d > DBL_EPSILON);
   ASSERT(W > DBL_EPSILON);
   ASSERT(Loss_A >= 0.);
   ASSERT(s_mx >= 0.);
    
   Cd = .611;                /* coefficiente di perdita */
   Mass = 1.175*valve_diameter*valve_diameter*valve_diameter*valve_density*M_PI; // massa della valvola 
   valve_area=0.25*valve_diameter*valve_diameter*M_PI;
}


Pressure_flow_control_valve::~Pressure_flow_control_valve(void)
{
   NO_OP;
}
   
/* Tipo di elemento idraulico (usato solo per debug ecc.) */
HydraulicType::Type Pressure_flow_control_valve::GetHydraulicType(void) const 
{
   return HydraulicType::PRESSURE_FLOW_CONTROL_VALVE;
}

/* Contributo al file di restart */
ostream& Pressure_flow_control_valve::Restart(ostream& out) const
{
   return out << "Pressure_flow_control_valve not implemented yet!" << endl;
}
   
unsigned int Pressure_flow_control_valve::iGetNumDof(void) const 
{
   return 2;
}
   
DofOrder::Order Pressure_flow_control_valve::SetDof(unsigned int i) const 
{
   ASSERT(i >= 0 && i <= 1);
   return DofOrder::DIFFERENTIAL;
}


void 
Pressure_flow_control_valve::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const 
{
   *piNumRows = 8; 
   *piNumCols = 8; 
}
      
VariableSubMatrixHandler& 
Pressure_flow_control_valve::AssJac(VariableSubMatrixHandler& WorkMat,
		  doublereal dCoef,
		  const VectorHandler& XCurr, 
		  const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering Control_valve::AssJac()" << endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.ResizeInit(8, 8, 0.);
   
   integer iNode1RowIndex = pNode1->iGetFirstRowIndex()+1;
   integer iNode2RowIndex = pNode2->iGetFirstRowIndex()+1;
   integer iNode3RowIndex = pNode3->iGetFirstRowIndex()+1;
   integer iNode4RowIndex = pNode4->iGetFirstRowIndex()+1;
   integer iNode5RowIndex = pNode5->iGetFirstRowIndex()+1;
   integer iNode6RowIndex = pNode6->iGetFirstRowIndex()+1;
   integer iNode1ColIndex = pNode1->iGetFirstColIndex()+1;
   integer iNode2ColIndex = pNode2->iGetFirstColIndex()+1;
   integer iNode3ColIndex = pNode3->iGetFirstColIndex()+1;
   integer iNode4ColIndex = pNode4->iGetFirstColIndex()+1;
   integer iNode5ColIndex = pNode5->iGetFirstColIndex()+1;
   integer iNode6ColIndex = pNode6->iGetFirstColIndex()+1;
   integer iFirstIndex = iGetFirstIndex();
   
   WM.fPutRowIndex(1, iNode1RowIndex);
   WM.fPutRowIndex(2, iNode2RowIndex);
   WM.fPutRowIndex(3, iNode3RowIndex);
   WM.fPutRowIndex(4, iNode4RowIndex);
   WM.fPutRowIndex(5, iNode5RowIndex);
   WM.fPutRowIndex(6, iNode6RowIndex);
   WM.fPutColIndex(1, iNode1ColIndex);
   WM.fPutColIndex(2, iNode2ColIndex);
   WM.fPutColIndex(3, iNode3ColIndex);
   WM.fPutColIndex(4, iNode4ColIndex);
   WM.fPutColIndex(5, iNode5ColIndex);
   WM.fPutColIndex(6, iNode6ColIndex);
 
   WM.fPutRowIndex(7, iFirstIndex+1);
   WM.fPutColIndex(7, iFirstIndex+1);
   WM.fPutRowIndex(8, iFirstIndex+2);
   WM.fPutColIndex(8, iFirstIndex+2);
   
   doublereal p1 = pNode1->dGetX();
   doublereal p2 = pNode2->dGetX();
   doublereal p3 = pNode3->dGetX();
   doublereal p4 = pNode4->dGetX();
   doublereal p5 = pNode5->dGetX();
   doublereal p6 = pNode6->dGetX();
   doublereal density = HF->dGetDensity();

   s = XCurr.dGetCoef(iFirstIndex+1); /* spostamento */
   v = XCurr.dGetCoef(iFirstIndex+2); /* velocita' */
   
   doublereal jumpPres12 = fabs(p1-p2); /* salto di pressione nodo1 & nodo2 */
   doublereal jumpPres13 = fabs(p1-p3); /* salto di pressione nodo1 & nodo3 */
   doublereal jumpPres24 = fabs(p2-p4); /* salto di pressione nodo2 & nodo4 */
   doublereal jumpPres34 = fabs(p3-p4); /* salto di pressione nodo3 & nodo4 */
   
    if (jumpPres12 == 0.) {
      jumpPres12 = DBL_EPSILON;
   }
   if (jumpPres13 == 0.) {
      jumpPres13 = DBL_EPSILON;
   }
   if (jumpPres24 == 0.) {
      jumpPres24 = DBL_EPSILON;
   }
   if (jumpPres34 == 0.) {
      jumpPres34 = DBL_EPSILON;
   }
   doublereal Jac17;
   doublereal Jac27;
   doublereal Jac37;
   doublereal Jac47;
    
   doublereal primo   = Cd*A1/sqrt(2.*jumpPres12/density);
   doublereal secondo = Cd*A2/sqrt(2.*jumpPres13/density);
   doublereal terzo   = Cd*A4/sqrt(2.*jumpPres24/density);
   doublereal quarto  = Cd*A3/sqrt(2.*jumpPres34/density);

   doublereal Jac11 = -primo-secondo;
   doublereal Jac12 = primo;
   doublereal Jac13 = secondo;
   doublereal Jac14 = 0.;
   doublereal Jac21 = primo;
   doublereal Jac22 = -primo-terzo;
   doublereal Jac23 = 0.;
   doublereal Jac24 = terzo;
   doublereal Jac31 = secondo;
   doublereal Jac32 = 0.;
   doublereal Jac33 = -secondo-quarto;
   doublereal Jac34 = quarto;
   doublereal Jac41 = 0.;
   doublereal Jac42 = terzo;
   doublereal Jac43 = quarto;
   doublereal Jac44 = -terzo-quarto;
  
   doublereal Jac57 = -valve_area;
   doublereal Jac67 = +valve_area;
     
   doublereal costante = density*Cd*width*valve_diameter;
   
   doublereal Jac71  = -.2*costante/sqrt(density*deltaP)*sp-.43*width*s; 
   doublereal Jac72  = 0.;
   doublereal Jac73  = 0.;
   doublereal Jac74  = 0.;
   doublereal Jac77  = -.4*valve_diameter*Cd*width*sqrt(density*deltaP)-dCoef*.43*width*deltaP-dCoef*c1-c2-dCoef*cf1-cf2;
   doublereal Jac78  = -Mass-c3-cf3;

   if (s > 0.) {  /* collegamento diretto 1->2   3->4 */
      Jac17 = dCoef*density*Cd*width*copysign(sqrt(2.*jumpPres12/density), p1-p2);
      Jac27 = -Jac17;
      Jac37 = dCoef*density*Cd*width*copysign(sqrt(2.*jumpPres34/density), p3-p4);
      Jac47 = -Jac37;
   } else {       /* collegamento inverso 1->3     2->4 */
      Jac17 = -dCoef*density*Cd*width*copysign(sqrt(2.*jumpPres13/density), p1-p3);
      Jac27 = -dCoef*density*Cd*width*copysign(sqrt(2.*jumpPres24/density), p2-p4);
      Jac37 = -Jac17;
      Jac47 = -Jac27;
   }
   
   doublereal Jac87  = -1;
   doublereal Jac88  = dCoef;  

   WM.fPutCoef(1, 1, Jac11);
   WM.fPutCoef(1, 2, Jac12);
   WM.fPutCoef(1, 3, Jac13);
   WM.fPutCoef(1, 4, Jac14);
   WM.fPutCoef(1, 7, Jac17);
   WM.fPutCoef(2, 1, Jac21);
   WM.fPutCoef(2, 2, Jac22);
   WM.fPutCoef(2, 3, Jac23);
   WM.fPutCoef(2, 4, Jac24);
   WM.fPutCoef(2, 7, Jac27);
   WM.fPutCoef(3, 1, Jac31);
   WM.fPutCoef(3, 2, Jac32);
   WM.fPutCoef(3, 3, Jac33);
   WM.fPutCoef(3, 4, Jac34);
   WM.fPutCoef(3, 7, Jac37);
   WM.fPutCoef(4, 1, Jac41);
   WM.fPutCoef(4, 2, Jac42);
   WM.fPutCoef(4, 3, Jac43);
   WM.fPutCoef(4, 4, Jac44);
   WM.fPutCoef(4, 7, Jac47);
   WM.fPutCoef(5, 7, Jac67);
   WM.fPutCoef(6, 7, Jac67);
   WM.fPutCoef(7, 1, Jac71);
   WM.fPutCoef(7, 2, Jac72);
   WM.fPutCoef(7, 3, Jac73);
   WM.fPutCoef(7, 7, Jac77);
   WM.fPutCoef(7, 8, Jac78);
   WM.fPutCoef(8, 7, Jac87);
   WM.fPutCoef(8, 8, Jac88);
 
   return WorkMat;
}

SubVectorHandler& 
Pressure_flow_control_valve::AssRes(SubVectorHandler& WorkVec,
			  doublereal dCoef,
			  const VectorHandler& XCurr, 
			  const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering Pressure_flow_control_valve::AssRes()" << endl);
   
   WorkVec.Resize(8);
   
   integer iNode1RowIndex = pNode1->iGetFirstRowIndex()+1;
   integer iNode2RowIndex = pNode2->iGetFirstRowIndex()+1;
   integer iNode3RowIndex = pNode3->iGetFirstRowIndex()+1;
   integer iNode4RowIndex = pNode4->iGetFirstRowIndex()+1;
   integer iNode5RowIndex = pNode5->iGetFirstRowIndex()+1;
   integer iNode6RowIndex = pNode6->iGetFirstRowIndex()+1;
   integer iFirstIndex = iGetFirstIndex();
   
   doublereal p1 = pNode1->dGetX();
   doublereal p2 = pNode2->dGetX();
   doublereal p3 = pNode3->dGetX();
   doublereal p4 = pNode4->dGetX();
   doublereal p5 = pNode5->dGetX();
   doublereal p6 = pNode6->dGetX();
   doublereal density = HF->dGetDensity();

   Force = pGetDriveCaller()->dGet();
   
   s = XCurr.dGetCoef(iFirstIndex+1);       /* spostamento */
   v = XCurr.dGetCoef(iFirstIndex+2);       /* velocita' */
   sp = XPrimeCurr.dGetCoef(iFirstIndex+1); /* velocita' */
   vp = XPrimeCurr.dGetCoef(iFirstIndex+2); /* accelerazione */
   
   doublereal jumpPres12 = fabs(p1-p2); /* salto di pressione nodo1 & nodo2 */
   doublereal jumpPres13 = fabs(p1-p3); /* salto di pressione nodo1 & nodo3 */
   doublereal jumpPres24 = fabs(p2-p4); /* salto di pressione nodo2 & nodo4 */
   doublereal jumpPres34 = fabs(p3-p4); /* salto di pressione nodo3 & nodo4 */

   /* qui decido se sono a fondo od ad inizio corsa e di conseguenza calcolo
    * i giusti coefficienti */

   doublereal x = s;  /* usato per calcolare le aree */
   
   doublereal area_max = width*s_max;         // Area massima
   doublereal area_min = area_max*loss_area;  // Area minimo 
   doublereal deltaA = x*width;

   if (x > 0.) { 
      A1 = area_min+deltaA;
      A2 = area_min;
      A3 = area_min+deltaA;
      A4 = area_min;
   } else { /* ho deltaA negativo */
      A1 = area_min;
      A2 = area_min-deltaA;
      A3 = area_min;
      A4 = area_min-deltaA;
   }
   
   doublereal Q12 = density*Cd*A1*copysign(sqrt(2.*jumpPres12/density), p1-p2);
   doublereal Q13 = density*Cd*A2*copysign(sqrt(2.*jumpPres13/density), p1-p3);
   doublereal Q24 = density*Cd*A4*copysign(sqrt(2.*jumpPres24/density), p2-p4);
   doublereal Q34 = density*Cd*A3*copysign(sqrt(2.*jumpPres34/density), p3-p4);

   doublereal Res_1 = Q12+Q13;
   doublereal Res_2 = -Q12+Q24;
   doublereal Res_3 = -Q13+Q34;
   doublereal Res_4 = -Q34-Q24;
   doublereal Res_5 = valve_area*sp;
   doublereal Res_6 = -Res_6;
 
   deltaP = p1;   
   if (deltaP == 0.) {
      deltaP = 10.; /* evito di dividere per zero nello jacobiano */
   }
   doublereal C = .4*valve_diameter*Cd*width*sqrt(density*deltaP);
   doublereal K = .43*width*deltaP;
   doublereal Res_7 = -Force+Mass*vp+C*sp+K*s+c1*s+c2*sp+c3*vp+cf1*(s-s_max)+cf2*sp+cf3*vp;
   doublereal Res_8 = sp-v;

   flow1 = -Res_1;
   flow2 = -Res_2;
   flow3 = -Res_3;
   flow4 = -Res_4;
   flow5 = -Res_5;
   flow6 = -Res_6;

#ifdef HYDR_DEVEL
   DEBUGCOUT("Force:         " << Force << endl);
   DEBUGCOUT("X equilibrio:  " << Force/K << endl);
   DEBUGCOUT("A1:            " << A1 << endl);
   DEBUGCOUT("A2:            " << A2  << endl);
   DEBUGCOUT("A3:            " << A3 << endl);
   DEBUGCOUT("A4:            " << A4 << endl);
   DEBUGCOUT("p1:            " << p1 << endl);
   DEBUGCOUT("p2:            " << p2 << endl);
   DEBUGCOUT("p3:            " << p3 << endl);
   DEBUGCOUT("p4:            " << p4 << endl);
   DEBUGCOUT("Cd:            " << Cd << endl);
   DEBUGCOUT("s_max :        " << s_max << endl);
   DEBUGCOUT("x :            " << x << endl);
   DEBUGCOUT("s :            " << s << endl);
   DEBUGCOUT("sp:            " << sp << endl);
   DEBUGCOUT("v :            " << v << endl);
   DEBUGCOUT("vp:            " << vp << endl);
   DEBUGCOUT("Valve_diameter:" << valve_diameter << endl);
   DEBUGCOUT("massa:         " << Mass << endl);
   DEBUGCOUT("smorzatore:    " << C << endl);
   DEBUGCOUT("molla:         " << K << endl);
   DEBUGCOUT("density:       " << density << endl);
   DEBUGCOUT("Valve_density: " << valve_density << endl);
   DEBUGCOUT("Area_max:      " << area_max << endl);
   DEBUGCOUT("Width:         " << width << endl);
   DEBUGCOUT("Loss_area:     " << loss_area << endl);
   DEBUGCOUT("c1:            " << c1 << endl);
   DEBUGCOUT("c2:            " << c2 << endl);
   DEBUGCOUT("c3:            " << c3 << endl);
   DEBUGCOUT("cf1:           " << cf1 << endl);
   DEBUGCOUT("cf2:           " << cf2 << endl);
   DEBUGCOUT("cf3:           " << cf3 << endl);
   DEBUGCOUT("Q12:           " << Q12 << endl);
   DEBUGCOUT("Q13:           " << Q13 << endl);
   DEBUGCOUT("Q24:           " << Q24 << endl);
   DEBUGCOUT("Q34:           " << Q34 << endl);
   DEBUGCOUT("PORTATE AI VARI NODI (positive se entranti)"<< endl);
   DEBUGCOUT("-Res_1 (portata nodo1): " << -Res_1 << endl); 
   DEBUGCOUT("-Res_2 (portata nodo2): " << -Res_2 << endl); 
   DEBUGCOUT("-Res_3 (portata nodo3): " << -Res_3 << endl); 
   DEBUGCOUT("-Res_4 (portata nodo4): " << -Res_4 << endl); 
   DEBUGCOUT("-Res_5 (portata nodo3): " << -Res_5 << endl); 
   DEBUGCOUT("-Res_6 (portata nodo4): " << -Res_6 << endl); 
   DEBUGCOUT("-Res_6   eq dinamica  : " << -Res_7 << endl); 
   DEBUGCOUT("-Res_7 sp-v           : " << -Res_8 << endl); 
#endif
   
   WorkVec.fPutItem(1, iNode1RowIndex, Res_1);
   WorkVec.fPutItem(2, iNode2RowIndex, Res_2);         
   WorkVec.fPutItem(3, iNode3RowIndex, Res_3);
   WorkVec.fPutItem(4, iNode4RowIndex, Res_4);    
   WorkVec.fPutItem(5, iNode5RowIndex, Res_5);
   WorkVec.fPutItem(6, iNode6RowIndex, Res_6);    
   WorkVec.fPutItem(7, iFirstIndex+1, Res_7);
   WorkVec.fPutItem(8, iFirstIndex+2, Res_8);    
    
   return WorkVec;
}


void Pressure_flow_control_valve::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) { 
      ostream& out = OH.Hydraulic();
      out << setw(8) << GetLabel()
	<< " " << s << " "  << sp << " "  << vp 
	<< " " << flow1 << " " << flow2 << " " << flow3 << " " << flow4
	<< " " << flow5 << " " << flow6
	<< " " << A1 << " "  << A2 << " "  << A3 << " " << A4 << " " << endl;
   }
}


void 
Pressure_flow_control_valve::SetValue(VectorHandler& X , VectorHandler& XP ) const 
{
   integer i = iGetFirstIndex();
   
   X.fPutCoef(i+1, start);
   X.fPutCoef(i+2, 0.);
   XP.fPutCoef(i+1, 0.);
   XP.fPutCoef(i+2, 0.);
}
 
/* Pressure_flow_control_valve - end */


/* Pressure_valve - begin */

Pressure_valve::Pressure_valve(unsigned int uL, const DofOwner* pDO,
			       HydraulicFluid* hf,
			       const PressureNode* p1, const PressureNode* p2, 
			       doublereal A_dia, 
			       doublereal mv, 
			       doublereal A_max, doublereal s_mx, doublereal K, 
			       doublereal F0, doublereal w, 
			       doublereal cs, doublereal cv, doublereal ca,
			       flag fOut)
: Elem(uL, ElemType::HYDRAULIC, fOut),
HydraulicElem(uL, pDO, hf, fOut),
pNode1(p1), pNode2(p2), area_diaf(A_dia), mass(mv),
area_max(A_max),s_max(s_mx), 
Kappa(K), force0(F0), width(w),c_spost(cs), c_vel(cv), c_acc(ca)
{
   ASSERT(pNode1 != NULL);
   ASSERT(pNode1->GetNodeType() == NodeType::HYDRAULIC);
   ASSERT(pNode2 != NULL);
   ASSERT(pNode2->GetNodeType() == NodeType::HYDRAULIC);
   ASSERT(A_dia > DBL_EPSILON);
   ASSERT(mv > DBL_EPSILON); 
   ASSERT(A_max > DBL_EPSILON);
   ASSERT(K >= 0. );
   ASSERT(F0 >= 0.);
   ASSERT(w > DBL_EPSILON); 
   ASSERT(s_mx>=0.);  // corsa massima della valvola

}

Pressure_valve::~Pressure_valve(void)
{
   NO_OP;
}
   
/* Tipo di elemento idraulico (usato solo per debug ecc.) */
HydraulicType::Type Pressure_valve::GetHydraulicType(void) const 
{
   return HydraulicType::PRESSURE_VALVE;
}


/* Contributo al file di restart */
ostream& Pressure_valve::Restart(ostream& out) const
{
   return out << "Pressure_valve not implemented yet!" << endl;
}


unsigned int Pressure_valve::iGetNumDof(void) const 
{
  return 2; 
}


DofOrder::Order Pressure_valve::SetDof(unsigned int i) const 
{   
   ASSERT(i >= 0 && i <= 1);
   return DofOrder::DIFFERENTIAL;
}


void 
Pressure_valve::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const 
{   
   *piNumRows = 4;
   *piNumCols = 4;   
    
}


VariableSubMatrixHandler& 
Pressure_valve::AssJac(VariableSubMatrixHandler& WorkMat,
		       doublereal dCoef,
		       const VectorHandler& XCurr, 
		       const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering Pressure_valve::AssJac()" << endl);
      
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.Resize(4, 4);
   
   integer iNode1RowIndex = pNode1->iGetFirstRowIndex()+1;
   integer iNode2RowIndex = pNode2->iGetFirstRowIndex()+1;
   integer iNode1ColIndex = pNode1->iGetFirstColIndex()+1;
   integer iNode2ColIndex = pNode2->iGetFirstColIndex()+1;
   integer iFirstIndex = iGetFirstIndex();
 
   WM.fPutRowIndex(1, iNode1RowIndex);
   WM.fPutRowIndex(2, iNode2RowIndex);
   WM.fPutColIndex(1, iNode1ColIndex);
   WM.fPutColIndex(2, iNode2ColIndex);

   WM.fPutRowIndex(3, iFirstIndex+1);
   WM.fPutColIndex(3, iFirstIndex+1);
   WM.fPutRowIndex(4, iFirstIndex+2);
   WM.fPutColIndex(4, iFirstIndex+2);

   doublereal p1 = pNode1->dGetX();
   doublereal p2 = pNode2->dGetX(); 
   /* p1 = XCurr.dGetCoef(pNode1->iGetFirstRowIndex()+1); */
   doublereal density = HF->dGetDensity();
  
   s = XCurr.dGetCoef(iFirstIndex+1); /* spostamento */
   v = XCurr.dGetCoef(iFirstIndex+2); /* velocita' */
   
   doublereal Cd = .6;       
   doublereal jumpPres = fabs(p1-p2);
 
   /* evito di dividere per un numero troppo piccolo */
   if (jumpPres < 1.e8*DBL_EPSILON) {
      jumpPres = 1.e8*DBL_EPSILON;
   }

   doublereal Jac11 = 0.;
   doublereal Jac12 = 0.;
   doublereal Jac13 = 0.;
   doublereal Jac14 = 0.;
   doublereal Jac21 = 0.;
   doublereal Jac22 = 0.;
   doublereal Jac23 = 0.;
   doublereal Jac24 = 0.;
   doublereal Jac31 = 0.;
   doublereal Jac32 = 0.;
   doublereal Jac33 = 0.;
   doublereal Jac34 = 0.;
   doublereal Jac41 = 0.;
   doublereal Jac42 = 0.;
   doublereal Jac43 = 0.;
   doublereal Jac44 = 0.;
   
#ifdef HYDR_DEVEL
   DEBUGCOUT("Valore di c1: " << c1 << endl);
   DEBUGCOUT("Valore di c2: " << c2 << endl);
   DEBUGCOUT("Valore di c3: " << c3 << endl);
#endif
   
   Jac11 = -density*width*s*Cd/sqrt(2*jumpPres/density)/density;
   Jac12 = -Jac11;
   Jac13 = -density*width*Cd*copysign(sqrt(2*jumpPres/density), p1-p2)*dCoef;
   Jac14 = 0.;
   
   Jac21 = -Jac11;
   Jac22 = Jac11;
   Jac23 = -Jac13; /* density*width*Cd*copysign(sqrt(2*jumpPres/density),(p1-p2))*dCoef; */
   Jac24 = 0.;
   
   Jac31 = area_max;
   Jac32 = 0.;
   Jac33 = -Kappa*dCoef-c1*dCoef-cf1*dCoef;
   Jac34 = -mass
     -density*area_max*pow(area_max/(area_diaf*Cd), 2)*fabs(v)*dCoef
     -c2*dCoef-c3-cf2*dCoef-cf3;
   
   Jac41 = 0.;
   Jac42 = 0.;
   Jac43 = -1.;
   Jac44 = dCoef;
   
#ifdef HYDR_DEVEL
   DEBUGCOUT("Jac smorzatore " << density*area_max*pow(area_max/(area_diaf*Cd), 2) << endl);
#endif
   
   WM.fPutCoef(1, 1, Jac11);
   WM.fPutCoef(1, 2, Jac12);
   WM.fPutCoef(1, 3, Jac13);
   WM.fPutCoef(1, 4, Jac14);
   WM.fPutCoef(2, 1, Jac21);
   WM.fPutCoef(2, 2, Jac22);
   WM.fPutCoef(2, 3, Jac23);
   WM.fPutCoef(2, 4, Jac24);
   WM.fPutCoef(3, 1, Jac31);
   WM.fPutCoef(3, 2, Jac32);
   WM.fPutCoef(3, 3, Jac33);
   WM.fPutCoef(3, 4, Jac34);
   WM.fPutCoef(4, 1, Jac41);
   WM.fPutCoef(4, 2, Jac42);
   WM.fPutCoef(4, 3, Jac43);
   WM.fPutCoef(4, 4, Jac44);

   return WorkMat;
}

SubVectorHandler& 
Pressure_valve::AssRes(SubVectorHandler& WorkVec,
		       doublereal dCoef,
		       const VectorHandler& XCurr, 
		       const VectorHandler& XPrimeCurr)

{   
   integer iNode1RowIndex = pNode1->iGetFirstRowIndex()+1;
   integer iNode2RowIndex = pNode2->iGetFirstRowIndex()+1;

   doublereal p1 = pNode1->dGetX();
   doublereal p2 = pNode2->dGetX();
   
   WorkVec.Resize(4);

   integer iFirstIndex = iGetFirstIndex();
   
   s = XCurr.dGetCoef(iFirstIndex+1);       /* spostamento */
   v = XCurr.dGetCoef(iFirstIndex+2);       /* velocita' */
   sp = XPrimeCurr.dGetCoef(iFirstIndex+1); /* velocita' */
   vp = XPrimeCurr.dGetCoef(iFirstIndex+2); /* accelerazione */
   doublereal density = HF->dGetDensity();
  
   doublereal Cd = .6;
   doublereal jumpPres = fabs(p1-p2);
   doublereal Res_1 = 0.;
   doublereal Res_2 = 0.;
   doublereal Res_3 = 0.;
   doublereal Res_4 = 0.;
   doublereal x0 = (p1*area_max-force0)/Kappa; /* punto di equilibrio */
   
#ifdef HYDR_DEVEL
   DEBUGCOUT("s:  " << s << endl);
   DEBUGCOUT("sp: " << sp << endl);
   DEBUGCOUT("v:  " << v << endl);
   DEBUGCOUT("vp: " << vp << endl);
#endif
   
   if (s <= 0.) { 
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

      c4 = 0.; /* (force0-p1*area_max); */
   } else {
       c1 = 0.;
       c2 = 0.;
       c3 = 0.;
       c4 = 0.;
   }
   
   if (s >= s_max) { 
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

      cf4 = 0.; /* (force0-p1*area_max); */
   } else {
      cf1 = 0.;
      cf2 = 0.;
      cf3 = 0.;
      cf4 = 0.;
   }

   Res_1 = density*width*s*Cd*copysign(sqrt(2.*jumpPres/density), p1-p2);
   Res_2 = -Res_1; 
   Res_3 = mass*vp-p1*area_max
     +copysign(.5*density*area_max*pow(sp*area_max/(area_diaf*Cd), 2), sp)
       +Kappa*s+force0+c1*s+c2*sp+c3*vp+c4+cf1*(s-s_max)+cf2*sp+cf3*vp+cf4;
   Res_4 = sp-v;
   
  
   flow1 = -Res_1;   /* per l'output */
   flow2 = -Res_2;   /* per l'output */

#ifdef HYDR_DEVEL
   DEBUGCOUT("width:    " << width << endl);
   DEBUGCOUT("Cd:       " << Cd << endl);
   DEBUGCOUT("jumpPres: " << jumpPres  << endl);
   DEBUGCOUT("density:  " << density << endl);
   DEBUGCOUT("p1:       " << p1 << endl);
   DEBUGCOUT("p2:       " << p2 << endl);
   DEBUGCOUT("s:        " << s << endl);
   DEBUGCOUT("sp:       " << sp << endl);
   DEBUGCOUT("v:        " << v << endl);
   DEBUGCOUT("vp:       " << vp << endl);
   DEBUGCOUT("area_max: " << area_max << endl);
   DEBUGCOUT("Kappa:    " << Kappa << endl);
   DEBUGCOUT("force0:   " << force0 << endl);
   DEBUGCOUT("mass:     " << mass  << endl);
   DEBUGCOUT("x0:       " << x0 << endl);
   DEBUGCOUT("s_max:    " << s_max << endl);
   DEBUGCOUT("Valore dello smorzatore: " 
	     << copysign(.5*density*area_max*pow(area_max/(area_diaf*Cd), 2), sp) << endl);
   DEBUGCOUT("c1:       " << c1 << endl);
   DEBUGCOUT("c2:       " << c2 << endl);
   DEBUGCOUT("c3:       " << c3 << endl);
   DEBUGCOUT("c4:       " << c4 << endl);
   DEBUGCOUT("cf1:      " << cf1 << endl);
   DEBUGCOUT("cf2:      " << cf2 << endl);
   DEBUGCOUT("cf3:      " << cf3 << endl);
   DEBUGCOUT("cf4:      " << cf4 << endl);
   DEBUGCOUT("PORTATE AI VARI NODI (positive se entranti)" << endl);
   DEBUGCOUT("-Res_1 (portata nodo1): " << -Res_1 << endl); 
   DEBUGCOUT("-Res_2 (portata nodo2): " << -Res_2 << endl);   
   DEBUGCOUT("Res_3:                  " << Res_3 << endl); 
   DEBUGCOUT("Res_4:                  " << Res_4 << endl); 
#endif
   
   WorkVec.fPutItem(1, iNode1RowIndex, Res_1);
   WorkVec.fPutItem(2, iNode2RowIndex, Res_2);         
   WorkVec.fPutItem(3, iFirstIndex+1, Res_3);
   WorkVec.fPutItem(4, iFirstIndex+2, Res_4);    
     	
   return WorkVec;
}


void Pressure_valve::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) { 
      ostream& out = OH.Hydraulic();
      out << setw(8) << GetLabel()
	<< " " << s << " " << v  << " "<< vp  
	<< " " << flow1  << " " << flow2 << endl;
   }  
}


void Pressure_valve::SetValue(VectorHandler& X, VectorHandler& XP) const 
{
   integer i = iGetFirstIndex();
   
   X.fPutCoef(i+1, 0.);
   X.fPutCoef(i+2, 0.);
   XP.fPutCoef(i+1, 0.);
   XP.fPutCoef(i+2, 0.);
}

/* Pressure_valve - end */
					

				  
 /* Flow_valve - begin */

Flow_valve:: Flow_valve(unsigned int uL, const DofOwner* pDO,
			HydraulicFluid* hf,
			const PressureNode* p1, const PressureNode* p2, 
			const PressureNode* p3, 
			doublereal A_dia, 
			doublereal mv, doublereal A_pipe, doublereal A_max, 
			doublereal K, doublereal F0, doublereal w,
			doublereal s_mx,  
			doublereal cs, doublereal cv, doublereal ca, 
			flag fOut)
: Elem(uL, ElemType::HYDRAULIC, fOut),
HydraulicElem(uL, pDO, hf, fOut),
pNode1(p1), pNode2(p2), pNode3(p3),
area_diaf(A_dia), mass(mv), area_pipe(A_pipe), area_max(A_max), 
Kappa(K), force0(F0), width(w), s_max(s_mx), c_spost(cs), c_vel(cv), c_acc(ca) 
{
   ASSERT(pNode1 != NULL);
   ASSERT(pNode1->GetNodeType() == NodeType::HYDRAULIC);
   ASSERT(pNode2 != NULL);
   ASSERT(pNode2->GetNodeType() == NodeType::HYDRAULIC);
   ASSERT(pNode3 != NULL);
   ASSERT(pNode3->GetNodeType() == NodeType::HYDRAULIC);
   ASSERT(A_dia > DBL_EPSILON);
   ASSERT(mv > DBL_EPSILON); 
   ASSERT(A_pipe > DBL_EPSILON);
   ASSERT(A_max > DBL_EPSILON);
   ASSERT(K > DBL_EPSILON );
   ASSERT(F0 >= 0.);
   ASSERT(w > DBL_EPSILON); 
   ASSERT(s_max >= 0.);
   
   h = .02; /* coefficiente di perdita di carico concentrata tra i nodi 1 e 2 (smorza il moto della valvola) */
}


Flow_valve::~Flow_valve(void)
{
   NO_OP;
}
 

/* Tipo di elemento idraulico (usato solo per debug ecc.) */
HydraulicType::Type  Flow_valve::GetHydraulicType(void) const 
{
   return HydraulicType:: FLOW_VALVE;
}


/* Contributo al file di restart */
ostream&  Flow_valve::Restart(ostream& out) const
{
   return out << " Flow_valve not implemented yet!" << endl;
}
 

unsigned int  Flow_valve::iGetNumDof(void) const 
{
   return 2;
}


DofOrder::Order Flow_valve::SetDof(unsigned int i) const 
{
   ASSERT(i >= 0 && i <= 1);
   return DofOrder::DIFFERENTIAL;
}


void 
Flow_valve::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const 
{
   *piNumRows = 5;
   *piNumCols = 5; 
}
      
VariableSubMatrixHandler& 
Flow_valve::AssJac(VariableSubMatrixHandler& WorkMat,
		   doublereal dCoef,
		   const VectorHandler& XCurr, 
		   const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering  Flow_valve::AssJac()" << endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.Resize(5, 5);
   
   integer iNode1RowIndex = pNode1->iGetFirstRowIndex()+1;
   integer iNode2RowIndex = pNode2->iGetFirstRowIndex()+1;
   integer iNode3RowIndex = pNode3->iGetFirstRowIndex()+1;
   integer iNode1ColIndex = pNode1->iGetFirstColIndex()+1;
   integer iNode2ColIndex = pNode2->iGetFirstColIndex()+1;
   integer iNode3ColIndex = pNode3->iGetFirstColIndex()+1;
   integer iFirstIndex = iGetFirstIndex();
   
   WM.fPutRowIndex(1, iNode1RowIndex);
   WM.fPutRowIndex(2, iNode2RowIndex);
   WM.fPutRowIndex(3, iNode3RowIndex);
   WM.fPutColIndex(1, iNode1ColIndex);
   WM.fPutColIndex(2, iNode2ColIndex);
   WM.fPutColIndex(3, iNode3ColIndex);

   WM.fPutRowIndex(4, iFirstIndex+1);
   WM.fPutColIndex(4, iFirstIndex+1);
   WM.fPutRowIndex(5, iFirstIndex+2);
   WM.fPutColIndex(5, iFirstIndex+2);

   doublereal p1 = pNode1->dGetX();
   doublereal p2 = pNode2->dGetX();
   doublereal p3 = pNode3->dGetX();
   s = XCurr.dGetCoef(iFirstIndex+1);  /* spostamento */
   v = XCurr.dGetCoef(iFirstIndex+2);  /* velocita' */
   
   doublereal Cd = .6;        /* verificare */
   doublereal jumpPres12 = fabs(p1-p2);
   doublereal jumpPres23 = fabs(p2-p3);
   doublereal jumpPres13 = fabs(p1-p3);
   doublereal density = HF->dGetDensity();
    
   /* evito di dividere per un numero troppo piccolo */
   if (jumpPres12 < 1.e8*DBL_EPSILON) {
      jumpPres12 = 1.e8*DBL_EPSILON;
   }
   if (jumpPres23 < 1.e8*DBL_EPSILON) {
      jumpPres23 = 1.e8*DBL_EPSILON;
   }
   if (jumpPres13 < 1.e8*DBL_EPSILON) {
      jumpPres13 = 1.e8*DBL_EPSILON;
   }
 
   doublereal Jac11 = 0.;
   doublereal Jac12 = 0.;
   doublereal Jac13 = 0.;
   doublereal Jac14 = 0.;
   doublereal Jac15 = 0.;
   doublereal Jac21 = 0.;
   doublereal Jac22 = 0.;
   doublereal Jac23 = 0.;
   doublereal Jac24 = 0.;
   doublereal Jac25 = 0.;
   doublereal Jac31 = 0.;
   doublereal Jac32 = 0.;
   doublereal Jac33 = 0.;
   doublereal Jac34 = 0.;
   doublereal Jac35 = 0.;
   doublereal Jac41 = 0.;
   doublereal Jac42 = 0.;
   doublereal Jac43 = 0.;
   doublereal Jac44 = 0.;
   doublereal Jac45 = 0.;
   doublereal Jac51 = 0.;
   doublereal Jac52 = 0.;
   doublereal Jac53 = 0.;
   doublereal Jac54 = 0.;
   doublereal Jac55 = 0.;

#ifdef HYDR_DEVEL
   DEBUGCOUT("Valore di p1: " << p1 << endl);
   DEBUGCOUT("Valore di p2: " << p2 << endl);
   DEBUGCOUT("Valore di p3: " << p3 << endl);
   DEBUGCOUT("Valore di s: " << s << endl);
#endif
   
   doublereal s12 = sqrt(2.*jumpPres12/density);
   doublereal s13 = sqrt(2.*jumpPres13/density);
   
   Jac11 = -width*s*Cd/s12-area_diaf*Cd/s13;
   Jac12 = width*s*Cd/s12;
   Jac13 = area_diaf*Cd/s13;
   Jac14 = -density*width*Cd*copysign(s12, p1-p2)*dCoef;
   Jac15 = 0.;
   
   Jac21 = width*s*Cd/s12;
   Jac22 = -width*s*Cd/s12;
   Jac23 = 0.;
   Jac24 = density*width*Cd*copysign(s12, p1-p2)*dCoef;
   Jac25 = 0.;
   
   Jac31 = area_diaf*Cd/s13;
   Jac32 = 0.;
   Jac33 = -area_diaf*Cd/s13;
   Jac34 = 0.;
   Jac35 = 0.;
   
   Jac41 = area_max;
   Jac42 = 0.;
   Jac43 =-area_max;
   Jac44 = -Kappa*dCoef-c1*dCoef-cf1*dCoef;
   
   
   doublereal Jac45old1;
   doublereal Jac45old2;
   Jac45old1 = -mass-c2*dCoef-c3-cf2*dCoef-cf3
     -density*area_max*pow(area_max/(area_diaf*Cd), 2.)*fabs(v)*dCoef;
   Jac45 = -mass-c2*dCoef-c3-cf2*dCoef-cf3
     -h*density*pow(area_max/area_pipe, 2.)*fabs(v)*dCoef;
   Jac45old2 = -mass-c2*dCoef-c3-cf2*dCoef-cf3-44.4*fabs(v)*dCoef;
   
   Jac51 = 0.;
   Jac52 = 0.;
   Jac53 = 0.;
   Jac54 = -1.;
   Jac55 = dCoef;
   
#ifdef HYDR_DEVEL
   DEBUGCOUT("Jac smorzatore " 
	     << density*area_max*pow(area_max/(area_diaf*Cd), 2) << endl);
   DEBUGCOUT("Jac smorzatore " << h*density << endl);
   DEBUGCOUT("Jac smorzatore Jac45old1 " << Jac45old1 << endl);
   DEBUGCOUT("Jac smorzatore Jac45old2 " << Jac45old2 << endl);
   DEBUGCOUT("Jac smorzatore Jac45     " << Jac45 << endl);
   
   DEBUGCOUT("Jac11: " << Jac11 << endl);
   DEBUGCOUT("Jac12: " << Jac12 << endl);
   DEBUGCOUT("Jac13: " << Jac13 << endl);
   DEBUGCOUT("Jac14: " << Jac14 << endl);
   DEBUGCOUT("Jac15: " << Jac15 << endl);
   DEBUGCOUT("Jac21: " << Jac21 << endl);
   DEBUGCOUT("Jac22: " << Jac22 << endl);
   DEBUGCOUT("Jac23: " << Jac23 << endl);
   DEBUGCOUT("Jac24: " << Jac24 << endl);
   DEBUGCOUT("Jac25: " << Jac25 << endl);
   DEBUGCOUT("Jac31: " << Jac31 << endl);
   DEBUGCOUT("Jac32: " << Jac32 << endl);
   DEBUGCOUT("Jac33: " << Jac33 << endl);
   DEBUGCOUT("Jac34: " << Jac34 << endl);
   DEBUGCOUT("Jac35: " << Jac35 << endl);
   DEBUGCOUT("Jac41: " << Jac41 << endl);
   DEBUGCOUT("Jac42: " << Jac42 << endl);
   DEBUGCOUT("Jac43: " << Jac43 << endl);
   DEBUGCOUT("Jac44: " << Jac44 << endl);
   DEBUGCOUT("Jac45: " << Jac45 << endl);
   DEBUGCOUT("Jac51: " << Jac51 << endl);
   DEBUGCOUT("Jac52: " << Jac52 << endl);
   DEBUGCOUT("Jac53: " << Jac53 << endl);
   DEBUGCOUT("Jac54: " << Jac54 << endl);
   DEBUGCOUT("Jac55: " << Jac55 << endl);
#endif
   
   WM.fPutCoef(1, 1, Jac11);
   WM.fPutCoef(1, 2, Jac12);
   WM.fPutCoef(1, 3, Jac13);
   WM.fPutCoef(1, 4, Jac14);
   WM.fPutCoef(1, 5, Jac15);
   WM.fPutCoef(2, 1, Jac21);
   WM.fPutCoef(2, 2, Jac22);
   WM.fPutCoef(2, 3, Jac23);
   WM.fPutCoef(2, 4, Jac24);
   WM.fPutCoef(2, 5, Jac25);
   WM.fPutCoef(3, 1, Jac31);
   WM.fPutCoef(3, 2, Jac32);
   WM.fPutCoef(3, 3, Jac33);
   WM.fPutCoef(3, 4, Jac34);
   WM.fPutCoef(3, 5, Jac35);
   WM.fPutCoef(4, 1, Jac41);
   WM.fPutCoef(4, 2, Jac42);
   WM.fPutCoef(4, 3, Jac43);
   WM.fPutCoef(4, 4, Jac44);
   WM.fPutCoef(4, 5, Jac45);
   WM.fPutCoef(5, 1, Jac51);
   WM.fPutCoef(5, 2, Jac52);
   WM.fPutCoef(5, 3, Jac53);
   WM.fPutCoef(5, 4, Jac54);
   WM.fPutCoef(5, 5, Jac55);
   
   return WorkMat;
}

SubVectorHandler& 
Flow_valve::AssRes(SubVectorHandler& WorkVec,
		   doublereal dCoef,
		   const VectorHandler& XCurr, 
		   const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering Flow_valve::AssRes()" << endl);
   
   WorkVec.Resize(5);
   
   integer iNode1RowIndex = pNode1->iGetFirstRowIndex()+1;
   integer iNode2RowIndex = pNode2->iGetFirstRowIndex()+1;
   integer iNode3RowIndex = pNode3->iGetFirstRowIndex()+1;
   
   doublereal p1 = pNode1->dGetX();
   doublereal p2 = pNode2->dGetX();
   doublereal p3 = pNode3->dGetX();
   integer iFirstIndex = iGetFirstIndex();

   s = XCurr.dGetCoef(iFirstIndex+1);        /* spostamento */
   v = XCurr.dGetCoef(iFirstIndex+2);        /* velocita' */
   sp = XPrimeCurr.dGetCoef(iFirstIndex+1);  /* velocita' */
   vp = XPrimeCurr.dGetCoef(iFirstIndex+2);  /* accelerazione */
  
   doublereal Cd = .6;
   doublereal jumpPres12 = fabs(p1-p2);
   doublereal jumpPres23 = fabs(p2-p3);
   doublereal jumpPres13 = fabs(p1-p3);
   
   doublereal Res_1 = 0.;
   doublereal Res_2 = 0.;
   doublereal Res_3 = 0.;  
   doublereal Res_4 = 0.;
   doublereal Res_5 = 0.;
   doublereal density = HF->dGetDensity();
   
   doublereal x0 = ((p1-p3)*area_max-force0)/Kappa;
   
   if (s <= 0.) { 
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

      c4 = 0.; /* (force0-p1*area_max); */
   } else {
      c1 = 0.;
      c2 = 0.;
      c3 = 0.;
      c4 = 0.;
   }
   
   if (s >= s_max) { 
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

      cf4 = 0.; /* (force0-p1*area_max); */
   } else {
      cf1 = 0.;
      cf2 = 0.;
      cf3 = 0.;
      cf4 = 0.;
   }

   doublereal s12 = sqrt(2.*jumpPres12/density);
   doublereal s13 = sqrt(2.*jumpPres13/density);
   
   Res_1 = density*width*s*Cd*copysign(s12, p1-p2)+
     density*area_diaf*Cd*copysign(s13, p1-p3);
   Res_2 = -density*width*s*Cd*copysign(s12, p1-p2);
   Res_3 = -density*area_diaf*Cd*copysign(s13, p1-p3);

   doublereal Res_4old;
   Res_4old = mass*vp-area_max*p1+force0+Kappa*s+area_max*p3+c1*s+c2*sp+c3*vp
     +c4+cf1*(s-s_max)+cf2*sp+cf3*vp+cf4
     +copysign(.5*density*area_max*pow(sp*area_max/(area_diaf*Cd), 2), sp);
   
   Res_4 = mass*vp-area_max*p1+force0+Kappa*s+area_max*p3+c1*s+c2*sp+c3*vp
     +c4+cf1*(s-s_max)+cf2*sp+cf3*vp+cf4
     +copysign(.5*h*density*pow(sp*area_max/area_pipe, 2), sp);
   
   Res_5 =sp-v;
   
   flow1=-Res_1;   /* per l'output */
   flow2=-Res_2;   /* per l'output */
   flow3=-Res_3;   /* per l'output */
   
#ifdef HYDR_DEVEL
   DEBUGCOUT("Res_4:         " <<  Res_4 << endl);
   DEBUGCOUT("Res_4old:      " <<  Res_4old << endl);
   DEBUGCOUT("smorazatore:   " 
	     << copysign(.5*density*area_max*pow(sp*area_max/(area_diaf*Cd), 2), sp) << endl);
   DEBUGCOUT("smorzatoreold: " 
	     << copysign(.5*h*density*pow(sp*area_max/area_pipe, 2), sp) << endl);
   
   DEBUGCOUT("width:         " << width << endl);
   DEBUGCOUT("Cd:            " << Cd << endl);
   DEBUGCOUT("jumpPres12:    " << jumpPres12  << endl);
   DEBUGCOUT("jumpPres23:    " << jumpPres23  << endl);
   DEBUGCOUT("jumpPres13:    " << jumpPres13  << endl);
   DEBUGCOUT("density:       " << density << endl);
   DEBUGCOUT("p1:            " << p1 << endl);
   DEBUGCOUT("p2:            " << p2 << endl);
   DEBUGCOUT("p3:            " << p3 << endl);
   DEBUGCOUT("s:             " << s << endl);
   DEBUGCOUT("sp:            " << sp << endl);
   DEBUGCOUT("v:             " << v << endl);
   DEBUGCOUT("vp:            " << vp << endl);
   DEBUGCOUT("area_max:      " << area_max << endl);
   DEBUGCOUT("Kappa:         " << Kappa << endl);
   DEBUGCOUT("force0:        " << force0 << endl);
   DEBUGCOUT("mass:          " << mass << endl);
   DEBUGCOUT("s_max:         " << s_max << endl);
   DEBUGCOUT("x0:            " << x0  << endl);
   DEBUGCOUT("c1:            " << c1 << endl);
   DEBUGCOUT("c2:            " << c2 << endl);
   DEBUGCOUT("c3:            " << c3 << endl);
   DEBUGCOUT("c4:            " << c4 << endl);
   DEBUGCOUT("cf1:           " << cf1 << endl);
   DEBUGCOUT("cf2:           " << cf2 << endl);
   DEBUGCOUT("cf3:           " << cf3 << endl);
   DEBUGCOUT("cf4:           " << cf4 << endl);
   DEBUGCOUT("PORTATE AI VARI NODI (positive se entranti)" << endl);
   DEBUGCOUT("-Res_1 (portata nodo1): " << -Res_1 << endl); 
   DEBUGCOUT("-Res_2 (portata nodo2): " << -Res_2 << endl); 
   DEBUGCOUT("-Res_3:(portata nodo3): " << -Res_3 << endl); 
   DEBUGCOUT("-Res_4:                 " << -Res_4 << endl); 
   DEBUGCOUT("-Res_5:                 " << -Res_5 << endl); 
#endif
   
   WorkVec.fPutItem(1, iNode1RowIndex, Res_1);
   WorkVec.fPutItem(2, iNode2RowIndex, Res_2);         
   WorkVec.fPutItem(3, iNode3RowIndex, Res_3);
   WorkVec.fPutItem(4, iFirstIndex+1, Res_4);    
   WorkVec.fPutItem(5, iFirstIndex+2, Res_5);    
     
   return WorkVec;
}
  
void Flow_valve::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) { 
      ostream& out = OH.Hydraulic();
      out << setw(8) << GetLabel()
	<< " " << s  << " " << v  << " "<< vp  
	<< " " << flow1  << " "<< flow2  << " "<< flow3  << endl;
   }  
}

void Flow_valve::SetValue(VectorHandler&  X  , VectorHandler&  XP  ) const 
{   
   integer i = iGetFirstIndex();
   X.fPutCoef(i+1, 0.);
   X.fPutCoef(i+2, 0.);
   XP.fPutCoef(i+1, 0.);
   XP.fPutCoef(i+2, 0.);  
}

/* Flow_valve - end */
