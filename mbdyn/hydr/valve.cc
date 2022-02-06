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

/* 
 * Copyright 1999-2000 Lamberto Puggelli <puggelli@tiscalinet.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cfloat>
#include <limits>

#include "valve.h"

/* Control_valve - begin */

Control_valve::Control_valve(unsigned int uL, const DofOwner* pDO,
			     HydraulicFluid* hf,
			     const PressureNode* p1, const PressureNode* p2, 
			     const PressureNode* p3, const PressureNode* p4,  
			     doublereal A_max, doublereal Loss_A, const DriveCaller* pDC,
			     flag fOut) 
: Elem(uL, fOut),
HydraulicElem(uL, pDO, hf, fOut),
DriveOwner(pDC),
pNode1(p1), pNode2(p2), pNode3(p3), pNode4(p4),
area_max(A_max), loss_area(Loss_A),
A1min(2.*area_max*loss_area),
A2min(area_max*loss_area),
A3min(2.*area_max*loss_area),
A4min(area_max*loss_area)
{
   ASSERT(pNode1 != NULL);
   ASSERT(pNode1->GetNodeType() == Node::HYDRAULIC);
   ASSERT(pNode2 != NULL);
   ASSERT(pNode2->GetNodeType() == Node::HYDRAULIC);
   ASSERT(pNode3 != NULL);
   ASSERT(pNode3->GetNodeType() == Node::HYDRAULIC);
   ASSERT(pNode4 != NULL);
   ASSERT(pNode4->GetNodeType() == Node::HYDRAULIC);
   ASSERT(A_max > std::numeric_limits<doublereal>::epsilon());
   ASSERT(loss_area >= 0.);
   
   /* 
    * Cd = pi / ( pi + 2 ) ~= .611
    * 
    * cfr. Merritt, Hydraulic Control Systems, Par. 3.4, pp. 39-45
    * John Wiley & Sons, New York, 1967
    */
   Cd = M_PI / ( M_PI + 2. ); 
}

Control_valve::~Control_valve(void)
{
   NO_OP;
}
   
/* Tipo di elemento idraulico (usato solo per debug ecc.) */
HydraulicElem::Type Control_valve::GetHydraulicType(void) const 
{
   return HydraulicElem::CONTROL_VALVE;
}

/* Contributo al file di restart */
std::ostream& Control_valve::Restart(std::ostream& out) const
{
   return out << "Control_valve not implemented yet!" << std::endl;
}
   
unsigned int Control_valve::iGetNumDof(void) const 
{
   return 0;
}
   
DofOrder::Order Control_valve::GetDofType(unsigned int i) const 
{
   silent_cerr("ControlValve(" << GetLabel() << ") has no dofs!" << std::endl);
   throw ErrGeneric(MBDYN_EXCEPT_ARGS);
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
   DEBUGCOUT("Entering Control_valve::AssJac()" << std::endl);
   
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
  
   WM.PutRowIndex(1, iNode1RowIndex);
   WM.PutRowIndex(2, iNode2RowIndex);
   WM.PutRowIndex(3, iNode3RowIndex);
   WM.PutRowIndex(4, iNode4RowIndex);
   WM.PutColIndex(1, iNode1ColIndex);
   WM.PutColIndex(2, iNode2ColIndex);
   WM.PutColIndex(3, iNode3ColIndex);
   WM.PutColIndex(4, iNode4ColIndex);   
   
   doublereal p1 = pNode1->dGetX();
   doublereal p2 = pNode2->dGetX();
   doublereal p3 = pNode3->dGetX();
   doublereal p4 = pNode4->dGetX();
   doublereal density = HF->dGetDensity();
 
   doublereal jumpPres12 = fabs(p1-p2); /* salto di pressione nodo1 & nodo3 */
   doublereal jumpPres13 = fabs(p1-p3); /* salto di pressione nodo1 & nodo4 */
   doublereal jumpPres24 = fabs(p2-p4); /* salto di pressione nodo2 & nodo3 */
   doublereal jumpPres34 = fabs(p3-p4); /* salto di pressione nodo2 & nodo4 */
   
   if (jumpPres12 < 1.e8*std::numeric_limits<doublereal>::epsilon()) {
      jumpPres12 = 1.e8*std::numeric_limits<doublereal>::epsilon();
   }
   if (jumpPres13 < 1.e8*std::numeric_limits<doublereal>::epsilon()) {
      jumpPres13 = 1.e8*std::numeric_limits<doublereal>::epsilon();
   }
   if (jumpPres24 < 1.e8*std::numeric_limits<doublereal>::epsilon()) {
      jumpPres24 = 1.e8*std::numeric_limits<doublereal>::epsilon();
   }
   if (jumpPres34 < 1.e8*std::numeric_limits<doublereal>::epsilon()) {
      jumpPres34 = 1.e8*std::numeric_limits<doublereal>::epsilon();
   }
   
   doublereal primo = Cd*A1/sqrt(2*jumpPres12/density);
   doublereal secondo = Cd*A2/sqrt(2*jumpPres13/density);
   doublereal terzo = Cd*A4/sqrt(2*jumpPres24/density);
   doublereal quarto = Cd*A3/sqrt(2*jumpPres34/density);

   doublereal Jac11 = -primo-secondo;
   doublereal Jac12 = primo;
   doublereal Jac13 = secondo;
   // doublereal Jac14 = 0.;
   doublereal Jac21 = primo;
   doublereal Jac22 = -primo-terzo;
   // doublereal Jac23 = 0.;
   doublereal Jac24 = terzo;
   doublereal Jac31 = secondo;
   // doublereal Jac32 = 0.;
   doublereal Jac33 = -secondo-quarto;
   doublereal Jac34 = quarto;
   // doublereal Jac41 = 0.;
   doublereal Jac42 = terzo;
   doublereal Jac43 = quarto;
   doublereal Jac44 = -terzo-quarto;

   WM.PutCoef(1, 1, Jac11);
   WM.PutCoef(1, 2, Jac12);
   WM.PutCoef(1, 3, Jac13);
   // WM.PutCoef(1, 4, Jac14);
   WM.PutCoef(2, 1, Jac21);
   WM.PutCoef(2, 2, Jac22);
   // WM.PutCoef(2, 3, Jac23);
   WM.PutCoef(2, 4, Jac24);
   WM.PutCoef(3, 1, Jac31);
   // WM.PutCoef(3, 2, Jac32);
   WM.PutCoef(3, 3, Jac33);
   WM.PutCoef(3, 4, Jac34);
   // WM.PutCoef(4, 1, Jac41);
   WM.PutCoef(4, 2, Jac42);
   WM.PutCoef(4, 3, Jac43);
   WM.PutCoef(4, 4, Jac44);

   return WorkMat;
}


SubVectorHandler& 
Control_valve::AssRes(SubVectorHandler& WorkVec,
		      doublereal dCoef,
		      const VectorHandler& XCurr, 
		      const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering Control_valve::AssRes()" << std::endl);
   
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

   doublereal dp12 = p1-p2;
   doublereal dp13 = p1-p3;
   doublereal dp24 = p2-p4;
   doublereal dp34 = p3-p4;
   
   doublereal jumpPres12 = fabs(dp12); /* salto di pressione nodo1 & nodo2 */
   doublereal jumpPres13 = fabs(dp13); /* salto di pressione nodo1 & nodo3 */
   doublereal jumpPres24 = fabs(dp24); /* salto di pressione nodo2 & nodo4 */
   doublereal jumpPres34 = fabs(dp34); /* salto di pressione nodo3 & nodo4 */

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

   doublereal Q12 = Cd*A1*copysign(sqrt(2.*jumpPres12*density), dp12);
   doublereal Q13 = Cd*A2*copysign(sqrt(2.*jumpPres13*density), dp13);
   doublereal Q24 = Cd*A4*copysign(sqrt(2.*jumpPres24*density), dp24);
   doublereal Q34 = Cd*A3*copysign(sqrt(2.*jumpPres34*density), dp34);

   doublereal Res_1 = Q12+Q13;
   doublereal Res_2 = -Q12+Q24;
   doublereal Res_3 = -Q13+Q34;
   doublereal Res_4 = -Q34-Q24;
        
   flow1 = -Res_1;
   flow2 = -Res_2;
   flow3 = -Res_3;
   flow4 = -Res_4;

#ifdef HYDR_DEVEL
   DEBUGCOUT("Stato:   " << Stato << std::endl);
   DEBUGCOUT("A1:      " << A1 << std::endl);
   DEBUGCOUT("A2:      " << A2  << std::endl);
   DEBUGCOUT("A3:      " << A3 << std::endl);
   DEBUGCOUT("A4:      " << A4 << std::endl);
   DEBUGCOUT("p1:      " << p1 << std::endl);
   DEBUGCOUT("p2:      " << p2 << std::endl);
   DEBUGCOUT("p3:      " << p3 << std::endl);
   DEBUGCOUT("p4:      " << p4 << std::endl);
   DEBUGCOUT("Cd:      " << Cd << std::endl);
   DEBUGCOUT("density: " << density << std::endl);
   DEBUGCOUT("Area_max:             " << area_max << std::endl);
   DEBUGCOUT("Loss_area in %:       " << loss_area << std::endl);
   DEBUGCOUT("Q12:     " << Q12 << std::endl);
   DEBUGCOUT("Q13:     " << Q13 << std::endl);
   DEBUGCOUT("Q24:     " << Q24 << std::endl);
   DEBUGCOUT("Q34:     " << Q34 << std::endl);
   DEBUGCOUT("PORTATE AI VARI NODI (positive se entranti)"<< std::endl);
   DEBUGCOUT("-Res_1 (portata nodo1): " << -Res_1 << std::endl); 
   DEBUGCOUT("-Res_2 (portata nodo2): " << -Res_2 << std::endl); 
   DEBUGCOUT("-Res_3 (portata nodo3): " << -Res_3 << std::endl); 
   DEBUGCOUT("-Res_4 (portata nodo4): " << -Res_4 << std::endl); 
#endif
   
   WorkVec.PutItem(1, iNode1RowIndex, Res_1);
   WorkVec.PutItem(2, iNode2RowIndex, Res_2);         
   WorkVec.PutItem(3, iNode3RowIndex, Res_3);
   WorkVec.PutItem(4, iNode4RowIndex, Res_4);    

   return WorkVec;
}
  
void Control_valve::Output(OutputHandler& OH) const
{
   if (bToBeOutput()) { 
      std::ostream& out = OH.Hydraulic();
      out << std::setw(8) << GetLabel()
	<< " " << Stato 
	<< " " << flow1 << " " << flow2 
	<< " " << flow3 << " " << flow4 << std::endl;
   }   
}

const OutputHandler::Dimensions 
Control_valve::GetEquationDimension(integer index) const {
   // DOF == 0
   return OutputHandler::Dimensions::UnknownDimension;
}

std::ostream&
Control_valve::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{

	out
		<< "It does not have any DOF" << std::endl;

	return out;
}
/* Control_valve - end */


/* Control_valve2 - begin */

Control_valve2::Control_valve2(unsigned int uL, const DofOwner* pDO,
			     HydraulicFluid* hf,
			     const PressureNode* p1, const PressureNode* p2, 
			     const PressureNode* p3, const PressureNode* p4,  
			     doublereal A_max, doublereal Loss_A, 
			     const DriveCaller* pDC,
			     flag fOut) 
: Elem(uL, fOut),
HydraulicElem(uL, pDO, hf, fOut),
DriveOwner(pDC),
area_max(A_max), loss_area(Loss_A), area_min(area_max*loss_area)
{
	pNode[N1] = p1;
	pNode[N2] = p2;
	pNode[N3] = p3;
	pNode[N4] = p4;
	
	ASSERT(pNode[N1] != NULL);
	ASSERT(pNode[N1]->GetNodeType() == Node::HYDRAULIC);
	ASSERT(pNode[N2] != NULL);
	ASSERT(pNode[N2]->GetNodeType() == Node::HYDRAULIC);
	ASSERT(pNode[N3] != NULL);
	ASSERT(pNode[N3]->GetNodeType() == Node::HYDRAULIC);
	ASSERT(pNode[N4] != NULL);
	ASSERT(pNode[N4]->GetNodeType() == Node::HYDRAULIC);

	ASSERT(area_max > std::numeric_limits<doublereal>::epsilon());
	ASSERT(loss_area > 1.e-9);	/* 
					 * se = 0. occorre fare un elemento
					 * apposta con solo 2 dof
					 */
   
	Cd = .611; /* coefficiente di perdita */
}

Control_valve2::~Control_valve2(void)
{
	NO_OP;
}
   
/* Tipo di elemento idraulico (usato solo per debug ecc.) */
HydraulicElem::Type 
Control_valve2::GetHydraulicType(void) const 
{
	return HydraulicElem::CONTROL_VALVE;
}

/* Contributo al file di restart */
std::ostream& 
Control_valve2::Restart(std::ostream& out) const
{
	return out << "Control_valve2 not implemented yet!" << std::endl;
}
   
unsigned int 
Control_valve2::iGetNumDof(void) const 
{
	return LAST_Q;
}
   
DofOrder::Order Control_valve2::GetDofType(unsigned int i) const 
{
	ASSERT(i >= 0 && i < iGetNumDof());
	return DofOrder::ALGEBRAIC;
}

void 
Control_valve2::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const 
{
	*piNumRows = 4+LAST_Q; 
	*piNumCols = 4+LAST_Q;
}

VariableSubMatrixHandler& 
Control_valve2::AssJac(VariableSubMatrixHandler& WorkMat,
		      doublereal dCoef,
		      const VectorHandler& XCurr, 
		      const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering Control_valve2::AssJac()" << std::endl);
	
	FullSubMatrixHandler& WM = WorkMat.SetFull();
	integer iNumRows, iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	integer iFirstIndex = iGetFirstIndex();
	
	for (int i = 0; i < LAST_N; i++) {
		WM.PutRowIndex(1+i, pNode[i]->iGetFirstRowIndex()+1);
		WM.PutColIndex(1+i, pNode[i]->iGetFirstColIndex()+1);
	}

	for (int i = 1; i <= LAST_Q; i++) {
		WM.PutRowIndex(4+i, iFirstIndex+i);
		WM.PutColIndex(4+i, iFirstIndex+i);
	}
	
	doublereal dKappa = 2.*Cd*Cd*density;

	/* Q12 */
	WM.PutCoef(1, 5, -1.);
	WM.PutCoef(2, 5,  1.);
	
	doublereal a = A[0]*A[0];
	WM.PutCoef(5, 1, -a);
	WM.PutCoef(5, 2,  a);

	/* Q34 */
	WM.PutCoef(3, 6, -1.);
	WM.PutCoef(4, 6,  1.);
	
	a = A[1]*A[1];
	WM.PutCoef(6, 3, -a);
	WM.PutCoef(6, 4,  a);
	
	/* Q13 */
	WM.PutCoef(1, 7, -1.);
	WM.PutCoef(3, 7,  1.);
	
	a = A[2]*A[2];
	WM.PutCoef(7, 1, -a);
	WM.PutCoef(7, 3,  a);
	
	/* Q24 */
	WM.PutCoef(2, 8, -1.);
	WM.PutCoef(4, 8,  1.);
	
	a = A[3]*A[3];
	WM.PutCoef(8, 2, -a);
	WM.PutCoef(8, 4,  a);

#ifdef VALVE_6
	/* Q14 */
	WM.PutCoef(1, 9, -1.);
	WM.PutCoef(4, 9,  1.);
	
	a = A[4]*A[4];
	WM.PutCoef(9, 1, -a);
	WM.PutCoef(9, 4,  a);
	
	/* Q23 */
	WM.PutCoef(2, 10, -1.);
	WM.PutCoef(3, 10,  1.);
	
	a = A[5]*A[5];
	WM.PutCoef(10, 2, -a);
	WM.PutCoef(10, 3,  a);
#endif /* VALVE_6 */

	for (int i = 0; i < LAST_Q; i++) {
		WM.PutCoef(5+i, 5+i, 2.*fabs(q[i])/dKappa);
	}
	
	return WorkMat;
}

void
Control_valve2::Prepare(void)
{
	doublereal p[LAST_N];
	doublereal pm = 0.;

	for (int i = 0; i < LAST_N; i++) {
		p[i] = pNode[i]->dGetX();
		pm += p[i];
	}
	pm /= 4.;

	density = HF->dGetDensity(pm);
	
	dp[Q12] = p[N1]-p[N2];
	dp[Q34] = p[N3]-p[N4];
	dp[Q13] = p[N1]-p[N3];
	dp[Q24] = p[N2]-p[N4];
#ifdef VALVE_6
	dp[Q14] = p[N1]-p[N4];
	dp[Q23] = p[N2]-p[N3];
#endif /* VALVE_6 */

	Stato = pGetDriveCaller()->dGet();
	if (Stato > 0.) { 
		if (Stato > 1.) {
			Stato = 1.;
		}

		A[Q12] = Stato*area_max+2.*area_min;
		A[Q34] = Stato*area_max+2.*area_min;
		A[Q13] = area_min;
		A[Q24] = area_min;
	} else {
		if (Stato < -1.) {
			Stato = -1.;
		}

		A[Q12] = 2.*area_min;
		A[Q34] = 2.*area_min;
		A[Q13] = -Stato*area_max+area_min;
		A[Q24] = -Stato*area_max+area_min;
	}
	
	/*
	 * vale sempre area_min (lo genero ogni volta perche' magari
	 * la si puo' far dipendere da qualche parametro
	 */
#ifdef VALVE_6
	A[Q14] = area_min;
	A[Q23] = area_min;
#endif /* VALVE_6 */
}

SubVectorHandler& 
Control_valve2::AssRes(SubVectorHandler& WorkVec,
		      doublereal dCoef,
		      const VectorHandler& XCurr, 
		      const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering Control_valve2::AssRes()" << std::endl);
	
	integer iNumRows, iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.Resize(iNumRows);
	
	integer iFirstIndex = iGetFirstIndex()+1;

	integer iNodeRowIndex[LAST_N];
	for (int i = 0; i < LAST_N; i++) {
		iNodeRowIndex[i] = pNode[i]->iGetFirstRowIndex()+1;
	}

	Prepare();
	doublereal dKappa = 2.*Cd*Cd*density;
	
	for (int i = 0; i < LAST_Q; i++) {
		q[i] = XCurr(iFirstIndex+i);
	}

#ifdef VALVE_6
	f[N1] = q[Q12]+q[Q13]+q[Q14];
	f[N2] = -q[Q12]+q[Q23]+q[Q24];
	f[N3] = -q[Q13]-q[Q23]+q[Q34];
	f[N4] = -q[Q14]-q[Q24]-q[Q34];
#else /* !VALVE_6 */
	f[N1] = q[Q12]+q[Q13];
	f[N2] = -q[Q12]+q[Q24];
	f[N3] = -q[Q13]+q[Q34];
	f[N4] = -q[Q24]-q[Q34];
#endif /* !VALVE_6 */

	for (int i = 0; i < LAST_N; i++) {
		WorkVec.PutItem(1+i, iNodeRowIndex[i], f[i]);
	}
	
	for (int i = 0; i < LAST_Q; i++) {
		WorkVec.PutItem(5+i, iFirstIndex+i, 
			A[i]*A[i]*dp[i]-q[i]*fabs(q[i])/dKappa);
	}

	return WorkVec;
}
  
void
Control_valve2::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) { 
		std::ostream& out = OH.Hydraulic();
		out << std::setw(8) << GetLabel()
			<< " " << Stato 
			<< " " << -f[N1] << " " << -f[N2] 
			<< " " << -f[N3] << " " << -f[N4] << std::endl;
	}   
}

void 
Control_valve2::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& /* XP */ ,
		SimulationEntity::Hints *ph)
{
	integer iFirstIndex = iGetFirstIndex()+1;
	
	const_cast<Control_valve2 *>(this)->Prepare();

	for (int i = 0; i < LAST_Q; i++) {

		/*
		 * q = sign(Dp)*Cd*A*sqrt(2.*rho*abs(Dp))
		 */
		
		X.PutCoef(iFirstIndex+i, 
			Cd*A[i]*copysign(sqrt(2.*density*fabs(dp[i])), dp[i]));
	}
}

const OutputHandler::Dimensions 
Control_valve2::GetEquationDimension(integer index) const {
   // DOF == LAST_Q = 6
   
   return OutputHandler::Dimensions::Force;
}

std::ostream&
Control_valve2::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{

	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << "->" << iIndex + LAST_Q << ": " <<
			"control valve 2 force balance" << std::endl;

	return out;
}

/* Control_valve2 - end */


/* Dynamic_control_valve - begin */

Dynamic_control_valve::Dynamic_control_valve(unsigned int uL, 
					     const DofOwner* pDO, 
					     HydraulicFluid* hf,
					     const PressureNode* p1, 
					     const PressureNode* p2, 
					     const PressureNode* p3, 
					     const PressureNode* p4,
					     const DriveCaller* pDC,
					     doublereal s0, doublereal s_mx,
					     doublereal W, doublereal Loss_A, 
					     doublereal Valve_d,
					     doublereal Valve_rho,
					     doublereal cs, doublereal cv, 
					     doublereal ca, flag fOut) 
: Elem(uL, fOut),
HydraulicElem(uL, pDO, hf, fOut), DriveOwner(pDC),
pNode1(p1), pNode2(p2), pNode3(p3), pNode4(p4),
start(s0), c_spost(cs), c_vel(cv), c_acc(ca),
width(W), loss_area(Loss_A),
valve_diameter(Valve_d), valve_density(Valve_rho), s_max(s_mx)
{
   ASSERT(pNode1 != NULL);
   ASSERT(pNode1->GetNodeType() == Node::HYDRAULIC);
   ASSERT(pNode2 != NULL);
   ASSERT(pNode2->GetNodeType() == Node::HYDRAULIC);
   ASSERT(pNode3 != NULL);
   ASSERT(pNode3->GetNodeType() == Node::HYDRAULIC);
   ASSERT(pNode4 != NULL);
   ASSERT(pNode4->GetNodeType() == Node::HYDRAULIC);
   ASSERT(s0 >= 0.);
   ASSERT(Valve_rho > std::numeric_limits<doublereal>::epsilon());
   ASSERT(Valve_d > std::numeric_limits<doublereal>::epsilon());
   ASSERT(W > std::numeric_limits<doublereal>::epsilon());
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
HydraulicElem::Type Dynamic_control_valve::GetHydraulicType(void) const 
{
   return HydraulicElem::DYNAMIC_CONTROL_VALVE;
}

/* Contributo al file di restart */
std::ostream& Dynamic_control_valve::Restart(std::ostream& out) const
{
   return out << "Dynamic_control_valve not implemented yet!" << std::endl;
}
   
unsigned int Dynamic_control_valve::iGetNumDof(void) const 
{
   return 2;
}
   
DofOrder::Order Dynamic_control_valve::GetDofType(unsigned int i) const 
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
   DEBUGCOUT("Entering Control_valve::AssJac()" << std::endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.ResizeReset(6, 6);
   
   integer iNode1RowIndex = pNode1->iGetFirstRowIndex()+1;
   integer iNode2RowIndex = pNode2->iGetFirstRowIndex()+1;
   integer iNode3RowIndex = pNode3->iGetFirstRowIndex()+1;
   integer iNode4RowIndex = pNode4->iGetFirstRowIndex()+1;
   integer iNode1ColIndex = pNode1->iGetFirstColIndex()+1;
   integer iNode2ColIndex = pNode2->iGetFirstColIndex()+1;
   integer iNode3ColIndex = pNode3->iGetFirstColIndex()+1;
   integer iNode4ColIndex = pNode4->iGetFirstColIndex()+1;
   integer iFirstIndex = iGetFirstIndex();
   
   WM.PutRowIndex(1, iNode1RowIndex);
   WM.PutRowIndex(2, iNode2RowIndex);
   WM.PutRowIndex(3, iNode3RowIndex);
   WM.PutRowIndex(4, iNode4RowIndex);
   WM.PutColIndex(1, iNode1ColIndex);
   WM.PutColIndex(2, iNode2ColIndex);
   WM.PutColIndex(3, iNode3ColIndex);
   WM.PutColIndex(4, iNode4ColIndex);
   
   WM.PutRowIndex(5, iFirstIndex+1);
   WM.PutColIndex(5, iFirstIndex+1);
   WM.PutRowIndex(6, iFirstIndex+2);
   WM.PutColIndex(6, iFirstIndex+2);
   
   doublereal p1 = pNode1->dGetX();
   doublereal p2 = pNode2->dGetX();
   doublereal p3 = pNode3->dGetX();
   doublereal p4 = pNode4->dGetX();
   doublereal density = HF->dGetDensity();

   s = XCurr(iFirstIndex+1); /* spostamento */
   v = XCurr(iFirstIndex+2); /* velocita' */
   
   doublereal jumpPres12 = fabs(p1-p2); /* salto di pressione nodo1 & nodo2 */
   doublereal jumpPres13 = fabs(p1-p3); /* salto di pressione nodo1 & nodo3 */
   doublereal jumpPres24 = fabs(p2-p4); /* salto di pressione nodo2 & nodo4 */
   doublereal jumpPres34 = fabs(p3-p4); /* salto di pressione nodo3 & nodo4 */
   
   if (jumpPres12 < 1.e8*std::numeric_limits<doublereal>::epsilon()) {
      jumpPres12 = 1.e8*std::numeric_limits<doublereal>::epsilon();
   }
   if (jumpPres13 < 1.e8*std::numeric_limits<doublereal>::epsilon()) {
      jumpPres13 = 1.e8*std::numeric_limits<doublereal>::epsilon();
   }
   if (jumpPres24 < 1.e8*std::numeric_limits<doublereal>::epsilon()) {
      jumpPres24 = 1.e8*std::numeric_limits<doublereal>::epsilon();
   }
   if (jumpPres34 < 1.e8*std::numeric_limits<doublereal>::epsilon()) {
      jumpPres34 = 1.e8*std::numeric_limits<doublereal>::epsilon();
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

#if 0
   doublereal Jac51  = -.2*costante*sp/sqrt(density*deltaP)-.43*width*s; 
   doublereal Jac52  = -Jac51;
   doublereal Jac53  = Jac51;
   doublereal Jac54  = -Jac51;
   doublereal Jac55  = -.4*valve_diameter*Cd*width*sqrt(density*deltaP)-dCoef*.43*width*deltaP-dCoef*c1-c2-dCoef*cf1-cf2;
   doublereal Jac56  = -Mass-c3-cf3;
#endif
   
   doublereal Jac51  = -.2*costante*sp/sqrt(density*deltaP)-.43*width*s; 
   doublereal Jac52  = 0.;
   doublereal Jac53  = 0.;
   //doublereal Jac54  = 0.;
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

   WM.PutCoef(1, 1, Jac11);
   WM.PutCoef(1, 2, Jac12);
   WM.PutCoef(1, 3, Jac13);
   WM.PutCoef(1, 4, Jac14);
   WM.PutCoef(1, 5, Jac15);
   WM.PutCoef(2, 1, Jac21);
   WM.PutCoef(2, 2, Jac22);
   WM.PutCoef(2, 3, Jac23);
   WM.PutCoef(2, 4, Jac24);
   WM.PutCoef(2, 5, Jac25);
   WM.PutCoef(3, 1, Jac31);
   WM.PutCoef(3, 2, Jac32);
   WM.PutCoef(3, 3, Jac33);
   WM.PutCoef(3, 4, Jac34);
   WM.PutCoef(3, 5, Jac35);
   WM.PutCoef(4, 1, Jac41);
   WM.PutCoef(4, 2, Jac42);
   WM.PutCoef(4, 3, Jac43);
   WM.PutCoef(4, 4, Jac44);
   WM.PutCoef(4, 5, Jac45);
   WM.PutCoef(5, 1, Jac51);
   WM.PutCoef(5, 2, Jac52);
   WM.PutCoef(5, 3, Jac53);
   WM.PutCoef(5, 5, Jac55);
   WM.PutCoef(5, 6, Jac56);
   WM.PutCoef(6, 5, Jac65);
   WM.PutCoef(6, 6, Jac66);
 
   return WorkMat;
}

SubVectorHandler& 
Dynamic_control_valve::AssRes(SubVectorHandler& WorkVec,
			  doublereal dCoef,
			  const VectorHandler& XCurr, 
			  const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering Dynamic_control_valve::AssRes()" << std::endl);
   
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
   
   s = XCurr(iFirstIndex+1);       /* spostamento */
   v = XCurr(iFirstIndex+2);       /* velocita' */
   sp = XPrimeCurr(iFirstIndex+1); /* velocita' */
   vp = XPrimeCurr(iFirstIndex+2); /* accelerazione */
   
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
   DEBUGCOUT("Force:         " << Force << std::endl);
   DEBUGCOUT("X equilibrio:  " << Force/K << std::endl);
   DEBUGCOUT("A1:            " << A1 << std::endl);
   DEBUGCOUT("A2:            " << A2  << std::endl);
   DEBUGCOUT("A3:            " << A3 << std::endl);
   DEBUGCOUT("A4:            " << A4 << std::endl);
   DEBUGCOUT("p1:            " << p1 << std::endl);
   DEBUGCOUT("p2:            " << p2 << std::endl);
   DEBUGCOUT("p3:            " << p3 << std::endl);
   DEBUGCOUT("p4:            " << p4 << std::endl);
   DEBUGCOUT("Cd:            " << Cd << std::endl);
   DEBUGCOUT("s_max :        " << s_max << std::endl);
   DEBUGCOUT("s :            " << s << std::endl);
   DEBUGCOUT("sp:            " << sp << std::endl);
   DEBUGCOUT("v :            " << v << std::endl);
   DEBUGCOUT("vp:            " << vp << std::endl);
   DEBUGCOUT("Valve_diameter:" << valve_diameter << std::endl);
   DEBUGCOUT("massa:         " << Mass << std::endl);
   DEBUGCOUT("smorzatore:    " << C << std::endl);
   DEBUGCOUT("molla:         " << K << std::endl);
   DEBUGCOUT("density:       " << density << std::endl);
   DEBUGCOUT("Valve_density: " << valve_density << std::endl);
   DEBUGCOUT("Area_max:      " << area_max << std::endl);
   DEBUGCOUT("Width:         " << width << std::endl);
   DEBUGCOUT("Loss_area:     " << loss_area << std::endl);
   DEBUGCOUT("c1:            " << c1 << std::endl);
   DEBUGCOUT("c2:            " << c2 << std::endl);
   DEBUGCOUT("c3:            " << c3 << std::endl);
   DEBUGCOUT("cf1:           " << cf1 << std::endl);
   DEBUGCOUT("cf2:           " << cf2 << std::endl);
   DEBUGCOUT("cf3:           " << cf3 << std::endl);
   DEBUGCOUT("Q12:           " << Q12 << std::endl);
   DEBUGCOUT("Q13:           " << Q13 << std::endl);
   DEBUGCOUT("Q24:           " << Q24 << std::endl);
   DEBUGCOUT("Q34:           " << Q34 << std::endl);
   DEBUGCOUT("PORTATE AI VARI NODI (positive se entranti)"<< std::endl);
   DEBUGCOUT("-Res_1 (portata nodo1): " << -Res_1 << std::endl); 
   DEBUGCOUT("-Res_2 (portata nodo2): " << -Res_2 << std::endl); 
   DEBUGCOUT("-Res_3 (portata nodo3): " << -Res_3 << std::endl); 
   DEBUGCOUT("-Res_4 (portata nodo4): " << -Res_4 << std::endl); 
   DEBUGCOUT("-Res_5  eq.dinamica   : " << -Res_5 << std::endl); 
   DEBUGCOUT("-Res_6 sp-v           : " << -Res_6 << std::endl); 
#endif
   
   WorkVec.PutItem(1, iNode1RowIndex, Res_1);
   WorkVec.PutItem(2, iNode2RowIndex, Res_2);         
   WorkVec.PutItem(3, iNode3RowIndex, Res_3);
   WorkVec.PutItem(4, iNode4RowIndex, Res_4);    
   WorkVec.PutItem(5, iFirstIndex+1, Res_5);
   WorkVec.PutItem(6, iFirstIndex+2, Res_6);    
    
   return WorkVec;
}


void Dynamic_control_valve::Output(OutputHandler& OH) const
{
   if (bToBeOutput()) { 
      std::ostream& out = OH.Hydraulic();
      out << std::setw(8) << GetLabel()
	<< " " << s << " "  << sp << " "  << vp 
	<< " " << flow1 << " " << flow2 << " " << flow3 << " " << flow4
	<< " " << A1 << " "  << A2 << " "  << A3 << " " << A4 << " " << std::endl;
   }
}


void 
Dynamic_control_valve::SetValue(DataManager *pDM,
		VectorHandler& X , VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
   integer i = iGetFirstIndex();
   
   X.PutCoef(i+1, start);
   X.PutCoef(i+2, 0.);
   XP.PutCoef(i+1, 0.);
   XP.PutCoef(i+2, 0.);
}

const OutputHandler::Dimensions 
Dynamic_control_valve::GetEquationDimension(integer index) const {
   // DOF == 2
   OutputHandler::Dimensions dimension = OutputHandler::Dimensions::UnknownDimension;

	switch (index)
	{
		case 1:
			dimension = OutputHandler::Dimensions::Force;
			break;
		case 2:
			dimension = OutputHandler::Dimensions::Velocity;
			break;
	}

	return dimension;
}

std::ostream&
Dynamic_control_valve::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{

	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << ": " <<
			"dynamic control valve force balance" << std::endl
      
      << prefix << iIndex + 2 << ": " <<
         "dynamic control valve velocity" << std::endl;

	return out;
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
: Elem(uL, fOut),
HydraulicElem(uL, pDO, hf, fOut), DriveOwner(pDC),
pNode1(p1), pNode2(p2), pNode3(p3), pNode4(p4), pNode5(p5), pNode6(p6),
start(s0), 
c_spost(cs), c_vel(cv), c_acc(ca),
width(W), 
loss_area(Loss_A),
valve_diameter(Valve_d),
valve_density(Valve_rho),
s_max(s_mx)
{
   ASSERT(pNode1 != NULL);
   ASSERT(pNode1->GetNodeType() == Node::HYDRAULIC);
   ASSERT(pNode2 != NULL);
   ASSERT(pNode2->GetNodeType() == Node::HYDRAULIC);
   ASSERT(pNode3 != NULL);
   ASSERT(pNode3->GetNodeType() == Node::HYDRAULIC);
   ASSERT(pNode4 != NULL);
   ASSERT(pNode4->GetNodeType() == Node::HYDRAULIC);
   ASSERT(pNode5 != NULL);
   ASSERT(pNode5->GetNodeType() == Node::HYDRAULIC);
   ASSERT(pNode6 != NULL);
   ASSERT(pNode6->GetNodeType() == Node::HYDRAULIC);
   ASSERT(s0 >= 0.);
   ASSERT(Valve_rho > std::numeric_limits<doublereal>::epsilon());
   ASSERT(Valve_d > std::numeric_limits<doublereal>::epsilon());
   ASSERT(W > std::numeric_limits<doublereal>::epsilon());
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
HydraulicElem::Type Pressure_flow_control_valve::GetHydraulicType(void) const 
{
   return HydraulicElem::PRESSURE_FLOW_CONTROL_VALVE;
}

/* Contributo al file di restart */
std::ostream& Pressure_flow_control_valve::Restart(std::ostream& out) const
{
   return out << "Pressure_flow_control_valve not implemented yet!" << std::endl;
}
   
unsigned int Pressure_flow_control_valve::iGetNumDof(void) const 
{
   return 2;
}
   
DofOrder::Order Pressure_flow_control_valve::GetDofType(unsigned int i) const 
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
   DEBUGCOUT("Entering Control_valve::AssJac()" << std::endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.ResizeReset(8, 8);
   
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
   
   WM.PutRowIndex(1, iNode1RowIndex);
   WM.PutRowIndex(2, iNode2RowIndex);
   WM.PutRowIndex(3, iNode3RowIndex);
   WM.PutRowIndex(4, iNode4RowIndex);
   WM.PutRowIndex(5, iNode5RowIndex);
   WM.PutRowIndex(6, iNode6RowIndex);
   WM.PutColIndex(1, iNode1ColIndex);
   WM.PutColIndex(2, iNode2ColIndex);
   WM.PutColIndex(3, iNode3ColIndex);
   WM.PutColIndex(4, iNode4ColIndex);
   WM.PutColIndex(5, iNode5ColIndex);
   WM.PutColIndex(6, iNode6ColIndex);
 
   WM.PutRowIndex(7, iFirstIndex+1);
   WM.PutColIndex(7, iFirstIndex+1);
   WM.PutRowIndex(8, iFirstIndex+2);
   WM.PutColIndex(8, iFirstIndex+2);
   
   doublereal p1 = pNode1->dGetX();
   doublereal p2 = pNode2->dGetX();
   doublereal p3 = pNode3->dGetX();
   doublereal p4 = pNode4->dGetX();
   //doublereal p5 = pNode5->dGetX();
   //doublereal p6 = pNode6->dGetX();
   doublereal density = HF->dGetDensity();

   s = XCurr(iFirstIndex+1); /* spostamento */
   v = XCurr(iFirstIndex+2); /* velocita' */
   
   doublereal jumpPres12 = fabs(p1-p2); /* salto di pressione nodo1 & nodo2 */
   doublereal jumpPres13 = fabs(p1-p3); /* salto di pressione nodo1 & nodo3 */
   doublereal jumpPres24 = fabs(p2-p4); /* salto di pressione nodo2 & nodo4 */
   doublereal jumpPres34 = fabs(p3-p4); /* salto di pressione nodo3 & nodo4 */
   
   if (jumpPres12 < 1.e8*std::numeric_limits<doublereal>::epsilon()) {
      jumpPres12 = 1.e8*std::numeric_limits<doublereal>::epsilon();
   }
   if (jumpPres13 < 1.e8*std::numeric_limits<doublereal>::epsilon()) {
      jumpPres13 = 1.e8*std::numeric_limits<doublereal>::epsilon();
   }
   if (jumpPres24 < 1.e8*std::numeric_limits<doublereal>::epsilon()) {
      jumpPres24 = 1.e8*std::numeric_limits<doublereal>::epsilon();
   }
   if (jumpPres34 < 1.e8*std::numeric_limits<doublereal>::epsilon()) {
      jumpPres34 = 1.e8*std::numeric_limits<doublereal>::epsilon();
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
   //doublereal Jac74  = 0.;
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

   WM.PutCoef(1, 1, Jac11);
   WM.PutCoef(1, 2, Jac12);
   WM.PutCoef(1, 3, Jac13);
   WM.PutCoef(1, 4, Jac14);
   WM.PutCoef(1, 7, Jac17);
   WM.PutCoef(2, 1, Jac21);
   WM.PutCoef(2, 2, Jac22);
   WM.PutCoef(2, 3, Jac23);
   WM.PutCoef(2, 4, Jac24);
   WM.PutCoef(2, 7, Jac27);
   WM.PutCoef(3, 1, Jac31);
   WM.PutCoef(3, 2, Jac32);
   WM.PutCoef(3, 3, Jac33);
   WM.PutCoef(3, 4, Jac34);
   WM.PutCoef(3, 7, Jac37);
   WM.PutCoef(4, 1, Jac41);
   WM.PutCoef(4, 2, Jac42);
   WM.PutCoef(4, 3, Jac43);
   WM.PutCoef(4, 4, Jac44);
   WM.PutCoef(4, 7, Jac47);
   WM.PutCoef(5, 7, Jac57); // (nearly) blind fix; was Jac67
   WM.PutCoef(6, 7, Jac67);
   WM.PutCoef(7, 1, Jac71);
   WM.PutCoef(7, 2, Jac72);
   WM.PutCoef(7, 3, Jac73);
   WM.PutCoef(7, 7, Jac77);
   WM.PutCoef(7, 8, Jac78);
   WM.PutCoef(8, 7, Jac87);
   WM.PutCoef(8, 8, Jac88);
 
   return WorkMat;
}

SubVectorHandler& 
Pressure_flow_control_valve::AssRes(SubVectorHandler& WorkVec,
			  doublereal dCoef,
			  const VectorHandler& XCurr, 
			  const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering Pressure_flow_control_valve::AssRes()" << std::endl);
   
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
   //doublereal p5 = pNode5->dGetX();
   //doublereal p6 = pNode6->dGetX();
   doublereal density = HF->dGetDensity();

   Force = pGetDriveCaller()->dGet();
   
   s = XCurr(iFirstIndex+1);       /* spostamento */
   v = XCurr(iFirstIndex+2);       /* velocita' */
   sp = XPrimeCurr(iFirstIndex+1); /* velocita' */
   vp = XPrimeCurr(iFirstIndex+2); /* accelerazione */
   
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
//#warning "????????????? Res_6 = -Res_6 ?"
//   doublereal Res_6 = -Res_6;
   doublereal Res_6 = -Res_5; // (nearly) blind fix
 
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
   DEBUGCOUT("Force:         " << Force << std::endl);
   DEBUGCOUT("X equilibrio:  " << Force/K << std::endl);
   DEBUGCOUT("A1:            " << A1 << std::endl);
   DEBUGCOUT("A2:            " << A2  << std::endl);
   DEBUGCOUT("A3:            " << A3 << std::endl);
   DEBUGCOUT("A4:            " << A4 << std::endl);
   DEBUGCOUT("p1:            " << p1 << std::endl);
   DEBUGCOUT("p2:            " << p2 << std::endl);
   DEBUGCOUT("p3:            " << p3 << std::endl);
   DEBUGCOUT("p4:            " << p4 << std::endl);
   DEBUGCOUT("Cd:            " << Cd << std::endl);
   DEBUGCOUT("s_max :        " << s_max << std::endl);
   DEBUGCOUT("x :            " << x << std::endl);
   DEBUGCOUT("s :            " << s << std::endl);
   DEBUGCOUT("sp:            " << sp << std::endl);
   DEBUGCOUT("v :            " << v << std::endl);
   DEBUGCOUT("vp:            " << vp << std::endl);
   DEBUGCOUT("Valve_diameter:" << valve_diameter << std::endl);
   DEBUGCOUT("massa:         " << Mass << std::endl);
   DEBUGCOUT("smorzatore:    " << C << std::endl);
   DEBUGCOUT("molla:         " << K << std::endl);
   DEBUGCOUT("density:       " << density << std::endl);
   DEBUGCOUT("Valve_density: " << valve_density << std::endl);
   DEBUGCOUT("Area_max:      " << area_max << std::endl);
   DEBUGCOUT("Width:         " << width << std::endl);
   DEBUGCOUT("Loss_area:     " << loss_area << std::endl);
   DEBUGCOUT("c1:            " << c1 << std::endl);
   DEBUGCOUT("c2:            " << c2 << std::endl);
   DEBUGCOUT("c3:            " << c3 << std::endl);
   DEBUGCOUT("cf1:           " << cf1 << std::endl);
   DEBUGCOUT("cf2:           " << cf2 << std::endl);
   DEBUGCOUT("cf3:           " << cf3 << std::endl);
   DEBUGCOUT("Q12:           " << Q12 << std::endl);
   DEBUGCOUT("Q13:           " << Q13 << std::endl);
   DEBUGCOUT("Q24:           " << Q24 << std::endl);
   DEBUGCOUT("Q34:           " << Q34 << std::endl);
   DEBUGCOUT("PORTATE AI VARI NODI (positive se entranti)"<< std::endl);
   DEBUGCOUT("-Res_1 (portata nodo1): " << -Res_1 << std::endl); 
   DEBUGCOUT("-Res_2 (portata nodo2): " << -Res_2 << std::endl); 
   DEBUGCOUT("-Res_3 (portata nodo3): " << -Res_3 << std::endl); 
   DEBUGCOUT("-Res_4 (portata nodo4): " << -Res_4 << std::endl); 
   DEBUGCOUT("-Res_5 (portata nodo3): " << -Res_5 << std::endl); 
   DEBUGCOUT("-Res_6 (portata nodo4): " << -Res_6 << std::endl); 
   DEBUGCOUT("-Res_6   eq dinamica  : " << -Res_7 << std::endl); 
   DEBUGCOUT("-Res_7 sp-v           : " << -Res_8 << std::endl); 
#endif
   
   WorkVec.PutItem(1, iNode1RowIndex, Res_1);
   WorkVec.PutItem(2, iNode2RowIndex, Res_2);         
   WorkVec.PutItem(3, iNode3RowIndex, Res_3);
   WorkVec.PutItem(4, iNode4RowIndex, Res_4);    
   WorkVec.PutItem(5, iNode5RowIndex, Res_5);
   WorkVec.PutItem(6, iNode6RowIndex, Res_6);    
   WorkVec.PutItem(7, iFirstIndex+1, Res_7);
   WorkVec.PutItem(8, iFirstIndex+2, Res_8);    
    
   return WorkVec;
}


void Pressure_flow_control_valve::Output(OutputHandler& OH) const
{
   if (bToBeOutput()) { 
      std::ostream& out = OH.Hydraulic();
      out << std::setw(8) << GetLabel()
	<< " " << s << " "  << sp << " "  << vp 
	<< " " << flow1 << " " << flow2 << " " << flow3 << " " << flow4
	<< " " << flow5 << " " << flow6
	<< " " << A1 << " "  << A2 << " "  << A3 << " " << A4 << " " << std::endl;
   }
}


void 
Pressure_flow_control_valve::SetValue(DataManager *pDM,
		VectorHandler& X , VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
   integer i = iGetFirstIndex();
   
   X.PutCoef(i+1, start);
   X.PutCoef(i+2, 0.);
   XP.PutCoef(i+1, 0.);
   XP.PutCoef(i+2, 0.);
}

const OutputHandler::Dimensions 
Pressure_flow_control_valve::GetEquationDimension(integer index) const {
   // DOF == 2
   OutputHandler::Dimensions dimension = OutputHandler::Dimensions::UnknownDimension;

	switch (index)
	{
		case 1:
			dimension = OutputHandler::Dimensions::Force;
			break;
		case 2:
			dimension = OutputHandler::Dimensions::Velocity;
			break;
	}

	return dimension;
}

std::ostream&
Pressure_flow_control_valve::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{

	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << ": " <<
			"pressure flow valve force balance" << std::endl
      
      << prefix << iIndex + 2 << ": " <<
         "pressure flow valve velocity" << std::endl;

	return out;
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
: Elem(uL, fOut),
HydraulicElem(uL, pDO, hf, fOut),
pNode1(p1), pNode2(p2), area_diaf(A_dia), mass(mv),
area_max(A_max),
Kappa(K), force0(F0), width(w),
s_max(s_mx),
c_spost(cs), 
c_vel(cv), 
c_acc(ca)
{
   ASSERT(pNode1 != NULL);
   ASSERT(pNode1->GetNodeType() == Node::HYDRAULIC);
   ASSERT(pNode2 != NULL);
   ASSERT(pNode2->GetNodeType() == Node::HYDRAULIC);
   ASSERT(A_dia > std::numeric_limits<doublereal>::epsilon());
   ASSERT(mv > std::numeric_limits<doublereal>::epsilon()); 
   ASSERT(A_max > std::numeric_limits<doublereal>::epsilon());
   ASSERT(K >= 0. );
   ASSERT(F0 >= 0.);
   ASSERT(w > std::numeric_limits<doublereal>::epsilon()); 
   ASSERT(s_mx>=0.);  // corsa massima della valvola

}

Pressure_valve::~Pressure_valve(void)
{
   NO_OP;
}
   
/* Tipo di elemento idraulico (usato solo per debug ecc.) */
HydraulicElem::Type Pressure_valve::GetHydraulicType(void) const 
{
   return HydraulicElem::PRESSURE_VALVE;
}


/* Contributo al file di restart */
std::ostream& Pressure_valve::Restart(std::ostream& out) const
{
   return out << "Pressure_valve not implemented yet!" << std::endl;
}


unsigned int Pressure_valve::iGetNumDof(void) const 
{
  return 2; 
}


DofOrder::Order Pressure_valve::GetDofType(unsigned int i) const 
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
   DEBUGCOUT("Entering Pressure_valve::AssJac()" << std::endl);
      
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.Resize(4, 4);
   
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

   doublereal p1 = pNode1->dGetX();
   doublereal p2 = pNode2->dGetX(); 
   /* p1 = XCurr(pNode1->iGetFirstRowIndex()+1); */
   doublereal density = HF->dGetDensity();
  
   s = XCurr(iFirstIndex+1); /* spostamento */
   v = XCurr(iFirstIndex+2); /* velocita' */
   
   doublereal Cd = .6;       
   doublereal jumpPres = fabs(p1-p2);
 
   /* evito di dividere per un numero troppo piccolo */
   if (jumpPres < 1.e8*std::numeric_limits<doublereal>::epsilon()) {
      jumpPres = 1.e8*std::numeric_limits<doublereal>::epsilon();
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
   DEBUGCOUT("Valore di c1: " << c1 << std::endl);
   DEBUGCOUT("Valore di c2: " << c2 << std::endl);
   DEBUGCOUT("Valore di c3: " << c3 << std::endl);
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
   DEBUGCOUT("Jac smorzatore " << density*area_max*pow(area_max/(area_diaf*Cd), 2) << std::endl);
#endif
   
   WM.PutCoef(1, 1, Jac11);
   WM.PutCoef(1, 2, Jac12);
   WM.PutCoef(1, 3, Jac13);
   WM.PutCoef(1, 4, Jac14);
   WM.PutCoef(2, 1, Jac21);
   WM.PutCoef(2, 2, Jac22);
   WM.PutCoef(2, 3, Jac23);
   WM.PutCoef(2, 4, Jac24);
   WM.PutCoef(3, 1, Jac31);
   WM.PutCoef(3, 2, Jac32);
   WM.PutCoef(3, 3, Jac33);
   WM.PutCoef(3, 4, Jac34);
   WM.PutCoef(4, 1, Jac41);
   WM.PutCoef(4, 2, Jac42);
   WM.PutCoef(4, 3, Jac43);
   WM.PutCoef(4, 4, Jac44);

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
   
   s = XCurr(iFirstIndex+1);       /* spostamento */
   v = XCurr(iFirstIndex+2);       /* velocita' */
   sp = XPrimeCurr(iFirstIndex+1); /* velocita' */
   vp = XPrimeCurr(iFirstIndex+2); /* accelerazione */
   doublereal density = HF->dGetDensity();
  
   doublereal Cd = .6;
   doublereal jumpPres = fabs(p1-p2);
   doublereal Res_1 = 0.;
   doublereal Res_2 = 0.;
   doublereal Res_3 = 0.;
   doublereal Res_4 = 0.;
   doublereal x0 = (p1*area_max-force0)/Kappa; /* punto di equilibrio */
   
#ifdef HYDR_DEVEL
   DEBUGCOUT("s:  " << s << std::endl);
   DEBUGCOUT("sp: " << sp << std::endl);
   DEBUGCOUT("v:  " << v << std::endl);
   DEBUGCOUT("vp: " << vp << std::endl);
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
   DEBUGCOUT("width:    " << width << std::endl);
   DEBUGCOUT("Cd:       " << Cd << std::endl);
   DEBUGCOUT("jumpPres: " << jumpPres  << std::endl);
   DEBUGCOUT("density:  " << density << std::endl);
   DEBUGCOUT("p1:       " << p1 << std::endl);
   DEBUGCOUT("p2:       " << p2 << std::endl);
   DEBUGCOUT("s:        " << s << std::endl);
   DEBUGCOUT("sp:       " << sp << std::endl);
   DEBUGCOUT("v:        " << v << std::endl);
   DEBUGCOUT("vp:       " << vp << std::endl);
   DEBUGCOUT("area_max: " << area_max << std::endl);
   DEBUGCOUT("Kappa:    " << Kappa << std::endl);
   DEBUGCOUT("force0:   " << force0 << std::endl);
   DEBUGCOUT("mass:     " << mass  << std::endl);
   DEBUGCOUT("x0:       " << x0 << std::endl);
   DEBUGCOUT("s_max:    " << s_max << std::endl);
   DEBUGCOUT("Valore dello smorzatore: " 
	     << copysign(.5*density*area_max*pow(area_max/(area_diaf*Cd), 2), sp) << std::endl);
   DEBUGCOUT("c1:       " << c1 << std::endl);
   DEBUGCOUT("c2:       " << c2 << std::endl);
   DEBUGCOUT("c3:       " << c3 << std::endl);
   DEBUGCOUT("c4:       " << c4 << std::endl);
   DEBUGCOUT("cf1:      " << cf1 << std::endl);
   DEBUGCOUT("cf2:      " << cf2 << std::endl);
   DEBUGCOUT("cf3:      " << cf3 << std::endl);
   DEBUGCOUT("cf4:      " << cf4 << std::endl);
   DEBUGCOUT("PORTATE AI VARI NODI (positive se entranti)" << std::endl);
   DEBUGCOUT("-Res_1 (portata nodo1): " << -Res_1 << std::endl); 
   DEBUGCOUT("-Res_2 (portata nodo2): " << -Res_2 << std::endl);   
   DEBUGCOUT("Res_3:                  " << Res_3 << std::endl); 
   DEBUGCOUT("Res_4:                  " << Res_4 << std::endl); 
#else
   // silence set but not used for variables used only with HYDR_DEVEL
   (void)x0;
#endif
   
   WorkVec.PutItem(1, iNode1RowIndex, Res_1);
   WorkVec.PutItem(2, iNode2RowIndex, Res_2);         
   WorkVec.PutItem(3, iFirstIndex+1, Res_3);
   WorkVec.PutItem(4, iFirstIndex+2, Res_4);    
     	
   return WorkVec;
}


void Pressure_valve::Output(OutputHandler& OH) const
{
   if (bToBeOutput()) { 
      std::ostream& out = OH.Hydraulic();
      out << std::setw(8) << GetLabel()
	<< " " << s << " " << v  << " "<< vp  
	<< " " << flow1  << " " << flow2 << std::endl;
   }  
}


void Pressure_valve::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
   integer i = iGetFirstIndex();
   
   X.PutCoef(i+1, 0.);
   X.PutCoef(i+2, 0.);
   XP.PutCoef(i+1, 0.);
   XP.PutCoef(i+2, 0.);
}

const OutputHandler::Dimensions 
Pressure_valve::GetEquationDimension(integer index) const {
   // DOF == 2
   OutputHandler::Dimensions dimension = OutputHandler::Dimensions::UnknownDimension;

	switch (index)
	{
		case 1:
			dimension = OutputHandler::Dimensions::Force;
			break;
		case 2:
			dimension = OutputHandler::Dimensions::Velocity;
			break;
	}

	return dimension;
}

std::ostream&
Pressure_valve::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{

	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << ": " <<
			"pressure valve force balance" << std::endl
      
      << prefix << iIndex + 2 << ": " <<
         "pressure valve velocity" << std::endl;

	return out;
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
: Elem(uL, fOut),
HydraulicElem(uL, pDO, hf, fOut),
pNode1(p1), pNode2(p2), pNode3(p3),
area_diaf(A_dia), mass(mv), area_pipe(A_pipe), area_max(A_max), 
Kappa(K), force0(F0), width(w), s_max(s_mx), c_spost(cs), c_vel(cv), c_acc(ca) 
{
   ASSERT(pNode1 != NULL);
   ASSERT(pNode1->GetNodeType() == Node::HYDRAULIC);
   ASSERT(pNode2 != NULL);
   ASSERT(pNode2->GetNodeType() == Node::HYDRAULIC);
   ASSERT(pNode3 != NULL);
   ASSERT(pNode3->GetNodeType() == Node::HYDRAULIC);
   ASSERT(A_dia > std::numeric_limits<doublereal>::epsilon());
   ASSERT(mv > std::numeric_limits<doublereal>::epsilon()); 
   ASSERT(A_pipe > std::numeric_limits<doublereal>::epsilon());
   ASSERT(A_max > std::numeric_limits<doublereal>::epsilon());
   ASSERT(K > std::numeric_limits<doublereal>::epsilon() );
   ASSERT(F0 >= 0.);
   ASSERT(w > std::numeric_limits<doublereal>::epsilon()); 
   ASSERT(s_max >= 0.);
   
   h = .02; /* coefficiente di perdita di carico concentrata tra i nodi 1 e 2 (smorza il moto della valvola) */
}


Flow_valve::~Flow_valve(void)
{
   NO_OP;
}
 

/* Tipo di elemento idraulico (usato solo per debug ecc.) */
HydraulicElem::Type  Flow_valve::GetHydraulicType(void) const 
{
   return HydraulicElem:: FLOW_VALVE;
}


/* Contributo al file di restart */
std::ostream&  Flow_valve::Restart(std::ostream& out) const
{
   return out << " Flow_valve not implemented yet!" << std::endl;
}
 

unsigned int  Flow_valve::iGetNumDof(void) const 
{
   return 2;
}


DofOrder::Order Flow_valve::GetDofType(unsigned int i) const 
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
   DEBUGCOUT("Entering  Flow_valve::AssJac()" << std::endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.Resize(5, 5);
   
   integer iNode1RowIndex = pNode1->iGetFirstRowIndex()+1;
   integer iNode2RowIndex = pNode2->iGetFirstRowIndex()+1;
   integer iNode3RowIndex = pNode3->iGetFirstRowIndex()+1;
   integer iNode1ColIndex = pNode1->iGetFirstColIndex()+1;
   integer iNode2ColIndex = pNode2->iGetFirstColIndex()+1;
   integer iNode3ColIndex = pNode3->iGetFirstColIndex()+1;
   integer iFirstIndex = iGetFirstIndex();
   
   WM.PutRowIndex(1, iNode1RowIndex);
   WM.PutRowIndex(2, iNode2RowIndex);
   WM.PutRowIndex(3, iNode3RowIndex);
   WM.PutColIndex(1, iNode1ColIndex);
   WM.PutColIndex(2, iNode2ColIndex);
   WM.PutColIndex(3, iNode3ColIndex);

   WM.PutRowIndex(4, iFirstIndex+1);
   WM.PutColIndex(4, iFirstIndex+1);
   WM.PutRowIndex(5, iFirstIndex+2);
   WM.PutColIndex(5, iFirstIndex+2);

   doublereal p1 = pNode1->dGetX();
   doublereal p2 = pNode2->dGetX();
   doublereal p3 = pNode3->dGetX();
   s = XCurr(iFirstIndex+1);  /* spostamento */
   v = XCurr(iFirstIndex+2);  /* velocita' */
   
   doublereal Cd = .6;        /* verificare */
   doublereal jumpPres12 = fabs(p1-p2);
   doublereal jumpPres23 = fabs(p2-p3);
   doublereal jumpPres13 = fabs(p1-p3);
   doublereal density = HF->dGetDensity();
    
   /* evito di dividere per un numero troppo piccolo */
   if (jumpPres12 < 1.e8*std::numeric_limits<doublereal>::epsilon()) {
      jumpPres12 = 1.e8*std::numeric_limits<doublereal>::epsilon();
   }
   if (jumpPres23 < 1.e8*std::numeric_limits<doublereal>::epsilon()) {
      jumpPres23 = 1.e8*std::numeric_limits<doublereal>::epsilon();
   }
   if (jumpPres13 < 1.e8*std::numeric_limits<doublereal>::epsilon()) {
      jumpPres13 = 1.e8*std::numeric_limits<doublereal>::epsilon();
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
   DEBUGCOUT("Valore di p1: " << p1 << std::endl);
   DEBUGCOUT("Valore di p2: " << p2 << std::endl);
   DEBUGCOUT("Valore di p3: " << p3 << std::endl);
   DEBUGCOUT("Valore di s: " << s << std::endl);
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
	     << density*area_max*pow(area_max/(area_diaf*Cd), 2) << std::endl);
   DEBUGCOUT("Jac smorzatore " << h*density << std::endl);
   DEBUGCOUT("Jac smorzatore Jac45old1 " << Jac45old1 << std::endl);
   DEBUGCOUT("Jac smorzatore Jac45old2 " << Jac45old2 << std::endl);
   DEBUGCOUT("Jac smorzatore Jac45     " << Jac45 << std::endl);
   
   DEBUGCOUT("Jac11: " << Jac11 << std::endl);
   DEBUGCOUT("Jac12: " << Jac12 << std::endl);
   DEBUGCOUT("Jac13: " << Jac13 << std::endl);
   DEBUGCOUT("Jac14: " << Jac14 << std::endl);
   DEBUGCOUT("Jac15: " << Jac15 << std::endl);
   DEBUGCOUT("Jac21: " << Jac21 << std::endl);
   DEBUGCOUT("Jac22: " << Jac22 << std::endl);
   DEBUGCOUT("Jac23: " << Jac23 << std::endl);
   DEBUGCOUT("Jac24: " << Jac24 << std::endl);
   DEBUGCOUT("Jac25: " << Jac25 << std::endl);
   DEBUGCOUT("Jac31: " << Jac31 << std::endl);
   DEBUGCOUT("Jac32: " << Jac32 << std::endl);
   DEBUGCOUT("Jac33: " << Jac33 << std::endl);
   DEBUGCOUT("Jac34: " << Jac34 << std::endl);
   DEBUGCOUT("Jac35: " << Jac35 << std::endl);
   DEBUGCOUT("Jac41: " << Jac41 << std::endl);
   DEBUGCOUT("Jac42: " << Jac42 << std::endl);
   DEBUGCOUT("Jac43: " << Jac43 << std::endl);
   DEBUGCOUT("Jac44: " << Jac44 << std::endl);
   DEBUGCOUT("Jac45: " << Jac45 << std::endl);
   DEBUGCOUT("Jac51: " << Jac51 << std::endl);
   DEBUGCOUT("Jac52: " << Jac52 << std::endl);
   DEBUGCOUT("Jac53: " << Jac53 << std::endl);
   DEBUGCOUT("Jac54: " << Jac54 << std::endl);
   DEBUGCOUT("Jac55: " << Jac55 << std::endl);
#else
   // silence set but not used warning for variables used only when HYDR_DEVEL
   (void)Jac45old1;
   (void)Jac45old2;
#endif
   
   WM.PutCoef(1, 1, Jac11);
   WM.PutCoef(1, 2, Jac12);
   WM.PutCoef(1, 3, Jac13);
   WM.PutCoef(1, 4, Jac14);
   WM.PutCoef(1, 5, Jac15);
   WM.PutCoef(2, 1, Jac21);
   WM.PutCoef(2, 2, Jac22);
   WM.PutCoef(2, 3, Jac23);
   WM.PutCoef(2, 4, Jac24);
   WM.PutCoef(2, 5, Jac25);
   WM.PutCoef(3, 1, Jac31);
   WM.PutCoef(3, 2, Jac32);
   WM.PutCoef(3, 3, Jac33);
   WM.PutCoef(3, 4, Jac34);
   WM.PutCoef(3, 5, Jac35);
   WM.PutCoef(4, 1, Jac41);
   WM.PutCoef(4, 2, Jac42);
   WM.PutCoef(4, 3, Jac43);
   WM.PutCoef(4, 4, Jac44);
   WM.PutCoef(4, 5, Jac45);
   WM.PutCoef(5, 1, Jac51);
   WM.PutCoef(5, 2, Jac52);
   WM.PutCoef(5, 3, Jac53);
   WM.PutCoef(5, 4, Jac54);
   WM.PutCoef(5, 5, Jac55);
   
   return WorkMat;
}

SubVectorHandler& 
Flow_valve::AssRes(SubVectorHandler& WorkVec,
		   doublereal dCoef,
		   const VectorHandler& XCurr, 
		   const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering Flow_valve::AssRes()" << std::endl);
   
   WorkVec.Resize(5);
   
   integer iNode1RowIndex = pNode1->iGetFirstRowIndex()+1;
   integer iNode2RowIndex = pNode2->iGetFirstRowIndex()+1;
   integer iNode3RowIndex = pNode3->iGetFirstRowIndex()+1;
   
   doublereal p1 = pNode1->dGetX();
   doublereal p2 = pNode2->dGetX();
   doublereal p3 = pNode3->dGetX();
   integer iFirstIndex = iGetFirstIndex();

   s = XCurr(iFirstIndex+1);        /* spostamento */
   v = XCurr(iFirstIndex+2);        /* velocita' */
   sp = XPrimeCurr(iFirstIndex+1);  /* velocita' */
   vp = XPrimeCurr(iFirstIndex+2);  /* accelerazione */
  
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
   DEBUGCOUT("Res_4:         " <<  Res_4 << std::endl);
   DEBUGCOUT("Res_4old:      " <<  Res_4old << std::endl);
   DEBUGCOUT("smorazatore:   " 
	     << copysign(.5*density*area_max*pow(sp*area_max/(area_diaf*Cd), 2), sp) << std::endl);
   DEBUGCOUT("smorzatoreold: " 
	     << copysign(.5*h*density*pow(sp*area_max/area_pipe, 2), sp) << std::endl);
   
   DEBUGCOUT("width:         " << width << std::endl);
   DEBUGCOUT("Cd:            " << Cd << std::endl);
   DEBUGCOUT("jumpPres12:    " << jumpPres12  << std::endl);
   DEBUGCOUT("jumpPres23:    " << jumpPres23  << std::endl);
   DEBUGCOUT("jumpPres13:    " << jumpPres13  << std::endl);
   DEBUGCOUT("density:       " << density << std::endl);
   DEBUGCOUT("p1:            " << p1 << std::endl);
   DEBUGCOUT("p2:            " << p2 << std::endl);
   DEBUGCOUT("p3:            " << p3 << std::endl);
   DEBUGCOUT("s:             " << s << std::endl);
   DEBUGCOUT("sp:            " << sp << std::endl);
   DEBUGCOUT("v:             " << v << std::endl);
   DEBUGCOUT("vp:            " << vp << std::endl);
   DEBUGCOUT("area_max:      " << area_max << std::endl);
   DEBUGCOUT("Kappa:         " << Kappa << std::endl);
   DEBUGCOUT("force0:        " << force0 << std::endl);
   DEBUGCOUT("mass:          " << mass << std::endl);
   DEBUGCOUT("s_max:         " << s_max << std::endl);
   DEBUGCOUT("x0:            " << x0  << std::endl);
   DEBUGCOUT("c1:            " << c1 << std::endl);
   DEBUGCOUT("c2:            " << c2 << std::endl);
   DEBUGCOUT("c3:            " << c3 << std::endl);
   DEBUGCOUT("c4:            " << c4 << std::endl);
   DEBUGCOUT("cf1:           " << cf1 << std::endl);
   DEBUGCOUT("cf2:           " << cf2 << std::endl);
   DEBUGCOUT("cf3:           " << cf3 << std::endl);
   DEBUGCOUT("cf4:           " << cf4 << std::endl);
   DEBUGCOUT("PORTATE AI VARI NODI (positive se entranti)" << std::endl);
   DEBUGCOUT("-Res_1 (portata nodo1): " << -Res_1 << std::endl); 
   DEBUGCOUT("-Res_2 (portata nodo2): " << -Res_2 << std::endl); 
   DEBUGCOUT("-Res_3:(portata nodo3): " << -Res_3 << std::endl); 
   DEBUGCOUT("-Res_4:                 " << -Res_4 << std::endl); 
   DEBUGCOUT("-Res_5:                 " << -Res_5 << std::endl); 
#else
   // silence set but not used warning for variables used only when HYDR_DEVEL
   (void)jumpPres23;
   (void)x0;
   (void)Res_4old;
#endif
   
   WorkVec.PutItem(1, iNode1RowIndex, Res_1);
   WorkVec.PutItem(2, iNode2RowIndex, Res_2);         
   WorkVec.PutItem(3, iNode3RowIndex, Res_3);
   WorkVec.PutItem(4, iFirstIndex+1, Res_4);    
   WorkVec.PutItem(5, iFirstIndex+2, Res_5);    
     
   return WorkVec;
}
  
void Flow_valve::Output(OutputHandler& OH) const
{
   if (bToBeOutput()) { 
      std::ostream& out = OH.Hydraulic();
      out << std::setw(8) << GetLabel()
	<< " " << s  << " " << v  << " "<< vp  
	<< " " << flow1  << " "<< flow2  << " "<< flow3  << std::endl;
   }  
}

void Flow_valve::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{   
   integer i = iGetFirstIndex();
   X.PutCoef(i+1, 0.);
   X.PutCoef(i+2, 0.);
   XP.PutCoef(i+1, 0.);
   XP.PutCoef(i+2, 0.);  
}

const OutputHandler::Dimensions 
Flow_valve::GetEquationDimension(integer index) const {
   // DOF == 2
   OutputHandler::Dimensions dimension = OutputHandler::Dimensions::UnknownDimension;

	switch (index)
	{
		case 1:
			dimension = OutputHandler::Dimensions::Force;
			break;
		case 2:
			dimension = OutputHandler::Dimensions::Velocity;
			break;
	}

	return dimension;
}

std::ostream&
Flow_valve::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{

	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << ": " <<
			"flow valve force balance" << std::endl
      
      << prefix << iIndex + 2 << ": " <<
         "flow valve velocity" << std::endl;

	return out;
}
/* Flow_valve - end */
