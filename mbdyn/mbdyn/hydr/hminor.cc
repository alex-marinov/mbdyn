/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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

#include "hminor.h"

/* Minor_losses - begin */

Minor_loss::Minor_loss(unsigned int uL, const DofOwner* pDO,
			   HydraulicFluid* hf,
			   const PressureNode* p1, const PressureNode* p2,
			   doublereal dK1, doublereal dK2, doublereal A,
			   flag fOut)
: Elem(uL, fOut),
HydraulicElem(uL, pDO, hf, fOut),
pNode1(p1), pNode2(p2),
dKappa1(dK1), dKappa2(dK2), area(A),
flow(0.),
vel(0.),
dKappa(0.)
{
   ASSERT(pNode1 != NULL);
   ASSERT(pNode1->GetNodeType() == Node::HYDRAULIC);
   ASSERT(pNode2 != NULL);
   ASSERT(pNode2->GetNodeType() == Node::HYDRAULIC);
   ASSERT(dK1 >= 0.);
   ASSERT(dK2 >= 0.);
   ASSERT(A > std::numeric_limits<doublereal>::epsilon());
}

Minor_loss::~Minor_loss(void)
{
   NO_OP;
}
   
/* Tipo di elemento idraulico (usato solo per debug ecc.) */
HydraulicElem::Type Minor_loss::GetHydraulicType(void) const {
   return HydraulicElem::MINOR_LOSS;
}

/* Contributo al file di restart */
std::ostream& Minor_loss::Restart(std::ostream& out) const
{
   return out << "Minor_loss not implemented yet!" << std::endl;
}
   
unsigned int Minor_loss::iGetNumDof(void) const { 
   return 0;
}
   
DofOrder::Order Minor_loss::GetDofType(unsigned int i) const {
   silent_cerr("Minor_loss has no dofs!" << std::endl);
   throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

void Minor_loss::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
   *piNumRows = 2; 
   *piNumCols = 2; 
}
      
VariableSubMatrixHandler& 
Minor_loss::AssJac(VariableSubMatrixHandler& WorkMat,
		  doublereal dCoef,
		  const VectorHandler& XCurr, 
		  const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering Minor_loss::AssJac()" << std::endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.Resize(2, 2);
   
   integer iNode1RowIndex = pNode1->iGetFirstRowIndex()+1;
   integer iNode2RowIndex = pNode2->iGetFirstRowIndex()+1;
   integer iNode1ColIndex = pNode1->iGetFirstColIndex()+1;
   integer iNode2ColIndex = pNode2->iGetFirstColIndex()+1;
   
   WM.PutRowIndex(1, iNode1RowIndex);
   WM.PutRowIndex(2, iNode2RowIndex);
   WM.PutColIndex(1, iNode1ColIndex);
   WM.PutColIndex(2, iNode2ColIndex);
   
   doublereal p1 = pNode1->dGetX();
   doublereal p2 = pNode2->dGetX();
   
   doublereal jumpPres = fabs(p1-p2);
   
   /* evito di dividere per un numero troppo piccolo */
   if (jumpPres < 1.e8*std::numeric_limits<doublereal>::epsilon()) {
      jumpPres = 1.e8*std::numeric_limits<doublereal>::epsilon();
   }
   /*
    * se voglio usare un fluido comprimibile, metto la pressione
    * media nel condotto:
   */
   
   doublereal density = HF->dGetDensity((p1+p2)/2.);
    
   
   /* altrimenti lascio la densita' di riferimento 
   doublereal density = HF->dGetDensity();
   */
   doublereal Jac = -density*.5*area*sqrt(2./(dKappa*density*jumpPres));
    
   WM.PutCoef(1, 1, Jac);
   WM.PutCoef(1, 2, -Jac);
   WM.PutCoef(2, 1, -Jac);
   WM.PutCoef(2, 2, Jac);
   
   return WorkMat;
}

SubVectorHandler& 
Minor_loss::AssRes(SubVectorHandler& WorkVec,
		     doublereal dCoef,
		     const VectorHandler& XCurr, 
		     const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering Minor_loss::AssRes()" << std::endl);
   
   WorkVec.Resize(2);
   
   integer iNode1RowIndex = pNode1->iGetFirstRowIndex()+1;
   integer iNode2RowIndex = pNode2->iGetFirstRowIndex()+1;
   
   doublereal p1 = pNode1->dGetX();
   doublereal p2 = pNode2->dGetX();
   
   doublereal jumpPres = fabs(p1-p2); 
  
   if (p1 > p2) {
      dKappa = dKappa1;  /* flusso diretto da 1 a 2 */
   } else { 
      dKappa = dKappa2;  /* flusso diretto da 2 a 1 */
   }

   doublereal density = HF->dGetDensity((p1+p2)/2.);
   flow = density*area*sqrt(2./(dKappa*density))*copysign(sqrt(jumpPres), p1-p2);
   vel= flow/(density*area);

#ifdef HYDR_DEVEL
   DEBUGCOUT("RES area :           " << area << std::endl);
   DEBUGCOUT("RES flow:            " << flow << std::endl);
   DEBUGCOUT("RES p1:              " << p1 << std::endl);
   DEBUGCOUT("RES p2:              " << p2 << std::endl);
   DEBUGCOUT("RES dKappa:          " << dKappa << std::endl);
   DEBUGCOUT("****************************************************" << std::endl);
   DEBUGCOUT("RES velocita':       " << vel << std::endl);
   DEBUGCOUT("    se positiva il fluido va dal nodo 1 al nodo 2 " << std::endl);
   DEBUGCOUT("RES portata (nodo2): " << flow << std::endl);
   DEBUGCOUT("****************************************************" << std::endl);
#endif
   WorkVec.PutItem(1, iNode1RowIndex, flow);
   WorkVec.PutItem(2, iNode2RowIndex, -flow);         

   return WorkVec;
}
   
void Minor_loss::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) { 
      std::ostream& out = OH.Hydraulic();
      out << std::setw(8) << GetLabel() 
          << " " << vel  << " " << flow << std::endl;
   }
}

/* Minor_loss - end */


/* ThreeWayMinorLoss - begin */

ThreeWayMinorLoss::ThreeWayMinorLoss(
		unsigned int uL, const DofOwner* pDO,
		HydraulicFluid* hf, const PressureNode* p0,
		const PressureNode* p1, const PressureNode* p2,
		doublereal dK1, doublereal dK2, 
		doublereal A1, doublereal A2, flag fOut)
: Elem(uL, fOut),
HydraulicElem(uL, pDO, hf, fOut),
pNode0(p0), pNode1(p1), pNode2(p2), pNodeN(NULL),
dKappa1(dK1), dKappa2(dK2), area1(A1), area2(A2),
area(0.),
flow(0.),
vel(0.),
dKappa(0.)
{
	ASSERT(pNode0 != NULL);
	ASSERT(pNode0->GetNodeType() == Node::HYDRAULIC);
	ASSERT(pNode1 != NULL);
	ASSERT(pNode1->GetNodeType() == Node::HYDRAULIC);
	ASSERT(pNode2 != NULL);
	ASSERT(pNode2->GetNodeType() == Node::HYDRAULIC);
	ASSERT(dK1 >= 0.);
	ASSERT(dK2 >= 0.);
	ASSERT(A1 > std::numeric_limits<doublereal>::epsilon());
	ASSERT(A2 > std::numeric_limits<doublereal>::epsilon());
}

ThreeWayMinorLoss::~ThreeWayMinorLoss(void)
{
	NO_OP;
}
   
/* Tipo di elemento idraulico (usato solo per debug ecc.) */
HydraulicElem::Type ThreeWayMinorLoss::GetHydraulicType(void) const {
	return HydraulicElem::THREEWAYMINORLOSS;
}

/* Contributo al file di restart */
std::ostream& ThreeWayMinorLoss::Restart(std::ostream& out) const
{
	return out << "ThreeWayMinorLoss not implemented yet!" << std::endl;
}
 
unsigned int ThreeWayMinorLoss::iGetNumDof(void) const { 
	return 0;
}
   
DofOrder::Order ThreeWayMinorLoss::GetDofType(unsigned int i) const {
	silent_cerr("Minor_loss has no dofs!" << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

void
ThreeWayMinorLoss::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 2;
	*piNumCols = 2;
}
      
VariableSubMatrixHandler& 
ThreeWayMinorLoss::AssJac(
		VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr
		)
{
	DEBUGCOUT("Entering Minor_loss::AssJac()" << std::endl);
   
	ASSERT(pNodeN != NULL);
	
	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.Resize(2, 2);
	
	integer iNode0RowIndex = pNode0->iGetFirstRowIndex()+1;
	integer iNodeNRowIndex = pNodeN->iGetFirstRowIndex()+1;
	integer iNode0ColIndex = pNode0->iGetFirstColIndex()+1;
	integer iNodeNColIndex = pNodeN->iGetFirstColIndex()+1;
	
	WM.PutRowIndex(1, iNode0RowIndex);
	WM.PutRowIndex(2, iNodeNRowIndex);
	WM.PutColIndex(1, iNode0ColIndex);
	WM.PutColIndex(2, iNodeNColIndex);

	doublereal p0 = pNode0->dGetX();
	doublereal p = pNodeN->dGetX();
	
	doublereal jumpPres = fabs(p0-p);
	
	/* evito di dividere per un numero troppo piccolo */
	if (jumpPres < 1.e8*std::numeric_limits<doublereal>::epsilon()) {
		jumpPres = 1.e8*std::numeric_limits<doublereal>::epsilon();
	}
	/*
	 * se voglio usare un fluido comprimibile, metto la pressione
	 * media nel condotto:
         */

	doublereal density = HF->dGetDensity((p0+p)/2.);
	
	/* 
	 * altrimenti lascio la densita' di riferimento 
	 * doublereal density = HF->dGetDensity();
	 */
	doublereal Jac = -density*.5*area*sqrt(2./(dKappa*density*jumpPres));
	
	WM.PutCoef(1, 1, Jac);
	WM.PutCoef(1, 2, -Jac);
	WM.PutCoef(2, 1, -Jac);
	WM.PutCoef(2, 2, Jac);
	
	return WorkMat;
}

SubVectorHandler& 
ThreeWayMinorLoss::AssRes(
		SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr
		)
{
	DEBUGCOUT("Entering ThreeWayMinorLoss::AssRes()" << std::endl);
	
	WorkVec.Resize(2);
	
	doublereal p0 = pNode0->dGetX();
	doublereal p1 = pNode1->dGetX();
	doublereal p2 = pNode2->dGetX();
	doublereal p;

	pNodeN = NULL;

	if (p1 > p2) {
		pNodeN = pNode1;
		p = p1;
		area = area1;
	} else {
		pNodeN = pNode2;
		p = p2;
		area = area2;
	}
	
	doublereal jumpPres = fabs(p0-p); 
	
	if (p0 > p) {
		dKappa = dKappa1;  /* flusso diretto da 0 a n */
	} else { 
		dKappa = dKappa2;  /* flusso diretto da n a 0 */
	}
	
	doublereal density = HF->dGetDensity((p0+p)/2.);
	flow = density*area*sqrt(2./(dKappa*density))*copysign(sqrt(jumpPres), p0-p);
	vel = flow/(density*area);
	
#ifdef HYDR_DEVEL
	DEBUGCOUT("RES area :           " << area << std::endl);
	DEBUGCOUT("RES flow:            " << flow << std::endl);
	DEBUGCOUT("RES p0:              " << p0 << std::endl);
	DEBUGCOUT("RES p:               " << p << std::endl);
	DEBUGCOUT("RES dKappa:          " << dKappa << std::endl);
	DEBUGCOUT("****************************************************" << std::endl);
	DEBUGCOUT("RES velocita':       " << vel << std::endl);
	DEBUGCOUT("    se positiva il fluido va dal nodo 0 al nodo n " << std::endl);
	DEBUGCOUT("RES portata (nodo n):" << flow << std::endl);
	DEBUGCOUT("****************************************************" << std::endl);
#endif
	integer iNode0RowIndex = pNode0->iGetFirstRowIndex()+1;
	integer iNodeNRowIndex = pNodeN->iGetFirstRowIndex()+1;
	
	WorkVec.PutItem(1, iNode0RowIndex, flow);
	WorkVec.PutItem(2, iNodeNRowIndex, -flow);         
	
	return WorkVec;
}

void 
ThreeWayMinorLoss::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) { 
		std::ostream& out = OH.Hydraulic();
		out << std::setw(8) << GetLabel() 
			<< " " << vel  << " " << flow << std::endl;
	}
}

/* ThreeWayMinorLoss - end */


/* Orifice - begin */

/* se Re < Rec avrò sicuramente moto laminare
 * se invece Re > Rec avrò sicuramente moto turbolento */

Orifice::Orifice(unsigned int uL, const DofOwner* pDO, 
		 HydraulicFluid* hf,
		 const PressureNode* p1, const PressureNode* p2,
		 doublereal Dh, doublereal A_diaf, doublereal A_pipe, doublereal ReCr, flag fOut)
: Elem(uL, fOut), 
HydraulicElem(uL, pDO, hf, fOut),
pNode1(p1), pNode2(p2),
diameter(Dh), area_diaf(A_diaf),
area_pipe(A_pipe), ReCr(ReCr),
flow(0.),
vel(0.),
Re(0.),
turbulent(false)
{
	ASSERT(pNode1 != NULL);
	ASSERT(pNode1->GetNodeType() == Node::HYDRAULIC);
	ASSERT(pNode2 != NULL);
	ASSERT(pNode2->GetNodeType() == Node::HYDRAULIC);
	ASSERT(Dh > std::numeric_limits<doublereal>::epsilon());
	ASSERT(A_diaf > std::numeric_limits<doublereal>::epsilon());
	ASSERT(A_pipe > std::numeric_limits<doublereal>::epsilon());
   
	/* se |p1-p2| < CriticJump avrò sicuramente moto laminare se no turbolento */
   
	viscosity = HF->dGetViscosity();
	/*
	 * Merritt, pp. 43-45:
	 *
	 * delta = 0.2 for a sharp-edged round orifice
	 * delta = 0.157 for a sharp-edged slit orifice
	 */
	delta = 0.2;
	// ReCr = pow(0.611/delta, 2);
	doublereal density = HF->dGetDensity((pNode2->dGetX()+pNode1->dGetX())/2.);
	CriticJump = area_pipe*ReCr*pow(viscosity/(diameter*delta), 2.)/(2.*density*area_diaf);   

	/* calcolo del Cd */
      
	/* Cc = funzione delle due aree; */
	doublereal rapp = area_diaf/area_pipe;
	/* parametro adimensionale in funzione del rapporto area_diaf/area_pipe */
	doublereal Cc = (((.4855*rapp-.4971)*rapp+.158)*rapp+.1707)*rapp+.6005;
      
	doublereal base = Cc*rapp;
	doublereal rad = 1.-base*base;
	if (rad < 1.e3*std::numeric_limits<doublereal>::epsilon()) {
		silent_cerr("Orifice(" << GetLabel() << "): error computing Cd" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	doublereal Cv = .98; /* costante adimensionale */
	Cd = Cv*Cc/sqrt(rad);
}

Orifice::~Orifice(void)
{
	NO_OP;
}
   
/* Tipo di elemento idraulico (usato solo per debug ecc.) */
HydraulicElem::Type
Orifice::GetHydraulicType(void) const 
{
	return HydraulicElem::ORIFICE;
}

/* Contributo al file di restart */
std::ostream&
Orifice::Restart(std::ostream& out) const
{
	return out << "Orifice not implemented yet!" << std::endl;
}
   
unsigned int
Orifice::iGetNumDof(void) const 
{
	return 0;
}
   
DofOrder::Order
Orifice::GetDofType(unsigned int i) const
{
	silent_cerr("Orifice has no dofs!" << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

void
Orifice::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 2; 
	*piNumCols = 2; 
}

VariableSubMatrixHandler& 
Orifice::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering Orifice::AssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.Resize(2, 2);

	integer iNode1RowIndex = pNode1->iGetFirstRowIndex()+1;
	integer iNode2RowIndex = pNode2->iGetFirstRowIndex()+1;
	integer iNode1ColIndex = pNode1->iGetFirstColIndex()+1;
	integer iNode2ColIndex = pNode2->iGetFirstColIndex()+1;

	WM.PutRowIndex(1, iNode1RowIndex);
	WM.PutRowIndex(2, iNode2RowIndex);
	WM.PutColIndex(1, iNode1ColIndex);
	WM.PutColIndex(2, iNode2ColIndex);

	doublereal p1 = pNode1->dGetX();
	doublereal p2 = pNode2->dGetX();
	doublereal jumpPres = fabs(p1-p2);

	doublereal Jac = 0.;
	doublereal density = HF->dGetDensity((p1+p2)/2.);

	if (jumpPres < CriticJump) {
		/*  moto sicuramente laminare (jacobiano) */
		Jac = -density*2.*(delta*delta)*diameter*area_diaf/viscosity;

	} else {
		/*  moto sicuramente turbolento  (jacobiano) */
		Jac = -Cd*area_diaf/sqrt(2.*jumpPres/density);
	}

	/* Jac *= dCoef; */

	WM.PutCoef(1, 1, Jac);
	WM.PutCoef(1, 2, -Jac);
	WM.PutCoef(2, 1, -Jac);
	WM.PutCoef(2, 2, Jac);

	return WorkMat;
}

SubVectorHandler&
Orifice::AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering Orifice::AssRes()" << std::endl);
 
	WorkVec.Resize(2);

	integer iNode1RowIndex = pNode1->iGetFirstRowIndex()+1;
	integer iNode2RowIndex = pNode2->iGetFirstRowIndex()+1;

	doublereal p1 = pNode1->dGetX();
	doublereal p2 = pNode2->dGetX();
   
	doublereal jumpPres = fabs(p1-p2);

	doublereal density = HF->dGetDensity((p1+p2)/2.);

	if (jumpPres < CriticJump) {
		/*  moto sicuramente laminare (residuo) */
#ifdef HYDR_DEVEL
		DEBUGCOUT("we are in orifice laminar" << std::endl);
#endif
		flow = density*2.*(delta*delta)*diameter*area_diaf*(p1-p2)/viscosity;
		turbulent = false;

	} else {
		/*  moto sicuramente turbolento  (residuo) */
#ifdef HYDR_DEVEL
		DEBUGCOUT("we are in orifice turbulent:" << std::endl);
#endif
		flow = density*Cd*area_diaf*copysign(sqrt(2.*jumpPres/density), p1-p2);
		turbulent = true;
	}

	vel = flow/(density*area_pipe);
	Re = fabs(density*vel*diameter/viscosity);

#ifdef HYDR_DEVEL
	DEBUGCOUT("jumpPres:       " << jumpPres << std::endl);
	DEBUGCOUT("density:        " << density << std::endl);
	DEBUGCOUT("Cd:             " << Cd << std::endl);
	DEBUGCOUT("p1:             " << p1 << std::endl);
	DEBUGCOUT("p2:             " << p2 << std::endl);
	DEBUGCOUT("Cc:             " << Cc << std::endl);
	DEBUGCOUT("CriticJump:     " << CriticJump << std::endl);
	DEBUGCOUT("delta:          " << delta << std::endl);
	DEBUGCOUT("viscosity:      " << viscosity << std::endl);
	DEBUGCOUT("rad:            " << rad << std::endl);
	DEBUGCOUT("area_pipe:      " << area_pipe << std::endl);
	DEBUGCOUT("area_diaf:      " << area_diaf << std::endl);
	DEBUGCOUT("RES area_pipe : " << area_pipe << std::endl);
	DEBUGCOUT("RES flow:       " << flow << std::endl);
	DEBUGCOUT("RES Reynolds:   " << Re << std::endl);
	DEBUGCOUT("******************************************" << std::endl);
	DEBUGCOUT("RES velocita': " << vel << std::endl);
	DEBUGCOUT("    se positiva il fluido va dal nodo 1 al nodo 2" << std::endl);
	DEBUGCOUT("*********************************************" << std::endl);
#endif

	WorkVec.PutItem(1, iNode1RowIndex, flow);	
	WorkVec.PutItem(2, iNode2RowIndex, -flow);

	return WorkVec;
}

void
Orifice::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) { 
		std::ostream& out = OH.Hydraulic();
		out << std::setw(8) << GetLabel()	/*  1 */
			<< " " << vel			/*  2 */
			<< " " << flow			/*  3 */
			<< " " << Re			/*  4 */
			<< " " << turbulent		/*  5 */
			<< std::endl;
	}
}

/* Orifice - end */

