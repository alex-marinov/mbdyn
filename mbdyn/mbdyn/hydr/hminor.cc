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

#include <hminor.h>

/* Minor_losses - begin */

Minor_loss::Minor_loss(unsigned int uL, const DofOwner* pDO,
			   HydraulicFluid* hf,
			   const PressureNode* p1, const PressureNode* p2,
			   doublereal dK1, doublereal dK2, doublereal A,
			   flag fOut)
: Elem(uL, Elem::HYDRAULIC, fOut),
HydraulicElem(uL, pDO, hf, fOut),
pNode1(p1), pNode2(p2),
dKappa1(dK1), dKappa2(dK2), area(A)
{
   ASSERT(pNode1 != NULL);
   ASSERT(pNode1->GetNodeType() == Node::HYDRAULIC);
   ASSERT(pNode2 != NULL);
   ASSERT(pNode2->GetNodeType() == Node::HYDRAULIC);
   ASSERT(dK1 >= 0.);
   ASSERT(dK2 >= 0.);
   ASSERT(A > DBL_EPSILON);
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
ostream& Minor_loss::Restart(ostream& out) const
{
   return out << "Minor_loss not implemented yet!" << endl;
}
   
unsigned int Minor_loss::iGetNumDof(void) const { 
   return 0;
}
   
DofOrder::Order Minor_loss::SetDof(unsigned int i) const {
   cerr << "Minor_loss has no dofs!" << endl;
   THROW(ErrGeneric());
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
   DEBUGCOUT("Entering Minor_loss::AssJac()" << endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.Resize(2, 2);
   
   integer iNode1RowIndex = pNode1->iGetFirstRowIndex()+1;
   integer iNode2RowIndex = pNode2->iGetFirstRowIndex()+1;
   integer iNode1ColIndex = pNode1->iGetFirstColIndex()+1;
   integer iNode2ColIndex = pNode2->iGetFirstColIndex()+1;
   
   WM.fPutRowIndex(1, iNode1RowIndex);
   WM.fPutRowIndex(2, iNode2RowIndex);
   WM.fPutColIndex(1, iNode1ColIndex);
   WM.fPutColIndex(2, iNode2ColIndex);
   
   doublereal p1 = pNode1->dGetX();
   doublereal p2 = pNode2->dGetX();
   
   doublereal jumpPres = fabs(p1-p2);
   
   /* evito di dividere per un numero troppo piccolo */
   if (jumpPres < 1.e8*DBL_EPSILON) {
      jumpPres = 1.e8*DBL_EPSILON;
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
    
   WM.fPutCoef(1, 1, Jac);
   WM.fPutCoef(1, 2, -Jac);
   WM.fPutCoef(2, 1, -Jac);
   WM.fPutCoef(2, 2, Jac);
   
   return WorkMat;
}

SubVectorHandler& 
Minor_loss::AssRes(SubVectorHandler& WorkVec,
		     doublereal dCoef,
		     const VectorHandler& XCurr, 
		     const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering Minor_loss::AssRes()" << endl);
   
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
   DEBUGCOUT("RES area :           " << area << endl);
   DEBUGCOUT("RES flow:            " << flow << endl);
   DEBUGCOUT("RES p1:              " << p1 << endl);
   DEBUGCOUT("RES p2:              " << p2 << endl);
   DEBUGCOUT("RES dKappa:          " << dKappa << endl);
   DEBUGCOUT("****************************************************" << endl);
   DEBUGCOUT("RES velocita':       " << vel << endl);
   DEBUGCOUT("    se positiva il fluido va dal nodo 1 al nodo 2 " << endl);
   DEBUGCOUT("RES portata (nodo2): " << flow << endl);
   DEBUGCOUT("****************************************************" << endl);
#endif
   WorkVec.fPutItem(1, iNode1RowIndex, flow);
   WorkVec.fPutItem(2, iNode2RowIndex, -flow);         

   return WorkVec;
}
   
void Minor_loss::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) { 
      ostream& out = OH.Hydraulic();
      out << setw(8) << GetLabel() 
          << " " << vel  << " " << flow << endl;
   }
}

void Minor_loss::SetValue(VectorHandler& /* X */ , VectorHandler& /* XP */ ) const {
   NO_OP;
}

/* Minor_loss - end */


/* ThreeWayMinorLoss - begin */

ThreeWayMinorLoss::ThreeWayMinorLoss(
		unsigned int uL, const DofOwner* pDO,
		HydraulicFluid* hf, const PressureNode* p0,
		const PressureNode* p1, const PressureNode* p2,
		doublereal dK1, doublereal dK2, 
		doublereal A1, doublereal A2, flag fOut)
: Elem(uL, Elem::HYDRAULIC, fOut),
HydraulicElem(uL, pDO, hf, fOut),
pNode0(p0), pNode1(p1), pNode2(p2), pNodeN(NULL),
dKappa1(dK1), dKappa2(dK2), area1(A1), area2(A2)
{
	ASSERT(pNode0 != NULL);
	ASSERT(pNode0->GetNodeType() == Node::HYDRAULIC);
	ASSERT(pNode1 != NULL);
	ASSERT(pNode1->GetNodeType() == Node::HYDRAULIC);
	ASSERT(pNode2 != NULL);
	ASSERT(pNode2->GetNodeType() == Node::HYDRAULIC);
	ASSERT(dK1 >= 0.);
	ASSERT(dK2 >= 0.);
	ASSERT(A1 > DBL_EPSILON);
	ASSERT(A2 > DBL_EPSILON);
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
ostream& ThreeWayMinorLoss::Restart(ostream& out) const
{
	return out << "ThreeWayMinorLoss not implemented yet!" << endl;
}
 
unsigned int ThreeWayMinorLoss::iGetNumDof(void) const { 
	return 0;
}
   
DofOrder::Order ThreeWayMinorLoss::SetDof(unsigned int i) const {
	cerr << "Minor_loss has no dofs!" << endl;
	THROW(ErrGeneric());
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
	DEBUGCOUT("Entering Minor_loss::AssJac()" << endl);
   
	ASSERT(pNodeN != NULL);
	
	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.Resize(2, 2);
	
	integer iNode0RowIndex = pNode0->iGetFirstRowIndex()+1;
	integer iNodeNRowIndex = pNodeN->iGetFirstRowIndex()+1;
	integer iNode0ColIndex = pNode0->iGetFirstColIndex()+1;
	integer iNodeNColIndex = pNodeN->iGetFirstColIndex()+1;
	
	WM.fPutRowIndex(1, iNode0RowIndex);
	WM.fPutRowIndex(2, iNodeNRowIndex);
	WM.fPutColIndex(1, iNode0ColIndex);
	WM.fPutColIndex(2, iNodeNColIndex);

	doublereal p0 = pNode0->dGetX();
	doublereal p = pNodeN->dGetX();
	
	doublereal jumpPres = fabs(p0-p);
	
	/* evito di dividere per un numero troppo piccolo */
	if (jumpPres < 1.e8*DBL_EPSILON) {
		jumpPres = 1.e8*DBL_EPSILON;
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
	
	WM.fPutCoef(1, 1, Jac);
	WM.fPutCoef(1, 2, -Jac);
	WM.fPutCoef(2, 1, -Jac);
	WM.fPutCoef(2, 2, Jac);
	
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
	DEBUGCOUT("Entering ThreeWayMinorLoss::AssRes()" << endl);
	
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
	DEBUGCOUT("RES area :           " << area << endl);
	DEBUGCOUT("RES flow:            " << flow << endl);
	DEBUGCOUT("RES p0:              " << p0 << endl);
	DEBUGCOUT("RES p:               " << p << endl);
	DEBUGCOUT("RES dKappa:          " << dKappa << endl);
	DEBUGCOUT("****************************************************" << endl);
	DEBUGCOUT("RES velocita':       " << vel << endl);
	DEBUGCOUT("    se positiva il fluido va dal nodo 0 al nodo n " << endl);
	DEBUGCOUT("RES portata (nodo n):" << flow << endl);
	DEBUGCOUT("****************************************************" << endl);
#endif
	integer iNode0RowIndex = pNode0->iGetFirstRowIndex()+1;
	integer iNodeNRowIndex = pNodeN->iGetFirstRowIndex()+1;
	
	WorkVec.fPutItem(1, iNode0RowIndex, flow);
	WorkVec.fPutItem(2, iNodeNRowIndex, -flow);         
	
	return WorkVec;
}

void 
ThreeWayMinorLoss::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) { 
		ostream& out = OH.Hydraulic();
		out << setw(8) << GetLabel() 
			<< " " << vel  << " " << flow << endl;
	}
}

void 
ThreeWayMinorLoss::SetValue(
		VectorHandler& /* X */ , 
		VectorHandler& /* XP */ 
		) const {
	NO_OP;
}

/* ThreeWayMinorLoss - end */


/* Orifice - begin */

/* se Re < Rec avrò sicuramente moto laminare
 * se invece Re > Rec avrò sicuramente moto turbolento */

Orifice::Orifice(unsigned int uL, const DofOwner* pDO, 
		 HydraulicFluid* hf,
		 const PressureNode* p1, const PressureNode* p2,
		 doublereal Dh, doublereal A_diaf, doublereal A_pipe, doublereal ReCr, flag fOut)
: Elem(uL, Elem::HYDRAULIC, fOut), 
HydraulicElem(uL, pDO, hf, fOut),
pNode1(p1), pNode2(p2),
diameter(Dh), area_diaf(A_diaf),
area_pipe(A_pipe),  ReCr(ReCr)
{
   ASSERT(pNode1 != NULL);
   ASSERT(pNode1->GetNodeType() == Node::HYDRAULIC);
   ASSERT(pNode2 != NULL);
   ASSERT(pNode2->GetNodeType() == Node::HYDRAULIC);
   ASSERT(Dh > DBL_EPSILON);
   ASSERT(A_diaf > DBL_EPSILON);
   ASSERT(A_pipe > DBL_EPSILON);
   
   /* se |p1-p2| < CriticJump avrò sicuramente moto laminare se no turbolento */
   
   viscosity = HF->dGetViscosity();
   delta = .2;            /* numero adimensionale Cd=delta*sqrt(Reynolds) */
   doublereal density = HF->dGetDensity((pNode2->dGetX()+pNode1->dGetX())/2.);
   CriticJump = area_pipe*ReCr*pow(viscosity/(diameter*delta), 2.)
     /(2.*density*area_diaf);   
}

Orifice::~Orifice(void)
{
   NO_OP;
}
   
/* Tipo di elemento idraulico (usato solo per debug ecc.) */
HydraulicElem::Type Orifice::GetHydraulicType(void) const 
{
   return HydraulicElem::ORIFICE;
}

/* Contributo al file di restart */
ostream& Orifice::Restart(ostream& out) const
{
   return out << "Orifice not implemented yet!" << endl;
}
   
unsigned int Orifice::iGetNumDof(void) const 
{
   return 0;
}
   
DofOrder::Order Orifice::SetDof(unsigned int i) const {
   cerr << "Orifice has no dofs!" << endl;
   THROW(ErrGeneric());
}

void Orifice::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
   *piNumRows = 2; 
   *piNumCols = 2; 
}
     
VariableSubMatrixHandler& 
Orifice::AssJac(VariableSubMatrixHandler& WorkMat,
		  doublereal dCoef,
		  const VectorHandler& XCurr, 
		  const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering Orifice::AssJac()" << endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.Resize(2, 2);
   
   integer iNode1RowIndex = pNode1->iGetFirstRowIndex()+1;
   integer iNode2RowIndex = pNode2->iGetFirstRowIndex()+1;
   integer iNode1ColIndex = pNode1->iGetFirstColIndex()+1;
   integer iNode2ColIndex = pNode2->iGetFirstColIndex()+1;
   
   WM.fPutRowIndex(1, iNode1RowIndex);
   WM.fPutRowIndex(2, iNode2RowIndex);
   WM.fPutColIndex(1, iNode1ColIndex);
   WM.fPutColIndex(2, iNode2ColIndex);
   
   doublereal p1 = pNode1->dGetX();
   doublereal p2 = pNode2->dGetX();
   doublereal jumpPres = fabs(p1-p2);

   doublereal Jac = 0.;
   doublereal Cv = .98; /* costante adimensionale */
   doublereal Cc;       /* parametro adimensionale in funzione del rapporto area_diaf/area_pipe */
   doublereal Cd;
   doublereal rad;
   doublereal density = HF->dGetDensity((p1+p2)/2.);
   
   if (jumpPres < CriticJump) {
      /*  moto sicuramente laminare (jacobiano) */
      Jac = -density*2.*(delta*delta)*diameter*area_diaf/viscosity;
   } else {
      /*  moto sicuramente turbolento  (jacobiano) */
      /* calcolo del Cd */
      /* Cc= funzione delle due aree; */
      doublereal rapp = area_diaf/area_pipe;
      doublereal rapp2 = rapp*rapp;
      doublereal rapp3 = rapp2*rapp;
      doublereal rapp4 = rapp3*rapp;
      Cc = .4855*rapp4-.4971*rapp3+.158*rapp2+.1707*rapp+.6005;

      rad = 1.-pow(Cc*area_diaf/area_pipe, 2.);
	   
      if (rad < 1.e3*DBL_EPSILON) {		  
	 cerr << "error in orifice: Cd" << endl;
	 THROW(ErrGeneric());
      }

      Cd = Cv*Cc/sqrt(rad);
      Jac = -density*Cd*area_diaf/(sqrt(2.*jumpPres/density)*density);
   }

   /* Jac *= dCoef; */
   
   WM.fPutCoef(1, 1, Jac);
   WM.fPutCoef(1, 2, -Jac);
   WM.fPutCoef(2, 1, -Jac);
   WM.fPutCoef(2, 2, Jac);
   
   return WorkMat;
}

SubVectorHandler& Orifice::AssRes(SubVectorHandler& WorkVec,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering Orifice::AssRes()" << endl);
   
   WorkVec.Resize(2);
   
   integer iNode1RowIndex = pNode1->iGetFirstRowIndex()+1;
   integer iNode2RowIndex = pNode2->iGetFirstRowIndex()+1;
   
   doublereal p1 = pNode1->dGetX();
   doublereal p2 = pNode2->dGetX();
   
   doublereal jumpPres = fabs(p1-p2);

   doublereal Cv = .98;
   doublereal Cc;  /* = .6;   valori inutili */
   doublereal Cd;  /* = .1;   valori inutili */
   doublereal rad; /* = 1.;   valori inutili */
   doublereal density = HF->dGetDensity((p1+p2)/2.);
   
   if (jumpPres < CriticJump) {
      /*  moto sicuramente laminare (residuo) */
#ifdef HYDR_DEVEL
      DEBUGCOUT("we are in orifice laminar" << endl);
#endif
      flow = density*2.*(delta*delta)*diameter*area_diaf*(p1-p2)/viscosity;
   } else {
      /*  moto sicuramente turbolento  (residuo) */
#ifdef HYDR_DEVEL
      DEBUGCOUT("we are in orifice turbulent:" << endl);
#endif
      /* calcolo del Cd */
      
      /* Cc= funzione delle due aree; */
      doublereal rapp = area_diaf/area_pipe;
      doublereal rapp2 = rapp*rapp;
      doublereal rapp3 = rapp2*rapp;
      doublereal rapp4 = rapp3*rapp;
      Cc = .4855*rapp4-.4971*rapp3+.158*rapp2+.1707*rapp+.6005;
      /* Cc = (((.4855*rapp-.4971)*rapp+.158)*rapp+.1707)*rapp+.6005; */
      
      doublereal base= Cc*area_diaf/area_pipe;
      rad = 1.-base*base;
      if (rad < 1.e3*DBL_EPSILON) {		  
	 cerr << "error in orifice: Cd" << endl;
	 THROW(ErrGeneric());
      }	   
		      
      Cd = Cv*Cc/sqrt(rad);
      flow = density*Cd*area_diaf*copysign(sqrt(2.*jumpPres/density), p1-p2);
   }
    
   vel = flow/(density*area_pipe);
   Re = fabs(density*vel*diameter/viscosity);
   
#ifdef HYDR_DEVEL
   DEBUGCOUT("jumpPres:       " << jumpPres << endl);
   DEBUGCOUT("density:        " << density << endl);
   DEBUGCOUT("Cd:             " << Cd << endl);
   DEBUGCOUT("p1:             " << p1 << endl);
   DEBUGCOUT("p2:             " << p2 << endl);
   DEBUGCOUT("Cc:             " << Cc << endl);
   DEBUGCOUT("CriticJump:     " << CriticJump << endl);
   DEBUGCOUT("delta:          " << delta << endl);
   DEBUGCOUT("viscosity:      " << viscosity << endl);
   DEBUGCOUT("rad:            " << rad << endl);
   DEBUGCOUT("area_pipe:      " << area_pipe << endl);
   DEBUGCOUT("area_diaf:      " << area_diaf << endl);
   DEBUGCOUT("RES area_pipe : " << area_pipe << endl);
   DEBUGCOUT("RES flow:       " << flow << endl);
   DEBUGCOUT("RES Reynolds:   " << Re << endl);
   DEBUGCOUT("******************************************" << endl);
   DEBUGCOUT("RES velocita': " << vel << endl);
   DEBUGCOUT("    se positiva il fluido va dal nodo 1 al nodo 2" << endl);
   DEBUGCOUT("*********************************************" << endl);
#endif
   
   WorkVec.fPutItem(1, iNode1RowIndex, flow);	
   WorkVec.fPutItem(2, iNode2RowIndex, -flow);

   return WorkVec;
}

void Orifice::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) { 
      ostream& out = OH.Hydraulic();
      out << setw(8) << GetLabel() << " " 
        << vel << " " 
	<< flow << " " 
	<< Re << endl;
   }
}

void Orifice::SetValue(VectorHandler& /* X */ , VectorHandler& /* XP */ ) const 
{
   NO_OP;
}

/* Orifice - end */

