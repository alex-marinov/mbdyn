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

/* Swash plate */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#ifdef USE_ELECTRIC_NODES

#include <swashpl.h>

extern "C" {
#include <mymath.h>
}

/* SwashPlate - begin */

SwashPlate::SwashPlate(unsigned int uL, const DofOwner* pDO,
		       const AbstractNode* pCollIn, // const DriveCaller* pColl, 
		       const AbstractNode* pLongIn, // const DriveCaller* pLong, 
		       const AbstractNode* pLatIn,  // const DriveCaller* pLat,
		       const AbstractNode* pN1,
		       const AbstractNode* pN2,
		       const AbstractNode* pN3,
		       doublereal dDynCoef,
		       doublereal dCyclFact,
		       doublereal dCollFact,
		       flag fCL,
		       doublereal dCMin,
		       doublereal dCMax,
		       flag fFL,
		       doublereal dFMin,
		       doublereal dFMax,
		       flag fLL,
		       doublereal dLMin,	      
		       doublereal dLMax,
		       flag fOut)
: Elem(uL, Elem::GENEL, fOut), 
Genel(uL, Genel::SWASHPLATE, pDO, fOut),
pCollectiveIn(pCollIn),   // Collective(pColl),
pLongitudinalIn(pLongIn), // Longitudinal(pLong),
pLateralIn(pLatIn),       // Lateral(pLat),
pNode1(pN1), pNode2(pN2), pNode3(pN3),
dDynamicCoef(dDynCoef), 
dCyclicFactor(dCyclFact), 
dCollectiveFactor(dCollFact),
fCollLimits(fCL), dCollMax(dCMax), dCollMin(dCMin),
fForeAftLimits(fFL), dForeAftMax(dFMax), dForeAftMin(dFMin),
fLatLimits(fLL), dLatMax(dLMax), dLatMin(dLMin)
{
   ASSERT(pCollectiveIn != NULL);
   ASSERT(pCollectiveIn->GetNodeType() == Node::ABSTRACT);
   ASSERT(pLongitudinalIn != NULL);
   ASSERT(pLongitudinalIn->GetNodeType() == Node::ABSTRACT);
   ASSERT(pLateralIn != NULL);
   ASSERT(pLateralIn->GetNodeType() == Node::ABSTRACT);
   
   ASSERT(pNode1 != NULL);
   ASSERT(pNode1->GetNodeType() == Node::ABSTRACT);
   ASSERT(pNode2 != NULL);
   ASSERT(pNode2->GetNodeType() == Node::ABSTRACT);
   ASSERT(pNode3 != NULL);
   ASSERT(pNode3->GetNodeType() == Node::ABSTRACT);
   
   ASSERT(dCyclicFactor != 0.);
   ASSERT(dCollectiveFactor != 0.);
   ASSERT(dDynamicCoef >= 0.);   
}


SwashPlate::~SwashPlate(void)
{
   NO_OP;
}


/* Scrive il contributo dell'elemento al file di restart */
ostream& SwashPlate::Restart(ostream& out) const
{
   Genel::Restart(out) << ", swash plate, "
     << pCollectiveIn->GetLabel() << ", ";
   // Collective.pGetDriveCaller()->Restart(out) << ", ";
   if(fCollLimits) {
      out << "limits, " << dCollMin << ", " << dCollMax << ", ";
   }
   out << pLongitudinalIn->GetLabel() << ", ";
   // Longitudinal.pGetDriveCaller()->Restart(out) << ", ";
   if(fForeAftLimits) {
      out << "limits, " << dForeAftMin << ", " << dForeAftMax << ", ";
   }   
   out << pLateralIn->GetLabel() << ", ";
   // Lateral.pGetDriveCaller()->Restart(out) << ", ";
   if(fLatLimits) {
      out << "limits, " << dLatMin << ", " << dLatMax << ", ";
   }
   out << pNode1->GetLabel() << ", "
     << pNode2->GetLabel() << ", "
     << pNode3->GetLabel() << ", "
     << dDynamicCoef << ", "
     << dCyclicFactor << ", " << dCollectiveFactor << ';' << endl;
   
   return out;
}


/* assemblaggio jacobiano */
VariableSubMatrixHandler& 
SwashPlate::AssJac(VariableSubMatrixHandler& WorkMat,
		   doublereal dCoef, 
		   const VectorHandler& /* XCurr */ ,
		   const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering SwashPlate::AssJac()" << endl);

   /* Casting di WorkMat */
   SparseSubMatrixHandler& WM = WorkMat.SetSparse();   
   WM.ResizeInit(6, 0, 0.);
      
   integer iCollFirstIndex = pCollectiveIn->iGetFirstIndex()+1;
   integer iLongFirstIndex = pLongitudinalIn->iGetFirstIndex()+1;
   integer iLatFirstIndex = pLateralIn->iGetFirstIndex()+1;

   integer iNode1FirstIndex = pNode1->iGetFirstIndex()+1;
   integer iNode2FirstIndex = pNode2->iGetFirstIndex()+1;
   integer iNode3FirstIndex = pNode3->iGetFirstIndex()+1;
   
   doublereal d = dDynamicCoef+dCoef;
   
   WM.fPutItem(1, iCollFirstIndex, iCollFirstIndex, dCoef);
   WM.fPutItem(2, iLongFirstIndex, iLongFirstIndex, dCoef);
   WM.fPutItem(3, iLatFirstIndex, iLatFirstIndex, dCoef);

   WM.fPutItem(4, iNode1FirstIndex, iNode1FirstIndex, d);
   WM.fPutItem(5, iNode2FirstIndex, iNode2FirstIndex, d);
   WM.fPutItem(6, iNode3FirstIndex, iNode3FirstIndex, d);   

   return WorkMat;
}


/* assemblaggio residuo */
SubVectorHandler& 
SwashPlate::AssRes(SubVectorHandler& WorkVec,
		   doublereal /* dCoef */ ,
		   const VectorHandler& XCurr, 
		   const VectorHandler& XPrimeCurr)
{   
   DEBUGCOUT("Entering SwashPlate::AssRes()" << endl);

   /* Dimensiona e resetta la matrice di lavoro */
   WorkVec.Resize(6);
   WorkVec.Reset(0.);
      
   integer iCollFirstIndex = pCollectiveIn->iGetFirstIndex()+1;
   integer iLongFirstIndex = pLongitudinalIn->iGetFirstIndex()+1;
   integer iLatFirstIndex = pLateralIn->iGetFirstIndex()+1;

   integer iNode1FirstIndex = pNode1->iGetFirstIndex()+1;
   integer iNode2FirstIndex = pNode2->iGetFirstIndex()+1;
   integer iNode3FirstIndex = pNode3->iGetFirstIndex()+1;
   
   WorkVec.fPutRowIndex(1, iCollFirstIndex);
   WorkVec.fPutRowIndex(2, iLongFirstIndex);
   WorkVec.fPutRowIndex(3, iLatFirstIndex);
         
   WorkVec.fPutRowIndex(4, iNode1FirstIndex);
   WorkVec.fPutRowIndex(5, iNode2FirstIndex);
   WorkVec.fPutRowIndex(6, iNode3FirstIndex);
         
   doublereal dXColl = XCurr.dGetCoef(iCollFirstIndex);
   doublereal dXLong = XCurr.dGetCoef(iLongFirstIndex);
   doublereal dXLat = XCurr.dGetCoef(iLatFirstIndex);
      
   WorkVec.fPutCoef(1, -dXColl);
   WorkVec.fPutCoef(2, -dXLong);
   WorkVec.fPutCoef(3, -dXLat);
   
   /* Limits on pitch angles */
   if(fCollLimits) {
      if(dXColl > dCollMax) {
	 dXColl = dCollMax;
      } else if(dXColl < dCollMin) {
	 dXColl = dCollMin;
      }
   }
   
   if(fForeAftLimits) {
      if(dXLong > dForeAftMax) {
	 dXLong = dForeAftMax;
      } else if(dXLong < dForeAftMin) {
	 dXLong = dForeAftMin;
      }
   }
   
   if(fLatLimits) {
      if(dXLat > dLatMax) {
	 dXLat = dLatMax;
      } else if(dXLat < dLatMin) {
	 dXLat = dLatMin;
      }
   }   
   
   WorkVec.fPutCoef(4, dCollectiveFactor*(dXColl-dCyclicFactor*dXLong)
		    -XCurr.dGetCoef(iNode1FirstIndex)
		    -dDynamicCoef*XPrimeCurr.dGetCoef(iNode1FirstIndex));
   WorkVec.fPutCoef(5, dCollectiveFactor*(dXColl+dCyclicFactor*(.5*dXLong-sqrt(3.)/2.*dXLat))
		    -XCurr.dGetCoef(iNode2FirstIndex)
		    -dDynamicCoef*XPrimeCurr.dGetCoef(iNode2FirstIndex));
   WorkVec.fPutCoef(6, dCollectiveFactor*(dXColl+dCyclicFactor*(.5*dXLong+sqrt(3.)/2.*dXLat))
		    -XCurr.dGetCoef(iNode3FirstIndex)
		    -dDynamicCoef*XPrimeCurr.dGetCoef(iNode3FirstIndex));

   return WorkVec;
}


void SwashPlate::Output(OutputHandler& /* OH */ ) const
{
   NO_OP;
}


void SwashPlate::SetInitialValue(VectorHandler& /* X */ ) const
{
   NO_OP;
}


void SwashPlate::SetValue(VectorHandler& /* X */ , 
			  VectorHandler& /* XP */ ) const
{
   NO_OP;
}

/* SwashPlate - end */

#endif /* USE_ELECTRIC_NODES */

