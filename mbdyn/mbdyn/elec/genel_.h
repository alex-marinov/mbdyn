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

#ifndef GENEL__H
#define GENEL__H

#include "genel.h"
#include "drive.h"
#include "constltp.h"

/* GenelClamp - begin */

class GenelClamp : virtual public Elem, public Genel, public DriveOwner {
 protected: 
   ScalarDof SD;
   doublereal dRct;
   
 public:
   GenelClamp(unsigned int uLabel, const DofOwner* pDO, const DriveCaller* pDC,
	      const ScalarDof& sd, flag fOutput)
     : Elem(uLabel, Elem::BULK, fOutput), 
     Genel(uLabel, Genel::CLAMP, pDO, fOutput),
     DriveOwner(pDC), SD(sd), dRct(0.) {
      NO_OP;
   };
   
   virtual ~GenelClamp(void) { 
      NO_OP;
   };
   
   virtual inline void* pGet(void) const { 
      return (void*)this;
   };
   
   virtual unsigned int iGetNumDof(void) const { 
      return 1;
   };   
      
   /* esegue operazioni sui dof di proprieta' dell'elemento */
#ifdef DEBUG   
   virtual DofOrder::Order SetDof(unsigned int i ) const {
#else /* DEBUG */
   virtual DofOrder::Order SetDof(unsigned int /* i */) const {
#endif /* DEBUG */
      ASSERT(i == 0);
      return DofOrder::ALGEBRAIC;
   };
   
   /* Scrive il contributo dell'elemento al file di restart */
   virtual ostream& Restart(ostream& out) const {
      return out; 
   };
   
   /* Tipo di Genel */
   virtual Genel::Type GetGenelType(void) const { 
      return Genel::CLAMP; 
   };
   
   /* Dimensioni del workspace */
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
      *piNumRows = 2;
      *piNumCols = 2;
   };
   
   void Output(OutputHandler& OH ) const {
      if (fToBeOutput()) {
	 ostream& out = OH.Genels();
	 out << setw(8) << GetLabel() << " " << dRct << endl;
      }
   };
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal /* dCoef */ ,
	    const VectorHandler& /* XCurr */ ,
	    const VectorHandler& /* XPrimeCurr */ ) {

	DEBUGCOUT("Entering GenelClamp::AssJac()" << endl);
	
	SparseSubMatrixHandler& WM = WorkMat.SetSparse();     
	WM.Resize(2, 0);
	
	integer iRowIndex = SD.pNode->iGetFirstRowIndex()+1;
	integer iColIndex = SD.pNode->iGetFirstColIndex()+1;
	integer iFirstReactionIndex = iGetFirstIndex()+1;
	
	WM.fPutItem(1, iRowIndex, iFirstReactionIndex, 1.);
	WM.fPutItem(2, iFirstReactionIndex, iColIndex, 1.);
	
	return WorkMat;
     };
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal dCoef,
				    const VectorHandler& XCurr,
				    const VectorHandler& /* XPrimeCurr */ ) {
      DEBUGCOUT("Entering GenelClamp::AssRes()" << endl);
      
      WorkVec.Resize(2);  

      integer iRowIndex = SD.pNode->iGetFirstRowIndex()+1;
      integer iFirstReactionIndex = iGetFirstIndex()+1;
      
      doublereal dVal = SD.pNode->dGetDofValue(1, SD.iOrder);
      dRct = XCurr.dGetCoef(iFirstReactionIndex);
      
      WorkVec.fPutItem(1, iRowIndex, -dRct);
      
      if (SD.iOrder == 1 || SD.pNode->SetDof(0) == DofOrder::ALGEBRAIC) {
	 WorkVec.fPutItem(2, iFirstReactionIndex, dGet()-dVal);
      } else {
	 if (dCoef != 0.) {
	    WorkVec.fPutItem(2, iFirstReactionIndex, (dGet()-dVal)/dCoef);
	 }
      }

      return WorkVec;
   };
   
   void SetValue(VectorHandler& X, VectorHandler& XP) const {
      if (SD.iOrder == 0) {
	 X.fPutCoef(SD.pNode->iGetFirstRowIndex()+1, dGet());
      } else if (SD.iOrder == 1) {
	 XP.fPutCoef(SD.pNode->iGetFirstRowIndex()+1, dGet());	 
      }
   };
   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(int& NumNodes, Node::Type* NdTyps, unsigned int* NdLabels) {
     NumNodes = 1;
     NdTyps[0] = SD.pNode->GetNodeType();
     NdLabels[0] = SD.pNode->GetLabel();
   };
   /* ************************************************ */
};

/* GenelClamp - end */


/* GenelDistance - begin */

class GenelDistance : virtual public Elem, public Genel, public DriveOwner {
 protected: 
   ScalarDof SD1;
   ScalarDof SD2;
   doublereal dRct;
   
 public:
   GenelDistance(unsigned int uLabel, const DofOwner* pDO, 
		 const DriveCaller* pDC,
		 const ScalarDof& sd1, const ScalarDof& sd2, flag fOutput)
     : Elem(uLabel, Elem::BULK, fOutput), 
     Genel(uLabel, Genel::DISTANCE, pDO, fOutput),
     DriveOwner(pDC), SD1(sd1), SD2(sd2), dRct(0.) {
      NO_OP;
   };
   
   virtual ~GenelDistance(void) { 
      NO_OP;
   };
   
   virtual inline void* pGet(void) const { 
      return (void*)this;
   };
   
   virtual unsigned int iGetNumDof(void) const { 
      return 1;
   };   
      
   /* esegue operazioni sui dof di proprieta' dell'elemento */
#ifdef DEBUG   
   virtual DofOrder::Order SetDof(unsigned int i ) const {
#else     
   virtual DofOrder::Order SetDof(unsigned int /* i */ ) const {
#endif    
      ASSERT(i == 0);
      return DofOrder::ALGEBRAIC;
   };
   
   /* Scrive il contributo dell'elemento al file di restart */
   virtual ostream& Restart(ostream& out) const {
      return out; 
   };
  
   /* Tipo di Genel */
   virtual Genel::Type GetGenelType(void) const { 
      return Genel::DISTANCE;
   };
   
   /* Dimensioni del workspace */
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
      *piNumRows = 3;
      *piNumCols = 3;
   };
   
   void Output(OutputHandler& OH) const {
      if (fToBeOutput()) {
	 ostream& out = OH.Genels();
	 out << setw(8) << GetLabel() << " " << dRct << endl;
      }
   };
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef, 
	    const VectorHandler& /* XCurr */ ,
	    const VectorHandler& /* XPrimeCurr */ ) {
	DEBUGCOUT("Entering GenelDistance::AssJac()" << endl);
	
	SparseSubMatrixHandler& WM = WorkMat.SetSparse();
	WM.ResizeInit(4, 0, 0.);
	
	integer iNode1RowIndex = SD1.pNode->iGetFirstRowIndex()+1;
	integer iNode1ColIndex = SD1.pNode->iGetFirstColIndex()+1;
	integer iNode2RowIndex = SD2.pNode->iGetFirstRowIndex()+1;
	integer iNode2ColIndex = SD2.pNode->iGetFirstColIndex()+1;
	integer iFirstReactionIndex = iGetFirstIndex()+1;
	
	WM.fPutItem(1, iNode1RowIndex, iFirstReactionIndex, -1.);
	WM.fPutItem(2, iNode2RowIndex, iFirstReactionIndex, 1.);
	
	doublereal d = dCoef;
	if ((SD1.iOrder == 0) && (SD2.iOrder == 0)) {
	   d = 1.;
	}
	
	if (SD1.iOrder == 1) {
	   WM.fPutItem(3, iFirstReactionIndex, iNode1ColIndex, -1.);
	} else {
	   WM.fPutItem(3, iFirstReactionIndex, iNode1ColIndex, -d);
	}      
	
	if (SD2.iOrder == 1) {
	   WM.fPutItem(4, iFirstReactionIndex, iNode2ColIndex, 1.);
	} else {
	   WM.fPutItem(4, iFirstReactionIndex, iNode2ColIndex, d);
	}      
	
	return WorkMat;
     };
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal dCoef,
				    const VectorHandler& XCurr,
				    const VectorHandler& /* XPrimeCurr */ ) {
      DEBUGCOUT("Entering GenelDistance::AssRes()" << endl);
      
      WorkVec.Resize(3);
      WorkVec.Reset(0.);

      integer iNode1RowIndex = SD1.pNode->iGetFirstRowIndex()+1;
      integer iNode2RowIndex = SD2.pNode->iGetFirstRowIndex()+1;
      integer iFirstReactionIndex = iGetFirstIndex()+1;
      
      doublereal dVal1 = SD1.pNode->dGetDofValue(1, SD1.iOrder);
      doublereal dVal2 = SD2.pNode->dGetDofValue(1, SD2.iOrder);
      dRct = XCurr.dGetCoef(iFirstReactionIndex);
      
      WorkVec.fPutItem(1, iNode1RowIndex, dRct);
      WorkVec.fPutItem(2, iNode2RowIndex, -dRct);
      
      if ((SD1.iOrder == 0) && (SD2.iOrder == 0)) {     
	 if (dCoef != 0.) {
	    WorkVec.fPutItem(3, iFirstReactionIndex, 
			     (dGet()-dVal2+dVal1)/dCoef);
	 }	 
      } else {
	 WorkVec.fPutItem(3, iFirstReactionIndex, dGet()-dVal2+dVal1);
      }

      return WorkVec;
   };
   
   void SetValue(VectorHandler& X, VectorHandler& XP) const {
      if (SD2.iOrder == 0) {
	 X.fPutCoef(SD2.pNode->iGetFirstRowIndex()+1,
		    dGet()-SD2.pNode->dGetX()
		    +SD1.pNode->dGetDofValue(1, SD1.iOrder));
      } else if (SD2.iOrder == 1) {
	 XP.fPutCoef(SD2.pNode->iGetFirstRowIndex()+1, 
		     dGet()-SD2.pNode->dGetXPrime()
		     +SD1.pNode->dGetDofValue(1, SD1.iOrder));
      }
   };
   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(int& NumNodes, Node::Type* NdTyps, unsigned int* NdLabels) {
     NumNodes = 2;
     NdTyps[0] = SD1.pNode->GetNodeType();
     NdLabels[0] = SD1.pNode->GetLabel();
     NdTyps[1] = SD2.pNode->GetNodeType();
     NdLabels[2] = SD2.pNode->GetLabel();

   };
   /* ************************************************ */
};

/* GenelDistance - end */

   
/* GenelSpring - begin */

class GenelSpring
: virtual public Elem, public Genel, public ConstitutiveLaw1DOwner {
 protected: 
   ScalarDof SD1;
   ScalarDof SD2;
   
 public:
   GenelSpring(unsigned int uLabel, const DofOwner* pDO, 
	       const ConstitutiveLaw1D* pCL,
	       const ScalarDof& sd1, const ScalarDof& sd2, flag fOutput)
     : Elem(uLabel, Elem::BULK, fOutput), 
     Genel(uLabel, Genel::SPRING, pDO, fOutput),
     ConstitutiveLaw1DOwner(pCL), SD1(sd1), SD2(sd2) { 
      NO_OP;
   };
   
   virtual ~GenelSpring(void) { 
      NO_OP;
   };
   
   virtual inline void* pGet(void) const { 
      return (void*)this;
   };   
   
   /* Scrive il contributo dell'elemento al file di restart */
   virtual ostream& Restart(ostream& out) const {
      return out; 
   };
  
   /* Tipo di Genel */
   virtual Genel::Type GetGenelType(void) const { 
      return Genel::SPRING; 
   };
   
   /* Dimensioni del workspace */
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
      *piNumRows = 2;
      *piNumCols = 2;
   };
   
   void Output(OutputHandler& /* OH */ ) const {
      NO_OP;
   };
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef, 
	    const VectorHandler& /* XCurr */ ,
	    const VectorHandler& /* XPrimeCurr */ ) {
	DEBUGCOUT("Entering GenelSpring::AssJac()" << endl);
	
	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeInit(2, 2, 0.);
	
	integer iNode1RowIndex = SD1.pNode->iGetFirstRowIndex()+1;
	integer iNode1ColIndex = SD1.pNode->iGetFirstColIndex()+1;
	integer iNode2RowIndex = SD2.pNode->iGetFirstRowIndex()+1;
	integer iNode2ColIndex = SD2.pNode->iGetFirstColIndex()+1;
	
	WM.fPutRowIndex(1, iNode1RowIndex);
	WM.fPutColIndex(1, iNode1ColIndex);
	WM.fPutRowIndex(2, iNode2RowIndex);
	WM.fPutColIndex(2, iNode2ColIndex);
	
	doublereal dFDE = GetFDE();
      	
	if (SD1.iOrder == 1) {
	   WM.fPutCoef(1, 1, dFDE);
	   WM.fPutCoef(2, 1, -dFDE);
	} else {
	   WM.fPutCoef(1, 1, dFDE*dCoef);
	   WM.fPutCoef(2, 1, -dFDE*dCoef);
	}
	
	if (SD2.iOrder == 1) {
	   WM.fPutCoef(1, 2, -dFDE);
	   WM.fPutCoef(2, 2, dFDE);
	} else {
	   WM.fPutCoef(1, 2, -dFDE*dCoef);
	   WM.fPutCoef(2, 2, dFDE*dCoef);
	}
	
	return WorkMat;
     };
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal /* dCoef */ ,
				    const VectorHandler& /* XCurr */ ,
				    const VectorHandler& /* XPrimeCurr */ ) {
      DEBUGCOUT("Entering GenelSpring::AssRes()" << endl);
      
      WorkVec.Resize(2);
      WorkVec.Reset(0.);

      integer iNode1RowIndex = SD1.pNode->iGetFirstRowIndex()+1;
      integer iNode2RowIndex = SD2.pNode->iGetFirstRowIndex()+1;    
      
      doublereal dVal1 = SD1.pNode->dGetDofValue(1, SD1.iOrder);
      doublereal dVal2 = SD2.pNode->dGetDofValue(1, SD2.iOrder);
      
      doublereal d = dVal2-dVal1;
      ConstitutiveLaw1DOwner::Update(d, 0.);
      
      d = GetF();
      
      WorkVec.fPutItem(1, iNode1RowIndex, d);
      WorkVec.fPutItem(2, iNode2RowIndex, -d);

      return WorkVec;
   };
   
   void SetValue(VectorHandler& /* X */ , VectorHandler& /* XP */ ) const {
      NO_OP;
   };
 /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(int& NumNodes, Node::Type* NdTyps, unsigned int* NdLabels) {
     NumNodes = 2;
     NdTyps[0] = SD1.pNode->GetNodeType();
     NdLabels[0] = SD1.pNode->GetLabel();
     NdTyps[1] = SD2.pNode->GetNodeType();
     NdLabels[2] = SD2.pNode->GetLabel();
     
   };
   /* ************************************************ */
};

/* GenelSpring - end */


/* GenelSpringSupport - begin */

class GenelSpringSupport
: virtual public Elem, public Genel, public ConstitutiveLaw1DOwner {
 protected: 
   ScalarDof SD;  
   
 public:
   GenelSpringSupport(unsigned int uLabel, const DofOwner* pDO, 
	       const ConstitutiveLaw1D* pCL,
	       const ScalarDof& sd, flag fOutput)
     : Elem(uLabel, Elem::BULK, fOutput), 
     Genel(uLabel, Genel::SPRINGSUPPORT, pDO, fOutput),
     ConstitutiveLaw1DOwner(pCL), SD(sd) {
      ASSERT(SD.iOrder == 0);
   };
   
   virtual ~GenelSpringSupport(void) { 
      NO_OP;
   };
   
   virtual inline void* pGet(void) const { 
      return (void*)this;
   };   
   
   /* Scrive il contributo dell'elemento al file di restart */
   virtual ostream& Restart(ostream& out) const {
      return out; 
   };

   /* Tipo di Genel */
   virtual Genel::Type GetGenelType(void) const { 
      return Genel::SPRINGSUPPORT;
   };
   
   /* Dimensioni del workspace */
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
      *piNumRows = 1;
      *piNumCols = 1;
   };
   
   void Output(OutputHandler& /* OH */ ) const {
      NO_OP;
   };
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef, 
	    const VectorHandler& /* XCurr */ ,
	    const VectorHandler& /* XPrimeCurr */ ) {
	DEBUGCOUT("Entering GenelSpringSupport::AssJac()" << endl);
	
	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeInit(1, 1, 0.);
	
	integer iNodeRowIndex = SD.pNode->iGetFirstRowIndex()+1;
	integer iNodeColIndex = SD.pNode->iGetFirstColIndex()+1;
	
	WM.fPutRowIndex(1, iNodeRowIndex);
	WM.fPutColIndex(1, iNodeColIndex);            
	     	              
	WM.fPutCoef(1, 1, GetFDE()*dCoef);
	
	return WorkMat;
     };
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal /* dCoef */ ,
				    const VectorHandler& /* XCurr */ ,
				    const VectorHandler& /* XPrimeCurr */ ) {
      DEBUGCOUT("Entering GenelSpringSupport::AssRes()" << endl);
      
      WorkVec.Resize(1);
      WorkVec.Reset(0.);
      
      integer iNodeRowIndex = SD.pNode->iGetFirstRowIndex()+1;
      
      doublereal dVal = SD.pNode->dGetX();    
      
      ConstitutiveLaw1DOwner::Update(dVal, 0.);
      
      WorkVec.fPutItem(1, iNodeRowIndex, -GetF());  

      return WorkVec;
   };
   
   void SetValue(VectorHandler& /* X */ , VectorHandler& /* XP */ ) const {
      NO_OP;
   };
   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(int& NumNodes, Node::Type* NdTyps, unsigned int* NdLabels) {
     NumNodes = 1;
     NdTyps[0] = SD.pNode->GetNodeType();
     NdLabels[0] = SD.pNode->GetLabel();
   };
   /* ************************************************ */

};

/* GenelSpringSupport - end */

   
/* GenelCrossSpringSupport - begin */

class GenelCrossSpringSupport
: virtual public Elem, public Genel, public ConstitutiveLaw1DOwner {
 protected: 
   ScalarDof SDRow;
   ScalarDof SDCol;
   
 public:
   GenelCrossSpringSupport(unsigned int uLabel, const DofOwner* pDO, 
			   const ConstitutiveLaw1D* pCL,
			   const ScalarDof& sdrow,
			   const ScalarDof& sdcol,
			   flag fOutput)
     : Elem(uLabel, Elem::BULK, fOutput), 
     Genel(uLabel, Genel::CROSSSPRINGSUPPORT, pDO, fOutput),
     ConstitutiveLaw1DOwner(pCL), SDRow(sdrow), SDCol(sdcol) {
      ASSERT(SDCol.iOrder == 0);
   };
   
   virtual ~GenelCrossSpringSupport(void) { 
      NO_OP;
   };
   
   virtual inline void* pGet(void) const { 
      return (void*)this;
   };   
   
   /* Scrive il contributo dell'elemento al file di restart */
   virtual ostream& Restart(ostream& out) const {
      return out; 
   };
   
   /* Tipo di Genel */
   virtual Genel::Type GetGenelType(void) const { 
      return Genel::CROSSSPRINGSUPPORT; 
   };
   
   /* Dimensioni del workspace */
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
      *piNumRows = 1;
      *piNumCols = 1;
   };
   
   void Output(OutputHandler& /* OH */ ) const {
      NO_OP;
   };
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef, 
	    const VectorHandler& /* XCurr */ ,
	    const VectorHandler& /* XPrimeCurr */ ) {
	DEBUGCOUT("Entering GenelCrossSpringSupport::AssJac()" << endl);
	
	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeInit(1, 1, 0.);
	
	integer iNodeRowIndex = SDRow.pNode->iGetFirstRowIndex()+1;
	integer iNodeColIndex = SDCol.pNode->iGetFirstColIndex()+1;
	
	WM.fPutRowIndex(1, iNodeRowIndex);
	WM.fPutColIndex(1, iNodeColIndex);            
	     	              
	WM.fPutCoef(1, 1, GetFDE()*dCoef);
	
	return WorkMat;
     };
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal /* dCoef */ ,
				    const VectorHandler& /* XCurr */ ,
				    const VectorHandler& /* XPrimeCurr */ ) {
      DEBUGCOUT("Entering GenelCrossSpringSupport::AssRes()" << endl);
      
      WorkVec.Resize(1);
      WorkVec.Reset(0.);
      
      integer iNodeRowIndex = SDRow.pNode->iGetFirstRowIndex()+1;
      
      doublereal dVal = SDCol.pNode->dGetX();    
      
      ConstitutiveLaw1DOwner::Update(dVal, 0.);
      
      WorkVec.fPutItem(1, iNodeRowIndex, -GetF());  

      return WorkVec;
   };
   
   void SetValue(VectorHandler& /* X */ , VectorHandler& /* XP */ ) const {
      NO_OP;
   };
   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(int& NumNodes, Node::Type* NdTyps, unsigned int* NdLabels) {
     NumNodes = 2;
     NdTyps[0] = SDRow.pNode->GetNodeType();
     NdLabels[0] = SDRow.pNode->GetLabel();
     NdTyps[1] = SDCol.pNode->GetNodeType();
     NdLabels[2] = SDCol.pNode->GetLabel();

   };
   /* ************************************************ */
};

/* GenelCrossSpringSupport - end */
   

/* GenelCrossSpringDamperSupport - begin */

class GenelCrossSpringDamperSupport
: virtual public Elem, public Genel, public ConstitutiveLaw1DOwner {
 protected: 
   ScalarDof SDRow;
   ScalarDof SDCol;
   
 public:
   GenelCrossSpringDamperSupport(unsigned int uLabel, const DofOwner* pDO, 
				 const ConstitutiveLaw1D* pCL,
				 const ScalarDof& sdrow,
				 const ScalarDof& sdcol,
				 flag fOutput)
     : Elem(uLabel, Elem::BULK, fOutput), 
     Genel(uLabel, Genel::CROSSSPRINGDAMPERSUPPORT, pDO, fOutput),
     ConstitutiveLaw1DOwner(pCL), SDRow(sdrow), SDCol(sdcol) {
      ASSERT(SDCol.iOrder == 0);
   };
   
   virtual ~GenelCrossSpringDamperSupport(void) { 
      NO_OP;
   };
   
   virtual inline void* pGet(void) const { 
      return (void*)this;
   };   
   
   /* Scrive il contributo dell'elemento al file di restart */
   virtual ostream& Restart(ostream& out) const {
      return out; 
   };
   
   /* Tipo di Genel */
   virtual Genel::Type GetGenelType(void) const { 
      return Genel::CROSSSPRINGDAMPERSUPPORT; 
   };
   
   /* Dimensioni del workspace */
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
      *piNumRows = 1;
      *piNumCols = 1;
   };
   
   void Output(OutputHandler& /* OH */ ) const {
      NO_OP;
   };
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef, 
	    const VectorHandler& /* XCurr */ ,
	    const VectorHandler& /* XPrimeCurr */ ) {
	DEBUGCOUT("Entering GenelCrossSpringDamperSupport::AssJac()" << endl);
	
	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeInit(1, 1, 0.);
	
	integer iNodeRowIndex = SDRow.pNode->iGetFirstRowIndex()+1;
	integer iNodeColIndex = SDCol.pNode->iGetFirstColIndex()+1;
	
	WM.fPutRowIndex(1, iNodeRowIndex);
	WM.fPutColIndex(1, iNodeColIndex);            
	     	              
	WM.fPutCoef(1, 1, GetFDE()*dCoef+GetFDEPrime());
	
	return WorkMat;
     };
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal /* dCoef */ ,
				    const VectorHandler& /* XCurr */ ,
				    const VectorHandler& /* XPrimeCurr */ ) {
      DEBUGCOUT("Entering GenelCrossSpringDamperSupport::AssRes()" << endl);
      
      WorkVec.Resize(1);
      WorkVec.Reset(0.);
      
      integer iNodeRowIndex = SDRow.pNode->iGetFirstRowIndex()+1;
      
      doublereal dVal = SDCol.pNode->dGetX();
      doublereal dValPrime = SDCol.pNode->dGetXPrime();
      
      ConstitutiveLaw1DOwner::Update(dVal, dValPrime);
      
      WorkVec.fPutItem(1, iNodeRowIndex, -GetF());  

      return WorkVec;
   };
   
   void SetValue(VectorHandler& /* X */ , VectorHandler& /* XP */ ) const {
      NO_OP;
   };
   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(int& NumNodes, Node::Type* NdTyps, unsigned int* NdLabels) {
     NumNodes = 2;
     NdTyps[0] = SDRow.pNode->GetNodeType();
     NdLabels[0] = SDRow.pNode->GetLabel();
     NdTyps[1] = SDCol.pNode->GetNodeType();
     NdLabels[2] = SDCol.pNode->GetLabel();
     
   };
   /* ************************************************ */
};

/* GenelCrossSpringDamperSupport - end */
   

/* GenelSpringDamperSupport - begin */

class GenelSpringDamperSupport
: virtual public Elem, public Genel, public ConstitutiveLaw1DOwner {
 protected: 
   ScalarDof SD;  
   
 public:
   GenelSpringDamperSupport(unsigned int uLabel, const DofOwner* pDO, 
			    const ConstitutiveLaw1D* pCL,
			    const ScalarDof& sd, flag fOutput)
     : Elem(uLabel, Elem::BULK, fOutput), 
     Genel(uLabel, Genel::SPRINGDAMPERSUPPORT, pDO, fOutput),
     ConstitutiveLaw1DOwner(pCL), SD(sd) {
	ASSERT(sd.pNode->SetDof(0) == DofOrder::DIFFERENTIAL);
	ASSERT(sd.iOrder == 0);
   };
   
   virtual ~GenelSpringDamperSupport(void) { 
      NO_OP;
   };
   
   virtual inline void* pGet(void) const { 
      return (void*)this;
   };   
   
   /* Scrive il contributo dell'elemento al file di restart */
   virtual ostream& Restart(ostream& out) const {
      return out; 
   };
   
   /* Tipo di Genel */
   virtual Genel::Type GetGenelType(void) const { 
      return Genel::SPRINGDAMPERSUPPORT; 
   };
   
   /* Dimensioni del workspace */
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
      *piNumRows = 1;
      *piNumCols = 1;
   };
   
   void Output(OutputHandler& /* OH */ ) const {
      NO_OP;
   };
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef, 
	    const VectorHandler& /* XCurr */ ,
	    const VectorHandler& /* XPrimeCurr */ ) {
	DEBUGCOUT("Entering GenelSpringDamperSupport::AssJac()" << endl);
	
	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeInit(1, 1, 0.);
	
	integer iNodeRowIndex = SD.pNode->iGetFirstRowIndex()+1;
	integer iNodeColIndex = SD.pNode->iGetFirstColIndex()+1;
	
	WM.fPutRowIndex(1, iNodeRowIndex);
	WM.fPutColIndex(1, iNodeColIndex);            
	     	              
	WM.fPutCoef(1, 1, GetFDE()*dCoef+GetFDEPrime());
	
	return WorkMat;
     };
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal /* dCoef */ ,
				    const VectorHandler& /* XCurr */ ,
				    const VectorHandler& /* XPrimeCurr */ ) {
      DEBUGCOUT("Entering GenelSpringDamperSupport::AssRes()" << endl);
      
      WorkVec.Resize(1);
      WorkVec.Reset(0.);
      
      integer iNodeRowIndex = SD.pNode->iGetFirstRowIndex()+1;
      
      doublereal dVal = SD.pNode->dGetX();
      doublereal dValPrime = SD.pNode->dGetXPrime();
      
      ConstitutiveLaw1DOwner::Update(dVal, dValPrime);    
      
      WorkVec.fPutItem(1, iNodeRowIndex, -GetF());  

      return WorkVec;
   };
   
   void SetValue(VectorHandler& /* X */ , VectorHandler& /* XP */ ) const {
      NO_OP;
   };
   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(int& NumNodes, Node::Type* NdTyps, unsigned int* NdLabels) {
     NumNodes = 1;
     NdTyps[0] = SD.pNode->GetNodeType();
     NdLabels[0] = SD.pNode->GetLabel();
   };
   /* ************************************************ */
};

/* GenelSpringDamperSupport - end */


/* GenelMass - begin */

class GenelMass : virtual public Elem, public Genel, public DriveOwner {
 protected: 
   ScalarDof SD;
   
 public:
   GenelMass(unsigned int uLabel, const DofOwner* pDO, const DriveCaller* pDC,
	      const ScalarDof& sd, flag fOutput)
     : Elem(uLabel, Elem::BULK, fOutput), 
     Genel(uLabel, Genel::MASS, pDO, fOutput),
     DriveOwner(pDC), SD(sd) { 
      NO_OP;
   };
   
   virtual ~GenelMass(void) { 
      NO_OP;
   };
   
   virtual inline void* pGet(void) const { 
      return (void*)this;
   };
   
   virtual unsigned int iGetNumDof(void) const { 
      return 1;
   };   
      
   /* esegue operazioni sui dof di proprieta' dell'elemento */
#ifdef DEBUG   
   virtual DofOrder::Order SetDof(unsigned int i ) const {
#else /* DEBUG */
   virtual DofOrder::Order SetDof(unsigned int /* i */ ) const {
#endif /* DEBUG */
      ASSERT(i == 0);
      return DofOrder::DIFFERENTIAL;
   };
   
   /* Scrive il contributo dell'elemento al file di restart */
   virtual ostream& Restart(ostream& out) const {
      return out; 
   };
   
   /* Tipo di Genel */
   virtual Genel::Type GetGenelType(void) const { 
      return Genel::MASS; 
   };
   
   /* Dimensioni del workspace */
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
      *piNumRows = 2;
      *piNumCols = 2;
   };
   
   void Output(OutputHandler& /* OH */ ) const {
      NO_OP;
   };
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef, 
	    const VectorHandler& /* XCurr */ ,
	    const VectorHandler& /* XPrimeCurr */ ) {

	DEBUGCOUT("Entering GenelMass::AssJac()" << endl);
	
	SparseSubMatrixHandler& WM = WorkMat.SetSparse();
	WM.ResizeInit(3, 0, 0.);
	
	integer iRowIndex = SD.pNode->iGetFirstRowIndex()+1;
	integer iColIndex = SD.pNode->iGetFirstColIndex()+1;
	integer iDerivativeIndex = iGetFirstIndex()+1;
	
	WM.fPutItem(1, iRowIndex, iDerivativeIndex, dGet());
	WM.fPutItem(2, iDerivativeIndex, iColIndex, -1.);
	WM.fPutItem(3, iDerivativeIndex, iDerivativeIndex, dCoef);
	
	return WorkMat;
     };
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal /* dCoef */ ,
				    const VectorHandler& XCurr,
				    const VectorHandler& XPrimeCurr) {
      DEBUGCOUT("Entering GenelMass::AssRes()" << endl);
      
      WorkVec.Resize(2);
      WorkVec.Reset(0.);

      integer iRowIndex = SD.pNode->iGetFirstRowIndex()+1;
      integer iDerivativeIndex = iGetFirstIndex()+1;
      
      doublereal dVal = SD.pNode->dGetXPrime();
      doublereal dDer = XCurr.dGetCoef(iDerivativeIndex);
      doublereal dDerPrime = XPrimeCurr.dGetCoef(iDerivativeIndex);
      
      WorkVec.fPutItem(1, iRowIndex, -dGet()*dDerPrime);      
      WorkVec.fPutItem(2, iDerivativeIndex, dVal-dDer);

      return WorkVec;
   };
   
   void SetValue(VectorHandler& X, VectorHandler& /* XP */ ) const {
      X.fPutCoef(iGetFirstIndex()+1, SD.pNode->dGetXPrime());
   };
   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(int& NumNodes, Node::Type* NdTyps, unsigned int* NdLabels) {
     NumNodes = 1;
     NdTyps[0] = SD.pNode->GetNodeType();
     NdLabels[0] = SD.pNode->GetLabel();
   };
   /* ************************************************ */
};

/* GenelMass - end */

#endif /* GENEL__H */

