/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2006
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

#ifndef GENEL__H
#define GENEL__H

#include "genel.h"
#include "drive.h"
#include "dataman.h"
#include "constltp.h"

/* GenelClamp - begin */

class GenelClamp : virtual public Elem, public Genel, public DriveOwner {
 protected: 
   ScalarDof SD;
   doublereal dRct;
   
 public:
   GenelClamp(unsigned int uLabel, const DofOwner* pDO, const DriveCaller* pDC,
	      const ScalarDof& sd, flag fOutput)
     : Elem(uLabel, fOutput), 
     Genel(uLabel, pDO, fOutput),
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
   virtual DofOrder::Order GetDofType(unsigned int i ) const {
      ASSERT(i == 0);
      return DofOrder::ALGEBRAIC;
   };
   
   /* esegue operazioni sui dof di proprieta' dell'elemento */
   virtual DofOrder::Order GetEqType(unsigned int i ) const {
      ASSERT(i == 0);
      return DofOrder::ALGEBRAIC;
   };
   
   /* Scrive il contributo dell'elemento al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const {
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
	 std::ostream& out = OH.Genels();
	 out << std::setw(8) << GetLabel() << " " << dRct << std::endl;
      }
   };
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal /* dCoef */ ,
	    const VectorHandler& /* XCurr */ ,
	    const VectorHandler& /* XPrimeCurr */ ) {

	DEBUGCOUT("Entering GenelClamp::AssJac()" << std::endl);
	
	SparseSubMatrixHandler& WM = WorkMat.SetSparse();     
	WM.Resize(2, 0);
	
	integer iRowIndex = SD.pNode->iGetFirstRowIndex()+1;
	integer iColIndex = SD.pNode->iGetFirstColIndex()+1;
	integer iFirstReactionIndex = iGetFirstIndex()+1;
	
	WM.PutItem(1, iRowIndex, iFirstReactionIndex, 1.);
	WM.PutItem(2, iFirstReactionIndex, iColIndex, 1.);
	
	return WorkMat;
     };
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal dCoef,
				    const VectorHandler& XCurr,
				    const VectorHandler& /* XPrimeCurr */ ) {
      DEBUGCOUT("Entering GenelClamp::AssRes()" << std::endl);
      
      WorkVec.Resize(2);  

      integer iRowIndex = SD.pNode->iGetFirstRowIndex()+1;
      integer iFirstReactionIndex = iGetFirstIndex()+1;
      
      doublereal dVal = SD.pNode->dGetDofValue(1, SD.iOrder);
      dRct = XCurr.dGetCoef(iFirstReactionIndex);
      
      WorkVec.PutItem(1, iRowIndex, -dRct);

      doublereal dConstr = dGet()-dVal;
      if (SD.iOrder == 0 
		      && SD.pNode->GetDofType(0) == DofOrder::DIFFERENTIAL
     		      && dCoef != 0.) {
	 dConstr /= dCoef;
      }
      WorkVec.PutItem(2, iFirstReactionIndex, dConstr);

      return WorkVec;
   };
   
   void SetValue(DataManager *pDM,
		   VectorHandler& X, VectorHandler& XP,
		   SimulationEntity::Hints *ph = 0)
   {
      if (SD.iOrder == 0) {
	 X.PutCoef(SD.pNode->iGetFirstRowIndex()+1, dGet());
      } else if (SD.iOrder == 1) {
	 XP.PutCoef(SD.pNode->iGetFirstRowIndex()+1, dGet());	 
      }
   };
   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) {
     connectedNodes.resize(1);
     connectedNodes[0] = SD.pNode;
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
     : Elem(uLabel, fOutput), 
     Genel(uLabel, pDO, fOutput),
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
   virtual DofOrder::Order GetDofType(unsigned int i ) const {
#else     
   virtual DofOrder::Order GetDofType(unsigned int /* i */ ) const {
#endif    
      ASSERT(i == 0);
      return DofOrder::ALGEBRAIC;
   };
   
   /* Scrive il contributo dell'elemento al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const {
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
	 std::ostream& out = OH.Genels();
	 out << std::setw(8) << GetLabel() << " " << dRct << std::endl;
      }
   };
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef, 
	    const VectorHandler& /* XCurr */ ,
	    const VectorHandler& /* XPrimeCurr */ ) {
	DEBUGCOUT("Entering GenelDistance::AssJac()" << std::endl);
	
	SparseSubMatrixHandler& WM = WorkMat.SetSparse();
	WM.ResizeReset(4, 0);
	
	integer iNode1RowIndex = SD1.pNode->iGetFirstRowIndex()+1;
	integer iNode1ColIndex = SD1.pNode->iGetFirstColIndex()+1;
	integer iNode2RowIndex = SD2.pNode->iGetFirstRowIndex()+1;
	integer iNode2ColIndex = SD2.pNode->iGetFirstColIndex()+1;
	integer iFirstReactionIndex = iGetFirstIndex()+1;
	
	WM.PutItem(1, iNode1RowIndex, iFirstReactionIndex, -1.);
	WM.PutItem(2, iNode2RowIndex, iFirstReactionIndex, 1.);
	
	doublereal d = dCoef;
	if ((SD1.iOrder == 0) && (SD2.iOrder == 0)) {
	   d = 1.;
	}
	
	if (SD1.iOrder == 1) {
	   WM.PutItem(3, iFirstReactionIndex, iNode1ColIndex, -1.);
	} else {
	   WM.PutItem(3, iFirstReactionIndex, iNode1ColIndex, -d);
	}      
	
	if (SD2.iOrder == 1) {
	   WM.PutItem(4, iFirstReactionIndex, iNode2ColIndex, 1.);
	} else {
	   WM.PutItem(4, iFirstReactionIndex, iNode2ColIndex, d);
	}      
	
	return WorkMat;
     };
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal dCoef,
				    const VectorHandler& XCurr,
				    const VectorHandler& /* XPrimeCurr */ ) {
      DEBUGCOUT("Entering GenelDistance::AssRes()" << std::endl);
      
      WorkVec.ResizeReset(3);

      integer iNode1RowIndex = SD1.pNode->iGetFirstRowIndex()+1;
      integer iNode2RowIndex = SD2.pNode->iGetFirstRowIndex()+1;
      integer iFirstReactionIndex = iGetFirstIndex()+1;
      
      doublereal dVal1 = SD1.pNode->dGetDofValue(1, SD1.iOrder);
      doublereal dVal2 = SD2.pNode->dGetDofValue(1, SD2.iOrder);
      dRct = XCurr.dGetCoef(iFirstReactionIndex);
      
      WorkVec.PutItem(1, iNode1RowIndex, dRct);
      WorkVec.PutItem(2, iNode2RowIndex, -dRct);
      
      if ((SD1.iOrder == 0) && (SD2.iOrder == 0)) {     
	 if (dCoef != 0.) {
	    WorkVec.PutItem(3, iFirstReactionIndex, 
			     (dGet()-dVal2+dVal1)/dCoef);
	 }	 
      } else {
	 WorkVec.PutItem(3, iFirstReactionIndex, dGet()-dVal2+dVal1);
      }

      return WorkVec;
   };
   
   void SetValue(DataManager *pDM,
		   VectorHandler& X, VectorHandler& XP,
		   SimulationEntity::Hints *ph = 0)
   {
      if (SD2.iOrder == 0) {
	 X.PutCoef(SD2.pNode->iGetFirstRowIndex()+1,
		    dGet()-SD2.pNode->dGetX()
		    +SD1.pNode->dGetDofValue(1, SD1.iOrder));
      } else if (SD2.iOrder == 1) {
	 XP.PutCoef(SD2.pNode->iGetFirstRowIndex()+1, 
		     dGet()-SD2.pNode->dGetXPrime()
		     +SD1.pNode->dGetDofValue(1, SD1.iOrder));
      }
   };
   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) {
     connectedNodes.resize(2);
     connectedNodes[0] = SD1.pNode;
     connectedNodes[1] = SD2.pNode;

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

   doublereal dVal;
   
 public:
   GenelSpring(unsigned int uLabel, const DofOwner* pDO, 
	       const ConstitutiveLaw1D* pCL,
	       const ScalarDof& sd1, const ScalarDof& sd2, flag fOutput)
     : Elem(uLabel, fOutput), 
     Genel(uLabel, pDO, fOutput),
     ConstitutiveLaw1DOwner(pCL), SD1(sd1), SD2(sd2), dVal(0.) { 
      NO_OP;
   };
   
   virtual ~GenelSpring(void) { 
      NO_OP;
   };
   
   virtual inline void* pGet(void) const { 
      return (void*)this;
   };   
   
   /* Scrive il contributo dell'elemento al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const {
      return out; 
   };
  
   virtual void AfterConvergence(const VectorHandler& X, 
			const VectorHandler& XP) {
      ConstitutiveLaw1DOwner::AfterConvergence(dVal, 0.);
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
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef, 
	    const VectorHandler& /* XCurr */ ,
	    const VectorHandler& /* XPrimeCurr */ ) {
	DEBUGCOUT("Entering GenelSpring::AssJac()" << std::endl);
	
	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeReset(2, 2);
	
	integer iNode1RowIndex = SD1.pNode->iGetFirstRowIndex()+1;
	integer iNode1ColIndex = SD1.pNode->iGetFirstColIndex()+1;
	integer iNode2RowIndex = SD2.pNode->iGetFirstRowIndex()+1;
	integer iNode2ColIndex = SD2.pNode->iGetFirstColIndex()+1;
	
	WM.PutRowIndex(1, iNode1RowIndex);
	WM.PutColIndex(1, iNode1ColIndex);
	WM.PutRowIndex(2, iNode2RowIndex);
	WM.PutColIndex(2, iNode2ColIndex);
	
	doublereal dFDE = GetFDE();
      	
	if (SD1.iOrder == 1) {
	   WM.PutCoef(1, 1, dFDE);
	   WM.PutCoef(2, 1, -dFDE);
	} else {
	   WM.PutCoef(1, 1, dFDE*dCoef);
	   WM.PutCoef(2, 1, -dFDE*dCoef);
	}
	
	if (SD2.iOrder == 1) {
	   WM.PutCoef(1, 2, -dFDE);
	   WM.PutCoef(2, 2, dFDE);
	} else {
	   WM.PutCoef(1, 2, -dFDE*dCoef);
	   WM.PutCoef(2, 2, dFDE*dCoef);
	}
	
	return WorkMat;
     };
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal /* dCoef */ ,
				    const VectorHandler& /* XCurr */ ,
				    const VectorHandler& /* XPrimeCurr */ ) {
      DEBUGCOUT("Entering GenelSpring::AssRes()" << std::endl);
      
      WorkVec.ResizeReset(2);

      integer iNode1RowIndex = SD1.pNode->iGetFirstRowIndex()+1;
      integer iNode2RowIndex = SD2.pNode->iGetFirstRowIndex()+1;    
      
      doublereal dVal1 = SD1.pNode->dGetDofValue(1, SD1.iOrder);
      doublereal dVal2 = SD2.pNode->dGetDofValue(1, SD2.iOrder);
      
      dVal = dVal2-dVal1;
      ConstitutiveLaw1DOwner::Update(dVal, 0.);
      
      doublereal d = GetF();
      
      WorkVec.PutItem(1, iNode1RowIndex, d);
      WorkVec.PutItem(2, iNode2RowIndex, -d);

      return WorkVec;
   };
   
 /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) {
     connectedNodes.resize(2);
     connectedNodes[0] = SD1.pNode;
     connectedNodes[1] = SD2.pNode;
   };
   /* ************************************************ */
};

/* GenelSpring - end */


/* GenelSpringSupport - begin */

class GenelSpringSupport
: virtual public Elem, public Genel, public ConstitutiveLaw1DOwner {
 protected: 
   ScalarDof SD;
   doublereal dVal;
   
 public:
   GenelSpringSupport(unsigned int uLabel, const DofOwner* pDO, 
	       const ConstitutiveLaw1D* pCL,
	       const ScalarDof& sd, flag fOutput)
     : Elem(uLabel, fOutput), 
     Genel(uLabel, pDO, fOutput),
     ConstitutiveLaw1DOwner(pCL), SD(sd), dVal(0.) {
      ASSERT(SD.iOrder == 0);
   };
   
   virtual ~GenelSpringSupport(void) { 
      NO_OP;
   };
   
   virtual inline void* pGet(void) const { 
      return (void*)this;
   };   
   
   /* Scrive il contributo dell'elemento al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const {
      return out; 
   };

   virtual void AfterConvergence(const VectorHandler& X, 
			const VectorHandler& XP) {
      ConstitutiveLaw1DOwner::AfterConvergence(dVal, 0.);
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
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef, 
	    const VectorHandler& /* XCurr */ ,
	    const VectorHandler& /* XPrimeCurr */ ) {
	DEBUGCOUT("Entering GenelSpringSupport::AssJac()" << std::endl);
	
	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeReset(1, 1);
	
	integer iNodeRowIndex = SD.pNode->iGetFirstRowIndex()+1;
	integer iNodeColIndex = SD.pNode->iGetFirstColIndex()+1;
	
	WM.PutRowIndex(1, iNodeRowIndex);
	WM.PutColIndex(1, iNodeColIndex);            
	     	              
	WM.PutCoef(1, 1, GetFDE()*dCoef);
	
	return WorkMat;
     };
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal /* dCoef */ ,
				    const VectorHandler& /* XCurr */ ,
				    const VectorHandler& /* XPrimeCurr */ ) {
      DEBUGCOUT("Entering GenelSpringSupport::AssRes()" << std::endl);
      
      WorkVec.ResizeReset(1);
      
      integer iNodeRowIndex = SD.pNode->iGetFirstRowIndex()+1;
      
      dVal = SD.pNode->dGetX();
      ConstitutiveLaw1DOwner::Update(dVal, 0.);
      
      WorkVec.PutItem(1, iNodeRowIndex, -GetF());  

      return WorkVec;
   };
   
   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) {
     connectedNodes.resize(1);
     connectedNodes[0] = SD.pNode;
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
   doublereal dVal;
   
 public:
   GenelCrossSpringSupport(unsigned int uLabel, const DofOwner* pDO, 
			   const ConstitutiveLaw1D* pCL,
			   const ScalarDof& sdrow,
			   const ScalarDof& sdcol,
			   flag fOutput)
     : Elem(uLabel, fOutput), 
     Genel(uLabel, pDO, fOutput),
     ConstitutiveLaw1DOwner(pCL), SDRow(sdrow), SDCol(sdcol), dVal(0.) {
      ASSERT(SDCol.iOrder == 0);
   };
   
   virtual ~GenelCrossSpringSupport(void) { 
      NO_OP;
   };
   
   virtual inline void* pGet(void) const { 
      return (void*)this;
   };   
   
   /* Scrive il contributo dell'elemento al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const {
      return out; 
   };
   
   virtual void AfterConvergence(const VectorHandler& X, 
			const VectorHandler& XP) {
      ConstitutiveLaw1DOwner::AfterConvergence(dVal, 0.);
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
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef, 
	    const VectorHandler& /* XCurr */ ,
	    const VectorHandler& /* XPrimeCurr */ ) {
	DEBUGCOUT("Entering GenelCrossSpringSupport::AssJac()" << std::endl);
	
	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeReset(1, 1);
	
	integer iNodeRowIndex = SDRow.pNode->iGetFirstRowIndex()+1;
	integer iNodeColIndex = SDCol.pNode->iGetFirstColIndex()+1;
	
	WM.PutRowIndex(1, iNodeRowIndex);
	WM.PutColIndex(1, iNodeColIndex);            
	     	              
	WM.PutCoef(1, 1, GetFDE()*dCoef);
	
	return WorkMat;
     };
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal /* dCoef */ ,
				    const VectorHandler& /* XCurr */ ,
				    const VectorHandler& /* XPrimeCurr */ ) {
      DEBUGCOUT("Entering GenelCrossSpringSupport::AssRes()" << std::endl);
      
      WorkVec.ResizeReset(1);
      
      integer iNodeRowIndex = SDRow.pNode->iGetFirstRowIndex()+1;
      
      dVal = SDCol.pNode->dGetX();    
      ConstitutiveLaw1DOwner::Update(dVal, 0.);
      
      WorkVec.PutItem(1, iNodeRowIndex, -GetF());  

      return WorkVec;
   };
   
   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) {
     connectedNodes.resize(2);
     connectedNodes[0] = SDRow.pNode;
     connectedNodes[1] = SDCol.pNode;
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
   doublereal dVal;
   doublereal dValPrime;
   
 public:
   GenelCrossSpringDamperSupport(unsigned int uLabel, const DofOwner* pDO, 
				 const ConstitutiveLaw1D* pCL,
				 const ScalarDof& sdrow,
				 const ScalarDof& sdcol,
				 flag fOutput)
     : Elem(uLabel, fOutput), 
     Genel(uLabel, pDO, fOutput),
     ConstitutiveLaw1DOwner(pCL), SDRow(sdrow), SDCol(sdcol),
     dVal(0.), dValPrime(0.) {
      ASSERT(SDCol.iOrder == 0);
   };
   
   virtual ~GenelCrossSpringDamperSupport(void) { 
      NO_OP;
   };
   
   virtual inline void* pGet(void) const { 
      return (void*)this;
   };   
   
   /* Scrive il contributo dell'elemento al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const {
      return out; 
   };
   
   virtual void AfterConvergence(const VectorHandler& X, 
			const VectorHandler& XP) {
      ConstitutiveLaw1DOwner::AfterConvergence(dVal, dValPrime);
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
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef, 
	    const VectorHandler& /* XCurr */ ,
	    const VectorHandler& /* XPrimeCurr */ ) {
	DEBUGCOUT("Entering GenelCrossSpringDamperSupport::AssJac()" << std::endl);
	
	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeReset(1, 1);
	
	integer iNodeRowIndex = SDRow.pNode->iGetFirstRowIndex()+1;
	integer iNodeColIndex = SDCol.pNode->iGetFirstColIndex()+1;
	
	WM.PutRowIndex(1, iNodeRowIndex);
	WM.PutColIndex(1, iNodeColIndex);            
	     	              
	WM.PutCoef(1, 1, GetFDE()*dCoef+GetFDEPrime());
	
	return WorkMat;
     };
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal /* dCoef */ ,
				    const VectorHandler& /* XCurr */ ,
				    const VectorHandler& /* XPrimeCurr */ ) {
      DEBUGCOUT("Entering GenelCrossSpringDamperSupport::AssRes()" << std::endl);
      
      WorkVec.ResizeReset(1);
      
      integer iNodeRowIndex = SDRow.pNode->iGetFirstRowIndex()+1;
      
      dVal = SDCol.pNode->dGetX();
      dValPrime = SDCol.pNode->dGetXPrime();
      ConstitutiveLaw1DOwner::Update(dVal, dValPrime);
      
      WorkVec.PutItem(1, iNodeRowIndex, -GetF());  

      return WorkVec;
   };
   
   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) {
     connectedNodes.resize(2);
     connectedNodes[0] = SDRow.pNode;
     connectedNodes[1] = SDCol.pNode;
   };
   /* ************************************************ */
};

/* GenelCrossSpringDamperSupport - end */
   

/* GenelSpringDamperSupport - begin */

class GenelSpringDamperSupport
: virtual public Elem, public Genel, public ConstitutiveLaw1DOwner {
 protected: 
   ScalarDof SD;  
   doublereal dVal;
   doublereal dValPrime;
   
 public:
   GenelSpringDamperSupport(unsigned int uLabel, const DofOwner* pDO, 
			    const ConstitutiveLaw1D* pCL,
			    const ScalarDof& sd, flag fOutput)
     : Elem(uLabel, fOutput), 
     Genel(uLabel, pDO, fOutput),
     ConstitutiveLaw1DOwner(pCL), SD(sd), dVal(0.),dValPrime(0.) {
	ASSERT(sd.pNode->GetDofType(0) == DofOrder::DIFFERENTIAL);
	ASSERT(sd.iOrder == 0);
   };
   
   virtual ~GenelSpringDamperSupport(void) { 
      NO_OP;
   };
   
   virtual inline void* pGet(void) const { 
      return (void*)this;
   };   
   
   /* Scrive il contributo dell'elemento al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const {
      return out; 
   };
   
   virtual void AfterConvergence(const VectorHandler& X, 
			const VectorHandler& XP) {
      ConstitutiveLaw1DOwner::AfterConvergence(dVal, dValPrime);
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
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef, 
	    const VectorHandler& /* XCurr */ ,
	    const VectorHandler& /* XPrimeCurr */ ) {
	DEBUGCOUT("Entering GenelSpringDamperSupport::AssJac()" << std::endl);
	
	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeReset(1, 1);
	
	integer iNodeRowIndex = SD.pNode->iGetFirstRowIndex()+1;
	integer iNodeColIndex = SD.pNode->iGetFirstColIndex()+1;
	
	WM.PutRowIndex(1, iNodeRowIndex);
	WM.PutColIndex(1, iNodeColIndex);            
	     	              
	WM.PutCoef(1, 1, GetFDE()*dCoef+GetFDEPrime());
	
	return WorkMat;
     };
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal /* dCoef */ ,
				    const VectorHandler& /* XCurr */ ,
				    const VectorHandler& /* XPrimeCurr */ ) {
      DEBUGCOUT("Entering GenelSpringDamperSupport::AssRes()" << std::endl);
      
      WorkVec.ResizeReset(1);
      
      integer iNodeRowIndex = SD.pNode->iGetFirstRowIndex()+1;
      
      dVal = SD.pNode->dGetX();
      dValPrime = SD.pNode->dGetXPrime();
      
      ConstitutiveLaw1DOwner::Update(dVal, dValPrime);    
      
      WorkVec.PutItem(1, iNodeRowIndex, -GetF());  

      return WorkVec;
   };
   
   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) {
     connectedNodes.resize(1);
     connectedNodes[0] = SD.pNode;
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
     : Elem(uLabel, fOutput), 
     Genel(uLabel, pDO, fOutput),
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
   virtual DofOrder::Order GetDofType(unsigned int i ) const {
#else /* DEBUG */
   virtual DofOrder::Order GetDofType(unsigned int /* i */ ) const {
#endif /* DEBUG */
      ASSERT(i == 0);
      return DofOrder::DIFFERENTIAL;
   };
   
   /* Scrive il contributo dell'elemento al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const {
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
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef, 
	    const VectorHandler& /* XCurr */ ,
	    const VectorHandler& /* XPrimeCurr */ ) {

	DEBUGCOUT("Entering GenelMass::AssJac()" << std::endl);
	
	SparseSubMatrixHandler& WM = WorkMat.SetSparse();
	WM.ResizeReset(3, 0);
	
	integer iRowIndex = SD.pNode->iGetFirstRowIndex()+1;
	integer iColIndex = SD.pNode->iGetFirstColIndex()+1;
	integer iDerivativeIndex = iGetFirstIndex()+1;
	
	WM.PutItem(1, iRowIndex, iDerivativeIndex, dGet());
	WM.PutItem(2, iDerivativeIndex, iColIndex, -1.);
	WM.PutItem(3, iDerivativeIndex, iDerivativeIndex, dCoef);
	
	return WorkMat;
     };
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal /* dCoef */ ,
				    const VectorHandler& XCurr,
				    const VectorHandler& XPrimeCurr) {
      DEBUGCOUT("Entering GenelMass::AssRes()" << std::endl);
      
      WorkVec.ResizeReset(2);

      integer iRowIndex = SD.pNode->iGetFirstRowIndex()+1;
      integer iDerivativeIndex = iGetFirstIndex()+1;
      
      doublereal dVal = SD.pNode->dGetXPrime();
      doublereal dDer = XCurr.dGetCoef(iDerivativeIndex);
      doublereal dDerPrime = XPrimeCurr.dGetCoef(iDerivativeIndex);
      
      WorkVec.PutItem(1, iRowIndex, -dGet()*dDerPrime);      
      WorkVec.PutItem(2, iDerivativeIndex, dVal-dDer);

      return WorkVec;
   };
   
   void SetValue(DataManager *pDM,
		   VectorHandler& X, VectorHandler& /* XP */ ,
		   SimulationEntity::Hints *ph = 0)
   {
      X.PutCoef(iGetFirstIndex()+1, SD.pNode->dGetXPrime());
   };
   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) {
     connectedNodes.resize(1);
     connectedNodes[0] = SD.pNode;
   };
   /* ************************************************ */
};

/* GenelMass - end */

#endif /* GENEL__H */

