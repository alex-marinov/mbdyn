/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

/* Bulk elements */

#ifndef BULK_H
#define BULK_H

#include "elem.h"
#include "node.h"
#include "drive.h"


/* Tipi di bulk */
class BulkType {
 public:
   enum Type {
      UNKNOWN = -1,
	SPRINGSUPPORT = 0,
	SPRING,
	
	LASTBULKTYPE
   };
};

extern const char* psBulkNames[];


/* Bulk - begin */

class Bulk : virtual public Elem {
 public:
   Bulk(unsigned int uLabel, flag fOutput)
     : Elem(uLabel, fOutput) { 
	NO_OP;
     };
   
   virtual ~Bulk(void) { 
      NO_OP;
   };
   
   virtual Elem::Type GetElemType(void) const {
      return Elem::BULK;
   };
};

/* Bulk - end */


/* BulkSpringSupport - begin */

class BulkSpringSupport 
: virtual public Elem, public Bulk, public DriveOwner {
 protected: 
   ScalarDof SD;
   
 public:
   BulkSpringSupport(unsigned int uLabel, const DriveCaller* pDC,
		     const ScalarDof& sd, flag fOutput)
     : Elem(uLabel, fOutput), Bulk(uLabel, fOutput),
     DriveOwner(pDC), SD(sd) { 
      NO_OP;
   };
   
   virtual ~BulkSpringSupport(void) { 
      NO_OP;
   };
   
   /* Scrive il contributo dell'elemento al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const {
      return out; 
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
	
	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeReset(1, 1);
	
	integer iRowIndex = SD.pNode->iGetFirstRowIndex()+1;
	integer iColIndex = SD.pNode->iGetFirstColIndex()+1;
	WM.PutRowIndex(1, iRowIndex);
	WM.PutColIndex(1, iColIndex);
	
	doublereal d = dGet();
	if (SD.iOrder == 0) {
	   d *= dCoef;
	}
	WM.PutCoef(1, 1, d);
	
	return WorkMat;
     };
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal /* dCoef */ ,
				    const VectorHandler& /* XCurr */ ,
				    const VectorHandler& /* XPrimeCurr */ ) {
      WorkVec.Resize(1);
      WorkVec.Reset();

      integer iRowIndex = SD.pNode->iGetFirstRowIndex()+1;
      doublereal dVal = SD.pNode->dGetDofValue(1, SD.iOrder);
      WorkVec.PutItem(1, iRowIndex, -dGet()*dVal);

      return WorkVec;
   };

   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
     connectedNodes.resize(1);
     connectedNodes[0] = SD.pNode;
   };
   /* ************************************************ */
};

/* BulkSpringSupport - end */

class DataManager;
class MBDynParser;

extern Elem* ReadBulk(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);

#endif
