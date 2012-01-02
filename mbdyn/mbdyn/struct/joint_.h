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

/* vincoli, tipo: Elem::Type JOINT */

#ifndef JOINT__H
#define JOINT__H

#include <cfloat>

#include <strnode.h>
#include <joint.h>


/* Classi di servizio per vincoli.
 * Queste classi sono intese per svincolare l'implementatore 
 * dalla gestione degli indici dei nodi e dei dof interni.
 */

/* Joint_ - begin */

class Joint_
: virtual public Elem, public Joint {
 protected:   
   
 public:
   Joint_(unsigned int uL,
	  Joint::Type T,
	  const DofOwner* pD,
	  flag fOut);
   virtual ~Joint_(void);
   
   virtual void AssRes_(SubVectorHandler& WorkVec,
			doublereal dCoef,
			const VectorHandler& XCurr,
			const VectorHandler& XPrimeCurr) = 0;
   
   virtual void AssJac_(VariableSubMatrixHandler& WorkMat,
			doublereal dCoef,
			const VectorHandler& XCurr,
			const VectorHandler& XPrimeCurr) = 0;

   virtual void InitialAssRes_(SubVectorHandler& WorkVec,
			       const VectorHandler& XCurr) = 0;
   
   virtual void InitialAssJac_(VariableSubMatrixHandler& WorkMat,
			       const VectorHandler& XCurr) = 0;
};

/* Joint_ - end */


/* Joint_1Node - begin */

class Joint_1Node
: virtual public Elem, public Joint_ {
 protected:
   StructNode* pNode;
   
 public:
   Joint_1Node(unsigned int uL,
	       Joint::Type T,
	       const DofOwner* pD,
	       flag fOut);
   virtual ~Joint_1Node(void);
   
   virtual SubVectorHandler& 
     AssRes(SubVectorHandler& WorkVec,
	    doublereal dCoef,
	    const VectorHandler& XCurr,
	    const VectorHandler& XPrimeCurr);
   
   virtual VariableSubMatrixHandler&
     AssJac(VariableSubMatrixHandler& WorkMat,
            doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr);
   
   virtual SubVectorHandler& 
     InitialAssRes(SubVectorHandler& WorkVec,
		   const VectorHandler& XCurr);

   virtual VariableSubMatrixHandler&
     InitialAssJac(VariableSubMatrixHandler& WorkMat,
            const VectorHandler& XCurr);
};

/* Joint_1Node - end */


/* Joint_2Nodes - begin */

class Joint_2Nodes
: virtual public Elem, public Joint_ {
 protected:
   StructNode* pNode1;
   StructNode* pNode2;
   
 public:
   Joint_2Nodes(unsigned int uL,
		Joint::Type T,
		const DofOwner* pD,
		flag fOut);
   virtual ~Joint_2Nodes(void);

   virtual SubVectorHandler& 
     AssRes(SubVectorHandler& WorkVec,
	    doublereal dCoef,
	    const VectorHandler& XCurr,
	    const VectorHandler& XPrimeCurr);

   virtual VariableSubMatrixHandler&
     AssJac(VariableSubMatrixHandler& WorkMat,
            doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr);
   
   virtual SubVectorHandler& 
     InitialAssRes(SubVectorHandler& WorkVec,
		   const VectorHandler& XCurr);
   
   virtual VariableSubMatrixHandler&
     InitialAssJac(VariableSubMatrixHandler& WorkMat,
            const VectorHandler& XCurr);
};

/* Joint_2Nodes - end */


/* Joint_NNodes - begin */

class Joint_NNodes
: virtual public Elem, public Joint {
 protected:
   StructNode** pNodes;      
   
 public:
   Joint_NNodes(unsigned int uL,
		Joint::Type T,
		const DofOwner* pD,
		flag fOut);
   virtual ~Joint_NNodes(void);
   
   virtual int iGetNNodes(void) const = 0;

   virtual SubVectorHandler& 
     AssRes(SubVectorHandler& WorkVec,
	    doublereal dCoef,
	    const VectorHandler& XCurr,
	    const VectorHandler& XPrimeCurr);

   virtual VariableSubMatrixHandler&
     AssJac(VariableSubMatrixHandler& WorkMat,
            doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr);

   virtual SubVectorHandler& 
     InitialAssRes(SubVectorHandler& WorkVec,
		   const VectorHandler& XCurr);

   virtual VariableSubMatrixHandler&
     InitialAssJac(VariableSubMatrixHandler& WorkMat,
            const VectorHandler& XCurr);
};

/* Joint_NNodes - end */

#endif // JOINT__H
