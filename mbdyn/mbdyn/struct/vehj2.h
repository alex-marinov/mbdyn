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

/* Deformable hinges */

#ifndef VEHJ2_H
#define VEHJ2_H

#include "joint.h"
#include "constltp.h"


/* DeformableDispHingeJoint - begin */

class DeformableDispHingeJoint : 
virtual public Elem, public Joint, public ConstitutiveLaw3DOwner {
 protected:
   const StructNode* pNode1;
   const StructNode* pNode2;
   const Vec3 f1;
   const Vec3 f2;
   const Mat3x3 R1h;
   const Mat3x3 R2h;      
   
   Vec3 k;
   Vec3 kPrime;

 public:
   /* Costruttore non banale */
   DeformableDispHingeJoint(unsigned int uL,
			    const DofOwner* pDO,
			    const ConstitutiveLaw3D* pCL,
			    const StructNode* pN1, 
			    const StructNode* pN2,
			    const Vec3& f1Tmp,
			    const Vec3& f2Tmp,
			    const Mat3x3& R1,
			    const Mat3x3& R2,
			    flag fOut);
   
   /* Distruttore */
   virtual ~DeformableDispHingeJoint(void);
   
   /* Tipo di Joint */
   virtual Joint::Type GetJointType(void) const {
      return Joint::DEFORMABLEHINGE; 
   };
   
   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;

   virtual void Output(OutputHandler& OH) const;
   
   /* Tipo di DeformableDispHinge */
   virtual DefHingeType::Type GetDefHingeType(void) const = 0;
   
   virtual unsigned int iGetNumDof(void) const { 
      return 0;
   };
   
   virtual DofOrder::Order SetDof(unsigned int /* i */ ) const {
      return DofOrder::UNKNOWN;
   };

   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
      *piNumRows = 12; 
      *piNumCols = 12; 
   };
         
   /* funzioni usate nell'assemblaggio iniziale */
   
   virtual unsigned int iGetInitialNumDof(void) const { 
      return 0;
   };
   
   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(int& NumNodes, Node::Type* NdTyps, unsigned int* NdLabels) {
     NumNodes = 2;
     NdTyps[0] = pNode1->GetNodeType();
     NdLabels[0] = pNode1->GetLabel();
     NdTyps[1] = pNode2->GetNodeType();
     NdLabels[1] = pNode2->GetLabel();
   };
   /* ************************************************ */  
};

/* DeformableDispHingeJoint - end */


/* ElasticDispHingeJoint - begin */

class ElasticDispHingeJoint : virtual public Elem, public DeformableDispHingeJoint {
 protected:
   void AssMat(FullSubMatrixHandler& WM, doublereal dCoef);
   void AssVec(SubVectorHandler& WorkVec);
   
 public:
   ElasticDispHingeJoint(unsigned int uL, 
			 const DofOwner* pDO, 
			 const ConstitutiveLaw3D* pCL,
			 const StructNode* pN1, 
			 const StructNode* pN2,
			 const Vec3& f1Tmp,
			 const Vec3& f2Tmp,
			 const Mat3x3& R1,
			 const Mat3x3& R2,
			 flag fOut);
   
   ~ElasticDispHingeJoint(void);
   
   virtual inline void* pGet(void) const { 
      return (void*)this;
   };

   /* Tipo di DeformableDispHinge */
   virtual DefHingeType::Type GetDefHingeType(void) const {
      return DefHingeType::ELASTIC; 
   };
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef, 
	    const VectorHandler& XCurr,
	    const VectorHandler& XPrimeCurr);
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& 
     AssRes(SubVectorHandler& WorkVec,
	    doublereal dCoef,
	    const VectorHandler& XCurr, 
	    const VectorHandler& XPrimeCurr);

   virtual void InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
      *piNumRows = 12; 
      *piNumCols = 12;
   };
   
   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   virtual VariableSubMatrixHandler& 
     InitialAssJac(VariableSubMatrixHandler& WorkMat, 
		   const VectorHandler& XCurr);
   
   /* Contributo al residuo durante l'assemblaggio iniziale */   
   virtual SubVectorHandler& 
     InitialAssRes(SubVectorHandler& WorkVec,
		   const VectorHandler& XCurr);   
};

/* ElasticDispHingeJoint - end */


/* ViscousDispHingeJoint - begin */

class ViscousDispHingeJoint : virtual public Elem, public DeformableDispHingeJoint {
 public:
   ViscousDispHingeJoint(unsigned int uL, 
			 const DofOwner* pDO, 
			 const ConstitutiveLaw3D* pCL,
			 const StructNode* pN1, 
			 const StructNode* pN2,
			 const Vec3& f1Tmp,
			 const Vec3& f2Tmp,
			 const Mat3x3& R1,
			 const Mat3x3& R2,
			 flag fOut);
   
   ~ViscousDispHingeJoint(void);
   
   virtual inline void* pGet(void) const { 
      return (void*)this;
   };

   /* Tipo di DeformableDispHinge */
   virtual DefHingeType::Type GetDefHingeType(void) const {
      return DefHingeType::VISCOUS; 
   };
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef, 
	    const VectorHandler& XCurr,
	    const VectorHandler& XPrimeCurr);
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& 
     AssRes(SubVectorHandler& WorkVec,
	    doublereal dCoef,
	    const VectorHandler& XCurr, 
	    const VectorHandler& XPrimeCurr);

   virtual void InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
      *piNumRows = 6; 
      *piNumCols = 12;
   };

   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   virtual VariableSubMatrixHandler& 
     InitialAssJac(VariableSubMatrixHandler& WorkMat, 
		   const VectorHandler& XCurr);
   
   /* Contributo al residuo durante l'assemblaggio iniziale */   
   virtual SubVectorHandler& 
     InitialAssRes(SubVectorHandler& WorkVec,
		   const VectorHandler& XCurr);   
};

/* ViscousDispHingeJoint - end */


/* ViscoElasticDispHingeJoint - begin */

class ViscoElasticDispHingeJoint 
: virtual public Elem, public DeformableDispHingeJoint {
 public:
   ViscoElasticDispHingeJoint(unsigned int uL, 
			      const DofOwner* pDO, 
			      const ConstitutiveLaw3D* pCL,
			      const StructNode* pN1, 
			      const StructNode* pN2,
			      const Vec3& f1Tmp,
			      const Vec3& f2Tmp,
			      const Mat3x3& R1,
			      const Mat3x3& R2, 
			      flag fOut);
   
   ~ViscoElasticDispHingeJoint(void);
   
   virtual inline void* pGet(void) const { return (void*)this; };   

   /* Tipo di DeformableDispHinge */
   virtual DefHingeType::Type GetDefHingeType(void) const {
      return DefHingeType::VISCOELASTIC; 
   };
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef, 
	    const VectorHandler& XCurr,
	    const VectorHandler& XPrimeCurr);
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& 
     AssRes(SubVectorHandler& WorkVec,
	    doublereal dCoef,
	    const VectorHandler& XCurr, 
	    const VectorHandler& XPrimeCurr);

   virtual void InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
      *piNumRows = 6; 
      *piNumCols = 12;
   };

   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   virtual VariableSubMatrixHandler& 
     InitialAssJac(VariableSubMatrixHandler& WorkMat, 
		   const VectorHandler& XCurr);
   
   /* Contributo al residuo durante l'assemblaggio iniziale */   
   virtual SubVectorHandler& 
     InitialAssRes(SubVectorHandler& WorkVec,
		   const VectorHandler& XCurr);   
};

/* ViscoElasticDispHingeJoint - end */

#endif /* VEHJ2_H */

