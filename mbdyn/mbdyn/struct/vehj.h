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

/* Deformable hinges */


#ifndef VEHJ_H
#define VEHJ_H

#include "joint.h"
#include "constltp.h"


/* Tipi di cerniere deformabili (ora sta in constltp.h)
class DefHingeType {
 public:
   enum Type {
      UNKNOWN = -1,
	ELASTIC = 0,
	VISCOUS,
	VISCOELASTIC,
	
	LASTDEFHINGETYPE
   };
};
 */

extern const char* psDefHingeNames[];


/* DeformableHingeJoint - begin */

class DeformableHingeJoint : 
virtual public Elem, public Joint, public ConstitutiveLaw3DOwner {
 private:   
   DefHingeType::Type DefHingeT;      
   
 protected:
   const StructNode* pNode1;
   const StructNode* pNode2;
   const Mat3x3 R1h;
   const Mat3x3 R2h;
   
   
   
   flag fFirstRes;
   
 public:
   /* Costruttore non banale */
   DeformableHingeJoint(unsigned int uL,
			DefHingeType::Type T,
			const DofOwner* pDO,
			const ConstitutiveLaw3D* pCL,
			const StructNode* pN1, 
			const StructNode* pN2,
			const Mat3x3& R1,
			const Mat3x3& R2, 
			flag fOut);

   /* Distruttore */
   virtual ~DeformableHingeJoint(void);
   
   /* Tipo di Joint 
   virtual JointType::Type GetJointType(void) const {
      return JointType::DEFORMABLEHINGE; 
   };
    */
   
   /* Contributo al file di restart */
   virtual ostream& Restart(ostream& out) const;

   virtual void Output(OutputHandler& OH) const;
   
   /* Tipo di DeformableHinge */
   virtual DefHingeType::Type GetDefHingeType(void) const {
      return DefHingeT;
   };
   
   virtual unsigned int iGetNumDof(void) const { 
      return 0;
   };
   
   virtual DofOrder::Order SetDof(unsigned int /* i */ ) const {
      return DofOrder::UNKNOWN;
   };

   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
      *piNumRows = 6; 
      *piNumCols = 6; 
   };
         
   /* funzioni usate nell'assemblaggio iniziale */
   
   virtual unsigned int iGetInitialNumDof(void) const { 
      return 0;
   };
   
#ifdef DEBUG
   virtual const char* sClassName(void) const { 
      return "DeformableHingeJoint";
   };
#endif   

   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(int& NumNodes, NodeType::Type* NdTyps, unsigned int* NdLabels) {
     NumNodes = 2;
     NdTyps[0] = pNode1->GetNodeType();
     NdLabels[0] = pNode1->GetLabel();
     NdTyps[1] = pNode2->GetNodeType();
     NdLabels[1] = pNode2->GetLabel();
   };
   /* ************************************************ */

};

/* DeformableHingeJoint - end */


/* ElasticHingeJoint - begin */

class ElasticHingeJoint : virtual public Elem, public DeformableHingeJoint {   
 protected:
   Vec3 ThetaRef;
   Vec3 ThetaCurr;
   
   Vec3 TaCurr;
   Vec3 TbCurr;
   
   Mat3x3 FDE;
   
   void AssMat(FullSubMatrixHandler& WM, doublereal dCoef);
   void AssVec(SubVectorHandler& WorkVec);
   
 public:
   ElasticHingeJoint(unsigned int uL, 
		     const DofOwner* pDO, 
		     const ConstitutiveLaw3D* pCL,
		     const StructNode* pN1, 
		     const StructNode* pN2,
		     const Mat3x3& R1,
		     const Mat3x3& R2,
		     flag fOut);
   
   ~ElasticHingeJoint(void);
   
   virtual inline void* pGet(void) const { 
      return (void*)this;
   };

   /* Tipo di DeformableHinge
   virtual DefHingeType::Type GetDefHingeType(void) const {
      return DefHingeType::ELASTIC; 
   };
    */
   
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

   /* Aggiorna le deformazioni ecc. */
   virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);           
   
   virtual void InitialWorkSpaceDim(integer* piNumRows,
				    integer* piNumCols) const { 
      *piNumRows = 6; 
      *piNumCols = 6;
   };

   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   virtual VariableSubMatrixHandler& 
     InitialAssJac(VariableSubMatrixHandler& WorkMat, 
		   const VectorHandler& XCurr);
   
   /* Contributo al residuo durante l'assemblaggio iniziale */   
   virtual SubVectorHandler& 
     InitialAssRes(SubVectorHandler& WorkVec,
		   const VectorHandler& XCurr);   

#ifdef DEBUG
   virtual const char* sClassName(void) const { 
      return "ElasticHingeJoint";
   };
#endif   
};

/* ElasticHingeJoint - end */


/* ViscousHingeJoint - begin */

class ViscousHingeJoint : virtual public Elem, public DeformableHingeJoint {
 protected:
   Vec3 ThetaRefPrime;
   Vec3 ThetaCurrPrime;
   
   Vec3 TaCurrPrime;
   Vec3 TbCurrPrime;
   
   Mat3x3 FDEPrime;
   
 public:
   ViscousHingeJoint(unsigned int uL, 
		     const DofOwner* pDO, 
		     const ConstitutiveLaw3D* pCL,
		     const StructNode* pN1, 
		     const StructNode* pN2,
		     const Mat3x3& R1,
		     const Mat3x3& R2,
		     flag fOut);
   
   ~ViscousHingeJoint(void);
   
   virtual inline void* pGet(void) const { 
      return (void*)this;
   };

   /* Tipo di DeformableHinge 
   virtual DefHingeType::Type GetDefHingeType(void) const {
      return DefHingeType::VISCOUS; 
   };
    */

   /* Aggiorna le deformazioni ecc. */
   virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);
   
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

   virtual void InitialWorkSpaceDim(integer* piNumRows,
				    integer* piNumCols) const  { 
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

#ifdef DEBUG
   virtual const char* sClassName(void) const { 
      return "ViscousHingeJoint";
   };
#endif   
};

/* ViscousHingeJoint - end */


/* ViscoElasticHingeJoint - begin */

class ViscoElasticHingeJoint 
: virtual public Elem, public DeformableHingeJoint {
 protected:
   Vec3 ThetaRef;
   Vec3 ThetaCurr;
   
   Vec3 ThetaRefPrime;
   Vec3 ThetaCurrPrime;
   
   Vec3 TaCurr;
   Vec3 TbCurr;
   
   Vec3 TaCurrPrime;
   Vec3 TbCurrPrime;
   
   Mat3x3 FDE;
   Mat3x3 FDEPrime;
   
 public:
   ViscoElasticHingeJoint(unsigned int uL, 
			  const DofOwner* pDO, 
			  const ConstitutiveLaw3D* pCL,
			  const StructNode* pN1, 
			  const StructNode* pN2,
			  const Mat3x3& R1,
			  const Mat3x3& R2, 
			  flag fOut);
   
   ~ViscoElasticHingeJoint(void);
   
   virtual inline void* pGet(void) const { 
      return (void*)this;
   };

   /* Tipo di DeformableHinge
   virtual DefHingeType::Type GetDefHingeType(void) const {
      return DefHingeType::VISCOELASTIC; 
   };
    */
   
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

   /* Aggiorna le deformazioni ecc. */
   virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);           

   virtual void InitialWorkSpaceDim(integer* piNumRows,
				    integer* piNumCols) const {
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

#ifdef DEBUG
   virtual const char* sClassName(void) const { 
      return "ViscoElasticHingeJoint";
   };
#endif   
};

/* ViscoElasticHingeJoint - end */

#endif
