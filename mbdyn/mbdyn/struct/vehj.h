/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2005
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


#ifndef VEHJ_H
#define VEHJ_H

#include "joint.h"
#include "constltp.h"

extern const char* psConstLawNames[];

/* DeformableHingeJoint - begin */

class DeformableHingeJoint : 
virtual public Elem, public Joint, public ConstitutiveLaw3DOwner {
protected:
	const StructNode* pNode1;
	const StructNode* pNode2;
	const Mat3x3 R1h;
	const Mat3x3 R2h;

	bool bFirstRes;

public:
	/* Costruttore non banale */
	DeformableHingeJoint(unsigned int uL,
			const DofOwner* pDO,
			const ConstitutiveLaw3D* pCL,
			const StructNode* pN1, 
			const StructNode* pN2,
			const Mat3x3& R1,
			const Mat3x3& R2, 
			flag fOut);

	/* Distruttore */
	virtual ~DeformableHingeJoint(void);

	/* Tipo di Joint */
	virtual Joint::Type GetJointType(void) const {
		return Joint::DEFORMABLEHINGE; 
	};

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	virtual void Output(OutputHandler& OH) const;

	/* Tipo di DeformableHinge */
	virtual ConstLawType::Type GetConstLawType(void) const = 0;

	virtual unsigned int iGetNumDof(void) const { 
		return 0;
	};

	virtual void
	WorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
		*piNumRows = 6; 
		*piNumCols = 6; 
	};
 
	/* funzioni usate nell'assemblaggio iniziale */
	
	virtual unsigned int iGetInitialNumDof(void) const { 
		return 0;
	};

	/* *******PER IL SOLUTORE PARALLELO******** */        
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void GetConnectedNodes(int& NumNodes, Node::Type* NdTyps,
			unsigned int* NdLabels) {
		NumNodes = 2;
		NdTyps[0] = pNode1->GetNodeType();
		NdLabels[0] = pNode1->GetLabel();
		NdTyps[1] = pNode2->GetNodeType();
		NdLabels[1] = pNode2->GetLabel();
	};
	/* ************************************************ */

	virtual unsigned int iGetNumPrivData(void) const;
	virtual unsigned int iGetPrivDataIdx(const char *s) const;
	virtual doublereal dGetPrivData(unsigned int i) const;
};

/* DeformableHingeJoint - end */


/* ElasticHingeJoint - begin */

class ElasticHingeJoint : virtual public Elem, public DeformableHingeJoint {   
protected:
	Vec3 ThetaRef;
	Vec3 ThetaCurr;

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

	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	/* Tipo di DeformableHinge */
	virtual ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::ELASTIC; 
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
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);   

	virtual unsigned int iGetNumPrivData(void) const {
		return DeformableHingeJoint::iGetNumPrivData();
	};
	virtual unsigned int iGetPrivDataIdx(const char *s) const {
		return DeformableHingeJoint::iGetPrivDataIdx(s);
	};
	virtual doublereal dGetPrivData(unsigned int i) const {
		return DeformableHingeJoint::dGetPrivData(i);
	};
};

/* ElasticHingeJoint - end */


/* ViscousHingeJoint - begin */

class ViscousHingeJoint : virtual public Elem, public DeformableHingeJoint {
protected:
	Vec3 ThetaCurrPrime;

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

	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	/* Tipo di DeformableHinge */
	virtual ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOUS; 
	};

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
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);   

	virtual unsigned int iGetNumPrivData(void) const {
		return DeformableHingeJoint::iGetNumPrivData();
	};
	virtual unsigned int iGetPrivDataIdx(const char *s) const {
		return DeformableHingeJoint::iGetPrivDataIdx(s);
	};
	virtual doublereal dGetPrivData(unsigned int i) const {
		return DeformableHingeJoint::dGetPrivData(i);
	};
};

/* ViscousHingeJoint - end */


/* ViscoElasticHingeJoint - begin */

class ViscoElasticHingeJoint 
: virtual public Elem, public DeformableHingeJoint {
protected:
	Vec3 ThetaRef;
	Vec3 ThetaCurr;

	Vec3 ThetaCurrPrime;

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

	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	/* Tipo di DeformableHinge */
	virtual ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOELASTIC; 
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
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);

	virtual unsigned int iGetNumPrivData(void) const {
		return DeformableHingeJoint::iGetNumPrivData();
	};
	virtual unsigned int iGetPrivDataIdx(const char *s) const {
		return DeformableHingeJoint::iGetPrivDataIdx(s);
	};
	virtual doublereal dGetPrivData(unsigned int i) const {
		return DeformableHingeJoint::dGetPrivData(i);
	};
};

/* ViscoElasticHingeJoint - end */

#endif /* VEHJ_H */

