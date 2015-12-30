/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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


/* DeformableDispJoint - begin */

class DeformableDispJoint :
virtual public Elem, public Joint, public ConstitutiveLaw3DOwner {
protected:
	const StructNode* pNode1;
	const StructNode* pNode2;
	mutable Vec3 tilde_f1;
	mutable Vec3 tilde_f2;
	mutable Mat3x3 tilde_R1h;
	mutable Mat3x3 tilde_R2h;
	mutable Vec3 tilde_R1hT_tilde_f1;

	Vec3 tilde_d;
	Vec3 tilde_dPrime;

	bool bFirstRes;

	Vec3 F;

	Mat3x3 FDE;
	Mat3x3 FDEPrime;

	void
	AssMatF(FullSubMatrixHandler& WMA,
		const Vec3& d1, const Vec3& d2, doublereal dCoef);
	void
	AssMatFDE(FullSubMatrixHandler& WMA,
		const Vec3& d1, const Vec3& d2, doublereal dCoef);
	void
	AssMatFDEPrime(FullSubMatrixHandler& WMA, FullSubMatrixHandler& WMB,
		const Vec3& d1, const Vec3& d2, doublereal dCoef);

	virtual void
	AssMats(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB,
		doublereal dCoef) = 0;

public:
	/* Costruttore non banale */
	DeformableDispJoint(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw3D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Vec3& tilde_f1,
		const Vec3& tilde_f2,
		const Mat3x3& tilde_R1,
		const Mat3x3& tilde_R2,
		flag fOut);

	/* Distruttore */
	virtual ~DeformableDispJoint(void);

	/* Tipo di Joint */
	virtual Joint::Type GetJointType(void) const {
		return Joint::DEFORMABLEDISPJOINT;
	};

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	void OutputPrepare(OutputHandler &OH);
	virtual void Output(OutputHandler& OH) const;

	void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0);

	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const;
	         
	/* Tipo di DeformableDispHinge */
	virtual ConstLawType::Type GetConstLawType(void) const = 0;

	virtual unsigned int iGetNumDof(void) const {
		return 0;
	};

	virtual DofOrder::Order GetDofType(unsigned int /* i */ ) const {
		return DofOrder::UNKNOWN;
	};

	virtual void
	WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 12;
		*piNumCols = 12;
	};

	/* inverse dynamics capable element */
	virtual bool bInverseDynamics(void) const;
 
	/* funzioni usate nell'assemblaggio iniziale */
	virtual unsigned int iGetInitialNumDof(void) const {
		return 0;
	};

	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* assemblaggio jacobiano */
	virtual void
	AssMats(VariableSubMatrixHandler& WorkMatA,
		VariableSubMatrixHandler& WorkMatB,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	virtual unsigned int iGetNumPrivData(void) const;
	virtual unsigned int iGetPrivDataIdx(const char *s) const;
	virtual doublereal dGetPrivData(unsigned int i) const;

	/* *******PER IL SOLUTORE PARALLELO******** */
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		connectedNodes.resize(2);
		connectedNodes[0] = pNode1;
		connectedNodes[1] = pNode2;
	};
	/* ************************************************ */
};

/* DeformableDispJoint - end */


/* ElasticDispJoint - begin */

class ElasticDispJoint : virtual public Elem, public DeformableDispJoint {
protected:
	void AssMats(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB,
		doublereal dCoef);
	void AssVec(SubVectorHandler& WorkVec);

public:
	ElasticDispJoint(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw3D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Vec3& tilde_f1,
		const Vec3& tilde_f2,
		const Mat3x3& tilde_R1,
		const Mat3x3& tilde_R2,
		flag fOut);

	~ElasticDispJoint(void);

	/* Tipo di DeformableDispHinge */
	virtual ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::ELASTIC;
	};

	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	/* Aggiorna le deformazioni ecc. */
	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);

	/* assemblaggio jacobiano */
	virtual void
	AssMats(VariableSubMatrixHandler& WorkMatA,
		VariableSubMatrixHandler& WorkMatB,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* assemblaggio residuo */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* Inverse Dynamics Jacobian matrix assembly */
	VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	/* Inverse Dynamics residual assembly */
	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr,
		const VectorHandler&  XPrimeCurr,
		const VectorHandler&  XPrimePrimeCurr,
		InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS);

	/* Inverse Dynamics update */
	void Update(const VectorHandler& XCurr, InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS);

	virtual void AfterConvergence(const VectorHandler& X,
			const VectorHandler& XP,
			const VectorHandler& XPP);

	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 12;
		*piNumCols = 12;
	};

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	/* Contributo al residuo durante l'assemblaggio iniziale */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);

#ifdef MBDYN_X_WORKAROUND_GCC_3_2
	virtual void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0)
	{
		DeformableDispJoint::SetValue(pDM, X, XP, ph);
	};

	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const
	{
		return DeformableDispJoint::ParseHint(pDM, s);
	};
	virtual unsigned int iGetNumPrivData(void) const {
		return DeformableDispJoint::iGetNumPrivData();
	};
	virtual unsigned int iGetPrivDataIdx(const char *s) const {
		return DeformableDispJoint::iGetPrivDataIdx(s);
	};
	virtual doublereal dGetPrivData(unsigned int i) const {
		return dGetPrivData(i);
	};
#endif /* MBDYN_X_WORKAROUND_GCC_3_2 */


};

/* ElasticDispJoint - end */


/* ElasticDispJointInv - begin */

class ElasticDispJointInv : virtual public Elem, public DeformableDispJoint {
protected:
	void AssMats(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB,
		doublereal dCoef);
	void AssVec(SubVectorHandler& WorkVec);

public:
	ElasticDispJointInv(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw3D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Vec3& tilde_f1,
		const Vec3& tilde_f2,
		const Mat3x3& tilde_R1,
		const Mat3x3& tilde_R2,
		flag fOut);

	~ElasticDispJointInv(void);

	/* Tipo di DeformableDispHinge */
	virtual ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::ELASTIC;
	};

	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	/* Aggiorna le deformazioni ecc. */
	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);

	/* assemblaggio jacobiano */
	virtual void
	AssMats(VariableSubMatrixHandler& WorkMatA,
		VariableSubMatrixHandler& WorkMatB,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* assemblaggio residuo */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* Inverse Dynamics residual assembly */
	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr,
		const VectorHandler&  XPrimeCurr,
		const VectorHandler&  XPrimePrimeCurr,
		InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS);

	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 12;
		*piNumCols = 12;
	};

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	/* Contributo al residuo durante l'assemblaggio iniziale */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);

#ifdef MBDYN_X_WORKAROUND_GCC_3_2
	virtual void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0)
	{
		DeformableDispJoint::SetValue(pDM, X, XP, ph);
	};

	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const
	{
		return DeformableDispJoint::ParseHint(pDM, s);
	};
	virtual unsigned int iGetNumPrivData(void) const {
		return DeformableDispJoint::iGetNumPrivData();
	};
	virtual unsigned int iGetPrivDataIdx(const char *s) const {
		return DeformableDispJoint::iGetPrivDataIdx(s);
	};
	virtual doublereal dGetPrivData(unsigned int i) const {
		return dGetPrivData(i);
	};
#endif /* MBDYN_X_WORKAROUND_GCC_3_2 */


};

/* ElasticDispJointInv - end */


/* ViscousDispJoint - begin */

class ViscousDispJoint : virtual public Elem, public DeformableDispJoint {
protected:
	void AssMats(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB,
		doublereal dCoef);
	void AssVec(SubVectorHandler& WorkVec);

public:
	ViscousDispJoint(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw3D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Vec3& tilde_f1,
		const Vec3& tilde_f2,
		const Mat3x3& tilde_R1,
		const Mat3x3& tilde_R2,
		flag fOut);

	~ViscousDispJoint(void);

	/* Tipo di DeformableDispHinge */
	virtual ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOUS;
	};

	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	/* Aggiorna le deformazioni ecc. */
	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);

	/* assemblaggio residuo */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* Inverse Dynamics residual assembly */
	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr,
		const VectorHandler&  XPrimeCurr,
		const VectorHandler&  XPrimePrimeCurr,
		InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS);

	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 12;
		*piNumCols = 12;
	};

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	/* Contributo al residuo durante l'assemblaggio iniziale */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);

#ifdef MBDYN_X_WORKAROUND_GCC_3_2
	virtual void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0)
	{
		DeformableDispJoint::SetValue(pDM, X, XP, ph);
	};

	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const
	{
		return DeformableDispJoint::ParseHint(pDM, s);
	};
	virtual unsigned int iGetNumPrivData(void) const {
		return DeformableDispJoint::iGetNumPrivData();
	};
	virtual unsigned int iGetPrivDataIdx(const char *s) const {
		return DeformableDispJoint::iGetPrivDataIdx(s);
	};
	virtual doublereal dGetPrivData(unsigned int i) const {
		return dGetPrivData(i);
	};
#endif /* MBDYN_X_WORKAROUND_GCC_3_2 */

};

/* ViscousDispJoint - end */


/* ViscoElasticDispJoint - begin */

class ViscoElasticDispJoint
: virtual public Elem, public DeformableDispJoint {
protected:
	void AssMats(FullSubMatrixHandler& WorkMatA,
		FullSubMatrixHandler& WorkMatB,
		doublereal dCoef);
	void AssVec(SubVectorHandler& WorkVec);

public:
	ViscoElasticDispJoint(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw3D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Vec3& tilde_f1,
		const Vec3& tilde_f2,
		const Mat3x3& tilde_R1,
		const Mat3x3& tilde_R2,
		flag fOut);

	~ViscoElasticDispJoint(void);

	/* Tipo di DeformableDispHinge */
	virtual ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOELASTIC;
	};

	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	/* Aggiorna le deformazioni ecc. */
	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);

	/* assemblaggio residuo */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* Inverse Dynamics residual assembly */
	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr,
		const VectorHandler&  XPrimeCurr,
		const VectorHandler&  XPrimePrimeCurr,
		InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS);

	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 12;
		*piNumCols = 24;
	};

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	/* Contributo al residuo durante l'assemblaggio iniziale */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);

#ifdef MBDYN_X_WORKAROUND_GCC_3_2
	virtual void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0)
	{
		DeformableDispJoint::SetValue(pDM, X, XP, ph);
	};

	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const
	{
		return DeformableDispJoint::ParseHint(pDM, s);
	};
	virtual unsigned int iGetNumPrivData(void) const {
		return DeformableDispJoint::iGetNumPrivData();
	};
	virtual unsigned int iGetPrivDataIdx(const char *s) const {
		return DeformableDispJoint::iGetPrivDataIdx(s);
	};
	virtual doublereal dGetPrivData(unsigned int i) const {
		return dGetPrivData(i);
	};
#endif /* MBDYN_X_WORKAROUND_GCC_3_2 */
};

/* ViscoElasticDispJoint - end */

#endif /* VEHJ2_H */

