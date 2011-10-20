/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2011
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


#ifndef VEHJ4_H
#define VEHJ4_H

#include "joint.h"
#include "constltp.h"

extern const char* psConstLawNames[];

/* DeformableAxialJoint - begin */

class DeformableAxialJoint :
virtual public Elem, public Joint, public ConstitutiveLaw1DOwner {
protected:
	const StructNode* pNode1;
	const StructNode* pNode2;
	mutable Mat3x3 tilde_R1h;
	mutable Mat3x3 tilde_R2h;

	bool bFirstRes;

	doublereal dTol;

	// from constitutive law
	doublereal dM;

	doublereal dMDE;
	doublereal dMDEPrime;

	// after aggregation
	Vec3 M;

	Mat3x3 MDE;
	Mat3x3 MDEPrime;

	/* Jacobian matrix helpers */
	virtual void
	AssMatM(FullSubMatrixHandler& WMA, doublereal dCoef);

	void
	AssMatMDE(FullSubMatrixHandler& WMA, doublereal dCoef);

	virtual void
	AssMatMDEPrime(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB, doublereal dCoef);

	virtual void AfterPredict(void) = 0;

public:
	/* Costruttore non banale */
	DeformableAxialJoint(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw1D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Mat3x3& tilde_R1h,
		const Mat3x3& tilde_R2h,
		flag fOut);

	/* Distruttore */
	virtual ~DeformableAxialJoint(void);

	/* Tipo di Joint */
	virtual Joint::Type GetJointType(void) const {
		return Joint::DEFORMABLEAXIALJOINT;
	};

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	virtual void Output(OutputHandler& OH) const;

	/* Aggiorna le deformazioni ecc. */
	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);

	void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0);
	virtual void SetInitialValue(VectorHandler& /* X */ );

	/* inverse dynamics capable element */
	virtual bool bInverseDynamics(void) const;
 
	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const;
	         
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
	virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		connectedNodes.resize(2);
		connectedNodes[0] = pNode1;
		connectedNodes[1] = pNode2;
	};
	/* ************************************************ */

	virtual unsigned int iGetNumPrivData(void) const;
	virtual unsigned int iGetPrivDataIdx(const char *s) const;
	virtual doublereal dGetPrivData(unsigned int i) const;
};

/* DeformableAxialJoint - end */


/* ElasticAxialJoint - begin */

class ElasticAxialJoint : virtual public Elem, public DeformableAxialJoint {
protected:
	doublereal dThetaRef;
	doublereal dThetaCurr;

	virtual void AfterPredict(void);
	virtual void AssMat(FullSubMatrixHandler& WM, doublereal dCoef);
	virtual void AssVec(SubVectorHandler& WorkVec);

public:
	ElasticAxialJoint(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw1D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Mat3x3& tilde_R1h,
		const Mat3x3& tilde_R2h,
		flag fOut);

	virtual ~ElasticAxialJoint(void);

	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	/* Tipo di DeformableHinge */
	virtual ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::ELASTIC;
	};

#ifdef MBDYN_X_WORKAROUND_GCC_3_2
	virtual void SetValue(DataManager *pDM,
			VectorHandler& X, VectorHandler& XP,
			SimulationEntity::Hints *ph = 0)
	{
		DeformableAxialJoint::SetValue(pDM, X, XP, ph);
	};

	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const
	{
		return DeformableAxialJoint::ParseHint(pDM, s);
	};
#endif /* MBDYN_X_WORKAROUND_GCC_3_2 */

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
		const VectorHandler& XPrimeCurr,
		const VectorHandler& XPrimePrimeCurr,
		InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS);

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

#if 0
	virtual unsigned int iGetNumPrivData(void) const {
		return DeformableAxialJoint::iGetNumPrivData();
	};

	virtual unsigned int iGetPrivDataIdx(const char *s) const {
		return DeformableAxialJoint::iGetPrivDataIdx(s);
	};

	virtual doublereal dGetPrivData(unsigned int i) const {
		return DeformableAxialJoint::dGetPrivData(i);
	};
#endif
};

/* ElasticAxialJoint - end */


/* ViscousAxialJoint - begin */

class ViscousAxialJoint : virtual public Elem, public DeformableAxialJoint {
protected:
	doublereal dOmega;

	virtual void AfterPredict(void);
	virtual void AssMats(FullSubMatrixHandler& WMA,
			FullSubMatrixHandler& WMB,
			doublereal dCoef);
	virtual void AssVec(SubVectorHandler& WorkVec);

public:
	ViscousAxialJoint(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw1D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Mat3x3& tilde_R1h,
		const Mat3x3& tilde_R2h,
		flag fOut);

	virtual ~ViscousAxialJoint(void);

	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	/* Tipo di DeformableHinge */
	virtual ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOUS;
	};

#ifdef MBDYN_X_WORKAROUND_GCC_3_2
	virtual void SetValue(DataManager *pDM,
			VectorHandler& X, VectorHandler& XP,
			SimulationEntity::Hints *ph = 0)
	{
		DeformableAxialJoint::SetValue(pDM, X, XP, ph);
	};

	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const
	{
		return DeformableAxialJoint::ParseHint(pDM, s);
	};
#endif /* MBDYN_X_WORKAROUND_GCC_3_2 */

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

#if 0
	virtual unsigned int iGetNumPrivData(void) const {
		return DeformableAxialJoint::iGetNumPrivData();
	};

	virtual unsigned int iGetPrivDataIdx(const char *s) const {
		return DeformableAxialJoint::iGetPrivDataIdx(s);
	};

	virtual doublereal dGetPrivData(unsigned int i) const {
		return DeformableAxialJoint::dGetPrivData(i);
	};
#endif
};

/* ViscousAxialJoint - end */


/* ViscoElasticAxialJoint - begin */

class ViscoElasticAxialJoint
: virtual public Elem, public DeformableAxialJoint {
protected:
	doublereal dThetaRef;
	doublereal dThetaCurr;

	doublereal dOmega;

	virtual void AssMats(FullSubMatrixHandler& WMA,
			FullSubMatrixHandler& WMB,
			doublereal dCoef);

	virtual void AfterPredict(void);
	virtual void AssVec(SubVectorHandler& WorkVec);

public:
	ViscoElasticAxialJoint(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw1D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Mat3x3& tilde_R1h,
		const Mat3x3& tilde_R2h,
		flag fOut);

	~ViscoElasticAxialJoint(void);

	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	/* Tipo di DeformableHinge */
	virtual ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOELASTIC;
	};

#ifdef MBDYN_X_WORKAROUND_GCC_3_2
	virtual void SetValue(DataManager *pDM,
			VectorHandler& X, VectorHandler& XP,
			SimulationEntity::Hints *ph = 0)
	{
		DeformableAxialJoint::SetValue(pDM, X, XP, ph);
	};

	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const
	{
		return DeformableAxialJoint::ParseHint(pDM, s);
	};
#endif /* MBDYN_X_WORKAROUND_GCC_3_2 */

	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
			doublereal dCoef,
			const VectorHandler& XCurr,
			const VectorHandler& XPrimeCurr);

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

#if 0
	virtual unsigned int iGetNumPrivData(void) const {
		return DeformableAxialJoint::iGetNumPrivData();
	};

	virtual unsigned int iGetPrivDataIdx(const char *s) const {
		return DeformableAxialJoint::iGetPrivDataIdx(s);
	};

	virtual doublereal dGetPrivData(unsigned int i) const {
		return DeformableAxialJoint::dGetPrivData(i);
	};
#endif
};

/* ViscoElasticAxialJoint - end */

#endif // VEHJ4_H

