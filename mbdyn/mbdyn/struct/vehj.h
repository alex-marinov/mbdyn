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
	mutable Mat3x3 tilde_R1h;
	mutable Mat3x3 tilde_R2h;

	OrientationDescription od;

private:
#ifdef USE_NETCDF
	NcVar *Var_Phi;
	NcVar *Var_Omega;
#endif // USE_NETCDF

protected:
	bool bFirstRes;

	Vec3 M;

	Mat3x3 MDE;
	Mat3x3 MDEPrime;

	/* for invariant stuff */
	Mat3x3	hat_I;
	Mat3x3	hat_IT;

	/* Jacobian matrix helpers */
	virtual void
	AssMatM(FullSubMatrixHandler& WMA, doublereal dCoef);

	void
	AssMatMInv(FullSubMatrixHandler& WMA, doublereal dCoef);

	void
	AssMatMDE(FullSubMatrixHandler& WMA, doublereal dCoef);

	virtual void
	AssMatMDEPrime(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB, doublereal dCoef);

	void
	AssMatMDEPrimeInv(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB, doublereal dCoef);

	/* output helper */
	void OutputInv(OutputHandler& OH) const;

	/* priv data helper */
	doublereal
	dGetPrivDataInv(unsigned int i) const;
	virtual void AfterPredict(void) = 0;

public:
	/* Costruttore non banale */
	DeformableHingeJoint(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw3D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Mat3x3& tilde_R1h,
		const Mat3x3& tilde_R2h,
		const OrientationDescription& od,
		flag fOut);

	/* Distruttore */
	virtual ~DeformableHingeJoint(void);

	/* Tipo di Joint */
	virtual Joint::Type GetJointType(void) const {
		return Joint::DEFORMABLEHINGE;
	};

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	void OutputPrepare(OutputHandler &OH);
	virtual void Output(OutputHandler& OH) const;

	/* Aggiorna le deformazioni ecc. */
	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);

	void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0);
	virtual void SetInitialValue(VectorHandler& /* X */ );

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

	/* inverse dynamics capable element */
	virtual bool bInverseDynamics(void) const;
 
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

/* DeformableHingeJoint - end */


/* ElasticHingeJoint - begin */

class ElasticHingeJoint : virtual public Elem, public DeformableHingeJoint {
protected:
	Vec3 ThetaRef;
	Vec3 ThetaCurr;

	virtual void AfterPredict(void);
	virtual void AssMat(FullSubMatrixHandler& WM, doublereal dCoef);
	virtual void AssVec(SubVectorHandler& WorkVec);

public:
	ElasticHingeJoint(unsigned int uL,
			const DofOwner* pDO,
			const ConstitutiveLaw3D* pCL,
			const StructNode* pN1,
			const StructNode* pN2,
			const Mat3x3& tilde_R1h,
			const Mat3x3& tilde_R2h,
			const OrientationDescription& od,
			flag fOut);

	virtual ~ElasticHingeJoint(void);

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
		DeformableHingeJoint::SetValue(pDM, X, XP, ph);
	};

	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const
	{
		return DeformableHingeJoint::ParseHint(pDM, s);
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

	/* Inverse Dynamics Jacobian matrix assembly */
	VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	/* Inverse Dynamics residual assembly */
	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr,
		const VectorHandler& XPrimePrimeCurr,
		InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS);

	/* Inverse Dynamics update */
	void Update(const VectorHandler& XCurr, InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS);

	virtual void AfterConvergence(const VectorHandler& X,
			const VectorHandler& XP,
			const VectorHandler& XPP);

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
		return DeformableHingeJoint::iGetNumPrivData();
	};

	virtual unsigned int iGetPrivDataIdx(const char *s) const {
		return DeformableHingeJoint::iGetPrivDataIdx(s);
	};

	virtual doublereal dGetPrivData(unsigned int i) const {
		return DeformableHingeJoint::dGetPrivData(i);
	};
#endif
};

/* ElasticHingeJoint - end */


/* ElasticHingeJointInv - begin */

class ElasticHingeJointInv : virtual public Elem, public ElasticHingeJoint {
protected:
	virtual void
	AssMatM(FullSubMatrixHandler& WMA, doublereal dCoef);

	/* AssMatMDE is OK as MDE is updated fine by AfterPredict();
	 * AssMatMDEPrime is not needed */

	virtual void AfterPredict(void);
	virtual void AssVec(SubVectorHandler& WorkVec);

public:
	ElasticHingeJointInv(unsigned int uL,
			const DofOwner* pDO,
			const ConstitutiveLaw3D* pCL,
			const StructNode* pN1,
			const StructNode* pN2,
			const Mat3x3& tilde_R1h,
			const Mat3x3& tilde_R2h,
			const OrientationDescription& od,
			flag fOut);

	virtual ~ElasticHingeJointInv(void);

	virtual void Output(OutputHandler& OH) const;

	virtual doublereal dGetPrivData(unsigned int i) const;

#ifdef MBDYN_X_WORKAROUND_GCC_3_2
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
	{
		ElasticHingeJoint::AfterConvergence(X, XP);
	};
	
	virtual void SetValue(DataManager *pDM,
			VectorHandler& X, VectorHandler& XP,
			SimulationEntity::Hints *ph = 0)
	{
		ElasticHingeJoint::SetValue(pDM, X, XP, ph);
	};

	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const
	{
		return ElasticHingeJoint::ParseHint(pDM, s);
	};

	virtual unsigned int iGetNumPrivData(void) const {
		return DeformableHingeJoint::iGetNumPrivData();
	};

	virtual unsigned int iGetPrivDataIdx(const char *s) const {
		return DeformableHingeJoint::iGetPrivDataIdx(s);
	};
#endif /* MBDYN_X_WORKAROUND_GCC_3_2 */
};

/* ElasticHingeJointInv - end */


/* ViscousHingeJoint - begin */

class ViscousHingeJoint : virtual public Elem, public DeformableHingeJoint {
protected:
	Vec3 Omega;

	virtual void AfterPredict(void);
	virtual void AssMats(FullSubMatrixHandler& WMA,
			FullSubMatrixHandler& WMB,
			doublereal dCoef);
	virtual void AssVec(SubVectorHandler& WorkVec);

public:
	ViscousHingeJoint(unsigned int uL,
			const DofOwner* pDO,
			const ConstitutiveLaw3D* pCL,
			const StructNode* pN1,
			const StructNode* pN2,
			const Mat3x3& tilde_R1h,
			const Mat3x3& tilde_R2h,
			const OrientationDescription& od,
			flag fOut);

	virtual ~ViscousHingeJoint(void);

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
		DeformableHingeJoint::SetValue(pDM, X, XP, ph);
	};

	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const
	{
		return DeformableHingeJoint::ParseHint(pDM, s);
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
		return DeformableHingeJoint::iGetNumPrivData();
	};

	virtual unsigned int iGetPrivDataIdx(const char *s) const {
		return DeformableHingeJoint::iGetPrivDataIdx(s);
	};

	virtual doublereal dGetPrivData(unsigned int i) const {
		return DeformableHingeJoint::dGetPrivData(i);
	};
#endif
};

/* ViscousHingeJoint - end */


/* ViscousHingeJointInv - begin */

class ViscousHingeJointInv : virtual public Elem, public ViscousHingeJoint {
protected:
	/* AssMatMDEPrime is not needed */
	virtual void
	AssMatM(FullSubMatrixHandler& WMA, doublereal dCoef);
	virtual void
	AssMatMDEPrime(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB, doublereal dCoef);

	virtual void AfterPredict(void);
	virtual void AssVec(SubVectorHandler& WorkVec);

public:
	ViscousHingeJointInv(unsigned int uL,
			const DofOwner* pDO,
			const ConstitutiveLaw3D* pCL,
			const StructNode* pN1,
			const StructNode* pN2,
			const Mat3x3& tilde_R1h,
			const Mat3x3& tilde_R2h,
			const OrientationDescription& od,
			flag fOut);

	virtual ~ViscousHingeJointInv(void);

	virtual void Output(OutputHandler& OH) const;

	virtual doublereal dGetPrivData(unsigned int i) const;

#ifdef MBDYN_X_WORKAROUND_GCC_3_2
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
	{
		ViscousHingeJoint::AfterConvergence(X, XP);
	};
	
	virtual void SetValue(DataManager *pDM,
			VectorHandler& X, VectorHandler& XP,
			SimulationEntity::Hints *ph = 0)
	{
		DeformableHingeJoint::SetValue(pDM, X, XP, ph);
	};

	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const
	{
		return DeformableHingeJoint::ParseHint(pDM, s);
	};

	virtual unsigned int iGetNumPrivData(void) const {
		return DeformableHingeJoint::iGetNumPrivData();
	};

	virtual unsigned int iGetPrivDataIdx(const char *s) const {
		return DeformableHingeJoint::iGetPrivDataIdx(s);
	};
#endif /* MBDYN_X_WORKAROUND_GCC_3_2 */
};

/* ViscousHingeJointInv - end */


/* ViscoElasticHingeJoint - begin */

class ViscoElasticHingeJoint
: virtual public Elem, public DeformableHingeJoint {
protected:
	Vec3 ThetaRef;
	Vec3 ThetaCurr;

	Vec3 Omega;

	virtual void AssMats(FullSubMatrixHandler& WMA,
			FullSubMatrixHandler& WMB,
			doublereal dCoef);

	virtual void AfterPredict(void);
	virtual void AssVec(SubVectorHandler& WorkVec);

public:
	ViscoElasticHingeJoint(unsigned int uL,
			const DofOwner* pDO,
			const ConstitutiveLaw3D* pCL,
			const StructNode* pN1,
			const StructNode* pN2,
			const Mat3x3& tilde_R1h,
			const Mat3x3& tilde_R2h,
			const OrientationDescription& od,
			flag fOut);

	~ViscoElasticHingeJoint(void);

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
		DeformableHingeJoint::SetValue(pDM, X, XP, ph);
	};

	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const
	{
		return DeformableHingeJoint::ParseHint(pDM, s);
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
		return DeformableHingeJoint::iGetNumPrivData();
	};

	virtual unsigned int iGetPrivDataIdx(const char *s) const {
		return DeformableHingeJoint::iGetPrivDataIdx(s);
	};

	virtual doublereal dGetPrivData(unsigned int i) const {
		return DeformableHingeJoint::dGetPrivData(i);
	};
#endif
};

/* ViscoElasticHingeJoint - end */

/* ViscoElasticHingeJointInv - begin */

class ViscoElasticHingeJointInv
: virtual public Elem, public ViscoElasticHingeJoint {
protected:
	/* AssMatMDEPrime is not needed */
	virtual void
	AssMatM(FullSubMatrixHandler& WMA, doublereal dCoef);
	virtual void
	AssMatMDEPrime(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB, doublereal dCoef);

	virtual void AfterPredict(void);
	virtual void AssVec(SubVectorHandler& WorkVec);

public:
	ViscoElasticHingeJointInv(unsigned int uL,
			const DofOwner* pDO,
			const ConstitutiveLaw3D* pCL,
			const StructNode* pN1,
			const StructNode* pN2,
			const Mat3x3& tilde_R1h,
			const Mat3x3& tilde_R2h,
			const OrientationDescription& od,
			flag fOut);

	~ViscoElasticHingeJointInv(void);

	virtual void Output(OutputHandler& OH) const;

	virtual doublereal dGetPrivData(unsigned int i) const;

#ifdef MBDYN_X_WORKAROUND_GCC_3_2
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
	{
		ViscoElasticHingeJoint::AfterConvergence(X, XP);
	};
	
	virtual void SetValue(DataManager *pDM,
			VectorHandler& X, VectorHandler& XP,
			SimulationEntity::Hints *ph = 0)
	{
		DeformableHingeJoint::SetValue(pDM, X, XP, ph);
	};

	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const
	{
		return DeformableHingeJoint::ParseHint(pDM, s);
	};

	virtual unsigned int iGetNumPrivData(void) const {
		return DeformableHingeJoint::iGetNumPrivData();
	};

	virtual unsigned int iGetPrivDataIdx(const char *s) const {
		return DeformableHingeJoint::iGetPrivDataIdx(s);
	};
#endif /* MBDYN_X_WORKAROUND_GCC_3_2 */
};

/* ViscoElasticHingeJointInv - end */


/* InvAngularCLR - begin */

struct InvAngularCLR : public ConstitutiveLawRead<Vec3, Mat3x3> {
	virtual ConstitutiveLaw<Vec3, Mat3x3> *
	Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType);
};

/* InvAngularCLR - end */

#endif /* VEHJ_H */

