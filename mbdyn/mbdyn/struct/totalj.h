/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2007
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

#ifndef TOTALJ_H
#define TOTALJ_H

#include "joint.h"
#include "constltp.h"

#if 0
// friction will be dealt with later
#include "friction.h"
#endif

/* TotalJoint - begin */

class TotalJoint :
virtual public Elem, public Joint {
private:
	const StructNode* pNode1;
	const StructNode* pNode2;
	Vec3 f1;
	Mat3x3 R1h;
	Mat3x3 R1hr;
	Vec3 f2;
	Mat3x3 R2h;
	Mat3x3 R2hr;
	bool bPosActive[3];
	bool bRotActive[3];
	
	TplDriveOwner<Vec3> XDrv;
	TplDriveOwner<Vec3> XPDrv;
	TplDriveOwner<Vec3> XPPDrv;
	
	TplDriveOwner<Vec3> ThetaDrv;
	TplDriveOwner<Vec3> OmegaDrv;
	TplDriveOwner<Vec3> OmegaPDrv;
	
	unsigned int nConstraints;
	unsigned int nPosConstraints;
	unsigned int nRotConstraints;
	
	unsigned int iPosIncid[3];
	unsigned int iRotIncid[3];
	
	Vec3 tilde_f1;

	mutable Vec3 M;
	mutable Vec3 F;
	mutable Vec3 ThetaDelta;
	mutable Vec3 ThetaDeltaPrev;

public:
	/* Constructor */
	TotalJoint(unsigned int uL, const DofOwner *pDO,
		bool bPos[3],
		TplDriveCaller<Vec3> *const pDCPos[3],
		bool bRot[3],
		TplDriveCaller<Vec3> *const pDCRot[3],
		const StructNode* pN1,
		const Vec3& f1Tmp, const Mat3x3& R1hTmp, const Mat3x3& R1hrTmp, 
		const StructNode* pN2,
		const Vec3& f2Tmp, const Mat3x3& R2hTmp, const Mat3x3& R2hrTmp, 
		flag fOut);

	/* Destructor */
	~TotalJoint(void);

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	/* Tipo di Joint */
	virtual Joint::Type
	GetJointType(void) const {
		return Joint::TOTALJOINT;
	};

	virtual unsigned int
	iGetNumDof(void) const {
		return nConstraints;
	};

	virtual std::ostream&
	DescribeDof(std::ostream& out,
		char *prefix = "",
		bool bInitial = false, int i = -1) const;

	virtual std::ostream&
	DescribeEq(std::ostream& out,
		char *prefix = "",
		bool bInitial = false, int i = -1) const;

	DofOrder::Order
	GetDofType(unsigned int i) const {
		ASSERT(i >= 0 && i < nConstraints);
		return DofOrder::ALGEBRAIC;
	};

	virtual void
	SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0);

	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const;

	virtual void
	AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP);

	void
	WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumCols = *piNumRows = 12 + nConstraints ;
	};

	VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* inverse dynamics Jacobian matrix assembly */
	VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	/* inverse dynamics residual assembly */
	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr,
		const VectorHandler&  XPrimeCurr,
		const VectorHandler&  XPrimePrimeCurr,
		int iOrder = -1);

	/* Inverse Dynamics update */
	virtual void Update(const VectorHandler& XCurr, int iOrder = -1);
	
	DofOrder::Order GetEqType(unsigned int i) const;

	void Output(OutputHandler& OH) const;

	/* funzioni usate nell'assemblaggio iniziale */

	virtual unsigned int
	iGetInitialNumDof(void) const {
		return  2*nConstraints;
	};
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumCols = *piNumRows = 24 + 2*nConstraints;
	};

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	/* Contributo al residuo durante l'assemblaggio iniziale */
	SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr);

	/* Dati privati */
	virtual unsigned int iGetNumPrivData(void) const;
	virtual unsigned int iGetPrivDataIdx(const char *s) const;
	virtual doublereal dGetPrivData(unsigned int i) const;

	/* *******PER IL SOLUTORE PARALLELO******** */
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) {
		connectedNodes.resize(2);
		connectedNodes[0] = pNode1;
		connectedNodes[1] = pNode2;
	};
	/* ************************************************ */
};

/* TotalJoint - end */

/* TotalPinJoint - begin */

class TotalPinJoint :
virtual public Elem, public Joint {
private:
	const StructNode* pNode;
	Vec3 Xc;
	Mat3x3 Rch;
	Mat3x3 Rchr;
	Vec3 tilde_fn;
	Mat3x3 tilde_Rnh;
	Mat3x3 tilde_Rnhr;
	bool bPosActive[3];
	bool bRotActive[3];
	
	Mat3x3 RchT;
	Vec3 tilde_Xc;
	Mat3x3 RchrT;

	TplDriveOwner<Vec3> XDrv;
	TplDriveOwner<Vec3> XPDrv;
	TplDriveOwner<Vec3> XPPDrv;
	
	TplDriveOwner<Vec3> ThetaDrv;
	TplDriveOwner<Vec3> OmegaDrv;
	TplDriveOwner<Vec3> OmegaPDrv;
	
	unsigned int nConstraints;
	unsigned int nPosConstraints;
	unsigned int nRotConstraints;
	
	unsigned int iPosIncid[3];
	unsigned int iRotIncid[3];
	
	mutable Vec3 M;
	mutable Vec3 F;
	mutable Vec3 ThetaDelta;
	mutable Vec3 ThetaDeltaPrev;

public:
	/* Constructor */
	TotalPinJoint(unsigned int uL, const DofOwner *pDO,
		bool bPos[3],
		TplDriveCaller<Vec3> *const pDCPos[3],
		bool bRot[3],
		TplDriveCaller<Vec3> *const pDCRot[3],
		const Vec3& XcTmp, const Mat3x3& RchTmp, const Mat3x3& RchrTmp,
		const StructNode* pN,
		const Vec3& fnTmp, const Mat3x3& RnhTmp, const Mat3x3& RnhrTmp, 
		flag fOut);

	/* Destructor */
	~TotalPinJoint(void);

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	/* Tipo di Joint */
	virtual Joint::Type
	GetJointType(void) const {
		return Joint::TOTALPINJOINT;
	};

	virtual unsigned int
	iGetNumDof(void) const {
		return nConstraints;
	};

	virtual std::ostream&
	DescribeDof(std::ostream& out,
		char *prefix = "",
		bool bInitial = false, int i = -1) const;

	virtual std::ostream&
	DescribeEq(std::ostream& out,
		char *prefix = "",
		bool bInitial = false, int i = -1) const;

	DofOrder::Order
	GetDofType(unsigned int i) const {
		ASSERT(i >= 0 && i < nConstraints);
		return DofOrder::ALGEBRAIC;
	};

	virtual void
	SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0);

	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const;

	virtual void
	AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP);

	void
	WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumCols = *piNumRows = 6 + nConstraints ;
	};

	VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

#if 0
	/* inverse dynamics Jacobian matrix assembly */
	VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	/* inverse dynamics residual assembly */
	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr,
		const VectorHandler&  XPrimeCurr,
		const VectorHandler&  XPrimePrimeCurr,
		int iOrder = -1);

	/* Inverse Dynamics update */
	virtual void Update(const VectorHandler& XCurr, int iOrder = -1);
#endif
	
	DofOrder::Order GetEqType(unsigned int i) const;

	void Output(OutputHandler& OH) const;

	/* funzioni usate nell'assemblaggio iniziale */

	virtual unsigned int
	iGetInitialNumDof(void) const {
		return  2*nConstraints;
	};
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumCols = *piNumRows = 12 + 2*nConstraints;
	};

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	/* Contributo al residuo durante l'assemblaggio iniziale */
	SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr);

	/* Dati privati */
	virtual unsigned int iGetNumPrivData(void) const;
	virtual unsigned int iGetPrivDataIdx(const char *s) const;
	virtual doublereal dGetPrivData(unsigned int i) const;

	/* *******PER IL SOLUTORE PARALLELO******** */
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) {
		connectedNodes.resize(1);
		connectedNodes[0] = pNode;
	};
	/* ************************************************ */
};

/* TotalPinJoint - end */

#endif // TOTALJ_H
