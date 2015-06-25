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

#ifndef TOTALEQUATION_H
#define TOTALEQUATION_H

#include "joint.h"

#if 0
// friction will be dealt with later
#include "friction.h"
#endif

// forward decls
class TotalReaction;


/* TotalEquation - begin */
class TotalEquation :
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
	
	bool bVelActive[3];
	bool bAgvActive[3];	/* Agv stands for AnGular Velocity */
	
	TplDriveOwner<Vec3> XDrv;
	TplDriveOwner<Vec3> XPDrv;
	TplDriveOwner<Vec3> XPPDrv;
	
	TplDriveOwner<Vec3> ThetaDrv;
	TplDriveOwner<Vec3> OmegaDrv;
	TplDriveOwner<Vec3> OmegaPDrv;
	
	unsigned int nConstraints;
	unsigned int nPosConstraints;
	unsigned int nRotConstraints;
	unsigned int nVelConstraints;
	unsigned int nAgvConstraints;
	
	unsigned int iPosIncid[3];
	unsigned int iRotIncid[3];
	unsigned int iVelIncid[3];
	unsigned int iAgvIncid[3];
	
	unsigned int iPosEqIndex[3];
	unsigned int iRotEqIndex[3];
	unsigned int iVelEqIndex[3];
	unsigned int iAgvEqIndex[3];
	
	Vec3 tilde_f1;

	mutable Vec3 M;
	mutable Vec3 F;
	mutable Vec3 ThetaDelta;
	mutable Vec3 ThetaDeltaPrev;
	
	friend class TotalReaction;

public:
	/* Constructor */
	TotalEquation(unsigned int uL, const DofOwner *pDO,
		bool bPos[3], bool bVel[3],
		TplDriveCaller<Vec3> *const pDCPos[3],
		bool bRot[3], bool bAgv[3],
		TplDriveCaller<Vec3> *const pDCRot[3],
		const StructNode* pN1,
		const Vec3& f1Tmp, const Mat3x3& R1hTmp, const Mat3x3& R1hrTmp, 
		const StructNode* pN2,
		const Vec3& f2Tmp, const Mat3x3& R2hTmp, const Mat3x3& R2hrTmp, 
		flag fOut);

	/* Destructor */
	~TotalEquation(void);

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	/* Tipo di Joint */
	virtual Joint::Type
	GetJointType(void) const {
		return Joint::TOTALEQUATION;
	};

	virtual unsigned int
	iGetNumDof(void) const {
		return nConstraints;
	};

	virtual std::ostream&
	DescribeDof(std::ostream& out,
		const char *prefix = "",
		bool bInitial = false) const;

	virtual void
	DescribeDof(std::vector<std::string>& desc,
		bool bInitial = false,
		int i = -1) const;

	virtual std::ostream&
	DescribeEq(std::ostream& out,
		const char *prefix = "",
		bool bInitial = false) const;

	virtual void
	DescribeEq(std::vector<std::string>& desc,
		bool bInitial = false,
		int i = -1) const;

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

	/* Inverse Dynamics: */
	virtual void
	AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP, const VectorHandler& XPP);

	void
	WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = nConstraints ;
		*piNumCols = 12;
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
	
	/* Inverse Dynamics Jacobian matrix assembly 
	VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	 Inverse Dynamics residual assembly 
	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr,
		const VectorHandler&  XPrimeCurr,
		const VectorHandler&  XPrimePrimeCurr,
		InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS);

	 Inverse Dynamics update 
	void Update(const VectorHandler& XCurr, InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS);	
	DofOrder::Order GetEqType(unsigned int i) const;
*/
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
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		connectedNodes.resize(2);
		connectedNodes[0] = pNode1;
		connectedNodes[1] = pNode2;
	};
	/* ************************************************ */
};


/* TotalEquation - end */

/* TotalReaction - begin */

class TotalReaction :
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
	
	TotalEquation * total_equation_element;

public:
	/* Constructor */
	TotalReaction(unsigned int uL, const DofOwner *pDO,
		bool bPos[3],
		bool bRot[3],
		const StructNode* pN1,
		const Vec3& f1Tmp, const Mat3x3& R1hTmp, const Mat3x3& R1hrTmp, 
		const StructNode* pN2,
		const Vec3& f2Tmp, const Mat3x3& R2hTmp, const Mat3x3& R2hrTmp,
		TotalEquation* t_elm, 
		flag fOut);

	/* Destructor */
	~TotalReaction(void);

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	/* Tipo di Joint */
	virtual Joint::Type
	GetJointType(void) const {
		return Joint::TOTALREACTION;
	};

	virtual unsigned int
	iGetNumDof(void) const {
		return 0;
	};

	virtual std::ostream&
	DescribeDof(std::ostream& out,
		const char *prefix = "",
		bool bInitial = false) const;

	virtual void
	DescribeDof(std::vector<std::string>& desc,
		bool bInitial = false,
		int i = -1) const;

#if 0
	virtual std::ostream&
	DescribeEq(std::ostream& out,
		char *prefix = "",
		bool bInitial = false) const;

	virtual void
	DescribeEq(std::vector<std::string>& desc,
		bool bInitial = false,
		int i = -1) const;
#endif

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

	virtual void
	AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP, const VectorHandler& XPP);

	void
	WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 12;
		*piNumCols = 12 + nConstraints ;
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

	DofOrder::Order GetEqType(unsigned int i) const;

	void Output(OutputHandler& OH) const;

	/* funzioni usate nell'assemblaggio iniziale */

	virtual unsigned int
	iGetInitialNumDof(void) const {
		return  2*nConstraints;
	};
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 24;
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
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		connectedNodes.resize(2);
		connectedNodes[0] = pNode1;
		connectedNodes[1] = pNode2;
	};
	/* ************************************************ */
};


/* TotalReaction - end */


#endif // TOTALEQUATION_H
