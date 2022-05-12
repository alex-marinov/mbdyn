/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

/* Forza */

#ifndef STRFORCE_IMPL_H
#define STRFORCE_IMPL_H

#include "force.h"

/* Force - begin */

/* AbsoluteDispForce - begin */

class AbsoluteDispForce : virtual public Elem, public Force {
protected:
	TplDriveOwner<Vec3> f;
	const StructDispNode* pNode;

#ifdef USE_NETCDF
	size_t Var_F;
#endif // USE_NETCDF

public:
	/* Costruttore non banale */
	AbsoluteDispForce(unsigned int uL, const StructDispNode* pN,
		const TplDriveCaller<Vec3>* pDC, flag fOut);

	~AbsoluteDispForce(void);

	/* Tipo di forza */
	virtual Force::Type GetForceType(void) const {
		return Force::ABSOLUTEDISPFORCE;
	};

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* Inverse Dynamics*/
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ ,
		const VectorHandler& /* XPrimePrimeCurr */ ,
		InverseDynamics::Order iOrder);

	virtual void OutputPrepare(OutputHandler &OH);
	virtual void Output(OutputHandler& OH) const;

	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	/* Contributo al residuo durante l'assemblaggio iniziale */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr);

	/* *******PER IL SOLUTORE PARALLELO******** */
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	   utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	/* ************************************************ */
};

/* AbsoluteDispForce - end */


/* AbsoluteInternalDispForce - begin */

class AbsoluteInternalDispForce : virtual public Elem, public Force {
protected:
	TplDriveOwner<Vec3> f;
	const StructDispNode* pNode1;
	const StructDispNode* pNode2;

#ifdef USE_NETCDF
	size_t Var_F;
#endif // USE_NETCDF

public:
	/* Costruttore non banale */
	AbsoluteInternalDispForce(unsigned int uL,
		const StructDispNode* pN1, const StructDispNode* pN2,
		const TplDriveCaller<Vec3>* pDC, flag fOut);

	~AbsoluteInternalDispForce(void);

	/* Tipo di forza */
	virtual Force::Type GetForceType(void) const {
		return Force::ABSOLUTEINTERNALDISPFORCE;
	};

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* Inverse Dynamics*/
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ ,
		const VectorHandler& /* XPrimePrimeCurr */ ,
		InverseDynamics::Order iOrder);

	virtual void OutputPrepare(OutputHandler &OH);
	virtual void Output(OutputHandler& OH) const;

	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	/* Contributo al residuo durante l'assemblaggio iniziale */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr);

	/* *******PER IL SOLUTORE PARALLELO******** */
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	   utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	/* ************************************************ */
};

/* AbsoluteInternalDispForce - end */


/* StructuralForce - begin */

class StructuralForce : virtual public Elem, public Force {
public:
	enum {
		OUTPUT_REL = ToBeOutput::OUTPUT_PRIVATE
	};

protected:
	TplDriveOwner<Vec3> f;
	const StructNode* pNode;

#ifdef USE_NETCDF
	size_t Var_F;
#endif // USE_NETCDF

public:
	/* Costruttore */
	StructuralForce(unsigned int uL,
		const StructNode* pN,
		const TplDriveCaller<Vec3> *pDC,
		flag fOut);

	virtual ~StructuralForce(void);

	/* *******PER IL SOLUTORE PARALLELO******** */
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	   utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	/* ************************************************ */
};

/* StructuralForce - end */


/* AbsoluteForce - begin */

class AbsoluteForce : virtual public Elem, public StructuralForce {
protected:
	const Vec3 Arm;

#ifdef USE_NETCDF
	size_t Var_A;
#endif // USE_NETCDF

public:
	/* Costruttore non banale */
	AbsoluteForce(unsigned int uL, const StructNode* pN,
		const TplDriveCaller<Vec3>* pDC,
		const Vec3& TmpArm,
		flag fOut);

	~AbsoluteForce(void);

	/* Tipo di forza */
	virtual Force::Type GetForceType(void) const {
		return Force::ABSOLUTEFORCE;
	};

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* Inverse Dynamics*/
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ ,
		const VectorHandler& /* XPrimePrimeCurr */ ,
		InverseDynamics::Order iOrder);

	virtual void OutputPrepare(OutputHandler &OH);
	virtual void Output(OutputHandler& OH) const;

	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	/* Contributo al residuo durante l'assemblaggio iniziale */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr);
};

/* AbsoluteForce - end */


/* FollowerForce - begin */

class FollowerForce : virtual public Elem, public StructuralForce {
protected:
	const Vec3 Arm;

#ifdef USE_NETCDF
	size_t Var_A;
#endif // USE_NETCDF

public:
	/* Costruttore banale */
	FollowerForce(unsigned int uL, const StructNode* pN,
		const TplDriveCaller<Vec3>* pDC,
		const Vec3& TmpArm,
		flag fOut);

	~FollowerForce(void);

	/* Tipo di forza */
	virtual Force::Type GetForceType(void) const {
		return Force::FOLLOWERFORCE;
	};

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* Inverse Dynamics*/
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ ,
		const VectorHandler& /* XPrimePrimeCurr */ ,
		InverseDynamics::Order iOrder);

	virtual void OutputPrepare(OutputHandler &OH);
	virtual void Output(OutputHandler& OH) const;

	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	/* Contributo al residuo durante l'assemblaggio iniziale */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr);
};

/* FollowerForce - end */


/* AbsoluteCouple - begin */

class AbsoluteCouple : virtual public Elem, public StructuralForce {
public:
	/* Costruttore banale */
	AbsoluteCouple(unsigned int uL, const StructNode* pN,
		const TplDriveCaller<Vec3>* pDC,
		flag fOut);

	~AbsoluteCouple(void);

	/* Tipo di forza */
	virtual Force::Type GetForceType(void) const {
		return Force::ABSOLUTECOUPLE;
	};

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* Inverse Dynamics*/
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ ,
		const VectorHandler& /* XPrimePrimeCurr */ ,
		InverseDynamics::Order iOrder);

	virtual void OutputPrepare(OutputHandler &OH);
	virtual void Output(OutputHandler& OH) const;

	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	/* Contributo al residuo durante l'assemblaggio iniziale */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr);
};

/* AbsoluteCouple - end */


/* FollowerCouple - begin */

class FollowerCouple : virtual public Elem, public StructuralForce {
public:
	/* Costruttore banale */
	FollowerCouple(unsigned int uL, const StructNode* pN,
		const TplDriveCaller<Vec3>* pDC,
		flag fOut);

	~FollowerCouple(void);

	/* Tipo di forza */
	virtual Force::Type GetForceType(void) const {
		return Force::FOLLOWERCOUPLE;
	};

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* Inverse Dynamics*/
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ ,
		const VectorHandler& /* XPrimePrimeCurr */ ,
		InverseDynamics::Order iOrder);

	virtual void OutputPrepare(OutputHandler &OH);
	virtual void Output(OutputHandler& OH) const;

	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	/* Contributo al residuo durante l'assemblaggio iniziale */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr);
};

/* FollowerCouple - end */


/* StructuralInternalForce - begin */

class StructuralInternalForce
: virtual public Elem, public Force {
protected:
	TplDriveOwner<Vec3> f;
	const StructNode* pNode1;
	const StructNode* pNode2;

#ifdef USE_NETCDF
	size_t Var_F;
#endif // USE_NETCDF

public:
	/* Costruttore */
	StructuralInternalForce(unsigned int uL,
		const StructNode* pN1, const StructNode* pN2,
		const TplDriveCaller<Vec3>* pDC,
		flag fOut);

	virtual ~StructuralInternalForce(void);

	/* *******PER IL SOLUTORE PARALLELO******** */
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	   utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	/* ************************************************ */
};

/* StructuralInternalForce - end */


/* AbsoluteInternalForce - begin */

class AbsoluteInternalForce
: virtual public Elem, public StructuralInternalForce {
protected:
	const Vec3 Arm1;
	const Vec3 Arm2;

#ifdef USE_NETCDF
	size_t Var_A1;
	size_t Var_A2;
#endif // USE_NETCDF

public:
	/* Costruttore non banale */
	AbsoluteInternalForce(unsigned int uL,
		const StructNode* pN1, const StructNode* pN2,
		const TplDriveCaller<Vec3>* pDC,
		const Vec3& TmpArm1, const Vec3& TmpArm2,
		flag fOut);

	~AbsoluteInternalForce(void);

	/* Tipo di forza */
	virtual Force::Type GetForceType(void) const {
		return Force::ABSOLUTEINTERNALFORCE;
	};

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* Inverse Dynamics*/
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ ,
		const VectorHandler& /* XPrimePrimeCurr */ ,
		InverseDynamics::Order iOrder);

	virtual void OutputPrepare(OutputHandler &OH);
	virtual void Output(OutputHandler& OH) const;

	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	/* Contributo al residuo durante l'assemblaggio iniziale */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr);
};

/* AbsoluteInternalForce - end */


/* FollowerInternalForce - begin */

class FollowerInternalForce
: virtual public Elem, public StructuralInternalForce {
protected:
	const Vec3 Arm1;
	const Vec3 Arm2;

#ifdef USE_NETCDF
	size_t Var_A1;
	size_t Var_A2;
#endif // USE_NETCDF

public:
	/* Costruttore banale */
	FollowerInternalForce(unsigned int uL,
		const StructNode* pN1, const StructNode* pN2,
		const TplDriveCaller<Vec3>* pDC,
		const Vec3& TmpArm1, const Vec3& TmpArm2,
		flag fOut);

	~FollowerInternalForce(void);

	/* Tipo di forza */
	virtual Force::Type GetForceType(void) const {
		return Force::FOLLOWERINTERNALFORCE;
	};

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* Inverse Dynamics*/
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ ,
		const VectorHandler& /* XPrimePrimeCurr */ ,
		InverseDynamics::Order iOrder);

	virtual void OutputPrepare(OutputHandler &OH);
	virtual void Output(OutputHandler& OH) const;

	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	/* Contributo al residuo durante l'assemblaggio iniziale */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr);
};

/* FollowerInternalForce - end */


/* AbsoluteInternalCouple - begin */

class AbsoluteInternalCouple
: virtual public Elem, public StructuralInternalForce {
public:
	/* Costruttore banale */
	AbsoluteInternalCouple(unsigned int uL,
		const StructNode* pN1, const StructNode* pN2,
		const TplDriveCaller<Vec3>* pDC,
		flag fOut);

	~AbsoluteInternalCouple(void);

	/* Tipo di forza */
	virtual Force::Type GetForceType(void) const {
		return Force::ABSOLUTEINTERNALCOUPLE;
	};

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* Inverse Dynamics*/
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ ,
		const VectorHandler& /* XPrimePrimeCurr */ ,
		InverseDynamics::Order iOrder);

	virtual void OutputPrepare(OutputHandler &OH);
	virtual void Output(OutputHandler& OH) const;

	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	/* Contributo al residuo durante l'assemblaggio iniziale */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr);
};

/* AbsoluteInternalCouple - end */


/* FollowerInternalCouple - begin */

class FollowerInternalCouple
: virtual public Elem, public StructuralInternalForce {
public:
	/* Costruttore banale */
	FollowerInternalCouple(unsigned int uL,
		const StructNode* pN1, const StructNode* pN2,
		const TplDriveCaller<Vec3>* pDC,
		flag fOut);

	~FollowerInternalCouple(void);

	/* Tipo di forza */
	virtual Force::Type GetForceType(void) const {
		return Force::FOLLOWERINTERNALCOUPLE;
	};

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* Inverse Dynamics*/
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ ,
		const VectorHandler& /* XPrimePrimeCurr */ ,
		InverseDynamics::Order iOrder);

	virtual void OutputPrepare(OutputHandler &OH);
	virtual void Output(OutputHandler& OH) const;

	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	/* Contributo al residuo durante l'assemblaggio iniziale */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr);
};

/* FollowerInternalCouple - end */

#endif // STRFORCE_H_IMPL

