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

/* Giunti universali */


#ifndef UNIVJ_H
#define UNIVJ_H

#include "joint.h"


/* UniversalHingeJoint - begin */

class UniversalHingeJoint : virtual public Elem, public Joint {
private:
	/* Giunto universale: l'asse 3 del primo nodo e l'asse 2 del secondo nodo
	 * rimangono ortogonali (giunto cardanico)
	 * I vettori F, M esprimono le reazioni vincolari di forza e coppia. */
	const StructNode* pNode1;
	const StructNode* pNode2;
	Vec3 d1;
	Mat3x3 R1h;
	Vec3 d2;
	Mat3x3 R2h;
	Vec3 F;
	doublereal dM;
#ifdef USE_NETCDF
	MBDynNcVar Var_Phi;
#endif // USE_NETCDF

public:
	/* Costruttore non banale */
	UniversalHingeJoint(unsigned int uL, const DofOwner* pDO,
		const StructNode* pN1, const StructNode* pN2,
		const Vec3& dTmp1, const Vec3& dTmp2,
		const Mat3x3& R1hTmp, const Mat3x3& R2hTmp, flag fOut);

	/* Distruttore */
	~UniversalHingeJoint(void);

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	/* Tipo di Joint */
	virtual Joint::Type GetJointType(void) const {
		return Joint::UNIVERSALHINGE;
	};

	virtual unsigned int iGetNumDof(void) const {
		return 4;
	};

	DofOrder::Order GetDofType(unsigned int i) const {
		ASSERT(i >= 0 && i < 4);
		return DofOrder::ALGEBRAIC;
	};

	void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 16;
		*piNumCols = 16;
	};

	VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
	SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	void OutputPrepare(OutputHandler& OH);
	void Output(OutputHandler& OH) const;

	/* funzioni usate nell'assemblaggio iniziale */

	virtual unsigned int iGetInitialNumDof(void) const {
		return 8;
	};
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 32;
		*piNumCols = 32;
	};

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	/* Contributo al residuo durante l'assemblaggio iniziale */
	SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr);

	/* *******PER IL SOLUTORE PARALLELO******** */
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	   utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		connectedNodes.resize(2);
		connectedNodes[0] = pNode1;
		connectedNodes[1] = pNode2;
	};
	/* ************************************************ */

	/* returns the dimension of the component */
	const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const;
};

/* UniversalHingeJoint - end */


/* UniversalRotationJoint - begin */

class UniversalRotationJoint : virtual public Elem, public Joint {
private:
	/* Giunto universale: l'asse 3 del primo nodo e l'asse 2 del secondo nodo
	 * rimangono ortogonali (giunto cardanico)
	 * I vettori F, M esprimono le reazioni vincolari di forza e coppia. */
	const StructNode* pNode1;
	const StructNode* pNode2;
#ifdef USE_NETCDF
	MBDynNcVar Var_Phi;
#endif // USE_NETCDF
	Mat3x3 R1h;
	Mat3x3 R2h;
	doublereal dM;

protected:
	OrientationDescription od;

public:
	/* Costruttore non banale */
	UniversalRotationJoint(unsigned int uL, const DofOwner* pDO,
		const StructNode* pN1, const StructNode* pN2,
		const Mat3x3& R1hTmp, const Mat3x3& R2hTmp,
		const OrientationDescription& od,
		flag fOut);

	/* Distruttore */
	~UniversalRotationJoint(void);

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	/* Tipo di Joint */
	virtual Joint::Type GetJointType(void) const {
		return Joint::UNIVERSALROTATION;
	};

	virtual unsigned int iGetNumDof(void) const {
		return 1;
	};

	DofOrder::Order GetDofType(unsigned int i) const {
		ASSERT(i >= 0 && i < 1);
		return DofOrder::ALGEBRAIC;
	};

	void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 7;
		*piNumCols = 7;
	};

	VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
	SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	void OutputPrepare(OutputHandler &OH);
	void Output(OutputHandler& OH) const;


	/* funzioni usate nell'assemblaggio iniziale */

	virtual unsigned int iGetInitialNumDof(void) const {
		return 2;
	};
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 14;
		*piNumCols = 14;
	};

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	/* Contributo al residuo durante l'assemblaggio iniziale */
	SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr);

	/* *******PER IL SOLUTORE PARALLELO******** */
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	   utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		connectedNodes.resize(2);
		connectedNodes[0] = pNode1;
		connectedNodes[1] = pNode2;
	};
	/* ************************************************ */

	/* returns the dimension of the component */
	const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const;
};

/* UniversalRotationJoint - end */


/* UniversalPinJoint - begin */

/* Incastro con liberta' di rotazione su un asse */

class UniversalPinJoint : virtual public Elem, public Joint {
private:
	const StructNode* pNode;
	Vec3 X0;
	Mat3x3 R0;
	Vec3 d;
	Mat3x3 Rh;
	Vec3 F;
	doublereal dM;
#ifdef USE_NETCDF
	MBDynNcVar Var_Phi;
#endif // USE_NETCDF
public:
	/* Costruttore non banale */
	UniversalPinJoint(unsigned int uL, const DofOwner* pDO,
		const StructNode* pN,
		const Vec3& X0Tmp, const Mat3x3& R0Tmp,
		const Vec3& dTmp, const Mat3x3& RhTmp, flag fOut);

	~UniversalPinJoint(void);

	/* Tipo di Joint */
	virtual Joint::Type GetJointType(void) const {
		return Joint::UNIVERSALPIN;
	};

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	virtual unsigned int iGetNumDof(void) const {
		return 4;
	};

	virtual DofOrder::Order
	GetDofType(unsigned int i) const {
		ASSERT(i >= 0 && i < 4);
		return DofOrder::ALGEBRAIC;
	};

	virtual void
	WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 10;
		*piNumCols = 10;
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

	void OutputPrepare(OutputHandler& OH);
	virtual void Output(OutputHandler& OH) const;

	/* funzioni usate nell'assemblaggio iniziale */

	virtual unsigned int iGetInitialNumDof(void) const {
		return 8;
	};
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 20;
		*piNumCols = 20;
	};

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	/* Contributo al residuo durante l'assemblaggio iniziale */
	SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr);

	/* *******PER IL SOLUTORE PARALLELO******** */
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	   utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		connectedNodes.resize(1);
		connectedNodes[0] = pNode;
	};
	/* ************************************************ */

<<<<<<< HEAD
	/* returns the dimension of the component */
	const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const;
=======
	/* test code for getting dimension of components */
	virtual void GetEquationDimension(int index);
>>>>>>> f4de6e740 (initial test of GetEquationDimension)
};

/* UniversalPinJoint - end */

#endif /* UNIVJ_H */

