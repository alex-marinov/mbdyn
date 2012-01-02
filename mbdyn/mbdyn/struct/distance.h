/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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

/* Vincoli generali */


#ifndef DISTANCE_H
#define DISTANCE_H

#ifdef MBDYN_X_DISTANCE_JOINT

#include <joint.h>
#include <drive.h>

/* DistanceJoint - begin */

class DistanceJoint : virtual public Elem, public Joint, public DriveOwner {
protected:
	const StructNode	*pNode1;
	const StructNode	*pNode2;
	Vec3			Vec;
	doublereal		dAlpha;
	doublereal		dDistance;

	void Abort(void);
   
	/*
	 * Assembla le due matrici
	 *
	 * A = dF/dx e B = dCoef * dF/dxp
	 */
	virtual void
	AssMat(FullSubMatrixHandler& WorkMatA,
			FullSubMatrixHandler& WorkMatB,
			doublereal dCoef,
			const VectorHandler& XCurr,
			const VectorHandler& XPrimeCurr);

public:
	/* Costruttore non banale */
	DistanceJoint(unsigned int uL, const DofOwner* pDO,
			const StructNode* pN1, const StructNode* pN2,
			const DriveCaller* pDC, flag fOut);
   
	~DistanceJoint(void);
   
	/* Tipo di Joint */
	virtual Joint::Type GetJointType(void) const { 
		return Joint::DISTANCE; 
	};
   
	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	virtual unsigned int iGetNumDof(void) const {
		return 1;
	};

	virtual DofOrder::Order GetDofType(unsigned int i) const {
		ASSERT(i == 0);
		return DofOrder::ALGEBRAIC;
	};
   
	virtual DofOrder::Order GetEqType(unsigned int i) const {
		ASSERT(i == 0);
		return DofOrder::ALGEBRAIC;
	};

	virtual void
	WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 7;
		*piNumCols = 7;
	};
   
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

	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
			doublereal dCoef,
			const VectorHandler& XCurr, 
			const VectorHandler& XPrimeCurr);
   
	virtual void Output(OutputHandler& OH) const;

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
	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
			const VectorHandler& XCurr);
   
	/* Contributo al residuo durante l'assemblaggio iniziale */   
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec,
			const VectorHandler& XCurr);
   
	virtual void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0);

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

	/* Adams output stuff */
	virtual unsigned int iGetNumDummyParts(void) const {
		return 1;
	};
	virtual void
	GetDummyPartPos(unsigned int part, Vec3& x, Mat3x3& R) const;
	virtual void
	GetDummyPartVel(unsigned int part, Vec3& v, Vec3& w) const;
#ifdef USE_ADAMS 
	virtual std::ostream&
	WriteAdamsDummyPartCmd(std::ostream& out, unsigned int part,
			unsigned int firstId) const;
#endif /* USE_ADAMS */
};

/* DistanceJoint - end */


/* DistanceJointWithOffset - begin */

class DistanceJointWithOffset : 
virtual public Elem, public DistanceJoint {
private:
	Vec3 f1;
	Vec3 f2;
 
	/*
	 * Assembla le due matrici
	 *
	 * A = dF/dx e B = dCoef * dF/dxp
	 */
	virtual void
	AssMat(FullSubMatrixHandler& WorkMatA,
			FullSubMatrixHandler& WorkMatB,
			doublereal dCoef,
			const VectorHandler& XCurr,
			const VectorHandler& XPrimeCurr);

public:
	/* Costruttore non banale */
	DistanceJointWithOffset(unsigned int uL, const DofOwner* pDO,
			const StructNode* pN1, const StructNode* pN2,
			const Vec3& f1Tmp, const Vec3& f2Tmp,
			const DriveCaller* pDC, flag fOut);

	~DistanceJointWithOffset(void);

	/* Tipo di Joint */
	virtual Joint::Type GetJointType(void) const { 
		return Joint::DISTANCEWITHOFFSET; 
	};

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	virtual void
	WorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
		*piNumRows = 13;
		*piNumCols = 13; 
	};
 
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

	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
			doublereal dCoef,
			const VectorHandler& XCurr, 
			const VectorHandler& XPrimeCurr);

	/* funzioni usate nell'assemblaggio iniziale */

	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 26;
		*piNumCols = 26;
	};
   
	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
			const VectorHandler& XCurr);

	/* Contributo al residuo durante l'assemblaggio iniziale */   
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec,
			const VectorHandler& XCurr);

	virtual void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0);

	/* Adams output stuff */
	virtual void
	GetDummyPartPos(unsigned int part, Vec3& x, Mat3x3& R) const;
	virtual void
	GetDummyPartVel(unsigned int part, Vec3& v, Vec3& w) const;
#ifdef USE_ADAMS 
	virtual std::ostream&
	WriteAdamsDummyPartCmd(std::ostream& out, unsigned int part,
			unsigned int firstId) const;
#endif /* USE_ADAMS */
};

/* DistanceJointWithOffset - end */

#endif

#endif /* DISTANCE_H */

