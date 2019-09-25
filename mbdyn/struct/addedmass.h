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

/* elementi di massa, tipo: Elem::Type BODY */

#ifndef ADDEDMASS_H
#define ADDEDMASS_H

/* include per derivazione della classe */

#include "elem.h"
#include "strnode.h"

/* AddedMass - begin */

class AddedMass :
virtual public Elem, public InitialAssemblyElem {
protected:
	const StructDispNode *pNode;
	Vec3 AddedMass;

	/* static moment */
	Vec3 GetS_int(void) const;

	/* moment of inertia */
	Mat3x3 GetJ_int(void) const;

	/* Scrive il contributo dell'elemento al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	void
	AssVecRBK_int(SubVectorHandler& WorkVec);

	void
	AssMatsRBK_int(
		FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB,
		const doublereal& dCoef);

public:
	/* Costruttore definitivo (da mettere a punto) */
	AddedMass(unsigned int uL, const StructDispNode *pNode,
		Vec3 AddedMassTmp, flag fOut);

	virtual ~AddedMass(void);

	/* massa totale */
	Vec3 GetAM(void) const;

	/* momento statico */
	Vec3 GetS(void) const;

	/* momento d'inerzia */
	Mat3x3 GetJ(void) const;

	/* nodo */
	const StructDispNode *pGetNode(void) const;

	/* Tipo dell'elemento (usato solo per debug ecc.) */
	virtual Elem::Type GetElemType(void) const {
		return Elem::BODY;
	};

	/* Numero gdl durante l'assemblaggio iniziale */
	virtual unsigned int iGetInitialNumDof(void) const {
		return 0;
	};

	/* Accesso ai dati privati */
	virtual unsigned int iGetNumPrivData(void) const;
	virtual unsigned int iGetPrivDataIdx(const char *s) const;
	virtual doublereal dGetPrivData(unsigned int i) const;

	/******** PER IL SOLUTORE PARALLELO *********/
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		connectedNodes.resize(1);
		connectedNodes[0] = pNode;
	};
	/**************************************************/
};

/* AddedMass - end */

/* DynamicAddedMass - begin */

class DynamicAddedMass :
virtual public Elem, public AddedMass {
private:

	Vec3 GetB_int(void) const;

	/* Assembla le due matrici necessarie per il calcolo degli
	 * autovalori e per lo jacobiano */
	void AssMats(FullSubMatrixHandler& WorkMatA,
		FullSubMatrixHandler& WorkMatB,
		doublereal dCoef);

public:
	/* Costruttore definitivo (da mettere a punto) */
	DynamicAddedMass(unsigned int uL, const DynamicStructDispNode* pNode,
		Vec3 AddedMassTmp, flag fOut);

	virtual ~DynamicAddedMass(void);

	void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 3;
		*piNumCols = 3;
	};

	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	void AssMats(VariableSubMatrixHandler& WorkMatA,
		VariableSubMatrixHandler& WorkMatB,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* Dimensione del workspace durante l'assemblaggio iniziale.
	 * Occorre tener conto del numero di dof che l'elemento definisce
	 * in questa fase e dei dof dei nodi che vengono utilizzati.
	 * Sono considerati dof indipendenti la posizione e la velocita'
	 * dei nodi */
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 3;
		*piNumCols = 3;
	};

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	/* Contributo al residuo durante l'assemblaggio iniziale */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);

	/* Usata per inizializzare la quantita' di moto */
	virtual void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0);
};

/* DynamicAddedMass - end */

/* StaticAddedMass - begin */

class StaticAddedMass :
virtual public Elem, public AddedMass {
private:
	/* Assembla le due matrici necessarie per il calcolo degli
	 * autovalori e per lo jacobiano */
	bool AssMats(FullSubMatrixHandler& WorkMatA,
		FullSubMatrixHandler& WorkMatB,
		doublereal dCoef);

public:
	/* Costruttore definitivo (da mettere a punto) */
	StaticAddedMass(unsigned int uL, const StaticStructDispNode* pNode,
		Vec3 AddedMass, flag fOut);

	virtual ~StaticAddedMass(void);

	void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 3;
		*piNumCols = 3;
	};

	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	void AssMats(VariableSubMatrixHandler& WorkMatA,
		VariableSubMatrixHandler& WorkMatB,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* inverse dynamics capable element */
	virtual bool bInverseDynamics(void) const;

	/* Inverse Dynamics: */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ ,
		const VectorHandler& /* XPrimePrimeCurr */ ,
		InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS);

	/* Dimensione del workspace durante l'assemblaggio iniziale.
	 * Occorre tener conto del numero di dof che l'elemento definisce
	 * in questa fase e dei dof dei nodi che vengono utilizzati.
	 * Sono considerati dof indipendenti la posizione e la velocita'
	 * dei nodi */
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 3;
		*piNumCols = 3;
	};

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	/* Contributo al residuo durante l'assemblaggio iniziale */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);

	/* Usata per inizializzare la quantita' di moto */
	virtual void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0);
};

/* StaticAddedMass - end */

/* Body - begin */

class Body :
virtual public Elem, public ElemGravityOwner, public InitialAssemblyElem {
protected:
	const StructNode *pNode;
	Vec3 AddedMass;
	Vec3 Xgc;
	Vec3 S0;
	Mat3x3 J0;

	mutable Vec3 STmp;
	mutable Mat3x3 JTmp;

	/* momento statico */
	Vec3 GetS_int(void) const;

	/* momento d'inerzia */
	Mat3x3 GetJ_int(void) const;

	void
	AssVecRBK_int(SubVectorHandler& WorkVec);

	void
	AssMatsRBK_int(
		FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB,
		const doublereal& dCoef,
		const Vec3& Sc);

public:
	/* Costruttore definitivo (da mettere a punto) */
	Body(unsigned int uL, const StructNode *pNode,
		Vec3 AddedMassTmp, const Vec3& XgcTmp, const Mat3x3& JTmp,
		flag fOut);

	virtual ~Body(void);

	/* Scrive il contributo dell'elemento al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	/* massa totale */
	Vec3 GetAM(void) const;

	/* momento statico */
	Vec3 GetS(void) const;

	/* momento d'inerzia */
	Mat3x3 GetJ(void) const;

	/* nodo */
	const StructNode *pGetNode(void) const;

	/* Tipo dell'elemento (usato solo per debug ecc.) */
	virtual Elem::Type GetElemType(void) const {
		return Elem::BODY;
	};

	/* Numero gdl durante l'assemblaggio iniziale */
	virtual unsigned int iGetInitialNumDof(void) const {
		return 0;
	};

	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);

	/* Accesso ai dati privati */
	virtual unsigned int iGetNumPrivData(void) const;
	virtual unsigned int iGetPrivDataIdx(const char *s) const;
	virtual doublereal dGetPrivData(unsigned int i) const;

	/******** PER IL SOLUTORE PARALLELO *********/
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		connectedNodes.resize(1);
		connectedNodes[0] = pNode;
	};
	/**************************************************/
};

/* Body - end */


/* DynamicBody - begin */

class DynamicBody :
virtual public Elem, public Body {
private:

	Vec3 GetB_int(void) const;
	Vec3 GetG_int(void) const;

	/* Assembla le due matrici necessarie per il calcolo degli
	 * autovalori e per lo jacobiano */
	void AssMats(FullSubMatrixHandler& WorkMatA,
		FullSubMatrixHandler& WorkMatB,
		doublereal dCoef,
		bool bGravity,
		const Vec3& GravityAcceleration);

public:
	/* Costruttore definitivo (da mettere a punto) */
	DynamicBody(unsigned int uL, const DynamicStructNode* pNodeTmp,
		Vec3 AddedMassTmp, const Vec3& XgcTmp, const Mat3x3& JTmp,
		flag fOut);

	virtual ~DynamicBody(void);

	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

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

	/* Dimensione del workspace durante l'assemblaggio iniziale.
	 * Occorre tener conto del numero di dof che l'elemento definisce
	 * in questa fase e dei dof dei nodi che vengono utilizzati.
	 * Sono considerati dof indipendenti la posizione e la velocita'
	 * dei nodi */
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	/* Contributo al residuo durante l'assemblaggio iniziale */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);

	/* Usata per inizializzare la quantita' di moto */
	virtual void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0);
};

/* DynamicBody - end */


/* ModalBody - begin */

class ModalBody :
virtual public Elem, public DynamicBody {
private:
	Vec3 XPP, WP;

	/* Assembla le due matrici necessarie per il calcolo degli
	 * autovalori e per lo jacobiano */
	void AssMats(FullSubMatrixHandler& WorkMatA,
		FullSubMatrixHandler& WorkMatB,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr,
		bool bGravity,
		const Vec3& GravityAcceleration);

public:
	/* Costruttore definitivo (da mettere a punto) */
	ModalBody(unsigned int uL, const ModalNode* pNodeTmp,
		Vec3 AddedMassTmp, const Vec3& XgcTmp, const Mat3x3& JTmp,
		flag fOut);

	virtual ~ModalBody(void);

	void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	void AssMats(VariableSubMatrixHandler& WorkMatA,
		VariableSubMatrixHandler& WorkMatB,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
};

/* ModalBody - end */


/* StaticBody - begin */

class StaticBody :
virtual public Elem, public Body {
private:
	/* Assembla le due matrici necessarie per il calcolo degli
	 * autovalori e per lo jacobiano */
	bool AssMats(FullSubMatrixHandler& WorkMatA,
		FullSubMatrixHandler& WorkMatB,
		doublereal dCoef);

public:
	/* Costruttore definitivo (da mettere a punto) */
	StaticBody(unsigned int uL, const StaticStructNode* pNode,
		Vec3 AddedMass, const Vec3& Xgc, const Mat3x3& J,
		flag fOut);

	virtual ~StaticBody(void);

	void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 6;
		*piNumCols = 6;
	};

	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	void AssMats(VariableSubMatrixHandler& WorkMatA,
		VariableSubMatrixHandler& WorkMatB,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* inverse dynamics capable element */
	virtual bool bInverseDynamics(void) const;

	/* Inverse Dynamics: */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ ,
		const VectorHandler& /* XPrimePrimeCurr */ ,
		InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS);

	/* Dimensione del workspace durante l'assemblaggio iniziale.
	 * Occorre tener conto del numero di dof che l'elemento definisce
	 * in questa fase e dei dof dei nodi che vengono utilizzati.
	 * Sono considerati dof indipendenti la posizione e la velocita'
	 * dei nodi */
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
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

	/* Usata per inizializzare la quantita' di moto */
	virtual void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0);
};

/* StaticBody - end */

class DataManager;
class MBDynParser;

extern Elem* ReadBody(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);

#endif /* ADDEDMASS_H */

