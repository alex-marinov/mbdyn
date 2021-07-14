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

#ifndef FORCE_H
#define FORCE_H

#include "elem.h"
#include "drive.h"
#include "strnode.h"

extern const char* psForceNames[];


/* Force - begin */

class Force 
: virtual public Elem, public InitialAssemblyElem {
public:
	/* Tipi di Force */
	enum Type {
		UNKNOWN = -1,
		ABSTRACTFORCE = 0,
		ABSTRACTINTERNALFORCE,

		// in mbdyn/struct/strforce[_impl].{h,cc}
		ABSOLUTEDISPFORCE,
		ABSOLUTEINTERNALDISPFORCE,

		ABSOLUTEFORCE,
		FOLLOWERFORCE,
		ABSOLUTECOUPLE,
		FOLLOWERCOUPLE,
	
		ABSOLUTEINTERNALFORCE,
		FOLLOWERINTERNALFORCE,
		ABSOLUTEINTERNALCOUPLE,
		FOLLOWERINTERNALCOUPLE,
	
		// in mbdyn/struct/totalj.{h,cc}
		TOTALFORCE,
		TOTALINTERNALFORCE,

		// in mbdyn/struct/strext.{h,cc}
		EXTERNALSTRUCTURAL,
	
		// in mbdyn/struct/modalforce.{h,cc}
		MODALFORCE,
	
		// in mbdyn/struct/modalext.{h,cc}
		EXTERNALMODAL,
	
		LASTFORCETYPE
	};

public:
	/* Costruttore banale */
	Force(unsigned int uL, flag fOut)
	: Elem(uL, fOut), 
	InitialAssemblyElem(uL, fOut)
	{ 
		NO_OP; 
	};
      
	virtual ~Force(void) { 
		NO_OP; 
	};

	/* Deformable element */
	virtual bool bIsDeformable() {
		return true;
	};
         
	/* Tipo dell'elemento (usato per debug ecc.) */
	virtual Elem::Type GetElemType(void) const { 
		return Elem::FORCE; 
	};   
   
	/* Tipo di forza */
	virtual Force::Type GetForceType(void) const = 0;
   
	virtual VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal /* dCoef */ ,
		const VectorHandler& /* XCurr */ , 
		const VectorHandler& /* XPrimeCurr */ )
	{
		DEBUGCOUT("Entering Force::AssJac()" << std::endl);
	
		WorkMat.SetNullMatrix();
		return WorkMat;
	};

	virtual std::ostream& Restart(std::ostream& out) const;

	virtual unsigned int iGetInitialNumDof(void) const { 
		return 0;
	};

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler& 
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& /* XCurr */ )
	{
		WorkMat.SetNullMatrix();
		return WorkMat;
	};
};

/* Force - end */


/* AbstractForce - begin */

class AbstractForce : virtual public Elem, public Force, public DriveOwner {
protected:
	const Node* pNode;
   
public:
	/* Costruttore banale */
	AbstractForce(unsigned int uL, const Node* pN, 
		const DriveCaller* pDC, flag fOut);
      
	virtual ~AbstractForce(void);
   
	/* Tipo di forza */
	virtual Force::Type GetForceType(void) const { 
		return Force::ABSTRACTFORCE;
	};

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
		*piNumRows = 1; 
		*piNumCols = 1; 
	};

	/* Contributo al residuo */
	SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr);

	virtual void Output(OutputHandler& OH) const;
   
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
	{
		*piNumRows = 1;
		*piNumCols = 1;
	};

	/* Contributo al residuo nell'assemblaggio iniziale */
	SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);

	/* *******PER IL SOLUTORE PARALLELO******** */        
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		connectedNodes.resize(1);
		connectedNodes[0] = pNode;
	};
	/* ************************************************ */
};

/* AbstractForce - end */


/* AbstractInternalForce - begin */

class AbstractInternalForce : virtual public Elem, public Force, public DriveOwner {
protected:
	const Node* pNode1;
	const Node* pNode2;
   
public:
	/* Costruttore banale */
	AbstractInternalForce(unsigned int uL, const Node* pN1, const Node* pN2, 
		const DriveCaller* pDC, flag fOut);
      
	virtual ~AbstractInternalForce(void);
   
	/* Tipo di forza */
	virtual Force::Type GetForceType(void) const { 
		return Force::ABSTRACTINTERNALFORCE;
	};

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
		*piNumRows = 2; 
		*piNumCols = 2; 
	};

	/* Contributo al residuo */
	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr);

	virtual void Output(OutputHandler& OH) const;
   
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
	{
		*piNumRows = 2;
		*piNumCols = 2;
	};

	/* Contributo al residuo nell'assemblaggio iniziale */
	SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);

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

/* AbstractInternalForce - end */

extern Elem* 
ReadForce(DataManager* pDM, 
	MBDynParser& HP, 
	unsigned int uLabel, 
	bool);

#endif /* FORCE_H */

