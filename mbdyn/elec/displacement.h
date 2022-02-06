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

/* elementi elettrici, tipo: Elem::ELECTRIC */

#ifndef DISPLACEMENT_H
#define DISPLACEMENT_H

#include "elec.h"

/* DispMeasure - begin */

class DispMeasure : virtual public Elem, public Electric {
private:
	const StructNode* pStrNode1;
	const StructNode* pStrNode2;
	mutable const ScalarDifferentialNode* pAbsNode;
	Vec3 f1;
	Vec3 f2;
   
public:
	DispMeasure(unsigned int uL, const DofOwner* pD,
		const StructNode* pS1, const StructNode* pS2, 
		const ScalarDifferentialNode* pA,
		const Vec3& Tmpf1, const Vec3& Tmpf2,
		flag fOut);
	~DispMeasure(void);

	virtual Electric::Type GetElectricType(void) const {
		return Electric::DISPLACEMENT;
	};
   
	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;
   
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
      
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
   
	/* Setta i valori iniziali delle variabili (e fa altre cose)
	 * prima di iniziare l'integrazione */
	virtual void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0);
   
	/* *******PER IL SOLUTORE PARALLELO******** */        
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		connectedNodes.resize(3);
		connectedNodes[0] = pStrNode1;
		connectedNodes[1] = pStrNode2;
		connectedNodes[2] = pAbsNode;
	};
	/* ************************************************ */

	/* returns the dimension of the component */
	const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const;

	/* describes the dimension of components of equation */
    virtual std::ostream& DescribeEq(std::ostream& out,
		  const char *prefix = "",
		  bool bInitial = false) const;
};

/* DispMeasure - end */

#endif /* DISPLACEMENT_H */

