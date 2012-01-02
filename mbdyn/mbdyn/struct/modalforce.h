/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2007-2012
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

#ifndef MODALFORCE_H
#define MODALFORCE_H

#include <vector>
#include <string>
#include "extforce.h"
#include "modal.h"

/* ModalForce - begin */

class ModalForce : virtual public Elem, public Force {
protected:
	const Modal *pModal;
	std::vector<unsigned int> modeList;
	std::vector<DriveCaller *> f;
	const Mat3xN *Mt;
	const Mat3xN *Mr;
	Vec3 F, M;

public:
	/* Costruttore */
	ModalForce(unsigned int uL,
		const Modal *pmodal,
		const std::vector<unsigned int>& modeList,
		std::vector<DriveCaller *>& f,
		const Mat3xN *Mt,
		const Mat3xN *Mr,
		flag fOut);

	virtual ~ModalForce(void);

	/* Tipo di forza */
	virtual Force::Type GetForceType(void) const { 
		return Force::MODALFORCE; 
	};
 
	void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
		*piNumRows = (pModal->pGetModalNode() ? 6 : 0) + modeList.size();
		*piNumCols = 1;
	};

	SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr);     

	virtual void Output(OutputHandler& OH) const;

	virtual void InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 0; 
		*piNumCols = 0; 
	};
   
	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler& 
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr)
	{
		WorkMat.SetNullMatrix();
		return WorkMat;
	};

	/* Contributo al residuo durante l'assemblaggio iniziale */   
	virtual SubVectorHandler& 
	InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr)
	{
		WorkVec.ResizeReset(0);
		return WorkVec;
	};

	/* *******PER IL SOLUTORE PARALLELO******** */        
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		if (pModal->pGetModalNode()) {
			connectedNodes.resize(1);
			connectedNodes[0] = pModal->pGetModalNode();
		} else {
			connectedNodes.resize(0);
		}
	};
	/* ************************************************ */
};

/* ModalForce - end */

class DataManager;
class MBDynParser;

extern Elem* ReadModalForce(DataManager* pDM, 
       MBDynParser& HP, 
       unsigned int uLabel);

#endif /* MODALFORCE_H */

