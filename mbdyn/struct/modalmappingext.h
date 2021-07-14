/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 2007-2017
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

#ifndef MODALMAPPINGEXT_H
#define MODALMAPPINGEXT_H

#include <vector>
#include <string>

#include "modalext.h"
#include "spmapmh.h"
#include "stlvh.h"

/* ModalMappingExt - begin */

class ModalMappingExt : virtual public Elem, public ExtForce {
protected:
	ExtModalForceBase *pEMF;
	unsigned uFlags;
	bool bOutputAccelerations;
	bool bUseReferenceNodeForces;
	bool bRotateReferenceNodeForces;

	const StructNode *pRefNode;

	// Moore-Penrose Generalized Inverse of MSD nodes to modes mapping
	SpMapMatrixHandler *pH;

	// Mapped nodes data
	struct NodeData {
		const StructNode* pNode;
		Vec3 X0;
		Mat3x3 R0;
		Vec3 F;
		Vec3 M;
	};

	std::vector<NodeData> Nodes;

	// rigid-body forces
	Vec3 F0, M0;
	Vec3 F1, M1;
	Vec3 F2, M2;

	STLVectorHandler f;
	STLVectorHandler p;

	STLVectorHandler x;
	STLVectorHandler xP;
	STLVectorHandler q;
	STLVectorHandler qP;

	bool Prepare(ExtFileHandlerBase *pEFH);
	void Send(ExtFileHandlerBase *pEFH, ExtFileHandlerBase::SendWhen when);
	void Recv(ExtFileHandlerBase *pEFH);

public:
	/* Costruttore */
	ModalMappingExt(unsigned int uL,
		DataManager *pDM,
		const StructNode *pRefNode,
		std::vector<const StructNode *>& n,
		SpMapMatrixHandler *pH,
		bool bOutputAccelerations,
		ExtFileHandlerBase *pEFH,
		ExtModalForceBase *pEMF,
		bool bSendAfterPredict,
		int iCoupling,
		ExtModalForceBase::BitMask bm,
		bool bUseReferenceNodeForces,
		bool bRotateReferenceNodeForces,
		flag fOut);

	virtual ~ModalMappingExt(void);

	/* Tipo di forza */
	virtual Force::Type GetForceType(void) const;

	void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	virtual void Output(OutputHandler& OH) const;

	/* *******PER IL SOLUTORE PARALLELO******** */
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	/* ************************************************ */
};

/* ModalMappingExt - end */

class DataManager;
class MBDynParser;

extern SpMapMatrixHandler *
ReadSparseMappingMatrix(MBDynParser& HP, integer& nRows, integer& nCols);

extern Elem*
ReadModalMappingExtForce(DataManager* pDM,
       MBDynParser& HP,
       unsigned int uLabel);

#endif // MODALMAPPINGEXT_H

