/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2007-2011
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

#ifndef STREXTMAPPING_H
#define STREXTMAPPING_H

#include <vector>
#include <string>

#include "extforce.h"
#include "spmapmh.h"
#include "stlvh.h"

/* StructMappingExtForce - begin */

class StructMappingExtForce : virtual public Elem, public ExtForce {
protected:
	const StructNode *pRefNode;
	bool bUseReferenceNodeForces;
	bool bRotateReferenceNodeForces;
	Vec3 F0, M0;
	Vec3 F1, M1;
	Vec3 F2, M2;

	// Mapping matrix
	SpMapMatrixHandler *pH;

	struct OffsetData {
		unsigned uLabel;
		Vec3 Offset;
		Vec3 F;
	};

	struct NodeData {
		const StructNode *pNode;
		std::vector<OffsetData> Offsets;
		Vec3 F;
		Vec3 M;
	};

	std::vector<NodeData> Nodes;
	unsigned uPoints;
	unsigned uMappedPoints;

	bool bLabels;
	bool bOutputAccelerations;

	unsigned uRRot;

	std::vector<uint32_t> m_qlabels;

	STLVectorHandler m_x;
	STLVectorHandler m_xP;
	STLVectorHandler m_xPP;
	STLVectorHandler m_q;
	STLVectorHandler m_qP;
	STLVectorHandler m_qPP;

	STLVectorHandler m_f;
	STLVectorHandler m_p;

	bool Prepare(ExtFileHandlerBase *pEFH);
	void Send(ExtFileHandlerBase *pEFH, ExtFileHandlerBase::SendWhen when);
	void Recv(ExtFileHandlerBase *pEFH);
   
	void SendToStream(std::ostream& outf, ExtFileHandlerBase::SendWhen when);
	void SendToFileDes(int outfd, ExtFileHandlerBase::SendWhen when);
	void RecvFromStream(std::istream& inf);
	void RecvFromFileDes(int infd);

public:
	/* Costruttore */
	StructMappingExtForce(unsigned int uL,
		DataManager *pDM,
		const StructNode *pRefNode,
		bool bUseReferenceNodeForces,
		bool bRotateReferenceNodeForces,
		std::vector<const StructNode *>& Nodes,
		std::vector<Vec3>& Offsets,
		std::vector<unsigned>& Labels,
		SpMapMatrixHandler *pH,
		std::vector<uint32_t>& MappedLabels,
		bool bLabels,
		bool bOutputAccelerations,
		unsigned uRRot,
		ExtFileHandlerBase *pEFH,
		bool bSendAfterPredict,
		int iCoupling,
		flag fOut);

	virtual ~StructMappingExtForce(void);

	/* Tipo di forza */
	virtual Force::Type GetForceType(void) const { 
		return Force::EXTERNALSTRUCTURAL; 
	};
 
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

/* StructMappingExtForce - end */

class DataManager;
class MBDynParser;

extern Elem*
ReadStructMappingExtForce(DataManager* pDM, 
       MBDynParser& HP, 
       unsigned int uLabel);

#endif // STREXTMAPPING_H

