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

#ifndef STREXT_H
#define STREXT_H

#include <vector>
#include <string>

#include "extforce.h"

/* StructExtForce - begin */

class StructExtForce : virtual public Elem, public ExtForce {
	// hack to recycle SendToStream for echo
	friend Elem* ReadStructExtForce(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);

protected:
	StructNode *pRefNode;
	bool bUseReferenceNodeForces;
	bool bRotateReferenceNodeForces;
	Vec3 F0, M0;
	Vec3 F1, M1;
	Vec3 F2, M2;
	unsigned uOutputFlags;

public:
	struct PointData {
		unsigned uLabel;
		StructNode *pNode;
		Vec3 Offset;
		Vec3 F;
		Vec3 M;
	};

protected:
	std::vector<PointData> m_Points;

	bool bLabels;
	bool bSorted;
	std::vector<bool> done;

	unsigned uRot;
	bool bOutputAccelerations;

	// buffer for filedes I/O
	unsigned node_kinematics_size;
	unsigned dynamics_size;
	std::vector<char> iobuf;
	uint32_t *iobuf_labels;
	doublereal *iobuf_x;
	doublereal *iobuf_R;
	doublereal *iobuf_theta;
	doublereal *iobuf_euler_123;
	doublereal *iobuf_xp;
	doublereal *iobuf_omega;
	doublereal *iobuf_xpp;
	doublereal *iobuf_omegap;
	doublereal *iobuf_f;
	doublereal *iobuf_m;

	bool Prepare(ExtFileHandlerBase *pEFH);
	void Send(ExtFileHandlerBase *pEFH, ExtFileHandlerBase::SendWhen when);
	void Recv(ExtFileHandlerBase *pEFH);
   
	virtual void SendToStream(std::ostream& outf, ExtFileHandlerBase::SendWhen when);
	virtual void SendToFileDes(int outfd, ExtFileHandlerBase::SendWhen when);
	virtual void RecvFromStream(std::istream& inf);
	virtual void RecvFromFileDes(int infd);
   
public:
	/* Costruttore */
	StructExtForce(unsigned int uL,
		DataManager *pDM,
		StructNode *pRefNode,
		bool bUseReferenceNodeForces,
		bool bRotateReferenceNodeForces,
		std::vector<unsigned>& Labels,
		std::vector<StructNode *>& Nodes,
		std::vector<Vec3>& Offsets,
		bool bSorted,
		bool bLabels,
		bool bOutputAccelerations,
		unsigned bRot,
		ExtFileHandlerBase *pEFH,
		bool bSendAfterPredict,
		int iCoupling,
		unsigned uOutputFlags,
		flag fOut);

	virtual ~StructExtForce(void);

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

/* StructExtForce - end */

class DataManager;
class MBDynParser;

extern Elem*
ReadStructExtForce(DataManager* pDM, 
       MBDynParser& HP, 
       unsigned int uLabel);

#endif // STREXT_H

