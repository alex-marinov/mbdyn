/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2007-2015
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

#ifndef STREDGE_H
#define STREDGE_H

#include "strext.h"

/* StructExtEDGEForce - begin */

class StructExtEDGEForce : virtual public Elem, public StructExtForce {
protected:
	// temporary, to avoid recomputing too much
	std::vector<Vec3> m_x;
	std::vector<Vec3> m_v;

	virtual void SendToStream(std::ostream& outf, ExtFileHandlerBase::SendWhen when);
	virtual void SendToFileDes(int outfd, ExtFileHandlerBase::SendWhen when);
	virtual void RecvFromStream(std::istream& inf);
	virtual void RecvFromFileDes(int infd);
   
public:
	/* Costruttore */
	StructExtEDGEForce(unsigned int uL,
		DataManager *pDM,
		const StructNode *pRefNode,
		bool bUseReferenceNodeForces,
		bool bRotateReferenceNodeForces,
		std::vector<unsigned>& Labels,
		std::vector<const StructNode *>& Nodes,
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

	virtual ~StructExtEDGEForce(void);
};

/* StructExtEDGEForce - end */

#endif // STREDGE_H

