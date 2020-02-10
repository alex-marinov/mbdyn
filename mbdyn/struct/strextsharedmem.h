/* $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/mbdyn/struct/stredge.h,v 1.9 2017/01/12 14:46:44 masarati Exp $ */
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

#ifndef STREXTSHAREDMEM_H
#define STREXTSHAREDMEM_H

#include "mbconfig.h"

#ifdef HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP

#include "strext.h"
#include "extsharedmem.h"
#include "sharedmem.h"

class StructExtSharedMemForce : virtual public Elem, public StructExtForce {
protected:

    mbdyn::shared_memory_buffer *buf;
    ExtSharedMemHandler* pESMH;

    bool Prepare(ExtFileHandlerBase *pEFH);
	virtual void SendToSharedMem(mbdyn::shared_memory_buffer *buf, ExtFileHandlerBase::SendWhen when);
	virtual void RecvFromSharedMem(mbdyn::shared_memory_buffer *buf);

	// override these from StructExtForce so we can throw error
	// if they are called
	virtual void SendToStream(std::ostream& outf, ExtFileHandlerBase::SendWhen when);
	virtual void SendToFileDes(int outfd, ExtFileHandlerBase::SendWhen when);
	virtual void RecvFromStream(std::istream& inf);
	virtual void RecvFromFileDes(int infd);

public:
	/* Constructor */
	StructExtSharedMemForce(unsigned int uL,
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

	virtual ~StructExtSharedMemForce(void);
};

#endif // HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP
#endif // STRSTREXTSHAREDMEM_H

