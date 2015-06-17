/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "dataman.h"
#include "socketstreammotionelem.h"
#include "geomdata.h"

/* StreamContentMotion - begin */

StreamContentMotion::StreamContentMotion(unsigned uFlags,
	std::vector<const StructNode *>& n,
	StreamContent::Modifier *pMod)
: StreamContent(0, pMod), uFlags(uFlags), nodes(n)
{
	/* FIXME: size depends on the type of the output signals */
	ASSERT(uFlags != 0);
	unsigned int size = 0;
	if (uFlags & GeometryData::X) {
		size += sizeof(doublereal)*3;
	}
	if (uFlags & GeometryData::R) {
		size += sizeof(doublereal)*9;
	}
	if (uFlags & GeometryData::V) {
		size += sizeof(doublereal)*3;
	}
	if (uFlags & GeometryData::W) {
		size += sizeof(doublereal)*3;
	}
	size *= nodes.size();

	buf.resize(size);
	memset(&buf[0], 0, size);
	m_pMod->Set(size, &buf[0]);
}

StreamContentMotion::~StreamContentMotion(void)
{
	NO_OP;
}

void
StreamContentMotion::Prepare(void)
{
	char *curbuf = &buf[0];
	for (unsigned int i = 0; i < nodes.size(); i++) {
		/* assign value somewhere into mailbox buffer */
		if (uFlags & GeometryData::X) {
			const Vec3& X = nodes[i]->GetXCurr();

			doublereal *dbuf = (doublereal *)curbuf;
			dbuf[0] = X(1);
			dbuf[1] = X(2);
			dbuf[2] = X(3);

			curbuf += 3*sizeof(doublereal);
		}

		if (uFlags & GeometryData::R) {
			const Mat3x3& R = nodes[i]->GetRCurr();

			doublereal *dbuf = (doublereal *)curbuf;
			dbuf[0] = R(1, 1);
			dbuf[1] = R(1, 2);
			dbuf[2] = R(1, 3);
			dbuf[3] = R(2, 1);
			dbuf[4] = R(2, 2);
			dbuf[5] = R(2, 3);
			dbuf[6] = R(3, 1);
			dbuf[7] = R(3, 2);
			dbuf[8] = R(3, 3);

			curbuf += 9*sizeof(doublereal);
		}

		if (uFlags & GeometryData::RT) {
			const Mat3x3& R = nodes[i]->GetRCurr();

			doublereal *dbuf = (doublereal *)curbuf;
			dbuf[0] = R(1, 1);
			dbuf[1] = R(2, 1);
			dbuf[2] = R(3, 1);
			dbuf[3] = R(1, 2);
			dbuf[4] = R(2, 2);
			dbuf[5] = R(3, 2);
			dbuf[6] = R(1, 3);
			dbuf[7] = R(2, 3);
			dbuf[8] = R(3, 3);

			curbuf += 9*sizeof(doublereal);
		}

		if (uFlags & GeometryData::V) {
			const Vec3& V = nodes[i]->GetVCurr();

			doublereal *dbuf = (doublereal *)curbuf;
			dbuf[0] = V(1);
			dbuf[1] = V(2);
			dbuf[2] = V(3);

			curbuf += 3*sizeof(doublereal);
		}

		if (uFlags & GeometryData::W) {
			const Vec3& W = nodes[i]->GetWCurr();

			doublereal *dbuf = (doublereal *)curbuf;
			dbuf[0] = W(1);
			dbuf[1] = W(2);
			dbuf[2] = W(3);

			curbuf += 3*sizeof(doublereal);
		}
	}

	m_pMod->Modify();

	ASSERT(curbuf == &buf[buf.size()]);
}

unsigned
StreamContentMotion::GetNumChannels(void) const
{
	return buf.size()/sizeof(doublereal);
}

/* StreamContentMotion - end */

