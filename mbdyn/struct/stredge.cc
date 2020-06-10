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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "dataman.h"
#include "extedge.h"
#include "stredge.h"
#include "Rot.hh"

#include <fstream>
#include <cerrno>
#include <algorithm>

/* StructExtEDGEForce - begin */

/* Costruttore */
StructExtEDGEForce::StructExtEDGEForce(unsigned int uL,
	DataManager *pDM,
	const StructNode *pRefNode,
	bool bUseReferenceNodeForces,
	bool bRotateReferenceNodeForces,
	std::vector<unsigned>& labels,
	std::vector<const StructNode *>& nodes,
	std::vector<Vec3>& offsets,
	bool bSorted,
	bool bLabels,
	bool bOutputAccelerations,
	unsigned uRot,
	ExtFileHandlerBase *pEFH,
	bool bSendAfterPredict,
	int iCoupling,
	unsigned uOutputFlags,
	flag fOut)
: Elem(uL, fOut), 
StructExtForce(uL, pDM, pRefNode, bUseReferenceNodeForces, bRotateReferenceNodeForces,
	labels, nodes, offsets, bSorted, bLabels, bOutputAccelerations, uRot,
	pEFH, bSendAfterPredict, iCoupling, uOutputFlags, fOut),
m_x(nodes.size()),
m_v(nodes.size())
{
	ASSERT(dynamic_cast<ExtFileHandlerEDGE *>(pEFH) != 0);
	ASSERT(bLabels);
	ASSERT(!bSorted);
}

StructExtEDGEForce::~StructExtEDGEForce(void)
{
	NO_OP;
}

void
StructExtEDGEForce::SendToStream(std::ostream& outf, ExtFileHandlerBase::SendWhen when)
{
	/* header:
ext_model, N, 0, 0, 3
*
grid_idents, IF, 1, 121, 0
1 2 3 4 5 6
...
115 116 117 118 119 120
121
str_coordinates, RF, 3, 121, 0
0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
...
0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
0.000000

0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
...
0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
0.000000

-0.000062 -0.000023 0.000000 0.000000 0.000000 0.000000
...
0.035062 0.036876 0.037783 0.039294 0.040654 0.042014
0.043526
velocity, RF, 3, 121, 0
0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
...
0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
0.000000

0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
...
0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
0.000000

0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
...
0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
0.000000
	 */

	if (pRefNode) {
		const Vec3& xRef = pRefNode->GetXCurr();
		const Mat3x3& RRef = pRefNode->GetRCurr();
		const Vec3& vRef = pRefNode->GetVCurr();
		const Vec3& wRef = pRefNode->GetWCurr();

		/*
			p = x + f
			R = R
			v = xp + w cross f
			w = w
			a = xpp + wp cross f + w cross w cross f
			wp = wp
		 */

		std::vector<PointData>::const_iterator point;
		std::vector<Vec3>::iterator ix;
		std::vector<Vec3>::iterator iv;
		for (point = m_Points.cbegin(), ix = m_x.begin(), iv = m_v.begin(); point != m_Points.cend(); ++point, ++ix, ++iv) {
			Vec3 f(point->pNode->GetRCurr()*point->Offset);
			Vec3 x(point->pNode->GetXCurr() + f);
			Vec3 v(point->pNode->GetVCurr() + point->pNode->GetWCurr().Cross(f));

			*ix = RRef.MulTV(x - xRef);
			*iv = RRef.MulTV(v - vRef - wRef.Cross(x - xRef));
		}

	} else {
		std::vector<PointData>::const_iterator point;
		std::vector<Vec3>::iterator ix;
		std::vector<Vec3>::iterator iv;
		for (point = m_Points.cbegin(), ix = m_x.begin(), iv = m_v.begin(); point != m_Points.cend(); ++point, ++ix, ++iv) {
			Vec3 f(point->pNode->GetRCurr()*point->Offset);

			*ix = point->pNode->GetXCurr() + f;
			*iv = point->pNode->GetVCurr() + point->pNode->GetWCurr().Cross(f);
		}
	}

	outf <<
		"ext_model, N, 0, 0, 3\n"
		"*\n";

	// labels - TODO: check
	outf <<
		"grid_idents, IF, 1, " << m_Points.size() << ", 0\n";

	int cnt = 0;
	for (std::vector<PointData>::const_iterator point = m_Points.cbegin(); point != m_Points.cend(); ++point, ++cnt) {
		if (cnt > 0) {
			if ((cnt%6) == 0) {
				outf << "\n";
			} else {
				outf << " ";
			}
		}
		outf << point->uLabel;
	}
	outf << "\n";

	// position
	outf <<
		"str_coordinates, RF, 3, " << m_Points.size() << ", 0\n";

	for (int c = 1; c <= 3; c++) {
		int cnt = 0;
		for (std::vector<Vec3>::const_iterator i = m_x.cbegin(); i != m_x.cend(); ++i, ++cnt) {
			if (cnt > 0) {
				if ((cnt%6) == 0) {
					outf << "\n";
				} else {
					outf << " ";
				}
			}

			outf << (*i)(c);
		}
		outf << "\n";
	}

	// velocity
	outf <<
		"velocity, RF, 3, " << m_Points.size() << ", 0\n";

	for (int c = 1; c <= 3; c++) {
		int cnt = 0;
		for (std::vector<Vec3>::const_iterator i = m_v.cbegin(); i != m_v.cend(); ++i, ++cnt) {
			if (cnt > 0) {
				if ((cnt%6) == 0) {
					outf << "\n";
				} else {
					outf << " ";
				}
			}

			outf << (*i)(c);
		}
		outf << "\n";
	}

	// TODO: other stuff?
}

void
StructExtEDGEForce::SendToFileDes(int outfd, ExtFileHandlerBase::SendWhen when)
{
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

void
StructExtEDGEForce::RecvFromStream(std::istream& inf)
{
	ASSERT(!bSorted);
	ASSERT(bLabels);

	/* header:
ext_model, N, 0, 0, 2
*
grid_idents, IF, 1, 121, 0
1 2 3 4 5 6
...
115 116 117 118 119 120
121
force, RF, 3, 121, 0
0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
...
0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
0.000000

0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
...
0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
0.000000

-0.000062 -0.000023 0.000000 0.000000 0.000000 0.000000
...
0.035062 0.036876 0.037783 0.039294 0.040654 0.042014
0.043526
	 */

	std::vector<unsigned> labels(m_Points.size());
	typedef std::vector<std::vector<PointData>::iterator> RIndexType;
	RIndexType points(m_Points.size());

	// cycle
	unsigned lineno = 0;
	while (true) {
		char buf[BUFSIZ], *p;
		if (mbedge_goto_eol(inf, buf, sizeof(buf))) {
			if (!inf) {
				break;
			}
		}

		lineno++;

		if (buf[0] == '*') {
			continue;
		}

		/*
ext_model, N, 0, 0, 2
		 */

		if (strncasecmp(buf, "ext_model", STRLENOF("ext_model")) == 0) {
			p = buf + STRLENOF("ext_model");

			size_t buflen = sizeof(buf) - STRLENOF("ext_model");
			p = &buf[0] + STRLENOF("ext_model");

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("StructExtEDGEForce(" << GetLabel() << "): unable to skip separator "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_field(p, buflen, "N");
			if (p == 0) {
				silent_cerr("StructExtEDGEForce(" << GetLabel() << "): unable to skip field \"N\" "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("StructExtEDGEForce(" << GetLabel() << "): unable to skip separator "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_field(p, buflen, "0");
			if (p == 0) {
				silent_cerr("StructExtEDGEForce(" << GetLabel() << "): unable to skip field \"0\" "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("StructExtEDGEForce(" << GetLabel() << "): unable to skip separator "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_field(p, buflen, "0");
			if (p == 0) {
				silent_cerr("StructExtEDGEForce(" << GetLabel() << "): unable to skip field \"0\" "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("StructExtEDGEForce(" << GetLabel() << "): unable to skip separator "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_field(p, buflen, "2");
			if (p == 0) {
				silent_cerr("StructExtEDGEForce(" << GetLabel() << "): unable to skip field \"2\" "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			continue;
		}

		/*
grid_idents, IF, 1, 121, 0
		 */

		if (strncasecmp(buf, "grid_idents", STRLENOF("grid_idents")) == 0) {
			p = buf + STRLENOF("grid_idents");

			size_t buflen = sizeof(buf) - STRLENOF("grid_idents");
			p = &buf[0] + STRLENOF("grid_idents");

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("StructExtEDGEForce(" << GetLabel() << "): unable to skip separator "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_field(p, buflen, "IF");
			if (p == 0) {
				silent_cerr("StructExtEDGEForce(" << GetLabel() << "): unable to skip field \"IF\" "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("StructExtEDGEForce(" << GetLabel() << "): unable to skip separator "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_field(p, buflen, "1");
			if (p == 0) {
				silent_cerr("StructExtEDGEForce(" << GetLabel() << "): unable to skip field \"1\" "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("StructExtEDGEForce(" << GetLabel() << "): unable to skip separator "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			char *next;
			errno = 0;
			long nids = strtol(p, &next, 10);
			int save_errno = errno;
			if (next == p) {
				silent_cerr("StructExtEDGEForce(" << GetLabel() << "): unable to read IDs number field "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);

			} else if (save_errno == ERANGE) {
				silent_cerr("StructExtEDGEForce(" << GetLabel() << "): IDs number field overflows "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (unsigned(nids) != m_Points.size()) {
				silent_cerr("StructExtEDGEForce(" << GetLabel() << "): IDs number " << nids << " does not match "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = next;

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("StructExtEDGEForce(" << GetLabel() << "): unable to skip separator "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_field(p, buflen, "0");
			if (p == 0) {
				silent_cerr("StructExtEDGEForce(" << GetLabel() << "): unable to skip field \"0\" "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			for (std::vector<unsigned>::iterator id = labels.begin(); id != labels.end(); ++id) {
				long idx;

				inf >> idx;

				if (idx < 0) {
					silent_cerr("StructExtEDGEForce(" << GetLabel() << "): invalid ID #" << id - labels.begin() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (std::find(labels.begin(), id, unsigned(idx)) < id) {
					silent_cerr("StructExtEDGEForce(" << GetLabel() << "): "
						"duplicate ID #" << id - labels.begin() << " \"" << idx << "\"" << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				*id = idx;
			}


			mbedge_goto_eol(inf, buf, sizeof(buf));

			for (std::vector<PointData>::iterator point = m_Points.begin(); point != m_Points.end(); ++point) {
				std::vector<unsigned>::iterator id = std::find(labels.begin(), labels.end(), point->uLabel);
				if (id == labels.end()) {
					silent_cerr("StructExtEDGEForce(" << GetLabel() << "): "
						"ID \"" << point->uLabel << "\" missing" << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				points[id - labels.begin()] = point;
			}

			continue;
		}

		/*
force, RF, 3, 121, 0
		 */

		if (strncasecmp(buf, "force", STRLENOF("force")) == 0) {
			p = buf + STRLENOF("force");

			size_t buflen = sizeof(buf) - STRLENOF("force");
			p = &buf[0] + STRLENOF("force");

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("StructExtEDGEForce(" << GetLabel() << "): unable to skip separator "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_field(p, buflen, "RF");
			if (p == 0) {
				silent_cerr("StructExtEDGEForce(" << GetLabel() << "): unable to skip field \"RF\" "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("StructExtEDGEForce(" << GetLabel() << "): unable to skip separator "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_field(p, buflen, "3");
			if (p == 0) {
				silent_cerr("StructExtEDGEForce(" << GetLabel() << "): unable to skip field \"3\" "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("StructExtEDGEForce(" << GetLabel() << "): unable to skip separator "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			char *next;
			errno = 0;
			long nids = strtol(p, &next, 10);
			int save_errno = errno;
			if (next == p) {
				silent_cerr("StructExtEDGEForce(" << GetLabel() << "): unable to read IDs number field "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);

			} else if (save_errno == ERANGE) {
				silent_cerr("StructExtEDGEForce(" << GetLabel() << "): IDs number field overflows "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			if (unsigned(nids) != m_Points.size()) {
				silent_cerr("StructExtEDGEForce(" << GetLabel() << "): IDs number " << nids << " does not match "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = next;

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("StructExtEDGEForce(" << GetLabel() << "): unable to skip separator "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_field(p, buflen, "0");
			if (p == 0) {
				silent_cerr("StructExtEDGEForce(" << GetLabel() << "): unable to skip field \"0\" "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			for (RIndexType::iterator p = points.begin(); p != points.end(); p++) {
				double d;

				inf >> d;

				(*p)->F(1) = d;
			}

			mbedge_goto_eol(inf, buf, sizeof(buf));

			for (RIndexType::const_iterator p = points.cbegin(); p != points.cend(); p++) {
				double d;

				inf >> d;

				(*p)->F(2) = d;
			}

			mbedge_goto_eol(inf, buf, sizeof(buf));

			for (RIndexType::const_iterator p = points.cbegin(); p != points.cend(); p++) {
				double d;

				inf >> d;

				(*p)->F(3) = d;
			}

			mbedge_goto_eol(inf, buf, sizeof(buf));

			continue;
		}

		silent_cerr("StructExtEDGEForce(" << GetLabel() << "): unexpected line=" << lineno << ", "
			"\"" << buf << "\"" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

void
StructExtEDGEForce::RecvFromFileDes(int infd)
{
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

