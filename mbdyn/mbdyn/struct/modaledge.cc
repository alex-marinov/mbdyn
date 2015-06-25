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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cerrno>

#include "extedge.h"
#include "modaledge.h"

/* ExtForceEDGE - begin */

ExtForceEDGE::ExtForceEDGE(DataManager *pDM)
: pAP(0)
{
	pAP = pDM->pFindElem<const AirProperties, Elem::AIRPROPERTIES>();
	if (pAP == 0) {
		silent_cerr("ExtForceEDGE::ExtForceEDGE: "
			"warning, unable to link to air properties (no airstream contribution to velocity)" << std::endl);
	}
}

ExtForceEDGE::~ExtForceEDGE(void)
{
	NO_OP;
}

bool
ExtForceEDGE::Prepare(ExtFileHandlerBase *pEFH, unsigned uLabel,
	bool bRigid, unsigned uModes)
{
	// consistency is checked later
	return true;
}

/* ExtForceEDGE - end */


/* ExtRigidForceEDGE - begin */

ExtRigidForceEDGE::ExtRigidForceEDGE(DataManager *pDM)
: ExtForceEDGE(pDM)
{
	NO_OP;
}

/*

* intestazione a scelta
body_forces,R,1,6,0
0.1 0.2 0.3 0.4 0.5 0.6

 * NOTE: force and moment orientation wrt/ the global reference frame,
 * the pole of the moment is the node
*/

unsigned
ExtRigidForceEDGE::Recv(ExtFileHandlerBase *pEFH,
	unsigned uFlags,
	unsigned& uLabel,
	Vec3& f, Vec3& m,
	std::vector<doublereal>& a)
{
	std::istream& fin = *pEFH->GetInStream();
	unsigned uOutFlags = ExtModalForceBase::EMF_NONE;

	ASSERT((uFlags & ExtModalForceBase::EMF_RIGID) != 0);
	ASSERT((uFlags & ~ExtModalForceBase::EMF_ALL) == 0);

	// fake
	uLabel = 0;

	// cycle
	unsigned lineno = 0;
	while (true) {
		char buf[BUFSIZ], *p;
		if (mbedge_goto_eol(fin, buf, sizeof(buf))) {
			if (!fin) {
				break;
			}
		}

		lineno++;

		if (buf[0] == '*') {
			continue;
		}

		if (strncasecmp(buf, "body_forces", STRLENOF("body_forces")) == 0) {
			size_t buflen = sizeof(buf) - STRLENOF("body_forces");
			p = &buf[0] + STRLENOF("body_forces");

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("ExtRigidForceEDGE: unable to skip separator "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_field(p, buflen, "R");
			if (p == 0) {
				silent_cerr("ExtRigidForceEDGE: unable to skip field \"R\" "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("ExtRigidForceEDGE: unable to skip separator "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_field(p, buflen, "1");
			if (p == 0) {
				silent_cerr("ExtRigidForceEDGE: unable to skip field \"1\" "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("ExtRigidForceEDGE: unable to skip separator "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_field(p, buflen, "6");
			if (p == 0) {
				silent_cerr("ExtRigidForceEDGE: unable to skip field \"6\" "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("ExtRigidForceEDGE: unable to skip separator "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_field(p, buflen, "0");
			if (p == 0) {
				silent_cerr("ExtRigidForceEDGE: unable to skip field \"0\" "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("ExtRigidForceEDGE: unable to skip separator "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (p[0] != '\0' && p[0] != '\n') {
				silent_cerr("ExtRigidForceEDGE: no line terminator "
					"at line=" << lineno << ", \"" << p << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (uFlags & ExtModalForceBase::EMF_RIGID_DETECT_MASK) {
				doublereal b[6];
				for (int i = 0; i < 6; i++) {
					fin >> b[i];
				}
				f = Vec3(&b[0]);
				m = Vec3(&b[3]);

				uOutFlags |= (uFlags & ExtModalForceBase::EMF_RIGID_DETECT_MASK);
			}

			// skip to eol
			if (mbedge_goto_eol(fin, buf, sizeof(buf))) {
				silent_cerr("ExtRigidForceEDGE: unable to skip to end of line "
					"at line=" << lineno << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			continue;
		}

		silent_cerr("ExtRigidForceEDGE: unexpected line=" << lineno << ", "
			"\"" << buf << "\"" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return uOutFlags;
}

/*

* intestazione a scelta
body_dynamics,N,0,0,3
* Body linear velocity in body axes
VRELV,R,1,3,0
100.0 0.0 0.0
* Body angular velocity in body axes
VRELM,R,1,3,0
0.0 0.1 0.0
* Body reference frame cosines (listed by columns)
OMGMAN,R,3,3,0
11 21 31
12 22 32
13 23 33

 * NOTE:
 *	- node velocity wrt/ the reference frame of the node,
 *	  including the wind velocity (from "air properties")
 *	- node angular velocity wrt/ the reference frame of the node
 *	- node orientation from global to node reference frame
 *	  (the transpose of node's R)
 */

void
ExtRigidForceEDGE::Send(
	ExtFileHandlerBase *pEFH, unsigned uFlags,
	unsigned uLabel,
	const Vec3& x, const Mat3x3& R, const Vec3& v, const Vec3& w,
	const std::vector<doublereal>& q,
	const std::vector<doublereal>& qP)
{
	std::ostream& fout = *pEFH->GetOutStream();
	ASSERT((uFlags & ExtModalForceBase::EMF_RIGID) != 0);
	ASSERT((uFlags & ~ExtModalForceBase::EMF_ALL) == 0);

	Vec3 vv = v;
	Vec3 v_infty;

	if (pAP && pAP->GetVelocity(x, v_infty)) {
		vv -= v_infty;
	}

	Vec3 vB = R.MulTV(vv);
	Vec3 wB = R.MulTV(w);

	fout << "* MBDyn to EDGE rigid body dynamics\n"
		"body_dynamics,N,0,0,3\n"
		"* Body linear velocity in body axes\n"
		"VRELV,R,1,3,0\n"
		<< vB(1) << " " << vB(2) << " " << vB(3) << "\n"
		"* Body angular velocity in body axes\n"
		"VRELM,R,1,3,0\n"
		<< wB(1) << " " << wB(2) << " " << wB(3) << "\n"
		"* Body reference frame cosines (listed by columns)\n"
		"OMGMAN,R,3,3,0\n"
		// NOTE: this is the transpose of the node's R
		<< R(1, 1) << " " << R(2, 1) << " " << R(3, 1) << "\n"
		<< R(1, 2) << " " << R(2, 2) << " " << R(3, 2) << "\n"
		<< R(1, 3) << " " << R(2, 3) << " " << R(3, 3) << "\n";
}

/* ExtRigidForceEDGE - end */

/* ExtModalForceEDGE - begin */

ExtModalForceEDGE::ExtModalForceEDGE(DataManager *pDM)
: ExtForceEDGE(pDM)
{
	NO_OP;
}

/*

* intestazione a scelta
modal_force_flow,R,5,1,0
0.1 0.2 0.3 0.4 0.5

 */

unsigned
ExtModalForceEDGE::Recv(
	ExtFileHandlerBase *pEFH,
	unsigned uFlags,
	unsigned& uLabel,
	Vec3& f, Vec3& m,
	std::vector<doublereal>& a)
{
	std::istream& fin = *pEFH->GetInStream();
	unsigned uOutFlags = ExtModalForceBase::EMF_NONE;

	ASSERT((uFlags & ExtModalForceBase::EMF_MODAL) != 0);
	ASSERT((uFlags & ~ExtModalForceBase::EMF_ALL) == 0);

	// fake
	uLabel = 0;

	// cycle
	unsigned lineno = 0;
	while (true) {
		char buf[BUFSIZ], *p;
		if (mbedge_goto_eol(fin, buf, sizeof(buf))) {
			if (!fin) {
				break;
			}
		}

		lineno++;

		if (buf[0] == '*') {
			continue;
		}

		if (strncasecmp(buf, "modal_force_flow", STRLENOF("modal_force_flow")) == 0) {
			p = buf + STRLENOF("modal_force_flow");

			size_t buflen = sizeof(buf) - STRLENOF("modal_force_flow");
			p = &buf[0] + STRLENOF("modal_force_flow");

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("ExtModalForceEDGE: unable to skip separator "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_field(p, buflen, "R");
			if (p == 0) {
				silent_cerr("ExtModalForceEDGE: unable to skip field \"R\" "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("ExtModalForceEDGE: unable to skip separator "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			char *next;
			errno = 0;
			long nmodes = strtol(p, &next, 10);
			int save_errno = errno;
			if (next == p) {
				silent_cerr("ExtModalForceEDGE: unable to skip modes number field "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);

			} else if (save_errno == ERANGE) {
				silent_cerr("ExtModalForceEDGE: modes number field overflows "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_sep(next, buflen);
			if (p == 0) {
				silent_cerr("ExtModalForceEDGE: unable to skip separator "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_field(p, buflen, "1");
			if (p == 0) {
				silent_cerr("ExtModalForceEDGE: unable to skip field \"1\" "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("ExtModalForceEDGE: unable to skip separator "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_field(p, buflen, "0");
			if (p == 0) {
				silent_cerr("ExtModalForceEDGE: unable to skip field \"0\" "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("ExtModalForceEDGE: unable to skip separator "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (p[0] != '\0' && p[0] != '\n') {
				silent_cerr("ExtModalForceEDGE: no line terminator "
					"at line=" << lineno << ", \"" << p << "\"" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (std::vector<doublereal>::size_type(nmodes) != a.size()) {
				silent_cerr("ExtModalForceEDGE: mode number mismatch, "
					"expected " << a.size() << " got " << nmodes << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (uFlags & ExtModalForceBase::EMF_MODAL_DETECT_MASK) {
				for (std::vector<doublereal>::iterator i = a.begin();
					i != a.end(); ++i)
				{
					fin >> *i;
				}

				uOutFlags |= (uFlags & ExtModalForceBase::EMF_MODAL_DETECT_MASK);
			}

			// skip to eol
			if (mbedge_goto_eol(fin, buf, sizeof(buf))) {
				silent_cerr("ExtModalForceEDGE: unable to skip to end of line "
					"at line=" << lineno << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			continue;
		}

		silent_cerr("ExtModalForceEDGE: unexpected line=" << lineno << ", "
			"\"" << buf << "\"" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return uOutFlags;
}

/*

* intestazione a scelta
modal_state,N,0,0,2
modal_coordinate,R,1,5,0
1.0 2.0 3.0 4.0 5.0
modal_velocity,R,1,5,0
0.1 0.2 0.3 0.4 0.5

 */

void
ExtModalForceEDGE::Send(
	ExtFileHandlerBase *pEFH, unsigned uFlags,
	unsigned uLabel,
	const Vec3& x, const Mat3x3& R, const Vec3& v, const Vec3& w,
	const std::vector<doublereal>& q,
	const std::vector<doublereal>& qP)
{
	std::ostream& fout = *pEFH->GetOutStream();
	ASSERT((uFlags & ExtModalForceBase::EMF_MODAL) != 0);
	ASSERT((uFlags & ~ExtModalForceBase::EMF_ALL) == 0);

	fout << "* MBDyn to EDGE modal dynamics\n"
		"modal_state,N,0,0,2\n"
		"modal_coordinate,R," << q.size() << ",1,0\n"
		<< q[0];
	for (std::vector<doublereal>::const_iterator i = q.begin() + 1;
		i < q.end(); ++i)
	{
		fout << " " << *i;
	}
	fout << "\n"
		"modal_velocity,R," << qP.size() << ",1,0\n"
		<< qP[0];
	for (std::vector<doublereal>::const_iterator i = qP.begin() + 1;
		i < qP.end(); ++i)
	{
		fout << " " << *i;
	}
	fout << "\n";
}

/* ExtModalForceEDGE - end */

