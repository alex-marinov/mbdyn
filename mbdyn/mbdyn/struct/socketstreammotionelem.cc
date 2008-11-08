/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2008
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#ifdef USE_SOCKET

#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/un.h>
#include <arpa/inet.h>

#include "dataman.h"
#include "socketstreammotionelem.h"
#include "sock.h"
#include "geomdata.h"

#define UNIX_PATH_MAX    108
#define DEFAULT_PORT	5501 /* FIXME: da definire meglio */
#define DEFAULT_HOST 	"127.0.0.1"

/* SocketStreamMotionElem - begin */

void
SocketStreamMotionElem::Init_int(void)
{
	/* FIXME: size depends on the type of the output signals */
	ASSERT(uFlags != 0);
	size = 0;
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
	SAFENEWARR(buf, char, size);
	memset(buf, 0, size);
}

SocketStreamMotionElem::SocketStreamMotionElem(unsigned int uL,
	unsigned uFlags,
	std::vector<StructNode *>& n,
	unsigned int oe,
	DataManager *pDM,
	const char *h, const char *m, unsigned short int p, bool c,
	int flags, bool bSF)
: Elem(uL, flag(0)),
uFlags(uFlags), nodes(n), size(-1), buf(0),
OutputEvery(oe), OutputCounter(0), 
pUS(0), name(m), send_flags(flags), bSendFirst(bSF)
{
	ASSERT(OutputEvery > 0);

	Init_int();

	SAFENEWWITHCONSTRUCTOR(pUS, UseInetSocket, UseInetSocket(h, p, c));
	if (c) {
		pDM->RegisterSocketUser(pUS);
	} else {
		pUS->Connect();
	}
}

SocketStreamMotionElem::SocketStreamMotionElem(unsigned int uL,
	unsigned uFlags,
	std::vector<StructNode *>& n,
	unsigned int oe,
	DataManager *pDM,
	const char *m, const char* const p, bool c,
	int flags, bool bSF)
: Elem(uL, flag(0)),
uFlags(uFlags), nodes(n), size(-1), buf(0),
OutputEvery(oe), OutputCounter(0), 
pUS(0), name(m), send_flags(flags), bSendFirst(bSF)
{
	ASSERT(OutputEvery > 0);

	Init_int();

	SAFENEWWITHCONSTRUCTOR(pUS, UseLocalSocket, UseLocalSocket(p, c));
	if (c) {
		pDM->RegisterSocketUser(pUS);
	} else {
		pUS->Connect();
	}
}

SocketStreamMotionElem::~SocketStreamMotionElem(void)
{
	if (pUS != 0) {
		SAFEDELETE(pUS);
	}

	if (buf != 0) {
		SAFEDELETEARR(buf);
	}
}

std::ostream&
SocketStreamMotionElem::Restart(std::ostream& out) const
{   	
	out << "  stream motion output: " << uLabel 
		<< ", stream name, \"" << name << "\"";
	pUS->Restart(out);

	out << ", output flags";
	if (uFlags & GeometryData::X) {
		out << ", position";
	}
	if (uFlags & GeometryData::R) {
		out << ", orientation matrix";
	}
	if (uFlags & GeometryData::V) {
		out << ", velocity";
	}
	if (uFlags & GeometryData::W) {
		out << ", angular velocity";
	}

	if (!bSendFirst) {
		out << ", no send first";
	}

	for (unsigned i = 0; i < nodes.size(); i++) {
		out << ", " << nodes[i]->GetLabel();
	}
	return out << ";" << std::endl;
}	

Elem::Type
SocketStreamMotionElem::GetElemType(void) const
{
	return Elem::SOCKETSTREAM_OUTPUT;
}

void
SocketStreamMotionElem::WorkSpaceDim(integer* piRows, integer* piCols) const
{
	*piRows = 0;
	*piCols = 0;
}

SubVectorHandler&
SocketStreamMotionElem::AssRes(SubVectorHandler& WorkVec, doublereal dCoef,
		const VectorHandler& X, const VectorHandler& XP)
{
	WorkVec.Resize(0);
	return WorkVec;
}

VariableSubMatrixHandler& 
SocketStreamMotionElem::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& X,
	const VectorHandler& XP)
{
	WorkMat.SetNullMatrix();
	return WorkMat;
}

void
SocketStreamMotionElem::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	if (bSendFirst) {
#if 1
		((SocketStreamMotionElem*)this)->AfterConvergence(X, XP);
#else
		if (send(pUS->GetSock(), (void *)buf, size, send_flags) == -1) {
			int save_errno = errno;
			char *msg = strerror(save_errno);
		
			silent_cerr("SocketStreamMotionElem(" << name << "): "
				"send() failed "
				"(" << save_errno << ": " << msg << ")"
				<< std::endl);

			pUS->Abandon();
		}
#endif
	}
}

void
SocketStreamMotionElem::AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP)
{
	/* by now, an abandoned element does not write any more;
	 * should we retry or what? */
	if (pUS->Abandoned()) {
		return;
	}

	ASSERT(pUS->Connected());

	/* output only every OutputEvery steps */
	OutputCounter++;
	if (OutputCounter != OutputEvery) {
		return;
	}
	OutputCounter = 0;

	char *curbuf = buf;
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
	
	if (send(pUS->GetSock(), (void *)buf, size, send_flags) == -1) {
		int save_errno = errno;
		char *msg = strerror(save_errno);
		
		silent_cerr("SocketStreamMotionElem(" << name << "): "
			"send() failed (" << save_errno << ": " << msg << ")"
			<< std::endl);

		pUS->Abandon();
	}
}

Elem *
ReadSocketStreamMotionElem(DataManager *pDM,
	MBDynParser& HP,
	unsigned int uLabel)
{
	bool create = false;
	unsigned short int port = DEFAULT_PORT;
	const char *name = 0;
	const char *host = 0;
	const char *path = 0;

	if (HP.IsKeyWord("name") || HP.IsKeyWord("stream" "name")) {
		const char *m = HP.GetStringWithDelims();
		if (m == 0) {
			silent_cerr("unable to read stream name "
				"for SocketStreamMotionElem(" << uLabel
				<< ") at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		} 
		
		SAFESTRDUP(name, m);

	} else {
		silent_cerr("missing stream name "
			"for SocketStreamMotionElem(" << uLabel
			<< ") at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (HP.IsKeyWord("create")) {
		if (HP.IsKeyWord("yes")) {
			create = true;

		} else if (HP.IsKeyWord("no")) {
			create = false;

		} else {
			silent_cerr("\"create\" must be either "
					"\"yes\" or \"no\" "
					"at line " << HP.GetLineData()
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
	
	if (HP.IsKeyWord("local") || HP.IsKeyWord("path")) {
		const char *m = HP.GetFileName();
		
		if (m == 0) {
			silent_cerr("unable to read local path for "
				<< psElemNames[Elem::SOCKETSTREAM_OUTPUT]
				<< "(" << uLabel << ") at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		
		SAFESTRDUP(path, m);	
	}

	if (HP.IsKeyWord("port")) {
		if (path != 0){
			silent_cerr("cannot specify a port "
					"for a local socket in "
				<< psElemNames[Elem::SOCKETSTREAM_OUTPUT]
				<< "(" << uLabel << ") at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);		
		}
		int p = HP.GetInt();
		/* Da sistemare da qui */

#ifdef IPPORT_USERRESERVED
		if (p <= IPPORT_USERRESERVED) {
			silent_cerr(psElemNames[Elem::SOCKETSTREAM_OUTPUT]
				<< "(" << uLabel << "): "
				"cannot listen on port " << port
				<< ": less than IPPORT_USERRESERVED=" 
				<< IPPORT_USERRESERVED
				<< " at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		/* if #undef'd, don't bother checking;
		 * the OS will do it for us */
#endif /* IPPORT_USERRESERVED */
		port = p;
	}

	if (HP.IsKeyWord("host")) {
		if (path != 0){
			silent_cerr("cannot specify an allowed host "
					"for a local socket in "
				<< psElemNames[Elem::SOCKETSTREAM_OUTPUT]
				<< "(" << uLabel << ") at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);		
		}

		const char *h;
		
		h = HP.GetStringWithDelims();
		if (h == 0) {
			silent_cerr("unable to read host for "
				<< psElemNames[Elem::SOCKETSTREAM_OUTPUT]
				<< "(" << uLabel << ") at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		SAFESTRDUP(host, h);

	} else if (!path && !create){
		/* INET sockets (!path) must be created if host is missing */
		silent_cerr("host undefined for "
			<< psElemNames[Elem::SOCKETSTREAM_OUTPUT]
			<< "(" << uLabel << ") at line "
			<< HP.GetLineData() << std::endl);
		silent_cerr("using default host: "
			<< DEFAULT_HOST << std::endl);
		SAFESTRDUP(host, DEFAULT_HOST);
	}

	int flags = 0;
	bool bSendFirst = true;
	while (HP.IsArg()) {
		if (HP.IsKeyWord("no" "signal")) {
			flags |= MSG_NOSIGNAL;

		} else if (HP.IsKeyWord("non" "blocking")) {
#ifdef MSG_DONTWAIT
			flags |= MSG_DONTWAIT;
#else /* !MSG_DONTWAIT */
			silent_cerr("SocketStreamMotionElem(" << uLabel << "): "
				"MSG_DONTWAIT undefined; "
				"your mileage may vary" << std::endl);
#endif /* !MSG_DONTWAIT */

		} else if (HP.IsKeyWord("no" "send" "first")) {
			bSendFirst = false;

		} else {
			break;
		}
	}

	unsigned int OutputEvery = 1;
	if (HP.IsKeyWord("output" "every")) {
		int i = HP.GetInt();
		if (i <= 0) {
			silent_cerr("invalid output every value " << i
					<< " at line " << HP.GetLineData()
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		OutputEvery = (unsigned int)i;
	}

	unsigned uFlags = GeometryData::X;
	if (HP.IsKeyWord("output" "flags")) {
		uFlags = 0;
		while (true) {
			if (HP.IsKeyWord("position")) {
				if (uFlags & GeometryData::X) {
					silent_cerr("SocketStreamMotionElem(" << uLabel << "): "
						"position flag already defined "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				uFlags |= GeometryData::X;

			} else if (HP.IsKeyWord("orientation" "matrix")) {
				if (uFlags & GeometryData::ORIENTATION_MASK) {
					silent_cerr("SocketStreamMotionElem(" << uLabel << "): "
						"orientation flag already defined "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				uFlags |= GeometryData::R;

			} else if (HP.IsKeyWord("orientation" "matrix" "transpose")) {
				if (uFlags & GeometryData::ORIENTATION_MASK) {
					silent_cerr("SocketStreamMotionElem(" << uLabel << "): "
						"orientation flag already defined "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				uFlags |= GeometryData::RT;

			} else if (HP.IsKeyWord("velocity")) {
				if (uFlags & GeometryData::V) {
					silent_cerr("SocketStreamMotionElem(" << uLabel << "): "
						"velocity flag already defined "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				uFlags |= GeometryData::V;

			} else if (HP.IsKeyWord("angular" "velocity")) {
				if (uFlags & GeometryData::W) {
					silent_cerr("SocketStreamMotionElem(" << uLabel << "): "
						"angular velocity flag already defined "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				uFlags |= GeometryData::W;

			} else {
				break;
			}
		}
	}

	std::vector<StructNode *> nodes;
	if (HP.IsKeyWord("all")) {
		/* FIXME: todo */

	} else {
		while (HP.IsArg()) {
			nodes.insert(nodes.end(), (StructNode *)pDM->ReadNode(HP, Node::STRUCTURAL));
		}
	}

	Elem *pEl = 0;

	if (path == 0){
		SAFENEWWITHCONSTRUCTOR(pEl, SocketStreamMotionElem,
				SocketStreamMotionElem(uLabel, uFlags, nodes,
					OutputEvery,
					pDM,
					host, name, port, create, flags,
					bSendFirst));
	} else {
		SAFENEWWITHCONSTRUCTOR(pEl, SocketStreamMotionElem,
				SocketStreamMotionElem(uLabel, uFlags, nodes,
					OutputEvery,
					pDM,
					name, path, create, flags,
					bSendFirst));
	}

	return pEl;
}

#endif // USE_SOCKET
