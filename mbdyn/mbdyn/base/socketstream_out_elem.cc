/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2007
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

/*
 * Michele Attolico <attolico@aero.polimi.it>
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
#include "socketstream_out_elem.h"
#include "sock.h"

#define UNIX_PATH_MAX    108
#define DEFAULT_PORT	5501 /* FIXME:da definire meglio */
#define DEFAULT_HOST 	"127.0.0.1"

/* SocketStreamElem - begin */

SocketStreamElem::SocketStreamElem(unsigned int uL,
		std::vector<ScalarValue *>& pn,
		unsigned int oe,
		DataManager *pDM,
		const char *h, const char *m, unsigned short int p, bool c,
		int flags, bool bSF)
: Elem(uL, flag(0)),
Values(pn), size(-1), buf(0),
OutputEvery(oe), OutputCounter(0), 
pUS(0), name(m), send_flags(flags), bSendFirst(bSF)
{
	ASSERT(OutputEvery > 0);

	/* FIXME: size depends on the type of the output signals */
	size = sizeof(doublereal)*Values.size();
	SAFENEWARR(buf, char, size);
	memset(buf, 0, size);

	SAFENEWWITHCONSTRUCTOR(pUS, UseInetSocket, UseInetSocket(h, p, c));
	if (c) {
		pDM->RegisterSocketUser(pUS);
	} else {
		pUS->Connect();
	}
}

SocketStreamElem::SocketStreamElem(unsigned int uL,
		std::vector<ScalarValue *>& pn,
		unsigned int oe,
		DataManager *pDM,
		const char *m, const char* const p, bool c,
		int flags, bool bSF)
: Elem(uL, flag(0)),
Values(pn), size(-1), buf(0),
OutputEvery(oe), OutputCounter(0), 
pUS(0), name(m), send_flags(flags), bSendFirst(bSF)
{
	ASSERT(OutputEvery > 0);

	/* FIXME: size depends on the type of the output signals */
	size = sizeof(doublereal)*Values.size();
	SAFENEWARR(buf, char, size);
	memset(buf, 0, size);
	
	SAFENEWWITHCONSTRUCTOR(pUS, UseLocalSocket, UseLocalSocket(p, c));
	if (c) {
		pDM->RegisterSocketUser(pUS);
	} else {
		pUS->Connect();
	}
}

SocketStreamElem::~SocketStreamElem(void)
{
	if (pUS != 0) {
		SAFEDELETE(pUS);
	}

	if (buf != 0) {
		SAFEDELETEARR(buf);
	}

	for (std::vector<ScalarValue *>::iterator i = Values.begin(); i != Values.end(); i++) {
		delete *i;
	}
}

std::ostream&
SocketStreamElem::Restart(std::ostream& out) const
{   	
	out << "  stream output: " << uLabel 
		<< ", stream name, \"" << name << "\"";
	pUS->Restart(out);
	if (!bSendFirst) {
		out << ", no send first";
	}
	out << ", " << Values.size();
	for (unsigned i = 0; i < Values.size(); i++) {
		out << ", " ;
		//Values[i].RestartScalarDofCaller(out);
	}
	return out << ";" << std::endl;
}	

Elem::Type
SocketStreamElem::GetElemType(void) const
{
	return Elem::SOCKETSTREAM_OUTPUT;
}

void
SocketStreamElem::WorkSpaceDim(integer* piRows, integer* piCols) const
{
	*piRows = 0;
	*piCols = 0;
}

SubVectorHandler&
SocketStreamElem::AssRes(SubVectorHandler& WorkVec, doublereal dCoef,
		const VectorHandler& X, const VectorHandler& XP)
{
	WorkVec.Resize(0);
	return WorkVec;
}

VariableSubMatrixHandler& 
SocketStreamElem::AssJac(VariableSubMatrixHandler& WorkMat, doublereal dCoef,
		const VectorHandler& X, const VectorHandler& XP)
{
	WorkMat.SetNullMatrix();
	return WorkMat;
}

void
SocketStreamElem::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
	if (bSendFirst) {
#if 1
		AfterConvergence(X, XP);
#else
		if (send(pUS->GetSock(), (void *)buf, size, send_flags) == -1) {
			int save_errno = errno;
			char *msg = strerror(save_errno);
		
			silent_cerr("SocketStreamElem(" << name << "): "
					"send() failed "
					"(" << save_errno << ": " << msg << ")"
					<< std::endl);

			pUS->Abandon();
		}
#endif
	}
}

void
SocketStreamElem::AfterConvergence(const VectorHandler& X, 
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
	for (std::vector<ScalarValue *>::iterator i = Values.begin(); i != Values.end(); i++) {
		/* assign value somewhere into mailbox buffer */
		doublereal v = (*i)->dGetValue();

		doublereal *dbuf = (doublereal *)curbuf;
		dbuf[0] = v;

		curbuf += sizeof(doublereal);
	}
	
	if (send(pUS->GetSock(), (void *)buf, size, send_flags) == -1) {
		int save_errno = errno;
		char *msg = strerror(save_errno);
		
		silent_cerr("SocketStreamElem(" << name << "): send() failed "
				"(" << save_errno << ": " << msg << ")"
				<< std::endl);

		pUS->Abandon();
	}
}

void
SocketStreamElem::AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP, const VectorHandler& XPP)
{
		AfterConvergence(X, XP);
}

Elem *
ReadSocketStreamElem(DataManager *pDM, MBDynParser& HP, unsigned int uLabel)
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
				"for SocketStreamElem(" << uLabel
				<< ") at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric();

		} 
		
		SAFESTRDUP(name, m);

	} else {
		silent_cerr("missing stream name "
			"for SocketStreamElem(" << uLabel
			<< ") at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric();
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
			throw ErrGeneric();
		}
	}
	
	if (HP.IsKeyWord("local") || HP.IsKeyWord("path")) {
		const char *m = HP.GetFileName();
		
		if (m == 0) {
			silent_cerr("unable to read local path for "
				<< psElemNames[Elem::SOCKETSTREAM_OUTPUT]
				<< "(" << uLabel << ") at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric();
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
			throw ErrGeneric();		
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
			throw ErrGeneric();
		}
		/* if #undef'd, don't bother checking;
		 * the OS will do it for us */
#endif /* IPPORT_USERRESERVED */
		port = p;
	}

	if (HP.IsKeyWord("host")) {
		if (path != 0) {
			silent_cerr("cannot specify an allowed host "
					"for a local socket in "
				<< psElemNames[Elem::SOCKETSTREAM_OUTPUT]
				<< "(" << uLabel << ") at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric();		
		}

		const char *h;
		
		h = HP.GetStringWithDelims();
		if (h == 0) {
			silent_cerr("unable to read host for "
				<< psElemNames[Elem::SOCKETSTREAM_OUTPUT]
				<< "(" << uLabel << ") at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric();
		}

		SAFESTRDUP(host, h);

	} else if (!path && !create) {
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
			silent_cerr("SocketStreamElem(" << uLabel << "): "
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
			throw ErrGeneric();
		}
		OutputEvery = (unsigned int)i;
	}

	int nch = HP.GetInt();
	if (nch <= 0) {
		silent_cerr("illegal number of channels for "
			<< psElemNames[Elem::SOCKETSTREAM_OUTPUT]
			<< "(" << uLabel << ") at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric();
	}

	std::vector<ScalarValue *> Values(nch);
	for (int i = 0; i < nch; i++) {
		if (HP.IsKeyWord("drive")) {
			Values[i] = new ScalarDriveValue(ReadDriveData(pDM, HP, false));
		} else {
			if (HP.IsKeyWord("node" "dof")) {
				NO_OP; // skip
			}
			Values[i] = new ScalarDofValue(ReadScalarDof(pDM, HP, 1));
		}
	}

   	(void)pDM->fReadOutput(HP, Elem::LOADABLE); 

	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		silent_cerr("semicolon expected at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric();
	}
      
	/* costruzione del nodo */
	Elem *pEl = 0;

	if (path == 0) {
		SAFENEWWITHCONSTRUCTOR(pEl, SocketStreamElem,
				SocketStreamElem(uLabel, Values,
					OutputEvery,
					pDM,
					host, name, port, create, flags,
					bSendFirst));
	} else {
		SAFENEWWITHCONSTRUCTOR(pEl, SocketStreamElem,
				SocketStreamElem(uLabel, Values,
					OutputEvery,
					pDM,
					name, path, create, flags,
					bSendFirst));
	}

	return pEl;
}

#endif // USE_SOCKET
