/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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


#include <socketstream_out_elem.h>
#include "sock.h"
#include <dataman.h>

#define UNIX_PATH_MAX    108
#define DEFAULT_PORT	5501 /*FIXME:da defineire meglio*/
#define DEFAULT_HOST 	"127.0.0.1"

/* SocketStreamElem - begin */

SocketStreamElem::SocketStreamElem(unsigned int uL, unsigned int nch, ScalarDof *& pn,
		unsigned int oe,
		const char *h, const char *m, unsigned short int p, bool c, int flags)
: Elem(uL, Elem::SOCKETSTREAM_OUTPUT, flag(0)),
NumChannels(nch), pNodes(pn), size(-1), buf(NULL),
OutputEvery(oe), OutputCounter(0), 
host(h), type(AF_INET), sock(0), name(m),
create(c), connected(false), abandoned(false), send_flags(flags)
{
	ASSERT(OutputEvery > 0);

	/* FIXME: size depends on the type of the output signals */
	size = sizeof(doublereal)*nch;
	SAFENEWARR(buf, char, size);
	
   	ASSERT(p > 0);
	data.Port = p;
	if (create) {
		struct sockaddr_in	addr_name;
		int			save_errno;

   		sock = make_inet_socket(&addr_name, host, data.Port, 1,
				&save_errno);
		
		if (sock == -1) {
			const char	*err_msg = strerror(save_errno);

      			silent_cerr("SocketStreamElem(" << name
				<< "): socket() failed "
				"(" << save_errno << ": " << err_msg << ")"
				<< std::endl);
      			throw ErrGeneric();

   		} else if (sock == -2) {
			const char	*err_msg = strerror(save_errno);

      			silent_cerr("SocketStreamElem(" << name
				<< "): bind() failed "
				"(" << save_errno << ": " << err_msg << ")"
				<< std::endl);
      			throw ErrGeneric();

   		} else if (sock == -3) {
      			silent_cerr("SocketStreamElem(" << name
				<< "): illegal host name "
				"(" << save_errno << ")"
				<< std::endl);
      			throw ErrGeneric();
   		}
		
   		if (listen(sock, 1) < 0) {
			save_errno = errno;
			const char	*err_msg = strerror(save_errno);

      			silent_cerr("SocketStreamElem(" << name
				<< "): listen() failed "
				"(" << save_errno << ": " << err_msg << ")"
				<< std::endl);
      			throw ErrGeneric();
   		}
	}
}

SocketStreamElem::SocketStreamElem(unsigned int uL, unsigned int nch, ScalarDof *& pn,
		unsigned int oe,
		const char *m, const char* const Path, bool c, int flags)
: Elem(uL, Elem::SOCKETSTREAM_OUTPUT, flag(0)),
NumChannels(nch), pNodes(pn), size(-1), buf(NULL),
OutputEvery(oe), OutputCounter(0), 
host(NULL), type(AF_LOCAL), sock(0), name(m),
create(c), connected(false), abandoned(false), send_flags(flags)
{
	ASSERT(OutputEvery > 0);

	/* FIXME: size depends on the type of the output signals */
	size = sizeof(doublereal)*nch;
	SAFENEWARR(buf, char, size);
	
	ASSERT(Path != NULL);

	SAFESTRDUP(data.Path, Path);
	
	if (create) {
		int			save_errno;

		sock = make_named_socket(data.Path, 1, &save_errno);
		
		if (sock == -1) {
			const char	*err_msg = strerror(save_errno);

      			silent_cerr("SocketStreamElem(" << name
				<< "): socket() failed "
				"(" << save_errno << ": " << err_msg << ")"
				<< std::endl);
      			throw ErrGeneric();

   		} else if (sock == -2) {
			const char	*err_msg = strerror(save_errno);

      			silent_cerr("SocketStreamElem(" << name
				<< "): bind() failed "
				"(" << save_errno << ": " << err_msg << ")"
				<< std::endl);
      			throw ErrGeneric();
   		}
		
   		if (listen(sock, 1) < 0) {
			save_errno = errno;
			const char	*err_msg = strerror(save_errno);

      			silent_cerr("SocketStreamElem(" << name
				<< "): listen failed "
				"(" << save_errno << ": " << err_msg << ")"
				<< std::endl);
      			throw ErrGeneric();
   		}
	}

}

SocketStreamElem::~SocketStreamElem(void)
{
	//std::cerr << "MBDyn: SocketStreamElem: sock clsed!" << std::endl;
	switch (type) { //connect
	case AF_LOCAL:
		if(create) {
			shutdown(sock,SHUT_RDWR);
			unlink(data.Path);
		}
		close(sock);
		if (data.Path) {
			unlink(data.Path);
			SAFEDELETEARR(data.Path);
			data.Path = 0;
		}
		break;
	case AF_INET:
		if(create) {
			shutdown(sock,SHUT_RDWR);
		}
		close(sock);
		if (host) {
			SAFEDELETEARR(host);
		}
		break;
	default:
		NO_OP;
		break;
	}
			
}

std::ostream&
SocketStreamElem::Restart(std::ostream& out) const
{   	
	//return out << "0. /* SocketStreamDrive not implemented yet */"
	//	<< std::endl;
	out << "  stream output: " << uLabel 
		<< ", stream name, \"" << name << "\",create , ";
	if(create) {
		out << "yes, ";
	} else	{
		out << "no, ";
	}
	/*if(nonblocking)	{
		out << "non blocking, ";
	}*/
	switch (type) {
	case AF_LOCAL:
		out << "path, " << "\"" << data.Path << "\", ";
		break;
	case AF_INET:
		out << "port, " << data.Port << ", ";
		if(host)
			out << "host, " << "\"" << host << "\", ";
		break;
	}
	out << NumChannels;
	for(int i = 0; i < NumChannels; i++) {
		out << ", " ;
		pNodes[i].RestartScalarDofCaller(out);
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
SocketStreamElem::AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP)
{
	/* output only every OutputEvery steps */
	OutputCounter++;
	if (OutputCounter != OutputEvery) {
		return;
	}
	OutputCounter = 0;

	if (!connected) {
		if (create) {
			int tmp_sock = sock;
		   	socklen_t socklen;

			switch (type) {
			case AF_LOCAL: {
   				struct sockaddr_un client_addr;

				pedantic_cout("accepting connection on local socket \""
						<< name << "\" (" << data.Path << ") ..." 
						<< std::endl);

				socklen = sizeof(struct sockaddr_un);
				client_addr.sun_path[0] = '\0';
				sock = accept(tmp_sock,
						(struct sockaddr *)&client_addr, &socklen);
				if (sock == -1) {
					int		save_errno = errno;
					const char	*err_msg = strerror(save_errno);

               				silent_cerr("SocketStreamElem(" << name << ") "
						"accept(" << tmp_sock << ") failed ("
						<< save_errno << ": " << err_msg << ")"
						<< std::endl);
					throw ErrGeneric();
       				}
				if (client_addr.sun_path[0] == '\0') {
					strncpy(client_addr.sun_path, data.Path,
							sizeof(client_addr.sun_path));
				}
      				silent_cout("SocketStreamElem(" << GetLabel()
	  				<< "): connect from " << client_addr.sun_path
					<< std::endl);
	
				break;
			}

			case AF_INET: {
				struct sockaddr_in client_addr;

				pedantic_cout("accepting connection on inet socket \""
					<< name << "\" (" 
					<< inet_ntoa(client_addr.sin_addr) 
					<< ":" << ntohs(client_addr.sin_port)
					<< ") ..." << std::endl);
						
				socklen = sizeof(struct sockaddr_in);
				sock = accept(tmp_sock,
						(struct sockaddr *)&client_addr, &socklen);
				if (sock == -1) {
					int		save_errno = errno;
					const char	*err_msg = strerror(save_errno);

               				silent_cerr("SocketStreamElem(" << name << ") "
						"accept(" << tmp_sock << ") failed ("
						<< save_errno << ": " << err_msg << ")"
						<< std::endl);
					throw ErrGeneric();
       				}
      				silent_cout("SocketStreamElem(" << GetLabel()
	  				<< "): connect from " << inet_ntoa(client_addr.sin_addr)
	  				<< ":" << ntohs(client_addr.sin_port) << std::endl);
				shutdown(tmp_sock,SHUT_RDWR);
				break;
			}

			default:
				break;
			}
			
		} else {
			//connect
			int flag;
			switch (type) {
			case AF_LOCAL: {
				struct sockaddr_un addr;
					
				sock = socket(PF_LOCAL, SOCK_STREAM, 0);
				if (sock < 0){
					int	save_errno = errno;
					const char	*err_msg = strerror(save_errno);
					
					silent_cerr("SocketStreamElem(" << name << ") "
						"socket(PF_LOCAL, SOCK_STREAM, 0) failed "
						"(" << save_errno << ": " << err_msg << ")"
						<< std::endl);
					throw ErrGeneric();
>>>>>>> 1.14
				}
				addr.sun_family = AF_UNIX;
				strncpy(addr.sun_path, data.Path, UNIX_PATH_MAX);

				pedantic_cout("connecting to local socket \""
					<< name << "\" (" << data.Path << ") ..." 
					<< std::endl);

				/* 1000 number of trials FIXME: make configurable */
				flag = MBDyn_connect(sock, (struct sockaddr *) &addr,
						sizeof(addr), 1000);
				break;
			}
	
			case AF_INET: {
				struct sockaddr_in addr;

				sock = socket(PF_INET, SOCK_STREAM, 0);
				if (sock < 0){
					int	save_errno = errno;
					const char	*err_msg = strerror(save_errno);
					
					silent_cerr("SocketStreamElem(" << name << ") "
						"socket(PF_INET, SOCK_STREAM, 0) failed "
						"(" << save_errno << ": " << err_msg << ")"
						<< std::endl);
					throw ErrGeneric();
				}				
				addr.sin_family = AF_INET;
				addr.sin_port = htons (data.Port);
				
				if (inet_aton(host, &(addr.sin_addr)) == 0) {
					silent_cerr("SocketStreamElem(" << name << ") "
						"unknow host " << host << std::endl);
					throw ErrGeneric();	
				}
				pedantic_cout("connecting to inet socket \""
					<< name << "\" (" 
					<< inet_ntoa(addr.sin_addr) 
					<< ":" << ntohs(addr.sin_port)
					<< ") ..." << std::endl);
						
				flag = MBDyn_connect(sock, (struct sockaddr *) &addr,
						sizeof(addr), 1000);
				break;
			}
				
			default:
				break;
			}
			
			switch(flag) {
			case -1:
				silent_cerr("SocketStreamDrive(" << name << ") "
					"connect() failed " << std::endl);
				throw ErrGeneric();
			case -2:
				silent_cerr("SocketStreamDrive(" << name << ") "
					"connection timeout reached " << std::endl);
				throw ErrGeneric();
			case -3:
				silent_cerr("SocketStreamDrive(" << name << ") "
					"poll() failed " << std::endl);
				throw ErrGeneric();
			case -4:
				silent_cerr("SocketStreamDrive(" << name << ") "
					"set socket option failed " << std::endl);
				throw ErrGeneric();
			default:
				NO_OP;
			}
	
		} /* create */
		struct linger lin;
		lin.l_onoff = 1;
		lin.l_linger = 0;
		
		if (setsockopt(sock, SOL_SOCKET, SO_LINGER, &lin, sizeof(lin))) {
			int	save_errno = errno;
			const char	*err_msg = strerror(save_errno);
						
      			silent_cerr("SocketStreamElem(" << name
				<< "): setsockopt failed "
				"(" << save_errno << ": " << err_msg << ")"
				<< std::endl);
      			throw ErrGeneric();
			
		}
		connected = true;
	} /* sock == NULL */

	/* by now, an abandoned element does not write any more;
	 * should we retry or what? */
	if (abandoned) {
		return;
	}

	char *curbuf = buf;
	for (unsigned int i = 0; i < NumChannels; i++) {
		/* assign value somewhere into mailbox buffer */
		doublereal v = pNodes[i].dGetValue();

		doublereal *dbuf = (doublereal *)curbuf;
		dbuf[0] = v;

		curbuf += sizeof(doublereal);
	}
	
	if (send(sock, (void *)buf, size, send_flags) == -1) {
		silent_cerr("SocketStreamElem(" << name << ") "
			<< "Communication closed by host" << std::endl);
		abandoned = true;
		/* FIXME: stop simulation? */
	}
}


Elem *
ReadSocketStreamElem(DataManager *pDM, MBDynParser& HP, unsigned int uLabel)
{
	bool create = false;
	unsigned short int port = DEFAULT_PORT;
	const char *name = NULL;
	const char *host = NULL;
	const char *path = NULL;

	if (HP.IsKeyWord("name") || HP.IsKeyWord("stream" "name")) {
		const char *m = HP.GetStringWithDelims();
		if (m == NULL) {
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
			silent_cerr("\"create\" must be \"yes\" or \"no\" "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric();
		}
	}
	
	if (HP.IsKeyWord("local") || HP.IsKeyWord("path")) {
		const char *m = HP.GetFileName();
		
		if (m == NULL) {
			silent_cerr("unable to read local path for "
				<< psElemNames[Elem::SOCKETSTREAM_OUTPUT]
				<< "(" << uLabel << ") at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric();
		}
		
		SAFESTRDUP(path, m);	
	}

	if (HP.IsKeyWord("port")) {
		if (path != NULL){
			silent_cerr("cannot specify a port "
					"for a local socket in "
				<< psElemNames[Elem::SOCKETSTREAM_OUTPUT]
				<< "(" << uLabel << ") at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric();		
		}
		int p = HP.GetInt();
		/*Da sistemare da qui*/
		
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
		port = p;
	}

	if (HP.IsKeyWord("host")) {
		if (path != NULL){
			silent_cerr("cannot specify an allowed host "
					"for a local socket in "
				<< psElemNames[Elem::SOCKETSTREAM_OUTPUT]
				<< "(" << uLabel << ") at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric();		
		}

		const char *h;
		
		h = HP.GetStringWithDelims();
		if (h == NULL) {
			silent_cerr("unable to read host for "
				<< psElemNames[Elem::SOCKETSTREAM_OUTPUT]
				<< "(" << uLabel << ") at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric();
		}

		SAFESTRDUP(host, h);

	} else if (!path && !create){
		silent_cerr("host undefined for "
			<< psElemNames[Elem::SOCKETSTREAM_OUTPUT]
			<< "(" << uLabel << ") at line "
			<< HP.GetLineData() << std::endl);
		silent_cerr("using default host: "
			<< DEFAULT_HOST << std::endl);
		SAFESTRDUP(host, DEFAULT_HOST);
	}

	int flags = 0;
	if (HP.IsKeyWord("no" "signal")) {
		flags = MSG_NOSIGNAL;
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

	ScalarDof *pNodes = NULL;
	SAFENEWARR(pNodes, ScalarDof, nch);
	for (int i = 0; i < nch; i++) {
		pNodes[i] = ReadScalarDof(pDM, HP, 1);
	}

   	(void)pDM->fReadOutput(HP, Elem::LOADABLE); 

	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		silent_cerr("semicolon expected at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric();
	}
      
      /* costruzione del nodo */
	Elem *pEl = NULL;

	if (path == NULL){
		SAFENEWWITHCONSTRUCTOR(pEl, SocketStreamElem,
				SocketStreamElem(uLabel, nch, pNodes,
					OutputEvery,
					host, name, port, create, flags));
	} else {
		SAFENEWWITHCONSTRUCTOR(pEl, SocketStreamElem,
				SocketStreamElem(uLabel, nch, pNodes,
					OutputEvery,
					name, path, create, flags));
	}

	return pEl;
}

