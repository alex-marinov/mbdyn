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
#define DEFAULT_PORT	5500 /*FIXME:da defineire meglio*/
#define SYSTEM_PORT	1000 /*FIXME:da defineire meglio*/
#define DEFAULT_HOST 	"127.0.0.1"

/* SocketStreamElem - begin */

SocketStreamElem::SocketStreamElem(unsigned int uL, unsigned int nch, ScalarDof *& pn,
		const char *h, const char *m, unsigned short int p, bool c)
: Elem(uL, Elem::SOCKETSTREAM_OUTPUT, flag(0)),
NumChannels(nch), pNodes(pn), size(-1), buf(NULL),
host(h), type(AF_INET), sock(0), name(m),
create(c), connection_flag(false)
{
	/* FIXME: size depends on the type of the output signals */
	size = sizeof(doublereal)*nch;
	SAFENEWARR(buf, char, size);
	
   	ASSERT(p > 0);
	data.Port = p;
	if(create){
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
				<< "): listen failed "
				"(" << save_errno << ": " << err_msg << ")"
				<< std::endl);
      			throw ErrGeneric();
   		}
	}
}

SocketStreamElem::SocketStreamElem(unsigned int uL, unsigned int nch, ScalarDof *& pn,
		const char *m, const char* const Path, bool c)
: Elem(uL, Elem::SOCKETSTREAM_OUTPUT, flag(0)),
NumChannels(nch), pNodes(pn), size(-1), buf(NULL),
host(NULL), type(AF_LOCAL), sock(0), name(m),
create(c), connection_flag(false)
{
	/* FIXME: size depends on the type of the output signals */
	size = sizeof(doublereal)*nch;
	SAFENEWARR(buf, char, size);
	
	ASSERT(Path != NULL);

	SAFESTRDUP(data.Path, Path);
	
	if(create){
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
	shutdown(sock,SHUT_WR);		
	switch (type) {
	case AF_LOCAL:
		if (data.Path) {
			SAFEDELETEARR(data.Path);
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
	return out << "# not implemented yet" << std::endl;
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
	if (!connection_flag) {
		if (create) {
			int tmp_sock = sock;
			//accept
#ifdef HAVE_SOCKLEN_T
		   	socklen_t socklen;
#else /* !HAVE_SOCKLEN_T */
   			int socklen;
#endif /* !HAVE_SOCKLEN_T */

			switch (type) {
			case AF_LOCAL:
				{
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
						strncpy(client_addr.sun_path, data.Path, sizeof(client_addr.sun_path));
					}
      					silent_cout("SocketStreamElem(" << GetLabel()
		  				<< "): connect from " << client_addr.sun_path
							<< std::endl);
					unlink(data.Path);
					
					break;
				}

			case AF_INET:
				{
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
				NO_OP;
				break;
			}
			
		} else {
			//connect
			switch (type) {
			case AF_LOCAL:
				{
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
					}
					addr.sun_family = AF_UNIX;
					strncpy(addr.sun_path, data.Path, UNIX_PATH_MAX);

					pedantic_cout("connecting to local socket \""
							<< name << "\" (" << data.Path << ") ..." 
							<< std::endl);

					if (connect(sock,(struct sockaddr *) &addr, sizeof (addr)) < 0){
						int	save_errno = errno;
						const char	*err_msg = strerror(save_errno);
						
						silent_cerr("SocketStreamElem(" << name << ") "
							"connect(" << sock << ", \"" << data.Path << "\""
							", " << sizeof(struct sockaddr_un) << ") "
							"failed (" << save_errno << ": " << err_msg << ")"
							<< std::endl);
						throw ErrGeneric();									
					} // da sistemare in modo da rendere non bloccante il connect					
					break;
				}
		
			case AF_INET:
				{
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
							
					if (connect(sock,(struct sockaddr *) &addr, sizeof (addr)) < 0){
						int	save_errno = errno;
						const char	*err_msg = strerror(save_errno);
						
						silent_cerr("SocketStreamElem(" << name << ") "
							"connect(" << sock << ", \"" << host << ":" << data.Port << "\""
							", " << sizeof(struct sockaddr_un) << ") "
							"failed (" << save_errno << ": " << err_msg << ")"
							<< std::endl);
						throw ErrGeneric();					
					} //da sistemare in modo da rendere non bloccante il connect
					break;
				}
				
			default:
				NO_OP;
				break;
			}
	
		} /* create */
		struct linger lin;
		lin.l_onoff = 1;
		lin.l_linger = 0;
		
		if (setsockopt(sock, SOL_SOCKET, SO_LINGER, &lin, sizeof(lin))){
			int	save_errno = errno;
			const char	*err_msg = strerror(save_errno);
						
      			silent_cerr("SocketStreamElem(" << name
				<< "): setsockopt failed "
				"(" << save_errno << ": " << err_msg << ")"
				<< std::endl);
      			throw ErrGeneric();
			
		}
		connection_flag = true;
	} /* sock == NULL */

	char *curbuf = buf;
	for (unsigned int i = 0; i < NumChannels; i++) {
		/* assign value somewhere into mailbox buffer */
		doublereal v = pNodes[i].dGetValue();

		doublereal *dbuf = (doublereal *)curbuf;
		dbuf[0] = v;

		curbuf += sizeof(doublereal);
	}
	
	if (send(sock, (void *)buf, size, 0) == -1) {
		silent_cerr("SocketStreamElem(" << name << ") "
			<< "Communication closed by host" << std::endl);
		// stop simulation?
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
		const char *m = HP.GetStringWithDelims();
		
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
					host, name, port, create));
	} else {
		SAFENEWWITHCONSTRUCTOR(pEl, SocketStreamElem,
				SocketStreamElem(uLabel, nch, pNodes,
					name, path, create));
	}

	return pEl;
}

