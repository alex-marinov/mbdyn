/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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

#ifdef USE_SOCKET_DRIVES

#include <dataman.h>
#include <sockdrv.h>

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <ctype.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/un.h>
#include <arpa/inet.h>

const size_t USERLEN = 8;
const size_t CREDLEN = 128;
const size_t BUFSIZE = 1024;
const char *MBDynSocketDrivePath = "/tmp/mbdyn.sock";

static int
make_socket(unsigned short int port)
{
   	int sock;
   	struct sockaddr_in name;
   
   	/* Create the socket. */
   	sock = socket (PF_LOCAL, SOCK_STREAM, 0);
   	if (sock < 0) {    
      		return -1;
   	}
   
   	/* Give the socket a name. */
   	name.sin_family = AF_INET;
   	name.sin_port = htons (port);
   	name.sin_addr.s_addr = htonl (INADDR_ANY);
   	if (bind (sock, (struct sockaddr *) &name, sizeof (name)) < 0) {
      		return -2;
   	}
   
   	return sock;
}

static int
make_named_socket(const char *path)
{
   	int sock;
   	struct sockaddr_un name;
	size_t size;
   
   	/* Create the socket. */
   	sock = socket (PF_INET, SOCK_STREAM, 0);
   	if (sock < 0) {    
      		return -1;
   	}
   
   	/* Give the socket a name. */
   	name.sun_family = AF_LOCAL;
   	strncpy(name.sun_path, path, sizeof(name.sun_path));
	size = (offsetof (struct sockaddr_un, sun_path)
			+strlen(name.sun_path)+1);
   	if (bind (sock, (struct sockaddr *) &name, size) < 0) {
      		return -2;
   	}
   
   	return sock;
}


SocketDrive::SocketDrive(unsigned int uL, const DriveHandler* pDH,
			 unsigned short int p,
			 AuthMethod* a,
			 integer nd)
: FileDrive(uL, pDH, "socket", nd), 
type(AF_INET),
auth(a), 
pFlags(NULL)
{
   	ASSERT(p > 0);
   	ASSERT(auth != NULL);
   	ASSERT(nd > 0);
   
   	SAFENEWARR(pFlags, int, nd+1);
   	for (int iCnt = 0; iCnt <= nd; iCnt++) {
      		pFlags[iCnt] = SocketDrive::DEFAULT;
   	}
   
   	/* Create the socket and set it up to accept connections. */
	data.Port = p;
   	sock = make_socket(data.Port);
   	if (sock == -1) {
      		std::cerr << "SocketDrive(" << GetLabel() 
			<< "): socket failed" << std::endl;
      		THROW(ErrGeneric());
   	} else if (sock == -2) {
      		std::cerr << "SocketDrive(" << GetLabel() 
			<< "): bind failed" << std::endl;
      		THROW(ErrGeneric());
   	}
   
   	/* non-blocking */
   	int oldflags = fcntl (sock, F_GETFL, 0);
   	if (oldflags == -1) {
      		/* FIXME: err ... */
   	} 
   	oldflags |= O_NONBLOCK;
   	if (fcntl(sock, F_SETFL, oldflags) == -1) {
      		/* FIXME: err ... */
   	}
   
   	if (listen(sock, 1) < 0) {
      		std::cerr << "SocketDrive(" << GetLabel() 
			<< "): listen failed" << std::endl;
      		THROW(ErrGeneric());
   	}
}

SocketDrive::SocketDrive(unsigned int uL, const DriveHandler* pDH,
		const char *path, integer nd)
: FileDrive(uL, pDH, "socket", nd), 
type(AF_LOCAL),
auth(NULL), 
pFlags(NULL)
{
   	ASSERT(path != NULL);
   	ASSERT(nd > 0);

	SAFENEW(auth, NoAuth);
   
   	SAFENEWARR(pFlags, int, nd+1);
   	for (int iCnt = 0; iCnt <= nd; iCnt++) {
      		pFlags[iCnt] = SocketDrive::DEFAULT;
   	}
   
   	/* Create the socket and set it up to accept connections. */
	SAFESTRDUP(data.Path, path);
   	sock = make_named_socket(data.Path);
   	if (sock == -1) {
      		std::cerr << "SocketDrive(" << GetLabel() 
			<< "): socket failed" << std::endl;
      		THROW(ErrGeneric());
   	} else if (sock == -2) {
      		std::cerr << "SocketDrive(" << GetLabel() 
			<< "): bind failed" << std::endl;
      		THROW(ErrGeneric());
   	}
   
   	/* non-blocking */
   	int oldflags = fcntl (sock, F_GETFL, 0);
   	if (oldflags == -1) {
      		/* FIXME: err ... */
   	} 
   	oldflags |= O_NONBLOCK;
   	if (fcntl(sock, F_SETFL, oldflags) == -1) {
      		/* FIXME: err ... */
   	}
   
   	if (listen(sock, 1) < 0) {
      		std::cerr << "SocketDrive(" << GetLabel() 
			<< "): listen failed" << std::endl;
      		THROW(ErrGeneric());
   	}
}


SocketDrive::~SocketDrive(void) 
{
   	/* some shutdown stuff ... */
   	shutdown(sock, SHUT_RDWR /* 2 */ );

	switch (type) {
	case AF_LOCAL:
		SAFEDELETEARR(data.Path);
		break;
	default:
		NO_OP;
		break;
	}
   
   	if (auth != NULL) {
      		SAFEDELETE(auth);
   	}
}

static char *
get_line(char *buf, size_t bufsize,  FILE *fd)
{   
   	int len;
   
   	if (fgets(buf, bufsize, fd) == NULL) {    
      		return NULL;
  	}
   
   	len = strlen(buf); 
   	if (len > 0 && buf[len-1] == '\n') {
      		buf[len-1] = '\0';
      		if (len > 1 && buf[len-2] == '\r') {
	 		buf[len-2] = '\0';
      		}
   	} else {
      		fprintf(stderr, "buffer overflow\n");
      		return NULL;
   	}
   
   	return buf;
}

/*
 * in/out: 
 *     user : puntatore a buffer di 9 bytes (8+'\0')
 *     cred : puntatore a buffer di 129 bytes (128+'\0')
 * 
 * out:
 *     buf : puntatore a buffer statico che contiene una copia della nuova 
 *           linea nel caso sia stata accidentalmente letta
 */
static int
get_auth_token(FILE *fd, char *user, char *cred, char **nextline)
{
   	char buf[BUFSIZE];
   
   	if (get_line(buf, BUFSIZE, fd) == NULL) {
      		return -1;
   	}
   
   	user[0] = '\0';
   	cred[0] = '\0';
   
   	if (strncasecmp(buf, "user:", 5) != 0) {
      		*nextline = (char *)malloc(sizeof(char)*(strlen(buf)+1));
		if ((*nextline) == NULL) {
	 		return -1;
      		}
      		strcpy(*nextline, buf);
      		return 0;
   	} else {
      		char *p;
      		unsigned int i;
      
      		p = buf+5;
      		while (isspace(*p)) {
	 		p++;
      		}
      		for (i = 0; i < USERLEN; i++) {
	 		if ((user[i] = p[i]) == '\0' || isspace(user[i])) {
	    			break;
	 		}
      		}
      		user[i] = '\0';
   	}
   
   	if (get_line(buf, BUFSIZE, fd) == NULL) {
      		return -1;
   	}
   
   	if (strncasecmp(buf, "password:", 9) != 0) {
      		*nextline = (char *)malloc(sizeof(char)*(strlen(buf)+1));
		if ((*nextline) == NULL) {
	 		return -1;
      		}
      		strcpy(*nextline, buf);
      		return 0;
   	} else {
      		char *p;
      		unsigned int i;
      
      		p = buf+9;
      		while (isspace(*p)) {
	 		p++;
      		}
      		for (i = 0; i < CREDLEN; i++) {
	 		if ((cred[i] = p[i]) == '\0' || isspace(cred[i])) {
	    			break;
	 		}
      		}
      		cred[i] = '\0';
   	}
   
   	return 1;
}

void 
SocketDrive::ServePending(const doublereal& /* t */ )
{
   	int cur_sock;
   	struct sockaddr_in client_name;
#ifdef HAVE_SOCKLEN_T
   	socklen_t socklen;
#else /* !HAVE_SOCKLEN_T */
   	int socklen;
#endif /* !HAVE_SOCKLEN_T */
   	FILE* fd;
   
   	/* prova */
   	for (integer iCnt = 1; iCnt <= iNumDrives; iCnt++) {
      		if (pFlags[iCnt] & SocketDrive::IMPULSIVE) {
	 		pdVal[iCnt] = 0.;
      		}
   	}

   	while (1) {
      		char user[USERLEN+1];
      		char cred[CREDLEN+1];

      		int got_value = 0;
      		char *nextline = NULL;
      		const size_t bufsize = BUFSIZE;
      		char buf[bufsize];
           
      		integer label;
      		doublereal value;
      
      		cur_sock = accept(sock, 
				  (struct sockaddr *)&client_name, 
				  &socklen);
      
      		if (cur_sock == -1) {
	 		if (errno != EWOULDBLOCK) {
	    			std::cerr << "SocketDrive(" << GetLabel() 
					<< "): accept failed" << std::endl;
	 		}
	 		return;
      		}
      
      		silent_cout("SocketDrive(" << GetLabel()
		  	<< "): connect from " << inet_ntoa(client_name.sin_addr)
		  	<< ":" << ntohs(client_name.sin_port) << std::endl);
      
      		fd = fdopen(cur_sock, "r");
      
      		if (get_auth_token(fd, user, cred, &nextline) == -1) {
	 		std::cerr << "SocketDrive(" << GetLabel() 
				<< "): corrupted stream" << std::endl;
	 		fclose(fd);
	 		continue;
      		}
      
      		DEBUGCOUT("got auth token: user=\"" << user 
			<< "\", cred=\"" << cred << "\"" << std::endl);
      
      		if (auth->Auth(user, cred) != AuthMethod::AUTH_OK) {
	 		std::cerr << "SocketDrive(" << GetLabel() 
				<< "): authentication failed" << std::endl;
	 		fclose(fd);
	 		continue;
      		}
      
      		DEBUGCOUT("authenticated" << std::endl);
      
      		/*
		 * la nuova linea puo' essere gia' stata letta
		 * da get_auth_token
		 */
      		if (nextline == NULL) {	 
	 		if (get_line(buf, bufsize, fd) == NULL) {
	    			std::cerr << "SocketDrive(" << GetLabel() 
					<< "): corrupted stream" << std::endl;
	    			fclose(fd);
	    			continue;
	 		}
      		} else {
	 		strncpy(buf, nextline, bufsize);
	 		free(nextline);
      		}
      		nextline = buf;      

      		/* legge la label */
      		if (strncasecmp(nextline, "label:", 6) != 0) {
	 		std::cerr << "SocketDrive(" << GetLabel() 
				<< "): missing label" << std::endl;
	 		fclose(fd);
	 		continue;
      		} else {
	 		char *p = nextline+6;
	 		while (isspace(p[0])) {
	    			p++;
	 		}
	 
	 		if (sscanf(p, "%ld", &label) != 1) {
	    			std::cerr << "SocketDrive(" << GetLabel() 
					<< "): unable to read label" << std::endl;
	    			fclose(fd);
	    			continue;
	 		}
	 
	 		if (label <= 0 || label > iNumDrives) {
	    			std::cerr << "SocketDrive(" << GetLabel() 
					<< "): illegal label " 
					<< label << std::endl;
	    			fclose(fd);
	    			continue;
	 		}
	 
	 		while (1) {
	    			if (get_line(buf, bufsize, fd) == NULL) {
	       				std::cerr << "SocketDrive(" << GetLabel()
						<< "): corrupted stream"
						<< std::endl;
	       				fclose(fd);
	       				break;
	    			}
	    
	    			nextline = buf;
	    
	    			if (nextline[0] == '.') {
	       				fclose(fd);
	       				break;
	    			}
	    
	    			if (strncasecmp(nextline, "value:", 6) == 0) {
	       				char *p = nextline+6;
	       				while (isspace(p[0])) {
		  				p++;
	       				}
	       
	       				if (sscanf(p, "%lf", &value) != 1) {
		  				std::cerr << "SocketDrive("
							<< GetLabel()
							<< "): unable to read"
							" value" << std::endl;
		  				fclose(fd);
		  				break;
	       				}
	       				got_value = 1;
	       
	    			} else if (strncasecmp(nextline, "inc:", 4) == 0) {
	       				char *p = nextline+4;
	       				while (isspace(p[0])) {
		  				p++;
	       				}
	       
	       				if (strncasecmp(p, "yes", 3) == 0) {
		  				pFlags[label] |= SocketDrive::INCREMENTAL;
	       				} else if (strncasecmp(p, "no", 2) == 0) {
		  				pFlags[label] &= !SocketDrive::INCREMENTAL;
	       				} else {
		  				std::cerr << "SocketDrive(" 
							<< GetLabel() 
		    					<< "): \"inc\" line"
							" in \"" << nextline 
		    					<< "\" looks corrupted"
							<< std::endl;
		  				fclose(fd);
		  				break;
	       				}
	       				nextline = NULL;
	       
	    			} else if (strncasecmp(nextline, "imp:", 4) == 0) {
	       				char *p = nextline+4;
	       				while (isspace(p[0])) {
		  				p++;
	       				}
	       
	       				if (strncasecmp(p, "yes", 3) == 0) {
		  				pFlags[label] |= SocketDrive::IMPULSIVE;
	       				} else if (strncasecmp(p, "no", 2) == 0) {
		  				pFlags[label] &= !SocketDrive::IMPULSIVE;
	       				} else {
		  				std::cerr << "SocketDrive(" 
							<< GetLabel() 
		    					<< "): \"imp\" line"
							" in \"" << nextline 
		    					<< "\" looks corrupted"
							<< std::endl;
		  				fclose(fd);
		  				break;
	       				}
	       				nextline = NULL;
	    			}
	 		}
          	 
	 		/* usa i valori */
	 		if (got_value) {
	    			if (pFlags[label] & SocketDrive::INCREMENTAL) {
	       				silent_cout("SocketDrive("
						<< GetLabel() << "): adding "
						<<  value << " to label "
						<< label << std::endl);
	       				pdVal[label] += value;
	    			} else {
	       				silent_cout("SocketDrive("
						<< GetLabel() 
			   			<< "): setting label "
						<< label << " to value "
						<<  value << std::endl);
	       				pdVal[label] = value;
	    			}
	 		}
      		}
   	}
}

FileDrive::Type 
SocketDrive::GetFileDriveType(void) const
{
   	return FileDrive::SOCKET;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
SocketDrive::Restart(std::ostream& out) const
{
   	return out << "SocketDrive not implemented yet" << std::endl;
}
   
#endif /* USE_SOCKET_DRIVES */

